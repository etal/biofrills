"""Compare the library with the outputs recorded under Python 2.7.

tests/data/README.md describes the fixtures under tests/data/py27/library and
documents every reason a current value may differ from a recorded one.
Values are compared exactly, with integers and floats kept distinct, except
where that file documents a floating-point difference; those comparisons use
the relative tolerances defined below.  A consensus case listed in
tests/data/differences.json is expected to give the outcome listed there
instead of the recorded one.  An exception is compared by its type alone,
since the message text is the interpreter's.
"""

import json
import math
from functools import cache

import pytest
from alignments import DATA, rows
from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment

from biofrills import alnutils, consensus, pairutils
from biofrills.stats import chisq

LIBRARY = DATA / "py27" / "library"
ALIGNMENT_FIXTURES = sorted(path.name for path in LIBRARY.glob("*.aln.json"))

# The summation order of the aa_frequencies normalizer and the compensated
# built-in sum() of Python 3.12 in sequence_weights move values by at most
# 3 units in the last place, a relative 4.2e-16.
SUM_REL_TOL = 1e-15
# math.lgamma depends on the platform's C math library; on arm64 macOS,
# chisqprob differs from the x86_64 capture by a relative 2.9e-14 at most.
LGAMMA_REL_TOL = 1e-13

# sequence_weights scalings whose normalizer is a sum of the raw weights.
SUMMED_SCALINGS = {"sum1", "avg1", "andy"}
# Weights names whose weighted aa_frequencies are normalized by a sum of
# non-integer counts.
FRACTIONAL_WEIGHTS = {True, "seqweights"}


@cache
def load(fname):
    with open(LIBRARY / fname, encoding="ascii") as handle:
        return json.load(handle)


@cache
def alignment(fname):
    return AlignIO.read(DATA / load(fname)["input"], "fasta")


@cache
def difference_entries():
    with open(DATA / "differences.json", encoding="ascii") as handle:
        return json.load(handle)["differences"]


@cache
def differences():
    """The expected current outcome of each case in differences.json."""
    return {
        (entry["file"], entry.get("alignment"), entry["case"]): entry["current"]
        for entry in difference_entries()
    }


def build_weights(aln, spec):
    """The sequence weights that a recorded weights name stands for."""
    if spec == "seqweights":
        return alnutils.sequence_weights(aln, "none")
    if spec == "unit":
        return [1] * len(aln)
    return spec


def as_records(records):
    return [[rec.id, rec.description, str(rec.seq)] for rec in records]


def outcome(func, *args, **kwargs):
    """Call func; return its result as JSON would record it, or its error."""
    try:
        value = func(*args, **kwargs)
    except Exception as exc:
        return {"error": {"type": type(exc).__name__}}
    # A round trip through JSON turns tuples into lists and dict subclasses
    # into dicts, and reads floats back exactly.
    return {"value": json.loads(json.dumps(value))}


def same(actual, expected, rel_tol):
    if type(actual) is not type(expected):
        return False
    if isinstance(expected, float) and rel_tol:
        return math.isclose(actual, expected, rel_tol=rel_tol, abs_tol=0.0)
    if isinstance(expected, dict):
        return actual.keys() == expected.keys() and all(
            same(actual[key], expected[key], rel_tol) for key in expected
        )
    if isinstance(expected, list):
        return len(actual) == len(expected) and all(
            same(a, e, rel_tol) for a, e in zip(actual, expected, strict=True)
        )
    return actual == expected


def comparable(recorded):
    """A recorded outcome without the message text of an exception."""
    if "error" in recorded:
        return {"error": {"type": recorded["error"]["type"]}}
    return recorded


def matches(actual, recorded, rel_tol=0.0):
    """Whether a current outcome agrees with a recorded one."""
    if "error" in recorded:
        return actual == comparable(recorded)
    return "value" in actual and same(actual["value"], recorded["value"], rel_tol)


def check(found):
    """Fail with the first few mismatches, if any."""
    assert not found, f"{len(found)} mismatches, first: {found[:3]}"


@pytest.mark.parametrize("fname", ALIGNMENT_FIXTURES)
def test_sequence_weights(fname):
    aln = alignment(fname)
    found = []
    for case in load(fname)["sequence_weights"]:
        rel_tol = SUM_REL_TOL if case["kwargs"]["scaling"] in SUMMED_SCALINGS else 0
        actual = outcome(alnutils.sequence_weights, aln, **case["kwargs"])
        if not matches(actual, case["result"], rel_tol):
            found.append((case["kwargs"], actual))
    check(found)


@pytest.mark.parametrize("fname", ALIGNMENT_FIXTURES)
@pytest.mark.parametrize("section", ["aa_counts", "aa_frequencies"])
def test_residue_counts(fname, section):
    aln = alignment(fname)
    func = getattr(alnutils, section)
    found = []
    for case in load(fname)[section]:
        fractional = section == "aa_frequencies" and (
            case["weights"] in FRACTIONAL_WEIGHTS
        )

        def call(spec=case["weights"], kwargs=case["kwargs"]):
            return func(aln, weights=build_weights(aln, spec), **kwargs)

        actual = outcome(call)
        if not matches(actual, case["result"], SUM_REL_TOL if fractional else 0):
            found.append((case["weights"], case["kwargs"], actual))
    check(found)


@pytest.mark.parametrize("fname", ALIGNMENT_FIXTURES)
def test_consensus(fname):
    aln = alignment(fname)
    expected_changes = differences()
    found = []
    for index, case in enumerate(load(fname)["consensus"]):

        def call(spec=case["weights"], kwargs=case["kwargs"]):
            return consensus.consensus(aln, weights=build_weights(aln, spec), **kwargs)

        actual = outcome(call)
        expected = expected_changes.get(("library/" + fname, None, index))
        if expected is None:
            expected = case["result"]
        if not matches(actual, expected):
            found.append((index, case["weights"], case["kwargs"], actual))
    check(found)


@pytest.mark.parametrize("fname", ALIGNMENT_FIXTURES)
def test_blocks(fname):
    aln = alignment(fname)
    found = []
    for case in load(fname)["blocks"]:
        actual = outcome(
            lambda case=case: as_records(
                alnutils.blocks(aln, weights=case["weights"], **case["kwargs"])
            )
        )
        if not matches(actual, case["result"]):
            found.append((case["weights"], case["kwargs"], actual))
    check(found)


@pytest.mark.parametrize("fname", ALIGNMENT_FIXTURES)
def test_remove_empty_cols(fname):
    aln = alignment(fname)
    actual = outcome(lambda: as_records(alnutils.remove_empty_cols(aln)))
    assert matches(actual, load(fname)["remove_empty_cols"])


@pytest.mark.parametrize("fname", ALIGNMENT_FIXTURES)
@pytest.mark.parametrize("section", ["score_pairwise", "identity_abs", "identity_avg"])
def test_pairwise(fname, section):
    seqs = rows(alignment(fname))
    func = getattr(pairutils, section)
    found = []
    for i, j, recorded in load(fname)[section]:
        actual = outcome(func, seqs[i], seqs[j])
        if not matches(actual, recorded):
            found.append((i, j, actual))
    check(found)


def tie_batch_alignments():
    data = load("tie-batch.json")
    groups = {}
    with open(DATA / data["input"], encoding="ascii") as handle:
        for rec in SeqIO.parse(handle, "fasta"):
            groups.setdefault(rec.id.split("_")[0], []).append(rec)
    for case in data["alignments"]:
        yield case, MultipleSeqAlignment(groups[case["id"]])


def test_tie_batch_sequence_weights():
    kwargs = load("tie-batch.json")["sequence_weights_kwargs"]
    found = []
    for case, aln in tie_batch_alignments():
        actual = outcome(alnutils.sequence_weights, aln, **kwargs)
        if not matches(actual, case["sequence_weights"], SUM_REL_TOL):
            found.append((case["id"], actual))
    check(found)


def test_tie_batch_consensus():
    all_kwargs = load("tie-batch.json")["consensus_kwargs"]
    expected_changes = differences()
    found = []
    for case, aln in tie_batch_alignments():
        for index, kwargs in enumerate(all_kwargs):
            actual = outcome(consensus.consensus, aln, **kwargs)
            expected = expected_changes.get(
                ("library/tie-batch.json", case["id"], index)
            )
            if expected is None:
                expected = case["consensus"][index]
            if not matches(actual, expected):
                found.append((case["id"], kwargs, actual))
    check(found)


def test_chisqprob():
    found = []
    for x, df, recorded in load("chisqprob.json")["rows"]:
        actual = outcome(chisq.chisqprob, x, df)
        if not matches(actual, recorded, LGAMMA_REL_TOL):
            found.append((x, df, actual))
    check(found)


def test_differences_name_changed_cases():
    """Each entry of differences.json names a recorded case that changed."""
    found = []
    for entry in difference_entries():
        data = load(entry["file"].removeprefix("library/"))
        if "alignment" in entry:
            (case,) = [c for c in data["alignments"] if c["id"] == entry["alignment"]]
            recorded = case["consensus"][entry["case"]]
        else:
            recorded = data["consensus"][entry["case"]]["result"]
        if entry["recorded"] != comparable(recorded) or matches(
            entry["current"], recorded
        ):
            found.append(entry)
    check(found)
