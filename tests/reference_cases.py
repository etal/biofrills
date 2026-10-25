"""The library cases recorded under Python 2.7, as calls on the current code.

tests/data/README.md describes the fixtures under tests/data/py27/library and
documents every reason a current value may differ from a recorded one.
tests/test_reference.py and tests/data/replay.py both enumerate the cases
here, so that the test and the differences.json that replay.py writes compare
outcomes in the same way: values exactly, with integers and floats kept
distinct, except where that file documents a floating-point difference, and
an exception by its type alone, since the message text is the interpreter's.
"""

import json
import math
from collections.abc import Callable
from dataclasses import dataclass
from functools import cache, partial

from alignments import DATA, rows
from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment

from biofrills import alnutils, consensus, pairutils
from biofrills.stats import chisq

LIBRARY = DATA / "py27" / "library"
ALIGNMENT_FIXTURES = sorted(path.name for path in LIBRARY.glob("*.aln.json"))
# The sections of each library/<input>.aln.json, in file order
ALIGNMENT_SECTIONS = [
    "sequence_weights",
    "aa_counts",
    "aa_frequencies",
    "consensus",
    "blocks",
    "remove_empty_cols",
    "score_pairwise",
    "identity_abs",
    "identity_avg",
]

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


@dataclass
class Case:
    """One recorded call: where it is recorded, how to repeat it, its outcome.

    ``alignment`` names the alignment of a tie-batch.json case and is None
    elsewhere; ``index`` is the position of the case in its section, or None
    for a section that records a single outcome.  ``details`` holds the
    fields that describe the call in a differences.json entry.
    """

    file: str
    alignment: str | None
    section: str
    index: int | None
    details: dict
    rel_tol: float
    call: Callable
    recorded: dict

    @property
    def key(self):
        return (self.file, self.alignment, self.section, self.index)


def entry_key(entry):
    """The Case.key of the case that a differences.json entry names."""
    return (entry["file"], entry.get("alignment"), entry["section"], entry["case"])


@cache
def load(fname):
    with open(LIBRARY / fname, encoding="ascii") as handle:
        return json.load(handle)


@cache
def alignment(fname):
    return AlignIO.read(DATA / load(fname)["input"], "fasta")


def build_weights(aln, spec):
    """The sequence weights that a recorded weights name stands for.

    This repeats build_weights in tests/data/capture_py27.py, which defines
    the names; the two must agree.
    """
    if spec == "seqweights":
        return alnutils.sequence_weights(aln, "none")
    if spec == "unit":
        return [1] * len(aln)
    return spec


def as_records(records):
    return [[rec.id, rec.description, str(rec.seq)] for rec in records]


def weighted_call(func, aln, spec, kwargs):
    return func(aln, weights=build_weights(aln, spec), **kwargs)


def records_call(func, *args, **kwargs):
    return as_records(func(*args, **kwargs))


def library_cases(fname):
    """The recorded cases of one library/<input>.aln.json, in file order."""
    data = load(fname)
    aln = alignment(fname)
    file = "library/" + fname

    def case(section, index, details, rel_tol, call, recorded):
        return Case(file, None, section, index, details, rel_tol, call, recorded)

    for index, rec in enumerate(data["sequence_weights"]):
        kwargs = rec["kwargs"]
        rel_tol = SUM_REL_TOL if kwargs["scaling"] in SUMMED_SCALINGS else 0.0
        call = partial(alnutils.sequence_weights, aln, **kwargs)
        yield case(
            "sequence_weights", index, {"kwargs": kwargs}, rel_tol, call, rec["result"]
        )
    for section in ("aa_counts", "aa_frequencies"):
        func = getattr(alnutils, section)
        for index, rec in enumerate(data[section]):
            fractional = (
                section == "aa_frequencies" and rec["weights"] in FRACTIONAL_WEIGHTS
            )
            details = {"weights": rec["weights"], "kwargs": rec["kwargs"]}
            call = partial(weighted_call, func, aln, rec["weights"], rec["kwargs"])
            rel_tol = SUM_REL_TOL if fractional else 0.0
            yield case(section, index, details, rel_tol, call, rec["result"])
    for index, rec in enumerate(data["consensus"]):
        details = {"weights": rec["weights"], "kwargs": rec["kwargs"]}
        call = partial(
            weighted_call, consensus.consensus, aln, rec["weights"], rec["kwargs"]
        )
        yield case("consensus", index, details, 0.0, call, rec["result"])
    for index, rec in enumerate(data["blocks"]):
        details = {"weights": rec["weights"], "kwargs": rec["kwargs"]}
        call = partial(
            records_call, alnutils.blocks, aln, weights=rec["weights"], **rec["kwargs"]
        )
        yield case("blocks", index, details, 0.0, call, rec["result"])
    call = partial(records_call, alnutils.remove_empty_cols, aln)
    yield case("remove_empty_cols", None, {}, 0.0, call, data["remove_empty_cols"])
    seqs = rows(aln)
    for section in ("score_pairwise", "identity_abs", "identity_avg"):
        func = getattr(pairutils, section)
        for index, (i, j, recorded) in enumerate(data[section]):
            call = partial(func, seqs[i], seqs[j])
            yield case(section, index, {"rows": [i, j]}, 0.0, call, recorded)


def tie_batch_cases():
    """The recorded cases of library/tie-batch.json, alignment by alignment."""
    file = "library/tie-batch.json"
    data = load("tie-batch.json")
    groups = {}
    with open(DATA / data["input"], encoding="ascii") as handle:
        for rec in SeqIO.parse(handle, "fasta"):
            groups.setdefault(rec.id.split("_")[0], []).append(rec)
    weights_kwargs = data["sequence_weights_kwargs"]
    for rec in data["alignments"]:
        name = rec["id"]
        aln = MultipleSeqAlignment(groups[name])
        call = partial(alnutils.sequence_weights, aln, **weights_kwargs)
        yield Case(
            file,
            name,
            "sequence_weights",
            None,
            {"kwargs": weights_kwargs},
            SUM_REL_TOL,
            call,
            rec["sequence_weights"],
        )
        for index, kwargs in enumerate(data["consensus_kwargs"]):
            call = partial(consensus.consensus, aln, **kwargs)
            yield Case(
                file,
                name,
                "consensus",
                index,
                {"kwargs": kwargs},
                0.0,
                call,
                rec["consensus"][index],
            )


def chisqprob_cases():
    """The recorded cases of library/chisqprob.json."""
    for index, (x, df, recorded) in enumerate(load("chisqprob.json")["rows"]):
        yield Case(
            "library/chisqprob.json",
            None,
            "chisqprob",
            index,
            {"x": x, "df": df},
            LGAMMA_REL_TOL,
            partial(chisq.chisqprob, x, df),
            recorded,
        )


def all_cases():
    """Every recorded library case."""
    for fname in ALIGNMENT_FIXTURES:
        yield from library_cases(fname)
    yield from tie_batch_cases()
    yield from chisqprob_cases()


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
