"""Compare the library with the outputs recorded under Python 2.7.

tests/reference_cases.py enumerates the recorded cases and defines how a
current outcome is compared with a recorded one.  A case listed in
tests/data/differences.json is expected to give the outcome listed there
instead of the recorded one; tests/data/README.md explains each entry.
"""

import json
from functools import cache

import pytest
from alignments import DATA
from reference_cases import (
    ALIGNMENT_FIXTURES,
    ALIGNMENT_SECTIONS,
    all_cases,
    chisqprob_cases,
    comparable,
    entry_key,
    library_cases,
    load,
    matches,
    outcome,
    tie_batch_cases,
)


@cache
def difference_entries():
    with open(DATA / "differences.json", encoding="ascii") as handle:
        return json.load(handle)["differences"]


@cache
def differences():
    """The expected current outcome of each case in differences.json."""
    return {entry_key(entry): entry["current"] for entry in difference_entries()}


def mismatches(cases):
    """The cases whose current outcome is not the expected one."""
    expected_changes = differences()
    found = []
    for case in cases:
        actual = outcome(case.call)
        expected = expected_changes.get(case.key, case.recorded)
        if not matches(actual, expected, case.rel_tol):
            found.append((case.key, case.details, actual))
    return found


def check(found):
    """Fail with the first few mismatches, if any."""
    assert not found, f"{len(found)} mismatches, first: {found[:3]}"


@pytest.mark.parametrize("fname", ALIGNMENT_FIXTURES)
@pytest.mark.parametrize("section", ALIGNMENT_SECTIONS)
def test_library(fname, section):
    check(mismatches(c for c in library_cases(fname) if c.section == section))


def test_every_recorded_section_is_compared():
    recorded = {key for fname in ALIGNMENT_FIXTURES for key in load(fname)}
    enumerated = {
        c.section for fname in ALIGNMENT_FIXTURES for c in library_cases(fname)
    }
    assert recorded - {"input", "ids"} == enumerated == set(ALIGNMENT_SECTIONS)


@pytest.mark.parametrize("section", ["sequence_weights", "consensus"])
def test_tie_batch(section):
    check(mismatches(c for c in tie_batch_cases() if c.section == section))


def test_chisqprob():
    check(mismatches(chisqprob_cases()))


def test_differences_name_changed_cases():
    """Each entry of differences.json names a recorded case that changed."""
    cases = {case.key: case for case in all_cases()}
    found = []
    for entry in difference_entries():
        case = cases.get(entry_key(entry))
        if (
            case is None
            or entry["recorded"] != comparable(case.recorded)
            or matches(entry["current"], case.recorded, case.rel_tol)
        ):
            found.append(entry)
    check(found)
