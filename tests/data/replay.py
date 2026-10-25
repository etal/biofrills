#!/usr/bin/env python3
"""Replay recorded reference cases against the current biofrills code.

Every library case recorded under py27/library is called again with the
current sources in this repository, and each case whose outcome differs from
the recorded one is written to differences.json beside this script.
README.md beside it explains each difference.

The cases and the comparison are those of tests/reference_cases.py, which
tests/test_reference.py uses as well: a value is compared exactly, except
for the floating-point differences that README.md documents, which are
compared with a relative tolerance and not written; an exception is compared
by its type only, since the message text is the interpreter's.  The script
requires Python 3.12 or later, whose built-in sum() the two t101 entries
depend on; the repository's default pixi environment provides a suitable
interpreter.

Usage, from the repository root:

    pixi run python tests/data/replay.py
"""

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TESTS = os.path.dirname(HERE)
REPO = os.path.dirname(TESTS)
OUT = os.path.join(HERE, "differences.json")

sys.path[:0] = [REPO, TESTS]

from reference_cases import all_cases, comparable, matches, outcome  # noqa: E402

from biofrills import consensus  # noqa: E402

assert os.path.dirname(os.path.abspath(consensus.__file__)) == os.path.join(
    REPO, "biofrills"
), consensus.__file__


def changed_columns(recorded, current):
    """1-based columns at which two equal-length strings differ."""
    old = recorded.get("value")
    new = current.get("value")
    if not (isinstance(old, str) and isinstance(new, str)) or len(old) != len(new):
        return None
    return [i + 1 for i, (a, b) in enumerate(zip(old, new, strict=True)) if a != b]


def difference(case, current):
    """One differences.json entry."""
    entry = {"file": case.file, "section": case.section}
    if case.alignment is not None:
        entry["alignment"] = case.alignment
    entry["case"] = case.index
    entry.update(case.details)
    entry["recorded"] = comparable(case.recorded)
    entry["current"] = current
    entry["columns"] = changed_columns(case.recorded, current)
    return entry


def main():
    """Replay every library case and rewrite differences.json."""
    found = []
    for case in all_cases():
        current = outcome(case.call)
        if not matches(current, case.recorded, case.rel_tol):
            found.append(difference(case, current))
    text = json.dumps(
        {"differences": found}, indent=1, separators=(",", ": "), ensure_ascii=True
    )
    with open(OUT, "w", encoding="ascii", newline="\n") as handle:
        handle.write(f"{text}\n")
    print(f"{len(found)} differing cases written to {OUT}")


if __name__ == "__main__":
    main()
