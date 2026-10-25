#!/usr/bin/env python3
"""Replay recorded reference cases against the current biofrills code.

Every consensus case recorded under py27/ (the ``consensus`` section of each
library/<input>.json and every alignment in library/tie-batch.json) is
called again with the current sources in this repository, and each case whose
outcome differs from the recorded one is written to differences.json beside
this script.  README.md beside it explains each difference.

Outcomes are compared as the capture records them: a value exactly, an
exception by its type only, since the message text is the interpreter's.
The script requires Python 3.12 or later, whose built-in sum() the two
t101 entries depend on; the repository's default pixi environment provides
a suitable interpreter.

Usage, from the repository root:

    pixi run python tests/data/replay.py
"""

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
FIXTURES = os.path.join(HERE, "py27", "library")
OUT = os.path.join(HERE, "differences.json")

sys.path.insert(0, REPO)

from Bio import AlignIO, SeqIO  # noqa: E402
from Bio.Align import MultipleSeqAlignment  # noqa: E402

from biofrills import alnutils, consensus  # noqa: E402

assert os.path.dirname(os.path.abspath(consensus.__file__)) == os.path.join(
    REPO, "biofrills"
), consensus.__file__


def build_weights(aln, spec):
    """The sequence weights that a recorded weights name stands for.

    This repeats build_weights in capture_py27.py, which defines the names,
    and in tests/test_reference.py; the three must agree.
    """
    if spec == "seqweights":
        return alnutils.sequence_weights(aln, "none")
    if spec == "unit":
        return [1] * len(aln)
    return spec


def outcome(func, *args, **kwargs):
    """Call func; return {"value": result} or {"error": {"type": ...}}."""
    try:
        return {"value": func(*args, **kwargs)}
    except Exception as exc:
        return {"error": {"type": type(exc).__name__}}


def comparable(result):
    """A recorded outcome, reduced to what the comparison uses."""
    if "error" in result:
        return {"error": {"type": result["error"]["type"]}}
    return {"value": result["value"]}


def changed_columns(recorded, current):
    """1-based columns at which two equal-length strings differ."""
    old = recorded.get("value")
    new = current.get("value")
    if not (isinstance(old, str) and isinstance(new, str)) or len(old) != len(new):
        return None
    return [i + 1 for i, (a, b) in enumerate(zip(old, new, strict=True)) if a != b]


def difference(case, recorded, current):
    """One differences.json entry."""
    entry = dict(case)
    entry["recorded"] = comparable(recorded)
    entry["current"] = current
    entry["columns"] = changed_columns(recorded, current)
    return entry


def replay_library(fname):
    """Differences in the consensus section of one library/<input>.json."""
    with open(os.path.join(FIXTURES, fname), encoding="ascii") as handle:
        data = json.load(handle)
    aln = AlignIO.read(os.path.join(HERE, data["input"]), "fasta")
    found = []
    for index, case in enumerate(data["consensus"]):

        def call(spec=case["weights"], kwargs=case["kwargs"]):
            return consensus.consensus(aln, weights=build_weights(aln, spec), **kwargs)

        current = outcome(call)
        if current != comparable(case["result"]):
            found.append(
                difference(
                    {
                        "file": "library/" + fname,
                        "section": "consensus",
                        "case": index,
                        "weights": case["weights"],
                        "kwargs": case["kwargs"],
                    },
                    case["result"],
                    current,
                )
            )
    return found


def replay_tie_batch():
    """Differences in library/tie-batch.json."""
    fname = "tie-batch.json"
    with open(os.path.join(FIXTURES, fname), encoding="ascii") as handle:
        data = json.load(handle)
    groups = {}
    for rec in SeqIO.parse(os.path.join(HERE, data["input"]), "fasta"):
        groups.setdefault(rec.id.split("_")[0], []).append(rec)
    found = []
    for alignment in data["alignments"]:
        aln = MultipleSeqAlignment(groups[alignment["id"]])
        for index, kw in enumerate(data["consensus_kwargs"]):
            current = outcome(consensus.consensus, aln, **kw)
            recorded = alignment["consensus"][index]
            if current != comparable(recorded):
                found.append(
                    difference(
                        {
                            "file": "library/" + fname,
                            "section": "consensus",
                            "alignment": alignment["id"],
                            "case": index,
                            "kwargs": kw,
                        },
                        recorded,
                        current,
                    )
                )
    return found


def main():
    """Replay every consensus case and rewrite differences.json."""
    found = []
    for fname in sorted(os.listdir(FIXTURES)):
        if fname.endswith(".aln.json"):
            found.extend(replay_library(fname))
    found.extend(replay_tie_batch())
    text = json.dumps(
        {"differences": found}, indent=1, separators=(",", ": "), ensure_ascii=True
    )
    with open(OUT, "w", encoding="ascii", newline="\n") as handle:
        handle.write(f"{text}\n")
    print(f"{len(found)} differing cases written to {OUT}")


if __name__ == "__main__":
    main()
