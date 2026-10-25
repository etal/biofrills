"""Run the installed console scripts and compare them with recorded runs.

Each case of tests/data/py27/scripts runs as the capture ran it: as a separate
process, in a scratch copy of tests/data/inputs, with PYTHONHASHSEED set to 0
and the recorded input file, if any, on standard input. Exit status, standard
output and the files the command wrote are compared exactly, and standard
error after the normalizations below. tests/data/README.md documents every
reason a current run may differ from a recorded one; the cases that differ
are listed in EXPECTED.
"""

import importlib.util
import json
import os
import re
import shutil
import subprocess
import sys
from functools import cache
from pathlib import Path

import pytest
from alignments import DATA, INPUTS

SCRIPT_CASES = DATA / "py27" / "scripts"
COMMANDS = [
    "seq-blocks",
    "seq-composition",
    "seq-convert",
    "seq-diff",
    "seq-dupes",
    "seq-extract",
    "seq-identity",
    "seq-intersect",
    "seq-orfs",
    "seq-pick",
    "seq-purge",
    "seq-sort",
    "seq-split",
]
# The commands that parse their arguments with argparse
ARGPARSE_COMMANDS = {"seq-blocks", "seq-orfs", "seq-pick", "seq-sort"}
HAVE_BIOCMA = importlib.util.find_spec("biocma") is not None


def command(name):
    """The path of an installed console script, beside the interpreter."""
    path = Path(sys.executable).with_name(name)
    assert path.exists(), f"{name} is not installed; run 'uv sync' or 'pixi install'"
    return str(path)


def docstring(name):
    """The docstring of the module that implements a command."""
    return importlib.import_module(f"biofrills.cli.{name.replace('-', '_')}").__doc__


def run(name, args, cwd, stdin=subprocess.DEVNULL, hashseed="0"):
    # argparse wraps help text to COLUMNS, and the width of a pipe is 80
    env = os.environ | {"PYTHONHASHSEED": hashseed, "COLUMNS": "80"}
    # The recorded runs used the default warning filters
    env.pop("PYTHONWARNINGS", None)
    return subprocess.run(
        [command(name), *args],
        cwd=cwd,
        env=env,
        stdin=stdin,
        capture_output=True,
        text=True,
        check=False,
    )


@cache
def load(fname):
    with open(SCRIPT_CASES / fname, encoding="utf-8") as handle:
        return json.load(handle)


def normalize_stderr(text):
    """Standard error without the parts that name the interpreter's files.

    - A log line's source line number, which moves with any edit to the
      logging module (tests/data/README.md, "Log line numbers").
    - The file, line and echoed source line of a warning, which belong to
      the Biopython installation.
    - The frames of a traceback; the first and last lines remain.
    """
    text = re.sub(r" \[@\d+\]:", " [@N]:", text)
    text = re.sub(r"^.*?:\d+: (\w*Warning: .*\n)(?:  .*\n)?", r"\1", text, flags=re.M)
    if text.startswith("Traceback (most recent call last):"):
        lines = text.splitlines(keepends=True)
        text = lines[0] + lines[-1]
    return text


def split_records(text):
    """The records of FASTA text, each with its header line."""
    return re.split(r"(?m)^(?=>)", text)


def fasta(records):
    """FASTA text of the given (identifier, sequence) pairs."""
    return "".join(f">{name}\n{seq}\n" for name, seq in records)


# Three records for seq-pick; c shares no residue with a or b.
THREE = [("a", "MKVLAYRGHQ"), ("b", "MKVLAYRGHE"), ("c", "WWWWWWWWWW")]


def keep_x_fragment(case):
    """seq-purge -c keeps a record that the recorded rule purged.

    tests/data/README.md, "Ambiguous residues in seq-purge -c": the record
    x_fragment_of_1 is kept, in length order before x_novel.
    """
    kept = [
        record
        for record in split_records(load("seq-purge.a.json")["stdout"])
        if record.startswith(">x_fragment_of_1\n")
    ]
    assert len(kept) == 1
    stdout = case["stdout"].replace(">x_novel\n", kept[0] + ">x_novel\n")
    stderr = case["stderr"].replace("Purged 6 sequences", "Purged 5 sequences")
    return case | {"stdout": stdout, "stderr": stderr}


def tied_medoid(case):
    """seq-pick -k chooses the other of two identical sequences as a medoid.

    CDC2-1_Sm and CDC2-2_Sm have the same aligned sequence, so the k-medoids
    cost ties exactly between them (tests/data/README.md, "Consensus ties").
    """
    stdout = case["stdout"].replace(">CDC2-1_Sm\n", ">CDC2-2_Sm\n")
    return case | {"stdout": stdout}


# Cases whose current run differs from the recorded one, for the reasons
# given in tests/data/README.md: the function gives the expected run from the
# recorded one.
EXPECTED = {
    "seq-purge.c.json": keep_x_fragment,
    "seq-pick.k3.json": tied_medoid,
}


def comparable(case):
    return {
        "returncode": case["returncode"],
        "stdout": case["stdout"],
        "stderr": normalize_stderr(case["stderr"]),
        "written_files": case["written_files"],
    }


@pytest.mark.parametrize("fname", sorted(path.name for path in SCRIPT_CASES.iterdir()))
def test_recorded_run(fname, tmp_path):
    recorded = load(fname)
    if not HAVE_BIOCMA and any(arg.endswith(".cma") for arg in recorded["argv"]):
        pytest.skip("reading CMA files requires biocma")
    expected = recorded
    if fname in EXPECTED:
        expected = EXPECTED[fname](recorded)
        assert expected != recorded, f"{fname} no longer differs as documented"

    shutil.copytree(INPUTS, tmp_path, dirs_exist_ok=True)
    before = set(os.listdir(tmp_path))
    if recorded["stdin"] is None:
        proc = run(recorded["script"], recorded["argv"], tmp_path)
    else:
        with open(tmp_path / recorded["stdin"], "rb") as stdin:
            proc = run(recorded["script"], recorded["argv"], tmp_path, stdin)
    written = {
        name: (tmp_path / name).read_text(encoding="utf-8")
        for name in sorted(set(os.listdir(tmp_path)) - before)
    }
    actual = {
        "returncode": proc.returncode,
        "stdout": proc.stdout,
        "stderr": proc.stderr,
        "written_files": written,
    }
    assert comparable(actual) == comparable(expected)


@pytest.mark.parametrize("name", COMMANDS)
def test_help(name, tmp_path):
    proc = run(name, ["--help"], tmp_path)
    assert (proc.returncode, proc.stderr) == (0, "")
    doc = docstring(name).strip()
    if name in ARGPARSE_COMMANDS:
        assert proc.stdout.startswith(f"usage: {name} ")
        assert doc.splitlines()[0] in proc.stdout
    else:
        assert proc.stdout == doc + "\n"


@pytest.mark.parametrize(
    "name, args, message",
    [
        ("seq-diff", ["-i", "a", "b", "c"], "Wrong number of arguments (4)"),
        ("seq-diff", ["-x", "a", "b"], "Invalid option '-x'"),
        ("seq-purge", ["-i", "a", "b"], "Wrong number of arguments (3)"),
        ("seq-purge", ["-x", "a"], "Invalid option '-x'"),
        ("seq-convert", ["fasta"], None),
        ("seq-split", ["-s"], None),
    ],
)
def test_usage_error(name, args, message, tmp_path):
    """A usage error prints the message, if any, and the usage; status 1."""
    proc = run(name, args, tmp_path)
    usage = docstring(name) + "\n"
    assert proc.returncode == 1
    assert proc.stderr == (f"{message}\n\n{usage}" if message else usage)


@pytest.mark.parametrize(
    "args, message",
    [
        (["-f", "0", "three.fa"], "argument nseqs: must be at least 1: 0"),
        (["3", "three.fa"], "one of the arguments -r/--random"),
        (["-r", "-f", "3", "three.fa"], "not allowed with argument"),
    ],
)
def test_pick_usage_error(args, message, tmp_path):
    proc = run("seq-pick", args, tmp_path)
    assert proc.returncode == 2
    assert message in proc.stderr


def test_split_keeps_every_record(tmp_path):
    """seq-split -s 2 writes 5 records as chunks of 2, 2 and 1."""
    ids = [f"r{i}" for i in range(5)]
    (tmp_path / "five.fa").write_text("".join(f">{i}\nACDE\n" for i in ids))
    proc = run("seq-split", ["-s", "2", "five.fa"], tmp_path)
    assert proc.returncode == 0
    chunks = sorted(tmp_path.glob("five.fa.*"))
    assert [path.name for path in chunks] == ["five.fa.0", "five.fa.1", "five.fa.2"]
    found = [
        line[1:]
        for path in chunks
        for line in path.read_text().splitlines()
        if line.startswith(">")
    ]
    assert found == ids
    assert [path.read_text().count(">") for path in chunks] == [2, 2, 1]


@pytest.mark.parametrize("nseqs", [2, 3, 10])
def test_pick_random_picks_distinct_records(nseqs, tmp_path):
    (tmp_path / "three.fa").write_text(fasta(THREE))
    proc = run("seq-pick", ["-r", str(nseqs), "three.fa"], tmp_path)
    assert (proc.returncode, proc.stderr) == (0, "")
    picked = split_records(proc.stdout)[1:]
    assert len(picked) == min(nseqs, 3)
    assert len(set(picked)) == len(picked)
    assert set(picked) <= set(split_records(fasta(THREE)))


def test_pick_first_stops_at_the_end_of_the_file(tmp_path):
    (tmp_path / "three.fa").write_text(fasta(THREE))
    proc = run("seq-pick", ["-f", "10", "three.fa"], tmp_path)
    assert (proc.returncode, proc.stderr, proc.stdout) == (0, "", fasta(THREE))


def test_pick_k_means_gives_an_unrelated_record_its_own_cluster(tmp_path):
    # c has identity 0 to every other record
    (tmp_path / "three.fa").write_text(fasta(THREE))
    proc = run("seq-pick", ["-k", "2", "three.fa"], tmp_path)
    assert (proc.returncode, proc.stderr) == (0, "")
    assert proc.stdout == fasta([THREE[2], THREE[1]])


@pytest.mark.parametrize("option", ["-c", "-k"])
def test_pick_by_clustering_from_a_single_record(option, tmp_path):
    (tmp_path / "one.fa").write_text(fasta(THREE[:1]))
    proc = run("seq-pick", [option, "1", "one.fa"], tmp_path)
    assert (proc.returncode, proc.stderr, proc.stdout) == (0, "", fasta(THREE[:1]))


@pytest.mark.parametrize("hashseed", ["1", "2", "3"])
def test_pick_cluster_order_does_not_depend_on_the_hash_seed(hashseed, tmp_path):
    # CDC2-1_Sm and CDC2-2_Sm have the same aligned sequence; with hash seeds
    # 1 and 2, networkx computes clustering coefficients for them that differ
    # in the last place.
    recorded = load("seq-pick.c5.json")
    shutil.copy(INPUTS / "cdc2-domain.aln.fasta", tmp_path)
    proc = run("seq-pick", recorded["argv"], tmp_path, hashseed=hashseed)
    assert (proc.returncode, proc.stdout) == (0, recorded["stdout"])


@pytest.mark.skipif(HAVE_BIOCMA, reason="biocma is installed")
@pytest.mark.parametrize(
    "name, args",
    [("seq-diff", ["-i", "fikk.cma", "fikk.cma"]), ("seq-purge", ["-i", "fikk.cma"])],
)
def test_cma_without_biocma(name, args, tmp_path):
    shutil.copy(INPUTS / "fikk.cma", tmp_path)
    proc = run(name, args, tmp_path)
    assert proc.returncode == 1
    assert proc.stderr == (
        "Reading CMA files requires biocma 0.3 or later: "
        "https://github.com/etal/biocma\n"
    )
