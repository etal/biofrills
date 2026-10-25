"""Run the installed console scripts and compare them with recorded runs.

Each case of tests/data/py27/scripts runs as the capture ran it: as a separate
process, in a scratch copy of tests/data/inputs, with PYTHONHASHSEED set to 0
and the recorded input file, if any, on standard input. Exit status, standard
output and the files the command wrote are compared exactly, and standard
error after the normalizations below. tests/data/README.md documents every
reason a current run may differ from a recorded one; the cases that differ
are listed in EXPECTED and TIED.
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
ARGPARSE_COMMANDS = {"seq-blocks", "seq-orfs", "seq-sort"}
HAVE_BIOCMA = importlib.util.find_spec("biocma") is not None


def command(name):
    """The path of an installed console script, beside the interpreter."""
    path = Path(sys.executable).with_name(name)
    assert path.exists(), f"{name} is not installed; run 'uv sync' or 'pixi install'"
    return str(path)


def docstring(name):
    """The docstring of the module that implements a command."""
    return importlib.import_module(f"biofrills.cli.{name.replace('-', '_')}").__doc__


def run(name, args, cwd, stdin=subprocess.DEVNULL):
    # argparse wraps help text to COLUMNS, and the width of a pipe is 80
    env = os.environ | {"PYTHONHASHSEED": "0", "COLUMNS": "80"}
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


def sort_tied(text, ids):
    """FASTA text with the records of the given ids sorted among themselves.

    The other records keep their positions.
    """
    records = split_records(text)
    headers = tuple(f">{id_}\n" for id_ in ids)
    places = [i for i, record in enumerate(records) if record.startswith(headers)]
    for i, record in zip(places, sorted(records[i] for i in places), strict=True):
        records[i] = record
    return "".join(records)


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
# Cases in which the given records may appear in either order. CDC2-1_Sm and
# CDC2-2_Sm have clustering coefficients that differ in the last place, and
# which is larger depends on the iteration order of sets in networkx, that is,
# on the hash seed (tests/data/README.md, "Consensus ties").
TIED = {"seq-pick.c5.json": ("CDC2-1_Sm", "CDC2-2_Sm")}


def comparable(case, fname):
    result = {
        "returncode": case["returncode"],
        "stdout": case["stdout"],
        "stderr": normalize_stderr(case["stderr"]),
        "written_files": case["written_files"],
    }
    if fname in TIED:
        result["stdout"] = sort_tied(result["stdout"], TIED[fname])
    return result


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
    assert comparable(actual, fname) == comparable(expected, fname)


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
        ("seq-pick", ["-f", "five"], None),
    ],
)
def test_usage_error(name, args, message, tmp_path):
    """A usage error prints the message, if any, and the usage; status 1."""
    proc = run(name, args, tmp_path)
    usage = docstring(name) + "\n"
    assert proc.returncode == 1
    assert proc.stderr == (f"{message}\n\n{usage}" if message else usage)


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
