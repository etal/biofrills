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
ARGPARSE_COMMANDS = {
    "seq-blocks",
    "seq-composition",
    "seq-dupes",
    "seq-intersect",
    "seq-orfs",
    "seq-pick",
    "seq-sort",
    "seq-split",
}
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
    """Standard error without the parts that vary between equivalent runs.

    - A log line's source line number, which moves with any edit to the
      logging module (tests/data/README.md, "Log line numbers").
    - The frames of a traceback; the first and last lines remain.
    """
    text = re.sub(r" \[@\d+\]:", " [@N]:", text)
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


# Three records; c shares no residue with a or b.
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


def no_partial_codon_warning(case):
    """seq-orfs translates each frame cut to whole codons.

    Biopython drops a trailing partial codon in either case, so standard output
    is unchanged and only its warning disappears (tests/data/README.md,
    "Command-line differences").
    """
    return case | {"stderr": ""}


def no_trailing_tab(case):
    """seq-identity no longer ends each row of its grid with a tab."""
    stdout = "\n".join(line.rstrip("\t") for line in case["stdout"].split("\n"))
    return case | {"stdout": stdout}


def removed_similar_option(case):
    """seq-sort -s, which sorted by increasing length, has been removed."""
    stderr = (
        "usage: seq-sort [-h] [-a | -l] [infile]\n"
        "seq-sort: error: unrecognized arguments: -s\n"
    )
    return case | {"returncode": 2, "stdout": "", "stderr": stderr}


# Cases whose current run differs from the recorded one, for the reasons
# given in tests/data/README.md: the function gives the expected run from the
# recorded one.
EXPECTED = {
    "seq-identity.cdc2.json": no_trailing_tab,
    "seq-identity.edge.json": no_trailing_tab,
    "seq-orfs.default.json": no_partial_codon_warning,
    "seq-orfs.m20.json": no_partial_codon_warning,
    "seq-pick.k3.json": tied_medoid,
    "seq-purge.c.json": keep_x_fragment,
    "seq-sort.s.json": removed_similar_option,
}


def expected_run(fname):
    """The expected run of a recorded case, by the name of its file."""
    recorded = load(fname)
    if fname not in EXPECTED:
        return recorded
    expected = EXPECTED[fname](recorded)
    assert expected != recorded, f"{fname} no longer differs as documented"
    return expected


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
    expected = expected_run(fname)

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
        ("seq-extract", [], None),
        ("seq-extract", ["a.fa", "2-4", "6"], None),
        (
            "seq-extract",
            ["a.fa", "5-3"],
            "Range must not end before it starts; got 5-3",
        ),
        ("seq-extract", ["0"], "Position must be at least 1; got 0"),
        ("seq-extract", ["-f"], None),
        ("seq-identity", [], None),
        ("seq-identity", ["a.fa", "fasta", "x"], None),
    ],
)
def test_usage_error(name, args, message, tmp_path):
    """A usage error prints the message, if any, and the usage; status 1."""
    proc = run(name, args, tmp_path)
    usage = docstring(name) + "\n"
    assert proc.returncode == 1
    assert proc.stderr == (f"{message}\n\n{usage}" if message else usage)


@pytest.mark.parametrize(
    "name, args, message",
    [
        ("seq-pick", ["-f", "0", "a.fa"], "argument nseqs: must be at least 1: 0"),
        ("seq-pick", ["3", "a.fa"], "one of the arguments -r/--random"),
        ("seq-pick", ["-r", "-f", "3", "a.fa"], "not allowed with argument"),
        (
            "seq-pick",
            ["-f", "3", "nope.fa"],
            "argument seqfname: no such file: nope.fa",
        ),
        ("seq-sort", ["-a", "-l"], "argument -l/--length: not allowed with"),
        ("seq-sort", ["nope.fa"], "argument infile: no such file: nope.fa"),
        ("seq-sort", ["."], "argument infile: not a file: ."),
        ("seq-composition", ["nope.fa"], "argument fname: no such file: nope.fa"),
        ("seq-split", ["-s"], "argument -s: expected one argument"),
        ("seq-split", ["-s", "0", "a.fa"], "argument -s: must be at least 1: 0"),
        ("seq-split", ["-s", "1.5", "a.fa"], "argument -s: not an integer: 1.5"),
        ("seq-split", ["-s", "2.5e0", "a.fa"], "argument -s: not an integer: 2.5e0"),
        ("seq-split", ["-s", "1e400", "a.fa"], "argument -s: out of range: 1e400"),
        ("seq-split", ["-b", "2", "-s", "2", "a.fa"], "not allowed with argument"),
        ("seq-split", ["a.fa"], "one of the arguments -b -s -nb -ns is required"),
        ("seq-split", ["-s", "2"], "the following arguments are required: fname"),
        ("seq-split", ["-s", "2", "nope.fa"], "argument fname: no such file: nope.fa"),
        ("seq-dupes", [], "the following arguments are required: fname"),
        ("seq-dupes", ["nope.fa"], "argument fname: no such file: nope.fa"),
        ("seq-intersect", [], "the following arguments are required: reference"),
        ("seq-intersect", ["a.fa", "r.fa"], "argument reference: no such file: r.fa"),
        ("seq-intersect", ["-c", "a.fa", "r.fa"], "argument reference: no such file"),
        # Standard input is /dev/null here, a character device, not a regular file
        ("seq-pick", ["-f", "3", "/dev/stdin"], "not a regular file: /dev/stdin"),
        ("seq-split", ["-s", "2", "/dev/stdin"], "not a regular file: /dev/stdin"),
        ("seq-intersect", ["-c", "a.fa", "/dev/stdin"], "not a regular file"),
    ],
)
def test_argparse_usage_error(name, args, message, tmp_path):
    """A usage error exits with status 2, names its cause, and writes no file.

    The scratch directory holds only a.fa, so a -c cache would show up.
    """
    (tmp_path / "a.fa").write_text(fasta(THREE))
    proc = run(name, args, tmp_path)
    assert proc.returncode == 2
    assert message in proc.stderr
    assert os.listdir(tmp_path) == ["a.fa"]


@pytest.mark.parametrize("fname", ["seq-dupes.dupes.json", "seq-sort.a.json"])
def test_reads_an_input_path_that_is_not_a_regular_file(fname, tmp_path):
    """A pipe reached through a path, here /dev/stdin, is read like a file.

    A shell's <(...) gives such a path.
    """
    expected = expected_run(fname)
    *options, infile = expected["argv"]
    read_end, write_end = os.pipe()
    # The inputs are smaller than a pipe's buffer, so the write does not block
    with open(write_end, "wb") as pipe:
        pipe.write((INPUTS / infile).read_bytes())
    with open(read_end, "rb") as stdin:
        proc = run(expected["script"], [*options, "/dev/stdin"], tmp_path, stdin)
    assert (proc.returncode, proc.stderr) == (0, expected["stderr"])
    assert proc.stdout == expected["stdout"]


def test_intersect_accepts_options_between_its_arguments(tmp_path):
    """seq-intersect q -r ref gives the recorded run of seq-intersect -r q ref."""
    expected = expected_run("seq-intersect.r.json")
    option, query, reference = expected["argv"]
    shutil.copytree(INPUTS, tmp_path, dirs_exist_ok=True)
    proc = run("seq-intersect", [query, option, reference], tmp_path)
    assert proc.returncode == expected["returncode"]
    assert proc.stdout == expected["stdout"]
    assert normalize_stderr(proc.stderr) == normalize_stderr(expected["stderr"])


@pytest.mark.parametrize(
    "fname",
    [
        "seq-blocks.cdc2.w.t09.json",
        "seq-composition.cdc2.json",
        "seq-orfs.m20.json",
        "seq-sort.a.json",
    ],
)
@pytest.mark.parametrize("dash", [False, True], ids=["implicit", "dash"])
def test_reads_standard_input(fname, dash, tmp_path):
    """Input on standard input, named '-' or not at all, gives the recorded run."""
    expected = expected_run(fname)
    *options, infile = expected["argv"]
    with open(INPUTS / infile, "rb") as stdin:
        proc = run(
            expected["script"], options + (["-"] if dash else []), tmp_path, stdin
        )
    assert proc.returncode == 0
    assert proc.stdout == expected["stdout"]
    assert normalize_stderr(proc.stderr) == normalize_stderr(expected["stderr"])


def test_identity_of_an_all_gap_sequence_is_zero(tmp_path):
    (tmp_path / "gaps.fa").write_text(
        fasta([("a", "AC-DE"), ("b", "-----"), ("c", "AC-DF")])
    )
    proc = run("seq-identity", ["gaps.fa"], tmp_path)
    assert (proc.returncode, proc.stderr) == (0, "")
    assert proc.stdout == (
        "a\n"
        "b\t0.00\n"
        "c\t75.00\t0.00\n"
        " \ta       b       c      \n"
        "\n"
        " Min: 0.00  ~  Mean: 25.00  ~  Max: 75.00\n"
    )


def test_identity_of_a_single_sequence_has_no_summary(tmp_path):
    proc = run("seq-identity", [str(INPUTS / "synth-single.aln.fasta")], tmp_path)
    assert (proc.returncode, proc.stderr) == (0, "")
    assert proc.stdout == "only\n    \tonly   \n"


def test_composition_of_an_input_without_sequence_characters(tmp_path):
    (tmp_path / "empty.fa").write_text(">a\n\n")
    proc = run("seq-composition", ["empty.fa"], tmp_path)
    assert (proc.returncode, proc.stdout) == (1, "")
    assert proc.stderr == "No sequence characters in the input.\n"


@pytest.mark.parametrize("parts", ["1000", "1E3"])
def test_split_into_more_parts_than_bytes(parts, tmp_path):
    """seq-split -nb with more parts than bytes writes a record per part.

    '1E3' checks that the exponent marker may be uppercase.
    """
    (tmp_path / "three.fa").write_text(fasta(THREE))
    proc = run("seq-split", ["-nb", parts, "three.fa"], tmp_path)
    assert (proc.returncode, proc.stderr) == (0, "")
    chunks = sorted(tmp_path.glob("three.fa.*"), key=lambda path: int(path.suffix[1:]))
    assert [path.read_text() for path in chunks] == [fasta([rec]) for rec in THREE]


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
