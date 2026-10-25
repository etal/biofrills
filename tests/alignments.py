"""Alignments and helpers shared by the test modules."""

from pathlib import Path

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

DATA = Path(__file__).parent / "data"
INPUTS = DATA / "inputs"

# Every aligned input, including the synthetic edge cases.
ALIGNED_INPUTS = sorted(path.name for path in INPUTS.glob("*.aln.fasta"))


def read_alignment(name):
    """Read a FASTA alignment from tests/data/inputs."""
    return AlignIO.read(INPUTS / name, "fasta")


def make_records(rows):
    """SeqRecords of the given sequence strings, named s0, s1, ...

    Unlike the rows of an alignment, these may differ in length.
    """
    return [SeqRecord(Seq(row), id=f"s{index}") for index, row in enumerate(rows)]


def make_alignment(rows):
    """An alignment of the given sequence strings, named s0, s1, ..."""
    return MultipleSeqAlignment(make_records(rows))


def rows(records):
    """The sequences of an alignment or an iterable of records, as strings."""
    return [str(rec.seq) for rec in records]


def columns(records):
    """The columns of an alignment or of aligned records, as strings."""
    return ["".join(col) for col in zip(*rows(records), strict=True)]
