"""Alignments and helpers shared by the test modules."""

from pathlib import Path

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

DATA = Path(__file__).parent / "data"
INPUTS = DATA / "inputs"

# The aligned inputs on which every library function succeeds: the three
# MAPGAPS alignments and the synthetic edge cases with more than one row
# and no residue type that occurs only in lowercase.
ALIGNED_INPUTS = [
    "cdc2-domain.aln.fasta",
    "scttlre-domain.aln.fasta",
    "fikk.aln.fasta",
    "synth-edge.aln.fasta",
    "synth-ties.aln.fasta",
    "synth-x.aln.fasta",
]


def read_alignment(name):
    """Read a FASTA alignment from tests/data/inputs."""
    return AlignIO.read(INPUTS / name, "fasta")


def make_alignment(rows):
    """An alignment of the given sequence strings, named s0, s1, ..."""
    return MultipleSeqAlignment(
        SeqRecord(Seq(row), id=f"s{index}") for index, row in enumerate(rows)
    )


def rows(records):
    """The sequences of an alignment or an iterable of records, as strings."""
    return [str(rec.seq) for rec in records]


def columns(records):
    """The columns of an alignment or of aligned records, as strings."""
    return ["".join(col) for col in zip(*rows(records), strict=True)]
