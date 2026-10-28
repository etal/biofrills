"""Trim an alignment to conserved "blocks", i.e. remove gappy columns."""

import argparse
import sys

from Bio import AlignIO

from biofrills import alnutils
from biofrills.cli import add_input_argument


def main(argv=None):
    AP = argparse.ArgumentParser(description=__doc__)
    add_input_argument(AP, "Input sequence alignment (FASTA format)")
    AP.add_argument(
        "-w", "--weight", action="store_true", help="Weight sequences by similarity."
    )
    AP.add_argument(
        "-t",
        "--threshold",
        type=float,
        default=0.5,
        help="Minimum proportion of non-gap characters to accept a column "
        "(default: %(default)s).",
    )

    args = AP.parse_args(argv)
    if not (0 < args.threshold <= 1):
        sys.exit("Non-gap threshold (-t) must be between 0 and 1.")

    aln = AlignIO.read(args.infile, "fasta")
    baln = alnutils.blocks(aln, args.threshold, args.weight)
    # Unwrapped FASTA
    for rec in baln:
        print(f">{rec.id} {rec.description}")
        print(rec.seq)
