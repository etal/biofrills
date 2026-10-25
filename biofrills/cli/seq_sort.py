"""Sort a set of sequences."""
# [ENH] sort by distance from center (sort-from-center)
#  see center_of.py ; reuse/refactor code
#      -> sort sequences by decreasing total score

import argparse
import sys

from Bio import SeqIO


def main(argv=None):
    AP = argparse.ArgumentParser(description=__doc__)

    AP.add_argument(
        "infile",
        default=sys.stdin,
        type=argparse.FileType("r"),
        # XXX make the stdin default work
        help="FASTA file name.",
    )
    AP.add_argument(
        "-a",
        "--alphabet",
        action="store_true",
        help="Sort alphabetically by sequence ID.",
    )
    AP.add_argument(
        "-l", "--length", action="store_true", help="Sort by sequence length."
    )
    AP.add_argument(
        "-s",
        "--similar",
        action="store_true",
        help="Sort by average similarity to other sequences.",
    )
    AP.add_argument(
        "-q",
        "--quiet",
        action="store_false",
        dest="verbose",
        default=True,
        help="Don't print status messages on stdout.",
    )

    args = AP.parse_args(argv)

    records = SeqIO.parse(args.infile, "fasta")
    if args.alphabet:

        def sortkey(rec):
            return str(rec.id)
    elif args.similar:
        sortkey = len  # TODO -- centroid calculation
    else:

        def sortkey(rec):
            return -len(rec.seq)  # Decreasing

    sortrecs = sorted(records, key=sortkey)
    SeqIO.write(sortrecs, sys.stdout, "fasta")
