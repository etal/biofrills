"""Sort a set of sequences, by decreasing length or by identifier."""

import argparse
import sys

from Bio import SeqIO

from biofrills.cli import add_input_argument


def main(argv=None):
    AP = argparse.ArgumentParser(description=__doc__)
    add_input_argument(AP, "Input sequences (FASTA format)")
    order = AP.add_mutually_exclusive_group()
    order.add_argument(
        "-a",
        "--alphabet",
        dest="order",
        action="store_const",
        const="alphabet",
        help="Sort alphabetically by sequence ID.",
    )
    order.add_argument(
        "-l",
        "--length",
        dest="order",
        action="store_const",
        const="length",
        help="Sort by decreasing sequence length (the default).",
    )
    AP.set_defaults(order="length")
    args = AP.parse_args(argv)

    records = SeqIO.parse(args.infile, "fasta")
    if args.order == "alphabet":

        def sortkey(rec):
            return str(rec.id)
    else:

        def sortkey(rec):
            return -len(rec.seq)

    sortrecs = sorted(records, key=sortkey)
    SeqIO.write(sortrecs, sys.stdout, "fasta")
