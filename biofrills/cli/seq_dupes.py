"""Find and print the sequences with identical IDs.

Each record whose ID repeats that of an earlier record is printed.
"""

import argparse

from Bio import SeqIO

from biofrills.cli import input_file


def main(argv=None):
    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument("fname", type=input_file, help="Input sequences (FASTA format).")
    args = AP.parse_args(argv)

    seen = set()
    for rec in SeqIO.parse(args.fname, "fasta"):
        acc = rec.id
        if acc in seen:
            print(rec.format("fasta"))
        else:
            seen.add(acc)
