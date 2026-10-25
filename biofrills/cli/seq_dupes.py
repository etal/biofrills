"""Find and print the sequences with identical IDs."""

import sys

from Bio import SeqIO

from biofrills.cli import help_requested


def main(argv=None):
    args = sys.argv[1:] if argv is None else list(argv)
    if help_requested(args, __doc__):
        return 0

    seen = set()
    for rec in SeqIO.parse(args[0], "fasta"):
        acc = rec.id
        if acc in seen:
            print(rec.format("fasta"))
        else:
            seen.add(acc)
