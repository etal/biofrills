"""Convert between sequence formats using Biopython's SeqIO.

Usage:
    seq-convert from_format to_format < infile > outfile

The formats are those that Bio.SeqIO reads and writes, named as in its
documentation: https://biopython.org/wiki/SeqIO#file-formats
"""

import sys

from Bio import SeqIO

from biofrills.cli import help_requested


def main(argv=None):
    args = sys.argv[1:] if argv is None else list(argv)
    if help_requested(args, __doc__):
        return 0
    if len(args) != 2:
        sys.exit(__doc__)

    count = SeqIO.convert(sys.stdin, args[0], sys.stdout, args[1])
    print(f"SeqIO: Converted {count} sequences", file=sys.stderr)
