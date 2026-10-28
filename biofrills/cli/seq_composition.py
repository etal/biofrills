"""Print the frequencies of each character in a sequence set.

Every character of the sequences is counted, including gaps; upper and lower
case are counted separately.
"""

import argparse
import fileinput
import itertools
import sys

from biofrills.cli import input_name
from biofrills.sequtils import aa_frequencies


def sequence_chars(lines):
    """The characters of the sequence lines of FASTA text."""
    for line in lines:
        if not line.startswith(">"):
            yield from line.strip()


def main(argv=None):
    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument(
        "fnames",
        metavar="fname",
        nargs="*",
        type=input_name,
        help="Input sequences (FASTA format); standard input if none or '-'.",
    )
    args = AP.parse_args(argv)

    with fileinput.FileInput(args.fnames) as lines:
        chars = sequence_chars(lines)
        first = next(chars, None)
        if first is None:
            sys.exit("No sequence characters in the input.")
        freqs = aa_frequencies(itertools.chain([first], chars), gap_chars="")

    for char, frq in sorted(freqs.items()):
        print(f"{char}: {frq:f}")
