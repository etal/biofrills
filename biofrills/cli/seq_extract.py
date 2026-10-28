"""Extract a range of columns from an alignment.

Usage:
    seq-extract [-f format] [infile] <range>

where <range> resembles: '1-5,10-20,34,35,40-'

Uses 1-based counting for column positions. The alignment is read from infile,
or from standard input if infile is omitted, in the Bio.AlignIO format given by
-f, "fasta" by default; the columns are written in the same format.
"""

import sys

from Bio import AlignIO

from biofrills.cli import help_requested


def _parse_position(text):
    if text.isdigit():
        pos = int(text) - 1
        if pos < 0:
            raise ValueError("Position must be at least 1; got " + text)
        return slice(pos, pos + 1)

    # ENH: allow ranges like "2:4"
    if not text.count("-") == 1:
        raise ValueError("String must resemble: '1', '2-4', '5-'")

    left, right = text.split("-")
    if left == "":
        left = 0
    else:
        left = int(left) - 1
        if left < 0:
            raise ValueError("Range must start from at least 1; got " + text)
    if right == "":
        return slice(left, None)
    else:
        right = int(right)
        if right <= left:
            raise ValueError("Range must not end before it starts; got " + text)
        return slice(left, right)


def parse_positions(text):
    """Convert a string of 1-based positions & ranges to Python.

    Return a list of slices (0-based indexing), one per position or range.

    Example:
    >>> parse_positions("6,7,10-20,25-")
    [slice(5, 6, None), slice(6, 7, None), slice(9, 20, None), slice(24, None, None)]
    """
    return [_parse_position(pos) for pos in text.split(",")]


def main(argv=None):
    args = sys.argv[1:] if argv is None else list(argv)
    if help_requested(args, __doc__):
        return 0

    fmt = "fasta"
    if args[:1] == ["-f"]:
        if len(args) < 2:
            sys.exit(__doc__)
        fmt = args[1]
        args = args[2:]

    if len(args) == 1:
        infile = sys.stdin
        positions = args[0]
    elif len(args) == 2:
        infile = args[0]
        positions = args[1]
    else:
        sys.exit(__doc__)

    try:
        indices = parse_positions(positions)
    except ValueError as exc:
        sys.exit(f"{exc}\n\n{__doc__}")
    aln = AlignIO.read(infile, fmt)
    first, rest = indices[0], indices[1:]
    out = aln[:, first]
    for index in rest:
        out += aln[:, index]

    AlignIO.write(out, sys.stdout, fmt)
