"""Print the frequencies of each letter in a sequence set."""
# TODO - move the calculation to module (take input stream, return freqs)

import fileinput
import sys
from collections import Counter

from biofrills.cli import help_requested


def main(argv=None):
    args = sys.argv[1:] if argv is None else list(argv)
    if help_requested(args, __doc__):
        return 0

    # Count character types
    counts = Counter()
    with fileinput.FileInput(args) as lines:
        for line in lines:
            if line.startswith(">") or not line.strip():
                continue
            counts.update(Counter(line.strip()))

    # Convert counts to frequencies
    scale = 1.0 / sum(counts.values())
    freqs = {char: cnt * scale for char, cnt in counts.items()}

    # Print a nice table
    for char, frq in sorted(freqs.items()):
        print(f"{char}: {frq:f}")
