"""Calculate sequence identity percentages across an alignment.

Usage:
    seq-identity <infile> [format]

where format names a Bio.AlignIO alignment format, "fasta" by default.

Per OrthoMCL:
percent_match = number of matched characters (in the shorter sequence)
                / length of the shorter sequence

Each row of the printed grid gives the identities of one sequence with each of
the preceding ones. A line of column labels follows the grid, then, for two or
more sequences, the minimum, mean and maximum over all pairs. A pair in which
one sequence is all gaps has identity 0.
"""

import sys

from Bio import AlignIO

from biofrills.cli import help_requested
from biofrills.pairutils import identity_abs


def ungap_len(seq):
    return len(seq) - seq.count("-") - seq.count(".")


def pct_identity(aseq, bseq):
    assert len(aseq) == len(bseq)
    shorter_len = min(ungap_len(aseq), ungap_len(bseq))
    if not shorter_len:
        return 0.0
    return 100.0 * identity_abs(aseq, bseq) / shorter_len


def main(argv=None):
    args = sys.argv[1:] if argv is None else list(argv)
    if help_requested(args, __doc__):
        return 0
    if not 1 <= len(args) <= 2:
        sys.exit(__doc__)

    infname = args[0]
    fmt = args[1] if len(args) == 2 else "fasta"

    records = list(AlignIO.read(infname, fmt))
    label_width = max(len(rec.id) for rec in records)

    identities = []
    for i, record in enumerate(records):
        aseq = str(record.seq)
        row = [record.id.ljust(label_width)]
        for j in range(i):
            pct_id = pct_identity(aseq, str(records[j].seq))
            row.append(f"{pct_id:2.2f}")
            identities.append(pct_id)
        print("\t".join(row))

    # Column labels
    print(" " * label_width + "\t" + " ".join([rec.id[:7].ljust(7) for rec in records]))

    if identities:
        print()
        print(
            f" Min: {min(identities):2.2f}  ~  "
            f"Mean: {sum(identities) / len(identities):2.2f}  ~  "
            f"Max: {max(identities):2.2f}"
        )
