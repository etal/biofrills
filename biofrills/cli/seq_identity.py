"""Calculate sequence identity percentages across an alignment.

Per OrthoMCL:
percent_match = number of matched characters (in the shorter sequence)
                / length of the shorter sequence

If a multiple sequence alignment is given, print a grid of pairwise identities.
"""

import sys

from Bio import AlignIO

from biofrills.cli import help_requested


def ungap_len(seq):
    return len(seq) - seq.count("-") - seq.count(".")


def pct_identity(aseq, bseq):
    assert len(aseq) == len(bseq)
    shorter_len = min(ungap_len(aseq), ungap_len(bseq))
    ident = sum(a == b for a, b in zip(aseq, bseq, strict=True) if a not in "-.")
    return 100.0 * ident / shorter_len


def main(argv=None):
    args = sys.argv[1:] if argv is None else list(argv)
    if help_requested(args, __doc__):
        return 0

    infname = args[0]
    fmt = args[1] if len(args) == 2 else "fasta"

    records = list(AlignIO.read(infname, fmt))
    label_width = max(len(rec.id) for rec in records)

    identities = []
    for i, record in enumerate(records):
        aseq = str(record.seq)
        print(record.id.ljust(label_width) + "\t", end="")
        for j in range(i):
            bseq = str(records[j].seq)
            pct_id = pct_identity(aseq, bseq)
            print(f"{pct_id:2.2f}\t", end="")
            identities.append(pct_id)
        print("\n", end="")

    # Column labels
    print(" " * label_width + "\t" + " ".join([rec.id[:7].ljust(7) for rec in records]))

    print()
    print(
        f" Min: {min(identities):2.2f}  ~  "
        f"Mean: {sum(identities) / len(identities):2.2f}  ~  "
        f"Max: {max(identities):2.2f}"
    )
