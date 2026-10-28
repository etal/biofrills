"""Translate the input sequences in 6 frames, split by ORFs.

Output: All of the translated ORFs, in FASTA format.
"""

import argparse

from Bio import SeqIO

from biofrills.cli import add_input_argument


def translate_six_frames(seq, genetic_code=1):
    """Generate the translations of a nucleotide sequence in all 6 frames.

    Each frame is cut to a whole number of codons before translation.
    """
    anti = seq.reverse_complement()
    for i in range(3):
        for frame in (seq[i:], anti[i:]):
            yield frame[: len(frame) - len(frame) % 3].translate(genetic_code)


def translate_orfs(record, long_enough):
    for frame in translate_six_frames(record.seq):
        for orf in filter(long_enough, map(str, frame.split("*"))):
            yield orf.strip("X")


def main(argv=None):
    AP = argparse.ArgumentParser(description=__doc__)
    add_input_argument(AP, "Input nucleotide sequences (FASTA format)")
    AP.add_argument(
        "-m",
        "--minlen",
        type=int,
        default=60,
        help="Minimum ORF length, in residues other than X (default: %(default)s).",
    )
    args = AP.parse_args(argv)

    def long_enough(seq):
        return len(seq.replace("X", "")) >= args.minlen

    for record in SeqIO.parse(args.infile, "fasta"):
        for idx, orf in enumerate(translate_orfs(record, long_enough)):
            print(f">{record.id}-{idx + 1}")
            print(orf)
