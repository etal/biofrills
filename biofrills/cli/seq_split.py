"""Split a FASTA file into multiple parts.

Each part of an input file fname is written to fname.K, where K increases with
each part; with -b and -nb the numbers K need not start at 0 or be consecutive.
"""

import argparse
import itertools
import os

from Bio import SeqIO

from biofrills.cli import regular_file


def sizeof(record):
    """Record size in bytes.

    = # chars in FASTA header + sequence, including newlines.

    Discards the description.
    """
    return len(record.id) + len(record.seq) + 3


def split_bytes(fname, byte_limit):
    """Split the sequence set into chunks of about ``byte_limit`` bytes."""
    curr_bytes = 0

    def byte_grouper(record):
        nonlocal curr_bytes
        curr_bytes += sizeof(record)
        # NB: keys are non-sequential if byte_limit is smaller than a record
        return curr_bytes // byte_limit

    return itertools.groupby(SeqIO.parse(fname, "fasta"), byte_grouper)


def split_count(fname, count_limit):
    """Split the sequence set into chunks of <= ``count_limit`` sequences."""
    curr_seqs = -1

    def count_grouper(record):
        nonlocal curr_seqs
        curr_seqs += 1
        return curr_seqs // count_limit

    return itertools.groupby(SeqIO.parse(fname, "fasta"), count_grouper)


def split_n_bytes(fname, num_chunks):
    tot_bytes = os.stat(fname).st_size
    return split_bytes(fname, max(1, tot_bytes // num_chunks))


def split_n_count(fname, num_chunks):
    with open(fname, "rb") as handle:
        tot_recs = sum(1 for line in handle if line.startswith(b">"))
    chunk_size, remainder = divmod(tot_recs, num_chunks)
    if remainder:
        chunk_size += 1
    return split_count(fname, chunk_size)


def write_chunk(fname, suffix, records):
    out_fname = f"{fname}.{suffix}"
    with open(out_fname, "w+") as outfile:
        outfile.writelines(f">{r.id}\n{r.seq}\n" for r in records)
    print("Wrote", out_fname)


def size_arg(text):
    """Parse a chunk size or count of at least 1, such as '100' or '1e3'."""
    try:
        number = float(text) if "e" in text.lower() else int(text)
        value = int(number)
    except ValueError:
        raise argparse.ArgumentTypeError(f"not an integer: {text}") from None
    except OverflowError:
        raise argparse.ArgumentTypeError(f"out of range: {text}") from None
    if value != number:
        raise argparse.ArgumentTypeError(f"not an integer: {text}")
    if value < 1:
        raise argparse.ArgumentTypeError(f"must be at least 1: {text}")
    return value


SPLITTERS = {
    "b": (split_bytes, "Split into chunks of about N bytes each."),
    "s": (split_count, "Split into chunks of at most N records each."),
    "nb": (split_n_bytes, "Split into about N chunks of similar size in bytes."),
    "ns": (split_n_count, "Split into at most N chunks of similar record counts."),
}


def main(argv=None):
    AP = argparse.ArgumentParser(
        description=__doc__, epilog="N may be written in exponent notation, as in 1e6."
    )
    mode = AP.add_mutually_exclusive_group(required=True)
    for name, (_func, help_text) in SPLITTERS.items():
        mode.add_argument(f"-{name}", metavar="N", type=size_arg, help=help_text)
    AP.add_argument(
        "fnames",
        metavar="fname",
        nargs="+",
        type=regular_file,
        help="FASTA file to split.",
    )
    args = AP.parse_args(argv)

    # The mode group is required and mutually exclusive, so exactly one is set
    (name,) = [name for name in SPLITTERS if getattr(args, name) is not None]
    split = SPLITTERS[name][0]
    size = getattr(args, name)
    for fname in args.fnames:
        for key, chunk in split(fname, size):
            write_chunk(fname, key, chunk)
