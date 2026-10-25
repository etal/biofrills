"""Split a FASTA file into multiple parts.

According to:
    -b: # bytes per chunk (approx.)
    -s: # records per chunk
    -nb: produce # similar-sized chunks, in terms of bytes
    -ns: produce # similar-sized chunks, in terms of records

Usage:
    seq-split [-b|-s|-nb|-ns] <size> <fname> [<fname2> ...]

"""

import itertools
import os
import subprocess
import sys

from Bio import SeqIO

from biofrills.cli import help_requested


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
    return split_bytes(fname, tot_bytes // num_chunks)


def split_n_count(fname, num_chunks):
    cnt = subprocess.Popen(["grep", "-c", r"^>", fname], stdout=subprocess.PIPE)
    tot_recs = int(cnt.stdout.readline())
    chunk_size, remainder = divmod(tot_recs, num_chunks)
    if remainder:
        chunk_size += 1
    return split_count(fname, chunk_size)


def write_chunk(fname, suffix, records):
    out_fname = f"{fname}.{suffix}"
    with open(out_fname, "w+") as outfile:
        outfile.writelines(f">{r.id}\n{r.seq}\n" for r in records)
    print("Wrote", out_fname)


def main(argv=None):
    args = sys.argv[1:] if argv is None else list(argv)
    if help_requested(args, __doc__):
        return 0

    chunk_func = {
        "-b": split_bytes,
        "-s": split_count,
        "-nb": split_n_bytes,
        "-ns": split_n_count,
    }

    try:
        opt, value, fnames = args[0], args[1], args[2:]
        assert opt[0] == "-"
        assert opt in chunk_func
        value = int(float(value)) if "e" in value else int(value)
    except Exception:
        sys.exit(__doc__)

    for fname in fnames:
        for key, chunk in chunk_func[opt](fname, value):
            write_chunk(fname, key, chunk)
