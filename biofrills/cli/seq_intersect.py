"""Extract the sequences found in both given FASTA files.

Shared sequences are written to standard output, taking sequences from the
second file.

Usage:
    seq-intersect query.fa reference.fa > output.fa
    cat query.fa | seq-intersect reference.fa > output.fa

"""

import sys

from Bio import SeqIO

from biofrills import cmdutils, sequtils
from biofrills.cli import help_requested


def main(argv=None):
    args = sys.argv[1:] if argv is None else list(argv)
    if help_requested(args, __doc__):
        return 0
    cmdutils.log_config()
    # TODO - argparse

    if "--cache" in args:
        do_cache = True
        args.remove("--cache")
    elif "-c" in args:
        do_cache = True
        args.remove("-c")
    else:
        do_cache = False

    if "-r" in args:
        do_clean_acc = True
        args.remove("-r")
    else:
        do_clean_acc = False

    if len(args) == 1:
        qfile = sys.stdin
        reffile = args[0]
        do_close_qfile = False
    elif len(args) == 2:
        qfile = open(args[0])
        reffile = args[1]
        do_close_qfile = True
    else:
        sys.exit(__doc__)

    try:
        accs = (record.id for record in SeqIO.parse(qfile, "fasta"))
        with open(reffile):
            matches = sequtils.intersect_keys(
                accs, reffile, cache=do_cache, clean_accs=do_clean_acc
            )
            SeqIO.write(matches, sys.stdout, "fasta")
    finally:
        if do_close_qfile:
            qfile.close()
