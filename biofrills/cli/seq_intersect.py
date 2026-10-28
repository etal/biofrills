"""Extract the sequences found in both given FASTA files.

The reference records whose identifiers occur in the query are written to
standard output.
"""

import argparse
import sys

from Bio import SeqIO

from biofrills import cmdutils, sequtils
from biofrills.cli import add_input_argument, regular_file


def main(argv=None):
    AP = argparse.ArgumentParser(description=__doc__)
    add_input_argument(
        AP, "Sequences whose identifiers are selected (FASTA format)", name="query"
    )
    # regular_file rejects a missing file or a pipe before anything is indexed;
    # with -c, indexing either would leave a cache file behind
    AP.add_argument(
        "reference", type=regular_file, help="Sequences to select from (FASTA format)."
    )
    AP.add_argument(
        "-c",
        "--cache",
        action="store_true",
        help="Save the index of the reference sequences as REFERENCE.sqlite and "
        "reuse it on later runs; delete that file after changing REFERENCE or "
        "the setting of -r.",
    )
    AP.add_argument(
        "-r",
        dest="clean_accs",
        action="store_true",
        help="Ignore everything after the last '/' of an identifier, such as "
        "HMMer's '/START-END' domain suffix, when matching.",
    )
    # Options may also come between the two positional arguments
    args = AP.parse_intermixed_args(argv)
    cmdutils.log_config()

    accs = (record.id for record in SeqIO.parse(args.query, "fasta"))
    matches = sequtils.intersect_keys(
        accs, args.reference, cache=args.cache, clean_accs=args.clean_accs
    )
    SeqIO.write(matches, sys.stdout, "fasta")
