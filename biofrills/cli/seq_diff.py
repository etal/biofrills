"""Show the difference in sequences between two FASTA or CMA files.

Usage:
    seq-diff [-h|-i|-s|-a|-b] seqfile1 seqfile2
    cat seqfile | seq-diff [-h|-i|-s|-a|-b]

Difference is actually the subtraction ``seqfile1 - seqfile2``, i.e. sequences
that appear in seqfile1 but not seqfile2. To get the sequences that are unique
to seqfile2, swap the arguments.

For CMA, the filename must be given and end with ".cma".

Options for how "duplicate" is determined:
    -h      Same header (sequence identifier and description)
    -i      Same identifiers (first header token)
    -s      Same sequences
    -a      Same headers and sequences ("all")
    -b      Same identifiers and sequences ("both")
"""

import sys

from Bio import SeqIO

from biofrills.cli import CMA_KEYS, RECORD_KEYS, help_requested, import_cma


def differ(infile1, infile2, key_func):
    """Diff between two sets of SeqRecords.

    Equivalents are defined according to mode (default: all fields the same).
    """
    # For each sequence file, create a dictionary mapping the 'key' (id,
    # description and/or sequence) to the full record.
    #   * SeqIO.to_dict raises an exception if records aren't unique
    seqrecs1 = SeqIO.to_dict(SeqIO.parse(infile1, "fasta"), key_func)
    seqrecs2 = SeqIO.to_dict(SeqIO.parse(infile2, "fasta"), key_func)
    # Use set subtraction to find the diff
    keys1 = set(seqrecs1.keys())
    keys2 = set(seqrecs2.keys())
    gain_recs = [seqrecs2[key] for key in sorted(keys2 - keys1)]
    loss_recs = [seqrecs1[key] for key in sorted(keys1 - keys2)]
    return gain_recs, loss_recs


def differ_cma(infile1, infile2, key_func):
    """Like differ(), but for CMA files."""
    cma = import_cma()
    # NB: only handles the first block of the CMA file
    block1 = cma.read(infile1)
    block2 = cma.read(infile2)

    seqs1 = {key_func(seq): seq for seq in block1["sequences"]}
    assert len(seqs1) == len(block1["sequences"]), (
        str(len(block1["sequences"]) - len(seqs1)) + " duplicated keys in file 1"
    )

    seqs2 = {key_func(seq): seq for seq in block2["sequences"]}
    assert len(seqs2) == len(block2["sequences"]), (
        str(len(block2["sequences"]) - len(seqs2)) + " duplicated keys in file 2"
    )

    # Use set subtraction to find the diff
    keys1 = set(seqs1.keys())
    keys2 = set(seqs2.keys())
    gain_recs = [seqs2[key] for key in sorted(keys2 - keys1)]
    loss_recs = [seqs1[key] for key in sorted(keys1 - keys2)]

    block2["sequences"] = gain_recs
    block1["sequences"] = loss_recs
    return block2, block1


def main(argv=None):
    args = sys.argv[1:] if argv is None else list(argv)
    if help_requested(args, __doc__):
        return 0

    if len(args) == 2:
        infile1 = sys.stdin
        infile2 = args[1]
    elif len(args) == 3:
        infile1 = args[1]
        infile2 = args[2]
    else:
        sys.exit(f"Wrong number of arguments ({len(args)})\n\n{__doc__}")

    opt = args[0]
    if not (len(opt) == 2 and opt[0] == "-" and opt[1] in "hisab"):
        sys.exit(f"Invalid option '{opt}'\n\n{__doc__}")
    opt = opt[1]

    if infile2.endswith(".cma"):
        cma = import_cma()
        _gain, loss = differ_cma(infile1, infile2, CMA_KEYS[opt])
        # XXX For now: just write the subtraction to stdout
        cma.write(loss, sys.stdout)
    else:
        _gain, loss = differ(infile1, infile2, RECORD_KEYS[opt])
        # XXX For now: just write the subtraction to stdout
        SeqIO.write(loss, sys.stdout, "fasta")
