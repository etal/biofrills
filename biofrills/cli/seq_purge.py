"""Remove duplicate sequences from a FASTA or CMA file.

For CMA, the filename must be given and end with ".cma".

Usage:
    seq-purge [-h|-i|-s|-a|-b|-c] seqfile
    cat seqfile | seq-purge [-h|-i|-s|-a|-b|-c]

Options for how "duplicate" is determined:
    -h      Same header (sequence identifier and description)
    -i      Same identifiers (first header token)
    -s      Same sequences
    -a      Same headers and sequences ("all")
    -b      Same identifiers and sequences ("both")
    -c      One sequence contains another ("contains")
"""

import re
import sys

from Bio import SeqIO

from biofrills.cli import CMA_KEYS, RECORD_KEYS, help_requested, import_cma


def purger(seqrecs, identify):
    """Filter out duplicate SeqRecords.

    Duplicates are defined according to mode (default: all fields the same).
    """
    seen = set()
    counter = 0
    for rec in seqrecs:
        try:
            this = identify(rec)
        except Exception:
            print(rec)
            raise
        if this not in seen:
            seen.add(this)
            yield rec
        else:
            counter += 1
    print(f"Purged {counter} sequence{'' if counter == 1 else 's'}", file=sys.stderr)


def container(seqrecs, is_cma):
    """Filter out SeqRecords that are completely contained by another.

    Records are visited in order of decreasing sequence length, and a record
    is dropped if its sequence occurs within a sequence already kept. An 'X'
    in the record's sequence matches any single character; an 'X' in a kept
    sequence matches only an 'X'. fammer's update_fasta.ungap_and_unique
    applies the same rule after removing gaps and converting to uppercase.
    """
    seen = set()
    counter = 0
    if is_cma:

        def getseq(rec):
            return rec["seq"]
    else:

        def getseq(rec):
            return str(rec.seq)

    for rec in sorted(seqrecs, key=lambda r: len(getseq(r)), reverse=True):
        recseq = getseq(rec)
        if "X" in recseq:
            # Each ambiguous residue matches any single character
            pattern = re.compile(re.escape(recseq).replace("X", "."))
            redundant = any(pattern.search(longseq) for longseq in seen)
        else:
            redundant = any(recseq in longseq for longseq in seen)
        if redundant:
            counter += 1
        else:
            seen.add(recseq)
            yield rec
    print(f"Purged {counter} sequence{'' if counter == 1 else 's'}", file=sys.stderr)


def main(argv=None):
    args = sys.argv[1:] if argv is None else list(argv)
    if help_requested(args, __doc__):
        return 0

    if len(args) == 1:
        infile = sys.stdin
    elif len(args) == 2:
        infile = args[1]
    else:
        sys.exit(f"Wrong number of arguments ({len(args)})\n\n{__doc__}")

    opt = args[0]
    if not (len(opt) == 2 and opt[0] == "-" and opt[1] in "hisabc"):
        sys.exit(f"Invalid option '{opt}'\n\n{__doc__}")

    if isinstance(infile, str) and infile.endswith(".cma"):
        # CMA file
        cma = import_cma()
        block = cma.read(infile)
        if opt == "-c":
            block["sequences"] = list(container(block["sequences"], True))
        else:
            block["sequences"] = list(purger(block["sequences"], CMA_KEYS[opt[1]]))
        cma.write(block, sys.stdout)
    else:
        # Assume FASTA
        inseqs = SeqIO.parse(infile, "fasta")
        if opt == "-c":
            results = container(inseqs, False)
        else:
            results = purger(inseqs, RECORD_KEYS[opt[1]])
        SeqIO.write(results, sys.stdout, "fasta")
