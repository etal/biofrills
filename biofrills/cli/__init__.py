"""Command-line programs installed with biofrills.

Each module provides the entry point of one console script, ``seq-NAME``, as a
function ``main(argv=None)``. The argument list ``argv`` holds the arguments
after the command name and defaults to ``sys.argv[1:]``.
"""

import sys

# The keys by which seq-diff and seq-purge compare records, by option letter,
# for Biopython SeqRecords and for biocma's CMA sequence dictionaries
RECORD_KEYS = {
    "h": lambda rec: (rec.id, rec.description),
    "i": lambda rec: rec.id,
    "s": lambda rec: str(rec.seq),
    "a": lambda rec: (rec.id, rec.description, str(rec.seq)),
    "b": lambda rec: (rec.id, str(rec.seq)),
}
CMA_KEYS = {
    "h": lambda rec: (rec["id"], rec["description"]),
    "i": lambda rec: rec["id"],
    "s": lambda rec: rec["seq"],
    "a": lambda rec: (rec["id"], rec["description"], rec["seq"]),
    "b": lambda rec: (rec["id"], rec["seq"]),
}


def help_requested(args, doc):
    """Print a command's usage and return True if args include '--help'.

    Commands that parse their arguments by hand call this first. Only the long
    spelling is recognized, because seq-diff and seq-purge use '-h' to select
    a mode.
    """
    if "--help" not in args:
        return False
    print(doc.strip())
    return True


def import_cma():
    """Import biocma's CMA module, or exit with a message if it is missing.

    Only the CMA input paths of seq-diff and seq-purge need biocma, and the
    Python 3 release it requires, 0.3, is not yet available from PyPI.
    """
    try:
        from biocma import cma
    except ImportError:
        sys.exit(
            "Reading CMA files requires biocma 0.3 or later: "
            "https://github.com/etal/biocma"
        )
    return cma
