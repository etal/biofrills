"""Command-line programs installed with biofrills.

Each module provides the entry point of one console script, ``seq-NAME``, as a
function ``main(argv=None)``. The argument list ``argv`` holds the arguments
after the command name and defaults to ``sys.argv[1:]``.
"""

import argparse
import os.path
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


def input_file(text):
    """Argument type for the path of an input file that is read once.

    A missing path or a directory is a usage error rather than a traceback.
    Pipes, such as a shell's process substitution, are accepted.
    """
    if os.path.isdir(text):
        raise argparse.ArgumentTypeError(f"not a file: {text}")
    if not os.path.exists(text):
        raise argparse.ArgumentTypeError(f"no such file: {text}")
    return text


def regular_file(text):
    """Argument type for the path of an input file that must be a regular file.

    For inputs that are indexed, read more than once, or used to name output
    files, as in seq-pick, seq-split and seq-intersect; pipes and devices do
    not support these uses.
    """
    if not os.path.isfile(input_file(text)):
        raise argparse.ArgumentTypeError(f"not a regular file: {text}")
    return text


def input_name(text):
    """Argument type for an input file name, where '-' means standard input."""
    return text if text == "-" else input_file(text)


def input_path(text):
    """Argument type for an input file path, or standard input for '-'."""
    return sys.stdin if text == "-" else input_file(text)


def add_input_argument(parser, help_text, name="infile"):
    """Add an optional input file argument that defaults to standard input.

    argparse ignores the default of a required positional argument, so
    ``nargs="?"`` is what lets an omitted file mean standard input. The value
    is a path, which Biopython's readers open after parsing, or ``sys.stdin``.
    """
    parser.add_argument(
        name,
        nargs="?",
        type=input_path,
        default=sys.stdin,
        help=f"{help_text}; standard input if omitted or '-'.",
    )


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
