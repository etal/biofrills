BioFrills
=========

Utilities for biomolecular sequence analysis, built on Biopython: a library of alignment and sequence-set functions, and 13 command-line programs for everyday work with FASTA files and alignments.
It supports the analyses of CladeCompare_ and Fammer_.

Requires Python 3.12 or later and Biopython 1.86 or later.

License: BSD 2-clause

Installation
------------

Install from a source checkout with pip::

    pip install .

Two optional extras add the dependencies of a few functions::

    pip install '.[graph,logo]'

``graph`` installs networkx, which ``seq-pick -c`` and ``-k`` and ``alnutils.to_graph`` need.
``logo`` installs weblogo, which the ``logoutils`` module needs.
Without networkx, ``seq-pick -c`` and ``-k`` exit with a message naming the extra.

Reading CMA files with ``seq-diff`` and ``seq-purge`` requires biocma 0.3 or later, which is not yet available from PyPI; install it from https://github.com/etal/biocma.

This version is not yet released.
The latest release on PyPI, 0.3.1, requires Python 2.7 and does not install the command-line programs.

Command-line programs
---------------------

Each program prints its usage with ``--help``.
Most programs read FASTA.
``seq-convert``, ``seq-extract`` and ``seq-identity`` also read other formats that Biopython supports, and ``seq-diff`` and ``seq-purge`` also read CMA.

=================== ===================================================================
Program             Purpose
=================== ===================================================================
``seq-blocks``      Remove gappy columns from an alignment, optionally weighting sequences by similarity.
``seq-composition`` Print the frequency of each character in a sequence set.
``seq-convert``     Convert sequences from one ``SeqIO`` format to another.
``seq-diff``        Print the records of one file that are absent from another, comparing headers, identifiers, sequences, or a combination.
``seq-dupes``       Print the records whose identifier repeats an earlier one.
``seq-extract``     Extract the alignment columns in a list of ranges, such as ``1-5,10-20,40-``.
``seq-identity``    Print the pairwise percent identities of an alignment as a grid, with their minimum, mean and maximum.
``seq-intersect``   Print the records of a reference file whose identifiers occur in a query file.
``seq-orfs``        Translate nucleotide sequences in six frames and print each stretch between stop codons of at least 60 residues by default.
``seq-pick``        Pick a number of records: the first, random, longest or shortest, or the most central of an alignment by clustering (``graph`` extra).
``seq-purge``       Remove duplicate records, comparing as ``seq-diff`` does or by one sequence containing another.
``seq-sort``        Sort records by decreasing length, or by identifier.
``seq-split``       Split a sequence file into parts of a given size, in bytes or records, or into a given number of similar parts.
=================== ===================================================================

Library
-------

``alnutils``
    Residue counts and frequencies of an alignment, optionally weighted; removal of gappy or empty columns; position-based sequence weights; a graph of pairwise identities.
``consensus``
    The consensus sequence of an alignment, by relative entropy or by plurality.
``pairutils``
    The score of two aligned sequences under BLOSUM62 with affine gap penalties, and their identity.
``sequtils``
    Extraction of records by identifier from an indexed FASTA file; residue frequencies of a sequence.
``logoutils``
    Per-column residue counts and entropy from weblogo (``logo`` extra).
``stats.chisq``
    The upper-tail probability of the chi-squared distribution, without SciPy.
``cmdutils``
    The logging configuration of the command-line programs.

Development
-----------

Two development workflows are defined in ``pyproject.toml``, and both lock files are committed.
Both install the package in editable mode with its ``graph`` and ``logo`` extras and the test and lint dependencies of the ``dev`` dependency group or pixi feature.
With uv_, which installs the dependencies from PyPI::

    uv sync
    uv run pytest
    uv run ruff check
    uv run ruff format

With pixi_, which installs them from conda-forge::

    pixi run test            # newest Python in the lock file
    pixi run -e py312 test   # Python 3.12, the oldest supported version
    pixi run lint            # ruff check
    pixi run fmt             # ruff format

Only weblogo comes from PyPI there, with three of its dependencies, because conda-forge lacks a recent enough build of it.

The tests compare the library and the command-line programs with reference outputs recorded under Python 2.7, and check properties that need no reference.
The programs are tested as installed: ``uv run`` reinstalls the package after a change to ``[project.scripts]``, and otherwise an existing environment needs ``uv sync`` or ``pixi install`` before the tests find them.
The ``py27-oracle`` pixi environment regenerates the reference outputs from a worktree of the pre-port sources, as described in ``tests/data/README.md``, which also lists every expected difference from them.

The GitHub Actions workflow in ``.github/workflows/ci.yml`` runs the lint checks once, and the tests on Linux and macOS at Python 3.12 and at the newest Python in ``pixi.lock``.
On Linux it also runs the tests with the versions in ``uv.lock`` and with the lowest declared versions of the dependencies, and checks the core modules and ``seq-pick`` without the extras.

Code is formatted with ``ruff format``, and the formatting-only commits are listed in ``.git-blame-ignore-revs``.
To skip them in ``git blame``::

    git config blame.ignoreRevsFile .git-blame-ignore-revs

.. _CladeCompare: https://github.com/etal/cladecompare
.. _Fammer: https://github.com/etal/fammer
.. _uv: https://docs.astral.sh/uv/
.. _pixi: https://pixi.prefix.dev/
