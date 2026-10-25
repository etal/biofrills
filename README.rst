BioFrills
=========

Utilities for biomolecular sequence analysis, underpinned by Biopython.

Requires Python 3.12 or later.

License: BSD 2-clause

Installation
------------

This is an ordinary Python package. You can install it from a source checkout
with pip::

    pip install .

This version is not yet released. It installs 13 command-line utilities,
``seq-blocks`` through ``seq-split``; each prints its usage with ``--help``.
``seq-pick -c`` and ``seq-pick -k`` require the ``graph`` extra, and the CMA
input of ``seq-diff`` and ``seq-purge`` requires biocma 0.3 or later, which is
not yet available from PyPI. The latest release on PyPI, 0.3.1, requires
Python 2.7 and does not install the utilities.

Development
-----------

Two development workflows are defined in ``pyproject.toml``, and both lock
files are committed. Both install the package in editable mode with its
``graph`` and ``logo`` extras and the test and lint dependencies of the
``dev`` dependency group or pixi feature. With uv_, which installs the
dependencies from PyPI::

    uv sync
    uv run pytest
    uv run ruff check
    uv run ruff format

With pixi_, which installs them from conda-forge::

    pixi run test            # newest Python in the lock file
    pixi run -e py312 test   # Python 3.12, the oldest supported version
    pixi run lint            # ruff check
    pixi run fmt             # ruff format

Only weblogo comes from PyPI there, with three of its dependencies, because
conda-forge lacks a recent enough build of it.

The tests compare the library and the command-line utilities with reference
outputs recorded under Python 2.7, and check properties that need no
reference. The utilities are tested as installed: ``uv run`` reinstalls the
package after a change to ``[project.scripts]``, and otherwise an existing
environment needs ``uv sync`` or ``pixi install`` before the tests find them.
The ``py27-oracle`` pixi environment regenerates the reference outputs from a
worktree of the pre-port sources, as described in ``tests/data/README.md``.

Code is formatted with ``ruff format``, and the formatting-only commits are
listed in ``.git-blame-ignore-revs``. To skip them in ``git blame``::

    git config blame.ignoreRevsFile .git-blame-ignore-revs

.. _uv: https://docs.astral.sh/uv/
.. _pixi: https://pixi.prefix.dev/
