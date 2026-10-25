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

This version is not yet released, and its port from Python 2 is in progress:
the library modules run on Python 3, but the command-line utilities in
scripts/ are still written for Python 2. The latest release on PyPI, 0.3.1,
requires Python 2.7. The scripts are not installed with the package.

Development
-----------

Two development workflows are defined in ``pyproject.toml``, and both lock
files are committed. With uv_, which installs the dependencies from PyPI::

    uv sync --all-extras
    uv run pytest

With pixi_, which installs them from conda-forge::

    pixi run test            # newest Python in the lock file
    pixi run -e py312 test   # Python 3.12, the oldest supported version

Only weblogo comes from PyPI there, with three of its dependencies, because
conda-forge lacks a recent enough build of it.

The tests compare the library with reference outputs recorded under Python
2.7, and check properties that need no reference. The tests of ``logoutils``
need the ``logo`` extra and are skipped without it. The ``py27-oracle`` pixi
environment regenerates the reference outputs from a worktree of the pre-port
sources, as described in ``tests/data/README.md``.

.. _uv: https://docs.astral.sh/uv/
.. _pixi: https://pixi.prefix.dev/
