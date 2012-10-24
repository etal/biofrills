BioFrills
=========

Utilities for biomolecular sequence analysis, underpinned by Biopython.

Requires Python version 2.7 or a recent PyPy. Optionally uses Cython to speed up
some CPU-intensive operations.

License: BSD

Installation
------------

This is an ordinary Python package. You can install it from source with the
setup.py script::

    python setup.py build
    python setup.py install

Some scripts that I find useful are in the scripts/ directory. By default these
are not installed, but you can include them by uncommenting the line in setup.py
that starts with ``scripts=glob``...

Alternatively, you can just copy those scripts into another directory in your
``$PATH``.
