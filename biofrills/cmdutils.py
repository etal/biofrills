"""Utilities for the command line (i.e. scripts)."""

import logging


def log_config(verbose=1):
    """Set up logging the way I like it."""
    # ENH:
    # - do print levelname before DEBUG and WARNING
    # - instead of %module, name the currently running script
    #   - make a subclass of logging.handlers.X instead?
    #   - tweak %root?
    #   - take __file__ as an argument?
    if verbose == 0:
        level = logging.WARNING
        fmt = "%(module)s: %(message)s"
    elif verbose == 1:
        level = logging.INFO
        fmt = "%(module)s [@%(lineno)s]: %(message)s"
    else:
        level = logging.DEBUG
        fmt = "%(module)s [%(lineno)s]: %(levelname)s: %(message)s"
    logging.basicConfig(format=fmt, level=level)

