"""Utilities for the command line (i.e. scripts)."""

import logging


def log_config():
    """Log INFO and above to standard error, prefixed by module and line."""
    logging.basicConfig(
        format="%(module)s [@%(lineno)s]: %(message)s", level=logging.INFO
    )
