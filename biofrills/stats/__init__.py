"""Statistical helpers: the chi-squared upper tail in the chisq module."""

import random


def uniqueof20(k, rep=10000):
    """Mean number of distinct amino acids among k drawn with replacement.

    Each of rep runs draws k residues uniformly from the 20 standard amino
    acids and counts the distinct types; the result is the mean over runs.
    A table of these Monte Carlo estimates for k from 0 to 100 is the
    expectk table of alnutils.sequence_weights. The exact expectation is
    20 * (1 - (19/20)**k).
    """
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    reps = [len(set(random.choice(alphabet) for i in range(k))) for j in range(rep)]
    return sum(reps) / len(reps)
