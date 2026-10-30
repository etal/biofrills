"""Statistical helpers: the chi-squared upper tail in the chisq module."""

import random


def uniqueof20(k, rep=10000):
    """Mean number of distinct amino acids among k drawn with replacement.

    Each of rep runs draws k residues uniformly from the 20 standard amino
    acids and counts the distinct types; the result is the mean over runs.
    The draws come from the unseeded module-level generator of random, so
    repeated calls give different estimates.  The exact expectation, which
    alnutils.sequence_weights inverts to estimate the effective number of
    sequences, is 20 * (1 - (19/20)**k).
    """
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    reps = [len(set(random.choice(alphabet) for i in range(k))) for j in range(rep)]
    return sum(reps) / len(reps)
