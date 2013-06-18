
from __future__ import division

import random

def uniqueof20(k, rep=10000):
    """Sample k times out of alphabet, how many different?"""
    alphabet = 'ACDEFGHIKLMNPQRSTVWY'
    reps = [len(set(random.choice(alphabet)
                    for i in range(k)))
            for j in range(rep)]
    return sum(reps) / len(reps) 

# ENH:
# B-H adjustment given a dict of {things: pvalues}
#   - if cutoff is given, return just the things that pass, otherwise all

# Medioidshift
# https://www.cs.cmu.edu/~yaser/new_medoid.htm

# Silhouette (to validate clusters -- use in Fammer?):
# https://en.wikipedia.org/wiki/Silhouette_%28clustering%29

