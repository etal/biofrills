"""Utilities for sequence logos.
"""
from cStringIO import StringIO

import weblogolib


def read_logodata(handle):
    """Get weblogo data for a sequence alignment.

    Returns a list of tuples: (posn, letter_counts, entropy, weight)
    """
    seqs = weblogolib.read_seq_data(handle)
    ldata = weblogolib.LogoData.from_seqs(seqs)
    letters = ldata.alphabet.letters()
    counts = ldata.counts.array
    logodata = []
    for i, coldata, entropy, weight in zip(range(len(counts)), counts,
                                           ldata.entropy, ldata.weight):
        cnts = dict((let, int(cnt))
                    for let, cnt in zip(letters, coldata))
        logodata.append((i + 1, cnts, entropy, weight))
    return logodata


def aln2logodata(aln):
    """Get weblogo data for an alignment object.

    Returns a list of tuples: (posn, letter_counts, entropy, weight)
    """
    handle = StringIO(aln.format('fasta'))
    logodata = read_logodata(handle)
    handle.close()
    return logodata


def logo_heights(aln):
    """
    Return a list of dicts: character to proportional height.
    """

def letter_scales(counts):
    """Convert letter counts to frequencies, sorted increasing.""" 
    try:
        scale = 1.0 / sum(counts.values())
    except ZeroDivisionError:
        # This logo is all gaps, nothing can be done
        return []
    freqs = [(aa, cnt*scale) for aa, cnt in counts.iteritems() if cnt]
    freqs.sort(key=lambda pair: pair[1])
    return freqs

