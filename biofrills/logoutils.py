"""Utilities for sequence logos."""

from io import StringIO

import weblogo
from weblogo.seq import unambiguous_protein_alphabet


def read_logodata(handle):
    """Get weblogo data for a sequence alignment.

    Returns a list of tuples: (posn, letter_counts, entropy, weight)
    """
    seqs = weblogo.read_seq_data(handle, alphabet=unambiguous_protein_alphabet)
    ldata = weblogo.LogoData.from_seqs(seqs)
    letters = ldata.alphabet.letters()
    counts = ldata.counts.array
    logodata = []
    for i, coldata, entropy, weight in zip(
        range(len(counts)), counts, ldata.entropy, ldata.weight, strict=True
    ):
        cnts = dict((let, int(cnt)) for let, cnt in zip(letters, coldata, strict=True))
        logodata.append((i + 1, cnts, entropy, weight))
    return logodata


def aln2logodata(aln):
    """Get weblogo data for an alignment object.

    Returns a list of tuples: (posn, letter_counts, entropy, weight)
    """
    handle = StringIO(format(aln, "fasta"))
    logodata = read_logodata(handle)
    handle.close()
    return logodata
