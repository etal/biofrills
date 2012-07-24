"""Utilities for manipulating sequence sets.

"""
import logging
import os.path
from collections import Counter

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def clean_accessions(records):
    """Extract accession keys from SeqRecords.

    The id of the given records is processed to remove domain location info
    added by HMMer. Most other records won't have a '/' character in the FASTA
    header key.
    """
    return (rec.id.rsplit('/', 1)[0] for rec in records)


# ENH: use key_function to get gi/dbref numbers from both sides
def intersect_keys(keys, reffile, cache=False):
    """Extract SeqRecords from the index by matching keys."""
    # Build/load the index of reference sequences
    index = None
    if cache:
        refcache = reffile + '.sqlite'
        if os.path.exists(refcache):
            if os.stat(refcache).st_mtime < os.stat(reffile).st_mtime:
                logging.warn("Outdated cache; rebuilding index")
            else:
                try:
                    index = SeqIO.index_db(refcache)
                except Exception:
                    logging.warn("Skipping corrupted cache; rebuilding index")
                    index = None
    else:
        refcache = ':memory:'
    if index is None:
        # Rebuild the index, for whatever reason
        index = SeqIO.index_db(refcache, [reffile], 'fasta')

    # Extract records by key
    for key in keys:
        try:
            record = index[key]
        except LookupError:
            # Missing keys are rare, so it's faster not to check every time
            logging.info("No match: %s", repr(key))
            continue
        yield record


def aa_frequencies(seq, gap_chars='-.'):
    """Calculate the amino acid frequencies in a sequence set."""
    aa_counts = Counter(seq)
    # Don't count gaps
    for gap_char in gap_chars:
        if gap_char in aa_counts:
            del aa_counts[gap_char]
    # Reduce to frequencies
    scale = 1.0 / sum(aa_counts.values())
    return dict((aa, cnt * scale) for aa, cnt in aa_counts.iteritems())

