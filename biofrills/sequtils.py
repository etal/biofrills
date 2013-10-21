"""Utilities for manipulating sequence sets.

"""
import logging
import os.path
from collections import Counter

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def clean_accession(rec_id):
    """Extract an accession key from one SeqRecord.

    The id of the given record is processed to remove domain location info
    added by HMMer. Most other records won't have a '/' character in the FASTA
    header key.
    """
    return rec_id.rsplit('/', 1)[0]


# Obsolete, but left here for backward compatibility for now
def clean_accessions(records):
    """Extract accession keys from an iterable of SeqRecords."""
    return (clean_accession(rec.id) for rec in records)


def intersect_keys(keys, reffile, cache=False, clean_accs=False):
    """Extract SeqRecords from the index by matching keys.

    keys - an iterable of sequence identifiers/accessions to select
    reffile - name of a FASTA file to extract the specified sequences from
    cache - save an index of the reference FASTA sequence offsets to disk?
    clean_accs - strip HMMer extensions from sequence accessions?
    """
    # Build/load the index of reference sequences
    index = None
    if cache:
        refcache = reffile + '.sqlite'
        if os.path.exists(refcache):
            if os.stat(refcache).st_mtime < os.stat(reffile).st_mtime:
                logging.warn("Outdated cache; rebuilding index")
            else:
                try:
                    index = (SeqIO.index_db(refcache,
                                            key_function=clean_accession)
                             if clean_accs
                             else SeqIO.index_db(refcache))

                except Exception:
                    logging.warn("Skipping corrupted cache; rebuilding index")
                    index = None
    else:
        refcache = ':memory:'
    if index is None:
        # Rebuild the index, for whatever reason
        index = (SeqIO.index_db(refcache, [reffile], 'fasta',
                                key_function=clean_accession)
                 if clean_accs
                 else SeqIO.index_db(refcache, [reffile], 'fasta'))

    # Extract records by key
    if clean_accs:
        keys = (clean_accession(k) for k in keys)
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

