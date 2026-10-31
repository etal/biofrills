"""Utilities for manipulating sequence sets."""

import logging
import os
import sqlite3
import uuid
from collections import Counter
from contextlib import closing

from Bio import SeqIO

logger = logging.getLogger(__name__)


def clean_accession(rec_id):
    """Extract an accession key from a sequence identifier.

    The given record ID is processed to remove domain location info added by
    HMMer. Most other records won't have a '/' character in the FASTA header
    key.
    """
    return rec_id.rsplit("/", 1)[0]


def intersect_keys(keys, reffile, cache=False, clean_accs=False):
    """Extract SeqRecords from the index by matching keys.

    keys - an iterable of sequence identifiers/accessions to select
    reffile - path of a FASTA file to extract the specified sequences from
    cache - save an index of the reference FASTA sequence offsets to disk?
    clean_accs - strip HMMer extensions from sequence accessions?

    reffile may be a str or an os.PathLike object. The cache is the reference
    file's path with '.sqlite' appended. It is reused only if it was built
    from a reference file of the same size and modification time and with the
    same setting of clean_accs; otherwise it is rebuilt.
    """
    key_function = clean_accession if clean_accs else None
    if cache:
        index = _cached_index(reffile, clean_accs, key_function)
    else:
        index = SeqIO.index_db(
            ":memory:", [reffile], "fasta", key_function=key_function
        )
    try:
        # Extract records by key
        if clean_accs:
            keys = (clean_accession(k) for k in keys)
        for key in keys:
            try:
                record = index[key]
            except LookupError:
                # Missing keys are rare, so it's faster not to check every time
                logger.info("No match: %s", repr(key))
                continue
            yield record
    finally:
        _close_index(index)


def _close_index(index):
    """Close an index from SeqIO.index_db and its database connection."""
    index.close()
    # close() closes only the sequence files (Biopython 1.88), and Biopython
    # offers no public way to close the private database connection
    index._con.close()


# Biopython looks up the rows of an index's meta_data table by key and ignores
# any other key, so this row does not affect loading.
_SIGNATURE_KEY = "biofrills_signature"


def _cached_index(reffile, clean_accs, key_function):
    """Return an index of reffile, reusing the cache beside it where possible.

    Biopython loads any existing file at the path given to index_db, so a
    rebuild is written to a new file that then replaces the cache. A matching
    cache of zero records is valid but cannot be loaded, so it is kept and the
    reference is indexed in memory instead.
    """
    refcache = os.fspath(reffile) + ".sqlite"
    # Taken before the reference is read, so that a change made while the
    # index is built makes the cache outdated
    stat = os.stat(reffile)
    signature = (
        f"clean_accs={bool(clean_accs)} size={stat.st_size} mtime_ns={stat.st_mtime_ns}"
    )
    if os.path.exists(refcache):
        try:
            cached_signature, count = _read_metadata(refcache)
            matches = cached_signature == signature
            if matches and int(count) != 0:
                return SeqIO.index_db(refcache, key_function=key_function)
        except Exception:
            logger.warning("Skipping corrupted cache; rebuilding index")
        else:
            if matches:
                # Biopython's _load_index raises TypeError on MAX(_ROWID_) of
                # an empty offset table and leaves its database connection
                # open (Biopython 1.74 through at least 1.88). This return is
                # in the else clause so that an error from indexing the
                # reference is not reported as a corrupt cache.
                return SeqIO.index_db(
                    ":memory:", [reffile], "fasta", key_function=key_function
                )
            logger.warning(
                "Cache does not record that it was built from this reference "
                "file with this setting of accession cleaning; rebuilding index"
            )

    # Built in the directory of the cache, since Biopython records the
    # reference path relative to the index file
    tmpcache = f"{refcache}.{uuid.uuid4().hex}.tmp"
    index = None
    try:
        index = SeqIO.index_db(tmpcache, [reffile], "fasta", key_function=key_function)
        with closing(sqlite3.connect(tmpcache)) as con, con:
            con.execute(
                "INSERT INTO meta_data (key, value) VALUES (?, ?);",
                (_SIGNATURE_KEY, signature),
            )
        os.replace(tmpcache, refcache)
    except BaseException:
        if index is not None:
            _close_index(index)
        if os.path.exists(tmpcache):
            os.remove(tmpcache)
        raise
    return index


def _read_metadata(refcache):
    """Return the signature and the record count recorded in a cached index.

    Either is None if the index lacks that row; an index written by an earlier
    version of biofrills records no signature.
    """
    with closing(sqlite3.connect(refcache)) as con:
        rows = dict(
            con.execute(
                "SELECT key, value FROM meta_data WHERE key IN (?, ?);",
                (_SIGNATURE_KEY, "count"),
            ).fetchall()
        )
    return rows.get(_SIGNATURE_KEY), rows.get("count")


def aa_frequencies(seq, gap_chars="-."):
    """Calculate the amino acid frequencies in a sequence set.

    Raises ValueError if the sequence contains only gap characters.
    """
    aa_counts = Counter(seq)
    # Don't count gaps
    for gap_char in gap_chars:
        if gap_char in aa_counts:
            del aa_counts[gap_char]
    if not aa_counts:
        raise ValueError("No residues to count: every character is a gap")
    # Reduce to frequencies
    scale = 1.0 / sum(aa_counts.values())
    return dict((aa, cnt * scale) for aa, cnt in aa_counts.items())
