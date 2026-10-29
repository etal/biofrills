import gc
import logging
import math
import os
import sqlite3
from contextlib import closing

import pytest
from alignments import rows

from biofrills import sequtils


def intersect(tmp_path, keys, **kwargs):
    """Return the identifiers that intersect_keys selects from ref.fa."""
    ref = tmp_path / "ref.fa"
    return [rec.id for rec in sequtils.intersect_keys(keys, ref, **kwargs)]


def test_intersect_keys_reuses_a_matching_cache(tmp_path, caplog):
    (tmp_path / "ref.fa").write_text(">a\nACDE\n>b\nFGHI\n")
    with caplog.at_level(logging.WARNING):
        assert intersect(tmp_path, ["b", "z"], cache=True) == ["b"]
        inode = os.stat(tmp_path / "ref.fa.sqlite").st_ino
        assert intersect(tmp_path, ["b", "z"], cache=True) == ["b"]
    assert not caplog.records
    assert os.stat(tmp_path / "ref.fa.sqlite").st_ino == inode
    assert sorted(os.listdir(tmp_path)) == ["ref.fa", "ref.fa.sqlite"]


def test_intersect_keys_caches_a_reference_given_as_a_string(tmp_path):
    """seq-intersect passes the reference path as a string."""
    (tmp_path / "ref.fa").write_text(">a\nACDE\n")
    ref = str(tmp_path / "ref.fa")
    assert [rec.id for rec in sequtils.intersect_keys(["a"], ref, cache=True)] == ["a"]
    inode = os.stat(tmp_path / "ref.fa.sqlite").st_ino
    assert [rec.id for rec in sequtils.intersect_keys(["a"], ref, cache=True)] == ["a"]
    assert os.stat(tmp_path / "ref.fa.sqlite").st_ino == inode


def test_intersect_keys_rebuilds_a_corrupt_cache(tmp_path):
    (tmp_path / "ref.fa").write_text(">a\nACDE\n>b\nFGHI\n")
    (tmp_path / "ref.fa.sqlite").write_bytes(b"not a database\n")
    assert intersect(tmp_path, ["b"], cache=True) == ["b"]
    assert intersect(tmp_path, ["a"], cache=True) == ["a"]
    assert sorted(os.listdir(tmp_path)) == ["ref.fa", "ref.fa.sqlite"]


def test_intersect_keys_rebuilds_a_cache_without_a_signature(tmp_path):
    """A cache written by an earlier version records no signature."""
    (tmp_path / "ref.fa").write_text(">a\nACDE\n")
    intersect(tmp_path, ["a"], cache=True)
    cache = tmp_path / "ref.fa.sqlite"
    with closing(sqlite3.connect(cache)) as con, con:
        con.execute("DELETE FROM meta_data WHERE key = ?;", (sequtils._SIGNATURE_KEY,))
    inode = os.stat(cache).st_ino
    assert intersect(tmp_path, ["a"], cache=True) == ["a"]
    # Replaced by a rebuild, although the old index was valid
    assert os.stat(cache).st_ino != inode


@pytest.mark.parametrize(
    "new_reference, mtime_shift",
    [
        # The same size, modified later
        (">a\nACDE\n>c\nFGHI\n", 1),
        # Larger, with the modification time unchanged, as a write within the
        # timestamp resolution of the file system can leave it
        (">a\nACDE\n>b\nFGHI\n>c\nKLMN\n", 0),
    ],
)
def test_intersect_keys_rebuilds_an_outdated_cache(
    tmp_path, new_reference, mtime_shift
):
    ref = tmp_path / "ref.fa"
    ref.write_text(">a\nACDE\n>b\nFGHI\n")
    assert intersect(tmp_path, ["c"], cache=True) == []
    mtime_ns = os.stat(ref).st_mtime_ns + mtime_shift * 10**9
    ref.write_text(new_reference)
    os.utime(ref, ns=(mtime_ns, mtime_ns))
    assert intersect(tmp_path, ["c"], cache=True) == ["c"]


@pytest.mark.parametrize("clean_accs", [False, True])
def test_intersect_keys_rebuilds_a_cache_built_with_the_other_clean_accs(
    tmp_path, clean_accs
):
    (tmp_path / "ref.fa").write_text(">x/1-5\nACDE\n")
    intersect(tmp_path, ["x/1-5"], cache=True, clean_accs=not clean_accs)
    assert intersect(tmp_path, ["x/1-5"], cache=True, clean_accs=clean_accs) == [
        "x/1-5"
    ]


def test_intersect_keys_leaves_no_cache_when_indexing_fails(tmp_path):
    (tmp_path / "ref.fa").write_text(">a\nACDE\n>a\nFGHI\n")
    with pytest.raises(ValueError, match="Duplicate key"):
        intersect(tmp_path, ["a"], cache=True)
    assert os.listdir(tmp_path) == ["ref.fa"]


def test_intersect_keys_closes_the_new_index_when_replacing_the_cache_fails(
    tmp_path, monkeypatch
):
    def fail(src, dst):
        raise OSError("replace failed")

    (tmp_path / "ref.fa").write_text(">a\nACDE\n")
    monkeypatch.setattr(os, "replace", fail)
    with pytest.raises(OSError, match="replace failed"):
        intersect(tmp_path, ["a"], cache=True)
    # The index holds the reference file open; an unclosed handle would be
    # reported here as a ResourceWarning
    gc.collect()
    assert os.listdir(tmp_path) == ["ref.fa"]


def test_aa_frequencies_of_a_sequence_sum_to_one_without_gaps(aligned_input):
    for seq in rows(aligned_input):
        freqs = sequtils.aa_frequencies(seq)
        assert math.fsum(freqs.values()) == pytest.approx(1.0, rel=1e-12)
        assert not {"-", "."} & freqs.keys()


def test_clean_accession_strips_the_hmmer_domain_location():
    assert sequtils.clean_accession("sp|P04637|P53_HUMAN/94-292") == (
        "sp|P04637|P53_HUMAN"
    )
    assert sequtils.clean_accession("P53_HUMAN") == "P53_HUMAN"


def test_aa_frequencies_of_gaps_alone_raise():
    with pytest.raises(ValueError, match="every character is a gap"):
        sequtils.aa_frequencies("-.-")
