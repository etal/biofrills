import math

import pytest
from alignments import columns, make_alignment, rows
from Bio.Align import MultipleSeqAlignment

from biofrills import alnutils


def test_sequence_weights_sum1_sums_to_one(aligned_input):
    weights = alnutils.sequence_weights(aligned_input, "sum1")
    assert math.fsum(weights) == pytest.approx(1.0, rel=1e-12)


def test_sequence_weights_max1_has_maximum_one(aligned_input):
    weights = alnutils.sequence_weights(aligned_input, "max1")
    assert max(weights) == pytest.approx(1.0, rel=1e-12)


def test_sequence_weights_avg1_has_mean_one(aligned_input):
    weights = alnutils.sequence_weights(aligned_input, "avg1")
    assert math.fsum(weights) / len(weights) == pytest.approx(1.0, rel=1e-12)


def test_sequence_weights_andy_caps_avg1_at_one(aligned_input):
    avg1 = alnutils.sequence_weights(aligned_input, "avg1")
    andy = alnutils.sequence_weights(aligned_input, "andy")
    assert andy == [min(weight, 1.0) for weight in avg1]


@pytest.mark.parametrize("scaling", ["none", "sum1", "max1", "avg1", "andy"])
def test_identical_sequences_receive_equal_weight(aligned_input, scaling):
    first, *rest = aligned_input
    aln = MultipleSeqAlignment([first, *rest, first])
    weights = alnutils.sequence_weights(aln, scaling)
    assert weights[0] == weights[-1]


def test_sequence_weights_rejects_unknown_scaling(aligned_input):
    with pytest.raises(ValueError, match="Unknown scaling"):
        alnutils.sequence_weights(aligned_input, "sum2")


@pytest.mark.parametrize("weighted", [False, True])
@pytest.mark.parametrize("gap_chars", ["-.", "-.X"])
def test_aa_frequencies_sum_to_one_without_gaps(aligned_input, weighted, gap_chars):
    freqs = alnutils.aa_frequencies(aligned_input, weighted or None, gap_chars)
    assert math.fsum(freqs.values()) == pytest.approx(1.0, rel=1e-12)
    assert not set(gap_chars) & freqs.keys()


def test_remove_empty_cols_keeps_every_column_with_a_residue():
    # Columns 1, 7 and 8 are all '-'; only those go.  A column of '-' and
    # '.' is kept, because remove_empty_cols treats only '-' as empty.
    aln = make_alignment(["--A-.C--", "-.-D-E--", "---F-G--"])
    cleaned = list(alnutils.remove_empty_cols(aln))
    assert columns(cleaned) == ["-.-", "A--", "-DF", ".--", "CEG"]
    assert [rec.id for rec in cleaned] == ["s0", "s1", "s2"]


def test_remove_empty_cols_removes_only_all_gap_columns(aligned_input):
    cleaned = list(alnutils.remove_empty_cols(aligned_input))
    assert columns(cleaned) == [
        col for col in columns(aligned_input) if set(col) != {"-"}
    ]


def test_remove_empty_cols_is_idempotent(aligned_input):
    once = list(alnutils.remove_empty_cols(aligned_input))
    twice = list(alnutils.remove_empty_cols(once))
    assert rows(twice) == rows(once)
    assert [(r.id, r.description) for r in twice] == [
        (r.id, r.description) for r in once
    ]


@pytest.mark.parametrize("threshold", [0.4, 0.5, 0.9])
def test_unweighted_blocks_is_idempotent(aligned_input, threshold):
    once = alnutils.blocks(aligned_input, threshold, weights=False)
    twice = alnutils.blocks(once, threshold, weights=False)
    assert rows(twice) == rows(once)


@pytest.mark.parametrize("threshold", [0.4, 0.5, 0.9])
def test_blocks_with_given_weights_is_idempotent(aligned_input, threshold):
    # Passing the weights keeps them fixed; by default, blocks computes them
    # again from its input, which after one pass is a different alignment.
    weights = alnutils.sequence_weights(aligned_input, "avg1")
    once = alnutils.blocks(aligned_input, threshold, weights=weights)
    twice = alnutils.blocks(once, threshold, weights=weights)
    assert rows(twice) == rows(once)


def test_blocks_keeps_columns_at_the_threshold():
    # Unweighted non-gap fractions by column: 1, 3/4, 1/2, 1/4, 0.
    aln = make_alignment(["AAAA-", "CCC--", "DD---", "E----"])
    assert columns(alnutils.blocks(aln, 0.5, weights=False)) == ["ACDE", "ACD-", "AC--"]
