import math
import random

import numpy as np
import pytest
from alignments import columns, make_alignment, make_records, rows
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


def test_sequence_weights_andy_halves_avg1_and_truncates_at_one(aligned_input):
    # Correctly rounded arithmetic commutes with scaling by a power of two,
    # so halving the avg1 scale gives the same result bit for bit.
    avg1 = alnutils.sequence_weights(aligned_input, "avg1")
    andy = alnutils.sequence_weights(aligned_input, "andy")
    assert andy == [min(weight / 2, 1.0) for weight in avg1]


AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"


def random_rows(n, length, seed):
    rng = random.Random(seed)
    return ["".join(rng.choices(AMINO_ACIDS, k=length)) for _ in range(n)]


@pytest.mark.parametrize("n", [2, 5, 20])
def test_none_weights_of_identical_residues_sum_to_one(n):
    # Gaps are not residue types, so where they fall does not matter.
    row = random_rows(1, 60, seed=n)[0]
    gapped = [row[:i] + "-" + row[i + 1 :] for i in range(n)]
    for seqs in ([row] * n, gapped):
        weights = alnutils.sequence_weights(make_alignment(seqs), "none")
        assert math.fsum(weights) == pytest.approx(1.0, rel=1e-12)


def test_none_weights_sum_to_the_number_of_independent_groups():
    # 60 rows in k groups, each holding copies of one of k random sequences.
    # The estimate assumes uniform residues, as here; over 40 seeds, its
    # relative error with 200 columns was at most 3.1%.
    sums = []
    for k in (1, 2, 5, 20, 60):
        founders = random_rows(k, 200, seed=k)
        aln = make_alignment([founders[i % k] for i in range(60)])
        sums.append(math.fsum(alnutils.sequence_weights(aln, "none")))
    assert sums == sorted(set(sums))
    assert sums == pytest.approx([1, 2, 5, 20, 60], rel=0.05)


@pytest.mark.parametrize("n", [3, 20])
def test_none_weights_are_capped_at_the_number_of_sequences(n):
    # Every column holds n distinct residues, more than n independent
    # sequences show on average; with n = 20 the mean reaches 20.
    seqs = [(AMINO_ACIDS[i:] + AMINO_ACIDS[:i]) * 3 for i in range(n)]
    weights = alnutils.sequence_weights(make_alignment(seqs), "none")
    assert math.fsum(weights) == pytest.approx(n, rel=1e-12)


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


def test_weights_may_be_a_numpy_array(aligned_input):
    # aa_counts and blocks also accept True, and blocks False, as flags; an
    # array of weights must not be compared with them elementwise.
    weights = alnutils.sequence_weights(aligned_input, "avg1")
    array = np.array(weights)
    assert alnutils.aa_counts(aligned_input, array) == alnutils.aa_counts(
        aligned_input, weights
    )
    assert rows(alnutils.blocks(aligned_input, 0.5, array)) == rows(
        alnutils.blocks(aligned_input, 0.5, weights)
    )


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


@pytest.mark.parametrize(("nrows", "threshold"), [(5, 0.2), (10, 0.1)])
def test_unweighted_blocks_keeps_a_column_of_one_residue_at_its_proportion(
    nrows, threshold
):
    # The non-gap proportion of column 2 is exactly 1/nrows, which equals the
    # threshold, although 1 - (nrows - 1) / nrows computes just below it.
    aln = make_alignment(["AA"] + ["A-"] * (nrows - 1))
    assert len(columns(alnutils.blocks(aln, threshold, weights=False))) == 2


def test_blocks_keeps_a_column_whose_weighted_proportion_equals_the_threshold():
    # The avg1 weights of these four rows are 11/12, 13/12, 13/12 and 11/12
    # in exact arithmetic, so the non-gap proportion of columns 2 and 4 is
    # exactly one half.  As stored, the two weights of 11/12 differ by one
    # unit in the last place, so column 4 computes just below one half.
    weights = [
        float.fromhex("0x1.d555555555556p-1"),
        float.fromhex("0x1.1555555555556p+0"),
        float.fromhex("0x1.1555555555556p+0"),
        float.fromhex("0x1.d555555555557p-1"),
    ]
    aln = make_alignment(["E-CC", "CDA-", "--FD", "EAE-"])
    expected = ["EC-E", "-D-A", "CAFE", "C-D-"]
    assert columns(alnutils.blocks(aln, 0.5, weights)) == expected


@pytest.mark.parametrize("scaling", ["none", "sum1", "max1", "avg1", "andy"])
def test_blocks_keeps_the_same_columns_under_every_weight_scaling(scaling):
    # A weighted proportion is a share of the total weight, whatever the
    # weights sum to.  The unweighted non-gap proportions of these columns are
    # 1, 3/4, 1/2, 1/4 and 0, and no row weighs more than 1.25 times another.
    aln = make_alignment(["AAAA-", "CCC--", "DD---", "E----"])
    weights = alnutils.sequence_weights(aln, scaling)
    assert columns(alnutils.blocks(aln, 0.5, weights)) == ["ACDE", "ACD-", "AC--"]


@pytest.mark.parametrize("heavy_row_first", [True, False])
def test_weighted_blocks_is_independent_of_row_order(heavy_row_first):
    # Column 2 holds one residue of weight 32769 and 2**16 residues of weight
    # 2**-39, under half a unit in the last place of 32769, so its non-gap
    # weight is 32769 + 2**-23, as is the weight of its one gap, and its exact
    # non-gap proportion is 1/2.  Added to 32769 one at a time, each small
    # weight rounds away, which puts the proportion 1.8e-12 below 1/2; added
    # to each other first, they sum exactly.  The weights are a NumPy array
    # because the built-in sum() compensates for rounding error only when
    # every term is exactly of type float, and numpy.float64 is a subclass.
    n_small = 2**16
    heavy = 32769.0
    small_rows = ["MA"] * n_small
    small_weights = [2.0**-39] * n_small
    if heavy_row_first:
        rows = ["MA", *small_rows, "M-"]
        weights = [heavy, *small_weights, heavy + n_small * 2.0**-39]
    else:
        rows = [*small_rows, "MA", "M-"]
        weights = [*small_weights, heavy, heavy + n_small * 2.0**-39]
    aln = make_alignment(rows)
    assert len(columns(alnutils.blocks(aln, 0.5, np.array(weights)))) == 2


# With no column, no per-column comparison of lengths takes place.
@pytest.mark.parametrize("seqs", [["AC-", "A-C"], ["", ""]])
def test_blocks_rejects_weights_of_another_length(seqs):
    with pytest.raises(ValueError, match="Length mismatch"):
        alnutils.blocks(make_alignment(seqs), 0.5, [1.0, 1.0, 1.0])


# No column holds two non-gap residues, so nothing tells the rows apart.
@pytest.mark.parametrize("seqs", [["MKV-LAYRGH"], ["AC--", "--DE"], ["A-", "--"]])
@pytest.mark.parametrize("scaling", ["none", "sum1", "max1", "avg1", "andy"])
def test_uninformative_alignment_gets_uniform_weights(seqs, scaling):
    weights = alnutils.sequence_weights(make_alignment(seqs), scaling)
    uniform = {"sum1": 1 / len(seqs), "andy": 0.5}.get(scaling, 1.0)
    assert weights == [uniform] * len(seqs)


def test_sequence_weights_rejects_an_empty_alignment():
    with pytest.raises(ValueError, match="No sequences"):
        alnutils.sequence_weights([])


def test_col_counts_default_to_unit_weights():
    assert alnutils.col_counts("AAC-.") == {"A": 2.0, "C": 1.0}
    assert alnutils.col_frequencies("AAC-.") == {"A": 2 / 3, "C": 1 / 3}


def test_col_counts_rejects_weights_of_another_length():
    with pytest.raises(ValueError):
        alnutils.col_counts("AAC", [1, 1])


def test_frequencies_of_gaps_alone_raise():
    with pytest.raises(ValueError, match="every character is a gap"):
        alnutils.aa_frequencies(make_alignment(["-.", ".-"]))
    with pytest.raises(ValueError, match="every character is a gap"):
        alnutils.col_frequencies("-.")


@pytest.mark.parametrize("weights", [False, None])
def test_blocks_keeps_every_record_when_no_column_is_kept(weights):
    # Each column is half gaps, below the threshold of 1.
    blocked = alnutils.blocks(make_alignment(["A-", "-C"]), 1.0, weights)
    assert [rec.id for rec in blocked] == ["s0", "s1"]
    assert rows(blocked) == ["", ""]


def test_remove_empty_cols_keeps_every_record_of_gaps_alone():
    cleaned = list(alnutils.remove_empty_cols(make_alignment(["--", "--"])))
    assert [rec.id for rec in cleaned] == ["s0", "s1"]
    assert rows(cleaned) == ["", ""]


@pytest.mark.parametrize(
    "call",
    [
        lambda recs: alnutils.aa_counts(recs, weights=[1, 1]),
        lambda recs: alnutils.blocks(recs, 0.5, weights=False),
        lambda recs: list(alnutils.remove_empty_cols(recs)),
        alnutils.sequence_weights,
    ],
)
def test_rows_of_different_lengths_raise(call):
    with pytest.raises(ValueError):
        call(make_records(["MKVL", "MKV"]))
