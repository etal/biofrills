import logging

import pytest
from alignments import make_alignment, make_records

from biofrills.consensus import consensus


@pytest.mark.parametrize("simple", [False, True])
@pytest.mark.parametrize("trim_ends", [True, False])
def test_consensus_of_identical_gapless_sequences_is_that_sequence(simple, trim_ends):
    aln = make_alignment(["MKVLAYRGHQ"] * 3)
    assert consensus(aln, simple=simple, trim_ends=trim_ends) == "MKVLAYRGHQ"


# The last column has no residue to report.  Before it, column 3 holds a
# single Q, a majority-gap column that trim_ends=False keeps as part of the
# C-terminal tail; the gap of the last column must follow it, not precede it.
@pytest.mark.parametrize("last", ["a", ".", "-"])
@pytest.mark.parametrize("simple", [False, True])
def test_consensus_keeps_column_order_behind_the_c_terminal_tail(last, simple):
    aln = make_alignment(["MKQ" + last] + ["MK-" + last] * 3)
    assert consensus(aln, simple=simple, trim_ends=False) == "MKQ-"
    assert consensus(aln, simple=simple, trim_ends=True) == "MK--"


@pytest.mark.parametrize("rows", [["AC", "CA"], ["CA", "AC"]])
def test_entropy_consensus_breaks_exact_ties_alphabetically(rows):
    # Every column holds A and C once, and so does the whole alignment, so
    # both residues have the same relative-entropy term in both columns.
    assert consensus(make_alignment(rows), weights=[1, 1]) == "AA"


@pytest.mark.parametrize("rows", [["D", "E"], ["E", "D"]])
def test_simple_consensus_breaks_remaining_ties_alphabetically(rows):
    # No preceding column and equal overall frequencies leave only the
    # alphabetical rule.
    assert consensus(make_alignment(rows), simple=True) == "D"


@pytest.mark.parametrize("rows", [["AC", "CA"], ["CA", "AC"]])
def test_simple_consensus_prefers_the_residue_following_the_previous_one(rows):
    # Column 1 goes to A alphabetically; in column 2, C is the residue that
    # follows an A.
    assert consensus(make_alignment(rows), simple=True) == "AC"


@pytest.mark.parametrize("simple", [False, True])
def test_consensus_uppercases_a_mixed_case_column(simple, caplog):
    aln = make_alignment(["AD", "Ad", "AD"])
    with caplog.at_level(logging.WARNING):
        assert consensus(aln, simple=simple) == "AD"
    assert "Mixed lowercase and uppercase letters in a column: DdD" in caplog.text


@pytest.mark.parametrize("simple", [False, True])
@pytest.mark.parametrize("trim_ends", [True, False])
def test_consensus_of_a_single_sequence_is_that_sequence(simple, trim_ends):
    aln = make_alignment(["MKV-LAYRGH"])
    assert consensus(aln, simple=simple, trim_ends=trim_ends) == "MKV-LAYRGH"


def test_entropy_consensus_of_a_residue_found_elsewhere_only_in_lowercase():
    assert consensus(make_alignment(["Ad", "AE", "AE"])) == "AE"


@pytest.mark.parametrize("simple", [False, True])
def test_background_frequencies_ignore_case(simple):
    # Column 1 is all lowercase and yields a gap.  In column 2, D and E tie;
    # counting the lowercase d of column 1, D is the more common residue in
    # the alignment, so its entropy term is smaller and the simple
    # algorithm's second rule prefers the less common E.
    aln = make_alignment(["dD", "dE"])
    assert consensus(aln, weights=[1, 1], simple=simple) == "-E"


@pytest.mark.parametrize("simple", [False, True])
def test_consensus_of_gaps_alone_raises(simple):
    with pytest.raises(ValueError, match="every character is a gap"):
        consensus(make_alignment(["--", "--"]), simple=simple)


@pytest.mark.parametrize("simple", [False, True])
def test_rows_of_different_lengths_raise(simple):
    with pytest.raises(ValueError):
        consensus(make_records(["MKVL", "MKV"]), simple=simple)
