import itertools

import pytest
from alignments import rows

from biofrills import pairutils


def test_blosum62_is_symmetric():
    table = pairutils.blosum62
    assert all(table[b, a] == score for (a, b), score in table.items())


def test_blosum62_covers_every_pair_of_its_residue_types():
    letters = {a for a, _ in pairutils.blosum62}
    assert pairutils.blosum62.keys() == set(itertools.product(letters, repeat=2))


@pytest.mark.parametrize("func", ["score_pairwise", "identity_abs", "identity_avg"])
def test_pairwise_measures_are_symmetric(aligned_input, func):
    measure = getattr(pairutils, func)
    seqs = rows(aligned_input)
    for aseq, bseq in itertools.combinations(seqs, 2):
        assert measure(aseq, bseq) == measure(bseq, aseq)


def test_score_pairwise_penalizes_only_internal_gaps():
    # Matches score by BLOSUM62 (A-A 4, C-C 9, D-D 6).  The internal gap of
    # two columns costs one opening (-10) and two extensions (-0.5 each);
    # terminal gaps cost their extensions only.
    assert pairutils.score_pairwise("AC--D", "ACEED") == 4 + 9 - 10 - 1 + 6
    assert pairutils.score_pairwise("--ACD", "EEACD") == -1 + 4 + 9 + 6
    assert pairutils.score_pairwise("ACD--", "ACDEE") == 4 + 9 + 6 - 1


def test_identity_avg_of_two_rows_without_residues_is_zero():
    assert pairutils.identity_avg("--.--", "..-..") == 0.0


def test_identity_ignores_letter_case():
    assert pairutils.identity_abs("ACDE", "acde") == 4
    assert pairutils.identity_abs("ACDE", "acdf") == 3
    assert pairutils.identity_avg("ACDE", "acd-") == 3 / 3.5


def test_identity_of_a_sequence_with_itself_is_its_residue_count():
    """The identity of a row with itself equals its ungap_len.

    seq-pick -k charges a medoid its ungap_len and every other sequence its
    identity_abs to the medoid; their agreement means that a medoid leaves
    none of its own residues unmatched, its distance to itself being zero.
    """
    seq = "MKa-VLw.GHX"
    assert pairutils.identity_abs(seq, seq) == pairutils.ungap_len(seq) == 9


def test_no_row_matches_another_at_more_sites_than_it_has_residues(aligned_input):
    """seq-pick -k relies on this bound to keep each medoid its own closest."""
    seqs = rows(aligned_input)
    for aseq, bseq in itertools.permutations(seqs, 2):
        assert pairutils.identity_abs(aseq, bseq) <= pairutils.ungap_len(aseq)
