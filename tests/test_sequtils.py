import math

import pytest
from alignments import rows

from biofrills import sequtils


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
