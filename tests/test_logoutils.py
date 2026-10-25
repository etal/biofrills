from collections import Counter

from alignments import columns, make_alignment

from biofrills import logoutils


def test_aln2logodata_counts_the_residues_of_each_column():
    aln = make_alignment(["ACD-", "ACE-", "AYE-"])
    logodata = logoutils.aln2logodata(aln)
    assert [posn for posn, _, _, _ in logodata] == [1, 2, 3, 4]
    for (_, counts, _, _), col in zip(logodata, columns(aln), strict=True):
        nonzero = {letter: count for letter, count in counts.items() if count}
        assert nonzero == Counter(col.replace("-", ""))
