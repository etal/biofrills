import numpy as np
import pytest
from scipy.stats import chi2

from biofrills.stats.chisq import chisqprob

# The largest relative difference from scipy.stats.chi2.sf and its two
# sources are given in the docstring of biofrills.stats.chisq.
REL_TOL = 5e-7


@pytest.mark.parametrize("df", range(1, 31))
def test_chisqprob_agrees_with_scipy(df):
    xs = np.arange(0.25, 80.01, 0.25)
    actual = [chisqprob(float(x), df) for x in xs]
    np.testing.assert_allclose(actual, chi2.sf(xs, df), rtol=REL_TOL, atol=0)


@pytest.mark.parametrize("x", [0, -1.5])
@pytest.mark.parametrize("df", [0, 1, 10])
def test_chisqprob_is_one_for_nonpositive_x(x, df):
    assert chisqprob(x, df) == 1.0


@pytest.mark.parametrize("df", [0, -1])
def test_chisqprob_rejects_nonpositive_df(df):
    with pytest.raises(ValueError, match="Domain error"):
        chisqprob(1.0, df)
