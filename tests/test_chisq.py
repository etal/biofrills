import numpy as np
import pytest
from scipy.stats import chi2

from biofrills.stats.chisq import chisqprob

# chisqprob stops its series and continued fraction at a relative change of
# MACHEP = 1e-7, which bounds its accuracy.  The largest relative error
# against scipy.stats.chi2.sf over x from 0.01 to 80 in steps of 0.01 and
# df from 1 to 30 is 3.6e-7, at x = 1 and df = 1.
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
