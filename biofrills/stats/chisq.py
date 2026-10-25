"""Python implementation of chisqprob, to avoid SciPy dependency.

Adapted from SciPy: scipy/special/cephes/{chdtr,igam}.c

The results are not those of SciPy's chi2.sf. The series and the continued
fraction stop at a relative change of MACHEP = 1e-7, where Cephes stops at
the machine epsilon, and chisqprob uses the series for x < 1 or x < df, where
Cephes uses it for x < 2 or x < df. The largest relative difference from
scipy.stats.chi2.sf over x from 0.01 to 80 and df from 1 to 30 is 3.6e-7, at
x = 1 and df = 1. Both choices are kept from earlier releases, so that
results do not change.
"""

import math

# Cephes Math Library Release 2.0:  April, 1987
# Copyright 1985, 1987 by Stephen L. Moshier
# Direct inquiries to 30 Frost Street, Cambridge, MA 02140

MACHEP = 0.0000001  # relative change at which the series and fraction stop
BIG = 4.503599627370496e15
BIGINV = 2.22044604925031308085e-16


def chisqprob(x, df):
    """Upper-tail probability of the chi-squared distribution.

    Parameters
    ----------
    x : float
        The chi-squared statistic. For x <= 0 the result is 1.0, whatever
        the value of df.
    df : float
        Degrees of freedom. For x > 0, ValueError is raised unless df is
        positive.

    Returns
    -------
    float
        The area from `x` to infinity under the chi-squared probability
        density with `df` degrees of freedom.
    """
    if x <= 0:
        return 1.0
    if df <= 0:
        raise ValueError("Domain error.")
    if x < 1.0 or x < df:
        return 1.0 - _igam(0.5 * df, 0.5 * x)
    return _igamc(0.5 * df, 0.5 * x)


def _igamc(a, x):
    """Complemented incomplete Gamma integral, by continued fraction.

    The function is defined by::

        igamc(a,x)   =   1 - igam(a,x)

                                inf.
                                   -
                          1       | |  -t  a-1
                    =   -----     |   e   t   dt.
                         -      | |
                        | (a)    -
                                    x

    Both arguments must be positive. chisqprob calls this function when its
    statistic is at least 1 and at least df, and the power series of _igam
    otherwise.
    """
    # Compute  x**a * exp(-x) / Gamma(a)
    ax = math.exp(a * math.log(x) - x - math.lgamma(a))

    # Continued fraction
    y = 1.0 - a
    z = x + y + 1.0
    c = 0.0
    pkm2 = 1.0
    qkm2 = x
    pkm1 = x + 1.0
    qkm1 = z * x
    ans = pkm1 / qkm1
    while True:
        c += 1.0
        y += 1.0
        z += 2.0
        yc = y * c
        pk = pkm1 * z - pkm2 * yc
        qk = qkm1 * z - qkm2 * yc
        if qk != 0:
            r = pk / qk
            t = abs((ans - r) / r)
            ans = r
        else:
            t = 1.0
        pkm2 = pkm1
        pkm1 = pk
        qkm2 = qkm1
        qkm1 = qk
        if abs(pk) > BIG:
            pkm2 *= BIGINV
            pkm1 *= BIGINV
            qkm2 *= BIGINV
            qkm1 *= BIGINV
        if t <= MACHEP:
            return ans * ax


def _igam(a, x):
    """Left tail of incomplete Gamma function, by power series::

            inf.      k
     a  -x   -       x
    x  e     >   ----------
             -     -
            k=0   | (a+k+1)
    """
    # Compute  x**a * exp(-x) / Gamma(a)
    ax = math.exp(a * math.log(x) - x - math.lgamma(a))

    # Power series
    r = a
    c = 1.0
    ans = 1.0

    while True:
        r += 1.0
        c *= x / r
        ans += c
        if c / ans <= MACHEP:
            return ans * ax / a
