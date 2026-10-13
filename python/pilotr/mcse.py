"""Monte Carlo uncertainty for the rates the power functions report, mirroring mcse.R.

Every power and decision probability is a proportion over a finite number of replicates, so it
carries its own uncertainty, and a bare number invites the reader to treat it as exact. At 60
replicates, a rate near 0.5 has a Monte Carlo standard error of 0.065, so a reported power of
0.80 is compatible with anything from roughly 0.70 to 0.90.

Each function performs the operations of its twin in the R package's mcse.R in the same order,
so that for the same counts both return the same doubles. python/tests/test_mcse.py and the R
package's test-mcse.R pin them. Sums are therefore written out as loops, as they are in R.
``statistics`` works in exact fractions, and from Python 3.12 ``sum()`` adds floats with
compensated summation, so either would round differently from the plain double folds of R.
"""

from __future__ import annotations
import math

_Z95 = 1.959963984540054  # qnorm(0.975), written out so the value cannot drift


def _mcse(p, n):
    """The plain Monte Carlo standard error of a proportion `p` over `n` replicates."""
    return math.nan if n <= 0 else math.sqrt(p * (1 - p) / n)


def _wilson(p, n, z=_Z95):
    """Wilson score interval for a proportion `p` over `n` replicates, as ``(lo, hi)``.

    The plain standard error collapses to zero at p = 0 and p = 1, which reads as certainty just
    where there is least of it: nought significant results out of 60 is not evidence that power
    is exactly zero. The Wilson interval stays inside [0, 1] and keeps a sensible width at the
    boundaries, and Brown, Cai and DasGupta (2001, doi:10.1214/ss/1009213286) recommend it in
    place of the Wald interval.
    """
    if n <= 0:
        return (math.nan, math.nan)
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    lo = centre - half
    hi = centre + half
    # At p = 0 and p = 1, the two terms cancel exactly in algebra but not in floating point, which
    # leaves a bound like 2.8e-17 where the answer is plainly 0. The endpoints are exact there, so
    # they are set by hand.
    if p <= 0:
        lo = 0.0
    if p >= 1:
        hi = 1.0
    return (max(0.0, lo), min(1.0, hi))


def _rate_with_error(count, n, prefix):
    """A rate with its Monte Carlo standard error and Wilson interval, as a flat dict.

    The keys are `prefix` and `prefix` followed by ``_mcse``, ``_lo`` and ``_hi``, ready to be
    spliced into a returned dict. With no replicates, all four are nan.
    """
    p = count / n if n > 0 else math.nan
    lo, hi = _wilson(0.0 if math.isnan(p) else p, n)
    return {prefix: p, prefix + "_mcse": _mcse(p, n), prefix + "_lo": lo, prefix + "_hi": hi}


def _type_m_mcse(ratios):
    """The Monte Carlo standard error of Type M over the significant replicates' `ratios`.

    The standard deviation of the exaggeration ratios, divided by the square root of their
    number. Type M averages over the significant replicates alone, often a handful, so it is the
    least settled number a design analysis reports. A standard deviation needs two values, so
    fewer give nan.
    """
    n = len(ratios)
    if n < 2:
        return math.nan
    s = 0.0
    for r in ratios:
        s = s + r
    m = s / n
    ss = 0.0
    for r in ratios:
        ss = ss + (r - m) * (r - m)
    return math.sqrt(ss / (n - 1) / n)


def _type_s_interval(n_wrong, n_sig, true_effect):
    """The Wilson interval for Type S, `n_wrong` wrong signs among `n_sig` significant replicates.

    Its plain standard error is zero whenever no sign is wrong. That is the usual result, and the
    standard error would then claim certainty from what may be a handful of replicates. Type S is
    defined only against a non-zero true effect and over at least one significant replicate, so
    otherwise both bounds are nan.
    """
    if math.isnan(true_effect) or true_effect == 0 or n_sig <= 0:
        return (math.nan, math.nan)
    return _wilson(n_wrong / n_sig, n_sig)
