"""Solving a simulated design curve for the swept value that meets a target.

The curve functions report a decision rate at each swept value and stop there, leaving the
reader to judge the crossing from a plot. That judgement is made over points whose Monte
Carlo intervals overlap, and what it yields is a bare number, the one that goes into a
preregistration, carrying no record of how wide the range of values compatible with the
simulation was. The functions here estimate the crossing instead, and report an interval
with it.

The model is a binomial regression of the rate on the swept value, weighted by the replicate
count behind each point, with a probit link. The link is not a convenience. A power curve is a
normal tail probability, and under the normal approximation to a two-group comparison the
probit of power is exactly linear in the square root of the sample size, so this
parameterisation estimates two coefficients of a curve the design analysis already implies
rather than bending a general sigmoid to fit.

The choice was checked rather than assumed, by ``tools/calibration/solve_curve_calibration.R``,
which produces every figure quoted here and writes them to a file beside itself. Against R's
``stats::power.t.test``, over twelve combinations of effect size and target power on each of three
grid shapes at 400 replicates a point, the probit solved to a mean absolute error of 2.75%, 2.18%
and 3.70%, the logit to 3.27%, 2.54% and 4.28%. The two separate most clearly on the grid that
sweeps the widest range: reaching a power of 0.995, the logit's Pearson chi-square per degree of
freedom rose to 1.65 against 1.15 for the probit, and its intervals covered the analytic answer
nine times in twelve where the probit covered twelve. The margin is not large, because at this
replicate count most of what separates a solved size from the analytic one is Monte Carlo noise in
the curve rather than the link fitted to it. Rerun at the seed base 77000, which the calibration
script takes as its third argument, and the probit still leads on the tight grid and by more on
the tall one, while the coarse grid swaps: four widely spaced points do not tell the two links
apart.

The inversion and its interval are the delta method that R's ``MASS::dose.p`` applies to a
fitted glm: the solved point on the fitted scale is ``(link(target) - intercept) / slope``, and
its standard error follows from the gradient of that expression in the two coefficients.

Fieller's exact interval for a ratio (Fieller, 1954) is the obvious alternative, and it was
tried. It does not separate from the delta method on this evidence. Over the same 36 checks it
covered the analytic answer 36 times against the delta method's 35, a difference of one check, at
a mean width of 9.48 subjects against 9.38, and it was bounded every time. These curves determine
the slope well enough that the two intervals nearly coincide, so nothing here argues for one over
the other on coverage. The delta method is kept because it is the interval R's ``MASS::dose.p``
reports, which gives the R test suite an independent implementation to check against, and because
the one condition under which Fieller's would differ, a slope too poorly determined for a bounded
interval, is the condition the slope refusal below already rejects outright. A one-check
difference in coverage is noise at this replicate count, and the rerun at seed base 77000 above
puts the two level at 34 each. Anyone reopening the question should rerun the calibration before
arguing from it.

One departure from ``dose.p``'s default: where the two-parameter model does not describe the
curve, the covariance is scaled by Pearson's chi-square over its degrees of freedom, the
heterogeneity factor of probit analysis (Finney, 1971). The factor is floored at one, so it only
ever widens the interval. Without it, a curve the model fits badly still reports the narrow
interval its replicate counts alone imply, which is the overconfidence this function exists to
remove.

The regression is written out rather than handed to a fitting library, because the R twin has
to reproduce it and because the package's core carries no numerical dependency. With an
intercept and one slope the weighted normal equations are a two-by-two solve in closed form,
small enough to mirror line for line, and the reductions are spelled out as loops for the
reason core.py gives: base R's ``sum`` and CPython's ``sum`` are each more accurate than a
plain double fold, and in different ways.

Extrapolation is the one thing these functions must not do. A curve that does not reach the
target within the range it swept is refused rather than extended, and a fit that solves
outside that range is refused as well.

Agreement with the R twin is checked by ``tools/parity/solve_cross.py``, at a relative
tolerance rather than bit for bit. The fit calls ``exp`` and the normal distribution function
at every iteration, and IEEE-754 fixes the rounding of neither, so the two languages' maths
libraries put the last bits in different places even though the arithmetic between them is
written to be identical.
"""

from __future__ import annotations

import math

from .core import as241

# Rate columns and replicate-count columns are looked for in this order. The names are those the
# curve functions already use, so a curve goes in unaltered.
_Y_COLUMNS = ("power", "p_meaningful")
_N_COLUMNS = ("n_returned", "n_converged", "n_sims")
_EFFECT_COLUMNS = ("effect", "param")
_TRANSFORMS = ("sqrt", "identity", "log")
_MAXIT = 100
_TOL = 1e-11
# 1 / sqrt(2 * pi). The normal density is written out with it rather than taken from a library,
# so that the R twin evaluates the same expression and the only difference left between the two
# is the rounding of exp().
_INV_SQRT_2PI = 0.3989422804014327


def _phi(x: float) -> float:
    """Standard normal CDF, the inverse of the probit link. R's twin calls stats::pnorm(); the
    two agree to a relative difference of about 1e-16 over the range a fit visits, which is why
    solve_curve() is held to a relative tolerance rather than to the bit-for-bit standard the
    generative core meets."""
    return 0.5 * math.erfc(-x / math.sqrt(2.0))


def _num(v: float) -> str:
    """Every number reaching a refusal message goes through this, so that the two engines
    format it the same way. Both call the C library's %g on a double."""
    return "%g" % v


def _forward(x: float, transform: str) -> float:
    """The scale the fit is linear on."""
    if transform == "sqrt":
        return math.sqrt(x)
    if transform == "log":
        return math.log(x)
    return x


def _back(u: float, transform: str) -> float:
    """The inverse of `_forward`, applied to the solved point and to both interval bounds, so it
    has to be monotone across everything the fit can return; under "sqrt" that means flooring the
    fitted scale at zero before squaring, since a negative square root is not a swept value the
    curve ever visited."""
    if transform == "sqrt":
        u = max(u, 0.0)
        return u * u
    if transform == "log":
        return math.exp(u)
    return u


def _wls(u, z, w):
    """Weighted least squares of `z` on `(1, u)`, with the unscaled covariance of the two
    coefficients. The binomial dispersion is one, so this covariance is what the delta method
    needs before any heterogeneity factor is applied to it. Returns None when the design is
    singular."""
    s0 = s1 = s2 = t0 = t1 = 0.0
    for i in range(len(u)):
        wi = w[i]
        s0 = s0 + wi
        s1 = s1 + wi * u[i]
        s2 = s2 + wi * u[i] * u[i]
        t0 = t0 + wi * z[i]
        t1 = t1 + wi * u[i] * z[i]
    det = s0 * s2 - s1 * s1
    if not math.isfinite(det) or det <= 0:
        return None
    return {"b0": (s2 * t0 - s1 * t1) / det, "b1": (s0 * t1 - s1 * t0) / det,
            "v00": s2 / det, "v01": -s1 / det, "v11": s0 / det}


def _irls(u, y, m):
    """Iteratively reweighted least squares for the binomial probit model with prior weights `m`,
    the number of replicates behind each rate. The starting rates are shrunk towards a half by
    half a replicate, which is what R's `glm` does, so that a rate of exactly zero or one does not
    start the iteration at an infinity. Returns None when the fit is singular or does not
    settle."""
    p0 = [(m[i] * y[i] + 0.5) / (m[i] + 1) for i in range(len(y))]
    eta = [as241(p) for p in p0]
    prev = None
    for _ in range(_MAXIT):
        mu = [_phi(e) for e in eta]
        dmu = [math.exp(-0.5 * e * e) * _INV_SQRT_2PI for e in eta]
        v = [p * (1 - p) for p in mu]
        # A fitted rate pinned at zero or one has neither a variance to divide by nor a
        # derivative, so the point carries no weight at that step rather than an infinity.
        w = [m[i] * dmu[i] * dmu[i] / v[i] if v[i] > 0 and dmu[i] > 0 else 0.0
             for i in range(len(mu))]
        z = [eta[i] + (y[i] - mu[i]) / dmu[i] if dmu[i] > 0 else eta[i] for i in range(len(mu))]
        fit = _wls(u, z, w)
        if fit is None:
            return None
        b = (fit["b0"], fit["b1"])
        eta = [fit["b0"] + fit["b1"] * ui for ui in u]
        if prev is not None and max(abs(b[0] - prev[0]), abs(b[1] - prev[1])) <= \
                _TOL * (1 + max(abs(b[0]), abs(b[1]))):
            fit["dispersion"] = _dispersion(eta, y, m)
            return fit
        prev = b
    return None


def _dispersion(eta, y, m):
    """Pearson's chi-square over the residual degrees of freedom, floored at one: the
    heterogeneity factor of probit analysis (Finney, 1971). A curve the two-parameter model does
    not describe leaves residuals larger than the replicate counts alone account for, and scaling
    the covariance by this factor widens the interval to say so. It never narrows one, because a
    factor below one would claim more precision than the replicates behind the curve support."""
    df = len(y) - 2
    if df < 1:
        return 1.0
    mu = [_phi(e) for e in eta]
    chi = 0.0
    for i in range(len(y)):
        v = mu[i] * (1 - mu[i])
        if v > 0:
            chi = chi + m[i] * (y[i] - mu[i]) * (y[i] - mu[i]) / v
    return max(chi / df, 1.0)


def _columns(curve):
    """The column names of a curve, taken from its first row."""
    return list(curve[0].keys())


def _column(curve, given, auto, label, arg):
    """Resolve a column named by the caller, or the first of `auto` the curve carries."""
    cols = _columns(curve)
    if given is not None:
        if given not in cols:
            raise ValueError("`curve` has no column named '%s'. Its columns are: %s."
                             % (given, ", ".join(cols)))
        return given
    hit = [c for c in auto if c in cols]
    if not hit:
        raise ValueError(
            "`curve` has no %s column. Name one with `%s`; the columns recognised "
            "automatically are %s." % (label, arg, ", ".join("'%s'" % c for c in auto)))
    return hit[0]


def _probability(v, arg):
    """A single probability strictly inside the unit interval. The two checks are separate so
    that a missing or non-numeric argument is never formatted into the message, which R would
    render as "NA" and Python as "nan"."""
    ok = isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v)
    if not ok:
        raise ValueError("`%s` must be a single number strictly between 0 and 1." % arg)
    if v <= 0 or v >= 1:
        raise ValueError("`%s` must be strictly between 0 and 1; got %s." % (arg, _num(v)))


def _rows(curve, effect):
    """The rows to fit, after selecting a focal effect if one was named."""
    if effect is None:
        return list(range(len(curve)))
    cols = _columns(curve)
    col = [c for c in _EFFECT_COLUMNS if c in cols]
    if not col:
        raise ValueError("`curve` has no 'effect' or 'param' column to select a focal effect from.")
    have = [str(r[col[0]]) for r in curve]
    keep = [i for i, v in enumerate(have) if v == str(effect)]
    if not keep:
        seen = list(dict.fromkeys(have))
        raise ValueError("`curve` holds no focal effect named '%s'. It holds: %s."
                         % (effect, ", ".join(seen)))
    return keep


def _numeric(curve, name, keep, what):
    """One column as a list of floats, refusing a column that is not numeric. A missing entry
    becomes a nan rather than a refusal, because that is what it is: the R twin reads the same
    absence as NA and drops the point for want of a rate, and a curve carried across the two
    languages as JSON arrives with its NAs written as nulls. What the refusal is for is a
    column of the wrong kind, a label or a string, which no amount of dropping would rescue."""
    out = []
    for i in keep:
        v = curve[i][name]
        if v is None:
            out.append(float("nan"))
            continue
        if isinstance(v, bool) or not isinstance(v, (int, float)):
            raise ValueError("the %s column '%s' is not numeric." % (what, name))
        out.append(float(v))
    return out


def _weights(curve, n, keep):
    """The replicate count behind each rate: a column name, a numeric value repeated over the
    rows, or the first count column the curve carries."""
    if isinstance(n, (int, float)) and not isinstance(n, bool):
        return [float(n)] * len(keep)
    if isinstance(n, (list, tuple)):
        if not len(n):
            raise ValueError("`n` must hold at least one replicate count.")
        return [float(n[i % len(n)]) for i in range(len(keep))]
    name = _column(curve, n, _N_COLUMNS, "replicate-count", "n")
    return _numeric(curve, name, keep, "replicate-count")


def solve_curve(curve, target, x=None, y=None, n=None, effect=None,
                transform="sqrt", level=0.95):
    """Solve a simulated design curve for the value that meets a target.

    Take the curve a sweep has already produced, fit the decision rate against the swept
    value, and solve for the value at which the rate meets a target. The solved value comes
    with a confidence interval, because a point read off a simulated curve without one
    repeats the overconfidence that design analysis exists to expose.

    The input is the list of records `power_curve` returns, used as it stands. The swept
    value is taken from the leading column, which is where the curve functions put it, and
    the rate from `power` or `p_meaningful`, whichever the curve carries. Each rate is a
    proportion over a known number of replicates, and that count, read from `n_returned`,
    `n_converged` or `n_sims`, weights the fit: a rate over 200 replicates should count for
    more than a rate over 20.

    The fit is a binomial regression with a probit link, and the solved value is the swept
    value at which the fitted rate equals `target`. The probit is chosen because power is a
    normal tail probability: under the normal approximation to a two-group comparison, the
    probit of power is linear in the square root of the sample size, so the model has the
    shape a design analysis already implies. Measured against R's `stats::power.t.test` across
    twelve combinations of effect size and target power on each of three grid shapes, at 400
    replicates a point, the solved sample size fell within 2.9% of the analytic answer on
    average, against 3.4% for a logit fitted the same way.

    The interval is the delta-method interval of R's `MASS::dose.p`, computed on the scale
    named by `transform` and mapped back, so it is symmetric on that scale rather than on the
    natural one. That asymmetry is the honest shape: at the top of a power curve a given
    change in rate costs far more sample size than the same change lower down. Where the
    two-parameter model does not describe the curve, the interval is widened by the
    heterogeneity factor of probit analysis, Pearson's chi-square over its degrees of freedom
    (Finney, 1971), reported as `dispersion`. It is floored at 1, so a well-fitting curve is
    left alone and a badly-fitting one cannot report a narrower interval than its own
    residuals justify.

    The default `transform` of `"sqrt"` suits a sample-size axis, where a rate rises with the
    square root of the sample size rather than with the sample size itself. Sweep something
    else, an effect size or a random-effect standard deviation, and `"identity"` is usually
    right.

    Nothing here extrapolates. A curve whose rates do not straddle the target is refused,
    with the range it did cover reported, and so is a fit that solves outside the swept
    range. A curve whose fitted slope cannot be told from zero is refused too: the crossing
    is then compatible with any value at all, and an interval that said otherwise would be
    false.

    What the interval covers is the Monte Carlo uncertainty of the fit, not the gap between
    the fitted shape and the true curve. Across 36 checks against R's `stats::power.t.test` at
    400 replicates a point, the solved size sat within 2.9% of the analytic answer on average
    and within 6.9% at worst, and the nominal 95% interval covered the analytic value 35 times
    out of 36. Nearly all of that error is the Monte Carlo noise the interval is describing,
    and it falls with the square root of the replicate count. Raise the count far enough and
    the interval narrows onto a fitted shape that is still slightly the wrong shape, so
    replicates alone do not make a solved size arbitrarily accurate. The remedies are a finer
    grid, more replicates, or a design with a closed form to check against.

    Parameters
    ----------
    curve : sequence of dict
        A curve, as returned by `power_curve`: one record per swept value, holding the swept
        value, a decision rate, and the number of replicates behind it.
    target : float
        The decision rate to solve for, strictly between 0 and 1.
    x : str, optional
        Name of the column holding the swept value. `None`, the default, takes the leading
        column.
    y : str, optional
        Name of the column holding the decision rate. `None`, the default, takes `power` or
        `p_meaningful`, whichever is present.
    n : str or float or sequence, optional
        The number of replicates behind each rate, either the name of a column or a numeric
        value. `None`, the default, takes `n_returned`, `n_converged` or `n_sims`, whichever
        is present.
    effect : str, optional
        Which focal effect to solve for, when the curve holds more than one. Matched against
        the `effect` or `param` column. `None`, the default, uses every record, which is
        correct only when the curve holds one effect.
    transform : str, optional
        The scale the swept value is fitted on: `"sqrt"` (the default, for a sample size),
        `"identity"` or `"log"`.
    level : float, optional
        Confidence level for the reported interval (default 0.95).

    Returns
    -------
    dict
        Keys: `value` (the solved swept value), `lo` and `hi` (its confidence bounds),
        `level`, `target`, `se` (the delta-method standard error on the fitted scale, the
        scale on which the interval is symmetric), `dispersion` (the heterogeneity factor
        applied, 1 where the model fits), `x` and `y` (the columns used), `transform`,
        `intercept` and `slope` (the fitted coefficients), `n_points` (the number of curve
        points the fit used), and `x_min` and `x_max` (the swept range). A bound falling
        outside that range is not an error but a message: the sweep was too narrow to pin the
        value down, and should be widened. A `dispersion` well above 1 says the curve is not
        the shape the model assumes, so the solve deserves a wider grid or more replicates
        rather than trust.

    Raises
    ------
    ValueError
        If `target` or `level` is not strictly between 0 and 1, if the curve carries fewer
        than three usable points, if its rate does not vary, if its rates do not straddle
        `target` within the swept range, if the fitted slope cannot be told from zero, or if
        the solve lands outside the swept range.

    References
    ----------
    Fieller, E. C. (1954). Some problems in interval estimation. *Journal of the Royal
    Statistical Society: Series B*, 16(2), 175-185. doi:10.1111/j.2517-6161.1954.tb00159.x

    Finney, D. J. (1971). *Probit analysis* (3rd ed.). Cambridge University Press.

    Examples
    --------
    >>> curve = [{"n_subject": 40, "power": 0.38, "n_sims": 100},
    ...          {"n_subject": 100, "power": 0.74, "n_sims": 100},
    ...          {"n_subject": 160, "power": 0.89, "n_sims": 100}]
    >>> round(solve_curve(curve, target=0.8)["value"])
    120
    """
    if not isinstance(curve, (list, tuple)) or not len(curve) or not _columns(curve):
        raise ValueError("`curve` must be a table of curve points, with one row per swept value.")
    _probability(target, "target")
    _probability(level, "level")
    if not isinstance(transform, str) or transform not in _TRANSFORMS:
        raise ValueError("`transform` must be one of %s."
                         % ", ".join("'%s'" % t for t in _TRANSFORMS))

    xname = _columns(curve)[0] if x is None else _column(curve, x, (), "swept-value", "x")
    yname = _column(curve, y, _Y_COLUMNS, "decision-rate", "y")

    keep = _rows(curve, effect)
    xv = _numeric(curve, xname, keep, "swept-value")
    yv = _numeric(curve, yname, keep, "decision-rate")
    mv = _weights(curve, n, keep)

    ok = [i for i in range(len(keep))
          if math.isfinite(xv[i]) and math.isfinite(yv[i]) and math.isfinite(mv[i]) and mv[i] > 0]
    xv = [xv[i] for i in ok]
    yv = [yv[i] for i in ok]
    mv = [mv[i] for i in ok]
    k = len(xv)
    if k < 3:
        raise ValueError("solving a curve needs at least 3 swept values with a rate and a "
                         "replicate count; this curve has %d." % k)
    if len(set(xv)) != k:
        raise ValueError("the curve has more than one row at the same swept value. Select a "
                         "single focal effect with `effect`, or subset the curve before solving.")
    bad = [v for v in yv if v < 0 or v > 1]
    if bad:
        raise ValueError("the column '%s' holds a value of %s, which is not a probability; "
                         "solve_curve() inverts a rate between 0 and 1." % (yname, _num(bad[0])))
    if max(yv) == min(yv):
        raise ValueError("the rate is %s at every swept value, so the curve has no trend to "
                         "invert." % _num(yv[0]))
    if transform == "sqrt" and min(xv) < 0:
        raise ValueError("the 'sqrt' transform is not defined for a swept value of %s. Use "
                         'transform = "identity" for an axis that is not a sample size.'
                         % _num(min(xv)))
    if transform == "log" and min(xv) <= 0:
        raise ValueError("the 'log' transform is not defined for a swept value of %s. Use "
                         'transform = "identity" for an axis that is not a sample size.'
                         % _num(min(xv)))
    if target < min(yv) or target > max(yv):
        raise ValueError(
            "the curve does not reach a %s of %s within its swept range. Over swept values from "
            "%s to %s the rate runs from %s to %s, and solving would extrapolate beyond the "
            "values simulated." % (yname, _num(target), _num(min(xv)), _num(max(xv)),
                                   _num(min(yv)), _num(max(yv))))

    u = [_forward(v, transform) for v in xv]
    fit = _irls(u, yv, mv)
    if fit is None:
        raise ValueError("the weighted fit to this curve did not settle within %d iterations, "
                         "so the curve cannot be inverted." % _MAXIT)

    z = as241((1 + level) / 2)
    se_slope = math.sqrt(fit["dispersion"] * fit["v11"])
    if not math.isfinite(se_slope) or abs(fit["b1"]) <= z * se_slope:
        raise ValueError(
            "the fitted slope (%s, standard error %s) cannot be told from zero at the %s level, "
            "so the curve does not determine a value; every swept value is compatible with the "
            "target." % (_num(fit["b1"]), _num(se_slope), _num(level)))

    et = as241(target)
    us = (et - fit["b0"]) / fit["b1"]
    g0 = -1 / fit["b1"]
    g1 = -us / fit["b1"]
    se = math.sqrt(fit["dispersion"]
                   * (g0 * g0 * fit["v00"] + 2 * g0 * g1 * fit["v01"] + g1 * g1 * fit["v11"]))
    value = _back(us, transform)
    if value < min(xv) or value > max(xv):
        raise ValueError("the solved value %s falls outside the swept range %s to %s, so "
                         "reporting it would extrapolate beyond the values simulated."
                         % (_num(value), _num(min(xv)), _num(max(xv))))

    return {"value": value, "lo": _back(us - z * se, transform),
            "hi": _back(us + z * se, transform),
            "level": level, "target": target, "se": se, "dispersion": fit["dispersion"],
            "x": xname, "y": yname, "transform": transform,
            "intercept": fit["b0"], "slope": fit["b1"], "n_points": k,
            "x_min": min(xv), "x_max": max(xv)}


def target_n(curve, target=0.8, **kwargs):
    """Solve a power curve for the sample size that reaches a target power.

    The sample-size case of `solve_curve`, and the number a power analysis is usually run to
    obtain. Takes the curve a sweep over sample size has produced and returns the size at
    which power reaches `target`, rounded up to a whole number of units alongside the exact
    solution.

    Everything `solve_curve` does applies here, including its refusals: a curve that never
    reaches the target within the sizes it swept is refused rather than extrapolated, and the
    reported interval can extend past the largest size simulated, which means the sweep was
    too narrow to settle the question.

    The whole-number fields round up rather than to nearest, because a design cannot recruit
    a fraction of a subject and rounding down would leave the study short of the target it
    was sized for.

    Parameters
    ----------
    curve : sequence of dict
        A power curve, as returned by `power_curve`.
    target : float, optional
        The power to reach (default 0.8, the convention this package's plots draw a line at).
    **kwargs
        Further arguments passed to `solve_curve`, such as `effect` to pick one focal effect
        out of a curve holding several, or `level` for the interval.

    Returns
    -------
    dict
        The dict `solve_curve` returns, with `n`, `n_lo` and `n_hi` added: `value`, `lo` and
        `hi` rounded up to whole numbers.

    Examples
    --------
    >>> curve = [{"n_subject": 40, "power": 0.38, "n_sims": 100},
    ...          {"n_subject": 100, "power": 0.74, "n_sims": 100},
    ...          {"n_subject": 160, "power": 0.89, "n_sims": 100}]
    >>> target_n(curve)["n"]
    120
    """
    out = solve_curve(curve, target=target, **kwargs)
    out["n"] = math.ceil(out["value"])
    out["n_lo"] = math.ceil(out["lo"])
    out["n_hi"] = math.ceil(out["hi"])
    return out
