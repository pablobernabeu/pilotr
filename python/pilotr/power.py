"""Simulation-based power and design analysis (Type S / Type M).

For each of the `n_sims` Monte Carlo replicates, we simulate a fresh data set from the
known ground-truth spec, fit the analysis model, and record the estimate and p-value. In
addition to classical power (the proportion of significant replicates), we report the
design-analysis quantities of Gelman and Carlin (2014), computed over the significant
replicates.

* Type S (sign) error: P(estimate has the wrong sign | significant)
* Type M (magnitude): E(|estimate| / |true effect| | significant)  (exaggeration ratio)

The two-group Gaussian design with one row per subject uses a two-sample t-test (`power`).
Crossed mixed-effects designs use a statsmodels MixedLM backend (`power_mixed`), which is
conservative for random-slope designs; the R package's lme4 backend is the reference there.

Every analysis takes a `workers` argument that spreads the replicates over local
processes. The replicate seeds are derived once from the specification's seed, so the
results are identical to a serial run whatever the worker count. The per-replicate
functions are module-level so that they pickle under the Windows spawn start method.
"""

from __future__ import annotations
import copy, functools, math, statistics, warnings
from .simulate import _simulate, _as_spec, load_spec
from .core import replicate_seeds
from .validate import _name_list


def _check_workers(workers):
    """Validate `workers` as a single positive whole number and return it as an int."""
    if isinstance(workers, bool) or not isinstance(workers, int) or workers < 1:
        raise ValueError("workers must be a positive whole number")
    return workers


def _map_replicates(rep, seeds, executor):
    """Run `rep` over the per-replicate seeds, serially or on the executor.

    `executor.map` returns results in input order, so the downstream reductions match the
    serial code exactly.
    """
    if executor is None:
        return [rep(seed) for seed in seeds]
    n_workers = getattr(executor, "_max_workers", 1) or 1
    return list(executor.map(rep, seeds, chunksize=max(1, len(seeds) // (4 * n_workers))))


# The refusal for a design whose rows are correlated, byte-identical in the R twin. It sends the
# user of either twin to R's power_mixed(), since power_mixed() here fits a single within factor
# crossed with items and refuses every between-subjects design.
_CORRELATED_ROWS = (
    "power_design() in R and power() in Python t-test every row as an independent observation, "
    "which is valid only when each subject contributes one row and no rows share a cluster. "
    "This design has %s, so its rows are correlated and the t-test would overstate power and "
    "understate Type M. Use power_mixed() in the R package, which fits the model the "
    "specification implies.")


def _two_group_refusal(spec):
    """Why the two-group backend cannot analyse a validated specification, or None when it can.

    A t-test of every row needs independent rows. It used to be applied to any Gaussian design
    with a two-level between factor, and under the null the test of 30 subjects crossed with 20
    items was then significant in 105 of 200 replicates. An item unit, a within factor or an
    extra grouping factor is therefore refused, the first found being named. Predictors and a
    ``subject`` entry stay allowed, since with one row per subject they vary independently from
    row to row. A between factor on items needs an item unit, so the item test covers it.
    Mirrors ``.two_group_refusal()`` in power.R, message for message.
    """
    if spec["response"]["family"] != "gaussian":
        return "The power backend currently handles only the gaussian two-group design."
    factors = spec.get("factors") or []
    between = [f for f in factors if f.get("between")]
    if len(between) != 1 or len(between[0]["levels"]) != 2:
        return "The power backend expects exactly one 2-level between factor."
    # A factor can set `vary_within` alongside `between`, and it then varies within subjects.
    within = [f for f in factors if f.get("vary_within")]
    # The grouping factors _simulate() draws: an `item` entry without an item unit groups
    # nothing, and with one the item unit is named first.
    extra = [g for g in (spec.get("random") or {}) if g not in ("subject", "item")]
    if "item" in spec["units"]:
        what = "an item unit"
    elif within:
        what = "the within factor '%s'" % within[0]["name"]
    elif extra:
        what = "the grouping factor '%s'" % extra[0]
    else:
        return None
    return _CORRELATED_ROWS % what


def _power_replicate(seed, spec, fname, lev0, lev1, yname):
    """One two-group replicate: simulate at `seed`, t-test, return (estimate, p-value)."""
    from scipy import stats  # lazy: imported in each worker process on first use

    # The specification was validated and normalised once, before the loop.
    d = _simulate(dict(spec, seed=seed))
    g0 = [r[yname] for r in d.rows if r[fname] == lev0]
    g1 = [r[yname] for r in d.rows if r[fname] == lev1]
    t = stats.ttest_ind(g1, g0, equal_var=True)
    return statistics.mean(g1) - statistics.mean(g0), t.pvalue


def power(spec, n_sims=1000, alpha=0.05, workers=1):
    """Simulation-based power and design analysis for a two-group Gaussian design.

    Repeatedly simulate from the specification, apply a two-sample t-test, and report power
    together with the Type S (sign) and Type M (magnitude) errors of Gelman and Carlin
    (2014), computed over the significant replicates.

    The t-test treats every row as an independent observation, which is valid only when each
    subject contributes one row and no rows share a cluster. `power` therefore takes a Gaussian
    design with exactly one two-level factor between subjects and no item unit, within factor or
    grouping factor besides ``subject``. Continuous predictors and by-subject random effects are
    allowed, since with one row per subject they vary independently from row to row. Any other
    design is refused with a message naming what makes its rows correlated. The R package's
    ``power_mixed()`` fits the model that such a design implies. Take 30 subjects crossed with
    20 items, with by-subject and residual standard deviations of 1 and a by-item one of 0.3.
    With no true effect, a t-test of every row is significant in about half of all replicates.

    Parameters
    ----------
    spec : dict or str
        A two-group Gaussian design specification with one row per subject (dict or path to a
        JSON file).
    n_sims : int, optional
        Number of Monte Carlo replicates (default 1000).
    alpha : float, optional
        Two-sided significance level (default 0.05).
    workers : int, optional
        Number of local worker processes over which to spread the replicates (default 1,
        serial). The replicate seeds are derived once from the specification's seed, so any worker
        count returns results identical to a serial run.

    Returns
    -------
    dict
        Keys: `n_sims`, `alpha`, `power`, `n_significant`, `true_effect`, `mean_estimate`,
        `type_s`, `type_m`. Both design-analysis quantities are `nan` when no replicate
        reached significance and when the true effect is zero, as in a null condition:
        neither is defined without a true value to compare against, and Type M divides by it.
        A specification with no coefficient for the factor's contrast has a true effect of 0.

    Raises
    ------
    NotImplementedError
        If the design is not Gaussian, does not have exactly one two-level between factor, or
        has rows that are correlated because of an item unit, a within factor or a grouping
        factor besides ``subject``.

    Notes
    -----
    Requires `scipy` (imported lazily, in each worker process when parallel); install the
    `power` or `dev` extra.
    """
    workers = _check_workers(workers)
    from scipy import stats as _stats  # noqa: F401  fail fast before simulating

    if workers == 1:
        return _power_impl(spec, n_sims, alpha, None)
    from concurrent.futures import ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=workers) as executor:
        return _power_impl(spec, n_sims, alpha, executor)


def _power_impl(spec, n_sims, alpha, executor):
    """The replicate loop behind `power`, taking an optional executor so that sweep
    functions can start one process pool and reuse it across grid points."""
    spec = _as_spec(spec)
    refusal = _two_group_refusal(spec)
    if refusal is not None:
        raise NotImplementedError(refusal)

    factor = next(f for f in spec["factors"] if f.get("between"))
    fname = factor["name"]
    lev0, lev1 = factor["levels"]
    col, vals = next(iter(factor["contrasts"].items()))
    # A specification may leave the contrast out of `coefficients`, `{}` included, and the
    # simulator then generates no effect, so the true effect is 0. The lookup used to raise
    # KeyError.
    true_effect = spec["fixed"]["coefficients"].get(col, 0.0) * (vals[1] - vals[0])
    yname = spec["response"]["name"]

    base_seed = spec["seed"]
    rep = functools.partial(_power_replicate, spec=spec, fname=fname,
                            lev0=lev0, lev1=lev1, yname=yname)
    results = _map_replicates(rep, replicate_seeds(base_seed, n_sims), executor)
    estimates = [r[0] for r in results]
    pvals = [r[1] for r in results]

    sig = [i for i, p in enumerate(pvals) if p < alpha]
    power_val = len(sig) / n_sims
    # Type S and Type M are defined relative to a true value, and Type M divides by it, so both
    # stay nan when the true effect is zero. The alternative was a division by zero, alongside a
    # sign-error rate that had quietly become "the estimate is positive". Same rule as
    # power_mixed(), and the null condition the sweep documentation recommends is exactly this
    # case.
    if sig and not math.isnan(true_effect) and true_effect != 0:
        type_s = sum(1 for i in sig if (estimates[i] > 0) != (true_effect > 0)) / len(sig)
        type_m = statistics.mean(abs(estimates[i]) / abs(true_effect) for i in sig)
    else:
        type_s = type_m = float("nan")

    return {
        "n_sims": n_sims, "alpha": alpha, "power": power_val,
        "n_significant": len(sig), "true_effect": true_effect,
        "mean_estimate": statistics.mean(estimates),
        "type_s": type_s, "type_m": type_m,
    }


# What power_mixed() fits, stated whenever a specification declares more. The model is written
# into _power_mixed_replicate() and ignores the rest of the specification, so whatever else the
# specification declared used to drop out of the analysis without a word. The R twin builds its
# model from the specification, through model_formula().
_LEFT_OUT = (
    "power_mixed() in Python fits 'yv ~ cc' with independent by-subject and by-item intercept "
    "and slope components, testing the first contrast of the within factor only. This "
    "specification also declares %s, which the fitted model leaves out. The R package's "
    "power_mixed() fits the model the specification implies.")


def _left_out_of_model(spec, col):
    """What a validated specification declares beyond the model ``yv ~ cc``, as phrases.

    `col` is the tested contrast. Each phrase names one kind of term and lists its members in
    the specification's order, and an empty list means that the model leaves nothing out. A
    coefficient counts only when it is not zero, since a zero effect left out of the fixed part
    is still the effect simulated. Correlations count as the R twin's model_formula() reads
    them, when ``correlated`` is true or is absent alongside ``correlations``, and only for a
    group with a slope to correlate with its intercept.
    """
    rs = spec["random"]

    def correlated(u):
        return bool(u.get("slopes")) and u.get("correlated", bool(u.get("correlations")))

    found = [
        ("another coefficient", "other coefficients",
         [k for k, v in spec["fixed"]["coefficients"].items() if k != col and v != 0]),
        ("a between factor", "between factors",
         [f["name"] for f in spec["factors"] if f.get("between")]),
        ("a predictor", "predictors", [p["name"] for p in spec["predictors"]]),
        ("a by-subject slope", "by-subject slopes",
         [k for k in rs.get("subject", {}).get("slopes", {}) if k != col]),
        ("a by-item slope", "by-item slopes",
         [k for k in rs.get("item", {}).get("slopes", {}) if k != col]),
        ("random-effect correlations", "random-effect correlations",
         [g for g in ("subject", "item") if g in rs and correlated(rs[g])]),
        ("an extra grouping factor", "extra grouping factors",
         [g for g in rs if g not in ("subject", "item")]),
    ]
    return ["%s (%s)" % (one if len(names) == 1 else many, ", ".join("'%s'" % n for n in names))
            for one, many, names in found if names]


def _power_mixed_replicate(seed, spec, fname, l2c, yname, fam, shift):
    """One mixed-model replicate: simulate at `seed`, fit MixedLM, return (estimate,
    p-value), or None when the fit fails."""
    import math, warnings
    import pandas as pd
    import statsmodels.formula.api as smf

    vcf = {"subj_i": "0 + C(subject)", "subj_s": "0 + C(subject):cc",
           "item_i": "0 + C(item)", "item_s": "0 + C(item):cc"}
    df = pd.DataFrame(_simulate(dict(spec, seed=seed)).rows)
    df["cc"] = df[fname].map(l2c)
    # Both lognormal families are analysed on the log scale, where their coefficients are, as
    # model_data() does in the R twin. The plain lognormal has no shift, so `shift` is 0 for it.
    # Only the shifted family used to be logged, so the estimate for a lognormal response came
    # out in the response's own units.
    if fam in ("lognormal", "shifted_lognormal"):
        # A response rounded down to the shift has no logarithm. R's log() gives -Inf for it and
        # lmer() refuses the fit, so the replicate returns no estimate here too. math.log()
        # raised ValueError on it, which ended the whole call.
        if any(v - shift <= 0 for v in df[yname]):
            return None
        df["yv"] = [math.log(v - shift) for v in df[yname]]
    else:
        df["yv"] = list(df[yname])
    df["grp"] = 1
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            m = smf.mixedlm("yv ~ cc", df, groups="grp", vc_formula=vcf).fit()  # REML (default)
        return float(m.fe_params["cc"]), float(m.pvalues["cc"])
    except Exception:
        return None


def power_mixed(spec, n_sims=50, alpha=0.05, workers=1):
    """Crossed mixed-effects simulation-based power in Python, via statsmodels MixedLM with
    by-subject and by-item random intercepts and slopes as (independent) variance components.

    The model is the same whatever the specification declares, and it tests one effect. It is
    ``yv ~ cc``, which regresses the response on ``cc``, the first contrast of the design's
    single within factor, with four independent variance components: by-subject and by-item
    intercepts and slopes on ``cc``. The response is analysed as ``log(y)`` for the
    ``lognormal`` family and as ``log(y - shift)`` for ``shifted_lognormal``, whose
    coefficients are on the log scale, as in the R twin's ``model_data()``. Every other family
    is analysed on its own scale. A replicate whose response rounds to 0, or to the shift, has
    no logarithm and returns no estimate, as in the R twin.

    A specification that declares anything this model leaves out draws one warning per call
    naming each such term. In the fixed part, that is a non-zero coefficient other than
    the tested contrast, a between factor or a predictor. In the random part, it is a
    by-subject or by-item slope on another term, a correlation between random effects or a
    grouping factor besides ``subject`` and ``item``. The R package's ``power_mixed()`` fits
    the model the specification implies and tests every coefficient.

    This is pilotr's own simulation loop over the portable design specification, not a
    wrapper around an existing power package. It covers territory pioneered by simr (Green
    and MacLeod, 2016, doi:10.1111/2041-210x.12504) and mixedpower (Kumle, Vo and Draschkow,
    2021, doi:10.3758/s13428-021-01546-0); pilotr differs in being driven by the portable
    cross-language spec, in reporting Type S and Type M errors and in built-in
    parallelisation via `workers` (default 1, serial), which returns results identical to a
    serial run for any worker count.

    Accuracy caveat (verified behaviour, not a bug). statsmodels fits crossed random effects
    as independent variance components and, in our tests, substantially overstates random-slope
    variance (for example, a by-subject slope SD of about 0.12 estimated against 0.04 true).
    This inflates the fixed-effect standard error. For designs with by-subject or by-item random
    slopes, the backend is therefore markedly conservative. On the crossed RT design it
    reports power around 0.48, against the R/lme4 reference of about 0.73. It still recovers the
    fixed effect (mean estimate about 0.048 against 0.05 true) and the Type S and Type M
    quantities correctly, and it is reliable for random-intercept designs. We recommend treating
    its output as a conservative lower bound and using the R/lme4 `power_mixed` as the reference
    whenever random slopes or random-effect correlations matter. Data generation is identical
    across R and Python. This discrepancy arises solely in the Python LMM estimator.

    Parameters
    ----------
    spec : dict or str
        A design specification (dict or path to a JSON file) with exactly one within-unit
        factor and a crossed design with an item unit.
    n_sims : int, optional
        Number of Monte Carlo replicates (default 50, smaller than `power`'s 1000 because
        each replicate fits a mixed model).
    alpha : float, optional
        Two-sided significance level (default 0.05).
    workers : int, optional
        Number of local worker processes over which to spread the replicates (default 1,
        serial). The replicate seeds are derived once from the specification's seed, so any worker
        count returns results identical to a serial run.

    Returns
    -------
    dict
        Keys: `backend` (the estimator used), `n_sims`, `n_converged` (how many replicates
        the model fit), `alpha`, `power`, `n_significant`, `true_effect`, `mean_estimate`,
        `type_s`, `type_m`. `power` is the proportion of significant results among the
        `n_converged` converged replicates, not among `n_sims`. `type_s` and `type_m` are
        `nan` when no replicate reached significance and when the true effect is zero.

    Raises
    ------
    ValueError
        If the spec has no item unit (`power_mixed` requires a crossed design).
    NotImplementedError
        If the design does not have exactly one within-unit factor.

    Warns
    -----
    UserWarning
        Once per call, when the specification declares a term that the fitted model leaves
        out. The message names every such term.

    Notes
    -----
    Requires the `mixed` extra (`statsmodels` and `pandas`, imported lazily in each worker
    process when parallel).
    """
    workers = _check_workers(workers)
    import pandas as _pd  # noqa: F401  fail fast before simulating
    import statsmodels.formula.api as _smf  # noqa: F401

    spec = _as_spec(spec)
    if "item" not in spec["units"]:
        raise ValueError("power_mixed() requires a crossed design with an item unit.")
    within = [f for f in spec["factors"] if f.get("vary_within")]
    if len(within) != 1:
        raise NotImplementedError(
            "The Python power_mixed backend expects exactly one within factor.")
    f = within[0]
    fname = f["name"]
    col, vals = next(iter(f["contrasts"].items()))
    l2c = {f["levels"][i]: vals[i] for i in range(len(f["levels"]))}
    # No coefficient for the contrast means no effect, as in _power_impl.
    beta = spec["fixed"]["coefficients"].get(col, 0.0)
    yname, fam = spec["response"]["name"], spec["response"]["family"]
    shift = spec["response"].get("shift", 0.0)
    base = spec["seed"]
    # The warning is raised here, once per call, since the replicates may run in other processes.
    left_out = _left_out_of_model(spec, col)
    if left_out:
        warnings.warn(_LEFT_OUT % _name_list(left_out), UserWarning, stacklevel=2)

    rep = functools.partial(_power_mixed_replicate, spec=spec, fname=fname, l2c=l2c,
                            yname=yname, fam=fam, shift=shift)
    seeds = replicate_seeds(base, n_sims)
    if workers == 1:
        results = _map_replicates(rep, seeds, None)
    else:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(max_workers=workers) as executor:
            results = _map_replicates(rep, seeds, executor)
    est = [r[0] for r in results if r is not None]
    pv = [r[1] for r in results if r is not None]

    sig = [i for i, p in enumerate(pv) if p < alpha]
    # See _power_impl: a zero true effect leaves both design-analysis quantities undefined.
    usable = bool(sig) and not math.isnan(beta) and beta != 0
    return {
        "backend": "statsmodels MixedLM (crossed variance components, REML)",
        "n_sims": n_sims, "n_converged": len(pv), "alpha": alpha,
        "power": len(sig) / len(pv) if pv else float("nan"),
        "n_significant": len(sig), "true_effect": beta,
        "mean_estimate": statistics.mean(est) if est else float("nan"),
        "type_s": (sum(1 for i in sig if (est[i] > 0) != (beta > 0)) / len(sig))
                  if usable else float("nan"),
        "type_m": statistics.mean(abs(est[i]) / abs(beta) for i in sig)
                  if usable else float("nan"),
    }


def power_curve(spec, subject_ns, n_sims=1000, alpha=0.05, workers=1):
    """Power curve over sample size for a two-group Gaussian design.

    Sweep the number of subjects and compute `power` at each grid point.

    Parameters
    ----------
    spec : dict or str
        A two-group Gaussian design specification with one row per subject, as `power`
        requires.
    subject_ns : iterable of int
        Subject counts to evaluate.
    n_sims : int, optional
        Replicates per grid point (default 1000).
    alpha : float, optional
        Significance level (default 0.05).
    workers : int, optional
        Number of local worker processes over which to spread the replicates at each grid
        point (default 1, serial). One process pool is started and reused across the whole
        sweep, and any worker count returns results identical to a serial run.

    Returns
    -------
    list of dict
        One dict per grid point, with keys `n_subject`, `power`, `type_m`, `n_sims` and
        `n_significant`. The two counts are what `solve_curve` weights the curve by, since a
        power estimate over 1000 replicates should count for more than one over 50.

    See Also
    --------
    solve_curve, target_n : solve this curve for the sample size that meets a target power.
    """
    workers = _check_workers(workers)
    if isinstance(spec, str):
        spec = load_spec(spec)

    def sweep(executor):
        out = []
        for n in subject_ns:
            s = copy.deepcopy(spec)
            s["units"]["subject"]["n"] = n
            r = _power_impl(s, n_sims, alpha, executor)
            out.append({"n_subject": n, "power": r["power"], "type_m": r["type_m"],
                        "n_sims": r["n_sims"], "n_significant": r["n_significant"]})
        return out

    if workers == 1:
        return sweep(None)
    from concurrent.futures import ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=workers) as executor:  # one pool for the sweep
        return sweep(executor)
