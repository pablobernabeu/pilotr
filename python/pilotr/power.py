"""Simulation-based power and design analysis (Type S / Type M).

For each of the `n_sims` Monte Carlo replicates, we simulate a fresh data set from the
known ground-truth spec, fit the analysis model, and record the estimate and p-value. In
addition to classical power (the proportion of significant replicates), we report the
design-analysis quantities of Gelman and Carlin (2014), computed over the significant
replicates.

* Type S (sign) error: P(estimate has the wrong sign | significant)
* Type M (magnitude): E(|estimate| / |true effect| | significant)  (exaggeration ratio)

The two-group Gaussian design with one row per subject uses a two-sample t-test (`power`).
Crossed mixed-effects designs use a statsmodels MixedLM backend (`power_mixed`), which fits
independent variance components by REML and tests the effect with a Wald z. The R package's
lme4 backend fits the model the specification implies and tests with Satterthwaite's
approximation. For accuracy and count responses, it fits a logistic or Poisson model and tests
with a Wald z.

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


# The families whose coefficients are on a link scale, with the scale's name. power_mixed() here
# fits a linear model to every response on its own scale, which for these families is not the
# scale the true effect is written on.
_LINK_SCALE = {"bernoulli": "logit", "poisson": "log", "ordinal": "logit", "beta": "logit"}

# Raised for those families, word for word as the R twin raises it for the ordinal and beta
# families, which it fits in the same way. The R twin fits bernoulli and poisson responses with
# glmer() on their link scale.
_LINEAR_ON_LINK = (
    "power_mixed() fits a linear model to the %s response on its own scale, while the "
    "specification's coefficients are on the %s scale, so the mean estimate and Type M are "
    "withheld (NA) and Type S compares signs only. For a model on the link scale, use "
    "generate_design_analysis() or brms_bridge() in the R package.")


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


# The ConvergenceWarnings that statsmodels raises about a MixedLM fit which may still be at the
# optimum. It flags a fit as possibly on the boundary whenever a variance component is below 0.01
# in absolute terms. On the log scale of a reaction time that is every fit, so singularity is
# judged from the estimates instead. It also reports each optimiser in its chain that gives way to
# the next, though the next may converge. A failure of the last one has a notice of its own.
# Every other ConvergenceWarning makes the replicate a fit with a warning, with the one exception
# below. Counting only a list of known failures would let a reworded failure pass unseen, whereas
# a reworded notice from this list shows up at once, as fits that never converge.
_NOTICES_THAT_PASS = ("The MLE may be on the boundary",
                      "Retrying MixedLM optimization",
                      "Maximum Likelihood optimization failed to converge")

# statsmodels also warns that the Hessian is not positive definite whenever an entry on its
# diagonal is not negative. It takes the Hessian with respect to the variances, and at a component
# held at zero the maximum is set by the boundary, so the curvature there need not be negative. In
# the 164 singular fits that drew this notice over five designs, only components at zero failed
# the check. In a singular fit the notice is therefore left to n_singular, as the R twin leaves
# lme4's notice of a singular fit to isSingular(). Given this model as its formula, the R twin's
# power_mixed() then reports the same five counts as this one on the 12 x 8 design of the Python
# power guide, with no fit warned. Counted there, the notice made 3 of those 12 fits warned.
_HESSIAN_NOTICE = "The Hessian matrix at the estimated parameter values is not positive definite"

# What a replicate fitted by the fallback records against its fit.
_FELL_BACK = "Powell's method raised, so statsmodels' default optimisers fitted the replicate."


def _notices(caught):
    """The messages of the warnings recorded during a MixedLM fit that may count against it."""
    from statsmodels.tools.sm_exceptions import ConvergenceWarning

    return [str(w.message) for w in caught if issubclass(w.category, ConvergenceWarning)
            and not str(w.message).startswith(_NOTICES_THAT_PASS)]


def _fit_mixed(model):
    """Fit a MixedLM by REML and return ``(result, notices)``, or None when no fit returned.

    ``notices`` lists the messages that may count against the fit, which the caller weighs once
    it knows whether the fit is singular. Powell's method runs first, and L-BFGS takes over if it
    does not converge. statsmodels' default chain of BFGS, L-BFGS and CG stopped short of the REML
    optimum in most crossed fits. On the first replicate of the crossed_mixed_rt example, it ended
    at a REML log-likelihood of -406.4, where lme4 reaches -389.8, and gave the effect a standard
    error of 0.0345 against lme4's 0.0186. That chain now serves only when Powell's method raises,
    and the replicate then counts as a fit with a warning. A fit is never repeated because of a
    warning or of its ``converged`` flag, since the boundary notice comes with every fit on the
    log scale, and a default fit flagged as converged could still fall short of the optimum.
    """
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            m = model.fit(reml=True, method=["powell", "lbfgs"])
        except Exception:
            m = None
    if m is not None:
        return m, _notices(caught)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # the fallback already makes this a fit with a warning
        try:
            return model.fit(reml=True), [_FELL_BACK]
        except Exception:
            return None


def _power_mixed_replicate(seed, spec, fname, l2c, yname, fam, shift):
    """One mixed-model replicate: simulate at `seed` and fit MixedLM.

    Returns a dict with the estimate and Wald z p-value of ``cc`` and the fit's flags,
    ``singular``, ``warned`` and ``converged``, or None when no fit returned a finite estimate and
    p-value.
    """
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
    # Building the model is kept out of the fits' error handling, so that a fault there stops the
    # call where it used to leave every replicate without a fit and power undefined.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model = smf.mixedlm("yv ~ cc", df, groups="grp", vc_formula=vcf)
    fit = _fit_mixed(model)
    if fit is None:
        return None
    m, notices = fit
    est, p = float(m.fe_params["cc"]), float(m.pvalues["cc"])
    # A Hessian that is not positive definite can leave the effect without a standard error, and
    # such a replicate has no test to count towards power.
    if not (math.isfinite(est) and math.isfinite(p)):
        return None
    # lme4's isSingular(): a variance component whose standard deviation is below 1e-4 of the
    # residual one. For independent components that ratio is lme4's theta, and on Powell fits the
    # test agreed with lme4 in each of 24 replicates. It is squared here, so that no variance needs
    # a square root.
    singular = any(v < 1e-8 * m.scale for v in m.vcomp)
    # See _HESSIAN_NOTICE: in a singular fit, the component at zero explains that notice.
    warned = any(not (singular and n.startswith(_HESSIAN_NOTICE)) for n in notices)
    return {"est": est, "p": p, "singular": singular, "warned": warned,
            "converged": bool(m.converged) and not singular and not warned}


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

    For ``bernoulli``, ``poisson``, ``ordinal`` and ``beta``, that scale is not the one the
    coefficients are written on, which is the logit scale, or the log scale for ``poisson``. The
    estimate of a linear model is then a difference on the response scale, such as a difference in
    the probability of a correct response, and comparing it with the true coefficient says nothing
    about exaggeration. For these families, `power_mixed` warns and returns `mean_estimate` and
    `type_m` as `nan`. `type_s` compares signs alone, which the link preserves, and power is the
    rate at which the linear model's test rejects. The R package's ``power_mixed()`` fits
    ``bernoulli`` and ``poisson`` responses with ``lme4::glmer()`` on their link scale, and its
    ``generate_design_analysis()`` writes a Bayesian design analysis for any family.

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

    Each replicate is fitted by REML with Powell's method, and L-BFGS takes over if Powell's does
    not converge. statsmodels' default chain of optimisers stopped short of the REML optimum in
    most crossed fits, which earlier versions described as statsmodels overstating random-slope
    variance. Over the first 12 replicates of the crossed_mixed_rt example, the standard error of
    the effect now matches the one from lme4's uncorrelated model, ``(1 + cc || subject) +
    (1 + cc || item)``, to four decimal places in 11. The same 6 replicates are significant in
    both. A replicate whose Powell fit raises is fitted again with the default chain and counts
    as a fit with a warning.

    The p-values are Wald z tests. A Wald z treats the estimate over its standard error as normal
    and so ignores the uncertainty in the variance components. With few subjects or items, it
    rejects more readily than a test with Satterthwaite's degrees of freedom (Luke, 2017,
    doi:10.3758/s13428-016-0809-y), which is the test of the R package's ``power_mixed()``. In
    small designs, the power here therefore tends to exceed R's on the same data. On the 12 x 8
    design of the documentation, 6 of 12 replicates are significant here and 4 of 12 in R. With
    8 subjects, 6 items, random intercepts alone and no true effect, the test here rejected the
    null hypothesis in 0.043 of 300 replicates (Monte Carlo standard error 0.012). On the same
    data, lme4's Satterthwaite test of the same model rejected in 0.030.

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
        Keys: `backend` (the estimator and its test), `n_sims`, `alpha`, the fit counts
        `n_attempted`, `n_returned`, `n_converged`, `n_singular` and `n_warning`, then `power`,
        `n_significant`, `true_effect`, `mean_estimate`, `type_s` and `type_m`.

        The fit counts mean what they mean in the R twin. `n_attempted` is `n_sims`, and
        `n_returned` counts the replicates that returned an estimate and a p-value. `power` is
        the proportion of significant results among those, not among `n_sims`. `n_singular`
        counts fits with a variance component whose standard deviation is below 1e-4 of the
        residual one, the test of lme4's ``isSingular()``. `n_warning` counts fits whose
        optimisers did not converge, whose Powell fit raised or whose Hessian was not positive
        definite. A component at zero can fail statsmodels' check of the Hessian by itself, so
        that notice counts only against a fit that is not singular. `n_converged` counts the fits
        that are neither singular nor warned. statsmodels' notice that a fit may be on the
        boundary counts towards none of them, since it comes with every fit that has a variance
        component below 0.01. Singular fits and fits with a warning stay in `power`, as in the R
        twin. `type_s` and `type_m` are `nan` when no replicate reached significance and when the
        true effect is zero. `mean_estimate` and `type_m` are also `nan` for a ``bernoulli``,
        ``poisson``, ``ordinal`` or ``beta`` response, whose true effect is on a link scale.

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
    UserWarning
        Once per call, for a ``bernoulli``, ``poisson``, ``ordinal`` or ``beta`` response,
        saying that the linear model is fitted on the response's own scale and that the mean
        estimate and Type M are withheld. The R twin raises the same words for the families it
        fits in that way.

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
    # The warnings are raised here, once per call, since the replicates may run in other
    # processes.
    left_out = _left_out_of_model(spec, col)
    if left_out:
        warnings.warn(_LEFT_OUT % _name_list(left_out), UserWarning, stacklevel=2)
    # A linear model of a 0/1, count, rating or proportion response estimates a difference on the
    # response scale, which used to be set against the true value on the link scale. On the R
    # twin's crossed accuracy design, such a model gave Type M 0.21 where a logistic one gives
    # 1.12. Type S survives, since a monotone link keeps the sign of an effect but not its size.
    withheld = fam in _LINK_SCALE
    if withheld:
        warnings.warn(_LINEAR_ON_LINK % (fam, _LINK_SCALE[fam]), UserWarning, stacklevel=2)

    rep = functools.partial(_power_mixed_replicate, spec=spec, fname=fname, l2c=l2c,
                            yname=yname, fam=fam, shift=shift)
    seeds = replicate_seeds(base, n_sims)
    if workers == 1:
        results = _map_replicates(rep, seeds, None)
    else:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(max_workers=workers) as executor:
            results = _map_replicates(rep, seeds, executor)
    fits = [r for r in results if r is not None]
    est = [r["est"] for r in fits]
    pv = [r["p"] for r in fits]

    sig = [i for i, p in enumerate(pv) if p < alpha]
    # See _power_impl: a zero true effect leaves both design-analysis quantities undefined.
    usable = bool(sig) and not math.isnan(beta) and beta != 0
    return {
        "backend": "statsmodels MixedLM (crossed variance components, REML; Wald z tests)",
        "n_sims": n_sims, "alpha": alpha,
        # The R twin's fit counts, with its meanings. n_converged used to count every fit that
        # returned, the count now called n_returned.
        "n_attempted": n_sims, "n_returned": len(fits),
        "n_converged": sum(r["converged"] for r in fits),
        "n_singular": sum(r["singular"] for r in fits),
        "n_warning": sum(r["warned"] for r in fits),
        "power": len(sig) / len(pv) if pv else float("nan"),
        "n_significant": len(sig), "true_effect": beta,
        "mean_estimate": statistics.mean(est) if est and not withheld else float("nan"),
        "type_s": (sum(1 for i in sig if (est[i] > 0) != (beta > 0)) / len(sig))
                  if usable else float("nan"),
        "type_m": statistics.mean(abs(est[i]) / abs(beta) for i in sig)
                  if usable and not withheld else float("nan"),
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
