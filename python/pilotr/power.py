"""Simulation-based power and design analysis (Type S / Type M).

For each of the `n_sims` Monte Carlo replicates, we simulate a fresh data set from the
known ground-truth spec, fit the analysis model, and record the estimate and p-value. In
addition to classical power (the proportion of significant replicates), we report the
design-analysis quantities of Gelman and Carlin (2014), computed over the significant
replicates.

* Type S (sign) error: P(estimate has the wrong sign | significant)
* Type M (magnitude): E(|estimate| / |true effect| | significant)  (exaggeration ratio)

Each of the three comes with its Monte Carlo error, from `pilotr.mcse`: power with its standard
error and Wilson interval, Type S with a Wilson interval and Type M with its standard error.

The two-group Gaussian design with one row per subject uses a two-sample t-test (`power`).
Mixed-effects designs use a statsmodels MixedLM backend (`power_mixed`), which fits the model the
specification implies, with each random effect as an independent variance component, by REML and
tests each effect with a Wald z. The R package's lme4 backend fits the same model with its
correlations and tests with Satterthwaite's approximation. For accuracy and count responses, it
fits a logistic or Poisson model and tests with a Wald z.

Every analysis takes a `workers` argument that spreads the replicates over local
processes. The replicate seeds are derived once from the specification's seed, so the
results are identical to a serial run whatever the worker count. The per-replicate
functions are module-level so that they pickle under the Windows spawn start method.
"""

from __future__ import annotations
import copy, functools, math, statistics, warnings
from collections.abc import Mapping
from .simulate import _simulate, _as_spec
from .core import replicate_seeds
from .autoformula import _interaction_keys, _model_data, _model_terms, _us
from .mcse import _rate_with_error, _type_m_mcse, _type_s_interval
from .validate import _is_num, _is_whole, _name_list


# The alternatives a test may take, spelt as in the R twin, whose default is the first.
_ALTERNATIVES = ("two.sided", "greater", "less")


def _check_alternative(alternative):
    """Validate `alternative` as one of the R twin's three values and return it.

    scipy spells the two-sided test "two-sided", and that spelling is refused here like any other,
    so that a call means the same in both twins and a result records one spelling.
    """
    if not isinstance(alternative, str) or alternative not in _ALTERNATIVES:
        raise ValueError("alternative must be 'two.sided', 'greater' or 'less'")
    return alternative


def _check_workers(workers):
    """Validate `workers` as a single positive whole number and return it as an int.

    Any whole number is accepted, a NumPy integer or a float such as 2.0 included, as the R twin's
    ``.check_workers()`` accepts one. Only a Python int used to pass.
    """
    if not _is_whole(workers) or workers < 1:
        raise ValueError("workers must be a positive whole number")
    return int(workers)


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
# user of either twin to R's power_mixed(), which fits every such design. power_mixed() here fits
# one only when it declares random effects.
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


def _power_replicate(seed, spec, fname, lev0, lev1, yname, alternative="two.sided"):
    """One two-group replicate: simulate at `seed`, t-test, return (estimate, p-value).

    scipy takes a one-sided p-value from the same t as the two-sided one, the upper tail for
    "greater" and the lower for "less", so it is half the two-sided p-value when t lies in the
    tested direction, as R's ``t.test()`` gives it.
    """
    from scipy import stats  # lazy: imported in each worker process on first use

    # The specification was validated and normalised once, before the loop.
    d = _simulate(dict(spec, seed=seed))
    g0 = [r[yname] for r in d.rows if r[fname] == lev0]
    g1 = [r[yname] for r in d.rows if r[fname] == lev1]
    t = stats.ttest_ind(g1, g0, equal_var=True,
                        alternative="two-sided" if alternative == "two.sided" else alternative)
    return statistics.mean(g1) - statistics.mean(g0), t.pvalue


def power(spec, n_sims=1000, alpha=0.05, workers=1, alternative="two.sided"):
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

    A directional prediction calls for a one-sided test, which `alternative` selects. The
    estimate is the mean of the second level minus that of the first, whose true value is
    `true_effect`. ``"greater"`` tests whether it is positive and ``"less"`` whether it is
    negative, by ``scipy.stats.ttest_ind`` with the same alternative. The one-sided p-value is
    half the two-sided one when the estimate lies in the tested direction. A one-sided test at
    `alpha` therefore rejects in exactly the replicates in which the two-sided test at
    ``2 * alpha`` rejects with the estimate in that direction. Type S and Type M are taken over
    the replicates significant in the tested direction. When that is the direction of the true
    effect, Type S is 0 by construction, and Type M is the exaggeration among the replicates that
    would support the prediction.

    Results from earlier versions were two-sided. A one-sided power at `alpha`, in the direction
    of the true effect, can be read off a two-sided result at ``2 * alpha``: it is
    ``power * (1 - type_s)``. Type M is recovered only when that result's Type S is 0, since the
    wrong-signed replicates cannot otherwise be told apart from the others.

    Parameters
    ----------
    spec : dict, str or os.PathLike
        A two-group Gaussian design specification with one row per subject (dict or path to a
        JSON file).
    n_sims : int, optional
        Number of Monte Carlo replicates (default 1000).
    alpha : float, optional
        Significance level of the test (default 0.05), two-sided unless `alternative` says
        otherwise.
    workers : int, optional
        Number of local worker processes over which to spread the replicates (default 1,
        serial). The replicate seeds are derived once from the specification's seed, so any worker
        count returns results identical to a serial run.
    alternative : {"two.sided", "greater", "less"}, optional
        The test's alternative hypothesis (default ``"two.sided"``), spelt as in the R twin.
        ``"greater"`` and ``"less"`` give a one-sided test of a directional prediction.

    Returns
    -------
    dict
        Keys: `n_sims`, `alpha`, `alternative`, `power`, `power_mcse` (its Monte Carlo
        standard error), `power_lo` and `power_hi` (its 95% Wilson interval), `n_significant`,
        `true_effect`, `mean_estimate`, `type_s`, `type_s_lo` and `type_s_hi` (its 95% Wilson
        interval), `type_m` and `type_m_mcse` (its Monte Carlo standard error). Both design-analysis
        quantities are `nan` when no replicate reached significance and when the true effect is
        zero, as in a null condition: neither is defined without a true value to compare
        against, and Type M divides by it. Their interval and standard error are `nan` then too,
        and `type_m_mcse` is also `nan` when a single replicate reached significance, since a
        standard deviation needs two. A specification with no coefficient for the factor's
        contrast has a true effect of 0.

        Only the significant replicates inform Type S and Type M, and they are often few. Type S
        comes with a Wilson interval, since its standard error is zero whenever no sign is wrong.
        Type M comes with its standard error, the standard deviation of the exaggeration ratios
        over the square root of their number. Each is computed as in the R twin's
        ``power_design()``, operation for operation, so the same counts give the same doubles.

    Raises
    ------
    NotImplementedError
        If the design is not Gaussian, does not have exactly one two-level between factor, or
        has rows that are correlated because of an item unit, a within factor or a grouping
        factor besides ``subject``.
    ValueError
        If `alternative` is not ``"two.sided"``, ``"greater"`` or ``"less"``.

    Notes
    -----
    Requires `scipy` (imported lazily, in each worker process when parallel); install the
    `power` or `dev` extra.
    """
    workers = _check_workers(workers)
    alternative = _check_alternative(alternative)
    from scipy import stats as _stats  # noqa: F401  fail fast before simulating

    if workers == 1:
        return _power_impl(spec, n_sims, alpha, None, alternative)
    from concurrent.futures import ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=workers) as executor:
        return _power_impl(spec, n_sims, alpha, executor, alternative)


def _power_impl(spec, n_sims, alpha, executor, alternative="two.sided"):
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
                            lev0=lev0, lev1=lev1, yname=yname, alternative=alternative)
    results = _map_replicates(rep, replicate_seeds(base_seed, n_sims), executor)
    estimates = [r[0] for r in results]
    pvals = [r[1] for r in results]

    # Under a one-sided test, the replicates significant in the tested direction, so Type S and
    # Type M below are taken over them alone.
    sig = [i for i, p in enumerate(pvals) if p < alpha]
    # Type S and Type M are defined relative to a true value, and Type M divides by it, so both
    # stay nan when the true effect is zero. The alternative was a division by zero, alongside a
    # sign-error rate that had quietly become "the estimate is positive". Same rule as
    # power_mixed(), and the null condition the sweep documentation recommends is exactly this
    # case. Each carries its Monte Carlo error, nan wherever the rate is.
    if sig and not math.isnan(true_effect) and true_effect != 0:
        n_wrong = sum(1 for i in sig if (estimates[i] > 0) != (true_effect > 0))
        ratios = [abs(estimates[i]) / abs(true_effect) for i in sig]
        type_s = n_wrong / len(sig)
        type_s_lo, type_s_hi = _type_s_interval(n_wrong, len(sig), true_effect)
        type_m = statistics.mean(ratios)
        type_m_mcse = _type_m_mcse(ratios)  # nan with a single significant replicate
    else:
        type_s = type_s_lo = type_s_hi = type_m = type_m_mcse = float("nan")

    return {
        "n_sims": n_sims, "alpha": alpha, "alternative": alternative,
        **_rate_with_error(len(sig), n_sims, "power"),
        "n_significant": len(sig), "true_effect": true_effect,
        "mean_estimate": statistics.mean(estimates),
        "type_s": type_s, "type_s_lo": type_s_lo, "type_s_hi": type_s_hi,
        "type_m": type_m, "type_m_mcse": type_m_mcse,
    }


# What power_mixed() leaves out of the model the specification implies. It fits each random effect
# as an independent variance component, so it estimates no correlation between a grouping factor's
# intercept and slopes, where the R twin's model_formula() writes `|` and lme4 estimates them. Up
# to 0.3, the model was written into the source as `yv ~ cc`, and every other term the
# specification declared was left out with it.
_UNCORRELATED = (
    "power_mixed() in Python fits each random effect as an independent variance component, so it "
    "leaves out the random-effect correlations this specification declares for %s. The R "
    "package's power_mixed() estimates them.")

# Raised for a design with no random effect, which leaves MixedLM no variance component to fit.
# The R twin fits such a design with lm() or glm().
_NO_RANDOM = (
    "power_mixed() in Python fits a mixed model, and this specification declares no random "
    "effects. For a two-group Gaussian design with one row per subject, use power(). The R "
    "package's power_mixed() fits a design without random effects with lm() or glm().")

# Word for word as the R twin's power_mixed() stops when there is no effect to test.
_NOTHING_TO_TEST = ("this specification has no fixed coefficients to test; name the focal effects "
                    "explicitly")

# Raised for a `focal` that is none of the three forms, such as the named sets of coefficients
# that the R twin tests jointly.
_FOCAL_FORMS = (
    "focal must be None, a list of effect names or a dict from effect names to true values. A "
    "joint test of a set of coefficients is planned only in the R package's power_mixed().")

# Word for word as the R twin's .warn_absent_focal() warns.
_ABSENT_FOCAL = ("focal effect%s %s never appeared in the fitted model, so %s results are empty by "
                 "construction; the fitted coefficients are %s")

# What each entry of `effects` holds, in this order, as the R twin's per-effect elements do.
_EFFECT_KEYS = ("power", "power_mcse", "power_lo", "power_hi", "n_significant", "true_effect",
                "mean_estimate", "type_s", "type_s_lo", "type_s_hi", "type_m", "type_m_mcse")

# Up to 0.3, the result held one effect's quantities under these keys, those of the first contrast
# of the single within factor. They are kept for one release beside `effects`, and reading one
# warns.
_SCALAR_DEPRECATED = (
    "result['%s'] from power_mixed() is deprecated and will be removed in the next release. Read "
    "result['effects']['%s']['%s'], which holds the same value for '%s', the first contrast of "
    "the within factor.")

# The size of the random-effects design, rows times random coefficients, from which its matrices
# are held sparse. Dense matrices suit a small design, and sparse ones a large one. With one BLAS
# thread on a Windows laptop, a replicate of the 12 x 8 design of the Python guide (7,680 entries)
# took 0.2 s dense and 1.0 s sparse, and one of the 24 x 20 design (84,480) 1.1 s and 1.4 s. One of
# crossed_mixed_rt (155,520) took 1.8 s and 1.3 s, of nested_clusters (364,800) 7 to 12 s and 2 to
# 3 s, and of reading_time_continuous (760,000) 18 to 20 s and 7 s. The two reached the same
# optimum, with standard errors that agreed to five significant digits or more.
_SPARSE_FROM = 100_000


class _MixedResult(dict):
    """The dict `power_mixed` returns, which warns when a deprecated scalar key is read.

    `within` names the contrast whose quantities the scalar keys repeat, or is None when the
    result has none. Iterating, comparing and serialising the dict read no key, so they draw no
    warning.
    """

    def __init__(self, *args, within=None, **kwargs):
        super().__init__(*args, **kwargs)
        self._within = within

    def _warn(self, key):
        within = getattr(self, "_within", None)
        if within is not None and key in _EFFECT_KEYS and dict.__contains__(self, key):
            warnings.warn(_SCALAR_DEPRECATED % (key, within, key, within), FutureWarning,
                          stacklevel=3)

    def __getitem__(self, key):
        self._warn(key)
        return dict.__getitem__(self, key)

    def get(self, key, default=None):
        self._warn(key)
        return dict.get(self, key, default)


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


def _resolve_focal(focal, spec):
    """The effects to test, as a dict from name to true value, as the R twin's .resolve_focal().

    None takes every coefficient, named by its analysis column (``a_b`` for ``"a:b"``), with its
    true value from the specification. A name or a list of names leaves the true values unknown
    (nan). A mapping gives them. A name may be any fixed effect of the model, ``"(Intercept)"``
    included.
    """
    if focal is None:
        return {_us(k): float(v) for k, v in spec["fixed"]["coefficients"].items()}
    if isinstance(focal, str):
        focal = [focal]
    if isinstance(focal, Mapping):
        if not all(isinstance(k, str) and _is_num(v) for k, v in focal.items()):
            raise ValueError(_FOCAL_FORMS)
        return {k: float(v) for k, v in focal.items()}
    try:
        names = list(focal)
    except TypeError:
        raise ValueError(_FOCAL_FORMS) from None
    if not all(isinstance(n, str) for n in names):
        raise ValueError(_FOCAL_FORMS)
    for n in names:
        if names.count(n) > 1:
            raise ValueError("focal names the effect '%s' twice" % n)
    return {n: math.nan for n in names}


def _variance_components(d, groups):
    """The variance components of the declared random terms, as a statsmodels VCSpec.

    Each grouping factor contributes a component for its intercept, ``0 + C(g)`` in statsmodels'
    formula notation, and one for each declared slope, ``0 + C(g):term``. All rows form one group,
    so that components of crossed factors can be fitted together. The components are named
    ``g`` and ``g:term`` and ordered by name, as MixedLM.from_formula() orders those of its
    `vc_formula`. Their matrices are dense or sparse by the size of the design (see
    _SPARSE_FROM).
    """
    import numpy as np
    import pandas as pd
    from scipy import sparse
    from statsmodels.regression.mixed_linear_model import VCSpec

    n = len(d)
    terms = []
    for g, slopes, _bar in groups:
        codes, levels = pd.factorize(d[g], sort=True)
        terms.append((g, codes, len(levels), np.ones(n)))
        for s in slopes:
            terms.append(("%s:%s" % (g, s), codes, len(levels), d[s].to_numpy(dtype=float)))
    held_sparse = n * sum(t[2] for t in terms) >= _SPARSE_FROM
    rows = np.arange(n)
    mats = {}
    for name, codes, k, values in terms:
        if held_sparse:
            mats[name] = sparse.csr_matrix((values, (rows, codes)), shape=(n, k))
        else:
            z = np.zeros((n, k))
            z[rows, codes] = values
            mats[name] = z
    names = sorted(mats)
    return VCSpec(names, [[[str(j) for j in range(mats[k].shape[1])]] for k in names],
                  [[mats[k]] for k in names])


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


def _power_mixed_replicate(seed, spec, fixed, groups, alternative="two.sided"):
    """One mixed-model replicate: simulate at `seed`, build the modelling data and fit MixedLM.

    `fixed` and `groups` are the model's terms, from _model_terms(). Returns a dict with the
    fitted coefficients' names, a dict of their estimates, a dict of their Wald z p-values and
    the fit's flags, ``singular``, ``warned`` and ``converged``, or None when the response has no
    logarithm or no fit returned. The p-values are two-sided, or one-sided for ``"greater"`` or
    ``"less"``. A coefficient that the fit leaves without a finite estimate or p-value is skipped
    by the caller.
    """
    import numpy as np
    import pandas as pd
    from statsmodels.regression.mixed_linear_model import MixedLM

    # model_data() takes both lognormal families to the log scale, where their coefficients are.
    d = _model_data(spec, _simulate(dict(spec, seed=seed)))
    y = d[".y"].to_numpy(dtype=float)
    # A response rounded down to 0, or to the shift, has no logarithm. R's log() gives -Inf for it
    # and lmer() refuses the fit, so the replicate returns no estimate here too.
    if not np.isfinite(y).all():
        return None
    n = len(d)
    exog = pd.DataFrame({"(Intercept)": np.ones(n), **{c: d[c].to_numpy(dtype=float)
                                                        for c in fixed}})
    # Building the model is kept out of the fits' error handling, so that a fault there stops the
    # call where it used to leave every replicate without a fit and power undefined.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model = MixedLM(y, exog, np.ones(n), exog_vc=_variance_components(d, groups))
    fit = _fit_mixed(model)
    if fit is None:
        return None
    m, notices = fit
    names = list(exog.columns)
    est = {k: float(m.fe_params[k]) for k in names}
    if alternative == "two.sided":
        p = {k: float(m.pvalues[k]) for k in names}
    else:
        from scipy.stats import norm
        # statsmodels gives the two-sided p-value as twice the normal tail beyond |z|. The
        # one-sided p-value is the tail on the tested side of the same z, the upper for "greater"
        # and the lower for "less", which is the upper tail beyond -z. Taking both from norm.sf
        # makes it exactly half the two-sided one when z lies in the tested direction.
        z = {k: float(m.tvalues[k]) for k in names}
        p = {k: float(norm.sf(z[k] if alternative == "greater" else -z[k])) for k in names}
    # lme4's isSingular(): a variance component whose standard deviation is below 1e-4 of the
    # residual one. For independent components that ratio is lme4's theta, and on Powell fits the
    # test agreed with lme4 in each of 24 replicates. It is squared here, so that no variance needs
    # a square root.
    singular = any(v < 1e-8 * m.scale for v in m.vcomp)
    # See _HESSIAN_NOTICE: in a singular fit, the component at zero explains that notice.
    warned = any(not (singular and note.startswith(_HESSIAN_NOTICE)) for note in notices)
    return {"names": names, "est": est, "p": p, "singular": singular, "warned": warned,
            "converged": bool(m.converged) and not singular and not warned}


def _effect_summary(fits, name, beta, alpha, withheld, no_sign):
    """One focal effect's entry in `effects`, over the replicates that returned a fit.

    Power is taken over the fits that give the effect a finite estimate and p-value, as the R twin
    takes it over the fits whose p-value is not NA. `withheld` withholds the mean estimate and
    Type M, which compare the response scale with a link scale, and `no_sign` withholds Type S as
    well, for an interaction on such a scale.
    """
    pairs = [(r["est"][name], r["p"][name]) for r in fits if name in r["est"]
             and math.isfinite(r["est"][name]) and math.isfinite(r["p"][name])]
    est = [e for e, _ in pairs]
    pv = [p for _, p in pairs]
    # Under a one-sided test, the replicates significant in the tested direction.
    sig = [i for i, p in enumerate(pv) if p < alpha]
    # See _power_impl: a zero or unknown true effect leaves both design-analysis quantities
    # undefined, and their Monte Carlo errors with them. Type M's is withheld with Type M.
    usable = bool(sig) and not math.isnan(beta) and beta != 0
    signs = usable and not no_sign
    n_wrong = sum(1 for i in sig if (est[i] > 0) != (beta > 0)) if signs else 0
    ratios = [abs(est[i]) / abs(beta) for i in sig] if usable and not withheld else []
    type_s_lo, type_s_hi = (_type_s_interval(n_wrong, len(sig), beta) if signs
                            else (math.nan, math.nan))
    return {
        **_rate_with_error(len(sig), len(pv), "power"),
        "n_significant": len(sig), "true_effect": beta,
        "mean_estimate": statistics.mean(est) if est and not withheld else math.nan,
        "type_s": n_wrong / len(sig) if signs else math.nan,
        "type_s_lo": type_s_lo, "type_s_hi": type_s_hi,
        "type_m": statistics.mean(ratios) if ratios else math.nan,
        "type_m_mcse": _type_m_mcse(ratios),
    }


def power_mixed(spec, n_sims=50, alpha=0.05, workers=1, alternative="two.sided", focal=None):
    """Mixed-effects simulation-based power in Python, for the model the specification implies,
    via statsmodels MixedLM.

    The model comes from the specification, as in the R twin: `model_formula` writes it and
    `model_data` builds its columns. Every coefficient is a fixed effect, an interaction ``"a:b"``
    entering as its product column ``a_b``. Every grouping factor that the specification gives
    random effects contributes the random intercept and the random slopes it declares, an extra
    grouping factor such as a site included, and nothing else. Each random effect is an
    independent variance component, ``0 + C(g)`` for an intercept and ``0 + C(g):term`` for a
    slope in statsmodels' notation, fitted over a single group so that crossed factors are fitted
    together. The response is analysed as ``log(y)`` for the ``lognormal`` family and as
    ``log(y - shift)`` for ``shifted_lognormal``, whose coefficients are on the log scale. Every
    other family is analysed on its own scale. A replicate whose response rounds to 0, or to the
    shift, has no logarithm and returns no estimate, as in the R twin.

    Up to 0.3, the model was the same whatever the specification declared, ``yv ~ cc`` with
    by-subject and by-item intercept and slope components, and it tested the first contrast of a
    single within factor crossed with items, refusing every other design. Any design that declares
    a random effect is now fitted, and one without random effects is refused, since MixedLM then
    has no variance component to fit. The R package's ``power_mixed()`` fits such a design with
    ``lm()`` or ``glm()``, and `power` analyses the two-group Gaussian design.

    Where the specification correlates a grouping factor's intercept and slopes, the independent
    components leave the correlations out, and one warning per call names the grouping factors
    concerned. The R package's ``power_mixed()`` estimates them, with ``|`` in its formula.

    `focal` chooses the effects to test, as in the R twin. By default it is every coefficient,
    with the true value the specification gives it. A list of names leaves the true values
    unknown, which leaves Type S and Type M undefined, and a dict from names to values gives the
    true values to compare with, which may differ from those simulated. A name is a column of the
    model, so ``"a_b"`` for the interaction ``"a:b"``, or ``"(Intercept)"``. A name that the model
    lacks draws a warning in the R twin's words, and its results are empty. A coefficient the
    specification does not list is not in the model, so a null effect is planned with a
    coefficient of 0.

    For ``bernoulli``, ``poisson``, ``ordinal`` and ``beta``, the response's own scale is not the
    one the coefficients are written on, which is the logit scale, or the log scale for
    ``poisson``. The estimate of a linear model is then a difference on the response scale, such
    as a difference in the probability of a correct response, and comparing it with the true
    coefficient says nothing about exaggeration. For these families, `power_mixed` warns and
    returns `mean_estimate` and `type_m` as `nan` for every effect. `type_s` compares signs alone,
    which the link preserves for a main effect. It does not preserve the sign of an interaction,
    since a bounded response compresses differences near the ends of its range, so the `type_s`
    of an interaction is `nan` too, as in the R twin. Power is the rate at which the linear
    model's test rejects. The R package's ``power_mixed()`` fits ``bernoulli`` and ``poisson``
    responses with ``lme4::glmer()`` on their link scale, and its ``generate_design_analysis()``
    writes a Bayesian design analysis for any family.

    Every replicate is analysed with this one model, singular or not. A preregistered
    simplification rule, which moves to a simpler model whenever a fit is singular, has power and
    a Type I error rate of its own. It is simulated only in the R package, whose ``power_mixed()``
    and ``precision_design()`` take a ladder of models as their ``formula``.

    This function has no joint test of several coefficients, such as the two contrasts of a
    three-level factor or the two coefficients of its interaction with a two-level factor. The
    power of such a multi-degree-of-freedom test differs from that of either coefficient. The R
    package's ``power_mixed()`` plans it when its ``focal`` names sets of coefficients, as in
    ``list(interaction = c("A1_B", "A2_B"))``.

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
    the effect now matches the one from lme4's uncorrelated model, ``(1 + cond || subject) +
    (1 + cond || item)``, to four decimal places in 11. The same 6 replicates are significant in
    both. A replicate whose Powell fit raises is fitted again with the default chain and counts
    as a fit with a warning. The random-effects design is held in dense matrices when it is small
    and in sparse ones from 100,000 entries, rows times random coefficients, which is quicker for
    a large design and reaches the same optimum.

    The fits take time, and more of it as the model grows. On a Windows laptop, with NumPy held to
    one thread, a replicate of the 12 x 8 design of the Python guide took about 0.3 s and one of
    crossed_mixed_rt (30 subjects, 24 items) about 2 s. One of reading_time_continuous, with 50
    subjects, 40 items, seven coefficients and four variance components, took about 7 s.
    `workers` spreads the replicates over processes.

    The p-values are Wald z tests. A Wald z treats the estimate over its standard error as normal
    and so ignores the uncertainty in the variance components. With few subjects or items, it
    rejects more readily than a test with Satterthwaite's degrees of freedom (Luke, 2017,
    doi:10.3758/s13428-016-0809-y), which is the test of the R package's ``power_mixed()``. In
    small designs, the power here therefore tends to exceed R's on the same data. On the 12 x 8
    design of the documentation, 6 of 12 replicates are significant here and 4 of 12 in R.
    With 8 subjects, 6 items, random intercepts alone and no true effect, the test here rejected
    the null hypothesis in 0.070 of 300 replicates (Monte Carlo standard error 0.015). On the same
    data, the R package's Satterthwaite test of the same model rejected in 0.057 (0.013).
    tools/calibration/wald_null_rate.py measures both.

    A directional prediction calls for a one-sided test, which `alternative` selects for every
    focal effect. ``"greater"`` tests whether the coefficient is positive and ``"less"`` whether
    it is negative. The one-sided p-value is the normal tail on the tested side of the same Wald
    z, so it is half the two-sided p-value when the estimate lies in the tested direction. A
    one-sided test at `alpha` therefore rejects in exactly the replicates in which the two-sided
    test at ``2 * alpha`` rejects with the estimate in that direction. Type S and Type M are taken
    over the replicates significant in the tested direction, and `backend` names the direction. A
    two-sided result at ``2 * alpha`` from an earlier version gives the one-sided power at
    `alpha`, in the direction of the true effect, as ``power * (1 - type_s)``.

    Parameters
    ----------
    spec : dict, str or os.PathLike
        A design specification (dict or path to a JSON file) that declares at least one random
        effect.
    n_sims : int, optional
        Number of Monte Carlo replicates (default 50, smaller than `power`'s 1000 because
        each replicate fits a mixed model).
    alpha : float, optional
        Significance level of the test (default 0.05), two-sided unless `alternative` says
        otherwise.
    workers : int, optional
        Number of local worker processes over which to spread the replicates (default 1,
        serial). The replicate seeds are derived once from the specification's seed, so any worker
        count returns results identical to a serial run.
    alternative : {"two.sided", "greater", "less"}, optional
        The alternative hypothesis of every focal effect's test (default ``"two.sided"``), as in
        `power`.
    focal : None, str, list of str or dict, optional
        The effects to test (default None, every coefficient, with the specification's true
        values). Names alone leave the true values unknown, and a dict from names to numbers gives
        them.

    Returns
    -------
    dict
        Keys: `backend` (the estimator and its test), `n_sims`, `alpha`, the fit counts
        `n_attempted`, `n_returned`, `n_converged`, `n_singular` and `n_warning`, and `effects`, a
        dict with one entry per focal effect, in the order of `focal`. Each entry is a dict with
        keys `power`, `power_mcse`, `power_lo`, `power_hi`, `n_significant`, `true_effect`,
        `mean_estimate`, `type_s`, `type_s_lo`, `type_s_hi`, `type_m` and `type_m_mcse`, as in
        `power`. ``pandas.DataFrame.from_dict(result["effects"], orient="index")`` makes a table
        of them.

        A design with a single within factor whose first contrast is a focal effect also has that
        effect's twelve quantities at the top level, where earlier versions returned them. Those
        keys are deprecated and will be removed in the next release, and reading one raises a
        FutureWarning that names its place in `effects`.

        The fit counts mean what they mean in the R twin. `n_attempted` is `n_sims`, and
        `n_returned` counts the replicates that returned a fit. Each effect's `power` is the
        proportion of significant results among the fits that give it a finite estimate and
        p-value, not among `n_sims`. `n_singular` counts fits with a variance component whose
        standard deviation is below 1e-4 of the residual one, the test of lme4's
        ``isSingular()``. `n_warning` counts fits whose optimisers did not converge, whose Powell
        fit raised or whose Hessian was not positive definite. A component at zero can fail
        statsmodels' check of the Hessian by itself, so that notice counts only against a fit that
        is not singular. `n_converged` counts the fits that are neither singular nor warned.
        statsmodels' notice that a fit may be on the boundary counts towards none of them, since
        it comes with every fit that has a variance component below 0.01. Singular fits and fits
        with a warning stay in `power`, as in the R twin. `power_mcse`, `power_lo` and `power_hi`
        are taken over the same fits. `type_s` and `type_m` are `nan` when no replicate reached
        significance and when the true effect is zero or unknown. `mean_estimate` and `type_m` are
        also `nan` for a ``bernoulli``, ``poisson``, ``ordinal`` or ``beta`` response, whose true
        effect is on a link scale, and so is the `type_s` of an interaction. The interval of
        `type_s` and the standard error of `type_m` are `nan` wherever their rate is, and
        `type_m_mcse` also when a single replicate reached significance.

    Raises
    ------
    ValueError
        If `alternative` is not ``"two.sided"``, ``"greater"`` or ``"less"``, if `focal` is none
        of its forms or names an effect twice, or if there is no effect to test, as for a
        specification without coefficients and no `focal`.
    NotImplementedError
        If the specification declares no random effect.

    Warns
    -----
    UserWarning
        Once per call, when the specification correlates the random effects of a grouping factor,
        which the independent components leave out.
    UserWarning
        Once per call, for a ``bernoulli``, ``poisson``, ``ordinal`` or ``beta`` response,
        saying that the linear model is fitted on the response's own scale and that the mean
        estimate and Type M are withheld. The R twin raises the same words for the families it
        fits in that way.
    UserWarning
        When a focal effect never appeared in a fitted model, in the R twin's words.
    FutureWarning
        When a deprecated top-level key of the result is read.

    Notes
    -----
    Requires the `mixed` extra (`statsmodels`, `pandas` and `scipy`, imported lazily in each
    worker process when parallel).
    """
    workers = _check_workers(workers)
    alternative = _check_alternative(alternative)
    import pandas as _pd  # noqa: F401  fail fast before simulating
    import scipy.sparse as _sparse  # noqa: F401
    import statsmodels.regression.mixed_linear_model as _mlm  # noqa: F401

    spec = _as_spec(spec)
    fixed, groups = _model_terms(spec)
    if not groups:
        raise NotImplementedError(_NO_RANDOM)
    true = _resolve_focal(focal, spec)
    if not true:
        raise ValueError(_NOTHING_TO_TEST)
    fam = spec["response"]["family"]
    # The warnings are raised here, once per call, since the replicates may run in other
    # processes. A group's random effects are correlated when the maximal structure writes `|`
    # between its intercept and slopes.
    correlated = ["'%s'" % g for g, slopes, bar in groups if slopes and bar == "|"]
    if correlated:
        warnings.warn(_UNCORRELATED % _name_list(correlated), UserWarning, stacklevel=2)
    # A linear model of a 0/1, count, rating or proportion response estimates a difference on the
    # response scale, which used to be set against the true value on the link scale. On the R
    # twin's crossed accuracy design, such a model gave Type M 0.21 where a logistic one gives
    # 1.12. Type S survives for a main effect, since a monotone link keeps its sign but not its
    # size. It does not survive for an interaction, a difference of differences, whose sign a
    # bounded response can reverse.
    withheld = fam in _LINK_SCALE
    if withheld:
        warnings.warn(_LINEAR_ON_LINK % (fam, _LINK_SCALE[fam]), UserWarning, stacklevel=2)
    interactions = {_us(k) for k in _interaction_keys(spec)} if withheld else set()

    rep = functools.partial(_power_mixed_replicate, spec=spec, fixed=fixed, groups=groups,
                            alternative=alternative)
    seeds = replicate_seeds(spec["seed"], n_sims)
    if workers == 1:
        results = _map_replicates(rep, seeds, None)
    else:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(max_workers=workers) as executor:
            results = _map_replicates(rep, seeds, executor)
    fits = [r for r in results if r is not None]
    effects = {name: _effect_summary(fits, name, beta, alpha, withheld, name in interactions)
               for name, beta in true.items()}

    # A focal name that no fit contains, most often one misspelt, would otherwise read as a design
    # that can decide nothing. The R twin warns in these words.
    absent = ["'%s'" % k for k in effects if fits and k not in fits[0]["est"]]
    if absent:
        warnings.warn(_ABSENT_FOCAL % (
            "s" if len(absent) > 1 else "", ", ".join(absent),
            "their" if len(absent) > 1 else "its",
            ", ".join("'%s'" % k for k in fits[0]["names"])), UserWarning, stacklevel=2)

    # The backend names the direction of a one-sided test, as the R twin's `fitter` does.
    tests = ("Wald z tests" if alternative == "two.sided"
             else "one-sided Wald z tests, %s" % alternative)
    result = {
        "backend": "statsmodels MixedLM (crossed variance components, REML; %s)" % tests,
        "n_sims": n_sims, "alpha": alpha,
        # The R twin's fit counts, with its meanings. n_converged used to count every fit that
        # returned, the count now called n_returned.
        "n_attempted": n_sims, "n_returned": len(fits),
        "n_converged": sum(r["converged"] for r in fits),
        "n_singular": sum(r["singular"] for r in fits),
        "n_warning": sum(r["warned"] for r in fits),
        "effects": effects,
    }
    # The quantities of the effect that earlier versions tested, the first contrast of the single
    # within factor, at the top level where they were, for one release.
    within = [f for f in spec["factors"] if f.get("vary_within")]
    col = next(iter(within[0]["contrasts"])) if len(within) == 1 else None
    if col in effects:
        result.update(effects[col])
        return _MixedResult(result, within=col)
    return _MixedResult(result)


# What each power_curve() record keeps of a power() result. A record used to hold power and Type M
# without their Monte Carlo errors, and no Type S.
_CURVE_KEYS = ("power", "power_mcse", "power_lo", "power_hi", "type_s", "type_s_lo", "type_s_hi",
               "type_m", "type_m_mcse", "n_sims", "n_significant")


def power_curve(spec, subject_ns, n_sims=1000, alpha=0.05, workers=1, seeds="independent",
                alternative="two.sided"):
    """Power curve over sample size for a two-group Gaussian design.

    Sweep the number of subjects and compute `power` at each grid point.

    Each grid point draws its own replicates. Point ``j`` of a sweep over ``J`` sizes takes
    ``replicate_seeds(seed, J)[j]`` as its seed, the rule of the R twin's ``sweep_spec()``. Its
    replicates take their seeds from that one, so no two points share their random numbers.
    `solve_curve` fits the points as independent binomial counts, and its interval keeps its
    nominal coverage only when they are. The rule is prefix-stable: appending sizes leaves the
    earlier points' results as they were. A point's seed follows its position in `subject_ns`,
    though, so reordering the sizes, or inserting one, changes the results at every point that
    moves. ``seeds="shared"`` gives every point the specification's own seed, as up to 0.3.0.
    The points then share their noise. Over 1,000 emulated curves for each of twelve settings on
    a tight grid of sample sizes, the nominal 95% interval of `solve_curve` covered the true size
    in 87.8% of solves with shared seeds and in 96.5% with independent ones.

    Parameters
    ----------
    spec : dict, str or os.PathLike
        A two-group Gaussian design specification with one row per subject, as `power`
        requires.
    subject_ns : iterable of int
        Subject counts to evaluate: any whole numbers, such as a range, a list or a NumPy
        array from ``np.arange`` or ``np.linspace``.
    n_sims : int, optional
        Replicates per grid point (default 1000).
    alpha : float, optional
        Significance level (default 0.05), two-sided unless `alternative` says otherwise.
    workers : int, optional
        Number of local worker processes over which to spread the replicates at each grid
        point (default 1, serial). One process pool is started and reused across the whole
        sweep, and any worker count returns results identical to a serial run.
    seeds : {"independent", "shared"}, optional
        How the grid points are seeded (default ``"independent"``). See above.
    alternative : {"two.sided", "greater", "less"}, optional
        The alternative hypothesis of the test at every grid point (default ``"two.sided"``), as
        in `power`.

    Returns
    -------
    list of dict
        One dict per grid point, with keys `n_subject`, `power`, `power_mcse`, `power_lo`,
        `power_hi`, `type_s`, `type_s_lo`, `type_s_hi`, `type_m`, `type_m_mcse`, `n_sims` and
        `n_significant`, each as `power` returns it. `n_subject` is a Python int, whatever type
        the grid held, so the records serialise as JSON. `power_lo` and `power_hi` give the band to
        draw around a plotted curve. The two counts are what `solve_curve` weights the curve by,
        since a power estimate over 1000 replicates should count for more than one over 50.

    Raises
    ------
    ValueError
        If `seeds` is neither ``"independent"`` nor ``"shared"``, or `alternative` is not
        ``"two.sided"``, ``"greater"`` or ``"less"``.

    See Also
    --------
    solve_curve, target_n : solve this curve for the sample size that meets a target power.
    """
    workers = _check_workers(workers)
    if seeds not in ("independent", "shared"):
        raise ValueError("seeds must be 'independent' or 'shared'")
    alternative = _check_alternative(alternative)
    spec = _as_spec(spec)
    subject_ns = list(subject_ns)
    # One seed per grid point, drawn from the specification's seed by the rule its replicates
    # follow. Every point used to keep the specification's seed, so replicate i drew the same
    # numbers at every point. The points of a curve then shared their noise, while solve_curve()
    # fits them as independent binomial counts, and its 95% intervals covered the true sample
    # size in 88% of solves on a tight grid. replicate_seeds() draws its seeds in turn, so
    # appending points leaves the earlier ones alone.
    if seeds == "independent":
        grid_seeds = replicate_seeds(spec["seed"], len(subject_ns))
    else:
        grid_seeds = [spec["seed"]] * len(subject_ns)

    def sweep(executor):
        out = []
        for n, seed in zip(subject_ns, grid_seeds):
            s = copy.deepcopy(spec)
            s["units"]["subject"]["n"] = n
            s["seed"] = seed
            r = _power_impl(s, n_sims, alpha, executor, alternative)
            # _power_impl() has validated n as a whole number, so int() is exact. A NumPy grid
            # would otherwise leave a NumPy scalar here, and json cannot write a NumPy integer.
            out.append({"n_subject": int(n), **{k: r[k] for k in _CURVE_KEYS}})
        return out

    if workers == 1:
        return sweep(None)
    from concurrent.futures import ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=workers) as executor:  # one pool for the sweep
        return sweep(executor)
