"""Simulation-based power and design analysis (Type S / Type M).

For each of the `n_sims` Monte Carlo replicates, we simulate a fresh data set from the
known ground-truth spec, fit the analysis model, and record the estimate and p-value. In
addition to classical power (the proportion of significant replicates), we report the
design-analysis quantities of Gelman and Carlin (2014), computed over the significant
replicates.

* Type S (sign) error: P(estimate has the wrong sign | significant)
* Type M (magnitude): E(|estimate| / |true effect| | significant)  (exaggeration ratio)

The two-group Gaussian design uses a two-sample t-test (`power`). Crossed mixed-effects
designs use a statsmodels MixedLM backend (`power_mixed`), which is conservative for
random-slope designs; the R package's lme4 backend is the reference there.

Every analysis takes a `workers` argument that spreads the replicates over local
processes. Each replicate seeds the shared cross-language RNG from its own index, so the
results are identical to a serial run whatever the worker count. The per-replicate
functions are module-level so that they pickle under the Windows spawn start method.
"""

from __future__ import annotations
import copy, functools, statistics
from .simulate import simulate, load_spec


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


def _power_replicate(seed, spec, fname, lev0, lev1, yname):
    """One two-group replicate: simulate at `seed`, t-test, return (estimate, p-value)."""
    from scipy import stats  # lazy: imported in each worker process on first use

    s = copy.deepcopy(spec)
    s["seed"] = seed
    d = simulate(s)
    g0 = [r[yname] for r in d.rows if r[fname] == lev0]
    g1 = [r[yname] for r in d.rows if r[fname] == lev1]
    t = stats.ttest_ind(g1, g0, equal_var=True)
    return statistics.mean(g1) - statistics.mean(g0), t.pvalue


def power(spec, n_sims=1000, alpha=0.05, workers=1):
    """Simulation-based power and design analysis for a two-group Gaussian design.

    Repeatedly simulate from the specification, apply a two-sample t-test, and report power
    together with the Type S (sign) and Type M (magnitude) errors of Gelman and Carlin
    (2014), computed over the significant replicates.

    Parameters
    ----------
    spec : dict or str
        A two-group Gaussian design specification (dict or path to a JSON file).
    n_sims : int, optional
        Number of Monte Carlo replicates (default 1000).
    alpha : float, optional
        Two-sided significance level (default 0.05).
    workers : int, optional
        Number of local worker processes over which to spread the replicates (default 1,
        serial). Each replicate seeds the shared RNG from its own index, so any worker
        count returns results identical to a serial run.

    Returns
    -------
    dict
        Keys: `n_sims`, `alpha`, `power`, `n_significant`, `true_effect`, `mean_estimate`,
        `type_s`, `type_m`.

    Raises
    ------
    NotImplementedError
        If the design is not a single two-level between-subjects Gaussian factor.

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
    if isinstance(spec, str):
        spec = load_spec(spec)
    if spec["response"]["family"] != "gaussian":
        raise NotImplementedError("The power backend currently handles only the gaussian two-group design.")

    between = [f for f in spec["factors"] if f.get("between")]
    if len(between) != 1 or len(between[0]["levels"]) != 2:
        raise NotImplementedError("The power backend expects exactly one 2-level between factor.")
    factor = between[0]
    fname = factor["name"]
    lev0, lev1 = factor["levels"]
    col, vals = next(iter(factor["contrasts"].items()))
    true_effect = spec["fixed"]["coefficients"][col] * (vals[1] - vals[0])
    yname = spec["response"]["name"]

    base_seed = spec["seed"]
    rep = functools.partial(_power_replicate, spec=spec, fname=fname,
                            lev0=lev0, lev1=lev1, yname=yname)
    results = _map_replicates(rep, [base_seed + i for i in range(n_sims)], executor)
    estimates = [r[0] for r in results]
    pvals = [r[1] for r in results]

    sig = [i for i, p in enumerate(pvals) if p < alpha]
    power_val = len(sig) / n_sims
    if sig:
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


def _power_mixed_replicate(seed, spec, fname, l2c, yname, fam, shift):
    """One mixed-model replicate: simulate at `seed`, fit MixedLM, return (estimate,
    p-value), or None when the fit fails."""
    import math, warnings
    import pandas as pd
    import statsmodels.formula.api as smf

    vcf = {"subj_i": "0 + C(subject)", "subj_s": "0 + C(subject):cc",
           "item_i": "0 + C(item)", "item_s": "0 + C(item):cc"}
    s = copy.deepcopy(spec)
    s["seed"] = seed
    df = pd.DataFrame(simulate(s).rows)
    df["cc"] = df[fname].map(l2c)
    df["yv"] = [math.log(v - shift) for v in df[yname]] if fam == "shifted_lognormal" else list(df[yname])
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

    This is pilotr's own simulation loop over the portable design specification, not a
    wrapper around an existing power package. It covers territory pioneered by simr (Green
    and MacLeod, 2016, doi:10.1111/2041-210x.12504) and mixedpower (Kumle, Vo and Draschkow,
    2021, doi:10.3758/s13428-021-01546-0); pilotr differs in being driven by the portable
    cross-language spec with bit-identical R and Python data, in reporting Type S and Type M
    errors and in built-in parallelisation via `workers` (default 1, serial), which returns
    results identical to a serial run for any worker count.

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
        serial). Each replicate seeds the shared RNG from its own index, so any worker
        count returns results identical to a serial run.

    Returns
    -------
    dict
        Keys: `backend` (the estimator used), `n_sims`, `n_converged` (how many replicates
        the model fit), `alpha`, `power`, `n_significant`, `true_effect`, `mean_estimate`,
        `type_s`, `type_m`. `power` is the proportion of significant results among the
        `n_converged` converged replicates, not among `n_sims`.

    Raises
    ------
    ValueError
        If the spec has no item unit (`power_mixed` requires a crossed design).
    NotImplementedError
        If the design does not have exactly one within-unit factor.

    Notes
    -----
    Requires the `mixed` extra (`statsmodels` and `pandas`, imported lazily in each worker
    process when parallel).
    """
    workers = _check_workers(workers)
    import pandas as _pd  # noqa: F401  fail fast before simulating
    import statsmodels.formula.api as _smf  # noqa: F401

    if isinstance(spec, str):
        spec = load_spec(spec)
    if "item" not in spec["units"]:
        raise ValueError("power_mixed() requires a crossed design with an item unit.")
    within = [f for f in spec["factors"] if f.get("vary_within")]
    if len(within) != 1:
        raise NotImplementedError("The Python power_mixed backend expects exactly one within factor.")
    f = within[0]
    fname = f["name"]
    col, vals = next(iter(f["contrasts"].items()))
    l2c = {f["levels"][i]: vals[i] for i in range(len(f["levels"]))}
    beta = spec["fixed"]["coefficients"][col]
    yname, fam = spec["response"]["name"], spec["response"]["family"]
    shift = spec["response"].get("shift", 0.0)
    base = spec["seed"]

    rep = functools.partial(_power_mixed_replicate, spec=spec, fname=fname, l2c=l2c,
                            yname=yname, fam=fam, shift=shift)
    seeds = [base + i for i in range(n_sims)]
    if workers == 1:
        results = _map_replicates(rep, seeds, None)
    else:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(max_workers=workers) as executor:
            results = _map_replicates(rep, seeds, executor)
    est = [r[0] for r in results if r is not None]
    pv = [r[1] for r in results if r is not None]

    sig = [i for i, p in enumerate(pv) if p < alpha]
    return {
        "backend": "statsmodels MixedLM (crossed variance components, REML)",
        "n_sims": n_sims, "n_converged": len(pv), "alpha": alpha,
        "power": len(sig) / len(pv) if pv else float("nan"),
        "n_significant": len(sig), "true_effect": beta,
        "mean_estimate": statistics.mean(est) if est else float("nan"),
        "type_s": (sum(1 for i in sig if (est[i] > 0) != (beta > 0)) / len(sig)) if sig else float("nan"),
        "type_m": statistics.mean(abs(est[i]) / abs(beta) for i in sig) if sig else float("nan"),
    }


def power_curve(spec, subject_ns, n_sims=1000, alpha=0.05, workers=1):
    """Power curve over sample size for a two-group Gaussian design.

    Sweep the number of subjects and compute `power` at each grid point.

    Parameters
    ----------
    spec : dict or str
        A two-group Gaussian design specification.
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
        One dict per grid point, with keys `n_subject`, `power`, `type_m`.
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
            out.append({"n_subject": n, "power": r["power"], "type_m": r["type_m"]})
        return out

    if workers == 1:
        return sweep(None)
    from concurrent.futures import ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=workers) as executor:  # one pool for the sweep
        return sweep(executor)
