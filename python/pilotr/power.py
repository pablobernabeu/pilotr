"""Simulation-based power and design analysis (Type S / Type M).

For each of `n_sims` Monte Carlo replicates we simulate a fresh data set from the *known*
ground-truth spec, fit the analysis model, and record the estimate and p-value. Beyond
classical power (proportion significant) we report Gelman & Carlin's (2014) design-analysis
quantities, computed over the *significant* replicates:

* Type S (sign) error  -- P(estimate has the wrong sign | significant)
* Type M (magnitude)   -- E(|estimate| / |true effect| | significant)  (exaggeration ratio)

v0.1 ships the analytic backend for a two-group Gaussian design (two-sample t-test). The R
package routes crossed mixed-effects designs to lme4/glmmTMB; the Python analysis backends
(statsmodels / pymer4) are on the roadmap.
"""

from __future__ import annotations
import copy, statistics
from .simulate import simulate, load_spec


def power(spec, n_sims=1000, alpha=0.05):
    from scipy import stats  # lazy: only the power demo needs scipy

    if isinstance(spec, str):
        spec = load_spec(spec)
    if spec["response"]["family"] != "gaussian":
        raise NotImplementedError("v0.1 power backend handles the gaussian two-group design.")

    between = [f for f in spec["factors"] if f.get("between")]
    if len(between) != 1 or len(between[0]["levels"]) != 2:
        raise NotImplementedError("v0.1 power backend expects exactly one 2-level between factor.")
    factor = between[0]
    fname = factor["name"]
    lev0, lev1 = factor["levels"]
    col, vals = next(iter(factor["contrasts"].items()))
    true_effect = spec["fixed"]["coefficients"][col] * (vals[1] - vals[0])
    yname = spec["response"]["name"]

    base_seed = spec["seed"]
    estimates, pvals = [], []
    for i in range(n_sims):
        s = copy.deepcopy(spec)
        s["seed"] = base_seed + i
        d = simulate(s)
        g0 = [r[yname] for r in d.rows if r[fname] == lev0]
        g1 = [r[yname] for r in d.rows if r[fname] == lev1]
        t = stats.ttest_ind(g1, g0, equal_var=True)
        estimates.append(statistics.mean(g1) - statistics.mean(g0))
        pvals.append(t.pvalue)

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


def power_mixed(spec, n_sims=50, alpha=0.05):
    """Crossed mixed-effects simulation-based power in Python, via statsmodels MixedLM with
    by-subject and by-item random intercepts and slopes as (independent) variance components.

    IMPORTANT — accuracy caveat (verified, not a bug). statsmodels fits crossed random
    effects as *independent* variance components and, in our tests, substantially overstates
    random-SLOPE variance (e.g. by-subject slope SD ~0.12 estimated vs 0.04 true). This
    inflates the fixed-effect standard error, so for designs with by-subject/by-item random
    SLOPES this backend is markedly CONSERVATIVE: on the flagship crossed RT design it reports
    power ~0.48 versus the R/lme4 reference ~0.73. It still recovers the fixed effect (mean
    estimate ~0.048 vs 0.05 true) and Type S/M correctly, and is reliable for random-INTERCEPT
    designs. Treat it as a conservative lower bound; use the R/lme4 `power_mixed` as the
    reference whenever random slopes (or random-effect correlations) matter. Data *generation*
    is identical across R and Python; this gap is purely in the Python LMM estimator.
    """
    import copy, math, statistics, warnings
    import pandas as pd
    import statsmodels.formula.api as smf

    if isinstance(spec, str):
        spec = load_spec(spec)
    within = [f for f in spec["factors"] if f.get("vary_within")]
    if len(within) != 1:
        raise NotImplementedError("v0.1 Python power_mixed expects exactly one within factor.")
    f = within[0]
    fname = f["name"]
    col, vals = next(iter(f["contrasts"].items()))
    l2c = {f["levels"][i]: vals[i] for i in range(len(f["levels"]))}
    beta = spec["fixed"]["coefficients"][col]
    yname, fam = spec["response"]["name"], spec["response"]["family"]
    shift = spec["response"].get("shift", 0.0)
    base = spec["seed"]
    vcf = {"subj_i": "0 + C(subject)", "subj_s": "0 + C(subject):cc",
           "item_i": "0 + C(item)", "item_s": "0 + C(item):cc"}

    est, pv = [], []
    for i in range(n_sims):
        s = copy.deepcopy(spec); s["seed"] = base + i
        df = pd.DataFrame(simulate(s).rows)
        df["cc"] = df[fname].map(l2c)
        df["yv"] = [math.log(v - shift) for v in df[yname]] if fam == "shifted_lognormal" else list(df[yname])
        df["grp"] = 1
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                m = smf.mixedlm("yv ~ cc", df, groups="grp", vc_formula=vcf).fit()  # REML (default)
            est.append(float(m.fe_params["cc"])); pv.append(float(m.pvalues["cc"]))
        except Exception:
            pass

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


def power_curve(spec, subject_ns, n_sims=1000, alpha=0.05):
    """Sweep the number of subjects and return power (and Type M) at each grid point."""
    import copy
    if isinstance(spec, str):
        spec = load_spec(spec)
    out = []
    for n in subject_ns:
        s = copy.deepcopy(spec)
        s["units"]["subject"]["n"] = n
        r = power(s, n_sims=n_sims, alpha=alpha)
        out.append({"n_subject": n, "power": r["power"], "type_m": r["type_m"]})
    return out
