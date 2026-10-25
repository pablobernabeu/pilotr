"""A sample-size justification, written in Markdown from a power result, mirroring power_report.R.

The documentation calls the value `target_n` returns the number a preregistration quotes. The
justification that goes with it names the analysis and its test, each effect size and where it came
from, the solved sample size with its interval and the power when the null hypothesis is true. A
simulation adds its Monte Carlo error and fit counts, and a reader needs what reproduces the run.
`power_report` writes all of these down.

It works in two steps. `_power_report_record` reads a result, a solved curve and a specification
file into a record of counts and strings, and `_render` writes the record as Markdown. Every rate
is computed in the renderer from two counts, the significant replicates and the replicates behind
them, with `_rate_with_error` and `_type_s_interval` from `pilotr.mcse`. The R twin renders a
record with the same operations and the same format strings, so one record gives the same bytes in
both languages, and ``tools/parity/report_cross.py`` checks this over the records in
``report_cases.json``. Numbers go through the % operator here and ``sprintf()`` there. Both round
a double correctly to the digits asked for, an exact tie going to the even digit, so they write
the same text for the same double.

The R twin also reports a ``precision_design()`` result, on the interval and the region of
practical equivalence, which has no Python counterpart.
"""

from __future__ import annotations

import hashlib, importlib, math, os, platform
from collections.abc import Mapping

from .mcse import _rate_with_error, _type_s_interval
from .validate import _name_list

# The refusals. Those about `solved`, `justification`, `null_result` and `spec_file` are word for
# word those of the R twin.
_X = ("`x` must be the result of power() or power_mixed() at one design. A curve goes to "
      "target_n(), and its result to `solved`.")
_SOLVED = "`solved` must be the result of target_n() or solve_curve()."
_JUSTIFICATION = ("`justification` must be one string, which applies to every focal effect, or "
                  "strings named by focal effect.")
_UNNAMED = "`x` has one effect, without a name, so `justification` must be a single string."
_NULL_KIND = "`null_result` must be a result of the same kind as `x`."
_SPEC_FILE = "`spec_file` must be the path of the specification file the analysis read."

# The headings of the one effect of a result that names none: that of `power`, and that of a
# `power_mixed` result from before its `effects`, which tested the within contrast alone.
_TWO_GROUP = "The difference between the two groups"
_WITHIN = "The within-factor contrast"

_NO_ESTIMATE = "not available, since no replicate returned an estimate for this effect"

_SOLVED_NEED = ("value", "lo", "hi", "target", "level", "transform", "x", "y", "dispersion")
_SOLVED_KEEP = ("target", "level", "value", "lo", "hi", "n", "n_lo", "n_hi", "transform",
                "dispersion", "x", "y", "n_points", "x_min", "x_max")


def power_report(x, solved=None, spec_file=None, justification=None, null_result=None):
    """Write a sample-size justification in Markdown from a power result.

    Turns the result of `power` or `power_mixed`, with the `target_n` result that solved a curve
    for the sample size, into a sample-size justification, ready for a preregistration, a
    registered report or a grant application.

    The sections follow the recommendations of Lakens (2022, Table 4,
    doi:10.1525/collabra.33267) for reporting an a-priori power analysis. The primary analysis
    comes first, with its model, test, alpha and sidedness, then the desired power and the sample
    size solved for it, with its interval. Each focal effect follows, with its true value, the
    justification given for that value and its power. Last comes the power when the null
    hypothesis is true. The simulation is then accounted for as Morris, White and Crowther (2019,
    doi:10.1002/sim.8086) recommend for simulation studies: the replicates simulated and those
    that returned an estimate, the fits that were singular or carried a warning, and the Monte
    Carlo standard error of every rate. The report ends with what a reader needs to reproduce the
    run: the software and its versions, the specification's seed and an MD5 checksum of the
    specification file's bytes. The versions are those loaded when the report is written, so
    write it in the session that ran the analysis.

    Each rate's Monte Carlo standard error and 95% Wilson interval are computed in the report from
    the number of significant replicates and the number of replicates behind the rate, with the
    functions that compute them in the results. Type S has its Wilson interval from the number of
    wrong signs among the significant replicates. The R twin's ``power_report()`` writes the same
    text, byte for byte, from the same counts, which ``tools/parity/report_cross.py`` in the
    source repository checks.

    The rates are those of the design `x` simulated. A sample size solved from a curve need not be
    that design's, so a justification usually pairs the solved size with a run of `power_mixed`
    at that size. `spec_file` names the specification of that run, and the report lists its
    units. A result carries neither its specification nor its seed, so without `spec_file` the
    report says that both are missing.

    `null_result` is the result of the same analysis with every focal effect set to zero. Its
    rejection rate is the Type I error rate of the test, which for a Wald z with few subjects or
    items can exceed alpha.

    The R twin also reports a ``precision_design()`` result, on the interval and the region of
    practical equivalence, which has no Python counterpart.

    Parameters
    ----------
    x : dict
        The result of `power` or `power_mixed` at one design. A curve goes to `target_n`, and its
        result to `solved`. A `power_mixed` result is read from its `effects` where it has them.
    solved : dict, optional
        The result of `target_n`, or of `solve_curve` for an axis other than sample size, which
        gives the desired power and the value solved for it.
    spec_file : str or os.PathLike, optional
        The path of the specification file the analysis read. The report gives the MD5 checksum
        of its bytes, its seed and its units.
    justification : str or dict, optional
        The justification of each focal effect's size: one string for every focal effect, or a
        dict from focal effect names to strings. An effect without one is reported as having
        none.
    null_result : dict, optional
        The result of the same analysis with every focal effect set to zero. A result whose focal
        effects are not zero is refused.

    Returns
    -------
    str
        The Markdown, ending in a newline. The same arguments always give the same text.

    Raises
    ------
    ValueError
        If `x` or `null_result` is not a result at one design, if `solved` is not the result of
        `target_n` or `solve_curve`, if `justification` is neither a string nor a dict of strings
        naming focal effects, if `null_result` has a focal effect that is not zero, or if
        `spec_file` names no file.

    Examples
    --------
    >>> from pilotr import power, power_report
    >>> spec = {"name": "d", "seed": 1, "units": {"subject": {"n": 64}},
    ...         "factors": [{"name": "group", "levels": ["a", "b"],
    ...                      "contrasts": {"effect": [-0.5, 0.5]}, "between": "subject"}],
    ...         "fixed": {"intercept": 100, "coefficients": {"effect": 5}},
    ...         "response": {"family": "gaussian", "name": "score", "sigma": 10}}
    >>> report = power_report(power(spec, n_sims=50),
    ...                       justification="Half a standard deviation.")  # doctest: +SKIP
    >>> print(report)  # doctest: +SKIP
    """
    return _render(_power_report_record(x, solved=solved, spec_file=spec_file,
                                        justification=justification, null_result=null_result))


def _missing(v):
    """A missing value: None, as JSON's null reads, or nan, as a result holds it."""
    return v is None or (isinstance(v, float) and math.isnan(v))


def _value(v):
    return None if _missing(v) else v


def _kind(x):
    """The kind of a result, "power" here, where the R twin also reads a precision result."""
    if (isinstance(x, Mapping) and "n_sims" in x and "alpha" in x
            and ("effects" in x or ("power" in x and "n_significant" in x))):
        return "power"
    raise ValueError(_X)


def _tested(p, k, base):
    """The number of replicates behind an effect's power, as in the R twin's .report_tested().

    `power_mixed` takes an effect's power over the fits that gave it an estimate, which can be
    fewer than `n_returned`, and the result keeps the rate and the count of significant
    replicates. The denominator is their ratio, which division recovers exactly. With no
    significant replicate it is `n_returned`, and `power` takes power over every replicate.
    """
    if _missing(p):
        return 0
    if p > 0:
        return round(k / p)
    return base


def _power_effects(x):
    """One record per focal effect. A result with `effects` is read from them alone, so that the
    deprecated scalar keys beside them, which warn when read, are left alone."""
    base = x.get("n_returned")
    if base is None:
        base = x["n_sims"]
    if "effects" in x:
        items = list(x["effects"].items())
        label = None
    else:
        items = [(None, x)]
        label = _WITHIN if "backend" in x else _TWO_GROUP
    out = []
    for name, e in items:
        k = e["n_significant"]
        out.append({"name": name, "label": label if name is None else None,
                    "coefficients": None,
                    "true_effect": _value(e.get("true_effect")),
                    "n_significant": k,
                    "n_tested": _tested(e.get("power"), k, base),
                    "type_s": _value(e.get("type_s")),
                    "type_m": _value(e.get("type_m")),
                    "type_m_mcse": _value(e.get("type_m_mcse"))})
    return out


def _justify(justification, effects):
    """Each focal effect's justification, or None where none was given."""
    if justification is None:
        return [None] * len(effects)
    if isinstance(justification, str):
        return [justification] * len(effects)
    if (not isinstance(justification, Mapping) or not justification
            or not all(isinstance(k, str) and k and isinstance(v, str)
                       for k, v in justification.items())):
        raise ValueError(_JUSTIFICATION)
    names = [e["name"] for e in effects]
    if any(n is None for n in names):
        raise ValueError(_UNNAMED)
    unknown = [k for k in justification if k not in names]
    if unknown:
        raise ValueError("`justification` names '%s', which is not a focal effect of `x`. Its "
                         "focal effects are %s." % (unknown[0],
                                                    _name_list("'%s'" % n for n in names)))
    return [justification.get(n) for n in names]


def _fitter(x):
    """The estimator and test a result names. `power` names none and always fits a pooled t."""
    backend = x.get("backend")
    if backend is not None:
        return backend
    if "effects" in x or "n_returned" in x:
        return None
    return "scipy.stats.ttest_ind, two-sample t with pooled variance"


def _alternative(x):
    """`power` records its alternative, and `power_mixed` names a one-sided test in `backend`.
    Earlier results were two-sided."""
    alternative = x.get("alternative")
    if alternative is not None:
        return alternative
    backend = x.get("backend")
    if not backend or "one-sided" not in backend:
        return "two.sided"
    return "greater" if "greater" in backend else "less"


def _fits(x):
    return {"n_sims": x["n_sims"], "n_attempted": x.get("n_attempted"),
            "n_returned": x.get("n_returned"), "n_converged": x.get("n_converged"),
            "n_singular": x.get("n_singular"), "n_warning": x.get("n_warning")}


def _null(null_result):
    """The null condition's counts, after checking that its focal effects are zero."""
    if null_result is None:
        return None
    try:
        _kind(null_result)
    except ValueError:
        raise ValueError(_NULL_KIND) from None
    effects = _power_effects(null_result)
    for e in effects:
        v = e["true_effect"]
        if v is not None and v != 0:
            raise ValueError("`null_result` must come from a design in which every focal effect "
                             "is zero. Its effect '%s' has a true value of %s."
                             % (e["label"] if e["name"] is None else e["name"], _g(v)))
    return {"alpha": null_result.get("alpha"), "alternative": _alternative(null_result),
            "n_sims": null_result["n_sims"], "n_returned": null_result.get("n_returned"),
            "effects": [{k: e[k] for k in ("name", "label", "n_significant", "n_tested")}
                        for e in effects]}


def _solved(solved):
    if solved is None:
        return None
    if not isinstance(solved, Mapping) or any(k not in solved for k in _SOLVED_NEED):
        raise ValueError(_SOLVED)
    return {k: solved.get(k) for k in _SOLVED_KEEP}


def _version(module, label):
    try:
        return "%s %s" % (label, importlib.import_module(module).__version__)
    except ImportError:
        return ""


def _software(fitter):
    """"pilotr 0.3.0, Python 3.13.1", with SciPy or statsmodels when the fitter used them."""
    from . import __version__
    f = fitter or ""
    parts = ["pilotr %s" % __version__, "Python %s" % platform.python_version()]
    if "statsmodels" in f:
        parts.append(_version("statsmodels", "statsmodels"))
    if "scipy" in f:
        parts.append(_version("scipy", "SciPy"))
    return ", ".join(p for p in parts if p)


def _reproducibility(spec_file, fitter):
    """The specification file's name, the MD5 checksum of its bytes, its seed and its units.

    A result records none of them. The checksum is taken over the file's bytes as they stand, with
    no decoding, as the R twin's ``tools::md5sum()`` takes it, so the same file gives the same
    checksum in both.
    """
    out = {"software": _software(fitter), "spec_file": None, "spec_md5": None, "seed": None,
           "units": None}
    if spec_file is None:
        return out
    if not isinstance(spec_file, (str, os.PathLike)):
        raise ValueError(_SPEC_FILE)
    path = os.fspath(spec_file)
    if not os.path.isfile(path):
        raise ValueError("%s No file is found at '%s'." % (_SPEC_FILE, path))
    from .simulate import load_spec
    spec = load_spec(path, validate=False)
    with open(path, "rb") as f:
        data = f.read()
    out.update(spec_file=os.path.basename(path),
               spec_md5=hashlib.md5(data, usedforsecurity=False).hexdigest(),
               seed=spec.get("seed"),
               units=[{"name": u, "n": v["n"]} for u, v in spec["units"].items()])
    return out


def _power_report_record(x, solved=None, spec_file=None, justification=None, null_result=None):
    """The record a report is rendered from: counts and strings, laid out as
    tools/parity/report_cases.json writes them, with None for anything missing."""
    kind = _kind(x)
    effects = _power_effects(x)
    for e, j in zip(effects, _justify(justification, effects)):
        e["justification"] = j
    fitter = _fitter(x)
    analysis = {"fitter": fitter, "alpha": x.get("alpha"), "alternative": _alternative(x),
                "joint": False, "rungs": None, "fallback": None}
    return {"kind": kind, "analysis": analysis, "effects": effects, "solved": _solved(solved),
            "null": _null(null_result), "fits": _fits(x),
            "reproducibility": _reproducibility(spec_file, fitter)}


# --- Rendering. -----------------------------------------------------------------------------------
#
# Every line below has a twin in power_report.R, which also renders a precision_design() result
# in lines that exist in R alone, as the top of this file says.

def _f3(v):
    return "%.3f" % v


def _g(v):
    return "%g" % v


def _int(v):
    return "%.0f" % v


def _name(e):
    return e.get("label") if _missing(e.get("name")) else "`%s`" % e["name"]


def _rate(count, n):
    """A rate with its Monte Carlo standard error and Wilson interval, computed from its counts."""
    if _missing(n) or n <= 0 or _missing(count):
        return _NO_ESTIMATE
    r = _rate_with_error(count, n, "p")
    return ("%s (Monte Carlo standard error %s; 95%% Wilson interval %s to %s), from %s of %s "
            "replicates" % (_f3(r["p"]), _f3(r["p_mcse"]), _f3(r["p_lo"]), _f3(r["p_hi"]),
                            _int(count), _int(n)))


def _sidedness(alternative, joint):
    if joint is True:
        return "a joint test of each set"
    if _missing(alternative) or alternative == "two.sided":
        return "two-sided"
    return "one-sided (%s)" % alternative


def _units(units):
    return ", ".join("%s %s" % (u["name"], _int(u["n"])) for u in units)


def _section_analysis(rec):
    a = rec["analysis"]
    fitter = "not recorded" if _missing(a.get("fitter")) else a["fitter"]
    level = ("not recorded" if _missing(a.get("alpha")) else
             "alpha = %s, %s" % (_g(a["alpha"]), _sidedness(a.get("alternative"), a.get("joint"))))
    out = ["## Primary analysis", "", "- Model and test: " + fitter,
           "- Significance level: " + level]
    rungs = a.get("rungs")
    if rungs:
        rule = ", then ".join(r["name"] for r in rungs)
        fallback = a.get("fallback")
        if not _missing(fallback):
            rule += ", moving on from " + ("singular fits and fits with warnings"
                                           if fallback == "warning" else "singular fits")
        out.append("- Simplification rule: " + rule)
    units = rec["reproducibility"].get("units")
    if units:
        out.append("- Simulated units: " + _units(units))
    return out


def _section_solved(rec):
    s = rec.get("solved")
    if not s:
        return []
    whole = not _missing(s.get("n"))
    rate = ("Desired probability of a meaningful decision" if s["y"] == "p_meaningful"
            else "Desired power")
    level = "%s%%" % _g(s["level"] * 100)
    out = ["## Desired power and sample size" if whole else "## Desired power and the solved value",
           "", "- %s: %s" % (rate, _g(s["target"]))]
    if whole:
        out.append("- Sample size: %s (%s confidence interval %s to %s), the exact solution %s "
                   "rounded up" % (_int(s["n"]), level, _int(s["n_lo"]), _int(s["n_hi"]),
                                   "%.2f" % s["value"]))
    else:
        out.append("- Solved value of `%s`: %s (%s confidence interval %s to %s)"
                   % (s["x"], _g(s["value"]), level, _g(s["lo"]), _g(s["hi"])))
    scale = {"sqrt": "the square root of `%s`" % s["x"],
             "log": "the logarithm of `%s`" % s["x"]}.get(s["transform"], "`%s`" % s["x"])
    span = ("" if _missing(s.get("n_points")) else ", over %s points from %s to %s"
            % (_int(s["n_points"]), _g(s["x_min"]), _g(s["x_max"])))
    out.append("- Fit: a probit regression of `%s` on %s%s" % (s["y"], scale, span))
    out.append("- Heterogeneity factor: %s, the dispersion that scales the variance behind the "
               "interval, which is 1 where the curve fits the probit model"
               % ("%.2f" % s["dispersion"]))
    return out


def _effect_power(e, joint):
    out = ["### " + _name(e), ""]
    if e.get("coefficients"):
        out.append("- Coefficients tested jointly: " + ", ".join(e["coefficients"]))
    known = not _missing(e.get("true_effect"))
    if joint is not True:
        out.append("- True value: " + (_g(e["true_effect"]) if known
                                       else "not given, so Type S and Type M are not defined"))
    out.append("- Justification: " + ("not given" if _missing(e.get("justification"))
                                      else e["justification"]))
    out.append("- Power: " + _rate(e.get("n_significant"), e.get("n_tested")))
    if joint is True or not known:
        return out
    if _missing(e.get("type_s")):
        type_s = "not available"
    else:
        k = e["n_significant"]
        lo, hi = _type_s_interval(round(e["type_s"] * k), k, e["true_effect"])
        type_s = (_f3(e["type_s"]) if math.isnan(lo) else
                  "%s (95%% Wilson interval %s to %s)" % (_f3(e["type_s"]), _f3(lo), _f3(hi)))
    if _missing(e.get("type_m")):
        type_m = "not available"
    elif _missing(e.get("type_m_mcse")):
        type_m = ("%s (too few significant replicates for a Monte Carlo standard error)"
                  % _f3(e["type_m"]))
    else:
        type_m = "%s (Monte Carlo standard error %s)" % (_f3(e["type_m"]), _f3(e["type_m_mcse"]))
    return out + ["- Type S: " + type_s, "- Type M: " + type_m]


def _section_effects(rec):
    out = ["## Focal effects", "",
           "The rates below are those of the simulated design, each over the replicates that "
           "returned an estimate for the effect."]
    for e in rec["effects"]:
        out += [""] + _effect_power(e, rec["analysis"].get("joint"))
    return out


def _section_null(rec):
    z = rec.get("null")
    if not z:
        return []
    if _missing(z.get("alpha")):
        intro = "With every focal effect set to zero:"
    else:
        intro = ("With every focal effect set to zero, at alpha = %s, %s:"
                 % (_g(z["alpha"]), _sidedness(z.get("alternative"),
                                               rec["analysis"].get("joint"))))
    lines = ["- %s: %s" % (_name(e), _rate(e.get("n_significant"), e.get("n_tested")))
             for e in z["effects"]]
    return ["## Power when the null hypothesis is true", "", intro, ""] + lines


def _replicates(n_sims, n_returned):
    if _missing(n_returned):
        return "%s simulated" % _int(n_sims)
    return "%s simulated, %s returned an estimate" % (_int(n_sims), _int(n_returned))


def _section_accounting(rec):
    f = rec["fits"]
    out = ["## Monte Carlo accounting", "",
           "- Replicates: " + _replicates(f["n_sims"], f.get("n_returned"))]
    if not _missing(f.get("n_converged")):
        out.append("- Fits: %s converged cleanly, %s singular, %s with a warning"
                   % (_int(f["n_converged"]), _int(f["n_singular"]), _int(f["n_warning"])))
    rungs = rec["analysis"].get("rungs")
    if rungs:
        out.append("- Analysed at each rung: "
                   + ", ".join("%s %s" % (r["name"], _f3(r["share"])) for r in rungs))
    z = rec.get("null")
    if z:
        out.append("- Null condition: " + _replicates(z["n_sims"], z.get("n_returned")))
    note = ("Each rate's Monte Carlo standard error is sqrt(p(1 - p) / n), over the n replicates "
            "that returned an estimate, and its interval is the 95% Wilson score interval.")
    if rec["analysis"].get("joint") is not True:
        note += (" The standard error of Type M is the standard deviation of the exaggeration "
                 "ratios over the square root of their number.")
    return out + ["", note]


def _section_reproducibility(rec):
    r = rec["reproducibility"]
    out = ["## Reproducibility", "", "- Software: " + r["software"]]
    if _missing(r.get("spec_file")):
        return out + ["- Specification: not given, so neither its seed nor its checksum is "
                      "recorded"]
    return out + ["- Specification: %s, MD5 checksum %s" % (r["spec_file"], r["spec_md5"]),
                  "- Seed: " + ("not given" if _missing(r.get("seed")) else _int(r["seed"]))]


def _render(rec):
    """The record as Markdown: one string, its lines joined by newlines and ending in one."""
    if rec.get("kind", "power") != "power":
        raise ValueError("a %s record is rendered by the R twin alone" % rec.get("kind"))
    lines = ["# Sample-size justification"]
    for section in (_section_analysis(rec), _section_solved(rec), _section_effects(rec),
                    _section_null(rec), _section_accounting(rec),
                    _section_reproducibility(rec)):
        if section:
            lines += [""] + section
    return "\n".join(lines) + "\n"
