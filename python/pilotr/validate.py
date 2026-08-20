"""Specification validation and version negotiation, mirroring pilotr/R/validate.R.

``load_spec`` was a bare ``json.load``: no validation and no version check. A strict draft-07
schema has shipped at spec/design.schema.json since 0.1, and no code path consulted it, so a
specification with a misspelled field loaded silently and kept the misspelling.

That matters more than a typo usually would, because several of the ways a specification can be
wrong produce plausible data and no error at all:

* A mistyped coefficient key resolves to no column, so that effect is silently set to zero. A
  spec whose focal effect is spelled "cnod" for "cond" generates exactly the data of a null
  design, and reports success.
* ``varies_by`` took anything other than "subject" as item-level, so a per-trial predictor was
  silently given one value per item.
* A response parameter belonging to another family is silently ignored.

It also matters for version negotiation. Once 0.3 features exist, a 0.3 specification opened by a
0.2 implementation, an un-upgraded twin, or a cached browser build produces different and wrong
data while reporting success. A 0.2 implementation has no version check and cannot be fixed
retrospectively, but a specification that uses a 0.3 feature can be made to say so, and every
implementation from 0.3 onwards refuses what it does not understand.
"""

from __future__ import annotations

import math

# The specification version this implementation writes and understands.
SPEC_VERSION = "0.3"

# Response families and the parameters each one uses. Anything else supplied under `response` is
# refused, because a leftover parameter from another family is usually a half-finished edit and
# silently dropping it would hide the mistake.
FAMILY_PARAMS = {
    "gaussian": ("sigma",),
    "lognormal": ("sigma",),
    "shifted_lognormal": ("sigma", "shift"),
    "exgaussian": ("sigma", "beta"),
    "bernoulli": (),
    "poisson": (),
    "ordinal": ("thresholds",),
    "beta": ("phi",),
}

# Families whose response value is rounded when `response.round` is set. For the others the
# outcome is an integer already, so `round` would do nothing and is refused as a likely mistake.
ROUNDING_FAMILIES = ("gaussian", "lognormal", "shifted_lognormal", "exgaussian", "beta")

_KNOWN_TOP = ("spec_version", "name", "seed", "units", "factors", "predictors",
              "fixed", "random", "response")


def _is_str(x) -> bool:
    return isinstance(x, str)


def _is_num(x) -> bool:
    # bool is a subclass of int in Python, and a boolean where a number belongs is a mistake.
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def _is_whole(x) -> bool:
    return _is_num(x) and float(x) == round(float(x))


def _parse_version(v):
    """Read a version as a (major, minor) pair, or None when it is not one.

    A version arrives as whatever JSON produced. A single part is read as a whole version, so
    "1" and the JSON number 1.0 both mean 1.0. The padding is what keeps the two engines agreeing
    about the same file: R renders the number 1.0 as "1" and Python as "1.0", so without it one
    implementation called the specification malformed while the other read it as version 1.
    """
    parts = str(v).split(".")
    if len(parts) == 1:
        parts.append("0")
    if len(parts) < 2:
        return None
    try:
        return int(parts[0]), int(parts[1])
    except ValueError:
        return None


def _features_0_3(spec) -> list:
    """Features introduced in 0.3. A specification using any of them is read differently by a
    0.2 implementation, so it has to declare 0.3 or later."""
    found = []
    for p in spec.get("predictors") or []:
        if not isinstance(p, dict):
            continue
        if p.get("varies_by") == "observation":
            found.append('predictors varies_by "observation"')
        if p.get("dist") is not None:
            found.append("predictors dist")
        rel = p.get("reliability")
        if rel is not None and rel != 1:
            found.append("predictors reliability")
    if (spec.get("response") or {}).get("family") == "exgaussian":
        found.append('the "exgaussian" family')
    for _g, re in (spec.get("random") or {}).items():
        if not isinstance(re, dict):
            continue
        if re.get("correlated") is not None:
            found.append("random correlated")
        if any(":" in k for k in (re.get("slopes") or {})):
            found.append("interaction random slopes")
    return list(dict.fromkeys(found))


def validate_spec(spec, strict: bool = True):
    """Validate a design specification, returning it so the call can be chained.

    Checks the specification against the portable schema and against the cross-field rules the
    schema cannot express, and checks that its declared ``spec_version`` is one this
    implementation understands. Called by ``load_spec`` by default.

    Parameters
    ----------
    spec : dict
        A parsed design specification.
    strict : bool
        Whether an unrecognised field is an error (the default) or a warning. Pass ``False`` to
        load a specification carrying private annotations, accepting that a misspelled field
        will then be ignored in silence.

    Returns
    -------
    dict
        The specification, unchanged.

    Raises
    ------
    ValueError
        If the specification is invalid. All problems found are reported together.
    """
    import warnings

    if not isinstance(spec, dict):
        raise ValueError("a design specification must be a JSON object")

    problems: list[str] = []
    soft: list[str] = []

    def bad(msg):
        problems.append(msg)

    def unknown(msg):
        (problems if strict else soft).append(msg)

    # ---- version ----
    declared = spec.get("spec_version", "0.2")
    dv, sv = _parse_version(declared), _parse_version(SPEC_VERSION)
    if dv is None:
        bad("spec_version '%s' is not of the form 'major.minor'" % declared)
    else:
        # The version as pilotr read it, so that the two engines report the same thing about a
        # JSON number they render differently.
        shown = "%d.%d" % dv
        if dv > sv:
            bad("this specification declares spec_version %s, which is newer than the %s this "
                "version of pilotr understands; please upgrade pilotr" % (shown, SPEC_VERSION))
        used = _features_0_3(spec)
        if used and dv < (0, 3):
            bad('this specification uses %s, which requires spec_version "0.3", but declares %s; '
                "a 0.2 implementation would read it differently and silently generate different "
                "data" % (", ".join(used), shown))

    # ---- top level ----
    for k in spec:
        if k not in _KNOWN_TOP:
            unknown("unknown top-level field '%s'; expected one of %s"
                    % (k, ", ".join(_KNOWN_TOP)))
    for k in ("name", "seed", "units", "fixed", "response"):
        if spec.get(k) is None:
            bad("required top-level field '%s' is missing" % k)
    if spec.get("name") is not None and not _is_str(spec["name"]):
        bad("'name' must be a single string")
    if spec.get("seed") is not None and not _is_whole(spec["seed"]):
        bad("'seed' must be a single whole number")

    # ---- units ----
    units = spec.get("units")
    has_item = False
    if units is not None:
        if not isinstance(units, dict):
            bad("'units' must be an object")
        else:
            for k in units:
                if k not in ("subject", "item"):
                    unknown("unknown unit '%s'; only 'subject' and 'item' exist" % k)
            if units.get("subject") is None:
                bad("'units.subject' is required")
            has_item = units.get("item") is not None
            for nm in ("subject", "item"):
                un = units.get(nm)
                if un is None:
                    continue
                if not isinstance(un, dict):
                    bad("'units.%s' must be an object" % nm)
                    continue
                for k in un:
                    if k not in ("n", "per_subject"):
                        unknown("unknown field 'units.%s.%s'" % (nm, k))
                if not _is_whole(un.get("n")) or un["n"] < 1:
                    bad("'units.%s.n' must be a whole number of at least 1" % nm)
                if un.get("per_subject") is not None:
                    if nm != "item":
                        bad("'per_subject' belongs to 'units.item', not 'units.%s'" % nm)
                    elif not _is_whole(un["per_subject"]) or un["per_subject"] < 1:
                        bad("'units.item.per_subject' must be a whole number of at least 1")
                    elif _is_whole(un.get("n")) and un["per_subject"] > un["n"]:
                        bad("'units.item.per_subject' (%s) cannot exceed the number of items (%s)"
                            % (un["per_subject"], un["n"]))

    # ---- factors ----
    contrast_cols: list[str] = []
    factors = spec.get("factors")
    if factors is not None:
        if not isinstance(factors, list):
            bad("'factors' must be an array of factor objects")
        else:
            for i, f in enumerate(factors, start=1):
                where = "factors[%d]" % i
                if not isinstance(f, dict):
                    bad("%s must be an object" % where)
                    continue
                for k in f:
                    if k not in ("name", "levels", "contrasts", "vary_within", "between"):
                        unknown("unknown field '%s.%s'" % (where, k))
                if not _is_str(f.get("name")):
                    bad("%s.name must be a single string" % where)
                levels = f.get("levels")
                nlev = len(levels) if isinstance(levels, list) else 0
                if not isinstance(levels, list) or nlev < 2 or not all(_is_str(v) for v in levels):
                    bad("%s.levels must be an array of at least two strings" % where)
                contrasts = f.get("contrasts")
                if not isinstance(contrasts, dict) or not contrasts:
                    bad("%s.contrasts must be a non-empty object mapping contrast columns to one "
                        "value per level" % where)
                else:
                    for cn, v in contrasts.items():
                        contrast_cols.append(cn)
                        if not isinstance(v, list) or not all(_is_num(x) for x in v):
                            bad("%s.contrasts.%s must be numeric" % (where, cn))
                        elif nlev >= 2 and len(v) != nlev:
                            bad("%s.contrasts.%s has %d value(s) but the factor has %d level(s)"
                                % (where, cn, len(v), nlev))
                vw = f.get("vary_within")
                if vw is not None:
                    # A single string is accepted where an array belongs. pilotr's own spec_json()
                    # emitted that form before 0.3, because a blanket auto_unbox collapsed every
                    # one-element array, so refusing it would mean refusing files pilotr itself
                    # wrote. The reading is unambiguous and both engines already treat them alike.
                    if isinstance(vw, str):
                        vw = [vw]
                    if not isinstance(vw, list) or not vw:
                        bad("%s.vary_within must be a unit name or an array of unit names" % where)
                    else:
                        for w in vw:
                            if w not in ("subject", "item"):
                                bad("%s.vary_within contains '%s'; only 'subject' and 'item' are "
                                    "allowed" % (where, w))
                            elif w == "item" and not has_item:
                                bad("%s.vary_within names 'item' but the design has no item unit"
                                    % where)
                bt = f.get("between")
                if bt is not None:
                    if bt not in ("subject", "item"):
                        bad("%s.between must be 'subject' or 'item'" % where)
                    elif bt == "item" and not has_item:
                        bad("%s.between is 'item' but the design has no item unit" % where)
                if vw is None and bt is None:
                    bad("%s must set either 'vary_within' or 'between'" % where)

    # ---- predictors ----
    pred_names: list[str] = []
    predictors = spec.get("predictors")
    if predictors is not None:
        if not isinstance(predictors, list):
            bad("'predictors' must be an array of predictor objects")
        else:
            for i, p in enumerate(predictors, start=1):
                where = "predictors[%d]" % i
                if not isinstance(p, dict):
                    bad("%s must be an object" % where)
                    continue
                for k in p:
                    if k not in ("name", "varies_by", "mean", "sd", "dist", "min", "max",
                                 "reliability"):
                        unknown("unknown field '%s.%s'" % (where, k))
                if not _is_str(p.get("name")):
                    bad("%s.name must be a single string" % where)
                else:
                    pred_names.append(p["name"])
                vb = p.get("varies_by")
                if vb not in ("subject", "item", "observation"):
                    bad("%s.varies_by must be 'subject', 'item' or 'observation'%s"
                        % (where, (", not '%s'" % vb) if vb is not None else ""))
                elif vb == "item" and not has_item:
                    bad("%s.varies_by is 'item' but the design has no item unit" % where)
                dist = p.get("dist", "normal")
                if dist not in ("normal", "uniform"):
                    bad("%s.dist must be 'normal' or 'uniform'" % where)
                elif dist == "uniform":
                    if not _is_num(p.get("min")) or not _is_num(p.get("max")):
                        bad("%s uses dist 'uniform' and so needs numeric 'min' and 'max'" % where)
                    elif p["min"] >= p["max"]:
                        bad("%s.min must be less than %s.max" % (where, where))
                    for k in ("mean", "sd"):
                        if k in p:
                            unknown("%s.%s is ignored when dist is 'uniform'" % (where, k))
                else:
                    for k in ("min", "max"):
                        if k in p:
                            unknown("%s.%s is ignored when dist is 'normal'" % (where, k))
                    if p.get("mean") is not None and not _is_num(p["mean"]):
                        bad("%s.mean must be a number" % where)
                    if p.get("sd") is not None and (not _is_num(p["sd"]) or p["sd"] < 0):
                        bad("%s.sd must be a number of at least 0" % where)
                rel = p.get("reliability")
                if rel is not None and (not _is_num(rel) or rel <= 0 or rel > 1):
                    bad("%s.reliability must be greater than 0 and at most 1" % where)
    dupes = [n for n in dict.fromkeys(pred_names) if pred_names.count(n) > 1]
    if dupes:
        bad("duplicated predictor name(s): %s" % ", ".join(dupes))
    known_cols = contrast_cols + pred_names

    def check_key(key, where):
        """Every coefficient and slope key must resolve to a contrast column or a predictor. An
        unresolved key contributes zero, so a typo silently removes the effect."""
        miss = [pp for pp in key.split(":") if pp not in known_cols]
        if miss:
            bad("%s '%s' names %s, which %s neither a contrast column nor a predictor; available "
                "columns are %s. An unresolved key contributes zero, so this would silently drop "
                "the term"
                % (where, key, ", ".join("'%s'" % m for m in miss),
                   "are" if len(miss) > 1 else "is",
                   ", ".join("'%s'" % c for c in known_cols) if known_cols else "(none)"))

    # ---- fixed ----
    fx = spec.get("fixed")
    if fx is not None:
        if not isinstance(fx, dict):
            bad("'fixed' must be an object")
        else:
            for k in fx:
                if k not in ("intercept", "coefficients"):
                    unknown("unknown field 'fixed.%s'" % k)
            if not _is_num(fx.get("intercept")):
                bad("'fixed.intercept' must be a single number")
            coeffs = fx.get("coefficients")
            if coeffs is None:
                bad("'fixed.coefficients' is required (use {} for none)")
            elif not isinstance(coeffs, dict):
                bad("'fixed.coefficients' must be an object")
            else:
                for k, v in coeffs.items():
                    if not _is_num(v):
                        bad("'fixed.coefficients.%s' must be a single number" % k)
                    check_key(k, "fixed.coefficients")

    # ---- random ----
    random_spec = spec.get("random")
    if random_spec:
        if not isinstance(random_spec, dict):
            bad("'random' must be an object keyed by grouping factor")
        else:
            for g, re in random_spec.items():
                where = "random.%s" % g
                if not isinstance(re, dict):
                    bad("%s must be an object" % where)
                    continue
                for k in re:
                    if k not in ("intercept_sd", "slopes", "correlations", "correlated",
                                 "over", "n"):
                        unknown("unknown field '%s.%s'" % (where, k))
                if not _is_num(re.get("intercept_sd")) or re["intercept_sd"] < 0:
                    bad("%s.intercept_sd is required and must be at least 0" % where)
                slopes = re.get("slopes")
                cols = ["intercept"] + list(slopes.keys() if isinstance(slopes, dict) else [])
                if slopes is not None:
                    if not isinstance(slopes, dict):
                        bad("%s.slopes must be an object" % where)
                    else:
                        for k, v in slopes.items():
                            if not _is_num(v) or v < 0:
                                bad("%s.slopes.%s must be a number of at least 0" % (where, k))
                            check_key(k, "%s.slopes" % where)
                cors = re.get("correlations")
                if cors is not None:
                    if not isinstance(cors, dict):
                        bad("%s.correlations must be an object" % where)
                    else:
                        for k, v in cors.items():
                            if not _is_num(v) or v < -1 or v > 1:
                                bad("%s.correlations.%s must be between -1 and 1" % (where, k))
                            parts = [s.strip() for s in k.replace("~", ",").split(",")]
                            if len(parts) != 2:
                                bad("%s.correlations key '%s' must name two terms, as 'a,b'"
                                    % (where, k))
                            else:
                                miss = [pp for pp in parts if pp not in cols]
                                if miss:
                                    bad("%s.correlations key '%s' names %s, which is not a "
                                        "random-effect term of %s; its terms are %s"
                                        % (where, k, ", ".join("'%s'" % m for m in miss), g,
                                           ", ".join("'%s'" % c for c in cols)))
                if re.get("correlated") is not None and not isinstance(re["correlated"], bool):
                    bad("%s.correlated must be true or false" % where)
                if re.get("correlated") is False and cors:
                    bad("%s sets correlated = false but also supplies correlations; one of the "
                        "two has to go" % where)
                if g in ("subject", "item"):
                    for k in ("over", "n"):
                        if k in re:
                            bad("%s.%s applies only to an extra grouping factor, not to '%s'"
                                % (where, k, g))
                else:
                    if re.get("over") not in ("subject", "item"):
                        bad("%s.over is required for an extra grouping factor and must be "
                            "'subject' or 'item'" % where)
                    elif re["over"] == "item" and not has_item:
                        bad("%s.over is 'item' but the design has no item unit" % where)
                    if not _is_whole(re.get("n")) or re["n"] < 1:
                        bad("%s.n is required for an extra grouping factor and must be a whole "
                            "number of at least 1" % where)

    # ---- response ----
    r = spec.get("response")
    if r is not None:
        if not isinstance(r, dict):
            bad("'response' must be an object")
        else:
            fam = r.get("family")
            if fam not in FAMILY_PARAMS:
                bad("'response.family' must be one of %s%s"
                    % (", ".join(FAMILY_PARAMS), (", not '%s'" % fam) if _is_str(fam) else ""))
            if not _is_str(r.get("name")) or not r["name"]:
                bad("'response.name' must be a non-empty string")
            if fam in FAMILY_PARAMS:
                needed = FAMILY_PARAMS[fam]
                allowed = ("family", "name", "round") + needed
                for k in r:
                    if k not in allowed:
                        unknown("'response.%s' is not used by the %s family; it would be silently "
                                "ignored" % (k, fam))
                for k in needed:
                    if r.get(k) is None:
                        bad("'response.%s' is required for the %s family" % (k, fam))
                for k in ("sigma", "beta", "phi"):
                    if r.get(k) is not None and (not _is_num(r[k]) or r[k] <= 0):
                        bad("'response.%s' must be greater than 0" % k)
                if r.get("shift") is not None and not _is_num(r["shift"]):
                    bad("'response.shift' must be a number")
                th = r.get("thresholds")
                if th is not None:
                    # As with vary_within, a single cut-point may arrive as a bare number, because
                    # pilotr's own spec_json() collapsed one-element arrays before 0.3.
                    if _is_num(th):
                        th = [th]
                    if not isinstance(th, list) or not th or not all(_is_num(x) for x in th):
                        bad("'response.thresholds' must be a number or a non-empty numeric array")
                    elif any(th[i + 1] <= th[i] for i in range(len(th) - 1)):
                        bad("'response.thresholds' must be strictly increasing")
                if r.get("round") is not None:
                    if not _is_whole(r["round"]) or r["round"] < 0:
                        bad("'response.round' must be a whole number of at least 0")
                    elif fam not in ROUNDING_FAMILIES:
                        unknown("'response.round' has no effect for the %s family, whose outcome "
                                "is already an integer" % fam)

    if soft:
        warnings.warn("in this design specification:\n  - " + "\n  - ".join(soft), stacklevel=2)
    if problems:
        raise ValueError("invalid design specification:\n  - " + "\n  - ".join(problems))
    return spec
