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

import copy
import math

# The specification version this implementation writes and understands.
SPEC_VERSION = "0.3"

# The largest seed both twins read exactly. R's jsonlite reads a JSON integer as a double, which
# holds every integer up to 2**53 - 1 and no longer every one beyond it, while Python's json reads
# an exact integer of any size. A larger seed was therefore two different numbers in the two
# readers, which then seeded different streams.
MAX_EXACT_SEED = 9007199254740991

_SEED_RANGE = ("'seed' must be a whole number between -9007199254740991 and 9007199254740991 "
               "(2^53 - 1), the range in which R and Python read a JSON integer exactly")


def _repeated_key_message(key) -> str:
    return ("the specification repeats the key '%s' within one object; JSON leaves a repeated "
            "key undefined, and R and Python read it differently" % key)

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


def _is_name(x) -> bool:
    """A name that reaches the data as a column name.

    A blank one passed ``_is_str``, and the spec then produced a column with no name at all,
    where the R twin failed inside its simulator with a message naming neither the field nor the
    emptied control.
    """
    return _is_str(x) and x.strip() != ""


def _name_list(xs) -> str:
    """Join names as 'a', 'a and b', 'a, b and c'."""
    xs = list(xs)
    return "".join(xs) if len(xs) < 2 else "%s and %s" % (", ".join(xs[:-1]), xs[-1])


def _is_num(x) -> bool:
    # bool is a subclass of int in Python, and a boolean where a number belongs is a mistake.
    if isinstance(x, bool) or not isinstance(x, (int, float)):
        return False
    # An integer too large for a double is infinite to R, which reads every JSON number as one,
    # and math.isfinite() raises OverflowError on it where a refusal belongs.
    try:
        return math.isfinite(x)
    except OverflowError:
        return False


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
    0.2 implementation, so it has to declare 0.3 or later.

    Only well-formed parts are inspected. This scan runs before the per-field checks, so a
    string-valued `response`, a list-valued `random` or a numeric `slopes` reached it first and
    raised AttributeError or TypeError; skipping them here leaves each field's own check to raise
    the ValueError the documentation promises.
    """
    found = []
    predictors = spec.get("predictors")
    for p in predictors if isinstance(predictors, list) else []:
        if not isinstance(p, dict):
            continue
        if p.get("varies_by") == "observation":
            found.append('predictors varies_by "observation"')
        if p.get("dist") is not None:
            found.append("predictors dist")
        rel = p.get("reliability")
        if rel is not None and not (_is_num(rel) and rel == 1):
            found.append("predictors reliability")
    response = spec.get("response")
    if isinstance(response, dict) and response.get("family") == "exgaussian":
        found.append('the "exgaussian" family')
    random_spec = spec.get("random")
    for re in random_spec.values() if isinstance(random_spec, dict) else []:
        if not isinstance(re, dict):
            continue
        if re.get("correlated") is not None:
            found.append("random correlated")
        slopes = re.get("slopes")
        if isinstance(slopes, dict) and any(":" in k for k in slopes):
            found.append("interaction random slopes")
    return list(dict.fromkeys(found))


def validate_spec(spec, strict: bool = True):
    """Validate a design specification, returning it so the call can be chained.

    Checks the specification against the portable schema and against the cross-field rules the
    schema cannot express, and checks that its declared ``spec_version`` is one this
    implementation understands. Called by ``load_spec`` by default.

    Names are checked together as well as one by one. Two columns with one name, a contrast
    column defined by two factors or named like another column, an interaction whose analysis
    column already exists, a level listed twice, a correlation that pairs a term with itself or
    gives one pair twice, a factor both between and within a unit, and a ``random.item`` entry in
    a design without items each validated and then changed the data in silence. The rules are set
    out under Names in the specification.

    At most one factor may be between each unit. pilotr assigns the levels of each between factor
    to blocks of its unit on its own, so two between factors over one unit fall into the same or
    overlapping blocks. A 2 x 2 between-subjects design over 40 subjects gave cells of 20, 0, 0
    and 20, and its effects could not be estimated apart. Such a design is written as one between
    factor whose levels are the cells, as the specification shows under Worked encodings.

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

    Warns
    -----
    UserWarning
        For an unrecognised field when ``strict`` is false, and in either mode for a within
        factor whose ``vary_within`` leaves out a unit of the design, which pilotr crosses with
        every unit anyway. From spec_version 0.4 the latter is an error.
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

    # A JSON null in an optional field means the field is absent, as it does to R's jsonlite, so
    # every optional field below is tested with `is None`, never by its presence or truthiness.

    # ---- version ----
    declared = spec.get("spec_version")
    if declared is None:
        declared = "0.2"
    # A version is one string or one number. R stopped with a base error on an empty array or
    # object and read a two-element array by its first element, so both twins refuse the lot.
    is_scalar = _is_str(declared) or _is_num(declared)
    dv = _parse_version(declared) if is_scalar else None
    sv = _parse_version(SPEC_VERSION)
    if not is_scalar:
        bad("'spec_version' must be a single string of the form 'major.minor'")
    elif dv is None:
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
    seed = spec.get("seed")
    if seed is not None and not (_is_whole(seed) and abs(seed) <= MAX_EXACT_SEED):
        bad(_SEED_RANGE)

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
            # In the order the file lists them, as R reports them, and a unit present as null
            # is refused as R refuses it. Skipping a null `item` let it through to simulate(),
            # which then failed on subscripting None.
            for nm in [k for k in units if k in ("subject", "item")]:
                un = units[nm]
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
                if not _is_name(f.get("name")):
                    bad("%s.name must be a non-empty string" % where)
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
                    if not isinstance(vw, list) or not vw or not all(_is_str(w) for w in vw):
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
                if not _is_name(p.get("name")):
                    bad("%s.name must be a non-empty string" % where)
                else:
                    pred_names.append(p["name"])
                vb = p.get("varies_by")
                if vb not in ("subject", "item", "observation"):
                    # Only a string is quoted back, as in R, which renders other values its own way.
                    bad("%s.varies_by must be 'subject', 'item' or 'observation'%s"
                        % (where, (", not '%s'" % vb) if _is_str(vb) else ""))
                elif vb == "item" and not has_item:
                    bad("%s.varies_by is 'item' but the design has no item unit" % where)
                dist = p.get("dist")
                if dist is None:
                    dist = "normal"
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
    # Only an empty object or array means "none", as in R, where `"random": 0` or `false` is a
    # malformed value and refused.
    if random_spec is not None and not (isinstance(random_spec, (dict, list)) and not random_spec):
        if not isinstance(random_spec, dict):
            bad("'random' must be an object keyed by grouping factor")
        else:
            for g, re in random_spec.items():
                if not _is_name(g):
                    bad("a 'random' grouping factor must have a non-empty name")
                    continue
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
            # The type test comes first: a list or an object is unhashable, and looking one up
            # in FAMILY_PARAMS raised TypeError where a ValueError is promised.
            known = _is_str(fam) and fam in FAMILY_PARAMS
            if not known:
                bad("'response.family' must be one of %s%s"
                    % (", ".join(FAMILY_PARAMS), (", not '%s'" % fam) if _is_str(fam) else ""))
            if not _is_name(r.get("name")):
                bad("'response.name' must be a non-empty string")
            if known:
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

    # ---- names and structure ----
    # An empty object is no object to R, which reads `{}` as an unnamed empty list.
    found, found_soft = _check_names(spec, has_item,
                                     units_ok=isinstance(units, dict) and len(units) > 0)
    problems.extend(found)
    soft.extend(found_soft)

    if soft:
        warnings.warn("in this design specification:\n  - " + "\n  - ".join(soft), stacklevel=2)
    if problems:
        raise ValueError("invalid design specification:\n  - " + "\n  - ".join(problems))
    return spec


def _check_names(spec, has_item, units_ok):
    """The rules that look at a specification's names together.

    They run after the per-field checks and mirror the R twin's ``.check_names()`` line for line.
    Each case they refuse used to validate and then move an effect, rescale a variance or
    overwrite a column without a word. Where a column was overwritten the twins also disagreed,
    R replacing the earlier column and Python appending a second one under the same name, so one
    specification gave two different tables. Only well-formed parts are inspected, since a
    malformed one is reported by its own check.

    Returns the refusals and the warnings, in the order both twins report them.
    """
    problems: list[str] = []
    soft: list[str] = []
    bad = problems.append
    factors = spec.get("factors") if isinstance(spec.get("factors"), list) else []
    predictors = spec.get("predictors") if isinstance(spec.get("predictors"), list) else []
    rs = spec.get("random") if isinstance(spec.get("random"), dict) else {}
    groups = [g for g in rs if _is_name(g)]
    r = spec.get("response")

    def repeats(xs):
        """Repeats in the order of their first appearance, as the R twin lists them."""
        return [x for x in dict.fromkeys(xs) if xs.count(x) > 1]

    # A level listed twice put the declared effect into the data, since the simulator works by
    # position, while R's model_data(), which matches labels, gave every row the first level's
    # value.
    for i, f in enumerate(factors, start=1):
        if not isinstance(f, dict):
            continue
        lv = f.get("levels")
        if isinstance(lv, list) and all(_is_str(v) for v in lv):
            for level in repeats(lv):
                bad("factors[%d].levels repeats '%s'; each level needs its own label" % (i, level))
        # With both fields each unit got one row per level, every one at the level that the
        # between assignment chose.
        if f.get("vary_within") is not None and f.get("between") is not None:
            bad("factors[%d] sets both 'vary_within' and 'between'; a factor has to set exactly "
                "one of them" % i)

    # The simulator assigns the levels of each between factor to blocks of its unit on its own,
    # so two such factors over one unit fell into the same or overlapping blocks. A 2 x 2
    # between-subjects design over 40 subjects gave cells of 20, 0, 0 and 20, and neither the
    # second effect nor the interaction could then be estimated. One factor whose levels are the
    # cells keeps every cell, and R's spec_from_model() writes that encoding for a fitted pilot.
    for unit in ["subject"] + (["item"] if has_item else []):
        on_unit = [f["name"] for f in factors
                   if isinstance(f, dict) and f.get("vary_within") is None
                   and _is_name(f.get("name")) and f.get("between") == unit]
        if len(on_unit) > 1:
            bad("the factors %s are %s between '%s'. pilotr assigns the levels of each between "
                "factor to blocks of %ss on its own, so the blocks of these factors coincide or "
                "overlap, which leaves some combinations of their levels without %ss and "
                "confounds their effects. Encode the design as one between factor whose levels "
                "are the cells, give it the contrast columns of these factors and key each "
                "interaction as 'a:b', as spec_from_model() in R does for a fitted pilot. "
                "Specification version 0.4 will allocate several between factors jointly."
                % (_name_list("'%s'" % n for n in on_unit),
                   "both" if len(on_unit) == 2 else "all", unit, unit, unit))

    # Every column of the simulated data, in the order the simulator writes them.
    cols: list[str] = []
    roles: list[str] = []

    def claim(nm, role):
        if _is_name(nm):
            cols.append(nm)
            roles.append(role)

    if units_ok:
        claim("subject", "the subject unit")
        if has_item:
            claim("item", "the item unit")
    for g in groups:
        if g not in ("subject", "item"):
            claim(g, "random.%s" % g)
    for i, f in enumerate(factors, start=1):
        if isinstance(f, dict):
            claim(f.get("name"), "factors[%d]" % i)
    for i, p in enumerate(predictors, start=1):
        if isinstance(p, dict):
            claim(p.get("name"), "predictors[%d]" % i)
    if isinstance(r, dict):
        claim(r.get("name"), "the response")

    def roles_of(nm):
        return [w for c, w in zip(cols, roles) if c == nm]

    for nm in repeats(cols):
        who = roles_of(nm)
        # Two predictors alone are already reported as a duplicated predictor name.
        if all(w.startswith("predictors[") for w in who):
            continue
        bad("the name '%s' is used for more than one column (%s); every unit, grouping factor, "
            "factor, predictor and the response needs a column of its own" % (nm, _name_list(who)))

    # Contrast columns are where the coefficients and slopes act. One defined by two factors
    # carried only the later factor's effect, and one named like a predictor took the
    # predictor's draws in place of the contrast. A contrast named after its own factor is
    # allowed: two-level specifications use it, and the factor's name is checked above.
    contrasts: list[tuple[str, int]] = []
    for i, f in enumerate(factors, start=1):
        cs = f.get("contrasts") if isinstance(f, dict) else None
        if isinstance(cs, dict):
            contrasts.extend((cn, i) for cn in cs)
    for cn in dict.fromkeys(c for c, _ in contrasts):
        fs = list(dict.fromkeys(i for c, i in contrasts if c == cn))
        if len(fs) > 1:
            bad("the contrast column '%s' is defined by more than one factor (%s); each contrast "
                "column has to belong to a single factor"
                % (cn, _name_list("factors[%d]" % i for i in fs)))
    for cn, i in contrasts:
        if cn == factors[i - 1].get("name"):
            continue
        who = roles_of(cn)
        if who:
            bad("the contrast column '%s' of factors[%d] is also the name of %s; a contrast "
                "column may share its own factor's name but no other column's"
                % (cn, i, _name_list(who)))

    # R's model_data() writes an interaction key "a:b" to a product column "a_b", over any column
    # that already has that name.
    fx = spec.get("fixed")
    keys: list[str] = []
    if isinstance(fx, dict) and isinstance(fx.get("coefficients"), dict):
        keys.extend(fx["coefficients"])
    for g in groups:
        re = rs[g]
        if isinstance(re, dict) and isinstance(re.get("slopes"), dict):
            keys.extend(re["slopes"])
    keys = [k for k in dict.fromkeys(keys) if ":" in k]
    for j, key in enumerate(keys):
        col = key.replace(":", "_")
        who = (roles_of(col)
               + ["a contrast column of factors[%d]" % i
                  for i in dict.fromkeys(i for c, i in contrasts if c == col)]
               + ["the analysis column of the interaction '%s'" % k
                  for k in keys[:j] if k.replace(":", "_") == col])
        if who:
            bad("the interaction '%s' becomes the analysis column '%s', which is already the name "
                "of %s" % (key, col, _name_list(who)))

    # A term paired with itself overwrote the unit diagonal of the correlation matrix, which
    # rescaled that term's variance, and a pair given twice kept whichever came later.
    for g in groups:
        re = rs[g]
        cors = re.get("correlations") if isinstance(re, dict) else None
        if not isinstance(cors, dict):
            continue
        seen: dict[tuple[str, str], str] = {}
        for k in cors:
            parts = [s.strip() for s in k.replace("~", ",").split(",")]
            if len(parts) != 2:
                continue
            a, b = parts
            if a == b:
                bad("random.%s.correlations key '%s' pairs '%s' with itself; a term's correlation "
                    "with itself is always 1" % (g, k, a))
                continue
            first = seen.get((a, b), seen.get((b, a)))
            if first is not None:
                bad("random.%s.correlations keys '%s' and '%s' name the same pair of terms; give "
                    "each pair once" % (g, first, k))
            else:
                seen[(a, b)] = k

    # The simulator drew no item effects without an item unit, while R's model_formula() still
    # emitted an item term for the analysis.
    if units_ok and not has_item and rs.get("item") is not None:
        bad("random.item describes an item unit the design does not have; add units.item or "
            "remove random.item")

    # pilotr crosses a within factor with every unit of the design, whatever vary_within lists,
    # so an incomplete list has never changed the data. It is a warning until spec_version 0.4.
    if units_ok:
        design_units = ["subject"] + (["item"] if has_item else [])
        for i, f in enumerate(factors, start=1):
            if not isinstance(f, dict) or f.get("between") is not None:
                continue
            vw = f.get("vary_within")
            if isinstance(vw, str):
                vw = [vw]
            if not isinstance(vw, list) or not vw or not all(w in design_units for w in vw):
                continue
            miss = [u for u in design_units if u not in vw]
            if miss:
                quoted = _name_list("'%s'" % u for u in miss)
                soft.append(
                    "factors[%d].vary_within lists %s but not %s. pilotr crosses a within factor "
                    "with every unit of the design, so it varies within %s as well; list every "
                    "unit, or make it between 'item' if items carry it. From spec_version 0.4 "
                    "this is an error."
                    % (i, _name_list("'%s'" % w for w in dict.fromkeys(vw)), quoted, quoted))

    return problems, soft


def _drop_nulls(x):
    """Remove every None-valued entry from the objects nested in `x`, in place."""
    if isinstance(x, dict):
        for k in [k for k, v in x.items() if v is None]:
            del x[k]
        for v in x.values():
            _drop_nulls(v)
    elif isinstance(x, list):
        for v in x:
            _drop_nulls(v)


def _as_int(d, key):
    """Turn a whole-number float at d[key] into the int the engine counts with."""
    v = d.get(key) if isinstance(d, dict) else None
    if isinstance(v, float) and v.is_integer():
        d[key] = int(v)


def _normalise(spec):
    """Return a copy of a validated specification in the one shape the engine reads.

    Validation accepts what the schema and R accept, and the engine was written for one form of
    each. Three readings bridge the two, applied once at every public entry point:

    * A JSON null in an optional field means the field is absent, as it does to R, so it is
      removed, and a missing `factors`, `predictors` or `random` becomes empty.
    * A count written as a whole-number float, such as ``"n": 8.0``, becomes an int. Draft-07
      counts 8.0 as an integer and R simulates it, while ``range()`` raised TypeError on it.
    * A single threshold written as a bare number, which pilotr's own spec_json() wrote before
      0.3, becomes a one-element list.
    """
    s = copy.deepcopy(spec)
    if not isinstance(s, dict):
        return s
    _drop_nulls(s)
    s.setdefault("factors", [])
    s.setdefault("predictors", [])
    s.setdefault("random", {})
    if isinstance(s["random"], list) and not s["random"]:
        s["random"] = {}
    _as_int(s, "seed")
    units = s.get("units")
    if isinstance(units, dict):
        for unit in units.values():
            _as_int(unit, "n")
            _as_int(unit, "per_subject")
    if isinstance(s["random"], dict):
        for entry in s["random"].values():
            _as_int(entry, "n")
    response = s.get("response")
    if isinstance(response, dict):
        _as_int(response, "round")
        th = response.get("thresholds")
        if _is_num(th):
            response["thresholds"] = [th]
    return s
