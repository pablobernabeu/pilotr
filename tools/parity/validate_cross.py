"""Cross-language agreement check for reading and validating a specification.

The two validators are the gate that decides whether a specification is usable, so a
specification accepted by one implementation and refused by the other is itself a parity bug: it
would mean a design that runs in R and fails in Python, or worse, one that runs in both but with
different meaning. This script runs four batteries through both twins.

The spec battery hands each twin's validate_spec() the same parsed specification. The file
battery writes raw JSON text to a file and has each twin read it with its own load_spec(), which
is where the readers used to differ: one-element arrays, `[]`, JSON null, whole-number floats,
repeated keys, seeds beyond 2^53 and a byte-order mark. Where both twins accept a file, each
simulates it and the two dumps, in the 17-digit format of run_r.R and run_py.py, must hash alike.
A refusal raised while simulating, such as the cap on a Poisson mean, is compared as text like
any other.

The power battery checks the two-group backend's refusals, which power_design() in R and power()
in Python raise before drawing any replicate. Each twin validates the specification and applies
its two-group checks, the first steps of both functions, so the battery needs neither scipy nor
any simulation. The checks refuse a design whose rows are correlated, naming the item unit,
within factor or grouping factor that makes them so, and a design with fewer than two subjects at
a level of the between factor.

The scale battery compares, for every response family, the warning power_mixed() raises when it
fits a linear model to a response whose coefficients are on a link scale. The R twin raises it for
ordinal and beta responses, and fits bernoulli and poisson ones with glmer() instead, while the
Python twin raises it for all four. The battery therefore compares the wording and the scale
each twin gives a family, not whether a given call raises it, and needs no model fit.

A case passes when both twins accept it, or both refuse it with the same message character for
character. In the spec battery the warnings the two validators raise, such as the deprecation of
an incomplete `vary_within`, must also match character for character. A Python exception other
than the one each battery expects (ValueError, and NotImplementedError from the two-group checks),
or an R error raised by base R rather than by pilotr (which always stops with call. = FALSE),
counts as a crash and fails the case whatever the other twin did.

The schema check holds spec/design.schema.json to strict validation. Every specification of the
spec battery that both twins accept with strict = TRUE must be valid under jsonschema's
Draft7Validator. The schema used to be looser than the validators in one direction and stricter
in the other. It accepted an unknown top-level field and a response parameter the family does
not use, and it refused the single-value forms of vary_within and thresholds that both
validators accept. The check needs the jsonschema package (pip install jsonschema).

Usage: python tools/parity/validate_cross.py
Exit status is 0 when the two agree on every case and the schema accepts every case they accept.
"""

from __future__ import annotations

import copy
import hashlib
import json
import math
import os
import subprocess
import sys
import tempfile
import warnings

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "python"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from pilotr.power import (  # noqa: E402
    _LINEAR_ON_LINK, _LINK_SCALE, _two_group_refusal, _two_group_size_refusal)
from pilotr.simulate import _as_spec, load_spec, simulate  # noqa: E402
from pilotr.validate import validate_spec  # noqa: E402
from run_py import _dump  # noqa: E402  the dump format the parity harness compares

try:
    import jsonschema
except ImportError:  # reported in main(), so that the message says what to install
    jsonschema = None

# The interpreter on PATH by default, so the script runs anywhere R is installed (CI
# included); PILOTR_RSCRIPT overrides it for a machine whose Rscript lives elsewhere.
RSCRIPT = os.environ.get("PILOTR_RSCRIPT", "Rscript")


def _read(*parts):
    with open(os.path.join(ROOT, *parts), encoding="utf-8") as f:
        return json.load(f)


def base_spec():
    return _read("spec", "examples", "crossed_mixed_rt.json")


def _mut(fn):
    s = base_spec()
    fn(s)
    return s


def _set_family(s, **kw):
    s["response"] = kw


# (label, spec) pairs. Each is fed to both validators; the verdicts and messages must match.
def cases():
    out = []
    for name in sorted(os.listdir(os.path.join(ROOT, "spec", "examples"))):
        if name.endswith(".json"):
            out.append(("shipped:" + name, _read("spec", "examples", name)))
    for name in sorted(os.listdir(os.path.join(ROOT, "tools", "parity", "cases"))):
        if name.endswith(".json"):
            out.append(("case:" + name, _read("tools", "parity", "cases", name)))

    out += [
        ("mistyped coefficient key",
         _mut(lambda s: s["fixed"].__setitem__("coefficients", {"cnod": 0.05}))),
        ("mistyped slope key",
         _mut(lambda s: s["random"]["subject"].__setitem__("slopes", {"conditon": 0.04}))),
        ("interaction coefficient with a bad component",
         _mut(lambda s: s["fixed"].__setitem__("coefficients", {"cond:nope": 0.05}))),
        ("foreign response parameter", _mut(lambda s: s["response"].__setitem__("phi", 8))),
        ("round on an integer family",
         _mut(lambda s: _set_family(s, family="poisson", name="count", round=2))),
        ("round 0 on beta",
         _mut(lambda s: _set_family(s, family="beta", name="p", phi=8, round=0))),
        ("round 1 on beta",
         _mut(lambda s: _set_family(s, family="beta", name="p", phi=8, round=1))),
        ("unknown top-level field", _mut(lambda s: s.__setitem__("notafield", 1))),
        ("unknown factor field",
         _mut(lambda s: s["factors"][0].__setitem__("levls", ["a", "b"]))),
        ("varies_by trial",
         _mut(lambda s: s.__setitem__("predictors", [{"name": "ISI", "varies_by": "trial"}]))),
        ("varies_by observation without declaring 0.3",
         _mut(lambda s: s.__setitem__("predictors",
                                      [{"name": "ISI", "varies_by": "observation"}]))),
        ("varies_by observation declaring 0.3",
         _mut(lambda s: (s.__setitem__("spec_version", "0.3"),
                         s.__setitem__("predictors", [{"name": "ISI", "varies_by": "observation"}]),
                         s["fixed"]["coefficients"].__setitem__("ISI", 0.1)))),
        ("contrast length mismatch",
         _mut(lambda s: s["factors"][0].__setitem__("contrasts", {"cond": [-0.5, 0, 0.5]}))),
        ("correlation naming a missing term",
         _mut(lambda s: s["random"]["subject"].__setitem__("correlations",
                                                           {"intercept,nope": 0.2}))),
        ("correlation out of range",
         _mut(lambda s: s["random"]["subject"].__setitem__("correlations",
                                                           {"intercept,cond": 1.4}))),
        ("correlation key with one term",
         _mut(lambda s: s["random"]["subject"].__setitem__("correlations", {"intercept": 0.2}))),
        ("per_subject over the item count",
         _mut(lambda s: s["units"]["item"].__setitem__("per_subject", 999))),
        ("per_subject on subject",
         _mut(lambda s: s["units"]["subject"].__setitem__("per_subject", 2))),
        ("extra group without over/n",
         _mut(lambda s: s["random"].__setitem__("site", {"intercept_sd": 0.5}))),
        ("extra group with a valid over/n",
         _mut(lambda s: s["random"].__setitem__("site", {"intercept_sd": 0.5, "over": "subject",
                                                        "n": 5}))),
        ("over/n on subject",
         _mut(lambda s: s["random"]["subject"].__setitem__("over", "subject"))),
        ("thresholds not increasing",
         _mut(lambda s: _set_family(s, family="ordinal", name="r", thresholds=[1, 0.5, 2]))),
        ("valid ordinal",
         _mut(lambda s: _set_family(s, family="ordinal", name="r", thresholds=[-1, 0, 1]))),
        ("a single threshold as a bare number",
         _mut(lambda s: _set_family(s, family="ordinal", name="r", thresholds=0.5))),
        ("valid bernoulli", _mut(lambda s: _set_family(s, family="bernoulli", name="acc"))),
        ("negative sigma", _mut(lambda s: s["response"].__setitem__("sigma", -1))),
        ("missing sigma", _mut(lambda s: s["response"].pop("sigma"))),
        ("missing response", _mut(lambda s: s.pop("response"))),
        ("missing units.subject", _mut(lambda s: s["units"].pop("subject"))),
        ("zero subjects", _mut(lambda s: s["units"]["subject"].__setitem__("n", 0))),
        ("non-integer seed", _mut(lambda s: s.__setitem__("seed", 1.5))),
        ("factor with neither vary_within nor between",
         _mut(lambda s: s["factors"][0].pop("vary_within"))),
        ("vary_within item with no item unit",
         _mut(lambda s: (s["units"].pop("item"), s["random"].pop("item")))),
        ("negative slope sd",
         _mut(lambda s: s["random"]["subject"]["slopes"].__setitem__("cond", -0.1))),
        ("correlated false with correlations",
         _mut(lambda s: s["random"]["subject"].__setitem__("correlated", False))),
        ("correlated flag without declaring 0.3",
         _mut(lambda s: (s["random"]["subject"].pop("correlations"),
                         s["random"]["subject"].__setitem__("correlated", False)))),
        ("exgaussian without declaring 0.3",
         _mut(lambda s: _set_family(s, family="exgaussian", name="RT", sigma=0.2, beta=0.3))),
        ("exgaussian declaring 0.3",
         _mut(lambda s: (s.__setitem__("spec_version", "0.3"),
                         _set_family(s, family="exgaussian", name="RT", sigma=0.2, beta=0.3)))),
        ("spec from the future", _mut(lambda s: s.__setitem__("spec_version", "9.9"))),
        ("malformed spec_version", _mut(lambda s: s.__setitem__("spec_version", "banana"))),
        ("unknown family", _mut(lambda s: s["response"].__setitem__("family", "weibull"))),
        ("empty response name", _mut(lambda s: s["response"].__setitem__("name", ""))),
        ("duplicate predictor names",
         _mut(lambda s: (s.__setitem__("predictors",
                                       [{"name": "z", "varies_by": "subject"},
                                        {"name": "z", "varies_by": "item"}]),
                         s["fixed"]["coefficients"].__setitem__("z", 0.1)))),
        ("uniform predictor without min/max",
         _mut(lambda s: (s.__setitem__("spec_version", "0.3"),
                         s.__setitem__("predictors",
                                       [{"name": "z", "varies_by": "subject", "dist": "uniform"}]),
                         s["fixed"]["coefficients"].__setitem__("z", 0.1)))),
        ("uniform predictor with min >= max",
         _mut(lambda s: (s.__setitem__("spec_version", "0.3"),
                         s.__setitem__("predictors",
                                       [{"name": "z", "varies_by": "subject", "dist": "uniform",
                                         "min": 1, "max": 0}]),
                         s["fixed"]["coefficients"].__setitem__("z", 0.1)))),
        ("valid uniform predictor",
         _mut(lambda s: (s.__setitem__("spec_version", "0.3"),
                         s.__setitem__("predictors",
                                       [{"name": "z", "varies_by": "subject", "dist": "uniform",
                                         "min": 0, "max": 1}]),
                         s["fixed"]["coefficients"].__setitem__("z", 0.1)))),
        ("reliability out of range",
         _mut(lambda s: (s.__setitem__("spec_version", "0.3"),
                         s.__setitem__("predictors",
                                       [{"name": "z", "varies_by": "subject", "reliability": 1.5}]),
                         s["fixed"]["coefficients"].__setitem__("z", 0.1)))),
        ("valid reliability",
         _mut(lambda s: (s.__setitem__("spec_version", "0.3"),
                         s.__setitem__("predictors",
                                       [{"name": "z", "varies_by": "subject", "reliability": 0.8}]),
                         s["fixed"]["coefficients"].__setitem__("z", 0.1)))),
    ]
    out += name_cases()
    out += allocation_cases()
    out += counterbalance_cases()
    return out


def _second_factor(s, **kw):
    f = {"name": "block", "levels": ["p", "q"], "contrasts": {"blk": [-0.5, 0.5]},
         "vary_within": ["subject", "item"]}
    f.update(kw)
    s["factors"].append(f)


def name_cases():
    """Names, columns and structure: each validated before and then changed the data in silence.

    The last few are accepted with a warning, which the spec battery compares as text too.
    """
    item_less = lambda s: (s["units"].pop("item"),  # noqa: E731
                           s["factors"][0].__setitem__("vary_within", ["subject"]))
    return [
        ("factor named like the subject column",
         _mut(lambda s: s["factors"][0].__setitem__("name", "subject"))),
        ("response named like the factor",
         _mut(lambda s: s["response"].__setitem__("name", s["factors"][0]["name"]))),
        ("grouping factor named like the response",
         _mut(lambda s: (s["random"].__setitem__("site", {"intercept_sd": 0.5, "over": "subject",
                                                          "n": 5}),
                         s["response"].__setitem__("name", "site")))),
        ("two factors sharing a contrast column",
         _mut(lambda s: _second_factor(s, contrasts={"cond": [-0.5, 0.5]}))),
        ("between factor reusing the within contrast",
         _mut(lambda s: (_second_factor(s, contrasts={"cond": [-0.5, 0.5]}, between="subject"),
                         s["factors"][1].pop("vary_within")))),
        ("predictor named like a contrast column",
         _mut(lambda s: s.__setitem__("predictors", [{"name": "cond", "varies_by": "subject"}]))),
        ("contrast column named item",
         _mut(lambda s: s["factors"][0]["contrasts"].__setitem__("item", [-0.5, 0.5]))),
        ("contrast column named like another factor",
         _mut(lambda s: _second_factor(s, contrasts={"condition": [-0.5, 0.5]}))),
        ("interaction column named like a predictor",
         _mut(lambda s: (s.__setitem__("predictors",
                                       [{"name": "freq", "varies_by": "item"},
                                        {"name": "cond_freq", "varies_by": "subject"}]),
                         s["fixed"]["coefficients"].__setitem__("cond:freq", 0.01)))),
        ("interaction column named like a contrast column",
         _mut(lambda s: (s["factors"][0]["contrasts"].update(cond_blk=[1, -1]),
                         _second_factor(s),
                         s["fixed"]["coefficients"].__setitem__("cond:blk", 0.01)))),
        ("a level listed twice",
         _mut(lambda s: s["factors"][0].__setitem__("levels", ["x", "x"]))),
        ("factor both between and within",
         _mut(lambda s: s["factors"][0].__setitem__("between", "subject"))),
        ("correlation pairing a term with itself",
         _mut(lambda s: s["random"]["subject"].__setitem__("correlations", {"cond,cond": 0.25}))),
        ("correlation pair given twice",
         _mut(lambda s: s["random"]["subject"].__setitem__(
             "correlations", {"intercept,cond": 0.9, "cond~intercept": -0.9}))),
        ("random.item without an item unit", _mut(item_less)),
        ("blank factor name", _mut(lambda s: s["factors"][0].__setitem__("name", "  "))),
        ("blank grouping factor name",
         _mut(lambda s: s["random"].__setitem__("", {"intercept_sd": 0.5, "over": "subject",
                                                    "n": 5}))),
        ("thresholds given as a string",
         _mut(lambda s: _set_family(s, family="ordinal", name="r", thresholds="low"))),
        ("correlated given as a string",
         _mut(lambda s: (s.__setitem__("spec_version", "0.3"),
                         s["random"]["subject"].pop("correlations"),
                         s["random"]["subject"].__setitem__("correlated", "yes")))),
        ("three levels, first contrast named after its factor",
         _mut(lambda s: (s["factors"][0].update(name="cond", levels=["a", "b", "c"],
                                                contrasts={"cond": [-1, 1, 0],
                                                           "cond2": [-1, 0, 1]}),
                         s["random"]["subject"]["slopes"].__setitem__("cond2", 0.01)))),
        ("vary_within omitting the item unit (warns)",
         _mut(lambda s: s["factors"][0].__setitem__("vary_within", ["subject"]))),
        ("vary_within omitting the subject unit (warns)",
         _mut(lambda s: s["factors"][0].__setitem__("vary_within", "item"))),
        ("vary_within omitting the item unit, declaring 0.3 (warns)",
         _mut(lambda s: (s.__setitem__("spec_version", "0.3"),
                         s["factors"][0].__setitem__("vary_within", ["subject"])))),
        ("vary_within omitting the item unit, declaring 0.4",
         _mut(lambda s: (s.__setitem__("spec_version", "0.4"),
                         s["factors"][0].__setitem__("vary_within", ["subject"])))),
        ("vary_within omitting the subject unit, declaring 0.4",
         _mut(lambda s: (s.__setitem__("spec_version", "0.4"),
                         s["factors"][0].__setitem__("vary_within", "item")))),
    ]


def _between_factor(name, unit, contrasts):
    """A between factor with one level for each value of its contrasts."""
    n_levels = len(next(iter(contrasts.values())))
    return {"name": name, "levels": ["%s%d" % (name.lower(), k) for k in range(1, n_levels + 1)],
            "contrasts": contrasts, "between": unit}


def allocation_cases():
    """Factors between one unit. Below specification version 0.4, pilotr assigned them to the
    same or overlapping blocks of it, and version 0.4 allocates them jointly.

    A specification declaring a version below 0.4 is refused for them. The rest are accepted: a
    factor between subjects and another between items do not alias, and version 0.4 allocates
    several factors on one unit jointly.
    """
    a = _between_factor("A", "subject", {"a": [-0.5, 0.5]})
    b = _between_factor("B", "subject", {"b": [-0.5, 0.5]})
    c = _between_factor("C", "subject", {"c1": [-1, 1, 0], "c2": [-1, 0, 1]})

    def declaring(version, factors):
        return _mut(lambda s: (s.__setitem__("spec_version", version),
                               s["factors"].extend(factors)))

    return [
        ("two factors between subject", _mut(lambda s: s["factors"].extend([a, b]))),
        ("three factors between subject", _mut(lambda s: s["factors"].extend([a, b, c]))),
        ("two factors between item",
         _mut(lambda s: s["factors"].extend([dict(a, between="item"), dict(b, between="item")]))),
        ("one factor between subject and one between item",
         _mut(lambda s: s["factors"].extend([a, dict(b, between="item")]))),
        ("two factors between subject, declaring 0.3", declaring("0.3", [a, b])),
        ("two factors between subject, declaring 0.4", declaring("0.4", [a, b])),
        ("three factors between subject, declaring 0.4", declaring("0.4", [a, b, c])),
        ("two factors between item, declaring 0.4",
         declaring("0.4", [dict(a, between="item"), dict(b, between="item")])),
    ]


_PRIME = {"name": "prime", "levels": ["p1", "p2"], "contrasts": {"pr": [-0.5, 0.5]},
          "vary_within": ["subject", "item"], "counterbalance": True}


def counterbalance_cases():
    """Counterbalanced lists, new in 0.4, which show each subject each item at one level.

    The field is refused in a specification declaring less than 0.4, whatever its value, and
    with a between factor, without items or with per_subject. Subjects or items that are not a
    multiple of the number of lists draw a warning, which the spec battery compares as text. The
    data of such a design name each subject's list in a column `list`, so another column of that
    name is refused, and accepted beside counterbalance false.
    """
    def counterbalanced(version="0.4", value=True, fn=None):
        def mutate(s):
            if version is not None:
                s["spec_version"] = version
            s["factors"][0]["counterbalance"] = value
            if fn is not None:
                fn(s)
        return _mut(mutate)

    return [
        ("counterbalanced condition, declaring 0.4", counterbalanced()),
        ("counterbalance false, declaring 0.4", counterbalanced(value=False)),
        ("counterbalanced condition, declaring 0.3", counterbalanced("0.3")),
        ("counterbalance false, declaring 0.3", counterbalanced("0.3", False)),
        ("counterbalanced condition, no spec_version", counterbalanced(None)),
        ("counterbalance given as a string", counterbalanced(value="yes")),
        ("counterbalance on a factor between items",
         counterbalanced(fn=lambda s: (s["factors"][0].pop("vary_within"),
                                       s["factors"][0].__setitem__("between", "item")))),
        ("counterbalanced without an item unit",
         counterbalanced(fn=lambda s: (s["units"].pop("item"), s["random"].pop("item"),
                                       s["factors"][0].__setitem__("vary_within", ["subject"])))),
        ("counterbalanced with per_subject",
         counterbalanced(fn=lambda s: s["units"]["item"].__setitem__("per_subject", 12))),
        ("counterbalanced over 31 subjects (warns)",
         counterbalanced(fn=lambda s: s["units"]["subject"].__setitem__("n", 31))),
        ("counterbalanced over 25 items (warns)",
         counterbalanced(fn=lambda s: s["units"]["item"].__setitem__("n", 25))),
        ("two counterbalanced factors over 32 subjects",
         counterbalanced(fn=lambda s: (s["factors"].append(dict(_PRIME)),
                                       s["units"]["subject"].__setitem__("n", 32),
                                       s["fixed"]["coefficients"].__setitem__("pr", 0.02)))),
        ("two counterbalanced factors over 30 subjects (warns)",
         counterbalanced(fn=lambda s: (s["factors"].append(dict(_PRIME)),
                                       s["fixed"]["coefficients"].__setitem__("pr", 0.02)))),
        ("counterbalanced beside a crossed within factor",
         counterbalanced(fn=_second_factor)),
        ("a factor named list beside a counterbalanced one",
         counterbalanced(fn=lambda s: _second_factor(s, name="list"))),
        ("a factor named list beside counterbalance false",
         counterbalanced(value=False, fn=lambda s: _second_factor(s, name="list"))),
        ("a predictor named list in a counterbalanced design",
         counterbalanced(fn=lambda s: s.__setitem__("predictors",
                                                    [{"name": "list", "varies_by": "subject"}]))),
        ("a grouping factor named list in a counterbalanced design",
         counterbalanced(fn=lambda s: s["random"].__setitem__(
             "list", {"intercept_sd": 0.5, "over": "subject", "n": 5}))),
        ("the response named list in a counterbalanced design",
         counterbalanced(fn=lambda s: s["response"].__setitem__("name", "list"))),
        ("a contrast column named list in a counterbalanced design",
         counterbalanced(fn=lambda s: _second_factor(s, contrasts={"list": [-0.5, 0.5]}))),
    ]


# ---- the file battery ---------------------------------------------------------------------

def _gaussian():
    """A crossed Gaussian design with slopes and correlations, so that a dump compares exactly."""
    return _read("tools", "parity", "cases", "control_rt_as_gaussian.json")


def _between():
    return _read("spec", "examples", "between_2group_gaussian.json")


def _poisson():
    """The shipped Poisson example with 200 subjects and no group effect, so that every row has
    the mean its intercept implies."""
    s = _read("spec", "examples", "poisson_counts_between.json")
    s["units"]["subject"]["n"] = 200
    s["fixed"]["coefficients"]["grp"] = 0.0
    return s


def _beta():
    """The shipped Beta example with no group effect, so that every row has the mean its
    intercept implies. It keeps the example's five decimals."""
    s = _read("spec", "examples", "beta_proportion.json")
    s["fixed"]["coefficients"]["grp"] = 0.0
    return s


def _lognormal(family, **extra):
    """The shipped two-group design with a lognormal response and no group effect."""
    s = _between()
    s["fixed"]["coefficients"]["grp"] = 0.0
    s["response"] = dict({"family": family, "name": "RT", "sigma": 0.3}, **extra)
    return s


def _text(spec):
    return json.dumps(spec, ensure_ascii=False, indent=1).encode("utf-8")


def _with(base, fn):
    s = base()
    fn(s)
    return _text(s)


def _counterbalance(value):
    """Set counterbalance on the first factor, in a specification declaring 0.4."""
    def fn(s):
        s["spec_version"] = "0.4"
        s["factors"][0]["counterbalance"] = value
    return fn


def _fixture(name):
    with open(os.path.join(ROOT, "python", "tests", "fixtures", name), "rb") as f:
        return f.read()


def _repeat_key(text, anchor, insert):
    """Insert `insert` right after the first occurrence of `anchor`, repeating a key."""
    s = text.decode("utf-8")
    if anchor not in s:
        raise SystemExit("file battery: anchor %r not found" % anchor)
    return s.replace(anchor, anchor + insert, 1).encode("utf-8")


def file_cases():
    """(label, raw bytes) pairs, each read by both twins' load_spec()."""
    g = _gaussian
    out = [
        ("write_json defaults (R fixture)", _fixture("write_json_defaults.json")),
        ("auto_unbox single threshold (R fixture)",
         _fixture("write_json_auto_unbox_ordinal.json")),
        ("seed [99]", _with(g, lambda s: s.__setitem__("seed", [99]))),
        ("spec_version [\"0.3\"]", _with(g, lambda s: s.__setitem__("spec_version", ["0.3"]))),
        ("between [\"subject\"]",
         _with(_between, lambda s: s["factors"][0].__setitem__("between", ["subject"]))),
        ("family and sigma as one-element arrays",
         _with(g, lambda s: s["response"].update(family=["gaussian"], sigma=[0.3]))),
        ("coefficients []",
         _with(_between, lambda s: s["fixed"].__setitem__("coefficients", []))),
        ("slopes [] and correlations []",
         _with(g, lambda s: s["random"]["item"].update(slopes=[], correlations=[]))),
        ("n 8.0", _with(_between, lambda s: s["units"]["subject"].__setitem__("n", 8.0))),
        ("per_subject 2.0",
         _with(lambda: _read("spec", "examples", "partial_crossing.json"),
               lambda s: s["units"]["item"].__setitem__("per_subject", 2.0))),
        ("random.site.n 2.0",
         _with(lambda: _read("spec", "examples", "nested_clusters.json"),
               lambda s: s["random"]["site"].__setitem__("n", 2.0))),
        ("round 2.0", _with(_between, lambda s: s["response"].__setitem__("round", 2.0))),
        ("factors null", _with(_between, lambda s: (s.__setitem__("factors", None),
                                                    s["fixed"].__setitem__("coefficients", {})))),
        ("predictors null", _with(g, lambda s: s.__setitem__("predictors", None))),
        ("slopes null and correlations null",
         _with(g, lambda s: s["random"]["item"].update(slopes=None, correlations=None))),
        ("predictor mean, sd, reliability and dist null",
         _with(g, lambda s: (s.__setitem__("predictors", [
             {"name": "z", "varies_by": "subject", "mean": None, "sd": None,
              "reliability": None, "dist": None}]),
             s["fixed"]["coefficients"].__setitem__("z", 0.1)))),
        ("spec_version null", _with(g, lambda s: s.__setitem__("spec_version", None))),
        ("units.item null",
         _with(_between, lambda s: s["units"].__setitem__("item", None))),
        ("repeated key in coefficients",
         _repeat_key(_text(_gaussian()), '"coefficients": {', '"cond": 0.5, ')),
        ("repeated key at the top level", _repeat_key(_text(_gaussian()), "{", '"seed": 1, ')),
        ("repeated key in a contrast",
         _repeat_key(_text(_gaussian()), '"contrasts": {', '"cond": [0.5, -0.5], ')),
        ("seed 2^53 - 1", _with(g, lambda s: s.__setitem__("seed", 2**53 - 1))),
        ("seed -(2^53 - 1)", _with(g, lambda s: s.__setitem__("seed", -(2**53 - 1)))),
        ("seed 2^53", _with(g, lambda s: s.__setitem__("seed", 2**53))),
        ("seed 9007199254740993", _with(g, lambda s: s.__setitem__("seed", 9007199254740993))),
        ("seed 123456789012345678", _with(g, lambda s: s.__setitem__("seed", 123456789012345678))),
        ("seed 10^400", _with(g, lambda s: s.__setitem__("seed", 10**400))),
        ("byte-order mark", b"\xef\xbb\xbf" + _text(_gaussian())),
        ("non-ASCII levels",
         _with(_between, lambda s: s["factors"][0].__setitem__("levels", ["fácil", "Łatwy"]))),
        ("levels [\"a\", null]",
         _with(_between, lambda s: s["factors"][0].__setitem__("levels", ["a", None]))),
        ("vary_within [\"subject\", null]",
         _with(g, lambda s: s["factors"][0].__setitem__("vary_within", ["subject", None]))),
        ("response \"gaussian\"", _with(g, lambda s: s.__setitem__("response", "gaussian"))),
        ("random [\"subject\"]", _with(g, lambda s: s.__setitem__("random", ["subject"]))),
        ("slopes 5", _with(g, lambda s: s["random"]["subject"].__setitem__("slopes", 5))),
        ("family [\"gaussian\", \"poisson\"]",
         _with(g, lambda s: s["response"].__setitem__("family", ["gaussian", "poisson"]))),
        ("family {\"name\": \"gaussian\"}",
         _with(g, lambda s: s["response"].__setitem__("family", {"name": "gaussian"}))),
        ("predictors 5", _with(g, lambda s: s.__setitem__("predictors", 5))),
        ("spec_version []", _with(g, lambda s: s.__setitem__("spec_version", []))),
        ("spec_version {}", _with(g, lambda s: s.__setitem__("spec_version", {}))),
        ("spec_version [\"0.3\", \"0.4\"]",
         _with(g, lambda s: s.__setitem__("spec_version", ["0.3", "0.4"]))),
        ("spec_version [null]", _with(g, lambda s: s.__setitem__("spec_version", [None]))),
        ("varies_by [\"subject\", \"item\"]",
         _with(g, lambda s: (s.__setitem__("predictors",
                                           [{"name": "z", "varies_by": ["subject", "item"]}]),
                             s["fixed"]["coefficients"].__setitem__("z", 0.1)))),
        ("random 0", _with(g, lambda s: s.__setitem__("random", 0))),
        ("random []", _with(g, lambda s: s.__setitem__("random", []))),
        # Two factors between subjects, allocated jointly from version 0.4.
        ("two between factors, declaring 0.4",
         _with(_between, lambda s: (s.__setitem__("spec_version", "0.4"),
                                    s["factors"].append(dict(_DOSE)),
                                    s["fixed"]["coefficients"].__setitem__("dose", 2.0)))),
        # Counterbalanced lists, new in 0.4. jsonlite reads [true] as true and [null] as NA, and a
        # null means the field is absent.
        ("counterbalance [true]", _with(g, _counterbalance([True]))),
        ("counterbalance null", _with(g, _counterbalance(None))),
        ("counterbalance [null]", _with(g, _counterbalance([None]))),
        ("counterbalanced over 31 subjects",
         _with(g, lambda s: (_counterbalance(True)(s),
                             s["units"]["subject"].__setitem__("n", 31)))),
        # Its data name each subject's list, in a column placed before the grouping factors.
        ("counterbalanced with subjects in sites",
         _with(g, lambda s: (_counterbalance(True)(s),
                             s["random"].__setitem__("site", {"intercept_sd": 0.5,
                                                              "over": "subject", "n": 4})))),
        # The Poisson cap is a refusal raised while simulating, below and above 708.4.
        ("poisson mean 700",
         _with(_poisson, lambda s: s["fixed"].__setitem__("intercept", math.log(700)))),
        ("poisson mean 742",
         _with(_poisson, lambda s: s["fixed"].__setitem__("intercept", math.log(742)))),
        ("poisson intercept 800 (exp overflows)",
         _with(_poisson, lambda s: s["fixed"].__setitem__("intercept", 800))),
        # The lognormal and Beta refusals, raised while simulating. At an intercept of 36.7, every
        # draw is exactly 1 and every row lands on 0.99999, so the data compare exactly too.
        ("lognormal intercept 750 (exp overflows)",
         _with(lambda: _lognormal("lognormal"),
               lambda s: s["fixed"].__setitem__("intercept", 750))),
        ("shifted_lognormal intercept 750 (exp overflows)",
         _with(lambda: _lognormal("shifted_lognormal", shift=200),
               lambda s: s["fixed"].__setitem__("intercept", 750))),
        ("beta intercept 36.7 (draws of 1 moved inside)",
         _with(_beta, lambda s: s["fixed"].__setitem__("intercept", 36.7))),
        ("beta intercept 40 (a shape of zero)",
         _with(_beta, lambda s: s["fixed"].__setitem__("intercept", 40))),
        ("beta intercept -750 (a shape of zero)",
         _with(_beta, lambda s: s["fixed"].__setitem__("intercept", -750))),
        ("beta intercept -700 and phi 1e-30 (mu * phi underflows)",
         _with(_beta, lambda s: (s["fixed"].__setitem__("intercept", -700),
                                 s["response"].__setitem__("phi", 1e-30)))),
        ("beta phi 0.001 (both Gamma variates underflow)",
         _with(_beta, lambda s: s["response"].__setitem__("phi", 0.001))),
    ]
    return out


# ---- the power battery --------------------------------------------------------------------

def _two_group(fn):
    """The shipped two-group design, one row per subject, changed by `fn`."""
    s = _between()
    fn(s)
    return s


_BLOCK = {"name": "block", "levels": ["x", "y"], "contrasts": {"blk": [-0.5, 0.5]},
          "vary_within": ["subject"]}
_DOSE = {"name": "dose", "levels": ["low", "high"], "contrasts": {"dose": [-0.5, 0.5]},
         "between": "subject"}


def power_cases():
    """(label, spec) pairs, each put through both twins' two-group check."""
    out = []
    for name in sorted(os.listdir(os.path.join(ROOT, "spec", "examples"))):
        if name.endswith(".json"):
            out.append(("shipped:" + name, _read("spec", "examples", name)))
    for name in sorted(os.listdir(os.path.join(ROOT, "tools", "parity", "cases"))):
        if name.endswith(".json"):
            out.append(("case:" + name, _read("tools", "parity", "cases", name)))

    out += [
        ("crossed with an item unit",
         _two_group(lambda s: (s["units"].__setitem__("item", {"n": 20}),
                               s.__setitem__("random", {"subject": {"intercept_sd": 1},
                                                        "item": {"intercept_sd": 0.3}})))),
        ("between factor on items",
         _two_group(lambda s: (s["units"].__setitem__("item", {"n": 10}),
                               s["factors"][0].__setitem__("between", "item")))),
        ("subjects nested in sites",
         _two_group(lambda s: s.__setitem__("random", {"site": {"intercept_sd": 0.5,
                                                                "over": "subject", "n": 8}}))),
        ("a within factor", _two_group(lambda s: s["factors"].append(dict(_BLOCK)))),
        ("a factor both between and within",
         _two_group(lambda s: s["factors"][0].__setitem__("vary_within", ["subject"]))),
        ("by-subject random effects and predictors",
         _two_group(lambda s: (s.__setitem__("spec_version", "0.3"),
                               s.__setitem__("random", {"subject": {"intercept_sd": 0.5}}),
                               s.__setitem__("predictors", [
                                   {"name": "age", "varies_by": "subject"},
                                   {"name": "noise", "varies_by": "observation"}]),
                               s["fixed"]["coefficients"].__setitem__("age", 0.2)))),
        ("an item entry without an item unit",
         _two_group(lambda s: s.__setitem__("random", {"item": {"intercept_sd": 0.3}}))),
        ("coefficients {}", _two_group(lambda s: s["fixed"].__setitem__("coefficients", {}))),
        ("three-level between factor",
         _two_group(lambda s: s["factors"][0].update(levels=["a", "b", "c"],
                                                     contrasts={"grp": [-1, 0, 1]}))),
        ("two between factors", _two_group(lambda s: s["factors"].append(dict(_DOSE)))),
        ("two between factors, declaring 0.4",
         _two_group(lambda s: (s.__setitem__("spec_version", "0.4"),
                               s["factors"].append(dict(_DOSE))))),
        ("no factor", _two_group(lambda s: (s.__setitem__("factors", []),
                                            s["fixed"].__setitem__("coefficients", {})))),
        ("lognormal response",
         _two_group(lambda s: s.__setitem__("response", {"family": "lognormal", "name": "RT",
                                                         "sigma": 0.3}))),
    ]
    # Too few subjects for the t-test. R stopped in t.test() at 1 and 2, where Python raised
    # StatisticsError at 1 and returned a power of 0 at 2. Both refuse up to 3 in one message.
    for n in (1, 2, 3, 4, 5):
        out.append(("%d subject%s" % (n, "" if n == 1 else "s"),
                    _two_group(lambda s, n=n: s["units"]["subject"].__setitem__("n", n))))
    return out


# ---- the scale battery --------------------------------------------------------------------

# Every family a specification may declare. A family whose coefficients are on the scale its
# response is fitted on has no warning in either twin.
FAMILIES = ["gaussian", "lognormal", "shifted_lognormal", "exgaussian",
            "bernoulli", "poisson", "ordinal", "beta"]


# ---- the R side ---------------------------------------------------------------------------

R_DRIVER = r'''
args <- commandArgs(trailingOnly = TRUE)
src <- args[1]; mode <- args[2]; payload <- args[3]; out <- args[4]
for (f in sort(list.files(src, pattern = "\\.R$", full.names = TRUE))) source(f)

# The dump of run_r.R, so that the hashes compare like the parity harness's.
.cell <- function(v) if (is.character(v)) v else sprintf("%.17g", as.numeric(v))
.dump <- function(d, path) {
  con <- file(path, open = "wb"); on.exit(close(con))
  lines <- character(nrow(d) + 1L)
  lines[1] <- paste(names(d), collapse = ",")
  cols <- lapply(d, function(col) vapply(col, .cell, character(1), USE.NAMES = FALSE))
  for (i in seq_len(nrow(d)))
    lines[i + 1L] <- paste(vapply(cols, `[`, character(1), i), collapse = ",")
  writeBin(charToRaw(enc2utf8(paste0(paste(lines, collapse = "\n"), "\n"))), con)
}

# A refusal from pilotr carries no call, since every one is raised with call. = FALSE; an
# error with a call came from base R or a package beneath, which is a crash. Warnings are
# collected, joined into one string, so that the spec battery can compare them as text.
run <- function(expr) {
  warns <- character(0)
  res <- tryCatch({
    withCallingHandlers(expr, warning = function(w) {
      warns <<- c(warns, enc2utf8(conditionMessage(w)))
      invokeRestart("muffleWarning")
    })
    list(verdict = "OK", message = "")
  }, error = function(e) list(verdict = if (is.null(conditionCall(e))) "ERROR" else "CRASH",
                              message = enc2utf8(conditionMessage(e))))
  res$warnings <- paste(warns, collapse = "\n")
  res
}

if (mode %in% c("spec", "power")) {
  specs <- jsonlite::fromJSON(payload, simplifyVector = TRUE, simplifyDataFrame = FALSE,
                              simplifyMatrix = FALSE)
  # The power mode takes the first steps of power_design(): validation, then the checks that
  # refuse a design the two-group t-test cannot analyse or has too few subjects for.
  res <- if (mode == "spec") lapply(specs, function(s) run(validate_spec(s, strict = TRUE)))
         else lapply(specs, function(s) run({
           s <- .as_spec(s)
           refusal <- .two_group_refusal(s)
           if (is.null(refusal)) refusal <- .two_group_size_refusal(s)
           if (!is.null(refusal)) stop(refusal, call. = FALSE)
         }))
} else if (mode == "scale") {
  # The wording power_mixed() gives the warning for a family fitted on a scale other than its
  # coefficients', compared as warning text.
  families <- jsonlite::fromJSON(payload)
  res <- lapply(families, function(f) list(
    verdict = "OK", message = "",
    warnings = if (f %in% names(.LINK_SCALE)) sprintf(.LINEAR_POWER, f, .LINK_SCALE[[f]]) else ""))
} else {
  paths <- readLines(payload, encoding = "UTF-8")
  res <- lapply(paths, function(p) run(.dump(simulate_design(load_spec(p)), paste0(p, ".r.txt"))))
}
con <- file(out, open = "wb")
writeBin(charToRaw(enc2utf8(as.character(jsonlite::toJSON(res, auto_unbox = TRUE)))), con)
close(con)
'''


def _run_r(td, mode, payload):
    driver = os.path.join(td, "driver.R")
    with open(driver, "w", encoding="utf-8") as f:
        f.write(R_DRIVER)
    out = os.path.join(td, "r_%s.json" % mode)
    proc = subprocess.run(
        [RSCRIPT, driver, os.path.join(ROOT, "r", "pilotr", "R"), mode, payload, out],
        capture_output=True, text=True, encoding="utf-8", errors="replace")
    if proc.returncode != 0 or not os.path.exists(out):
        raise SystemExit("R driver failed:\n" + proc.stdout + proc.stderr)
    with open(out, encoding="utf-8") as f:
        return json.load(f)


# ---- the Python side ----------------------------------------------------------------------

def _py(fn, refusals=(ValueError,)):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            fn()
            res = {"verdict": "OK", "message": ""}
        except refusals as e:
            res = {"verdict": "ERROR", "message": str(e)}
        except Exception as e:  # noqa: BLE001  anything else is the crash looked for
            res = {"verdict": "CRASH", "message": "%s: %s" % (type(e).__name__, e)}
    res["warnings"] = "\n".join(str(w.message) for w in caught)
    return res


def _py_power(spec):
    """The first steps of power(): validation, then the two-group checks, which raise
    NotImplementedError for a design the t-test cannot analyse and ValueError for one with too
    few subjects."""
    def check():
        s = _as_spec(copy.deepcopy(spec))
        refusal = _two_group_refusal(s)
        if refusal is not None:
            raise NotImplementedError(refusal)
        refusal = _two_group_size_refusal(s)
        if refusal is not None:
            raise ValueError(refusal)
    return _py(check, refusals=(ValueError, NotImplementedError))


def _py_scale(family):
    """The wording of power_mixed()'s warning for `family`, empty for a family it fits on the
    scale of its coefficients."""
    text = _LINEAR_ON_LINK % (family, _LINK_SCALE[family]) if family in _LINK_SCALE else ""
    return {"verdict": "OK", "message": "", "warnings": text}


def _sha256(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


# ---- comparison ---------------------------------------------------------------------------

def _compare(rv, pv, data=None, warned=False):
    """Return a list of problems with one case; empty when the twins agree.

    With `warned`, the warnings each twin raised must also match character for character.
    """
    problems = []
    if "CRASH" in (rv["verdict"], pv["verdict"]):
        problems.append("crash")
    elif rv["verdict"] != pv["verdict"]:
        problems.append("verdict")
    elif rv["verdict"] == "ERROR" and rv["message"] != pv["message"]:
        problems.append("text")
    if warned and not problems and rv.get("warnings", "") != pv.get("warnings", ""):
        problems.append("warning text")
    if data is not None and not problems and rv["verdict"] == "OK":
        r_dump, p_dump = data
        if not os.path.exists(r_dump) or not os.path.exists(p_dump):
            problems.append("data missing")
        elif _sha256(r_dump) != _sha256(p_dump):
            problems.append("data")
    return problems


def _report(title, labels, r_res, p_res, data=None, warned=False):
    print("\n" + title)
    print("%-52s %-6s %-6s %s" % ("case", "R", "Python", ""))
    print("-" * 96)
    failed = 0
    for i, label in enumerate(labels):
        rv, pv = r_res[i], p_res[i]
        problems = _compare(rv, pv, data[i] if data else None, warned)
        flag = "<== " + ", ".join(problems) if problems else ""
        print("%-52s %-6s %-6s %s" % (label[:52], rv["verdict"], pv["verdict"], flag))
        if problems:
            failed += 1
            for who, v in (("R", rv), ("Python", pv)):
                for kind in ("message", "warnings"):
                    if v.get(kind):
                        print("    %s %s: %s" % (who, kind, v[kind].replace("\n", "\n      ")))
    n_ok = sum(1 for v in p_res if v["verdict"] == "OK")
    print("%d cases: %d accepted by Python, %d refused, %d disagreements"
          % (len(labels), n_ok, len(labels) - n_ok, failed))
    return failed


def _schema_check(labels, specs, r_res, p_res):
    """Count the specifications both twins accept with strict = TRUE that the schema refuses."""
    with open(os.path.join(ROOT, "spec", "design.schema.json"), encoding="utf-8") as f:
        schema = json.load(f)
    jsonschema.Draft7Validator.check_schema(schema)
    validator = jsonschema.Draft7Validator(schema)
    print("\nSchema check: spec/design.schema.json on every case both validators accept")
    print("-" * 96)
    failed = checked = 0
    for label, spec, rv, pv in zip(labels, specs, r_res, p_res):
        if rv["verdict"] != "OK" or pv["verdict"] != "OK":
            continue
        checked += 1
        errors = sorted(e.message for e in validator.iter_errors(spec))
        if errors:
            failed += 1
            print("%-52s <== refused by the schema" % label[:52])
            for message in errors:
                print("    " + message)
    print("%d cases accepted by both validators: %d valid under the schema, %d refused"
          % (checked, checked - failed, failed))
    return failed


def main() -> int:
    if jsonschema is None:
        print("the schema check needs the jsonschema package: pip install jsonschema")
        return 1
    battery = cases()
    files = file_cases()
    powers = power_cases()
    py_spec = [_py(lambda s=s: validate_spec(copy.deepcopy(s), strict=True)) for _l, s in battery]
    py_power = [_py_power(s) for _l, s in powers]

    with tempfile.TemporaryDirectory() as td:
        payload = os.path.join(td, "specs.json")
        with open(payload, "w", encoding="utf-8") as f:
            json.dump([s for _l, s in battery], f)
        r_spec = _run_r(td, "spec", payload)
        power_payload = os.path.join(td, "power_specs.json")
        with open(power_payload, "w", encoding="utf-8") as f:
            json.dump([s for _l, s in powers], f)
        r_power = _run_r(td, "power", power_payload)
        scale_payload = os.path.join(td, "families.json")
        with open(scale_payload, "w", encoding="utf-8") as f:
            json.dump(FAMILIES, f)
        r_scale = _run_r(td, "scale", scale_payload)

        paths, data = [], []
        for i, (_label, raw) in enumerate(files, start=1):
            path = os.path.join(td, "case%02d.json" % i)
            with open(path, "wb") as f:
                f.write(raw)
            paths.append(path)
            data.append((path + ".r.txt", path + ".py.txt"))
        py_files = [_py(lambda p=p: _dump(simulate(load_spec(p)), p + ".py.txt")) for p in paths]
        listing = os.path.join(td, "files.txt")
        with open(listing, "w", encoding="utf-8") as f:
            f.write("\n".join(p.replace("\\", "/") for p in paths) + "\n")
        r_files = _run_r(td, "file", listing)

        if (len(r_spec) != len(battery) or len(r_files) != len(files)
                or len(r_power) != len(powers) or len(r_scale) != len(FAMILIES)):
            print("R returned %d, %d, %d and %d results for %d, %d, %d and %d cases"
                  % (len(r_spec), len(r_files), len(r_power), len(r_scale),
                     len(battery), len(files), len(powers), len(FAMILIES)))
            return 1
        failed = _report("Spec battery: validate_spec() on the same parsed specification",
                         [label for label, _s in battery], r_spec, py_spec, warned=True)
        failed += _report("File battery: load_spec() on the same bytes, then simulate()",
                          [label for label, _r in files], r_files, py_files, data)
        failed += _report("Power battery: the two-group check of power_design() and power()",
                          [label for label, _s in powers], r_power, py_power)
        failed += _report("Scale battery: power_mixed()'s warning for a family on a link scale",
                          FAMILIES, r_scale, [_py_scale(f) for f in FAMILIES], warned=True)
        failed += _schema_check([label for label, _s in battery], [s for _l, s in battery],
                                r_spec, py_spec)

    print("\n%d disagreement%s in all" % (failed, "" if failed == 1 else "s"))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
