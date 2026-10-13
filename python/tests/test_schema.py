"""The portable schema against strict validation.

spec/design.schema.json was looser than validate_spec(). It accepted an unknown top-level field,
sigma or round on a Bernoulli response and an ordinal or Gaussian response without its
parameters, all of which strict validation refuses. It also refused the single-value forms of
vary_within and thresholds that strict validation accepts. These tests read the schema from the
repository around the package and skip outside it. The structural checks compare the schema with
the validator's own tables and need nothing else. The behavioural ones put specifications through
jsonschema's Draft7Validator and skip without it. Twinned with test-contract.R, which compares the
schema with the R validator's tables.
"""

import copy
import json
import os
import sys
import warnings

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pytest

from pilotr import pilotr_example, validate_spec
from pilotr.validate import FAMILY_PARAMS, ROUNDING_FAMILIES

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, "..", ".."))
SCHEMA = os.path.join(REPO, "spec", "design.schema.json")

# Every response parameter some family reads. A family refuses the ones it does not read.
PARAMS = ("sigma", "shift", "beta", "phi", "thresholds")


def _schema():
    if not os.path.exists(SCHEMA):
        pytest.skip("spec/design.schema.json not available outside the repository")
    with open(SCHEMA, encoding="utf-8") as f:
        return json.load(f)


def _family_rules(schema):
    """Each family's required and refused response fields, read off the schema's if/then blocks."""
    rules = {}
    for block in schema["properties"]["response"].get("allOf", []):
        cond = block.get("if", {}).get("properties", {}).get("family", {})
        families = [cond["const"]] if "const" in cond else list(cond.get("enum", []))
        then = block.get("then", {})
        for fam in families:
            r = rules.setdefault(fam, {"required": set(), "refused": set(), "round_min": 0})
            r["required"].update(then.get("required", []))
            for key, sub in then.get("properties", {}).items():
                if sub is False:
                    r["refused"].add(key)
                elif key == "round":
                    r["round_min"] = sub.get("minimum", 0)
    return rules


def test_the_schema_refuses_unknown_top_level_fields():
    assert _schema()["additionalProperties"] is False


def test_each_family_requires_and_refuses_what_the_validator_does():
    rules = _family_rules(_schema())
    assert set(rules) == set(FAMILY_PARAMS)
    for fam, needed in FAMILY_PARAMS.items():
        refused = set(PARAMS) - set(needed)
        if fam not in ROUNDING_FAMILIES:
            refused.add("round")
        assert rules[fam]["required"] == set(needed), fam
        assert rules[fam]["refused"] == refused, fam
    # A proportion rounded to whole numbers is 0 or 1, outside the family's support.
    assert rules["beta"]["round_min"] == 1


# ---- behaviour, through jsonschema ----------------------------------------------------------

def _example(name):
    with open(pilotr_example(name), encoding="utf-8") as f:
        return json.load(f)


def _with(name, fn):
    s = _example(name)
    fn(s)
    return s


def _response(**kw):
    return lambda s: s.__setitem__("response", kw)


def _predictor(**kw):
    def fn(s):
        s["spec_version"] = "0.3"
        s["predictors"] = [dict({"name": "z", "varies_by": "subject"}, **kw)]
        s["fixed"]["coefficients"]["z"] = 0.1
    return fn


def _site(**kw):
    return lambda s: s["random"].__setitem__("site", dict({"intercept_sd": 0.5}, **kw))


RT, TWO = "crossed_mixed_rt", "between_2group_gaussian"

# (label, specification, whether strict validation accepts it)
CASES = [(name, _example(name), True) for name in pilotr_example()] + [
    ("a single string for vary_within",
     _with(RT, lambda s: s["factors"][0].__setitem__("vary_within", "item")), True),
    ("a single threshold as a bare number",
     _with(TWO, _response(family="ordinal", name="r", thresholds=0.5)), True),
    ("bernoulli", _with(TWO, _response(family="bernoulli", name="acc")), True),
    ("exgaussian declaring 0.3",
     _with(TWO, lambda s: (s.__setitem__("spec_version", "0.3"),
                           _response(family="exgaussian", name="RT", sigma=40, beta=80)(s))),
     True),
    ("beta rounded to one place", _with(TWO, _response(family="beta", name="p", phi=8, round=1)),
     True),
    ("a grouping factor with over and n", _with(RT, _site(over="subject", n=4)), True),
    ("a uniform predictor with min and max", _with(RT, _predictor(dist="uniform", min=0, max=1)),
     True),
    ("an unknown top-level field", _with(TWO, lambda s: s.__setitem__("notes", "x")), False),
    ("sigma on bernoulli", _with(TWO, _response(family="bernoulli", name="acc", sigma=1)), False),
    ("round on bernoulli", _with(TWO, _response(family="bernoulli", name="acc", round=2)), False),
    ("round on poisson", _with(TWO, _response(family="poisson", name="n", round=0)), False),
    ("ordinal without thresholds", _with(TWO, _response(family="ordinal", name="r")), False),
    ("gaussian without sigma", _with(TWO, _response(family="gaussian", name="y")), False),
    ("shifted_lognormal without shift",
     _with(TWO, _response(family="shifted_lognormal", name="RT", sigma=0.3)), False),
    ("exgaussian without beta",
     _with(TWO, lambda s: (s.__setitem__("spec_version", "0.3"),
                           _response(family="exgaussian", name="RT", sigma=40)(s))), False),
    ("beta without phi", _with(TWO, _response(family="beta", name="p")), False),
    ("shift on gaussian", _with(TWO, _response(family="gaussian", name="y", sigma=1, shift=200)),
     False),
    ("phi on lognormal", _with(TWO, _response(family="lognormal", name="RT", sigma=0.3, phi=8)),
     False),
    ("beta rounded to whole numbers",
     _with(TWO, _response(family="beta", name="p", phi=8, round=0)), False),
    ("per_subject on the subject unit",
     _with(TWO, lambda s: s["units"]["subject"].__setitem__("per_subject", 2)), False),
    ("a grouping factor without over and n", _with(RT, _site()), False),
    ("over on random.subject", _with(RT, lambda s: s["random"]["subject"].__setitem__("over",
                                                                                     "subject")),
     False),
    ("a uniform predictor without min and max", _with(RT, _predictor(dist="uniform")), False),
    ("a mean on a uniform predictor",
     _with(RT, _predictor(dist="uniform", min=0, max=1, mean=0.5)), False),
    ("a min on a normal predictor", _with(RT, _predictor(min=0)), False),
    ("correlated false with correlations",
     _with(RT, lambda s: (s.__setitem__("spec_version", "0.3"),
                          s["random"]["subject"].__setitem__("correlated", False))), False),
    ("an empty vary_within", _with(RT, lambda s: s["factors"][0].__setitem__("vary_within", [])),
     False),
]


def _strict_accepts(spec):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # the single-string vary_within draws a deprecation warning
        try:
            validate_spec(copy.deepcopy(spec), strict=True)
        except ValueError:
            return False
    return True


@pytest.fixture(scope="module")
def validator():
    jsonschema = pytest.importorskip("jsonschema")
    schema = _schema()
    jsonschema.Draft7Validator.check_schema(schema)
    return jsonschema.Draft7Validator(schema)


@pytest.mark.parametrize("label, spec, accepted", CASES, ids=[c[0] for c in CASES])
def test_the_schema_accepts_what_strict_validation_accepts(validator, label, spec, accepted):
    assert _strict_accepts(spec) is accepted, label  # the case is what its label says
    errors = sorted(e.message for e in validator.iter_errors(spec))
    assert validator.is_valid(spec) is accepted, errors


def test_every_parity_case_is_valid_under_the_schema(validator):
    folder = os.path.join(REPO, "tools", "parity", "cases")
    if not os.path.isdir(folder):
        pytest.skip("tools/parity/cases not available outside the repository")
    for name in sorted(os.listdir(folder)):
        if name.endswith(".json"):
            with open(os.path.join(folder, name), encoding="utf-8") as f:
                spec = json.load(f)
            errors = [e.message for e in validator.iter_errors(spec)]
            assert not errors, (name, errors)
