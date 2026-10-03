"""The two-group backend, `power`, and the coefficient lookup it shares with `power_mixed`.

`power` t-tests every row, and that test holds its level only with one row per subject and no
shared cluster. It used to accept any Gaussian design with a two-level between factor, and report
power for a test that did not hold its level. Each refusal below is byte-identical to the R
twin's, which tools/parity/validate_cross.py checks.
"""
import math, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pytest

pytest.importorskip("scipy")  # power() t-tests with scipy, and imports it before anything else

from pilotr import power, power_curve, power_mixed


def _one_row_per_subject(n=40, effect=0.5):
    """A two-group Gaussian design with one row per subject, the design the t-test analyses."""
    return {
        "spec_version": "0.3", "name": "tg", "seed": 11,
        "units": {"subject": {"n": n}},
        "factors": [{"name": "group", "levels": ["a", "b"],
                     "contrasts": {"effect": [-0.5, 0.5]}, "between": "subject"}],
        "fixed": {"intercept": 0, "coefficients": {"effect": effect}},
        "random": {},
        "response": {"family": "gaussian", "name": "y", "sigma": 1, "round": 4},
    }


def _correlated_rows(what):
    return ("power_design() in R and power() in Python t-test every row as an independent "
            "observation, which is valid only when each subject contributes one row and no rows "
            "share a cluster. This design has %s, so its rows are correlated and the t-test would "
            "overstate power and understate Type M. Use power_mixed() in the R package, which "
            "fits the model the specification implies." % what)


def _refusal(fn):
    with pytest.raises(NotImplementedError) as caught:
        fn()
    return str(caught.value)


def test_a_design_crossed_with_items_is_refused():
    # 30 subjects by 20 items with no true effect. Analysed row by row, 105 of 200 replicates
    # were significant, against 16 for a t-test of the subject means on the same data.
    spec = _one_row_per_subject(n=30, effect=0.0)
    spec["units"]["item"] = {"n": 20}
    spec["random"] = {"subject": {"intercept_sd": 1}, "item": {"intercept_sd": 0.3}}
    assert _refusal(lambda: power(spec, n_sims=2)) == _correlated_rows("an item unit")
    # power_curve() runs the same check at its first grid point.
    assert (_refusal(lambda: power_curve(spec, subject_ns=[20, 40], n_sims=2))
            == _correlated_rows("an item unit"))


def test_subjects_nested_in_sites_are_refused_and_the_grouping_factor_is_named():
    # One row per subject, but the rows of a site share its effect: 0.200 under the null.
    spec = _one_row_per_subject(n=120, effect=0.0)
    spec["random"] = {"site": {"intercept_sd": 0.5, "over": "subject", "n": 12}}
    assert (_refusal(lambda: power(spec, n_sims=2))
            == _correlated_rows("the grouping factor 'site'"))


def test_a_within_factor_is_refused_and_named():
    spec = _one_row_per_subject()
    spec["factors"].append({"name": "block", "levels": ["x", "y"],
                            "contrasts": {"blk": [-0.5, 0.5]}, "vary_within": "subject"})
    assert (_refusal(lambda: power(spec, n_sims=2))
            == _correlated_rows("the within factor 'block'"))


def test_predictors_and_by_subject_random_effects_keep_one_row_per_subject():
    spec = _one_row_per_subject()
    spec["random"] = {"subject": {"intercept_sd": 0.5}}
    spec["predictors"] = [{"name": "age", "varies_by": "subject", "mean": 0, "sd": 1},
                          {"name": "noise", "varies_by": "observation", "mean": 0, "sd": 1}]
    spec["fixed"]["coefficients"]["age"] = 0.2
    r = power(spec, n_sims=20)
    assert r["true_effect"] == 0.5
    assert 0.0 <= r["power"] <= 1.0


def test_designs_the_backend_never_covered_keep_their_messages():
    counts = _one_row_per_subject()
    counts["response"] = {"family": "poisson", "name": "count"}
    assert (_refusal(lambda: power(counts, n_sims=2))
            == "The power backend currently handles only the gaussian two-group design.")
    three = _one_row_per_subject()
    three["factors"][0]["levels"] = ["a", "b", "c"]
    three["factors"][0]["contrasts"] = {"effect": [-1, 0, 1]}
    assert (_refusal(lambda: power(three, n_sims=2))
            == "The power backend expects exactly one 2-level between factor.")


def test_a_specification_with_no_coefficient_for_the_contrast_has_a_true_effect_of_0():
    # `"coefficients": {}` is a valid way to write a null design. power() raised KeyError.
    spec = _one_row_per_subject()
    spec["fixed"]["coefficients"] = {}
    r = power(spec, n_sims=20)
    assert r["true_effect"] == 0
    assert math.isnan(r["type_s"])
    assert math.isnan(r["type_m"])
    assert 0.0 <= r["power"] <= 1.0


def test_power_mixed_reads_a_missing_coefficient_as_a_true_effect_of_0():
    pytest.importorskip("statsmodels")
    pytest.importorskip("pandas")
    spec = {
        "name": "null", "seed": 3,
        "units": {"subject": {"n": 10}, "item": {"n": 6}},
        "factors": [{"name": "cond", "levels": ["a", "b"],
                     "contrasts": {"cond": [-0.5, 0.5]}, "vary_within": ["subject", "item"]}],
        "fixed": {"intercept": 6, "coefficients": {}},
        "random": {"subject": {"intercept_sd": 0.12}, "item": {"intercept_sd": 0.08}},
        "response": {"family": "gaussian", "name": "y", "sigma": 0.3},
    }
    r = power_mixed(spec, n_sims=2)
    assert r["true_effect"] == 0
    assert math.isnan(r["type_s"])
    assert math.isnan(r["type_m"])


def test_power_validates_the_specification_once(monkeypatch):
    # The replicates simulate through _simulate(), which skips validation, where each of them used
    # to validate the specification again.
    module = sys.modules["pilotr.simulate"]  # the attribute pilotr.simulate is the function
    real = module.validate_spec
    calls = []

    def counting(*args, **kwargs):
        calls.append(1)
        return real(*args, **kwargs)

    monkeypatch.setattr(module, "validate_spec", counting)
    power(_one_row_per_subject(), n_sims=10)
    assert len(calls) == 1
