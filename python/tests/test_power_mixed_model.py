"""`power_mixed` fits the model the specification implies, and tests every coefficient.

It used to fit one model written into its source, `yv ~ cc` with by-subject and by-item intercept
and slope components, whatever the specification declared. It tested the first contrast of a
single within factor and refused any other design. The model now comes from `model_formula` and
`model_data`, as in the R twin: every coefficient is a fixed effect, and every declared random
intercept and slope, an extra grouping factor's included, is a variance component. Each focal
effect has its own results under `effects`, and the scalar keys of the within contrast are kept for
one release with a FutureWarning.
"""
import math, os, re, sys, warnings
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pytest

pytest.importorskip("statsmodels")
pytest.importorskip("pandas")

from statsmodels.regression.mixed_linear_model import MixedLM

from pilotr import load_spec, pilotr_example, power_mixed

_UNCORRELATED = ("power_mixed() in Python fits each random effect as an independent variance "
                 "component, so it leaves out the random-effect correlations this specification "
                 "declares for %s. The R package's power_mixed() estimates them.")

_EFFECT_KEYS = ["power", "power_mcse", "power_lo", "power_hi", "n_significant", "true_effect",
                "mean_estimate", "type_s", "type_s_lo", "type_s_hi", "type_m", "type_m_mcse"]


def _quiet(fn):
    """Call `fn` with every warning recorded, and return its result with the warnings."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = fn()
    return result, caught


def _fitted_components(monkeypatch):
    """Record the names of the variance components of every model that power_mixed() fits."""
    real, seen = MixedLM.fit, []

    def recording(self, *args, **kwargs):
        seen.append(list(self.exog_vc.names))
        return real(self, *args, **kwargs)

    monkeypatch.setattr(MixedLM, "fit", recording)
    return seen


def test_reading_time_continuous_returns_all_seven_coefficients():
    # The model used to be `yv ~ cc`, so the six other coefficients of this example, its
    # predictors, their interactions and the by-subject slopes on SyntaxPC and CoherencePC were
    # left out. The unit counts are cut to keep the fits quick, and the declared structure is kept.
    spec = load_spec(pilotr_example("reading_time_continuous"))
    spec["units"] = {"subject": {"n": 16}, "item": {"n": 12}}
    r, caught = _quiet(lambda: power_mixed(spec, n_sims=2))
    # Nothing declared is left out, since the example has no correlations.
    assert [str(w.message) for w in caught if issubclass(w.category, UserWarning)] == []
    assert list(r["effects"]) == ["SyntaxPC", "CoherencePC", "age", "cond", "SyntaxPC_age",
                                  "CoherencePC_age", "cond_age"]
    for name, e in r["effects"].items():
        assert list(e) == _EFFECT_KEYS
        key = name.replace("_", ":")
        assert e["true_effect"] == spec["fixed"]["coefficients"][key]
        assert 0 <= e["power"] <= 1 and math.isfinite(e["mean_estimate"])
    assert (r["n_attempted"], r["n_returned"]) == (2, 2)


def test_a_two_by_two_within_design_runs_and_tests_both_effects_and_their_interaction():
    # Two within factors used to raise NotImplementedError.
    spec = {
        "name": "twobytwo", "seed": 7,
        "units": {"subject": {"n": 12}, "item": {"n": 8}},
        "factors": [
            {"name": "A", "levels": ["a1", "a2"], "contrasts": {"a": [-0.5, 0.5]},
             "vary_within": ["subject", "item"]},
            {"name": "B", "levels": ["b1", "b2"], "contrasts": {"b": [-0.5, 0.5]},
             "vary_within": ["subject", "item"]}],
        "fixed": {"intercept": 0, "coefficients": {"a": 0.4, "b": 0.2, "a:b": 0.3}},
        "random": {"subject": {"intercept_sd": 0.5, "slopes": {"a": 0.2}},
                   "item": {"intercept_sd": 0.3}},
        "response": {"family": "gaussian", "name": "y", "sigma": 1},
    }
    r = power_mixed(spec, n_sims=3)
    assert list(r["effects"]) == ["a", "b", "a_b"]
    assert [r["effects"][k]["true_effect"] for k in ("a", "b", "a_b")] == [0.4, 0.2, 0.3]
    assert r["n_returned"] == 3
    # Neither factor is the single within factor, so there are no scalar keys to keep.
    assert "power" not in r


def test_only_the_declared_random_terms_are_fitted(monkeypatch):
    # Random intercepts alone were fitted with slope components as well, and an extra grouping
    # factor such as nested_clusters' site was left out.
    seen = _fitted_components(monkeypatch)
    spec = load_spec(pilotr_example("nested_clusters"))
    spec["units"] = {"subject": {"n": 24}, "item": {"n": 6}}
    spec["random"]["site"]["n"] = 4
    power_mixed(spec, n_sims=1)
    assert seen == [["item", "site", "subject"]]

    seen.clear()
    spec = load_spec(pilotr_example("reading_time_continuous"))
    spec["units"] = {"subject": {"n": 8}, "item": {"n": 6}}
    power_mixed(spec, n_sims=1)
    assert seen == [["item", "subject", "subject:CoherencePC", "subject:SyntaxPC"]]


def _small():
    return {
        "name": "small", "seed": 5,
        "units": {"subject": {"n": 10}, "item": {"n": 6}},
        "factors": [{"name": "condition", "levels": ["a", "b"],
                     "contrasts": {"cond": [-0.5, 0.5]}, "vary_within": ["subject", "item"]}],
        "fixed": {"intercept": 0, "coefficients": {"cond": 0.5}},
        "random": {"subject": {"intercept_sd": 0.5, "slopes": {"cond": 0.2}},
                   "item": {"intercept_sd": 0.3, "slopes": {"cond": 0.1}}},
        "response": {"family": "gaussian", "name": "y", "sigma": 1},
    }


def test_declared_correlations_are_fitted_as_independent_with_a_warning():
    spec = _small()
    spec["random"]["subject"]["correlations"] = {"intercept,cond": 0.3}
    spec["random"]["item"]["correlations"] = {"intercept,cond": -0.2}
    with pytest.warns(UserWarning) as seen:
        power_mixed(spec, n_sims=1)
    assert [str(w.message) for w in seen] == [_UNCORRELATED % "'subject' and 'item'"]


def test_a_design_without_items_is_fitted_from_its_by_subject_terms(monkeypatch):
    # A design without an item unit used to be refused outright.
    seen = _fitted_components(monkeypatch)
    spec = _small()
    spec["units"] = {"subject": {"n": 20}}
    spec["factors"][0]["vary_within"] = "subject"
    del spec["random"]["item"]
    r = power_mixed(spec, n_sims=2)
    assert seen == [["subject", "subject:cond"]] * 2
    assert r["n_returned"] == 2 and list(r["effects"]) == ["cond"]


def test_a_design_without_random_effects_is_refused():
    spec = load_spec(pilotr_example("between_2group_gaussian"))
    with pytest.raises(NotImplementedError, match=re.escape(
            "power_mixed() in Python fits a mixed model, and this specification declares no "
            "random effects. For a two-group Gaussian design with one row per subject, use "
            "power(). The R package's power_mixed() fits a design without random effects with "
            "lm() or glm().")):
        power_mixed(spec, n_sims=1)


def test_focal_names_the_effects_and_may_give_their_true_values():
    spec = _small()
    spec["factors"][0]["levels"] = ["a", "b", "c"]
    spec["factors"][0]["contrasts"] = {"cond": [-1, 0, 1], "quad": [1, -2, 1]}
    spec["fixed"]["coefficients"] = {"cond": 0.5, "quad": 0.2}
    # By default every coefficient, with its true value from the specification.
    default = power_mixed(spec, n_sims=2)
    assert {k: e["true_effect"] for k, e in default["effects"].items()} == {
        "cond": 0.5, "quad": 0.2}
    # Names alone leave the true values unknown, so Type S and Type M are undefined.
    named = power_mixed(spec, n_sims=2, focal=["quad"])
    assert list(named["effects"]) == ["quad"]
    assert math.isnan(named["effects"]["quad"]["true_effect"])
    assert math.isnan(named["effects"]["quad"]["type_s"])
    assert named["effects"]["quad"]["power"] == default["effects"]["quad"]["power"]
    # A dict gives the true values to compare with, here other than the one simulated.
    valued = power_mixed(spec, n_sims=2, focal={"cond": 1.0, "(Intercept)": 0.0})
    assert valued["effects"]["cond"]["true_effect"] == 1.0
    assert valued["effects"]["cond"]["mean_estimate"] == default["effects"]["cond"][
        "mean_estimate"]
    assert "(Intercept)" in valued["effects"]


def test_a_focal_name_the_model_lacks_is_reported_in_the_r_twins_words():
    with pytest.warns(UserWarning, match=re.escape(
            "focal effect 'Cond' never appeared in the fitted model, so its results are empty by "
            "construction; the fitted coefficients are '(Intercept)', 'cond'")):
        r = power_mixed(_small(), n_sims=1, focal=["cond", "Cond"])
    assert math.isnan(r["effects"]["Cond"]["power"]) and r["effects"]["Cond"]["n_significant"] == 0


@pytest.mark.parametrize("focal, message", [
    ([], "this specification has no fixed coefficients to test; name the focal effects explicitly"),
    (["cond", "cond"], "focal names the effect 'cond' twice"),
    ({"cond": ["cond", "quad"]}, "focal must be None, a list of effect names or a dict from "
                                 "effect names to true values"),
    ([1], "focal must be None"),
])
def test_a_malformed_focal_is_refused(focal, message):
    with pytest.raises(ValueError, match=re.escape(message)):
        power_mixed(_small(), n_sims=1, focal=focal)


def test_a_specification_without_coefficients_has_nothing_to_test():
    spec = _small()
    spec["fixed"]["coefficients"] = {}
    with pytest.raises(ValueError, match=re.escape(
            "this specification has no fixed coefficients to test; name the focal effects "
            "explicitly")):
        power_mixed(spec, n_sims=1)


def test_the_scalar_keys_of_the_within_contrast_warn_that_they_are_deprecated():
    r = power_mixed(_small(), n_sims=2)
    effect = r["effects"]["cond"]
    for key in _EFFECT_KEYS:
        with pytest.warns(FutureWarning, match=re.escape(
                "result['%s'] from power_mixed() is deprecated and will be removed in the next "
                "release. Read result['effects']['cond']['%s'], which holds the same value for "
                "'cond', the first contrast of the within factor." % (key, key))):
            value = r[key]
        assert value == effect[key] or (math.isnan(value) and math.isnan(effect[key]))
    with pytest.warns(FutureWarning):
        assert r.get("power") == effect["power"]
    # The keys that are not deprecated are read without a word.
    _, caught = _quiet(lambda: (r["effects"], r["n_returned"], r.get("backend")))
    assert caught == []


def test_an_interaction_on_a_link_scale_has_no_type_s():
    # A bounded response compresses differences near the ends of its range, so an interaction can
    # change sign between the response scale, where the linear model estimates it, and the logit
    # scale of its true value. The R twin withholds its Type S for that reason.
    spec = _small()
    spec["factors"].append({"name": "B", "levels": ["b1", "b2"], "contrasts": {"b": [-0.5, 0.5]},
                            "vary_within": ["subject", "item"]})
    spec["fixed"] = {"intercept": 0.3, "coefficients": {"cond": 0.8, "b": 0.5, "cond:b": 0.6}}
    spec["response"] = {"family": "bernoulli", "name": "y"}
    r, _ = _quiet(lambda: power_mixed(spec, n_sims=2, alpha=1))
    assert all(r["effects"][k]["n_significant"] == 2 for k in ("cond", "b", "cond_b"))
    assert not math.isnan(r["effects"]["cond"]["type_s"])
    assert math.isnan(r["effects"]["cond_b"]["type_s"])
    assert math.isnan(r["effects"]["cond_b"]["type_s_lo"])


def test_a_large_random_effects_design_is_held_sparse_and_gives_the_same_fit(monkeypatch):
    # Dense matrices are quicker for a small design and sparse ones for a large one. The two give
    # the same optimum, so the switch changes the running time and not the result.
    from scipy import sparse
    real, held = MixedLM.fit, []

    def recording(self, *args, **kwargs):
        held.append(sparse.issparse(self.exog_vc.mats[0][0]))
        return real(self, *args, **kwargs)

    monkeypatch.setattr(MixedLM, "fit", recording)
    spec = _small()
    dense = power_mixed(spec, n_sims=2)
    # The module, which the function of the same name hides as an attribute of the package.
    monkeypatch.setattr(sys.modules["pilotr.power"], "_SPARSE_FROM", 0)
    held_sparse = power_mixed(spec, n_sims=2)
    assert held == [False, False, True, True]
    for k in ("mean_estimate", "power"):
        assert held_sparse["effects"]["cond"][k] == pytest.approx(dense["effects"]["cond"][k],
                                                                   rel=1e-6)
