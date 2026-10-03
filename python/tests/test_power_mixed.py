"""`power_mixed`: the scale it analyses a response on, and what its fixed model leaves out.

The statsmodels backend fits one model whatever the specification declares, `yv ~ cc` with
independent by-subject and by-item intercept and slope components, and tests the first contrast of
the single within factor. It fitted a `lognormal` response on its raw scale, while the R twin and
its own `shifted_lognormal` path analyse the log, and it said nothing when a specification declared
terms that its model leaves out.
"""
import os, sys, warnings
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pytest

pytest.importorskip("statsmodels")
pytest.importorskip("pandas")

from pilotr import load_spec, pilotr_example, power_mixed, replicate_seeds, simulate

_LEFT_OUT = ("power_mixed() in Python fits 'yv ~ cc' with independent by-subject and by-item "
             "intercept and slope components, testing the first contrast of the within factor "
             "only. This specification also declares %s, which the fitted model leaves out. The "
             "R package's power_mixed() fits the model the specification implies.")


def _left_out_warnings(fn):
    """Call `fn` and return its result with the messages of the scope warnings it raised."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = fn()
    return result, [str(w.message) for w in caught
                    if issubclass(w.category, UserWarning)
                    and str(w.message).startswith("power_mixed() in Python fits")]


def _priming(family, **response):
    """The 12 x 8 crossed design of the Python documentation, with random intercepts only."""
    return {
        "name": "priming", "seed": 1,
        "units": {"subject": {"n": 12}, "item": {"n": 8}},
        "factors": [{"name": "condition", "levels": ["related", "unrelated"],
                     "contrasts": {"cond": [-0.5, 0.5]}, "vary_within": ["subject", "item"]}],
        "fixed": {"intercept": 6.0, "coefficients": {"cond": 0.1}},
        "random": {"subject": {"intercept_sd": 0.12}, "item": {"intercept_sd": 0.08}},
        "response": dict({"family": family, "name": "RT", "sigma": 0.3, "round": 4}, **response),
    }


def test_a_lognormal_response_is_analysed_on_the_log_scale():
    # Both families draw the same eta + sigma * z, so on the log scale, log(y) for one and
    # log(y - 200) for the other, they give one estimate. The lognormal was fitted on its raw
    # scale, in milliseconds: a mean estimate of 42.3 against a true 0.1, and Type M 423.
    ln, ln_warned = _left_out_warnings(lambda: power_mixed(_priming("lognormal"), n_sims=2))
    sh, sh_warned = _left_out_warnings(
        lambda: power_mixed(_priming("shifted_lognormal", shift=200), n_sims=2))
    assert ln["mean_estimate"] == pytest.approx(sh["mean_estimate"], abs=1e-6)
    assert ln["n_significant"] == sh["n_significant"]
    assert ln["type_m"] == pytest.approx(sh["type_m"], abs=1e-5)
    assert abs(ln["mean_estimate"] - ln["true_effect"]) < 0.1
    # Random intercepts and nothing else, so the fitted model leaves nothing out.
    assert ln_warned == [] and sh_warned == []


@pytest.mark.parametrize("family, response", [
    ("lognormal", {}),
    ("shifted_lognormal", {"shift": 100}),
])
def test_a_response_with_no_logarithm_leaves_its_replicate_without_an_estimate(family, response):
    # Rounded to whole units, a response near exp(1.5) is sometimes 0, or the shift itself, and
    # has no logarithm. The R twin's log() gives -Inf there, which lmer() refuses, so R counts
    # the replicate as returning no fit. math.log() raised ValueError here and ended the call.
    spec = _priming(family, **response)
    spec.update(seed=4, units={"subject": {"n": 6}, "item": {"n": 4}})
    spec["fixed"]["intercept"] = 1.5
    spec["response"].update(sigma=1, round=0)
    shift = response.get("shift", 0)
    no_log = [any(v - shift <= 0 for v in simulate(dict(spec, seed=s)).column("RT"))
              for s in replicate_seeds(spec["seed"], 4)]
    assert no_log == [False, False, True, True]
    r = power_mixed(spec, n_sims=4)
    assert r["n_sims"] == 4 and r["n_converged"] == 2


def test_the_warning_names_what_reading_time_continuous_declares_beyond_the_model():
    # The shipped example declares six further coefficients, three predictors and two by-subject
    # slopes, none of which `yv ~ cc` contains. Its unit counts are cut to keep the fits quick,
    # and the declared structure, which is all the warning reads, is untouched. Two replicates
    # still give one warning.
    spec = load_spec(pilotr_example("reading_time_continuous"))
    spec["units"] = {"subject": {"n": 8}, "item": {"n": 6}}
    _, warned = _left_out_warnings(lambda: power_mixed(spec, n_sims=2))
    assert warned == [_LEFT_OUT % (
        "other coefficients ('SyntaxPC', 'CoherencePC', 'age', 'SyntaxPC:age', "
        "'CoherencePC:age', 'cond:age'), predictors ('SyntaxPC', 'CoherencePC', 'age') and "
        "by-subject slopes ('SyntaxPC', 'CoherencePC')")]


def _small():
    """A 6 x 4 crossed design whose specification matches the fitted model term for term."""
    return {
        "spec_version": "0.3", "name": "small", "seed": 5,  # 0.3 for the `correlated` field
        "units": {"subject": {"n": 6}, "item": {"n": 4}},
        "factors": [{"name": "condition", "levels": ["a", "b", "c"],
                     "contrasts": {"cond": [-1, 0, 1], "quad": [1, -2, 1]},
                     "vary_within": ["subject", "item"]}],
        "fixed": {"intercept": 0, "coefficients": {"cond": 0.3}},
        "random": {"subject": {"intercept_sd": 0.3, "slopes": {"cond": 0.1}},
                   "item": {"intercept_sd": 0.2, "slopes": {"cond": 0.1}}},
        "response": {"family": "gaussian", "name": "y", "sigma": 1},
    }


def _second_contrast(s, value=0.2):
    s["fixed"]["coefficients"]["quad"] = value


def _zero_coefficient(s):
    _second_contrast(s, 0)


def _between_factors(s):
    s["factors"] += [{"name": "group", "levels": ["x", "y"], "contrasts": {"grp": [-0.5, 0.5]},
                      "between": "subject"},
                     {"name": "set", "levels": ["p", "q"], "contrasts": {"st": [-0.5, 0.5]},
                      "between": "item"}]


def _predictor_with_effect(s):
    s["predictors"] = [{"name": "age", "varies_by": "subject"}]
    s["fixed"]["coefficients"]["age"] = 0.1


def _by_item_slope(s):
    s["predictors"] = [{"name": "age", "varies_by": "subject"}]
    s["random"]["item"]["slopes"]["age"] = 0.05


def _correlations(s):
    s["random"]["subject"]["correlations"] = {"intercept,cond": 0.3}


def _correlated_flag(s):
    s["random"]["item"]["correlated"] = True


def _uncorrelated_flag(s):
    s["random"]["subject"]["correlated"] = False


def _extra_groups(s):
    s["random"]["site"] = {"over": "subject", "n": 3, "intercept_sd": 0.2}
    s["random"]["block"] = {"over": "item", "n": 2, "intercept_sd": 0.1}


@pytest.mark.parametrize("change, left_out", [
    (_second_contrast, "another coefficient ('quad')"),
    (_zero_coefficient, None),
    (_between_factors, "between factors ('group', 'set')"),
    (_predictor_with_effect, "another coefficient ('age') and a predictor ('age')"),
    (_by_item_slope, "a predictor ('age') and a by-item slope ('age')"),
    (_correlations, "random-effect correlations ('subject')"),
    (_correlated_flag, "random-effect correlations ('item')"),
    (_uncorrelated_flag, None),
    (_extra_groups, "extra grouping factors ('site', 'block')"),
])
def test_each_term_the_fitted_model_leaves_out_is_named(change, left_out):
    # A zero coefficient changes no data, and `correlated: false` is what the independent
    # components assume, so neither is reported.
    spec = _small()
    change(spec)
    _, warned = _left_out_warnings(lambda: power_mixed(spec, n_sims=1))
    assert warned == ([] if left_out is None else [_LEFT_OUT % left_out])
