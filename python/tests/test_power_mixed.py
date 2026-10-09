"""`power_mixed`: the scale it analyses a response on, what its fixed model leaves out, whether
its fits reach the optimum and how it counts them.

The statsmodels backend fits one model whatever the specification declares, `yv ~ cc` with
independent by-subject and by-item intercept and slope components, and tests the first contrast of
the single within factor. It fitted a `lognormal` response on its raw scale, while the R twin and
its own `shifted_lognormal` path analyse the log, and it said nothing when a specification declared
terms that its model leaves out. Its fits stopped short of the REML optimum in most crossed
replicates, and it counted every fit that returned as converged. For an accuracy, count, ordinal or
proportion response, it compared an estimate on the response scale with a true value on the link
scale, and reported the mean estimate and Type M of that comparison without a word.
"""
import math, os, re, sys, warnings
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pytest

pytest.importorskip("statsmodels")
pytest.importorskip("pandas")

import numpy as np
from statsmodels.regression.mixed_linear_model import MixedLM
from statsmodels.tools.sm_exceptions import ConvergenceWarning

from pilotr import load_spec, pilotr_example, power_mixed, replicate_seeds, simulate

_LEFT_OUT = ("power_mixed() in Python fits 'yv ~ cc' with independent by-subject and by-item "
             "intercept and slope components, testing the first contrast of the within factor "
             "only. This specification also declares %s, which the fitted model leaves out. The "
             "R package's power_mixed() fits the model the specification implies.")

# lme4's REML fit of (1 + cond || subject) + (1 + cond || item), the uncorrelated model the backend
# fits, to replicate 0 of the shipped crossed_mixed_rt example, as tools/parity/lme4_reference.R
# prints them (lme4 2.0.1). Rerun that script and update both values whenever a change to the
# simulator or to the example moves that replicate's data.
_LME4_SE = 0.018597
_LME4_REML_LL = -389.7568


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
    assert r["n_attempted"] == 4 and r["n_returned"] == 2


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


@pytest.fixture(scope="module")
def crossed_rt():
    """`power_mixed` over the first replicate of crossed_mixed_rt, with the fits it made."""
    real, fits = MixedLM.fit, []

    def recording(self, *args, **kwargs):
        m = real(self, *args, **kwargs)
        fits.append(m)
        return m

    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(MixedLM, "fit", recording)
        # The example declares random-effect correlations, which the fitted model leaves out.
        with pytest.warns(UserWarning, match="random-effect correlations"):
            result = power_mixed(load_spec(pilotr_example("crossed_mixed_rt")), n_sims=1)
    return result, fits


def test_the_fit_reaches_the_reml_optimum_that_lme4_reaches(crossed_rt):
    # statsmodels' default chain of optimisers, BFGS, L-BFGS and CG, stopped short of the optimum
    # on this replicate, at a REML log-likelihood of -406.4 and a standard error of 0.0345. The
    # documentation had put the gap down to statsmodels overstating random-slope variance.
    _, fits = crossed_rt
    assert fits[-1].llf == pytest.approx(_LME4_REML_LL, abs=0.01)
    assert fits[-1].bse_fe["cc"] == pytest.approx(_LME4_SE, rel=0.01)


def test_the_result_names_the_test_behind_its_p_values(crossed_rt):
    result, _ = crossed_rt
    assert result["backend"] == (
        "statsmodels MixedLM (crossed variance components, REML; Wald z tests)")


def test_a_fit_with_a_variance_component_at_zero_is_singular_and_not_converged():
    # The specification has random intercepts alone, so the slope components the backend always
    # fits are zero in truth, and most fits put one of them on the boundary. Every fit that returned
    # used to count as converged, where the R twin counts singular fits apart from clean ones.
    r = power_mixed(_priming("lognormal"), n_sims=6)
    assert (r["n_attempted"], r["n_returned"]) == (6, 6)
    assert r["n_singular"] > 0
    # statsmodels flags all six fits as possibly on the boundary, as it does any fit with a
    # variance component below 0.01, so its notice cannot be what decides convergence.
    assert 0 < r["n_converged"] < r["n_returned"]


def _clean():
    """A 12 x 8 Gaussian design whose fits converge without a notice of any kind."""
    return {
        "name": "clean", "seed": 1,
        "units": {"subject": {"n": 12}, "item": {"n": 8}},
        "factors": [{"name": "condition", "levels": ["a", "b"],
                     "contrasts": {"cond": [-0.5, 0.5]}, "vary_within": ["subject", "item"]}],
        "fixed": {"intercept": 0, "coefficients": {"cond": 0.5}},
        "random": {"subject": {"intercept_sd": 1, "slopes": {"cond": 0.8}},
                   "item": {"intercept_sd": 1, "slopes": {"cond": 0.8}}},
        "response": {"family": "gaussian", "name": "y", "sigma": 1},
    }


_HESSIAN = "The Hessian matrix at the estimated parameter values is not positive definite."


def _notice_after_each_fit(monkeypatch, notice):
    """Make every MixedLM fit raise `notice` as a ConvergenceWarning once it has returned."""
    real = MixedLM.fit

    def fit_and_notice(self, *args, **kwargs):
        m = real(self, *args, **kwargs)
        warnings.warn(notice, ConvergenceWarning, stacklevel=2)  # where statsmodels raises it
        return m

    monkeypatch.setattr(MixedLM, "fit", fit_and_notice)


@pytest.mark.parametrize("notice, counts", [
    ("The MLE may be on the boundary of the parameter space.", False),
    ("Retrying MixedLM optimization with lbfgs", False),
    ("Maximum Likelihood optimization failed to converge. Check mle_retvals", False),
    ("MixedLM optimization failed, trying a different optimizer may help.", True),
    ("Gradient optimization failed, |grad| = 0.104300", True),
    (_HESSIAN, True),
])
def test_only_a_notice_against_the_fit_returned_makes_it_a_fit_with_a_warning(
        monkeypatch, notice, counts):
    # statsmodels raises the first three about fits that may be at the optimum: the first whenever
    # a variance component is below 0.01, the other two whenever an optimiser in its chain gives
    # way to the next, which may then converge. The last three say that the fit returned is not at
    # a proper optimum. With nothing added, both fits of this design converge and neither is
    # singular or warned, so the Hessian notice counts here (see the next test).
    _notice_after_each_fit(monkeypatch, notice)
    r = power_mixed(_clean(), n_sims=2)
    assert (r["n_returned"], r["n_warning"], r["n_converged"]) == (
        (2, 2, 0) if counts else (2, 0, 2))


def test_the_hessian_notice_counts_only_against_a_fit_that_is_not_singular(monkeypatch):
    # statsmodels checks the diagonal of the Hessian with respect to the variances, and a
    # component held at zero can fail that check by itself. On the power guide's example, 3 of the
    # 12 fits, all singular, drew the notice in that way. Counted there, it made them fits with a
    # warning, while the R twin's power_mixed(), given the same model, warns about none of them.
    # Added to every fit of a design whose fits are partly singular, the notice counts against
    # each fit that is not singular and against no other.
    _notice_after_each_fit(monkeypatch, _HESSIAN)
    r = power_mixed(_priming("lognormal"), n_sims=6)
    assert 0 < r["n_singular"] < r["n_returned"]
    assert r["n_warning"] == r["n_returned"] - r["n_singular"]
    assert r["n_converged"] == 0


def test_a_replicate_whose_powell_fit_raises_is_refitted_and_counts_as_a_warning(monkeypatch):
    # Powell's method can raise where the default chain succeeds, as LinAlgError did in one of
    # 300 null fits at 8 x 6. The replicate is fitted again with the default chain, and that
    # fallback makes it a fit with a warning.
    real = MixedLM.fit

    def powell_raises(self, *args, method=None, **kwargs):
        if method is not None and "powell" in method:
            raise np.linalg.LinAlgError("Singular matrix")
        return real(self, *args, method=method, **kwargs)

    monkeypatch.setattr(MixedLM, "fit", powell_raises)
    r = power_mixed(_clean(), n_sims=2)
    assert (r["n_returned"], r["n_warning"], r["n_converged"]) == (2, 2, 0)

    # When the default chain raises too, the replicate returns no fit.
    def always_raises(self, *args, **kwargs):
        raise np.linalg.LinAlgError("Singular matrix")

    monkeypatch.setattr(MixedLM, "fit", always_raises)
    r = power_mixed(_clean(), n_sims=2)
    assert (r["n_attempted"], r["n_returned"]) == (2, 0)
    assert math.isnan(r["power"])


# The warning for a family whose coefficients are on a link scale, word for word as the R twin's
# power_mixed() raises it for the ordinal and Beta families.
_LINEAR = ("power_mixed() fits a linear model to the %s response on its own scale, while the "
           "specification's coefficients are on the %s scale, so the mean estimate and Type M are "
           "withheld (NA) and Type S compares signs only. For a model on the link scale, use "
           "generate_design_analysis() or brms_bridge() in the R package.")


def _one_within(family, **response):
    """A 12 x 8 crossed design with one within factor and random intercepts, in `family`."""
    return {
        "name": family, "seed": 3,
        "units": {"subject": {"n": 12}, "item": {"n": 8}},
        "factors": [{"name": "condition", "levels": ["a", "b"],
                     "contrasts": {"cond": [-0.5, 0.5]}, "vary_within": ["subject", "item"]}],
        "fixed": {"intercept": 0.5, "coefficients": {"cond": 0.8}},
        "random": {"subject": {"intercept_sd": 0.5}, "item": {"intercept_sd": 0.3}},
        "response": dict({"family": family, "name": "y"}, **response),
    }


def _scale_warnings(fn):
    """Call `fn` and return its result with the messages of the scale warnings it raised."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = fn()
    return result, [str(w.message) for w in caught
                    if issubclass(w.category, UserWarning)
                    and str(w.message).startswith("power_mixed() fits a linear model")]


def test_a_bernoulli_response_warns_and_withholds_the_mean_estimate_and_type_m():
    # A linear model of the 0/1 response estimates a difference in probability, which was set
    # against a true value on the logit scale, so the mean estimate and Type M compared two
    # scales. On the R twin's crossed accuracy design, a linear mixed model gave a Type M of 0.21
    # where a logistic one gives 1.12, so a design that exaggerates read as one that underestimates.
    with pytest.warns(UserWarning, match=re.escape(_LINEAR % ("bernoulli", "logit"))):
        r = power_mixed(_one_within("bernoulli"), n_sims=2)
    assert math.isnan(r["mean_estimate"]) and math.isnan(r["type_m"])
    assert r["n_returned"] == 2 and not math.isnan(r["power"])


@pytest.mark.parametrize("family, scale, response", [
    ("bernoulli", "logit", {}),
    ("poisson", "log", {}),
    ("ordinal", "logit", {"thresholds": [-1, 0, 1]}),
    ("beta", "logit", {"phi": 8}),
    ("gaussian", None, {"sigma": 1}),
])
def test_each_family_on_a_link_scale_warns_once_with_its_scale(family, scale, response):
    r, warned = _scale_warnings(lambda: power_mixed(_one_within(family, **response), n_sims=1))
    if scale is None:
        # The response's own scale is the coefficients' scale, so nothing is withheld.
        assert warned == [] and not math.isnan(r["mean_estimate"])
    else:
        assert warned == [_LINEAR % (family, scale)]
        assert math.isnan(r["mean_estimate"]) and math.isnan(r["type_m"])
