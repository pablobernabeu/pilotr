"""Directional tests: `alternative` in `power`, `power_curve` and `power_mixed`.

A directional prediction is tested one-sided, and each power function took only a two-sided test.
The one-sided p-value comes from the statistic behind the two-sided one, so it is half the
two-sided p-value when the estimate lies in the tested direction. On the same replicates, a
one-sided test at alpha therefore rejects exactly where the two-sided test at 2 alpha rejects with
the estimate in that direction. One-sided power at alpha is then two-sided power at 2 alpha times
one minus Type S at 2 alpha, which the tests below hold in counts, exactly, as the R twin's
tests/testthat/test-alternative.R does.
"""
import math, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pytest

pytest.importorskip("scipy")  # power() t-tests with scipy

from pilotr import power, power_curve, replicate_seeds
from pilotr.power import _power_replicate
from pilotr.simulate import _as_spec


def _small_two_group():
    """A two-group design with a small effect, so that some replicates significant at 2 alpha
    have the wrong sign and the identity is tested in both directions."""
    return {
        "name": "dir", "seed": 5,
        "units": {"subject": {"n": 20}},
        "factors": [{"name": "group", "levels": ["a", "b"],
                     "contrasts": {"effect": [-0.5, 0.5]}, "between": "subject"}],
        "fixed": {"intercept": 0, "coefficients": {"effect": 0.2}},
        "response": {"family": "gaussian", "name": "y", "sigma": 1},
    }


def test_one_sided_power_at_alpha_is_the_two_sided_test_at_twice_alpha_in_one_direction():
    spec = _small_two_group()
    two = power(spec, n_sims=400, alpha=0.1)
    greater = power(spec, n_sims=400, alpha=0.05, alternative="greater")
    less = power(spec, n_sims=400, alpha=0.05, alternative="less")
    n_wrong = round(two["type_s"] * two["n_significant"])
    assert n_wrong > 0
    assert greater["n_significant"] == two["n_significant"] - n_wrong
    assert less["n_significant"] == n_wrong
    assert greater["power"] == pytest.approx(two["power"] * (1 - two["type_s"]), rel=1e-12)
    # Type S and Type M are taken over the replicates significant in the tested direction.
    # Against a positive effect, every one of them has the right sign under "greater" and the
    # wrong one under "less", and the exaggeration ratios of the two directions make up the
    # two-sided set.
    assert greater["type_s"] == 0 and less["type_s"] == 1
    assert (greater["type_m"] * greater["n_significant"] + less["type_m"] * less["n_significant"]
            == pytest.approx(two["type_m"] * two["n_significant"], rel=1e-12))
    assert greater["mean_estimate"] == two["mean_estimate"]


def test_power_records_its_alternative_and_two_sided_stays_the_default():
    spec = _small_two_group()
    default = power(spec, n_sims=30)
    assert default["alternative"] == "two.sided"
    assert default == power(spec, n_sims=30, alternative="two.sided")
    assert power(spec, n_sims=30, alternative="less")["alternative"] == "less"


def test_each_replicate_one_sided_p_is_half_its_two_sided_p_in_its_direction():
    spec = _as_spec(_small_two_group())  # the replicates take a validated specification

    def at(seed, alternative):
        return _power_replicate(seed, spec, "group", "a", "b", "y", alternative)

    for seed in replicate_seeds(spec["seed"], 6):
        est, p = at(seed, "two.sided")
        toward, away = ("greater", "less") if est > 0 else ("less", "greater")
        assert at(seed, toward)[1] == pytest.approx(p / 2, rel=1e-12)
        assert at(seed, away)[1] == pytest.approx(1 - p / 2, rel=1e-12)
        assert at(seed, toward)[0] == est


def test_power_curve_passes_the_alternative_to_every_point():
    spec = _small_two_group()
    curve = power_curve(spec, subject_ns=[20, 40], n_sims=30, alternative="greater")
    two = power_curve(spec, subject_ns=[20, 40], n_sims=30, alpha=0.1)
    for rec, ref in zip(curve, two):
        n_wrong = round(ref["type_s"] * ref["n_significant"])
        assert rec["n_significant"] == ref["n_significant"] - n_wrong


@pytest.mark.parametrize("value", ["both", "two-sided", "Greater", None])
def test_an_unknown_alternative_is_refused(value):
    message = "alternative must be 'two.sided', 'greater' or 'less'"
    spec = _small_two_group()
    with pytest.raises(ValueError, match=message):
        power(spec, n_sims=2, alternative=value)
    with pytest.raises(ValueError, match=message):
        power_curve(spec, subject_ns=[20], n_sims=2, alternative=value)


def _crossed(effect=0.08, seed=4):
    """A 10 x 8 crossed Gaussian design with one within factor and random intercepts."""
    return {
        "name": "dirmix", "seed": seed,
        "units": {"subject": {"n": 10}, "item": {"n": 8}},
        "factors": [{"name": "cond", "levels": ["a", "b"], "contrasts": {"cond": [-0.5, 0.5]},
                     "vary_within": ["subject", "item"]}],
        "fixed": {"intercept": 6, "coefficients": {"cond": effect}},
        "random": {"subject": {"intercept_sd": 0.12}, "item": {"intercept_sd": 0.08}},
        "response": {"family": "gaussian", "name": "y", "sigma": 0.3},
    }


def _mixed():
    pytest.importorskip("statsmodels")
    pytest.importorskip("pandas")
    from pilotr import power_mixed
    return power_mixed


def test_one_sided_power_mixed_at_alpha_is_the_two_sided_test_at_twice_alpha():
    power_mixed = _mixed()
    # An effect small enough that over 12 replicates each direction reaches significance once.
    spec = _crossed(effect=0.005)
    two = power_mixed(spec, n_sims=12, alpha=0.2)
    greater = power_mixed(spec, n_sims=12, alpha=0.1, alternative="greater")
    less = power_mixed(spec, n_sims=12, alpha=0.1, alternative="less")
    assert greater["n_significant"] > 0 and less["n_significant"] > 0
    assert greater["n_significant"] + less["n_significant"] == two["n_significant"]
    n_wrong = round(two["type_s"] * two["n_significant"])
    assert less["n_significant"] == n_wrong
    assert (greater["type_s"], less["type_s"]) == (0, 1)
    for key in ("n_returned", "n_converged", "n_singular", "n_warning", "mean_estimate"):
        assert greater[key] == two[key], key


def test_each_mixed_replicate_takes_its_one_sided_p_from_the_wald_z():
    _mixed()
    from pilotr.power import _power_mixed_replicate
    spec = _as_spec(_crossed(effect=0.3, seed=1))

    def at(seed, alternative):
        return _power_mixed_replicate(seed, spec, "cond", {"a": -0.5, "b": 0.5}, "y",
                                      "gaussian", 0.0, alternative)

    for seed in replicate_seeds(spec["seed"], 2):
        two = at(seed, "two.sided")
        toward, away = ("greater", "less") if two["est"] > 0 else ("less", "greater")
        # statsmodels gives the two-sided p as twice the normal tail beyond |z|, and the one-sided
        # p is that tail, so the two agree to the last bit.
        assert at(seed, toward)["p"] == two["p"] / 2
        assert at(seed, away)["p"] == pytest.approx(1 - two["p"] / 2, rel=1e-12)
        assert at(seed, toward)["est"] == two["est"]


def test_the_backend_names_a_one_sided_test():
    power_mixed = _mixed()
    spec = _crossed()
    assert power_mixed(spec, n_sims=1)["backend"] == (
        "statsmodels MixedLM (crossed variance components, REML; Wald z tests)")
    assert power_mixed(spec, n_sims=1, alternative="greater")["backend"] == (
        "statsmodels MixedLM (crossed variance components, REML; one-sided Wald z tests, "
        "greater)")
    with pytest.raises(ValueError, match="alternative must be 'two.sided', 'greater' or 'less'"):
        power_mixed(spec, n_sims=1, alternative="both")


def test_a_null_effect_leaves_type_s_and_type_m_undefined_in_either_direction():
    spec = _small_two_group()
    spec["fixed"]["coefficients"] = {}
    r = power(spec, n_sims=50, alternative="greater")
    assert math.isnan(r["type_s"]) and math.isnan(r["type_m"])
    assert 0 <= r["power"] <= 1
