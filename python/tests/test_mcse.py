"""The Monte Carlo error of every rate the power functions report, held bit for bit against R.

`pilotr.mcse` performs the operations of the R package's mcse.R in the same order, so for the
same counts the twins return the same doubles. The pins below are those doubles as
`float.hex()` writes them. The R package's tests/testthat/test-mcse.R holds the same strings,
which the last test here checks whenever the R package sits beside this one.
"""
import copy, math, os, re, statistics, sys, warnings
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pytest

from pilotr import power, power_curve, power_mixed
from pilotr.mcse import _rate_with_error, _type_m_mcse, _type_s_interval

# count, n, then the rate, its standard error and its Wilson bounds.
RATE_PINS = [
    (0, 60, "0x0.0p+0", "0x0.0p+0", "0x0.0p+0", "0x1.eced8522e178fp-5"),
    (1, 60, "0x1.1111111111111p-6", "0x1.0ec813a58caffp-6", "0x1.826bd3cdd1810p-9",
     "0x1.6ab4c3dd18080p-4"),
    (8, 50, "0x1.47ae147ae147bp-3", "0x1.a8b8c95a80272p-5", "0x1.55803137393f1p-4",
     "0x1.23fc4e7a78f78p-2"),
    (16, 20, "0x1.999999999999ap-1", "0x1.6e5b7d16657e0p-4", "0x1.2affc35bd30ebp-1",
     "0x1.d6b409f2188edp-1"),
    (30, 60, "0x1.0000000000000p-1", "0x1.08654a2d4f6dap-4", "0x1.82681a193ffe0p-2",
     "0x1.3ecbf2f360010p-1"),
    (59, 60, "0x1.f777777777777p-1", "0x1.0ec813a58cb06p-6", "0x1.d2a967845cfefp-1",
     "0x1.fe7d942c322e7p-1"),
    (60, 60, "0x1.0000000000000p+0", "0x0.0p+0", "0x1.e13127add1e86p-1",
     "0x1.0000000000000p+0"),
]

# Seven exaggeration ratios, and the standard error and Type S interval the twins pin for them.
TYPE_M_RATIOS = [1.82, 2.41, 1.37, 3.05, 1.96, 2.73, 1.51]
TYPE_M_PIN = "0x1.e6a0d0061eecfp-3"
TYPE_S_PIN = ("0x1.6f49c628ba020p-6", "0x1.e23082fab5585p-2")  # 1 wrong sign in 8

_R_TESTS = os.path.join(os.path.dirname(__file__), "..", "..", "r", "pilotr", "tests", "testthat",
                        "test-mcse.R")

_POWER_KEYS = ("power", "power_mcse", "power_lo", "power_hi")


@pytest.mark.parametrize("count, n, p, mcse, lo, hi", RATE_PINS)
def test_each_rate_its_standard_error_and_its_wilson_bounds_match_the_pinned_doubles(
        count, n, p, mcse, lo, hi):
    got = _rate_with_error(count, n, "p")
    assert list(got) == ["p", "p_mcse", "p_lo", "p_hi"]
    assert [v.hex() for v in got.values()] == [p, mcse, lo, hi]


def test_a_rate_over_no_replicates_is_undefined_with_its_error():
    assert all(math.isnan(v) for v in _rate_with_error(0, 0, "power").values())


def test_type_m_mcse_is_the_standard_deviation_of_the_ratios_over_the_root_of_their_number():
    assert _type_m_mcse(TYPE_M_RATIOS).hex() == TYPE_M_PIN
    # statistics.stdev() works in exact fractions and rounds once, so it can differ from the
    # plain double fold in the last bits, and it agrees to within that.
    assert _type_m_mcse(TYPE_M_RATIOS) == pytest.approx(
        statistics.stdev(TYPE_M_RATIOS) / math.sqrt(len(TYPE_M_RATIOS)), rel=1e-14)
    # One ratio has no standard deviation.
    assert math.isnan(_type_m_mcse([1.4]))
    assert math.isnan(_type_m_mcse([]))


def test_the_type_s_interval_is_a_wilson_interval_withheld_without_a_true_effect():
    assert tuple(v.hex() for v in _type_s_interval(1, 8, 0.4)) == TYPE_S_PIN
    assert tuple(v.hex() for v in _type_s_interval(1, 8, -0.4)) == TYPE_S_PIN
    for args in ((0, 8, 0.0), (0, 8, math.nan), (0, 0, 0.4)):
        assert all(math.isnan(v) for v in _type_s_interval(*args))


def _two_group(effect=5.0):
    """A two-group Gaussian design with one row per subject."""
    return {
        "name": "mc", "seed": 2024,
        "units": {"subject": {"n": 64}},
        "factors": [{"name": "group", "levels": ["control", "treatment"],
                     "contrasts": {"effect": [-0.5, 0.5]}, "between": "subject"}],
        "fixed": {"intercept": 100, "coefficients": {"effect": effect}},
        "response": {"family": "gaussian", "name": "score", "sigma": 10},
    }


def _same(a, b):
    return a == b or (math.isnan(a) and math.isnan(b))


def test_power_reports_the_monte_carlo_error_of_power_type_s_and_type_m():
    pytest.importorskip("scipy")
    r = power(_two_group(), n_sims=100)
    assert list(r) == ["n_sims", "alpha", "power", "power_mcse", "power_lo", "power_hi",
                       "n_significant", "true_effect", "mean_estimate", "type_s", "type_s_lo",
                       "type_s_hi", "type_m", "type_m_mcse"]
    assert ({k: r[k] for k in _POWER_KEYS}
            == _rate_with_error(r["n_significant"], r["n_sims"], "power"))
    wrong = round(r["type_s"] * r["n_significant"])
    assert ((r["type_s_lo"], r["type_s_hi"])
            == _type_s_interval(wrong, r["n_significant"], r["true_effect"]))
    assert 0 < r["type_m_mcse"] < r["type_m"]


def test_type_m_needs_two_significant_replicates_for_its_error_and_type_s_an_effect():
    pytest.importorskip("scipy")
    # With alpha = 1, every replicate is significant, so n_sims sets the count exactly.
    one = power(_two_group(), n_sims=1, alpha=1)
    assert one["n_significant"] == 1
    assert not math.isnan(one["type_m"]) and math.isnan(one["type_m_mcse"])
    assert not (math.isnan(one["type_s_lo"]) or math.isnan(one["type_s_hi"]))
    two = power(_two_group(), n_sims=2, alpha=1)
    assert two["n_significant"] == 2 and math.isfinite(two["type_m_mcse"])
    # Type S and Type M are undefined at a zero effect, and so is their uncertainty. Power keeps
    # its own, since it is defined there.
    null = power(_two_group(effect=0.0), n_sims=50)
    assert all(math.isnan(null[k])
               for k in ("type_s", "type_s_lo", "type_s_hi", "type_m", "type_m_mcse"))
    assert not any(math.isnan(null[k]) for k in _POWER_KEYS)


def test_each_power_curve_record_carries_the_errors_power_reports():
    pytest.importorskip("scipy")
    spec = _two_group()
    curve = power_curve(spec, subject_ns=[32, 64], n_sims=50)
    for rec in curve:
        assert list(rec) == ["n_subject", "power", "power_mcse", "power_lo", "power_hi", "type_s",
                             "type_s_lo", "type_s_hi", "type_m", "type_m_mcse", "n_sims",
                             "n_significant"]
        s = copy.deepcopy(spec)
        s["units"]["subject"]["n"] = rec["n_subject"]
        r = power(s, n_sims=50)
        assert all(_same(rec[k], r[k]) for k in list(rec)[1:])


def _crossed(family="gaussian", effect=0.3, **response):
    """A 10 x 8 crossed design with one within factor and random intercepts."""
    return {
        "name": "mc", "seed": 1,
        "units": {"subject": {"n": 10}, "item": {"n": 8}},
        "factors": [{"name": "cond", "levels": ["a", "b"], "contrasts": {"cond": [-0.5, 0.5]},
                     "vary_within": ["subject", "item"]}],
        "fixed": {"intercept": 6, "coefficients": {"cond": effect}},
        "random": {"subject": {"intercept_sd": 0.12}, "item": {"intercept_sd": 0.08}},
        "response": dict({"family": family, "name": "y"}, **response),
    }


def test_power_mixed_reports_the_monte_carlo_error_of_power_type_s_and_type_m():
    pytest.importorskip("statsmodels")
    pytest.importorskip("pandas")
    r = power_mixed(_crossed(sigma=0.3), n_sims=6)
    assert r["n_significant"] >= 2
    # power is taken over the replicates that returned a fit, and so is its error.
    assert ({k: r[k] for k in _POWER_KEYS}
            == _rate_with_error(r["n_significant"], r["n_returned"], "power"))
    wrong = round(r["type_s"] * r["n_significant"])
    assert ((r["type_s_lo"], r["type_s_hi"])
            == _type_s_interval(wrong, r["n_significant"], r["true_effect"]))
    assert r["type_m_mcse"] > 0

    null = power_mixed(_crossed(effect=0.0, sigma=0.3), n_sims=3)
    assert all(math.isnan(null[k])
               for k in ("type_s", "type_s_lo", "type_s_hi", "type_m", "type_m_mcse"))


def test_power_mixed_withholds_the_error_of_a_type_m_it_withholds():
    pytest.importorskip("statsmodels")
    pytest.importorskip("pandas")
    # A Bernoulli response is fitted on the response scale while its effect is on the logit
    # scale, so Type M is withheld, and with it its standard error. Type S compares signs, which
    # the link keeps, so its interval stays. With alpha = 1, both replicates are significant.
    spec = _crossed("bernoulli", effect=0.8)
    spec["fixed"]["intercept"] = 0.5
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        r = power_mixed(spec, n_sims=2, alpha=1)
    assert r["n_significant"] == 2
    assert math.isnan(r["type_m"]) and math.isnan(r["type_m_mcse"])
    assert not (math.isnan(r["type_s_lo"]) or math.isnan(r["type_s_hi"]))


def test_the_r_twin_pins_the_same_doubles():
    if not os.path.exists(_R_TESTS):
        pytest.skip("the R package's tests are not beside this checkout")
    with open(_R_TESTS, encoding="utf-8") as f:
        text = f.read()
    rows = re.findall(r"^\s*(\d+)\s+(\d+)\s+(0x\S+)\s+(0x\S+)\s+(0x\S+)\s+(0x\S+)\s*$", text,
                      flags=re.MULTILINE)
    assert [(int(count), int(n), *pins) for count, n, *pins in rows] == RATE_PINS
    ratios = re.search(r"type_m_ratios <- c\(([^)]*)\)", text).group(1)
    assert [float(v) for v in ratios.split(",")] == TYPE_M_RATIOS
    assert re.search(r'type_m_pin <- "([^"]+)"', text).group(1) == TYPE_M_PIN
    assert re.search(r'type_s_pin <- c\("([^"]+)", "([^"]+)"\)', text).groups() == TYPE_S_PIN
