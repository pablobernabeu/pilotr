"""Minimal validation suite (run with: pytest python/tests).

Covers the three things a methods-package reviewer checks first: the RNG is
deterministic, the inverse-normal is numerically correct, and the engine recovers the
ground-truth parameters it was given.
"""
import os, re, sys, statistics, math
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pytest

from pilotr import RNG, as241, replicate_seeds, simulate, load_spec, validate_spec
from pilotr.core import gamma_mt

SPEC = os.path.join(os.path.dirname(__file__), "..", "..", "spec", "examples")


def test_rng_is_deterministic():
    a = [RNG(123).uniform() for _ in range(1)]
    b = RNG(123).uniform()
    assert a[0] == b
    stream = RNG(7)
    xs = [stream.uniform() for _ in range(1000)]
    assert all(0.0 < x < 1.0 for x in xs)
    assert len(set(xs)) == len(xs)  # no immediate repeats


def test_as241_matches_known_quantiles():
    # qnorm(0.975) = 1.959963985; qnorm(0.5) = 0; symmetric
    assert abs(as241(0.975) - 1.959963984540054) < 1e-9
    assert abs(as241(0.5)) < 1e-12
    assert abs(as241(0.975) + as241(0.025)) < 1e-9


def test_gaussian_recovers_effect_at_large_n():
    s = load_spec(os.path.join(SPEC, "between_2group_gaussian.json"))
    s["units"]["subject"]["n"] = 40000
    d = simulate(s)
    g0 = [r["score"] for r in d.rows if r["group"] == "control"]
    g1 = [r["score"] for r in d.rows if r["group"] == "treatment"]
    assert len(g0) == len(g1) == 20000
    assert abs((statistics.mean(g1) - statistics.mean(g0)) - 5.0) < 0.3


def test_mixed_rt_has_expected_structure():
    d = simulate(os.path.join(SPEC, "crossed_mixed_rt.json"))
    assert len(d) == 30 * 24 * 2
    rel = [r["RT"] for r in d.rows if r["condition"] == "related"]
    unr = [r["RT"] for r in d.rows if r["condition"] == "unrelated"]
    assert min(rel) > 200  # shifted-lognormal floor = shift (200 ms non-decision time)
    assert statistics.mean(unr) > statistics.mean(rel)  # positive priming effect


def test_simulation_is_reproducible():
    d1 = simulate(os.path.join(SPEC, "crossed_mixed_rt.json"))
    d2 = simulate(os.path.join(SPEC, "crossed_mixed_rt.json"))
    assert [r["RT"] for r in d1.rows] == [r["RT"] for r in d2.rows]


def test_beta_family_in_unit_interval_and_ordered():
    d = simulate(os.path.join(SPEC, "beta_proportion.json"))
    vals = [r["prop"] for r in d.rows]
    assert all(0.0 < v < 1.0 for v in vals)
    ctrl = [r["prop"] for r in d.rows if r["group"] == "control"]
    trt = [r["prop"] for r in d.rows if r["group"] == "treatment"]
    assert statistics.mean(trt) > statistics.mean(ctrl)  # grp +0.8 on the logit scale


def test_partial_crossing_subset_size():
    d = simulate(os.path.join(SPEC, "partial_crossing.json"))
    items_by_subject = {}
    for r in d.rows:
        items_by_subject.setdefault(r["subject"], set()).add(r["item"])
    assert len(items_by_subject) == 60
    assert all(len(v) == 12 for v in items_by_subject.values())  # per_subject = 12


def test_per_subject_must_lie_between_1_and_n_items():
    # The same inputs must raise in the R twin, preserving cross-language parity.
    spec = {
        "name": "pc", "seed": 1,
        "units": {"subject": {"n": 2}, "item": {"n": 3, "per_subject": 5}},
        "factors": [{"name": "cond", "levels": ["a", "b"],
                     "contrasts": {"cond": [-0.5, 0.5]}, "vary_within": ["subject", "item"]}],
        "fixed": {"intercept": 6, "coefficients": {"cond": 0.05}},
        "random": {"subject": {"intercept_sd": 0.1}, "item": {"intercept_sd": 0.1}},
        "response": {"family": "gaussian", "name": "y", "sigma": 0.3},
    }
    with pytest.raises(ValueError, match="cannot exceed the number of items"):
        simulate(spec)
    spec["units"]["item"]["per_subject"] = 0
    with pytest.raises(ValueError, match="at least 1"):
        simulate(spec)


# The refusal is byte-identical to the R twin's; Python's "inf" is normalised to R's "Inf" when
# the linear predictor overflows exp().
def _poisson_refusal(mean):
    return ("the poisson mean exp(eta) = %s is too large for the inverse-CDF sampler: beyond a "
            "mean of about 708.4, exp(-mean) is subnormal and the walk loses accuracy. Lower the "
            "poisson intercept or coefficients until the implied mean is at most 708." % mean)


def test_a_poisson_mean_beyond_the_samplers_accurate_range_is_refused():
    def pois(intercept, n=4):
        return {
            "name": "p", "seed": 1,
            "units": {"subject": {"n": n}},
            "factors": [{"name": "group", "levels": ["a", "b"],
                         "contrasts": {"grp": [-0.5, 0.5]}, "between": "subject"}],
            "fixed": {"intercept": intercept, "coefficients": {"grp": 0.0}},
            "random": {},
            "response": {"family": "poisson", "name": "count"},
        }

    # From a mean of about 738, the walk's total mass departs from 1. Over 40,000 draws, a mean
    # of 742 was refused for 4% of them. A mean of 744 drew counts with a mean of 733.5 and a
    # variance of 407.
    with pytest.raises(ValueError, match=re.escape(_poisson_refusal("742"))):
        simulate(pois(math.log(742)))
    # An intercept of 7 implies a mean of exp(7), about 1097, where exp(-mean) underflows to
    # zero; every count then used to come back as the sampler's 1e6 iteration cap.
    with pytest.raises(ValueError, match=re.escape(_poisson_refusal("1096.63"))):
        simulate(pois(7))
    with pytest.raises(ValueError, match=re.escape(_poisson_refusal("Inf"))):
        simulate(pois(800))  # overflows math.exp itself; R's exp() returns Inf
    # A mean of 700 is served, and its counts centre on it.
    counts = [r["count"] for r in simulate(pois(math.log(700), n=200)).rows]
    assert all(isinstance(c, int) for c in counts)
    assert abs(statistics.fmean(counts) - 700) < 5 * math.sqrt(700 / 200)
    # Feasible means are untouched: the shipped poisson example still draws real counts.
    d = simulate(os.path.join(SPEC, "poisson_counts_between.json"))
    assert all(0 <= r["count"] < 1_000_000 for r in d.rows)


def test_the_poisson_cap_is_the_double_nearest_minus_log_of_the_smallest_normal_double():
    from pilotr.core import poisson_inv

    # 1022 * log(2) rounds to the literal both twins compare against. exp(-literal) is still a
    # normal double, so that mean is served. The next double up gives a subnormal and is refused.
    assert poisson_inv(708.3964185322641, 0.5) == 708
    with pytest.raises(ValueError, match=re.escape(_poisson_refusal("708.396"))):
        poisson_inv(float.fromhex("0x1.6232bdd7abcd3p+9"), 0.5)


def test_the_poisson_walk_evaluates_each_term_as_the_r_twin_does():
    from pilotr.core import poisson_inv

    # From exp(-3.7), the cumulative sum at k = 3 is 0.4941532441504185 when each term is
    # (p * lam) / k, as in R, and 0.49415324415041845 when it is p * (lam / k), as here before.
    # This uniform therefore drew 3 in R and 4 in Python. The case was derived from the
    # correctly rounded exp(-3.7).
    if math.exp(-3.7) != float.fromhex("0x1.9511fc6871044p-6"):
        pytest.skip("this platform's exp(-3.7) is not the correctly rounded value")
    assert poisson_inv(3.7, 0.4941532441504185) == 3


# The edges of the beta and lognormal families. Each count below is the one the R twin pins, and
# each refusal is byte-identical to the twin's.
_BELOW_ONE = 1.0 - sys.float_info.epsilon / 2          # 1 - 2^-53, the largest double below 1
_SMALLEST = sys.float_info.min * sys.float_info.epsilon  # 2^-1074, the smallest positive double


def _beta_spec(intercept, phi, n=4000, seed=11, grp=0.5, rnd=None):
    spec = {
        "name": "b", "seed": seed,
        "units": {"subject": {"n": n}},
        "factors": [{"name": "group", "levels": ["a", "b"],
                     "contrasts": {"grp": [-0.5, 0.5]}, "between": "subject"}],
        "fixed": {"intercept": intercept, "coefficients": {"grp": grp}},
        "response": {"family": "beta", "name": "p", "phi": phi},
    }
    if rnd is not None:
        spec["response"]["round"] = rnd
    return spec


def _beta_shape_refusal(eta, phi):
    return ("the beta mean invlogit(eta) at a linear predictor of %s and precision phi = %s leaves "
            "a shape parameter of zero, so no proportion can be drawn. Keep eta between about -745 "
            "and 36.7 and phi large enough that mu * phi and (1 - mu) * phi stay positive."
            % (eta, phi))


def _beta_underflow_refusal(eta, phi):
    return ("both Gamma variates underflowed to zero at a linear predictor of %s and precision "
            "phi = %s, so the proportion is undefined. Raise phi." % (eta, phi))


def _exp_overflow_refusal(family, eta):
    return ("the %s response exp(eta + sigma * z) overflows at a linear predictor of %s, so no "
            "response can be drawn. Lower the intercept or the coefficients." % (family, eta))


def test_a_beta_response_never_equals_1_and_moving_it_inside_draws_nothing():
    # At an intercept of 3 and phi 2, (1 - mu) * phi is about 0.1. The ratio x / (x + y) then
    # rounded to exactly 1 in 161 of these 4,000 rows, a value brms's Beta() refuses. Each becomes
    # 1 - 2^-53, the largest double below 1, where 9 draws already lay.
    p = simulate(_beta_spec(3, 2)).column("p")
    assert not any(v in (0.0, 1.0) for v in p)
    assert p.count(_BELOW_ONE) == 170
    # The move draws nothing, so every later row is as it was. These are the sum and the last
    # value recorded before the change: the 161 moves of 2^-53 lie far inside the tolerance, and a
    # stream shifted by one draw far outside it.
    assert math.isclose(math.fsum(p), 3811.6303419786809, rel_tol=1e-12)
    assert math.isclose(p[-1], 0.99719118581120192, rel_tol=1e-12)


def test_a_beta_response_never_equals_0():
    # At an intercept of -7, mu * phi is about 0.002, and the first Gamma variate underflowed to
    # zero in 97 of these 400 rows, each of which was then exactly 0. Each is now 2^-1074, the
    # smallest positive double.
    p = simulate(_beta_spec(-7, 2, n=400, grp=0.0)).column("p")
    assert not any(v in (0.0, 1.0) for v in p)
    assert p.count(_SMALLEST) == 97


def test_a_rounded_beta_response_is_moved_inside_the_unit_interval_at_its_decimals():
    # A low mean rounded to three decimals: 694 of these 4,000 values, 17% of the rows, rounded to
    # exactly 0. Each becomes 0.001, where 218 values already were.
    p = simulate(_beta_spec(-2.5, 4, seed=909, rnd=3)).column("p")
    assert not any(v in (0.0, 1.0) for v in p)
    assert min(p) == 0.001
    assert p.count(0.001) == 912
    # A high mean rounded to three decimals: 2,088 values rounded to exactly 1 and become 0.999.
    p = simulate(_beta_spec(3, 2, rnd=3)).column("p")
    assert not any(v in (0.0, 1.0) for v in p)
    assert max(p) == 0.999
    assert p.count(0.999) == 2306


def test_the_bounds_of_a_rounded_beta_response_follow_the_specifications_rule():
    from pilotr.core import beta_round

    def bounds(d):
        return beta_round(0.0, d), beta_round(1.0, d)

    assert bounds(1) == (0.1, 0.9)
    assert bounds(3) == (0.001, 0.999)
    assert bounds(15) == (1e-15, 1 - 1e-15)
    # 1 - 1e-16 rounds to 1 - 2^-53, and from 17 places 1 - lo rounds to 1 itself, so the upper
    # bound stays at the largest double below 1.
    assert bounds(16) == (1e-16, _BELOW_ONE)
    assert bounds(17) == (1e-17, _BELOW_ONE)
    assert bounds(22) == (1e-22, _BELOW_ONE)
    # 10^23 is not a double, so the lower bound falls back to the smallest positive one.
    assert bounds(23) == (_SMALLEST, _BELOW_ONE)


def test_round_0_is_refused_for_the_beta_family():
    with pytest.raises(ValueError, match=re.escape(
            "'response.round' must be at least 1 for the beta family, because a proportion "
            "rounded to whole numbers is 0 or 1, outside the family's support")):
        validate_spec(_beta_spec(0, 8, rnd=0))
    validate_spec(_beta_spec(0, 8, rnd=1))


def test_a_beta_response_with_a_shape_of_zero_or_two_gamma_variates_of_zero_is_refused():
    # invlogit(40) rounds to 1, so (1 - mu) * phi is 0: R returned 1 in every row, and Python
    # raised ZeroDivisionError. invlogit(-750) underflows to 0, the mirror case, and at -700 with
    # phi 1e-30 mu * phi underflows to 0 although mu does not.
    with pytest.raises(ValueError, match=re.escape(_beta_shape_refusal("40", "2"))):
        simulate(_beta_spec(40, 2, n=8, grp=0.0))
    with pytest.raises(ValueError, match=re.escape(_beta_shape_refusal("-750", "2"))):
        simulate(_beta_spec(-750, 2, n=8, grp=0.0))
    with pytest.raises(ValueError, match=re.escape(_beta_shape_refusal("-700", "1e-30"))):
        simulate(_beta_spec(-700, 1e-30, n=8, grp=0.0))
    # At phi 0.001 both shapes are 0.0005, and both Gamma variates of a row can underflow to zero.
    # R returned NaN for 96 of 200 rows, and Python raised ZeroDivisionError.
    with pytest.raises(ValueError, match=re.escape(_beta_underflow_refusal("0", "0.001"))):
        simulate(_beta_spec(0, 0.001, n=200, grp=0.0))
    # Within the bounds the message gives, a proportion is drawn.
    p = simulate(_beta_spec(36.7, 2, n=8, grp=0.0)).column("p")
    assert all(0.0 < v < 1.0 for v in p)


class _FixedDraws:
    """Normals of 0 and one fixed uniform. gamma_mt() then accepts its first proposal and returns
    shape - 1/3 for a shape of at least 1. For a shape below 1, it boosts the variate for
    shape + 1 with the same uniform."""

    def __init__(self, u):
        self.u = u

    def normal(self):
        return 0.0

    def uniform(self):
        return self.u


def test_a_gamma_shape_of_one_half_is_boosted_by_u_times_u_as_in_the_r_twin():
    # R's ^ computes a square as u * u, where ** calls pow(). For this uniform, the 860th that
    # RNG(1) draws, glibc's pow() returned the double below u * u, so a Beta design with eta = 0
    # and phi = 1, whose Gamma shapes are both 1/2, drew different data in the twins.
    u = 0.6979060118654794
    assert u * u == float.fromhex("0x1.f2c336631a082p-2")
    g = gamma_mt(_FixedDraws(u), 1.5)
    assert gamma_mt(_FixedDraws(u), 0.5) == g * (u * u)
    # Any other exponent goes through pow(), as in R before 4.6.0.
    for shape in (1 / 3, 0.25, 0.4):
        g = gamma_mt(_FixedDraws(u), shape + 1.0)
        assert gamma_mt(_FixedDraws(u), shape) == g * u ** (1.0 / shape)
    # The boost's uniform is still drawn after the Gamma(3/2) variate.
    a, b = RNG(2026), RNG(2026)
    for _ in range(2000):
        x = gamma_mt(a, 0.5)
        g = gamma_mt(b, 1.5)
        v = b.uniform()
        assert x == g * (v * v)


def test_a_lognormal_response_that_overflows_is_refused():
    def rt(family):
        spec = _beta_spec(750, 1, n=8, grp=0.0)
        spec["response"] = {"family": family, "name": "RT", "sigma": 0.3}
        if family == "shifted_lognormal":
            spec["response"]["shift"] = 200
        return spec

    # exp(750) overflows, so R returned Inf in every row, and Python raised a bare OverflowError.
    with pytest.raises(ValueError, match=re.escape(_exp_overflow_refusal("lognormal", "750"))):
        simulate(rt("lognormal"))
    with pytest.raises(ValueError,
                       match=re.escape(_exp_overflow_refusal("shifted_lognormal", "750"))):
        simulate(rt("shifted_lognormal"))


def test_the_shipped_beta_example_simulates_the_data_it_always_has():
    # Moving a value inside (0, 1) touches only a value at 0 or 1, and this example has none.
    p = simulate(os.path.join(SPEC, "beta_proportion.json")).column("prop")
    assert len(p) == 200
    assert p[:3] == pytest.approx([0.29928, 0.54044, 0.25637])
    assert p[-3:] == pytest.approx([0.60073, 0.91115, 0.43263])
    assert math.fsum(p) == pytest.approx(102.45347)


# One specification file is run through both engines, so the same mistake has to be reported the
# same way by both. Each message below is byte-identical to the R twin's.
def test_a_whole_version_is_read_the_same_however_it_was_written():
    # R renders the JSON number 1.0 as "1", Python as "1.0". Neither engine may call it malformed.
    for declared in (1.0, 1, "1"):
        s = load_spec(os.path.join(SPEC, "between_2group_gaussian.json"))
        s["spec_version"] = declared
        with pytest.raises(ValueError,
                           match=re.escape("declares spec_version 1.0, which is newer")):
            validate_spec(s)


def test_a_non_object_unit_is_reported_rather_than_crashing():
    s = load_spec(os.path.join(SPEC, "between_2group_gaussian.json"))
    s["units"]["subject"] = 5
    with pytest.raises(ValueError, match=re.escape("'units.subject' must be an object")):
        validate_spec(s)
    # An empty object is a missing n, not a wrong shape, which is what the twin says too.
    s["units"]["subject"] = {}
    with pytest.raises(ValueError, match=re.escape("'units.subject.n' must be a whole number")):
        validate_spec(s)


def test_a_non_whole_seed_truncates_as_in_the_r_twin():
    # Reachable only outside validation, the fast path the replicate loops use.
    assert RNG(2.7).uniform() == RNG(2).uniform()
    assert RNG(3.5).uniform() == RNG(3).uniform()
    assert RNG(-2.7).uniform() == RNG(2).uniform()


def test_replicate_seeds_are_whole_numbers():
    # The annotation says list[int]; math.floor() + 1 has to keep it true.
    assert all(isinstance(s, int) for s in replicate_seeds(90210, 5))


def test_zero_true_effect_leaves_type_s_and_type_m_undefined():
    # A null condition (every effect zero) is a recommended input, not a hypothetical one, and
    # dividing by it used to raise ZeroDivisionError. The R twin returns NaN for the same spec.
    pytest.importorskip("scipy")
    from pilotr import power

    s = load_spec(os.path.join(SPEC, "between_2group_gaussian.json"))
    s["fixed"]["coefficients"]["grp"] = 0.0
    r = power(s, n_sims=100)
    assert r["true_effect"] == 0
    assert math.isnan(r["type_s"])
    assert math.isnan(r["type_m"])
    assert 0.0 <= r["power"] <= 1.0  # power itself is still reported


def test_power_mixed_also_refuses_to_divide_by_a_zero_true_effect():
    pytest.importorskip("statsmodels")
    pytest.importorskip("pandas")
    from pilotr import power_mixed

    spec = {
        "name": "null", "seed": 3,
        "units": {"subject": {"n": 10}, "item": {"n": 6}},
        "factors": [{"name": "cond", "levels": ["a", "b"],
                     "contrasts": {"cond": [-0.5, 0.5]}, "vary_within": ["subject", "item"]}],
        "fixed": {"intercept": 6, "coefficients": {"cond": 0.0}},
        "random": {"subject": {"intercept_sd": 0.12}, "item": {"intercept_sd": 0.08}},
        "response": {"family": "gaussian", "name": "y", "sigma": 0.3},
    }
    run = power_mixed(spec, n_sims=4)
    r = run["effects"]["cond"]
    assert run["n_returned"] > 0  # so the NaNs below are the guard, not a failed fit
    assert math.isnan(r["type_s"])
    assert math.isnan(r["type_m"])
    # Power is the significant share of the replicates that returned a fit, as in the R twin.
    assert r["n_significant"] == round(r["power"] * run["n_returned"])


def test_additional_grouping_column_present():
    d = simulate(os.path.join(SPEC, "nested_clusters.json"))
    assert "site" in d.columns
    sites = set(r["site"] for r in d.rows)
    assert sites == set(range(1, 13))  # 12 clusters
