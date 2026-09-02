"""Minimal validation suite (run with: pytest python/tests).

Covers the three things a methods-package reviewer checks first: the RNG is
deterministic, the inverse-normal is numerically correct, and the engine recovers the
ground-truth parameters it was given.
"""
import json, os, re, sys, statistics, math, warnings
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pytest

from pilotr import RNG, as241, replicate_seeds, simulate, load_spec, validate_spec

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


def test_a_poisson_mean_past_the_sampler_is_an_error_not_the_iteration_cap():
    # The refusal is byte-identical to the R twin's; Python's "inf" is normalised to R's
    # "Inf" when the linear predictor overflows exp(). An intercept of 7 implies a mean of
    # exp(7), about 1097, where exp(-mean) underflows to zero; every count then used to
    # come back as the sampler's 1e6 iteration cap.
    spec = {
        "name": "p", "seed": 1,
        "units": {"subject": {"n": 4}},
        "factors": [{"name": "group", "levels": ["a", "b"],
                     "contrasts": {"grp": [-0.5, 0.5]}, "between": "subject"}],
        "fixed": {"intercept": 7, "coefficients": {"grp": 0.0}},
        "random": {},
        "response": {"family": "poisson", "name": "count"},
    }
    with pytest.raises(ValueError, match=re.escape(
            "the poisson mean exp(eta) = 1096.63 is too large for the inverse-CDF "
            "sampler: the cumulative distribution cannot reach the drawn uniform, so no "
            "count can be drawn. Lower the poisson intercept or coefficients until the "
            "implied mean is simulable.")):
        simulate(spec)
    spec["fixed"]["intercept"] = 800  # overflows math.exp itself; R's exp() returns Inf
    with pytest.raises(ValueError, match=re.escape("exp(eta) = Inf")):
        simulate(spec)
    # Feasible means are untouched: the shipped poisson example still draws real counts.
    d = simulate(os.path.join(SPEC, "poisson_counts_between.json"))
    assert all(0 <= r["count"] < 1_000_000 for r in d.rows)


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


def test_a_name_that_would_take_another_column_s_place_is_refused():
    # A response named after the factor appended a second column of the same name here and
    # overwrote the factor column in the R twin, so one specification exported two different
    # tables.
    s = load_spec(os.path.join(SPEC, "between_2group_gaussian.json"))
    s["response"]["name"] = "group"
    with pytest.raises(ValueError, match=re.escape(
            "the name 'group' is used by factors[1].name and response.name")):
        validate_spec(s)
    s = load_spec(os.path.join(SPEC, "between_2group_gaussian.json"))
    s["factors"][0]["name"] = "subject"
    with pytest.raises(ValueError, match=re.escape(
            "the name 'subject' is used by units.subject and factors[1].name")):
        validate_spec(s)
    s = load_spec(os.path.join(SPEC, "nested_clusters.json"))
    s["response"]["name"] = "site"
    with pytest.raises(ValueError, match=re.escape(
            "the name 'site' is used by random.site and response.name")):
        validate_spec(s)


def test_a_blank_name_is_refused_where_it_used_to_reach_the_data():
    # An emptied 'Factor name' box built a spec that validated and then produced a column with no
    # name at all; the R twin stopped inside its simulator with "replacement has length zero".
    s = load_spec(os.path.join(SPEC, "between_2group_gaussian.json"))
    for blank in ("", "  "):
        s["factors"][0]["name"] = blank
        with pytest.raises(ValueError,
                           match=re.escape("factors[1].name must be a non-empty string")):
            validate_spec(s)
    s = load_spec(os.path.join(SPEC, "reading_time_continuous.json"))
    s["predictors"][0]["name"] = ""
    with pytest.raises(ValueError,
                       match=re.escape("predictors[1].name must be a non-empty string")):
        validate_spec(s)
    s = load_spec(os.path.join(SPEC, "between_2group_gaussian.json"))
    s["response"]["name"] = " "
    with pytest.raises(ValueError,
                       match=re.escape("'response.name' must be a non-empty string")):
        validate_spec(s)
    s = load_spec(os.path.join(SPEC, "nested_clusters.json"))
    s["random"][""] = s["random"].pop("site")
    with pytest.raises(ValueError, match=re.escape(
            "a 'random' grouping factor must have a non-empty name")):
        validate_spec(s)


def test_only_false_skips_validation_and_other_values_are_lenient(tmp_path):
    # `validate` is not a plain flag: False skips validation, True is strict and anything else
    # validates leniently, which is what both twins' docstrings promise.
    s = load_spec(os.path.join(SPEC, "between_2group_gaussian.json"))
    s["extra_field"] = "x"
    path = tmp_path / "spec.json"
    path.write_text(json.dumps(s))
    with pytest.raises(ValueError, match="unknown top-level field 'extra_field'"):
        load_spec(str(path))
    with pytest.warns(UserWarning, match="unknown top-level field 'extra_field'"):
        load_spec(str(path), validate=None)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        load_spec(str(path), validate=False)
    with pytest.warns(UserWarning, match="unknown top-level field 'extra_field'"):
        simulate(s, validate=None)


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


def test_power_refuses_a_between_level_with_fewer_than_two_units():
    # scipy answered one subject per group with a nan p-value and a power of 0, where R's
    # stats::t.test() stopped; both languages now refuse the design in the same words.
    pytest.importorskip("scipy")
    from pilotr import power

    s = load_spec(os.path.join(SPEC, "between_2group_gaussian.json"))
    for n in (2, 3):
        s["units"]["subject"]["n"] = n
        with pytest.raises(ValueError, match="at least 2 subjects at each level"):
            power(s, n_sims=5)
    s["units"]["subject"]["n"] = 4
    assert 0.0 <= power(s, n_sims=5)["power"] <= 1.0


def test_power_mixed_also_refuses_to_divide_by_a_zero_true_effect():
    pytest.importorskip("statsmodels")
    pytest.importorskip("pandas")
    from pilotr import power_mixed

    spec = {
        "name": "null", "seed": 3,
        "units": {"subject": {"n": 10}, "item": {"n": 6}},
        "factors": [{"name": "cond", "levels": ["a", "b"],
                     "contrasts": {"cond": [-0.5, 0.5]}, "vary_within": "subject"}],
        "fixed": {"intercept": 6, "coefficients": {"cond": 0.0}},
        "random": {"subject": {"intercept_sd": 0.12}, "item": {"intercept_sd": 0.08}},
        "response": {"family": "gaussian", "name": "y", "sigma": 0.3},
    }
    r = power_mixed(spec, n_sims=4)
    assert r["n_converged"] > 0  # so the NaNs below are the guard, not a failed fit
    assert math.isnan(r["type_s"])
    assert math.isnan(r["type_m"])
    # This n_converged counts the replicates that were fit, which the R twin calls n_returned,
    # so it is the denominator of power here.
    assert r["n_significant"] == round(r["power"] * r["n_converged"])


def test_additional_grouping_column_present():
    d = simulate(os.path.join(SPEC, "nested_clusters.json"))
    assert "site" in d.columns
    sites = set(r["site"] for r in d.rows)
    assert sites == set(range(1, 13))  # 12 clusters
