"""Tests for solving a simulated design curve (run with: pytest python/tests).

The curves are written out rather than simulated, so that a change to the generative core
cannot move a number here, and so that the R suite can hold exactly the same tables in
test-solve-curve.R.

`analytic_curve` is the power of a two-sample t-test at a standardised effect of 0.5, taken
from R's `stats::power.t.test` at each total sample size and rounded to four places. It is the
external reference the solver is checked against: the same function reports that 0.80 power
needs 127.5315 subjects in total and 0.90 needs 170.0626, and the solve has to land near those
without ever having been told them.
"""
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pytest

from pilotr import solve_curve, target_n

ANALYTIC_N80 = 127.5315
ANALYTIC_N90 = 170.0626

# What the R engine returns for `solve_curve(analytic_curve(), target=0.8)`, recorded at full
# double precision. The two engines are not required to agree bit for bit here, because the fit
# calls exp() and log() and the languages do not share a maths library; they are required to
# agree far more closely than any reading of the result could depend on. The full cross-engine
# comparison, including every refusal message, is tools/parity/solve_cross.py.
R_REFERENCE = {"value": 127.57452392848577, "lo": 125.23529527081689,
               "hi": 129.9353977416252, "se": 0.053078331919020047}
PARITY_TOLERANCE = 1e-9


def analytic_curve(n_sims=2000):
    ns = [40, 60, 80, 100, 120, 140, 160, 180]
    pw = [0.3377, 0.4778, 0.5981, 0.6969, 0.7753, 0.8358, 0.8816, 0.9156]
    return [{"n_subject": n, "power": p, "n_sims": n_sims} for n, p in zip(ns, pw)]


def test_solved_size_matches_the_analytic_answer():
    s80 = solve_curve(analytic_curve(), target=0.8)
    s90 = solve_curve(analytic_curve(), target=0.9)
    assert abs(s80["value"] - ANALYTIC_N80) / ANALYTIC_N80 < 0.02
    assert abs(s90["value"] - ANALYTIC_N90) / ANALYTIC_N90 < 0.02
    # The interval is the point of the exercise, so it has to cover the value it is estimating.
    assert s80["lo"] <= ANALYTIC_N80 <= s80["hi"]
    assert s90["lo"] <= ANALYTIC_N90 <= s90["hi"]


def test_solve_agrees_with_the_r_engine():
    s = solve_curve(analytic_curve(), target=0.8)
    for k, expected in R_REFERENCE.items():
        assert abs(s[k] - expected) / abs(expected) < PARITY_TOLERANCE


def test_the_returned_solve_describes_itself():
    s = solve_curve(analytic_curve(), target=0.8)
    assert set(s) == {"value", "lo", "hi", "level", "target", "se", "dispersion", "x", "y",
                      "transform", "intercept", "slope", "n_points", "x_min", "x_max"}
    assert s["x"] == "n_subject"
    assert s["y"] == "power"
    assert s["transform"] == "sqrt"
    assert s["n_points"] == 8
    assert s["x_min"] == 40 and s["x_max"] == 180
    assert s["lo"] < s["value"] < s["hi"]
    assert s["slope"] > 0


def test_a_wider_level_gives_a_wider_interval():
    narrow = solve_curve(analytic_curve(), target=0.8, level=0.8)
    wide = solve_curve(analytic_curve(), target=0.8, level=0.99)
    assert narrow["value"] == wide["value"]
    assert wide["lo"] < narrow["lo"] and wide["hi"] > narrow["hi"]


def test_more_replicates_narrow_the_interval_but_leave_the_point_alone():
    few = solve_curve(analytic_curve(n_sims=50), target=0.8)
    many = solve_curve(analytic_curve(n_sims=5000), target=0.8)
    # A weight common to every point cancels out of the fitted line, though not out of the last
    # bits of it, so the two solves agree to rounding rather than exactly.
    assert abs(few["value"] - many["value"]) / many["value"] < 1e-12
    assert many["hi"] - many["lo"] < few["hi"] - few["lo"]


def test_target_n_rounds_up():
    s = target_n(analytic_curve())
    assert s["target"] == 0.8
    assert s["n"] == math.ceil(s["value"])
    assert s["n_lo"] == math.ceil(s["lo"])
    assert s["n_hi"] == math.ceil(s["hi"])
    assert 0 <= s["n"] - s["value"] < 1


def test_a_precision_curve_is_solved_through_its_own_rate_column():
    ns = [15, 30, 60, 100, 140, 180, 220, 260]
    pm = [0.06, 0.14, 0.33, 0.55, 0.71, 0.83, 0.90, 0.94]
    curve = [{"n_subject": n, "p_meaningful": p, "n_returned": 200} for n, p in zip(ns, pm)]
    s = solve_curve(curve, target=0.9)
    assert s["y"] == "p_meaningful"
    assert 180 < s["value"] < 260


def test_a_non_sample_size_axis_is_solved_on_the_identity_scale():
    curve = [{"effect_size": d, "power": p, "n_sims": 400} for d, p in
             zip([0.1, 0.2, 0.3, 0.4, 0.5], [0.11, 0.31, 0.58, 0.79, 0.92])]
    s = solve_curve(curve, target=0.8, transform="identity")
    assert s["x"] == "effect_size"
    assert 0.3 < s["value"] < 0.5


def multi_effect_curve():
    ns = [20, 20, 40, 40, 60, 60, 80, 80, 100, 100]
    names = ["a", "b"] * 5
    pw = [0.20, 0.10, 0.40, 0.18, 0.62, 0.27, 0.78, 0.35, 0.88, 0.44]
    return [{"n_subject": n, "effect": e, "power": p, "n_sims": 300}
            for n, e, p in zip(ns, names, pw)]


def test_a_focal_effect_is_selected_out_of_a_curve_holding_several():
    curve = multi_effect_curve()
    s = solve_curve(curve, target=0.7, effect="a")
    assert s["n_points"] == 5
    assert 40 < s["value"] < 100
    with pytest.raises(ValueError) as e:
        solve_curve(curve, target=0.7)
    assert str(e.value) == (
        "the curve has more than one row at the same swept value. Select a single focal "
        "effect with `effect`, or subset the curve before solving.")
    with pytest.raises(ValueError) as e:
        solve_curve(curve, target=0.7, effect="z")
    assert str(e.value) == "`curve` holds no focal effect named 'z'. It holds: a, b."


def test_the_replicate_count_may_be_given_rather_than_read_off_a_column():
    curve = [{"n_subject": n, "power": p} for n, p in
             zip([40, 80, 120, 160], [0.34, 0.60, 0.78, 0.88])]
    with pytest.raises(ValueError) as e:
        solve_curve(curve, target=0.7)
    assert str(e.value) == (
        "`curve` has no replicate-count column. Name one with `n`; the columns recognised "
        "automatically are 'n_returned', 'n_converged', 'n_sims'.")
    s = solve_curve(curve, target=0.7, n=500)
    assert 40 < s["value"] < 160


# --- refusals ---------------------------------------------------------------------------


@pytest.mark.parametrize("bad", [0, 1, 1.5, -0.2])
def test_a_target_outside_the_unit_interval_is_refused(bad):
    with pytest.raises(ValueError) as e:
        solve_curve(analytic_curve(), target=bad)
    assert str(e.value) == "`target` must be strictly between 0 and 1; got %g." % bad


def test_a_target_that_is_not_a_single_number_is_refused():
    for bad in (float("nan"), None, [0.5, 0.8]):
        with pytest.raises(ValueError) as e:
            solve_curve(analytic_curve(), target=bad)
        assert str(e.value) == "`target` must be a single number strictly between 0 and 1."
    with pytest.raises(ValueError) as e:
        solve_curve(analytic_curve(), target=0.8, level=1)
    assert str(e.value) == "`level` must be strictly between 0 and 1; got 1."


def test_a_curve_that_never_reaches_the_target_is_refused_with_its_range():
    with pytest.raises(ValueError) as e:
        solve_curve(analytic_curve(), target=0.95)
    assert str(e.value) == (
        "the curve does not reach a power of 0.95 within its swept range. Over swept values "
        "from 40 to 180 the rate runs from 0.3377 to 0.9156, and solving would extrapolate "
        "beyond the values simulated.")


def test_a_solve_landing_outside_the_swept_range_is_refused():
    # The target is exactly the lowest rate simulated, so the bracket test passes; the fitted
    # curve sits a shade above that point, so the solve falls just short of the smallest size
    # swept, which is outside the range whether it misses by one subject or a hundred.
    with pytest.raises(ValueError) as e:
        solve_curve(analytic_curve(), target=0.3377)
    assert str(e.value).endswith(
        "falls outside the swept range 40 to 180, so reporting it would extrapolate beyond "
        "the values simulated.")


def test_a_target_exactly_on_the_highest_rate_simulated_is_solved():
    s = solve_curve(analytic_curve(), target=0.9156)
    assert 40 <= s["value"] <= 180


def test_a_curve_the_model_does_not_describe_widens_its_interval():
    clean = solve_curve(analytic_curve(), target=0.8)
    assert clean["dispersion"] == 1
    # The same curve with its middle points pushed about, which no two-parameter model can pass
    # through. The point estimate barely moves; the interval has to say that it is less certain.
    rough = [dict(r) for r in analytic_curve()]
    for r, p in zip(rough, [0.3377, 0.5600, 0.5400, 0.7600, 0.7100, 0.8700, 0.8500, 0.9156]):
        r["power"] = p
    ragged = solve_curve(rough, target=0.8)
    assert ragged["dispersion"] > 1
    assert ragged["hi"] - ragged["lo"] > 3 * (clean["hi"] - clean["lo"])


def test_a_curve_with_too_few_points_is_refused():
    curve = analytic_curve()
    for rows, k in ((curve[:2], 2), (curve[:1], 1)):
        with pytest.raises(ValueError) as e:
            solve_curve(rows, target=0.5)
        assert str(e.value) == (
            "solving a curve needs at least 3 swept values with a rate and a replicate "
            "count; this curve has %d." % k)
    with pytest.raises(ValueError) as e:
        solve_curve([], target=0.5)
    assert str(e.value) == "`curve` must be a table of curve points, with one row per swept value."
    # A point with no converged fits has no rate, and is dropped before the count is taken.
    # None is the same absence written the other way: it is what R's NA becomes when a curve
    # crosses to this side as JSON, and the R twin drops such a point rather than refusing.
    for missing in (float("nan"), None):
        thin = [dict(r) for r in curve[:4]]
        thin[1]["power"] = missing
        thin[2]["power"] = missing
        with pytest.raises(ValueError) as e:
            solve_curve(thin, target=0.5)
        assert str(e.value) == (
            "solving a curve needs at least 3 swept values with a rate and a replicate count; "
            "this curve has 2.")


def test_a_curve_with_no_trend_is_refused():
    flat = [{"n_subject": n, "power": 0.5, "n_sims": 100} for n in (10, 20, 30)]
    with pytest.raises(ValueError) as e:
        solve_curve(flat, target=0.5)
    assert str(e.value) == (
        "the rate is 0.5 at every swept value, so the curve has no trend to invert.")


def test_a_curve_whose_slope_cannot_be_told_from_zero_is_refused():
    noise = [{"n_subject": n, "power": p, "n_sims": 30} for n, p in
             zip([10, 20, 30, 40, 50], [0.50, 0.52, 0.49, 0.51, 0.50])]
    with pytest.raises(ValueError) as e:
        solve_curve(noise, target=0.5)
    assert str(e.value).endswith(
        "cannot be told from zero at the 0.95 level, so the curve does not determine a "
        "value; every swept value is compatible with the target.")


def test_a_transform_undefined_on_the_swept_values_is_refused():
    curve = [{"delta": d, "power": p, "n_sims": 200} for d, p in
             zip([-0.2, 0.0, 0.2], [0.20, 0.50, 0.85])]
    for transform in ("sqrt", "log"):
        with pytest.raises(ValueError) as e:
            solve_curve(curve, target=0.6, transform=transform)
        assert str(e.value) == (
            "the '%s' transform is not defined for a swept value of -0.2. Use "
            'transform = "identity" for an axis that is not a sample size.' % transform)
    s = solve_curve(curve, target=0.6, transform="identity")
    assert -0.2 < s["value"] < 0.2


def test_a_column_that_is_not_a_rate_or_not_there_is_refused():
    with pytest.raises(ValueError) as e:
        solve_curve(analytic_curve(), target=0.8, y="nope")
    assert str(e.value) == (
        "`curve` has no column named 'nope'. Its columns are: n_subject, power, n_sims.")
    widths = [{"n": n, "width": w} for n, w in zip(range(1, 6), [0.9, 0.7, 0.5, 0.4, 0.3])]
    with pytest.raises(ValueError) as e:
        solve_curve(widths, target=0.5)
    assert str(e.value) == (
        "`curve` has no decision-rate column. Name one with `y`; the columns recognised "
        "automatically are 'power', 'p_meaningful'.")
    wide = [{"n": n, "power": p, "n_sims": 100} for n, p in
            zip([1, 2, 3, 4], [0.1, 0.5, 0.9, 1.4])]
    with pytest.raises(ValueError) as e:
        solve_curve(wide, target=0.5)
    assert str(e.value) == (
        "the column 'power' holds a value of 1.4, which is not a probability; solve_curve() "
        "inverts a rate between 0 and 1.")
    with pytest.raises(ValueError) as e:
        solve_curve(analytic_curve(), target=0.8, transform="cube")
    assert str(e.value) == "`transform` must be one of 'sqrt', 'identity', 'log'."


def test_solve_curve_consumes_what_power_curve_produces():
    scipy = pytest.importorskip("scipy")  # noqa: F841  power() needs it for the t-test
    from pilotr import power_curve
    spec = {
        "name": "s", "seed": 11,
        "units": {"subject": {"n": 40}},
        "factors": [{"name": "group", "levels": ["a", "b"],
                     "contrasts": {"effect": [-0.5, 0.5]}, "between": "subject"}],
        "fixed": {"intercept": 0, "coefficients": {"effect": 0.6}},
        "response": {"family": "gaussian", "name": "score", "sigma": 1},
    }
    curve = power_curve(spec, subject_ns=[20, 40, 60, 80, 100], n_sims=200)
    s = target_n(curve)
    assert s["x"] == "n_subject"
    assert s["y"] == "power"
    assert 20 <= s["n"] <= 100
