"""Minimal validation suite (run with: pytest toolkit/python/tests).

Covers the three things a methods-package reviewer checks first: the RNG is
deterministic, the inverse-normal is numerically correct, and the engine recovers the
ground-truth parameters it was given.
"""
import os, sys, statistics, math
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from simdgp import RNG, as241, simulate, load_spec

SPEC = os.path.join(os.path.dirname(__file__), "..", "..", "spec", "examples")


def test_rng_is_deterministic():
    a = [RNG(123).uniform() for _ in range(1)]
    b = RNG(123).uniform()
    assert a[0] == b
    stream = RNG(7); xs = [stream.uniform() for _ in range(1000)]
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
