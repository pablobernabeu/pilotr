"""Parallel execution must reproduce the serial stream exactly. Every replicate seeds
the shared cross-language RNG from its own index, so the worker count can change only
the wall-clock time, never a number.
"""
import math, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pytest


def assert_identical(a, b):
    """Assert two result dicts are identical, treating NaN as equal to NaN (the Type S and
    Type M entries are NaN when no replicate reaches significance)."""
    assert a.keys() == b.keys()
    for k in a:
        va, vb = a[k], b[k]
        if isinstance(va, float) and isinstance(vb, float) and math.isnan(va) and math.isnan(vb):
            continue
        assert va == vb, f"{k}: {va!r} != {vb!r}"

pytest.importorskip("scipy")  # the power backends need scipy, an optional extra

from pilotr import power, power_curve


def _spec():
    return {
        "name": "par", "seed": 11,
        "units": {"subject": {"n": 24}},
        "factors": [{"name": "group", "levels": ["control", "treatment"],
                     "contrasts": {"effect": [-0.5, 0.5]}, "between": "subject"}],
        "fixed": {"intercept": 100, "coefficients": {"effect": 5}},
        "response": {"family": "gaussian", "name": "score", "sigma": 10},
    }


def test_workers_is_validated():
    for bad in (0, -1, 1.5, "two", True):
        with pytest.raises(ValueError):
            power(_spec(), n_sims=2, workers=bad)


def test_power_parallel_matches_serial():
    serial = power(_spec(), n_sims=8, workers=1)
    parallel2 = power(_spec(), n_sims=8, workers=2)
    assert serial == parallel2


def test_power_curve_parallel_matches_serial():
    serial = power_curve(_spec(), subject_ns=[16, 24], n_sims=8, workers=1)
    parallel2 = power_curve(_spec(), subject_ns=[16, 24], n_sims=8, workers=2)
    assert serial == parallel2


def test_power_mixed_rejects_spec_without_items():
    pytest.importorskip("statsmodels")
    pytest.importorskip("pandas")
    from pilotr import power_mixed

    spec = {
        "name": "w", "seed": 1,
        "units": {"subject": {"n": 20}},
        "factors": [{"name": "cond", "levels": ["a", "b"],
                     "contrasts": {"cond": [-0.5, 0.5]}, "vary_within": "subject"}],
        "fixed": {"intercept": 6, "coefficients": {"cond": 0.05}},
        "random": {"subject": {"intercept_sd": 0.12, "slopes": {"cond": 0.04},
                               "correlations": {"intercept~cond": 0.2}}},
        "response": {"family": "gaussian", "name": "y", "sigma": 0.3},
    }
    with pytest.raises(ValueError, match="requires a crossed design with an item unit"):
        power_mixed(spec, n_sims=2)


def test_power_mixed_parallel_matches_serial():
    pytest.importorskip("statsmodels")
    pytest.importorskip("pandas")
    from pilotr import power_mixed

    spec = {
        "name": "parmix", "seed": 3,
        "units": {"subject": {"n": 10}, "item": {"n": 6}},
        "factors": [{"name": "cond", "levels": ["a", "b"],
                     "contrasts": {"cond": [-0.5, 0.5]}, "vary_within": "subject"}],
        "fixed": {"intercept": 6, "coefficients": {"cond": 0.05}},
        "random": {
            "subject": {"intercept_sd": 0.12, "slopes": {"cond": 0.04},
                        "correlations": {"intercept~cond": 0.2}},
            "item": {"intercept_sd": 0.08, "slopes": {"cond": 0.02},
                     "correlations": {"intercept~cond": -0.1}},
        },
        "response": {"family": "shifted_lognormal", "name": "RT", "sigma": 0.3, "shift": 200},
    }
    serial = power_mixed(spec, n_sims=4, workers=1)
    parallel2 = power_mixed(spec, n_sims=4, workers=2)
    assert_identical(serial, parallel2)
