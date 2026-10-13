"""Each grid point of `power_curve` draws its own replicates.

Every point used to reuse the specification's seed, so neighbouring points shared their random
numbers, while `solve_curve` fits them as independent binomial counts. Its nominal 95% intervals
then covered the true sample size 88% of the time on a tight grid. Point j of a sweep over J values
now takes ``replicate_seeds(seed, J)[j]`` as its seed, the rule of the R twin's ``sweep_spec()``.
"""
import copy, math, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pytest

pytest.importorskip("scipy")  # power() t-tests with scipy

from pilotr import power, power_curve, replicate_seeds, simulate


def _two_group(effect=0.6):
    return {
        "name": "seeds", "seed": 7,
        "units": {"subject": {"n": 20}},
        "factors": [{"name": "group", "levels": ["a", "b"],
                     "contrasts": {"effect": [-0.5, 0.5]}, "between": "subject"}],
        "fixed": {"intercept": 0, "coefficients": {"effect": effect}},
        "response": {"family": "gaussian", "name": "score", "sigma": 1},
    }


def _at(spec, n, seed, n_sims):
    """What `power` returns at subject count `n` and seed `seed`, keyed as a curve record."""
    s = copy.deepcopy(spec)
    s["units"]["subject"]["n"] = n
    s["seed"] = seed
    return {"n_subject": n, **power(s, n_sims=n_sims)}


def _same(record, full):
    """Whether a curve record holds the values `power` returned, NaN matching NaN."""
    def eq(a, b):
        return a == b or (isinstance(a, float) and isinstance(b, float)
                          and math.isnan(a) and math.isnan(b))
    return all(eq(record[k], full[k]) for k in record)


def test_grid_point_j_of_j_takes_replicate_seeds_seed_j_as_its_seed():
    spec = _two_group()
    sizes = [20, 40, 60]
    curve = power_curve(spec, subject_ns=sizes, n_sims=30)
    seeds = replicate_seeds(spec["seed"], len(sizes))
    for rec, n, seed in zip(curve, sizes, seeds):
        assert _same(rec, _at(spec, n, seed, 30)), n


def test_the_twins_give_each_grid_point_the_same_replicates():
    # The R package's tests/testthat/test-sweep.R pins the same counts for sweep_spec() over
    # power_design(). Shared seeds gave 8, 18 and 19.
    curve = power_curve(_two_group(), subject_ns=[20, 40, 60], n_sims=30)
    assert [r["n_significant"] for r in curve] == [4, 14, 19]


def test_adding_grid_points_leaves_the_earlier_points_unchanged():
    spec = _two_group()
    two = power_curve(spec, subject_ns=[20, 40], n_sims=30)
    three = power_curve(spec, subject_ns=[20, 40, 60], n_sims=30)
    assert len(three) == 3
    assert all(_same(a, b) for a, b in zip(three[:2], two))


def test_seeds_shared_gives_every_point_the_specifications_own_seed_as_before():
    spec = _two_group()
    sizes = [20, 40]
    curve = power_curve(spec, subject_ns=sizes, n_sims=30, seeds="shared")
    for rec, n in zip(curve, sizes):
        assert _same(rec, _at(spec, n, spec["seed"], 30)), n


def test_a_point_seed_follows_its_position_in_the_grid():
    # The data make the rule plain. Under the shared seed, the rows of the smaller design were a
    # prefix of the larger one's, and under the grid rule they are not.
    spec = _two_group()
    seeds = replicate_seeds(spec["seed"], 2)

    def rows(n, seed):
        s = copy.deepcopy(spec)
        s["units"]["subject"]["n"] = n
        s["seed"] = seed
        return [r["score"] for r in simulate(s).rows]

    assert rows(20, spec["seed"])[:10] == rows(40, spec["seed"])[:10]
    assert rows(20, seeds[0])[:10] != rows(40, seeds[1])[:10]
    # Reversing the grid hands each size the other seed.
    up = power_curve(spec, subject_ns=[20, 40], n_sims=30)
    down = power_curve(spec, subject_ns=[40, 20], n_sims=30)
    assert _same(up[0], _at(spec, 20, seeds[0], 30))
    assert _same(down[1], _at(spec, 20, seeds[1], 30))


def test_subject_counts_may_come_from_any_iterable():
    spec = _two_group()
    lazy = power_curve(spec, subject_ns=(n for n in [20, 40]), n_sims=10)
    listed = power_curve(spec, subject_ns=[20, 40], n_sims=10)
    assert len(lazy) == 2
    assert all(_same(a, b) for a, b in zip(lazy, listed))


def test_an_unknown_seeds_value_is_refused():
    with pytest.raises(ValueError, match="seeds must be 'independent' or 'shared'"):
        power_curve(_two_group(), subject_ns=[20, 40], n_sims=2, seeds="common")
