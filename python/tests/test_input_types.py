"""The number, path and table types users pass: NumPy scalars, pathlib paths and data frames.

The type checks tested ``isinstance(x, (int, float))``. A NumPy integer was therefore refused as
not a whole number and a NumPy float32 as not a number. A NumPy float64, which subclasses float,
was accepted and then carried into the simulated data as a NumPy scalar. A pathlib path was
refused as "not a JSON object" where the same path as a string was read, and a pandas data frame
was refused as "not a table of curve points". NumPy registers its scalar types with the
``numbers`` ABCs, so the checks now test those. Before anything is drawn, every number becomes a
Python int or float, on a copy of the caller's specification. NumPy and pandas are optional, so
the tests that need them skip without them.
"""
import copy, json, os, re, sys
from pathlib import Path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pytest

from pilotr import load_spec, simulate, solve_curve, target_n, validate_spec

SPEC = os.path.join(os.path.dirname(__file__), "..", "..", "spec", "examples")
DUPLICATE = ("the specification repeats the key '%s' within one object; JSON leaves a repeated "
             "key undefined, and R and Python read it differently")


def _two_group(n=16):
    return {
        "name": "types", "seed": 2024,
        "units": {"subject": {"n": n}},
        "factors": [{"name": "group", "levels": ["a", "b"],
                     "contrasts": {"effect": [-0.5, 0.5]}, "between": "subject"}],
        "fixed": {"intercept": 100, "coefficients": {"effect": 5}},
        "response": {"family": "gaussian", "name": "score", "sigma": 10},
    }


def _crossed():
    """A small crossed design, the kind power_mixed() takes."""
    return {
        "name": "crossed", "seed": 3,
        "units": {"subject": {"n": 10}, "item": {"n": 6}},
        "factors": [{"name": "cond", "levels": ["a", "b"],
                     "contrasts": {"cond": [-0.5, 0.5]}, "vary_within": ["subject", "item"]}],
        "fixed": {"intercept": 6, "coefficients": {"cond": 0.2}},
        "random": {"subject": {"intercept_sd": 0.12}, "item": {"intercept_sd": 0.08}},
        "response": {"family": "gaussian", "name": "y", "sigma": 0.3},
    }


def _as_numpy(x, np):
    """`x` with every int made a NumPy int64 and every float a NumPy float64."""
    if isinstance(x, dict):
        return {k: _as_numpy(v, np) for k, v in x.items()}
    if isinstance(x, list):
        return [_as_numpy(v, np) for v in x]
    if isinstance(x, bool):
        return x
    if isinstance(x, int):
        return np.int64(x)
    if isinstance(x, float):
        return np.float64(x)
    return x


def _same(a, b):
    """Equal as JSON text, where NaN equals NaN and a NumPy integer, which json cannot write,
    raises."""
    return json.dumps(a, sort_keys=True) == json.dumps(b, sort_keys=True)


def _plain_types(rows):
    return {type(v) for r in rows for v in r.values()} <= {int, float, str}


# ---- NumPy numbers ------------------------------------------------------------------------


def test_every_shipped_example_simulates_the_same_from_numpy_numbers():
    np = pytest.importorskip("numpy")
    for name in sorted(os.listdir(SPEC)):
        if not name.endswith(".json"):
            continue
        spec = load_spec(os.path.join(SPEC, name))
        rows = simulate(_as_numpy(spec, np)).rows
        assert rows == simulate(spec).rows, name
        # A NumPy float64 reached the data as one, and to_csv() then wrote "np.float64(...)".
        assert _plain_types(rows), name


def test_numpy_float32_and_small_integer_types_are_read_as_their_values():
    np = pytest.importorskip("numpy")
    spec = _two_group()
    spec["seed"] = np.uint16(2024)
    spec["units"]["subject"]["n"] = np.int32(16)
    spec["factors"][0]["contrasts"]["effect"] = [np.float32(-0.5), np.float32(0.5)]
    spec["fixed"]["intercept"] = np.float32(100)
    spec["fixed"]["coefficients"]["effect"] = np.float32(5)
    spec["response"]["sigma"] = np.float32(10)
    # Each value is exact in float32, so the data have to be those of the plain specification,
    # computed in double precision throughout.
    rows = simulate(spec).rows
    assert rows == simulate(_two_group()).rows
    assert _plain_types(rows)


def test_the_callers_specification_is_left_as_given():
    np = pytest.importorskip("numpy")
    spec = _as_numpy(_two_group(), np)
    before = copy.deepcopy(spec)
    simulate(spec)
    validate_spec(spec)
    assert spec == before
    assert type(spec["units"]["subject"]["n"]) is np.int64


def test_a_numpy_boolean_or_a_fractional_count_is_refused_as_before():
    np = pytest.importorskip("numpy")
    for bad in (np.bool_(True), np.float64(20.5), np.float32(20.5)):
        spec = _two_group()
        spec["units"]["subject"]["n"] = bad
        with pytest.raises(ValueError, match=re.escape(
                "'units.subject.n' must be a whole number of at least 1")):
            simulate(spec)
    spec = _two_group()
    spec["fixed"]["intercept"] = np.bool_(False)
    with pytest.raises(ValueError, match=re.escape("'fixed.intercept' must be a single number")):
        validate_spec(spec)


def test_power_curve_takes_a_numpy_grid_and_records_whole_counts():
    np = pytest.importorskip("numpy")
    pytest.importorskip("scipy")
    from pilotr import power_curve

    ref = power_curve(_two_group(), [16, 32, 48], n_sims=8)
    for grid in (np.arange(16, 49, 16), np.linspace(16, 48, 3), [16.0, 32.0, 48.0]):
        curve = power_curve(_two_group(), grid, n_sims=8)
        assert [type(r["n_subject"]) for r in curve] == [int, int, int]
        assert _same(curve, ref)


def test_numpy_numbers_in_the_specification_give_the_same_power():
    np = pytest.importorskip("numpy")
    pytest.importorskip("scipy")
    from pilotr import power

    assert _same(power(_as_numpy(_two_group(), np), n_sims=20), power(_two_group(), n_sims=20))


def test_workers_may_be_any_whole_number_as_in_r():
    np = pytest.importorskip("numpy")
    pytest.importorskip("scipy")
    from pilotr import power

    ref = power(_two_group(), n_sims=6)
    for workers in (np.int64(1), np.uint8(1), 1.0):
        assert _same(power(_two_group(), n_sims=6, workers=workers), ref)
    for bad in (np.bool_(True), np.float64(1.5), np.int64(0)):
        with pytest.raises(ValueError, match="workers must be a positive whole number"):
            power(_two_group(), n_sims=2, workers=bad)


# ---- paths --------------------------------------------------------------------------------


def test_a_path_object_is_read_wherever_a_path_string_is():
    name = os.path.join(SPEC, "between_2group_gaussian.json")
    path = Path(name)
    assert load_spec(path) == load_spec(name)
    assert simulate(path).rows == simulate(name).rows


def test_the_power_functions_read_a_path_object():
    pytest.importorskip("scipy")
    from pilotr import power, power_curve

    name = os.path.join(SPEC, "between_2group_gaussian.json")
    assert _same(power(Path(name), n_sims=10), power(name, n_sims=10))
    assert _same(power_curve(Path(name), [20, 40], n_sims=5),
                 power_curve(name, [20, 40], n_sims=5))


def test_power_mixed_reads_a_path_object(tmp_path):
    pytest.importorskip("statsmodels")
    pytest.importorskip("pandas")
    from pilotr import power_mixed

    path = tmp_path / "crossed.json"
    path.write_text(json.dumps(_crossed()), encoding="utf-8")
    assert _same(power_mixed(path, n_sims=2), power_mixed(_crossed(), n_sims=2))


def test_validate_spec_reads_a_path_as_the_r_twin_does(tmp_path):
    name = os.path.join(SPEC, "between_2group_gaussian.json")
    expected = load_spec(name)
    for given in (name, Path(name)):
        assert validate_spec(given) == expected

    bad = load_spec(name)
    bad["fixed"]["coefficients"] = {"efect": 0.5}
    path = tmp_path / "bad.json"
    path.write_text(json.dumps(bad), encoding="utf-8")
    with pytest.raises(ValueError, match="names 'efect', which is neither a contrast column"):
        validate_spec(path)

    noted = dict(load_spec(name), note="private")
    path = tmp_path / "noted.json"
    path.write_text(json.dumps(noted), encoding="utf-8")
    with pytest.raises(ValueError, match="unknown top-level field 'note'"):
        validate_spec(path)
    with pytest.warns(UserWarning, match="unknown top-level field 'note'"):
        assert validate_spec(path, strict=False) == noted

    # The file is read as load_spec() reads it, so a repeated key is refused as R refuses it.
    text = json.dumps(expected).replace('"coefficients": {', '"coefficients": {"grp": 9, ', 1)
    path = tmp_path / "repeated.json"
    path.write_text(text, encoding="utf-8")
    with pytest.raises(ValueError) as info:
        validate_spec(str(path))
    assert str(info.value) == DUPLICATE % "grp"

    with pytest.raises(ValueError, match="a design specification must be a JSON object"):
        validate_spec(5)


# ---- curves -------------------------------------------------------------------------------


def _curve():
    ns = [40, 60, 80, 100, 120, 140, 160, 180]
    pw = [0.3377, 0.4778, 0.5981, 0.6969, 0.7753, 0.8358, 0.8816, 0.9156]
    return [{"n_subject": n, "power": p, "n_sims": 400} for n, p in zip(ns, pw)]


def test_a_data_frame_curve_is_solved_as_its_records():
    pd = pytest.importorskip("pandas")
    assert target_n(pd.DataFrame(_curve())) == target_n(_curve())
    assert solve_curve(pd.DataFrame(_curve()), target=0.9) == solve_curve(_curve(), target=0.9)


def test_numpy_numbers_in_a_curve_and_in_the_solver_arguments():
    np = pytest.importorskip("numpy")
    plain = solve_curve(_curve(), target=0.75)
    records = [{"n_subject": np.int64(r["n_subject"]), "power": np.float64(r["power"]),
                "n_sims": np.int32(r["n_sims"])} for r in _curve()]
    assert solve_curve(records, target=0.75) == plain
    # 0.75 is exact in float32, so the solve is the same one.
    assert solve_curve(_curve(), target=np.float32(0.75)) == plain
    solved = solve_curve(_curve(), target=np.float64(0.75), level=np.float64(0.95))
    assert solved == plain
    assert all(type(solved[k]) is float for k in ("value", "lo", "hi", "target", "level"))
    # A NumPy count is a replicate count, which used to be looked up as a column name.
    bare = [{"n_subject": r["n_subject"], "power": r["power"]} for r in _curve()]
    assert solve_curve(bare, target=0.75, n=np.int64(400)) == solve_curve(bare, 0.75, n=400)
