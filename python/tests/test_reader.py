"""How a specification is read: as R's jsonlite reads it, in UTF-8, with repeated keys refused,
JSON null read as an absent field, whole-number floats accepted, and every malformed shape
reported as a ValueError in the R validator's words.

The two fixtures under tests/fixtures were written by the R twin:

    spec <- build_spec(list(name = "write_json defaults", seed = 2024, design_kind = "within",
      include_items = TRUE, n_subject = 6, n_item = 4, factor_name = "cond", lev1 = "a",
      lev2 = "b", intercept = 6, effect = 0.25, subj_int_sd = 0.5, subj_slope_sd = 0.25,
      subj_corr = 0.5, item_int_sd = 0.5, item_slope_sd = 0.125, item_corr = -0.25,
      family = "gaussian", resp_name = "", sigma = 1))
    jsonlite::write_json(spec, "write_json_defaults.json")

    ord <- build_spec(list(name = "auto_unbox ordinal", seed = 7, design_kind = "between",
      factor_name = "group", lev1 = "a", lev2 = "b", n_subject = 20, intercept = 0,
      effect = 1, family = "ordinal", resp_name = "", thresholds = "0"))
    jsonlite::write_json(ord, "write_json_auto_unbox_ordinal.json", auto_unbox = TRUE)

The first writes every scalar as a one-element array, the second a single threshold as a bare
number. R reads and simulates both.
"""
import copy, json, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pytest

from pilotr import load_spec, simulate, validate_spec

HERE = os.path.dirname(__file__)
FIXTURES = os.path.join(HERE, "fixtures")
SPEC = os.path.join(HERE, "..", "..", "spec", "examples")

# R's first three responses from write_json_defaults.json, as simulate_design(load_spec(f))
# returns them, written with sprintf("%a").
R_WRITE_JSON_HEAD = ["0x1.a154c985f06f7p+2", "0x1.96339c0ebedfap+2", "0x1.806f694467382p+2"]
R_ORDINAL = [2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 2, 1, 2, 2, 2, 2, 2]

DUPLICATE = ("the specification repeats the key '%s' within one object; JSON leaves a repeated "
             "key undefined, and R and Python read it differently")
SEED_RANGE = ("'seed' must be a whole number between -9007199254740991 and 9007199254740991 "
              "(2^53 - 1), the range in which R and Python read a JSON integer exactly")
SPEC_VERSION_SHAPE = "'spec_version' must be a single string of the form 'major.minor'"
FAMILIES = ("'response.family' must be one of gaussian, lognormal, shifted_lognormal, exgaussian, "
            "bernoulli, poisson, ordinal, beta")


def _example(name):
    with open(os.path.join(SPEC, name + ".json"), encoding="utf-8") as f:
        return json.load(f)


def _write(tmp_path, text, name="spec.json", bom=False):
    path = tmp_path / name
    path.write_bytes((b"\xef\xbb\xbf" if bom else b"") + text.encode("utf-8"))
    return str(path)


def _refusal(spec):
    with pytest.raises(ValueError) as info:
        validate_spec(spec)
    return str(info.value)


def test_a_write_json_file_reads_as_r_reads_it():
    spec = load_spec(os.path.join(FIXTURES, "write_json_defaults.json"))
    assert spec["seed"] == 2024 and spec["spec_version"] == "0.3"
    assert spec["factors"][0]["name"] == "cond"
    scores = simulate(spec).column("score")
    assert scores[:3] == [float.fromhex(h) for h in R_WRITE_JSON_HEAD]


def test_a_bare_single_threshold_simulates():
    path = os.path.join(FIXTURES, "write_json_auto_unbox_ordinal.json")
    assert simulate(load_spec(path)).column("rating") == R_ORDINAL
    assert simulate(path).column("rating") == R_ORDINAL


def test_one_element_arrays_and_empty_arrays_read_as_jsonlite_reads_them(tmp_path):
    spec = _example("crossed_mixed_rt")
    boxed = copy.deepcopy(spec)
    boxed["seed"] = [spec["seed"]]
    boxed["spec_version"] = ["0.3"]
    boxed["response"]["family"] = [spec["response"]["family"]]
    boxed["random"]["item"]["intercept_sd"] = [spec["random"]["item"]["intercept_sd"]]
    path = _write(tmp_path, json.dumps(boxed))
    assert simulate(load_spec(path)).rows == simulate(dict(spec, spec_version="0.3")).rows

    between = _example("between_2group_gaussian")
    between["factors"][0]["between"] = ["subject"]
    between["fixed"]["coefficients"] = []
    path = _write(tmp_path, json.dumps(between), "empty.json")
    loaded = load_spec(path)
    assert loaded["factors"][0]["between"] == "subject"
    assert loaded["fixed"]["coefficients"] == {}

    nested = _example("crossed_mixed_rt")
    nested["random"]["subject"]["slopes"] = []
    nested["random"]["subject"]["correlations"] = []
    loaded = load_spec(_write(tmp_path, json.dumps(nested), "slopes.json"))
    assert loaded["random"]["subject"]["slopes"] == {}
    assert loaded["random"]["subject"]["correlations"] == {}


def test_whole_number_floats_simulate_as_their_integers():
    cases = [
        ("between_2group_gaussian", lambda s: s["units"]["subject"].__setitem__("n", 8.0),
         lambda s: s["units"]["subject"].__setitem__("n", 8)),
        ("partial_crossing", lambda s: s["units"]["item"].__setitem__("per_subject", 2.0),
         lambda s: s["units"]["item"].__setitem__("per_subject", 2)),
        ("between_2group_gaussian", lambda s: s["response"].__setitem__("round", 2.0),
         lambda s: s["response"].__setitem__("round", 2)),
        ("nested_clusters", lambda s: s["random"]["site"].__setitem__("n", 2.0),
         lambda s: s["random"]["site"].__setitem__("n", 2)),
    ]
    for name, as_float, as_int in cases:
        a, b = _example(name), _example(name)
        as_float(a)
        as_int(b)
        assert simulate(a).rows == simulate(b).rows, name


def test_validate_spec_raises_only_value_error_for_malformed_shapes():
    spec = _example("crossed_mixed_rt")
    s = copy.deepcopy(spec)
    s["response"]["family"] = ["gaussian"]
    assert FAMILIES in _refusal(s)
    s = copy.deepcopy(spec)
    s["response"]["family"] = {"name": "gaussian"}
    assert FAMILIES in _refusal(s)
    s = copy.deepcopy(spec)
    s["response"] = "gaussian"
    assert "'response' must be an object" in _refusal(s)
    s = copy.deepcopy(spec)
    s["random"] = ["subject"]
    assert "'random' must be an object keyed by grouping factor" in _refusal(s)
    s = copy.deepcopy(spec)
    s["random"]["subject"]["slopes"] = 5
    assert "random.subject.slopes must be an object" in _refusal(s)
    s = copy.deepcopy(spec)
    s["predictors"] = 5
    assert "'predictors' must be an array of predictor objects" in _refusal(s)
    # Falsy but not empty: R refuses both, and Python used to read them as no random effects.
    for value in (0, False):
        s = copy.deepcopy(spec)
        s["random"] = value
        assert "'random' must be an object keyed by grouping factor" in _refusal(s)


def test_a_spec_version_that_is_not_one_string_or_number_is_refused():
    # R stopped with a base error on an empty array or object and read a two-element array by
    # its first element; both twins now refuse every such value in the same words.
    spec = _example("between_2group_gaussian")
    for value in ([], {}, ["0.3", "0.4"], [None], True):
        assert SPEC_VERSION_SHAPE in _refusal(dict(spec, spec_version=value))


def test_a_list_valued_varies_by_is_reported_without_its_value():
    # Only a string is quoted back, as R does: R pasted a two-element array in as two problems.
    spec = _example("between_2group_gaussian")
    spec["predictors"] = [{"name": "x", "varies_by": ["subject", "item"]}]
    spec["fixed"]["coefficients"]["x"] = 0.1
    lines = [line.strip() for line in _refusal(spec).splitlines()]
    assert "- predictors[1].varies_by must be 'subject', 'item' or 'observation'" in lines
    assert not any(", not '" in line for line in lines)


def test_a_repeated_key_is_refused_at_parse_time(tmp_path):
    text = json.dumps(_example("between_2group_gaussian"))
    dup = text.replace('"coefficients": {', '"coefficients": {"grp": 50, ', 1)
    path = _write(tmp_path, dup)
    for validate in (True, False):
        with pytest.raises(ValueError) as info:
            load_spec(path, validate=validate)
        assert str(info.value) == DUPLICATE % "grp"
    top = _write(tmp_path, '{"seed": 7, ' + text[1:], "top.json")
    with pytest.raises(ValueError, match="repeats the key 'seed'"):
        load_spec(top, validate=False)


def test_utf8_with_or_without_a_byte_order_mark(tmp_path):
    spec = _example("between_2group_gaussian")
    spec["factors"][0]["levels"] = ["fácil", "Łatwy"]
    text = json.dumps(spec, ensure_ascii=False)
    for bom in (False, True):
        path = _write(tmp_path, text, "bom.json" if bom else "plain.json", bom=bom)
        loaded = load_spec(path)
        assert loaded["factors"][0]["levels"] == ["fácil", "Łatwy"]
        assert set(simulate(loaded).column("group")) == {"fácil", "Łatwy"}


def test_null_in_an_optional_field_reads_as_absent():
    spec = _example("crossed_mixed_rt")
    spec["spec_version"] = "0.3"
    spec["predictors"] = [{"name": "z", "varies_by": "subject"}]
    spec["fixed"]["coefficients"]["z"] = 0.1
    ref = simulate(spec).rows

    nulls = copy.deepcopy(spec)
    nulls["predictors"][0].update(mean=None, sd=None, dist=None, reliability=None)
    assert simulate(nulls).rows == ref

    nulls = copy.deepcopy(spec)
    nulls["spec_version"] = None
    assert simulate(nulls).rows == ref

    bare = _example("crossed_mixed_rt")
    bare["random"]["item"] = {"intercept_sd": 0.08}
    nulls = copy.deepcopy(bare)
    nulls["random"]["item"].update(slopes=None, correlations=None)
    assert simulate(nulls).rows == simulate(bare).rows

    between = _example("between_2group_gaussian")
    nulls = copy.deepcopy(between)
    nulls["predictors"] = None
    nulls["random"] = None
    assert simulate(nulls).rows == simulate(between).rows

    empty = copy.deepcopy(between)
    empty["factors"] = None
    empty["fixed"]["coefficients"] = {}
    assert len(simulate(empty).rows) == between["units"]["subject"]["n"]


def test_a_null_item_unit_is_refused_in_r_words():
    spec = _example("between_2group_gaussian")
    spec["units"]["item"] = None
    assert "'units.item' must be an object" in _refusal(spec)


def test_a_seed_must_lie_within_the_exact_range():
    spec = _example("between_2group_gaussian")
    for seed in (2**53, -(2**53), 9007199254740993, 123456789012345678, 1.5):
        spec["seed"] = seed
        assert SEED_RANGE in _refusal(spec)
    for seed in (2**53 - 1, -(2**53 - 1), 9007199254740991.0):
        spec["seed"] = seed
        validate_spec(spec)
