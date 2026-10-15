"""`model_formula` and `model_data`: the analysis model a specification implies, as the R twin
builds it.

The R twin derives its mixed model from the specification, while the Python `power_mixed` used to
fit one model written into its source. These are ports of the R functions. The formula strings
below are the R twin's, as `deparse(model_formula(spec, structure))` writes them for the shipped
examples, and tools/parity/formula_cross.py checks every example and parity case against R
itself, the modelling columns included.
"""
import math, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pytest

from pilotr import load_spec, model_formula, pilotr_example, simulate

# R 4.6.1, pilotr from this checkout: maximal, zcp and intercepts for each shipped example.
_R_FORMULAS = {
    "beta_proportion": [".y ~ grp"] * 3,
    "between_2group_gaussian": [".y ~ grp"] * 3,
    "crossed_mixed_rt": [
        ".y ~ cond + (1 + cond | subject) + (1 + cond | item)",
        ".y ~ cond + (1 + cond || subject) + (1 + cond || item)",
        ".y ~ cond + (1 | subject) + (1 | item)"],
    "nested_clusters": [".y ~ cond + (1 | subject) + (1 | item) + (1 | site)"] * 3,
    "ordinal_likert_between": [".y ~ grp"] * 3,
    "partial_crossing": [".y ~ cond + (1 | subject) + (1 | item)"] * 3,
    "poisson_counts_between": [".y ~ grp"] * 3,
    "reading_time_continuous": [
        ".y ~ SyntaxPC + CoherencePC + age + cond + SyntaxPC_age + CoherencePC_age + cond_age + "
        "(1 + SyntaxPC + CoherencePC || subject) + (1 | item)"] * 2 + [
        ".y ~ SyntaxPC + CoherencePC + age + cond + SyntaxPC_age + CoherencePC_age + cond_age + "
        "(1 | subject) + (1 | item)"],
}


@pytest.mark.parametrize("name", sorted(_R_FORMULAS))
def test_each_example_has_the_r_twins_formula_under_each_structure(name):
    spec = pilotr_example(name)
    got = [model_formula(spec, s) for s in ("maximal", "zcp", "intercepts")]
    assert got == _R_FORMULAS[name]
    assert model_formula(spec) == _R_FORMULAS[name][0]


def test_a_flag_decides_whether_a_group_is_correlated_as_in_the_r_twin():
    spec = load_spec(pilotr_example("crossed_mixed_rt"))
    spec["spec_version"] = "0.3"  # the version that introduced `correlated`
    for entry in spec["random"].values():
        del entry["correlations"]
    spec["random"]["subject"]["correlated"] = False
    spec["random"]["item"]["correlated"] = True
    assert model_formula(spec) == ".y ~ cond + (1 + cond || subject) + (1 + cond | item)"
    del spec["random"]["subject"]["correlated"]
    assert model_formula(spec) == ".y ~ cond + (1 + cond || subject) + (1 + cond | item)"


def test_a_random_slope_on_an_interaction_takes_its_product_column():
    spec = load_spec(pilotr_example("reading_time_continuous"))
    spec["spec_version"] = "0.3"  # the version that introduced interaction random slopes
    spec["random"]["item"]["slopes"] = {"cond:age": 0.02}
    assert model_formula(spec).endswith("(1 + SyntaxPC + CoherencePC || subject) + "
                                        "(1 + cond_age || item)")


def test_a_specification_without_any_term_is_refused_where_r_fails_to_parse():
    spec = load_spec(pilotr_example("between_2group_gaussian"))
    spec["fixed"]["coefficients"] = {}
    with pytest.raises(ValueError, match="no fixed coefficients and no random effects"):
        model_formula(spec)
    with pytest.raises(ValueError, match="structure must be 'maximal', 'zcp' or 'intercepts'"):
        model_formula(pilotr_example("crossed_mixed_rt"), "uncorrelated")


# ---- model_data ---------------------------------------------------------------------------


@pytest.fixture
def pd():
    return pytest.importorskip("pandas")


def test_model_data_adds_the_columns_in_the_r_twins_order(pd):
    from pilotr import model_data
    spec = pilotr_example("reading_time_continuous")
    d = model_data(spec, simulate(spec))
    # The column names and order R's model_data() gives the same data.
    assert list(d.columns) == [
        "subject", "item", "narration", "SyntaxPC", "CoherencePC", "age",
        "reading_time_per_word", ".y", "cond", "SyntaxPC_age", "CoherencePC_age", "cond_age"]
    # A lognormal response is analysed on the log scale, where its coefficients are.
    assert d[".y"].tolist() == [math.log(v) for v in d["reading_time_per_word"]]
    assert set(d["cond"]) == {-0.5, 0.5}
    # Products are taken in the key's order, as R's Reduce() takes them.
    assert d["cond_age"].tolist() == [c * a for c, a in zip(d["cond"], d["age"])]


def test_model_data_logs_a_shifted_response_after_removing_the_shift(pd):
    from pilotr import model_data
    spec = pilotr_example("crossed_mixed_rt")
    data = simulate(spec)
    d = model_data(spec, data)
    assert d[".y"].tolist() == [math.log(v - 200) for v in data.column("RT")]


def test_a_contrast_named_after_its_own_factor_leaves_no_later_column_empty(pd):
    # Every factor's level index is read before any contrast column is written. Written first, a
    # contrast named after its factor replaces the labels, and the factor's second contrast then
    # matches nothing. The R twin's model_data() left that column NA in every row until it took
    # the same rule.
    from pilotr import model_data
    spec = {
        "name": "own", "seed": 2,
        "units": {"subject": {"n": 4}, "item": {"n": 3}},
        "factors": [{"name": "cond", "levels": ["a", "b", "c"],
                     "contrasts": {"cond": [-1, 0, 1], "cond2": [1, -2, 1]},
                     "vary_within": ["subject", "item"]}],
        "fixed": {"intercept": 0, "coefficients": {"cond": 0.3, "cond2": 0.1}},
        "random": {"subject": {"intercept_sd": 0.5}},
        "response": {"family": "gaussian", "name": "y", "sigma": 1},
    }
    data = simulate(spec)
    d = model_data(spec, data)
    assert list(d.columns) == ["subject", "item", "cond", "y", ".y", "cond2"]
    position = {"a": 0, "b": 1, "c": 2}
    labels = data.column("cond")
    assert d["cond"].tolist() == [[-1.0, 0.0, 1.0][position[v]] for v in labels]
    assert d["cond2"].tolist() == [[1.0, -2.0, 1.0][position[v]] for v in labels]


def test_model_data_takes_a_data_frame_and_leaves_it_as_it_was(pd):
    from pilotr import model_data
    spec = pilotr_example("crossed_mixed_rt")
    data = simulate(spec)
    frame = pd.DataFrame(data.rows, columns=data.columns)
    before = frame.copy()
    d = model_data(spec, frame)
    assert d.equals(model_data(spec, data))
    assert frame.equals(before)
    with pytest.raises(TypeError, match="pilotr Dataset or a pandas DataFrame"):
        model_data(spec, data.rows)
