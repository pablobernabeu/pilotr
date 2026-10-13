"""How units are allocated to levels and clusters.

Up to specification version 0.3, each between factor was assigned to blocks of its unit on its
own, so two between factors over one unit fell into the same or overlapping blocks. A 2 x 2
between-subjects design over 40 subjects gave cells of 20, 0, 0 and 20, and its effects could not
be estimated apart, so a specification declaring a version below 0.4 is refused for them. Version
0.4 allocates them jointly, one block of units per combination of their levels. The R twin's
tests/testthat/test-allocation.R runs the same cases and expects the same text, and
tools/parity/validate_cross.py compares the two validators on them.
"""
import copy, os, sys, warnings
from collections import Counter

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pytest

from pilotr import load_spec, simulate, validate_spec

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, "..", ".."))
SPEC = os.path.join(REPO, "spec", "examples")


def _ex(name):
    return load_spec(os.path.join(SPEC, name + ".json"))


def between_2x2(n=40):
    return {
        "name": "between_2x2", "seed": 1,
        "units": {"subject": {"n": n}},
        "factors": [
            {"name": "A", "levels": ["a1", "a2"], "contrasts": {"a": [-0.5, 0.5]},
             "between": "subject"},
            {"name": "B", "levels": ["b1", "b2"], "contrasts": {"b": [-0.5, 0.5]},
             "between": "subject"}],
        "fixed": {"intercept": 0, "coefficients": {"a": 0.5, "b": 0.3, "a:b": 0.2}},
        "random": {},
        "response": {"family": "gaussian", "name": "y", "sigma": 1}}


def aliased(names, quantifier, unit):
    """The refusal for factors between one unit, given the quoted names, "both" or "all", and
    the unit."""
    return ("the factors %s are %s between '%s'. Below specification version 0.4, pilotr assigns "
            "the levels of each between factor to blocks of %ss on its own. The blocks of these "
            "factors then coincide or overlap, which leaves some combinations of their levels "
            "without %ss and confounds their effects. Declare spec_version \"0.4\", which "
            "allocates several between factors jointly. Alternatively, encode the design as one "
            "between factor whose levels are the cells, give it the contrast columns of these "
            "factors and key each interaction as 'a:b', as spec_from_model() in R does for a "
            "fitted pilot. Earlier releases of pilotr read that encoding too."
            % (names, quantifier, unit, unit, unit))


def cells_of(d, sep=" "):
    """The labels of each row's levels of A and B, for reading off a joint allocation."""
    return [r["A"] + sep + r["B"] for r in d.rows]


def declare(spec, version):
    """The specification, declaring `version`."""
    spec["spec_version"] = version
    return spec


def refused(spec, msg):
    with pytest.raises(ValueError) as e:
        validate_spec(spec)
    assert msg in str(e.value)


def silent(spec):
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        validate_spec(spec)


def first_rows(ds):
    """One row per subject, for counting what each subject was assigned."""
    seen, out = set(), []
    for r in ds.rows:
        if r["subject"] not in seen:
            seen.add(r["subject"])
            out.append(r)
    return out


def parity_case(name):
    """A worked encoding from the specification, which lives in the repository beside the
    package."""
    path = os.path.join(REPO, "tools", "parity", "cases", name + ".json")
    if not os.path.exists(path):
        pytest.skip("the parity cases are in the repository, not the package")
    return load_spec(path)


def contrast(spec, d, factor, column):
    """The contrast value of each row, as R's model_data() writes it from the level labels."""
    f = next(f for f in spec["factors"] if f["name"] == factor)
    return [f["contrasts"][column][f["levels"].index(r[factor])] for r in d.rows]


def test_two_factors_between_one_unit_are_refused_below_0_4_in_the_words_the_r_twin_uses():
    s = between_2x2()
    refused(s, "the factors 'A' and 'B' are both between 'subject'. Below specification version "
               "0.4, pilotr assigns the levels of each between factor to blocks of subjects on its "
               "own. The blocks of these factors then coincide or overlap, which leaves some "
               "combinations of their levels without subjects and confounds their effects. "
               "Declare spec_version \"0.4\", which allocates several between factors jointly. "
               "Alternatively, encode the design as one between factor whose levels are the "
               "cells, give it the contrast columns of these factors and key each interaction as "
               "'a:b', as spec_from_model() in R does for a fitted pilot. Earlier releases of "
               "pilotr read that encoding too.")
    # Declaring 0.3 changes nothing: that version's rule put these factors in the same blocks,
    # leaving cells of 20, 0, 0 and 20 subjects.
    refused(declare(s, "0.3"), aliased("'A' and 'B'", "both", "subject"))

    s2 = _ex("between_2group_gaussian")
    s2["factors"].append(s["factors"][1])
    refused(s2, aliased("'group' and 'B'", "both", "subject"))


def test_three_factors_between_one_unit_and_two_between_items_are_refused_alike():
    s = between_2x2(12)
    s["factors"].append({"name": "C", "levels": ["c1", "c2", "c3"],
                         "contrasts": {"c1": [-1, 1, 0], "c2": [-1, 0, 1]}, "between": "subject"})
    refused(s, aliased("'A', 'B' and 'C'", "all", "subject"))

    s = _ex("crossed_mixed_rt")
    s["factors"].append({"name": "frequency", "levels": ["low", "high"],
                         "contrasts": {"freq": [-0.5, 0.5]}, "between": "item"})
    s["factors"].append({"name": "length", "levels": ["short", "long"],
                         "contrasts": {"len": [-0.5, 0.5]}, "between": "item"})
    refused(s, aliased("'frequency' and 'length'", "both", "item"))


def test_one_between_factor_per_unit_stays_valid():
    # One factor between subjects and one between items: the two-list encoding relies on it.
    s = between_2x2()
    s["units"]["item"] = {"n": 4}
    s["factors"][1]["between"] = "item"
    silent(s)

    # A between factor beside a within factor.
    s = _ex("crossed_mixed_rt")
    s["factors"].append({"name": "group", "levels": ["x", "y"],
                         "contrasts": {"grp": [-0.5, 0.5]}, "between": "subject"})
    silent(s)


def test_a_between_factors_blocks_follow_the_rule_the_specification_states():
    # Unit u of N takes level floor((u - 1) * L / N), so 3 levels over 10 subjects give blocks of
    # 4, 3 and 3. Equal blocks need N to be a multiple of L.
    s = _ex("between_2group_gaussian")
    s["units"]["subject"]["n"] = 10
    s["factors"][0] = {"name": "group", "levels": ["x", "y", "z"],
                       "contrasts": {"g1": [-1, 1, 0], "g2": [-1, 0, 1]}, "between": "subject"}
    s["fixed"]["coefficients"] = {"g1": 1}
    assert simulate(s).column("group") == ["x"] * 4 + ["y"] * 3 + ["z"] * 3


def test_a_between_factor_and_a_grouping_factor_over_one_unit_nest_as_the_specification_says():
    s = _ex("nested_clusters")
    s["factors"] = [{"name": "grp", "levels": ["control", "treatment"],
                     "contrasts": {"g": [-0.5, 0.5]}, "between": "subject"}]
    s["fixed"]["coefficients"] = {"g": 0.3}

    def levels_per_site(spec):
        rows = first_rows(simulate(spec))
        sites = sorted({r["site"] for r in rows})
        return [len({r["grp"] for r in rows if r["site"] == k}) for k in sites]

    # 12 sites, a multiple of the 2 levels: each site lies wholly in one condition.
    assert levels_per_site(s) == [1] * 12
    # 5 sites: the middle site straddles the two conditions.
    s["random"]["site"]["n"] = 5
    assert levels_per_site(s) == [1, 1, 2, 1, 1]


def test_a_2x2_between_design_written_as_one_factor_of_cells_keeps_every_cell():
    s = parity_case("between_cells_2x2")
    assert Counter(r["cell"] for r in first_rows(simulate(s))) == {
        "a1.b1": 10, "a1.b2": 10, "a2.b1": 10, "a2.b2": 10}
    s["units"]["subject"]["n"] = 80
    assert Counter(r["cell"] for r in first_rows(simulate(s))) == {
        "a1.b1": 20, "a1.b2": 20, "a2.b1": 20, "a2.b2": 20}


def test_a_two_list_counterbalanced_design_shows_each_subject_each_item_once():
    s = parity_case("two_list_counterbalanced")
    d = simulate(s)
    assert len(d) == 24 * 18
    assert len({(r["subject"], r["item"]) for r in d.rows}) == len(d)
    # The condition is the product of the two contrasts, the column l_g of R's model_data().
    # Every item appears in both conditions across subjects and every subject sees both. Each
    # list meets each condition 108 times.
    cond = [lv * gv for lv, gv in zip(contrast(s, d, "list", "l"),
                                      contrast(s, d, "item_set", "g"))]
    for unit, n in (("item", 18), ("subject", 24)):
        per_unit = {}
        for r, c in zip(d.rows, cond):
            per_unit.setdefault(r[unit], set()).add(c)
        assert len(per_unit) == n and all(len(v) == 2 for v in per_unit.values())
    assert Counter(zip(d.column("list"), cond)) == {
        ("list1", -0.5): 108, ("list1", 0.5): 108, ("list2", -0.5): 108, ("list2", 0.5): 108}


def test_randomisation_within_clusters_puts_both_arms_in_every_site():
    s = parity_case("within_cluster_randomised")
    d = simulate(s)
    assert len(d) == 120
    arms = Counter(zip(d.column("site"), contrast(s, d, "arm", "trt")))
    assert arms == {(k, t): 5 for k in range(1, 13) for t in (-0.5, 0.5)}
    # The level labels name the site each level falls in.
    assert all(int(r["arm"][1:3]) == r["site"] for r in d.rows)


# Version 0.4 allocates the between factors of one unit jointly. The units fall into one block per
# combination of levels, a cell, in order, with the first-listed factor varying slowest.

def test_a_0_4_specification_allocates_two_factors_between_subjects_jointly():
    s = declare(between_2x2(), "0.4")
    silent(s)
    # Subjects 1 to 10 take the first cell, 11 to 20 the second and so on, B varying fastest.
    assert cells_of(simulate(s)) == (["a1 b1"] * 10 + ["a1 b2"] * 10 + ["a2 b1"] * 10
                                     + ["a2 b2"] * 10)


def test_a_3x2_design_over_12_subjects_puts_two_subjects_in_each_cell():
    s = declare(between_2x2(12), "0.4")
    s["factors"][0] = {"name": "A", "levels": ["a1", "a2", "a3"],
                       "contrasts": {"a1": [-0.5, 0.5, 0], "a2": [-0.5, 0, 0.5]},
                       "between": "subject"}
    s["fixed"]["coefficients"] = {"a1": 0.4, "a2": 0.7, "b": 0.3}
    assert cells_of(simulate(s)) == [c for c in ("a1 b1", "a1 b2", "a2 b1", "a2 b2", "a3 b1",
                                                 "a3 b2") for _ in range(2)]


def test_the_cells_are_equal_only_when_the_units_are_a_multiple_of_the_cells():
    # Unit u of 10 takes cell floor((u - 1) * 4 / 10), so the four cells hold 3, 2, 3 and 2.
    s = declare(between_2x2(10), "0.4")
    assert Counter(cells_of(simulate(s))) == {"a1 b1": 3, "a1 b2": 2, "a2 b1": 3, "a2 b2": 2}


def test_two_factors_between_items_are_allocated_jointly_too():
    s = declare(_ex("crossed_mixed_rt"), "0.4")
    s["factors"].append({"name": "frequency", "levels": ["low", "high"],
                         "contrasts": {"freq": [-0.5, 0.5]}, "between": "item"})
    s["factors"].append({"name": "length", "levels": ["short", "long"],
                         "contrasts": {"len": [-0.5, 0.5]}, "between": "item"})
    s["fixed"]["coefficients"] = {"cond": 0.05, "freq": 0.03, "len": 0.02}
    d = simulate(s)
    items = {}
    for r in d.rows:
        items.setdefault(r["item"], r["frequency"] + " " + r["length"])
    # The 24 items fall into four cells of six, in item order.
    assert [items[t] for t in range(1, 25)] == [
        c for c in ("low short", "low long", "high short", "high long") for _ in range(6)]


def test_joint_allocation_draws_nothing_so_it_gives_the_data_of_one_factor_of_cells():
    joint = declare(between_2x2(), "0.4")
    joint["units"]["item"] = {"n": 3}
    joint["random"] = {"subject": {"intercept_sd": 1}, "item": {"intercept_sd": 0.5}}
    # The worked encoding that every version reads: one factor whose levels are the cells, in the
    # same order, carrying both contrast columns.
    cells = copy.deepcopy(joint)
    cells["spec_version"] = "0.3"
    cells["factors"] = [{"name": "cell", "levels": ["a1.b1", "a1.b2", "a2.b1", "a2.b2"],
                         "contrasts": {"a": [-0.5, -0.5, 0.5, 0.5], "b": [-0.5, 0.5, -0.5, 0.5]},
                         "between": "subject"}]
    dj, dc = simulate(joint), simulate(cells)
    assert cells_of(dj, ".") == dc.column("cell")
    assert dj.column("y") == dc.column("y")


def test_the_parity_cases_for_joint_allocation_hold_the_cells_they_describe():
    joint = simulate(parity_case("two_between_joint"))
    cells = simulate(parity_case("between_cells_2x2"))
    assert cells_of(joint, ".") == cells.column("cell")
    assert joint.column("score") == cells.column("score")

    d = simulate(parity_case("three_by_two_joint"))
    assert [r["dose"] + " " + r["group"] for r in first_rows(d)] == [
        c for c in ("d0 g1", "d0 g2", "d1 g1", "d1 g2", "d2 g1", "d2 g2") for _ in range(2)]
    items = [r["freq"] + " " + r["len"] for r in d.rows if r["subject"] == 1 and r["cond"] == "x"]
    assert items == [c for c in ("lo short", "lo long", "hi short", "hi long") for _ in range(2)]


def test_the_cells_of_a_joint_allocation_nest_in_a_grouping_factor_over_the_same_unit():
    s = declare(_ex("nested_clusters"), "0.4")
    s["factors"] = between_2x2()["factors"]
    s["fixed"]["coefficients"] = {"a": 0.3, "b": 0.2}

    def cells_per_site(spec):
        rows = first_rows(simulate(spec))
        sites = sorted({r["site"] for r in rows})
        return [len({r["A"] + " " + r["B"] for r in rows if r["site"] == k}) for k in sites]

    # 12 sites, a multiple of the 4 cells: each site lies wholly in one cell.
    assert cells_per_site(s) == [1] * 12
    # 6 sites: sites 2 and 5 straddle two cells each.
    s["random"]["site"]["n"] = 6
    assert cells_per_site(s) == [1, 2, 1, 1, 2, 1]


def test_declaring_0_4_changes_nothing_for_one_between_factor_per_unit():
    # With one factor, the joint rule is the block rule of earlier versions, and it draws nothing.
    # tools/parity/golden.json holds every parity case to the same.
    for name in ("between_2group_gaussian", "nested_clusters"):
        before = simulate(_ex(name))
        after = simulate(declare(_ex(name), "0.4"))
        assert (after.columns, after.rows) == (before.columns, before.rows)
    s = parity_case("two_list_counterbalanced")
    before = simulate(s)
    after = simulate(declare(s, "0.4"))
    assert (after.columns, after.rows) == (before.columns, before.rows)
