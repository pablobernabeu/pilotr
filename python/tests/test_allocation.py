"""How units are allocated to levels and clusters, and the one allocation a specification may not
ask for.

Each between factor is assigned to blocks of its unit on its own, so two between factors over one
unit fell into the same or overlapping blocks. A 2 x 2 between-subjects design over 40 subjects
gave cells of 20, 0, 0 and 20, and its effects could not be estimated apart. The R twin's
tests/testthat/test-allocation.R runs the same cases and expects the same text, and
tools/parity/validate_cross.py compares the two validators on them.
"""
import os, sys, warnings
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
    return ("the factors %s are %s between '%s'. pilotr assigns the levels of each between factor "
            "to blocks of %ss on its own, so the blocks of these factors coincide or overlap, "
            "which leaves some combinations of their levels without %ss and confounds their "
            "effects. Encode the design as one between factor whose levels are the cells, give it "
            "the contrast columns of these factors and key each interaction as 'a:b', as "
            "spec_from_model() in R does for a fitted pilot. Specification version 0.4 will "
            "allocate several between factors jointly." % (names, quantifier, unit, unit, unit))


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


def test_two_factors_between_one_unit_are_refused_in_the_words_the_r_twin_uses():
    s = between_2x2()
    # What the refusal prevents: every subject in A's first block is in B's first block too.
    cells = Counter((r["A"], r["B"]) for r in simulate(s, validate=False).rows)
    assert cells == Counter({("a1", "b1"): 20, ("a2", "b2"): 20})

    refused(s, "the factors 'A' and 'B' are both between 'subject'. pilotr assigns the levels of "
               "each between factor to blocks of subjects on its own, so the blocks of these "
               "factors coincide or overlap, which leaves some combinations of their levels "
               "without subjects and confounds their effects. Encode the design as one between "
               "factor whose levels are the cells, give it the contrast columns of these factors "
               "and key each interaction as 'a:b', as spec_from_model() in R does for a fitted "
               "pilot. Specification version 0.4 will allocate several between factors jointly.")

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
