"""Counterbalanced lists, new in specification version 0.4.

A within factor that sets ``counterbalance`` shows each subject each item at one level, rotated
across subjects and items as in a Latin square: subject s meets item t at combination
((s - 1) + (t - 1)) mod C of the counterbalanced factors' levels, read with the first-listed factor
varying slowest. Until 0.4, every within factor was crossed with every subject-item pair, so
planning a counterbalanced study doubled its observations per subject with two conditions. The R
twin's tests/testthat/test-counterbalance.R runs the same cases and expects the same text, and
tools/parity/validate_cross.py compares the two validators on them.
"""
import copy, json, os, sys, warnings
from collections import Counter

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pytest

from pilotr import load_spec, simulate, validate_spec
from pilotr.validate import _features_0_4, _required_spec_version

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, "..", ".."))


def cb_spec(n_subject=24, n_item=18):
    """A Gaussian version of the priming design, counterbalanced: 24 subjects and 18 items, the
    dimensions of the two-list encoding in tools/parity/cases."""
    return {
        "spec_version": "0.4", "name": "counterbalanced", "seed": 2417,
        "units": {"subject": {"n": n_subject}, "item": {"n": n_item}},
        "factors": [{"name": "cond", "levels": ["related", "unrelated"],
                     "contrasts": {"cond": [-0.5, 0.5]},
                     "vary_within": ["subject", "item"], "counterbalance": True}],
        "fixed": {"intercept": 6, "coefficients": {"cond": 0.06}},
        "random": {"subject": {"intercept_sd": 0.12, "slopes": {"cond": 0.04}},
                   "item": {"intercept_sd": 0.08, "slopes": {"cond": 0.02}}},
        "response": {"family": "gaussian", "name": "logRT", "sigma": 0.3}}


def cb3_spec(n_subject=18, n_item=12):
    """The same design with three conditions, which the two-list encoding cannot write."""
    s = cb_spec(n_subject, n_item)
    s["factors"][0] = {"name": "cond", "levels": ["a", "b", "c"],
                       "contrasts": {"c1": [-0.5, 0.5, 0], "c2": [-0.5, 0, 0.5]},
                       "vary_within": ["subject", "item"], "counterbalance": True}
    s["fixed"]["coefficients"] = {"c1": 0.05, "c2": 0.08}
    s["random"]["subject"]["slopes"] = {"c1": 0.04, "c2": 0.03}
    s["random"]["item"]["slopes"] = {"c1": 0.02, "c2": 0.02}
    return s


def rotated(d, levels):
    """The level each row takes by the rule: combination ((s - 1) + (t - 1)) mod C."""
    return [levels[(r["subject"] - 1 + r["item"] - 1) % len(levels)] for r in d.rows]


def refused(spec, msg, strict=True):
    with pytest.raises(ValueError) as e:
        validate_spec(spec, strict=strict)
    assert msg in str(e.value)


def silent(spec):
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        validate_spec(spec)


def warned(spec, msg, strict=True):
    with pytest.warns(UserWarning) as caught:
        validate_spec(spec, strict=strict)
    assert any(msg in str(w.message) for w in caught), [str(w.message) for w in caught]


def quiet_simulate(spec):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return simulate(spec)


def parity_case(name):
    path = os.path.join(REPO, "tools", "parity", "cases", name + ".json")
    if not os.path.exists(path):
        pytest.skip("the parity cases are in the repository, not the package")
    return load_spec(path)


def unequal(lead, subjects=None, items=None, lists=None):
    """The warning for unequal lists, given its lead and the counts that are not multiples."""
    out = lead + ", and the lists are unequal."
    if subjects is not None:
        out += (" The %d subjects are not a multiple of %d, so the lists differ in size by one "
                "subject, and an item is seen in some conditions by one subject more than in "
                "others." % (subjects, lists))
    if items is not None:
        out += (" The %d items are not a multiple of %d, so each list shows some conditions on "
                "one item more than others." % (items, lists))
    return out


def counts(pairs):
    """The number of rows for each (unit, condition) pair, in a fixed order."""
    return sorted(Counter(pairs).values())


def test_a_counterbalanced_factor_gives_each_subject_item_pair_one_row():
    s = cb_spec()
    silent(s)
    d = simulate(s)
    assert len(d) == 24 * 18
    assert len({(r["subject"], r["item"]) for r in d.rows}) == len(d)
    # The canonical order: subjects slowest, then items.
    assert d.column("subject") == [s for s in range(1, 25) for _ in range(18)]
    assert d.column("item") == list(range(1, 19)) * 24
    assert d.columns == ["subject", "item", "list", "cond", "logRT"]


def test_the_data_of_a_counterbalanced_design_name_each_subjects_list():
    # Subject s follows list ((s - 1) mod C) + 1, in a column after the units. A list is a set of
    # subjects who meet every item at the same level, and list k meets the first item at level k.
    d = simulate(cb_spec())
    assert all(type(v) is int for v in d.column("list"))
    assert d.column("list") == [(r["subject"] - 1) % 2 + 1 for r in d.rows]
    assert ([r["list"] for r in d.rows if r["item"] == 1]
            == [["related", "unrelated"].index(r["cond"]) + 1 for r in d.rows if r["item"] == 1])
    pairs = {(r["list"], r["item"], r["cond"]) for r in d.rows}
    assert len(pairs) == 2 * 18
    assert len({(lst, t) for lst, t, _ in pairs}) == len(pairs)
    d3 = simulate(cb3_spec())
    assert d3.column("list") == [(r["subject"] - 1) % 3 + 1 for r in d3.rows]
    # Two counterbalanced factors rotate over the four combinations of their levels.
    s = cb_spec(8, 8)
    s["factors"].append({"name": "prime", "levels": ["p1", "p2"],
                         "contrasts": {"pr": [-0.5, 0.5]},
                         "vary_within": ["subject", "item"], "counterbalance": True})
    s["fixed"]["coefficients"]["pr"] = 0.03
    assert simulate(s).column("list") == ([k for k in range(1, 5) for _ in range(8)] * 2)
    # 25 subjects over two lists put 13 in the first and 12 in the second.
    d = quiet_simulate(cb_spec(25, 18))
    assert Counter(r["list"] for r in d.rows if r["item"] == 1) == {1: 13, 2: 12}
    # The column comes before an additional grouping factor, as the units do.
    s = cb_spec()
    s["random"]["site"] = {"over": "subject", "n": 4, "intercept_sd": 0.1}
    assert simulate(s).columns == ["subject", "item", "list", "site", "cond", "logRT"]
    assert list(simulate(s).rows[0]) == ["subject", "item", "list", "site", "cond", "logRT"]
    # Without a counterbalanced factor, the data have no such column.
    crossed = cb_spec()
    del crossed["factors"][0]["counterbalance"]
    assert simulate(crossed).columns == ["subject", "item", "cond", "logRT"]


def test_a_counterbalanced_design_keeps_the_name_list_for_its_list_column():
    def taken(who):
        return ("the name 'list' is used for more than one column (the list column of a "
                "counterbalanced design and %s); every unit, grouping factor, factor, predictor "
                "and the response needs a column of its own" % who)

    s = cb_spec()
    s["factors"].append({"name": "list", "levels": ["l1", "l2"], "contrasts": {"lst": [-1, 1]},
                         "between": "subject"})
    refused(s, taken("factors[2]"))
    s = cb_spec()
    s["predictors"] = [{"name": "list", "varies_by": "subject"}]
    refused(s, taken("predictors[1]"))
    s = cb_spec()
    s["random"]["list"] = {"over": "subject", "n": 4, "intercept_sd": 0.1}
    refused(s, taken("random.list"))
    s = cb_spec()
    s["response"]["name"] = "list"
    refused(s, taken("the response"))
    s = cb_spec()
    s["factors"][0]["contrasts"] = {"list": [-0.5, 0.5]}
    s["fixed"]["coefficients"] = {"list": 0.06}
    s["random"]["subject"]["slopes"] = {"list": 0.04}
    s["random"]["item"]["slopes"] = {"list": 0.02}
    refused(s, "the contrast column 'list' of factors[1] is also the name of the list column of a "
               "counterbalanced design; a contrast column may share its own factor's name but no "
               "other column's")
    # Without a counterbalanced factor the name is free, as the two-list encoding's factor shows.
    s["factors"][0]["counterbalance"] = False
    silent(s)
    silent(parity_case("two_list_counterbalanced"))


def test_subject_s_meets_item_t_at_level_s_plus_t_mod_c():
    d = simulate(cb_spec())
    assert d.column("cond") == rotated(d, ["related", "unrelated"])
    d3 = simulate(cb3_spec())
    assert d3.column("cond") == rotated(d3, ["a", "b", "c"])
    # Subject 1 follows the first list, item 1 at the first level, and subject 2 the second list.
    assert [r["cond"] for r in d3.rows if r["subject"] == 1][:4] == ["a", "b", "c", "a"]
    assert [r["cond"] for r in d3.rows if r["subject"] == 2][:4] == ["b", "c", "a", "b"]


def test_every_item_appears_in_every_condition_and_every_subject_sees_each_equally_often():
    # 18 subjects and 12 items, both multiples of the three lists.
    d = simulate(cb3_spec())
    assert Counter((r["item"], r["cond"]) for r in d.rows) == {
        (t, c): 6 for t in range(1, 13) for c in "abc"}
    assert Counter((r["subject"], r["cond"]) for r in d.rows) == {
        (s, c): 4 for s in range(1, 19) for c in "abc"}
    # With 14 items, a subject sees some conditions on one item more than others.
    d = quiet_simulate(cb3_spec(18, 14))
    per_subject = Counter((r["subject"], r["cond"]) for r in d.rows)
    assert len(per_subject) == 18 * 3
    assert sorted(set(per_subject.values())) == [4, 5]
    assert Counter(r["subject"] for r in d.rows) == {s: 14 for s in range(1, 19)}


def test_several_counterbalanced_factors_rotate_jointly_the_first_listed_slowest():
    s = cb_spec(8, 8)
    s["factors"].append({"name": "prime", "levels": ["p1", "p2"],
                         "contrasts": {"pr": [-0.5, 0.5]},
                         "vary_within": ["subject", "item"], "counterbalance": True})
    s["fixed"]["coefficients"]["pr"] = 0.03
    silent(s)
    d = simulate(s)
    assert len(d) == 64
    cells = [(r["subject"] - 1 + r["item"] - 1) % 4 for r in d.rows]
    assert d.column("cond") == [["related", "unrelated"][c // 2] for c in cells]
    assert d.column("prime") == [["p1", "p2"][c % 2] for c in cells]


def test_beside_a_counterbalanced_factor_only_the_crossed_within_factors_expand_rows():
    s = cb_spec(4, 6)
    s["factors"].append({"name": "task", "levels": ["read", "listen"],
                         "contrasts": {"tk": [-0.5, 0.5]}, "vary_within": ["subject", "item"]})
    s["fixed"]["coefficients"]["tk"] = 0.1
    d = simulate(s)
    assert len(d) == 4 * 6 * 2
    # Each subject-item pair has one row per level of the crossed factor, in level order, and
    # keeps the one condition its rotation gives it.
    assert d.column("subject") == [s for s in range(1, 5) for _ in range(12)]
    assert d.column("item") == [t for t in range(1, 7) for _ in range(2)] * 4
    assert d.column("task") == ["read", "listen"] * 24
    assert d.column("cond") == rotated(d, ["related", "unrelated"])


def test_counterbalancing_draws_nothing_so_it_gives_the_responses_of_an_encoding_by_unit():
    # The rotation of two conditions, written with a factor between subjects and one between
    # items whose contrasts multiply to the condition, each subject and each item a level of its
    # own. The two specifications have the same rows in the same order, the same draws and the
    # same linear predictor term by term, so their responses agree to the last bit.
    native = cb_spec()
    encoded = copy.deepcopy(native)
    encoded["spec_version"] = "0.3"
    encoded["factors"] = [
        {"name": "subject_list", "levels": ["s%02d" % k for k in range(1, 25)],
         "contrasts": {"l": [1, -1] * 12}, "between": "subject"},
        {"name": "item_list", "levels": ["t%02d" % k for k in range(1, 19)],
         "contrasts": {"g": [-0.5, 0.5] * 9}, "between": "item"}]
    encoded["fixed"]["coefficients"] = {"l:g": 0.06}
    encoded["random"]["subject"]["slopes"] = {"l:g": 0.04}
    encoded["random"]["item"]["slopes"] = {"l:g": 0.02}
    assert simulate(native).column("logRT") == simulate(encoded).column("logRT")


def test_two_counterbalanced_conditions_have_the_cell_counts_of_the_two_list_encoding():
    native = simulate(parity_case("counterbalanced_two_conditions"))
    two_list = parity_case("two_list_counterbalanced")
    encoded = simulate(two_list)
    assert len(native) == len(encoded)
    # The condition of the two-list encoding is the product of its two contrasts, the column l_g
    # of R's model_data().
    lists, sets = two_list["factors"]
    by_list = dict(zip(lists["levels"], lists["contrasts"]["l"]))
    by_set = dict(zip(sets["levels"], sets["contrasts"]["g"]))
    cond = [by_list[r["list"]] * by_set[r["item_set"]] for r in encoded.rows]
    # Nine items in each condition for every subject, and 12 subjects in each for every item.
    for unit in ("subject", "item"):
        assert (counts((r[unit], r["cond"]) for r in native.rows)
                == counts(zip(encoded.column(unit), cond)))
    # Subject s follows list (s - 1) mod 2, and each list meets each condition 108 times.
    assert counts(((r["subject"] - 1) % 2, r["cond"]) for r in native.rows) == [108] * 4
    assert counts(zip(encoded.column("list"), cond)) == [108] * 4


def test_power_mixed_analyses_a_counterbalanced_design_and_its_two_list_encoding():
    # The condition stays a within factor, so the statsmodels backend fits a counterbalanced
    # design as it fits the crossed one. The two-list encoding has no within factor, which the
    # backend used to refuse. It now fits the model the specification implies, in which the
    # condition is the product column l_g of the encoding's interaction.
    pytest.importorskip("statsmodels")
    pytest.importorskip("pandas")
    from pilotr import power_mixed

    res = power_mixed(cb_spec(12, 8), n_sims=2)
    assert (res["n_attempted"], res["effects"]["cond"]["true_effect"]) == (2, 0.06)
    with pytest.warns(UserWarning, match="random-effect correlations"):
        res = power_mixed(parity_case("two_list_counterbalanced"), n_sims=1)
    assert list(res["effects"]) == ["l_g"] and res["effects"]["l_g"]["true_effect"] == 0.06
    assert res["n_returned"] == 1


def test_the_parity_cases_hold_the_designs_they_describe():
    d = simulate(parity_case("counterbalanced_two_conditions"))
    assert len(d) == 432
    assert len({(r["subject"], r["item"]) for r in d.rows}) == 432
    assert d.column("cond") == rotated(d, ["related", "unrelated"])

    d = simulate(parity_case("counterbalanced_three_conditions"))
    assert len(d) == 18 * 12 * 2
    assert d.column("cond") == rotated(d, ["a", "b", "c"])
    assert d.column("task") == ["read", "listen"] * (18 * 12)
    assert Counter((r["item"], r["cond"]) for r in d.rows) == {
        (t, c): 12 for t in range(1, 13) for c in "abc"}


def test_counterbalance_false_and_an_absent_field_cross_the_factor_alike():
    s = cb_spec()
    s["factors"][0]["counterbalance"] = False
    silent(s)
    crossed = copy.deepcopy(s)
    del crossed["factors"][0]["counterbalance"]
    a, b = simulate(s), simulate(crossed)
    assert (a.columns, a.rows) == (b.columns, b.rows)
    assert len(a) == 24 * 18 * 2


def test_counterbalance_needs_specification_version_0_4_whatever_its_value():
    msg = ('this specification uses factors counterbalance, which requires spec_version "0.4", '
           "but declares 0.3. An implementation that reads no version above 0.3 refuses such a "
           "field or ignores it, and ignoring it generates different data.")
    s = cb_spec()
    s["spec_version"] = "0.3"
    refused(s, msg)
    refused(s, msg, strict=False)
    # Releases that read 0.3 do not know the field, whatever its value.
    s["factors"][0]["counterbalance"] = False
    refused(s, msg)
    s["factors"][0]["counterbalance"] = True
    del s["spec_version"]
    refused(s, msg.replace("declares 0.3", "declares 0.2"))


def test_counterbalance_applies_to_a_within_factor_of_a_fully_crossed_design_with_items():
    s = cb_spec()
    s["factors"].append({"name": "group", "levels": ["x", "y"], "contrasts": {"grp": [-0.5, 0.5]},
                         "between": "subject", "counterbalance": True})
    refused(s, "factors[2].counterbalance applies only to a within factor, one that sets "
               "vary_within")

    s = cb_spec()
    del s["units"]["item"]
    del s["random"]["item"]
    s["factors"][0]["vary_within"] = "subject"
    refused(s, "factors[1].counterbalance rotates the levels over items, but the design has no "
               "item unit")

    s = cb_spec()
    s["units"]["item"]["per_subject"] = 6
    refused(s, "factors[1].counterbalance needs every subject to see every item, so it cannot be "
               "combined with units.item.per_subject")

    s = cb_spec()
    s["factors"][0]["counterbalance"] = "yes"
    refused(s, "factors[1].counterbalance must be true or false")
    s["factors"][0]["counterbalance"] = 1
    refused(s, "factors[1].counterbalance must be true or false")


def test_a_design_whose_subjects_or_items_are_not_a_multiple_of_the_lists_warns():
    lead = "the counterbalanced factor 'cond' rotates its 3 levels over 3 lists"
    s = cb3_spec(18, 12)
    silent(s)
    s["units"]["subject"]["n"] = 20
    warned(s, unequal(lead, subjects=20, lists=3))
    s["units"]["subject"]["n"] = 18
    s["units"]["item"]["n"] = 16
    warned(s, unequal(lead, items=16, lists=3))
    s["units"]["subject"]["n"] = 20
    warned(s, unequal(lead, 20, 16, 3), strict=False)
    # The data are still the rotation, with lists one subject or one item apart.
    d = quiet_simulate(s)
    assert d.column("cond") == rotated(d, ["a", "b", "c"])

    # Two counterbalanced factors rotate the combinations of their levels.
    s = cb_spec(6, 8)
    s["factors"].append({"name": "prime", "levels": ["p1", "p2"],
                         "contrasts": {"pr": [-0.5, 0.5]},
                         "vary_within": ["subject", "item"], "counterbalance": True})
    warned(s, unequal("the counterbalanced factors 'cond' and 'prime' rotate the 4 combinations "
                      "of their levels over 4 lists", subjects=6, lists=4))


def test_the_version_a_counterbalanced_specification_needs_is_0_4():
    s = cb_spec()
    assert _features_0_4(s) == ["factors counterbalance"]
    assert _required_spec_version(s) == "0.4"
    s["factors"][0]["counterbalance"] = False
    assert _required_spec_version(s) == "0.4"
    del s["factors"][0]["counterbalance"]
    assert _features_0_4(s) == []
    assert _required_spec_version(s) == "0.3"


def test_a_counterbalance_written_as_jsonlite_writes_it_reads_as_its_value(tmp_path):
    # jsonlite::write_json() writes every scalar as a one-element array by default, and R reads
    # [true] as true, so the Python reader does too. A null means the field is absent.
    s = cb_spec()
    s["factors"][0]["counterbalance"] = [True]
    path = tmp_path / "boxed.json"
    path.write_text(json.dumps(s), encoding="utf-8")
    assert simulate(load_spec(path)).column("logRT") == simulate(cb_spec()).column("logRT")
    s["factors"][0]["counterbalance"] = None
    path.write_text(json.dumps(s), encoding="utf-8")
    assert len(simulate(load_spec(path))) == 24 * 18 * 2
