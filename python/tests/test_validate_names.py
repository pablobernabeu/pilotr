"""Names, columns and structure.

Each specification below used to validate and then move an effect, rescale a variance or
overwrite a column without a word. The R twin's tests/testthat/test-validate-names.R runs the
same cases and expects the same text, and tools/parity/validate_cross.py compares the two
validators on them.
"""
import os, sys, warnings
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pytest

from pilotr import load_spec, simulate, validate_spec

SPEC = os.path.join(os.path.dirname(__file__), "..", "..", "spec", "examples")

COLUMN_OF_ITS_OWN = ("; every unit, grouping factor, factor, predictor and the response needs a "
                     "column of its own")
OWN_FACTOR_ONLY = "; a contrast column may share its own factor's name but no other column's"


def _ex(name):
    return load_spec(os.path.join(SPEC, name + ".json"))


def between():
    return _ex("between_2group_gaussian")


def crossed():
    return _ex("crossed_mixed_rt")


def refused(spec, msg):
    with pytest.raises(ValueError) as e:
        validate_spec(spec)
    assert msg in str(e.value)


def test_two_columns_with_one_name_are_refused():
    s = between()
    s["factors"][0]["name"] = "subject"
    refused(s, "the name 'subject' is used for more than one column (the subject unit and "
               "factors[1])" + COLUMN_OF_ITS_OWN)

    s = between()
    s["response"]["name"] = "group"
    refused(s, "the name 'group' is used for more than one column (factors[1] and the response)"
            + COLUMN_OF_ITS_OWN)

    s = _ex("nested_clusters")
    s["predictors"] = [{"name": "site", "varies_by": "subject"}]
    s["response"]["name"] = "site"
    refused(s, "the name 'site' is used for more than one column (random.site, predictors[1] and "
               "the response)" + COLUMN_OF_ITS_OWN)


def test_a_contrast_column_belongs_to_one_factor_and_shares_no_other_columns_name():
    s = between()
    s["factors"][0]["contrasts"] = {"x": [-0.5, 0.5]}
    s["factors"].append({"name": "block", "levels": ["p", "q"],
                         "contrasts": {"x": [-0.5, 0.5]}, "between": "subject"})
    s["fixed"]["coefficients"] = {"x": 5}
    refused(s, "the contrast column 'x' is defined by more than one factor (factors[1] and "
               "factors[2]); each contrast column has to belong to a single factor")

    # A within factor reusing the between factor's column put the whole effect on the block.
    s = between()
    s["factors"][0]["contrasts"] = {"effect": [-0.5, 0.5]}
    s["factors"].append({"name": "block", "levels": ["p", "q"],
                         "contrasts": {"effect": [-0.5, 0.5]}, "vary_within": "subject"})
    s["fixed"]["coefficients"] = {"effect": 5}
    refused(s, "the contrast column 'effect' is defined by more than one factor (factors[1] and "
               "factors[2]); each contrast column has to belong to a single factor")

    # A predictor named like the contrast replaced the group difference with its own draws.
    s = between()
    s["factors"][0]["contrasts"] = {"effect": [-0.5, 0.5]}
    s["predictors"] = [{"name": "effect", "varies_by": "subject"}]
    s["fixed"]["coefficients"] = {"effect": 5}
    refused(s, "the contrast column 'effect' of factors[1] is also the name of predictors[1]"
            + OWN_FACTOR_ONLY)

    # A contrast named 'item' replaced the item identifiers in R's model_data(), so the
    # formula's (1 | item) grouped by condition.
    s = crossed()
    s["factors"][0]["contrasts"]["item"] = [-0.5, 0.5]
    refused(s, "the contrast column 'item' of factors[1] is also the name of the item unit"
            + OWN_FACTOR_ONLY)

    s = between()
    s["factors"].append({"name": "block", "levels": ["p", "q"],
                         "contrasts": {"group": [-0.5, 0.5]}, "between": "subject"})
    refused(s, "the contrast column 'group' of factors[2] is also the name of factors[1]"
            + OWN_FACTOR_ONLY)

    # A contrast named after its own factor stays allowed: two-level specifications use it.
    s = between()
    s["factors"][0]["contrasts"] = {"group": [-0.5, 0.5]}
    s["fixed"]["coefficients"] = {"group": 5}
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        validate_spec(s)


def test_an_interactions_analysis_column_may_not_take_another_columns_name():
    s = crossed()
    s["predictors"] = [{"name": "freq", "varies_by": "item"},
                       {"name": "cond_freq", "varies_by": "subject"}]
    s["fixed"]["coefficients"]["cond:freq"] = 0.01
    refused(s, "the interaction 'cond:freq' becomes the analysis column 'cond_freq', which is "
               "already the name of predictors[2]")

    s = crossed()
    s["factors"][0]["contrasts"]["cond_cond2"] = [1, -1]
    s["factors"][0]["contrasts"]["cond2"] = [-1, 1]
    s["fixed"]["coefficients"]["cond:cond2"] = 0.01
    refused(s, "the interaction 'cond:cond2' becomes the analysis column 'cond_cond2', which is "
               "already the name of a contrast column of factors[1]")


def test_a_level_listed_twice_is_refused():
    s = between()
    s["factors"][0]["levels"] = ["x", "x"]
    refused(s, "factors[1].levels repeats 'x'; each level needs its own label")


def test_a_factor_both_between_and_within_a_unit_is_refused():
    s = crossed()
    s["factors"][0]["between"] = "subject"
    refused(s, "factors[1] sets both 'vary_within' and 'between'; a factor has to set exactly one "
               "of them")


def test_a_correlation_pairs_two_different_terms_and_each_pair_once():
    s = crossed()
    s["random"]["subject"]["correlations"] = {"cond,cond": 0.25}
    refused(s, "random.subject.correlations key 'cond,cond' pairs 'cond' with itself; a term's "
               "correlation with itself is always 1")

    s = crossed()
    s["random"]["subject"]["correlations"] = {"intercept,cond": 0.9, "cond~intercept": -0.9}
    refused(s, "random.subject.correlations keys 'intercept,cond' and 'cond~intercept' name the "
               "same pair of terms; give each pair once")


def test_random_item_is_refused_in_a_design_without_items():
    s = between()
    s["random"] = {"item": {"intercept_sd": 0.5}}
    refused(s, "random.item describes an item unit the design does not have; add units.item or "
               "remove random.item")


def test_a_blank_name_is_refused_where_it_used_to_reach_the_data():
    # Ported from 7693ab2 on claude/meridian-packages-apps-review-g2a74i. An emptied name built a
    # spec that validated and then produced a column with no name at all, where the R twin
    # stopped inside its simulator with "replacement has length zero".
    s = between()
    for blank in ("", "  "):
        s["factors"][0]["name"] = blank
        refused(s, "factors[1].name must be a non-empty string")
    s = _ex("reading_time_continuous")
    s["predictors"][0]["name"] = ""
    refused(s, "predictors[1].name must be a non-empty string")
    s = between()
    s["response"]["name"] = " "
    refused(s, "'response.name' must be a non-empty string")
    s = _ex("nested_clusters")
    s["random"][""] = s["random"].pop("site")
    refused(s, "a 'random' grouping factor must have a non-empty name")


def test_a_within_factor_that_omits_a_unit_of_the_design_warns_in_either_mode():
    s = crossed()
    full = simulate(s)
    s["factors"][0]["vary_within"] = "subject"
    msg = ("factors[1].vary_within lists 'subject' but not 'item'. pilotr crosses a within factor "
           "with every unit of the design, so it varies within 'item' as well; list every unit, or "
           "make it between 'item' if items carry it. From spec_version 0.4 this is an error.")
    for strict in (True, False):
        with pytest.warns(UserWarning) as rec:
            validate_spec(s, strict=strict)
        assert any(msg in str(w.message) for w in rec)
    # The list never changed the data, which is why it is a warning and not yet an error.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        part = simulate(s)
    assert part.columns == full.columns
    assert part.rows == full.rows

    # A list naming every unit, in any order, and 'subject' alone in a design without items are
    # complete.
    s["factors"][0]["vary_within"] = ["item", "subject"]
    b = between()
    b["factors"][0].pop("between")
    b["factors"][0]["vary_within"] = "subject"
    for spec in (s, b):
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            validate_spec(spec)
