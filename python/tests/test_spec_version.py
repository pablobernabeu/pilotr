"""The specification version this release reads, and the one a specification needs.

The R twin's tests/testthat/test-spec-version.R runs the same cases and expects the same text, and
also checks that build_spec() and spec_from_model(), which the Python twin lacks, declare the
version the specification needs.
"""
import os, sys, warnings

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pytest

from pilotr import SPEC_VERSION, load_spec, validate_spec
from pilotr.validate import _features_0_4, _required_spec_version

SPEC = os.path.join(os.path.dirname(__file__), "..", "..", "spec", "examples")


def _ex(name):
    return load_spec(os.path.join(SPEC, name + ".json"))


def _factor_b(unit="subject"):
    return {"name": "B", "levels": ["b1", "b2"], "contrasts": {"b": [-0.5, 0.5]}, "between": unit}


def silent(spec):
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        validate_spec(spec)


def test_this_release_reads_specification_version_0_4_and_refuses_a_later_one():
    assert SPEC_VERSION == "0.4"
    s = _ex("between_2group_gaussian")
    s["spec_version"] = "0.4"
    silent(s)
    s["spec_version"] = "0.5"
    with pytest.raises(ValueError) as e:
        validate_spec(s)
    assert ("this specification declares spec_version 0.5, which is newer than the 0.4 this "
            "version of pilotr understands; please upgrade pilotr") in str(e.value)


def test_a_specification_needs_0_4_for_several_between_factors_on_one_unit():
    # The other 0.4 feature, a counterbalanced factor, is tested in test_counterbalance.py.
    s = _ex("between_2group_gaussian")
    assert _features_0_4(s) == []
    assert _required_spec_version(s) == "0.3"
    # A factor between subjects beside one between items needs nothing new.
    s["units"]["item"] = {"n": 4}
    s["factors"].append(_factor_b("item"))
    assert _required_spec_version(s) == "0.3"
    s["factors"][1]["between"] = "subject"
    assert _features_0_4(s) == ["several between factors on one unit"]
    assert _required_spec_version(s) == "0.4"


def test_the_version_a_specification_needs_is_never_below_0_3():
    # A specification without a 0.3 feature still needs 0.3, the version the R twin's
    # build_spec() has written since 0.3.0. Releases before 0.3 ignore the field.
    s = _ex("between_2group_gaussian")
    s.pop("spec_version", None)
    assert _required_spec_version(s) == "0.3"
    # Malformed parts are left to the validator, as for the 0.3 features.
    assert _required_spec_version({"factors": "none"}) == "0.3"
    assert _required_spec_version({"factors": [5, {"name": "A"}]}) == "0.3"
