"""The bundled example specifications load and simulate from the installed package."""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pytest

from pilotr import load_spec, pilotr_example, simulate


def test_pilotr_example_lists_and_resolves_each_spec():
    names = pilotr_example()
    assert isinstance(names, list)
    assert "between_2group_gaussian" in names
    for name in names:
        path = pilotr_example(name)
        assert os.path.exists(path)
        # Every shipped example loads and simulates without error.
        data = simulate(load_spec(path))
        assert len(data) > 0
    # The .json extension is optional.
    assert pilotr_example("between_2group_gaussian") == pilotr_example(
        "between_2group_gaussian.json"
    )


def test_pilotr_example_rejects_unknown_names():
    with pytest.raises(ValueError, match="Unknown example"):
        pilotr_example("no_such_example")
