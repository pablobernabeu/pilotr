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


def _repo_examples_dir():
    here = os.path.dirname(os.path.abspath(__file__))
    return os.path.normpath(os.path.join(here, "..", "..", "spec", "examples"))


def test_packaged_examples_are_byte_identical_to_the_repository_copies():
    """The specifications under spec/examples/ are canonical; both packages carry
    mirrors so an installed copy can reach them. Nothing enforced the mirror, and
    the load-and-simulate test above cannot: a stale packaged copy still loads and
    simulates perfectly well, it simply describes a different design from the one
    the repository documents. Twinned with test-examples.R. Skips when the package
    is tested in isolation from the repository."""
    canonical = _repo_examples_dir()
    if not os.path.isdir(canonical):
        pytest.skip("repository specifications not available")
    expected = sorted(f for f in os.listdir(canonical) if f.endswith(".json"))
    assert expected, canonical
    packaged = sorted(name + ".json" for name in pilotr_example())
    assert packaged == expected, "the packaged mirror and spec/examples/ differ in membership"
    for name in expected:
        with open(os.path.join(canonical, name), "rb") as fh:
            want = fh.read()
        with open(pilotr_example(name), "rb") as fh:
            got = fh.read()
        assert got == want, f"packaged {name} differs from spec/examples/{name}"
