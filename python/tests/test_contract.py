"""The cross-language contract as the documentation states it.

The guarantee was stated as bit-identity for the Gaussian family and for every family with
rounding set, and neither held. AS 241's tail branch passes about 15% of normal draws through
log(), whose last bit each platform's maths library decides, and the validators refuse round for
the discrete families. The replicate_seeds docstring also gave the wrong reason why consecutive
seeds correlate. The document checks read the files around the package and skip outside the
repository. Twinned with test-contract.R.
"""

import os
import re
import statistics
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pytest

from pilotr import RNG, replicate_seeds

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, "..", ".."))


def _text(*parts):
    path = os.path.join(REPO, *parts)
    if not os.path.exists(path):
        pytest.skip(f"{os.path.join(*parts)} not available outside the repository")
    with open(path, encoding="utf-8") as f:
        return re.sub(r"\s+", " ", f.read())  # a sentence may be wrapped anywhere


def test_the_specification_drops_the_guarantees_that_do_not_hold():
    spec = _text("spec", "SPEC.md")
    for claim in (
        # round is refused for bernoulli, poisson and ordinal, and rounding guarantees nothing.
        "also holds for every family when `response.round` is set",
        "either set `round` or stay with `gaussian`",
        # AS 241's tail branch calls log(), so normals agree only where the libraries agree.
        "so the two languages agree to full double precision",
        # A Gamma proposal draws a uniform only when v > 0.
        "each consuming normal–uniform pairs until acceptance",
    ):
        assert claim not in spec, claim
    assert "## Sampling algorithms" in spec


@pytest.mark.parametrize("parts", [
    ("python", "README.md"),
    ("python", "docs", "index.md"),
    ("python", "docs", "getting-started.md"),
    ("python", "docs", "cross-language.md"),
])
def test_the_python_documents_state_the_platform_condition(parts):
    text = _text(*parts)
    for claim in ("bit for bit apart from a documented tolerance",
                  "to full floating-point precision",
                  "agree to full double precision"):
        assert claim not in text, claim
    assert "maths library" in text


def test_replicate_seeds_gives_the_affine_mechanism():
    doc = re.sub(r"\s+", " ", replicate_seeds.__doc__)
    assert "a few steps apart" not in doc
    assert "affine" in doc
    assert "Matsumoto" in doc and "10.1145/1276927.1276928" in doc


def test_power_mixed_says_that_a_simplification_rule_is_simulated_in_r_only():
    # The R twin's power_mixed() and precision_design() simulate a preregistered rule that
    # simplifies the model when a fit is singular. This backend fits one model, so its guide and
    # docstring say where the rule is.
    from pilotr import power_mixed

    for text in (re.sub(r"\s+", " ", power_mixed.__doc__), _text("python", "docs", "power.md")):
        assert "simplification rule" in text
        assert "only in the R package" in text


def _lag1(xs):
    """The lag-1 Pearson correlation, written out because statistics.correlation needs 3.10."""
    a, b = xs[:-1], xs[1:]
    ma, mb = statistics.fmean(a), statistics.fmean(b)
    cov = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    return cov / (sum((x - ma) ** 2 for x in a) * sum((y - mb) ** 2 for y in b)) ** 0.5


def test_consecutive_seeds_step_by_a_constant_as_the_docstring_says():
    first = [RNG(s).uniform() for s in range(1000, 3000)]
    steps = [((b - a + 0.5) % 1.0) - 0.5 for a, b in zip(first[:-1], first[1:])]
    # The step is the same for every pair of consecutive seeds, to within 1e-6.
    assert max(abs(d + 0.009056) for d in steps) < 1e-6
    d = abs(statistics.mean(steps))
    assert _lag1(first) == pytest.approx(1 - 6 * d * (1 - d), abs=0.01)
    assert _lag1(first) > 0.94
    # Seeds drawn from the generator do not step, so their first draws barely correlate.
    seeds = replicate_seeds(90210, 2000)
    assert abs(_lag1([RNG(s).uniform() for s in seeds])) < 0.05
