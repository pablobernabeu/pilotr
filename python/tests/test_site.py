"""The documents and site files the two packages share.

The contributing guide named a CSV comparison with a tolerance as the check that keeps the twins
identical, where CI gates on the harness in tools/parity. The note on advanced designs predated
0.3. The hub's theme toggle stopped working wherever the browser blocks site storage, because
its pages read localStorage outside a try block. These checks read the repository around the
package and skip outside it. Twinned with test-site.R, which also checks the apps, the vignette
on them, the R site's favicons and the browser build's page.
"""

import os
import re

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, "..", ".."))

GATES = (
    "tools/parity/run_r.R",
    "tools/parity/run_py.py",
    "tools/parity/compare.py",
    "tools/parity/validate_cross.py",
    "tools/parity/solve_cross.py",
)


def _read(*parts):
    path = os.path.join(REPO, *parts)
    if not os.path.exists(path):
        pytest.skip(f"{os.path.join(*parts)} not available outside the repository")
    with open(path, encoding="utf-8") as f:
        return f.read()


def _one_line(*parts):
    return re.sub(r"\s+", " ", _read(*parts))  # a sentence may be wrapped anywhere


def _unguarded_script(html):
    """The page's JavaScript with every try block removed: the code that runs unguarded."""
    js = "\n".join(re.findall(r"<script[^>]*>.*?</script>", html, flags=re.S))
    js = re.sub(r"(^|\s)//[^\n]*", r"\1", js)  # line comments
    while True:
        m = re.search(r"\btry\s*\{", js)
        if m is None:
            return js
        depth = 0
        for i in range(m.end() - 1, len(js)):
            if js[i] == "{":
                depth += 1
            elif js[i] == "}":
                depth -= 1
                if depth == 0:
                    break
        js = js[: m.start()] + js[i + 1:]


def test_the_contributing_guide_names_the_parity_gates_ci_runs():
    text = _one_line(".github", "CONTRIBUTING.md")
    for gate in GATES:
        assert gate in text, gate
    assert "exercises that contract and must stay green" not in text
    assert "# confirm cross-language agreement" not in text
    assert "maths library" in text


def test_the_pull_request_template_names_the_parity_harness():
    text = _one_line(".github", "PULL_REQUEST_TEMPLATE.md")
    assert "parity_check.py` confirms the outputs stay bit-identical" not in text
    assert "tools/parity" in text


def test_the_note_on_advanced_designs_describes_specification_0_3():
    text = _one_line("docs", "mixed_models_and_design_analysis.md")
    for term in ('"observation"', '"uniform"', "reliability", "exgaussian"):
        assert term in text, term
    assert "5e-15" not in text


@pytest.mark.parametrize("page", ["landing.html", "404.html"])
def test_the_hub_pages_read_and_write_site_storage_inside_try_blocks(page):
    html = _read(".github", page)
    assert "localStorage" in html
    assert "localStorage" not in _unguarded_script(html)
