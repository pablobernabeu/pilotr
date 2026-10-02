"""Release metadata: citations name one package's version, and every version string agrees.

The R and Python packages are released separately, so their version numbers can differ
(R 0.3.1 reached CRAN while PyPI stayed at 0.3.0). A citation that reads "R and Python
package version" names a release of one of them that may not exist. These checks read the
repository around the package and skip when the package is tested in isolation from it.
Twinned with test-citation.R, which checks the R side.
"""

import importlib.util
import os
import shutil
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pytest

import pilotr

HERE = os.path.dirname(os.path.abspath(__file__))
PYTHON_DIR = os.path.normpath(os.path.join(HERE, ".."))
REPO = os.path.normpath(os.path.join(HERE, "..", ".."))
CHECKER = os.path.join(REPO, "tools", "check_versions.py")


def _repo_file(*parts):
    path = os.path.join(REPO, *parts)
    if not os.path.exists(path):
        pytest.skip(f"{os.path.join(*parts)} not available outside the repository")
    return path


def _read(path):
    with open(path, encoding="utf-8") as fh:
        return fh.read()


def _checker():
    # Skip only outside the repository; inside it, a missing checker is a failure.
    _repo_file("r", "pilotr", "DESCRIPTION")
    spec = importlib.util.spec_from_file_location("check_versions", CHECKER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_about_page_cites_the_python_package_version_alone():
    text = _read(_repo_file("python", "docs", "about.md"))
    assert "R and Python package version" not in text
    assert "Python package version {version}" in text


def test_committed_bibtex_names_the_installed_python_version():
    text = _read(_repo_file("python", "docs", "pilotr.bib"))
    assert f"note   = {{Python package version {pilotr.__version__}}}," in text


def test_repository_versions_agree():
    checker = _checker()
    assert checker.check(REPO) == []
    assert checker.main(["--root", REPO]) == 0


# The files the checker reads, copied into a scratch tree so each test can break one.
CHECKED = [
    ("r", "pilotr", "DESCRIPTION"),
    ("python", "pyproject.toml"),
    ("python", "pilotr", "__init__.py"),
    ("python", "mkdocs.yml"),
    ("python", "docs", "pilotr.bib"),
    ("CITATION.cff",),
    ("README.md",),
    ("python", "README.md"),
    ("r", "pilotr", "README.md"),
]


@pytest.fixture
def tree(tmp_path):
    _checker()
    for parts in CHECKED:
        target = tmp_path.joinpath(*parts)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(_repo_file(*parts), target)
    return tmp_path


def _replace(path, old, new):
    text = path.read_text(encoding="utf-8")
    assert old in text, (path, old)
    path.write_text(text.replace(old, new, 1), encoding="utf-8")


def test_checker_passes_an_untouched_copy(tree):
    assert _checker().check(str(tree)) == []


def test_checker_reports_python_versions_that_disagree(tree):
    _replace(tree / "python" / "pyproject.toml", f'version = "{pilotr.__version__}"',
             'version = "9.9.9"')
    problems = _checker().check(str(tree))
    assert any("python/pyproject.toml" in p and "9.9.9" in p for p in problems), problems
    assert _checker().main(["--root", str(tree)]) == 1


def test_checker_wants_citation_cff_at_the_newer_release(tree):
    checker = _checker()
    cff = tree / "CITATION.cff"
    newer = checker.newer(checker.r_version(str(tree)), checker.python_version(str(tree)))
    _replace(cff, f'version: "{newer}"', 'version: "0.0.1"')
    problems = checker.check(str(tree))
    assert any("CITATION.cff" in p and "0.0.1" in p and newer in p for p in problems), problems


def test_checker_refuses_a_version_in_a_readme_citation(tree):
    readme = tree / "python" / "README.md"
    _replace(readme, "## Citation\n", "## Citation\n\nPython package version 0.3.0.\n")
    problems = _checker().check(str(tree))
    assert any("python/README.md" in p and "0.3.0" in p for p in problems), problems


def test_checker_reports_a_stale_bibtex_file(tree):
    _replace(tree / "python" / "docs" / "pilotr.bib", pilotr.__version__, "0.0.1")
    problems = _checker().check(str(tree))
    assert any("python/docs/pilotr.bib" in p for p in problems), problems


def test_newer_compares_numerically():
    checker = _checker()
    assert checker.newer("0.3.1", "0.3.0") == "0.3.1"
    assert checker.newer("0.9.0", "0.10.0") == "0.10.0"
    assert checker.newer("1.0.0", "1.0.0") == "1.0.0"
