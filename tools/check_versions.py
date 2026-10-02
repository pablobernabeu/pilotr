"""Check that the version strings set by hand at release agree with the packages.

The R and Python packages are released separately, so their version numbers can differ
(R 0.3.1 reached CRAN while PyPI stayed at 0.3.0), and the check compares each language
with itself instead of requiring every file to agree:

- R: DESCRIPTION's Version is the only R-side copy. inst/CITATION and the About article
  read it at run time.
- Python: pyproject.toml, pilotr.__version__, the 'extra.version' chip in mkdocs.yml and
  the note in docs/pilotr.bib must all give the same version. The bib file is written by
  the docs build and committed, so a bump without a docs build leaves it behind.
- CITATION.cff, which GitHub and Zenodo read, carries the newer of the two releases.
- The citation sections of the three READMEs carry no version number, since a single
  number cannot be right for both packages. They point to citation("pilotr") and the
  About pages, which give each package's own.

Uses the standard library only, so CI can run it before installing anything.
Run from anywhere: python tools/check_versions.py [--root PATH]. Exits 1 on a mismatch.
"""

import argparse
import os
import re
import sys

README_FILES = ("README.md", "python/README.md", "r/pilotr/README.md")
VERSION_IN_TEXT = re.compile(r"\b\d+\.\d+\.\d+\b")


def _read(root, rel):
    with open(os.path.join(root, *rel.split("/")), encoding="utf-8") as fh:
        return fh.read()


def _first(pattern, text, what):
    match = re.search(pattern, text, flags=re.MULTILINE)
    if match is None:
        raise ValueError(f"no version found in {what}")
    return match.group(1)


def _section(text, opener):
    """The lines from `opener` to the next line that starts a section of the same kind."""
    lines = text.splitlines()
    for i, line in enumerate(lines):
        if line.rstrip() == opener:
            end = next((j for j in range(i + 1, len(lines)) if _same_level(lines[j], opener)),
                       len(lines))
            return "\n".join(lines[i + 1:end])
    return None


def _same_level(line, opener):
    if opener.startswith("["):  # TOML table
        return line.startswith("[")
    if opener.startswith("#"):  # Markdown heading
        level = len(opener) - len(opener.lstrip("#"))
        return re.match(rf"#{{1,{level}}} ", line) is not None
    return bool(line) and not line[0].isspace()  # YAML top-level key


def r_version(root):
    return _first(r"^Version:\s*(\S+)\s*$", _read(root, "r/pilotr/DESCRIPTION"),
                  "r/pilotr/DESCRIPTION")


# Each Python-side copy of the version: the file, the section it sits in (None for the
# whole file) and the pattern that captures it. pyproject.toml comes first because it is
# the one the others are compared with.
PYTHON_COPIES = (
    ("python/pyproject.toml", "[project]", r'^version\s*=\s*"([^"]+)"'),
    ("python/pilotr/__init__.py", None, r'^__version__\s*=\s*"([^"]+)"'),
    ("python/mkdocs.yml", "extra:", r'^\s+version:\s*"?([^"\s]+)"?'),
    ("python/docs/pilotr.bib", None, r"note\s*=\s*\{Python package version ([^}]+)\}"),
)


def _python_copy(root, path, opener, pattern):
    text = _read(root, path)
    if opener is not None:
        text = _section(text, opener) or ""
        return _first(pattern, text, f"{path} {opener}")
    return _first(pattern, text, path)


def python_versions(root):
    """Every Python-side copy of the version, keyed by the file it was read from."""
    return {path: _python_copy(root, path, opener, pattern)
            for path, opener, pattern in PYTHON_COPIES}


def python_version(root):
    path, opener, pattern = PYTHON_COPIES[0]
    return _python_copy(root, path, opener, pattern)


def cff_version(root):
    return _first(r'^version:\s*"?([^"\s]+)"?', _read(root, "CITATION.cff"), "CITATION.cff")


def _release_key(version):
    """The release a version string belongs to, as integers for comparison.

    Only the leading dotted integers count, so a Python pre-release such as 0.4.0.dev1
    compares as 0.4.0. R marks a development version with a fourth component of 9000 or
    more (0.3.1.9000), which still belongs to the 0.3.1 release.
    """
    parts = []
    for piece in version.split("."):
        digits = re.match(r"\d+", piece)
        if digits is None:
            break
        parts.append(int(digits.group(0)))
        if digits.group(0) != piece:
            break
    if len(parts) == 4 and parts[3] >= 9000:
        parts = parts[:3]
    return tuple(parts)


def newer(a, b):
    """The newer of two version strings (the first when they belong to the same release)."""
    return b if _release_key(b) > _release_key(a) else a


def check(root):
    """Every mismatch, as a list of messages; empty when the metadata agree."""
    problems = []

    def read(reader, *args):
        try:
            return reader(root, *args)
        except (OSError, ValueError) as exc:
            problems.append(str(exc))
            return None

    r_ver = read(r_version)
    cff = read(cff_version)
    py_vers = {path: read(_python_copy, path, opener, pattern)
               for path, opener, pattern in PYTHON_COPIES}

    py_ver = py_vers[PYTHON_COPIES[0][0]]
    if py_ver is not None:
        for path, found in py_vers.items():
            if found is not None and found != py_ver:
                problems.append(f"{path} gives Python version {found}, "
                                f"but python/pyproject.toml gives {py_ver}")

    if None not in (r_ver, py_ver, cff):
        want = newer(r_ver, py_ver)
        if _release_key(cff) != _release_key(want):
            problems.append(
                f"CITATION.cff gives version {cff}, but the newer release is {want} "
                f"(R {r_ver}, Python {py_ver})"
            )

    for path in README_FILES:
        try:
            section = _section(_read(root, path), "## Citation")
        except OSError as exc:
            problems.append(str(exc))
            continue
        if section is None:
            problems.append(f"{path} has no '## Citation' section")
            continue
        found = VERSION_IN_TEXT.findall(section)
        if found:
            problems.append(
                f"{path} gives version {', '.join(found)} in its citation. The R and Python "
                "packages are released separately, so the citation should point to "
                "citation(\"pilotr\") and the About pages for the version"
            )
    return problems


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--root",
        default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        help="repository root (default: the parent of tools/)",
    )
    args = parser.parse_args(argv)
    problems = check(args.root)
    if problems:
        for problem in problems:
            print(f"check_versions: {problem}", file=sys.stderr)
        return 1
    print(
        f"Release metadata agree: R {r_version(args.root)}, Python "
        f"{python_version(args.root)}, CITATION.cff {cff_version(args.root)}."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
