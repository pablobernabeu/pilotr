# About

## Citing pilotr

If you use pilotr in published work, please cite it:

````python exec="1"
# The citation, the BibTeX entry and the downloadable .bib file are all built here
# from one string, using the version the package itself reports, so none of them can
# drift away from the code. This reaches none of the copies that sit where no code can
# run, and those still have to be bumped by hand when a release is cut. They are the
# 'extra.version' chip in mkdocs.yml, the 'version' field in CITATION.cff, and the
# citation block in the repository README and in the Python package's README.
import pilotr

version = pilotr.__version__

bibtex = (
    "@Manual{pilotr,\n"
    "  title  = {{pilotr}: Simulate experimental and behavioural data from a portable design specification},\n"
    "  author = {Pablo Bernabeu},\n"
    "  year   = {2026},\n"
    f"  note   = {{R and Python package version {version}}},\n"
    "  doi    = {10.5281/zenodo.21266313},\n"
    "  url    = {https://doi.org/10.5281/zenodo.21266313},\n"
    "}\n"
)

# docs/pilotr.bib is the file the download link below serves, and it is written from
# the same string, so nobody has to maintain it by hand. MkDocs renders the pages before it
# copies the static files, so the copy that reaches the built site is this one. The
# file is only rewritten when its contents have actually changed, because 'mkdocs
# serve' watches the docs directory and an unconditional write would set off a fresh
# rebuild after every pass. A build that rewrites it leaves a change to commit. The
# read is deliberately unguarded, since MkDocs lists the static files before any page
# runs, so a file created here would never be copied and the link would quietly 404.
with open("docs/pilotr.bib", encoding="utf-8") as f:
    on_disk = f.read()
if on_disk != bibtex:
    with open("docs/pilotr.bib", "w", encoding="utf-8") as f:
        f.write(bibtex)

print("> Bernabeu, P. (2026). *pilotr: Simulate experimental and behavioural data from a portable")
print(f"> design specification* (R and Python package version {version}).")
print("> https://doi.org/10.5281/zenodo.21266313")
print()
print("```bibtex")
print(bibtex, end="")
print("```")

# The download link is left as ordinary page content below, outside this block. MkDocs
# rewrites relative links only in the page's own Markdown, so a link printed from
# executed code would keep the path as written, resolving against about/ when it should
# resolve against the site root. It names no version, so it cannot drift.
````

[Download the BibTeX entry](pilotr.bib){ download="pilotr.bib" }

In R, `citation("pilotr")` returns the same reference.

## The developer

pilotr is developed by [Pablo Bernabeu](https://pablobernabeu.github.io), a researcher in
the Department of Education at the University of Oxford. His work spans cognitive psychology,
neuroscience, linguistics, education and research methods, with hands-on experience of
behavioural experiments, EEG, corpus analysis, computational modelling and statistics. He
develops open, reproducible research software in R and Python, and is a Fellow of the Software
Sustainability Institute. pilotr and its
[R twin](https://pablobernabeu.github.io/pilotr/r/) share one design specification, keeping a
simulation reproducible across both languages. His
[ORCID record](https://orcid.org/0000-0003-1083-2460) lists his other work.

## Licence

pilotr is released under the MIT licence, reproduced in full on the
[licence page](licence.md). The licence covers the Python and R packages, the no-code app
and the design specification alike.

## Versioning and archival

Each release is tagged on GitHub and archived on Zenodo. The concept DOI,
[10.5281/zenodo.21266313](https://doi.org/10.5281/zenodo.21266313), always resolves to the
latest archived version, so a citation stays current without naming a version. The
[changelog](changelog.md) records what changed in each release.

## Contributing and support

Bugs and feature requests are best raised on the
[GitHub issues page](https://github.com/pablobernabeu/pilotr/issues). The
[contributing guide](https://github.com/pablobernabeu/pilotr/blob/main/.github/CONTRIBUTING.md)
describes the development setup and the conventions the repository follows, including the
cross-language random-number contract that any change to the generative core must keep
intact.
