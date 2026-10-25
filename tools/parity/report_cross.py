"""Cross-language agreement check for power_report().

power_report() reads a result into a record of counts and strings, and then writes the record as
Markdown. The two twins must write the same bytes for the same record: a justification quoted in
a preregistration should not depend on the language that produced it. Every rate in the text, its
Monte Carlo standard error and its Wilson interval, is computed while rendering, from the count of
significant replicates and the count behind it, with the functions of mcse.R and mcse.py, which
return the same doubles for the same counts. Both twins format those doubles with the same format
strings, through sprintf() in R and the % operator in Python, and both round correctly, so the
text is required to match exactly, with no tolerance.

The records are stored in report_cases.json, not produced by live analyses, so that the check
tests the renderer alone. A live result depends on its fitter, lme4 in R and statsmodels in
Python, and its counts would differ between the languages for reasons that have nothing to do with
the report. The cases cover every section the two twins share and these variants of a line: a
ladder of models, a joint test, a one-sided test, an effect with no true value, an effect no
replicate estimated, a Type M with too few significant replicates for a standard error, a solve
over an axis other than sample size, a solve over the precision analysis's rate, a null
condition, a missing specification, a negative seed and a justification with a character outside
ASCII. A precision_design() record, which the R twin alone renders, is not checked here.

Usage: python tools/parity/report_cross.py
Exit status is 0 when the two twins write the same text for every case.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "python"))

from pilotr.report import _render  # noqa: E402

# The interpreter on PATH by default, so the script runs anywhere R is installed (CI
# included); PILOTR_RSCRIPT overrides it for a machine whose Rscript lives elsewhere.
RSCRIPT = os.environ.get("PILOTR_RSCRIPT", "Rscript")

CASES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "report_cases.json")

# Each case's text is written to a file of its own as UTF-8 bytes, so that neither the console's
# encoding nor its line endings come between R and the comparison.
R_DRIVER = r"""
args <- commandArgs(trailingOnly = TRUE)
for (f in sort(list.files(args[1], pattern = "[.]R$", full.names = TRUE))) source(f)
cases <- jsonlite::fromJSON(args[2], simplifyVector = FALSE)
for (nm in names(cases)) {
  con <- file(file.path(args[3], paste0(nm, ".md")), open = "wb")
  writeBin(charToRaw(enc2utf8(.report_render(cases[[nm]]))), con)
  close(con)
}
"""


def _run_r(out_dir):
    """Render every case in R, into one file per case in `out_dir`."""
    driver = os.path.join(out_dir, "driver.R")
    with open(driver, "w", encoding="utf-8") as f:
        f.write(R_DRIVER)
    src = os.path.join(ROOT, "r", "pilotr", "R")
    proc = subprocess.run([RSCRIPT, "--vanilla", driver, src, CASES, out_dir],
                          capture_output=True, text=True, encoding="utf-8", errors="replace")
    if proc.returncode != 0:
        raise SystemExit("the R driver failed:\n" + proc.stderr)


def _first_difference(a, b):
    """The first line at which two texts differ, as a short description."""
    la, lb = a.split("\n"), b.split("\n")
    for i, (x, y) in enumerate(zip(la, lb)):
        if x != y:
            return "line %d: R=%r Py=%r" % (i + 1, x, y)
    return "R has %d lines, Python %d" % (len(la), len(lb))


def main() -> int:
    with open(CASES, encoding="utf-8") as f:
        cases = json.load(f)
    failures = 0
    print("%-28s %-10s %s" % ("case", "parity", "detail"))
    print("-" * 92)
    with tempfile.TemporaryDirectory() as tmp:
        _run_r(tmp)
        for name, record in cases.items():
            with open(os.path.join(tmp, name + ".md"), "rb") as f:
                r_bytes = f.read()
            py_bytes = _render(record).encode("utf-8")
            if r_bytes == py_bytes:
                print("%-28s %-10s %d bytes identical" % (name, "ok", len(r_bytes)))
                continue
            failures += 1
            print("%-28s %-10s %s" % (name, "DIFFERS",
                                      _first_difference(r_bytes.decode("utf-8", "replace"),
                                                        py_bytes.decode("utf-8"))))
    print("\n%d of %d cases failed" % (failures, len(cases)))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
