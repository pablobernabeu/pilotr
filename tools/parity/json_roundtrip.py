"""Cross-language check that the numbers R writes with spec_json() read back exactly in Python.

A specification file is the one artefact the two packages share, so each number in it has to
name the same double in both. spec_json() writes a number at the shortest width, 15, 16 or 17
significant digits, that a correctly rounded reader brings back to the double R holds. It checks
each candidate with jsonlite, whose strtod() rounds correctly, as Python's json does. A check with
as.numeric() alone, which old_rule() in the R driver reproduces, is not enough: on a build of R
with long doubles, such as R on x86 Windows and Linux, as.numeric() rounds a decimal twice, first
to long double and then to double, so it passes some numbers at a width that jsonlite and Python
read one unit in the last place away, and a specification saved and reloaded would simulate
different data.

The R side sources r/pilotr/R and writes spec_json() for 20,000 doubles drawn with a fixed
seed, after four fixed ones, together with each double's exact hexadecimal form from
sprintf("%a"). Half of the random doubles come from runif(), with its 32-bit resolution, and
half from rnorm() scaled over 31 orders of magnitude. Python reads the document with json and
compares every number with float.fromhex(). The R side also counts how many of the doubles that
as.numeric() check would have written in a form jsonlite reads differently. On a build with long
doubles, that count shows the batch is large enough to catch the fault. On a build without them,
it is 0.

Usage: python tools/parity/json_roundtrip.py
Exit status is 0 when every number reads back exactly.
"""

from __future__ import annotations

import json
import os
import re
import subprocess
import sys
import tempfile

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

# The interpreter on PATH by default, so that the script runs wherever R is installed, CI
# included. PILOTR_RSCRIPT overrides it on a machine whose Rscript lives elsewhere.
RSCRIPT = os.environ.get("PILOTR_RSCRIPT", "Rscript")

R_DRIVER = r"""
args <- commandArgs(trailingOnly = TRUE)
for (f in sort(list.files(args[1], pattern = "[.]R$", full.names = TRUE))) source(f)
set.seed(20261013)
n <- 10000
vals <- c(-0x1.0efdc5258p+3, 0x1.e2b06cc2p+2, 0x1.1fe271f9p+3, 1 / 3,
          runif(n, -10, 10), rnorm(n) * 10^sample(-15:15, n, replace = TRUE))
con <- file(args[2], open = "wb")
writeLines(spec_json(list(values = vals)), con, sep = "\n")
close(con)
con <- file(args[3], open = "wb")
writeLines(sprintf("%a", vals), con, sep = "\n")
close(con)
# The 0.3.1 rule: the first of 15 and 16 digits that as.numeric() reads back, else 17.
old_rule <- function(z) {
  for (d in 15:16) {
    s <- sprintf(paste0("%.", d, "g"), z)
    if (as.numeric(s) == z) return(s)
  }
  sprintf("%.17g", z)
}
old_wrong <- sum(vapply(vals, function(z) jsonlite::fromJSON(old_rule(z)) != z, logical(1)))
cat(jsonlite::toJSON(list(r = R.version.string, long_double = capabilities("long.double"),
                          old_wrong = old_wrong), auto_unbox = TRUE))
"""


def _significant_digits(text):
    """The significant digits a number is written with, as spec_json() counts its widths."""
    mantissa = re.split("[eE]", text.lstrip("-"))[0]
    digits = mantissa.replace(".", "").lstrip("0")
    return len(digits.rstrip("0")) if "." in mantissa else len(digits)


def _run_r(json_path, hex_path):
    driver = os.path.join(os.path.dirname(json_path), "driver.R")
    with open(driver, "w") as f:
        f.write(R_DRIVER)
    src = os.path.join(ROOT, "r", "pilotr", "R")
    proc = subprocess.run([RSCRIPT, "--vanilla", driver, src, json_path, hex_path],
                          capture_output=True, text=True)
    if proc.returncode != 0:
        raise SystemExit("the R driver failed:\n" + proc.stderr)
    return json.loads(proc.stdout)


def main() -> int:
    with tempfile.TemporaryDirectory() as tmp:
        json_path = os.path.join(tmp, "values.json")
        hex_path = os.path.join(tmp, "hex.txt")
        meta = _run_r(json_path, hex_path)
        with open(json_path, encoding="utf-8") as f:
            text = f.read()
        with open(hex_path, encoding="utf-8") as f:
            held = [float.fromhex(h) for h in f.read().split()]
    read = json.loads(text)["values"]
    # The same document again, keeping each number's text, to report the widths written.
    written = json.loads(text, parse_float=str, parse_int=str)["values"]
    if len(read) != len(held):
        print("R wrote %d numbers and %d hexadecimal forms" % (len(read), len(held)))
        return 1

    widths = [_significant_digits(t) for t in written]
    misread = [i for i, (r, h) in enumerate(zip(read, held)) if float(r) != h]
    print("%d numbers written by spec_json() in %s (long double: %s)"
          % (len(read), meta["r"], meta["long_double"]))
    print("  widths: %d at 15 significant digits or fewer, %d at 16, %d at 17"
          % (sum(w <= 15 for w in widths), sum(w == 16 for w in widths),
             sum(w == 17 for w in widths)))
    print("  checked with as.numeric() alone, as in 0.3.1, %d would have been written in a form"
          " jsonlite reads differently" % meta["old_wrong"])
    for i in misread[:10]:
        print("  MISREAD %s: written %s, read %s, held %s"
              % (i, written[i], float(read[i]).hex(), held[i].hex()))
    print("%d of %d numbers read back differently in Python" % (len(misread), len(read)))
    return 1 if misread else 0


if __name__ == "__main__":
    sys.exit(main())
