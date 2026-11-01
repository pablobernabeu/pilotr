"""Cross-language check of the analysis model: model_formula() and model_data().

The Python power_mixed() fits the model the specification implies, which ports of the R twin's
model_formula() and model_data() build (python/pilotr/autoformula.py). This script holds the ports
to the R functions on every shipped example and parity case, and on three cases written below that
exercise what those do not: a contrast named after its own three-level factor, the `correlated`
flag in both directions, and a three-way interaction with a random slope on an interaction.

For each specification, R writes its formula under each of the three structures, the data set
simulate_design() draws and the data model_data() builds from it. Python reads R's simulated data,
builds its own modelling data from it and writes its own formulas. The formulas must be identical,
character for character, and where R refuses, Python must refuse in the same words, as both do for
a specification with neither a coefficient nor a random effect. The modelling data must have the
same columns in the same order, and every cell must agree bit for bit, compared through its text at
17 significant digits, which reads back as the same double. The one exception is the .y column of
the two lognormal families. It comes from the C library's log(), which can differ by an ulp between
R and Python on Windows, as tools/parity/tolerance.json records, so it is allowed the 8 ulp that
file gives the unrounded dumps of the lognormal, shifted-lognormal and ex-Gaussian cases. The
contrast and product columns are copies and products of doubles, which IEEE 754 rounds correctly,
so they are allowed nothing.

Python is given R's simulated data, so that the simulator stays out of the comparison. The dump
harness beside this file checks the simulator, and a normal draw that differed by an ulp between
the twins would otherwise show here as a product column that differs.

Usage: python tools/parity/formula_cross.py
Exit status is 0 when the two twins agree on every case.
"""

from __future__ import annotations

import json
import math
import os
import struct
import subprocess
import sys
import tempfile

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "python"))

from pilotr import load_spec, model_data, model_formula  # noqa: E402

# The interpreter on PATH by default, so the script runs anywhere R is installed, CI included.
# PILOTR_RSCRIPT overrides it for a machine whose Rscript lives elsewhere.
RSCRIPT = os.environ.get("PILOTR_RSCRIPT", "Rscript")

STRUCTURES = ("maximal", "zcp", "intercepts")

# The allowance for .y under a lognormal family, in units in the last place, that tolerance.json
# gives the unrounded dumps of crossed_mixed_rt, reading_time_continuous and stress_exgaussian.
LOG_ULP = 8

# Cases written here, for what no example or parity case exercises.
EXTRA = [
    ("own_name_three_levels", {
        "name": "own_name_three_levels", "seed": 12,
        "units": {"subject": {"n": 6}, "item": {"n": 4}},
        "factors": [{"name": "cond", "levels": ["a", "b", "c"],
                     "contrasts": {"cond": [-1, 0, 1], "cond2": [1, -2, 1]},
                     "vary_within": ["subject", "item"]}],
        "fixed": {"intercept": 0, "coefficients": {"cond": 0.3, "cond2": 0.1}},
        "random": {"subject": {"intercept_sd": 0.5, "slopes": {"cond": 0.2, "cond2": 0.1}},
                   "item": {"intercept_sd": 0.3}},
        "response": {"family": "gaussian", "name": "y", "sigma": 1},
    }),
    ("correlated_flags", {
        "spec_version": "0.3", "name": "correlated_flags", "seed": 13,
        "units": {"subject": {"n": 8}, "item": {"n": 6}},
        "factors": [{"name": "condition", "levels": ["a", "b"],
                     "contrasts": {"cond": [-0.5, 0.5]}, "vary_within": ["subject", "item"]}],
        "fixed": {"intercept": 1, "coefficients": {"cond": 0.4}},
        "random": {"subject": {"intercept_sd": 0.4, "slopes": {"cond": 0.2},
                               "correlated": False},
                   "item": {"intercept_sd": 0.3, "slopes": {"cond": 0.1}, "correlated": True}},
        "response": {"family": "shifted_lognormal", "name": "rt", "sigma": 0.3, "shift": 150,
                     "round": 2},
    }),
    ("three_way_interaction", {
        "spec_version": "0.3", "name": "three_way_interaction", "seed": 14,
        "units": {"subject": {"n": 10}, "item": {"n": 6}},
        "predictors": [{"name": "freq", "varies_by": "item", "mean": 0, "sd": 1}],
        "factors": [
            {"name": "A", "levels": ["a1", "a2"], "contrasts": {"a": [-0.5, 0.5]},
             "vary_within": ["subject", "item"]},
            {"name": "G", "levels": ["g1", "g2"], "contrasts": {"g": [-0.5, 0.5]},
             "between": "subject"}],
        "fixed": {"intercept": -1, "coefficients": {"a": 0.2, "g": 0.1, "freq": 0.15,
                                                    "a:g:freq": 0.05}},
        "random": {"subject": {"intercept_sd": 0.2, "slopes": {"a:freq": 0.05}},
                   "item": {"intercept_sd": 0.1, "slopes": {"a:g": 0.04}}},
        "response": {"family": "lognormal", "name": "time", "sigma": 0.25},
    }),
]

R_DRIVER = r"""
args <- commandArgs(trailingOnly = TRUE)
for (f in sort(list.files(args[1], pattern = "[.]R$", full.names = TRUE))) source(f)
paths <- jsonlite::fromJSON(args[2])
# Every cell as text: a label as it is, a number at 17 significant digits, which reads back as the
# same double, and a missing or undefined number as NA.
column <- function(v) {
  if (is.character(v) || is.factor(v)) return(as.character(v))
  s <- sprintf("%.17g", as.numeric(v))
  s[is.na(v)] <- "NA"
  s
}
table <- function(d) list(names = names(d), columns = unname(lapply(d, column)))
out <- lapply(paths, function(p) {
  spec <- load_spec(p)
  formulas <- lapply(c(maximal = "maximal", zcp = "zcp", intercepts = "intercepts"), function(s)
    tryCatch(paste(trimws(deparse(model_formula(spec, s), width.cutoff = 500L)), collapse = " "),
             error = function(e) list(error = conditionMessage(e))))
  d <- simulate_design(spec)
  list(formulas = formulas, data = table(d), model = table(model_data(spec, d)))
})
cat(jsonlite::toJSON(out, auto_unbox = TRUE, null = "null"))
"""


def _cell(v) -> str:
    """A cell as the R driver writes it: a label as it is, a number at 17 significant digits."""
    if isinstance(v, str):
        return v
    v = float(v)
    if math.isnan(v):
        return "NA"
    if math.isinf(v):
        return "Inf" if v > 0 else "-Inf"
    return "%.17g" % v


_SPECIAL = {"NA": math.nan, "NaN": math.nan, "Inf": math.inf, "-Inf": -math.inf}


def _number(s: str) -> float:
    """The double a cell's text stands for, R's "NA", "Inf" and "-Inf" included."""
    return _SPECIAL[s] if s in _SPECIAL else float(s)


def _ulps(a: float, b: float) -> float:
    """Signed-magnitude ordinal distance between two doubles, as compare.py measures it."""
    if a == b or (math.isnan(a) and math.isnan(b)):
        return 0.0
    if math.isnan(a) or math.isnan(b) or math.isinf(a) or math.isinf(b):
        return math.inf

    def ordinal(x: float) -> int:
        (n,) = struct.unpack("<q", struct.pack("<d", x))
        return n if n >= 0 else -(n & 0x7FFFFFFFFFFFFFFF) - (1 << 63)

    return float(abs(ordinal(a) - ordinal(b)))


def _cases():
    """Every specification to check, as (name, path, spec) with a dict spec for the extras."""
    out = []
    for folder in (os.path.join(ROOT, "spec", "examples"),
                   os.path.join(ROOT, "tools", "parity", "cases")):
        for name in sorted(os.listdir(folder)):
            if name.endswith(".json"):
                path = os.path.join(folder, name)
                out.append((name[:-len(".json")], path, load_spec(path)))
    return out


def _run_r(paths):
    """Run the R driver over the specification files and return its output, case by case."""
    with tempfile.TemporaryDirectory() as tmp:
        paths_file = os.path.join(tmp, "paths.json")
        driver_path = os.path.join(tmp, "driver.R")
        with open(paths_file, "w", encoding="utf-8") as f:
            json.dump(paths, f)
        with open(driver_path, "w", encoding="utf-8") as f:
            f.write(R_DRIVER)
        proc = subprocess.run([RSCRIPT, "--vanilla", driver_path,
                               os.path.join(ROOT, "r", "pilotr", "R"), paths_file],
                              capture_output=True, text=True, encoding="utf-8", errors="replace")
    if proc.returncode != 0:
        raise SystemExit("the R driver failed:\n" + proc.stderr)
    return json.loads(proc.stdout)


def _check(spec, r_out):
    """Python's formulas and modelling data against R's.

    Returns the differences as messages and the largest gap in ulp among the .y cells of a
    lognormal family that differ within the allowance.
    """
    import pandas as pd

    problems = []
    for s in STRUCTURES:
        expected = r_out["formulas"][s]
        try:
            got = model_formula(spec, s)
        except ValueError as e:
            got = {"error": str(e)}
        # A formula or a refusal, {"error": message}, each compared character for character.
        if got != expected:
            problems.append("%s: R=%r Py=%r" % (s, expected, got))

    # R's simulated data, with the factors' labels as text and every other column as numbers.
    labels = {f["name"] for f in spec.get("factors") or []}
    data = r_out["data"]
    frame = pd.DataFrame({name: (col if name in labels else [_number(c) for c in col])
                          for name, col in zip(data["names"], data["columns"])},
                         columns=data["names"])
    got = model_data(spec, frame)
    expected = r_out["model"]
    if list(got.columns) != expected["names"]:
        return problems + ["columns: R=%s Py=%s" % (expected["names"], list(got.columns))], 0.0
    logged = spec["response"]["family"] in ("lognormal", "shifted_lognormal")
    worst = 0.0
    for name, r_col in zip(expected["names"], expected["columns"]):
        py_col = [_cell(v) for v in got[name]]
        for i, (a, b) in enumerate(zip(r_col, py_col)):
            if a == b:
                continue
            gap = _ulps(_number(a), _number(b)) if name == ".y" and logged else math.inf
            if gap <= LOG_ULP:
                worst = max(worst, gap)
                continue
            problems.append("%s row %d: R=%s Py=%s" % (name, i + 1, a, b))
            break
    return problems, worst


def main() -> int:
    cases = _cases()
    with tempfile.TemporaryDirectory() as tmp:
        for name, spec in EXTRA:
            path = os.path.join(tmp, name + ".json")
            with open(path, "w", encoding="utf-8") as f:
                json.dump(spec, f)
            cases.append((name, path, spec))
        r_out = _run_r([path for _, path, _ in cases])

    failures = 0
    print("%-36s %-8s %s" % ("case", "parity", "detail"))
    print("-" * 92)
    for (name, _path, spec), rr in zip(cases, r_out):
        problems, worst = _check(spec, rr)
        if problems:
            failures += 1
            print("%-36s %-8s %s" % (name, "DIFFERS", problems[0]))
            for p in problems[1:5]:
                print("%46s%s" % ("", p))
        else:
            refused = all(isinstance(rr["formulas"][s], dict) for s in STRUCTURES)
            detail = ("formulas refused alike; " if refused else "formulas identical; ")
            detail += ("columns bit for bit" if not worst
                       else "columns bit for bit, .y within %d ulp" % worst)
            print("%-36s %-8s %s" % (name, "ok", detail))
    print("\n%d of %d cases failed" % (failures, len(cases)))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
