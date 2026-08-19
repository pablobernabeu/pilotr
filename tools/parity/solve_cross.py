"""Cross-language agreement check for solve_curve().

The dump harness beside this file compares simulated data, which is where the bit-for-bit
guarantee lives. solve_curve() is not part of that guarantee and cannot be: its fit calls
exp() and the normal distribution function at every iteration, for neither of which IEEE-754
fixes a rounding, and R for Windows (MinGW-w64) and CPython (UCRT) do not share a maths
library. A golden hash over its output would pin an
accident of one platform's libm, which is exactly what tolerance.json's classification
exists to prevent, so this feature is deliberately kept out of golden.json and out of
tolerance.json's per-case dump list, and is gated here instead.

What is required is that the two engines agree far more closely than any statistically
meaningful difference. The allowance below is a relative difference, not an ulp count,
because the quantity compared is the end of an iterative fit rather than a single arithmetic
result: a hundred IRLS steps, each rounding a few libm calls, accumulate more than one ulp
even though nothing has gone wrong. TOLERANCE is set at 1e-9, which on a solved sample size
of 220 subjects is two ten-millionths of a subject. The largest disagreement measured across
the cases below is 5.0e-15, so the allowance sits some five orders of magnitude above it,
which leaves headroom for a third platform's libm without admitting anything a reader could
act on. A logic error would show as a differing refusal message, a differing column, or a
discrepancy many orders of magnitude larger, all of which still fail.

The curves are fixed tables written into this file rather than simulated, so that the check
tests the solver alone and cannot be disturbed by a change to the generative core.

Usage: python tools/parity/solve_cross.py
Exit status is 0 when the two engines agree on every case.
"""

from __future__ import annotations

import json
import math
import os
import subprocess
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "python"))

from pilotr.solve import target_n  # noqa: E402

# The interpreter on PATH by default, so the script runs anywhere R is installed (CI
# included); PILOTR_RSCRIPT overrides it for a machine whose Rscript lives elsewhere.
RSCRIPT = os.environ.get("PILOTR_RSCRIPT", "Rscript")

TOLERANCE = 1e-9

# One shallow curve, one steep one, one over a wide range, and one whose rate column is the
# precision analysis's rather than the power analysis's, so both auto-detected names are
# exercised. The counts vary across points, which is what makes the weighting visible.
CASES = [
    ("power_shallow", {
        "curve": [{"n_subject": 20, "power": 0.21, "n_returned": 200},
                  {"n_subject": 40, "power": 0.38, "n_returned": 200},
                  {"n_subject": 60, "power": 0.52, "n_returned": 190},
                  {"n_subject": 80, "power": 0.63, "n_returned": 200},
                  {"n_subject": 120, "power": 0.78, "n_returned": 180},
                  {"n_subject": 160, "power": 0.87, "n_returned": 200}],
        "target": 0.8, "transform": "sqrt", "level": 0.95}),
    ("power_steep", {
        "curve": [{"n_subject": 10, "power": 0.08, "n_sims": 500},
                  {"n_subject": 20, "power": 0.34, "n_sims": 500},
                  {"n_subject": 30, "power": 0.66, "n_sims": 500},
                  {"n_subject": 40, "power": 0.85, "n_sims": 500},
                  {"n_subject": 50, "power": 0.94, "n_sims": 500}],
        "target": 0.9, "transform": "sqrt", "level": 0.99}),
    ("precision_rope", {
        "curve": [{"n_subject": 15, "p_meaningful": 0.06, "n_returned": 200},
                  {"n_subject": 30, "p_meaningful": 0.14, "n_returned": 200},
                  {"n_subject": 60, "p_meaningful": 0.33, "n_returned": 200},
                  {"n_subject": 100, "p_meaningful": 0.55, "n_returned": 200},
                  {"n_subject": 140, "p_meaningful": 0.71, "n_returned": 199},
                  {"n_subject": 180, "p_meaningful": 0.83, "n_returned": 200},
                  {"n_subject": 220, "p_meaningful": 0.90, "n_returned": 200},
                  {"n_subject": 260, "p_meaningful": 0.94, "n_returned": 200}],
        "target": 0.9, "transform": "sqrt", "level": 0.95}),
    ("effect_axis_identity", {
        "curve": [{"effect_size": 0.10, "power": 0.11, "n_sims": 400},
                  {"effect_size": 0.20, "power": 0.31, "n_sims": 400},
                  {"effect_size": 0.30, "power": 0.58, "n_sims": 400},
                  {"effect_size": 0.40, "power": 0.79, "n_sims": 400},
                  {"effect_size": 0.50, "power": 0.92, "n_sims": 400}],
        "target": 0.8, "transform": "identity", "level": 0.95}),
    ("log_axis", {
        "curve": [{"items": 4, "power": 0.18, "n_sims": 300},
                  {"items": 8, "power": 0.36, "n_sims": 300},
                  {"items": 16, "power": 0.59, "n_sims": 300},
                  {"items": 32, "power": 0.80, "n_sims": 300},
                  {"items": 64, "power": 0.93, "n_sims": 300}],
        "target": 0.75, "transform": "log", "level": 0.9}),
]

# Inputs every refusal path should reject, so that the two engines are checked to refuse the
# same things with the same words rather than only to agree when they succeed.
REFUSALS = [
    ("target_out_of_range", {
        "curve": CASES[0][1]["curve"], "target": 1.2, "transform": "sqrt", "level": 0.95}),
    ("never_reaches_target", {
        "curve": CASES[0][1]["curve"], "target": 0.95, "transform": "sqrt", "level": 0.95}),
    ("too_few_points", {
        "curve": CASES[0][1]["curve"][:2], "target": 0.3, "transform": "sqrt", "level": 0.95}),
    ("flat_curve", {
        "curve": [{"n_subject": 10, "power": 0.5, "n_sims": 100},
                  {"n_subject": 20, "power": 0.5, "n_sims": 100},
                  {"n_subject": 30, "power": 0.5, "n_sims": 100}],
        "target": 0.5, "transform": "sqrt", "level": 0.95}),
    ("slope_not_distinguishable", {
        "curve": [{"n_subject": 10, "power": 0.50, "n_sims": 30},
                  {"n_subject": 20, "power": 0.52, "n_sims": 30},
                  {"n_subject": 30, "power": 0.49, "n_sims": 30},
                  {"n_subject": 40, "power": 0.51, "n_sims": 30},
                  {"n_subject": 50, "power": 0.50, "n_sims": 30}],
        "target": 0.5, "transform": "sqrt", "level": 0.95}),
    ("solve_outside_range", {
        "curve": [{"n_subject": 10, "power": 0.50, "n_sims": 400},
                  {"n_subject": 20, "power": 0.50, "n_sims": 400},
                  {"n_subject": 30, "power": 0.50, "n_sims": 400},
                  {"n_subject": 40, "power": 0.50, "n_sims": 400},
                  {"n_subject": 50, "power": 0.90, "n_sims": 400}],
        "target": 0.85, "transform": "sqrt", "level": 0.95}),
    ("negative_swept_value", {
        "curve": [{"delta": -0.2, "power": 0.20, "n_sims": 200},
                  {"delta": 0.0, "power": 0.50, "n_sims": 200},
                  {"delta": 0.2, "power": 0.85, "n_sims": 200}],
        "target": 0.6, "transform": "sqrt", "level": 0.95}),
    ("repeated_swept_value", {
        "curve": [{"n_subject": 10, "effect": "a", "power": 0.2, "n_sims": 100},
                  {"n_subject": 10, "effect": "b", "power": 0.4, "n_sims": 100},
                  {"n_subject": 20, "effect": "a", "power": 0.5, "n_sims": 100},
                  {"n_subject": 20, "effect": "b", "power": 0.7, "n_sims": 100}],
        "target": 0.5, "transform": "sqrt", "level": 0.95}),
]

# The fields compared. The strings have to match exactly; the numbers are compared relatively.
NUMERIC = ("value", "lo", "hi", "se", "dispersion", "intercept", "slope", "x_min", "x_max",
           "n", "n_lo", "n_hi")
EXACT = ("x", "y", "transform", "n_points", "target", "level")

R_DRIVER = r"""
args <- commandArgs(trailingOnly = TRUE)
src <- args[1]
for (f in sort(list.files(src, pattern = "[.]R$", full.names = TRUE))) source(f)
cases <- jsonlite::fromJSON(args[2], simplifyDataFrame = FALSE)
out <- lapply(cases, function(cs) {
  curve <- do.call(rbind, lapply(cs$curve, function(row)
    as.data.frame(row, stringsAsFactors = FALSE)))
  res <- tryCatch(
    target_n(curve, target = cs$target, transform = cs$transform, level = cs$level),
    error = function(e) list(error = conditionMessage(e)))
  lapply(res, function(v) if (is.numeric(v)) as.numeric(v) else v)
})
cat(jsonlite::toJSON(out, auto_unbox = TRUE, digits = NA, null = "null"))
"""


def _run_r(cases):
    """Run target_n() in R over the same cases and return its results."""
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        cases_path = os.path.join(tmp, "cases.json")
        driver_path = os.path.join(tmp, "driver.R")
        with open(cases_path, "w") as f:
            json.dump([c[1] for c in cases], f)
        with open(driver_path, "w") as f:
            f.write(R_DRIVER)
        src = os.path.join(ROOT, "r", "pilotr", "R")
        proc = subprocess.run([RSCRIPT, "--vanilla", driver_path, src, cases_path],
                              capture_output=True, text=True)
    if proc.returncode != 0:
        raise SystemExit("the R driver failed:\n" + proc.stderr)
    return json.loads(proc.stdout)


def _python(case):
    try:
        return target_n(case["curve"], target=case["target"],
                        transform=case["transform"], level=case["level"])
    except ValueError as e:
        return {"error": str(e)}


def main() -> int:
    cases = CASES + REFUSALS
    r_out = _run_r(cases)
    failures = 0
    worst = 0.0
    print("%-28s %-10s %s" % ("case", "parity", "detail"))
    print("-" * 92)
    for (name, case), rr in zip(cases, r_out):
        pr = _python(case)
        if ("error" in rr) != ("error" in pr):
            print("%-28s %-10s %s" % (name, "DIFFERS", "one engine refused and the other did not"))
            failures += 1
            continue
        if "error" in rr:
            same = rr["error"] == pr["error"]
            print("%-28s %-10s %s" % (name, "ok" if same else "DIFFERS",
                                      "refusal text identical" if same else
                                      "R=%s\n%39sPy=%s" % (rr["error"], "", pr["error"])))
            failures += 0 if same else 1
            continue
        worst_here = 0.0
        detail = []
        for k in EXACT:
            if rr[k] != pr[k]:
                detail.append("%s: R=%r Py=%r" % (k, rr[k], pr[k]))
        for k in NUMERIC:
            a, b = float(rr[k]), float(pr[k])
            rel = 0.0 if a == b else abs(a - b) / max(abs(a), abs(b), 1e-300)
            worst_here = max(worst_here, rel)
            if not math.isfinite(rel) or rel > TOLERANCE:
                detail.append("%s: R=%.17g Py=%.17g (rel %.3g)" % (k, a, b, rel))
        worst = max(worst, worst_here)
        if detail:
            print("%-28s %-10s %s" % (name, "DIFFERS", "; ".join(detail)))
            failures += 1
        else:
            print("%-28s %-10s worst relative difference %.3g" % (name, "ok", worst_here))

    print("\nworst relative difference over all solved cases: %.3g (tolerance %g)"
          % (worst, TOLERANCE))
    print("%d of %d cases failed" % (failures, len(cases)))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
