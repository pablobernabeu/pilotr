"""The Type I error rate of the Python twin's Wald z test against the R twin's Satterthwaite test.

The Python power_mixed() tests each effect with a Wald z, which treats the estimate over its
standard error as normal and so ignores the uncertainty in the variance components. With few
subjects or items, such a test rejects a true null hypothesis more often than a test with
Satterthwaite's degrees of freedom, the test of the R twin's power_mixed() (Luke, 2017,
doi:10.3758/s13428-016-0809-y). This script measures both rejection rates on the same replicates of
a small null design, which the Python documentation quotes.

The design is the 12 x 8 priming design of the Python power guide cut to 8 subjects and 6 items,
with random intercepts alone and no true effect. Both twins simulate the same data, which the
parity harness checks, and both fit the model the specification implies,
.y ~ cond + (1 | subject) + (1 | item), by REML. Each runs its own power_mixed() from this source
tree, the R twin through Rscript, as tools/parity/solve_cross.py runs it.

Results are written to wald_null_rate.txt beside this file, with LF endings. Every replicate takes
its seed from the specification, so a rerun with the same statsmodels and lme4 reproduces the file.

Usage: python tools/calibration/wald_null_rate.py [replicates] [output-file]

The default of 300 replicates takes about a minute in each twin on one core.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "python"))

import pilotr  # noqa: E402
from pilotr import power_mixed  # noqa: E402

# The interpreter on PATH by default, as in tools/parity. PILOTR_RSCRIPT overrides it.
RSCRIPT = os.environ.get("PILOTR_RSCRIPT", "Rscript")

SPEC = {
    "name": "null_8x6", "seed": 1,
    "units": {"subject": {"n": 8}, "item": {"n": 6}},
    "factors": [{"name": "condition", "levels": ["related", "unrelated"],
                 "contrasts": {"cond": [-0.5, 0.5]}, "vary_within": ["subject", "item"]}],
    "fixed": {"intercept": 6, "coefficients": {"cond": 0}},
    "random": {"subject": {"intercept_sd": 0.12}, "item": {"intercept_sd": 0.08}},
    "response": {"family": "shifted_lognormal", "name": "RT", "sigma": 0.3, "shift": 200},
}

R_DRIVER = r"""
args <- commandArgs(trailingOnly = TRUE)
for (f in sort(list.files(args[1], pattern = "[.]R$", full.names = TRUE))) source(f)
r <- power_mixed(load_spec(args[2]), n_sims = as.integer(args[3]))
cat(jsonlite::toJSON(list(
  power = unname(r$power[["cond"]]), power_mcse = unname(r$power_mcse[["cond"]]),
  power_lo = unname(r$power_lo[["cond"]]), power_hi = unname(r$power_hi[["cond"]]),
  n_significant = unname(r$n_significant[["cond"]]), n_returned = r$n_returned,
  fitter = r$fitter, versions = sprintf("%s with lme4 %s and lmerTest %s", R.version.string,
    as.character(utils::packageVersion("lme4")), as.character(utils::packageVersion("lmerTest")))),
  auto_unbox = TRUE, digits = NA))
"""


def _run_r(n_sims):
    """Run the R twin's power_mixed() over the same specification and return its rates."""
    with tempfile.TemporaryDirectory() as tmp:
        spec_path = os.path.join(tmp, "spec.json")
        driver_path = os.path.join(tmp, "driver.R")
        with open(spec_path, "w", encoding="utf-8") as f:
            json.dump(SPEC, f)
        with open(driver_path, "w", encoding="utf-8") as f:
            f.write(R_DRIVER)
        proc = subprocess.run([RSCRIPT, "--vanilla", driver_path,
                               os.path.join(ROOT, "r", "pilotr", "R"), spec_path, str(n_sims)],
                              capture_output=True, text=True, encoding="utf-8", errors="replace")
    if proc.returncode != 0:
        raise SystemExit("the R driver failed:\n" + proc.stderr)
    return json.loads(proc.stdout)


def main() -> int:
    n_sims = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    out_path = sys.argv[2] if len(sys.argv) > 2 else os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "wald_null_rate.txt")
    import statsmodels

    run = power_mixed(SPEC, n_sims=n_sims)
    py = dict(run["effects"]["cond"], n_returned=run["n_returned"], fitter=run["backend"])
    r = _run_r(n_sims)

    def row(label, x):
        return "  %-8s %6.3f  %6.3f  [%.3f, %.3f]  %4d of %d" % (
            label, x["power"], x["power_mcse"], x["power_lo"], x["power_hi"],
            x["n_significant"], x["n_returned"])

    lines = [
        "The Type I error rate of a Wald z test and of Satterthwaite's test on a small null design",
        "",
        "Produced by tools/calibration/wald_null_rate.py with %d replicates, under" % n_sims,
        "Python %s with pilotr %s and statsmodels %s, and" % (
            sys.version.split()[0], pilotr.__version__, statsmodels.__version__),
        "%s." % r["versions"],
        "Rerun it to regenerate this file. Nothing here is timed or dated, so a rerun that "
        "changes a",
        "line has changed a result.",
        "",
        "The Python power guide's priming design cut to 8 subjects and 6 items, with random",
        "intercepts alone (standard deviations 0.12 by subject and 0.08 by item), a residual",
        "standard deviation of 0.3 on the log scale and no true effect. Both twins fit",
        ".y ~ cond + (1 | subject) + (1 | item) by REML to the same replicates and test cond",
        "two-sided at alpha 0.05.",
        "",
        "  Python   %s" % py["fitter"],
        "  R        %s" % r["fitter"],
        "",
        "  twin     rate    MCSE    Wilson 95%       significant",
        row("Python", py),
        row("R", r),
        "",
    ]
    with open(out_path, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(lines))
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
