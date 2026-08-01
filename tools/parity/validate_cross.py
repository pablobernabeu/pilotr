"""Cross-language agreement check for validate_spec.

The two validators are the gate that decides whether a specification is usable, so a
specification accepted by one implementation and refused by the other is itself a parity bug: it
would mean a design that runs in R and fails in Python, or worse, one that runs in both but with
different meaning. This script feeds an identical battery of valid and invalid specifications to
both and reports any verdict that differs.

Usage: python tools/parity/validate_cross.py
Exit status is 0 when the two agree on every case.
"""

from __future__ import annotations

import copy
import json
import os
import subprocess
import sys
import tempfile

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "python"))

from pilotr.validate import validate_spec  # noqa: E402

RSCRIPT = os.environ.get("PILOTR_RSCRIPT", r"C:\Program Files\R\R-4.6.1\bin\x64\Rscript.exe")


def base_spec():
    with open(os.path.join(ROOT, "spec", "examples", "crossed_mixed_rt.json")) as f:
        return json.load(f)


def _mut(fn):
    s = base_spec()
    fn(s)
    return s


def _set_family(s, **kw):
    s["response"] = kw


# (label, spec) pairs. Each is fed to both validators; the verdicts must match.
def cases():
    out = []
    for name in sorted(os.listdir(os.path.join(ROOT, "spec", "examples"))):
        if name.endswith(".json"):
            with open(os.path.join(ROOT, "spec", "examples", name)) as f:
                out.append(("shipped:" + name, json.load(f)))
    for name in sorted(os.listdir(os.path.join(ROOT, "tools", "parity", "cases"))):
        if name.endswith(".json"):
            with open(os.path.join(ROOT, "tools", "parity", "cases", name)) as f:
                out.append(("case:" + name, json.load(f)))

    out += [
        ("mistyped coefficient key",
         _mut(lambda s: s["fixed"].__setitem__("coefficients", {"cnod": 0.05}))),
        ("mistyped slope key",
         _mut(lambda s: s["random"]["subject"].__setitem__("slopes", {"conditon": 0.04}))),
        ("interaction coefficient with a bad component",
         _mut(lambda s: s["fixed"].__setitem__("coefficients", {"cond:nope": 0.05}))),
        ("foreign response parameter", _mut(lambda s: s["response"].__setitem__("phi", 8))),
        ("round on an integer family",
         _mut(lambda s: _set_family(s, family="poisson", name="count", round=2))),
        ("unknown top-level field", _mut(lambda s: s.__setitem__("notafield", 1))),
        ("unknown factor field",
         _mut(lambda s: s["factors"][0].__setitem__("levls", ["a", "b"]))),
        ("varies_by trial",
         _mut(lambda s: s.__setitem__("predictors", [{"name": "ISI", "varies_by": "trial"}]))),
        ("varies_by observation without declaring 0.3",
         _mut(lambda s: s.__setitem__("predictors",
                                      [{"name": "ISI", "varies_by": "observation"}]))),
        ("varies_by observation declaring 0.3",
         _mut(lambda s: (s.__setitem__("spec_version", "0.3"),
                         s.__setitem__("predictors", [{"name": "ISI", "varies_by": "observation"}]),
                         s["fixed"]["coefficients"].__setitem__("ISI", 0.1)))),
        ("contrast length mismatch",
         _mut(lambda s: s["factors"][0].__setitem__("contrasts", {"cond": [-0.5, 0, 0.5]}))),
        ("correlation naming a missing term",
         _mut(lambda s: s["random"]["subject"].__setitem__("correlations",
                                                           {"intercept,nope": 0.2}))),
        ("correlation out of range",
         _mut(lambda s: s["random"]["subject"].__setitem__("correlations",
                                                           {"intercept,cond": 1.4}))),
        ("correlation key with one term",
         _mut(lambda s: s["random"]["subject"].__setitem__("correlations", {"intercept": 0.2}))),
        ("per_subject over the item count",
         _mut(lambda s: s["units"]["item"].__setitem__("per_subject", 999))),
        ("per_subject on subject",
         _mut(lambda s: s["units"]["subject"].__setitem__("per_subject", 2))),
        ("extra group without over/n",
         _mut(lambda s: s["random"].__setitem__("site", {"intercept_sd": 0.5}))),
        ("extra group with a valid over/n",
         _mut(lambda s: s["random"].__setitem__("site", {"intercept_sd": 0.5, "over": "subject",
                                                        "n": 5}))),
        ("over/n on subject",
         _mut(lambda s: s["random"]["subject"].__setitem__("over", "subject"))),
        ("thresholds not increasing",
         _mut(lambda s: _set_family(s, family="ordinal", name="r", thresholds=[1, 0.5, 2]))),
        ("valid ordinal",
         _mut(lambda s: _set_family(s, family="ordinal", name="r", thresholds=[-1, 0, 1]))),
        ("negative sigma", _mut(lambda s: s["response"].__setitem__("sigma", -1))),
        ("missing sigma", _mut(lambda s: s["response"].pop("sigma"))),
        ("missing response", _mut(lambda s: s.pop("response"))),
        ("missing units.subject", _mut(lambda s: s["units"].pop("subject"))),
        ("zero subjects", _mut(lambda s: s["units"]["subject"].__setitem__("n", 0))),
        ("non-integer seed", _mut(lambda s: s.__setitem__("seed", 1.5))),
        ("factor with neither vary_within nor between",
         _mut(lambda s: s["factors"][0].pop("vary_within"))),
        ("vary_within item with no item unit",
         _mut(lambda s: (s["units"].pop("item"), s["random"].pop("item")))),
        ("negative slope sd",
         _mut(lambda s: s["random"]["subject"]["slopes"].__setitem__("cond", -0.1))),
        ("correlated false with correlations",
         _mut(lambda s: s["random"]["subject"].__setitem__("correlated", False))),
        ("correlated flag without declaring 0.3",
         _mut(lambda s: (s["random"]["subject"].pop("correlations"),
                         s["random"]["subject"].__setitem__("correlated", False)))),
        ("exgaussian without declaring 0.3",
         _mut(lambda s: _set_family(s, family="exgaussian", name="RT", sigma=0.2, beta=0.3))),
        ("exgaussian declaring 0.3",
         _mut(lambda s: (s.__setitem__("spec_version", "0.3"),
                         _set_family(s, family="exgaussian", name="RT", sigma=0.2, beta=0.3)))),
        ("spec from the future", _mut(lambda s: s.__setitem__("spec_version", "9.9"))),
        ("malformed spec_version", _mut(lambda s: s.__setitem__("spec_version", "banana"))),
        ("unknown family", _mut(lambda s: s["response"].__setitem__("family", "weibull"))),
        ("empty response name", _mut(lambda s: s["response"].__setitem__("name", ""))),
        ("duplicate predictor names",
         _mut(lambda s: (s.__setitem__("predictors",
                                       [{"name": "z", "varies_by": "subject"},
                                        {"name": "z", "varies_by": "item"}]),
                         s["fixed"]["coefficients"].__setitem__("z", 0.1)))),
        ("uniform predictor without min/max",
         _mut(lambda s: (s.__setitem__("spec_version", "0.3"),
                         s.__setitem__("predictors",
                                       [{"name": "z", "varies_by": "subject", "dist": "uniform"}]),
                         s["fixed"]["coefficients"].__setitem__("z", 0.1)))),
        ("uniform predictor with min >= max",
         _mut(lambda s: (s.__setitem__("spec_version", "0.3"),
                         s.__setitem__("predictors",
                                       [{"name": "z", "varies_by": "subject", "dist": "uniform",
                                         "min": 1, "max": 0}]),
                         s["fixed"]["coefficients"].__setitem__("z", 0.1)))),
        ("valid uniform predictor",
         _mut(lambda s: (s.__setitem__("spec_version", "0.3"),
                         s.__setitem__("predictors",
                                       [{"name": "z", "varies_by": "subject", "dist": "uniform",
                                         "min": 0, "max": 1}]),
                         s["fixed"]["coefficients"].__setitem__("z", 0.1)))),
        ("reliability out of range",
         _mut(lambda s: (s.__setitem__("spec_version", "0.3"),
                         s.__setitem__("predictors",
                                       [{"name": "z", "varies_by": "subject", "reliability": 1.5}]),
                         s["fixed"]["coefficients"].__setitem__("z", 0.1)))),
        ("valid reliability",
         _mut(lambda s: (s.__setitem__("spec_version", "0.3"),
                         s.__setitem__("predictors",
                                       [{"name": "z", "varies_by": "subject", "reliability": 0.8}]),
                         s["fixed"]["coefficients"].__setitem__("z", 0.1)))),
    ]
    return out


R_DRIVER = r'''
args <- commandArgs(trailingOnly = TRUE)
src <- args[1]; payload <- args[2]
for (f in sort(list.files(src, pattern = "\\.R$", full.names = TRUE))) source(f)
specs <- jsonlite::fromJSON(payload, simplifyVector = TRUE, simplifyDataFrame = FALSE,
                            simplifyMatrix = FALSE)
out <- character(length(specs))
for (i in seq_along(specs)) {
  out[i] <- tryCatch({
    suppressWarnings(validate_spec(specs[[i]], strict = TRUE)); "OK"
  }, error = function(e) "ERROR")
}
cat(paste(out, collapse = "\n"))
'''


def main() -> int:
    battery = cases()
    py_verdicts = []
    for _label, spec in battery:
        try:
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                validate_spec(copy.deepcopy(spec), strict=True)
            py_verdicts.append("OK")
        except ValueError:
            py_verdicts.append("ERROR")

    with tempfile.TemporaryDirectory() as td:
        payload = os.path.join(td, "specs.json")
        with open(payload, "w") as f:
            json.dump([s for _l, s in battery], f)
        driver = os.path.join(td, "driver.R")
        with open(driver, "w") as f:
            f.write(R_DRIVER)
        proc = subprocess.run(
            [RSCRIPT, driver, os.path.join(ROOT, "r", "pilotr", "R"), payload],
            capture_output=True, text=True)
        if proc.returncode != 0:
            print("R driver failed:\n" + proc.stdout + proc.stderr)
            return 1
        r_verdicts = proc.stdout.strip().splitlines()

    if len(r_verdicts) != len(battery):
        print("R returned %d verdicts for %d cases" % (len(r_verdicts), len(battery)))
        return 1

    print("%-56s %-8s %-8s %s" % ("case", "R", "Python", ""))
    print("-" * 88)
    mismatches = 0
    for (label, _s), rv, pv in zip(battery, r_verdicts, py_verdicts):
        flag = ""
        if rv != pv:
            flag = "<== DISAGREE"
            mismatches += 1
        print("%-56s %-8s %-8s %s" % (label[:56], rv, pv, flag))
    n_ok = sum(1 for v in py_verdicts if v == "OK")
    print("\n%d cases: %d accepted, %d refused, %d disagreements"
          % (len(battery), n_ok, len(battery) - n_ok, mismatches))
    return 1 if mismatches else 0


if __name__ == "__main__":
    sys.exit(main())
