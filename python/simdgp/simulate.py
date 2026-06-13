"""Spec parsing + the generative simulation engine.

Implements the data-generating process the prototype app lacks: a linear predictor with
user-specified fixed effect sizes plus crossed by-subject / by-item random intercepts and
slopes, pushed through a link + response family. RNG draw order follows spec/SPEC.md.
"""

from __future__ import annotations
import json, math, itertools
from .core import RNG, cholesky, matvec, inv_logit, poisson_inv, ordinal_inv


class Dataset:
    def __init__(self, columns, rows):
        self.columns = columns
        self.rows = rows  # list of dict

    def __len__(self):
        return len(self.rows)

    def column(self, name):
        return [r[name] for r in self.rows]

    def head(self, n=6):
        return self.rows[:n]

    def to_csv(self, path):
        with open(path, "w", newline="") as f:
            f.write(",".join(self.columns) + "\n")
            for r in self.rows:
                f.write(",".join(_fmt(r[c]) for c in self.columns) + "\n")


def _fmt(v):
    if isinstance(v, float):
        return repr(v)
    return str(v)


def load_spec(path):
    with open(path, "r") as f:
        return json.load(f)


def _ranef(unit_spec):
    """Return (columns, lower-Cholesky L) for a unit's random-effect covariance."""
    cols = ["intercept"] + list(unit_spec.get("slopes", {}).keys())
    sds = [unit_spec["intercept_sd"]] + [unit_spec.get("slopes", {})[c] for c in cols[1:]]
    n = len(cols)
    R = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    for key, val in unit_spec.get("correlations", {}).items():
        a, b = [s.strip() for s in key.replace("~", ",").split(",")]
        i, j = cols.index(a), cols.index(b)
        R[i][j] = R[j][i] = val
    cov = [[sds[i] * R[i][j] * sds[j] for j in range(n)] for i in range(n)]
    return cols, cholesky(cov)


def simulate(spec) -> Dataset:
    if isinstance(spec, str):
        spec = load_spec(spec)

    S = spec["units"]["subject"]["n"]
    has_item = "item" in spec["units"]
    I = spec["units"]["item"]["n"] if has_item else 1

    factors = spec["factors"]
    within = [f for f in factors if f.get("vary_within")]
    between = [f for f in factors if f.get("between")]

    # ---- build design rows in canonical order (consumes no randomness) ----
    within_level_ranges = [range(len(f["levels"])) for f in within]
    rows = []
    for s in range(1, S + 1):
        for t in range(1, I + 1):
            for combo in (itertools.product(*within_level_ranges) if within else [()]):
                level_idx = {}
                for f, li in zip(within, combo):
                    level_idx[f["name"]] = li
                for f in between:
                    n_lev = len(f["levels"])
                    unit = s if f["between"] == "subject" else t
                    n_unit = S if f["between"] == "subject" else I
                    level_idx[f["name"]] = ((unit - 1) * n_lev) // n_unit
                cvals = {}
                labels = {}
                for f in factors:
                    li = level_idx[f["name"]]
                    labels[f["name"]] = f["levels"][li]
                    for col, vals in f["contrasts"].items():
                        cvals[col] = vals[li]
                rows.append({"subject": s, "item": t if has_item else None,
                             "labels": labels, "cvals": cvals})

    # ---- draw random effects in the documented order ----
    rng = RNG(spec["seed"])
    random_spec = spec.get("random", {}) or {}
    b_subject, subj_cols = {}, []
    if "subject" in random_spec:
        subj_cols, L = _ranef(random_spec["subject"])
        for s in range(1, S + 1):
            b_subject[s] = matvec(L, rng.normals(len(subj_cols)))
    b_item, item_cols = {}, []
    if has_item and "item" in random_spec:
        item_cols, L = _ranef(random_spec["item"])
        for t in range(1, I + 1):
            b_item[t] = matvec(L, rng.normals(len(item_cols)))

    # ---- linear predictor + response per row (residual draws happen here) ----
    intercept = spec["fixed"]["intercept"]
    coeffs = spec["fixed"]["coefficients"]
    resp = spec["response"]
    family, yname = resp["family"], resp["name"]
    sigma = resp.get("sigma")
    shift = resp.get("shift", 0.0)
    thresholds = resp.get("thresholds")
    ndp = resp.get("round")

    out_cols = (["subject"] + (["item"] if has_item else []) +
                [f["name"] for f in factors] + [yname])
    out_rows = []
    for r in rows:
        eta = intercept + sum(beta * r["cvals"].get(col, 0.0) for col, beta in coeffs.items())
        if r["subject"] in b_subject:
            b = b_subject[r["subject"]]
            eta += b[0] + sum(b[j] * r["cvals"].get(subj_cols[j], 0.0) for j in range(1, len(subj_cols)))
        if has_item and r["item"] in b_item:
            b = b_item[r["item"]]
            eta += b[0] + sum(b[j] * r["cvals"].get(item_cols[j], 0.0) for j in range(1, len(item_cols)))

        if family == "gaussian":
            y = eta + sigma * rng.normal()
        elif family == "shifted_lognormal":
            y = shift + math.exp(eta + sigma * rng.normal())
        elif family == "bernoulli":
            y = 1 if rng.uniform() < inv_logit(eta) else 0
        elif family == "poisson":
            y = poisson_inv(math.exp(eta), rng.uniform())
        elif family == "ordinal":
            y = ordinal_inv(eta, thresholds, rng.uniform())
        else:
            raise ValueError(f"unknown family: {family}")
        if ndp is not None and isinstance(y, float):
            y = round(y, ndp)

        out = {"subject": r["subject"]}
        if has_item:
            out["item"] = r["item"]
        out.update(r["labels"])
        out[yname] = y
        out_rows.append(out)

    return Dataset(out_cols, out_rows)
