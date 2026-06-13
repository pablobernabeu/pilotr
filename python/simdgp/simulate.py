"""Spec parsing + the generative simulation engine.

Implements the data-generating process the prototype app lacks: a linear predictor with
user-specified fixed effect sizes (for categorical contrasts AND continuous predictors,
including their interactions) plus crossed by-subject / by-item random intercepts and slopes
(on contrasts or continuous predictors), pushed through a link + response family. RNG draw
order follows spec/SPEC.md: continuous predictors -> subject random effects -> item random
effects -> per-row residuals. Specs without a `predictors` block keep the original stream.
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
    return repr(v) if isinstance(v, float) else str(v)


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


def _design_value(cvals, key):
    """A coefficient/slope key is a column name or an 'a:b' interaction (product of columns)."""
    if ":" in key:
        v = 1.0
        for part in key.split(":"):
            v *= cvals.get(part, 0.0)
        return v
    return cvals.get(key, 0.0)


def simulate(spec) -> Dataset:
    if isinstance(spec, str):
        spec = load_spec(spec)

    S = spec["units"]["subject"]["n"]
    has_item = "item" in spec["units"]
    I = spec["units"]["item"]["n"] if has_item else 1

    factors = spec.get("factors", [])
    predictors = spec.get("predictors", [])
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
                cvals, labels = {}, {}
                for f in factors:
                    li = level_idx[f["name"]]
                    labels[f["name"]] = f["levels"][li]
                    for col, vals in f["contrasts"].items():
                        cvals[col] = vals[li]
                rows.append({"subject": s, "item": t if has_item else None,
                             "labels": labels, "cvals": cvals})

    rng = RNG(spec["seed"])

    # ---- continuous predictors: one draw per unit (subject- or item-level) ----
    pred_values = {}
    for p in predictors:
        unit, n_unit = p["varies_by"], (S if p["varies_by"] == "subject" else I)
        if unit == "item" and not has_item:
            raise ValueError(f"predictor '{p['name']}' varies_by item but design has no items")
        mean, sd = p.get("mean", 0.0), p.get("sd", 1.0)
        pred_values[p["name"]] = {u: mean + sd * rng.normal() for u in range(1, n_unit + 1)}
    for r in rows:
        for p in predictors:
            r["cvals"][p["name"]] = pred_values[p["name"]][r["subject"] if p["varies_by"] == "subject" else r["item"]]

    # ---- random effects (subject then item) ----
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

    # ---- additional grouping factors (e.g. units nested in higher-level clusters) ----
    # Any random-effect entry other than subject/item declares `over` (which unit it groups)
    # and `n` (number of groups); units are assigned to groups in equal blocks.
    extra = [(k, v) for k, v in random_spec.items() if k not in ("subject", "item")]
    b_group, group_meta = {}, {}
    for gname, gspec in extra:
        over, K = gspec["over"], gspec["n"]
        n_over = S if over == "subject" else I
        cols, L = _ranef(gspec)
        group_of = {u: ((u - 1) * K) // n_over for u in range(1, n_over + 1)}
        group_meta[gname] = (over, cols, group_of)
        b_group[gname] = {g: matvec(L, rng.normals(len(cols))) for g in range(K)}

    # ---- linear predictor + response per row (residual draws happen here) ----
    intercept = spec["fixed"]["intercept"]
    coeffs = spec["fixed"]["coefficients"]
    resp = spec["response"]
    family, yname = resp["family"], resp["name"]
    sigma, shift = resp.get("sigma"), resp.get("shift", 0.0)
    thresholds, ndp = resp.get("thresholds"), resp.get("round")

    out_cols = (["subject"] + (["item"] if has_item else []) + [g for g, _ in extra] +
                [f["name"] for f in factors] + [p["name"] for p in predictors] + [yname])
    out_rows = []
    for r in rows:
        cv = r["cvals"]
        eta = intercept + sum(beta * _design_value(cv, col) for col, beta in coeffs.items())
        if r["subject"] in b_subject:
            b = b_subject[r["subject"]]
            eta += b[0] + sum(b[j] * cv.get(subj_cols[j], 0.0) for j in range(1, len(subj_cols)))
        if has_item and r["item"] in b_item:
            b = b_item[r["item"]]
            eta += b[0] + sum(b[j] * cv.get(item_cols[j], 0.0) for j in range(1, len(item_cols)))
        for gname, (over, cols, group_of) in group_meta.items():
            b = b_group[gname][group_of[r[over]]]
            eta += b[0] + sum(b[j] * cv.get(cols[j], 0.0) for j in range(1, len(cols)))

        if family == "gaussian":
            y = eta + sigma * rng.normal()
        elif family == "shifted_lognormal":
            y = shift + math.exp(eta + sigma * rng.normal())
        elif family == "lognormal":
            y = math.exp(eta + sigma * rng.normal())
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
        for gname, (over, cols, group_of) in group_meta.items():
            out[gname] = group_of[r[over]] + 1
        out.update(r["labels"])
        for p in predictors:
            out[p["name"]] = cv[p["name"]]
        out[yname] = y
        out_rows.append(out)

    return Dataset(out_cols, out_rows)
