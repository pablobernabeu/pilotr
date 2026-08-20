"""Spec parsing and the generative simulation engine.

The model is a linear predictor with user-specified fixed effect sizes, covering categorical
contrasts and continuous predictors as well as their interactions. To this it adds crossed
by-subject and by-item random intercepts and slopes (on contrasts or continuous predictors), and
passes the result through a link and a response family. The RNG draw order follows spec/SPEC.md,
namely per-subject item subsets (partial crossing only) -> continuous predictors -> subject random
effects -> item random effects -> extra grouping-factor random effects -> per-row response
draws (one deviate per row, except the beta family's rejection sampler, which consumes a
variable number). Specs without partial crossing or a `predictors` block keep the original
stream.
"""

from __future__ import annotations
import json, math, itertools
from .core import RNG, cholesky, matvec, inv_logit, poisson_inv, ordinal_inv, beta_draw
from .validate import validate_spec


class Dataset:
    """A simulated data set: named columns and a list of row dicts.

    Returned by `simulate`. Lightweight and dependency-free; build a pandas `DataFrame`
    from `dataset.rows` when you need one.

    Attributes
    ----------
    columns : list of str
        Column names, in order.
    rows : list of dict
        One dict per observation, keyed by column name.
    """

    def __init__(self, columns, rows):
        self.columns = columns
        self.rows = rows  # list of dict

    def __len__(self):
        return len(self.rows)

    def column(self, name):
        """Return the values in column `name` as a list."""
        return [r[name] for r in self.rows]

    def head(self, n=6):
        """Return the first `n` rows (default 6) as a list of dicts."""
        return self.rows[:n]

    def to_csv(self, path):
        """Write the data set to `path` as CSV: a header row then one row per observation."""
        with open(path, "w", newline="") as f:
            f.write(",".join(self.columns) + "\n")
            for r in self.rows:
                f.write(",".join(_fmt(r[c]) for c in self.columns) + "\n")


def _fmt(v):
    return repr(v) if isinstance(v, float) else str(v)


def load_spec(path, validate=True):
    """Load a JSON design specification from a file.

    The specification is validated by default, because several ways of getting one wrong produce
    plausible data and no error at all: a mistyped coefficient key resolves to no column and so
    silently sets that effect to zero, and a response parameter left over from another family is
    ignored. Validation also refuses a specification declaring a ``spec_version`` newer than this
    implementation understands, and never reads such a file in part.

    Parameters
    ----------
    path : str
        Path to a JSON design-specification file.
    validate : bool
        Whether to validate the specification after reading it. ``True`` (the default) applies
        `validate_spec` strictly; ``False`` skips validation; ``"lenient"`` or any other truthy
        non-``True`` value validates with ``strict=False``.

    Returns
    -------
    dict
        The parsed specification, ready for `simulate` or `power`.
    """
    with open(path) as f:
        spec = json.load(f)
    if validate is not False:
        validate_spec(spec, strict=validate is True)
    return spec


def _ranef(unit_spec, label=None):
    """Return (columns, lower-Cholesky L) for a unit's random-effect covariance.

    `label` names the grouping factor, so that a covariance which is not positive definite can
    be reported against the entry the user wrote, where it used to appear as an anonymous
    matrix failure.
    """
    cols = ["intercept"] + list(unit_spec.get("slopes", {}).keys())
    sds = [unit_spec["intercept_sd"]] + [unit_spec.get("slopes", {})[c] for c in cols[1:]]
    n = len(cols)
    R = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    for key, val in unit_spec.get("correlations", {}).items():
        a, b = [s.strip() for s in key.replace("~", ",").split(",")]
        if a not in cols or b not in cols:
            raise ValueError(
                "correlation '%s' for '%s' names a random-effect term that does not exist; "
                "available terms are %s"
                % (key, label if label is not None else "unknown group",
                   ", ".join("'%s'" % c for c in cols)))
        i, j = cols.index(a), cols.index(b)
        R[i][j] = R[j][i] = val
    # Σ = D R D, bracketed as (sd_i · sd_j) · r_ij to match R's `outer(sds, sds) * R`.
    # Multiplication is commutative but not associative in floating point, so the alternative
    # grouping (sd_i · r_ij) · sd_j lands on a different double for some inputs, and the
    # difference propagates through the Cholesky factor into every random effect drawn.
    cov = [[(sds[i] * sds[j]) * R[i][j] for j in range(n)] for i in range(n)]
    return cols, cholesky(cov, label=label, cols=cols)


def _attenuate(latent, pmean, psd, rho, z):
    """Contaminate a latent predictor value down to a stated reliability.

        observed = mean + (true - mean + sd * sqrt((1 - rho) / rho) * z) * sqrt(rho)

    The observed variable then has the same variance as the latent one and correlates sqrt(rho)
    with it. Reliability in the classical sense is that squared correlation, which is why the
    field is rho, and the correlation itself is its square root.

    The attenuation is sqrt(rho), where the textbook regression-dilution result gives rho,
    because both variables are put on the same variance here. Standardising the observed
    variable back to the latent one's variance absorbs the 1 / sqrt(rho) factor that result
    carries.

    The moments used are the population ones. R's mean() and sd() accumulate in long double and
    Python's do not, so standardising against the sample mean and standard deviation of the
    values drawn would reintroduce exactly the cross-language divergence this release removes.
    """
    return pmean + (latent - pmean + psd * math.sqrt((1.0 - rho) / rho) * z) * math.sqrt(rho)


def _design_value(cvals, key):
    """A coefficient/slope key is a column name or an 'a:b' interaction (product of columns)."""
    if ":" in key:
        v = 1.0
        for part in key.split(":"):
            v *= cvals.get(part, 0.0)
        return v
    return cvals.get(key, 0.0)


def _sample_items(rng, n_items, m):
    """Sample m distinct items from 1..n_items via partial Fisher-Yates on the shared RNG.
    This supports partial crossing, where each subject sees a self-selected subset of items."""
    pool = list(range(1, n_items + 1))
    for k in range(m):
        j = k + int(rng.uniform() * (n_items - k))
        pool[k], pool[j] = pool[j], pool[k]
    return sorted(pool[:m])


def simulate(spec, validate=True) -> Dataset:
    """Simulate a data set from a design specification.

    Build a linear predictor from the fixed effect sizes (categorical contrasts, continuous
    predictors and their interactions) plus the crossed by-subject and by-item random
    intercepts and slopes, then map it through the chosen response family.

    Parameters
    ----------
    spec : dict or str
        A design specification, either an already-parsed `dict` or a path to a JSON spec
        file. The format is documented at
        https://github.com/pablobernabeu/pilotr/blob/main/spec/SPEC.md.
    validate : bool
        Whether to validate the specification first. The default ``True`` catches the errors
        that would otherwise pass silently, such as a mistyped coefficient key, which resolves
        to no column and so sets that effect to zero. Validation costs a few milliseconds, so
        replicate loops validate once and then pass ``False``.

    Returns
    -------
    Dataset
        A table with one row per observation: a `subject` column, an optional `item` column,
        any grouping, factor, and continuous-predictor columns, and the response column named
        by ``spec["response"]["name"]``.
    """
    if isinstance(spec, str):
        spec = load_spec(spec, validate=validate)
    elif validate is not False:
        validate_spec(spec, strict=validate is True)

    S = spec["units"]["subject"]["n"]
    has_item = "item" in spec["units"]
    I = spec["units"]["item"]["n"] if has_item else 1  # noqa: E741  S and I as in simulate.R

    factors = spec.get("factors", [])
    predictors = spec.get("predictors", [])
    within = [f for f in factors if f.get("vary_within")]
    between = [f for f in factors if f.get("between")]

    rng = RNG(spec["seed"])
    per_subject = spec["units"]["item"].get("per_subject") if has_item else None
    if per_subject is not None:
        if per_subject < 1:
            raise ValueError(f"per_subject ({per_subject}) must be at least 1")
        if per_subject > I:
            raise ValueError(f"per_subject ({per_subject}) cannot exceed the number of items ({I})")

    # ---- build design rows in canonical order ----
    # When per_subject is set, each subject's item subset is sampled here, which is the first
    # RNG consumption (see spec/SPEC.md). Full-crossing specs draw nothing here, keeping the stream.
    within_level_ranges = [range(len(f["levels"])) for f in within]
    rows = []
    for s in range(1, S + 1):
        items_s = _sample_items(rng, I, per_subject) if per_subject else range(1, I + 1)
        for t in items_s:
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

    # ---- continuous predictors ----
    # One draw per unit for a subject- or item-level predictor, and one per row for an
    # observation-level one. The rows are already enumerated by this point, so the third case is
    # one more branch of the same dispatch, and no restructuring. Before 0.3 `varies_by` was
    # read as "subject or else item", so a predictor declared to vary by "trial" was silently
    # given one value per item; it is now validated against the three names that exist.
    #
    # `dist` chooses the draw. A normal consumes one uniform through the inverse CDF and a
    # uniform consumes one directly, so the two cost the same number of draws and switching
    # between them does not move the stream.
    #
    # `reliability` draws one further normal per value, and only when it is present and below
    # one, so a specification that does not use it keeps the original stream. The latent value
    # drives the linear predictor and any random slope keyed on the predictor, while the
    # returned data carries the contaminated observed value, which is what an analyst would
    # have measured.
    n_rows = len(rows)
    pred_latent, pred_observed, pred_unit = {}, {}, {}
    for p in predictors:
        unit = p["varies_by"]
        if unit == "item" and not has_item:
            raise ValueError(f"predictor '{p['name']}' varies_by item but design has no items")
        if unit == "subject":
            n_unit = S
        elif unit == "item":
            n_unit = I
        elif unit == "observation":
            n_unit = n_rows
        else:
            raise ValueError(f"predictor '{p['name']}' has varies_by '{unit}'; "
                             "expected 'subject', 'item' or 'observation'")
        dist = p.get("dist", "normal")
        # Population moments of the latent variable, needed by the reliability contamination.
        if dist == "uniform":
            pmean = (p["min"] + p["max"]) / 2.0
            psd = (p["max"] - p["min"]) / math.sqrt(12.0)
        else:
            pmean, psd = p.get("mean", 0.0), p.get("sd", 1.0)
        rho = p.get("reliability", 1.0)
        lat, obs = {}, {}
        for u in range(1, n_unit + 1):
            lat[u] = (p["min"] + (p["max"] - p["min"]) * rng.uniform() if dist == "uniform"
                      else pmean + psd * rng.normal())
            obs[u] = _attenuate(lat[u], pmean, psd, rho, rng.normal()) if rho < 1 else lat[u]
        pred_latent[p["name"]], pred_observed[p["name"]], pred_unit[p["name"]] = lat, obs, unit
    for i, r in enumerate(rows, start=1):
        r["obs"] = {}
        for p in predictors:
            pu = pred_unit[p["name"]]
            u = r["subject"] if pu == "subject" else (r["item"] if pu == "item" else i)
            r["cvals"][p["name"]] = pred_latent[p["name"]][u]
            r["obs"][p["name"]] = pred_observed[p["name"]][u]

    # ---- random effects (subject then item) ----
    random_spec = spec.get("random", {}) or {}
    b_subject, subj_cols = {}, []
    if "subject" in random_spec:
        subj_cols, L = _ranef(random_spec["subject"], "subject")
        for s in range(1, S + 1):
            b_subject[s] = matvec(L, rng.normals(len(subj_cols)))
    b_item, item_cols = {}, []
    if has_item and "item" in random_spec:
        item_cols, L = _ranef(random_spec["item"], "item")
        for t in range(1, I + 1):
            b_item[t] = matvec(L, rng.normals(len(item_cols)))

    # ---- additional grouping factors (e.g. units nested in higher-level clusters) ----
    # Any random-effect entry other than subject/item declares `over` (which unit it groups)
    # and `n` (number of groups). Units are assigned to groups in equal blocks.
    extra = [(k, v) for k, v in random_spec.items() if k not in ("subject", "item")]
    b_group, group_meta = {}, {}
    for gname, gspec in extra:
        over, K = gspec["over"], gspec["n"]
        n_over = S if over == "subject" else I
        cols, L = _ranef(gspec, gname)
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
    beta_exg = resp.get("beta")

    out_cols = (["subject"] + (["item"] if has_item else []) + [g for g, _ in extra] +
                [f["name"] for f in factors] + [p["name"] for p in predictors] + [yname])
    out_rows = []
    for r in rows:
        cv = r["cvals"]
        # Accumulate eta as a strict left fold, one term at a time, in the order set out in
        # spec/SPEC.md. Floating-point addition is not associative, so `intercept + sum(...)`
        # and a running `eta += ...` give answers that differ in the last few ulps; only an
        # accumulation order fixed across both languages keeps the two ports bit-identical.
        eta = intercept
        for col, beta in coeffs.items():
            eta += beta * _design_value(cv, col)
        if r["subject"] in b_subject:
            b = b_subject[r["subject"]]
            eta += b[0]
            for j in range(1, len(subj_cols)):
                eta += b[j] * _design_value(cv, subj_cols[j])
        if has_item and r["item"] in b_item:
            b = b_item[r["item"]]
            eta += b[0]
            for j in range(1, len(item_cols)):
                eta += b[j] * _design_value(cv, item_cols[j])
        for gname, (over, cols, group_of) in group_meta.items():
            b = b_group[gname][group_of[r[over]]]
            eta += b[0]
            for j in range(1, len(cols)):
                eta += b[j] * _design_value(cv, cols[j])

        if family == "gaussian":
            y = eta + sigma * rng.normal()
        elif family == "shifted_lognormal":
            y = shift + math.exp(eta + sigma * rng.normal())
        elif family == "lognormal":
            y = math.exp(eta + sigma * rng.normal())
        elif family == "exgaussian":
            # A normal plus an exponential, mean-centred by subtracting the exponential's own
            # mean so that eta stays the mean of the response. That is brms's exgaussian(mu,
            # sigma, beta) parameterisation, in which mu is the mean, so a spec and the model
            # fitted to it agree on what the intercept means. -log(u) is a unit exponential,
            # hence the two draws.
            y = eta + sigma * rng.normal() - beta_exg * (math.log(rng.uniform()) + 1.0)
        elif family == "bernoulli":
            y = 1 if rng.uniform() < inv_logit(eta) else 0
        elif family == "poisson":
            # math.exp raises OverflowError past about eta = 710 where R's exp() returns
            # Inf; the refusal has to be poisson_inv's, shared byte for byte with the R
            # twin, so the overflow becomes the infinite mean it stands for.
            try:
                lam = math.exp(eta)
            except OverflowError:
                lam = math.inf
            y = poisson_inv(lam, rng.uniform())
        elif family == "ordinal":
            y = ordinal_inv(eta, thresholds, rng.uniform())
        elif family == "beta":
            mu, phi = inv_logit(eta), resp.get("phi", 10.0)
            y = beta_draw(rng, mu * phi, (1.0 - mu) * phi)
        else:
            raise ValueError(f"unknown family: {family}")
        if ndp is not None and isinstance(y, float):
            y = round(y, ndp)

        out = {"subject": r["subject"]}
        if has_item:
            out["item"] = r["item"]
        for gname, (over, _cols, group_of) in group_meta.items():
            out[gname] = group_of[r[over]] + 1
        out.update(r["labels"])
        # The observed (contaminated) value goes into the data; the latent one, in cv, drove eta.
        for p in predictors:
            out[p["name"]] = r["obs"][p["name"]]
        out[yname] = y
        out_rows.append(out)

    return Dataset(out_cols, out_rows)
