"""The analysis model a design specification implies: its modelling data and its formula.

A port of the R twin's autoformula.R. `model_data` adds to a simulated data set the columns the
model reads, and `model_formula` writes the model as the R twin writes it, term for term, in lme4's
notation. `power_mixed` fits the model the two describe, so it analyses every coefficient and every
declared random effect, where it used to fit one model written into its source whatever the
specification declared.

Interactions ``"a:b"`` become product columns named ``"a_b"``. Categorical factors become their
numeric contrast columns. The response is logged for the lognormal families, whose coefficients
are on the log scale. tools/parity/formula_cross.py checks every shipped example and parity case
against the R twin: the formula strings must be identical, and so must the contrast and product
columns, bit for bit.
"""

from __future__ import annotations
import math
from .simulate import Dataset, _as_spec

# The random-effects structures of model_formula(), named as in the R twin, whose default is the
# first.
_STRUCTURES = ("maximal", "zcp", "intercepts")

# The families whose coefficients are on the log scale of the response, which model_data() logs.
_LOG_FAMILIES = ("lognormal", "shifted_lognormal")

# Raised when a specification has neither a coefficient nor a random effect. The R twin's
# model_formula() then writes ".y ~ " with nothing after it, and as.formula() stops on it with a
# parse error.
_NO_MODEL = ("this specification has no fixed coefficients and no random effects, so it implies "
             "no model to fit")


def _us(key):
    """The analysis column of a coefficient or slope key: ``"a:b"`` becomes ``"a_b"``."""
    return key.replace(":", "_")


def _interaction_keys(spec):
    """Every interaction key the model needs a product column for.

    These are the keys of the fixed coefficients and of every grouping factor's random slopes that
    contain ``":"``, de-duplicated in first-seen order, as ``.interaction_keys()`` in the R twin
    takes them. A random slope on an interaction needs its column as much as a coefficient does.
    """
    keys = list(spec["fixed"]["coefficients"])
    for entry in (spec.get("random") or {}).values():
        keys += list(entry.get("slopes") or {})
    out = []
    for k in keys:
        if ":" in k and k not in out:
            out.append(k)
    return out


def _re_correlated(entry):
    """Whether a grouping factor's random effects are correlated, by the R twin's rule.

    An explicit ``correlated`` decides. Without it, the effects are correlated when the entry
    gives ``correlations``, which is when the simulator draws them correlated.
    """
    flag = entry.get("correlated")
    if flag is not None:
        return flag is True
    return bool(entry.get("correlations"))


def _model_terms(spec, structure="maximal"):
    """The terms of the model a validated specification implies.

    Returns ``(fixed, groups)``. `fixed` lists the analysis columns of the coefficients, in the
    specification's order. `groups` holds one ``(name, slopes, bar)`` triple per grouping factor,
    in the specification's order, where `slopes` lists the analysis columns of its random slopes
    and `bar` is ``"|"`` or ``"||"``, as `model_formula` writes them.
    """
    fixed = [_us(k) for k in spec["fixed"]["coefficients"]]
    groups = []
    for g, entry in (spec.get("random") or {}).items():
        slopes = ([] if structure == "intercepts"
                  else [_us(k) for k in (entry.get("slopes") or {})])
        # `|` only where the process correlates the terms, as in the R twin. The two simpler
        # structures estimate no correlation, and a group left with its intercept alone takes `|`,
        # since lme4 refuses `(1 || g)`.
        if structure == "maximal":
            bar = "||" if entry.get("slopes") and not _re_correlated(entry) else "|"
        else:
            bar = "||" if slopes else "|"
        groups.append((g, slopes, bar))
    return fixed, groups


def model_formula(spec, structure="maximal"):
    """The lme4 formula a specification implies, as the R twin's ``model_formula()`` writes it.

    The formula's response is ``.y``, the column `model_data` adds. Each coefficient is a fixed
    term, an interaction ``"a:b"`` as its product column ``a_b``. Each grouping factor has a random
    intercept and the random slopes the specification declares.

    The three structures are the usual steps of a rule that simplifies a model when its fit is
    singular or fails to converge (Barr et al., 2013, doi:10.1016/j.jml.2012.11.001; Matuschek et
    al., 2017, doi:10.1016/j.jml.2017.01.001). ``"maximal"`` has every random slope the
    specification declares, correlated with the intercept where the specification correlates them,
    and uncorrelated (``||``) where it does not. ``"zcp"``, the zero-correlation-parameter model,
    keeps those slopes and writes ``||`` for every grouping factor with a slope. A grouping factor
    without slopes keeps ``(1 | g)``, since lme4 refuses ``(1 || g)``. ``"intercepts"`` keeps a
    random intercept for each grouping factor and drops every slope. The R package's
    ``power_mixed()`` takes these names as the rungs of such a rule, which this package does not
    simulate.

    The string is the R twin's, character for character, which tools/parity/formula_cross.py
    checks on every shipped example and parity case. ``power_mixed`` fits the ``"maximal"`` model,
    with each random effect as an independent variance component.

    Parameters
    ----------
    spec : dict, str or os.PathLike
        A design specification (dict or path to a JSON file). It is validated first.
    structure : {"maximal", "zcp", "intercepts"}, optional
        The random-effects structure (default ``"maximal"``).

    Returns
    -------
    str
        The formula, such as ``".y ~ cond + (1 + cond | subject) + (1 + cond | item)"``.

    Raises
    ------
    ValueError
        If `structure` is not one of the three names, or if the specification has neither a
        coefficient nor a random effect, so that the formula would have no terms. The R twin then
        stops with a parse error.
    """
    if structure not in _STRUCTURES:
        raise ValueError("structure must be 'maximal', 'zcp' or 'intercepts'")
    fixed, groups = _model_terms(_as_spec(spec), structure)
    terms = fixed + ["(%s %s %s)" % (" + ".join(["1"] + slopes), bar, g)
                     for g, slopes, bar in groups]
    if not terms:
        raise ValueError(_NO_MODEL)
    return ".y ~ " + " + ".join(terms)


def _log(x):
    """R's log() of one double: -inf at 0 and nan below it, where math.log() raises.

    math.log() and R's log() call the platform's C library, which on Linux is the same one, so the
    two give the same double there. NumPy's log may take a vectorised path of its own.
    """
    if x > 0:
        return math.log(x)
    return -math.inf if x == 0 else math.nan


def model_data(spec, dataset):
    """The modelling data for a simulated data set, as the R twin's ``model_data()`` builds it.

    Three kinds of column are added to the data, in this order. ``.y`` is the analysis response:
    ``log(y)`` for the ``lognormal`` family, ``log(y - shift)`` for ``shifted_lognormal`` and the
    response itself for every other family. Each factor's contrast columns follow, read from its
    labels. Every factor's level index is taken before any contrast column is written, since a
    contrast column may carry its own factor's name and then replaces the labels. Last comes one
    product column per interaction key among the coefficients and the random slopes, ``a_b`` for
    ``"a:b"``, the product of its parts taken from left to right.

    A column that already exists is replaced where it stands, and a new one goes at the end, as in
    the R twin. The contrast and product columns equal the R twin's bit for bit, given the same
    data. ``.y`` passes through the C library's ``log()``, which can differ in the last place
    between R and Python on Windows.

    Parameters
    ----------
    spec : dict, str or os.PathLike
        A design specification (dict or path to a JSON file). It is validated first.
    dataset : Dataset or pandas.DataFrame
        A data set simulated from `spec`, as `simulate` returns it, or a DataFrame with its columns.
        A DataFrame is copied, so the one given is left as it was.

    Returns
    -------
    pandas.DataFrame
        The data with the ``.y``, contrast and product columns that `model_formula` reads.

    Raises
    ------
    TypeError
        If `dataset` is neither a Dataset nor a DataFrame.

    Notes
    -----
    Requires `pandas`, from the `mixed` extra.
    """
    return _model_data(_as_spec(spec), dataset)


def _model_data(spec, dataset):
    """`model_data` for a specification already validated and normalised, as in the replicates."""
    import pandas as pd

    if isinstance(dataset, Dataset):
        d = pd.DataFrame(dataset.rows, columns=dataset.columns)
    elif isinstance(dataset, pd.DataFrame):
        d = dataset.copy()
    else:
        raise TypeError("dataset must be a pilotr Dataset or a pandas DataFrame")
    resp = spec["response"]
    y = d[resp["name"]]
    if resp["family"] in _LOG_FAMILIES:
        shift = resp.get("shift", 0.0)
        d[".y"] = [_log(float(v) - shift) for v in y]
    else:
        d[".y"] = [float(v) for v in y]
    # The level index of every factor first, then the columns. Writing a contrast column named
    # after its own factor replaces the labels that the factor's later columns are matched against.
    factors = spec["factors"]
    index = []
    for f in factors:
        position = {level: i for i, level in enumerate(f["levels"])}
        index.append([position.get(v) for v in d[f["name"]]])
    for f, idx in zip(factors, index):
        for col, values in f["contrasts"].items():
            d[col] = [math.nan if i is None else float(values[i]) for i in idx]
    # A left fold, ((a * b) * c), as R's Reduce(`*`, ...) takes it. Each product of two doubles is
    # correctly rounded, so the two twins agree bit for bit.
    for key in _interaction_keys(spec):
        parts = key.split(":")
        product = d[parts[0]].to_numpy(dtype=float)
        for p in parts[1:]:
            product = product * d[p].to_numpy(dtype=float)
        d[_us(key)] = product
    return d
