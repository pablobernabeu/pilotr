"""Portable numerical core for pilotr.

Everything here is implemented to be bit-identical with the R port. The
components are a shared L'Ecuyer (1988) combined LCG for uniforms, Wichura's
AS 241 inverse-normal, a hand-rolled Cholesky factorisation, and inverse-CDF
transforms for the response families. The routines avoid NumPy randomness, so
the random stream remains aligned across the two languages.
"""

from __future__ import annotations
import math
import sys

# --- L'Ecuyer (1988) combined LCG -------------------------------------------------
# All products stay below 2**53, so integer arithmetic is exact in IEEE-754 doubles.

_M1 = 2147483563
_M2 = 2147483399
_A1 = 40014
_A2 = 40692


class RNG:
    """Shared cross-language random-number generator.

    The combined linear congruential generator of L'Ecuyer (1988), implemented identically
    in R and Python so that a given seed yields the same stream in both. All intermediate
    products stay below `2**53`, so the arithmetic is exact in IEEE-754 doubles. The
    draw-order contract is documented in the specification at
    https://github.com/pablobernabeu/pilotr/blob/main/spec/SPEC.md.

    Parameters
    ----------
    seed : int
        Seed for the generator (coerced to a non-negative integer).
    """

    __slots__ = ("s1", "s2")

    def __init__(self, seed: int):
        seed = int(abs(seed))
        self.s1 = 1 + (seed % (_M1 - 1))
        self.s2 = 1 + ((_A2 * self.s1) % (_M2 - 1))
        for _ in range(10):  # warm-up
            self.uniform()

    def uniform(self) -> float:
        """Return one draw from the standard uniform distribution on (0, 1)."""
        self.s1 = (_A1 * self.s1) % _M1
        self.s2 = (_A2 * self.s2) % _M2
        d = self.s1 - self.s2
        if d < 1:
            d += _M1 - 1
        return d / _M1

    def normal(self) -> float:
        """Return one draw from the standard normal distribution."""
        return as241(self.uniform())

    def normals(self, k: int) -> list[float]:
        """Return a list of `k` standard-normal draws."""
        return [self.normal() for _ in range(k)]


def replicate_seeds(base, n: int) -> list[int]:
    """The seeds pilotr's replicate loops give to their replicates, from a specification's seed.

    Until 0.3 the rule was ``base + i``. Consecutive seeds are not independent streams in this
    generator: seeding sets ``s1`` to ``1 + (seed mod 2147483562)`` and ``s2`` from ``s1``, and only
    ten warm-up draws are discarded, so replicate ``i`` and replicate ``i + 1`` begin a few steps
    apart in the same sequence rather than in unrelated parts of it. Measured over 2,000
    replicates, the first draw of replicate ``i`` correlated 0.95 with the first draw of ``i + 1``.

    An arithmetic scramble does not fix that, because the seeding rule is itself linear in the
    seed: adding a Weyl increment and applying a Lehmer step left the first draws correlated at
    -0.27. Drawing the seeds from the shared generator does work, since successive outputs of the
    combined generator are what it exists to make look independent. Duplicates are skipped, so no
    two replicates are handed the same seed and silently produce identical data, and the skipping
    is deterministic, so this matches the R implementation exactly.

    Parameters
    ----------
    base : int
        The specification's seed.
    n : int
        How many replicate seeds to return.

    Returns
    -------
    list of int
        `n` distinct seeds.
    """
    rng = RNG(base)
    out, seen = [], set()
    while len(out) < n:
        s = math.floor(rng.uniform() * (_M1 - 1)) + 1
        if s not in seen:
            seen.add(s)
            out.append(s)
    return out


# --- Wichura (1988) Algorithm AS 241: inverse normal CDF (PPND16) -----------------

def as241(p: float) -> float:
    """Inverse standard-normal CDF (quantile function), Wichura's (1988) Algorithm AS 241.

    The PPND16 routine underlying R's `qnorm`, accurate to full double precision.

    Parameters
    ----------
    p : float
        A probability in the open interval (0, 1).

    Returns
    -------
    float
        The standard-normal quantile for `p`.
    """
    q = p - 0.5
    if abs(q) <= 0.425:
        r = 0.180625 - q * q
        num = (((((((2509.0809287301226727 * r + 33430.575583588128105) * r +
                    67265.770927008700853) * r + 45921.953931549871457) * r +
                  13731.693765509461125) * r + 1971.5909503065514427) * r +
                133.14166789178437745) * r + 3.387132872796366608)
        den = (((((((5226.495278852854561 * r + 28729.085735721942674) * r +
                    39307.89580009271061) * r + 21213.794301586595867) * r +
                  5394.1960214247511077) * r + 687.1870074920579083) * r +
                42.313330701600911252) * r + 1.0)
        return q * num / den
    r = p if q < 0 else 1.0 - p
    r = math.sqrt(-math.log(r))
    if r <= 5.0:
        r -= 1.6
        num = (((((((7.7454501427834140764e-4 * r + 0.0227238449892691845833) * r +
                    0.24178072517745061177) * r + 1.27045825245236838258) * r +
                  3.64784832476320460504) * r + 5.7694972214606914055) * r +
                4.6303378461565452959) * r + 1.42343711074968357734)
        den = (((((((1.05075007164441684324e-9 * r + 5.475938084995344946e-4) * r +
                    0.0151986665636164571966) * r + 0.14810397642748007459) * r +
                  0.68976733498510000455) * r + 1.6763848301838038494) * r +
                2.05319162663775882187) * r + 1.0)
    else:
        r -= 5.0
        num = (((((((2.01033439929228813265e-7 * r + 2.71155556874348757815e-5) * r +
                    0.0012426609473880784386) * r + 0.026532189526576123093) * r +
                  0.29656057182850489123) * r + 1.7848265399172913358) * r +
                5.4637849111641143699) * r + 6.6579046435011037772)
        den = (((((((2.04426310338993978564e-15 * r + 1.4215117583164458887e-7) * r +
                    1.8463183175100546818e-5) * r + 7.868691311456132591e-4) * r +
                  0.0148753612908506148525) * r + 0.13692988092273580531) * r +
                0.59983220655588793769) * r + 1.0)
    z = num / den
    return -z if q < 0 else z


# --- Linear algebra (kept hand-rolled so rounding matches the R port exactly) -----

def dot(a, b, m: int) -> float:
    """Inner product of the first `m` elements, as an explicit IEEE-754 double fold.

    Neither language's built-in reduction is a plain double fold, and the two disagree with
    each other: CPython's `sum` has applied Neumaier compensation to float sequences since
    version 3.12, while base R's `sum` accumulates in 80-bit long double on x86. Both are more
    accurate than a naive fold, but they are more accurate in different ways, so an inner
    product of length three or more can land on different doubles in the two ports. Spelling
    the loop out keeps the arithmetic identical, which is what the cross-language guarantee
    actually needs.
    """
    s = 0.0
    for k in range(m):
        s += a[k] * b[k]
    return s


def _stop_not_pd(label, cols, i: int, d: float):
    """Report which grouping factor failed, which random-effect column the factorisation
    reached, and how negative the pivot was. The offending column is more actionable than the
    smallest eigenvalue of the whole matrix, because it points at the correlations the user has
    to change, and it costs no extra numerics, so the R and Python ports report it identically.
    """
    where = ((" at column '%s'" % cols[i]) if cols is not None and len(cols) > i
             else " at position %d" % (i + 1))
    involved = ("" if cols is None else " Check the correlations among %s."
                % ", ".join("'%s'" % c for c in cols[:i + 1]))
    raise ValueError(
        "the random-effect covariance for '%s' is not positive definite: the Cholesky "
        "factorisation failed%s (pivot %g).%s No random-effect distribution has the requested "
        "standard deviations and correlations."
        % (label if label is not None else "unknown group", where, d, involved))


def cholesky(cov: list[list[float]], label=None, cols=None) -> list[list[float]]:
    """Lower Cholesky factor L with L Lᵀ = cov (Cholesky–Banachiewicz).

    A negative pivot means the covariance implied by the specified standard deviations and
    correlations is not positive definite, so no random-effect distribution has those moments.
    Clamping the pivot at zero and carrying on silently did not fail, but it did not honour the
    specification either: the factor returned generated random effects whose standard
    deviations were several times the requested ones. A correlation set that a standard
    Cholesky routine rejects has to be an error, because the alternative is plausible-looking
    data from a process the user never described.

    A pivot of exactly zero is left alone. That is a genuinely useful case, not a failure: it is
    what a slope with a standard deviation of zero produces, which is how a term is held fixed
    while the rest of the structure is kept intact. The tolerance absorbs the rounding of an
    exactly-singular matrix, so a term that is only numerically rather than truly negative is
    clamped as before rather than rejected.
    """
    n = len(cov)
    L = [[0.0] * n for _ in range(n)]
    tol = sys.float_info.epsilon * n * max([cov[i][i] for i in range(n)] + [1.0])
    for i in range(n):
        for j in range(i + 1):
            s = dot(L[i], L[j], j)
            if i == j:
                d = cov[i][i] - s
                if d < -tol:
                    _stop_not_pd(label, cols, i, d)
                L[i][j] = math.sqrt(max(d, 0.0))
            else:
                L[i][j] = (cov[i][j] - s) / L[j][j] if L[j][j] != 0.0 else 0.0
    return L


def matvec(L: list[list[float]], z: list[float]) -> list[float]:
    return [dot(L[i], z, len(z)) for i in range(len(L))]


# --- Response transforms (inverse-CDF / link functions) ---------------------------

def inv_logit(x: float) -> float:
    """Logistic (inverse-logit) link `1 / (1 + exp(-x))`, evaluated stably for large `|x|`."""
    if x >= 0:
        return 1.0 / (1.0 + math.exp(-x))
    e = math.exp(x)
    return e / (1.0 + e)


def _stop_poisson_mean(lam: float):
    """Refuse a mean the inverse-CDF walk cannot serve. exp(-lam) underflows to exactly
    zero near lam = 746, after which every term of the cumulative sum is zero and the walk
    would return its iteration cap as though it were a drawn count (a poisson intercept of
    7 implies lam = exp(7), about 1097, already past that point); a huge but representable
    mean exhausts the cap the same way. Both twins raise this text byte for byte; Python
    renders an infinite value as "inf" where R prints "Inf", hence the normalisation.
    """
    raise ValueError(
        "the poisson mean exp(eta) = %s is too large for the inverse-CDF sampler: the "
        "cumulative distribution cannot reach the drawn uniform, so no count can be drawn. "
        "Lower the poisson intercept or coefficients until the implied mean is simulable."
        % ("Inf" if math.isinf(lam) else "%g" % lam))


def poisson_inv(lam: float, u: float) -> int:
    """Inverse-CDF Poisson draw from a uniform u."""
    p = math.exp(-lam)
    if p == 0.0:
        _stop_poisson_mean(lam)
    cum = p
    k = 0
    while u > cum and k < 1_000_000:
        k += 1
        p *= lam / k
        cum += p
    if u > cum:
        _stop_poisson_mean(lam)
    return k


def ordinal_inv(eta: float, thresholds: list[float], u: float) -> int:
    """Cumulative-logit ordinal draw: P(Y<=k) = invlogit(theta_k - eta). 1-indexed category."""
    for k, th in enumerate(thresholds):
        if u <= inv_logit(th - eta):
            return k + 1
    return len(thresholds) + 1


def gamma_mt(rng, shape: float) -> float:
    """Marsaglia-Tsang Gamma(shape, scale=1) draw from the shared RNG (normals and uniforms in a
    rejection loop). The floating-point operations are identical in R and Python, so the
    accept/reject decisions also agree, which preserves parity across the two languages.

    The cube is written as two multiplications rather than ``** 3``. Python's ``**`` calls the
    library pow(), whereas R's ``^`` special-cases small integer exponents into repeated
    multiplication; measured over 200,000 draws in this generator's range the two disagree on
    a third of inputs, by up to 6 ulp. Two explicit multiplies are exactly rounded in both
    languages, and they also decide the rejection step, so aligning them keeps the
    accept/reject sequence (and hence the number of draws consumed) identical.
    """
    if shape < 1.0:
        g = gamma_mt(rng, shape + 1.0)
        return g * rng.uniform() ** (1.0 / shape)
    d = shape - 1.0 / 3.0
    c = 1.0 / math.sqrt(9.0 * d)
    while True:
        x = rng.normal()
        t = 1.0 + c * x
        v = t * t * t
        if v <= 0.0:
            continue
        if math.log(rng.uniform()) < 0.5 * x * x + d - d * v + d * math.log(v):
            return d * v


def beta_draw(rng, a: float, b: float) -> float:
    """Beta(a, b) via two Gamma draws: X/(X+Y), X~Gamma(a), Y~Gamma(b)."""
    x = gamma_mt(rng, a)
    return x / (x + gamma_mt(rng, b))
