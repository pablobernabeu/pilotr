"""Portable numerical core for pilotr.

Everything here is implemented to be bit-identical with the R port. The
components are a shared L'Ecuyer (1988) combined LCG for uniforms, Wichura's
AS 241 inverse-normal, a hand-rolled Cholesky factorisation, inverse-CDF
transforms for the response families and a port of R's round(). The routines
avoid NumPy randomness, so the random stream remains aligned across the two
languages.
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
    needs.
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

    A pivot of exactly zero is left alone. It is a useful case, and no kind of failure: it is
    what a slope with a standard deviation of zero produces, which is how a term is held fixed
    while the rest of the structure is kept intact. The tolerance absorbs the rounding of an
    exactly-singular matrix, so a term that is negative only numerically, with no true
    negativity behind it, is still clamped as before.
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
    """Logistic (inverse-logit) link, `1 / (1 + exp(-x))`.

    The two branches keep the exponent negative on either side of zero, so a large `|x|`
    underflows towards 0 or 1 where a single expression would overflow `math.exp`. The R twin,
    `.inv_logit`, branches at the same value, so both languages evaluate the same expression
    for a given `x`.

    Parameters
    ----------
    x : float
        A value on the logit scale.

    Returns
    -------
    float
        The corresponding probability, in (0, 1).
    """
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


# --- Rounding a response: R's round(), src/nmath/fround.c -------------------------
# log10(2) as Rmath.h writes it, the decimal digits a double always carries (float.h's DBL_DIG)
# and the largest power of ten a double holds (DBL_MAX_10_EXP).

_M_LOG10_2 = 0.301029995663981195213738894724493027
_DBL_DIG = 15
_MAX10E = 308


def _pow_di(x: float, n: int) -> float:
    """`x` to the integer power `n` by binary powering, as R's R_pow_di() computes it.

    10^n is exact up to n = 22 however it is computed. From n = 25, repeated multiplication can
    land on a different double, so the multiplications have to follow R's order.
    """
    xn = 1.0
    if n != 0:
        negative = n < 0
        if negative:
            n = -n
        while True:
            if n & 1:
                xn *= x
            n >>= 1
            if not n:
                break
            x *= x
        if negative:
            xn = 1.0 / xn
    return xn


def fround(x: float, digits) -> float:
    """Round `x` to `digits` decimal places exactly as R's ``round(x, digits)`` does.

    A port of R's ``fround()`` (src/nmath/fround.c), whose algorithm dates from R 4.0.0. It
    returns whichever of ``floor(x * 10^d) / 10^d`` and ``ceil(x * 10^d) / 10^d`` lies nearer to
    `x` when the distances are computed in double precision, with ``d = digits``. A tie goes to
    the candidate whose last digit is even. CPython's ``round()`` rounds the exact binary value
    of `x` instead, so near a decimal tie the two can take different candidates. For instance,
    ``round(1.05, 1)`` is 1 in R and 1.1 in CPython. See "Rounding" in spec/SPEC.md.

    Parameters
    ----------
    x : float
        The value to round.
    digits : int or float
        Decimal places. A value that is not whole counts as ``floor(digits + 0.5)``, as in R.

    Returns
    -------
    float
        The value R's ``round()`` returns, to the bit.
    """
    if math.isnan(x) or math.isnan(digits):
        return x + digits
    if math.isinf(x) or digits > _MAX10E + _DBL_DIG or x == 0.0:
        return x
    if digits < -_MAX10E:
        return 0.0
    if digits == 0:
        # C's nearbyint(): the nearest whole number, a tie going to the even one. CPython's
        # round() computes exactly that with 0 digits, and keeps the sign of a zero result.
        return round(x, 0)
    dig = math.floor(digits + 0.5)
    sign = 1.0
    if x < 0.0:
        sign, x = -1.0, -x
    # The first term approximates log10(x) from the binary exponent alone. frexp() gives
    # x = m * 2^e with m in [0.5, 1), so e - 1 is C's logb(x). A rounding position beyond the 15
    # significant digits a double always carries leaves x as it is.
    if _M_LOG10_2 * (0.5 + (math.frexp(x)[1] - 1)) + dig > _DBL_DIG:
        return sign * x
    if dig <= _MAX10E:
        pow10 = _pow_di(10.0, dig)
        x10 = x * pow10
        i10 = float(math.floor(x10))
        xd = i10 / pow10
        xu = float(math.ceil(x10)) / pow10
    else:
        # 10^dig would overflow, so R applies 10^308 and then the rest, in that order.
        p10 = _pow_di(10.0, dig - _MAX10E)
        pow10 = _pow_di(10.0, _MAX10E)
        x10 = (x * pow10) * p10
        i10 = float(math.floor(x10))
        xd = i10 / pow10 / p10
        xu = float(math.ceil(x10)) / pow10 / p10
    du, dd = xu - x, x - xd
    return sign * (xu if du < dd or (math.fmod(i10, 2.0) == 1.0 and du == dd) else xd)
