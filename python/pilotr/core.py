"""Portable numerical core for pilotr.

Everything here is implemented to be *bit-identical* with the R port: a shared
L'Ecuyer (1988) combined LCG for uniforms, Wichura's AS 241 inverse-normal, a hand-rolled
Cholesky factorisation, and inverse-CDF transforms for the response families. No NumPy
randomness is used, so the stream cannot drift between languages.
"""

from __future__ import annotations
import math

# --- L'Ecuyer (1988) combined LCG -------------------------------------------------
# All products stay below 2**53, so integer arithmetic is exact in IEEE-754 doubles.

_M1 = 2147483563
_M2 = 2147483399
_A1 = 40014
_A2 = 40692


class RNG:
    """Combined LCG producing uniforms in (0, 1). See spec/SPEC.md for the contract."""

    __slots__ = ("s1", "s2")

    def __init__(self, seed: int):
        seed = int(abs(seed))
        self.s1 = 1 + (seed % (_M1 - 1))
        self.s2 = 1 + ((_A2 * self.s1) % (_M2 - 1))
        for _ in range(10):  # warm-up
            self.uniform()

    def uniform(self) -> float:
        self.s1 = (_A1 * self.s1) % _M1
        self.s2 = (_A2 * self.s2) % _M2
        d = self.s1 - self.s2
        if d < 1:
            d += _M1 - 1
        return d / _M1

    def normal(self) -> float:
        return as241(self.uniform())

    def normals(self, k: int) -> list[float]:
        return [self.normal() for _ in range(k)]


# --- Wichura (1988) Algorithm AS 241: inverse normal CDF (PPND16) -----------------

def as241(p: float) -> float:
    """Inverse standard-normal CDF; full double precision. Same algorithm as R's qnorm."""
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

def cholesky(cov: list[list[float]]) -> list[list[float]]:
    """Lower Cholesky factor L with L Lᵀ = cov (Cholesky–Banachiewicz)."""
    n = len(cov)
    L = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1):
            s = sum(L[i][k] * L[j][k] for k in range(j))
            if i == j:
                L[i][j] = math.sqrt(max(cov[i][i] - s, 0.0))
            else:
                L[i][j] = (cov[i][j] - s) / L[j][j] if L[j][j] != 0.0 else 0.0
    return L


def matvec(L: list[list[float]], z: list[float]) -> list[float]:
    return [sum(L[i][k] * z[k] for k in range(len(z))) for i in range(len(L))]


# --- Response transforms (inverse-CDF / link functions) ---------------------------

def inv_logit(x: float) -> float:
    if x >= 0:
        return 1.0 / (1.0 + math.exp(-x))
    e = math.exp(x)
    return e / (1.0 + e)


def poisson_inv(lam: float, u: float) -> int:
    """Inverse-CDF Poisson draw from a uniform u."""
    p = math.exp(-lam)
    cum = p
    k = 0
    while u > cum and k < 1_000_000:
        k += 1
        p *= lam / k
        cum += p
    return k


def ordinal_inv(eta: float, thresholds: list[float], u: float) -> int:
    """Cumulative-logit ordinal draw: P(Y<=k) = invlogit(theta_k - eta). 1-indexed category."""
    for k, th in enumerate(thresholds):
        if u <= inv_logit(th - eta):
            return k + 1
    return len(thresholds) + 1


def gamma_mt(rng, shape: float) -> float:
    """Marsaglia-Tsang Gamma(shape, scale=1) draw from the shared RNG (normals + uniforms in a
    rejection loop). Identical float ops in R and Python => identical accept/reject => parity."""
    if shape < 1.0:
        g = gamma_mt(rng, shape + 1.0)
        return g * rng.uniform() ** (1.0 / shape)
    d = shape - 1.0 / 3.0
    c = 1.0 / math.sqrt(9.0 * d)
    while True:
        x = rng.normal()
        v = (1.0 + c * x) ** 3
        if v <= 0.0:
            continue
        if math.log(rng.uniform()) < 0.5 * x * x + d - d * v + d * math.log(v):
            return d * v


def beta_draw(rng, a: float, b: float) -> float:
    """Beta(a, b) via two Gamma draws: X/(X+Y), X~Gamma(a), Y~Gamma(b)."""
    x = gamma_mt(rng, a)
    return x / (x + gamma_mt(rng, b))
