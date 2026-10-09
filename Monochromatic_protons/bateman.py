"""
bateman.py
==========
Solution of the Bateman equations for a linear decay chain

    dN_0/dt = -λ_0 N_0
    dN_k/dt = -λ_k N_k + λ_{k-1} N_{k-1}        k = 1 .. n-1

with N_0(0) = 1 and N_k(0) = 0 (Campana et al. 2026, Eq. A7).

The analytic solution (Eq. A8–A9),

    N_k(t) = (λ_0 ... λ_{k-1}) Σ_{j=0..k} exp(-λ_j t) / Π_{m=0..k, m≠j} (λ_m - λ_j),

is exact but loses precision when decay constants are close, or span many
orders of magnitude, or at times short compared with the half-lives, because
the terms of the sum cancel (by up to hundreds of orders of magnitude).
Here it is evaluated in double precision, and every value whose terms cancel
by more than CANCEL_LIMIT (or that is not finite, or not positive) is
recomputed with mpmath, with the working precision raised until two
successive evaluations agree to MP_RTOL.  Exactly equal decay constants
(different nuclides with the same tabulated half-life) are separated by a
relative 1e-12, which changes the result by a similar relative amount.

Only numpy and mpmath (a sympy dependency) are needed.

    from bateman import chain_populations, chain_activities
    N = chain_populations([l0, l1, l2], times)      # shape (3, len(times))
    A = chain_activities([l0, l1, l2], times)       # A_k = λ_k N_k  [1/s]
"""

from __future__ import annotations

import math

import numpy as np
import mpmath

CANCEL_LIMIT = 1e3     # recompute in multiple precision if Σ|terms| > CANCEL_LIMIT·|N|
MP_DPS = 40            # initial mpmath working precision (decimal digits)
MP_RTOL = 1e-13        # agreement required between successive precisions
MP_MAX_DPS = 5000
DEGENERATE_REL = 1e-12


def _separate(lams: np.ndarray) -> np.ndarray:
    """Separate exactly (or almost exactly) equal decay constants."""
    lams = lams.astype(float).copy()
    n = len(lams)
    for i in range(n):
        for j in range(i):
            if abs(lams[i] - lams[j]) <= DEGENERATE_REL * max(lams[i], lams[j]):
                lams[i] = lams[j] * (1.0 + 10 * DEGENERATE_REL * (i - j))
    return lams


def _mp_eval(lams: np.ndarray, k: int, t: float, dps: int) -> tuple:
    """(N_k(t), Σ|terms|) at the given precision."""
    with mpmath.workdps(dps):
        L = [mpmath.mpf(float(x)) for x in lams[:k + 1]]
        T = mpmath.mpf(float(t))
        pref = mpmath.fprod(L[:k]) if k else mpmath.mpf(1)
        terms = [mpmath.exp(-L[j] * T) / mpmath.fprod(L[m] - L[j] for m in range(k + 1) if m != j)
                 for j in range(k + 1)]
        return pref * mpmath.fsum(terms), pref * mpmath.fsum(abs(x) for x in terms)


def _mp_value(lams: np.ndarray, k: int, t: float) -> float:
    """N_k(t) evaluated with mpmath at adaptive precision."""
    dps = MP_DPS
    val, mag = _mp_eval(lams, k, t, dps)
    while True:
        # digits lost to cancellation, plus margin
        lost = int(mpmath.log10(mag / abs(val))) if val != 0 else dps
        new_dps = min(max(2 * dps, lost + 30), MP_MAX_DPS)
        new_val, mag = _mp_eval(lams, k, t, new_dps)
        if new_val != 0 and abs(new_val - val) <= MP_RTOL * abs(new_val):
            return float(new_val)
        if new_dps >= MP_MAX_DPS:
            raise ArithmeticError(f"Bateman sum did not converge (k={k}, t={t})")
        val, dps = new_val, new_dps


def chain_populations(lambdas, times) -> np.ndarray:
    """
    Normalised populations N_k(t) of a linear chain, N_0(0) = 1.

    Parameters
    ----------
    lambdas : decay constants [1/s] of the chain members, all > 0
    times   : evaluation times [s], >= 0

    Returns
    -------
    array (n, len(times))
    """
    lams = _separate(np.asarray(lambdas, dtype=float))
    ts = np.asarray(times, dtype=float)
    n = len(lams)
    if n == 0:
        return np.zeros((0, len(ts)))
    if np.any(lams <= 0):
        raise ValueError("all decay constants must be positive")

    out = np.zeros((n, len(ts)))
    expo = np.exp(-np.outer(lams, ts))                  # (n, T)
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        for k in range(n):
            pref = math.prod(lams[:k])
            coef = np.empty(k + 1)
            for j in range(k + 1):
                d = lams[:k + 1] - lams[j]
                d[j] = 1.0
                coef[j] = pref / np.prod(d)
            terms = coef[:, None] * expo[:k + 1]        # (k+1, T)
            val = terms.sum(axis=0)
            mag = np.abs(terms).sum(axis=0)
            bad = ~np.isfinite(val) | ~np.isfinite(mag) | (mag > CANCEL_LIMIT * np.abs(val)) \
                | ((val <= 0) & (mag > 0))
            if k == 0:
                bad[:] = False                          # pure exponential: always exact
            for i in np.flatnonzero(bad):
                val[i] = _mp_value(lams, k, ts[i])
            out[k] = np.maximum(val, 0.0)
    return out


def chain_activities(lambdas, times) -> np.ndarray:
    """Activities A_k(t) = λ_k N_k(t) [decays/s per initial root nucleus]."""
    lams = np.asarray(lambdas, dtype=float)
    return lams[:, None] * chain_populations(lams, times)
