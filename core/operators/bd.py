"""
core.operators.bd
=================

The discrete smeared Benincasa-Dowker (BD) d'Alembertian on a sprinkled
causal order, in its static reduction (field equation A5'). Ported from
GTS Lab 3.x (fields.laplacians, BD part only; the slice and kernel
Laplacians of 3.x are not carried over).

    B_eps phi(x) = -a_d (eps rho)^(2/d) [ phi(x) + sum_{y < x} K_eps(n_xy) phi(y) ]

K_eps(n) = eps (1-eps)^n sum_i C_i binom(n, i) (eps/(1-eps))^i with the
action coefficients C (d = 2: -2, 4, -2; d = 3: -1, 27/8, -9/4;
d = 4: -1, 9, -16, 8), n_xy = number of elements strictly between y and x.
Its expectation tends to (-d_t^2 + lap) phi, so for a STATIC field it is
lap phi. The nonlocality scale is l_k = (eps rho)^(-1/d).

Normalisation a_d: 2 (d = 2), 4/sqrt(6) (d = 4), and
a_3 = (pi / (3 sqrt 2))^(2/3) / Gamma(5/3) = 0.9067 (d = 3); all three were
checked by integrating the expectation of the operator on phi = t^2 at the
tip of a causal diamond (GTS Lab 3.x).

Static reduction: the field is represented on a radial basis (hat
functions on nodes r_k), so every sampled element x gives one linear row
in the nodal values: (B Pi)(x) = sum_k A[x, k] Pi_k. Only interval counts
(order) and the radii of the elements (where the static field is
evaluated) enter.

Null option (test F1.2): with `null_rng`, a second set of rows is built in
which the kernel weights K_eps(n_xy) of each row are randomly permuted
among that row's candidate ancestors. The weights keep their distribution,
but the link between the causal structure (which y has which interval
count) and the field is destroyed.
"""

from __future__ import annotations

import math
from typing import Any, Dict, Optional, Tuple

import numpy as np

from ..order.relation import prec_matrix, tau2_matrix
from ..order.sprinkle import interval_volume

BD_COEFFS = {2: (-2.0, 4.0, -2.0), 3: (-1.0, 27.0 / 8.0, -9.0 / 4.0), 4: (-1.0, 9.0, -16.0, 8.0)}
BD_NORM = {2: 2.0, 3: (math.pi / (3.0 * math.sqrt(2.0))) ** (2.0 / 3.0) / math.gamma(5.0 / 3.0),
           4: 4.0 / math.sqrt(6.0)}


def smeared_kernel(n: np.ndarray, d: int, eps: float) -> np.ndarray:
    """K_eps(n) for an array of interval counts (vectorised)."""
    C = BD_COEFFS[d]
    n = np.asarray(n, float)
    if eps >= 1.0:
        out = np.zeros_like(n)
        for i, c in enumerate(C):
            out += np.where(n == i, c, 0.0)
        return out
    q, r = 1.0 - eps, eps / (1.0 - eps)
    s = np.zeros_like(n)
    binom = np.ones_like(n)                 # binom(n, i), built up from i = 0
    for i, c in enumerate(C):
        if i > 0:
            binom = binom * (n - i + 1) / i
        s += c * np.where(n >= i, binom, 0.0) * r ** i
    return eps * np.exp(n * math.log(q)) * s


def kernel_cutoff_n(eps: float, d: int, tiny: float = 1e-13) -> int:
    """Smallest n beyond which |K_eps| stays below `tiny` (scanned)."""
    n = np.arange(0, 200000)
    k = np.abs(smeared_kernel(n, d, eps))
    big = np.nonzero(k > tiny)[0]
    return int(big[-1]) + 1 if big.size else 1


def hat_weights(r: np.ndarray, nodes: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Piecewise-linear interpolation: (node index lo, weight lo, weight hi) for every r."""
    K = len(nodes) - 1
    k = np.clip(np.searchsorted(nodes, r, side="right") - 1, 0, K - 1)
    f = np.clip((r - nodes[k]) / (nodes[k + 1] - nodes[k]), 0.0, 1.0)
    return k, 1.0 - f, f


def bd_static_rows(X: np.ndarray, rho: float, eps: float, rows_idx: np.ndarray, nodes: np.ndarray,
                   R: float, lam_extra: float = 50.0,
                   null_rng: Optional[np.random.Generator] = None) -> Dict[str, Any]:
    """
    For each sampled element x (indices rows_idx into the time-sorted X), the
    coefficients of the static BD operator on the radial basis:
        (B_eps Pi)(x) = sum_k A[x, k] Pi_k
    (elements beyond the outer node are clamped to it). Also the operator
    applied to the static quadratic centred at x, |y - x|^2, whose local
    Laplacian is 2n (reported). With null_rng, also "A_null": the same rows
    with the kernel weights permuted among each row's candidates (module notes).
    """
    N, D = X.shape
    pref = -BD_NORM[D] * (eps * rho) ** (2.0 / D)
    nmax = kernel_cutoff_n(eps, D)
    lam = 1.5 * nmax + lam_extra              # rho V cutoff for candidates
    tc2 = (lam / (rho * interval_volume(D, 1.0))) ** (2.0 / D)
    K = len(nodes)
    A = np.zeros((len(rows_idx), K))
    A_null = np.zeros((len(rows_idx), K)) if null_rng is not None else None
    r2_val = np.zeros(len(rows_idx))
    ncand = []
    rr = np.sqrt((X[:, 1:] ** 2).sum(axis=1))

    def basis_row(i, cand, kw):
        coef = np.zeros(K)
        for idxs, wts in ((np.array([i]), np.array([1.0])), (cand, kw)):
            if idxs.size == 0:
                continue
            k, w0, w1 = hat_weights(np.minimum(rr[idxs], nodes[-1]), nodes)
            np.add.at(coef, k, wts * w0)
            np.add.at(coef, k + 1, wts * w1)
        return pref * coef

    for row, i in enumerate(rows_idx):
        dt, t2 = tau2_matrix(X[:i], X[i:i + 1])
        cand = np.nonzero((dt[:, 0] > 0) & (t2[:, 0] >= 0) & (t2[:, 0] <= tc2))[0]
        ncand.append(cand.size)
        if cand.size:
            Xc = X[cand]
            n_between = np.zeros(cand.size)
            for c0 in range(0, cand.size, 2048):
                n_between[c0:c0 + 2048] = prec_matrix(Xc[c0:c0 + 2048], Xc).sum(axis=1)
            kw = smeared_kernel(n_between, D, eps)
        else:
            kw = np.zeros(0)
        A[row] = basis_row(i, cand, kw)
        if A_null is not None:
            A_null[row] = basis_row(i, cand, null_rng.permutation(kw) if kw.size else kw)
        if cand.size:
            dq = ((X[cand, 1:] - X[i, 1:]) ** 2).sum(axis=1)
            r2_val[row] = pref * float((kw * dq).sum())
    out = {"A": A, "lap_r2": r2_val, "candidates_mean": float(np.mean(ncand)) if ncand else 0.0, "n_cut": nmax}
    if A_null is not None:
        out["A_null"] = A_null
    return out


def stratified_rows(X: np.ndarray, r: np.ndarray, t0: float, t1: float, nodes: np.ndarray,
                    r_rows: float, per_element: int, rng: np.random.Generator) -> np.ndarray:
    """
    Elements of the time band [t0, t1], per_element of them (at random) in
    every radial element [nodes[k], nodes[k+1]) below r_rows (all of them
    if fewer are available). A fixed radial layout gives a fixed source charge.
    """
    band = (X[:, 0] >= t0) & (X[:, 0] <= t1)
    out = []
    for lo, hi in zip(nodes[:-1], nodes[1:]):
        if lo >= r_rows:
            break
        c = np.nonzero(band & (r >= lo) & (r < min(hi, r_rows)))[0]
        if c.size:
            out.extend(rng.choice(c, size=min(per_element, c.size), replace=False).tolist())
    return np.sort(np.array(out, dtype=np.int64))
