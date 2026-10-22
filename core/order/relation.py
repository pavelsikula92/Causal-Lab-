"""
core.order.relation
===================

Causal relations computed from coordinates, never stored as an N x N
matrix or as bitsets (the 2.x bitsets need memory ~ N^2: unusable at
N = 128 000 in 3+1).

    FLAT        x < y  iff  t_y > t_x  and  (t_y - t_x)^2 >= |x_y - x_x|^2
                (Minkowski; also the relation of a conformal sprinkling,
                whose metric Pi^2 eta has the same light cones)
    H1          the light cones of ds^2 = Pi^2 dt^2 - Pi^-2 dx^2 (the
                isotropic weak-field form with PPN gamma = 1): coordinate
                light speed Pi^2, so  x < y  iff  t_y - t_x >= the optical
                length  int Pi^-2 dl  along the straight segment. A
                straight segment is exact to first order in the field
                (Fermat: the first-order change of the path does not change
                the arrival time), and it is how the relation is DEFINED
                for K.2 (E47)(c); the core test checks the flat limit.

Chains
------
`chain_dp` is the dynamic programme behind every longest or heaviest
chain in the lab: W(v) = w_v + max over u < v of W(u), processed in time
order in blocks; the maximum over all predecessors equals the maximum
over link predecessors because a heaviest chain with positive weights is
saturated. The same routine counts maximal chains (for uniform
tie-breaking) in log space. Memory is O(block * N), time O(N^2).
"""

from __future__ import annotations

import math
from typing import Callable, Optional, Tuple

import numpy as np


# =============================================================================
# Flat relation
# =============================================================================

def tau2_matrix(A: np.ndarray, B: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    For points A (m x D) and B (k x D): dt = t_B - t_A (m x k) and the
    squared proper time tau^2 = dt^2 - |dx|^2 (m x k). The spatial sum is
    accumulated one dimension at a time (no m x k x D temporary).
    """
    dt = B[None, :, 0] - A[:, None, 0]
    t2 = dt * dt
    for d in range(1, A.shape[1]):
        dd = B[None, :, d] - A[:, None, d]
        t2 -= dd * dd
    return dt, t2


def prec_matrix(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Boolean (m x k): A_i < B_j in the flat relation (strict)."""
    dt, t2 = tau2_matrix(A, B)
    return (dt > 0) & (t2 >= 0)


def prec_vector(x: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Boolean (k,): x < B_j (flat)."""
    dt = B[:, 0] - x[0]
    dx2 = ((B[:, 1:] - x[1:]) ** 2).sum(axis=1) if B.shape[1] > 1 else 0.0
    return (dt > 0) & (dt * dt >= dx2)


def succ_vector(x: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Boolean (k,): B_j < x (flat)."""
    dt = x[0] - B[:, 0]
    dx2 = ((B[:, 1:] - x[1:]) ** 2).sum(axis=1) if B.shape[1] > 1 else 0.0
    return (dt > 0) & (dt * dt >= dx2)


def interval_members(X: np.ndarray, p: int, q: int) -> np.ndarray:
    """Indices strictly between p and q (flat), using the time order."""
    lo, hi = p + 1, q
    if hi <= lo:
        return np.zeros(0, dtype=np.int64)
    B = X[lo:hi]
    m = prec_vector(X[p], B) & succ_vector(X[q], B)
    return lo + np.nonzero(m)[0]


# =============================================================================
# H1 relation (straight-line optical length)
# =============================================================================

_GL_NODES, _GL_WEIGHTS = np.polynomial.legendre.leggauss(16)


def optical_length(A: np.ndarray, B: np.ndarray, inv_speed: Callable[[np.ndarray], np.ndarray]
                   ) -> np.ndarray:
    """
    Optical length int (1/c(x)) dl along straight spatial segments from the
    spatial positions of A (m x n) to those of B (m x n), pairwise row by
    row (16-point Gauss-Legendre). inv_speed(P) takes an (..., n) array
    of spatial points and returns 1/c there.
    """
    a, b = A, B
    d = b - a
    L = np.sqrt((d ** 2).sum(axis=-1))
    s = 0.5 * (_GL_NODES + 1.0)                       # nodes on [0, 1]
    pts = a[..., None, :] + s[:, None] * d[..., None, :]
    val = inv_speed(pts)                              # (..., 16)
    return L * 0.5 * (val * _GL_WEIGHTS).sum(axis=-1)


def h1_prec_pairs(A: np.ndarray, B: np.ndarray, Pi: Callable[[np.ndarray], np.ndarray]
                  ) -> np.ndarray:
    """
    Row-by-row H1 relation for paired arrays A, B (m x D): A_i < B_i iff
    t_B - t_A >= int Pi^-2 dl. Pi takes spatial points (..., n).
    """
    inv = lambda P: Pi(P) ** -2.0
    dt = B[:, 0] - A[:, 0]
    opt = optical_length(A[:, 1:], B[:, 1:], inv)
    return (dt > 0) & (dt >= opt)


# =============================================================================
# Chains
# =============================================================================

def chain_dp(X: np.ndarray, w: Optional[np.ndarray] = None, *, start: Optional[int] = None,
             block: int = 256) -> dict:
    """
    Heaviest chains in the flat order of time-sorted points X.

    W[v]    the largest total weight of a chain ending at v (w[v] included);
            with w = 1 the number of elements of the longest chain ending
            at v (height + 1)
    pred[v] an argmax predecessor (-1 if none): one heaviest chain

    If `start` is given, only chains that begin at element `start` count;
    elements outside its future get W = -inf.
    Blocked: predecessors in earlier blocks by one (b0 x B) matrix, inside
    the block one element at a time.
    """
    N = len(X)
    w = np.ones(N) if w is None else np.asarray(w, float)
    W = np.full(N, -np.inf)
    pred = np.full(N, -1, dtype=np.int64)
    allowed = np.ones(N, dtype=bool)
    if start is not None:
        allowed[:] = False
        allowed[start] = True
        allowed[start + 1:] = prec_vector(X[start], X[start + 1:])
    for b0 in range(0, N, block):
        b1 = min(N, b0 + block)
        Xb = X[b0:b1]
        if b0 > 0:
            R = prec_matrix(X[:b0], Xb)
            arg = np.where(R, W[:b0, None], -np.inf).argmax(axis=0)
            cols = np.arange(b1 - b0)
            best = np.where(R[arg, cols], W[arg], -np.inf)
        else:
            best = np.full(b1 - b0, -np.inf)
            arg = np.full(b1 - b0, -1)
        Rin = prec_matrix(Xb, Xb)
        for k in range(b1 - b0):
            v = b0 + k
            if not allowed[v]:
                continue
            if start is not None and v == start:
                W[v] = w[v]
                continue
            bv, av = best[k], int(arg[k])
            inner = np.nonzero(Rin[:k, k])[0]
            if inner.size:
                wi = W[b0 + inner]
                j = int(np.argmax(wi))
                if wi[j] > bv:
                    bv, av = wi[j], b0 + int(inner[j])
            if np.isfinite(bv):
                W[v], pred[v] = bv + w[v], av
            elif start is None:
                W[v] = w[v]
    return {"W": W, "pred": pred}


def count_heaviest(X: np.ndarray, w: np.ndarray, W: np.ndarray, *, start: Optional[int] = None,
                   only: Optional[np.ndarray] = None, tol: float = 1e-9) -> np.ndarray:
    """
    log of the number of heaviest chains ending at each element, given W
    from chain_dp with the same w and start: the log-sum of the counts of
    the predecessors u with W[u] + w[v] = W[v] (relative tolerance tol).
    A chain start (no such predecessor, or v == start) counts 1.
    `only` (boolean mask) restricts the work to a subset that is closed
    under tight predecessors (e.g. the elements on heaviest chains to a
    fixed end); everything else gets -inf.
    """
    N = len(X)
    logc = np.full(N, -np.inf)
    todo = np.isfinite(W) if only is None else (np.asarray(only, bool) & np.isfinite(W))
    for v in np.nonzero(todo)[0]:
        if start is not None and v == start:
            logc[v] = 0.0
            continue
        target = W[v] - w[v]
        ok = np.zeros(v, dtype=bool)
        if v > 0:
            ok = succ_vector(X[v], X[:v]) & todo[:v] & \
                (np.abs(W[:v] - target) <= tol * max(1.0, abs(target)))
        if ok.any():
            logc[v] = float(np.logaddexp.reduce(logc[:v][ok]))
        elif start is None:
            logc[v] = 0.0
    return logc


def trace_chain(pred: np.ndarray, end: int) -> np.ndarray:
    """The chain ending at `end`, following pred back to its start."""
    out = [end]
    while pred[out[-1]] >= 0:
        out.append(int(pred[out[-1]]))
    return np.array(out[::-1], dtype=np.int64)


def reversed_points(X: np.ndarray) -> np.ndarray:
    """Time reversal t -> -t, re-sorted; returns (Xr, perm) with Xr = X'[perm]."""
    Xr = X.copy()
    Xr[:, 0] = -Xr[:, 0]
    perm = np.argsort(Xr[:, 0], kind="stable")
    return Xr[perm], perm


# =============================================================================
# Membership in the causal future of a SET (source slab -> screen), K.2 (E47)
# =============================================================================

def flat_future_of_set(S: np.ndarray, Q: np.ndarray, chunk: int = 2048) -> np.ndarray:
    """Boolean (len Q): Q_i is in the flat causal future of some point of S."""
    out = np.zeros(len(Q), dtype=bool)
    for c0 in range(0, len(Q), chunk):
        Qc = Q[c0:c0 + chunk]
        dt, t2 = tau2_matrix(S, Qc)
        out[c0:c0 + chunk] = ((dt > 0) & (t2 >= 0)).any(axis=0)
    return out


def h1_future_of_set(S: np.ndarray, Q: np.ndarray, Pi: Callable[[np.ndarray], np.ndarray],
                     margin: float, chunk: int = 512) -> np.ndarray:
    """
    Boolean (len Q): Q_i is in the H1 causal future of some point of S, i.e.
    t_Q - t_S >= int Pi^-2 dl along the straight segment, for some S.
    Pi <= 1 makes every H1 pair a flat pair with t_Q - t_S - |dx| in
    [0, margin] or larger; candidates are the flat pairs whose flat slack
    is at most `margin` (an upper bound of the optical excess), plus any
    pair with slack > margin, which is related for sure.
    """
    out = np.zeros(len(Q), dtype=bool)
    n = S.shape[1] - 1
    for c0 in range(0, len(Q), chunk):
        Qc = Q[c0:c0 + chunk]
        dt = Qc[None, :, 0] - S[:, None, 0]
        dx = np.zeros_like(dt)
        for d in range(1, n + 1):
            dd = Qc[None, :, d] - S[:, None, d]
            dx += dd * dd
        dx = np.sqrt(dx)
        slack = dt - dx
        sure = (slack > margin).any(axis=0)
        cand = (slack >= 0) & (slack <= margin)
        res = sure.copy()
        si, qi = np.nonzero(cand & ~sure[None, :])
        if si.size:
            ok = h1_prec_pairs(S[si], Qc[qi], Pi)
            hit = np.zeros(len(Qc), dtype=bool)
            np.logical_or.at(hit, qi[ok], True)
            res |= hit
        out[c0:c0 + chunk] = res
    return out
