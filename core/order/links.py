"""
core.order.links
================

Links (the Hasse diagram, the transitive reduction) of a flat sprinkling,
from coordinates, without the N x N relation.

The near-cone shell
-------------------
k < j is a link iff the Alexandrov interval I(k, j) contains no point.
If rho is the density, an interval of volume V is empty with probability
exp(-rho V). Every candidate link therefore lies (with probability
1 - exp(-lam) per pair) in the SHELL

    C_j = { k < j : rho_min * V(k, j) <= lam }

hugging the past light cone of j. Every point inside I(k, j) for k in C_j
is itself in C_j (its interval to j is contained in I(k, j)), so the links
of j are exactly the MAXIMAL elements of C_j: no candidate needs a search
outside the shell. lam = 25 by default: an interval of rho V > 25 is
empty with probability < 1.4e-11; the expected number of links missed
in a whole sprinkling of 1e5 points is far below one (the core test
compares with brute force on every pair).

For a non-uniform (conformal) sprinkling pass rho_min, the smallest
density in the region; the shell is then wider, never too narrow.

Maximal elements, vectorized
----------------------------
Candidates are taken in decreasing time. k is maximal iff it precedes no
other candidate; it suffices to test it against the maximal elements
found so far (any later candidate above k lies below one of them) and
against the other candidates of its own block. Blocks of `block`
candidates make this a few matrix operations per element.

Output
------
CSR arrays of PARENTS (past links): parents of j are
idx[ptr[j]:ptr[j+1]], sorted increasingly.
"""

from __future__ import annotations

import math
from typing import Optional, Tuple

import numpy as np

from .relation import prec_matrix, tau2_matrix
from .sprinkle import interval_volume


def _tau_cut(D: int, rho_min: float, lam: float) -> float:
    """Proper time at which rho_min * V(tau) = lam."""
    unit = interval_volume(D, 1.0)
    return (lam / (rho_min * unit)) ** (1.0 / D)


def _maximal(X: np.ndarray, cand: np.ndarray, block: int) -> np.ndarray:
    """
    Maximal elements of the candidate set `cand` (indices, any order).
    Candidates are taken by decreasing time in blocks; a candidate is dead
    if it precedes a maximal element already found or another candidate of
    its own block (see the module notes for why this is complete).
    """
    cand = cand[np.argsort(-X[cand, 0], kind="stable")]
    D = X.shape[1]
    M = np.zeros((0, D))
    keep_all = []
    for b0 in range(0, cand.size, block):
        cb = cand[b0:b0 + block]
        Xc = X[cb]
        dead = prec_matrix(Xc, Xc).any(axis=1)
        if M.shape[0]:
            dead |= prec_matrix(Xc, M).any(axis=1)
        keep = cb[~dead]
        if keep.size:
            keep_all.append(keep)
            M = np.vstack([M, X[keep]])
    return np.concatenate(keep_all) if keep_all else np.zeros(0, dtype=np.int64)


def links_csr(X: np.ndarray, rho_min: float, *, lam: float = 25.0, rows: int = 64,
              block: int = 128, targets: Optional[np.ndarray] = None,
              progress=None) -> Tuple[np.ndarray, np.ndarray]:
    """
    Past links of every point of the time-sorted flat sprinkling X.
    Returns (ptr, idx) in CSR form. rows = how many j are scanned at once
    against all earlier points (memory ~ rows * N * 8 bytes * 3).
    targets (optional): only these elements get their links computed; the
    others get none (e.g. the link graph of a time band).
    """
    N, D = X.shape
    tc2 = _tau_cut(D, rho_min, lam) ** 2
    ptr = np.zeros(N + 1, dtype=np.int64)
    chunks = []
    total = 0
    want = None
    if targets is not None:
        want = np.zeros(N, dtype=bool)
        want[np.asarray(targets)] = True
    for j0 in range(0, N, rows):
        j1 = min(N, j0 + rows)
        if want is not None and not want[j0:j1].any():
            ptr[j0 + 1:j1 + 1] = total
            continue
        if j1 <= 1:
            for j in range(j0, j1):
                ptr[j + 1] = total
            continue
        dt, t2 = tau2_matrix(X[:j1], X[j0:j1])            # (j1 x R)
        shell = (dt > 0) & (t2 >= 0) & (t2 <= tc2)
        for r, j in enumerate(range(j0, j1)):
            cand = np.nonzero(shell[:j, r])[0] if (want is None or want[j]) else np.zeros(0, dtype=np.int64)
            if cand.size == 0:
                ptr[j + 1] = total
                continue
            maxi = _maximal(X, cand, block)
            par = np.sort(maxi)
            chunks.append(par)
            total += par.size
            ptr[j + 1] = total
        if progress and (j0 // rows) % 50 == 0:
            progress(f"links {j1}/{N}")
    idx = np.concatenate(chunks) if chunks else np.zeros(0, dtype=np.int64)
    return ptr, idx.astype(np.int64)


def links_bruteforce(X: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Reference implementation (small N only): k < j is a link iff no m with
    k < m < j, tested on the full relation matrix.
    """
    N = len(X)
    R = prec_matrix(X, X)
    ptr = np.zeros(N + 1, dtype=np.int64)
    out = []
    for j in range(N):
        anc = np.nonzero(R[:, j])[0]
        if anc.size:
            sub = R[np.ix_(anc, anc)]
            par = anc[~sub.any(axis=1)]
        else:
            par = anc
        out.append(np.sort(par))
        ptr[j + 1] = ptr[j] + par.size
    idx = np.concatenate(out) if out else np.zeros(0, dtype=np.int64)
    return ptr, idx.astype(np.int64)


def children_counts(ptr: np.ndarray, idx: np.ndarray, N: int) -> np.ndarray:
    """Number of future links (children) of every element."""
    return np.bincount(idx, minlength=N).astype(np.int64)


def heights_from_links(ptr: np.ndarray, idx: np.ndarray) -> np.ndarray:
    """Height (longest chain below, 0 for minimal elements) from the links."""
    N = len(ptr) - 1
    h = np.zeros(N, dtype=np.int64)
    for j in range(N):
        a, b = ptr[j], ptr[j + 1]
        if b > a:
            h[j] = 1 + h[idx[a:b]].max()
    return h
