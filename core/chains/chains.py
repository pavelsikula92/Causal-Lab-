"""
core.chains.chains
==================

The discrete motion of a body between two fixed events p < q: the chain
that maximises W = sum_v w_v (w = Pi for a body in a clock-rate field on
a flat order; w = 1, the longest chain, on a conformal sprinkling).

Ties
----
On a conformal sprinkling (w = 1) the longest chain is far from unique.
The 2.x lab picked the straightest one, which biases toward a straight
path. Here the reading is the AVERAGE over ALL heaviest chains, uniform:
forward counts from p and backward counts from q (log space), and the
mid-time position is averaged over the consecutive pairs (u, v) of
heaviest chains that straddle the middle time, each weighted by
(#chains p..u) (#chains v..q) / (#chains p..q). With continuous weights
the heaviest chain is unique and this is just that chain.

Null: the SHORTEST maximal chain (fewest links, F0.2 (E44a)), averaged the same
way over all shortest chains, by breadth-first search on the links.
"""

from __future__ import annotations

from typing import Any, Dict, Optional, Sequence, Tuple

import numpy as np

from ..order.relation import chain_dp, count_heaviest, reversed_points, prec_vector


def _mid_position(xu: np.ndarray, xv: np.ndarray, tm: float) -> np.ndarray:
    """Spatial position at time tm on the straight segment between events xu and xv."""
    f = (tm - xu[..., 0]) / (xv[..., 0] - xu[..., 0])
    return xu[..., 1:] + f[..., None] * (xv[..., 1:] - xu[..., 1:])


def heaviest_between(X: np.ndarray, w: np.ndarray, p: int, q: int, e_hat: Sequence[float],
                     tol: float = 1e-9) -> Dict[str, Any]:
    """
    Heaviest chains from p to q in the flat order of time-sorted X.
    Returns W_total, the number of heaviest chains (log), the number of
    elements on them, and `disp`: the mid-time displacement from the
    straight line along e_hat, averaged uniformly over all heaviest chains.
    """
    N = len(X)
    w = np.asarray(w, float)
    e = np.asarray(e_hat, float)
    fw = chain_dp(X, w, start=p)
    Wf = fw["W"]
    Xr, perm = reversed_points(X)
    inv = np.empty(N, dtype=np.int64)
    inv[perm] = np.arange(N)
    bw = chain_dp(Xr, w[perm], start=int(inv[q]))
    Wb = np.full(N, -np.inf)
    Wb[perm] = bw["W"]
    Wtot = Wf[q]
    if not np.isfinite(Wtot):
        return {"W_total": float("nan"), "disp": float("nan"), "n_opt": 0, "log_count": float("nan")}
    thr = tol * max(1.0, abs(Wtot))
    opt = np.isfinite(Wf) & np.isfinite(Wb) & (np.abs(Wf + Wb - w - Wtot) <= thr)
    lcf = count_heaviest(X, w, Wf, start=p, only=opt, tol=tol)
    lcbr = count_heaviest(Xr, w[perm], bw["W"], start=int(inv[q]), only=opt[perm], tol=tol)
    lcb = np.full(N, -np.inf)
    lcb[perm] = lcbr
    ltot = lcf[q]
    tm = 0.5 * (X[p, 0] + X[q, 0])
    x0 = 0.5 * (X[p, 1:] + X[q, 1:])
    U = np.nonzero(opt & (X[:, 0] <= tm))[0]
    V = np.nonzero(opt & (X[:, 0] > tm))[0]
    num, den = 0.0, 0.0
    for u in U:
        rel = prec_vector(X[u], X[V])
        tight = rel & (np.abs(Wf[u] + Wb[V] - Wtot) <= thr)
        if not tight.any():
            continue
        vv = V[tight]
        lw = lcf[u] + lcb[vv] - ltot
        wt = np.exp(lw)
        pos = _mid_position(np.repeat(X[u][None, :], vv.size, axis=0), X[vv], tm)
        num += float((wt * ((pos - x0) @ e)).sum())
        den += float(wt.sum())
    disp = num / den if den > 0 else float("nan")
    return {"W_total": float(Wtot), "disp": disp, "n_opt": int(opt.sum()),
            "log_count": float(ltot), "weight_check": den}


def shortest_between(X: np.ndarray, ptr: np.ndarray, idx: np.ndarray, p: int, q: int,
                     e_hat: Sequence[float]) -> Dict[str, Any]:
    """
    Shortest maximal chains (fewest links) from p to q by breadth-first
    search on the links (parents CSR ptr/idx), averaged uniformly like
    heaviest_between.
    """
    N = len(X)
    e = np.asarray(e_hat, float)
    # children lists from the parents CSR
    v_of = np.repeat(np.arange(N), np.diff(ptr))
    order = np.argsort(idx, kind="stable")
    cptr = np.zeros(N + 1, dtype=np.int64)
    np.add.at(cptr, idx + 1, 1)
    cptr = np.cumsum(cptr)
    cidx = v_of[order]

    def bfs(start, nbr_ptr, nbr_idx):
        d = np.full(N, -1, dtype=np.int64)
        lc = np.full(N, -np.inf)
        d[start] = 0
        lc[start] = 0.0
        frontier = [start]
        while frontier:
            nxt = {}
            for a in frontier:
                for b in nbr_idx[nbr_ptr[a]:nbr_ptr[a + 1]]:
                    b = int(b)
                    if d[b] == -1 or d[b] == d[a] + 1:
                        if d[b] == -1:
                            d[b] = d[a] + 1
                            nxt[b] = True
                        lc[b] = np.logaddexp(lc[b], lc[a])
            frontier = list(nxt.keys())
        return d, lc

    df, lcf = bfs(p, cptr, cidx)
    db, lcb = bfs(q, ptr, idx)
    if df[q] < 0:
        return {"length": -1, "disp": float("nan")}
    L = df[q]
    tm = 0.5 * (X[p, 0] + X[q, 0])
    x0 = 0.5 * (X[p, 1:] + X[q, 1:])
    u_all = idx
    v_all = v_of
    ok = (df[u_all] >= 0) & (db[v_all] >= 0) & (df[u_all] + 1 + db[v_all] == L) & \
         (X[u_all, 0] <= tm) & (X[v_all, 0] > tm)
    uu, vv = u_all[ok], v_all[ok]
    if uu.size == 0:
        return {"length": int(L), "disp": float("nan")}
    lw = lcf[uu] + lcb[vv] - lcf[q]
    wt = np.exp(lw)
    pos = _mid_position(X[uu], X[vv], tm)
    disp = float((wt * ((pos - x0) @ e)).sum() / wt.sum())
    return {"length": int(L), "disp": disp, "weight_check": float(wt.sum())}


def _children_csr(ptr: np.ndarray, idx: np.ndarray, N: int):
    """Children lists (CSR) from the parents CSR ptr/idx."""
    v_of = np.repeat(np.arange(N), np.diff(ptr))
    order = np.argsort(idx, kind="stable")
    cptr = np.zeros(N + 1, dtype=np.int64)
    np.add.at(cptr, idx + 1, 1)
    cptr = np.cumsum(cptr)
    return cptr, v_of[order], v_of


def allpaths_between(X: np.ndarray, ptr: np.ndarray, idx: np.ndarray, p: int, q: int,
                     e_hat: Sequence[float]) -> Dict[str, Any]:
    """
    F2.4 (E46b) functional U: the uniform average over ALL maximal chains from p to
    q, i.e. over all paths p -> q along links (a chain is maximal in [p, q]
    exactly when consecutive elements are links). No selection principle.

    Exact, no sampling: lcf[v] = log #(link paths p -> v) forward in time
    order, lcb[v] = log #(link paths v -> q) backward; every path crosses
    the mid time t_m on exactly one link (u, v) with t_u <= t_m < t_v, so
    the weight of that link is #(p -> u) #(v -> q) / #(p -> q) and the
    weights sum to 1 (weight_check). X must be time-sorted; ptr/idx are
    the parents CSR of links_csr.
    """
    N = len(X)
    e = np.asarray(e_hat, float)
    cptr, cidx, v_of = _children_csr(ptr, idx, N)
    lcf = np.full(N, -np.inf)
    lcf[p] = 0.0
    for v in range(p + 1, N):                              # time order = linear extension
        par = idx[ptr[v]:ptr[v + 1]]
        if par.size:
            vals = lcf[par]
            if np.isfinite(vals).any():
                lcf[v] = float(np.logaddexp.reduce(vals))
    lcb = np.full(N, -np.inf)
    lcb[q] = 0.0
    for u in range(q - 1, -1, -1):
        ch = cidx[cptr[u]:cptr[u + 1]]
        if ch.size:
            vals = lcb[ch]
            if np.isfinite(vals).any():
                lcb[u] = float(np.logaddexp.reduce(vals))
    if not np.isfinite(lcf[q]):
        return {"disp": float("nan"), "log_count": float("nan"), "weight_check": 0.0}
    tm = 0.5 * (X[p, 0] + X[q, 0])
    x0 = 0.5 * (X[p, 1:] + X[q, 1:])
    uu, vv = idx, v_of
    ok = np.isfinite(lcf[uu]) & np.isfinite(lcb[vv]) & (X[uu, 0] <= tm) & (X[vv, 0] > tm)
    uu, vv = uu[ok], vv[ok]
    wt = np.exp(lcf[uu] + lcb[vv] - lcf[q])
    pos = _mid_position(X[uu], X[vv], tm)
    disp = float((wt * ((pos - x0) @ e)).sum() / wt.sum()) if wt.size else float("nan")
    return {"disp": disp, "log_count": float(lcf[q]), "weight_check": float(wt.sum())}


def greedy_between(X: np.ndarray, ptr: np.ndarray, idx: np.ndarray, w: np.ndarray, p: int, q: int,
                   e_hat: Sequence[float]) -> Dict[str, Any]:
    """
    F2.4 (E46b) diagnostic G (REPORTED, no decision role): a local greedy walk from
    p along links. At element v the next element is the child c of v with
    c = q or c < q (inside the interval) that has the largest weight w[c];
    ties go to the lowest index (earliest in time). Stops at q. The mid-time
    position is read on the link that crosses t_m.
    """
    N = len(X)
    e = np.asarray(e_hat, float)
    cptr, cidx, _ = _children_csr(ptr, idx, N)
    tm = 0.5 * (X[p, 0] + X[q, 0])
    x0 = 0.5 * (X[p, 1:] + X[q, 1:])
    v, steps, disp = p, 0, float("nan")
    while v != q and steps < N:
        ch = cidx[cptr[v]:cptr[v + 1]]
        if ch.size == 0:
            break
        inside = (ch == q) | _prec_to(X, ch, q)
        ch = ch[inside]
        if ch.size == 0:
            break
        best = ch[np.argmax(w[ch])]          # argmax returns the first (lowest index) on ties
        if X[v, 0] <= tm < X[best, 0]:
            pos = _mid_position(X[v], X[best], tm)
            disp = float((pos - x0) @ e)
        v, steps = int(best), steps + 1
    return {"disp": disp if v == q else float("nan"), "steps": steps, "reached": bool(v == q)}


def _prec_to(X: np.ndarray, cand: np.ndarray, q: int) -> np.ndarray:
    """Boolean: X[cand] < X[q] in the flat relation (strict)."""
    dt = X[q, 0] - X[cand, 0]
    dx2 = ((X[cand, 1:] - X[q, 1:]) ** 2).sum(axis=1)
    return (dt > 0) & (dt * dt >= dx2)
