"""
Computational helpers for rules of motion read on a window of the order (no test here)
=====================================================================================

Shared by the F5.1f drafts and by the diagnostics diag/band.py: the kept links of the order among a set of events
(Window), the relaxation of LocalLinks for longest chains, counts of the events in the past of several seeds at once
(bit masks), the spheroid of a diamond against the faces of a box, a hash of event identities, and the longest chain
between two events used to build launch histories. None of these is a rule of the theory; they compute what a rule
asks for.
"""

from __future__ import annotations

import numpy as np

from phases.F5_gravity.F5_01c_free_fall_longest import LocalLinks

_SWEEP = LocalLinks._sweep                               # the relaxation of LocalLinks (static)
_ONE = np.uint64(1)


# ============================================================================================== the order

def sprinkle_flat_box(T, hx, hyz, rng):
    """Poisson events of density 1 in t [0, T], |x| <= hx, |y|, |z| <= hyz, sorted by time."""
    n = rng.poisson(T * (2 * hx) * (2 * hyz) ** 2)
    X = np.column_stack([rng.uniform(0, T, n), rng.uniform(-hx, hx, n),
                         rng.uniform(-hyz, hyz, n), rng.uniform(-hyz, hyz, n)])
    return X[np.argsort(X[:, 0], kind="stable")]


class Window:
    """The kept links of the order among the events sel (global indices, ascending = time order), in local indices:
    f_u -> f_w sorted by the later end (forward sweeps), b_u -> b_w sorted by the earlier end, and the children of
    every event (CSR)."""

    def __init__(self, order, sel):
        self.ids, self.n = sel, sel.size
        lo = order.pptr[sel]
        cnt = order.pptr[sel + 1] - lo
        tot = int(cnt.sum())
        child = np.repeat(np.arange(self.n), cnt)
        offs = np.repeat(lo - np.concatenate([[0], np.cumsum(cnt)[:-1]]), cnt) + np.arange(tot)
        par_g = order.pidx[offs]                          # the parents of every event of sel
        par = np.minimum(np.searchsorted(sel, par_g), max(self.n - 1, 0))
        ok = sel[par] == par_g                            # links with both ends in sel
        self.f_u, self.f_w = par[ok], child[ok]           # sorted by the child
        o = np.argsort(self.f_u, kind="stable")
        self.b_u, self.b_w = self.f_u[o], self.f_w[o]
        self.c_ptr = np.concatenate([[0], np.cumsum(np.bincount(self.b_u, minlength=self.n))]).astype(np.int64)
        self.c_idx = self.b_w

    def loc(self, g):
        """Local index of the global event g, -1 if it is not in the window."""
        j = int(np.searchsorted(self.ids, g))
        return j if j < self.n and self.ids[j] == g else -1


def back_counts(win, seeds, alive):
    """For every seed (local index) the number of alive events in its past (the seed included), by one bit per seed
    propagated backwards along the links among alive events."""
    W = (len(seeds) + 63) // 64
    M = np.zeros((win.n, W), dtype=np.uint64)
    for k, s in enumerate(seeds):
        M[s, k >> 6] |= _ONE << np.uint64(k & 63)
    ok = alive[win.b_u] & alive[win.b_w]
    a, b = win.b_w[ok], win.b_u[ok]                       # child -> parent, sorted by the parent
    if b.size:
        starts = np.flatnonzero(np.r_[True, b[1:] != b[:-1]])
        tg = b[starts]
        while True:
            new = M[tg] | np.bitwise_or.reduceat(M[a], starts, axis=0)
            if np.array_equal(new, M[tg]):
                break
            M[tg] = new
    nodes = np.flatnonzero(alive)
    out = np.zeros(len(seeds), dtype=np.int64)
    sh = np.arange(64, dtype=np.uint64)
    for w in range(W):
        c = ((M[nodes, w][:, None] >> sh[None, :]) & _ONE).sum(axis=0)
        k = min(64, len(seeds) - 64 * w)
        out[64 * w:64 * w + k] = c[:k]
    return out


def diamond_margin(xp, tp, Xq, half):
    """Smallest distance of the spatial projection of the diamonds I(p, q) (spheroids with foci x_p, x_q and the sum
    of distances t_q - t_p) from the faces of the box |x_i| <= half_i; negative: the diamond leaves the box."""
    xq, tq = Xq[:, 1:], Xq[:, 0]
    a = 0.5 * (tq - tp)
    dv = xq - xp[None, :]
    c = 0.5 * np.linalg.norm(dv, axis=1)
    f = np.where(c[:, None] > 0, dv / np.maximum(2 * c[:, None], 1e-300), 0.0)
    b2 = np.maximum(a * a - c * c, 0.0)
    ext = np.sqrt(a[:, None] ** 2 * f ** 2 + b2[:, None] * (1.0 - f ** 2))
    cen = 0.5 * (xq + xp[None, :])
    return np.min(half[None, :] - np.abs(cen) - ext, axis=1)


def id_hash(g):
    """A pseudo-random order of events by their identity (splitmix64), independent of coordinates and of time."""
    h = np.asarray(g, dtype=np.uint64) * np.uint64(0x9E3779B97F4A7C15)
    h ^= h >> np.uint64(30)
    h *= np.uint64(0xBF58476D1CE4E5B9)
    h ^= h >> np.uint64(27)
    return h


def chain_between(order, s_, q, region):
    """The longest chain s_ -> q inside a cylinder (stepping to the earliest child along the tube), as in the
    diagnostics; None if q is not reachable."""
    X = order.X
    ll = LocalLinks(order, X[s_, 1:], region, float(X[s_, 0]) - 1.0, float(X[q, 0]) + 1.0)
    js, jq = ll.pos[s_], ll.pos[q]
    if js < 0 or jq < 0:
        return None
    Lf, Lb = ll.forward(js), ll.backward(jq)
    if Lf[jq] < 2:
        return None
    path, j = [js], js
    while j != jq:
        kids = ll.c_idx[ll.c_ptr[j]:ll.c_ptr[j + 1]]
        kids = kids[(Lf[kids] == Lf[j] + 1) & (Lb[kids] == Lb[j] - 1)]
        if kids.size == 0:
            return None
        j = int(kids.min())
        path.append(j)
    return [int(ll.ids[k]) for k in path]


