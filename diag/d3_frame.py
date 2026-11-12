"""
D3, development stage 1: the rest frame of matter read from the order (decision 19; decides nothing)
====================================================================================================
Decision 19: the rest frame of the matter that loads an event comes from coarse-grained matter, not from the pointwise
T_ab (degenerate wherever grad phi is spacelike). Candidates, at a target event x:
  F1a  world tube from x: the core events p of the defect (V(phi(p)) above a threshold, a scalar) in the past of x at
       longest-chain distance L_px in a window [k, k + 2); u ~ mean of the unit directions (x - p) / tau_px.
  F1b  world tube from pairs: core events q in the past of x close to it (L_qx <= 3), and for each of them core events p
       in the past of q with L_pq in [k, k + 2); u ~ mean of (q - p) / tau_pq. Does not point from the source to x when x
       lies off the core.
  F2r  coarse-grained T_ab: the sum of T_ab over the loaded events p in the past of x with L_px <= k (a Lorentz
       invariant region), with the gradient of the E1 fit of diag/d2_gradient.py (h = 3 rho^(-1/3)); u = its timelike
       eigenvector. (With the exact gradient this is exact on rigidly moving fields and tests nothing.)
The chain window is set by a proper time tau: k = 1.3 rho^(1/3) tau (the direction of a world tube is defined over a
proper time longer than its width).
The selection of events uses only the order and the scalar V; the averaging of directions and the comparison with the
true velocity use coordinates (development readout). Exact fields on the events of a pyramid sprinkling:
  wall v      the sine-Gordon kink moving at v along x (the boost along y is free: physical degeneracy);
  lump (vx, vy)  a moving Gaussian, round in its rest frame (localised, not a solution; only its world tube matters);
  pair        the target lump moving at (0, 0.3) and a second lump at rest ~4.5 away (does the frame stay local?).
Targets: on the core, and off the core at distance d from the centre line (wall) or the centre (lump).
Error: the relative rapidity eta between the estimated and the true four-velocity (cosh eta = u_est . u_true).

    python -m diag.d3_frame
"""

from __future__ import annotations

import math
import time

import numpy as np

from phases.F3_structure import field1d as F
from diag.b0_green2p1 import sprinkle_pyramid
from diag.b0_chain_alg import links3
from diag.b0_weight_valid import children

_nb = F._nb

if _nb is not None:
    @_nb.njit(cache=True)
    def back_chain(xi, cptr, cidx):
        """L(u) = longest chain u -> xi for all u <= xi (-1 if u is not in the past of xi)."""
        L = np.full(xi + 1, -1, np.int64)
        L[xi] = 0
        for u in range(xi - 1, -1, -1):
            b = -1
            for q in range(cptr[u], cptr[u + 1]):
                v = cidx[q]
                if v <= xi and L[v] >= 0 and L[v] + 1 > b:
                    b = L[v] + 1
            L[u] = b
        return L


def minkowski_unit(d):
    n2 = d[..., 0] ** 2 - d[..., 1] ** 2 - d[..., 2] ** 2
    return d / np.sqrt(np.maximum(n2, 1e-300))[..., None]


def rapidity(u, w):
    c = u[0] * w[0] - u[1] * w[1] - u[2] * w[2]
    return float(np.arccosh(max(c, 1.0)))


def config(name, t, x, y):
    if name.startswith("wall"):
        v = float(name.split()[1])
        g = 1 / math.sqrt(1 - v * v)
        z = g * (x - v * t)
        phi = 4 * np.arctan(np.exp(z))
        d = 2.0 / np.cosh(z)
        grad = np.column_stack([-g * v * d, g * d, 0 * d])
        V = 1 - np.cos(phi)
        core = V > 1.0                                         # |phi - pi| < pi / 2
        u = np.array([g, g * v, 0.0])
        dist = lambda tt, xx, yy: np.abs(g * (xx - v * tt))     # rest-frame distance from the centre line
        return phi, grad, V, core, u, dist
    if name.startswith("pair"):
        lumps = [(-2.5, 0.0, 0.0, 0.3), (2.0, 0.0, 0.0, 0.0)]      # target lump first; the other at rest, ~4.5 away
    else:
        vx, vy = (float(q) for q in name.split()[1:])
        lumps = [(0.0, 0.0, vx, vy)]
    phi = np.zeros_like(t)
    grad = np.zeros((t.size, 3))
    for cx, cy, vx, vy in lumps:
        p_, g_ = _lump(t, x - cx, y - cy, vx, vy)
        phi += p_
        grad += g_
    cx, cy, vx, vy = lumps[0]
    g = 1 / math.sqrt(1 - vx * vx - vy * vy)
    V = 0.5 * phi * phi
    core = phi > 0.5
    u = np.array([g, g * vx, g * vy])
    dist = lambda tt, xx, yy: np.hypot(*_rest(tt, xx - cx, yy - cy, vx, vy, g))
    return phi, grad, V, core, u, dist


def _lump(t, x, y, vx, vy):
    """A Gaussian exp(-|r_rest|^2), round in its rest frame, moving at (vx, vy); returns phi and its gradient."""
    g = 1 / math.sqrt(1 - vx * vx - vy * vy)
    vv = vx * vx + vy * vy
    Xr, Yr = _rest(t, x, y, vx, vy, g)
    phi = np.exp(-(Xr * Xr + Yr * Yr))
    dpx = -2 * Xr * phi
    dpy = -2 * Yr * phi
    if vv > 0:
        dXr_dx = 1 + (g - 1) * vx * vx / vv
        dXr_dy = (g - 1) * vx * vy / vv
        dYr_dx = (g - 1) * vy * vx / vv
        dYr_dy = 1 + (g - 1) * vy * vy / vv
        dXr_dt = -vx * dXr_dx - vy * dXr_dy
        dYr_dt = -vx * dYr_dx - vy * dYr_dy
    else:
        dXr_dx, dXr_dy, dYr_dx, dYr_dy, dXr_dt, dYr_dt = 1, 0, 0, 1, 0, 0
    grad = np.column_stack([dpx * dXr_dt + dpy * dYr_dt, dpx * dXr_dx + dpy * dYr_dx, dpx * dXr_dy + dpy * dYr_dy])
    return phi, grad


def _rest(tt, xx, yy, vx, vy, g):
    X, Y = xx - vx * tt, yy - vy * tt
    vv = vx * vx + vy * vy
    if vv == 0:
        return X, Y
    par = (X * vx + Y * vy) / vv
    return X + (g - 1) * par * vx, Y + (g - 1) * par * vy


def tensor_up(grad, V):
    """T^{ab} (upper indices, metric diag(1, -1, -1)) per event, shape (n, 3, 3)."""
    gu = grad * np.array([1.0, -1.0, -1.0])
    d2 = grad[:, 0] ** 2 - grad[:, 1] ** 2 - grad[:, 2] ** 2
    Lg = 0.5 * d2 - V
    eta = np.diag([1.0, -1.0, -1.0])
    return gu[:, :, None] * gu[:, None, :] - Lg[:, None, None] * eta[None]


def frame_from_T(Tsum):
    """Timelike eigenvector of T^a_b = T^{ac} eta_cb (future pointing, normalised)."""
    M = Tsum @ np.diag([1.0, -1.0, -1.0])
    w, vec = np.linalg.eig(M)
    best = None
    for k in range(3):
        e = np.real(vec[:, k])
        n2 = e[0] ** 2 - e[1] ** 2 - e[2] ** 2
        if abs(np.imag(w[k])) < 1e-9 and n2 > 1e-12:
            e = e / math.sqrt(n2)
            if e[0] < 0:
                e = -e
            if best is None or np.real(w[k]) > best[0]:
                best = (np.real(w[k]), e)
    return None if best is None else best[1]


def run(rho, name, taus, dists, seed=5, ntarg=25):
    from diag.d2_gradient import fit
    T = 5.0
    t, x, y = sprinkle_pyramid(rho, T + 4, T + 4, T, seed)
    P = np.column_stack([t, x, y])
    ptr, idx = links3(t, x, y, (12.0 * 25.0 / (math.pi * rho)) ** (1.0 / 3.0))
    cptr, cidx = children(ptr, idx, t.size)
    phi, grad, V, core, u_true, dist = config(name, t, x, y)
    loaded = V > 0.05 * V.max()
    h = 3.0 * rho ** (-1.0 / 3.0)                      # E1 neighbourhood: ~28 past events
    cache = {}

    def grad_e1(v):
        if v not in cache:
            lo = np.searchsorted(t, t[v] - h)
            u = np.arange(lo, v)
            d = P[v] - P[u]
            u = u[d[:, 0] ** 2 > d[:, 1] ** 2 + d[:, 2] ** 2]
            cache[v] = fit(P[u] - P[v], phi[u] - phi[v])
        return cache[v]

    ks = {tau: int(round(1.3 * rho ** (1.0 / 3.0) * tau)) for tau in taus}
    rng = np.random.default_rng(seed + 7)
    out = {}
    D = dist(t, x, y)
    for d in dists:
        sel = np.nonzero((t > T - 1.2) & (t < T - 0.8) & (np.abs(D - d) < 0.15) & (np.abs(x) < 3) & (np.abs(y) < 3))[0]
        if sel.size == 0:
            continue
        tg = rng.choice(sel, min(ntarg, sel.size), replace=False)
        res = {f"{m} tau{tau:g}": [] for m in ("F1a", "F1b", "F2r") for tau in taus}
        for xi in tg:
            L = back_chain(int(xi), cptr, cidx)
            past = np.nonzero(L > 0)[0]
            cp = past[core[past]]
            for tau, k in ks.items():
                m = cp[(L[cp] >= k) & (L[cp] < k + 2)]
                if m.size >= 3:
                    ue = minkowski_unit(P[xi] - P[m]).mean(axis=0)
                    res[f"F1a tau{tau:g}"].append(rapidity(minkowski_unit(ue), u_true))
                r = past[(L[past] <= k) & loaded[past]]
                if r.size:
                    r = rng.choice(r, min(150, r.size), replace=False)
                    G = [grad_e1(int(v)) for v in r]
                    ok = np.array([gg is not None for gg in G])
                    if ok.sum() >= 5:
                        Ge = np.array([gg for gg in G if gg is not None])
                        ue = frame_from_T(tensor_up(Ge, V[r[ok]]).sum(axis=0))
                        if ue is not None:
                            res[f"F2r tau{tau:g}"].append(rapidity(ue, u_true))
            qs = cp[L[cp] <= 3]
            if qs.size:
                qs = rng.choice(qs, min(6, qs.size), replace=False)
                Lq = [back_chain(int(q), cptr, cidx) for q in qs]
                for tau, k in ks.items():
                    dirs = []
                    for q, Lqq in zip(qs, Lq):
                        pp = np.nonzero((Lqq >= k) & (Lqq < k + 2))[0]
                        pp = pp[core[pp]]
                        if pp.size:
                            dirs.append(minkowski_unit(P[q] - P[pp]))
                    if dirs:
                        ue = np.concatenate(dirs).mean(axis=0)
                        res[f"F1b tau{tau:g}"].append(rapidity(minkowski_unit(ue), u_true))
        out[d] = (tg.size, res)
    return out, int(t.size), ks


def main():
    t0 = time.time()
    taus = [1.0, 2.0, 3.0]
    lines = ["D3 stage 1: rest frame of matter from the order; mean (rms) relative rapidity eta between the estimated and "
             "the true four-velocity; d = rest-frame distance of the target from the core centre; the chain window "
             "[k, k + 2) with k = 1.3 rho^(1/3) tau; F2r = summed T_ab with the E1 gradient (h = 3 rho^(-1/3)) over the "
             "loaded events of the past with L <= k. (F2 with the exact gradient is exact on rigidly moving fields: "
             "u is then an eigenvector of the pointwise T, so it is not a test.)"]
    for rho in (50.0, 150.0):
        for name, dists in (("wall 0.3", [0.0, 1.0, 2.0]), ("lump 0.3 0", [0.0, 1.0]), ("lump 0.4 0.3", [0.0, 1.0]),
                            ("lump 0 0", [0.0, 1.0]), ("pair", [0.0, 1.0])):
            out, N, ks = run(rho, name, taus, dists)
            u0 = config(name, np.zeros(1), np.zeros(1), np.zeros(1))[4]
            lines.append(f"=== rho {rho:g}, {name} (N {N}; true rapidity {math.acosh(u0[0]):.3f}; k {ks})")
            for d, (n, res) in out.items():
                parts = []
                for key, v in res.items():
                    parts.append(f"{key} {np.mean(v):.3f} ({np.sqrt(np.mean(np.square(v))):.3f}, n {len(v)})" if v
                                 else f"{key} -")
                lines.append(f"  d {d:g} ({n} targets): " + "; ".join(parts))
    print("\n".join(lines))
    print(f"({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
