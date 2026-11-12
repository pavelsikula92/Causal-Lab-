"""
D2, development stage 1: the gradient of the field from the causal past of an event, and T_ab (decides nothing)
==============================================================================================================
Decision 16, D2: order of work  grad phi from the causal neighbourhood -> T_ab -> its components against the continuum ->
its trace -> only then the gradient in the load. Convention (fixed here for development, to be fixed before freezing):
signature (+, -, -), L = 1/2 (d phi)^2 - V, T_ab = d_a phi d_b phi - g_ab (1/2 (d phi)^2 - V),
trace T = -1/2 (d phi)^2 + 3 V in 2+1 (in 1+1 the same convention gives T = 2V, the load of F5.2b).

To separate the estimator from the evolution, phi is an EXACT field evaluated on the events of a sprinkling:
  wall  the sine-Gordon kink along y moving at v = 0.3 (gradient spacelike);
  wave  a plane wave of the Klein-Gordon equation, k = (1, 0.5), amplitude 0.5 (gradient changes character);
  lump  a moving Gaussian lump exp(-|r - v t|^2), v = (0.2, 0.1) (a localised, genuinely 2+1 gradient).
Estimators at a target event v, from its causal PAST only (the load must be causal):
  E1(h)  all past events u with 0 < t_v - t_u <= h: least squares of phi(u) - phi(v) on the coordinate displacement,
         linear + quadratic terms (10 parameters) -- a coordinate neighbourhood, development proxy;
  E2     the past links of v and the past links of those (two steps): the same fit. The neighbour set is intrinsic;
         the displacements are read in coordinates (development).
Reported per density: bias and rms of the gradient components relative to |grad phi| scale, of (d phi)^2, of the
components T_tt, T_tx, T_xx, T_yy and of the trace T, against the exact values; the number of neighbours used.
Coordinates enter for the exact fields, the displacements of the fit and the readout.

    python -m diag.d2_gradient
"""

from __future__ import annotations

import math
import time

import numpy as np

from phases.F3_structure import field1d as F
from diag.b0_green2p1 import sprinkle_pyramid
from diag.b0_chain_alg import links3


def fields(name, t, x, y):
    """phi and its exact gradient (phi_t, phi_x, phi_y) and V (mu = 1, V = 1 - cos phi for the wall, V = phi^2/2 else)."""
    if name == "wall":
        v = 0.3
        g = 1.0 / math.sqrt(1 - v * v)
        z = g * (x - v * t)
        phi = 4 * np.arctan(np.exp(z))
        d = 2.0 / np.cosh(z)                         # d phi / dz
        return phi, np.column_stack([-g * v * d, g * d, 0 * d]), 1 - np.cos(phi)
    if name == "wave":
        kx, ky, A = 1.0, 0.5, 0.5
        om = math.sqrt(kx * kx + ky * ky + 1.0)
        ph = kx * x + ky * y - om * t
        phi = A * np.cos(ph)
        s = -A * np.sin(ph)
        return phi, np.column_stack([-om * s, kx * s, ky * s]), 0.5 * phi ** 2
    if name == "lump":
        vx, vy = 0.2, 0.1
        X, Y = x - vx * t, y - vy * t
        phi = np.exp(-(X * X + Y * Y))
        return phi, np.column_stack([phi * 2 * (X * vx + Y * vy), -2 * X * phi, -2 * Y * phi]), 0.5 * phi ** 2
    raise ValueError(name)


def tensor(grad, V):
    """T_ab (lower indices) with metric diag(1, -1, -1); returns T_tt, T_tx, T_xx, T_yy, trace."""
    gt, gx, gy = grad[:, 0], grad[:, 1], grad[:, 2]
    d2 = gt * gt - gx * gx - gy * gy
    Lg = 0.5 * d2 - V
    return np.column_stack([gt * gt - Lg, gt * gx, gx * gx + Lg, gy * gy + Lg, -0.5 * d2 + 3 * V])


def fit(du, dphi):
    dt, dx, dy = du[:, 0], du[:, 1], du[:, 2]
    M = np.column_stack([dt, dx, dy, dt * dt, dx * dx, dy * dy, dt * dx, dt * dy, dx * dy])
    if M.shape[0] < 12:
        return None
    sc = np.abs(M).max(axis=0) + 1e-300
    c, *_ = np.linalg.lstsq(M / sc, dphi, rcond=None)
    return (c / sc)[:3]


def run(rho, seed, names, hs, ntarg=300):
    T = 4.0
    t, x, y = sprinkle_pyramid(rho, T + 3, T + 3, T, seed)
    ptr, idx = links3(t, x, y, (12.0 * 25.0 / (math.pi * rho)) ** (1.0 / 3.0))
    rng = np.random.default_rng(seed + 1)
    cand = np.nonzero((t > T - 1.5) & (t < T - 0.5) & (np.abs(x) < 1.5) & (np.abs(y) < 1.5))[0]
    tg = rng.choice(cand, min(ntarg, cand.size), replace=False)
    out = {}
    for name in names:
        phi, grad, V = fields(name, t, x, y)
        P = np.column_stack([t, x, y])
        ex_g = grad[tg]
        ex_T = tensor(ex_g, V[tg])
        est = {}
        for h in hs:
            G = []
            nn = []
            for v in tg:
                lo = np.searchsorted(t, t[v] - h)
                u = np.arange(lo, v)
                d = P[v] - P[u]
                inside = d[:, 0] ** 2 > d[:, 1] ** 2 + d[:, 2] ** 2
                u = u[inside]
                nn.append(u.size)
                G.append(fit(P[u] - P[v], phi[u] - phi[v]))
            est[f"E1 h {h:g}"] = (G, nn)
        G, nn = [], []
        for v in tg:
            p1 = idx[ptr[v]:ptr[v + 1]]
            p2 = np.concatenate([idx[ptr[p]:ptr[p + 1]] for p in p1]) if p1.size else np.empty(0, np.int64)
            u = np.unique(np.concatenate([p1, p2]))
            nn.append(u.size)
            G.append(fit(P[u] - P[v], phi[u] - phi[v]))
        est["E2 links (2 steps)"] = (G, nn)
        out[name] = (ex_g, ex_T, V[tg], est)
    return out, int(t.size)


def main():
    t0 = time.time()
    names = ["wall", "wave", "lump"]
    hs = [0.5, 0.8, 1.2]
    lines = ["D2 stage 1: gradient from the causal past (exact fields on the events); errors relative to the rms of the "
             "exact quantity over the targets; bias = mean error / rms exact"]
    for rho in (50.0, 200.0):
        res, N = run(rho, 11, names, hs)
        lines.append(f"=== rho {rho:g} (N {N})")
        for name in names:
            ex_g, ex_T, V, est = res[name]
            sg = np.sqrt(np.mean(np.sum(ex_g ** 2, axis=1)))
            lines.append(f"  {name}: |grad| rms {sg:.3f}")
            for lab, (G, nn) in est.items():
                ok = np.array([g is not None for g in G])
                if ok.sum() < 10:
                    lines.append(f"    {lab:20s} n~{int(np.median(nn)):4d} | too few neighbours for the 9-term fit")
                    continue
                Ge = np.array([g for g in G if g is not None])
                eg = Ge - ex_g[ok]
                Te = tensor(Ge, V[ok])
                eT = Te - ex_T[ok]
                sT = np.sqrt(np.mean(ex_T[ok] ** 2, axis=0)) + 1e-300
                d2e = Ge[:, 0] ** 2 - Ge[:, 1] ** 2 - Ge[:, 2] ** 2
                d2x = ex_g[ok, 0] ** 2 - ex_g[ok, 1] ** 2 - ex_g[ok, 2] ** 2
                sd = np.sqrt(np.mean(d2x ** 2)) + 1e-300
                lines.append(
                    f"    {lab:20s} n~{int(np.median(nn)):4d} | grad rms err t/x/y "
                    + "/".join(f"{np.sqrt(np.mean(eg[:, k] ** 2)) / sg:.3f}" for k in range(3))
                    + f" bias " + "/".join(f"{np.mean(eg[:, k]) / sg:+.3f}" for k in range(3))
                    + f" | (dphi)^2 rms {np.sqrt(np.mean((d2e - d2x) ** 2)) / sd:.3f} bias {np.mean(d2e - d2x) / sd:+.3f}"
                    + " | T_tt,T_tx,T_xx,T_yy,trace rms " + ",".join(f"{np.sqrt(np.mean(eT[:, k] ** 2)) / sT[k]:.2f}"
                                                                     for k in range(5))
                    + " bias " + ",".join(f"{np.mean(eT[:, k]) / sT[k]:+.2f}" for k in range(5)))
    print("\n".join(lines))
    print(f"({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
