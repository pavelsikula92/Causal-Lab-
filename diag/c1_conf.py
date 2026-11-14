"""
C1, development: the C0 instrument on a curved 3+1 order, first step (14 Nov; decides nothing)
==============================================================================================
The simplest curved order with an exactly known causal structure: a static, planar, conformally flat metric
    g = Omega(x)^2 (-dt^2 + dx^2 + dy^2 + dz^2),   Omega(x) = 1 - a exp(-x^2 / (2 L^2))   (a well; Omega <= 1).
Its light cones are the flat ones, so the order is the flat order of the points. Only the density changes: a sprinkling
of proper density rho has the coordinate density rho Omega^4. It is made by thinning a flat sprinkling of density rho
with probability Omega(x)^4. The flat parent and its curved child are paired: they share every surviving point, so
their difference has little noise (as kappa against kappa = 0 in F5.2d). The pyramid domain stays past-complete,
because the cones are the flat ones (the author's point of 14 Nov, 16:23: a curved domain must follow the true cones).

The instrument is that of C0 unchanged (MV8+ct, c and delta at the proper density rho; diag/c0_fast.evolve_fast). It
reads only the order and the counts in intervals. The counts now measure the curved volume.

Question (development): does the field of the rule follow the field of the curved continuum, and with which coupling
to the curvature? The continuum candidates solve phi = phi0 - G_xi * (mu phi) with G_xi the retarded Green function of
L_xi = -box_g + xi R. For data phi0 = A (constant), L_xi phi + mu phi = xi R A, in the planar case
    phi_tt = Omega^-2 (Omega^2 phi_x)_x - Omega^2 (xi R + mu) phi + Omega^2 xi R A,   R = -6 Omega^-3 Omega''.
xi = 0 (minimal) and xi = 1/6 (conformal) are compared. The readout is the paired difference
D(t, x) = <phi_curved> - <phi_flat> in bins of (t, x) over the inner events, against the same difference of the
continuum solutions (the finite-density effects of C0 largely cancel in the pair).

    python -m diag.c1_conf --rho 20 --seeds 0,1 --workers 2

First rows (14 Nov, container, rho 20, 4 paired sprinklings; decide nothing): the paired difference follows the
continuum at early times best with xi ~ 1/6 (chi2/bin 7.8 for t 1-3.5; xi 0: 15), at later times with xi ~ 0.25-0.33
(t 3.5-6). A constant xi does not fit: at rho 20 the finite-density effects of C0 (up to 0.25 by t 6) do not cancel in
the pair, because the curved region runs on a different proper time. The density dependence decides: rho 30 and 60.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import time

import numpy as np

from core.stats.parallel import parallel_map, default_workers
from diag.c0_apex import delta
from diag.c0_fast import evolve_fast
from diag.c0_links4 import sprinkle4
from diag.c0_minvar import NMAX_ALL, minvar

A_WELL, L_WELL = 0.15, 1.5
MU, AMP = 1.0, 0.05
TB_W, XB_W = 0.5, 0.5


def omega(x, a=A_WELL, L=L_WELL):
    return 1.0 - a * np.exp(-x * x / (2 * L * L))


def omega_xx(x, a=A_WELL, L=L_WELL):
    return -a * np.exp(-x * x / (2 * L * L)) * (x * x / L ** 4 - 1.0 / L ** 2)


def ricci(x):
    return -6.0 * omega_xx(x) / omega(x) ** 3


def task(rho, seed, X=7.0, T=6.0):
    t, x, y, z = sprinkle4(rho, X, T, 21000 + seed)
    keep = np.random.default_rng(31000 + seed).uniform(size=t.size) < omega(x) ** 4
    w = minvar(8)
    ct = delta(w, rho)
    c = 1.0 / (2.0 * math.pi * math.sqrt(6.0 * rho))
    tc = (24.0 * (NMAX_ALL + 6) / (math.pi * rho)) ** 0.25
    out = {"rho": rho, "seed": seed, "X": X, "T": T, "a": A_WELL, "L": L_WELL}
    for lab, m in (("flat", np.ones(t.size, bool)), ("curved", keep)):
        tt, xx, yy, zz = t[m], x[m], y[m], z[m]
        phi0 = AMP * np.ones((1, tt.size))
        a0 = time.time()
        phi, _, _ = evolve_fast(tt, xx, yy, zz, tc, NMAX_ALL, w[None, :], np.array([ct]), phi0, np.array([True]), c, MU,
                                np.zeros(tt.size, bool), 0.25, 1)
        inner = (np.abs(xx) < X - tt - 0.3) & (np.abs(yy) < X - tt - 0.3) & (np.abs(zz) < X - tt - 0.3)
        ib = np.floor(tt[inner] / TB_W).astype(int)
        jb = np.floor((xx[inner] + X) / XB_W).astype(int)
        nt, nx = int(T / TB_W) + 1, int(2 * X / XB_W) + 1
        S = np.zeros((nt, nx))
        N = np.zeros((nt, nx))
        np.add.at(S, (ib, jb), phi[0, inner] / AMP)
        np.add.at(N, (ib, jb), 1)
        out[lab] = {"S": S.tolist(), "N": N.tolist(), "sec": time.time() - a0, "n": int(tt.size),
                    "tx": np.c_[tt[inner], xx[inner]].tolist() if lab == "curved" and False else None}
        # event lists for the continuum averaging (t, x of the inner events, both sets)
        out[lab]["ev"] = np.c_[tt[inner], xx[inner]].astype(np.float32).tolist()
    return out


def cached(rho, seed, X=7.0, T=6.0):
    import hashlib
    key = hashlib.sha256(json.dumps([str(rho), seed, X, T, A_WELL, L_WELL, "c1-v1"]).encode()).hexdigest()[:20]
    fn = os.path.join("diag", "out", "c1_conf_cache", key + ".json")
    if os.path.exists(fn):
        return json.load(open(fn))
    o = task(rho, seed, X, T)
    os.makedirs(os.path.dirname(fn), exist_ok=True)
    json.dump(o, open(fn + ".tmp", "w"))
    os.replace(fn + ".tmp", fn)
    return o


def continuum(xi, T, X, a=A_WELL, L=L_WELL, nx=1601, dt=0.002, flat=False):
    """phi(t, x)/A of L_xi phi + mu phi = xi R A, phi(0) = A, phi_t(0) = 0, on x in [-X - T, X + T] (far edges flat)."""
    xs = np.linspace(-X - T - 2, X + T + 2, nx)
    dx = xs[1] - xs[0]
    if flat:
        Om = np.ones_like(xs)
        R = np.zeros_like(xs)
    else:
        Om = omega(xs, a, L)
        R = ricci(xs)
    Om2 = Om * Om
    Omh2 = 0.5 * (Om2[1:] + Om2[:-1])

    def acc(p):
        flux = Omh2 * (p[1:] - p[:-1]) / dx
        div = np.zeros_like(p)
        div[1:-1] = (flux[1:] - flux[:-1]) / dx
        return div / Om2 - Om2 * (xi * R + MU) * p + Om2 * xi * R

    p0 = np.ones_like(xs)
    p1 = p0 + 0.5 * dt * dt * acc(p0)
    nt = int(round(T / dt)) + 1
    ts = [0.0]
    fs = [p0.copy()]
    for n in range(1, nt):
        if n % 5 == 0:
            ts.append(n * dt)
            fs.append(p1.copy())
        p2 = 2 * p1 - p0 + dt * dt * acc(p1)
        p0, p1 = p1, p2
    return np.array(ts), xs, np.array(fs)


def at_events(sol, ev):
    ts, xs, f = sol
    ev = np.asarray(ev, float)
    tt, xx = ev[:, 0], ev[:, 1]
    ft = tt / (ts[1] - ts[0])
    j0 = np.clip(np.floor(ft).astype(int), 0, ts.size - 2)
    at = ft - j0
    fx = (xx - xs[0]) / (xs[1] - xs[0])
    i0 = np.clip(np.floor(fx).astype(int), 0, xs.size - 2)
    ax = fx - i0
    v0 = f[j0, i0] * (1 - ax) + f[j0, i0 + 1] * ax
    v1 = f[j0 + 1, i0] * (1 - ax) + f[j0 + 1, i0 + 1] * ax
    return v0 * (1 - at) + v1 * at


XIS = (0.0, 1.0 / 12.0, 1.0 / 6.0, 0.25, 1.0 / 3.0, 0.5)


def summary(R, xis=XIS):
    """The paired difference D against the continuum differences for a scan of xi (descriptive), by time window."""
    X, T = R[0]["X"], R[0]["T"]
    nt, nx = np.array(R[0]["flat"]["S"]).shape
    sols = {"flat": continuum(0.0, T, X, flat=True)}
    for xi in xis:
        sols[xi] = continuum(xi, T, X)
    Dm = []
    Dc = {xi: [] for xi in xis}
    for r in R:
        Sf, Nf = np.array(r["flat"]["S"]), np.array(r["flat"]["N"])
        Sc, Nc = np.array(r["curved"]["S"]), np.array(r["curved"]["N"])
        with np.errstate(invalid="ignore", divide="ignore"):
            Dm.append(Sc / Nc - Sf / Nf)
        evf, evc = np.array(r["flat"]["ev"]), np.array(r["curved"]["ev"])
        bins = {}
        for lab, ev in (("c", evc), ("f", evf)):
            bins[lab] = (np.floor(ev[:, 0] / TB_W).astype(int), np.floor((ev[:, 1] + X) / XB_W).astype(int))
        Nn = {lab: np.zeros((nt, nx)) for lab in bins}
        for lab in bins:
            np.add.at(Nn[lab], bins[lab], 1)
        vf = np.zeros((nt, nx))
        np.add.at(vf, bins["f"], at_events(sols["flat"], evf))
        for xi in xis:
            vc = np.zeros((nt, nx))
            np.add.at(vc, bins["c"], at_events(sols[xi], evc))
            with np.errstate(invalid="ignore", divide="ignore"):
                Dc[xi].append(vc / Nn["c"] - vf / Nn["f"])
    Dm = np.array(Dm)
    n = len(R)
    mean = np.nanmean(Dm, 0)
    se = np.nanstd(Dm, 0, ddof=1) / math.sqrt(n)
    C = {xi: np.nanmean(np.array(Dc[xi]), 0) for xi in xis}
    xc = (np.arange(nx) + 0.5) * XB_W - X
    tcs = (np.arange(nt) + 0.5) * TB_W
    print(f"C1 conformal well a {A_WELL}, L {L_WELL}; rho {R[0]['rho']:g}; {n} paired sprinklings; "
          f"D = <phi>_curved - <phi>_flat at |x| <= 0.5 (centre) and 2.0-2.5 (edge), in A")
    print("  t     meas(centre)   " + " ".join(f"xi {xi:.3f}" for xi in xis))
    for i in range(nt):
        cen = [j for j in range(nx) if abs(xc[j]) < 0.5 and np.isfinite(mean[i, j])]
        if not cen or tcs[i] < 0.5:
            continue
        m = np.mean(mean[i, cen])
        e = math.sqrt(np.sum(se[i, cen] ** 2)) / len(cen)
        print(f"  {tcs[i]:.2f}  {m:+.3f} +- {e:.3f}  " + " ".join(f"{np.mean(C[xi][i, cen]):+8.3f}" for xi in xis))
    for lo, hi in ((1.0, 3.5), (3.5, 6.0)):
        chi = {}
        for xi in xis:
            c2, nn = 0.0, 0
            for i in range(nt):
                if not (lo <= tcs[i] < hi):
                    continue
                for jx in range(nx):
                    if abs(xc[jx]) <= 2.5 and np.isfinite(mean[i, jx]) and se[i, jx] > 0:
                        e2 = se[i, jx] ** 2 + 0.005 ** 2
                        c2 += (mean[i, jx] - C[xi][i, jx]) ** 2 / e2
                        nn += 1
            chi[xi] = c2
        print(f"  chi2/bin, t {lo}-{hi} ({nn} bins, SE floor 0.005): " + ", ".join(f"xi {xi:.3f}: {chi[xi] / nn:.1f}" for xi in xis))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rho", type=float, default=20.0)
    ap.add_argument("--seeds", type=str, default="0,1")
    ap.add_argument("--X", type=float, default=7.0)
    ap.add_argument("--T", type=float, default=6.0)
    ap.add_argument("--workers", type=int, default=0)
    a = ap.parse_args()
    jobs = [(a.rho, int(s), a.X, a.T) for s in a.seeds.split(",")]
    t0 = time.time()
    res = parallel_map(cached, jobs, default_workers(a.workers), print,
                       lambda i, r, d, n: f"[C1 conf] {d}/{n}: rho {r['rho']:g} seed {r['seed']} "
                                          f"{r['flat']['sec'] + r['curved']['sec']:.0f}s ({time.time() - t0:.0f}s)")
    summary(res)


if __name__ == "__main__":
    main()
