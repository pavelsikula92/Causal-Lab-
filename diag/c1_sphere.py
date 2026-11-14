"""
C1, development: a spherical conformally flat well, the onion in 3D (14 Nov, the author's choice of 19:58; decides nothing)
=========================================================================================================================
g = Omega(r)^2 eta with r = |x| (3D), Omega = 1 - a exp(-r^2 / (2 L^2)). The cones stay flat, so the order is exact and the
pyramid past-complete. Paired flat / curved sprinklings (thinning with Omega(r)^4); the C0 rule unchanged.

The exact mean (as diag/c1_apex, Mecke for the inhomogeneous sprinkling), now with spherical symmetry. f(t, r) obeys
    f(t, r_v) (1 - delta mu) = A - c mu Int d^4u lambda(u) K(Lambda(u, v)) f(t_u, r_u).
With u - v = (-dt, l n), n at polar angle theta to the radial direction of v, sigma = dt^2 - l^2 (flat):
    d^4u = pi l d(dt) d(sigma) d(cos theta),   r_u = sqrt(r_v^2 + l^2 + 2 r_v l cos theta).
The mean of Omega^4 over the diamond uses the full 3D point of the unit diamond: rest-frame time u (weight R^3) and a
point of the unit ball (z1 along the pair, z2, z3 across); in the lab
    w = x_mid + n (l u + dt R z1) + tau R (z2 e2 + z3 e3),   tau = sqrt(sigma), R = 1/2 - |u|.
mode 0: the exact mean; mode 1: the conformally coupled Green function Omega_u^-1 Omega_v^-1 G_eta with the same flat
smearing (as F5.3).

    python -m diag.c1_sphere --predict --rho 60
    python -m diag.c1_sphere --rho 60 --seeds 0-7 --workers 8
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
from phases.F3_structure import field1d as F

_nb = F._nb
A_S, L_S = 0.15, 1.5
MU, AMP = 1.0, 0.05
TB_W, RB_W = 0.5, 0.5


def omega(r, a=A_S, L=L_S):
    return 1.0 - a * np.exp(-r * r / (2 * L * L))


def diamond3_nodes(nu=6, nb=26):
    """(u, z1, z2, z3, weight) for the unit 4D diamond: u with weight R^3 on each half (Gauss-Legendre), the ball by a
    fixed low-discrepancy set of nb points (radius^3 uniform, directions by a Fibonacci sphere)."""
    gu, wu = np.polynomial.legendre.leggauss(nu)
    U, W = [], []
    for sgn in (-1.0, 1.0):
        u = sgn * (0.25 + 0.25 * gu)
        U.append(u)
        W.append(0.25 * wu * (0.5 - np.abs(u)) ** 3)
    U, W = np.concatenate(U), np.concatenate(W)
    k = np.arange(nb) + 0.5
    rad = (k / nb) ** (1.0 / 3.0)
    phi = math.pi * (1 + 5 ** 0.5) * k
    cz = 1 - 2 * ((k * 0.6180339887) % 1.0)
    sz = np.sqrt(1 - cz * cz)
    Z = np.c_[rad * sz * np.cos(phi), rad * sz * np.sin(phi), rad * cz]
    nodes = []
    for u, w in zip(U, W):
        R = 0.5 - abs(u)
        for z in Z:
            nodes.append((u, R * z[0], R * z[1], R * z[2], w / nb))
    A = np.array(nodes)
    A[:, 4] /= A[:, 4].sum()
    return A


if _nb is not None:
    @_nb.njit(cache=True)
    def _om(r, a, L):
        return 1.0 - a * math.exp(-r * r / (2.0 * L * L))

    @_nb.njit(cache=True)
    def kernel_table(rs, dts, a, L, rho, tc2, w, c, mu, D, gs, gws, gc, gwc, mode):
        """W[ir, it, jr]: the weight of f(t - dt_it, r_jr) in the mean sum at (t, r_ir); linear spreading in r."""
        nr = rs.size
        nt = dts.size
        DR = rs[1] - rs[0]
        out = np.zeros((nr, nt, nr))
        nn = w.size
        for ir in range(nr):
            rv = rs[ir]
            ov = _om(rv, a, L)
            for it in range(1, nt):
                dt = dts[it]
                smax = dt * dt
                if smax > tc2:
                    smax = tc2
                for p in range(gs.size):
                    sig = 0.5 * smax * (gs[p] + 1.0)
                    wsig = 0.5 * smax * gws[p]
                    l = math.sqrt(dt * dt - sig)
                    tau = math.sqrt(sig)
                    for q in range(gc.size):
                        ct = gc[q]
                        wq = gwc[q]
                        st = math.sqrt(max(0.0, 1.0 - ct * ct))
                        # v at (0, 0, rv); n = (st, 0, ct); u at v + l n
                        xu, zu = l * st, rv + l * ct
                        ru = math.sqrt(xu * xu + zu * zu)
                        ou = _om(ru, a, L)
                        lam_u = rho * ou ** 4
                        if mode == 1:
                            lam_u = lam_u / (ou * ov)
                        e = 1.0
                        if mode == 0:
                            e = 0.0
                            mx, mz = 0.5 * xu, rv + 0.5 * l * ct
                            # e2 = (ct, 0, -st), e3 = (0, 1, 0)
                            for k in range(D.shape[0]):
                                al = l * D[k, 0] + dt * D[k, 1]
                                px = mx + st * al + tau * ct * D[k, 2]
                                py = tau * D[k, 3]
                                pz = mz + ct * al - tau * st * D[k, 2]
                                e += D[k, 4] * _om(math.sqrt(px * px + py * py + pz * pz), a, L) ** 4
                        lam = rho * math.pi / 24.0 * sig * sig * e
                        pk = math.exp(-lam)
                        kk = 0.0
                        for n in range(nn):
                            if n > 0:
                                pk *= lam / n
                            kk += w[n] * pk
                        val = c * mu * math.pi * l * lam_u * kk * wsig * wq
                        fj = ru / DR
                        j0 = int(math.floor(fj))
                        fr = fj - j0
                        if j0 < nr:
                            out[ir, it, j0] += val * (1.0 - fr)
                        if j0 + 1 < nr:
                            out[ir, it, j0 + 1] += val * fr
        return out

    @_nb.njit(cache=True)
    def solve_r(W, nr, nt, dtau, A, ct, mu):
        f = np.zeros((nt, nr))
        for ir in range(nr):
            f[0, ir] = A / (1.0 - ct * mu)
        for it in range(1, nt):
            for ir in range(nr):
                s = 0.0
                for k in range(1, it + 1):
                    wk = 0.5 if k == it else 1.0
                    for jr in range(nr):
                        wv = W[ir, k, jr]
                        if wv != 0.0:
                            s += wk * wv * f[it - k, jr]
                f[it, ir] = (A - dtau * s) / (1.0 - ct * mu)
        return f


def mean_field(rho, mode=0, a=A_S, T=6.0, Rg=13.0, DR=0.1, DT=0.05, ns=24, nc=24, lam_cut=100.0, nu=10, nb=80):
    """f(t, r) on (ts, rs); Rg must exceed the largest r of interest + T (the far field is then never reached)."""
    w = minvar(8)
    ct = delta(w, rho)
    rs = np.arange(0.0, Rg + 1e-9, DR)
    ts = np.arange(0.0, T + 1e-9, DT)
    tc2 = math.sqrt(24.0 * lam_cut / (math.pi * rho))
    c = 1.0 / (2.0 * math.pi * math.sqrt(6.0 * rho))
    D = diamond3_nodes(nu, nb)
    gs, gws = np.polynomial.legendre.leggauss(ns)
    gc, gwc = np.polynomial.legendre.leggauss(nc)
    Wt = kernel_table(rs, ts, a, L_S, rho, tc2, w, c, MU, D, gs, gws, gc, gwc, mode)
    return ts, rs, solve_r(Wt, rs.size, ts.size, DT, 1.0, ct, MU)


def predict(rho):
    t0 = time.time()
    ts, rs, f0 = mean_field(rho, 0)
    _, _, f1 = mean_field(rho, 1)
    _, _, ff = mean_field(rho, 0, a=0.0)
    print(f"spherical well a {A_S}, L {L_S}; rho {rho:g}; ({time.time() - t0:.0f}s)")
    rq = [0.0, 1.0, 2.0, 3.0, 4.5]
    for lab, g in (("exact - flat (D)", f0 - ff), ("exact - local", f0 - f1)):
        print(f"  {lab} by t (rows) and r = " + " ".join(f"{r:.1f}" for r in rq))
        for tq in (1.0, 2.0, 3.0, 4.0, 5.0, 6.0):
            i = int(round(tq / (ts[1] - ts[0])))
            print(f"   t {tq:.1f}: " + " ".join(f"{g[i, int(round(r / (rs[1] - rs[0])))]:+.4f}" for r in rq))
    return ts, rs, f0, f1, ff


def task(rho, seed, X=7.0, T=6.0):
    t, x, y, z = sprinkle4(rho, X, T, 61000 + seed)
    r = np.sqrt(x * x + y * y + z * z)
    keep = np.random.default_rng(71000 + seed).uniform(size=t.size) < omega(r) ** 4
    w = minvar(8)
    ct = delta(w, rho)
    c = 1.0 / (2.0 * math.pi * math.sqrt(6.0 * rho))
    tc = (24.0 * 100.0 / (math.pi * rho)) ** 0.25          # the cut never binds (proper rho V <= 100 Omega^4 >= 52)
    out = {"rho": rho, "seed": seed, "X": X, "T": T, "shape": [A_S, L_S], "lam_cut": 100.0}
    for lab, m in (("flat", np.ones(t.size, bool)), ("curved", keep)):
        tt, xx, yy, zz = t[m], x[m], y[m], z[m]
        a0 = time.time()
        phi, _, _ = evolve_fast(tt, xx, yy, zz, tc, NMAX_ALL, w[None, :], np.array([ct]), AMP * np.ones((1, tt.size)),
                                np.array([True]), c, MU, np.zeros(tt.size, bool), 0.25, 1)
        inner = (np.abs(xx) < X - tt - 0.3) & (np.abs(yy) < X - tt - 0.3) & (np.abs(zz) < X - tt - 0.3)
        rr = np.sqrt(xx * xx + yy * yy + zz * zz)
        out[lab] = {"sec": time.time() - a0, "n": int(tt.size),
                    "ev": np.c_[tt[inner], rr[inner], phi[0, inner] / AMP].astype(np.float32).tolist()}
    return out


def cached(rho, seed, X=7.0, T=6.0):
    import hashlib
    key = hashlib.sha256(json.dumps([str(rho), seed, X, T, A_S, L_S, "c1sphere-v1"]).encode()).hexdigest()[:20]
    fn = os.path.join("diag", "out", "c1_sphere_cache", key + ".json")
    if os.path.exists(fn):
        return json.load(open(fn))
    o = task(rho, seed, X, T)
    os.makedirs(os.path.dirname(fn), exist_ok=True)
    json.dump(o, open(fn + ".tmp", "w"))
    os.replace(fn + ".tmp", fn)
    return o


def compare(R):
    rho = R[0]["rho"]
    T = R[0]["T"]
    ts, rs, f0 = mean_field(rho, 0, T=T)
    _, _, f1 = mean_field(rho, 1, T=T)
    _, _, ff = mean_field(rho, 0, a=0.0, T=T)
    nt, nr = int(T / TB_W) + 1, 8

    def interp(f, ev):
        it = np.clip(ev[:, 0] / (ts[1] - ts[0]), 0, ts.size - 1.0001)
        ir = np.clip(ev[:, 1] / (rs[1] - rs[0]), 0, rs.size - 1.0001)
        i0, j0 = it.astype(int), ir.astype(int)
        a_, b_ = it - i0, ir - j0
        return ((f[i0, j0] * (1 - b_) + f[i0, j0 + 1] * b_) * (1 - a_) + (f[i0 + 1, j0] * (1 - b_) + f[i0 + 1, j0 + 1] * b_) * a_)

    Dm, P0, P1 = [], [], []
    for r in R:
        cur, fl = np.array(r["curved"]["ev"]), np.array(r["flat"]["ev"])
        res = {}
        for lab, ev in (("c", cur), ("f", fl)):
            ib = (ev[:, 0] / TB_W).astype(int)
            jb = (ev[:, 1] / RB_W).astype(int)
            m = jb < nr
            N = np.zeros((nt, nr))
            S = np.zeros((nt, nr))
            np.add.at(N, (ib[m], jb[m]), 1)
            np.add.at(S, (ib[m], jb[m]), ev[m, 2])
            p = {}
            for k, f in (("0", f0), ("1", f1), ("f", ff)):
                Pp = np.zeros((nt, nr))
                np.add.at(Pp, (ib[m], jb[m]), interp(f, ev[m]))
                p[k] = Pp
            res[lab] = (S, N, p)
        Sc, Nc, pc = res["c"]
        Sf, Nf, pf = res["f"]
        with np.errstate(invalid="ignore", divide="ignore"):
            Dm.append(Sc / Nc - Sf / Nf)
            P0.append(pc["0"] / Nc - pf["f"] / Nf)
            P1.append(pc["1"] / Nc - pf["f"] / Nf)
    Dm, P0, P1 = np.array(Dm), np.nanmean(np.array(P0), 0), np.nanmean(np.array(P1), 0)
    n = len(R)
    mean = np.nanmean(Dm, 0)
    se = np.nanstd(Dm, 0, ddof=1) / math.sqrt(n)
    print(f"######## spherical well a {A_S}, L {L_S}; rho {rho:g}; {n} pairs; D = <phi>_curved - <phi>_flat by (t, r)")
    for i in range(nt):
        tc_ = (i + 0.5) * TB_W
        if tc_ < 1.0:
            continue
        cells = []
        for j in range(nr):
            if np.isfinite(mean[i, j]) and np.isfinite(se[i, j]):
                cells.append(f"r {(j + 0.5) * RB_W:.2f}: {mean[i, j]:+.3f}+-{se[i, j]:.3f} (ex {P0[i, j]:+.3f}, loc {P1[i, j]:+.3f})")
        print(f"  t {tc_:.2f}: " + "; ".join(cells))
    sel = np.isfinite(mean) & np.isfinite(se) & (se > 0)
    sel[: int(1.0 / TB_W)] = False
    z0 = np.abs(mean - P0)[sel] / np.hypot(se[sel], 0.005)
    z1 = np.abs(mean - P1)[sel] / np.hypot(se[sel], 0.005)
    print(f"  {sel.sum()} bins (t >= 1): chi2/bin exact {np.mean(z0 ** 2):.2f} (max z {z0.max():.2f}); local "
          f"{np.mean(z1 ** 2):.2f} (max z {z1.max():.2f})")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--predict", action="store_true")
    ap.add_argument("--rho", type=float, default=60.0)
    ap.add_argument("--seeds", type=str, default="0-7")
    ap.add_argument("--workers", type=int, default=0)
    a = ap.parse_args()
    if a.predict:
        predict(a.rho)
        return
    lo, hi = (int(v) for v in a.seeds.split("-"))
    t0 = time.time()
    res = parallel_map(cached, [(a.rho, s) for s in range(lo, hi + 1)], default_workers(a.workers), print,
                       lambda i, r, d, n: f"[C1 sphere] {d}/{n}: rho {r['rho']:g} seed {r['seed']} "
                                          f"{r['flat']['sec'] + r['curved']['sec']:.0f}s ({time.time() - t0:.0f}s)")
    import contextlib
    import io
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        compare(res)
    print(buf.getvalue())
    fn = os.path.join("diag", "out", f"c1_sphere_rho{a.rho:g}_{lo}-{hi}.txt")
    open(fn, "w").write(buf.getvalue())
    print("written", fn)


if __name__ == "__main__":
    main()
