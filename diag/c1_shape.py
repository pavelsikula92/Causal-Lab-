"""
C1, development: where the nonlocality lives - a plateau with walls (14 Nov, after F5.3; decides nothing)
========================================================================================================
F5.3 (PASSED) showed that on a curved order the C0 rule weighs a near-null pair by the curved volume along it,
rho Omega_u^4 <Omega^4>_I^(-1/2). The local conformally coupled field needs Omega_u Omega_v instead. The two are equal
wherever Omega is constant along the pair. The author's picture of 14 Nov, 19:27: local inside a uniformly
"saturated" region, nonlocal where a signal crosses from one region to another.

Predictions (derived, before any run):
  N1  a plateau Omega = 1 - a S(x), S a smooth top hat of half-width W and wall width w: deep inside the plateau, and far
      outside it, P_exact = P_local as long as the past cone of the event has not reached a wall;
  N2  the difference P_exact - P_local starts when the past light cone of the event first reaches the wall (onset at
      dt ~ W - |x| - w for an event inside), and grows with the length of the pairs that cross the wall.
Both are computed exactly as in F5.3 (diag/c1_apex: Mecke, the diamond average), with the profile replaced, and
measured on paired sprinklings (flat parent, curved child thinned with Omega^4) with the C0 rule unchanged.

    python -m diag.c1_shape --predict                      (the maps, before the run)
    python -m diag.c1_shape --rho 60 --seeds 0-7 --workers 8
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
from diag.c1_apex import diamond_nodes, solve, _flat_selfconsistent, grid_at
from phases.F3_structure import field1d as F

_nb = F._nb
A_P, W_P, WW_P = 0.15, 2.5, 0.4          # depth, half-width, wall width of the plateau
MU, AMP = 1.0, 0.05
TB_W, XB_W = 0.5, 0.5


def omega(x, a=A_P, W=W_P, ww=WW_P):
    return 1.0 - a * 0.5 * (np.tanh((x + W) / ww) - np.tanh((x - W) / ww))


if _nb is not None:
    @_nb.njit(cache=True)
    def _om(x, a, W, ww):
        return 1.0 - a * 0.5 * (math.tanh((x + W) / ww) - math.tanh((x - W) / ww))

    @_nb.njit(cache=True)
    def kernel_table_p(xs, dts, J, DX, a, W, ww, rho, tc2, w, c, mu, nU, nZ, nW, gs, gws, gx, gwx, mode):
        """As diag/c1_apex.kernel_table for the plateau profile (mode 0 exact mean, mode 1 local conformal)."""
        nx = xs.size
        nt = dts.size
        out = np.zeros((nx, nt, 2 * J + 1))
        nn = w.size
        for ix in range(nx):
            xv = xs[ix]
            ov = _om(xv, a, W, ww)
            for it in range(1, nt):
                dt = dts[it]
                for q in range(gx.size):
                    ddx = dt * gx[q]
                    wq = dt * gwx[q]
                    smax = dt * dt - ddx * ddx
                    if smax > tc2:
                        smax = tc2
                    if smax <= 0.0:
                        continue
                    xu = xv + ddx
                    ou = _om(xu, a, W, ww)
                    lam_u = rho * ou ** 4
                    if mode == 1:
                        lam_u = lam_u / (ou * ov)
                    acc = 0.0
                    for p in range(gs.size):
                        sig = 0.5 * smax * (gs[p] + 1.0)
                        wsig = 0.5 * smax * gws[p]
                        e = 1.0
                        if mode == 0:
                            e = 0.0
                            xm = xv + 0.5 * ddx
                            sp = math.sqrt(sig + ddx * ddx)
                            for k in range(nW.size):
                                e += nW[k] * _om(xm + ddx * nU[k] + sp * nZ[k], a, W, ww) ** 4
                        lam = rho * math.pi / 24.0 * sig * sig * e
                        pk = math.exp(-lam)
                        kk = 0.0
                        for n in range(nn):
                            if n > 0:
                                pk *= lam / n
                            kk += w[n] * pk
                        acc += wsig * kk
                    val = c * mu * math.pi * lam_u * acc * wq
                    fj = ddx / DX
                    j0 = int(math.floor(fj))
                    fr = fj - j0
                    if -J <= j0 <= J:
                        out[ix, it, J + j0] += val * (1.0 - fr)
                    if -J <= j0 + 1 <= J:
                        out[ix, it, J + j0 + 1] += val * fr
        return out


def mean_field(rho, mode=0, a=A_P, T=6.0, Xg=9.0, DX=0.1, DT=0.05, ns=10, nfine=48, lam_cut=22.0):
    w = minvar(8)
    ct = delta(w, rho)
    xs = np.arange(-Xg, Xg + 1e-9, DX)
    ts = np.arange(0.0, T + 1e-9, DT)
    J = int(math.ceil(T / DX)) + 1
    tc2 = math.sqrt(24.0 * lam_cut / (math.pi * rho))
    c = 1.0 / (2.0 * math.pi * math.sqrt(6.0 * rho))
    U, Z, Wd = diamond_nodes()
    gs, gws = np.polynomial.legendre.leggauss(ns)
    gx, gwx = np.polynomial.legendre.leggauss(nfine)
    W = kernel_table_p(xs, ts, J, DX, a, W_P, WW_P, rho, tc2, w, c, MU, U, Z, Wd, gs, gws, gx, gwx, mode)
    Wf = kernel_table_p(np.array([0.0]), ts, J, DX, 0.0, W_P, WW_P, rho, tc2, w, c, MU, U, Z, Wd, gs, gws, gx, gwx, 0)
    far = _flat_selfconsistent(Wf, ts.size, J, DT, 1.0, ct, MU)
    f = solve(W, xs.size, ts.size, J, DT, 1.0, ct, MU, far)
    return ts, xs, f, far


def predict(rho=60.0):
    ts, xs, f0, far = mean_field(rho, 0)
    _, _, f1, _ = mean_field(rho, 1)
    print(f"plateau a {A_P}, half-width {W_P}, wall {WW_P}; rho {rho:g}; P_exact - P_local (in A) by t (rows) and x:")
    xq = [0.0, 1.0, 1.75, 2.5, 3.25, 4.5]
    print("   t    " + " ".join(f"x {x:+.2f}" for x in xq))
    for tq in (1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0):
        i = int(round(tq / (ts[1] - ts[0])))
        row = [f0[i, int(round((x - xs[0]) / (xs[1] - xs[0])))] - f1[i, int(round((x - xs[0]) / (xs[1] - xs[0])))]
               for x in xq]
        print(f"  {tq:.1f}  " + " ".join(f"{v:+8.4f}" for v in row))
    # onset at the centre: first t where |diff| > 0.005
    i0 = int(round((0.0 - xs[0]) / (xs[1] - xs[0])))
    d = np.abs(f0[:, i0] - f1[:, i0])
    on = ts[np.argmax(d > 0.005)] if np.any(d > 0.005) else float("nan")
    print(f"onset at the centre (|P_exact - P_local| > 0.005): t {on:.2f}; the past cone reaches the wall at "
          f"t ~ {W_P - WW_P:.2f} (inner edge) .. {W_P:.2f} (mid-wall)")
    return ts, xs, f0, f1


def task(rho, seed, X=7.0, T=6.0):
    t, x, y, z = sprinkle4(rho, X, T, 41000 + seed)
    keep = np.random.default_rng(51000 + seed).uniform(size=t.size) < omega(x) ** 4
    w = minvar(8)
    ct = delta(w, rho)
    c = 1.0 / (2.0 * math.pi * math.sqrt(6.0 * rho))
    tc = (24.0 * (NMAX_ALL + 6) / (math.pi * rho)) ** 0.25
    out = {"rho": rho, "seed": seed, "X": X, "T": T, "shape": [A_P, W_P, WW_P]}
    for lab, m in (("flat", np.ones(t.size, bool)), ("curved", keep)):
        tt, xx, yy, zz = t[m], x[m], y[m], z[m]
        a0 = time.time()
        phi, _, _ = evolve_fast(tt, xx, yy, zz, tc, NMAX_ALL, w[None, :], np.array([ct]), AMP * np.ones((1, tt.size)),
                                np.array([True]), c, MU, np.zeros(tt.size, bool), 0.25, 1)
        inner = (np.abs(xx) < X - tt - 0.3) & (np.abs(yy) < X - tt - 0.3) & (np.abs(zz) < X - tt - 0.3)
        ib = np.floor(tt[inner] / TB_W).astype(int)
        jb = np.floor((xx[inner] + X) / XB_W).astype(int)
        nt, nx = int(T / TB_W) + 1, int(2 * X / XB_W) + 1
        S = np.zeros((nt, nx))
        N = np.zeros((nt, nx))
        np.add.at(S, (ib, jb), phi[0, inner] / AMP)
        np.add.at(N, (ib, jb), 1)
        out[lab] = {"S": S.tolist(), "N": N.tolist(), "sec": time.time() - a0, "n": int(tt.size),
                    "ev": np.c_[tt[inner], xx[inner]].astype(np.float32).tolist()}
    return out


def cached(rho, seed, X=7.0, T=6.0):
    import hashlib
    key = hashlib.sha256(json.dumps([str(rho), seed, X, T, A_P, W_P, WW_P, "c1shape-v1"]).encode()).hexdigest()[:20]
    fn = os.path.join("diag", "out", "c1_shape_cache", key + ".json")
    if os.path.exists(fn):
        return json.load(open(fn))
    o = task(rho, seed, X, T)
    os.makedirs(os.path.dirname(fn), exist_ok=True)
    json.dump(o, open(fn + ".tmp", "w"))
    os.replace(fn + ".tmp", fn)
    return o


def compare(R):
    """D = <phi>_curved - <phi>_flat by (t, x) against P_exact and P_local; the N1/N2 readout."""
    rho = R[0]["rho"]
    X, T = R[0]["X"], R[0]["T"]
    ts, xs, f0, far = mean_field(rho, 0, T=T)
    _, _, f1, _ = mean_field(rho, 1, T=T)
    nt, nx = np.array(R[0]["flat"]["S"]).shape
    Dm, P0, P1 = [], [], []
    for r in R:
        Sf, Nf = np.array(r["flat"]["S"]), np.array(r["flat"]["N"])
        Sc, Nc = np.array(r["curved"]["S"]), np.array(r["curved"]["N"])
        evc, evf = np.array(r["curved"]["ev"]), np.array(r["flat"]["ev"])
        bc = (np.floor(evc[:, 0] / TB_W).astype(int), np.floor((evc[:, 1] + X) / XB_W).astype(int))
        bf = (np.floor(evf[:, 0] / TB_W).astype(int), np.floor((evf[:, 1] + X) / XB_W).astype(int))
        pf = np.zeros((nt, nx))
        np.add.at(pf, bf, np.interp(evf[:, 0], ts, far))
        out = []
        for f in (f0, f1):
            pc = np.zeros((nt, nx))
            np.add.at(pc, bc, grid_at(ts, xs, f, evc))
            with np.errstate(invalid="ignore", divide="ignore"):
                out.append(pc / Nc - pf / Nf)
        with np.errstate(invalid="ignore", divide="ignore"):
            Dm.append(Sc / Nc - Sf / Nf)
        P0.append(out[0])
        P1.append(out[1])
    Dm = np.array(Dm)
    n = len(R)
    mean = np.nanmean(Dm, 0)
    se = np.nanstd(Dm, 0, ddof=1) / math.sqrt(n)
    P0, P1 = np.nanmean(np.array(P0), 0), np.nanmean(np.array(P1), 0)
    xc = (np.arange(nx) + 0.5) * XB_W - X
    tcs = (np.arange(nt) + 0.5) * TB_W
    print(f"######## plateau a {A_P}, W {W_P}, wall {WW_P}; rho {rho:g}; {n} pairs; per bin: D - local (D - exact) +- SE")
    cols = [j for j in range(nx) if xc[j] >= -0.25 and xc[j] <= 4.75]
    print("   t    " + " ".join(f"   x {xc[j]:+.2f}      " for j in cols))
    for i in range(nt):
        if tcs[i] < 0.5:
            continue
        cells = []
        for j in cols:
            if np.isfinite(mean[i, j]) and np.isfinite(se[i, j]):
                cells.append(f"{mean[i, j] - P1[i, j]:+.3f}({mean[i, j] - P0[i, j]:+.3f})+-{se[i, j]:.3f}")
            else:
                cells.append(" " * 21)
        print(f"  {tcs[i]:.2f} " + " ".join(cells))
    print("   predicted P_exact - P_local per bin:")
    for i in range(nt):
        if tcs[i] < 0.5:
            continue
        print(f"  {tcs[i]:.2f} " + " ".join(f"{P0[i, j] - P1[i, j]:+8.4f}" if np.isfinite(P0[i, j]) else "        "
                                          for j in cols))


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
    lo, hi = (int(v) for v in a.seeds.split("-")) if "-" in a.seeds else (int(a.seeds), int(a.seeds))
    jobs = [(a.rho, s) for s in range(lo, hi + 1)]
    t0 = time.time()
    res = parallel_map(cached, jobs, default_workers(a.workers), print,
                       lambda i, r, d, n: f"[C1 plateau] {d}/{n}: rho {r['rho']:g} seed {r['seed']} "
                                          f"{r['flat']['sec'] + r['curved']['sec']:.0f}s ({time.time() - t0:.0f}s)")
    compare(res)


if __name__ == "__main__":
    main()
