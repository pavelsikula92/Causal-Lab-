"""
F6.8f-dev, step 1b: the cut loci in theta and their contribution at finite q (development, decides nothing)
=========================================================================================================

docs/derivations/f68f_dev.md, step 1b. At a cut locus theta_c the first arrival changes from ray a to ray b; their
actions cross linearly, S_b - S_a = alpha (theta - theta_c). At finite q both needles lie in the interval while
|S_b - S_a| < c = q^2/2, the second one with the weight (1 - |S_b - S_a| / c)^2 (its height is reduced by the action
difference; the needle volume goes as the square of the height). Integrated over theta, the extra volume per cut
locus in the theta-mean is

    dA(q) = (1 / 2 pi) (A_a + A_b) c / (3 alpha),     c = q^2 / 2,

with A_a, A_b the plane-wave volumes of the two rays at theta_c. A grid in theta (spacing 2 pi / 128 = 0.049) does not
resolve windows of width c / alpha; a grid point exactly at a tie (p = 0: theta = 0, pi) counts both needles with
full weight instead. This file locates the cut loci, computes alpha = |d(S_b - S_a)/dtheta| (envelope theorem:
dS/dtheta = -eps int sin(p s + theta + psi) ds along a fixed-end stationary ray), A_a and A_b, and reports dA(q).

Run from the repository root (one job = one point (eps, p), cached in diag/out/f68f_1b_cut_cache/):
    python -m diag.f68f_1b_cut --points 9.87:0,12:0,40:0,12:1,16:1,40:1,16:2,40:2,40:4 --workers 8
Output diag/out/f68f_1b_cut_<UTC time>.json and .txt. The partial run of 16 Nov before this form (p = 0, eps 9.95 -
16, stopped by a restart of the container) is diag/out/f68f_1b_cut_partial.txt.
"""
from __future__ import annotations

import argparse
import math
import os
import time

import numpy as np

from core.stats.parallel import default_workers, parallel_map
from diag import f68f_1b as B
from diag import f68f_common as C

N = 4000
QS = [0.5, 0.25, 0.125]


def rays(eps, p, th):
    r = B.branches(eps, p, np.array([th]), 1000)[0]
    vv = B.newton(r, eps, p, np.full(r.size, th), N)
    vv = np.unique(np.round(vv, 10))
    _, _, S = B.end_state(vv, eps, p, np.full(vv.size, th), N)
    return vv, S


def dS_dth(v, eps, p, th):
    out = B.integrate(np.atleast_1d(v), eps, p, np.full(np.size(v), th), N)
    s = np.linspace(0, 1, N + 1)[:, None]
    ph = p * s + th + out[:, 0]
    w = np.ones(N + 1)
    w[1:-1:2] = 4
    w[2:-1:2] = 2
    return -eps * (w[:, None] * np.sin(ph)).sum(0) / (3 * N)


def follow(v, eps, p, th):
    return float(B.newton(np.array([v]), eps, p, np.array([th]), N)[0])


def S_of(v, eps, p, th):
    return float(B.end_state(np.array([v]), eps, p, np.array([th]), N)[2][0])


def cut_loci(eps, p, n_th=B.N_TH):
    ths = 2 * math.pi * np.arange(n_th) / n_th
    roots = B.branches(eps, p, ths, 1000)
    if max(len(r) for r in roots) < 2:
        return []
    data = []
    for i, th in enumerate(ths):
        vv = np.unique(np.round(B.newton(roots[i], eps, p, np.full(roots[i].size, th), N), 10))
        _, _, S = B.end_state(vv, eps, p, np.full(vv.size, th), N)
        data.append((vv, S))
    vmin = np.array([d[0][np.argmin(d[1])] for d in data])
    found = []
    for i in range(n_th):
        j = (i + 1) % n_th
        th0, th1 = ths[i], ths[i] + 2 * math.pi / n_th
        va, vb = vmin[i], vmin[j]
        # continue the left minimiser to the right point and compare
        try:
            va1 = follow(va, eps, p, th1)
        except Exception:
            va1 = np.nan
        if np.isfinite(va1) and abs(va1 - vb) < 1e-6:
            continue                                   # same branch: no switch in this interval
        # a switch, or a tie exactly at a grid point: check ties at th0 first
        vv0, S0 = data[i]
        order = np.argsort(S0)
        if len(S0) > 1 and S0[order[1]] - S0[order[0]] < 1e-10:
            a, b = vv0[order[0]], vv0[order[1]]
            thc = th0
        else:
            vb0 = follow(vb, eps, p, th0)           # right minimiser continued to the left point
            lo, hi = th0, th1
            a, b = va, vb0
            f_lo = S_of(a, eps, p, lo) - S_of(b, eps, p, lo)
            for _ in range(40):
                m = 0.5 * (lo + hi)
                a_m, b_m = follow(a, eps, p, m), follow(b, eps, p, m)
                f_m = S_of(a_m, eps, p, m) - S_of(b_m, eps, p, m)
                if f_m * f_lo > 0:
                    lo, f_lo, a, b = m, f_m, a_m, b_m
                else:
                    hi = m
                    a, b = a_m, b_m
            thc = 0.5 * (lo + hi)
            a, b = follow(a, eps, p, thc), follow(b, eps, p, thc)
        if any(abs(((thc - f[0] + math.pi) % (2 * math.pi)) - math.pi) < 1e-9 for f in found):
            continue
        da, db = dS_dth(np.array([a, b]), eps, p, thc)
        alpha = abs(db - da)
        A, m1, m2, W = B.volume_A(np.array([a, b]), eps, p, np.array([thc, thc]), N)
        gap = abs(S_of(a, eps, p, thc) - S_of(b, eps, p, thc))
        found.append((thc, alpha, float(A[0]), float(A[1]), gap, bool(np.all(m1 > 0) and np.all(m2 > 0))))
    return found


NAME = "f68f_1b_cut"
CODE = C.code_md5([os.path.abspath(__file__), B.__file__])


def job(eps, p, n_th):
    return C.cached(NAME, CODE, cut_loci, (eps, p, n_th))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--points", type=str, required=True, help="eps:p,eps:p,...")
    ap.add_argument("--n-theta", type=int, default=B.N_TH)
    ap.add_argument("--workers", type=int, default=0)
    a = ap.parse_args()
    t0 = time.time()
    pts = C.pairs(a.points)
    jobs = [(e, p, a.n_theta) for (e, p) in pts]
    res = parallel_map(job, jobs, default_workers(a.workers), print,
                       lambda i, r, d, n: f"[f68f_1b_cut] {d}/{n}: eps {jobs[i][0]:g} p {jobs[i][1]:g}, {len(r)} cut loci "
                                          f"({time.time() - t0:.0f}s)")
    lines = ["F6.8f-dev step 1b: cut loci in theta (switch of the first-arrival ray) and their contribution to the "
             "theta-mean at finite q, dA = (A_a + A_b) c / (6 pi alpha), c = q^2/2; development, decides nothing", ""]
    data = []
    for (e, p), f in zip(pts, res):
        if not f:
            lines.append(f"  p {p:3.1f} eps {e:6.2f}: no cut locus (one ray at every phase)")
        for (thc, alpha, Aa, Ab, gap, ok) in f:
            d = {q: (Aa + Ab) * (q * q / 2) / (6 * math.pi * alpha) for q in QS}
            lines.append(f"  p {p:3.1f} eps {e:6.2f}: theta_c {thc:.6f}  alpha {alpha:.4e}  A_a {Aa:.5f}  A_b {Ab:.5f}  "
                         f"|dS| {gap:.1e}  no-conj {ok}  dA(1/2) {d[0.5]:.2e}  dA(1/4) {d[0.25]:.2e}  dA(1/8) {d[0.125]:.2e}")
        data.append({"eps": e, "p": p, "cut_loci": [{"theta_c": c[0], "alpha": c[1], "A_a": c[2], "A_b": c[3],
                                                     "dS": c[4], "no_conj": c[5]} for c in f]})
    lines.append(f"cost {time.time() - t0:.0f} s")
    stem = C.save(NAME, vars(a), CODE, data, lines, t0)
    print("\n".join(lines))
    print("saved", stem + ".json")


if __name__ == "__main__":
    main()
