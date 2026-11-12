"""
B1, development rows: inertia of a domain wall in 2+1 on the intrinsic instrument (variant B; decides nothing)
=============================================================================================================
Decision 20: the frozen tests of direction B use variant B (the longest chain with the weight of decision 17). These
rows measure, before a card is written, the spread of the quantities B1 would decide on, so that the tolerances and the
sample size follow the rule of 10 Nov (expected 2 SE + the largest development offset + a margin).

Model as F3.14a in 2+1: phi = phi_0 - G*(mu sin phi + f chi(t)), G from the chain weight, mu = 1; the wall is the
sine-Gordon kink along y; pyramid |x| <= T + 8 - t, |y| <= T + 2.5 - t; readouts on the band |y| < 2 (coordinates,
declared): per slab of width 1 the kink fit (centre, width), the charge Q from the plateaus, the energy per unit length
from local regressions (field1d.energy on the band; the wall does not depend on y). Arms on the same sprinkling:
  rest; v0.3; v0.6 (launched); kick (at rest, force f in t in [1, 3)).
Reported: velocity, contraction C = w(v) gamma / w(rest), energy per unit length against 8, Q, drift of the rest wall,
the kick Delta(gamma v) and r_M = 2 pi Q f tau / Delta(gamma v) / M_rest, each with its spread per sprinkling (sd), which
fixes the sample size of the card.

    python -m diag.b1_dev --rho 30 --seeds 8 --T 8 --workers 2
"""

from __future__ import annotations

import argparse
import json
import math
import os
import time

import numpy as np

from core.order.sprinkle import derive_seed
from core.stats.parallel import parallel_map, default_workers
from phases.F3_structure import field1d as F
from diag.b0_green2p1 import sprinkle_pyramid
from diag.b0_chain_alg import links3
from diag.b0_chain_field import weight_table, newest_calib

_nb = F._nb

if _nb is not None:
    @_nb.njit(cache=True)
    def evolve_force(ptr, idx, t, phi0, rho, mu, wL, f, t1, t2):
        N = phi0.size
        acc = np.zeros(N)
        phi = np.empty(N)
        L = np.zeros(N, np.int64)
        stamp = np.full(N, -1, np.int64)
        for v in range(N):
            ph = phi0[v] - acc[v]
            phi[v] = ph
            J = mu * math.sin(ph)
            if t[v] >= t1 and t[v] < t2:
                J += f
            coef = J / rho
            L[v] = 0
            stamp[v] = v
            for u in range(v + 1, N):
                b = -1
                for q in range(ptr[u], ptr[u + 1]):
                    p = idx[q]
                    if p >= v and stamp[p] == v and L[p] + 1 > b:
                        b = L[p] + 1
                if b > 0:
                    L[u] = b
                    stamp[u] = v
                    acc[u] += wL[b] * coef
        return phi


def energy_db(t, x, phi, tm, xc, mu=1.0, half=5.0, cell=0.5, step=0.25):
    """As field1d.energy (local linear regressions in cells), with the noise bias of the squared estimates removed:
    E[c^2] = c^2 + Var(c) and E[cos(c0 + d)] = cos(c0) exp(-Var(c0)/2), Var from the residuals of each regression."""
    s = np.abs(t - tm) < 0.5 * cell
    ts, xs, ps = t[s], x[s], phi[s]
    o = np.argsort(xs)
    ts, xs, ps = ts[o], xs[o], ps[o]
    E, used = 0.0, 0
    for xg in np.arange(xc - half, xc + half + 1e-9, step):
        i0, i1 = np.searchsorted(xs, [xg - 0.5 * cell, xg + 0.5 * cell])
        n = i1 - i0
        if n < 8:
            continue
        A = np.column_stack([np.ones(n), ts[i0:i1] - tm, xs[i0:i1] - xg])
        c, *_ = np.linalg.lstsq(A, ps[i0:i1], rcond=None)
        r = ps[i0:i1] - A @ c
        s2 = float(r @ r) / (n - 3)
        C = s2 * np.linalg.inv(A.T @ A)
        E += (0.5 * (c[1] ** 2 - C[1, 1]) + 0.5 * (c[2] ** 2 - C[2, 2])
              + mu * (1.0 - math.cos(c[0]) * math.exp(0.5 * C[0, 0]))) * step
        used += 1
    n_all = int(round(2 * half / step)) + 1
    return E * n_all / used if used else float("nan")


ARMS = {"rest": (0.0, 0.0), "v0.3": (0.3, 0.0), "v0.6": (0.6, 0.0), "kick": (0.0, 0.05)}


def readout(t, x, y, phi, v, T):
    band = np.abs(y) < 2.0
    tb, xb, pb = t[band], x[band], phi[band]
    xc, w, Q, E, Ed = [], [], [], [], []
    for a in np.arange(0.0, T - 1.0 + 1e-9, 1.0):
        g = xc[-1] if (xc and np.isfinite(xc[-1])) else v * (a + 0.5)
        s = (tb >= a) & (tb < a + 1.0)
        r = F.fit_kink(tb[s], xb[s], pb[s], a, a + 1.0, g, 1.0)
        xc.append(r["xc"] if r else float("nan"))
        w.append(r["w"] if r else float("nan"))
        xg = v * (a + 0.5)
        hi = s & (xb > xg + 5) & (xb < xg + 7)
        lo = s & (xb < xg - 5) & (xb > xg - 7)
        Q.append(float((np.mean(pb[hi]) - np.mean(pb[lo])) / (2 * math.pi)) if hi.any() and lo.any() else float("nan"))
        xe = xc[-1] if np.isfinite(xc[-1]) else xg
        E.append(F.energy(tb, xb, pb, a + 0.5, xe, 1.0))
        Ed.append(energy_db(tb, xb, pb, a + 0.5, xe, 1.0))
    return {"xc": xc, "w": w, "Q": Q, "E": E, "Ed": Ed}


CACHE = os.path.join("diag", "out", "b1_dev_cache")


def task(rho, seed, arm, T, calib):
    fn = os.path.join(CACHE, f"r{rho:g}_T{T:g}_{seed}_{arm}.json")
    if os.path.exists(fn):
        with open(fn) as fh:
            return json.load(fh)
    out = _task(rho, seed, arm, T, calib)
    os.makedirs(CACHE, exist_ok=True)
    with open(fn + ".tmp", "w") as fh:
        json.dump(out, fh)
    os.replace(fn + ".tmp", fn)
    return out


def _task(rho, seed, arm, T, calib):
    v, f = ARMS[arm]
    t, x, y = sprinkle_pyramid(rho, T + 8.0, T + 2.5, T, seed)
    t0 = time.time()
    ptr, idx = links3(t, x, y, (12.0 * 25.0 / (math.pi * rho)) ** (1.0 / 3.0))
    phi = evolve_force(ptr, idx, t, F.kink_free(t, x, v), float(rho), 1.0, weight_table("B", rho, calib), f, 1.0, 3.0)
    out = {"rho": rho, "seed": seed, "arm": arm, "N": int(t.size), "sec": time.time() - t0, "v": v, "f": f}
    out.update(readout(t, x, y, phi, v, T))
    if f == 0.0:
        ex = readout(t, x, y, F.kink_exact(t, x, v), v, T)
        out.update({"E_ex": ex["E"], "w_ex": ex["w"], "Q_ex": ex["Q"]})
    return out


def slope(tm, xc, m):
    return float(np.polyfit(tm[m], np.array(xc)[m], 1)[0])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rho", type=float, default=30.0)
    ap.add_argument("--seeds", type=int, default=8)
    ap.add_argument("--T", type=float, default=8.0)
    ap.add_argument("--workers", type=int, default=0)
    a = ap.parse_args()
    if _nb is None:
        raise SystemExit("numba required")
    calib = newest_calib()[:3]
    jobs = [(a.rho, derive_seed(20261111, "B1-dev", a.rho, k), arm, a.T, calib) for k in range(a.seeds) for arm in ARMS]
    t0 = time.time()
    res = parallel_map(task, jobs, default_workers(a.workers), print,
                       lambda i, r, d, n: f"[B1 dev] {d}/{n}: {r['arm']} N {r['N']} {r['sec']:.0f}s ({time.time() - t0:.0f}s)")
    T = a.T
    tm = np.arange(0.5, T, 1.0)
    allm = tm < T - 1.0            # slab rule: the last slab is not read (at v 0.6 the plateau window x_c + 7 leaves the
    late = (tm >= 4.0) & allm      # pyramid |x| <= T + 8 - t there; development rows of 12 Nov)

    def mm(a):
        return float(np.nanmean(np.array(a, float)[allm]))
    by = {arm: {r["seed"]: r for r in res if r["arm"] == arm} for arm in ARMS}
    seeds = sorted(by["rest"])
    L_ = [f"B1 development rows: rho {a.rho:g}, T {T:g}, {len(seeds)} sprinklings, variant B; mean +- SE (sd per sprinkling)"]

    def ms(vals):
        vals = np.array([q for q in vals if np.isfinite(q)])
        return f"{np.mean(vals):+.4f} +- {np.std(vals, ddof=1) / math.sqrt(vals.size):.4f} (sd {np.std(vals, ddof=1):.4f}, n {vals.size})"

    rest = by["rest"]
    Mrest = [mm(rest[s]["Ed"]) for s in seeds]
    wrest = [mm(rest[s]["w"]) for s in seeds]
    L_.append(f"  rest: drift x_c(last read) - x_c(first) {ms([np.array(rest[s]['xc'])[allm][-1] - rest[s]['xc'][0] for s in seeds])}; "
              f"E (noise-debiased) / 8 - 1 {ms([m / 8 - 1 for m in Mrest])}; plain E - E_ex (/8) "
              f"{ms([(mm(rest[s]['E']) - mm(rest[s]['E_ex'])) / 8 for s in seeds])}; Q - Q_ex {ms([mm(rest[s]['Q']) - mm(rest[s]['Q_ex']) for s in seeds])}; "
              f"w / w_ex - 1 {ms([mm(rest[s]['w']) / mm(rest[s]['w_ex']) - 1 for s in seeds])}")
    for arm in ("v0.3", "v0.6"):
        v = ARMS[arm][0]
        g = 1 / math.sqrt(1 - v * v)
        R = by[arm]
        L_.append(f"  {arm}: v_meas - v {ms([slope(tm, R[s]['xc'], allm) - v for s in seeds])}; "
                  f"C - 1 {ms([mm(R[s]['w']) * g / wr - 1 for s, wr in zip(seeds, wrest)])}; "
                  f"w / w_ex - 1 {ms([mm(R[s]['w']) / mm(R[s]['w_ex']) - 1 for s in seeds])}; "
                  f"e_B (debiased E) {ms([(mm(R[s]['Ed']) - 8 * g) / (8 * g) for s in seeds])}; "
                  f"plain (E - E_ex)/(8 gamma) {ms([(mm(R[s]['E']) - mm(R[s]['E_ex'])) / (8 * g) for s in seeds])}; "
                  f"Q - Q_ex {ms([mm(R[s]['Q']) - mm(R[s]['Q_ex']) for s in seeds])}")
    K = by["kick"]
    f, tau = ARMS["kick"][1], 2.0
    dp = []
    for s in seeds:
        vk = slope(tm, K[s]["xc"], late)
        vr = slope(tm, rest[s]["xc"], late)
        dp.append(vk / math.sqrt(1 - vk * vk) - vr / math.sqrt(1 - vr * vr))
    dp = np.array(dp)
    Qk = np.nanmean([mm(K[s]["Q"]) for s in seeds])
    Min = 2 * math.pi * Qk * f * tau / np.mean(dp)
    boot = []
    rng = np.random.default_rng(1)
    for _ in range(2000):
        i = rng.integers(0, len(seeds), len(seeds))
        boot.append(2 * math.pi * Qk * f * tau / np.mean(dp[i]) / np.mean(np.array(Mrest)[i]))
    L_.append(f"  kick (f {f}, tau {tau:g}): Delta(gamma v) {ms(dp)} (continuum {2 * math.pi * f * tau / 8:.4f}); "
              f"r_M = M_in / M_rest {Min / np.mean(Mrest):.4f} +- {np.std(boot, ddof=1):.4f} "
              f"(per-sprinkling sd of r_M {np.std(2 * math.pi * Qk * f * tau / dp / np.array(Mrest), ddof=1):.4f})")
    L_.append(f"  cost: N ~{int(np.mean([r['N'] for r in res]))}, {np.mean([r['sec'] for r in res]):.0f} s per run")
    print("\n".join(L_))
    os.makedirs(os.path.join("diag", "out"), exist_ok=True)
    fn = os.path.join("diag", "out", f"b1_dev_{time.strftime('%Y%m%dT%H%M%S', time.gmtime())}.json")
    with open(fn, "w") as fh:
        json.dump({"args": vars(a), "summary": L_, "results": res}, fh)
    print(f"raw records: {fn} ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
