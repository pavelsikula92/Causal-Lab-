"""
B0, development stage 1: the retarded Green function of the order in 2+1 as a field evolution, with COORDINATE tau
==================================================================================================================
(development; decides nothing; the coordinate stage approved by the author, decision 16 D1, 11 Nov 2026)

In 2+1 the retarded Green function of d_t^2 - d_x^2 - d_y^2 is G = 1 / (2 pi tau) inside the past cone (tau the proper
time). The field rule of F3.14a becomes
    phi(x) = phi_0(x) - (1 / rho) * sum over y < x of G(x, y) mu sin(phi(y)),
phi_0 the free wave of the initial data. Here tau is computed from coordinates (a declared development proxy), with a cap
tau >= tau_c = 0.68 rho^(-1/3), the proper time of one link for a longest-chain estimator with alpha_3 ~ 1.47 (the
variance of the sum near the cone diverges logarithmically without it; the intrinsic estimator has the cap built in).

Known answers used as instrument checks:
  - a y-independent source: int dy / (2 pi sqrt(a^2 - y^2)) over |y| < a = 1/2, so a domain wall (the 1+1 kink extended
    along y) must reproduce the 1+1 kink exactly (static and moving, F3.14a);
  - plane waves k = (k_x, k_y) of the linearised equation: omega = sqrt(k^2 + mu), also for k_y != 0 (genuinely 2+1);
  - vacuum stability of small random modes.
Domain: the pyramid of dependence |x| <= X - t, |y| <= Y - t of the slice t = 0 (no wall enters the dynamics).
Readout: the kink fit of F3.14a on the events of a slab with |y| < 2 (the wall's centre line), phi against x.

Run from the root of the laboratory (numba required):
    python -m diag.b0_green2p1 --quick
    python -m diag.b0_green2p1 --seeds 8 --rhos 50,100 --workers 6
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

_nb = F._nb

if _nb is not None:
    @_nb.njit(cache=True, parallel=False)
    def _evolve3(t, x, y, phi0, rho, mu, tau_c):
        N = t.size
        phi = np.empty(N)
        w = np.empty(N)
        inv2pi = 1.0 / (2.0 * math.pi)
        for i in range(N):                       # sorted by t
            s = 0.0
            ti, xi, yi = t[i], x[i], y[i]
            for j in range(i):
                dt = ti - t[j]
                dx = xi - x[j]
                dy = yi - y[j]
                q = dt * dt - dx * dx - dy * dy
                if q > 0.0:
                    tau = math.sqrt(q)
                    if tau < tau_c:
                        tau = tau_c
                    s += w[j] * inv2pi / tau
            ph = phi0[i] - s / rho
            phi[i] = ph
            w[i] = mu * math.sin(ph)
        return phi


def sprinkle_pyramid(rho, X, Y, T, seed):
    rng = np.random.default_rng(seed)
    vol = 4.0 * (X * Y * T - 0.5 * (X + Y) * T * T + T ** 3 / 3.0)
    n = int(rng.poisson(rho * vol))
    t, x, y = np.empty(0), np.empty(0), np.empty(0)
    while t.size < n:
        m = int(1.5 * (n - t.size)) + 100
        tt = rng.uniform(0, T, m)
        xx = rng.uniform(-X, X, m)
        yy = rng.uniform(-Y, Y, m)
        ok = (np.abs(xx) <= X - tt) & (np.abs(yy) <= Y - tt)
        t, x, y = np.concatenate([t, tt[ok]]), np.concatenate([x, xx[ok]]), np.concatenate([y, yy[ok]])
    t, x, y = t[:n], x[:n], y[:n]
    o = np.argsort(t, kind="stable")
    return t[o], x[o], y[o]


def modes_free_2d(t, x, y, modes, mu=1.0):
    out = np.zeros_like(x)
    for a, kx, ky, th in modes:
        k = math.hypot(kx, ky)
        om = math.sqrt(k * k + mu)
        ph = kx * x + ky * y + th
        out += a * (np.cos(ph) * np.cos(k * t) + (om / k) * np.sin(ph) * np.sin(k * t))
    return out


def task(rho, seed, arm, T, quick):
    X, Y = T + 8.0, T + 3.0
    t, x, y = sprinkle_pyramid(rho, X, Y, T, seed)
    tau_c = 0.68 * rho ** (-1.0 / 3.0)
    t0 = time.time()
    out = {"rho": rho, "seed": seed, "arm": arm, "N": int(t.size)}
    if arm in ("rest", "v0.3"):
        v = 0.0 if arm == "rest" else 0.3
        phi = _evolve3(t, x, y, F.kink_free(t, x, v), float(rho), 1.0, tau_c)
        ex = F.kink_exact(t, x, v)
        c = np.abs(y) < 2.0
        xc, rms = [], []
        for a in np.arange(0.0, T - 1.0 + 1e-9, 1.0):
            g = xc[-1] if (xc and np.isfinite(xc[-1])) else v * (a + 0.5)
            s = c & (t >= a) & (t < a + 1.0)
            r = F.fit_kink(t[s], x[s], phi[s], a, a + 1.0, g, 1.0)
            xc.append(r["xc"] if r else float("nan"))
            sn = s & (np.abs(x - v * (a + 0.5)) < 4.0)
            rms.append(float(np.sqrt(np.mean((phi[sn] - ex[sn]) ** 2))))
        out.update({"xc": xc, "rms": rms, "v": v})
    elif arm.startswith("k"):
        kx, ky = (float(u) for u in arm[1:].split("_"))
        modes = [(0.01, kx, ky, 0.0)]
        phi = _evolve3(t, x, y, modes_free_2d(t, x, y, modes), float(rho), 1.0, tau_c)
        ph, tm = [], []
        c = (np.abs(x) < 6.0) & (np.abs(y) < 4.0)
        for a in np.arange(0.0, T - 0.25 + 1e-9, 0.25):
            s = c & (t >= a) & (t < a + 0.25)
            arg = kx * x[s] + ky * y[s]
            M = np.column_stack([np.cos(arg), np.sin(arg), np.ones(s.sum())])
            cf, *_ = np.linalg.lstsq(M, phi[s], rcond=None)
            ph.append(math.atan2(cf[1], cf[0]))
            tm.append(a + 0.125)
        om = float(np.polyfit(np.array(tm), np.unwrap(np.array(ph)), 1)[0])
        out.update({"omega": om, "omega_c": math.sqrt(kx * kx + ky * ky + 1.0)})
    elif arm == "vacuum":
        rng = np.random.default_rng(derive_seed(seed, "modes"))
        modes = []
        for _ in range(8):
            k = rng.uniform(0.3, 3.0)
            th = rng.uniform(0, 2 * math.pi)
            modes.append((0.01 / 2.0, k * math.cos(th), k * math.sin(th), rng.uniform(0, 2 * math.pi)))
        phi = _evolve3(t, x, y, modes_free_2d(t, x, y, modes), float(rho), 1.0, tau_c)
        ex = np.zeros_like(x)
        for a_, kx, ky, th in modes:
            ex += a_ * np.cos(kx * x + ky * y - math.sqrt(kx * kx + ky * ky + 1.0) * t + th)
        c = (np.abs(x) < 6.0) & (np.abs(y) < 4.0)
        f0 = c & (t < 1.0)
        f1 = c & (t >= T - 1.0)
        out.update({"R": float(np.sqrt(np.mean(phi[f1] ** 2)) / np.sqrt(np.mean(phi[f0] ** 2))),
                    "R_c": float(np.sqrt(np.mean(ex[f1] ** 2)) / np.sqrt(np.mean(ex[f0] ** 2))),
                    "err": float(np.sqrt(np.mean((phi[f1] - ex[f1]) ** 2)) / np.sqrt(np.mean(ex[f1] ** 2)))})
    out["sec"] = time.time() - t0
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=8)
    ap.add_argument("--rhos", type=str, default="50,100")
    ap.add_argument("--T", type=float, default=6.0)
    ap.add_argument("--workers", type=int, default=0)
    ap.add_argument("--quick", action="store_true")
    a = ap.parse_args()
    if _nb is None:
        raise SystemExit("numba required")
    rhos = [float(r) for r in a.rhos.split(",")]
    ns, T = a.seeds, a.T
    if a.quick:
        rhos, ns, T = [20.0], 2, 4.0
    arms = ["rest", "v0.3", "vacuum", "k1_0", "k0_1", "k1_1", "k2_1"]
    jobs = [(rho, derive_seed(20261111, "B0", rho, k), arm, T, a.quick) for rho in rhos for k in range(ns) for arm in arms]
    jobs.sort(key=lambda j: -j[0])
    t0 = time.time()
    res = parallel_map(task, jobs, default_workers(a.workers), print,
                       lambda i, r, d, n: f"[B0] {d}/{n}: {r['arm']} rho {r['rho']:g} N {r['N']} ({time.time() - t0:.0f}s)")
    L = [f"2+1, coordinate tau (development proxy), cap 0.68 rho^-1/3; pyramid |x| <= T+8-t, |y| <= T+3-t, T {T:g}"]
    for rho in rhos:
        L.append(f"--- rho {rho:g}")
        for arm in arms:
            R = [r for r in res if r["rho"] == rho and r["arm"] == arm]
            if arm in ("rest", "v0.3"):
                X = np.array([r["xc"] for r in R])
                tm = np.arange(0.5, T, 1.0)
                v = R[0]["v"]
                sl = [np.polyfit(tm, xc, 1)[0] for xc in X if np.all(np.isfinite(xc))]
                rms = np.mean([r["rms"] for r in R], axis=0)
                L.append(f"  wall {arm}: velocity {np.mean(sl):+.4f} +- {np.std(sl, ddof=1) / math.sqrt(len(sl)):.4f} (launch {v}); "
                         f"centre at the end {np.nanmean(X[:, -1]):+.3f} (exact {v * tm[-1]:+.3f}); rms against the 1+1 kink "
                         + " ".join(f"{q:.3f}" for q in rms))
            elif arm == "vacuum":
                L.append(f"  vacuum: R {np.mean([r['R'] for r in R]):.3f} (exact {np.mean([r['R_c'] for r in R]):.3f}); "
                         f"rel. error at the end {np.mean([r['err'] for r in R]):.3f}")
            else:
                om = [r["omega"] for r in R]
                oc = R[0]["omega_c"]
                L.append(f"  plane wave {arm}: omega {np.mean(om):.4f} +- {np.std(om, ddof=1) / math.sqrt(len(om)):.4f} "
                         f"(continuum {oc:.4f}; ratio {np.mean(om) / oc:.4f})")
        L.append(f"  N per run ~{int(np.mean([r['N'] for r in res if r['rho'] == rho])):d}; seconds per run "
                 f"~{np.mean([r['sec'] for r in res if r['rho'] == rho]):.0f}")
    print("\n".join(L))
    os.makedirs(os.path.join("diag", "out"), exist_ok=True)
    fn = os.path.join("diag", "out", f"b0_green2p1_{time.strftime('%Y%m%dT%H%M%S', time.gmtime())}.json")
    with open(fn, "w") as fh:
        json.dump({"args": vars(a), "summary": L, "results": res}, fh)
    print(f"raw records: {fn} ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
