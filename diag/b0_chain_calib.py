"""
B0, development stage 2: calibration of the longest-chain estimator of proper time in 2+1 (decides nothing)
===========================================================================================================
Decision 16, D1 (11 Nov 2026): tau_hat(x, y) = L_max(x, y) / (alpha_3 rho^(1/3)); alpha_3 to be calibrated before
freezing, also at several levels of thinning, so that a change of density is not taken for a change of proper time.
The interval count N ~ (pi/12) rho tau^3 is a diagnostic only.

For random timelike pairs y < x of a flat 2+1 Poisson sprinkling of density rho (thinned uniformly by a factor p:
realised density p rho), the longest chain L between them is computed from the order alone: dynamic programming over the
events of the causal interval, sorted by time, with the relation as the only input (coordinates only to find the
interval and to read the true tau for the calibration). Reported per (rho, p):
  alpha = L / ((p rho)^(1/3) tau) in bins of (p rho)^(1/3) tau (finite-size dependence);
  the spread of tau_hat / tau;
  E[1 / tau_hat] against 1 / tau (the weight of the Green function), and the interval count against (pi/12) p rho tau^3.
Convention: L counts the links of the chain (a link has L = 1).

Run from the root of the laboratory (numba required):
    python -m diag.b0_chain_calib --quick
    python -m diag.b0_chain_calib --pairs 400 --workers 6
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
    @_nb.njit(cache=True)
    def _longest(t, x, y):
        """Longest chain (number of links) from event 0 to event m-1 among events sorted by t; relation by coordinates
        (Minkowski 2+1) -- the relation of the sprinkling, the only input of the order."""
        m = t.size
        best = np.full(m, -1, np.int64)
        best[0] = 0
        for i in range(1, m):
            b = -1
            for j in range(i):
                if best[j] >= 0:
                    dt = t[i] - t[j]
                    if dt > 0.0:
                        dx = x[i] - x[j]
                        dy = y[i] - y[j]
                        if dt * dt > dx * dx + dy * dy and best[j] + 1 > b:
                            b = best[j] + 1
            best[i] = b
        return best[m - 1]


def task(rho, p, seed, n_pairs, tau_max):
    rng = np.random.default_rng(seed)
    out = []
    lam = p * rho
    for _ in range(n_pairs):
        tau = rng.uniform(0.05, tau_max) if rng.uniform() < 0.7 else rng.uniform(0.05, 0.6)
        eta = rng.uniform(0, 1.5)                       # rapidity of the pair
        T = tau * math.cosh(eta)
        R = tau * math.sinh(eta)
        th = rng.uniform(0, 2 * math.pi)
        ex, ey = R * math.cos(th), R * math.sin(th)
        # sprinkle the interval (in the rest frame of the pair, then boost: the sprinkling is Lorentz invariant, so
        # sprinkling the rest-frame interval is the same process)
        vol = math.pi * tau ** 3 / 12.0
        n = rng.poisson(lam * vol)
        pts = []
        while len(pts) < n:
            m = 4 * (n - len(pts)) + 10
            tt = rng.uniform(0, tau, m)
            xx = rng.uniform(-tau / 2, tau / 2, m)
            yy = rng.uniform(-tau / 2, tau / 2, m)
            r = np.hypot(xx, yy)
            ok = (r < tt) & (r < tau - tt)
            for a, b, c in zip(tt[ok], xx[ok], yy[ok]):
                pts.append((a, b, c))
        pts = np.array(pts[:n]).reshape(-1, 3)
        tt = np.concatenate([[0.0], pts[:, 0], [tau]])
        xx = np.concatenate([[0.0], pts[:, 1], [0.0]])
        yy = np.concatenate([[0.0], pts[:, 2], [0.0]])
        o = np.argsort(tt, kind="stable")
        L = int(_longest(tt[o], xx[o], yy[o]))
        out.append((tau, L, n))
    return {"rho": rho, "p": p, "pairs": out}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", type=int, default=400)
    ap.add_argument("--workers", type=int, default=0)
    ap.add_argument("--quick", action="store_true")
    a = ap.parse_args()
    if _nb is None:
        raise SystemExit("numba required")
    cells = [(25.0, 1.0), (50.0, 1.0), (100.0, 1.0), (100.0, 0.5), (100.0, 0.25), (200.0, 0.25)]
    npairs, tau_max, nj = a.pairs, 4.0, 8
    if a.quick:
        cells, npairs, tau_max, nj = [(50.0, 1.0), (100.0, 0.5)], 40, 2.5, 2
    jobs = [(rho, p, derive_seed(20261111, "chain", rho, p, k), npairs // nj, tau_max) for rho, p in cells for k in range(nj)]
    t0 = time.time()
    res = parallel_map(task, jobs, default_workers(a.workers), print,
                       lambda i, r, d, n: f"[chain] {d}/{n} ({time.time() - t0:.0f}s)")
    L_ = ["longest chain L (links) between y < x against the true tau; realised density p rho"]
    for rho, p in cells:
        P = np.array([q for r in res if r["rho"] == rho and r["p"] == p for q in r["pairs"]], float)
        tau, L, n = P[:, 0], P[:, 1], P[:, 2]
        s = (p * rho) ** (1.0 / 3.0) * tau
        bins = [0, 1, 2, 4, 8, 100]
        parts = []
        for lo, hi in zip(bins[:-1], bins[1:]):
            m = (s >= lo) & (s < hi)
            if m.sum() > 3:
                parts.append(f"[{lo},{hi}): alpha {np.mean(L[m] / s[m]):.3f} (n {m.sum()})")
        big = s >= 4
        al = float(np.mean(L[big] / s[big])) if big.sum() > 3 else float("nan")
        th = L / (al * (p * rho) ** (1.0 / 3.0))
        Ninv = np.mean(n / (math.pi / 12 * p * rho * tau ** 3))
        L_.append(f"--- rho {rho:g}, p {p:g} (realised density {p * rho:g}): alpha (s >= 4) {al:.3f}; " + "; ".join(parts))
        L_.append(f"  tau_hat/tau for s >= 4: mean {np.mean(th[big] / tau[big]):.3f}, sd {np.std(th[big] / tau[big]):.3f}; "
                  f"links (L = 1) {np.mean(L == 1):.3f} of the pairs; mean 1/tau_hat / mean 1/tau for s < 2: "
                  f"{np.mean(1 / th[s < 2]) / np.mean(1 / tau[s < 2]):.3f}; interval count / (pi/12 p rho tau^3) {Ninv:.3f}")
    print("\n".join(L_))
    os.makedirs(os.path.join("diag", "out"), exist_ok=True)
    fn = os.path.join("diag", "out", f"b0_chain_calib_{time.strftime('%Y%m%dT%H%M%S', time.gmtime())}.json")
    with open(fn, "w") as fh:
        json.dump({"args": vars(a), "summary": L_, "results": res}, fh)
    print(f"raw records: {fn} ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
