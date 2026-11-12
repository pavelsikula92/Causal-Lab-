"""
Wandering of a field defect on a 1+1 order (development; decides nothing)
=========================================================================

Why (10 Nov 2026; findings F3.14a and F5.2a2; agreed with the author at 18:01)
------------------------------------------------------------------------------
The kink of F3.14a at rest, with no force, wanders. Its displacement spread grows faster than linearly in time, falls
with the density, and is larger for narrower (heavier) kinks. In the F5.2a2 flat arms (2 x 128 sprinklings at rho 800,
2 x 32 at 400) the spread fits a single dimensionless form, sd/w ~ 0.2 (rho w^2)^(-1/2) (t/w)^1.4. The wandering is the
main noise of every test of this branch: both replications (F3.14a, F5.2a2) failed on resolving power set from 16
development sprinklings. It is also a candidate physical effect: a random walk of a bound object in momentum from
discreteness, the "swerves" of Dowker, Henson and Sorkin (2004), never computed for a bound object.

Questions:
  1. Over longer times (T = 30), is the spread a random walk in momentum (Var x ~ t^3), in position (~ t), a spread of
     velocities set early (~ t^2), or a mixture? Fit Var(x(t) - x(0)) = a t + b t^2 + c t^3 (+ a readout floor f).
  2. How do the coefficients scale with the density rho and the width w = 1/sqrt(mu) (mass M = 8 sqrt(mu))? Does the
     spread depend only on rho w^2 and t/w (one dimensionless law)?
  3. Does the kink's momentum diffusion follow the noise of the field around it? The rms of the field in the vacuum
     away from the kink (against the exact vacuum, 0 or 2 pi) is read per slab as the "noise background".
  4. Is there a mean drift (it should be zero by symmetry)?

Model: the rule of F3.14a unchanged (phases/F3_structure/field1d.py), flat 1+1 trapezoid |x| <= X - t with X = T + 10,
kink at rest at x = 0 (exact initial data), mu in {0.25, 1, 4}. Readout: the kink fit per slab of width 1 (tracking the
centre), the local velocity v_s of the fit, the noise background on |x| in [8, 8 + (X - t - 9)] per slab.

Run from the root of the laboratory (numba required):
    python -m diag.wander --quick
    python -m diag.wander --seeds 32 --rhos 200,400,800 --T 30 --workers 8
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


def task(rho, mu, seed, T):
    X = T + 10.0
    t, x = F.sprinkle_trapezoid(rho, X, T, seed)
    t0 = time.time()
    phi = F.evolve(t, x, F.kink_free(t, x, 0.0, mu), rho, mu=mu)
    sec = time.time() - t0
    xc, vs, w, bg = [], [], [], []
    for a in np.arange(0.0, T - 1.0 + 1e-9, 1.0):
        guess = xc[-1] if (xc and np.isfinite(xc[-1])) else 0.0
        r = F.fit_kink(t, x, phi, a, a + 1.0, guess, mu)
        xc.append(r["xc"] if r else float("nan"))
        vs.append(r["vs"] if r else float("nan"))
        w.append(r["w"] if r else float("nan"))
        s = (t >= a) & (t < a + 1.0) & (np.abs(x) > 8.0) & (np.abs(x) < X - a - 1.0)
        if s.sum() > 20:
            vac = np.where(x[s] > 0, 2.0 * math.pi, 0.0)
            bg.append(float(np.sqrt(np.mean((phi[s] - vac) ** 2))))
        else:
            bg.append(float("nan"))
    return {"rho": rho, "mu": mu, "seed": seed, "N": int(t.size), "sec": sec, "xc": xc, "vs": vs, "w": w, "bg": bg}


def fit_var(tt, var):
    """Var = f + a t + b t^2 + c t^3, least squares with weights 1/var (relative errors), non-negative check reported."""
    A = np.column_stack([np.ones_like(tt), tt, tt ** 2, tt ** 3])
    wts = 1.0 / np.maximum(var, 1e-12)
    coef, *_ = np.linalg.lstsq(A * wts[:, None], var * wts, rcond=None)
    return coef


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=32)
    ap.add_argument("--rhos", type=str, default="200,400,800")
    ap.add_argument("--mus", type=str, default="0.25,1,4")
    ap.add_argument("--T", type=float, default=30.0)
    ap.add_argument("--workers", type=int, default=0)
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--base", type=int, default=20261110)
    a = ap.parse_args()
    if F._nb is None:
        raise SystemExit("numba required")
    rhos = [float(r) for r in a.rhos.split(",")]
    mus = [float(m) for m in a.mus.split(",")]
    ns, T = a.seeds, a.T
    if a.quick:
        rhos, mus, ns, T = [100.0], [1.0], 4, 12.0
    jobs = [(rho, mu, derive_seed(a.base, "wander", rho, mu, k), T) for rho in rhos for mu in mus for k in range(ns)]
    jobs.sort(key=lambda j: -j[0])
    t0 = time.time()
    res = parallel_map(task, jobs, default_workers(a.workers), print,
                       lambda i, r, d, n: f"[wander] {d}/{n}: rho {r['rho']:g} mu {r['mu']:g} N {r['N']} "
                                          f"({time.time() - t0:.0f}s)")
    L = [f"flat 1+1, kink at rest, T {T:g}; {ns} sprinklings per (rho, mu); Var(x(t) - x(0.5)) over sprinklings"]
    tm = np.arange(0.5, T, 1.0)
    dt = tm - tm[0]
    for mu in mus:
        w0 = 1.0 / math.sqrt(mu)
        for rho in rhos:
            R = [r for r in res if r["rho"] == rho and r["mu"] == mu]
            X = np.array([r["xc"] for r in R])
            good = np.all(np.isfinite(X), axis=1)
            X = X[good]
            D = X - X[:, [0]]
            var = D.var(axis=0, ddof=1)
            sd = np.sqrt(var)
            idx = [int(k) for k in (2, 5, 10, 20, len(dt) - 1) if k < len(dt)]
            coef = fit_var(dt[1:], var[1:])
            i_a, i_b = np.argmin(np.abs(dt - 0.3 * dt[-1])), len(dt) - 1
            ex = math.log(sd[i_b] / sd[i_a]) / math.log(dt[i_b] / dt[i_a])
            V = np.array([r["vs"] for r in R])[good]
            vv = V.var(axis=0, ddof=1)
            bg = np.nanmean(np.array([r["bg"] for r in R])[good], axis=0)
            L.append(f"--- mu {mu:g} (M {8 * math.sqrt(mu):g}, w {w0:.2f}), rho {rho:g} (rho w^2 {rho * w0 * w0:g}), n {good.sum()}")
            L.append("  sd at dt " + ", ".join(f"{dt[k]:g}: {sd[k]:.4f}" for k in idx) +
                     f"; exponent (dt {dt[i_a]:g} -> {dt[i_b]:g}) {ex:.2f}; mean drift at the end {D[:, -1].mean():+.4f}")
            L.append(f"  Var fit: floor {coef[0]:+.2e}, a (t) {coef[1]:+.2e}, b (t^2) {coef[2]:+.2e}, c (t^3) {coef[3]:+.2e}")
            L.append("  Var of the fitted local velocity v_s at dt " + ", ".join(f"{dt[k]:g}: {vv[k]:.2e}" for k in idx))
            L.append("  noise background rms at dt " + ", ".join(f"{dt[k]:g}: {bg[k]:.4f}" for k in idx))
            L.append(f"  dimensionless: sd/w at t/w = {dt[-1] / w0:.1f}: {sd[-1] / w0:.4f}; "
                     f"0.2 (rho w^2)^(-1/2) (t/w)^1.4 = {0.2 * (rho * w0 * w0) ** -0.5 * (dt[-1] / w0) ** 1.4:.4f}")
    print("\n".join(L))
    os.makedirs(os.path.join("diag", "out"), exist_ok=True)
    fn = os.path.join("diag", "out", f"wander_{time.strftime('%Y%m%dT%H%M%S', time.gmtime())}.json")
    with open(fn, "w") as fh:
        json.dump({"args": vars(a), "summary": L, "results": res}, fh)
    print(f"raw records: {fn} ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
