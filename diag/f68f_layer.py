"""
F6.8f-dev, step 0' numerical check (development, decides nothing): the layer of the age cut at second order
========================================================================================================

docs/derivations/f68f_dev.md, step 0'. With the flat age cut s <= l2 (T^2 - r^2), exactly (per direction)

    dJ/dT = int_0^{l2 T^2} ds (T / R) R^2 m(s, R) F(s, R(s)),   R(s) = sqrt(T^2 - s / l2),   m = 1 / (2 l2 dt),

and at r = R(s) the pair has dt = T, so m = 1 / (2 l2 T) there. F6.8e/e3 computed only the bulk term
int ds T^2 m(s, T) F(s, T) (all pairs at r = T, full s range). This file computes both with the stable second-order
instrument of F6.8e3 (same integrand F as its pair_density, same s nodes, same p nodes) and reports the age slope
-a4 rho / k^5 (units of k^2 per unit n) for each, at k = 0.12 and n = 24, 96, 384.

Not included: the explicit b-dependence of the cut (the boundary term of F6.8c), as in F6.8e3.

Uses the 80-bit instrument of F6.8e3: refuses to run without an 80-bit long double (Linux or WSL).

Run from the repository root (one job = one direction p at one n, cached in diag/out/f68f_layer_cache/):
    python -m diag.f68f_layer --ns 24,96,384 --workers 8           bulk and exact-with-layer slopes (16 Nov run)
    python -m diag.f68f_layer --ns 16,20,24 --workers 8            the secant over n = 16 ... 24 (Simpson)
    python -m diag.f68f_layer --ns 24 --bulk-only --p-panels 4x12  direction quadrature check (48 nodes)
    python -m diag.f68f_layer --ns 24 --bulk-only --p-list 0.5,3,10 --n-s 16,24,32,48   s quadrature check
Output diag/out/f68f_layer_<UTC time>.json and .txt. The runs of 16 Nov before this form are
diag/out/f68f_layer.txt, f68f_layer_secant.txt, f68f_pconv.txt and f68f_sconv.txt (the last two by the removed
scripts diag/f68f_pconv.py and diag/f68f_sconv.py, whose checks are the last two commands above).
"""
from __future__ import annotations

import argparse
import math
import os
import time

import numpy as np

from core.stats.parallel import default_workers, parallel_map
from diag import f68f_common as Cm
from phases.F6_order_dynamics.F6_08e3_bd_action_optical_asymptotic_rev2 import (
    second_order_volume_stable, _Ef, _bend, _p_nodes)
from phases.F6_order_dynamics.F6_08a_bd_action_conformal import _A4, _L2, kernel_derivs
from phases.F6_order_dynamics.F6_08b2_bd_action_mixed_rev import _s_nodes

LEV = (24, 16, 16)          # n_sig, n_te, n_phi (lev1 of F6.8e3)
N_S = 16
P_PANELS = [[0.0, 12.0, 10], [12.0, 24.0, 6]]


def p_panels(n):
    """p = n cos(k, pair) runs over [0, n]. F6.8e3 truncates at p = 24 (its tail item I4); for n < 24 the panels end
    at p = n (cos <= 1). The first version used [0, 24] for every n, which for n = 16, 20 added directions with
    cos > 1 (kperp clipped to 0): those rows were wrong."""
    if n >= 24.0:
        return P_PANELS
    return [[0.0, n / 2.0, 10], [n / 2.0, n, 6]]


def integrand(k, CO, r, s):
    """F(s, r) of F6.8e3.pair_density (the optical bb part), per s node, without the measure."""
    s = np.atleast_1d(np.asarray(s, float))
    kpar = k * CO
    tau2 = s / _L2
    dt = np.sqrt(r * r + tau2)
    g0, g1, g2 = kernel_derivs(s * s)
    beta, gam = second_order_volume_stable(k, CO, r, tau2, *LEV)
    E_xy, _ = _Ef(kpar * r)
    Dxy = complex(r * E_xy)
    Bxy = float(_bend(k, r, kpar * r))
    one = 1.0 + np.exp(1j * kpar * r)
    mu = -r * Dxy / (dt * dt)
    mup = -2.0 * r * (0.25 * r + Bxy) / dt ** 2 - 0.5 * abs(Dxy) ** 2 * (dt * dt - 2.0 * r * r) / dt ** 4
    X1 = 3.0 * one + mu
    return ((0.5 * np.abs(X1) ** 2 + mup) * g0 + np.real(X1 * np.conj(beta)) * g1 + gam * g1
            + 0.5 * np.abs(beta) ** 2 * g2)


def task(k, n, p, n_s, bulk_only):
    """Bulk and layer densities of one direction p = n cos at n (ages T = n / k), with n_s s nodes per panel."""
    T = n / k
    CO = p / (k * T)
    s, ws = _s_nodes(np.array([1e12]), n_s, True)
    s, ws = s[0], ws[0]
    keep = ws > 0
    s, ws = s[keep], ws[keep]
    Fb = integrand(k, CO, T, s)
    bulk = float((ws * T * T / (2.0 * _L2 * np.sqrt(T * T + s / _L2)) * Fb).sum())
    if bulk_only:
        return {"bulk": bulk, "layer": None, "n_nodes": int(s.size)}
    lay = 0.0
    for si, wi in zip(s, ws):
        R = math.sqrt(T * T - si / _L2)
        lay += wi * R / (2.0 * _L2) * float(integrand(k, CO, R, [si])[0])
    return {"bulk": bulk, "layer": lay, "n_nodes": int(s.size)}


NAME = "f68f_layer"
CODE = Cm.code_md5([os.path.abspath(__file__)])


def job(*args):
    return Cm.cached(NAME, CODE, task, args)


def directions(n, spec):
    """p nodes and weights: 'auto' = the panels of F6.8e3 (ending at p = n for n < 24); 'MxN' = M equal panels of N
    Gauss nodes over [0, min(n, 24)]."""
    if spec == "auto":
        return _p_nodes(p_panels(n))
    m, q = (int(v) for v in spec.split("x"))
    top = min(n, 24.0)
    e = np.linspace(0.0, top, m + 1)
    return _p_nodes([[e[i], e[i + 1], q] for i in range(m)])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=float, default=0.12)
    ap.add_argument("--ns", type=str, default="24,96,384", help="ages n = kT")
    ap.add_argument("--n-s", type=str, default="16", help="s nodes per panel (a list only with --p-list)")
    ap.add_argument("--p-panels", type=str, default="auto")
    ap.add_argument("--p-list", type=str, default="", help="single directions p, no direction integral")
    ap.add_argument("--bulk-only", action="store_true")
    ap.add_argument("--workers", type=int, default=0)
    a = ap.parse_args()
    Cm.ld_guard()
    t0 = time.time()
    k, ns, nss = a.k, Cm.floats(a.ns), [int(v) for v in a.n_s.split(",")]
    lines = ["F6.8f-dev step 0': the layer of the age cut at second order (flat cut); development, decides nothing", ""]
    if a.p_list:
        jobs = [(k, n, p, m, True) for n in ns for p in Cm.floats(a.p_list) for m in nss]
        res = parallel_map(job, jobs, default_workers(a.workers), print,
                           lambda i, r, d, nn: f"[f68f_layer] {d}/{nn} ({time.time() - t0:.0f}s)")
        lines.append(f"bulk pair density at r = T, single directions, k = {k:g}")
        for (kk, n, p, m, _), r in zip(jobs, res):
            lines.append(f"  n {n:g}  p {p:5.1f}  N_S {m:3d} ({r['n_nodes']} nodes): {r['bulk']:+.10e}")
        data = [{"n": j[1], "p": j[2], "n_s": j[3], **r} for j, r in zip(jobs, res)]
    else:
        if len(nss) != 1:
            raise SystemExit("--n-s takes one value without --p-list")
        nodes = {n: directions(n, a.p_panels) for n in ns}
        jobs = [(k, n, float(p), nss[0], a.bulk_only) for n in ns for p in nodes[n][0]]
        res = parallel_map(job, jobs, default_workers(a.workers), print,
                           lambda i, r, d, nn: f"[f68f_layer] {d}/{nn} ({time.time() - t0:.0f}s)")
        lines += [f"k = {k:g}; slope = -a4 rho / k^5 (units of k^2 per unit n); p panels {a.p_panels}, N_S {nss[0]}", "",
                  "  n      bulk (r = T, F6.8e3)   exact with layer      layer / bulk      (layer/bulk) * n^2"]
        data = []
        for n in ns:
            p, wp = nodes[n]
            R = [r for j, r in zip(jobs, res) if j[1] == n]
            rho_b = 2.0 * math.pi / n * 2.0 * float((wp * np.array([r["bulk"] for r in R])).sum())
            sb = -_A4 * rho_b / k ** 5
            if a.bulk_only:
                lines.append(f"  {n:5g}  {sb:+.10e}")
                data.append({"n": n, "p": p, "w": wp, "bulk": [r["bulk"] for r in R], "slope_bulk": sb})
                continue
            rho_l = 2.0 * math.pi / n * 2.0 * float((wp * np.array([r["layer"] for r in R])).sum())
            sl = -_A4 * rho_l / k ** 5
            d = (sl - sb) / sb
            lines.append(f"  {n:5g}  {sb:+.6e}        {sl:+.6e}      {d:+.3e}       {d * n * n:+.3e}")
            data.append({"n": n, "p": p, "w": wp, "bulk": [r["bulk"] for r in R], "layer": [r["layer"] for r in R],
                         "slope_bulk": sb, "slope_layer": sl})
    lines.append(f"cost {time.time() - t0:.0f} s")
    stem = Cm.save(NAME, vars(a), CODE, data, lines, t0)
    print("\n".join(lines))
    print("saved", stem + ".json")


if __name__ == "__main__":
    main()
