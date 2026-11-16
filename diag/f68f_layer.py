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

Run: python diag/f68f_layer.py [workers] [n list, e.g. 16,20,24] [output tag]   (output diag/out/f68f_layer<tag>.txt)
"""
import math
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from phases.F6_order_dynamics.F6_08e3_bd_action_optical_asymptotic_rev2 import (  # noqa: E402
    second_order_volume_stable, _Ef, _bend, _p_nodes)
from phases.F6_order_dynamics.F6_08a_bd_action_conformal import _A4, _L2, kernel_derivs  # noqa: E402
from phases.F6_order_dynamics.F6_08b2_bd_action_mixed_rev import _s_nodes  # noqa: E402

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


def task(args):
    k, n, p = args
    T = n / k
    CO = p / (k * T)
    s, ws = _s_nodes(np.array([1e12]), N_S, True)
    s, ws = s[0], ws[0]
    keep = ws > 0
    s, ws = s[keep], ws[keep]
    Fb = integrand(k, CO, T, s)
    bulk = float((ws * T * T / (2.0 * _L2 * np.sqrt(T * T + s / _L2)) * Fb).sum())
    lay = 0.0
    for si, wi in zip(s, ws):
        R = math.sqrt(T * T - si / _L2)
        lay += wi * R / (2.0 * _L2) * float(integrand(k, CO, R, [si])[0])
    return (k, n, p, bulk, lay)


def main():
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    k = 0.12
    ns = [float(x) for x in sys.argv[2].split(",")] if len(sys.argv) > 2 else [24.0, 96.0, 384.0]
    tag = sys.argv[3] if len(sys.argv) > 3 else ""
    nodes = {n: _p_nodes(p_panels(n)) for n in ns}
    jobs = [(k, n, float(pp)) for n in ns for pp in nodes[n][0]]
    t0 = time.time()
    cdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "f68f_layer_cache")
    os.makedirs(cdir, exist_ok=True)
    done = {}
    for fn in os.listdir(cdir):
        with open(os.path.join(cdir, fn)) as fh:
            kk, n_, p_, b_, l_ = fh.read().replace("np.float64(", "").replace(")", "").split()
            done[(float(kk), float(n_), float(p_))] = (float(b_), float(l_))
    todo = [a for a in jobs if a not in done]
    print(f"jobs {len(jobs)}, cached {len(jobs) - len(todo)}", flush=True)
    with Pool(workers) as pool:
        for kk, n_, p_, b_, l_ in pool.imap_unordered(task, todo, chunksize=1):
            with open(os.path.join(cdir, f"{n_:.0f}_{p_!r}.txt"), "w") as fh:
                fh.write(f"{float(kk)!r} {float(n_)!r} {float(p_)!r} {float(b_)!r} {float(l_)!r}\n")
            done[(kk, n_, p_)] = (b_, l_)
    res = [a + done[a] for a in jobs]
    lines = ["F6.8f-dev step 0': the layer of the age cut at second order (flat cut), k = 0.12; slope = -a4 rho / k^5 "
             "(units of k^2 per unit n); development, decides nothing", "",
             "  n      bulk (r = T, F6.8e3)   exact with layer      layer / bulk      (layer/bulk) * n^2"]
    for n in ns:
        wp = nodes[n][1]
        rb = np.array([x[3] for x in res if x[1] == n])
        rl = np.array([x[4] for x in res if x[1] == n])
        rho_b = 2.0 * math.pi / n * 2.0 * float((wp * rb).sum())
        rho_l = 2.0 * math.pi / n * 2.0 * float((wp * rl).sum())
        sb, sl = -_A4 * rho_b / k ** 5, -_A4 * rho_l / k ** 5
        d = (sl - sb) / sb
        lines.append(f"  {n:5.0f}  {sb:+.6e}        {sl:+.6e}      {d:+.3e}       {d * n * n:+.3e}")
    lines.append("")
    lines.append("  F6.8c2 full rows (bulk + layer + explicit boundary), k = 0.12, n 16 -> 24: 1.244; F6.8e3 reported bulk "
                 "at n = 24: 1.109")
    lines.append(f"  cost {time.time() - t0:.0f} s, {workers} processes")
    txt = "\n".join(lines)
    print(txt)
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", f"f68f_layer{tag}.txt")
    with open(out, "w") as fh:
        fh.write(txt + "\n")


if __name__ == "__main__":
    main()
