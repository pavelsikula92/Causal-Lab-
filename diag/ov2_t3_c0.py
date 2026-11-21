"""
OV2 sec. 14, T3 pilot stage 2 (development; identity, sign, normalisation and readout fixed by decision 42).
==========================================================================================================
On each sprinkling (the F3.14b domain of dependence {0 <= t <= T, |x| <= X - t}, density rho, the C0 rule of F3.14b:
phi = phi_0 - (1/2 rho) sum over the causal past of sin phi, mu = 1, launched kink of velocity v):
  C0 arm:     E_x and the ramp terms on the C0 values; (L_xi phi)_x = d_t Phi_fit(x), Phi_fit the kink fitted to the
              run (centre x_c(t) by the F3.14b slab fits, linear in t over the slabs covering the window; velocity
              v_f; gamma_f = 1/sqrt(1 - v_f^2); decision 42);
  smooth arm: the exact kink and its exact Phi_t on the SAME sprinkling (the identity is exact in mean there), so the
              paired differences C0 - smooth carry less of the common noise of the order.
Per sprinkling: B, R = Q2 + Q1, Q2 = Q(t2), Q1 = -Q(t1) for both arms; the continuum-readout energy of F3.14b
(F.energy) at t1 and t2; F-b: the winding (mean phi on [x_c + 5, x_c + 7] minus on [x_c - 7, x_c - 5]) / 2 pi in the
slabs at t1 and t2 (a separate readout).
Expected for a conserved energy: <B_C0> = <R_C0> = 0 within the scatter; -<Q2>/H and <Q1>/H as in the smooth arm.

    python -m diag.ov2_t3_c0 --rho 100 --n 48 --workers 2
Progress lines with an estimate of the remaining time; with --out, the rows so far are saved at each progress line.
"""
import argparse
import json
import math
import time
from multiprocessing import Pool

import numpy as np

from diag.ov2_t3_kernels import noether_multi
from phases.F3_structure import field1d as F


def _winding(t, x, phi, tm, xc):
    s = np.abs(t - tm) < 0.5
    lo = s & (x > xc - 7) & (x < xc - 5)
    hi = s & (x > xc + 5) & (x < xc + 7)
    if lo.sum() < 5 or hi.sum() < 5:
        return float("nan")
    return float((phi[hi].mean() - phi[lo].mean()) / (2 * math.pi))


def one(args):
    seed, rho, X, T, vel, eps, t1, t2, d = args
    t0 = time.time()
    t, x = F.sprinkle_trapezoid(rho, X, T, 1000003 * seed + 17)
    phi = F.evolve(t, x, F.kink_free(t, x, vel), rho)
    # slab fits over the window
    tms, xcs = [], []
    for a in np.arange(math.floor(t1 - d), math.ceil(t2 + d), 1.0):
        r = F.fit_kink(t, x, phi, a, a + 1.0, vel * (a + 0.5))
        if r is not None:
            tms.append(a + 0.5)
            xcs.append(r["xc"])
    if len(tms) < 3:
        return {"seed": seed, "fail": True}
    vf, x0 = np.polyfit(np.array(tms), np.array(xcs), 1)
    gf = 1 / math.sqrt(1 - vf * vf)
    phit_fit = -2 * gf * vf / np.cosh(gf * (x - x0 - vf * t))
    g = 1 / math.sqrt(1 - vel * vel)
    phis = F.kink_exact(t, x, vel)
    phits = -2 * g * vel / np.cosh(g * (x - vel * t))
    # both arms in one pass over the pairs (fast 2D path; equal to the bitset kernel to rounding, diag/ov2_t3_exact.py)
    B, _, Q2, Q1, _ = noether_multi(t, x, np.vstack([phi, phis]), np.vstack([phit_fit, phits]), rho, eps, 1.0,
                                    t1, t2, d)
    Bc, Q2c, Q1c, Bs, Q2s, Q1s = B[0], Q2[0], Q1[0], B[1], Q2[1], Q1[1]
    H1 = F.energy(t, x, phi, t1, x0 + vf * t1)
    H2 = F.energy(t, x, phi, t2, x0 + vf * t2)
    return {"seed": seed, "N": int(t.size), "vf": float(vf), "x0": float(x0),
            "Bc": Bc, "Q2c": Q2c, "Q1c": Q1c, "Bs": Bs, "Q2s": Q2s, "Q1s": Q1s, "H1": H1, "H2": H2,
            "W1": _winding(t, x, phi, t1, x0 + vf * t1), "W2": _winding(t, x, phi, t2, x0 + vf * t2),
            "sec": time.time() - t0}


def summary(rows, vel, header):
    rows = [r for r in rows if not r.get("fail")]
    H = 8 / math.sqrt(1 - vel * vel)
    a = {k: np.array([r[k] for r in rows]) for k in rows[0] if k not in ("seed", "fail")}
    Rc = a["Q2c"] + a["Q1c"]
    Rs = a["Q2s"] + a["Q1s"]
    se = lambda z: z.std(ddof=1) / math.sqrt(z.size)
    nreq = lambda s: (1.96 * s / (0.01 * H)) ** 2
    out = [header + f"; {len(rows)} sprinklings, <N> {a['N'].mean():.0f}; H = 8 gamma = {H:.4f}"]
    for name, z in (("B C0", a["Bc"]), ("R C0", Rc), ("B-R C0", a["Bc"] - Rc), ("B smooth", a["Bs"]),
                    ("R smooth", Rs), ("dB C0-sm", a["Bc"] - a["Bs"]), ("dR C0-sm", Rc - Rs),
                    ("-Q2 C0", -a["Q2c"]), ("Q1 C0", a["Q1c"]), ("-Q2 sm", -a["Q2s"]), ("Q1 sm", a["Q1s"]),
                    ("H(t1) rd", a["H1"]), ("H(t2) rd", a["H2"]), ("H2-H1 rd", a["H2"] - a["H1"]),
                    ("v_fit", a["vf"]), ("W(t1)", a["W1"]), ("W(t2)", a["W2"])):
        out.append(f"  {name:9s} mean {z.mean():+9.4f} +- {se(z):7.4f}   SD {z.std(ddof=1):8.4f}  "
                   f"(SD/H {z.std(ddof=1) / H:.3f}; n(1 % H) {nreq(z.std(ddof=1)):.0f})")
    return "\n".join(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rho", type=float, default=100.0)
    ap.add_argument("--X", type=float, default=14.0)
    ap.add_argument("--T", type=float, default=8.0)
    ap.add_argument("--v", type=float, default=0.3)
    ap.add_argument("--eps", type=float, default=0.12)
    ap.add_argument("--win", default="2.5,5.5,1.0")
    ap.add_argument("--n", type=int, default=48)
    ap.add_argument("--seed0", type=int, default=0)
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--out", default="")
    q = ap.parse_args()
    t1, t2, d = [float(z) for z in q.win.split(",")]
    jobs = [(q.seed0 + i, q.rho, q.X, q.T, q.v, q.eps, t1, t2, d) for i in range(q.n)]
    t0 = time.time()
    rows = []
    with Pool(q.workers) as pool:
        for r in pool.imap_unordered(one, jobs):
            rows.append(r)
            k = len(rows)
            if k % max(1, q.workers) == 0 or k == len(jobs):
                el = time.time() - t0
                print(f"  [{k}/{len(jobs)}] {el:.0f} s elapsed, ~{el / k * (len(jobs) - k):.0f} s left", flush=True)
                if q.out:
                    json.dump({"args": vars(q), "rows": rows, "partial": k < len(jobs)}, open(q.out, "w"), indent=1)
    rows.sort(key=lambda r: r["seed"])
    print(summary(rows, q.v, f"stage 2: rho {q.rho:g} X {q.X:g} T {q.T:g} v {q.v:g} eps {q.eps:g} window {q.win}; "
                             f"{time.time() - t0:.0f} s"), flush=True)
    if q.out:
        json.dump({"args": vars(q), "rows": rows}, open(q.out, "w"), indent=1)


if __name__ == "__main__":
    main()
