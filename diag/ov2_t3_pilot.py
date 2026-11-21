"""
OV2 sec. 14, T3 pilot stage 1 (development, decides nothing): the windowed Noether balance on single sprinklings with a
smooth moving kink (the identity <B> = <R> is exact in expectation here; the scatter is the feasibility number).
================================================================================================================
2D rectangle t in [0, T], z in [-Z, Z], Poisson density rho; the pairwise action of T1 with the exact kernel
f(n, eps); field phi_x = Phi(t_x, z_x) = 4 arctan(exp(m gamma (z - v t))), Phi_t analytic; window chi(t) with cos^2
ramps of half-width d at t1 and t2. Per sprinkling (diag/ov2_t3_kernels.noether): B, the pair term P, the ramp parts
Q2 = Q(t2) and Q1 = -Q(t1), R = Q2 + Q1. The continuum energy H = 8 m gamma.
Reported: means and SDs over sprinklings (the sprinkling is the unit), <B> - <R> with its SE, -<Q2>/H and <Q1>/H, and
n = (1.96 s / delta)^2 for delta = 1 % of H, for s = SD(B), SD(R) and SD(B - R).

    python -m diag.ov2_t3_pilot --rho 100 --n 24 --workers 2
"""
import argparse
import json
import math
import time
from multiprocessing import Pool

import numpy as np

from diag.ov2_t3_kernels import future_bits, noether, past_bits


def one(args):
    seed, rho, T, Z, m, vel, eps, t1, t2, d = args
    rng = np.random.default_rng(seed)
    n = rng.poisson(rho * T * 2 * Z)
    t = np.sort(rng.uniform(0, T, n))
    z = rng.uniform(-Z, Z, n)
    gam = 1 / math.sqrt(1 - vel * vel)
    X = m * gam * (z - vel * t)
    phi = 4 * np.arctan(np.exp(X))
    phit = -2 * m * gam * vel / np.cosh(X)
    t0 = time.time()
    P = past_bits(t, z)
    F = future_bits(P)
    B, Pp, Q2, Q1, S, E = noether(t, z, phi, phit, rho, eps, m, t1, t2, d, P, F)
    return {"seed": seed, "N": int(n), "B": B, "P": Pp, "Q2": Q2, "Q1": Q1, "S": S, "sec": time.time() - t0}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rho", type=float, default=100.0)
    ap.add_argument("--T", type=float, default=6.0)
    ap.add_argument("--Z", type=float, default=8.0)
    ap.add_argument("--m", type=float, default=1.0)
    ap.add_argument("--v", type=float, default=0.3)
    ap.add_argument("--eps", type=float, default=0.12)
    ap.add_argument("--win", default="2,4,0.5")
    ap.add_argument("--n", type=int, default=24)
    ap.add_argument("--seed0", type=int, default=0)
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--out", default="")
    q = ap.parse_args()
    t1, t2, d = [float(x) for x in q.win.split(",")]
    jobs = [(q.seed0 + i, q.rho, q.T, q.Z, q.m, q.v, q.eps, t1, t2, d) for i in range(q.n)]
    t0 = time.time()
    with Pool(q.workers) as pool:
        res = pool.map(one, jobs)
    H = 8 * q.m / math.sqrt(1 - q.v ** 2)
    a = {k: np.array([r[k] for r in res]) for k in ("B", "P", "Q2", "Q1", "N")}
    R = a["Q2"] + a["Q1"]
    se = lambda x: x.std(ddof=1) / math.sqrt(x.size)
    nreq = lambda s: (1.96 * s / (0.01 * H)) ** 2
    print(f"rho {q.rho:g} T {q.T:g} Z {q.Z:g} m {q.m:g} v {q.v:g} eps {q.eps:g} window {q.win}; {q.n} sprinklings, "
          f"<N> {a['N'].mean():.0f}; H = {H:.4f}; {time.time() - t0:.0f} s")
    for name, x in (("B", a["B"]), ("R", R), ("B - R", a["B"] - R), ("-Q2", -a["Q2"]), ("Q1", a["Q1"])):
        print(f"  {name:6s} mean {x.mean():+10.4f} +- {se(x):8.4f}   SD {x.std(ddof=1):9.4f}   "
              f"(SD/H {x.std(ddof=1) / H:.3f}; n for 1 % of H: {nreq(x.std(ddof=1)):.0f})")
    print(f"  -<Q2>/H {-a['Q2'].mean() / H:.4f} +- {se(a['Q2']) / H:.4f};  <Q1>/H {a['Q1'].mean() / H:.4f} +- "
          f"{se(a['Q1']) / H:.4f}")
    if q.out:
        json.dump({"args": vars(q), "rows": res}, open(q.out, "w"), indent=1)


if __name__ == "__main__":
    main()
