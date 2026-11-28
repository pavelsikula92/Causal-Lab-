"""
Q.1-dev -- Johnston's sum-over-chains propagator on a 1+1 sprinkling (development, decides nothing)
====================================================================================================
The first instrument of the quantum layer (docs/quantum/Q1_card_draft.md). A particle on a causal set propagates
from x to y as a sum over chains x = z_0 < z_1 < ... < z_n = y with an amplitude a per relation and b per intermediate
element (S. Johnston, arXiv:0806.3083, 1+1: the causal matrix C):

    K = a C + a^2 b C^2 + ... = a C (1 - a b C)^(-1),   a = 1/2,  b = -m^2 / rho.

For two fixed points with proper time tau and a Poisson sprinkling of density rho between them, the expected number
of chains with k intermediate elements is (rho V)^k / (k!)^2 with V = tau^2 / 2, so E K = (1/2) J_0(m tau) exactly at
every density: the continuum retarded propagator of a scalar of mass m in 1+1. The mean is therefore a code check;
the physics question is the single sprinkling: does K(x, y) of one order approach (1/2) J_0(m tau) as rho grows
(self-averaging), and how fast.

Model: the light-cone square u, v in [0, 1] (u = t + x, v = t - x, area element dt dx = du dv / 2), x = (0, 0) and
y = (1, 1) added as fixed elements (tau = 1, V = 1/2), N ~ Poisson(rho / 2) interior elements. The forward recursion
k_j = a (C_0j + b sum_{i<j} k_i C_ij) in the order of u gives K(x, z) for every z; K(x, y) is the last entry.
Memory O(N): column j of C is (u_i < u_j) & (v_i < v_j) for i < j.

Usage:  python -m diag.q1_dev --n 250,1000,4000 --m 2,5,10 --reps 400 --workers 6
"""

from __future__ import annotations

import argparse
import json
import math
import time
from datetime import datetime, timezone
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy.special import j0

NAME = "q1_dev"


def k_row(u, v, a, b):
    """K(x, z) for the source x = (0, 0) below every element; u sorted increasingly."""
    n = u.size
    k = np.empty(n)
    for j in range(n):
        rel = v[:j] < v[j]                     # u[:j] < u[j] by the sort (ties have probability 0)
        k[j] = a * (1.0 + b * float(k[:j][rel].sum()))
    return k


def one(args):
    n_exp, m, seed = args
    rng = np.random.default_rng(seed)
    rho = 2.0 * n_exp
    n = rng.poisson(n_exp)
    u = rng.random(n); v = rng.random(n)
    o = np.argsort(u); u = np.append(u[o], 1.0); v = np.append(v[o], 1.0)   # y = (1, 1) last
    a, b = 0.5, -m * m / rho
    k = k_row(u, v, a, b)
    tau = np.sqrt(u[:-1] * v[:-1])
    inner = k[:-1] - 0.5 * j0(m * tau)
    return {"K_y": float(k[-1]), "n": int(n), "inner_rms": float(np.sqrt(np.mean(inner ** 2))) if n else 0.0}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", default="250,1000,4000")
    ap.add_argument("--m", default="2,5,10")
    ap.add_argument("--reps", type=int, default=400)
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--seed", type=int, default=20261128)
    a = ap.parse_args()
    t0 = time.time()
    rows = []
    L = [f"Q.1-dev: Johnston's chain propagator on a 1+1 sprinkling (development, decides nothing); tau = 1, "
         f"exact mean (1/2) J0(m tau); {a.reps} sprinklings per point"]
    with Pool(a.workers) as P:
        for m in (float(x) for x in a.m.split(",")):
            for n_exp in (int(x) for x in a.n.split(",")):
                jobs = [(n_exp, m, a.seed + 1000003 * int(m * 10) + 7919 * n_exp + r) for r in range(a.reps)]
                res = P.map(one, jobs)
                K = np.array([r["K_y"] for r in res]); ex = 0.5 * float(j0(m))
                mean, sd = float(K.mean()), float(K.std(ddof=1))
                inner = float(np.mean([r["inner_rms"] for r in res]))
                rows.append({"m": m, "n_exp": n_exp, "mean": mean, "se": sd / math.sqrt(a.reps), "sd": sd,
                             "exact": ex, "inner_rms": inner})
                L.append(f"  m {m:g}  N_exp {n_exp:6d}  (m^2/rho {m*m/(2*n_exp):.4f}): K(x,y) mean {mean:+.5f} +- "
                         f"{sd/math.sqrt(a.reps):.5f}  exact {ex:+.5f}  (mean-exact)/SE {(mean-ex)/(sd/math.sqrt(a.reps)):+.2f}"
                         f" | single-sprinkling sd {sd:.5f}  rms over the interior {inner:.5f}")
    for m in sorted({r["m"] for r in rows}):
        sel = [r for r in rows if r["m"] == m]
        if len(sel) > 1:
            ex = [math.log(sel[i + 1]["sd"] / sel[i]["sd"]) / math.log(sel[i + 1]["n_exp"] / sel[i]["n_exp"])
                  for i in range(len(sel) - 1)]
            L.append(f"  m {m:g}: d ln sd / d ln N = " + ", ".join(f"{e:+.3f}" for e in ex))
    L.append(f"total {time.time() - t0:.0f} s")
    print("\n".join(L))
    out = Path(__file__).resolve().parent / "out"
    out.mkdir(exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    (out / f"{NAME}_{stamp}.txt").write_text("\n".join(L) + "\n")
    (out / f"{NAME}_{stamp}.json").write_text(json.dumps({"args": vars(a), "rows": rows}, indent=1))
    print(f"saved diag/out/{NAME}_{stamp}.json")


if __name__ == "__main__":
    main()
