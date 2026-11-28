"""
Q.1c-dev2 -- Im K_F of one sprinkled order at fixed probe pairs against the continuum SJ state of the same diamond
==================================================================================================================
(development, decides nothing). The probe points of diag/q1c_sj.py (centre +- k/24, timelike and spacelike pairs; the
null pairs are left out) are added as fixed elements to every Poisson sprinkling of the light-cone square, as x and y
in Q.1. K_F = K_R + iQ is built from the order alone (Johnston, arXiv:0909.0944): K_R = aC(1 - abC)^(-1), a = 1/2,
b = -m^2/rho, Delta = K_R - K_R^T, Q = Pos(i Delta) = (|Delta| + i Delta)/2, Im K_F = Re Q = |Delta|/2, no truncation
of the spectrum. The reference W_SJ(p, q) is the continuum SJ state of record diag/out/q1c_sj_20261128T105641Z,
extrapolated in n with the order estimated from the three grids; its uncertainty is the size of the last step.

Readout per (m, N_exp, pair): d = (Re Q - W_SJ) / (1/4), mean and sd over realisations.

Probes in general position (revision 2, 28 Nov). Inserting all probe points at once put some of them exactly on the
light cones of others (equal u or equal v: for example (-1, 1) is null to both (-1, -1) and (1, 1)). Such
coincidences have probability zero in a sprinkling, but two of them moved Re Q of the pair (-1,-1)-(1,1) by -0.32 (in
units of 1/4) at m 10, N_exp 2000 (alone: -0.010 +- 0.018; with (0, 0): -0.021; with (-1, 1), (1, -1): -0.316). The
pairs are therefore inserted in five groups whose points share no u and no v, one group per sprinkling (record
20261128T115328Z is the first, flawed version).

Usage: python -m diag.q1c_dev2 --n 1000,2000,4000 --m 5,10,20 --reps 20 --workers 2
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

NAME = "q1c_dev2"
SJ_REC = Path(__file__).resolve().parent / "out" / "q1c_sj_20261128T105641Z.json"
BASE = 24


def reference():
    d = json.loads(SJ_REC.read_text())
    pr = [tuple(map(tuple, p[:2])) + (p[2],) for p in d["probes"]]
    ref = {}
    for blk in d["data"]:
        m = blk["m"]; ns = sorted(blk["rows"], key=int)
        for i, p in enumerate(pr):
            if p[2] == "null":
                continue
            v = [blk["rows"][n][i]["ReW"] for n in ns]
            d1, d2 = v[1] - v[0], v[2] - v[1]
            r = d2 / d1 if d1 != 0 and 0 < d2 / d1 < 1 else 0.5
            lim = v[2] + d2 * r / (1 - r)
            ref[(m, i)] = {"W": lim, "unc": abs(lim - v[2]), "kind": p[2], "ms": blk["rows"][ns[-1]][i]["ms"],
                           "p": p[0], "q": p[1]}
    return pr, ref


GROUPS = [  # (kind, k) per group; in each group all u and all v of the probe points are distinct
    [("skew", 1), ("time", 3), ("space", 4)],
    [("skew", 2), ("time", 1), ("space", 3)],
    [("skew", 3), ("time", 2), ("space", 1)],
    [("skew", 4), ("time", 3), ("space", 2)],
    [("time", 4)],
]


def pair_index(kind, k):
    """Index of the probe pair in diag/q1c_sj.probes(): per k the order is time, space, null, skew."""
    return 4 * (k - 1) + {"time": 0, "space": 1, "null": 2, "skew": 3}[kind]


def one(args):
    n_exp, m, seed, pts = args
    rng = np.random.default_rng(seed)
    n = rng.poisson(n_exp)
    u = np.concatenate([rng.random(n), [0.5 + a / BASE for a, b in pts]])
    v = np.concatenate([rng.random(n), [0.5 + b / BASE for a, b in pts]])
    o = np.argsort(u, kind="stable"); u = u[o]; v = v[o]
    pos = {k: int(np.where(o == n + k)[0][0]) for k in range(len(pts))}
    N = u.size; rho = 2.0 * n_exp
    a, b = 0.5, -m * m / rho
    C = ((u[:, None] < u[None, :]) & (v[:, None] < v[None, :])).astype(float)
    KR = a * np.linalg.solve((np.eye(N) - a * b * C).T, C.T).T
    D = KR - KR.T
    lam, E = np.linalg.eigh(D.T @ D)
    absD = (E * np.sqrt(np.clip(lam, 0.0, None))) @ E.T
    return {"ReQ": {f"{i},{j}": float(0.5 * absD[pos[i], pos[j]]) for i in range(len(pts)) for j in range(len(pts))},
            "n": int(n)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", default="1000,2000,4000")
    ap.add_argument("--m", default="5,10,20")
    ap.add_argument("--reps", type=int, default=20)
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--seed", type=int, default=20261128)
    a = ap.parse_args()
    t0 = time.time()
    pr, ref = reference()
    L = ["Q.1c-dev2: Re Q (= Im K_F) of one sprinkled order at fixed probe pairs against the continuum SJ state of the "
         f"same diamond (development, decides nothing); {a.reps} realisations per point; d = (Re Q - W_SJ)/(1/4)"]
    rows = []
    with Pool(a.workers) as P:
        for m in (float(x) for x in a.m.split(",")):
            for n_exp in (int(x) for x in a.n.split(",")):
                L.append(f"  m {m:g} N_exp {n_exp} (m^2/rho {m*m/(2*n_exp):.4f}):")
                got = {}
                for g, grp in enumerate(GROUPS):
                    ids = [pair_index(kd, k) for kd, k in grp]
                    pts = sorted({p for i in ids for p in pr[i][:2]})
                    assert len({p[0] for p in pts}) == len(pts) and len({p[1] for p in pts}) == len(pts)
                    idx = {p: j for j, p in enumerate(pts)}
                    jobs = [(n_exp, m, a.seed + 1000003 * int(m) + 7919 * n_exp + 104729 * g + r, pts)
                            for r in range(a.reps)]
                    res = P.map(one, jobs)
                    for i in ids:
                        key = f"{idx[pr[i][0]]},{idx[pr[i][1]]}"
                        got[i] = [r["ReQ"][key] for r in res]
                for i in sorted(got):
                    r0 = ref[(m, i)]
                    dd = np.array([(x - r0["W"]) / 0.25 for x in got[i]])
                    rows.append({"m": m, "n_exp": n_exp, "pair": i, "kind": r0["kind"], "ms": r0["ms"],
                                 "W": r0["W"], "unc": r0["unc"] / 0.25, "d": dd.tolist()})
                    L.append(f"    {r0['kind']:5s} m s {r0['ms']:6.3f}: d {dd.mean():+.4f} +- {dd.std(ddof=1)/math.sqrt(dd.size):.4f}"
                             f" (sd one order {dd.std(ddof=1):.4f}; reference unc {r0['unc']/0.25:.4f}; W_SJ/(1/4) {r0['W']/0.25:+.4f})")
                print("\n".join(L[-1 - len(got):]), flush=True)
    L.append(f"total {time.time() - t0:.0f} s")
    out = Path(__file__).resolve().parent / "out"
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    (out / f"{NAME}_{stamp}.txt").write_text("\n".join(L) + "\n")
    (out / f"{NAME}_{stamp}.json").write_text(json.dumps({"args": vars(a), "rows": rows}))
    print(L[-1]); print(f"saved diag/out/{NAME}_{stamp}.json")


if __name__ == "__main__":
    main()
