"""
Q.1c-dev -- Johnston's Feynman propagator on a 1+1 sprinkling (development, decides nothing)
============================================================================================
Johnston (arXiv:0909.0944): K_R = aC(1 - abC)^(-1) (Q.1), K_A = K_R^T, Delta = K_R - K_A; i Delta is Hermitian with
eigenvalues in +- pairs; Q = the positive part (sum over lambda > 0 of lambda u u^dagger); K_F = K_R + i Q. Then
Re K_F = (K_R + K_A)/2 (exact in the mean, from Q.1) and Im K_F = Re Q, the non-linear (spectral) part, whose mean
is NOT exact. For a real antisymmetric Delta, Q = (|Delta| + i Delta)/2 with |Delta| = (Delta^T Delta)^(1/2), so
Re Q = |Delta|/2: one real symmetric eigendecomposition.

Continuum (Minkowski vacuum, 2D): G_F = (1/4) H0^(2)(m s); Im G_F = -(1/4) Y0(m tau) for timelike pairs,
+(1/2 pi) K0(m sigma) for spacelike pairs. The SJ state of a finite diamond differs from the Minkowski vacuum near
the boundary and when m L is not large; the comparison is made for pairs in the centre of the diamond.

Model: the light-cone square u, v in [0, 1] (dt dx = du dv / 2), N ~ Poisson(N_exp); the centre: u, v in
[0.5 - h, 0.5 + h]; pairs binned in m s (s = sqrt|tau^2|), timelike and spacelike apart.

Usage: python -m diag.q1c_dev --n 1000,2000,4000 --m 5,10,20 --reps 20 --workers 2
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
from scipy.special import k0, y0

NAME = "q1c_dev"
BINS = np.array([0.25, 0.35, 0.5, 0.7, 1.0, 1.4, 2.0, 2.8, 4.0])


def one(args):
    n_exp, m, seed, h = args
    rng = np.random.default_rng(seed)
    n = rng.poisson(n_exp)
    u = rng.random(n); v = rng.random(n)
    o = np.argsort(u); u = u[o]; v = v[o]
    rho = 2.0 * n_exp
    a, b = 0.5, -m * m / rho
    C = ((u[:, None] < u[None, :]) & (v[:, None] < v[None, :])).astype(float)   # C[i, j] = 1 if i < j
    KR = a * np.linalg.solve((np.eye(n) - a * b * C).T, C.T).T                 # K_R = a C (1 - abC)^(-1)
    D = KR - KR.T
    w, V = np.linalg.eigh(D.T @ D)
    absD = (V * np.sqrt(np.clip(w, 0.0, None))) @ V.T
    ReQ = 0.5 * absD
    c = np.where((np.abs(u - 0.5) < h) & (np.abs(v - 0.5) < h))[0]
    du = u[c][:, None] - u[c][None, :]; dv = v[c][:, None] - v[c][None, :]
    t2 = du * dv                                                               # tau^2 (>0 timelike)
    s = np.sqrt(np.abs(t2)); ms = m * s
    iu = np.triu_indices(c.size, 1)
    ms, t2, q = ms[iu], t2[iu], ReQ[np.ix_(c, c)][iu]
    out = {}
    for kind, sel in (("time", t2 > 0), ("space", t2 < 0)):
        idx = np.digitize(ms[sel], BINS) - 1
        qq, xx = q[sel], ms[sel]
        cont = -0.25 * y0(xx) if kind == "time" else k0(xx) / (2 * math.pi)
        for k in range(len(BINS) - 1):
            mk = idx == k
            if mk.sum() >= 5:
                out[f"{kind}{k}"] = (float(qq[mk].mean()), float(cont[mk].mean()), int(mk.sum()))
    diag = float(np.mean(np.diag(ReQ)[c]))
    return {"bins": out, "n": int(n), "diagQ": diag}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", default="1000,2000,4000")
    ap.add_argument("--m", default="5,10,20")
    ap.add_argument("--reps", type=int, default=20)
    ap.add_argument("--h", type=float, default=0.15)
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--seed", type=int, default=20261128)
    a = ap.parse_args()
    t0 = time.time()
    L = [f"Q.1c-dev: Johnston's Feynman propagator, Im K_F = Re Q against the Minkowski vacuum (development, decides "
         f"nothing); centre half-width {a.h}; {a.reps} sprinklings per point; deviation (mean Re Q - continuum) / (1/4) per bin of m s"]
    rows = []
    with Pool(a.workers) as P:
        for m in (float(x) for x in a.m.split(",")):
            for n_exp in (int(x) for x in a.n.split(",")):
                jobs = [(n_exp, m, a.seed + 1000003 * int(m) + 7919 * n_exp + r, a.h) for r in range(a.reps)]
                res = P.map(one, jobs)
                line = f"  m {m:g} N_exp {n_exp} (m^2/rho {m*m/(2*n_exp):.4f}):"
                row = {"m": m, "n_exp": n_exp, "bins": {}}
                for kind in ("time", "space"):
                    parts = []
                    for k in range(len(BINS) - 1):
                        key = f"{kind}{k}"
                        vals = [r["bins"][key] for r in res if key in r["bins"]]
                        if len(vals) < 3:
                            continue
                        rq = np.array([v[0] for v in vals]); ct = np.array([v[1] for v in vals])
                        d = (rq - ct) / 0.25                  # deviation in units of the amplitude 1/4 of G_F
                        dm, se, sd1 = d.mean(), d.std(ddof=1) / math.sqrt(len(d)), d.std(ddof=1)
                        row["bins"][key] = {"ms_lo": float(BINS[k]), "dev": float(dm), "se": float(se),
                                            "sd1": float(sd1), "cont": float(ct.mean())}
                        parts.append(f"{BINS[k]:.2g}-{BINS[k+1]:.2g}: {dm:+.4f}+-{se:.4f} (sd1 {sd1:.4f}; cont {ct.mean()/0.25:+.3f})")
                    line += f"\n      {kind:5s} " + "  ".join(parts)
                rows.append(row)
                L.append(line)
                print(line, flush=True)
    L.append(f"total {time.time() - t0:.0f} s")
    out = Path(__file__).resolve().parent / "out"; out.mkdir(exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    (out / f"{NAME}_{stamp}.txt").write_text("\n".join(L) + "\n")
    (out / f"{NAME}_{stamp}.json").write_text(json.dumps({"args": vars(a), "rows": rows}))
    print(L[-1]); print(f"saved diag/out/{NAME}_{stamp}.json")


if __name__ == "__main__":
    main()
