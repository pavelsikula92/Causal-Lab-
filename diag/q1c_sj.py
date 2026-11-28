"""
Q.1c-sj -- the continuum Sorkin-Johnston state of the 1+1 diamond, independent of any causal set (development)
==============================================================================================================
The reference for Q.1c: the positive part of the continuum Pauli-Jordan operator of a free scalar of mass m on the
light-cone square D = [0,1]^2 (u, v; dV = du dv / 2), by quadrature of the exact continuum kernel. No sprinkling and no
chain sum enters: the kernel values are the continuum G_R itself.

Kernel: Delta(x, y) = G_R(x, y) - G_R(y, x), G_R(x, y) = (1/2) J0(m tau_xy) for x < y; on a null line (du = 0 or
dv = 0, not both) half the jump. Quadrature: the product trapezoid rule on the nested grids u, v = i/n (weights w).
With B = W^(1/2) Delta W^(1/2) real antisymmetric, Pos(iB) = (|B| + iB)/2, |B| = (B^T B)^(1/2), and the SJ function
at the nodes is W^(-1/2) Pos(iB) W^(-1/2): Re W_SJ = W^(-1/2) |B| W^(-1/2) / 2 (and Im W_SJ = Delta/2 exactly).
Probe pairs are node pairs (offsets k/24 around the centre), common to the nested grids n = 24, 48, 96; the trapezoid
rule with a kernel that jumps on the light cone converges as 1/n, so the limit is taken by Richardson extrapolation.
A first version with Gauss nodes and the Nystrom extension of the eigenfunctions did not converge (division by small
eigenvalues; record 20261128T104703Z) and was replaced.

Usage: python -m diag.q1c_sj --m 5,10,20 --ns 24,48,96
"""

from __future__ import annotations

import argparse
import json
import math
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.special import j0, k0, y0

NAME = "q1c_sj"
BASE = 24                                     # probe offsets are multiples of 1/BASE


def probes():
    """Probe pairs (u, v in units of 1/BASE from the centre index): timelike tau = d, d/sqrt2 shapes, spacelike."""
    out = []
    for k in (1, 2, 3, 4):
        out.append(((-k, -k), (k, k), "time"))            # du = dv = 2k/BASE
        out.append(((-k, k), (k, -k), "space"))
        out.append(((-k, 0), (k, 0), "null"))             # null separated (du = 2k, dv = 0): reported
        out.append(((-2 * k, -k), (0, 0), "time"))        # du = 2k, dv = k
    return out


def kernel(du, dv, m):
    t2 = du * dv
    g = 0.5 * j0(m * np.sqrt(np.abs(t2)))
    fwd = (du > 0) & (dv > 0); bwd = (du < 0) & (dv < 0)
    nf = ((du == 0) & (dv > 0)) | ((dv == 0) & (du > 0))
    nb = ((du == 0) & (dv < 0)) | ((dv == 0) & (du < 0))
    return np.where(fwd, g, 0.0) - np.where(bwd, g, 0.0) + 0.5 * (np.where(nf, g, 0.0) - np.where(nb, g, 0.0))


def sj(m, n, pr):
    x = np.arange(n + 1) / n
    w1 = np.full(n + 1, 1.0 / n); w1[0] = w1[-1] = 0.5 / n
    U, V = np.meshgrid(x, x, indexing="ij"); U = U.ravel(); V = V.ravel()
    w = np.outer(w1, w1).ravel() / 2.0
    B = kernel(U[None, :] - U[:, None], V[None, :] - V[:, None], m)     # Delta(x_i, x_j)
    sw = np.sqrt(w)
    B *= sw[:, None]; B *= sw[None, :]
    lam, E = np.linalg.eigh(B.T @ B)
    del B
    absB = (E * np.sqrt(np.clip(lam, 0.0, None))) @ E.T
    del E
    c = n // 2; f = n // BASE
    res = []
    for (p, q, kind) in pr:
        ip = (c + p[0] * f) * (n + 1) + (c + p[1] * f); iq = (c + q[0] * f) * (n + 1) + (c + q[1] * f)
        reW = 0.5 * absB[ip, iq] / (sw[ip] * sw[iq])
        du = (q[0] - p[0]) / BASE; dv = (q[1] - p[1]) / BASE; s = math.sqrt(abs(du * dv))
        mink = (-0.25 * y0(m * s) if kind == "time" else k0(m * s) / (2 * math.pi) if kind == "space" else float("nan"))
        res.append({"kind": kind, "du": du, "dv": dv, "s": s, "ms": m * s, "ReW": float(reW), "mink": float(mink)})
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--m", default="5,10,20")
    ap.add_argument("--ns", default="24,48,96")
    a = ap.parse_args()
    t0 = time.time()
    pr = probes()
    ns = [int(z) for z in a.ns.split(",")]
    L = ["Q.1c-sj: the continuum SJ state of the diamond (trapezoid quadrature of the exact kernel; development); "
         "Re W_SJ at probe node pairs for n = " + ", ".join(map(str, ns)) + ", Richardson limit (1/n), and "
         "(limit - Im G_F Minkowski) / (1/4)"]
    data = []
    for m in (float(z) for z in a.m.split(",")):
        rows = {}
        for n in ns:
            t1 = time.time(); rows[n] = sj(m, n, pr)
            print(f"  m {m:g} n {n}: {time.time() - t1:.0f} s", flush=True)
        L.append(f"  m {m:g}:")
        for i, p in enumerate(pr):
            vals = [rows[n][i]["ReW"] for n in ns]
            lim = 2 * vals[-1] - vals[-2]                         # Richardson for an O(1/n) error
            r = rows[ns[-1]][i]
            L.append(f"    {r['kind']:5s} du {r['du']:.4f} dv {r['dv']:+.4f} m s {r['ms']:6.3f}: Re W_SJ "
                     + " / ".join(f"{v:+.6f}" for v in vals) + f"  -> limit {lim:+.6f}  Minkowski {r['mink']:+.6f}  "
                     f"offset/(1/4) {(lim - r['mink'])/0.25:+.4f}")
        data.append({"m": m, "rows": {str(n): rows[n] for n in ns}})
        print("\n".join(L[-len(pr) - 1:]), flush=True)
    L.append(f"total {time.time() - t0:.0f} s")
    out = Path(__file__).resolve().parent / "out"; out.mkdir(exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    (out / f"{NAME}_{stamp}.txt").write_text("\n".join(L) + "\n")
    (out / f"{NAME}_{stamp}.json").write_text(json.dumps({"args": vars(a), "probes": pr, "data": data}))
    print(L[-1]); print(f"saved diag/out/{NAME}_{stamp}.json")


if __name__ == "__main__":
    main()
