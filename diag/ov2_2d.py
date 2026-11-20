"""
OV2 in 2D, development (decides nothing): tests T0 and T1 of docs/derivations/order_variation_OV2.md.
=====================================================================================================
Order: the light-cone square [0, 1]^2, ds^2 = -du dv (t = (u + v)/2, x = (v - u)/2), sqrt|g| = 1/2; a Poisson
sprinkling of proper density rho = 2 N_exp, so N_exp du dv elements per coordinate area. A conformal factor
g -> e^{2 sigma} g does not change the order, only the density: rho e^{2 sigma}.

Matter at FIXED configuration (decision 39, item 3): a static sine-Gordon kink phi(x) = 4 arctan(e^{m x}),
V = m^2 (1 - cos phi), centred on the diagonal of the diamond.

Matter action S_m = S_kin + S_pot, S_pot = -(1/rho) sum_v V(phi_v). The 2D smeared BD operator (Surya 2012 convention,
S_BD = -(2/rho) sum B1 = 4 eps (N - 2 eps sum N_n f)): (B phi)(x) = 4 eps rho [-phi(x)/2 + eps sum_{y < x} f(n) phi(y)].
Three kinetic forms, all S_kin = (1/rho) sum_x (1/2) phi_x (K phi)_x in the bulk:
  plain       K = B                      -> S_kin = 2 eps sum_x [-phi_x^2/2] + 2 eps^2 P
  row-diff    K = B - diag(B1)           -> S_kin = 2 eps^2 R      (R = sum f (phi_y - phi_x) phi_x, x later)
  Dirichlet   pairwise                   -> S_kin = -eps^2 D       (D = sum f (phi_x - phi_y)^2)
The Dirichlet form is the symmetric form of the row-difference operator without its time-asymmetric diagonal
2 eps^2 rho (F_future - F_past); only it is exactly invariant under phi -> phi + 2 pi k (T0).
Continuum: S_kin = -(1/2) int (d phi)^2 sqrt|g| = int phi_u phi_v du dv (conformally invariant in 2D);
S_pot = -(1/2) int e^{2 sigma} V du dv.

T0: S(phi + 2 pi) - S(phi) for the three forms, on the flat diamond and with a conformal bump.
T1: sigma = a exp(-((u - 1/2)^2 + (v - 1/2)^2)/(2 s^2)); coupled sprinklings (one master at the largest density,
    thinned): delta S_kin and delta S_pot between +a, 0, -a. Expected at first order: delta S_kin = 0 (conformal
    invariance in 2D), delta S_pot = -(1/2) int (e^{2 sigma} - 1) V du dv, i.e. delta S_m = int sigma T^a_a sqrt|g|
    with T^a_a = -2V.

    python -m diag.ov2_2d --eps 0.12 --nexp 16640 --seeds 0-15 --workers 2
"""
from __future__ import annotations

import argparse
import json
import math
import os
import time

import numpy as np

from diag.ov1_2d import fker, n_cut_of, ranks
from diag.ov2_kernels import matter_sums

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "ov2_2d")


def kink(u, v, m):
    x = 0.5 * (v - u)
    return 4.0 * np.arctan(np.exp(m * x))


def sigma_of(u, v, a, s):
    return a * np.exp(-((u - 0.5) ** 2 + (v - 0.5) ** 2) / (2 * s * s))


def actions(u, v, phi, eps, n_exp, f, nc):
    o = np.argsort(u)
    u, v, phi = u[o], v[o], phi[o]
    vr, iv = ranks(u, v)
    D, P, R, Fp = matter_sums(vr, iv, nc, f, phi)
    rho = 2.0 * n_exp
    return {"plain": 2 * eps * float(np.sum(-0.5 * phi ** 2)) + 2 * eps * eps * P,
            "rowdiff": 2 * eps * eps * R, "dirichlet": -eps * eps * D, "N": int(u.size)}


def continuum(m, a, s, n=1200):
    g = (np.arange(n) + 0.5) / n
    U, Vv = np.meshgrid(g, g, indexing="ij")
    x = 0.5 * (Vv - U)
    ph_x = 4 * m * np.exp(m * x) / (1 + np.exp(2 * m * x))           # d phi / dx
    phu, phv = -0.5 * ph_x, 0.5 * ph_x
    Vp = m * m * (1 - np.cos(kink(U, Vv, m)))
    dA = 1.0 / n ** 2
    out = {"kin": float(np.sum(phu * phv) * dA), "pot0": float(-0.5 * np.sum(Vp) * dA)}
    for sgn in (1, -1):
        sg = sigma_of(U, Vv, sgn * a, s)
        out[f"dpot{sgn:+d}"] = float(-0.5 * np.sum((np.exp(2 * sg) - 1) * Vp) * dA)
        out[f"first{sgn:+d}"] = float(np.sum(sg * (-2 * Vp) * 0.5) * dA)     # int sigma T^a_a sqrt|g|
    return out


def task(eps, n_exp, seed, m, a, s):
    t0 = time.time()
    nc = n_cut_of(eps)
    f = fker(eps, nc + 2)
    rng = np.random.default_rng(seed)
    smax = a                                                   # sigma <= a
    n_master = rng.poisson(n_exp * math.exp(2 * smax))
    u, v = rng.uniform(0, 1, n_master), rng.uniform(0, 1, n_master)
    r = rng.uniform(0, 1, n_master)
    sg = sigma_of(u, v, a, s)
    keep = {+1: r < np.exp(2 * sg - 2 * smax), 0: r < np.exp(-2 * smax), -1: r < np.exp(-2 * sg - 2 * smax)}
    rho = 2.0 * n_exp
    out = {"seed": seed}
    for k, sel in keep.items():
        uu, vv = u[sel], v[sel]
        ph = kink(uu, vv, m)
        A = actions(uu, vv, ph, eps, n_exp, f, nc)
        A2 = actions(uu, vv, ph + 2 * math.pi, eps, n_exp, f, nc)
        Vp = m * m * (1 - np.cos(ph))
        out[str(k)] = {"N": A["N"], "pot": float(-np.sum(Vp) / rho),
                       **{f: A[f] for f in ("plain", "rowdiff", "dirichlet")},
                       **{f"shift_{f}": A2[f] - A[f] for f in ("plain", "rowdiff", "dirichlet")}}
    out["sec"] = time.time() - t0
    return out


def _run(args):
    eps, n_exp, seed, m, a, s = args
    fn = os.path.join(OUT, f"ov2_{eps:g}_{n_exp}_{m:g}_{a:g}_{s:g}_{seed}.json")
    if os.path.exists(fn):
        return json.load(open(fn))
    r = task(eps, n_exp, seed, m, a, s)
    json.dump(r, open(fn + ".tmp", "w"))
    os.replace(fn + ".tmp", fn)
    return r


def summary(eps, n_exp, m, a, s):
    import glob
    R = [json.load(open(f)) for f in glob.glob(os.path.join(OUT, f"ov2_{eps:g}_{n_exp}_{m:g}_{a:g}_{s:g}_*.json"))]
    C = continuum(m, a, s)
    ms = lambda x: (np.mean(x), np.std(x, ddof=1) / np.sqrt(len(x)) if len(x) > 1 else float("nan"))
    print(f"eps {eps}, N_exp {n_exp}, kink m {m}, bump a {a} s {s}; {len(R)} sprinklings, {np.mean([r['sec'] for r in R]):.0f} s each")
    print(f"  continuum: S_kin {C['kin']:+.4f}  S_pot {C['pot0']:+.4f}  dS_pot(+a) {C['dpot+1']:+.4f}  dS_pot(-a) {C['dpot-1']:+.4f}"
          f"  first order +-{C['first+1']:+.4f}")
    print("  flat (k = 0):")
    for fm in ("plain", "rowdiff", "dirichlet"):
        x = ms([r["0"][fm] for r in R])
        sh = ms([r["0"][f"shift_{fm}"] for r in R])
        print(f"    S_kin {fm:9}: {x[0]:+.4f} +- {x[1]:.4f}   T0 shift 2 pi: {sh[0]:+.3e} +- {sh[1]:.1e}")
    x = ms([r["0"]["pot"] for r in R])
    print(f"    S_pot          : {x[0]:+.4f} +- {x[1]:.4f}")
    print("  conformal bump (T1): differences against the flat thinning of the same master sprinkling")
    for k in ("1", "-1"):
        for fm in ("plain", "rowdiff", "dirichlet"):
            d = ms([r[k][fm] - r["0"][fm] for r in R])
            sh = ms([r[k][f"shift_{fm}"] for r in R])
            print(f"    {'+a' if k == '1' else '-a'} dS_kin {fm:9}: {d[0]:+.4f} +- {d[1]:.4f}   T0 shift: {sh[0]:+.3e} +- {sh[1]:.1e}")
        d = ms([r[k]["pot"] - r["0"]["pot"] for r in R])
        print(f"    {'+a' if k == '1' else '-a'} dS_pot          : {d[0]:+.4f} +- {d[1]:.4f}   (continuum {C['dpot+1' if k == '1' else 'dpot-1']:+.4f})")
    for fm in ("dirichlet",):
        d1 = ms([(r["1"][fm] + r["1"]["pot"] - r["-1"][fm] - r["-1"]["pot"]) / 2 for r in R])
        print(f"  first order (S(+a) - S(-a))/2, Dirichlet + pot: {d1[0]:+.4f} +- {d1[1]:.4f}  (continuum {C['first+1']:+.4f})")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--eps", type=float, default=0.12)
    ap.add_argument("--nexp", type=int, default=16640)
    ap.add_argument("--m", type=float, default=12.0)
    ap.add_argument("--a", type=float, default=0.2)
    ap.add_argument("--s", type=float, default=0.12)
    ap.add_argument("--seeds", default="0-3")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--summary", action="store_true")
    q = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    if not q.summary:
        lo, _, hi = q.seeds.partition("-")
        jobs = [(q.eps, q.nexp, s, q.m, q.a, q.s) for s in range(int(lo), int(hi or lo) + 1)]
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(max_workers=q.workers) as ex:
            for r in ex.map(_run, jobs):
                print(f"seed {r['seed']}: {r['sec']:.0f} s", flush=True)
    summary(q.eps, q.nexp, q.m, q.a, q.s)


if __name__ == "__main__":
    main()
