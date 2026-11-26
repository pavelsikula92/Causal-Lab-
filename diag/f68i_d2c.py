"""
F6.8i-dev, step D2c (development, decides nothing): the cusp constants along eps. docs/derivations/f68i_dev.md
section 6.4.

Revision 2 (26 Nov). The ray count of f68f_1b_cut.rays uses a v grid of B.N_V = 400 points; near the cusp the three
rays fall into one grid cell and the count drops to 1 before p_c, so revision 1 placed p_c too low (by 1.6e-4 at eps 16,
4.0e-4 at eps 20; f68i_dev section 6.5). Revision 2 sets B.N_V = NV_FINE in every job, measures the geometry at
Delta in {4D, 2D, D, D/2}, and fits sig = sig' (Delta + d) + b (Delta + d)^2 on all four with a free shift d
of the origin (d is the remaining error of p_c; it should be ~0), alongside the two-point extrapolation of revision 1.
The cut-locus count (128 phases) is optional (--loci); revision 1 recorded it.

Section 6.2 found that each cusp of the first arrival adds C Psi (2 pi)^-1 gamma^(11/7) to the age term per unit p
(C = -8.83e-3, the normal form V = a^4/4 + sig a^2/2 + tau a), with Psi = (sig' tau')^-1 N'^(-9/7):
sig' = d sig / dp, tau' = d tau / d theta and A = N' w(c) / c^2 at the cusp. D2a gave Psi at eps 12, 16, 20. This file
maps Psi(eps) from just above pi^2 (where the cut locus is born at p = 0) to just below 4 pi^2 (the next caustic),
and counts the cusps.

For each eps:
- ray counts on both symmetric lines theta = -p/2 and -p/2 + pi on a p list; p_c by bisection where the count on a
  line changes (the change and the counts on both sides are reported; the normal form below is applied to 3 <-> 1
  changes where the 3-ray side is below);
- the cut loci over 128 phases (f68f_1b_cut) at p = p_c / 2 and the number of them off the symmetric lines;
- at p = p_c - Delta, Delta in {D, D/2}, D = min(0.05, p_c / 20): sig = -2 sqrt(third-ray gap); the tie-gap slope
  alpha at theta_c from d(delta) at delta = w/10 and w/20 (d odd in delta: alpha = (8 d(delta/2) - d(delta)) / (3 delta));
  tau' = alpha / (2 sqrt|sig|); the half-width w (bisection) against 2 (|sig|/3)^(3/2) / tau';
- sig' and tau' extrapolated linearly to Delta -> 0 from the two Delta;
- N' from the exact needle (f68f_1b_exact) at p_c, theta_c: N' = A c^(1/4) / (4^(1/4) B), B = int (1 - x^4)^(3/2) dx,
  at q = 0.0625, 0.044, 0.031 (the drift and the local slopes of ln A against ln q, -1/2 for a cusp, are reported);
- Psi = (sig' tau')^-1 N'^(-9/7) (N' at q = 0.0625).

    python -m diag.f68i_d2c --workers 8
Output diag/out/f68i_d2c_<UTC time>.json and .txt.
"""
from __future__ import annotations

import argparse
import math
import os
import time
from multiprocessing import Pool

import numpy as np
from scipy.integrate import quad
from scipy.optimize import least_squares

from diag import f68f_1b as B
from diag import f68f_1b_cut as CUT
from diag import f68f_1b_exact as EX
from diag import f68f_common as C
from diag.f68i_d2a import nrays

EPS_DEFAULT = "9.9,10,10.25,10.5,11,12,14,16,20,25,30,35,39"
P_SCAN = [0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0]
QN = [0.0625, 0.044, 0.031]
NV_FINE = 6400
LOCI = False
BQ = quad(lambda x: (1 - x ** 4) ** 1.5, -1, 1, epsabs=1e-14, epsrel=1e-14)[0]
NAME = "f68i_d2c"


def lowest_two(eps, p, th):
    vv, S = CUT.rays(eps, p, th % (2 * math.pi))
    S = np.sort(S)
    return S, vv.size


def bisect_change(eps, line, lo, hi, n_lo, it=30):
    for _ in range(it):
        m = 0.5 * (lo + hi)
        if nrays(eps, m, -m / 2 + line) == n_lo:
            lo = m
        else:
            hi = m
    return 0.5 * (lo + hi), 0.5 * (hi - lo)


def half_width(eps, p, th, hi=0.6, it=24):
    lo = 0.0
    if nrays(eps, p, th + hi) >= 3:
        return float("inf")
    for _ in range(it):
        m = 0.5 * (lo + hi)
        if nrays(eps, p, th + m) >= 3:
            lo = m
        else:
            hi = m
    return 0.5 * (lo + hi)


def geometry(eps, p):
    th = -p / 2
    S, n = lowest_two(eps, p, th)
    out = {"p": p, "nray": int(n)}
    if n < 3:
        return out
    sig = -2 * math.sqrt(max(S[2] - S[0], 0.0))
    w = half_width(eps, p, th)
    d = {}
    for f in (10.0, 20.0):
        dl = w / f
        Sp, _ = lowest_two(eps, p, th + dl)
        d[f] = float(Sp[1] - Sp[0])
    dl = w / 10.0
    alpha = (8 * d[20.0] - d[10.0]) / (3 * dl)
    taup = alpha / (2 * math.sqrt(-sig))
    out.update(tie_gap=float(S[1] - S[0]), third_gap=float(S[2] - S[0]), sig=sig, w=w, alpha=alpha,
               alpha_lin=d[20.0] / (w / 20.0), tau_p=taup, w_pred=2 * (-sig / 3) ** 1.5 / taup)
    return out


def needle(eps, p):
    th = (-p / 2) % (2 * math.pi)
    roots = B.branches(eps, p, np.array([th]), 1000)[0]
    vv = B.newton(roots, eps, p, np.full(roots.size, th), 4000)
    _, _, S = B.end_state(vv, eps, p, np.full(vv.size, th), 4000)
    Smin = float(S.min())
    near = vv[S - Smin < max(QN) ** 2 / 2]
    paths = B.integrate(near, eps, p, np.full(near.size, th), EX.N_RK)[:, 0, :]
    xs, _ = np.polynomial.legendre.leggauss(EX.N_S)
    sn = 0.5 * (xs + 1.0)
    sg = np.linspace(0, 1, EX.N_RK + 1)
    psi_rays = [np.interp(sn, sg, paths[:, j]) for j in range(paths.shape[1])]
    A, worst = EX.phase_A(eps, p, th, Smin, psi_rays, qs=QN)
    q = np.array(QN)
    Av = np.array([A[x] for x in QN])
    c = q * q / 2
    Np = Av * c ** 0.25 / (4 ** 0.25 * BQ)
    return {"A": Av.tolist(), "N": Np.tolist(), "slopes": (np.diff(np.log(Av)) / np.diff(np.log(q))).tolist(),
            "worst": float(worst), "nray": int(vv.size)}


def sig_fit(rows):
    """sig = sig' (Delta + d) + b (Delta + d)^2, least squares on the four Delta; returns sig', d, b, max residual."""
    g = sorted([r for r in rows if "sig" in r], key=lambda r: r["delta"])
    if len(g) < 3:
        return None
    D = np.array([r["delta"] for r in g])
    s = -np.array([r["sig"] for r in g])
    f = lambda a: (a[0] * (D + a[1]) + a[2] * (D + a[1]) ** 2) / s - 1
    a = least_squares(f, [s[0] / D[0], 0.0, 0.0], xtol=1e-15, ftol=1e-15, gtol=1e-15).x
    return {"sig_p": float(a[0]), "d": float(a[1]), "b": float(a[2]), "resid": float(np.abs(f(a)).max())}


def eps_job(eps):
    B.N_V = NV_FINE
    t0 = time.time()
    out = {"eps": eps, "scan": {}, "changes": [], "cusps": []}
    for line, lname in ((0.0, "0"), (math.pi, "pi")):
        cnt = [nrays(eps, p, -p / 2 + line) for p in P_SCAN]
        out["scan"][lname] = cnt
        for i in range(len(P_SCAN) - 1):
            if cnt[i] != cnt[i + 1]:
                pc, err = bisect_change(eps, line, P_SCAN[i], P_SCAN[i + 1], cnt[i])
                out["changes"].append({"line": lname, "p_c": pc, "err": err, "n_below": cnt[i], "n_above": cnt[i + 1]})
    for ch in out["changes"]:
        if not (ch["line"] == "0" and ch["n_below"] == 3 and ch["n_above"] == 1):
            continue
        pc = ch["p_c"]
        D = min(0.05, pc / 20)
        rows = [geometry(eps, pc - dl) | {"delta": dl} for dl in (4 * D, 2 * D, D, D / 2)]
        cu = {"p_c": pc, "rows": rows}
        if LOCI:
            loci = CUT.job(eps, pc / 2, B.N_TH)
            off = [l for l in loci if abs((l[0] + pc / 4 + math.pi) % (2 * math.pi) - math.pi) > 1e-3]
            cu |= {"loci_at_half": len(loci), "loci_off_symmetric": len(off)}
        cu["fit"] = sig_fit(rows)
        if all("sig" in r for r in rows[2:]):
            s1, s2 = rows[2]["sig"] / D, rows[3]["sig"] / (D / 2)
            t1, t2 = rows[2]["tau_p"], rows[3]["tau_p"]
            cu["sig_p"] = -(2 * s2 - s1)
            cu["tau_p"] = 2 * t2 - t1
            nd = needle(eps, pc)
            cu["needle"] = nd
            Np = nd["N"][0]
            cu["Psi"] = 1.0 / (cu["sig_p"] * cu["tau_p"]) * Np ** (-9 / 7)
        out["cusps"].append(cu)
    out["seconds"] = time.time() - t0
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--eps", type=str, default=EPS_DEFAULT)
    ap.add_argument("--workers", type=int, default=0)
    ap.add_argument("--loci", action="store_true")
    a = ap.parse_args()
    global LOCI
    LOCI = a.loci
    t0 = time.time()
    code = C.code_md5([os.path.abspath(__file__), CUT.__file__, B.__file__, EX.__file__,
                       os.path.join(os.path.dirname(os.path.abspath(__file__)), "f68i_d2a.py")])
    eps_list = [float(x) for x in a.eps.split(",")]
    lines = ["F6.8i-dev D2c: the cusp constants along eps (development, decides nothing; "
             "docs/derivations/f68i_dev.md section 6.4)", "",
             f"revision 2: ray counts with B.N_V = {NV_FINE}; p scan: {P_SCAN}; "
             f"normal form at p_c - Delta, Delta in (4D, 2D, D, D/2), D = min(0.05, p_c/20); "
             f"N' at q = {QN}", ""]
    res = []
    with Pool(a.workers or os.cpu_count()) as pool:
        for r in pool.imap(eps_job, eps_list):
            res.append(r)
            L = [f"eps {r['eps']:g} ({r['seconds']:.0f} s): ray counts on theta = -p/2: {r['scan']['0']}; "
                 f"on -p/2 + pi: {r['scan']['pi']}"]
            for ch in r["changes"]:
                L.append(f"  change on line {ch['line']}: {ch['n_below']} -> {ch['n_above']} at p = {ch['p_c']:.6f}")
            for cu in r["cusps"]:
                L.append(f"  cusp p_c {cu['p_c']:.7f}" + (f": cut loci at p_c/2: {cu['loci_at_half']} "
                         f"({cu['loci_off_symmetric']} off the symmetric line)" if "loci_at_half" in cu else ""))
                if cu["fit"]:
                    fz = cu["fit"]
                    L.append(f"    fit on the 4 Delta: sig' {fz['sig_p']:.5f}  shift of p_c {fz['d']:+.2e}  "
                             f"b {fz['b']:+.4f}  residual {fz['resid']:.1e}")
                for g in cu["rows"]:
                    if "sig" in g:
                        L.append(f"    Delta {g['delta']:.4f}: sig {g['sig']:+.5f} (sig/Delta {g['sig'] / g['delta']:+.4f})  "
                                 f"alpha {g['alpha']:.4f} (linear {g['alpha_lin']:.4f})  tau' {g['tau_p']:.4f}  "
                                 f"w {g['w']:.4e} pred {g['w_pred']:.4e} (ratio {g['w'] / g['w_pred']:.5f})")
                    else:
                        L.append(f"    Delta {g['delta']:.4f}: {g['nray']} rays (no 3-ray geometry)")
                if "Psi" in cu:
                    nd = cu["needle"]
                    L.append(f"    sig' {cu['sig_p']:.4f}  tau' {cu['tau_p']:.4f}  N' at q {QN}: "
                             f"{', '.join(f'{x:.4f}' for x in nd['N'])}  slopes {', '.join(f'{x:+.3f}' for x in nd['slopes'])}  "
                             f"=> Psi {cu['Psi']:.5f}")
            lines += L
            print("\n".join(L), flush=True)
    res.sort(key=lambda r: r["eps"])
    lines += ["", "Psi(eps), cusps with a 3 -> 1 change on theta = -p/2 (sig' two-point; fit in brackets):"]
    for r in res:
        for cu in r["cusps"]:
            if "Psi" in cu:
                fs = f" [{cu['fit']['sig_p']:.4f}, shift {cu['fit']['d']:+.1e}]" if cu.get("fit") else ""
                lines.append(f"  eps {r['eps']:6.2f}  p_c {cu['p_c']:.6f}  sig' {cu['sig_p']:.4f}{fs}  tau' {cu['tau_p']:.4f}  "
                             f"N' {cu['needle']['N'][0]:.4f}  Psi {cu['Psi']:.5f}")
    lines.append(f"total {time.time() - t0:.0f} s")
    stem = C.save(NAME, vars(a), code, res, lines, t0)
    print("\n".join(lines[-len(res) - 2:]))
    print("saved", stem + ".json")


if __name__ == "__main__":
    main()
