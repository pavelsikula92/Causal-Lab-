"""
F6.8i-dev, step D2a (development, decides nothing): a map of the cusp in p at fixed eps. docs/derivations/f68i_dev.md
section 6.1.

At fixed eps the symmetric cut locus theta_c = -p/2 exists for p < p_c(eps); at p_c the two mirror rays merge with the
middle ray (a cusp of the caustic) and for p > p_c a single ray remains, whose transverse stiffness K dips towards 0 near
the symmetric phase. The D3a p-scan (diag/out/f68i_d3a_20261125T002324Z) found the largest values of <s4hat A0^-2>_theta
there, unconverged or flagged. This file maps the neighbourhood of p_c:

- p_c by bisection on the number of rays at the symmetric phase (3 below, 1 above);
- for p = p_c -/+ Delta: the half-width w of the 3-ray region around theta* = -p/2 (bisection on the ray count), the
  window slope alpha (p < p_c), min K over the phases;
- the smooth part (1/2pi) int dtheta s4hat A0^-2 (first arrival), on one piece (theta*, theta* + 2pi) graded
  geometrically towards theta* from both sides (first subinterval --d0, ratio --ratio), Gauss-Legendre at 6 and 12
  nodes per subinterval (their difference is the quadrature error);
- the window term W3/(4 mu_3 8) of D1 (p < p_c) where the window model fits inside the 3-ray region at the gammas
  used: the model needs delta up to 1.25 c_max/alpha with c_max = gamma_max s_max; it is computed with the D1 gammas
  scaled down to fit (scale reported), and skipped ("window not uniform") if the scale would fall below 0.05;
- power-law fits of |smooth|, min K, w and alpha against Delta on each side.

    python -m diag.f68i_d2a --eps 12,16,20 --deltas 0.4,0.2,0.1,0.05,0.025 --workers 8
Output diag/out/f68i_d2a_<UTC time>.json and .txt.
"""
from __future__ import annotations

import argparse
import math
import os
import time
from multiprocessing import Pool

import numpy as np

from diag import f68f_1b as B
from diag import f68f_1b_cut as CUT
from diag import f68f_1c as SER
from diag import f68f_common as C
from diag import f68i_d1_checks as D1
from diag import f68i_d3a as D3

S_MAX_WIN = 7.0          # s up to which the window integrand matters (g0(s^2 A) ~ 0 beyond, A ~ 1 to 1.5)


def nrays(eps, p, th):
    vv, S = CUT.rays(eps, p, th % (2 * math.pi))
    return int(vv.size)


def pc_job(args):
    eps, lo, hi, it = args
    if nrays(eps, lo, -lo / 2) < 3 or nrays(eps, hi, -hi / 2) >= 3:
        return {"eps": eps, "pc": float("nan"), "bracket_ok": False}
    for _ in range(it):
        m = 0.5 * (lo + hi)
        if nrays(eps, m, -m / 2) >= 3:
            lo = m
        else:
            hi = m
    return {"eps": eps, "pc": 0.5 * (lo + hi), "bracket_ok": True, "err": 0.5 * (hi - lo)}


def geom_job(args):
    """3-ray half-width w around theta* (p < p_c) and the tie gap at theta*."""
    eps, p = args
    th = -p / 2
    vv, S = CUT.rays(eps, p, th % (2 * math.pi))
    o = np.argsort(S)
    out = {"eps": eps, "p": p, "nray": int(vv.size)}
    if vv.size >= 3:
        out["tie_gap"] = float(S[o[1]] - S[o[0]])
        out["third_gap"] = float(S[o[2]] - S[o[0]])
        lo, hi = 0.0, 0.6
        if nrays(eps, p, th + hi) >= 3:
            out["w"] = float("inf")
        else:
            for _ in range(30):
                m = 0.5 * (lo + hi)
                if nrays(eps, p, th + m) >= 3:
                    lo = m
                else:
                    hi = m
            out["w"] = 0.5 * (lo + hi)
    return out


def graded_sym(th0, d0, ratio):
    return D3.graded(th0, th0 + 2 * math.pi, d0=d0, r=ratio)


def window_scaled(pool, eps, p, thc, w, n_nodes):
    """W3/(4 mu_3 8) with the D1 model, the gammas scaled so that the model stays inside the 3-ray region."""
    if not (w > 3e-4):
        return {"theta_c": thc, "alpha": float("nan"), "scale": 0.0, "term": None, "note": "window not uniform (w < 3e-4)"}
    probe = pool.map(D1.node_job, [(eps, p, thc - 1e-4), (eps, p, thc + 1e-4)])
    alpha0 = 0.5 * (probe[0]["d"] + probe[1]["d"]) / 1e-4
    c_need = max(D1.GAMMAS) * S_MAX_WIN
    scale = min(1.0, 0.7 * w * alpha0 / (1.25 * c_need))
    if scale < 0.05:
        return {"theta_c": thc, "alpha": alpha0, "scale": scale, "term": None, "note": "window not uniform"}
    gam = [g * scale for g in D1.GAMMAS]
    cmax = max(gam) * S_MAX_WIN
    dmax = 1.25 * cmax / alpha0
    u = 0.5 * (1 - np.cos(math.pi * (np.arange(n_nodes) + 0.5) / n_nodes))
    deltas = list(dmax * u / u.max())
    res = pool.map(D1.node_job, [(eps, p, thc + sg * dl) for sg in (-1, 1) for dl in deltas], chunksize=1)
    sides = []
    for k in range(2):
        rows = []
        for i, dl in enumerate(deltas):
            r = dict(res[k * n_nodes + i])
            r["delta"] = dl
            rows.append(r)
        sides.append(D1.fit_side(rows))
    saved = D1.C_MAX
    D1.C_MAX = cmax * 1.0001
    try:
        W = np.array([D1.window_W(sides, g) for g in gam])
        N = np.array([D1.window_W(sides, g, absval=True) for g in gam[:3]])
    finally:
        D1.C_MAX = saved
    G = np.array(gam)
    coef = np.linalg.lstsq(np.stack([G ** m for m in (1, 2, 3, 4, 5, 6)], 1), W, rcond=None)[0]
    coef3 = np.linalg.lstsq(np.stack([G ** m for m in (3, 4, 5, 6)], 1), W, rcond=None)[0]
    N1 = float(np.median(N / G[:3]))
    return {"theta_c": thc, "alpha": float(0.5 * (sides[0]["alpha"] + sides[1]["alpha"])), "scale": scale,
            "W3": float(coef[2]), "W3_3to6": float(coef3[0]), "W1_rel": float(coef[0] / N1),
            "W2_rel": float(coef[1] / N1), "term": float(coef[2] / (4 * D1.MU3 * 8)),
            "res_a": max(M["res_a"] for M in sides), "note": ""}


def smooth_graded(pool, eps, p, th0, gls, d0, ratio):
    subs = graded_sym(th0, d0, ratio)
    jobs = [(eps, p, a, b, n, False) for n in gls for a, b in subs]
    res = pool.map(D3.piece_job, jobs, chunksize=1)
    out = {}
    k = 0
    for n in gls:
        rr = res[k:k + len(subs)]
        k += len(subs)
        out[n] = {"value": sum(r["sum"] for r in rr) / (2 * math.pi), "minK": min(r["minK"] for r in rr),
                  "fmin": min(r["fmin"] for r in rr), "fmax": max(r["fmax"] for r in rr)}
    return out, len(subs)


def powfit(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = (x > 0) & (np.abs(y) > 0) & np.isfinite(y)
    if ok.sum() < 2:
        return float("nan"), float("nan")
    b, a = np.polyfit(np.log(x[ok]), np.log(np.abs(y[ok])), 1)
    return float(b), float(math.exp(a))


NAME = "f68i_d2a"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--eps", default="12,16,20")
    ap.add_argument("--deltas", default="0.4,0.2,0.1,0.05,0.025")
    ap.add_argument("--gl", default="6,12")
    ap.add_argument("--d0", type=float, default=1e-4)
    ap.add_argument("--ratio", type=float, default=2.0)
    ap.add_argument("--nodes", type=int, default=9)
    ap.add_argument("--workers", type=int, default=2)
    a = ap.parse_args()
    t0 = time.time()
    epss = [float(x) for x in a.eps.split(",")]
    deltas = [float(x) for x in a.deltas.split(",")]
    gls = [int(x) for x in a.gl.split(",")]
    here = os.path.dirname(os.path.abspath(__file__))
    code = C.code_md5([os.path.abspath(__file__), B.__file__, CUT.__file__, SER.__file__,
                       os.path.join(here, "f68i_d1_checks.py"), os.path.join(here, "f68i_d1.py"),
                       os.path.join(here, "f68i_d3a.py")])
    lines = ["F6.8i-dev D2a: a map of the cusp in p at fixed eps (development, decides nothing; docs/derivations/"
             "f68i_dev.md section 6.1)", ""]
    data = []
    with Pool(a.workers) as pool:
        pcs = pool.map(pc_job, [(e, 0.5, 12.0, 30) for e in epss], chunksize=1)
        for e, pcr in zip(epss, pcs):
            n0 = len(lines)
            pc = pcr["pc"]
            lines.append(f"eps {e:g}: p_c = {pc:.6f} (bisection on the ray count at the symmetric phase)")
            if not pcr["bracket_ok"]:
                lines.append("  bracket failed (no 3-ray region at p = 0.5 or still 3 rays at p = 12): skipped")
                print("\n".join(lines[n0:]), flush=True)
                continue
            rows = []
            ps = [pc - d for d in deltas] + [pc + d for d in deltas]
            geo = pool.map(geom_job, [(e, p) for p in ps], chunksize=1)
            for p, g in zip(ps, geo):
                below = p < pc
                th0 = -p / 2
                sm, nsub = smooth_graded(pool, e, p, th0, gls, a.d0, a.ratio)
                win = None
                if below and g["nray"] >= 3 and g.get("tie_gap", 1.0) < 1e-6:
                    win = window_scaled(pool, e, p, D1.refine_tie(e, p, th0 % (2 * math.pi)), g["w"], a.nodes)
                v6, v12 = sm[gls[0]]["value"], sm[gls[-1]]["value"]
                row = {"p": p, "delta": abs(p - pc), "side": "below" if below else "above", "geom": g,
                       "smooth": {str(n): sm[n] for n in gls}, "nsub": nsub, "window": win,
                       "minK": sm[gls[-1]]["minK"], "qerr": abs(v12 - v6)}
                rows.append(row)
                wtxt = ""
                if win is not None:
                    if win["term"] is None:
                        wtxt = f"  window: {win['note']} (scale {win['scale']:.3f}, alpha {win['alpha']:.3f})"
                    else:
                        wtxt = (f"  window {win['term']:+.4e} (alpha {win['alpha']:.3f}, gamma scale {win['scale']:.3f}, "
                                f"W1/N1 {win['W1_rel']:+.0e} W2/N1 {win['W2_rel']:+.0e})")
                gtxt = (f"w {g['w']:.4e}  tie {g['tie_gap']:.0e}" if "w" in g else f"rays {g['nray']}")
                lines.append(f"  p {p:.5f} ({row['side']}, Delta {row['delta']:.3f}): {gtxt}  min K {row['minK']:.4f}  "
                             f"smooth {v12:+.6e} (6/12 differ {row['qerr']:.1e}, {nsub} subintervals)  "
                             f"f in [{sm[gls[-1]]['fmin']:+.2e}, {sm[gls[-1]]['fmax']:+.2e}]" + wtxt)
            for side in ("below", "above"):
                rr = [r for r in rows if r["side"] == side]
                dd = [r["delta"] for r in rr]
                bs, _ = powfit(dd, [r["smooth"][str(gls[-1])]["value"] for r in rr])
                bk, _ = powfit(dd, [r["minK"] for r in rr])
                txt = f"  fits on the {side} side (log-log over Delta): |smooth| ~ Delta^{bs:+.2f}, min K ~ Delta^{bk:+.2f}"
                if side == "below":
                    ww = [r["geom"].get("w", float("nan")) for r in rr]
                    bw, _ = powfit(dd, ww)
                    al = [r["window"]["alpha"] for r in rr if r["window"] is not None]
                    da = [r["delta"] for r in rr if r["window"] is not None]
                    ba, _ = powfit(da, al)
                    txt += f", w ~ Delta^{bw:+.2f}, alpha ~ Delta^{ba:+.2f}"
                lines.append(txt)
            data.append({"eps": e, "pc": pcr, "rows": rows})
            print("\n".join(lines[n0:]), flush=True)
    lines.append(f"total {time.time() - t0:.0f} s")
    stem = C.save(NAME, vars(a), code, data, lines, t0)
    print(lines[-1])
    print("saved", stem + ".json")


if __name__ == "__main__":
    main()
