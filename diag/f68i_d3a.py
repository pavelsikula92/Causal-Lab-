"""
F6.8i-dev, step D3a (development, decides nothing): <s4hat A0^-2>_theta at fixed (eps, p) with the cut loci treated
properly. docs/derivations/f68i_dev.md section 5.3.

    <s4hat A0^-2>_theta = (1/2pi) sum_pieces int_piece dtheta s4hat A0^-2 (first arrival)  +  sum_loci W3 / (4 mu_3 8)

- The cut loci theta_c in [0, 2pi) are those of diag/f68f_1b_cut (grid of 128 phases, bisection; cached there).
- Between two neighbouring cut loci the first arrival is one smooth branch, but its s4hat A0^-2 peaks sharply at the
  cut loci (width ~0.05 at eps 12): each piece is cut into subintervals graded geometrically towards both ends (first
  0.002 from theta_c, ratio 2.5), Gauss-Legendre on each at two node counts (--gl, default 6 and 12; the difference is
  the quadrature error). Without a cut locus the function is smooth and periodic: equal phases, 8 x the count.
- Each cut locus adds the window term of D1 (f68i_dev 5.2): W3 / (4 mu_3 8), W3 the gamma^3 coefficient of the window
  part of the age term, from the two-well model fitted around theta_c (diag/f68i_d1_checks: the same model, nodes, fit
  and W(gamma) quadrature; check 2 with the exact needle is not repeated here).
- Reported per piece: min K (the transverse stiffness along the first arrival over the piece's nodes). K -> 0 marks a
  conjugate point near the first arrival (a cusp of the cut locus), where the q-series of A diverges (D2, open): such a
  piece is flagged and its value is not trusted.
- For scale, the equal-phase mean on 64 phases (diag/f68f_1c.point, as in E1) is reported too.

    python -m diag.f68i_d3a --points 12:1,16:1,12:0,4:1 --workers 8
Output diag/out/f68i_d3a_<UTC time>.json and .txt.
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

KMIN_FLAG = 0.5          # min K below this along a piece: conjugate point near (series not trusted)
GRADE_D0, GRADE_R = 0.002, 2.5   # graded subintervals: first at 0.002 from a cut locus, ratio 2.5


# ------------------------------------------------------------------------------------------------ jobs (in workers)

def loci_job(args):
    eps, p = args
    return CUT.job(eps, p, B.N_TH)


def piece_job(args):
    """sum_i w_i s4hat A0^-2 over Gauss-Legendre nodes on (a, b), or equal phases if periodic."""
    eps, p, a, b, n, periodic = args
    if periodic:
        ths = a + (b - a) * np.arange(n) / n
        w = np.full(n, (b - a) / n)
    else:
        x, wg = np.polynomial.legendre.leggauss(n)
        ths = 0.5 * (a + b) + 0.5 * (b - a) * x
        w = 0.5 * (b - a) * wg
    r = SER.series_coefficients(eps, p, ths)
    A = r["A"]
    f = np.array([SER.s4hat(row) for row in A]) * A[:, 0] ** -2.0
    return {"a": a, "b": b, "n": n, "sum": float((w * f).sum()), "minK": r["minK"], "nray": r["nray"],
            "fmin": float(f.min()), "fmax": float(f.max())}


def equal_job(args):
    eps, p, n = args
    return SER.point(eps, p, n)["phi"]


def graded(a, b, d0=GRADE_D0, r=GRADE_R):
    """Subintervals of (a, b) graded geometrically towards both ends (the integrand peaks at the cut loci)."""
    half = 0.5 * (b - a)
    ds = [0.0]
    d = d0
    while d < 0.7 * half:
        ds.append(d)
        d *= r
    ds.append(half)
    left = [a + x for x in ds]
    right = [b - x for x in reversed(ds[:-1])]
    edges = left + right
    return list(zip(edges[:-1], edges[1:]))


# ------------------------------------------------------------------------------------------------ window term

def window_term(pool, eps, p, thc0, n_nodes):
    """W3 / (4 mu_3 8) at one cut locus, with the D1 model (f68i_d1_checks.run_point without check 2)."""
    thc = D1.refine_tie(eps, p, thc0)
    probe = pool.map(D1.node_job, [(eps, p, thc - 1e-3), (eps, p, thc + 1e-3)])
    alpha0 = 0.5 * (probe[0]["d"] + probe[1]["d"]) / 1e-3
    dmax = 1.25 * D1.C_MAX / alpha0
    u = 0.5 * (1 - np.cos(math.pi * (np.arange(n_nodes) + 0.5) / n_nodes))
    deltas = list(dmax * u / u.max())
    jobs = [(eps, p, thc + sg * dl) for sg in (-1, 1) for dl in deltas]
    res = pool.map(D1.node_job, jobs, chunksize=1)
    sides = []
    for k, sg in enumerate((-1, 1)):
        rows = []
        for i, dl in enumerate(deltas):
            r = dict(res[k * n_nodes + i])
            r["delta"] = dl
            rows.append(r)
        M = D1.fit_side(rows)
        M["min_d3"] = min(r["d3"] for r in rows)
        sides.append(M)
    cont = max(abs(float(sides[0]["P1"][j](0.0) - sides[1]["P2"][j](0.0))) for j in range(4))
    G = np.array(D1.GAMMAS)
    W = np.array([D1.window_W(sides, g) for g in D1.GAMMAS])
    V = np.stack([G ** m for m in (1, 2, 3, 4, 5, 6)], 1)
    coef = np.linalg.lstsq(V, W, rcond=None)[0]
    V3 = np.stack([G ** m for m in (3, 4, 5, 6)], 1)
    coef3 = np.linalg.lstsq(V3, W, rcond=None)[0]
    N1 = float(np.median(np.array([D1.window_W(sides, g, absval=True) for g in D1.GAMMAS[:3]]) / G[:3]))
    return {"theta_c": thc, "alpha": [sides[0]["alpha"], sides[1]["alpha"]], "W3": float(coef[2]),
            "W3_3to6": float(coef3[0]), "W1_rel": float(coef[0] / N1), "W2_rel": float(coef[1] / N1),
            "term": float(coef[2] / (4 * D1.MU3 * 8)), "continuity": cont,
            "res_a": max(M["res_a"] for M in sides), "min_gap_third": min(M["min_d3"] for M in sides)}


# ------------------------------------------------------------------------------------------------ one point

def run_point(pool, eps, p, loci, gls, n_nodes, lines):
    t0 = time.time()
    thcs = sorted(float(f[0]) % (2 * math.pi) for f in loci)
    if thcs:
        edges = thcs + [thcs[0] + 2 * math.pi]
        pieces = [(edges[i], edges[i + 1], False) for i in range(len(thcs))]
    else:
        pieces = [(0.0, 2 * math.pi, True)]
    jobs, tags = [], []
    for k, n in enumerate(gls):
        for i, (a, b, per) in enumerate(pieces):
            if per:
                jobs.append((eps, p, a, b, 8 * n, True))
                tags.append((k, i))
            else:
                for u, v in graded(a, b):
                    jobs.append((eps, p, u, v, n, False))
                    tags.append((k, i))
    res = pool.map(piece_job, jobs, chunksize=1)
    eq = pool.map(equal_job, [(eps, p, 64)])[0]
    smooth = {}
    per_piece = []
    for k, n in enumerate(gls):
        rr = []
        for i, (a, b, per) in enumerate(pieces):
            sub = [r for r, t in zip(res, tags) if t == (k, i)]
            rr.append({"a": a, "b": b, "n": n, "nsub": len(sub), "sum": sum(r["sum"] for r in sub),
                       "minK": min(r["minK"] for r in sub), "nray": max(r["nray"] for r in sub),
                       "fmin": min(r["fmin"] for r in sub), "fmax": max(r["fmax"] for r in sub)})
        smooth[n] = sum(r["sum"] for r in rr) / (2 * math.pi)
        per_piece.append(rr)
    wins = [window_term(pool, eps, p, thc, n_nodes) for thc in thcs]
    wsum = sum(w["term"] for w in wins)
    nhi = gls[-1]
    total = smooth[nhi] + wsum
    minK = min(r["minK"] for r in per_piece[-1])
    flag = minK < KMIN_FLAG
    lines.append(f"eps {eps:g} p {p:g}: {len(thcs)} cut loci" + (" at theta_c " + ", ".join(f"{t:.6f}" for t in thcs)
                                                                if thcs else "") + f"; rays up to {max(r['nray'] for r in per_piece[-1])}")
    for i, (a, b, per) in enumerate(pieces):
        vals = "  ".join(f"n {n}: {per_piece[k][i]['sum'] / (2 * math.pi):+.6e}" for k, n in enumerate(gls))
        r = per_piece[-1][i]
        label = " periodic" if per else f" in {r['nsub']} graded subintervals"
        lines.append(f"  piece ({a:.4f}, {b:.4f}){label}: {vals}  min K {r['minK']:.3f}  "
                     f"f in [{r['fmin']:+.2e}, {r['fmax']:+.2e}]" + ("  CONJUGATE POINT NEAR (D2)" if r['minK'] < KMIN_FLAG else ""))
    sm = "  ".join(f"{n}/piece {smooth[n]:+.8e}" for n in gls)
    lines.append(f"  smooth part: {sm}  (difference {abs(smooth[gls[-1]] - smooth[gls[0]]):.1e})")
    for w in wins:
        lines.append(f"  window at {w['theta_c']:.6f}: alpha {w['alpha'][0]:.4f}/{w['alpha'][1]:.4f}  W1/N1 {w['W1_rel']:+.0e} "
                     f"W2/N1 {w['W2_rel']:+.0e}  W3 {w['W3']:+.6e} (m 3-6: {w['W3_3to6']:+.6e})  term {w['term']:+.6e}  "
                     f"continuity {w['continuity']:.0e}  fit res A {w['res_a']:.0e}  third-ray gap {w['min_gap_third']:.3f}")
    lines.append(f"  TOTAL <s4hat A0^-2>_theta = smooth {smooth[nhi]:+.6e} + windows {wsum:+.6e} = {total:+.6e}"
                 f"   [64 equal phases, as E1: {eq:+.6e}]" + ("   (NOT TRUSTED: conjugate point, D2)" if flag else ""))
    lines.append(f"  ({time.time() - t0:.0f} s)")
    return {"eps": eps, "p": p, "loci": [list(map(float, f)) for f in loci], "pieces": per_piece,
            "smooth": {str(n): v for n, v in smooth.items()}, "windows": wins, "window_sum": wsum, "total": total,
            "equal64": eq, "minK": minK, "conjugate_flag": bool(flag)}


NAME = "f68i_d3a"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--points", default="12:1,16:1,12:0,4:1", help="eps:p,...")
    ap.add_argument("--gl", default="6,12", help="Gauss-Legendre nodes per graded subinterval (two counts)")
    ap.add_argument("--nodes", type=int, default=9, help="window-model nodes per side")
    ap.add_argument("--workers", type=int, default=2)
    a = ap.parse_args()
    t0 = time.time()
    gls = [int(x) for x in a.gl.split(",")]
    pts = [tuple(float(x) for x in it.split(":")) for it in a.points.split(",")]
    here = os.path.dirname(os.path.abspath(__file__))
    code = C.code_md5([os.path.abspath(__file__), B.__file__, CUT.__file__, SER.__file__,
                       os.path.join(here, "f68i_d1_checks.py"), os.path.join(here, "f68i_d1.py")])
    lines = ["F6.8i-dev D3a: <s4hat A0^-2>_theta at fixed (eps, p), piecewise between the cut loci plus the window terms "
             "(development, decides nothing; docs/derivations/f68i_dev.md section 5.3)", ""]
    data = []
    with Pool(a.workers) as pool:
        loci_all = pool.map(loci_job, pts, chunksize=1)
        for (eps, p), loci in zip(pts, loci_all):
            n0 = len(lines)
            data.append(run_point(pool, eps, p, loci, gls, a.nodes, lines))
            print("\n".join(lines[n0:]), flush=True)
    lines.append(f"total {time.time() - t0:.0f} s")
    stem = C.save(NAME, vars(a), code, data, lines, t0)
    print(lines[-1])
    print("saved", stem + ".json")


if __name__ == "__main__":
    main()
