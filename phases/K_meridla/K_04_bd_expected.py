"""
K.4 -- The expected smeared BD operator (quadrature)
===================================================

The expected Benincasa-Dowker operator on a static field
(core.operators.bd_expect) is the reference against which the field
equation A5' is tested in one sprinkling (F1.2) and from which the
asymptotic shape is derived (F1.3). This test calibrates it.

Part "quadrature":
  (a) 1+1: against an independent light-cone integral (null coordinates,
      a separate code path), Gaussian test field;
  (b) local limit: E[B f] / lap f -> 1 as l_k / w -> 0 in 1+1, 2+1, 3+1;
  (c) the hat-basis rows reproduce the direct quadrature of the interpolant.
Parts "shape_d3", "shape_d4": the nonlocal field equation E[B] u = q 1(r < a)
  solved in units of l_k (one computation per dimension) at two radial
  resolutions must agree (convergence of quadrature and basis). The shape g
  is only REPORTED here; the physics question belongs to F1.3.
"""

from __future__ import annotations

import math
import time

import numpy as np

from core.operators import bd_expect as BE
from core.operators.continuum import ball_source_solution
from core.operators.solve import fit_profile
from core.stats.stats import fmt
from core.testing import Result

TEST = {
    "id": "K.4",
    "phase": "K",
    "title": "Očekávaný BD operátor (kvadratura)",
    "theory_part": "měřidlo pro A5′ (F1.2, F1.3); bez fyzikálního závěru",
    "question": "Počítá kvadratura očekávanou hodnotu rozmazaného BD operátoru správně, včetně lokální limity a konvergence?",
    "method": "(a) nezávislý integrál přes světelný kužel v 1+1; (b) poměr k Laplaceovu operátoru pro l_k/w -> 0; "
              "(c) řádky na bázi proti přímé kvadratuře; tvar nelokální rovnice ve dvou rozlišeních",
    "prediction": "shoda (a) v 0,005; (b) poměr 1 v 0,015; (c) shoda na 1e-9; dvě rozlišení shodná",
    "criterion": "quadrature: |ref - quad| <= tol_lightcone, |poměr - 1| <= tol_local ve všech dimenzích, "
                 "relativní rozdíl řádků <= 1e-9. shape: |g(res1) - g(res2)| <= conv_g a |dev(res1) - dev(res2)| <= conv_dev.",
    "null_model": "—",
    "level": "calibration",
    "meters": [],
    "depends_on": [],
    "old_id": "E45b",
    "params": {"tol_lightcone": 0.005, "tol_local": 0.015, "conv_g": 0.05, "conv_dev": 0.02,
               "shape_R": 40.0, "shape_a": 1.0, "shape_lo": 5.0, "shape_hi_frac": 0.5,
               "res1": [0.25, 1.06, 3], "res2": [0.2, 1.04, 5]},
    "parts": {"quadrature": {"mode": "quadrature"}, "shape_d3": {"mode": "shape", "D": 3},
              "shape_d4": {"mode": "shape", "D": 4}},
    "quick": {"shape_R": 24.0, "res1": [0.5, 1.15, 2], "res2": [0.4, 1.12, 3]},
    "estimate": "quadrature ~1 min; shape ~5-15 min na dimenzi",
}


def _lightcone(eps, rho, tx, w, n=2000):
    """1+1, x on the axis, null coordinates a = dt - dx, b = dt + dx (dV = da db / 2), difference form."""
    al = np.linspace(math.log(1e-12), math.log(4.0), n)
    A = np.exp(al)
    dal = al[1] - al[0]
    a, b = A[:, None], A[None, :]
    dt, dx = (a + b) / 2, (b - a) / 2
    ok = (tx - dt >= 0) & (np.abs(dx) <= 0.5)
    uu = eps * rho * a * b / 2
    K = (-2 + 4 * uu - uu * uu) * np.exp(-uu) * a * b / 2 * dal * dal * ok
    one, diff = K.sum(), (K * (np.exp(-dx ** 2 / (2 * w * w)) - 1)).sum()
    return -2 * eps * rho * (1 + rho * eps * (one + diff)) / (-1 / w ** 2)


def _quadrature(p):
    w = 0.2
    gauss = lambda r: np.exp(-np.asarray(r) ** 2 / (2 * w * w))
    ref = float(np.mean([_lightcone(0.05, 2000, tx, w) for tx in (0.85, 0.95)]))
    quad = float(np.mean([BE.expected_apply(2, 2000, 0.05, tx, 0.0, 0.5, gauss, 1.0) / (-1 / w ** 2) for tx in (0.85, 0.95)]))
    lim = {}
    for D, er in ((2, 1e7), (3, 1e9), (4, 1e9)):
        lim[D] = BE.expected_apply(D, er / 0.05, 0.05, 3.0, 0.0, 2.0, gauss, 1.0) / (-(D - 1) / w ** 2)
    nodes = np.linspace(0.0, 0.5, 51)
    f = np.exp(-nodes ** 2 / (2 * 0.15 ** 2))
    rows = BE.expected_radial_rows(3, 4000, 0.1, 0.9, [0.2], 0.5, nodes)[0]
    direct = BE.expected_apply(3, 4000, 0.1, 0.9, 0.2, 0.5, lambda r: np.interp(r, nodes, f), float(np.interp(0.2, nodes, f)))
    rel = abs(rows @ f - direct) / max(1.0, abs(direct))
    ok_a = abs(ref - quad) <= p["tol_lightcone"]
    ok_b = all(abs(v - 1.0) <= p["tol_local"] for v in lim.values())
    ok_c = rel <= 1e-9
    lines = [f"  (a) 1+1 light cone {ref:.4f} vs quadrature {quad:.4f} (tol {p['tol_lightcone']}) -> {'OK' if ok_a else 'FAILS'}",
             "  (b) local limit E[B f] / lap f: " + ", ".join(f"{D-1}+1 {v:.4f}" for D, v in lim.items())
             + f" (tol {p['tol_local']}) -> {'OK' if ok_b else 'FAILS'}",
             f"  (c) hat rows {rows @ f:.8g} vs direct {direct:.8g} (rel {rel:.1e}) -> {'OK' if ok_c else 'FAILS'}"]
    return ok_a and ok_b and ok_c, lines, {"lightcone": ref, "quadrature": quad, "local_limit": lim, "rows_rel": rel}


def _shape(p, ctx):
    D = int(p["D"])
    R, a = float(p["shape_R"]), float(p["shape_a"])
    window = [float(p["shape_lo"]), float(p["shape_hi_frac"]) * R]
    out = {}
    for label in ("res1", "res2"):
        dr, ratio, per = p[label]
        t0 = time.time()
        prof = BE.expected_green_profile(D, R, a, dr, ratio, int(per))
        nodes, u = prof["nodes"], prof["u"]
        fit = fit_profile(nodes, u, window)
        uc = ball_source_solution(D - 1, 1.0, a, R)(nodes)
        wn = (nodes >= window[0]) & (nodes <= window[1])
        dev = float(np.sqrt(((u[wn] - uc[wn]) ** 2).sum() / (uc[wn] ** 2).sum()))
        out[label] = {"g": fit["g"], "dev": dev, "cond": prof["cond"], "nodes": nodes, "u": u}
        ctx.progress(f"[K.4] {D-1}+1 shape {label}: g = {fit['g']:+.2f}, deviation {dev:.3f} ({time.time()-t0:.0f}s)")
    dg = abs(out["res1"]["g"] - out["res2"]["g"])
    dd = abs(out["res1"]["dev"] - out["res2"]["dev"])
    ok = dg <= p["conv_g"] and dd <= p["conv_dev"]
    lines = [f"  {D-1}+1 nonlocal equation in units of l_k: R = {R:g}, source radius {a:g}, window [{window[0]:g}, {window[1]:g}]",
             f"     res1 g = {out['res1']['g']:+.3f} dev {out['res1']['dev']:.3f}; res2 g = {out['res2']['g']:+.3f} dev {out['res2']['dev']:.3f}",
             f"     convergence |dg| = {dg:.3f} (tol {p['conv_g']}), |d dev| = {dd:.3f} (tol {p['conv_dev']}) -> {'OK' if ok else 'FAILS'}",
             "     (g is reported only; whether the profile is 1/r is the question of F1.3)"]
    return ok, lines, out


def run(p, ctx):
    if p["mode"] == "quadrature":
        ok, lines, data = _quadrature(p)
    else:
        ok, lines, data = _shape(p, ctx)
    return Result("PASSED" if ok else "FAILED", lines, data)
