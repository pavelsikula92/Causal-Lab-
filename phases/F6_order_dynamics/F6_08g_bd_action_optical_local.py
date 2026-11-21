"""
F6.8g -- The young-field stiffness of the order at fixed n = kT = 4 (optical, conformal, mixed), k l_k = 0.045-0.12
===================================================================================================================

FROZEN by the author, 21 Nov 2026 15:48 (decision 46; card docs/F6_8g_card_draft.md rev. 2, decisions 44 (option B)
and 46). Every judged point was exercised at full resolution before freezing: the rows k >= 0.045 at n = 4 in two
development runs of revision 1 (Claude's workspace and the author's machine), and the new I3 point (0.045, n 4) at
level 2 (K_hat_bb 0.961251 against 0.961281 at level 1, 3.0e-5).

Revision 2 (decision 46): the judged rows are k = 0.045 - 0.12. The development runs of revision 1 (Claude's workspace and
the author's machine, sha 6bfb3471, not verdicts) showed that the optical rows at large age are not numerically robust at
small k: the two platforms differ by 2.7e-4 at (0.03, T 133), 1.4e-4 at (0.045, T 133) and 5e-5 at (0.06, T 100), but
by <= 6e-6 at T <= 89; and the row (0.03, T 133) lies ~3.5e-3 below the n-dependence of every other k (its fit residual
1e-3 against 5e-5). k = 0.03 (T = 133 at n = 4) is therefore not used. On k >= 0.045 the two runs give K0_bb(n 4) =
0.99944 / 0.99943 (residual 4.7e-5).

Why this test exists
--------------------
OV3-lin-a (docs/derivations/order_variation_OV3.md, sections 1, 5, 8; decision 41) needs the conformal, mixed and
optical stiffness at the same k and the same age. The young field is the field of age T = n_y / k with n_y = 4
(decision 44). Revision 0 of this card (decision 43) tried to separate a local part from an age part; its quick run
and a probe showed that this is not defined at order k^2 (the conformal sector at small n differs from F6.8a's
asymptotic-intercept local part by ~1 k^2, growing with n) and that the k = 0 null carries a young-age term
J_bb ~ 13 / T^4 (a row is clean only for n T >~ 100).

What is computed
----------------
The code of F6.8c (optical_row: K_hat_bb and K_hat_c on the same nodes) and of F6.8b2 (mixed_row: K_hat_ab),
unchanged; eps rho = 1, l_k = 1. Judged rows: k = 0.045, 0.06, 0.08, 0.10, 0.12 at n = 4 (T = 4/k = 89 ... 33), optical
and mixed. Reported rows: optical at n = 3 and n = 6 (T <= 160; the row (0.045, n 6) has T = 133, above the age where
the platforms agree to 1e-5, reported only).

Reading (stated before the run): at fixed n, per sector, K_hat_s(k, n/k) = K0_s - k^2 (c1 ln(1/k) + c0) over the k.

Instrument (INSTRUMENT_FAILED if any fails)
-------------------------------------------
I1 the conformal intercept at n = 4: |K0_c - 1| <= 1e-3 (F6.8c2: 0.99976 at n = 6, 8, 10);
I2 the geometry of F6.8c (<T2> vs Fermat ray <= 2e-3; d2V/db2 vs numerical roots <= 2e-3; uniform null <= 1e-9) and of
   F6.8b2 (dV/db, d2V/da db <= 2e-3; uniform null <= 1e-9);
I3 optical rows (0.045, n 4) and (0.12, n 4) at level 2 vs level 1: |diff K_hat_bb| <= 2e-3;
I4 the k = 0 optical null at T = 4/k for each judged k: a_4 |J_bb| / k^2 <= 2e-3.

Criterion
---------
From the optical fit at n = 4: PASSED if |K0_bb - 1| <= 0.01; FAILED if |K0_bb - 1| > 0.03; UNDECIDED otherwise.
Reported, not judged: every row; the fits at n = 3, 4, 6 (optical, conformal) and n = 4 (mixed); the inputs of
OV3-lin-a at k = 0.06, 0.12 (K_hat_c, K_hat_ab, K_hat_bb at n = 4); phi_2^(C) = 2 K_hat_bb - 2 - 4 phi_0,
phi_0 = (2/3)(K_hat_c - 1) (to be checked by an independent phi_2 later, never a check here).
"""

from __future__ import annotations

import math
import time

import numpy as np

from core.stats.parallel import parallel_map
from core.testing import Result
from phases.F6_order_dynamics.F6_08a_bd_action_conformal import _A4
from phases.F6_order_dynamics.F6_08b2_bd_action_mixed_rev import check_geometry as mixed_geometry
from phases.F6_order_dynamics.F6_08b2_bd_action_mixed_rev import mixed_row
from phases.F6_order_dynamics.F6_08c_bd_action_optical import check_geometry, check_t2, optical_row

TEST = {
    "id": "F6.8g",
    "phase": "F6",
    "title": "Young-field stiffness of the order at fixed n = kT = 4 (optical, conformal, mixed), k l_k = 0.045-0.12",
    "theory_part": "the static propagator of the order; inputs of OV3-lin-a in the young-field reading (decisions 41, "
                   "44)",
    "question": "At the fixed age n = kT = 4, does the optical stiffness of the expected BD action tend to the EH value "
                "(C -> 2) at long wavelength, and what are the three stiffnesses at k l_k = 0.045-0.12?",
    "method": "optical_row of F6.8c and mixed_row of F6.8b2 (unchanged) at T = 4/k for six k; fixed-n fit K0 - k^2 "
              "(c1 ln(1/k) + c0) per sector; reported rows at n = 3 and 6",
    "prediction": "K0_bb within 0.01 of 1 (F6.8c2: 0.99971 / 1.00118 at n = 6 / 8)",
    "criterion": "PASSED if |K0_bb - 1| <= 0.01 at n = 4; FAILED if > 0.03; else UNDECIDED; INSTRUMENT_FAILED "
                 "if I1 (K0_c at n = 4 within 1e-3 of 1), I2 (geometry of F6.8c and F6.8b2), I3 (level 2 within "
                 "2e-3) or I4 (k = 0 null at the judged ages <= 2e-3 in K_hat_bb units) fails",
    "null_model": "constant index (k = 0): the response vanishes",
    "level": "chain",
    "meters": [],
    "depends_on": ["F6.8a", "F6.8b2", "F6.8c2", "F6.10b"],
    "params": {"ks": [0.045, 0.06, 0.08, 0.10, 0.12], "n_judge": 4.0, "n_report": [3.0, 6.0], "T_max": 160.0,
               "mixed": True, "check": [0.045, 0.12], "ov3_ks": [0.06, 0.12],
               "quad": [4, 16, 0.8, 8, 16], "bnd_nodes": 32,
               "tol_I1": 1e-3, "tol_t2": 2e-3, "tol_geom": 2e-3, "tol_unif": 1e-9, "tol_level": 2e-3,
               "tol_null": 2e-3, "tol_pass": 0.01, "tol_fail": 0.03},
    "parts": {"main": {}},
    "quick": {"ks": [0.06, 0.12], "n_report": [], "check": [0.12], "ov3_ks": [0.06, 0.12]},
    "estimate": "~2.5 CPU-h (14 optical rows of 70-170 s, 2 level-2 rows, 5 null rows, 5 mixed rows); ~1 h wall with "
                "16 processes (set by the longest rows)",
}


def _dispatch(kind, args):
    if kind == "orow":
        return optical_row(*args)
    if kind == "mrow":
        return mixed_row(*args)
    if kind == "t2":
        return check_t2()
    if kind == "geom":
        return check_geometry()
    if kind == "mgeom":
        return mixed_geometry()
    raise ValueError(kind)


def _kfit(ks, K):
    """K = K0 - k^2 (c1 ln(1/k) + c0); returns (K0, c1, c0, max |res|)."""
    ks = np.asarray(ks, float)
    A = np.c_[np.ones_like(ks), -ks ** 2 * np.log(1 / ks), -ks ** 2]
    c, *_ = np.linalg.lstsq(A, np.asarray(K, float), rcond=None)
    return float(c[0]), float(c[1]), float(c[2]), float(np.max(np.abs(A @ c - K)))


def run(p, ctx):
    t0 = time.time()
    quad = tuple(p["quad"])
    bn = int(p["bnd_nodes"])
    ks = [float(k) for k in p["ks"]]
    nj = float(p["n_judge"])
    Tmax = float(p["T_max"])
    orows = [(k, nj, nj / k) for k in ks]
    orows += [(k, float(n), float(n) / k) for n in p["n_report"] for k in ks if float(n) / k <= Tmax + 1e-9]
    jobs = [("orow", (k, T, 1.0, quad, bn)) for k, n, T in orows]
    n_o = len(jobs)
    mrows = [(k, nj, nj / k) for k in ks] if p["mixed"] else []
    jobs += [("mrow", (k, T, 1.0, quad, bn)) for k, n, T in mrows]
    checks = [float(k) for k in p["check"]]
    jobs += [("orow", (k, nj / k, 2.0, quad, 2 * bn)) for k in checks]
    jobs += [("orow", (0.0, nj / k, 1.0, quad, bn)) for k in ks]
    jobs += [("t2", ()), ("geom", ()), ("mgeom", ())]
    res = parallel_map(_dispatch, jobs, ctx.workers, ctx.progress,
                       lambda i, r_, d, n: f"[F6.8g] task {d}/{n} ({time.time() - t0:.0f}s)")
    i = 0
    ro = res[i:i + n_o]; i += n_o
    rm = res[i:i + len(mrows)]; i += len(mrows)
    lev2 = res[i:i + len(checks)]; i += len(checks)
    nulls = res[i:i + len(ks)]; i += len(ks)
    t2c, geom, mgeom = res[i], res[i + 1], res[i + 2]
    O = {(k, n): r_ for (k, n, T), r_ in zip(orows, ro)}
    M = {k: r_ for (k, n, T), r_ in zip(mrows, rm)}

    lines = [f"3+1, eps rho = 1 (l_k = 1); young-field stiffness at fixed n = kT (judged n = {nj:g}); K_hat = second "
             "derivative of the expected BD action / EH (6 k^2 conformal, 8 k^2 mixed, 2 k^2 optical)"]
    # ---------------- instruments
    ok_t2 = all(c_["rel"] <= float(p["tol_t2"]) for c_ in t2c)
    ok_g = geom["rel_gam"] <= float(p["tol_geom"]) and geom["unif"] <= float(p["tol_unif"])
    ok_mg = (mgeom["rel_beta"] <= float(p["tol_geom"]) and mgeom["rel_gam"] <= float(p["tol_geom"])
             and mgeom["unif"] <= float(p["tol_unif"]))
    lines.append("  I2 optical: <T2> vs Fermat " + "; ".join(f"{c_['rel']:.1e}" for c_ in t2c)
                 + f"; d2V/db2 rel {geom['rel_gam']:.1e}; null {geom['unif']:.1e}. Mixed: rel {mgeom['rel_beta']:.1e} / "
                 f"{mgeom['rel_gam']:.1e}; null {mgeom['unif']:.1e} -> {'OK' if ok_t2 and ok_g and ok_mg else 'FAILS'}")
    ok3 = True
    for k, r2 in zip(checks, lev2):
        r1 = O[(k, nj)]
        d = abs(r2["Kbb"] - r1["Kbb"])
        ok3 &= d <= float(p["tol_level"])
        lines.append(f"  I3 k {k:g} T {nj / k:.4g}: K_hat_bb level 1 {r1['Kbb']:+.6f}, level 2 {r2['Kbb']:+.6f}, "
                     f"|diff| {d:.1e} -> {'OK' if d <= float(p['tol_level']) else 'FAILS'}")
    nu = {k: _A4 * abs(n_["J_bb"]) / k ** 2 for k, n_ in zip(ks, nulls)}
    ok4 = all(v <= float(p["tol_null"]) for v in nu.values())
    lines.append("  I4 null k = 0 at T = 4/k: " + "; ".join(f"k {k:g} (T {nj / k:.4g}): {v:.1e}" for k, v in nu.items())
                 + f" (K_hat_bb units; tol {p['tol_null']:g}) -> {'OK' if ok4 else 'FAILS'}")
    # ---------------- fits
    fits = {}
    for n in [nj] + [float(x) for x in p["n_report"]]:
        kk = [k for k in ks if (k, n) in O]
        if len(kk) >= 4:
            fits[("bb", n)] = _kfit(kk, [O[(k, n)]["Kbb"] for k in kk])
            fits[("c", n)] = _kfit(kk, [O[(k, n)]["Kc"] for k in kk])
    if len(M) >= 4:
        fits[("ab", nj)] = _kfit(sorted(M), [M[k]["Kab"] for k in sorted(M)])
    if ("c", nj) in fits:
        K0c = fits[("c", nj)][0]
        ok1 = abs(K0c - 1) <= float(p["tol_I1"])
        lines.append(f"  I1 conformal intercept at n = {nj:g}: K0_c = {K0c:.5f} (tol {p['tol_I1']:g}) -> "
                     f"{'OK' if ok1 else 'FAILS'}")
    else:
        ok1 = True
        lines.append("  I1 not evaluated (fewer than four wavenumbers: quick run)")
    # ---------------- table
    lines.append("  k       n   T        K_hat_bb    K_hat_c     K_hat_ab    phi_2^(C) (reported)")
    for (k, n) in sorted(O, key=lambda kn: (kn[1], kn[0])):
        r_ = O[(k, n)]
        kab = M[k]["Kab"] if (n == nj and k in M) else float("nan")
        ph0 = (2.0 / 3.0) * (r_["Kc"] - 1.0)
        lines.append(f"  {k:<7g} {n:<3g} {r_['T']:<8.4g} {r_['Kbb']:+.6f}   {r_['Kc']:+.6f}   {kab:+.6f}   "
                     f"{2 * r_['Kbb'] - 2 - 4 * ph0:+.5f}")
    for (s, n), f in sorted(fits.items(), key=lambda x: (x[0][1], x[0][0])):
        lines.append(f"  fit {s:2s} at n = {n:g}: K0 = {f[0]:.5f}, c1 {f[1]:.3f}, c0 {f[2]:.3f}, max res {f[3]:.1e}")
    ov3 = {}
    for k in [float(x) for x in p["ov3_ks"]]:
        if (k, nj) in O:
            ov3[k] = {"K_c": O[(k, nj)]["Kc"], "K_ab": M[k]["Kab"] if k in M else float("nan"),
                      "K_bb": O[(k, nj)]["Kbb"]}
    lines.append("  inputs of OV3-lin-a at n = 4: " + "; ".join(
        f"k {k:g}: K {v['K_c']:.5f}, K_ab {v['K_ab']:.5f}, K_bb {v['K_bb']:.5f}" for k, v in ov3.items()))
    inst = {"I1": ok1, "I2": ok_t2 and ok_g and ok_mg, "I3": ok3, "I4": ok4}
    K0 = fits[("bb", nj)][0] if ("bb", nj) in fits else float("nan")
    if not all(inst.values()):
        status = "INSTRUMENT_FAILED"
    elif math.isnan(K0):
        status = "UNDECIDED"
    elif abs(K0 - 1) <= float(p["tol_pass"]):
        status = "PASSED"
    elif abs(K0 - 1) > float(p["tol_fail"]):
        status = "FAILED"
    else:
        status = "UNDECIDED"
    lines.append(f"  optical K0_bb at n = {nj:g}: {K0:.5f} (PASSED <= {p['tol_pass']:g}, FAILED > {p['tol_fail']:g})")
    lines.append(f"cost: {time.time() - t0:.0f} s")
    return Result(status, lines, {
        "optical": [dict(r_, n=n) for (k, n, T), r_ in zip(orows, ro)],
        "mixed": rm, "nulls": [{"k": k, "T": n_["T"], "J_bb": n_["J_bb"], "units": nu[k]} for k, n_ in zip(ks, nulls)],
        "fits": {f"{s}@{n:g}": f for (s, n), f in fits.items()}, "ov3": ov3, "instrument": inst, "K0_bb": K0})
