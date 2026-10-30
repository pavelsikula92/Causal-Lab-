"""
F6.8c2 -- Optical (nu^2) static sector of the BD action, read at fixed k T (revision of F6.8c)
==============================================================================================

Why this revision exists (29 Oct 2026)
--------------------------------------
F6.8c (full run 29 Oct) ended FAILED and keeps that verdict and its criteria: its reading fitted, per k,
K_hat_bb(k, T) = K_loc,bb(k) - s(k) T over three ages and extrapolated to T -> 0. For the optical sector this
model is wrong: K_hat_bb falls faster than linearly with T (k = 0.14 and 0.2: K_hat_bb ~ K_loc - q k^4 T^2 with
q ~ 0.019, no linear term), so the extrapolated K_loc,bb and the fit K0_bb = 1.559 (rms 4e-2) carry no
information. The rows with the largest ages (k = 0.05 at T = 240, 320; k = 0.07 at T = 229) also break the
scaling that the other rows obey, which points to a numerical loss at T >~ 200 in the bending terms (not
seen by the k = 0 null, which has no bending); this is not settled.

The reading used here does not assume any form of the age dependence. In units l_k = 1 the age-dependent
part of every sector scales as k^2 F(k T) (conformal: F = c_T kT; the conformal and mixed data obey this to
<1e-4), so at FIXED n = k T

    K_hat(k, T = n/k) = K0 - k^2 [c1 ln(1/k) + c0(n)],

and the intercept K0 is the long-wavelength local stiffness for every n. Applied to the F6.8c and F6.8b2 data
at n = 8 (a description of data already seen, NOT used for this criterion): K0 = 0.99993 (conformal; known
1.00005), 0.99988 (mixed), 0.99914 (optical; residual 1.4e-4). Because those data were seen, this revision is
judged on FRESH wavenumbers only.

What is computed
----------------
As F6.8c (same code: phases/F6_order_dynamics/F6_08c_bd_action_optical.py, optical_row), on new k with ages
T = n/k for n = 6, 8, 10 (T <= 167, below the ages where F6.8c lost the scaling), plus reported rows that
probe the growth of the optical age term at larger n with moderate T.

Criteria (stated before the full run)
-------------------------------------
Rows: k = 0.06, 0.085, 0.12, 0.17 at n = k T = 6, 8, 10.
Instrument (INSTRUMENT_FAILED if any fails):
  I1 geometry as F6.8c (a) <T2> vs numerical Fermat ray <= 2e-3, (b) d2V/db2 vs numerical roots + finite
     differences <= 2e-3, (c) uniform null <= 1e-9;
  I2 the conformal K_hat on this grid against F6.8a's own row (k = 0.2, T = 40) <= 1e-3;
  I3 null k = 0 at T = 40, 100, 170 in K_hat_bb units at the smallest k whose rows reach down to that age <= 2e-3;
  I4 level 2 vs level 1 on the check rows <= 2e-3;
  I5 the reading method on the conformal sector of the same rows: fixed-n intercept at n = 8,
     |K0_conf - 1| <= 2e-3 (F6.8a established K0 = 1.00005).
Status (fit at n = 8 over the four new k, K_hat_bb = K0_bb - k^2 (c1 ln(1/k) + c0)):
  PASSED     |K0_bb - 1| <= 0.01: the third coefficient of the EH form emerges, C(k -> 0) = 2;
  FAILED     |K0_bb - 1| > 0.03;
  UNDECIDED  otherwise.
Reported, not judged: K0_bb at n = 6 and 10 (the intercept should not depend on n); the same fits for the
conformal sector; per k the shape of the age dependence F_bb(n) = (K0_bb - K_hat_bb)/k^2 - [c1 ln(1/k) + ...];
the rows at k = 0.12, n = 16 and 24 (T = 133, 200); gamma(0) = 2 K0_ab / K0_bb - 1 with K0_ab = 0.99988 (F6.8b2 data
read at n = 8, descriptive).
"""

from __future__ import annotations

import math
import time

import numpy as np

from core.stats.parallel import parallel_map
from core.testing import Result
from phases.F6_order_dynamics.F6_08a_bd_action_conformal import _A4
from phases.F6_order_dynamics.F6_08c_bd_action_optical import _dispatch

TEST = {
    "id": "F6.8c2",
    "phase": "F6",
    "title": "Optical (nu^2) static sector of the BD action read at fixed kT, revision of F6.8c: does C -> 2 emerge?",
    "theory_part": "derivation of the substrate rule from the order's action; third stage of F6.8; input of F8.1",
    "question": "Read without assuming the form of its age dependence (fixed n = kT, fresh wavenumbers), does the nu^2 "
                "static coefficient C of the expected discrete BD action tend to the Einstein-Hilbert value 2 at long "
                "wavelength?",
    "method": "as F6.8c (same row computation); k = 0.06, 0.085, 0.12, 0.17 at T = n/k, n = 6, 8, 10; at fixed n the fit "
              "K_hat = K0 - k^2 (c1 ln(1/k) + c0) gives the long-wavelength intercept K0 for any age dependence of the "
              "form k^2 F(kT); the method is checked on the conformal sector of the same rows (I5)",
    "prediction": "EH: K0_bb = 1 (C -> 2), the same intercept at every n",
    "criterion": "INSTRUMENT: I1a/b/c, I2, I3, I4 as F6.8c (null ages 40, 100, 170); I5 conformal fixed-n intercept at n = 8 "
                 "|K0_conf - 1| <= 2e-3. Fit at n = 8 over the four k: PASSED |K0_bb - 1| <= 0.01; FAILED > 0.03; "
                 "UNDECIDED otherwise.",
    "null_model": "constant index (k = 0): flat space rescaled, the response vanishes (finite-age residual ~ T^-4)",
    "level": "chain",
    "meters": [],
    "depends_on": ["F6.8a", "F6.8b2", "F6.8c"],
    "params": {"ks": [0.06, 0.085, 0.12, 0.17], "n_kT": [6.0, 8.0, 10.0], "judge_n": 8.0,
               "extra": [[0.12, 16.0], [0.12, 24.0]], "check": [[0.06, 8.0], [0.17, 8.0]],
               "null_ages": [40.0, 100.0, 170.0], "quad": [4, 16, 0.8, 8, 16], "bnd_nodes": 32,
               "tol_t2": 2e-3, "tol_geom": 2e-3, "tol_unif": 1e-9, "tol_conf": 1e-3, "tol_null": 2e-3,
               "tol_level": 2e-3, "tol_method": 2e-3, "tol_pass": 0.01, "tol_fail": 0.03,
               "K0_ab_descriptive": 0.99988},
    "parts": {"main": {}},
    "quick": {"ks": [0.12, 0.17], "n_kT": [8.0], "extra": [], "check": [[0.17, 8.0]], "null_ages": [40.0]},
    "estimate": "like F6.8c or less (12 rows with T <= 167, 2 extra rows, 2 level-2 rows, 3 nulls): ~1-2 h wall with 10 "
                "processes",
}


def _fit(ks, K):
    """K = K0 - k^2 (c1 ln(1/k) + c0); returns K0, c1, c0, max residual (0 if exactly determined)."""
    k = np.asarray(ks, float)
    A = np.c_[np.ones_like(k), -k ** 2 * np.log(1.0 / k), -k ** 2]
    c, *_ = np.linalg.lstsq(A, np.asarray(K, float), rcond=None)
    res = float(np.max(np.abs(np.asarray(K) - A @ c))) if len(k) > 3 else 0.0
    return float(c[0]), float(c[1]), float(c[2]), res


def run(p, ctx):
    t0 = time.time()
    ks = [float(k) for k in p["ks"]]
    ns = [float(n) for n in p["n_kT"]]
    quad = tuple(p["quad"])
    bn = int(p["bnd_nodes"])
    grid = [(k, n, n / k) for k in ks for n in ns]
    jobs = [("row", (k, T, 1.0, quad, bn)) for k, n, T in grid]
    n1 = len(jobs)
    extra = [(float(e[0]), float(e[1])) for e in p["extra"]]
    jobs += [("row", (k, n / k, 1.0, quad, bn)) for k, n in extra]
    checks = [(float(c[0]), float(c[1])) for c in p["check"]]
    jobs += [("row", (k, n / k, 2.0, quad, 2 * bn)) for k, n in checks]
    jobs += [("row", (0.0, float(T), 1.0, quad, bn)) for T in p["null_ages"]]
    jobs += [("row", (0.2, 40.0, 1.0, quad, bn))]                     # the conformal row compared with F6.8a (I2)
    jobs += [("t2", ()), ("geom", ()), ("conf", ())]
    res = parallel_map(_dispatch, jobs, ctx.workers, ctx.progress,
                       lambda i, r_, d, n: f"[F6.8c2] task {d}/{n} ({time.time() - t0:.0f}s)")
    rows = res[:n1]
    i = n1
    rx = res[i:i + len(extra)]; i += len(extra)
    lev2 = res[i:i + len(checks)]; i += len(checks)
    nulls = res[i:i + len(p["null_ages"])]
    rc = res[-4]
    t2c, geom, conf = res[-3], res[-2], res[-1]
    for r_, (k, n, T) in zip(rows, grid):
        r_["n"] = n

    lines = ["3+1, eps rho = 1 (l_k = 1); optical sector as F6.8c, read at fixed n = k T; K_hat_bb = d2E[S]/db2 / EH 2k^2"]
    ok1a = all(c_["rel"] <= float(p["tol_t2"]) for c_ in t2c)
    lines.append("  I1a <T2> vs numerical Fermat ray: " + "; ".join(f"rel {c_['rel']:.1e}" for c_ in t2c)
                 + f" (tol {p['tol_t2']:g})")
    ok1b = geom["rel_gam"] <= float(p["tol_geom"]) and geom["unif"] <= float(p["tol_unif"])
    lines.append(f"  I1b d2V/db2 {geom['gam']:+.7f} vs {geom['gam_sph']:+.7f} (rel {geom['rel_gam']:.1e}); uniform null "
                 f"{geom['unif']:.1e} (tol {p['tol_geom']:g}, {p['tol_unif']:g})")
    dconf = abs(rc["Kc"] - conf["Kc_8a"])
    ok2 = dconf <= float(p["tol_conf"])
    lines.append(f"  I2 conformal K_hat {rc['Kc']:+.6f} vs F6.8a row {conf['Kc_8a']:+.6f} (k = {conf['k']:g}, "
                 f"T = {conf['T']:g}): |diff| {dconf:.1e} (tol {p['tol_conf']:g})")

    def k_ref(Tn):
        cand = [k for k in ks if min(n / k for n in ns) <= Tn * (1 + 1e-9)]
        return min(cand) if cand else max(ks)
    nunits = [(n_["T"], k_ref(n_["T"]), abs(n_["J_bb"]) * _A4 / k_ref(n_["T"]) ** 2) for n_ in nulls]
    ok3 = max(u_[2] for u_ in nunits) <= float(p["tol_null"])
    lines.append("  I3 null k = 0: " + "; ".join(f"T = {n_['T']:g}: J_bb {n_['J_bb']:+.1e}, {u_[2]:.1e} at k = {u_[1]:g}"
                                                  for n_, u_ in zip(nulls, nunits)) + f" (tol {p['tol_null']:g})")
    ok4 = True
    for r2 in lev2:
        r1 = next((r_ for r_ in rows if abs(r_["k"] - r2["k"]) < 1e-12 and abs(r_["T"] - r2["T"]) < 1e-9), None)
        if r1 is None:
            continue
        d = abs(r2["Kbb"] - r1["Kbb"])
        ok4 &= d <= float(p["tol_level"])
        lines.append(f"  I4 k = {r2['k']:g}, T = {r2['T']:g}: K_hat_bb level 1 {r1['Kbb']:+.6f}, level 2 {r2['Kbb']:+.6f}, "
                     f"|diff| {d:.1e} (tol {p['tol_level']:g})")

    by_n = {n: sorted([r_ for r_ in rows if r_["n"] == n], key=lambda r_: r_["k"]) for n in ns}
    fits = {}
    for n in ns:
        kk = [r_["k"] for r_ in by_n[n]]
        fits[n] = {"bb": _fit(kk, [r_["Kbb"] for r_ in by_n[n]]), "conf": _fit(kk, [r_["Kc"] for r_ in by_n[n]])}
    jn = float(p["judge_n"])
    if len(ks) >= 4 and jn in fits:
        k0c = fits[jn]["conf"][0]
        ok5 = abs(k0c - 1.0) <= float(p["tol_method"])
        lines.append(f"  I5 reading method on the conformal sector, n = {jn:g}: K0_conf = {k0c:.5f} (F6.8a 1.00005; "
                     f"tol {p['tol_method']:g})")
    else:
        ok5 = True
        lines.append("  I5 not evaluated (fewer than four wavenumbers: quick run)")
    inst_ok = ok1a and ok1b and ok2 and ok3 and ok4 and ok5

    lines.append("  k       n    T          K_hat_bb   K_hat (conformal)")
    for r_ in sorted(rows, key=lambda r_: (r_["k"], r_["n"])):
        lines.append(f"  {r_['k']:<7g} {r_['n']:<4g} {r_['T']:<10.4g} {r_['Kbb']:+.6f}  {r_['Kc']:+.6f}")
    for n in ns:
        fb, fc = fits[n]["bb"], fits[n]["conf"]
        lines.append(f"  fit at n = {n:g}: optical K0_bb = {fb[0]:.5f} (c1 {fb[1]:.3f}, c0 {fb[2]:.3f}, max res {fb[3]:.1e}); "
                     f"conformal K0 = {fc[0]:.5f} (max res {fc[3]:.1e})")
    for r_ in rx:
        lines.append(f"  reported row k = {r_['k']:g}, T = {r_['T']:g} (n = {r_['k'] * r_['T']:g}): K_hat_bb {r_['Kbb']:+.6f}, "
                     f"K_hat {r_['Kc']:+.6f}")
    # shape of the optical age dependence per k: (K0_bb(judge) - K_hat_bb) / k^2 - [c1 ln(1/k) + c0(judge)] + c0(n) ...
    if jn in fits:
        k0b = fits[jn]["bb"][0]
        g0 = 2.0 * float(p["K0_ab_descriptive"]) / k0b - 1.0
        lines.append(f"  reported: gamma(0) = 2 K0_ab / K0_bb - 1 = 2 x {p['K0_ab_descriptive']:g} / {k0b:.5f} - 1 = {g0:+.5f}")
        lines.append("  reported: age dependence per k, D(n) = (K0_bb - K_hat_bb(k, n/k)) / k^2: " + "; ".join(
            f"k = {k:g}: " + ", ".join(f"{(k0b - r_['Kbb']) / k ** 2:.3f}" for r_ in sorted(
                [r_ for r_ in rows if r_["k"] == k], key=lambda r_: r_["n"])) for k in ks))
    if not inst_ok:
        status = "INSTRUMENT_FAILED"
    elif jn not in fits or len(ks) < 4:
        status = "UNDECIDED"
    elif abs(fits[jn]["bb"][0] - 1.0) <= float(p["tol_pass"]):
        status = "PASSED"
    elif abs(fits[jn]["bb"][0] - 1.0) > float(p["tol_fail"]):
        status = "FAILED"
    else:
        status = "UNDECIDED"
    data = {"rows": rows, "extra": rx, "level2": lev2, "nulls": nulls, "t2": t2c, "geom": geom, "conf": conf,
            "fits": {str(n): v for n, v in fits.items()}}
    return Result(status, lines, data)
