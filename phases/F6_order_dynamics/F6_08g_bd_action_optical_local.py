"""
F6.8g -- The local optical stiffness K_c(k) of the BD action at k l_k = 0.02-0.12, separated from its age component
===================================================================================================================

DRAFT, not frozen (21 Nov 2026; card docs/F6_8g_card_draft.md rev. 0, approved by the author as drafted, decision 43;
frozen by the author only).

Why this test exists
--------------------
OV3-lin-a (docs/derivations/order_variation_OV3.md, sections 1, 5 and 8; decision 41, item 4: option A) needs the
local part of the optical (nu^2) stiffness K_c(k) = C(k)/2 at the k of the source form factors (F6.10b) and of the
conformal and mixed stiffness (F6.8a, F6.8b2). F6.8c2's C(0) = 2.002 at n = kT = 8 contains the age contribution;
it is not a separate measurement of the local K_c(k).

What is computed
----------------
The code of F6.8c, unchanged: optical_row(k, T) gives, on the same nodes, the phase-averaged second-order densities of
the optical sector K_hat_bb(k, T) and of the conformal sector K_hat_c(k, T) (eps rho = 1, l_k = 1, age T of the static
field). Rows: k = 0.02, 0.03, 0.045, 0.06, 0.08, 0.10, 0.12 at the ages T = 20, 28, 40, 57, 80, 113, 160 with
n = kT <= 12; rows at T = 10 and 14 are reported only (young-age regime: at k = 0.06, T = 10 the optical row is
dominated by the boundary term of the age cut, K_hat_bb = 1.600).

The reading (stated before the run)
-----------------------------------
For each sector s, one joint linear least-squares fit over the rows with T >= 20 and n <= 8:
    K_hat_s(k, T) = K_loc,s(k) - k^2 F_s(kT),     F_s(n) = a_1 n + a_2 n^2,   F_s(0) = 0,
K_loc,s(k) free per k; the age component F_s one function of n for all k (the scaling the conformal and mixed data
obey to < 1e-4, F6.8c2). Reported separately for every k: (1) the total K_hat_bb(k, T) of every row; (2) the age
component k^2 F_bb(kT) of every row; (3) K_loc,bb(k); (4) the sensitivity: the spread (max - min) of K_loc,bb(k) over
the variants n <= 4 / 8 / 12 x T >= 20 / 28 x J = 2 / 3, and the youngest row taken directly (n <= 1, where it exists).

Instrument (INSTRUMENT_FAILED if any fails)
-------------------------------------------
I1 the method on the known sector: the same separation of K_hat_c reproduces F6.8a's local part (its fit
   K0 - k^2 (c1 ln(1/k) + c0), K0 = 1.00005, c1 = 3.98, c0 = 0.91) within 1e-3 at every k <= 0.1;
I2 the geometry of F6.8c: <T2> against the numerical Fermat ray <= 2e-3; d2V/db2 against numerical roots <= 2e-3; the
   uniform null <= 1e-9;
I3 level 2 against level 1 on the rows (0.03, 28), (0.06, 57), (0.12, 40): |diff K_hat_bb| <= 2e-3;
I4 the k = 0 null at T = 20, 40, 160: a_4 |J_bb| / k_ref^2 <= 2e-3, k_ref the smallest k of the test.

Criterion
---------
Fit K_loc,bb(k) = K0_bb - k^2 (c1 ln(1/k) + c0) over k <= 0.1 (main reading):
  PASSED     |K0_bb - 1| <= 0.01 and the sensitivity spread <= 0.005 at every k <= 0.1;
  FAILED     |K0_bb - 1| > 0.03;
  UNDECIDED  otherwise (including a larger spread: the separation is then not stable enough for OV3-lin-a).
Reported, not judged: phi_2^(C) = C_loc - 2 - 4 phi_0 with C_loc = 2 K_loc,bb and phi_0 = (2/3)(K_loc,c - 1) (to be
checked by an independent phi_2 later; decision 43: F6.8h not developed now); the rows at T = 10 and 14.
"""

from __future__ import annotations

import itertools
import math
import time

import numpy as np

from core.stats.parallel import parallel_map
from core.testing import Result
from phases.F6_order_dynamics.F6_08a_bd_action_conformal import _A4
from phases.F6_order_dynamics.F6_08c_bd_action_optical import _dispatch

TEST = {
    "id": "F6.8g",
    "phase": "F6",
    "title": "Local optical stiffness K_c(k) of the BD action at k l_k = 0.02-0.12, separated from its age component "
             "(draft)",
    "theory_part": "derivation of the substrate rule from the order's action; the static propagator of the order; "
                   "input of OV3-lin-a (decision 41)",
    "question": "Separated from its age component by a joint fit at fixed scaling k^2 F(kT), is the local optical "
                "stiffness K_loc,bb(k) at k l_k = 0.02-0.12 stable against the choice of the fitted range, and does it "
                "tend to the EH value (C -> 2) as k -> 0?",
    "method": "optical_row of F6.8c (unchanged) on 7 k x the age ladder T = 20-160 (n = kT <= 12); joint linear fit "
              "K_hat = K_loc(k) - k^2 (a1 n + a2 n^2) per sector; sensitivity over the fitted range; the method "
              "validated on the conformal sector against F6.8a",
    "prediction": "K0_bb within 0.01 of 1 (F6.8c2: 1.00118 at n = 8, age included); K_loc,bb(k) close to the conformal "
                  "and mixed local parts at k <= 0.05, below them at larger k",
    "criterion": "PROPOSED: PASSED if |K0_bb - 1| <= 0.01 and the sensitivity spread <= 0.005 at every k <= 0.1; FAILED if "
                 "|K0_bb - 1| > 0.03; else UNDECIDED; INSTRUMENT_FAILED if I1 (conformal separation vs F6.8a within "
                 "1e-3 at k <= 0.1), I2 (geometry of F6.8c), I3 (level 2 within 2e-3) or I4 (k = 0 null <= 2e-3) fails",
    "null_model": "constant index (k = 0): the response vanishes",
    "level": "chain",
    "meters": [],
    "depends_on": ["F6.8a", "F6.8c2", "F6.8e3"],
    "params": {"ks": [0.02, 0.03, 0.045, 0.06, 0.08, 0.10, 0.12], "ladder": [20, 28, 40, 57, 80, 113, 160],
               "report_T": [10, 14], "n_max": 12.0, "fit_n": 8.0, "fit_T": 20.0, "J": 2,
               "variants_n": [4.0, 8.0, 12.0], "variants_T": [20.0, 28.0], "variants_J": [2, 3], "young_n": 1.0,
               "check": [[0.03, 28], [0.06, 57], [0.12, 40]], "null_ages": [20, 40, 160],
               "quad": [4, 16, 0.8, 8, 16], "bnd_nodes": 32, "k_fit_max": 0.1,
               "f68a": [1.00005, 3.98, 0.91], "tol_I1": 1e-3, "tol_t2": 2e-3, "tol_geom": 2e-3, "tol_unif": 1e-9,
               "tol_level": 2e-3, "tol_null": 2e-3, "tol_pass": 0.01, "tol_fail": 0.03, "tol_spread": 0.005},
    "parts": {"main": {}},
    "quick": {"ks": [0.06, 0.12], "ladder": [20, 40], "report_T": [], "check": [[0.12, 40]], "null_ages": [20],
              "variants_n": [8.0], "variants_T": [20.0], "variants_J": [1, 2], "J": 1},
    "estimate": "~2.5 CPU-h (about 60 rows of 70-170 s, three level-2 rows); ~15 min wall with 16 processes",
}


def _rows_grid(p):
    out = []
    for k in p["ks"]:
        for T in list(p["ladder"]) + list(p["report_T"]):
            if float(k) * float(T) <= float(p["n_max"]) + 1e-9:
                out.append((float(k), float(T)))
    return out


def separate(rows, key, n_max, T_min, J):
    """Joint least squares K = K_loc(k) - k^2 sum_j a_j n^j over rows with T >= T_min, n <= n_max.
    Returns ({k: K_loc}, [a_j], max |residual|)."""
    use = [r for r in rows if r["T"] >= T_min - 1e-9 and r["k"] * r["T"] <= n_max + 1e-9]
    ks = sorted({r["k"] for r in use})
    if not use or len(use) < len(ks) + J:
        return {}, [], float("nan")
    A = np.zeros((len(use), len(ks) + J))
    y = np.zeros(len(use))
    for i, r in enumerate(use):
        n = r["k"] * r["T"]
        A[i, ks.index(r["k"])] = 1.0
        for j in range(1, J + 1):
            A[i, len(ks) + j - 1] = -r["k"] ** 2 * n ** j
        y[i] = r[key]
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    res = float(np.max(np.abs(A @ c - y)))
    return {k: float(c[i]) for i, k in enumerate(ks)}, [float(a) for a in c[len(ks):]], res


def _kfit(ks, K):
    """K = K0 - k^2 (c1 ln(1/k) + c0)."""
    ks = np.asarray(ks, float)
    A = np.c_[np.ones_like(ks), -ks ** 2 * np.log(1 / ks), -ks ** 2]
    c, *_ = np.linalg.lstsq(A, np.asarray(K, float), rcond=None)
    return float(c[0]), float(c[1]), float(c[2]), float(np.max(np.abs(A @ c - K)))


def run(p, ctx):
    t0 = time.time()
    quad = tuple(p["quad"])
    bn = int(p["bnd_nodes"])
    grid = _rows_grid(p)
    jobs = [("row", (k, T, 1.0, quad, bn)) for k, T in grid]
    n1 = len(jobs)
    checks = [(float(c[0]), float(c[1])) for c in p["check"]]
    jobs += [("row", (k, T, 2.0, quad, 2 * bn)) for k, T in checks]
    jobs += [("row", (0.0, float(T), 1.0, quad, bn)) for T in p["null_ages"]]
    jobs += [("t2", ()), ("geom", ())]
    res = parallel_map(_dispatch, jobs, ctx.workers, ctx.progress,
                       lambda i, r_, d, n: f"[F6.8g] task {d}/{n} ({time.time() - t0:.0f}s)")
    rows = [dict(r_) for r_ in res[:n1]]
    lev2 = res[n1:n1 + len(checks)]
    nulls = res[n1 + len(checks):n1 + len(checks) + len(p["null_ages"])]
    t2c, geom = res[-2], res[-1]
    for r_ in rows:
        r_["n"] = r_["k"] * r_["T"]

    lines = ["3+1, eps rho = 1 (l_k = 1); rows of F6.8c optical_row; K_hat_bb (optical, EH 2k^2) and K_hat_c (conformal, "
             "EH 6k^2) on the same nodes; separation K_hat = K_loc(k) - k^2 F(kT), F(0) = 0"]
    # ---------------- instruments
    ok_t2 = all(c_["rel"] <= float(p["tol_t2"]) for c_ in t2c)
    ok_g = geom["rel_gam"] <= float(p["tol_geom"]) and geom["unif"] <= float(p["tol_unif"])
    lines.append("  I2 <T2> vs Fermat ray: " + "; ".join(f"{c_['rel']:.1e}" for c_ in t2c)
                 + f"; d2V/db2 rel {geom['rel_gam']:.1e}; uniform null {geom['unif']:.1e} -> "
                 + ("OK" if ok_t2 and ok_g else "FAILS"))
    ok3 = True
    for (k, T), r2 in zip(checks, lev2):
        r1 = next((r_ for r_ in rows if abs(r_["k"] - k) < 1e-12 and abs(r_["T"] - T) < 1e-9), None)
        if r1 is None:
            continue
        d = abs(r2["Kbb"] - r1["Kbb"])
        ok3 &= d <= float(p["tol_level"])
        lines.append(f"  I3 k {k:g} T {T:g}: K_hat_bb level 1 {r1['Kbb']:+.6f}, level 2 {r2['Kbb']:+.6f}, |diff| {d:.1e}"
                     f" -> {'OK' if d <= float(p['tol_level']) else 'FAILS'}")
    kref = min(float(k) for k in p["ks"])
    nu = [_A4 * abs(n_["J_bb"]) / kref ** 2 for n_ in nulls]
    ok4 = max(nu) <= float(p["tol_null"])
    lines.append("  I4 null k = 0: " + "; ".join(f"T {n_['T']:g}: {u:.1e}" for n_, u in zip(nulls, nu))
                 + f" (units of K_hat_bb at k = {kref:g}; tol {p['tol_null']:g}) -> {'OK' if ok4 else 'FAILS'}")
    # ---------------- separation, main reading
    fn, fT, J = float(p["fit_n"]), float(p["fit_T"]), int(p["J"])
    Lbb, abb, rbb = separate(rows, "Kbb", fn, fT, J)
    Lc, ac, rc = separate(rows, "Kc", fn, fT, J)
    K0a, c1a, c0a = [float(x) for x in p["f68a"]]
    kmax = float(p["k_fit_max"])
    dI1 = {k: Lc[k] - (K0a - k * k * (c1a * math.log(1 / k) + c0a)) for k in Lc}
    ok1 = all(abs(v) <= float(p["tol_I1"]) for k, v in dI1.items() if k <= kmax + 1e-12)
    lines.append(f"  I1 conformal separation (a = {', '.join(f'{a:+.4f}' for a in ac)}, max res {rc:.1e}) against "
                 "F6.8a's local part: " + "; ".join(f"k {k:g}: {v:+.1e}" for k, v in sorted(dI1.items()))
                 + f" (tol {p['tol_I1']:g} at k <= {kmax:g}) -> {'OK' if ok1 else 'FAILS'}")
    # ---------------- sensitivity
    variants = {}
    for vn, vT, vJ in itertools.product(p["variants_n"], p["variants_T"], p["variants_J"]):
        L, _, _ = separate(rows, "Kbb", float(vn), float(vT), int(vJ))
        if L:
            variants[f"n<={vn:g},T>={vT:g},J={vJ}"] = L
    young = {}
    for k in sorted({r_["k"] for r_ in rows}):
        cand = [r_ for r_ in rows if r_["k"] == k and r_["T"] >= fT - 1e-9 and r_["n"] <= float(p["young_n"]) + 1e-9]
        if cand:
            young[k] = min(cand, key=lambda r_: r_["n"])["Kbb"]
    variants["youngest row (n <= %g)" % float(p["young_n"])] = young
    spread = {}
    for k in Lbb:
        vals = [v[k] for v in variants.values() if k in v]
        spread[k] = (max(vals) - min(vals)) if len(vals) > 1 else 0.0
    ok_spread = all(spread[k] <= float(p["tol_spread"]) for k in spread if k <= kmax + 1e-12)
    # ---------------- table
    lines.append("  k       T      n       (1) K_hat_bb   (2) k^2 F_bb(n)   K_hat_c")
    for r_ in sorted(rows, key=lambda r_: (r_["k"], r_["T"])):
        n = r_["n"]
        Fb = sum(a * n ** (j + 1) for j, a in enumerate(abb))
        tag = "" if r_["T"] >= fT - 1e-9 and n <= fn + 1e-9 else ("  (reported)" if r_["T"] < fT else "  (n > fit range)")
        lines.append(f"  {r_['k']:<7g} {r_['T']:<6g} {n:<7.3g} {r_['Kbb']:+.6f}      {r_['k'] ** 2 * Fb:+.6f}        "
                     f"{r_['Kc']:+.6f}{tag}")
    lines.append(f"  age component, optical: F_bb(n) = " + " + ".join(f"{a:+.5f} n^{j + 1}" for j, a in enumerate(abb))
                 + f" (max res {rbb:.1e}); conformal: " + " + ".join(f"{a:+.5f} n^{j + 1}" for j, a in enumerate(ac)))
    lines.append("  k       (3) K_loc,bb   (4) spread   K_loc,c   phi_0     phi_2^(C) (reported)")
    phi2 = {}
    for k in sorted(Lbb):
        ph0 = (2.0 / 3.0) * (Lc[k] - 1.0) if k in Lc else float("nan")
        phi2[k] = 2 * Lbb[k] - 2 - 4 * ph0
        lines.append(f"  {k:<7g} {Lbb[k]:+.6f}      {spread[k]:.1e}      {Lc.get(k, float('nan')):+.6f}  {ph0:+.5f}  "
                     f"{phi2[k]:+.5f}")
    for name, v in variants.items():
        lines.append(f"    variant {name}: " + ", ".join(f"{k:g}: {x:.5f}" for k, x in sorted(v.items())))
    kk = [k for k in sorted(Lbb) if k <= kmax + 1e-12]
    if len(kk) >= 4:
        K0, c1, c0, mres = _kfit(kk, [Lbb[k] for k in kk])
        lines.append(f"  fit K_loc,bb = K0 - k^2 (c1 ln(1/k) + c0) over k <= {kmax:g}: K0_bb = {K0:.5f}, c1 {c1:.3f}, "
                     f"c0 {c0:.3f}, max res {mres:.1e}; max spread {max(spread[k] for k in kk):.1e} "
                     f"(tol {p['tol_spread']:g})")
    else:
        K0 = float("nan")
        lines.append("  fit K0_bb not evaluated (fewer than four wavenumbers: quick run)")
    inst = {"I1": ok1, "I2": ok_t2 and ok_g, "I3": ok3, "I4": ok4}
    if not all(inst.values()):
        status = "INSTRUMENT_FAILED"
    elif math.isnan(K0):
        status = "UNDECIDED"
    elif abs(K0 - 1) <= float(p["tol_pass"]) and ok_spread:
        status = "PASSED"
    elif abs(K0 - 1) > float(p["tol_fail"]):
        status = "FAILED"
    else:
        status = "UNDECIDED"
    lines.append(f"cost: {time.time() - t0:.0f} s")
    return Result(status, lines, {"rows": rows, "K_loc_bb": Lbb, "K_loc_c": Lc, "F_bb": abb, "F_c": ac,
                                  "spread": spread, "variants": variants, "phi2_C": phi2, "K0_bb": K0,
                                  "I1": dI1, "instrument": inst,
                                  "nulls": [{"T": n_["T"], "J_bb": n_["J_bb"]} for n_ in nulls]})
