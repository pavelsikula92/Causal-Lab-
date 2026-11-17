"""
F6.8f -- The static response K_eff of the optical sector at finite amplitude (DRAFT, NOT FROZEN)
===============================================================================================

Why this test exists (the author's direction of 16 Nov 2026, 09:59, and the order of 22:49, point 3)
---------------------------------------------------------------------------------------------------
F6.8e3 (PASSED) found the optical age term of the second variation growing as (kT)^5. The second variation holds only
for eps = b (kT)^2 << 1; static fields of physical interest are old (eps huge). The author's question of 09:59: does the
expected BD action give a usable static response at finite b and large T? Quantities: E_opt(k, T, b) = dE[S_BD]/db and
K_eff = (E_opt(b) - E_opt(0)) / b, with the projection of F6.8e3.

What the development found (docs/derivations/f68f_dev.md, steps 1c, 2; decides nothing)
-------------------------------------------------------------------------------------
  * In the scaling limit (kappa = kT -> infinity, b -> 0 at fixed eps, p) the age density is dE/dn = c Phi(b n^2),
    Phi(eps) = int dp <s4hat A_0^-2>_theta from the needle volumes of K.18; c C2 is the b^2 coefficient of F6.8e3
    (checked: the p-profile and the constant agree with F6.8e3 to 5.5e-4 at n = 1536).
  * Hence, at fixed k and age n = kT, with M(eps) = int_0^eps sqrt(u) Phi'(u) du:
        E_opt = (c / 2 b^{3/2}) M(eps),  K_eff = (c / 2 b^{5/2}) M(eps),  K_2 = lim_{b->0} K_eff = (2/5) c C2 n^5,
        R(eps) = K_eff / K_2 = M(eps) / ((4/5) C2 eps^{5/2}).
    The K_eff curves of all (b, k, T) collapse onto one function R(eps); R = 1 is the second-order extrapolation from
    F6.8e3 (the null model). At fixed b, K_eff(T) = K_2(T) R(b k^2 T^2), with the local exponent
        gamma_T = d ln|K_eff| / d ln T = 2 d ln|M| / d ln eps       (5 for the second order; 0 for a static limit).
  * M = sqrt(eps) Phi(eps) - int_0^{sqrt eps} Phi(t^2) dt, Gauss-Legendre on panels in t (diag/f68f_keff.py).

What this card can and cannot decide
------------------------------------
The instrument (K.18) is calibrated for eps <= 4 (one ray per phase; the first caustic is at eps = pi^2 for p = 0).
The card reads the K_eff curve on that range: whether it keeps the second-order growth, bends towards saturation, or
changes sign. Whether an old field (eps -> infinity) has a finite static response is a statement about Phi at large
eps; it needs the multi-ray range (an extension of K.18 with the cut loci, and a card F6.8g). Reading A of this card
is therefore an indication on eps <= 4, not a proof of the static limit.

Question, quantities, criteria (PROPOSALS for the author's decision; written 16 Nov 23:40 UTC, before any K_eff row of
the development was seen; nothing is frozen)
-------------------------------------------------------------------------------------------------------------------
  Judged eps = 0.25, 1, 2, 4 (the ends of the t panels [0, 1/2], [1/2, 1], [1, sqrt 2], [sqrt 2, 2], G = 4 Gauss
  nodes each); C2 from eps = 0.1, 0.2 by Richardson. Development used eps 0.1, 0.2, 0.5, 3 and the nodes of panels
  ending at 0.5 and 3 at reduced resolution (16 phases, 4 directions); the judged rows at this resolution are not seen.
  gamma_T on the judged intervals 0.25 -> 1, 1 -> 2, 2 -> 4: gamma = 2 ln(M_b / M_a) / ln(eps_b / eps_a).
  INSTRUMENT_FAILED: K.18 not PASSED (separate record); more than one ray at some phase; stationarity of the ray
    above 1e-8; |R(G) - R(G - 1)| above 1e-3 at a judged eps; Phi with 128 against 64 phases above 1e-3 (relative)
    at a judged eps.
  Reading A (the K_eff curve bends towards saturation on eps <= 4): M keeps the sign of C2 at every judged eps, and
    gamma_T(2 -> 4) <= 2.5 (at most half the second-order exponent), and gamma_T decreases over the three intervals.
  Reading B (K_eff grows or changes sign): M changes sign at a judged eps (K_eff changes sign), or
    gamma_T(2 -> 4) >= 5 (no slowdown against the second order).
  Reading C (only young fields: the curve departs from the second order without a bend that reaches A): otherwise.
  Status: PASSED = reading A, FAILED = reading B, UNDECIDED = reading C.
  Recorded with every status: R(eps) and N(eps) = Phi / (C2 eps^2) at the judged eps (the null model R = N = 1), and
  |R - 1| against 0.05 (whether the second order still holds at eps = 4).
  Checks of the author's list (09:59): 1 b -> 0 (R, N -> 1 as eps -> 0; F6.8e3 in K.18 part e3); 2 zero age
  (M ~ eps^{5/2} -> 0); 3 k = 0 (Phi = 0, computed at eps = 0); 4 both signs of b (exact: n(-b, theta) =
  n(b, theta + pi), an even number of phases maps the grid onto itself, so K_eff(-b) = K_eff(b)); 5 conformal and
  mixed sectors (no code path here, unchanged); 6 the stable instrument (K.18, and the G and phase checks above).
Level: mechanism (the age term of the expected BD action in a declared model: a static plane-wave optical index).
"""

from __future__ import annotations

import math
import time

import numpy as np

from core.stats.parallel import parallel_map
from core.testing import Result

TEST = {
    "id": "F6.8f",
    "phase": "F6",
    "title": "Static response K_eff of the optical sector of the expected BD action at finite amplitude eps = b (kT)^2 "
             "(draft, not frozen)",
    "theory_part": "the static response of the order for old fields (B6-full); input of the static limit and of C.1",
    "question": "Does K_eff = (dE[S_BD]/db (b) - dE[S_BD]/db (0)) / b of the optical sector keep the second-order growth "
                "(kT)^5, bend towards saturation, or change sign on eps = b (kT)^2 <= 4?",
    "method": "scaling limit of the age derivative: needle volumes as a q^2 series to q^6 (K.18), the BD moments select "
              "the q^6 survivor Phi(eps); R(eps) = K_eff / K_2 = M(eps) / ((4/5) C2 eps^{5/2}), M = int sqrt(u) Phi'(u) du "
              "by parts on Gauss panels in t = sqrt(u); p on the panels of F6.8e3, 128 phases",
    "prediction": "none fixed (development saw eps 0.1, 0.2, 0.5, 3 and reduced-resolution nodes up to 3); the null model "
                  "is R = 1",
    "criterion": "INSTRUMENT: K.18 PASSED; one ray per phase; stationarity <= 1e-8; |R(G) - R(G-1)| <= 1e-3; phases "
                 "64 vs 128 <= 1e-3. PASSED (reading A): sign kept, gamma_T(2->4) <= 2.5 and decreasing; FAILED "
                 "(reading B): sign change or gamma_T(2->4) >= 5; UNDECIDED (reading C) otherwise (proposal, not frozen)",
    "null_model": "the second-order extrapolation from F6.8e3: R(eps) = 1, gamma_T = 5",
    "level": "mechanism",
    "meters": ["K.18"],
    "depends_on": ["F6.8e3", "K.18"],
    "params": {"ends": [0.25, 1.0, 2.0, 4.0], "g": 4, "c2_eps": [0.1, 0.2],
               "p_panels": [[0.0, 12.0, 10], [12.0, 24.0, 6]], "n_theta": 128,
               "gamma_A": 2.5, "gamma_B": 5.0, "tol_lin": 1e-8, "tol_g": 1e-3, "tol_theta": 1e-3, "tol_null": 0.05},
    "parts": {"main": {}},
    "quick": {"ends": [0.5, 3.0], "p_panels": [[0.0, 12.0, 2], [12.0, 24.0, 1]], "n_theta": 8, "g": 3},  # development eps only
    "estimate": "~9 h with 2 processes (22 + 12 eps x 16 directions x 128 phases of the series, ~115 s each; the "
                "phase check at 64 phases); ~2-3 h with 8 processes; float64 only (no 80-bit long double), not yet run on Windows",
}


def _point(eps, p, n_th):
    from diag import f68f_1c as S
    r = S.point(eps, p, n_th)
    return {"phi": r["phi"], "nray": r["nray"], "lin": r["lin"]}


def run(p, ctx):
    t0 = time.time()
    from diag import f68f_keff as KF
    from phases.F6_order_dynamics.F6_08e3_bd_action_optical_asymptotic_rev2 import _p_nodes
    pn, wp = _p_nodes(p["p_panels"])
    nt, g = int(p["n_theta"]), int(p["g"])
    ends = [float(e) for e in p["ends"]]
    e1, e2 = (float(v) for v in p["c2_eps"])
    pan, pan_c = KF.t_nodes(ends, g), KF.t_nodes(ends, g - 1)
    nodes = sorted(set([e1, e2] + ends + [float(x * x) for t, _ in pan for x in t]))
    nodes_c = sorted(set(float(x * x) for t, _ in pan_c for x in t) - set(nodes))
    jobs = [(e, float(x), nt) for e in nodes + nodes_c for x in pn]
    jobs += [(e, float(x), nt // 2) for e in ends for x in pn]
    jobs += [(0.0, float(x), nt) for x in pn]
    res = parallel_map(_point, jobs, ctx.workers, ctx.progress, lambda i, r, d, n: f"[F6.8f] {d}/{n}")
    tab = {j: r for j, r in zip(jobs, res)}

    def Phi(e, n=nt):
        return float((wp * np.array([tab[(e, float(x), n)]["phi"] for x in pn])).sum())

    C2 = (4 * Phi(e1) / e1 ** 2 - Phi(e2) / e2 ** 2) / 3 if abs(e2 - 2 * e1) < 1e-12 else Phi(e1) / e1 ** 2
    M = KF.assemble(ends, pan, Phi)
    Mc = KF.assemble(ends, pan_c, Phi)
    R = [m / (0.8 * C2 * e ** 2.5) for m, e in zip(M, ends)]
    Rc = [m / (0.8 * C2 * e ** 2.5) for m, e in zip(Mc, ends)]
    N = [Phi(e) / (C2 * e * e) for e in ends]
    dth = [abs(Phi(e, nt // 2) / Phi(e) - 1) for e in ends]
    nray = max(r["nray"] for r in res)
    lin = max(r["lin"] for r in res)
    k0 = max(abs(tab[(0.0, float(x), nt)]["phi"]) for x in pn)
    lines = [f"F6.8f (draft): K_eff / K_2 = R(eps); C2 = {C2:+.10e} (from eps {e1:g}, {e2:g}); {len(pn)} directions, "
             f"{nt} phases, G = {g}",
             "     eps      Phi(eps)          N(eps)      R(eps) = K_eff/K_2   R with G-1"]
    for e, n_, r, rc in zip(ends, N, R, Rc):
        lines.append(f"  {e:6.3f}  {Phi(e):+.10e}  {n_:+.6f}   {r:+.6f}            {rc:+.6f}")
    gam = []
    for (ea, ma), (eb, mb) in zip(zip(ends[:-1], M[:-1]), zip(ends[1:], M[1:])):
        gam.append(2 * math.log(mb / ma) / math.log(eb / ea) if ma * mb > 0 else float("nan"))
        lines.append(f"  gamma_T {ea:g} -> {eb:g}: {gam[-1]:+.4f}  (second order 5)")
    sign_kept = all(m * C2 > 0 for m in M)
    g_last = gam[-1]
    decreasing = all(b_ < a_ for a_, b_ in zip(gam[:-1], gam[1:]))
    dg = max(abs(a_ - b_) for a_, b_ in zip(R, Rc))
    inst_ok = nray == 1 and lin <= p["tol_lin"] and dg <= p["tol_g"] and max(dth) <= p["tol_theta"]
    lines += [f"  instrument: rays per phase {nray}; stationarity {lin:.1e}; |R(G) - R(G-1)| {dg:.1e} "
              f"(<= {p['tol_g']:g}); phases {nt // 2} vs {nt}: {max(dth):.1e} (<= {p['tol_theta']:g}); "
              "K.18 and the F6.8e3 check are separate records",
              f"  checks: k = 0 (eps = 0) largest |Phi| {k0:.1e}; both signs of b exact; zero age M ~ eps^5/2; "
              f"null model: largest |R - 1| {max(abs(r - 1) for r in R):.4f} (second order holds within "
              f"{p['tol_null']:g}: {max(abs(r - 1) for r in R) <= p['tol_null']})"]
    if not inst_ok:
        status, reading = "INSTRUMENT_FAILED", "-"
    elif not sign_kept or (math.isfinite(g_last) and g_last >= p["gamma_B"]):
        status, reading = "FAILED", "B"
    elif sign_kept and g_last <= p["gamma_A"] and decreasing:
        status, reading = "PASSED", "A"
    else:
        status, reading = "UNDECIDED", "C"
    lines.append(f"  reading {reading}")
    lines.append(f"  cost {time.time() - t0:.0f} s")
    return Result(status, lines, {"ends": ends, "Phi": [Phi(e) for e in ends], "C2": C2, "N": N, "M": M, "R": R,
                                  "R_g_minus_1": Rc, "gamma_T": gam, "reading": reading, "nray": nray, "lin": lin,
                                  "dphi_theta": dth, "k0": k0})
