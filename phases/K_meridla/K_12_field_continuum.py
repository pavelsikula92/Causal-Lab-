"""
K.12 -- Does the self-consistent causal-set field of (a) in 2+1 follow the continuum, or does it feel the wall?
==============================================================================================================

Why this test exists
--------------------
F2.7 (21-22 Oct 2026, UNDECIDED) found the motion exactly on the geodesic
of the field the causal set made (H1 = 2.027 +- 0.071, geodesic in that
field P_cs = 2.023), but that field is steeper in the motion region than
the continuum solution of the same equation (P_cs / P_a = 1.21). The saved
profile (R = 100) showed why: in 2+1 equation (a) reads
    lap Pi = (kappa / 3) Pi^2 rho_s,
so outside the source the flux r dPi/dr must be constant (Gauss's law; the
continuum keeps 0.1197 from r ~ 30 to the wall). The causal-set field does
not: r dPi/dr peaks at 0.129 near r = 22 and then falls towards the wall
(0.111 at 56, 0.073 at 79, 0.006 at 99), while the core depth agrees
(1 - Pi(0) 0.2887 against 0.2901); the slope ratio cs / continuum falls
smoothly from 1.28 at r = 4 to 1 near r = 43 and below 1 further out. A
different constant coupling would scale the slope by one factor
everywhere, so the profile is not harmonic outside the source. Two
candidates remain:
  (1) a boundary effect of the instrument: rows near the wall R have pasts
      that reach beyond it, where the geometry is the vacuum (Pi = 1); the
      vacuum subtraction E[B_flat] 1 cancels the cut-off exactly only in
      flat geometry. In 2+1 the logarithmic field is still well below 1
      close to the wall (continuum 0.97 at r = 80), so the effect is larger
      than in 3+1, where 1/r has nearly died out there;
  (2) screening by the nonlocal operator on the scale l_k (a property of
      the operator, i.e. physics of the model).
They differ in how the deviation depends on the wall.

Development runs before the full test (22 Oct 2026, coarse quadrature
[40, 6, 12, 3, 8, 1], decide nothing; at sigma = 8 the coarse level gives
P_cs / P_a = 1.211 against 1.209 of the full F2.7 field, so the quadrature
is not the cause there):
  wall:  sigma = 8, R = 100 -> 150: P_cs / P_a 1.211 -> 1.214, slope ratios
         in the motion region unchanged; the flux falls near the wall in both.
         The wall is not what makes the motion region steeper.
  size:  the self-similar family (total source charge q sigma^2 fixed,
         R = 12.5 sigma, b = 1.25 sigma, T = 5 sigma, node step sigma / 4:
         the continuum solution is the same function of r / sigma, Pi(0) =
         0.7099 and P_a / sigma = 0.2096 for every sigma, so any change of
         P_cs / P_a is the finite scale l_k alone):
         sigma = 4, 8, 16, 32 -> P_cs / P_a = 1.406, 1.211, 1.128, 1.137.
         The excess falls from 4 to 16 and then stops in this coarse run; the
         stop may be the coarse quadrature at large R (rows at t_x = 3R).
         At sigma = 32 a finer quadrature [120, 8, 48, 6, 24, 1] gives 1.108
         (Pi(0) 0.7284 against 0.7140 coarse, continuum 0.7099): at large R
         the quadrature is NOT converged at those levels, so every part
         checks its own quadrature (below).
So the question became: does the field link converge to the continuum when
the source is large against the nonlocality scale (the continuum limit of
the theory), and from which sigma is an F2.7-type test fair?

Parts (one per sigma, the self-similar family; the instrument of F2.7:
grid ratio 1.12, quadrature [240, 8, 96, 8, 48, 1], C^2 spline, start from
the vacuum):
  s08        sigma = 8,  q = 0.01,     R = 100, b = 10, T = 40,  step 2  (the F2.7 field; optional)
  s16        sigma = 16, q = 0.0025,   R = 200, b = 20, T = 80,  step 4
  s32        sigma = 32, q = 0.000625, R = 400, b = 40, T = 160, step 8
  s08_R200   sigma = 8,  q = 0.01,     R = 200, b = 10, T = 40,  step 2  (the wall check at the full level)
Every part is compared with the continuum solution of (a) with the same
wall, and reports the F2.7 field and the earlier K.12 parts beside its own
(re-read from the saved records), with P_cs / P_a against sigma.

Criterion (calibration of the field link for an F2.7-type test, stated
before the full runs)
---------------------------------------------------------------------
The quantity F2.7 needs from the field is the geodesic displacement in the
motion region. F2.7 decides with the band d_alt / 4 = P_a / 4 (the
alternative "no field" lies at distance P_a); by the quarter rule of the
laboratory the instrument may use at most a quarter of that band:
    PASSED   |P_cs - P_a| <= P_a / 16     (the field link is calibrated at this sigma)
    FAILED   otherwise
INSTRUMENT_FAILED if a solve did not converge, the RMS residual of the
main solve exceeds res_tol * q (as F2.7), or the quadrature is not
converged at this sigma: every part solves the field twice, at the
check level [160, 8, 48, 6, 24, 1] (K.11c level 1, from the vacuum) and at
the main level [240, 8, 96, 8, 48, 1] (K.11c level 2, the F2.7 instrument,
started from the check solution), and requires
|P_cs(main) - P_cs(check)| <= P_a / 64 (a quarter of the tolerance, the
quarter rule again). K.11c calibrated the levels only at R = 100; this
extends the check to every R of the scan. Reported, not decided: P_cs / P_a against
sigma over the parts available; the slope ratio cs / continuum at
r = 1.25, 2.5, 3.75 sigma; the flux r dPi/dr of both fields at r = 3 sigma
and at 1/4, 1/2, 3/4, 9/10 of R; the Gauss deficit; Pi(0). Reading
declared now: P_cs / P_a falling towards 1 with sigma = a finite-l_k
correction that vanishes in the continuum limit (then F2.7 is repeated at
the smallest sigma that PASSED); a floor that does not fall = a property of
the operator that the continuum chain of the theory does not contain.

Seen before this test was written: the F2.7 field and the development runs
above. The criterion uses only the F2.7 band and the quarter rule.
"""

from __future__ import annotations

import time

import numpy as np

from core import provenance
from core.chains.geodesics import geodesic_1p1
from core.operators.bd_expect import shape_nodes
from core.operators.bd_expect_curved import continuum_radial, pi_geometry, solve_selfconsistent
from core.stats.parallel import parallel_map
from core.testing import Result

TEST = {
    "id": "K.12",
    "phase": "K",
    "title": "Pole kauzální množiny (a) ve 2+1 proti spojitému řešení: konverguje s velikostí zdroje vůči ℓ_k?",
    "theory_part": "měřidlo článku zdroj → pole pro F2.7 (samokonzistentní rovnice (a) ve 2+1) a jeho spojitá limita; "
                   "bez fyzikálního závěru",
    "question": "Blíží se geodetika v samokonzistentním poli (a) kauzální množiny ve 2+1 geodetice spojitého řešení, když "
                "zdroj roste vůči nelokální škále ℓ_k (σ = 8, 16, 32 v sobě podobné rodině), a od kterého σ je odchylka "
                "v toleranci testu typu F2.7? Mění ji poloha stěny?",
    "method": "pole přesně jako F2.7 (síť s krokem σ/4 a poměrem 1,12, kvadratura [240, 8, 96, 8, 48, 1], C² splajn, start "
              "z vakua) v sobě podobné rodině (náboj qσ² pevný, R = 12,5σ, b = 1,25σ, T = 5σ) a jedna kontrola stěny "
              "(σ = 8, R = 200); srovnání se spojitým řešením (a) se stejným R: geodetika v oblasti pohybu, poměr sklonů, "
              "tok r·dΠ/dr, Gaussův deficit; P_cs/P_a proti σ přes dostupné části",
    "prediction": "bez předpovědi; vývojové běhy s hrubou kvadraturou: 1,41 / 1,21 / 1,13 / 1,14 pro σ = 4 / 8 / 16 / 32, "
                  "stěna 100 → 150 bez vlivu",
    "criterion": "MĚŘIDLO: obě řešení zkonvergovala, reziduum <= res_tol·q, kvadratura zkonvergovala: |P_cs(hlavní "
                 "úroveň) − P_cs(kontrolní úroveň K.11c L1)| <= P_a/64 (čtvrtina tolerance). PROŠEL: |P_cs − P_a| <= P_a/16 (čtvrtina pásma F2.7 "
                 "d_alt/4 = P_a/4). NEPROŠEL: jinak. P_cs/P_a proti σ, poměr sklonů, tok a Gaussův deficit hlášeny.",
    "null_model": "—",
    "level": "calibration",
    "meters": ["K.11", "K.11c"],
    "depends_on": ["K.11c", "F2.7"],
    "params": {"D": 3, "ratio": 1.12, "rows_per_element": 2, "quad": [240, 8, 96, 8, 48, 1], "quad_check": [160, 8, 48, 6, 24, 1],
               "geom": "spline",
               "tol": 1e-5, "max_iter": 30, "res_tol": 0.05, "chunk": 4, "motion_fracs": [1.25, 2.5, 3.75],
               "flux_fracs": [0.25, 0.5, 0.75, 0.9], "h_frac": 0.0625, "ref_grid": 40001},
    "parts": {"s08": {"sigma": 8.0, "q": 0.01, "R": 100.0, "b": 10.0, "T": 40.0, "dr": 2.0},
              "s16": {"sigma": 16.0, "q": 0.0025, "R": 200.0, "b": 20.0, "T": 80.0, "dr": 4.0},
              "s32": {"sigma": 32.0, "q": 0.000625, "R": 400.0, "b": 40.0, "T": 160.0, "dr": 8.0},
              "s08_R200": {"sigma": 8.0, "q": 0.01, "R": 200.0, "b": 10.0, "T": 40.0, "dr": 2.0}},
    "quick": {"quad": [40, 6, 12, 3, 8, 1], "quad_check": [24, 4, 8, 2, 6, 1], "tol": 1e-4, "max_iter": 12,
              "ref_grid": 6001},
    "quick_parts": {"s08": {"sigma": 3.0, "q": 0.0711, "R": 37.5, "b": 3.75, "T": 15.0, "dr": 0.75},
                    "s16": {"sigma": 6.0, "q": 0.0178, "R": 75.0, "b": 7.5, "T": 30.0, "dr": 1.5},
                    "s32": {"sigma": 12.0, "q": 0.00444, "R": 150.0, "b": 15.0, "T": 60.0, "dr": 3.0},
                    "s08_R200": {"sigma": 3.0, "q": 0.0711, "R": 75.0, "b": 3.75, "T": 15.0, "dr": 0.75}},
    "estimate": "každá část ~4–6 h s 10 procesy (kontrolní úroveň z vakua + hlavní úroveň z ní; stejný počet uzlů jako "
                "pole F2.7; s32 může trvat déle); s08 je volitelná (hlavní úroveň se shoduje s polem F2.7)",
}

_EPS, _RHO = 0.05, 20.0          # the field computation of F2.7: eps * rho = 1, lengths in units of l_k


def geodesic_mid(Pi_of_r, b, T, power=1):
    """Mid-time displacement away from the source of the geodesic of Pi^(2 power) eta between events at distance b (F2.7)."""
    return float(geodesic_1p1(lambda y: float(Pi_of_r(b + y)) ** power, T, 0.0, fracs=(0.5,))["disp"][0.5])


def flux(f, r, h):
    """r dPi/dr by a central difference of step h (the Gauss-law flux per unit angle in 2+1)."""
    r = np.asarray(r, float)
    return r * (np.asarray(f(r + h), float) - np.asarray(f(np.maximum(r - h, 0.0)), float)) / (2.0 * h)


def profile_metrics(Pi_cs_f, R, p):
    """Everything this test reports for one causal-set field (callable) with wall R, against the continuum of (a)."""
    sig, qs = float(p["sigma"]), float(p["q"])
    src = lambda r: qs * np.exp(-np.asarray(r, float) ** 2 / (2.0 * sig ** 2))
    c = continuum_radial("a", int(p["D"]), R, src, n_grid=int(p["ref_grid"]))
    Pi_a_f = lambda r, c=c: np.interp(r, c["r"], c["Pi"])
    b, T, h = float(p["b"]), float(p["T"]), float(p["h_frac"]) * sig
    rm = np.asarray(p["motion_fracs"], float) * sig
    rf = np.r_[3.0 * sig, np.asarray(p["flux_fracs"], float) * R]
    P_cs, P_a = geodesic_mid(Pi_cs_f, b, T), geodesic_mid(Pi_a_f, b, T)
    fl_cs, fl_a = flux(Pi_cs_f, rf, h), flux(Pi_a_f, rf, h)
    return {"R": R, "sigma": sig, "P_cs": P_cs, "P_a": P_a, "ratio_P": P_cs / P_a,
            "r_motion": rm, "slope_ratio": flux(Pi_cs_f, rm, h) / flux(Pi_a_f, rm, h),
            "r_flux": rf, "flux_cs": fl_cs, "flux_a": fl_a,
            "gauss_deficit_cs": 1.0 - fl_cs[2] / fl_cs[0], "gauss_deficit_a": 1.0 - fl_a[2] / fl_a[0],
            "Pi0_cs": float(Pi_cs_f(0.0)), "Pi0_a": float(Pi_a_f(0.0))}


def metric_lines(m, tag):
    """The report lines of one field."""
    return [f"  {tag}: sigma = {m['sigma']:g}, R = {m['R']:g}; Pi(0) cs {m['Pi0_cs']:.4f} continuum {m['Pi0_a']:.4f}; P_cs = {m['P_cs']:+.4f}, "
            f"P_a = {m['P_a']:+.4f}, P_cs/P_a = {m['ratio_P']:.3f}",
            "     slope ratio cs/continuum at r = " + ", ".join(f"{r:g}" for r in m["r_motion"]) + ": "
            + " ".join(f"{x:.3f}" for x in m["slope_ratio"]),
            "     flux r dPi/dr at r = " + ", ".join(f"{r:.4g}" for r in m["r_flux"]) + ":  cs "
            + " ".join(f"{x:.4f}" for x in m["flux_cs"]) + "   continuum " + " ".join(f"{x:.4f}" for x in m["flux_a"]),
            f"     Gauss deficit 1 - flux(R/2)/flux(3 sigma): cs {m['gauss_deficit_cs']:+.3f}, continuum "
            f"{m['gauss_deficit_a']:+.3f}"]


def earlier_fields(p, ctx):
    """The F2.7 field and earlier full K.12 parts, re-read from their saved records and analysed the same way."""
    out = []
    if ctx.quick:
        return out
    for rec in provenance.load_all("F2.7"):
        d = rec.get("data") or {}
        if d.get("nodes") and d.get("Pi_cs"):
            rp = dict(p, **{k: rec["params"][k] for k in ("sigma", "q", "b", "T")})
            f = pi_geometry(np.asarray(d["nodes"], float), np.asarray(d["Pi_cs"], float), str(rec["params"]["geom"]))
            out.append((f"F2.7 field ({rec['time_utc']})", profile_metrics(f, float(rec["params"]["R"]), rp)))
    for rec in provenance.load_all("K.12"):
        d = rec.get("data") or {}
        if d.get("nodes") and d.get("Pi_cs") and rec.get("part") != ctx.part and "sigma" in rec["params"]:
            rp = dict(p, **{k: rec["params"][k] for k in ("sigma", "q", "b", "T")})
            f = pi_geometry(np.asarray(d["nodes"], float), np.asarray(d["Pi_cs"], float), str(rec["params"]["geom"]))
            out.append((f"K.12 {rec['part']} ({rec['time_utc']}, {rec['status']})",
                        profile_metrics(f, float(rec["params"]["R"]), rp)))
    return out


def run(p, ctx):
    D, R = int(p["D"]), float(p["R"])
    sig, qs = float(p["sigma"]), float(p["q"])
    t0 = time.time()
    nodes = shape_nodes(R, 2.0 * sig, float(p["dr"]), float(p["ratio"]))
    per = int(p["rows_per_element"])
    fr = (np.arange(per) + 0.5) / per
    rows = np.sort((nodes[:-1, None] + fr[None, :] * np.diff(nodes)[:, None]).ravel())
    src = lambda r: qs * np.exp(-np.asarray(r, float) ** 2 / (2.0 * sig ** 2))
    kw = dict(rho=_RHO, eps=_EPS, tol=float(p["tol"]), max_iter=int(p["max_iter"]), chunk=int(p["chunk"]),
              geom=str(p["geom"]))
    # the coarser check level first (from the vacuum), then the instrument level of F2.7 started from it (the fixed
    # point does not depend on the start; starting near it saves iterations)
    chk = solve_selfconsistent(D, R, nodes, rows, src(rows), [int(x) for x in p["quad_check"]], parallel_map,
                               ctx.workers, ctx.progress, label=f"[K.12 sigma{sig:g} R{R:g} check level]", **kw)
    ctx.progress(f"[K.12] check level: Pi(0) = {chk['Pi'][0]:.4f}, {chk['iterations']} iterations ({time.time()-t0:.0f}s)")
    sol = solve_selfconsistent(D, R, nodes, rows, src(rows), [int(x) for x in p["quad"]], parallel_map, ctx.workers,
                               ctx.progress, label=f"[K.12 sigma{sig:g} R{R:g} main level]",
                               Pi0=np.asarray(chk["Pi"], float), **kw)
    Pi_cs = np.asarray(sol["Pi"], float)
    ctx.progress(f"[K.12] field sigma = {sig:g}, R = {R:g}: Pi(0) = {Pi_cs[0]:.4f}, {sol['iterations']} iterations, "
                 f"converged {sol['converged']} ({time.time()-t0:.0f}s)")
    m = profile_metrics(pi_geometry(nodes, Pi_cs, str(p["geom"])), R, p)
    m_chk = profile_metrics(pi_geometry(nodes, np.asarray(chk["Pi"], float), str(p["geom"])), R, p)
    tol_P = m["P_a"] / 16.0
    dq = abs(m["P_cs"] - m_chk["P_cs"])
    inst = {"field solve converged": bool(sol["converged"]) and bool(chk["converged"]),
            "field residual": float(sol["residual"]) <= float(p["res_tol"]) * qs,
            "quadrature converged": dq <= tol_P / 4.0}
    dev = abs(m["P_cs"] - m["P_a"])
    if not all(inst.values()):
        status = "INSTRUMENT_FAILED"
    elif dev <= tol_P:
        status = "PASSED"
    else:
        status = "FAILED"
    lines = [f"2+1, lengths in l_k: source sigma = {sig:g}, q = {qs:g}; field (a) self-consistent, R = {R:g}, {len(nodes)} nodes, "
             f"{len(rows)} rows, quadrature {list(p['quad'])}, geometry {p['geom']}; {sol['iterations']} iterations, "
             f"residual {float(sol['residual']):.1e}; motion region: b = {p['b']:g}, T = {p['T']:g}"]
    lines += metric_lines(m, "this part")
    lines.append(f"  check level {list(p['quad_check'])}: P_cs = {m_chk['P_cs']:+.4f}, Pi(0) {m_chk['Pi0_cs']:.4f}; "
                 f"|P_cs(main) - P_cs(check)| = {dq:.4f} vs P_a/64 = {tol_P / 4:.4f}")
    lines.append(f"  decision: |P_cs - P_a| = {dev:.4f} vs P_a/16 = {tol_P:.4f} (a quarter of the F2.7 band P_a/4)")
    prev = earlier_fields(p, ctx)
    if prev:
        lines.append("  for comparison (re-read from saved records, analysed the same way):")
        for tag, mm in prev:
            lines += metric_lines(mm, tag)
    allm = [("this part", m)] + prev
    lines.append("  P_cs / P_a against sigma (all fields above): " + "; ".join(
        f"sigma {mm['sigma']:g} R {mm['R']:g}: {mm['ratio_P']:.3f}" for _, mm in sorted(allm, key=lambda x: (x[1]['sigma'], x[1]['R']))))
    lines.append("instrument: " + "; ".join(f"{k}: {'OK' if v else 'FAILS'}" for k, v in inst.items()))
    data = {"nodes": nodes, "Pi_cs": Pi_cs, "Pi_check": np.asarray(chk["Pi"], float), "metrics_check": m_chk,
            "solve": {k: sol[k] for k in ("iterations", "converged", "residual", "steps")},
            "solve_check": {k: chk[k] for k in ("iterations", "converged", "residual", "steps")},
            "metrics": m, "tolerance": tol_P, "instrument": inst,
            "earlier": [{"tag": tag, "metrics": mm} for tag, mm in prev]}
    return Result(status, lines, data)
