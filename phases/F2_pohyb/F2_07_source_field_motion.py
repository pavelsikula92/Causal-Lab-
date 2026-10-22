"""
F2.7 -- Closing chain: source -> field -> motion (2+1)
======================================================

Why this test exists
--------------------
The motion tests so far (F2.2, F2.4, F2.5; E46, E46b) put the field in by
hand: Pi = 1 - eps / sqrt(r^2 + soft^2). The field tests (F1.3, F1.4, F1.5)
stop at the field. Here nothing is put in by hand except the source: the
field is the one the causal set itself produces from the source through the
self-consistent field equation (a) (author's choice, 18 Oct 2026),

    B_g Pi = kappa rho_s,     g = Pi^2 eta,

and the motion is read from the heaviest chains of the action W = sum Pi
on a sprinkling, with that field as the clock rate of every element.

What is new compared with F1.5 + F2.4, stated honestly: each link was
tested on its own (F1.5: the field of (a) in 2+1, preliminary; F2.4: chains
follow the geodesic of an inserted field). This test checks that the two
links compose: that the causal-set field (nonlocal, with the wall offset and
the coupling kappa / D of (a)) moves bodies as the continuum chain of the
theory predicts, in the same units, and that the result is told apart from
the competing field equations and from no field at all. A failure here
without a failure of F1.5 or F2.4 would mean the errors of the two links
add up beyond the tolerance, or that the units do not connect.

Setup (lengths in units of the nonlocality scale l_k = (eps rho)^(-1/d))
------------------------------------------------------------------------
Source: Gaussian q exp(-r^2 / 2 sigma^2) at the origin, sigma = 8 l_k,
q = 0.01 (the core source of F1.4; central Pi ~0.71 in (a)).
Field: the self-consistent solution of E[B_g(Pi)] Pi = q rho_hat in a
cylinder of radius R = 100 l_k (core.operators.bd_expect_curved.
solve_selfconsistent; the flat operator on the vacuum subtracted), node
grid g2 (step 2 l_k, ratio 1.12), quadrature [240, 8, 96, 8, 48, 1], the
geometry a C^2 spline through the nodal values (geom = 'spline'). This is
the instrument K.11c far_d3_g2 calibrates; the test may run only after
K.11c far_d3_g2 PASSED with this level within tolerance (dependency).
Starting field: the vacuum (not the prediction).
Motion: end events p = (0, x0) and q = (T, x0) at distance b = 10 l_k
from the source, T = 40 l_k (the field gradient over T is then comparable
with E46: predicted displacement ~0.04 T). N0 points uniform in I(p, q);
weights w = Pi(|x - m|) (H1, the action of the theory) and Pi^2 (H2,
reported), Pi the causal-set field interpolated with the same spline.
Reading: the mid-time displacement away from the source, uniform average
over all heaviest chains, mirror-antisymmetrised with common points
(source at -b e and +b e; F2.5).
Units: the expected operator depends on eps and rho only through l_k
(core.operators.bd_expect). The sprinkling of the motion has density
rho_c = N0 / V(T) per l_k^3; the same causal set has the field of this
test when its smearing is eps = 1 / rho_c (0.52 for N0 = 32 000, T = 40;
eps <= 1 is checked).

Predictions (continuum, the whole chain of the theory without the causal set)
-----------------------------------------------------------------------------
  P_a    geodesic of Pi_a^2 eta, Pi_a the continuum solution of (a)
         (local form of B_g: D Pi^-2 lap Pi + c_D Pi^-3 |grad Pi|^2)
Alternatives (the same source, the same geodesic rule):
  box_g  field of box_g Pi = kappa rho_s (no -R/2: coupling kappa, not kappa/D)
  flat   field of lap Pi = kappa rho_s (no back-reaction; at this source
         Pi < 0 near the centre: no static solution, then not an alternative)
  none   no field: displacement 0
Reported: P_cs, the geodesic of the causal-set field itself (separates the
field link from the motion link), and P2_a, the geodesic of Pi_a^4 eta (H2).

Criterion (stated before the run; the quarter rule of F1.3 with 2 SE)
---------------------------------------------------------------------
d_alt = min over the existing alternatives of |P_alt - P_a|; band = d_alt / 4.
  PASSED     |mean(H1) - P_a| + 2 SE <= band
  FAILED     |mean(H1) - P_a| - 2 SE >  band
  UNDECIDED  otherwise
SE: bootstrap over whole seeds (B resamples).
Instrument (any failing -> INSTRUMENT_FAILED):
  (a) the self-consistent solve converged and RMS residual <= res_tol * q;
  (b) mirror identity at Pi = 1 at N0: |d(P) + d(MP)| <= 1e-12 for H1 and
      H2 on mirror_seeds seeds, readings non-trivial (sd > 0) (K.7, E46);
  (c) path weights sum to 1 within 1e-6 on every seed;
  (d) seeds with an undefined reading <= 5 %;
  (e) eps = 1 / rho_c <= 1.

Seen before this test was written (development, 19 Oct): with the
CONTINUUM (a) field as weights, N0 = 16 000 and 8 seeds, H1 = 1.94 +- 1.54
(sd per seed) against P_a = 1.68, H2 = 2.91 +- 2.36 against P2_a = 3.32;
this fixed T = 40 and the number of seeds (SE ~0.1 l_k), nothing else.
Known: in 2+1 H1 reads ~10 % above the geodesic in the E46 field
(docs/nalezy.md; F2.5 d3 measures whether it is a finite-N effect). The
criterion is not widened for it.
"""

from __future__ import annotations

import math
import time

import numpy as np

from core.chains.chains import heaviest_between
from core.chains.geodesics import geodesic_1p1
from core.chains.readings import with_ends
from core.operators.bd_expect import shape_nodes
from core.operators.bd_expect_curved import continuum_radial, pi_geometry, solve_selfconsistent
from core.order.sprinkle import Interval, derive_seed
from core.stats.parallel import parallel_map
from core.testing import Result

TEST = {
    "id": "F2.7",
    "phase": "F2",
    "title": "Uzavírací řetězec zdroj → pole → pohyb (2+1)",
    "theory_part": "A5′ jako samokonzistentní rovnice (a) B_g Π = κρ_s a akce W = Σ Π: pole vytvořené zdrojem řídí pohyb",
    "question": "Pohybují se nejtěžší řetězce Σ Π v poli, které kauzální množina sama vytvoří ze zdroje rovnicí (a), "
                "po geodetice, kterou předpovídá spojitý řetězec teorie, a odliší se to od konkurenčních rovnic pole a od "
                "nulového pole?",
    "method": "gaussovský zdroj σ = 8 l_k, q = 0,01; samokonzistentní pole (a) (K.11c, C² geometrie); sprinkling intervalu "
              "T = 40 l_k ve vzdálenosti b = 10 l_k od zdroje, váhy Π z pole kauzální množiny; zrcadlově antisymetrizovaný "
              "posun v polovině času; bootstrap přes seedy",
    "prediction": "H1 → P_a, geodetika Π_a²η, Π_a spojité řešení (a)",
    "criterion": "MĚŘIDLO: řešení pole zkonvergovalo a reziduum <= res_tol·q; zrcadlová identita při Π = 1; váhy 1 ± 1e-6; "
                 "neplatné seedy <= 5 %; eps = 1/ρ_c <= 1. d_alt = min |P_alt − P_a| (□_g bez −R/2, bez zpětné vazby, "
                 "pokud existuje statické řešení, a nulové pole). PROŠEL: |H1 − P_a| + 2 SE <= d_alt/4. NEPROŠEL: "
                 "|H1 − P_a| − 2 SE > d_alt/4. NEROZHODNUTO: jinak.",
    "null_model": "nulové pole (posun 0) a konkurenční rovnice pole (□_g bez −R/2; bez zpětné vazby)",
    "level": "chain",
    "meters": ["K.7", "K.11", "K.11c"],
    "depends_on": ["K.11c", "F1.5", "F2.4", "F2.5"],
    "params": {"D": 3, "sigma": 8.0, "q": 0.01, "R": 100.0, "dr": 2.0, "ratio": 1.12, "rows_per_element": 2,
               "quad": [240, 8, 96, 8, 48, 1], "geom": "spline", "tol": 1e-5, "max_iter": 30, "res_tol": 0.05,
               "chunk": 4, "b": 10.0, "T": 40.0, "N0": 32000, "n_seeds": 512, "boot_B": 4000, "n_se": 2.0,
               "invalid_max": 0.05, "mirror_tol": 1e-12, "weight_tol": 1e-6, "mirror_seeds": 4},
    "parts": {},
    "quick": {"sigma": 3.0, "R": 30.0, "dr": 1.0, "ratio": 1.2, "quad": [40, 6, 12, 3, 8, 1], "tol": 1e-4,
              "max_iter": 12, "b": 4.0, "T": 12.0, "N0": 1500, "n_seeds": 6, "boot_B": 200, "mirror_seeds": 2},
    "estimate": "pole: ~10 iterací po stavbě ~45 řádků jemnou kvadraturou (~1–2 h s 10 procesy); pohyb: ~250 s na seed "
                "při N0 = 32 000, 512 seedů ~35 h CPU, ~3–4 h s 10 procesy",
}

_EPS, _RHO = 0.05, 20.0          # the field computation: eps * rho = 1, lengths in units of l_k


def one_seed(D, N0, seed, nodes, Pi_nodes, geom, b, T):
    """Mirror-antisymmetrised H1 and H2 readings of one sprinkling in the causal-set field, and the worst weight error."""
    Pi = pi_geometry(nodes, Pi_nodes, geom)
    rng = np.random.default_rng(seed)
    X = with_ends(Interval(D, T).propose(rng, N0), D, T)
    q = len(X) - 1
    e = np.zeros(D - 1)
    e[0] = 1.0
    d = {"H1": {}, "H2": {}}
    wchk = []
    for sgn in (1, -1):
        m = -sgn * b * e                      # source on the -e side for sgn = +1: "away" = +e
        w = Pi(np.sqrt(((X[:, 1:] - m) ** 2).sum(axis=1)))
        h1 = heaviest_between(X, w, 0, q, sgn * e)
        h2 = heaviest_between(X, w * w, 0, q, sgn * e)
        d["H1"][sgn], d["H2"][sgn] = h1["disp"], h2["disp"]
        wchk += [h1.get("weight_check", 0.0), h2.get("weight_check", 0.0)]
    return {"H1": 0.5 * (d["H1"][1] + d["H1"][-1]), "H2": 0.5 * (d["H2"][1] + d["H2"][-1]),
            "weight_worst": max(abs(c - 1.0) for c in wchk)}


def mirror_task(D, N0, seed, T):
    """Pi = 1: one-sided readings on P and on its mirror MP; |d(P) + d(MP)| (H1 = H2 here) and d(P)."""
    rng = np.random.default_rng(seed)
    P = Interval(D, T).propose(rng, N0)
    PM = P.copy()
    PM[:, 1] = -PM[:, 1]
    e = np.zeros(D - 1)
    e[0] = 1.0
    vals = []
    for pts in (P, PM):
        X = with_ends(pts, D, T)
        vals.append(heaviest_between(X, np.ones(len(X)), 0, len(X) - 1, e)["disp"])
    return abs(vals[0] + vals[1]), vals[0]


def geodesic_mid(Pi_of_r, b, T, power=1):
    """Mid-time displacement away from the source of the geodesic of Pi^(2 power) eta between events at distance b."""
    return float(geodesic_1p1(lambda y: float(Pi_of_r(b + y)) ** power, T, 0.0, fracs=(0.5,))["disp"][0.5])


def run(p, ctx):
    D = int(p["D"])
    R, sig, qs = float(p["R"]), float(p["sigma"]), float(p["q"])
    b, T, N0, ns = float(p["b"]), float(p["T"]), int(p["N0"]), int(p["n_seeds"])
    geom = str(p["geom"])
    quad = [int(x) for x in p["quad"]]
    t0 = time.time()

    # ---- the field the causal set produces from the source --------------------------------------------------------
    nodes = shape_nodes(R, 2.0 * sig, float(p["dr"]), float(p["ratio"]))
    per = int(p["rows_per_element"])
    fr = (np.arange(per) + 0.5) / per
    rows = np.sort((nodes[:-1, None] + fr[None, :] * np.diff(nodes)[:, None]).ravel())
    src = lambda r: qs * np.exp(-np.asarray(r, float) ** 2 / (2.0 * sig ** 2))
    sol = solve_selfconsistent(D, R, nodes, rows, src(rows), quad, parallel_map, ctx.workers, ctx.progress,
                               label=f"[F2.7] {D-1}+1 field", rho=_RHO, eps=_EPS, tol=float(p["tol"]),
                               max_iter=int(p["max_iter"]), chunk=int(p["chunk"]), geom=geom)
    Pi_cs = sol["Pi"]
    Pi_cs_f = pi_geometry(nodes, Pi_cs, geom)
    ctx.progress(f"[F2.7] field: Pi(0) = {Pi_cs[0]:.4f}, {sol['iterations']} iterations, converged {sol['converged']} "
                 f"({time.time()-t0:.0f}s)")

    # ---- continuum predictions -------------------------------------------------------------------------------------
    def cont(model):
        try:
            c = continuum_radial(model, D, R, src)
        except ValueError:
            return None                                     # no static solution of that strength
        if float(np.min(c["Pi"])) <= 0.0:
            return None
        return lambda r, c=c: np.interp(r, c["r"], c["Pi"])
    fields = {m: cont(m) for m in ("a", "box_g", "flat")}
    P_a = geodesic_mid(fields["a"], b, T)
    P2_a = geodesic_mid(fields["a"], b, T, power=2)
    P_cs = geodesic_mid(Pi_cs_f, b, T)
    alts = {"none": 0.0}
    for m in ("box_g", "flat"):
        alts[m] = geodesic_mid(fields[m], b, T) if fields[m] is not None else float("nan")
    d_alt = min(abs(v - P_a) for v in alts.values() if np.isfinite(v))
    band = d_alt / 4.0

    # ---- motion: heaviest chains in the causal-set field -----------------------------------------------------------
    rho_c = N0 / Interval(D, T).volume()
    eps_c = 1.0 / rho_c
    args = [(D, N0, derive_seed(ctx.seed, "F2.7", D, N0, k), nodes, Pi_cs, geom, b, T) for k in range(ns)]
    res = parallel_map(one_seed, args, ctx.workers, ctx.progress,
                       lambda i, r_, dn, nn: f"[F2.7] seed {dn}/{nn} H1 {r_['H1']:+.4f} H2 {r_['H2']:+.4f} "
                                             f"({time.time()-t0:.0f}s)")
    valid = [r_ for r_ in res if np.isfinite(r_["H1"]) and np.isfinite(r_["H2"])]
    h1 = np.array([r_["H1"] for r_ in valid])
    h2 = np.array([r_["H2"] for r_ in valid])
    wworst = max(r_["weight_worst"] for r_ in res)
    B = int(p["boot_B"])
    idx = np.random.default_rng(derive_seed(ctx.seed, "F2.7-boot", D, N0)).integers(0, h1.size, size=(B, h1.size))
    st = {k: {"mean": float(v.mean()), "se": float(v[idx].mean(axis=1).std(ddof=1)),
              "sd": float(v.std(ddof=1))} for k, v in (("H1", h1), ("H2", h2))}
    ms = [derive_seed(ctx.seed, "F2.7-mirror", D, N0, k) for k in range(int(p["mirror_seeds"]))]
    mres = parallel_map(mirror_task, [(D, N0, s, T) for s in ms], ctx.workers, ctx.progress,
                        lambda i, r_, dn, nn: f"[F2.7] mirror {dn}/{nn} ({time.time()-t0:.0f}s)")
    mirror_worst = max(m_[0] for m_ in mres)
    mirror_sd = float(np.std([m_[1] for m_ in mres])) if len(mres) > 1 else float("nan")

    inst = {"field solve converged": bool(sol["converged"]),
            "field residual": sol["residual"] <= float(p["res_tol"]) * qs,
            "mirror identity": mirror_worst <= float(p["mirror_tol"]) and (len(mres) < 2 or mirror_sd > 0),
            "weights": wworst <= float(p["weight_tol"]),
            "invalid seeds": (ns - len(valid)) <= float(p["invalid_max"]) * ns,
            "eps <= 1": eps_c <= 1.0}
    nse = float(p["n_se"])
    dev = abs(st["H1"]["mean"] - P_a)
    se = st["H1"]["se"]
    if not all(inst.values()):
        status = "INSTRUMENT_FAILED"
    elif dev + nse * se <= band:
        status = "PASSED"
    elif dev - nse * se > band:
        status = "FAILED"
    else:
        status = "UNDECIDED"

    f = lambda x: "no static solution" if not np.isfinite(x) else f"{x:+.4f}"
    lines = [f"{D-1}+1, lengths in l_k: source sigma = {sig:g}, q = {qs:g}; field (a) self-consistent, R = {R:g}, "
             f"{nodes.size} nodes, {rows.size} rows, quadrature {quad}, geometry {geom}; Pi_cs(0) = {Pi_cs[0]:.4f} "
             f"(continuum (a) {float(fields['a'](0.0)):.4f}), {sol['iterations']} iterations, residual {sol['residual']:.1e}",
             f"  motion: b = {b:g}, T = {T:g}, N0 = {N0}, rho_c = {rho_c:.3f} per l_k^{D} (eps = {eps_c:.3f}); "
             f"{len(valid)}/{ns} valid seeds",
             "  Pi along the path r = b ... b + T/2:  causal set " +
             " ".join(f"{float(Pi_cs_f(r)):.4f}" for r in np.linspace(b, b + T / 2, 5)) +
             "   continuum (a) " + " ".join(f"{float(fields['a'](r)):.4f}" for r in np.linspace(b, b + T / 2, 5)),
             f"  predictions (mid-time displacement): P_a = {P_a:+.4f} (geodesic of Pi_a^2 eta); alternatives: "
             + ", ".join(f"{k} {f(v)}" for k, v in alts.items()) + f"; d_alt = {d_alt:.4f}, band d_alt/4 = {band:.4f}",
             f"  reported: P_cs = {P_cs:+.4f} (geodesic in the causal-set field), P2_a = {P2_a:+.4f} (Pi_a^4 eta, for H2)",
             f"  H1 = {st['H1']['mean']:+.4f} +- {se:.4f} (sd {st['H1']['sd']:.3f}); |H1 - P_a| = {dev:.4f}, "
             f"relative {dev / abs(P_a):.3f}; H1 / P_a = {st['H1']['mean'] / P_a:.3f}",
             f"  H2 = {st['H2']['mean']:+.4f} +- {st['H2']['se']:.4f} (P2_a {P2_a:+.4f}, H2 / P2_a = "
             f"{st['H2']['mean'] / P2_a:.3f}; reported only)",
             f"  decision: |H1 - P_a| + {nse:g} SE = {dev + nse * se:.4f}, |H1 - P_a| - {nse:g} SE = {dev - nse * se:.4f} "
             f"vs band {band:.4f}",
             "instrument: " + "; ".join(f"{k}: {'OK' if v else 'FAILS'}" for k, v in inst.items())
             + f" (mirror worst {mirror_worst:.1e}, weights worst {wworst:.1e})"]
    data = {"nodes": nodes, "Pi_cs": Pi_cs, "solve": {k: sol[k] for k in ("iterations", "converged", "residual", "steps")},
            "P_a": P_a, "P2_a": P2_a, "P_cs": P_cs, "alternatives": alts, "d_alt": d_alt, "band": band,
            "stats": st, "H1": h1, "H2": h2, "rho_c": rho_c, "eps_c": eps_c, "instrument": inst,
            "mirror_worst": mirror_worst, "weight_worst": wworst}
    return Result(status, lines, data)
