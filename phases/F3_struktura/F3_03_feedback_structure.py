"""
F3.3 -- Does feedback of the load create a localized, persistent structure?
===========================================================================

Why this test exists
--------------------
F3.2 (E55b, archive) found that capacity-limited transport WITHOUT
feedback creates no structure: under homogeneous driving the saturated
events neither repeat across time slabs nor cluster beyond the common
envelope (Delta_pers +0.006 +- 0.006, E_loc +0.004 +- 0.004, margin 0.05).
Matter in Causal Physics is a persistent structure of high impedance; if
it is to arise from the substrate, some feedback of the load on the
transport must create what F3.2 does not show. This test switches feedback
on, one channel at a time, on the SAME causal sets as the feedback-off
routing (author's decision, 21 Oct 2026: test the channels separately,
the combination D only after A-C).

Model (core.flux.feedback)
--------------------------
Response of axiom 5: sigma_v = F_v / C_v (capped at 1), Pi_v = f(sigma_v)
= 1 - (1 - Pi_sat) sigma_v^m with Pi_sat = 0.5 and m = 2 (author, 21 Oct:
experimental parameters, not constants of the theory; F3.4 scans them).
One channel per part:
  capacity  (A) C_v = C Pi_v: slow clocks process less (positive feedback;
                the largest local solution, i.e. relaxation from Pi = 1);
  routing   (B) an event forwards to its children with weights Pi_k / I_k
                (I_k = 1 until F0.5 defines rigidity): flux prefers low
                impedance. NOT a causal update rule (the children are in the
                future of the parent): computed as the stationary
                self-consistent state (damped passes until Pi changes by
                < tol); a candidate transport rule, not a consequence of the
                current axioms;
  source    (C) internal events run at the clock rate: s_v Pi_v.

Setup (F3.2 / E55b, with one change)
------------------------------------
2+1 cylinder |x| <= R = 0.5, homogeneous driving s0 = 1, outflow boundary
|x| > R_abs = 0.4; analysis events |x| < R_abs, t in [3.5, 5.5); cells 0.1,
slabs 0.25, min_count 3, n_rbins 16, n_perm 20; C = the 0.99-quantile of
the uncapped offered flux of the analysis events of 8 independent
calibration sprinklings (stream "F3.3-cal", no feedback). The one change:
the cylinder is T = 7 instead of 6 (N = 56 000 instead of 48 000, the same
density), so that the survival of a structure after the source is
switched off can be followed for 1.25 time units (author's request). The
analysis window and its statistics are those of F3.2 (the routing of an
event depends only on its past; for 'routing' see below).

Per sprinkling and arm (off = no feedback, on = the part's channel):
  * main routing, source on for all t: saturated chi_v = 1[dQ_v > 0];
    E_loc (one value per sprinkling) and the cell sequence for Delta_pers
    (one value per disjoint pair), exactly as F3.2;
  * switch-off routing, source on only for t < t_off = 5.5: saturated
    fraction of the analysis region in five slabs after t_off and the
    Spearman r_s of the cell values of each such slab with the last slab
    before t_off (REPORTED: how long a structure survives without driving).
For 'routing' the two routings are separate stationary states (its
weights look into the future, so the switch-off would otherwise reach
back into the window).

Criteria (stated before the run)
--------------------------------
Paired differences on the same causal sets: d_loc = E_loc(on) - E_loc(off)
per sprinkling, d_pers = Delta_pers(on) - Delta_pers(off) per pair; mean
+- 2 SE (bootstrap over whole sprinklings / pairs), margin 0.05 (the F3.2
margins): EQUIVALENT if |mean| + 2 SE <= 0.05, EXCESS if |mean| - 2 SE >
0.05, else INCONCLUSIVE. "positive" = EXCESS with mean > 0.
  PASSED     d_pers and d_loc both positive: the channel creates a
             localized AND persistent structure beyond the F3.2 baseline
  FAILED     neither positive and both decided (EQUIVALENT or a negative
             EXCESS); or the channel removes saturation altogether (mean
             saturated fraction of the analysis events, on, < sat_min =
             0.001, a tenth of the F3.2 level): no saturated events, no
             structure
  UNDECIDED  otherwise (one statistic positive, or INCONCLUSIVE)
Instrument (any failing -> INSTRUMENT_FAILED): the off arm equals
core.flux.capacity.route_capacity exactly on every sprinkling; balance of
axiom 3 within 1e-12 on every routing; 'routing' converged on every
sprinkling; undefined statistics in at most 5 % of the sprinklings / pairs
(counted only where the on arm saturates).
Reported: the absolute Delta_pers and E_loc of both arms (the off arm is a
replication of F3.2 on this geometry), saturated fractions, the Pi of the
saturated events, the survival after switch-off, iterations ('routing').

Null model: the feedback-off routing on the same causal sets (the F3.2
baseline); the F3.2 degree/time DAG null is not repeated here.

Seen before this test was written (development, 21 Oct, one sprinkling of
N = 6000, T = 7, C from that sprinkling): saturated fraction of the
analysis events off 0.0098, capacity 0.428, routing 0.000 (21 passes),
source 0.000; minimum Pi 0.50, 0.60, 0.57. This fixed nothing in the
criteria except that the 'no saturation' outcome needed a rule (sat_min).
"""

from __future__ import annotations

import time
import warnings

import numpy as np

from core.flux.capacity import route_capacity
from core.flux.feedback import route_feedback, balance_residual_fb
from core.flux import structure as ST
from core.order.links import links_csr
from core.order.sprinkle import Cylinder, derive_seed, sprinkle
from core.stats.parallel import parallel_map
from core.stats.stats import boot_se, equivalence
from core.testing import Result

TEST = {
    "id": "F3.3",
    "phase": "F3",
    "title": "Vytváří zpětná vazba zatížení lokalizovanou trvalou strukturu?",
    "theory_part": "axiom 4 a 5 (kapacita, odezva Π = f(σ)); hmota jako trvalá struktura vysoké impedance",
    "question": "Vytvoří jeden kanál zpětné vazby zatížení na transport (kapacita CΠ, směrování k nízké impedanci, zdroj sΠ) "
                "lokalizovanou a trvalou strukturu saturovaných událostí nad úrovní F3.2 bez zpětné vazby?",
    "method": "2+1 válec jako F3.2 (T = 7, N = 56 000), homogenní buzení, C = Q_0,99 z nezávislých kalibrací; na stejných "
              "kauzálních množinách směrování bez zpětné vazby a s jedním kanálem (Π_sat = 0,5, m = 2); párové rozdíly "
              "Δ_pers a E_loc; přežití struktury po vypnutí zdroje (hlášeno)",
    "prediction": "bez předpovědi směru (kanál A je kladná zpětná vazba, B a C záporná); rozhoduje párový rozdíl",
    "criterion": "MĚŘIDLO: větev bez zpětné vazby přesně rovna route_capacity, bilance <= 1e-12, směrování (B) "
                 "zkonvergovalo, neplatné statistiky <= 5 %. PROŠEL: d_pers i d_loc kladně EXCESS (mez 0,05, 2 SE). "
                 "NEPROŠEL: žádný kladný a oba rozhodnuté, nebo saturace zmizí (podíl < 0,001). NEROZHODNUTO: jinak.",
    "null_model": "směrování bez zpětné vazby na stejných kauzálních množinách (základ F3.2)",
    "level": "mechanism",
    "meters": ["K.9"],
    "depends_on": ["F3.2"],
    "params": {"D": 3, "T": 7.0, "R": 0.5, "R_abs": 0.4, "s0": 1.0, "t1": 3.5, "t2": 5.5, "t_off": 5.5, "post_end": 6.75,
               "cell": 0.1, "slab": 0.25, "min_count": 3, "n_rbins": 16, "n_perm": 20, "q": 0.99, "N": 56000,
               "n_seeds": 128, "n_cal": 8, "Pi_sat": 0.5, "m": 2.0, "margin_pers": 0.05, "margin_loc": 0.05,
               "n_se": 2.0, "boot_B": 4000, "sat_min": 0.001, "balance_tol": 1e-12, "invalid_max": 0.05,
               "route_tol": 1e-6, "route_max_iter": 200, "route_damp": 0.5},
    "parts": {"capacity": {"mode": "capacity"}, "routing": {"mode": "routing"}, "source": {"mode": "source"}},
    "quick": {"N": 6000, "n_seeds": 4, "n_cal": 1, "n_perm": 4, "boot_B": 200},
    "estimate": "capacity a source ~20–40 min s 8 procesy (~1 min na sprinkling: 4 směrování a 2 lokalizace); routing "
                "~1–2 h (asi 20 průchodů na stacionární stav)",
}

F32_REFERENCE = {"Delta_pers": (0.006, 0.006), "E_loc": (0.004, 0.004)}     # F3.2 (E55b), archive


def geometry(p, seed):
    """One cylinder sprinkling, its links, radii, times and the masks of F3.2."""
    S = sprinkle(Cylinder(int(p["D"]), 0.0, float(p["T"]), float(p["R"])), int(p["N"]), seed)
    ptr, idx = links_csr(S.X, S.rho)
    r, t = S.r(), S.X[:, 0]
    absorb = r > float(p["R_abs"])
    win = (t >= float(p["t1"])) & (t < float(p["t2"])) & (r < float(p["R_abs"]))
    return S, ptr, idx, r, t, absorb, win


def calibration_task(p, seed):
    """Uncapped offered flux of the analysis events of one calibration sprinkling (no feedback)."""
    S, ptr, idx, r, t, absorb, win = geometry(p, seed)
    return route_capacity(ptr, idx, np.full(S.n, float(p["s0"])), np.inf, absorb)["offered"][win]


def survival(S, r, t, chi, p):
    """Saturated fraction in the post-switch-off slabs and r_s of each with the last slab before t_off."""
    pp = dict(p, t1=float(p["t_off"]) - float(p["slab"]), t2=float(p["post_end"]))
    reg = (r < float(p["R_abs"])) & (t >= pp["t1"]) & (t < pp["t2"])
    seq = ST.cells(S.X, chi, reg, pp)
    edges = np.arange(pp["t1"], pp["t2"] + 1e-9, float(p["slab"]))
    frac = [float(chi[reg & (t >= edges[k]) & (t < edges[k + 1])].mean()) for k in range(len(edges) - 1)]
    rs = [ST.rs_successive([seq[0], seq[k]], int(p["min_count"])) for k in range(1, len(seq))]
    return {"frac": frac, "rs_with_last_on": rs}


def seed_task(p, seed, C, mode):
    """Both arms on one sprinkling: main routing (statistics) and switch-off routing (survival)."""
    S, ptr, idx, r, t, absorb, win = geometry(p, seed)
    s_on = np.full(S.n, float(p["s0"]))
    s_off = np.where(t < float(p["t_off"]), float(p["s0"]), 0.0)
    kw = dict(Pi_sat=float(p["Pi_sat"]), m=float(p["m"]), tol=float(p["route_tol"]),
              max_iter=int(p["route_max_iter"]), damp=float(p["route_damp"]))
    out = {}
    ref = route_capacity(ptr, idx, s_on, C, absorb)
    for arm, md in (("off", "off"), ("on", mode)):
        main = route_feedback(ptr, idx, s_on, C, absorb, md, **kw)
        swo = route_feedback(ptr, idx, s_off, C, absorb, md, **kw)
        chi = (main["stored"] > 0).astype(float)
        rng = np.random.default_rng(derive_seed(seed, "F3.3-perm", arm))
        loc = ST.localization(S.X, r, chi, win, p, rng)
        sat = chi[win] > 0
        out[arm] = {"sat_fraction": float(chi[win].mean()), "E_loc": loc["E_loc"], "I_obs": loc["I_obs"],
                    "cells": ST.cells(S.X, chi, win, p),
                    "Pi_sat_median": float(np.median(main["Pi"][win][sat])) if sat.any() else float("nan"),
                    "Pi_min": float(main["Pi"][win].min()),
                    "balance": max(balance_residual_fb(main), balance_residual_fb(swo)),
                    "iterations": [int(main["iterations"]), int(swo["iterations"])],
                    "converged": bool(main["converged"] and swo["converged"]),
                    "survival": survival(S, r, t, (swo["stored"] > 0).astype(float), p)}
        if arm == "off":
            out["off_equals_capacity"] = all(np.array_equal(main[k], ref[k]) for k in ("offered", "accepted", "stored", "Phi"))
    return out


def run(p, ctx):
    mode = str(p["mode"])
    n, B, nse = int(p["n_seeds"]), int(p["boot_B"]), float(p["n_se"])
    t0 = time.time()
    cal = parallel_map(calibration_task, [(p, derive_seed(ctx.seed, "F3.3-cal", int(p["N"]), k)) for k in range(int(p["n_cal"]))],
                       ctx.workers, ctx.progress, lambda i, r_, d, nn: f"[F3.3] calibration {d}/{nn} ({time.time()-t0:.0f}s)")
    C = float(np.quantile(np.concatenate(cal), float(p["q"])))
    seeds = [derive_seed(ctx.seed, "F3.3", int(p["N"]), k) for k in range(n)]
    rows = parallel_map(seed_task, [(p, sd, C, mode) for sd in seeds], ctx.workers, ctx.progress,
                        lambda i, r_, d, nn: f"[F3.3] {mode} seed {d}/{nn}: saturated off {r_['off']['sat_fraction']:.4f} "
                                             f"on {r_['on']['sat_fraction']:.4f}; E_loc off {r_['off']['E_loc']:+.4f} on "
                                             f"{r_['on']['E_loc']:+.4f} ({time.time()-t0:.0f}s)")
    mc = int(p["min_count"])
    pers = {arm: ST.pairs_persistence([r_[arm]["cells"] for r_ in rows], mc) for arm in ("off", "on")}
    d_pers = np.array([a["delta"] - b["delta"] for a, b in zip(pers["on"], pers["off"])])
    d_loc = np.array([r_["on"]["E_loc"] - r_["off"]["E_loc"] for r_ in rows])
    bs = lambda tag: derive_seed(ctx.seed, "F3.3-boot", mode, tag)
    dp, dl = boot_se(d_pers, B, bs("pers")), boot_se(d_loc, B, bs("loc"))
    absval = {arm: {"pers": boot_se(np.array([x["delta"] for x in pers[arm]]), B, bs("p" + arm)),
                    "loc": boot_se(np.array([r_[arm]["E_loc"] for r_ in rows]), B, bs("l" + arm))} for arm in ("off", "on")}
    dec_p = equivalence(dp[0], dp[1], float(p["margin_pers"]), nse)
    dec_l = equivalence(dl[0], dl[1], float(p["margin_loc"]), nse)
    pos = lambda dec, st: dec == "EXCESS" and st[0] > 0
    decided = lambda dec, st: dec == "EQUIVALENT" or (dec == "EXCESS" and st[0] < 0)
    sat_on = float(np.mean([r_["on"]["sat_fraction"] for r_ in rows]))
    sat_off = float(np.mean([r_["off"]["sat_fraction"] for r_ in rows]))
    no_saturation = sat_on < float(p["sat_min"])
    saturating = [r_["on"]["sat_fraction"] > 0 for r_ in rows]
    invalid = sum(1 for r_, ok in zip(rows, saturating) if ok and not np.isfinite(r_["on"]["E_loc"])) \
        + sum(1 for j, x in enumerate(pers["on"]) if saturating[2 * j] and saturating[2 * j + 1] and not np.isfinite(x["delta"]))
    n_stat = n + len(pers["on"])
    inst = {"off arm = route_capacity": all(r_["off_equals_capacity"] for r_ in rows),
            "balance": max(max(r_["off"]["balance"], r_["on"]["balance"]) for r_ in rows) <= float(p["balance_tol"]),
            "stationary routing converged": all(r_["on"]["converged"] for r_ in rows),
            "invalid statistics": no_saturation or invalid <= float(p["invalid_max"]) * n_stat}
    if not all(inst.values()):
        status = "INSTRUMENT_FAILED"
    elif no_saturation:
        status = "FAILED"
    elif pos(dec_p, dp) and pos(dec_l, dl):
        status = "PASSED"
    elif not pos(dec_p, dp) and not pos(dec_l, dl) and decided(dec_p, dp) and decided(dec_l, dl):
        status = "FAILED"
    else:
        status = "UNDECIDED"
    f = lambda st: f"{st[0]:+.4f} +- {st[1]:.4f}" if np.isfinite(st[0]) else "undefined"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")              # all-nan columns: no saturated cells left after switch-off
        surv = {arm: np.nanmean([r_[arm]["survival"]["frac"] for r_ in rows], axis=0) for arm in ("off", "on")}
        surv_rs = {arm: np.nanmean([r_[arm]["survival"]["rs_with_last_on"] for r_ in rows], axis=0) for arm in ("off", "on")}
    lines = [f"2+1 cylinder T = {p['T']}, N = {p['N']}, s0 = {p['s0']}; analysis |x| < {p['R_abs']}, t in [{p['t1']}, {p['t2']}); "
             f"C = Q_{p['q']} = {C:.4f} ({p['n_cal']} calibration sprinklings); channel '{mode}', Pi_sat = {p['Pi_sat']}, m = {p['m']}; "
             f"{n} sprinklings = {len(d_pers)} pairs",
             f"  saturated fraction (analysis events): off {sat_off:.4f}, on {sat_on:.4f}; Pi of saturated events (median) "
             f"off {np.nanmedian([r_['off']['Pi_sat_median'] for r_ in rows]):.3f}, on "
             f"{np.nanmedian([r_['on']['Pi_sat_median'] for r_ in rows]):.3f}; min Pi on {min(r_['on']['Pi_min'] for r_ in rows):.3f}",
             f"  Delta_pers: off {f(absval['off']['pers'])}, on {f(absval['on']['pers'])}   (F3.2: "
             f"{F32_REFERENCE['Delta_pers'][0]:+.3f} +- {F32_REFERENCE['Delta_pers'][1]:.3f})",
             f"  E_loc:      off {f(absval['off']['loc'])}, on {f(absval['on']['loc'])}   (F3.2: "
             f"{F32_REFERENCE['E_loc'][0]:+.3f} +- {F32_REFERENCE['E_loc'][1]:.3f})",
             f"  paired d_pers = {f(dp)} (margin {p['margin_pers']}) -> {dec_p}; d_loc = {f(dl)} (margin {p['margin_loc']}) -> {dec_l}"
             + ("; no saturation with feedback" if no_saturation else ""),
             "  switch-off at t = %g: saturated fraction in the last slab before it, then in the slabs of %g after it: off "
             % (p["t_off"], p["slab"])
             + " ".join(f"{x:.4f}" for x in surv["off"]) + "; on " + " ".join(f"{x:.4f}" for x in surv["on"]),
             "     r_s of each slab after with the last slab before (nan: no saturated cells left): off "
             + " ".join(f"{x:+.3f}" for x in surv_rs["off"]) + "; on " + " ".join(f"{x:+.3f}" for x in surv_rs["on"])
             + "  (reported)"]
    if mode == "routing":
        its = [it for r_ in rows for it in r_["on"]["iterations"]]
        lines.append(f"  stationary routing passes: median {int(np.median(its))}, max {max(its)}")
    lines.append("instrument: " + "; ".join(f"{k}: {'OK' if v else 'FAILS'}" for k, v in inst.items())
                 + f" (invalid {invalid}/{n_stat})")
    data = {"C": C, "mode": mode, "d_pers": dp, "d_loc": dl, "decisions": {"pers": dec_p, "loc": dec_l},
            "absolute": absval, "sat_fraction": {"off": sat_off, "on": sat_on}, "survival": surv, "survival_rs": surv_rs,
            "pairs": pers, "instrument": inst,
            "seeds": [{arm: {k: v for k, v in r_[arm].items() if k != "cells"} for arm in ("off", "on")} for r_ in rows]}
    return Result(status, lines, data)
