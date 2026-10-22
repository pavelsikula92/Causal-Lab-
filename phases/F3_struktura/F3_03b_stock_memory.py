"""
F3.3b -- Does a carried stock (memory) create a localized, persistent structure?
================================================================================

Why this test exists
--------------------
F3.2 (no feedback) and F3.3 capacity (feedback through C Pi) created no
structure. In both, the excess of an event over its capacity is a local
increment dQ_v that never acts again: the model has no memory, and a
persistent structure (a worldtube) has nothing that carries it forward.
The author decided on 21 Oct 2026 to add memory in the smallest way the
axioms allow: axiom 3 already contains the stock, so the stock is made a
carried state that moves forward along the causal links. Nothing else is
added: no inheritance of Pi, no rigidity I, no feedback on Pi, the
capacity C is constant. F3.2 and F3.3 refute only the variants they
tested; this is a new hypothesis after a failure and changes no earlier
criterion.

Model (core.flux.stock)
-----------------------
At every event, in any linear extension of the causal order:
    offered F~ = s + flux from the parents;  arriving stock Q~ from the parents;
    accepted A = min(F~ + Q~, C) (new flux first, the stock fills the rest);
    held stock Q = F~ + Q~ - A;  A and Q both go to the children uniformly.
Boundary events forward nothing: flux and stock leave the region. Exact
balance: sum s = sum over non-forwarding events of (A + Q). Saturated:
chi_v = 1[Q_v > 0] (the event holds stock). Control arm: the same graphs,
source and capacity, the stock discarded (exactly route_capacity, F3.2).

Lorentz invariance and foliation. The state of an event depends only on
its parents, so every linear extension of the causal order gives the
identical result; the rule has no global time step and no preferred
foliation. Instrument (f) checks it on every sprinkling: the routing
evaluated in the time order of a boosted frame (beta = 0.6 along x1)
must give bit-identical results. (Statistical independence of the
results from the sprinkling is covered by the bootstrap over sprinklings;
this does not make the analysis window, a fixed region in one frame,
frame-independent.)

Two capacity levels, both declared before any run: part q99 (C = Q_0.99,
the F3.2 / F3.3 level) and part q90 (C = Q_0.9, ~5-9 % saturated in the
F3.2 pilot). Reason: at Q_0.99 the excess of an event is small and its
children mostly have spare capacity, so a carried stock may be absorbed
within one link (development run below); the q90 part asks the same
question where the stock is not negligible. Each part is decided on its
own by the same criterion.

Setup and statistics: exactly F3.3 (2+1 cylinder T = 7, N = 56 000, s0 = 1,
R_abs = 0.4, analysis window t in [3.5, 5.5), cells 0.1, slabs 0.25,
min_count 3, n_rbins 16, n_perm 20, C = Q_0.99 of the uncapped offered flux
of 8 independent calibration sprinklings, stream "F3.3b-cal"; q per part).
Per sprinkling and arm: E_loc; cell sequence for Delta_pers (pairs); the
switch-off routing (source only for t < t_off = 5.5): saturated fraction
and mean held stock per event in the slabs after t_off, and the r_s of each
slab's cell pattern with the last slab before t_off (REPORTED: lifetime of
a structure without driving and whether it is carried by the stock).

Criteria (stated before the run; those of F3.3)
-----------------------------------------------
Paired differences on the same causal sets: d_loc = E_loc(on) - E_loc(off)
per sprinkling, d_pers = Delta_pers(on) - Delta_pers(off) per pair; mean
+- 2 SE (bootstrap), margin 0.05: EQUIVALENT / EXCESS / INCONCLUSIVE;
positive = EXCESS with mean > 0.
  PASSED     d_pers and d_loc both positive
  FAILED     neither positive and both decided (EQUIVALENT or negative
             EXCESS); or no saturation with the stock (fraction < 0.001)
  UNDECIDED  otherwise
Instrument (any failing -> INSTRUMENT_FAILED): (a) the control arm equals
route_capacity exactly; (b) the stock arm with keep_stock = False equals it
too (the same code path); (c) balance within 1e-12 on every routing;
(d) undefined statistics <= 5 %; (f) boosted-order evaluation bit-identical
on every sprinkling.

Null model: the control arm (no memory) on the same causal sets.

Seen before this test was written (development, 21 Oct, one sprinkling of
N = 6000, T = 7, C from that sprinkling): saturated fraction of the
analysis events without memory 0.0098, with the carried stock 0.0107
(C = Q_0.99); quick runs (N = 6000, 4 sprinklings): q99 0.0109 vs 0.0111,
q90 0.091 vs 0.116 (paired differences of the quick runs inconclusive);
boosted order bit-identical; balance exact. Nothing in the criteria was
set from it.
"""

from __future__ import annotations

import time
import warnings

import numpy as np

from core.flux.capacity import route_capacity
from core.flux.stock import route_stock, stock_balance
from core.flux import structure as ST
from core.order.links import links_csr
from core.order.sprinkle import Cylinder, derive_seed, sprinkle
from core.stats.parallel import parallel_map
from core.stats.stats import boot_se, equivalence
from core.testing import Result

TEST = {
    "id": "F3.3b",
    "phase": "F3",
    "title": "Vytváří přenášená zásoba (paměť) lokalizovanou trvalou strukturu?",
    "theory_part": "axiom 3 (zásoba v bilanci toku) jako přenášený stav; hmota jako trvalá struktura",
    "question": "Vytvoří zásoba, která se přenáší po kauzálních linkách a znovu nabízí následníkům, lokalizovanou a trvalou "
                "strukturu saturovaných událostí nad úrovní téhož modelu bez paměti?",
    "method": "geometrie a statistiky F3.3 (2+1 válec, N = 56 000, C = Q_0,99); na stejných kauzálních množinách směrování se "
              "zásobou přenášenou k potomkům a kontrola se zásobou zahozenou (F3.2); párové rozdíly Δ_pers a E_loc; "
              "přežití a zásoba po vypnutí zdroje (hlášeno); vyhodnocení v pořadí boostované soustavy",
    "prediction": "bez předpovědi směru; rozhoduje párový rozdíl",
    "criterion": "MĚŘIDLO: kontrola = route_capacity, bilance <= 1e-12, neplatné statistiky <= 5 %, pořadí boostované soustavy "
                 "dá bitově stejný výsledek. PROŠEL: d_pers i d_loc kladně EXCESS (mez 0,05, 2 SE). NEPROŠEL: žádný kladný a "
                 "oba rozhodnuté, nebo saturace zmizí (< 0,001). NEROZHODNUTO: jinak.",
    "null_model": "stejné kauzální množiny bez paměti (zásoba zahozena, F3.2)",
    "level": "mechanism",
    "meters": ["K.9"],
    "depends_on": ["F3.2", "F3.3"],
    "params": {"D": 3, "T": 7.0, "R": 0.5, "R_abs": 0.4, "s0": 1.0, "t1": 3.5, "t2": 5.5, "t_off": 5.5, "post_end": 6.75,
               "cell": 0.1, "slab": 0.25, "min_count": 3, "n_rbins": 16, "n_perm": 20, "q": 0.99, "N": 56000,
               "n_seeds": 128, "n_cal": 8, "margin_pers": 0.05, "margin_loc": 0.05, "n_se": 2.0, "boot_B": 4000,
               "sat_min": 0.001, "balance_tol": 1e-12, "invalid_max": 0.05, "boost_beta": 0.6},
    "parts": {"q99": {"q": 0.99}, "q90": {"q": 0.9}},
    "quick": {"N": 6000, "n_seeds": 4, "n_cal": 1, "n_perm": 4, "boot_B": 200},
    "estimate": "každá část ~20–40 min s 8 procesy (5 směrování a 2 lokalizace na sprinkling)",
}


def geometry(p, seed):
    """One cylinder sprinkling, its links, radii, times and the masks of F3.2."""
    S = sprinkle(Cylinder(int(p["D"]), 0.0, float(p["T"]), float(p["R"])), int(p["N"]), seed)
    ptr, idx = links_csr(S.X, S.rho)
    r, t = S.r(), S.X[:, 0]
    absorb = r > float(p["R_abs"])
    win = (t >= float(p["t1"])) & (t < float(p["t2"])) & (r < float(p["R_abs"]))
    return S, ptr, idx, r, t, absorb, win


def calibration_task(p, seed):
    """Uncapped offered flux of the analysis events of one calibration sprinkling."""
    S, ptr, idx, r, t, absorb, win = geometry(p, seed)
    return route_capacity(ptr, idx, np.full(S.n, float(p["s0"])), np.inf, absorb)["offered"][win]


def survival(S, r, t, res, p):
    """After switch-off: saturated fraction and mean held stock per slab, r_s of each slab with the last one before."""
    pp = dict(p, t1=float(p["t_off"]) - float(p["slab"]), t2=float(p["post_end"]))
    reg = (r < float(p["R_abs"])) & (t >= pp["t1"]) & (t < pp["t2"])
    chi = (res["stock"] > 0).astype(float)
    seq = ST.cells(S.X, chi, reg, pp)
    edges = np.arange(pp["t1"], pp["t2"] + 1e-9, float(p["slab"]))
    sl = [reg & (t >= edges[k]) & (t < edges[k + 1]) for k in range(len(edges) - 1)]
    return {"frac": [float(chi[m].mean()) for m in sl], "stock": [float(res["stock"][m].mean()) for m in sl],
            "rs_with_last_on": [ST.rs_successive([seq[0], seq[k]], int(p["min_count"])) for k in range(1, len(seq))]}


def seed_task(p, seed, C):
    """Both arms on one sprinkling: statistics, switch-off survival, and the boosted-order check."""
    S, ptr, idx, r, t, absorb, win = geometry(p, seed)
    s_on = np.full(S.n, float(p["s0"]))
    s_off = np.where(t < float(p["t_off"]), float(p["s0"]), 0.0)
    ref = route_capacity(ptr, idx, s_on, C, absorb)
    out = {}
    for arm, keep in (("off", False), ("on", True)):
        main = route_stock(ptr, idx, s_on, C, absorb, keep_stock=keep)
        swo = route_stock(ptr, idx, s_off, C, absorb, keep_stock=keep)
        chi = (main["stock"] > 0).astype(float)
        loc = ST.localization(S.X, r, chi, win, p, np.random.default_rng(derive_seed(seed, "F3.3b-perm", arm)))
        out[arm] = {"sat_fraction": float(chi[win].mean()), "E_loc": loc["E_loc"], "I_obs": loc["I_obs"],
                    "cells": ST.cells(S.X, chi, win, p), "stock_mean": float(main["stock"][win].mean()),
                    "balance": max(stock_balance(s_on, main, keep), stock_balance(s_off, swo, keep)),
                    "survival": survival(S, r, t, swo, p)}
        if arm == "off":
            out["off_equals_capacity"] = (np.array_equal(main["accepted"], ref["accepted"])
                                          and np.array_equal(main["stock"], ref["stored"])
                                          and np.array_equal(main["Phi"], ref["Phi"]))
        else:
            beta = float(p["boost_beta"])
            tb = (t - beta * S.X[:, 1]) / np.sqrt(1.0 - beta * beta)
            boosted = route_stock(ptr, idx, s_on, C, absorb, keep_stock=True, order=np.argsort(tb, kind="stable"))
            out["boost_identical"] = all(np.array_equal(main[k], boosted[k]) for k in ("accepted", "stock", "Phi", "Psi"))
    return out


def run(p, ctx):
    n, B, nse = int(p["n_seeds"]), int(p["boot_B"]), float(p["n_se"])
    t0 = time.time()
    cal = parallel_map(calibration_task, [(p, derive_seed(ctx.seed, "F3.3b-cal", int(p["N"]), k)) for k in range(int(p["n_cal"]))],
                       ctx.workers, ctx.progress, lambda i, r_, d, nn: f"[F3.3b] calibration {d}/{nn} ({time.time()-t0:.0f}s)")
    C = float(np.quantile(np.concatenate(cal), float(p["q"])))
    seeds = [derive_seed(ctx.seed, "F3.3b", int(p["N"]), k) for k in range(n)]
    rows = parallel_map(seed_task, [(p, sd, C) for sd in seeds], ctx.workers, ctx.progress,
                        lambda i, r_, d, nn: f"[F3.3b] seed {d}/{nn}: saturated off {r_['off']['sat_fraction']:.4f} "
                                             f"on {r_['on']['sat_fraction']:.4f}; E_loc off {r_['off']['E_loc']:+.4f} on "
                                             f"{r_['on']['E_loc']:+.4f} ({time.time()-t0:.0f}s)")
    mc = int(p["min_count"])
    pers = {arm: ST.pairs_persistence([r_[arm]["cells"] for r_ in rows], mc) for arm in ("off", "on")}
    d_pers = np.array([a["delta"] - b["delta"] for a, b in zip(pers["on"], pers["off"])])
    d_loc = np.array([r_["on"]["E_loc"] - r_["off"]["E_loc"] for r_ in rows])
    bs = lambda tag: derive_seed(ctx.seed, "F3.3b-boot", tag)
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
    inst = {"control arm = route_capacity": all(r_["off_equals_capacity"] for r_ in rows),
            "balance": max(max(r_["off"]["balance"], r_["on"]["balance"]) for r_ in rows) <= float(p["balance_tol"]),
            "invalid statistics": no_saturation or invalid <= float(p["invalid_max"]) * n_stat,
            "boosted order bit-identical": all(r_["boost_identical"] for r_ in rows)}
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
        surv = {arm: {k: np.nanmean([r_[arm]["survival"][k] for r_ in rows], axis=0) for k in ("frac", "stock", "rs_with_last_on")}
                for arm in ("off", "on")}
    lines = [f"2+1 cylinder T = {p['T']}, N = {p['N']}, s0 = {p['s0']}; analysis |x| < {p['R_abs']}, t in [{p['t1']}, {p['t2']}); "
             f"C = Q_{p['q']} = {C:.4f} ({p['n_cal']} calibration sprinklings); carried stock vs discarded; "
             f"{n} sprinklings = {len(d_pers)} pairs",
             f"  saturated (holding stock) fraction of the analysis events: no memory {sat_off:.4f}, carried stock {sat_on:.4f}; "
             f"mean held stock per event {np.mean([r_['off']['stock_mean'] for r_ in rows]):.4f} vs "
             f"{np.mean([r_['on']['stock_mean'] for r_ in rows]):.4f}",
             f"  Delta_pers: no memory {f(absval['off']['pers'])}, carried {f(absval['on']['pers'])}   (F3.2: +0.006 +- 0.006)",
             f"  E_loc:      no memory {f(absval['off']['loc'])}, carried {f(absval['on']['loc'])}   (F3.2: +0.004 +- 0.004)",
             f"  paired d_pers = {f(dp)} (margin {p['margin_pers']}) -> {dec_p}; d_loc = {f(dl)} (margin {p['margin_loc']}) -> {dec_l}"
             + ("; no saturation with the stock" if no_saturation else ""),
             "  switch-off at t = %g; the last slab before it, then the slabs of %g after it:" % (p["t_off"], p["slab"]),
             "     saturated fraction  no memory " + " ".join(f"{x:.4f}" for x in surv["off"]["frac"])
             + "; carried " + " ".join(f"{x:.4f}" for x in surv["on"]["frac"]),
             "     mean held stock     no memory " + " ".join(f"{x:.4f}" for x in surv["off"]["stock"])
             + "; carried " + " ".join(f"{x:.4f}" for x in surv["on"]["stock"]),
             "     r_s with the last slab before (nan: no saturated cells left): no memory "
             + " ".join(f"{x:+.3f}" for x in surv["off"]["rs_with_last_on"]) + "; carried "
             + " ".join(f"{x:+.3f}" for x in surv["on"]["rs_with_last_on"]) + "  (reported)",
             "instrument: " + "; ".join(f"{k}: {'OK' if v else 'FAILS'}" for k, v in inst.items()) + f" (invalid {invalid}/{n_stat})"]
    data = {"C": C, "d_pers": dp, "d_loc": dl, "decisions": {"pers": dec_p, "loc": dec_l}, "absolute": absval,
            "sat_fraction": {"off": sat_off, "on": sat_on}, "survival": surv, "pairs": pers, "instrument": inst,
            "seeds": [{arm: {k: v for k, v in r_[arm].items() if k != "cells"} for arm in ("off", "on")} for r_ in rows]}
    return Result(status, lines, data)
