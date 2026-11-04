"""
F3.3c -- Does a copied (undiluted) pulse create a localized, persistent structure?
==================================================================================

Why this test exists
--------------------
F3.3b: a stock carried to the children is real memory but dilutes; it is
split over the children like the flux and spreads over the future cone
within a step or two (paired differences EQUIVALENT to zero with SE ~0.0015).
The author's hypothesis (21 Oct 2026): the state of a persistent structure
(a particle, a node) is not an extensive quantity divided along the links
(sum out = in) but an intensive state of the event, like a frequency or a
phase, which every child receives in full ("the causal pulse is not
diluted"). Two ways of combining the pulses of several parents, both
proposed by the author, are tested as separate parts:
  max   P_in(v) = max over the parents of P(u)       (a dominant excitation)
  main  P_in(v) = P(u*) of the parent that sends v the most flux (the pulse
        follows the carrying flux: a worldline-like lineage)
The flux-weighted mean was rejected by the author before the run (it is a
discrete diffusion, the dilution again).

Choices made to keep the test from being a tautology (stated before the run)
-----------------------------------------------------------------------------
  * The pulse has an effect: it acts as internal events, s_eff = s0 + P_in.
    A label copied forward would persist by construction; here structure is
    judged on the FLUX (the saturated events), as in F3.2-F3.3b.
  * The pulse lives only at the capacity: a saturated event (dQ > 0) keeps
    the inherited pulse or creates one of value omega = s0 = 1; an
    unsaturated event has P = 0. Without such a limit 'max' fills the whole
    future cone.
  * Nothing else changes: constant capacity C, the stock is local (not
    carried), no Pi feedback, no rigidity. One change against F3.2.
omega = 1 is an experimental parameter (a later scan belongs to F3.4).

Model: core.flux.pulse.route_pulse. The state of an event depends only on
its parents, so every linear extension of the causal order gives the same
result (checked on every sprinkling in a boosted-frame order, beta = 0.6).
The flux balance of axiom 3 is exact (the pulse injection is a source).

Setup, statistics and criteria: exactly F3.3b (2+1 cylinder T = 7,
N = 56 000, analysis window t in [3.5, 5.5), |x| < 0.4; Delta_pers per pair,
E_loc per sprinkling; paired differences pulse - no pulse on the same
causal sets; margin 0.05, 2 SE; PASSED both positive EXCESS, FAILED neither
positive and both decided or no saturation, UNDECIDED otherwise; switch-off
at t = 5.5 reported, now also the fraction of events carrying a pulse).
Capacity levels as F3.3b: q99 (C = Q_0.99) and q90 (C = Q_0.9), declared
before the run. Four parts: max_q99, main_q99, max_q90, main_q90.
A flood is not a structure: if the pulse saturates most of the region,
E_loc (permutation within radial bins) shows no localization; the
saturated and pulse-carrying fractions are reported.

Seen before this test was written (development, 21 Oct, one sprinkling of
N = 6000, C from that sprinkling): saturated (= pulse-carrying) fraction of
the analysis events: q99 max 0.56, main 0.038; q90 max 0.77, main 0.61
(without pulse ~0.01 and ~0.09); after switch-off saturated events in the
analysis region in four slabs: q99 max 27, 17, 5, 1, main 1, 0, 0, 0;
q90 max 48, 27, 15, 6, main 28, 8, 2, 0. Nothing in the criteria was set
from it.
"""

from __future__ import annotations

import time
import warnings

import numpy as np

from core.flux.capacity import route_capacity
from core.flux.pulse import route_pulse, pulse_balance
from core.flux import structure as ST
from core.order.links import links_csr
from core.order.sprinkle import Cylinder, derive_seed, sprinkle
from core.stats.parallel import parallel_map
from core.stats.stats import boot_se, equivalence
from core.testing import Result

TEST = {
    "id": "F3.3c",
    "phase": "F3",
    "title": "Does a copied (undiluted) pulse create a localised persistent structure?",
    "theory_part": "the state of a persistent structure as an intensive quantity (frequency, phase) that is copied, "
                   "not divided",
    "question": "Does a pulse that every descendant receives whole (maximum over the parents, or from the main "
                "parent), that lives only in saturated events and acts as a source, create a localised and "
                "persistent structure of saturated events above the level of the same model without a pulse?",
    "method": "geometry and statistics of F3.3b (2+1 cylinder, N = 56 000); on the same causal sets routing with a "
              "pulse (ω = 1, source s0 + P) and without a pulse (F3.2); rules max and main, capacities Q_0.99 and "
              "Q_0.9; paired differences Δ_pers and E_loc; survival after the source is switched off (reported); "
              "evaluation in the order of a boosted frame",
    "prediction": "no prediction of the direction; the paired difference decides (flooding of the region does not "
                  "count as localisation)",
    "criterion": "INSTRUMENT: control = route_capacity, flow balance <= 1e-12, invalid statistics <= 5 %, the order "
                 "of a boosted frame bitwise identical. PASSED: d_pers and d_loc both positively EXCESS (margin "
                 "0.05, 2 SE). FAILED: none positive and both decided, or the saturation disappears (< 0.001). "
                 "UNDECIDED: otherwise.",
    "null_model": "the same causal sets without a pulse (F3.2)",
    "level": "mechanism",
    "meters": ["K.9"],
    "depends_on": ["F3.2", "F3.3b"],
    "params": {"D": 3, "T": 7.0, "R": 0.5, "R_abs": 0.4, "s0": 1.0, "t1": 3.5, "t2": 5.5, "t_off": 5.5, "post_end": 6.75,
               "cell": 0.1, "slab": 0.25, "min_count": 3, "n_rbins": 16, "n_perm": 20, "q": 0.99, "N": 56000,
               "n_seeds": 128, "n_cal": 8, "omega": 1.0, "margin_pers": 0.05, "margin_loc": 0.05, "n_se": 2.0,
               "boot_B": 4000, "sat_min": 0.001, "balance_tol": 1e-12, "invalid_max": 0.05, "boost_beta": 0.6},
    "parts": {"max_q99": {"rule": "max", "q": 0.99}, "main_q99": {"rule": "main", "q": 0.99},
              "max_q90": {"rule": "max", "q": 0.9}, "main_q90": {"rule": "main", "q": 0.9}},
    "quick": {"N": 6000, "n_seeds": 4, "n_cal": 1, "n_perm": 4, "boot_B": 200},
    "estimate": "each part ~20–40 min with 8 processes (5 routings and 2 localisations per sprinkling)",
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
    """After switch-off: saturated and pulse-carrying fraction per slab, r_s of each slab with the last one before."""
    pp = dict(p, t1=float(p["t_off"]) - float(p["slab"]), t2=float(p["post_end"]))
    reg = (r < float(p["R_abs"])) & (t >= pp["t1"]) & (t < pp["t2"])
    chi = (res["stored"] > 0).astype(float)
    seq = ST.cells(S.X, chi, reg, pp)
    edges = np.arange(pp["t1"], pp["t2"] + 1e-9, float(p["slab"]))
    sl = [reg & (t >= edges[k]) & (t < edges[k + 1]) for k in range(len(edges) - 1)]
    return {"frac": [float(chi[m].mean()) for m in sl], "pulse": [float((res["P"][m] > 0).mean()) for m in sl],
            "rs_with_last_on": [ST.rs_successive([seq[0], seq[k]], int(p["min_count"])) for k in range(1, len(seq))]}


def seed_task(p, seed, C):
    """Both arms on one sprinkling: statistics, switch-off survival, and the boosted-order check."""
    S, ptr, idx, r, t, absorb, win = geometry(p, seed)
    s_on = np.full(S.n, float(p["s0"]))
    s_off = np.where(t < float(p["t_off"]), float(p["s0"]), 0.0)
    ref = route_capacity(ptr, idx, s_on, C, absorb)
    om = float(p["omega"])
    out = {}
    for arm, rule in (("off", "none"), ("on", str(p["rule"]))):
        main = route_pulse(ptr, idx, s_on, C, absorb, rule, om)
        swo = route_pulse(ptr, idx, s_off, C, absorb, rule, om)
        chi = (main["stored"] > 0).astype(float)
        loc = ST.localization(S.X, r, chi, win, p, np.random.default_rng(derive_seed(seed, "F3.3c-perm", arm)))
        out[arm] = {"sat_fraction": float(chi[win].mean()), "E_loc": loc["E_loc"], "I_obs": loc["I_obs"],
                    "cells": ST.cells(S.X, chi, win, p), "pulse_fraction": float((main["P"][win] > 0).mean()),
                    "balance": max(pulse_balance(main), pulse_balance(swo)),
                    "survival": survival(S, r, t, swo, p)}
        if arm == "off":
            out["off_equals_capacity"] = (np.array_equal(main["accepted"], ref["accepted"])
                                          and np.array_equal(main["stored"], ref["stored"])
                                          and np.array_equal(main["Phi"], ref["Phi"]))
        else:
            beta = float(p["boost_beta"])
            tb = (t - beta * S.X[:, 1]) / np.sqrt(1.0 - beta * beta)
            boosted = route_pulse(ptr, idx, s_on, C, absorb, rule, om, order=np.argsort(tb, kind="stable"))
            out["boost_identical"] = all(np.array_equal(main[k], boosted[k]) for k in ("accepted", "stored", "Phi", "P"))
    return out


def run(p, ctx):
    n, B, nse = int(p["n_seeds"]), int(p["boot_B"]), float(p["n_se"])
    t0 = time.time()
    cal = parallel_map(calibration_task, [(p, derive_seed(ctx.seed, "F3.3c-cal", int(p["N"]), k)) for k in range(int(p["n_cal"]))],
                       ctx.workers, ctx.progress, lambda i, r_, d, nn: f"[F3.3c] calibration {d}/{nn} ({time.time()-t0:.0f}s)")
    C = float(np.quantile(np.concatenate(cal), float(p["q"])))
    seeds = [derive_seed(ctx.seed, "F3.3c", int(p["N"]), k) for k in range(n)]
    rows = parallel_map(seed_task, [(p, sd, C) for sd in seeds], ctx.workers, ctx.progress,
                        lambda i, r_, d, nn: f"[F3.3c] seed {d}/{nn}: saturated off {r_['off']['sat_fraction']:.4f} "
                                             f"on {r_['on']['sat_fraction']:.4f}; E_loc off {r_['off']['E_loc']:+.4f} on "
                                             f"{r_['on']['E_loc']:+.4f} ({time.time()-t0:.0f}s)")
    mc = int(p["min_count"])
    pers = {arm: ST.pairs_persistence([r_[arm]["cells"] for r_ in rows], mc) for arm in ("off", "on")}
    d_pers = np.array([a["delta"] - b["delta"] for a, b in zip(pers["on"], pers["off"])])
    d_loc = np.array([r_["on"]["E_loc"] - r_["off"]["E_loc"] for r_ in rows])
    bs = lambda tag: derive_seed(ctx.seed, "F3.3c-boot", tag)
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
        surv = {arm: {k: np.nanmean([r_[arm]["survival"][k] for r_ in rows], axis=0) for k in ("frac", "pulse", "rs_with_last_on")}
                for arm in ("off", "on")}
    lines = [f"2+1 cylinder T = {p['T']}, N = {p['N']}, s0 = {p['s0']}; analysis |x| < {p['R_abs']}, t in [{p['t1']}, {p['t2']}); "
             f"C = Q_{p['q']} = {C:.4f} ({p['n_cal']} calibration sprinklings); copied pulse, rule '{p['rule']}', omega = {p['omega']}; "
             f"{n} sprinklings = {len(d_pers)} pairs",
             f"  saturated fraction of the analysis events: no pulse {sat_off:.4f}, pulse {sat_on:.4f}; events carrying a "
             f"pulse {np.mean([r_['on']['pulse_fraction'] for r_ in rows]):.4f}",
             f"  Delta_pers: no pulse {f(absval['off']['pers'])}, pulse {f(absval['on']['pers'])}   (F3.2: +0.006 +- 0.006)",
             f"  E_loc:      no pulse {f(absval['off']['loc'])}, pulse {f(absval['on']['loc'])}   (F3.2: +0.004 +- 0.004)",
             f"  paired d_pers = {f(dp)} (margin {p['margin_pers']}) -> {dec_p}; d_loc = {f(dl)} (margin {p['margin_loc']}) -> {dec_l}"
             + ("; no saturation with the pulse" if no_saturation else ""),
             "  switch-off at t = %g; the last slab before it, then the slabs of %g after it:" % (p["t_off"], p["slab"]),
             "     saturated fraction  no pulse " + " ".join(f"{x:.4f}" for x in surv["off"]["frac"])
             + "; pulse " + " ".join(f"{x:.4f}" for x in surv["on"]["frac"]),
             "     carrying a pulse    " + " ".join(f"{x:.4f}" for x in surv["on"]["pulse"]),
             "     r_s with the last slab before (nan: no saturated cells left): no pulse "
             + " ".join(f"{x:+.3f}" for x in surv["off"]["rs_with_last_on"]) + "; pulse "
             + " ".join(f"{x:+.3f}" for x in surv["on"]["rs_with_last_on"]) + "  (reported)",
             "instrument: " + "; ".join(f"{k}: {'OK' if v else 'FAILS'}" for k, v in inst.items()) + f" (invalid {invalid}/{n_stat})"]
    data = {"C": C, "d_pers": dp, "d_loc": dl, "decisions": {"pers": dec_p, "loc": dec_l}, "absolute": absval,
            "sat_fraction": {"off": sat_off, "on": sat_on}, "survival": surv, "pairs": pers, "instrument": inst,
            "seeds": [{arm: {k: v for k, v in r_[arm].items() if k != "cells"} for arm in ("off", "on")} for r_ in rows]}
    return Result(status, lines, data)
