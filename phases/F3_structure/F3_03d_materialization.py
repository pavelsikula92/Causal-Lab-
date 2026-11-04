"""
F3.3d -- Does a materialization threshold turn the copied pulse into a localized node?
=====================================================================================

Why this test exists
--------------------
F3.3c main_q99: the copied pulse (main-parent rule, alive wherever an event
is saturated, acting as a source) is memory that does not dilute and keeps
the position (after switch-off the saturated fraction decays 44 % -> 8 % ->
1.6 % -> 0.05 % over three slabs, the pattern correlated with the one before,
r_s +0.66, +0.41, +0.12), but while driven it floods the region (30 %
saturated) evenly: no localization, no persistence beyond the envelope
(paired differences EQUIVALENT). The author's hypothesis (21 Oct 2026):
materialization is a threshold. The pulse stays bound (copied in full
along the main lineage) only where the excess load reaches a critical value
dQ >= dQ_crit; below it the pulse is lost. Between the subcritical regime
(the pulse fades, F3.3b) and the flood without bounds (F3.3c) a bound,
localized state is expected.

The second half of the hypothesis (self-focusing: a slower clock narrows
the causal cone of the children and pulls the surrounding links into the
core) cannot be tested on a FIXED causal set: narrowing J+(v) changes which
pairs are related, i.e. the dynamics of the order (axiom 7). It is not part
of this test.

Model: F3.3c with one change (core.flux.pulse.route_pulse, dq_crit): the
pulse is alive at v only if dQ_v >= dQ_crit (dQ_v = offered - accepted, in
the units of the driving s0 = omega = 1); alive: v keeps the inherited
pulse or creates omega; not alive: P = 0 at once (no decay parameter is
added). Rule 'main' (the author's choice for the lineage), C = Q_0.99.
dQ_crit = 0 is F3.3c main_q99.

Scan, declared before the run: dQ_crit = 0.1, 0.25, 0.5, 1.0 (one part
each). The roadmap gate of F3 asks for a structure in a pre-declared wide
range, not at one tuned value; every part is decided on its own and the
scan is reported as a whole. Range set from the distribution of dQ of the
saturated events with the pulse in a development run (one sprinkling,
N = 12 000, dQ_crit = 0: quantiles 10/50/90/99 % = 0.03, 0.36, 0.95, 1.48;
without the pulse 0.008, 0.09, 0.48, 0.78).

Setup, statistics and criteria: exactly F3.3c (2+1 cylinder T = 7,
N = 56 000, analysis window t in [3.5, 5.5), |x| < 0.4; paired differences
pulse - no pulse on the same causal sets; margin 0.05, 2 SE; PASSED both
positive EXCESS, FAILED neither positive and both decided or no
saturation, UNDECIDED otherwise; switch-off reported; boosted-order
identity on every sprinkling).

Seen before this test was written: F3.3c main_q99 (above), the
development run named in the scan paragraph, and the quick runs (N = 6000,
4 sprinklings): crit025 saturated 0.0185 vs 0.0095 without the pulse,
pulse on 1.0 % of the events; crit100 pulse on 0.02 %, saturation as
without it. Nothing in the criteria was set from them.
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
    "id": "F3.3d",
    "phase": "F3",
    "title": "Does a materialisation threshold turn a copied pulse into a localised node?",
    "theory_part": "materialisation as a threshold: a bound state of the pulse only where the excess load reaches "
                   "ΔQ_crit",
    "question": "Does a copied pulse (from the main parent) that survives only where the excess load reaches ΔQ_crit "
                "create a localised and persistent structure above the level of the same model without a pulse?",
    "method": "F3.3c main_q99 with a single change: the pulse lives only at ΔQ >= ΔQ_crit; a pre-set scan ΔQ_crit = "
              "0.1; 0.25; 0.5; 1.0; paired differences Δ_pers and E_loc against the branch without a pulse on the "
              "same causal sets",
    "prediction": "according to the author's hypothesis a bound localised state between extinction and flooding; the "
                  "paired difference decides",
    "criterion": "INSTRUMENT: control = route_capacity, flow balance <= 1e-12, invalid statistics <= 5 %, the order "
                 "of a boosted frame bitwise identical. PASSED: d_pers and d_loc both positively EXCESS (margin "
                 "0.05, 2 SE). FAILED: none positive and both decided, or the saturation disappears (< 0.001). "
                 "UNDECIDED: otherwise. Each value of ΔQ_crit separately.",
    "null_model": "the same causal sets without a pulse (F3.2)",
    "level": "mechanism",
    "meters": ["K.9"],
    "depends_on": ["F3.3c"],
    "params": {"D": 3, "T": 7.0, "R": 0.5, "R_abs": 0.4, "s0": 1.0, "t1": 3.5, "t2": 5.5, "t_off": 5.5, "post_end": 6.75,
               "cell": 0.1, "slab": 0.25, "min_count": 3, "n_rbins": 16, "n_perm": 20, "q": 0.99, "N": 56000,
               "n_seeds": 128, "n_cal": 8, "omega": 1.0, "rule": "main", "margin_pers": 0.05, "margin_loc": 0.05,
               "n_se": 2.0, "boot_B": 4000, "sat_min": 0.001, "balance_tol": 1e-12, "invalid_max": 0.05, "boost_beta": 0.6},
    "parts": {"crit010": {"dq_crit": 0.1}, "crit025": {"dq_crit": 0.25}, "crit050": {"dq_crit": 0.5},
              "crit100": {"dq_crit": 1.0}},
    "quick": {"N": 6000, "n_seeds": 4, "n_cal": 1, "n_perm": 4, "boot_B": 200},
    "estimate": "each part ~20–40 min with 8–10 processes",
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
    dqc = float(p["dq_crit"])
    out = {}
    for arm, rule in (("off", "none"), ("on", str(p["rule"]))):
        main = route_pulse(ptr, idx, s_on, C, absorb, rule, om, dq_crit=dqc)
        swo = route_pulse(ptr, idx, s_off, C, absorb, rule, om, dq_crit=dqc)
        chi = (main["stored"] > 0).astype(float)
        loc = ST.localization(S.X, r, chi, win, p, np.random.default_rng(derive_seed(seed, "F3.3d-perm", arm)))
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
            boosted = route_pulse(ptr, idx, s_on, C, absorb, rule, om, order=np.argsort(tb, kind="stable"), dq_crit=dqc)
            out["boost_identical"] = all(np.array_equal(main[k], boosted[k]) for k in ("accepted", "stored", "Phi", "P"))
    return out


def run(p, ctx):
    n, B, nse = int(p["n_seeds"]), int(p["boot_B"]), float(p["n_se"])
    t0 = time.time()
    cal = parallel_map(calibration_task, [(p, derive_seed(ctx.seed, "F3.3d-cal", int(p["N"]), k)) for k in range(int(p["n_cal"]))],
                       ctx.workers, ctx.progress, lambda i, r_, d, nn: f"[F3.3d] calibration {d}/{nn} ({time.time()-t0:.0f}s)")
    C = float(np.quantile(np.concatenate(cal), float(p["q"])))
    seeds = [derive_seed(ctx.seed, "F3.3d", int(p["N"]), k) for k in range(n)]
    rows = parallel_map(seed_task, [(p, sd, C) for sd in seeds], ctx.workers, ctx.progress,
                        lambda i, r_, d, nn: f"[F3.3d] seed {d}/{nn}: saturated off {r_['off']['sat_fraction']:.4f} "
                                             f"on {r_['on']['sat_fraction']:.4f}; E_loc off {r_['off']['E_loc']:+.4f} on "
                                             f"{r_['on']['E_loc']:+.4f} ({time.time()-t0:.0f}s)")
    mc = int(p["min_count"])
    pers = {arm: ST.pairs_persistence([r_[arm]["cells"] for r_ in rows], mc) for arm in ("off", "on")}
    d_pers = np.array([a["delta"] - b["delta"] for a, b in zip(pers["on"], pers["off"])])
    d_loc = np.array([r_["on"]["E_loc"] - r_["off"]["E_loc"] for r_ in rows])
    bs = lambda tag: derive_seed(ctx.seed, "F3.3d-boot", tag)
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
             f"C = Q_{p['q']} = {C:.4f} ({p['n_cal']} calibration sprinklings); copied pulse, rule '{p['rule']}', omega = {p['omega']}, dQ_crit = {p['dq_crit']}; "
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
