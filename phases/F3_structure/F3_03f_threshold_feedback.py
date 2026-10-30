"""
F3.3f -- Does the clock-rate feedback turn the thresholded pulse into a localized, persistent node?
==================================================================================================

Why this test exists
--------------------
Step D of the author's F3.3 plan (21 Oct 2026): the best partial mechanism
combined with ONE feedback channel. No channel of F3.3 won on its own
(capacity C Pi floods 54 %, source s Pi and routing Pi_k / I_k remove the
saturation). The best partial mechanism of the series is the copied pulse
with the materialization threshold (F3.3d, F3.3e): weak but real
localization, no persistence, memory gone a slab or two after switch-off;
F3.3e confirm at dQ_crit = 0.2 (new sprinklings) d_loc = +0.0095 +- 0.0043,
d_pers = -0.009 +- 0.007, FAILED. The author's hypothesis (22 Oct 2026):
once the flux above the threshold slows the local clock (Pi < 1), the
slowing holds the load in place and the node sustains itself.

Model: core.flux.pulse.route_pulse with rule 'main', omega = 1,
dq_crit = 0.2 (the value chosen from the F3.3e search: CHOSEN FROM DATA,
stated in every record) and one channel of the clock rate (axiom 5
response Pi = 1 - (1 - Pi_sat) min(sigma, 1)^m, Pi_sat = 0.5, m = 2, the
author's F3.3 values), one part per channel:
  capacity  C_v = C Pi_v (slow clocks process less; positive feedback, the
            channel that flooded alone in F3.3);
  source    s_eff = (s + P_in) Pi_v (internal events and the pulse run at
            the clock rate; negative feedback, removed the saturation alone).
Routing (Pi_k / I_k) is not combined: it is a stationary state, not a
causal pass, and the pulse is a causal pass.

What is not tested: the second half of the author's picture, "the other
causal lines bend into the well". On a fixed causal set the links do not
move; only chains (motion, F2) respond to Pi. Pulling links into the node
is the dynamics of the order (axiom 7).

Arms, all on the same sprinklings:
  base       no pulse, no channel (F3.2; route_capacity, the null of the series)
  threshold  pulse with dQ_crit, no channel (F3.3e confirm's model)
  combined   pulse with dQ_crit and the channel
Decision on combined - base (the criterion of the whole F3 series:
structure above the level without the mechanism). Reported: combined -
threshold (what the channel adds), switch-off survival of all arms, Pi of
the saturated events.

Criteria exactly F3.3d / F3.3e confirm: paired d_pers, d_loc, margin 0.05,
2 SE. PASSED both positive EXCESS; FAILED neither positive and both
decided, or no saturation (< 0.001); UNDECIDED otherwise. Instrument: base
arm = route_capacity bit for bit, flux balance <= 1e-12 in every arm,
invalid statistics <= 5 %, boosted-order identity of the combined arm.

Setup: 2+1 cylinder T = 7, N = 56 000, window t in [3.5, 5.5), |x| < 0.4,
C = Q_0.99 from 8 calibration sprinklings of this test's own stream,
128 sprinklings of its own stream (independent of F3.3d and F3.3e).

Seen before this test was written: F3.3, F3.3c, F3.3d, F3.3e (numbers
above) and the core check on a small sprinkling (1 500 events: the
channels with rule 'none' equal core.flux.feedback.route_feedback, balance
and order-independence with the pulse). No run of the combined model on
the test geometry; nothing in the criteria was set from data.
Quick runs made after the criteria were written (N = 6 000, 8 sprinklings,
decide nothing): capacity: combined arm saturated 98 % of the last slab
before switch-off and the pulse kept 52 %, 42 %, 24 %, 10 %, 8 % of the
slabs after it (threshold arm alone 2.6 %, 0.7 %, 0); source: no saturated
event in the combined arm.
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
    "id": "F3.3f",
    "phase": "F3",
    "title": "Does a clock (Π) feedback turn a thresholded pulse into a localised persistent node?",
    "theory_part": "step D of the F3.3 plan: materialisation threshold + one channel of Π feedback (axiom 5)",
    "question": "Does a copied pulse with the threshold ΔQ_crit = 0.2 together with one channel of clock feedback "
                "(capacity C·Π, or source s·Π) create a localised and persistent structure above the level of F3.2?",
    "method": "route_pulse (main, ω = 1, ΔQ_crit = 0.2 selected from the data of F3.3e) with a Π channel (Π_sat = "
              "0.5, m = 2), one part per channel; three branches on the same causal sets (no pulse, threshold only, "
              "threshold + channel); threshold + channel against the branch without a pulse decides, the "
              "contribution of the channel (against the branch with the threshold only) is reported",
    "prediction": "according to the author's hypothesis the slowing of the clocks above the threshold keeps the load "
                  "and the node sustains itself; the paired difference decides",
    "criterion": "INSTRUMENT: branch without a pulse = route_capacity, flow balance <= 1e-12 in all branches, "
                 "invalid statistics <= 5 %, the order of a boosted frame bitwise identical. PASSED: d_pers and "
                 "d_loc both positively EXCESS (margin 0.05, 2 SE). FAILED: none positive and both decided, or the "
                 "saturation disappears (< 0.001). UNDECIDED: otherwise.",
    "null_model": "the same causal sets without a pulse and without a channel (F3.2)",
    "level": "mechanism",
    "meters": ["K.9"],
    "depends_on": ["F3.3", "F3.3e"],
    "params": {"D": 3, "T": 7.0, "R": 0.5, "R_abs": 0.4, "s0": 1.0, "t1": 3.5, "t2": 5.5, "t_off": 5.5, "post_end": 6.75,
               "cell": 0.1, "slab": 0.25, "min_count": 3, "n_rbins": 16, "n_perm": 20, "q": 0.99, "N": 56000,
               "n_seeds": 128, "n_cal": 8, "omega": 1.0, "rule": "main", "dq_crit": 0.2, "Pi_sat": 0.5, "m": 2.0,
               "margin_pers": 0.05, "margin_loc": 0.05, "n_se": 2.0, "boot_B": 4000, "sat_min": 0.001,
               "balance_tol": 1e-12, "invalid_max": 0.05, "boost_beta": 0.6},
    "parts": {"capacity": {"channel": "capacity"}, "source": {"channel": "source"}},
    "quick": {"N": 6000, "n_seeds": 4, "n_cal": 1, "n_perm": 4, "boot_B": 200},
    "estimate": "each part ~40 min with 8–10 processes (like F3.3d)",
}

ARMS = ("base", "threshold", "combined")


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


def arm_kwargs(p, arm):
    """route_pulse arguments of one arm."""
    if arm == "base":
        return {"rule": "none", "dq_crit": 0.0, "channel": "none"}
    ch = str(p["channel"]) if arm == "combined" else "none"
    return {"rule": str(p["rule"]), "dq_crit": float(p["dq_crit"]), "channel": ch}


def seed_task(p, seed, C):
    """The three arms on one sprinkling: statistics, switch-off survival, Pi of saturated events, instrument checks."""
    S, ptr, idx, r, t, absorb, win = geometry(p, seed)
    s_on = np.full(S.n, float(p["s0"]))
    s_off = np.where(t < float(p["t_off"]), float(p["s0"]), 0.0)
    om, ps, mm = float(p["omega"]), float(p["Pi_sat"]), float(p["m"])
    out = {}
    for arm in ARMS:
        kw = arm_kwargs(p, arm)
        main = route_pulse(ptr, idx, s_on, C, absorb, kw["rule"], om, dq_crit=kw["dq_crit"], channel=kw["channel"],
                           Pi_sat=ps, m=mm)
        swo = route_pulse(ptr, idx, s_off, C, absorb, kw["rule"], om, dq_crit=kw["dq_crit"], channel=kw["channel"],
                          Pi_sat=ps, m=mm)
        chi = (main["stored"] > 0).astype(float)
        loc = ST.localization(S.X, r, chi, win, p, np.random.default_rng(derive_seed(seed, "F3.3f-perm", arm)))
        sat = (chi > 0) & win
        out[arm] = {"sat_fraction": float(chi[win].mean()), "E_loc": loc["E_loc"], "cells": ST.cells(S.X, chi, win, p),
                    "pulse_fraction": float((main["P"][win] > 0).mean()),
                    "Pi_sat_median": float(np.median(main["Pi"][sat])) if sat.any() else float("nan"),
                    "Pi_min": float(main["Pi"][win].min()),
                    "balance": max(pulse_balance(main), pulse_balance(swo)), "survival": survival(S, r, t, swo, p)}
        if arm == "base":
            ref = route_capacity(ptr, idx, s_on, C, absorb)
            out["base_equals_capacity"] = all(np.array_equal(main[k], ref[k]) for k in ("accepted", "stored", "Phi"))
        if arm == "combined":
            beta = float(p["boost_beta"])
            order = np.argsort((t - beta * S.X[:, 1]) / np.sqrt(1.0 - beta * beta), kind="stable")
            b = route_pulse(ptr, idx, s_on, C, absorb, kw["rule"], om, order=order, dq_crit=kw["dq_crit"],
                            channel=kw["channel"], Pi_sat=ps, m=mm)
            out["boost_identical"] = all(np.array_equal(main[k], b[k]) for k in ("accepted", "stored", "Phi", "P", "Pi"))
    return out


def fmt(st):
    return f"{st[0]:+.4f} +- {st[1]:.4f}" if np.isfinite(st[0]) else "undefined"


def run(p, ctx):
    n, B, nse, mc = int(p["n_seeds"]), int(p["boot_B"]), float(p["n_se"]), int(p["min_count"])
    t0 = time.time()
    cal = parallel_map(calibration_task, [(p, derive_seed(ctx.seed, "F3.3f-cal", int(p["N"]), k)) for k in range(int(p["n_cal"]))],
                       ctx.workers, ctx.progress, lambda i, r_, d, nn: f"[F3.3f] calibration {d}/{nn} ({time.time()-t0:.0f}s)")
    C = float(np.quantile(np.concatenate(cal), float(p["q"])))
    seeds = [derive_seed(ctx.seed, "F3.3f", int(p["N"]), k) for k in range(n)]
    rows = parallel_map(seed_task, [(p, sd, C) for sd in seeds], ctx.workers, ctx.progress,
                        lambda i, r_, d, nn: f"[F3.3f] seed {d}/{nn}: saturated base {r_['base']['sat_fraction']:.4f} "
                                             f"threshold {r_['threshold']['sat_fraction']:.4f} combined "
                                             f"{r_['combined']['sat_fraction']:.4f}; E_loc base {r_['base']['E_loc']:+.4f} "
                                             f"combined {r_['combined']['E_loc']:+.4f} ({time.time()-t0:.0f}s)")
    pers = {a: ST.pairs_persistence([r_[a]["cells"] for r_ in rows], mc) for a in ARMS}
    bs = lambda *tag: derive_seed(ctx.seed, "F3.3f-boot", ctx.part, *tag)

    def paired(a, b):
        dpv = np.array([x["delta"] - y["delta"] for x, y in zip(pers[a], pers[b])])
        dlv = np.array([r_[a]["E_loc"] - r_[b]["E_loc"] for r_ in rows])
        return boot_se(dpv, B, bs("pers", a, b)), boot_se(dlv, B, bs("loc", a, b))

    dp, dl = paired("combined", "base")
    dp_ch, dl_ch = paired("combined", "threshold")
    absval = {a: {"pers": boot_se(np.array([x["delta"] for x in pers[a]]), B, bs("p", a)),
                  "loc": boot_se(np.array([r_[a]["E_loc"] for r_ in rows]), B, bs("l", a))} for a in ARMS}
    dec_p = equivalence(dp[0], dp[1], float(p["margin_pers"]), nse)
    dec_l = equivalence(dl[0], dl[1], float(p["margin_loc"]), nse)
    pos = lambda dec, st: dec == "EXCESS" and st[0] > 0
    decided = lambda dec, st: dec == "EQUIVALENT" or (dec == "EXCESS" and st[0] < 0)
    sat = {a: float(np.mean([r_[a]["sat_fraction"] for r_ in rows])) for a in ARMS}
    no_saturation = sat["combined"] < float(p["sat_min"])
    saturating = [r_["combined"]["sat_fraction"] > 0 for r_ in rows]
    invalid = sum(1 for r_, ok in zip(rows, saturating) if ok and not np.isfinite(r_["combined"]["E_loc"])) \
        + sum(1 for j, x in enumerate(pers["combined"]) if saturating[2 * j] and saturating[2 * j + 1] and not np.isfinite(x["delta"]))
    n_stat = n + len(pers["combined"])
    inst = {"base arm = route_capacity": all(r_["base_equals_capacity"] for r_ in rows),
            "balance": max(max(r_[a]["balance"] for a in ARMS) for r_ in rows) <= float(p["balance_tol"]),
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
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")              # all-nan columns: nothing saturated left / no saturated events
        surv = {a: {k: np.nanmean([r_[a]["survival"][k] for r_ in rows], axis=0) for k in ("frac", "pulse", "rs_with_last_on")}
                for a in ARMS}
        pi_med = float(np.nanmedian([r_["combined"]["Pi_sat_median"] for r_ in rows]))
    pi_min = min(r_["combined"]["Pi_min"] for r_ in rows)
    lines = [f"2+1 cylinder T = {p['T']}, N = {p['N']}, s0 = {p['s0']}; analysis |x| < {p['R_abs']}, t in [{p['t1']}, {p['t2']}); "
             f"C = Q_{p['q']} = {C:.4f} ({p['n_cal']} calibration sprinklings, own stream); pulse rule '{p['rule']}', omega = "
             f"{p['omega']}, dQ_crit = {p['dq_crit']} (CHOSEN FROM DATA, F3.3e search); channel '{p['channel']}', Pi_sat = "
             f"{p['Pi_sat']}, m = {p['m']}; {n} sprinklings = {len(pers['base'])} pairs",
             "  saturated fraction: base %.4f, threshold %.4f, combined %.4f; pulse carried: threshold %.4f, combined %.4f"
             % (sat["base"], sat["threshold"], sat["combined"],
                np.mean([r_["threshold"]["pulse_fraction"] for r_ in rows]), np.mean([r_["combined"]["pulse_fraction"] for r_ in rows])),
             f"  combined: Pi of saturated events (median) {pi_med:.3f}; min Pi {pi_min:.3f}",
             "  Delta_pers: " + "; ".join(f"{a} {fmt(absval[a]['pers'])}" for a in ARMS) + "   (F3.2: +0.006 +- 0.006)",
             "  E_loc:      " + "; ".join(f"{a} {fmt(absval[a]['loc'])}" for a in ARMS) + "   (F3.2: +0.004 +- 0.004)",
             f"  paired combined - base: d_pers = {fmt(dp)} (margin {p['margin_pers']}) -> {dec_p}; d_loc = {fmt(dl)} "
             f"(margin {p['margin_loc']}) -> {dec_l}" + ("; no saturation" if no_saturation else ""),
             f"  reported, combined - threshold (what the channel adds): d_pers = {fmt(dp_ch)}; d_loc = {fmt(dl_ch)}",
             "  switch-off at t = %g; the last slab before it, then the slabs of %g after it:" % (p["t_off"], p["slab"])]
    for a in ARMS:
        lines.append(f"     {a:9s} saturated " + " ".join(f"{x:.4f}" for x in surv[a]["frac"])
                     + ("; pulse " + " ".join(f"{x:.4f}" for x in surv[a]["pulse"]) if a != "base" else "")
                     + "; r_s " + " ".join(f"{x:+.3f}" for x in surv[a]["rs_with_last_on"]))
    lines.append("instrument: " + "; ".join(f"{k}: {'OK' if v else 'FAILS'}" for k, v in inst.items())
                 + f" (invalid {invalid}/{n_stat})")
    data = {"C": C, "d_pers": dp, "d_loc": dl, "d_pers_channel": dp_ch, "d_loc_channel": dl_ch,
            "decisions": {"pers": dec_p, "loc": dec_l}, "absolute": absval, "sat_fraction": sat, "survival": surv,
            "Pi_sat_median": pi_med, "Pi_min": pi_min, "instrument": inst,
            "seeds": [{a: {k: v for k, v in r_[a].items() if k != "cells"} for a in ARMS} for r_ in rows]}
    return Result(status, lines, data)
