"""
F3.3e -- Where does the materialization threshold work best, and does it hold on fresh data?
============================================================================================

Why this test exists
--------------------
F3.3d (pre-declared scan, rule 'main', C = Q_0.99) put the paired
localization of the copied pulse with the threshold dQ_crit above zero
for the first time in the F3 series, but below the margin 0.05:
    dQ_crit = 0     (F3.3c)  d_loc -0.004 +- 0.004   flood (30 % saturated)
    dQ_crit = 0.25           d_loc +0.019 +- 0.004   (4.7 SE above zero)
    dQ_crit = 0.5            d_loc +0.007 +- 0.002
so the localization has a maximum somewhere in (0, 0.5). The author's
estimate (21 Oct 2026): the window of materialization lies near
dQ_crit 0.15-0.25. This test looks for the best value and checks it on
data that took no part in the choice.

Two parts, two roles
--------------------
search   (exploratory, status RECORDED, no verdict) A grid of dQ_crit
         values, all on the SAME causal sets (the sprinkling is the
         expensive step; one sprinkling serves every value), so the curve
         is compared value against value without the noise of different
         graphs. The grid is set at the command line (settable lo, hi,
         step, n_seeds), so a coarse grid can be followed by a finer one
         around its maximum. Every value is analysed exactly as in F3.3d
         (paired d_pers and d_loc against the arm without the pulse).
         The quantity searched for is the one the criterion decides on:
         PASSED needs both d_pers and d_loc above the margin by 2 SE, so
         the objective is  min(d_pers - 2 SE, d_loc - 2 SE)  and a value
         would pass if it exceeds 0.05. No new weight or rule is added.
         The per-value verdict printed by the search is for orientation
         only: the maximum of a noisy curve is biased upwards (the winner
         of many draws), which is why it is not a result.
confirm  (decides) One value dQ_crit, given at the command line (settable
         dq_crit; there is no default, the value comes from a search),
         128 NEW sprinklings from their own seed stream, independent of
         the search and of F3.3d. Criteria exactly F3.3d. The record says
         that the value was chosen from data.

Guard against choosing on the confirmation data: a full confirm run
refuses a second value on the same seed (the same confirmation
sprinklings); another value needs another base seed (--seed K), and the
output lists every earlier confirmation, so a series of attempts is
visible in the record.

What a PASSED confirm would and would not mean: a localized, persistent
structure at ONE value selected from data. The roadmap gate of F3 asks for
a structure without tuning over a wide range; whether the window is wide
is what the shape of the search curve shows (a broad plateau or a sharp
peak), and F3.4 decides it. A FAILED confirm at the optimum closes the
threshold model on a fixed causal set, not only within the F3.3d scan.

Model, setup and statistics: exactly F3.3d (core.flux.pulse.route_pulse
with dq_crit, rule 'main', omega = 1, C = Q_0.99 from 8 independent
calibration sprinklings of this test's own stream, shared by both parts;
2+1 cylinder T = 7, N = 56 000, window t in [3.5, 5.5), |x| < 0.4;
switch-off reported; boosted-order identity for every value).

Seen before this test was written: F3.3c main_q99, F3.3d crit025 and
crit050 (numbers above), the author's estimate of the window. Nothing in
the criteria was set from them; the default search grid 0.05-0.40 by
0.025 spans the region between the flood and the F3.3d value 0.5 on both
sides of the estimate.
"""

from __future__ import annotations

import time
import warnings

import numpy as np

from core import provenance
from core.flux.capacity import route_capacity
from core.flux.pulse import route_pulse, pulse_balance
from core.flux import structure as ST
from core.order.links import links_csr
from core.order.sprinkle import Cylinder, derive_seed, sprinkle
from core.stats.parallel import parallel_map
from core.stats.stats import boot_se, equivalence
from core.testing import Result

TEST = {
    "id": "F3.3e",
    "phase": "F3",
    "title": "Search for the best materialisation threshold ΔQ_crit and its confirmation on new data",
    "theory_part": "materialisation as a threshold: a bound state of the pulse only where the excess load reaches "
                   "ΔQ_crit",
    "question": "At which ΔQ_crit does a copied pulse create the strongest localised and persistent structure, and "
                "does that value pass the criterion of F3.3d on causal sets that took no part in the selection?",
    "method": "search: a grid of ΔQ_crit (lo, hi, step from the command line) on the same causal sets, paired d_pers "
              "and d_loc as F3.3d, the quantity searched is min(d_pers − 2 SE, d_loc − 2 SE); confirm: one value "
              "from the command line, 128 new sprinklings from their own seed stream, the criteria of F3.3d",
    "prediction": "according to the author's estimate a materialisation window around ΔQ_crit 0.15–0.25; confirm "
                  "decides",
    "criterion": "search: exploration, RECORDED, no verdict (the verdict per value only indicative). confirm: "
                 "INSTRUMENT as F3.3d (control = route_capacity, balance <= 1e-12, invalid statistics <= 5 %, the "
                 "order of a boosted frame bitwise identical). PASSED: d_pers and d_loc both positively EXCESS "
                 "(margin 0.05, 2 SE). FAILED: none positive and both decided, or the saturation disappears (< "
                 "0.001). UNDECIDED: otherwise. The value was selected from the data, and the record says so.",
    "null_model": "the same causal sets without a pulse (F3.2)",
    "level": "mechanism",
    "meters": ["K.9"],
    "depends_on": ["F3.3d"],
    "params": {"D": 3, "T": 7.0, "R": 0.5, "R_abs": 0.4, "s0": 1.0, "t1": 3.5, "t2": 5.5, "t_off": 5.5, "post_end": 6.75,
               "cell": 0.1, "slab": 0.25, "min_count": 3, "n_rbins": 16, "n_perm": 20, "q": 0.99, "N": 56000,
               "n_seeds": 128, "n_cal": 8, "omega": 1.0, "rule": "main", "margin_pers": 0.05, "margin_loc": 0.05,
               "n_se": 2.0, "boot_B": 4000, "sat_min": 0.001, "balance_tol": 1e-12, "invalid_max": 0.05, "boost_beta": 0.6},
    "parts": {"search": {"mode": "search", "lo": 0.05, "hi": 0.40, "step": 0.025, "n_seeds": 64},
              "confirm": {"mode": "confirm", "dq_crit": None}},
    "settable": {"search": ["lo", "hi", "step", "n_seeds"], "confirm": ["dq_crit"]},
    "quick": {"N": 6000, "n_seeds": 4, "n_cal": 1, "n_perm": 4, "boot_B": 200},
    "estimate": "search: one causal set serves all values; 64 seeds and 15 values ~25–30 min, every further 15 "
                "values ~+8 min; confirm ~40 min (like one part of F3.3d); with 8–10 processes",
}


# ----------------------------------------------------------------------------------------------- the model (as F3.3d)

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


def vtag(dq):
    """Stable text form of a threshold value: the permutation and bootstrap streams depend on the value, not the grid."""
    return f"{float(dq):.6g}"


def arm_stats(S, r, t, win, p, main, swo, seed, tag):
    """Statistics of one arm on one sprinkling (as F3.3d)."""
    chi = (main["stored"] > 0).astype(float)
    loc = ST.localization(S.X, r, chi, win, p, np.random.default_rng(derive_seed(seed, "F3.3e-perm", tag)))
    return {"sat_fraction": float(chi[win].mean()), "E_loc": loc["E_loc"], "I_obs": loc["I_obs"],
            "cells": ST.cells(S.X, chi, win, p), "pulse_fraction": float((main["P"][win] > 0).mean()),
            "balance": max(pulse_balance(main), pulse_balance(swo)), "survival": survival(S, r, t, swo, p)}


def seed_task(p, seed, C, grid):
    """One sprinkling: the arm without the pulse once, then the pulse arm for every threshold value of the grid."""
    S, ptr, idx, r, t, absorb, win = geometry(p, seed)
    s_on = np.full(S.n, float(p["s0"]))
    s_off = np.where(t < float(p["t_off"]), float(p["s0"]), 0.0)
    om, rule = float(p["omega"]), str(p["rule"])
    ref = route_capacity(ptr, idx, s_on, C, absorb)
    # without the pulse the threshold plays no role (P stays 0): one control arm serves every value
    main = route_pulse(ptr, idx, s_on, C, absorb, "none", om)
    swo = route_pulse(ptr, idx, s_off, C, absorb, "none", om)
    out = {"off": arm_stats(S, r, t, win, p, main, swo, seed, "off"),
           "off_equals_capacity": (np.array_equal(main["accepted"], ref["accepted"])
                                   and np.array_equal(main["stored"], ref["stored"])
                                   and np.array_equal(main["Phi"], ref["Phi"])),
           "on": []}
    beta = float(p["boost_beta"])
    border = np.argsort((t - beta * S.X[:, 1]) / np.sqrt(1.0 - beta * beta), kind="stable")
    for dq in grid:
        main = route_pulse(ptr, idx, s_on, C, absorb, rule, om, dq_crit=dq)
        swo = route_pulse(ptr, idx, s_off, C, absorb, rule, om, dq_crit=dq)
        st = arm_stats(S, r, t, win, p, main, swo, seed, "on-" + vtag(dq))
        boosted = route_pulse(ptr, idx, s_on, C, absorb, rule, om, order=border, dq_crit=dq)
        st["boost_identical"] = all(np.array_equal(main[k], boosted[k]) for k in ("accepted", "stored", "Phi", "P"))
        out["on"].append(st)
    return out


# ----------------------------------------------------------------------------------------------- analysis of one value

def analyse(rows, j, dq, p, ctx, pers_off):
    """Paired statistics, instruments and the F3.3d decision rule for grid value j (dq) over all sprinklings."""
    n, B, nse, mc = len(rows), int(p["boot_B"]), float(p["n_se"]), int(p["min_count"])
    on = [r_["on"][j] for r_ in rows]
    off = [r_["off"] for r_ in rows]
    pers_on = ST.pairs_persistence([x["cells"] for x in on], mc)
    d_pers = np.array([a["delta"] - b["delta"] for a, b in zip(pers_on, pers_off)])
    d_loc = np.array([a["E_loc"] - b["E_loc"] for a, b in zip(on, off)])
    bs = lambda tag: derive_seed(ctx.seed, "F3.3e-boot", ctx.part, tag, vtag(dq))
    dp, dl = boot_se(d_pers, B, bs("pers")), boot_se(d_loc, B, bs("loc"))
    absval = {"on": {"pers": boot_se(np.array([x["delta"] for x in pers_on]), B, bs("pon")),
                     "loc": boot_se(np.array([x["E_loc"] for x in on]), B, bs("lon"))}}
    dec_p = equivalence(dp[0], dp[1], float(p["margin_pers"]), nse)
    dec_l = equivalence(dl[0], dl[1], float(p["margin_loc"]), nse)
    pos = lambda dec, st: dec == "EXCESS" and st[0] > 0
    decided = lambda dec, st: dec == "EQUIVALENT" or (dec == "EXCESS" and st[0] < 0)
    sat_on = float(np.mean([x["sat_fraction"] for x in on]))
    no_saturation = sat_on < float(p["sat_min"])
    saturating = [x["sat_fraction"] > 0 for x in on]
    invalid = sum(1 for x, ok in zip(on, saturating) if ok and not np.isfinite(x["E_loc"])) \
        + sum(1 for k, x in enumerate(pers_on) if saturating[2 * k] and saturating[2 * k + 1] and not np.isfinite(x["delta"]))
    n_stat = n + len(pers_on)
    inst = {"control arm = route_capacity": all(r_["off_equals_capacity"] for r_ in rows),
            "balance": max(max(r_["off"]["balance"], r_["on"][j]["balance"]) for r_ in rows) <= float(p["balance_tol"]),
            "invalid statistics": no_saturation or invalid <= float(p["invalid_max"]) * n_stat,
            "boosted order bit-identical": all(x["boost_identical"] for x in on)}
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
    lower = [st[0] - nse * st[1] for st in (dp, dl)]
    objective = float(min(lower)) if all(np.isfinite(lower)) else float("nan")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")              # all-nan columns: no saturated cells left after switch-off
        surv = {k: np.nanmean([x["survival"][k] for x in on], axis=0) for k in ("frac", "pulse", "rs_with_last_on")}
    return {"dq_crit": float(dq), "status": status, "d_pers": dp, "d_loc": dl, "decisions": {"pers": dec_p, "loc": dec_l},
            "objective": objective, "absolute": absval, "sat_on": sat_on, "no_saturation": no_saturation,
            "pulse_fraction": float(np.mean([x["pulse_fraction"] for x in on])), "survival": surv,
            "instrument": inst, "invalid": invalid, "n_stat": n_stat, "pairs": pers_on}


# ----------------------------------------------------------------------------------------------- run

def grid_of(p):
    """The search grid lo, lo + step, ..., <= hi (validated)."""
    lo, hi, step = float(p["lo"]), float(p["hi"]), float(p["step"])
    if not (lo >= 0.0 and hi >= lo and step > 0.0):
        raise SystemExit(f"F3.3e search: need 0 <= lo <= hi and step > 0 (got lo={lo}, hi={hi}, step={step})")
    k = int(np.floor((hi - lo) / step + 1e-9))
    return [round(lo + i * step, 10) for i in range(k + 1)]


def confirm_guard(p, ctx):
    """Refuse a second confirmation value on the same seed; list earlier confirmations and the latest search."""
    dq = p.get("dq_crit")
    if dq is None:
        raise SystemExit("F3.3e confirm: give the value chosen from a search, e.g. --set dq_crit=0.2")
    if not (isinstance(dq, (int, float)) and float(dq) >= 0.0):
        raise SystemExit(f"F3.3e confirm: dq_crit must be a number >= 0 (got {dq!r})")
    recs = provenance.load_all("F3.3e")
    conf = [r_ for r_ in recs if r_.get("part") == "confirm"]
    if not ctx.quick:
        clash = [r_ for r_ in conf if r_["seed"] == ctx.seed and r_["params"].get("dq_crit") != float(dq)]
        if clash:
            used = sorted({r_["params"].get("dq_crit") for r_ in clash})
            raise SystemExit(f"F3.3e confirm: seed {ctx.seed} already confirmed dQ_crit = {used}; another value on the same "
                             f"confirmation sprinklings would choose on them. Use another base seed (--seed K).")
    lines = ["  earlier confirmations: " + ("none" if not conf else ", ".join(
        f"dQ_crit {r_['params'].get('dq_crit')} seed {r_['seed']} {r_['status']}" for r_ in conf))]
    srch = [r_ for r_ in recs if r_.get("part") == "search" and (r_.get("data") or {}).get("best")]
    if srch:
        b = srch[-1]["data"]["best"]
        lines.append(f"  latest search ({srch[-1]['time_utc']}, seed {srch[-1]['seed']}): best dQ_crit = {b.get('dq_crit')}, "
                     f"min(d - 2 SE) = {b.get('objective')}")
    return float(dq), lines


def fmt(st):
    return f"{st[0]:+.4f} +- {st[1]:.4f}" if np.isfinite(st[0]) else "undefined"


def run(p, ctx):
    mode = str(p["mode"])
    n = int(p["n_seeds"])
    if n < 4 or n % 2:
        raise SystemExit(f"F3.3e: n_seeds must be even and >= 4 (pairs of sprinklings; got {n})")
    pre = []
    if mode == "search":
        grid = grid_of(p)
    else:
        dq, pre = confirm_guard(p, ctx)
        grid = [dq]
    t0 = time.time()
    cal = parallel_map(calibration_task, [(p, derive_seed(ctx.seed, "F3.3e-cal", int(p["N"]), k)) for k in range(int(p["n_cal"]))],
                       ctx.workers, ctx.progress, lambda i, r_, d, nn: f"[F3.3e] calibration {d}/{nn} ({time.time()-t0:.0f}s)")
    C = float(np.quantile(np.concatenate(cal), float(p["q"])))
    seeds = [derive_seed(ctx.seed, "F3.3e-" + mode, int(p["N"]), k) for k in range(n)]

    def prog(i, r_, d, nn):
        sat = [x["sat_fraction"] for x in r_["on"]]
        el = [x["E_loc"] for x in r_["on"]]
        return (f"[F3.3e] seed {d}/{nn}: saturated off {r_['off']['sat_fraction']:.4f} on {min(sat):.4f}-{max(sat):.4f}; "
                f"E_loc off {r_['off']['E_loc']:+.4f} on {np.nanmin(el):+.4f}..{np.nanmax(el):+.4f} ({time.time()-t0:.0f}s)")

    rows = parallel_map(seed_task, [(p, sd, C, grid) for sd in seeds], ctx.workers, ctx.progress, prog)
    mc, B = int(p["min_count"]), int(p["boot_B"])
    pers_off = ST.pairs_persistence([r_["off"]["cells"] for r_ in rows], mc)
    bs = lambda tag: derive_seed(ctx.seed, "F3.3e-boot", ctx.part, tag)
    off_abs = {"pers": boot_se(np.array([x["delta"] for x in pers_off]), B, bs("poff")),
               "loc": boot_se(np.array([r_["off"]["E_loc"] for r_ in rows]), B, bs("loff"))}
    sat_off = float(np.mean([r_["off"]["sat_fraction"] for r_ in rows]))
    res = [analyse(rows, j, dq, p, ctx, pers_off) for j, dq in enumerate(grid)]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        surv_off = {k: np.nanmean([r_["off"]["survival"][k] for r_ in rows], axis=0) for k in ("frac", "rs_with_last_on")}
    head = (f"2+1 cylinder T = {p['T']}, N = {p['N']}, s0 = {p['s0']}; analysis |x| < {p['R_abs']}, t in [{p['t1']}, {p['t2']}); "
            f"C = Q_{p['q']} = {C:.4f} ({p['n_cal']} calibration sprinklings, own stream); copied pulse, rule '{p['rule']}', "
            f"omega = {p['omega']}; {n} sprinklings = {len(pers_off)} pairs")
    off_line = (f"  no pulse: saturated {sat_off:.4f}; Delta_pers {fmt(off_abs['pers'])}; E_loc {fmt(off_abs['loc'])}   "
                f"(F3.2: +0.006 +- 0.006, +0.004 +- 0.004)")
    inst_all = {k: all(x["instrument"][k] for x in res) for k in res[0]["instrument"]}
    inv = f"invalid {max(x['invalid'] for x in res)}/{res[0]['n_stat']} (worst value)"
    inst_line = "instrument: " + "; ".join(f"{k}: {'OK' if v else 'FAILS'}" for k, v in inst_all.items()) + f" ({inv})"
    data = {"C": C, "grid": grid, "sat_off": sat_off, "absolute_off": off_abs, "survival_off": surv_off,
            "values": [{k: v for k, v in x.items() if k != "pairs"} for x in res], "instrument": inst_all,
            "seeds": [{"off": {k: v for k, v in r_["off"].items() if k != "cells"},
                       "on": [{k: v for k, v in x.items() if k != "cells"} for x in r_["on"]]} for r_ in rows]}

    if mode == "search":
        lines = [head + f"; search grid dQ_crit = {grid[0]:g} .. {grid[-1]:g} by {float(p['step']):g} ({len(grid)} values, "
                 f"the same sprinklings for all)", off_line,
                 "  dQ_crit  saturated  pulse    d_pers             d_loc              min(d - 2SE)  orientation   "
                 "after switch-off: saturated / pulse in slab 1"]
        for x in res:
            lines.append(f"  {x['dq_crit']:7.4g}  {x['sat_on']:9.4f}  {x['pulse_fraction']:.4f}  {fmt(x['d_pers']):17s}  "
                         f"{fmt(x['d_loc']):17s}  {x['objective']:+12.4f}  {x['status']:12s}  "
                         f"{x['survival']['frac'][1]:.4f} / {x['survival']['pulse'][1]:.4f}")
        ok = [x for x in res if np.isfinite(x["objective"])]
        if ok:
            best = max(ok, key=lambda x: x["objective"])
            above = [f"{x['dq_crit']:g}" for x in ok if x["objective"] > float(p["margin_loc"])]
            dl_best = max(ok, key=lambda x: x["d_loc"][0])
            lines += [f"  best by the criterion quantity: dQ_crit = {best['dq_crit']:g}, min(d - 2 SE) = {best['objective']:+.4f} "
                      f"(a pass needs > {p['margin_loc']}); largest d_loc {fmt(dl_best['d_loc'])} at dQ_crit = {dl_best['dq_crit']:g}",
                      f"  values with min(d - 2 SE) > margin: {', '.join(above) if above else 'none'}",
                      "  orientation only: the maximum of a noisy curve is biased upwards; the value is confirmed on new "
                      f"sprinklings: python run.py F3.3e --part confirm --set dq_crit={best['dq_crit']:g}"]
            data["best"] = {"dq_crit": best["dq_crit"], "objective": best["objective"], "d_pers": best["d_pers"],
                            "d_loc": best["d_loc"]}
        else:
            lines.append("  no value has defined statistics")
        lines.append(inst_line)
        status = "RECORDED" if all(inst_all.values()) else "INSTRUMENT_FAILED"
        return Result(status, lines, data)

    x = res[0]
    lines = [head + f"; dQ_crit = {x['dq_crit']:g} CHOSEN FROM DATA (command line), confirmation sprinklings independent of "
             "the search and of F3.3d", *pre, off_line,
             f"  pulse: saturated {x['sat_on']:.4f}; events carrying a pulse {x['pulse_fraction']:.4f}; "
             f"Delta_pers {fmt(x['absolute']['on']['pers'])}; E_loc {fmt(x['absolute']['on']['loc'])}",
             f"  paired d_pers = {fmt(x['d_pers'])} (margin {p['margin_pers']}) -> {x['decisions']['pers']}; d_loc = "
             f"{fmt(x['d_loc'])} (margin {p['margin_loc']}) -> {x['decisions']['loc']}"
             + ("; no saturation with the pulse" if x["no_saturation"] else ""),
             "  switch-off at t = %g; the last slab before it, then the slabs of %g after it:" % (p["t_off"], p["slab"]),
             "     saturated fraction  no pulse " + " ".join(f"{v:.4f}" for v in surv_off["frac"])
             + "; pulse " + " ".join(f"{v:.4f}" for v in x["survival"]["frac"]),
             "     carrying a pulse    " + " ".join(f"{v:.4f}" for v in x["survival"]["pulse"]),
             "     r_s with the last slab before (nan: no saturated cells left): no pulse "
             + " ".join(f"{v:+.3f}" for v in surv_off["rs_with_last_on"]) + "; pulse "
             + " ".join(f"{v:+.3f}" for v in x["survival"]["rs_with_last_on"]) + "  (reported)",
             inst_line]
    return Result(x["status"], lines, data)
