"""
F3.3g -- Can the bistable saturated phase exist as a bounded droplet, and does it survive perturbations?
=======================================================================================================

Why this test exists
--------------------
F3.3f capacity: the copied pulse with the materialization threshold
(dQ_crit = 0.2) and the clock-rate feedback on the capacity (C_v = C Pi_v,
Pi_sat = 0.5, m = 2) is a BISTABLE medium. Under homogeneous driving the
whole region switches to a saturated phase (99 %, every saturated event at
Pi = Pi_sat) that outlives the source (99 % -> 51 % -> 34 % -> 17 % ->
6.8 % -> 3.4 % per slab of 0.25 after switch-off), but it fills
everything: no boundary. In a bistable medium what decides whether a
bounded object exists is the boundary between the two phases: it advances
(flood), retreats (the droplet dies) or stays (a bounded droplet). Under
homogeneous driving a droplet cannot form; here the saturated phase is
ignited only in a disk on a weak background and left to itself. The
author's plan (22 Oct 2026): first the boundary without the field (this
test), then the droplet as a source of the field Pi of (a), then the field
acting on the routing, then the self-consistent state with the order updated.

What is and is not in the model
-------------------------------
Model: core.flux.pulse.route_pulse, rule 'main', omega = 1, dq_crit = 0.2,
channel 'capacity' (the local response Pi = f(sigma) of axiom 5 lowers
the event's own capacity; this is what makes the medium bistable). NOT in
the model: the field Pi of equation (a) (no source -> field step), any
feedback of Pi on the routing (the flux is split uniformly among the
children), any change of the causal order. Everything runs on a fixed
sprinkling.

Quantities (distinct, never identified with each other)
-------------------------------------------------------
  Phi      flux on a link (uniform share of the parent's accepted flux)
  C_v      capacity of event v (C Pi_v in the bistable arm, C otherwise)
  dQ_v     local overflow offered - accepted at v (not carried in this model)
  chi_v    saturation indicator, chi = 1 if dQ_v > 0
  P_v      the copied pulse (omega or 0)
Per time slab of 0.25 (events with r < R_abs, the absorbing ring excluded):
  N_drop   number of saturated events
  R_rms    sqrt(2 <r^2>) over the saturated events, r from the ignition centre;
           a uniform disk of radius a gives R_rms = a
  R_N      sqrt(N_drop / (pi n_A)), n_A = events per unit area of the slab:
           the radius of a disk holding N_drop events
  B        R_N / R_rms: compactness of the boundary (1 for a filled disk,
           < 1 for a sparse or fragmented set)
  C_drop   saturated fraction of the events with r < R_rms
  Q_drop   sum of dQ over the saturated events (the overflow, not a stock)
  pulse    number of events carrying a pulse
  ring     saturated fraction in the ring 0.3 <= r < R_abs (spontaneous
           ignition of the background, away from the droplet)
The radii are read from the coordinates of the sprinkling: a readout of
where the saturated events are, independent of the rule that produces
them (the rule never uses coordinates).

Setup
-----
2+1 cylinder T = 9, radius 0.5, N = 72 000 (the density of the F3 series),
outflow absorbed beyond R_abs = 0.4. C = Q_0.99 of the uncapped offered
flux at unit drive in the F3 analysis window (the capacity definition of the
whole series, own calibration stream). Drive: background s_bg everywhere
at all times, plus the ignition s_ign = 5 in the ignition region during
t in [3.0, 4.0) (the flux near its steady state, F3 series; strong enough
to switch the region to the saturated phase in every arm, development runs
below).
Observation: slabs from t = 3.0 (the ignition included) to 8.5; early slabs t in [4.5, 5.0), late
slabs t in [8.0, 8.5).

Arms, all on the same sprinklings and the same drive:
  control    no pulse, no channel (route_capacity, F3.2): ignition only
  threshold  pulse with dQ_crit, no channel
  bistable   pulse with dQ_crit and the capacity channel (decides)

Part A: does a boundary exist without the field? (decides)
-----------------------------------------------------------
Ignition disk radius r0 = 0.15; background grid declared now, all values
on the same sprinklings in one part (the sprinkling is the expensive step):
s_bg = 0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8 (no drive ... near the
homogeneous drive that flooded in F3.3f). In a bistable medium the front
between the phases retreats on one side of a critical background and
advances on the other (development runs below: dies at 0, floods at 0.5);
a droplet can only stand near that point, so the grid is dense there.
Physics to keep in mind: in a continuous bistable medium a standing
droplet is usually the critical nucleus, which is unstable (a larger one
grows, a smaller one shrinks); a stable bounded droplet needs something
more, for example fronts pinned by the discreteness of the medium. Part B
(B_small, B_large) is what tells the two apart.
Every grid value gets the decision below; part A PASSES if at least one
value passes (the value is then confirmed on new sprinklings by part B,
whose seed stream is its own). Per sprinkling of the bistable arm: died if N_drop(late) < 0.05
N_drop(early); otherwise growth g = R_rms(late) - R_rms(early). Over the
sprinklings (bootstrap SE):
  PASSED      at most 5 % died, |mean g| + 2 SE <= r0 / 4, and the ring
              stays below 2 x control (no spontaneous ignition)
  FAILED      more than 50 % died (the droplet dies), or mean g - 2 SE > r0/4
              (it grows: flood), or mean g + 2 SE < -r0/4 (it shrinks), or
              the ring >= 2 x control + 0.01 (the background ignites)
  UNDECIDED   otherwise, or the control arm meets PASSED too (not specific)
The tolerance r0 / 4 is the quarter rule of the laboratory.

Part B: does a droplet return after a perturbation? (decides, after A)
-----------------------------------------------------------------------
Only meaningful at a background where part A found a bounded droplet, so
s_bg is a settable parameter of the B parts with no default (it comes from
part A: CHOSEN FROM DATA, stated in the record). Each B part runs two
bistable arms on the same sprinklings, the reference droplet (r0 = 0.15)
and one perturbation, plus the control arm with the perturbation:
  B_small   ignition radius 0.10
  B_large   ignition radius 0.22
  B_asym    elliptic ignition, semi-axes 0.21 x 0.107 (the area of r0 = 0.15)
  B_cut     reference ignition, then the events with x > 0, r < 0.25,
            t in [5.0, 5.25) forward nothing (part of the flux removed)
  (seed: every part already runs 128 independent sprinklings)
Per sprinkling: the perturbed droplet died as in A; otherwise
d = R_rms(late, perturbed) - R_rms(late, reference).
  PASSED      the reference meets the A criterion, at most 5 % of the
              perturbed died, and |mean d| + 2 SE <= r0 / 4 (it returns to
              the radius of the reference)
  FAILED      the reference meets A, and more than 50 % of the perturbed
              died or |mean d| - 2 SE > r0 / 4
  UNDECIDED   otherwise (including: no reference droplet at this s_bg)

Instrument (any failing -> INSTRUMENT_FAILED): control arm = route_capacity
bit for bit, flux balance <= 1e-12 in every arm, the bistable arm
bit-identical in the order of a boosted frame.

Seen before this test was written: F3.3f (above) and development quick
runs of this file (N = 7 000 and 14 000, 6 sprinklings, decide nothing):
with s_ign = 2 in [3.5, 4.0) and no background the disk did not ignite;
with s_ign = 2 from t = 0 or s_ign = 5 in [3.0, 4.0) and no background the
bistable droplet ignited (N ~40 events per slab, R_rms ~0.15-0.19) and
died within ~0.75 after the ignition (the front first ran outwards, R_rms
-> 0.25-0.39, then vanished); with background 0.5 it ran out to the
absorbing ring within ~0.5 (R_rms 0.21 -> 0.37) and filled the region
(ring 48 % saturated). A quick run of the grid (N = 14 000, 6 sprinklings):
dies at s_bg <= 0.2, fills the region (R_rms 0.35-0.40, ring 16-91 %) at
s_bg >= 0.3. These fixed the ignition (s_ign = 5 in [3.0, 4.0)) and the
density of the grid (0.15 and 0.25 added around the change); nothing in
the criteria was set from them.
"""

from __future__ import annotations

import math
import time
import warnings

import numpy as np

from core.flux.capacity import route_capacity
from core.flux.pulse import route_pulse, pulse_balance
from core.order.links import links_csr
from core.order.sprinkle import Cylinder, derive_seed, sprinkle
from core.stats.parallel import parallel_map
from core.stats.stats import boot_se
from core.testing import Result

TEST = {
    "id": "F3.3g",
    "phase": "F3",
    "title": "Umí bistabilní saturovaná fáze existovat jako ohraničená kapka a přežije poruchy?",
    "theory_part": "vznik hmoty: hranice bistabilní saturované fáze (pulz s prahem + kapacita C·Π) na pevném kauzálním pozadí",
    "question": "Udrží se saturovaná fáze zapálená v disku na slabém pozadí jako ohraničená kapka stálého poloměru (A), a vrátí "
                "se po poruše k témuž poloměru (B)?",
    "method": "route_pulse (main, ω = 1, ΔQ_crit = 0,2, kanál capacity, Π_sat = 0,5, m = 2), bez pole (a) a bez zpětné vazby na "
              "směrování; zapálení s_ign = 5 v disku r0 = 0,15 v čase [3,0; 4,0) na pozadí s_bg (A: mřížka 0–0,8 na stejných "
              "sprinklinzích; B: s_bg z A, nové sprinklingy); poloměr R_rms a R_N, kompaktnost, "
              "N, saturace, přetečení dQ po vrstvách 0,25; tři větve na stejných sprinklinzích (bez mechanismu, jen práh, bistabilní)",
    "prediction": "bez předpovědi; bistabilní médium dává buď šířící se, nebo ustupující, nebo stojící hranici",
    "criterion": "MĚŘIDLO: kontrola = route_capacity, bilance <= 1e-12, boostované pořadí bitově shodné. A (každá hodnota s_bg; "
                 "část projde, projde-li aspoň jedna): PROŠEL <= 5 % kapek "
                 "zaniklo, |růst R| + 2 SE <= r0/4 a pozadí se nezapálí; NEPROŠEL > 50 % zaniklo, nebo růst/pokles nad r0/4, "
                 "nebo se zapálí pozadí; jinak NEROZHODNUTO (i když kontrola projde také). B (s_bg z části A): PROŠEL reference "
                 "splní A, <= 5 % porušených zaniklo a |R_porušená − R_reference| + 2 SE <= r0/4; NEPROŠEL reference splní A a "
                 "> 50 % zaniklo nebo rozdíl − 2 SE > r0/4; jinak NEROZHODNUTO.",
    "null_model": "stejné sprinklingy a buzení bez pulzu a bez kanálu (F3.2)",
    "level": "mechanism",
    "meters": ["K.9"],
    "depends_on": ["F3.3f"],
    "params": {"D": 3, "T": 9.0, "R": 0.5, "R_abs": 0.4, "N": 72000, "n_seeds": 128, "n_cal": 8, "q": 0.99,
               "cal_t1": 3.5, "cal_t2": 5.5, "s0": 1.0, "s_ign": 5.0, "t_ign1": 3.0, "t_ign2": 4.0, "r0": 0.15,
               "ign_axes": [0.15, 0.15], "slab": 0.25, "t_obs1": 3.0, "t_obs2": 8.5,
               "early": [4.5, 5.0], "late": [8.0, 8.5], "ring": 0.3, "die_frac": 0.05,
               "cut": None, "omega": 1.0, "rule": "main", "dq_crit": 0.2, "channel": "capacity", "Pi_sat": 0.5, "m": 2.0,
               "max_died_pass": 0.05, "max_died_fail": 0.5, "n_se": 2.0, "boot_B": 4000, "balance_tol": 1e-12,
               "boost_beta": 0.6, "s_bg": None},
    "parts": {"A": {"mode": "A", "bg_grid": [0.0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8]},
              "B_small": {"mode": "B", "pert": {"ign_axes": [0.10, 0.10]}},
              "B_large": {"mode": "B", "pert": {"ign_axes": [0.22, 0.22]}},
              "B_asym": {"mode": "B", "pert": {"ign_axes": [0.21, 0.107]}},
              "B_cut": {"mode": "B", "pert": {"cut": [0.0, 0.25, 5.0, 5.25]}}},
    "settable": {"B_small": ["s_bg"], "B_large": ["s_bg"], "B_asym": ["s_bg"], "B_cut": ["s_bg"]},
    "quick": {"N": 7000, "n_seeds": 4, "n_cal": 1, "boot_B": 200},
    "estimate": "každá část ~50–70 min s 10 procesy (geometrie N = 72 000 převládá; A: 8 hodnot pozadí × 3 větve na "
                "sprinkling, B: 3 větve)",
}


def geometry(p, seed):
    """One cylinder sprinkling, its links, radii and times."""
    S = sprinkle(Cylinder(int(p["D"]), 0.0, float(p["T"]), float(p["R"])), int(p["N"]), seed)
    ptr, idx = links_csr(S.X, S.rho)
    return S, ptr, idx, S.r(), S.X[:, 0]


def calibration_task(p, seed):
    """Uncapped offered flux at unit drive of the F3 analysis window (the capacity definition of the series)."""
    S, ptr, idx, r, t = geometry(p, seed)
    absorb = r > float(p["R_abs"])
    win = (t >= float(p["cal_t1"])) & (t < float(p["cal_t2"])) & (r < float(p["R_abs"]))
    return route_capacity(ptr, idx, np.full(S.n, float(p["s0"])), np.inf, absorb)["offered"][win]


def drive(S, t, p, axes, s_bg):
    """Background s_bg everywhere plus the ignition in an ellipse with the given semi-axes during the ignition window."""
    x, y = S.X[:, 1], S.X[:, 2]
    a, b = float(axes[0]), float(axes[1])
    inside = (x / a) ** 2 + (y / b) ** 2 < 1.0
    on = (t >= float(p["t_ign1"])) & (t < float(p["t_ign2"]))
    return float(s_bg) + float(p["s_ign"]) * (inside & on)


def absorb_mask(S, r, t, p, cut):
    """Outflow beyond R_abs; optionally a cut: events with x > x0, r < rc, t in [tc1, tc2) forward nothing."""
    ab = r > float(p["R_abs"])
    if cut:
        x0, rc, tc1, tc2 = (float(v) for v in cut)
        ab = ab | ((S.X[:, 1] > x0) & (r < rc) & (t >= tc1) & (t < tc2))
    return ab


def slab_series(S, r, t, res, p):
    """Per slab of the observation window: N_drop, R_rms, R_N, B, C_drop, Q_drop, pulse count, ring fraction."""
    Rabs, ring = float(p["R_abs"]), float(p["ring"])
    edges = np.arange(float(p["t_obs1"]), float(p["t_obs2"]) + 1e-9, float(p["slab"]))
    sat = res["stored"] > 0
    out = {k: [] for k in ("t", "N", "R_rms", "R_N", "B", "C_drop", "Q_drop", "pulse", "ring")}
    for k in range(len(edges) - 1):
        m = (r < Rabs) & (t >= edges[k]) & (t < edges[k + 1])
        ms = m & sat
        n_all = int(m.sum())
        nA = n_all / (math.pi * Rabs ** 2)                     # events per unit area in this slab
        Nd = int(ms.sum())
        Rrms = float(np.sqrt(2.0 * np.mean(r[ms] ** 2))) if Nd else float("nan")
        RN = float(np.sqrt(Nd / (math.pi * nA))) if nA > 0 else float("nan")
        inside = m & (r < Rrms) if Nd else np.zeros_like(m)
        ringm = m & (r >= ring)
        out["t"].append(0.5 * (edges[k] + edges[k + 1]))
        out["N"].append(Nd)
        out["R_rms"].append(Rrms)
        out["R_N"].append(RN)
        out["B"].append(RN / Rrms if Nd and Rrms > 0 else float("nan"))
        out["C_drop"].append(float(sat[inside].mean()) if inside.any() else float("nan"))
        out["Q_drop"].append(float(res["stored"][ms].sum()))
        out["pulse"].append(int((res["P"][m] > 0).sum()))
        out["ring"].append(float(sat[ringm].mean()) if ringm.any() else float("nan"))
    return out


def window_stats(ser, p):
    """Early / late means of one arm on one sprinkling: N, R_rms, ring; died flag; growth."""
    t = np.asarray(ser["t"])
    e = (t >= p["early"][0]) & (t < p["early"][1])
    l_ = (t >= p["late"][0]) & (t < p["late"][1])
    N = np.asarray(ser["N"], float)
    R = np.asarray(ser["R_rms"], float)
    Ne, Nl = float(N[e].mean()), float(N[l_].mean())
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        Re, Rl = _nanmean(R[e]), _nanmean(R[l_])
    died = (Ne == 0.0) or (Nl < float(p["die_frac"]) * Ne)
    if died:
        Rl = float("nan")                            # a few stray saturated events are no droplet radius
    return {"N_early": Ne, "N_late": Nl, "R_early": Re, "R_late": Rl, "died": bool(died),
            "growth": (Rl - Re) if not died else float("nan"),
            "ring_late": _nanmean(np.asarray(ser["ring"], float)[l_])}


def _nanmean(v):
    """nanmean without the empty-slice warning (nan if nothing finite)."""
    v = np.asarray(v, float)
    return float(v[np.isfinite(v)].mean()) if np.isfinite(v).any() else float("nan")


def run_arm(ptr, idx, s, C, absorb, p, arm, order=None):
    """route_pulse for one arm: control (no pulse, no channel), threshold (pulse), bistable (pulse + channel)."""
    if arm == "control":
        return route_pulse(ptr, idx, s, C, absorb, "none", float(p["omega"]), order=order)
    ch = str(p["channel"]) if arm.startswith("bistable") else "none"
    return route_pulse(ptr, idx, s, C, absorb, str(p["rule"]), float(p["omega"]), order=order,
                       dq_crit=float(p["dq_crit"]), channel=ch, Pi_sat=float(p["Pi_sat"]), m=float(p["m"]))


def configs(p):
    """(label, s_bg, {arm: (ignition axes, cut)}) for every background of the part."""
    ref_axes = p["ign_axes"]
    if p["mode"] == "A":
        arms = {"control": (ref_axes, None), "threshold": (ref_axes, None), "bistable": (ref_axes, None)}
        return [(f"s_bg {float(b):g}", float(b), arms) for b in p["bg_grid"]]
    pert = p["pert"]
    pax, pcut = pert.get("ign_axes", ref_axes), pert.get("cut")
    return [(f"s_bg {float(p['s_bg']):g}", float(p["s_bg"]),
             {"control": (pax, pcut), "bistable_ref": (ref_axes, None), "bistable": (pax, pcut)})]


def seed_task(p, seed, C):
    """All configurations and arms of one sprinkling; per configuration {arm: stats}, and the instrument checks."""
    S, ptr, idx, r, t = geometry(p, seed)
    beta = float(p["boost_beta"])
    border = np.argsort((t - beta * S.X[:, 1]) / np.sqrt(1.0 - beta * beta), kind="stable")
    out = {"cfg": [], "control_equals_capacity": True, "boost_identical": True, "balance": 0.0}
    for label, s_bg, arms in configs(p):
        one = {}
        for arm, (axes, cut) in arms.items():
            s = drive(S, t, p, axes, s_bg)
            ab = absorb_mask(S, r, t, p, cut)
            res = run_arm(ptr, idx, s, C, ab, p, arm)
            out["balance"] = max(out["balance"], pulse_balance(res))
            ser = slab_series(S, r, t, res, p)
            one[arm] = {"series": ser, **window_stats(ser, p)}
            if arm == "control":
                ref = route_capacity(ptr, idx, s, C, ab)
                out["control_equals_capacity"] &= all(np.array_equal(res[k], ref[k]) for k in ("accepted", "stored", "Phi"))
            if arm == "bistable":
                b = run_arm(ptr, idx, s, C, ab, p, arm, order=border)
                out["boost_identical"] &= all(np.array_equal(res[k], b[k]) for k in ("accepted", "stored", "Phi", "P", "Pi"))
        out["cfg"].append(one)
    return out


def summarise(rows, arm, p, B, bs):
    """Fraction died, mean growth with bootstrap SE, ring, and the mean slab series of one arm."""
    died = np.array([r_[arm]["died"] for r_ in rows])
    g = np.array([r_[arm]["growth"] for r_ in rows])
    gs = boot_se(g[np.isfinite(g)], B, bs("g", arm)) if np.isfinite(g).sum() >= 2 else (float("nan"), float("nan"), None)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")              # all-nan columns: no saturated event left in a slab
        ring = float(np.nanmean([r_[arm]["ring_late"] for r_ in rows]))
        mean_ser = {k: np.nanmean(np.array([r_[arm]["series"][k] for r_ in rows], float), axis=0)
                    for k in ("t", "N", "R_rms", "R_N", "B", "C_drop", "Q_drop", "pulse", "ring")}
        R_early = float(np.nanmean([r_[arm]["R_early"] for r_ in rows]))
        R_late = float(np.nanmean([r_[arm]["R_late"] for r_ in rows]))
    return {"died_frac": float(died.mean()), "growth": gs, "ring_late": ring, "series": mean_ser,
            "R_early": R_early, "R_late": R_late}


def verdict_A(st, ctrl, p):
    """The part-A decision for one arm (st) against the control arm (ctrl)."""
    r0, nse = float(p["r0"]), float(p["n_se"])
    g, se = st["growth"][0], st["growth"][1]
    ring_bad = np.isfinite(st["ring_late"]) and st["ring_late"] >= 2.0 * ctrl["ring_late"] + 0.01
    ring_ok = np.isfinite(st["ring_late"]) and st["ring_late"] < 2.0 * ctrl["ring_late"] + 0.01
    if st["died_frac"] > float(p["max_died_fail"]) or ring_bad:
        return "FAILED"
    if np.isfinite(g) and np.isfinite(se):
        if g - nse * se > r0 / 4 or g + nse * se < -r0 / 4:
            return "FAILED"
        if st["died_frac"] <= float(p["max_died_pass"]) and abs(g) + nse * se <= r0 / 4 and ring_ok:
            return "PASSED"
    return "UNDECIDED"


def series_lines(name, st):
    s = st["series"]
    f = lambda v, fmt: " ".join(fmt.format(x) if np.isfinite(x) else "  nan" for x in v)
    return [f"     {name}: died {st['died_frac']:.3f}; R_rms early {st['R_early']:.3f} late {st['R_late']:.3f}; growth "
            + (f"{st['growth'][0]:+.4f} +- {st['growth'][1]:.4f}" if np.isfinite(st["growth"][0]) else "undefined")
            + f"; ring (late) {st['ring_late']:.4f}",
            "        N      " + f(s["N"], "{:6.1f}"),
            "        R_rms  " + f(s["R_rms"], "{:6.3f}"),
            "        B      " + f(s["B"], "{:6.3f}"),
            "        C_drop " + f(s["C_drop"], "{:6.3f}"),
            "        Q_drop " + f(s["Q_drop"], "{:6.2f}"),
            "        pulse  " + f(s["pulse"], "{:6.1f}")]


def run(p, ctx):
    p = dict(p)
    p["mode"] = str(p.get("mode", "A"))
    if p["mode"] == "B" and p.get("s_bg") is None:
        raise SystemExit("F3.3g part B: give the background where part A found a bounded droplet, e.g. --set s_bg=0.3")
    n, B = int(p["n_seeds"]), int(p["boot_B"])
    t0 = time.time()
    cal = parallel_map(calibration_task, [(p, derive_seed(ctx.seed, "F3.3g-cal", int(p["N"]), k)) for k in range(int(p["n_cal"]))],
                       ctx.workers, ctx.progress, lambda i, r_, d, nn: f"[F3.3g] calibration {d}/{nn} ({time.time()-t0:.0f}s)")
    C = float(np.quantile(np.concatenate(cal), float(p["q"])))
    stream = "F3.3g-" + p["mode"]                        # part B confirms on sprinklings part A never saw
    seeds = [derive_seed(ctx.seed, stream, int(p["N"]), k) for k in range(n)]
    cf = configs(p)

    def prog(i, r_, d, nn):
        bits = " ".join(f"{lab.split()[-1]}:{c['bistable']['N_late']:.0f}" for (lab, _, _), c in zip(cf, r_["cfg"]))
        return f"[F3.3g] seed {d}/{nn}: bistable N late per s_bg {bits} ({time.time()-t0:.0f}s)"

    rows = parallel_map(seed_task, [(p, sd, C) for sd in seeds], ctx.workers, ctx.progress, prog)
    inst = {"control arm = route_capacity": all(r_["control_equals_capacity"] for r_ in rows),
            "balance": max(r_["balance"] for r_ in rows) <= float(p["balance_tol"]),
            "boosted order bit-identical": all(r_["boost_identical"] for r_ in rows)}
    r0, nse = float(p["r0"]), float(p["n_se"])
    head = (f"2+1 cylinder T = {p['T']}, N = {p['N']}, R_abs = {p['R_abs']}; C = Q_{p['q']} = {C:.4f} ({p['n_cal']} calibration "
            f"sprinklings, own stream); ignition s_ign = {p['s_ign']} in t [{p['t_ign1']}, {p['t_ign2']}); pulse '{p['rule']}', "
            f"dQ_crit = {p['dq_crit']}, channel '{p['channel']}', Pi_sat = {p['Pi_sat']}, m = {p['m']}; {n} sprinklings "
            f"(stream {stream}); slabs of {p['slab']} from t = {p['t_obs1']} (early {p['early']}, late {p['late']}); "
            f"tolerance r0/4 = {r0 / 4:.4f}")
    lines, summaries, verdicts = [head], [], []
    for j, (label, s_bg, arms) in enumerate(cf):
        sub = [{**r_["cfg"][j]} for r_ in rows]
        bs = lambda *tag, j=j: derive_seed(ctx.seed, "F3.3g-boot", ctx.part, j, *tag)
        st = {a: summarise(sub, a, p, B, bs) for a in arms}
        lines.append(f"  {label}:")
        for a in arms:
            lines += series_lines(a, st[a])
        if p["mode"] == "A":
            v_b, v_c = verdict_A(st["bistable"], st["control"], p), verdict_A(st["control"], st["control"], p)
            v = "UNDECIDED" if (v_b == "PASSED" and v_c == "PASSED") else v_b
            lines.append(f"     -> bistable arm {v_b}; control arm alone {v_c}; value {v}")
        else:
            ref_ok = verdict_A(st["bistable_ref"], st["control"], p) == "PASSED"
            d = np.array([c["bistable"]["R_late"] - c["bistable_ref"]["R_late"]
                          if not (c["bistable"]["died"] or c["bistable_ref"]["died"]) else float("nan") for c in sub])
            ds = boot_se(d[np.isfinite(d)], B, bs("d")) if np.isfinite(d).sum() >= 2 else (float("nan"), float("nan"), None)
            died = st["bistable"]["died_frac"]
            if not ref_ok:
                v = "UNDECIDED"
            elif died > float(p["max_died_fail"]) or (np.isfinite(ds[0]) and abs(ds[0]) - nse * ds[1] > r0 / 4):
                v = "FAILED"
            elif died <= float(p["max_died_pass"]) and np.isfinite(ds[0]) and abs(ds[0]) + nse * ds[1] <= r0 / 4:
                v = "PASSED"
            else:
                v = "UNDECIDED"
            lines[0] += "; s_bg CHOSEN FROM DATA (part A)"
            lines.append(f"     perturbation {p['pert']}: reference droplet meets A: {ref_ok}; perturbed died {died:.3f}; "
                         "R_late perturbed - reference = "
                         + (f"{ds[0]:+.4f} +- {ds[1]:.4f}" if np.isfinite(ds[0]) else "undefined") + f" -> {v}")
            st["difference"] = ds
        summaries.append({"label": label, "s_bg": s_bg, "summary": st, "verdict": v})
        verdicts.append(v)
    if p["mode"] == "A":
        passing = [s["label"] for s in summaries if s["verdict"] == "PASSED"]
        if passing:
            status = "PASSED"
        elif all(v == "FAILED" for v in verdicts):
            status = "FAILED"
        else:
            status = "UNDECIDED"
        lines.append("  part A: values passing: " + (", ".join(passing) if passing else "none")
                     + "; confirm a passing value with the B parts (--set s_bg=...)")
    else:
        status = verdicts[0]
    if not all(inst.values()):
        status = "INSTRUMENT_FAILED"
    lines.append("instrument: " + "; ".join(f"{k}: {'OK' if v else 'FAILS'}" for k, v in inst.items()))
    data = {"C": C, "configs": summaries, "instrument": inst,
            "seeds": [[{a: {k: v for k, v in c[a].items() if k != "series"} for a in c} for c in r_["cfg"]] for r_ in rows]}
    return Result(status, lines, data)
