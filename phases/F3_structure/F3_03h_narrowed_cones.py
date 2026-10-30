"""
F3.3h -- Do narrowed future cones (c_eff = c Pi^2) hold the saturated phase together as a bounded droplet?
=========================================================================================================

Why this test exists
--------------------
F3.3g A (22 Oct 2026, FAILED at every background 0-0.8): the bistable
saturated phase (copied pulse with the threshold dQ_crit = 0.2 and the
capacity channel C Pi, F3.3f) has NO boundary. Right after the ignition
the droplet dissolves into a sparse cloud: R_rms 0.15 -> ~0.31 within
0.25-0.5 (about the speed of light), compactness B ~0.8 -> ~0.45, the
saturation inside ~0.5 -> 0.13-0.25, whatever the background; the cloud
then fades (s_bg <= 0.25) or fills the region (s_bg >= 0.3). The reason
is the routing: every event splits its flux uniformly over the children in
its whole future cone, so the influence of a saturated region spreads at
the speed of light and nothing holds a boundary.
The author's hypothesis (21-22 Oct 2026): a slow clock narrows the causal
cone of the event's children (self-focusing, "a causal lens"); gravity
compresses matter. The laboratory has the size of the narrowing from an
independent result: full light bending needs c_eff = c Pi^2 (F7). Here
the narrowing is put into the routing in its smallest, local form.

Model
-----
F3.3g with one change (core.flux.pulse.route_pulse, narrow_power = 2): an
event v forwards its accepted flux only to the linked children k inside
its narrowed future cone |x_k - x_v| <= Pi_v^2 (t_k - t_v), uniformly
among them; Pi_v is v's own clock rate, known when v forwards (the rule is
a causal pass, no stationary state). With Pi_v = 1 nothing changes. If no
linked child lies inside (a sprinkling has few nearly timelike links), v
forwards to its single most timelike child: the flux stays in the medium.
The share of such events is reported.
What this is and is not: a stand-in for axiom 7. The cone is drawn in the
coordinates of the sprinkling, so the rule is frame-dependent; in the
theory the metric would narrow the cones itself and the ORDER would change
(the true children of v in a narrowed order need not be links of the
original order). A PASSED result says that cone narrowing of this size can
hold the phase together; the order dynamics that would produce it remains
to be built (F6). Not in the model: the field Pi of equation (a), any
other feedback on the routing.

Arms (all on the same sprinklings and the same drive)
-----------------------------------------------------
  control   no pulse, no channel, no narrowing (route_capacity, F3.2)
  bistable  pulse + threshold + capacity channel, full cones (F3.3g)
  narrow    pulse + threshold + narrowed cones, Pi from the response of
            axiom 5 without acting on the capacity (the lens alone)
  focused   pulse + threshold + capacity channel + narrowed cones (decides)

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

Part A: does a boundary exist without the field? (decides)
-----------------------------------------------------------
Ignition disk radius r0 = 0.15; background grid declared now, all values
on the same sprinklings in one part (the sprinkling is the expensive step):
s_bg = 0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8 (no drive ... near the
homogeneous drive that flooded in F3.3f). In a bistable medium the front
between the phases retreats on one side of a critical background and
advances on the other (F3.3g without narrowing: dies at <= 0.25, floods at >= 0.3);
a droplet can only stand near that point, so the grid is dense there.
Physics to keep in mind: in a continuous bistable medium a standing
droplet is usually the critical nucleus, which is unstable (a larger one
grows, a smaller one shrinks); a stable bounded droplet needs something
more, for example fronts pinned by the discreteness of the medium. Part B
(B_small, B_large) is what tells the two apart.
Every grid value gets the decision below; part A PASSES if at least one
value passes (the value is then confirmed on new sprinklings by part B,
whose seed stream is its own). Per sprinkling of the focused arm: died if N_drop(late) < 0.05
N_drop(early); otherwise growth g = R_rms(late) - R_rms(early). Over the
sprinklings (bootstrap SE):
  PASSED      at most 5 % died, |mean g| + 2 SE <= r0 / 4, the ring stays
              below 2 x control (no spontaneous ignition), and the droplet is a
              DOMAIN of the saturated phase: the saturated fraction inside R_rms
              in the late slabs C_drop(late) - 2 SE >= 0.5 (most of its events in
              that phase; added for F3.3h before its first full run, because
              F3.3g A showed that the phase can survive as a sparse cloud)
  FAILED      more than 50 % died (the droplet dies), or C_drop(late) + 2 SE
              < 0.5 (a sparse cloud, not a domain), or mean g - 2 SE > r0/4
              (it grows: flood), or mean g + 2 SE < -r0/4 (it shrinks), or
              the ring >= 2 x control + 0.01 (the background ignites)
  UNDECIDED   otherwise, or the control or the bistable arm (full cones) meets PASSED
              too (then the narrowing is not what holds it)
The tolerance r0 / 4 is the quarter rule of the laboratory.

Part B: does a droplet return after a perturbation? (decides, after A)
-----------------------------------------------------------------------
Only meaningful at a background where part A found a bounded droplet, so
s_bg is a settable parameter of the B parts with no default (it comes from
part A: CHOSEN FROM DATA, stated in the record). Each B part runs two
focused arms on the same sprinklings, the reference droplet (r0 = 0.15)
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
bit for bit, flux balance <= 1e-12 in every arm, the focused arm
bit-identical in the order of a boosted frame (a computational check: the
rule is a causal pass; the narrowed cone itself is frame-dependent, above).

Seen before this test was written: F3.3g A (above), the core check of the
narrowed routing on a small sprinkling (1 500 events: Pi = 1 changes
nothing, balance exact, order-free, 330 fallback events), and the quick
runs of this file (N = 14 000, 6 sprinklings, decide nothing): with the
narrowed cones the ignited phase survived at backgrounds 0.1-0.25 where it
died with full cones, but as a sparse, slowly thinning cloud (R_rms ~0.3,
B 1.0 -> ~0.3, C_drop 0.9 -> 0.06-0.10, ~10 saturated events per slab
late); about 6-27 % of the narrowed events had no linked child inside the
narrowed cone. The criteria are those of F3.3g plus the domain condition
C_drop(late) >= 0.5, added after these quick runs and before any full run,
so that a sparse surviving cloud cannot pass as a droplet; nothing else
was set from data.
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
    "id": "F3.3h",
    "phase": "F3",
    "title": "Do narrowed future cones (c_eff = cΠ²) hold the saturated phase together as a bounded droplet?",
    "theory_part": "self-focusing (the author's hypothesis): slow clocks narrow the causal cone of the descendants; "
                   "the measure of narrowing c_eff = cΠ² from F7",
    "question": "Does the bistable saturated phase (F3.3g), when an event sends flow only to descendants in its "
                "narrowed cone |Δx| <= Π²Δt, persist as a bounded droplet of constant radius (A) and return after a "
                "perturbation (B)?",
    "method": "F3.3g with a single change: route_pulse with narrow_power = 2 (a cone narrowed by the event's own Π; "
              "with no descendant inside, the most timelike descendant); four branches (control, bistable with full "
              "cones, lens only, bistable with narrowed cones = decides); ignition of a disc r0 = 0.15 on a "
              "background s_bg (A: grid 0–0.8, B: s_bg from A)",
    "prediction": "according to the author's hypothesis the narrowing of the cones in the saturated region keeps the "
                  "flow at its world line and the droplet keeps its radius",
    "criterion": "INSTRUMENT: control = route_capacity, balance <= 1e-12, boosted order bitwise identical. A (every "
                 "value of s_bg; the part passes if at least one passes): PASSED <= 5 % of droplets died, |growth of "
                 "R| + 2 SE <= r0/4, the background does not ignite and the droplet is a domain of the phase (late "
                 "C_drop − 2 SE >= 0.5); FAILED > 50 % died, or a sparse cloud (C_drop + 2 SE < 0.5), or "
                 "growth/decline above r0/4, or the background ignites; otherwise UNDECIDED (also when the control "
                 "or the bistable branch with full cones passes). B (s_bg from part A): PASSED the reference meets "
                 "A, <= 5 % of the perturbed died and |R_perturbed − R_reference| + 2 SE <= r0/4; FAILED the "
                 "reference meets A and > 50 % died or the difference − 2 SE > r0/4; otherwise UNDECIDED.",
    "null_model": "the same sprinklings and driving without a pulse, without a channel and without narrowing (F3.2); "
                  "the bistable branch with full cones (F3.3g)",
    "level": "mechanism",
    "meters": ["K.9"],
    "depends_on": ["F3.3g", "F7.1"],
    "params": {"D": 3, "T": 9.0, "R": 0.5, "R_abs": 0.4, "N": 72000, "n_seeds": 128, "n_cal": 8, "q": 0.99,
               "cal_t1": 3.5, "cal_t2": 5.5, "s0": 1.0, "s_ign": 5.0, "t_ign1": 3.0, "t_ign2": 4.0, "r0": 0.15,
               "ign_axes": [0.15, 0.15], "slab": 0.25, "t_obs1": 3.0, "t_obs2": 8.5,
               "early": [4.5, 5.0], "late": [8.0, 8.5], "ring": 0.3, "die_frac": 0.05,
               "cut": None, "omega": 1.0, "rule": "main", "dq_crit": 0.2, "channel": "capacity", "Pi_sat": 0.5, "m": 2.0,
               "narrow_power": 2.0, "c_domain": 0.5, "max_died_pass": 0.05, "max_died_fail": 0.5, "n_se": 2.0, "boot_B": 4000,
               "balance_tol": 1e-12, "boost_beta": 0.6, "s_bg": None},
    "parts": {"A": {"mode": "A", "bg_grid": [0.0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8]},
              "B_small": {"mode": "B", "pert": {"ign_axes": [0.10, 0.10]}},
              "B_large": {"mode": "B", "pert": {"ign_axes": [0.22, 0.22]}},
              "B_asym": {"mode": "B", "pert": {"ign_axes": [0.21, 0.107]}},
              "B_cut": {"mode": "B", "pert": {"cut": [0.0, 0.25, 5.0, 5.25]}}},
    "settable": {"B_small": ["s_bg"], "B_large": ["s_bg"], "B_asym": ["s_bg"], "B_cut": ["s_bg"]},
    "quick": {"N": 7000, "n_seeds": 4, "n_cal": 1, "boot_B": 200},
    "estimate": "A ~1.5–2 h with 10 processes (10 background values × 4 branches per sprinkling, geometry N = 72 "
                "000), B ~50–60 min",
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
    C_late = _nanmean(np.asarray(ser["C_drop"], float)[l_]) if not died else float("nan")
    return {"N_early": Ne, "N_late": Nl, "R_early": Re, "R_late": Rl, "died": bool(died), "C_late": C_late,
            "growth": (Rl - Re) if not died else float("nan"),
            "ring_late": _nanmean(np.asarray(ser["ring"], float)[l_])}


def _nanmean(v):
    """nanmean without the empty-slice warning (nan if nothing finite)."""
    v = np.asarray(v, float)
    return float(v[np.isfinite(v)].mean()) if np.isfinite(v).any() else float("nan")


def run_arm(ptr, idx, s, C, absorb, p, arm, X, order=None):
    """route_pulse for one arm: control, bistable (full cones), narrow (lens alone), focused[_ref] (bistable + lens)."""
    om = float(p["omega"])
    if arm == "control":
        return route_pulse(ptr, idx, s, C, absorb, "none", om, order=order)
    ch = "capacity" if (arm == "bistable" or arm.startswith("focused")) else "none"
    npow = float(p["narrow_power"]) if (arm == "narrow" or arm.startswith("focused")) else 0.0
    return route_pulse(ptr, idx, s, C, absorb, str(p["rule"]), om, order=order, dq_crit=float(p["dq_crit"]), channel=ch,
                       Pi_sat=float(p["Pi_sat"]), m=float(p["m"]), narrow_power=npow, X=X if npow > 0 else None)


def configs(p):
    """(label, s_bg, {arm: (ignition axes, cut)}) for every background of the part."""
    ref_axes = p["ign_axes"]
    if p["mode"] == "A":
        arms = {"control": (ref_axes, None), "bistable": (ref_axes, None), "narrow": (ref_axes, None),
                "focused": (ref_axes, None)}
        return [(f"s_bg {float(b):g}", float(b), arms) for b in p["bg_grid"]]
    pert = p["pert"]
    pax, pcut = pert.get("ign_axes", ref_axes), pert.get("cut")
    return [(f"s_bg {float(p['s_bg']):g}", float(p["s_bg"]),
             {"control": (pax, pcut), "focused_ref": (ref_axes, None), "focused": (pax, pcut)})]


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
            res = run_arm(ptr, idx, s, C, ab, p, arm, S.X)
            out["balance"] = max(out["balance"], pulse_balance(res))
            ser = slab_series(S, r, t, res, p)
            one[arm] = {"series": ser, **window_stats(ser, p)}
            ev = (res["n_allowed"] >= 0) & (r < float(p["R_abs"])) & (t >= float(p["t_obs1"])) & (t < float(p["t_obs2"]))
            one[arm]["fallback_frac"] = float(res["narrow_fallback"][ev].mean()) if ev.any() else float("nan")
            one[arm]["narrowed_frac"] = float((res["Pi"][ev] < 1.0).mean()) if ev.any() else float("nan")
            if arm == "control":
                ref = route_capacity(ptr, idx, s, C, ab)
                out["control_equals_capacity"] &= all(np.array_equal(res[k], ref[k]) for k in ("accepted", "stored", "Phi"))
            if arm == "focused":
                b = run_arm(ptr, idx, s, C, ab, p, arm, S.X, order=border)
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
        fb = float(np.nanmean([r_[arm]["fallback_frac"] for r_ in rows]))
        nf = float(np.nanmean([r_[arm]["narrowed_frac"] for r_ in rows]))
    cl = np.array([r_[arm]["C_late"] for r_ in rows])
    cs = boot_se(cl[np.isfinite(cl)], B, bs("c", arm)) if np.isfinite(cl).sum() >= 2 else (float("nan"), float("nan"), None)
    return {"died_frac": float(died.mean()), "growth": gs, "ring_late": ring, "series": mean_ser,
            "R_early": R_early, "R_late": R_late, "fallback_frac": fb, "narrowed_frac": nf, "C_late": cs}


def verdict_A(st, ctrl, p):
    """The part-A decision for one arm (st) against the control arm (ctrl)."""
    r0, nse = float(p["r0"]), float(p["n_se"])
    g, se = st["growth"][0], st["growth"][1]
    ring_bad = np.isfinite(st["ring_late"]) and st["ring_late"] >= 2.0 * ctrl["ring_late"] + 0.01
    ring_ok = np.isfinite(st["ring_late"]) and st["ring_late"] < 2.0 * ctrl["ring_late"] + 0.01
    c, cse = st["C_late"][0], st["C_late"][1]
    cmin = float(p["c_domain"])
    if st["died_frac"] > float(p["max_died_fail"]) or ring_bad:
        return "FAILED"
    if np.isfinite(c) and np.isfinite(cse) and c + nse * cse < cmin:
        return "FAILED"                              # a sparse cloud, not a domain of the saturated phase
    if np.isfinite(g) and np.isfinite(se):
        if g - nse * se > r0 / 4 or g + nse * se < -r0 / 4:
            return "FAILED"
        if (st["died_frac"] <= float(p["max_died_pass"]) and abs(g) + nse * se <= r0 / 4 and ring_ok
                and np.isfinite(c) and c - nse * cse >= cmin):
            return "PASSED"
    return "UNDECIDED"


def series_lines(name, st):
    s = st["series"]
    f = lambda v, fmt: " ".join(fmt.format(x) if np.isfinite(x) else "  nan" for x in v)
    return [f"     {name}: died {st['died_frac']:.3f}; R_rms early {st['R_early']:.3f} late {st['R_late']:.3f}; growth "
            + (f"{st['growth'][0]:+.4f} +- {st['growth'][1]:.4f}" if np.isfinite(st["growth"][0]) else "undefined")
            + f"; ring (late) {st['ring_late']:.4f}"
            + (f"; C_drop (late) {st['C_late'][0]:.3f} +- {st['C_late'][1]:.3f}" if np.isfinite(st["C_late"][0]) else "")
            + (f"; narrowed cones: events with Pi < 1 {st['narrowed_frac']:.3f}, of which forwarding to the most timelike "
               f"child {st['fallback_frac']:.3f}" if np.isfinite(st["fallback_frac"]) else ""),
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
        raise SystemExit("F3.3h part B: give the background where part A found a bounded droplet, e.g. --set s_bg=0.3")
    n, B = int(p["n_seeds"]), int(p["boot_B"])
    t0 = time.time()
    cal = parallel_map(calibration_task, [(p, derive_seed(ctx.seed, "F3.3h-cal", int(p["N"]), k)) for k in range(int(p["n_cal"]))],
                       ctx.workers, ctx.progress, lambda i, r_, d, nn: f"[F3.3h] calibration {d}/{nn} ({time.time()-t0:.0f}s)")
    C = float(np.quantile(np.concatenate(cal), float(p["q"])))
    stream = "F3.3h-" + p["mode"]                        # part B confirms on sprinklings part A never saw
    seeds = [derive_seed(ctx.seed, stream, int(p["N"]), k) for k in range(n)]
    cf = configs(p)

    def prog(i, r_, d, nn):
        bits = " ".join(f"{lab.split()[-1]}:{c['focused']['N_late']:.0f}" for (lab, _, _), c in zip(cf, r_["cfg"]))
        return f"[F3.3h] seed {d}/{nn}: focused N late per s_bg {bits} ({time.time()-t0:.0f}s)"

    rows = parallel_map(seed_task, [(p, sd, C) for sd in seeds], ctx.workers, ctx.progress, prog)
    inst = {"control arm = route_capacity": all(r_["control_equals_capacity"] for r_ in rows),
            "balance": max(r_["balance"] for r_ in rows) <= float(p["balance_tol"]),
            "boosted order bit-identical": all(r_["boost_identical"] for r_ in rows)}
    r0, nse = float(p["r0"]), float(p["n_se"])
    head = (f"2+1 cylinder T = {p['T']}, N = {p['N']}, R_abs = {p['R_abs']}; C = Q_{p['q']} = {C:.4f} ({p['n_cal']} calibration "
            f"sprinklings, own stream); ignition s_ign = {p['s_ign']} in t [{p['t_ign1']}, {p['t_ign2']}); pulse '{p['rule']}', "
            f"dQ_crit = {p['dq_crit']}, channel '{p['channel']}', Pi_sat = {p['Pi_sat']}, m = {p['m']}; narrowed cones "
            f"|dx| <= Pi^{p['narrow_power']:g} dt (frame of the sprinkling); {n} sprinklings "
            f"(stream {stream}); slabs of {p['slab']} from t = {p['t_obs1']} (early {p['early']}, late {p['late']}); "
            f"tolerance r0/4 = {r0 / 4:.4f}")
    lines, summaries, verdicts = [head], [], []
    for j, (label, s_bg, arms) in enumerate(cf):
        sub = [{**r_["cfg"][j]} for r_ in rows]
        bs = lambda *tag, j=j: derive_seed(ctx.seed, "F3.3h-boot", ctx.part, j, *tag)
        st = {a: summarise(sub, a, p, B, bs) for a in arms}
        lines.append(f"  {label}:")
        for a in arms:
            lines += series_lines(a, st[a])
        if p["mode"] == "A":
            v_f = verdict_A(st["focused"], st["control"], p)
            v_c = verdict_A(st["control"], st["control"], p)
            v_b = verdict_A(st["bistable"], st["control"], p)
            v_n = verdict_A(st["narrow"], st["control"], p)
            v = "UNDECIDED" if (v_f == "PASSED" and (v_c == "PASSED" or v_b == "PASSED")) else v_f
            lines.append(f"     -> focused arm {v_f}; control {v_c}; bistable with full cones {v_b}; lens alone {v_n} "
                         f"(reported); value {v}")
        else:
            ref_ok = verdict_A(st["focused_ref"], st["control"], p) == "PASSED"
            d = np.array([c["focused"]["R_late"] - c["focused_ref"]["R_late"]
                          if not (c["focused"]["died"] or c["focused_ref"]["died"]) else float("nan") for c in sub])
            ds = boot_se(d[np.isfinite(d)], B, bs("d")) if np.isfinite(d).sum() >= 2 else (float("nan"), float("nan"), None)
            died = st["focused"]["died_frac"]
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
