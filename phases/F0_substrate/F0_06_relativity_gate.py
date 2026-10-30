"""
F0.6 -- The relativity gate: do the local rules of the laboratory pick a preferred frame?
========================================================================================

Why this test exists (22 Oct 2026)
----------------------------------
Lorentz invariance is tested in nature to many orders of magnitude. A
discrete model most often dies on it: a rule evaluated "slab by slab" in
the coordinates of the computation quietly picks a preferred frame. The
author asked for a gate that every rule of the substrate must pass before
its consequences are interpreted as physics (the first gate on the way to
the universality-of-free-fall test: relativity, existence, variety, fall).

Principle
---------
The physics of a moving object must be the Lorentz transform of the
physics of the same object at rest. In the laboratory this can be tested
directly: take ONE sprinkling (the same events and the same causal order)
and let the rule draw its geometry in the coordinates of a frame S'
moving with velocity beta along x (t' = g(t - beta x), x' = g(x - beta t),
y' = y). Everything else stays in the frame S of the source: the drive
(ignition disk, background), the absorbing ring, the cylinder, and the
readout of the observables. Run (S) is "source at rest in the frame of
the rule"; run (S') is "source moving with -beta in the frame of the
rule", read back in the source frame. A covariant rule gives the same
physics in both (bit for bit, since the events are the same); a rule with
a preferred frame does not.

Rules that use only the order (links, counts, flux from parents: capacity
routing, the copied pulse, the threshold, the capacity channel) cannot see
the frame. Rules that use coordinates can:
  * the narrowed future cones of F3.3h, |x_k - x_v| <= Pi_v^2 (t_k - t_v),
    and their fallback "most timelike child" (smallest |dx|/dt);
  * the local relation of the dynamic order of F6.3b,
    z R x iff |x_x - x_z| <= Pi_z^2 (t_x - t_z)  (a different ORDER in S').
The laboratory already calls them "frame-dependent stand-ins". This test
measures HOW frame-dependent they are, against the size of what they do.

Parts
-----
fixed    the fixed order of a sprinkling (core.order.links), the F3.3h model:
         2+1 cylinder T = 9, radius 0.5, absorbed beyond R_abs = 0.4,
         N = 72 000, C = Q_0.99 at unit drive (the F3 capacity definition),
         ignition s_ign = 5 in the disk r0 = 0.15 during t [3, 4).
         Arms:
           bistable   pulse + threshold dQ_crit 0.2 + capacity channel C Pi,
                      full cones (F3.3g): the reference, no coordinates;
           focused    the same with narrowed cones (F3.3h, power 2):
                      run in S and in S' for each beta;
           light      positive control: focused with Pi_sat = 1, so every
                      narrowed cone is the light cone, an invariant set;
                      S and S' (largest beta) must agree bit for bit.
         Backgrounds s_bg = 0 (decision 3 of the author: the base) and 0.2
         (the validation background of F6.3b, where the rules act).
dynamic  the dynamic order of F6.3b (self-consistent order and flux, exact
         construction core.order.dynamic, burn-in 8, average over 4),
         N = 36 000, s_bg = 0.2. Arms: fixed (F3.3g bistable on the vacuum
         order: the reference) and dynamic, built in S and in S' for each
         beta (the local relation drawn in S'; the time window of a link
         scaled by g (1 + beta), the largest dilation of a link's duration).
         Instrument: the vacuum order (Pi = 1) built in S' equals the one
         built in S link for link (the flat relation is invariant).

Observables (per sprinkling and arm, read in S; the definitions of F3.3g/h):
N_late (saturated events per late slab), R_late (R_rms of the saturated
events, late slabs), C_late (saturated fraction inside R_rms, late),
X_late (mean x of the saturated events, late slabs: a drift along the
boost axis is the plainest signature of a preferred frame). Late slabs
t in [8.0, 8.5); early [4.5, 5.0) reported. Died: N_late < 0.05 N_early
(then R, C, X are undefined for that sprinkling).

Criterion (stated before any full run)
--------------------------------------
For each part, background, beta and observable O, over the sprinklings
(paired, bootstrap SE, 2 SE):
  frame effect      Delta = O(S') - O(S) of the coordinate arm (focused / dynamic);
  mechanism effect  d_mech = O(coordinate arm, S) - O(reference arm, S).
Tolerance (the quarter rule of the laboratory):
  N_late, R_late, C_late   |d_mech| / 4, used only where the mechanism is
                           resolved (|d_mech| - 2 SE > 0); a frame effect as
                           large as a quarter of what the rule does means the
                           result of the rule is not physics;
  X_late                   r0 / 4 (the radius tolerance of the F3 series;
                           by symmetry the mechanism does not move the
                           centroid, so d_mech is no scale for it).
Each decisive comparison is EQUIVALENT (|Delta| + 2 SE <= tol), EXCESS
(|Delta| - 2 SE > tol) or INCONCLUSIVE. Per part:
  PASSED     every decisive comparison EQUIVALENT (and at least one exists):
             the rule is covariant within the resolution;
  FAILED     at least one EXCESS: the rule picks a preferred frame at a size
             that matters against what it does;
  UNDECIDED  otherwise (no decisive comparison, or some INCONCLUSIVE).
INSTRUMENT_FAILED if a positive control fails: fixed: the light arm is not
bit-identical in S and S', or the flux balance exceeds 1e-12; dynamic: the
vacuum order in S' differs from S, a link goes backwards in t, or the
balance exceeds 1e-12.

Reading declared now. FAILED for the narrowed cones does not refute the
idea that a slow clock narrows the cones; it says that this coordinate
stand-in cannot carry it, and that results built on it (F3.3h, F6.3b) are
not yet physics. A covariant version exists in principle (the cone
narrowed relative to a local rest frame of the flux, not of the
laboratory), but it is a new hypothesis that may only be introduced after
this result, with its own criteria (R132).

Seen before this test was written: F3.3g A, the quick runs of F3.3h and
the development runs of F6.3 / F6.3b (in the frame of the sprinkling
only). No run in a moving frame had been made. After the criteria above
were written, the quick runs of this file (decide nothing; fixed N = 7 000,
dynamic N = 4 000, 4 sprinklings, burn-in 2): with the rule drawn in a
moving frame the saturated cloud that survives in S mostly dies
(N_late 5.6 -> 0.2 fixed, 16.9 -> 2.8 / 0.4 dynamic at s_bg = 0.2) and what
remains drifts along the boost axis (X_late +0.30 / +0.38 fixed); the light
control and the vacuum order were identical in both frames. The criteria
were not changed after these runs.

Note on R_late and C_late: their tolerance needs the reference arm alive
(d_mech); where the reference dies on every sprinkling (the full-cone
bistable arm does at these backgrounds, F3.3g), they are reported, not
decisive, and N_late and X_late decide.
"""

from __future__ import annotations

import math
import time
import warnings

import numpy as np

from core.flux.capacity import route_capacity
from core.flux.pulse import pulse_balance, route_pulse
from core.order.dynamic import link_set, links_increase_time, order_links
from core.order.links import links_csr
from core.order.sprinkle import Cylinder, derive_seed, sprinkle
from core.stats.parallel import parallel_map
from core.stats.stats import boot_se
from core.testing import Result

TEST = {
    "id": "F0.6",
    "phase": "F0",
    "title": "Relativity gate: do the local rules of the laboratory select a preferred frame?",
    "theory_part": "Lorentz covariance of the substrate rules (narrowed cones F3.3h, local relation of the dynamic "
                   "order F6.3b); a gate before their results are read as physics",
    "question": "Does a rule that draws cones in coordinates give the same physics when it is evaluated in a frame "
                "moving at velocity β (the same events, the same order, driving and readout in the rest frame of the "
                "source)? Is any difference small compared with what the rule itself does?",
    "method": "fixed order (F3.3h, N = 72 000): the bistable branch with full cones (reference), narrowed cones in S "
              "and in S' (β = 0.3; 0.6), positive control Π_sat = 1; dynamic order (F6.3b, N = 36 000): fixed "
              "(reference) and dynamic built in S and in S'; observables N, R_rms, C_drop and centroid x in the late "
              "layers, read in S; paired differences per sprinkling",
    "prediction": "rules from the order alone are covariant by construction; rules with coordinates probably not",
    "criterion": "For every part, background, β and quantity: frame effect Δ = O(S') − O(S), mechanism effect d_mech "
                 "= O(branch with coordinates, S) − O(reference, S). Tolerance: N, R, C |d_mech|/4 (only where the "
                 "mechanism is resolved, |d_mech| − 2 SE > 0); centroid x r0/4. EQUIVALENT |Δ| + 2 SE <= tol, "
                 "EXCEEDED |Δ| − 2 SE > tol. PASSED: all deciding comparisons equivalent (at least one exists); "
                 "FAILED: at least one exceeded; otherwise UNDECIDED. INSTRUMENT: positive control (Π_sat = 1 in S "
                 "and S' bitwise identical; vacuum order in S' = in S), balance <= 1e-12, links forward in t.",
    "null_model": "positive control: a cone narrowed to the light cone (Π_sat = 1) is an invariant set and must give "
                  "bitwise the same",
    "level": "mechanism",
    "meters": ["K.9"],
    "depends_on": ["F3.3h", "F6.3b"],
    "params": {"D": 3, "T": 9.0, "R": 0.5, "R_abs": 0.4, "n_cal": 8, "q": 0.99, "cal_t1": 3.5, "cal_t2": 5.5, "s0": 1.0,
               "s_ign": 5.0, "t_ign1": 3.0, "t_ign2": 4.0, "r0": 0.15, "ign_axes": [0.15, 0.15], "slab": 0.25,
               "t_obs1": 3.0, "t_obs2": 8.5, "early": [4.5, 5.0], "late": [8.0, 8.5], "ring": 0.3, "die_frac": 0.05,
               "omega": 1.0, "rule": "main", "dq_crit": 0.2, "Pi_sat": 0.5, "m": 2.0, "power": 2.0, "lam": 25.0,
               "window": 1.5, "damp": 0.5, "burn_in": 8, "n_avg": 4, "betas": [0.3, 0.6], "n_se": 2.0, "quarter": 4.0,
               "boot_B": 4000, "balance_tol": 1e-12},
    "parts": {"fixed": {"mode": "fixed", "N": 72000, "n_seeds": 128, "bg_grid": [0.0, 0.2]},
              "dynamic": {"mode": "dynamic", "N": 36000, "n_seeds": 24, "bg_grid": [0.2]}},
    "quick": {"n_seeds": 4, "n_cal": 1, "boot_B": 200, "burn_in": 2, "n_avg": 2},
    "quick_parts": {"fixed": {"N": 7000}, "dynamic": {"N": 4000}},
    "estimate": "fixed ~30–40 min with 10 processes (12 flow passes per sprinkling); dynamic ~4–6 h with 10 "
                "processes (3 frames × 12 builds of the order per sprinkling, ~100–200 s per build at N = 36 000; "
                "memory ~200 MB per process)",
}

OBS = ("N_late", "R_late", "C_late", "X_late")
SER_KEYS = ("t", "N", "R_rms", "C_drop", "X", "ring")


# ============================================================================================== frames

def boost(X: np.ndarray, beta: float) -> np.ndarray:
    """Coordinates of the events in the frame moving with velocity beta along axis 1 (c = 1)."""
    if beta == 0.0:
        return X.copy()
    g = 1.0 / math.sqrt(1.0 - beta * beta)
    Y = X.copy()
    Y[:, 0] = g * (X[:, 0] - beta * X[:, 1])
    Y[:, 1] = g * (X[:, 1] - beta * X[:, 0])
    return Y


def order_links_in_frame(X: np.ndarray, Pi: np.ndarray, beta: float, power: float, window: float):
    """
    Links of the dynamic order whose local relation is drawn in the frame beta (core.order.dynamic.order_links on
    the boosted coordinates, sorted by the new time), returned as a CSR of parents per child in the ORIGINAL event
    numbering. Every pair of the relation is inside the flat light cone (Pi <= 1), so the original time order is
    still a linear extension. The window (the longest link duration seen) is scaled by g (1 + beta), the largest
    dilation of a duration under the boost.
    """
    if beta == 0.0:
        return order_links(X, Pi, power, window)
    g = 1.0 / math.sqrt(1.0 - beta * beta)
    Xb = boost(X, beta)
    perm = np.argsort(Xb[:, 0], kind="stable")
    pb, ib = order_links(Xb[perm], Pi[perm], power, window * g * (1.0 + beta))
    child_new = np.repeat(np.arange(len(X)), np.diff(pb))
    child = perm[child_new]
    parent = perm[ib]
    o = np.lexsort((parent, child))                          # by child, parents ascending (as order_links)
    child, parent = child[o], parent[o]
    ptr = np.concatenate([[0], np.cumsum(np.bincount(child, minlength=len(X)))]).astype(np.int64)
    return ptr, parent.astype(np.int64)


# ============================================================================================== model (F3.3g/h, F6.3b)

def calibration_task(p, seed):
    """Uncapped offered flux at unit drive in the F3 analysis window (the capacity definition of the series)."""
    S = sprinkle(Cylinder(int(p["D"]), 0.0, float(p["T"]), float(p["R"])), int(p["N"]), seed)
    if p["mode"] == "dynamic":
        ptr, idx = order_links(S.X, np.ones(S.n), float(p["power"]), float(p["window"]))
    else:
        ptr, idx = links_csr(S.X, S.rho)
    r, t = S.r(), S.X[:, 0]
    win = (t >= float(p["cal_t1"])) & (t < float(p["cal_t2"])) & (r < float(p["R_abs"]))
    return route_capacity(ptr, idx, np.full(S.n, float(p["s0"])), np.inf, r > float(p["R_abs"]))["offered"][win]


def drive(S, t, p, s_bg):
    """Background everywhere plus the ignition in the disk during the ignition window (F3.3g), in the source frame."""
    a, b = (float(v) for v in p["ign_axes"])
    inside = (S.X[:, 1] / a) ** 2 + (S.X[:, 2] / b) ** 2 < 1.0
    on = (t >= float(p["t_ign1"])) & (t < float(p["t_ign2"]))
    return float(s_bg) + float(p["s_ign"]) * (inside & on)


def bistable(ptr, idx, s, C, ab, p, Pi_sat=None, narrow=0.0, X=None):
    """The F3.3g bistable flux (pulse + threshold + capacity channel), optionally with narrowed cones drawn in X."""
    ps = float(p["Pi_sat"]) if Pi_sat is None else float(Pi_sat)
    return route_pulse(ptr, idx, s, C, ab, str(p["rule"]), float(p["omega"]), dq_crit=float(p["dq_crit"]),
                       channel="capacity", Pi_sat=ps, m=float(p["m"]), narrow_power=narrow,
                       X=X if narrow > 0 else None)


def slab_series(S, r, t, res, p):
    """Per slab (read in the source frame S): N_drop, R_rms, C_drop, mean x of the saturated events, ring fraction."""
    Rabs, ring = float(p["R_abs"]), float(p["ring"])
    edges = np.arange(float(p["t_obs1"]), float(p["t_obs2"]) + 1e-9, float(p["slab"]))
    sat = res["stored"] > 0
    out = {k: [] for k in SER_KEYS}
    for k in range(len(edges) - 1):
        m = (r < Rabs) & (t >= edges[k]) & (t < edges[k + 1])
        ms = m & sat
        Nd = int(ms.sum())
        Rrms = float(np.sqrt(2.0 * np.mean(r[ms] ** 2))) if Nd else float("nan")
        inside = m & (r < Rrms) if Nd else np.zeros_like(m)
        ringm = m & (r >= ring)
        out["t"].append(0.5 * (edges[k] + edges[k + 1]))
        out["N"].append(Nd)
        out["R_rms"].append(Rrms)
        out["C_drop"].append(float(sat[inside].mean()) if inside.any() else float("nan"))
        out["X"].append(float(S.X[ms, 1].mean()) if Nd else float("nan"))
        out["ring"].append(float(sat[ringm].mean()) if ringm.any() else float("nan"))
    return out


def _nanmean(v):
    v = np.asarray(v, float)
    return float(v[np.isfinite(v)].mean()) if np.isfinite(v).any() else float("nan")


def average_series(series_list):
    """Element-wise mean over the iterations of the self-consistent ensemble (F6.3b)."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return {k: list(np.nanmean(np.array([s[k] for s in series_list], float), axis=0)) for k in SER_KEYS}


def observables(ser, p):
    """N_early, N_late, R_late, C_late, X_late, died (F3.3g definitions; R, C, X undefined if died)."""
    t = np.asarray(ser["t"])
    e = (t >= p["early"][0]) & (t < p["early"][1])
    l_ = (t >= p["late"][0]) & (t < p["late"][1])
    N = np.asarray(ser["N"], float)
    Ne, Nl = float(N[e].mean()), float(N[l_].mean())
    died = (Ne == 0.0) or (Nl < float(p["die_frac"]) * Ne)
    nan = float("nan")
    return {"N_early": Ne, "N_late": Nl, "died": bool(died),
            "R_late": nan if died else _nanmean(np.asarray(ser["R_rms"], float)[l_]),
            "C_late": nan if died else _nanmean(np.asarray(ser["C_drop"], float)[l_]),
            "X_late": nan if died else _nanmean(np.asarray(ser["X"], float)[l_])}


# ============================================================================================== tasks

def fixed_task(p, seed, C):
    """Fixed order: per background the reference, focused in S and S' (every beta), light control in S and S'."""
    S = sprinkle(Cylinder(int(p["D"]), 0.0, float(p["T"]), float(p["R"])), int(p["N"]), seed)
    ptr, idx = links_csr(S.X, S.rho)
    r, t = S.r(), S.X[:, 0]
    ab = r > float(p["R_abs"])
    betas = [float(b) for b in p["betas"]]
    frames = {0.0: S.X, **{b: boost(S.X, b) for b in betas}}
    npow = float(p["power"])
    out = {"cfg": [], "light_identical": True, "balance": 0.0, "frame_changed_kids": {}}
    for s_bg in p["bg_grid"]:
        s = drive(S, t, p, s_bg)
        one = {}
        ref = bistable(ptr, idx, s, C, ab, p)
        out["balance"] = max(out["balance"], pulse_balance(ref))
        one["reference"] = observables(slab_series(S, r, t, ref, p), p)
        base = None
        for b, Xf in frames.items():
            res = bistable(ptr, idx, s, C, ab, p, narrow=npow, X=Xf)
            out["balance"] = max(out["balance"], pulse_balance(res))
            one[f"coord_{b:g}"] = observables(slab_series(S, r, t, res, p), p)
            if b == 0.0:
                base = res
            else:                                            # how many events have a different number of allowed children
                ev = (base["n_allowed"] >= 0)
                one[f"diff_events_{b:g}"] = float((base["n_allowed"][ev] != res["n_allowed"][ev]).mean()) if ev.any() else 0.0
        bmax = max(betas)
        l0 = bistable(ptr, idx, s, C, ab, p, Pi_sat=1.0, narrow=npow, X=S.X)
        l1 = bistable(ptr, idx, s, C, ab, p, Pi_sat=1.0, narrow=npow, X=frames[bmax])
        out["light_identical"] &= all(np.array_equal(l0[k], l1[k]) for k in ("accepted", "stored", "Phi", "P", "Pi"))
        out["cfg"].append(one)
    return out


def solve_dynamic_in_frame(S, r, t, s, C, ab, p, beta):
    """F6.3b self-consistent order and flux with the local relation drawn in the frame beta (average after burn-in)."""
    n_burn, n_avg, d = int(p["burn_in"]), int(p["n_avg"]), float(p["damp"])
    Pi = np.ones(S.n)
    series, bal, acyc = [], 0.0, True
    for it in range(n_burn + n_avg):
        ptr, idx = order_links_in_frame(S.X, Pi.copy(), beta, float(p["power"]), float(p["window"]))
        acyc &= links_increase_time(S.X, ptr, idx)
        res = bistable(ptr, idx, s, C, ab, p)
        bal = max(bal, pulse_balance(res))
        if it >= n_burn:
            series.append(slab_series(S, r, t, res, p))
        change = np.abs(res["Pi"] - Pi)
        if float(change.max()) == 0.0:
            if not series:
                series.append(slab_series(S, r, t, res, p))
            break
        Pi = (1.0 - d) * Pi + d * res["Pi"]
    return average_series(series), bal, acyc


def dynamic_task(p, seed, C):
    """Dynamic order: vacuum invariance check, the fixed reference, the dynamic arm in S and in S' (every beta)."""
    S = sprinkle(Cylinder(int(p["D"]), 0.0, float(p["T"]), float(p["R"])), int(p["N"]), seed)
    r, t = S.r(), S.X[:, 0]
    ab = r > float(p["R_abs"])
    betas = [float(b) for b in p["betas"]]
    ptr0, idx0 = order_links(S.X, np.ones(S.n), float(p["power"]), float(p["window"]))
    pv, iv = order_links_in_frame(S.X, np.ones(S.n), max(betas), float(p["power"]), float(p["window"]))
    out = {"cfg": [], "vacuum_identical": link_set(ptr0, idx0) == link_set(pv, iv), "balance": 0.0, "acyclic": True}
    for s_bg in p["bg_grid"]:
        s = drive(S, t, p, s_bg)
        one = {}
        ref = bistable(ptr0, idx0, s, C, ab, p)
        out["balance"] = max(out["balance"], pulse_balance(ref))
        one["reference"] = observables(slab_series(S, r, t, ref, p), p)
        for b in [0.0] + betas:
            ser, bal, acyc = solve_dynamic_in_frame(S, r, t, s, C, ab, p, b)
            out["balance"] = max(out["balance"], bal)
            out["acyclic"] &= acyc
            one[f"coord_{b:g}"] = observables(ser, p)
        out["cfg"].append(one)
    return out


# ============================================================================================== statistics

def compare(rows, j, beta, key, p, bs):
    """Paired frame effect and mechanism effect of one observable at one background and beta."""
    def pairs(a, b):
        d = np.array([r_["cfg"][j][a][key] - r_["cfg"][j][b][key] for r_ in rows], float)
        return d[np.isfinite(d)]
    B, nse = int(p["boot_B"]), float(p["n_se"])
    fr = pairs(f"coord_{beta:g}", "coord_0")
    me = pairs("coord_0", "reference")
    nan2 = (float("nan"), float("nan"))
    f_m, f_se = boot_se(fr, B, bs("frame", key, beta))[:2] if fr.size >= 2 else nan2
    m_m, m_se = boot_se(me, B, bs("mech", key))[:2] if me.size >= 2 else nan2
    if key == "X_late":
        tol, decisive = float(p["r0"]) / float(p["quarter"]), np.isfinite(f_m) and np.isfinite(f_se)
    else:
        resolved = np.isfinite(m_m) and np.isfinite(m_se) and abs(m_m) - nse * m_se > 0
        tol = abs(m_m) / float(p["quarter"]) if resolved else float("nan")
        decisive = resolved and np.isfinite(f_m) and np.isfinite(f_se)
    if not decisive:
        verdict = "NOT DECISIVE"
    elif abs(f_m) + nse * f_se <= tol:
        verdict = "EQUIVALENT"
    elif abs(f_m) - nse * f_se > tol:
        verdict = "EXCESS"
    else:
        verdict = "INCONCLUSIVE"
    return {"key": key, "beta": beta, "frame": (f_m, f_se), "n_frame": int(fr.size), "mech": (m_m, m_se),
            "n_mech": int(me.size), "tol": tol, "verdict": verdict}


def arm_means(rows, j, arm):
    """Mean observables and died fraction of one arm (report)."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = {k: float(np.nanmean([r_["cfg"][j][arm][k] for r_ in rows])) for k in ("N_early",) + OBS}
    m["died"] = float(np.mean([r_["cfg"][j][arm]["died"] for r_ in rows]))
    return m


# ============================================================================================== run

def run(p, ctx):
    mode = str(p["mode"])
    t0 = time.time()
    cal = parallel_map(calibration_task, [(p, derive_seed(ctx.seed, "F0.6-cal", mode, int(p["N"]), k))
                                          for k in range(int(p["n_cal"]))],
                       ctx.workers, ctx.progress, lambda i, r_, d, nn: f"[F0.6] calibration {d}/{nn} ({time.time()-t0:.0f}s)")
    C = float(np.quantile(np.concatenate(cal), float(p["q"])))
    seeds = [derive_seed(ctx.seed, "F0.6", mode, int(p["N"]), k) for k in range(int(p["n_seeds"]))]
    task = fixed_task if mode == "fixed" else dynamic_task
    rows = parallel_map(task, [(p, sd, C) for sd in seeds], ctx.workers, ctx.progress,
                        lambda i, r_, d, nn: f"[F0.6] {mode} sprinkling {d}/{nn} ({time.time()-t0:.0f}s)")
    betas = [float(b) for b in p["betas"]]
    if mode == "fixed":
        inst = {"light control bit-identical in S and S'": all(r_["light_identical"] for r_ in rows),
                "balance": max(r_["balance"] for r_ in rows) <= float(p["balance_tol"])}
        rule = (f"narrowed cones |dx| <= Pi^{p['power']:g} dt (F3.3h) on the fixed order; reference: bistable with full cones")
    else:
        inst = {"vacuum order in S' = in S": all(r_["vacuum_identical"] for r_ in rows),
                "links forward in t": all(r_["acyclic"] for r_ in rows),
                "balance": max(r_["balance"] for r_ in rows) <= float(p["balance_tol"])}
        rule = (f"dynamic order z R x iff |dx| <= Pi_z^{p['power']:g} dt (F6.3b, self-consistent, burn-in {p['burn_in']}, "
                f"average over {p['n_avg']}); reference: bistable on the vacuum order")
    lines = [f"{mode}: 2+1 cylinder T = {p['T']}, N = {p['N']}, R_abs = {p['R_abs']}; C = Q_{p['q']} = {C:.4f}; ignition "
             f"s_ign = {p['s_ign']} in the disk r0 = {p['r0']}, t [{p['t_ign1']}, {p['t_ign2']}); {len(rows)} sprinklings; "
             f"rule frames beta = 0 (source frame) and {betas}; readout in the source frame",
             f"  rule under test: {rule}"]
    comps = []
    for j, s_bg in enumerate(p["bg_grid"]):
        bs = lambda *tag, j=j: derive_seed(ctx.seed, "F0.6-boot", mode, j, *tag)
        lines.append(f"  s_bg {float(s_bg):g}:")
        for arm in ["reference", "coord_0"] + [f"coord_{b:g}" for b in betas]:
            m = arm_means(rows, j, arm)
            lines.append(f"     {arm:10s}: died {m['died']:.3f}; N early {m['N_early']:.1f} late {m['N_late']:.1f}; "
                         f"R_late {m['R_late']:.4f}; C_late {m['C_late']:.3f}; X_late {m['X_late']:+.4f}")
        if mode == "fixed":
            for b in betas:
                fr = float(np.mean([r_["cfg"][j][f"diff_events_{b:g}"] for r_ in rows]))
                lines.append(f"     beta {b:g}: events with a different number of allowed children than in S: {fr:.4f}")
        for b in betas:
            for key in OBS:
                c = compare(rows, j, b, key, p, bs)
                c["s_bg"] = float(s_bg)
                comps.append(c)
                lines.append(f"     beta {b:g} {key:7s}: frame effect {c['frame'][0]:+.4f} +- {c['frame'][1]:.4f} "
                             f"(n {c['n_frame']}); mechanism {c['mech'][0]:+.4f} +- {c['mech'][1]:.4f} (n {c['n_mech']}); "
                             f"tol {c['tol']:.4f} -> {c['verdict']}")
    dec = [c for c in comps if c["verdict"] != "NOT DECISIVE"]
    if not all(inst.values()):
        status = "INSTRUMENT_FAILED"
    elif any(c["verdict"] == "EXCESS" for c in dec):
        status = "FAILED"
    elif dec and all(c["verdict"] == "EQUIVALENT" for c in dec):
        status = "PASSED"
    else:
        status = "UNDECIDED"
    lines.append(f"decision: {len(dec)} decisive comparisons: "
                 + ", ".join(f"{v} {sum(c['verdict'] == v for c in dec)}" for v in ("EQUIVALENT", "EXCESS", "INCONCLUSIVE")))
    lines.append("instrument: " + "; ".join(f"{k}: {'OK' if v else 'FAILS'}" for k, v in inst.items()))
    return Result(status, lines, {"C": C, "comparisons": comps, "instrument": inst, "rows": rows})
