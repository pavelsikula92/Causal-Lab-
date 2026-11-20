"""
F5.1e -- Free fall from the apex in 3+1: does a longest-chain history that starts at the top of a thrown ball fall?
==================================================================================================================

Why this test exists (30 Oct 2026, after F5.1d PASSED in 2+1; the author's priority)
-----------------------------------------------------------------------------------
F5.1d: in the 2+1 throughput order of F7.1b (sigma0 = 0.45) a history that starts at the apex of the longest chain
between two events at the same place, and is then continued by the longest-chain rule of F5.1c (m = H = 8), falls with
R = a / a_pred = 1.069 +- 0.051 (criterion 1/2 < R -+ 2 SE < 2; flat R = 0.146 +- 0.098; longest chains R_max = 0.742 +-
0.075). Reported there: a velocity jump of ~0.09 at the apex and R / R_max = 1.44. 2+1 is not our world: its potential
is not 1/r far away and general relativity has no Newtonian limit in it. This test is the same experiment in 3+1, the
dimension of the throughput order with gamma = 1 (F7.1c raw 1.000 +- 0.045 and 0.972 +- 0.074).

Model (locked before any run)
-----------------------------
The 3+1 throughput order of F7.1c: events of coordinate density exp(sigma), sigma = sigma0 exp(-r^2 / 2 w^2), w = 4;
links of the flat sprinkling (core.order.links, rho_min = 1) kept iff their speed is <= exp(-sigma) at the parent event
(CapacityOrder3 of F7.1c). sigma0 = 0.3 (F7.1c reported 0.924 +- 0.027 and 0.901 +- 0.045 for this strength: the order
is still metric; 0.2 would halve the signal). Flat order: sigma0 = 0.
Everything else is F5.1d scaled to the smaller box that 3+1 can carry:
  launch ring r_launch = 4 (= w, the largest gradient) at t_launch = 1, in the plane z = 0, n_ring = 60 chains;
  apex: the longest chain from the launch event to the event nearest to the same place dT_launch = 6 later; the history
    is its rising half, the apex its middle event (F5.1d, apex_history);
  rule: F5.1c / F5.1d unchanged (walk_history of F5.1d: p = m = 8 steps back; Q = {q : L(v, q) = H = 8};
    q* = argmin L(p, q), ties min n(p, q); one step along a longest chain to q*; rolling horizon), computed on the kept
    links inside a cylinder of radius region = 4 around the start;
  box t in [0, 24], |x|, |y|, |z| <= 8; fit over span = 12 from the apex. The box is sized so that the horizon of H = 8
    links (~7 units of time in this order) still fits above the end of the readout (apex ~4 + 12 + 7 < 24), the lesson of
    F5.1c, and so that one order (~100 000 events) fits in the memory of one process (~1.5 GB).
  16 sprinklings per order (loaded and flat).
Readout and prediction as F5.1d: the displacement from the apex projected on the outward radial unit vector (3D) at
the apex, mean over chains, fitted by v0 tau + a tau^2 / 2. Prediction: with the 3+1 dictionary (K.16c: density
phi - 3 psi, cones phi + psi) the order has phi = psi = -sigma/2, so a slow body falls from rest at r_apex with
a = -dphi/dr (RK4, fall_curve3); a_pred per sprinkling from the same fit; R = a / a_pred. The flat order uses the
a_pred of the loaded geometry for the same apex radii (the null: no fall).
Control (as F5.1): longest chains between events at r0 = 4, t0 = 2 and t0 + dT, dT = 16; their outward sag at
mid-time against the sag of the proper-time-maximising path of phi = -sigma/2 (sag_prediction3): R_max.

Criteria (stated before any run; those of F5.1d unchanged)
---------------------------------------------------------
Classification by R over sprinklings:
  falls as predicted   1/2 < R - 2 SE  and  R + 2 SE < 2;
  falls too strongly   R - 2 SE >= 2;
  does not fall        R + 2 SE < 1/2;
  inconclusive         otherwise.
Instrument (INSTRUMENT_FAILED if any fails): flat |R| + 2 SE <= 1/2 for the walk and for the longest chains; loaded
longest chains 1/2 < R_max -+ 2 SE < 2.
PASSED if the walk falls as predicted; FAILED if it does not fall or falls too strongly; UNDECIDED otherwise.
Reported: R / R_max, v0 and the radial velocity of the rising half (the jump at the apex of F5.1d), the fraction of
chains that reach the end, the displacement series, <v^2> of the links of the walk.

Reading declared now. PASSED: the result of F5.1d is not a property of 2+1; in the order of our dimension a history
continued as a longest chain from a free state falls with the Newtonian acceleration of phi = -sigma/2. FAILED (does
not fall): the rule that falls in 2+1 loses the gravity of the order in 3+1 (shorter chains, more directions; a
property of the rule at m = H = 8). FAILED (falls too strongly): the magnitude is not Newtonian in 3+1.
INSTRUMENT_FAILED: the longest chains do not read the geometry in this box; no conclusion about the rule.
Limits known in advance: the field is Gaussian, not 1/r, so this tests the local acceleration, not Kepler's law; the
source is put in (the order is shaped by sigma, not by flows); one strength, one launch radius.

Development (30 Oct, before any full run; decides nothing)
----------------------------------------------------------
A probe at the first size (box t [0, 20], |x| <= 9, ~123 000 events): one order took 29 min to build in the
container (4.7e7 links, 27 % kept), each walk ~2 s. Too large for 10 processes; the box was reduced.
Two rows at the second size (box t [0, 20], |x| <= 8, span 8, control dT 12; 60 chains; one sprinkling each):
  loaded: N 87 988, links kept 20 %, a = -0.0090 against a_pred = -0.0220 (R 0.41), v0 -0.079 (rising half +0.111),
    displacement from the apex (tau = 0, 1, ..., 8) 0.00 -0.08 -0.20 -0.27 -0.39 -0.49 -0.62 -0.78 -0.92;
    control sag +0.190 against +0.406 (R_max 0.47); <v^2> 0.46; 995 s;
  flat:   N 81 597, a = -0.0022 (R 0.10), v0 -0.064 (rising half +0.066), displacement 0.00 -0.09 -0.20 -0.21 -0.24
    -0.29 -0.44 -0.50 -0.59; control sag +0.098 (R_max 0.24); <v^2> 0.62; 840 s.
One sprinkling each, so nothing is read from R. What they show is that the signal of that size is small against the
discreteness: the predicted fall over 8 is -0.73 and the predicted sag 0.41, while a single sprinkling scatters by ~0.2.
The box was therefore made taller before any full run: span 12 (predicted fall -1.62, 2.2x) and control dT 16
(predicted sag 0.71, 1.75x). Strength, launch radius, rule and criteria are unchanged. Both rows also show an inward
start velocity at the apex (v0 < 0 while the rising half moved outwards), as F5.1d reported in 2+1; it is reported, the
criterion reads the acceleration.
"""

from __future__ import annotations

import math
import time

import numpy as np

from core.order.sprinkle import derive_seed
from core.stats.parallel import parallel_map
from core.stats.stats import boot_se
from core.testing import Result
from phases.F5_gravity.F5_01_free_fall_order import fit_a, longest_chain
from phases.F5_gravity.F5_01d_free_fall_apex import classify_band, walk_history
from phases.F7_light.F7_01c_throughput_cones_3d import CapacityOrder3, sprinkle_density

TEST = {
    "id": "F5.1e",
    "phase": "F5",
    "title": "Free fall from the apex in 3+1: does a longest-chain history that starts at the top of a thrown ball fall?",
    "theory_part": "F5 gravity and equivalence; free fall as inertia plus the geodesic of the order, in 3+1",
    "question": "In the 3+1 throughput order of F7.1c (sigma0 = 0.3), does a chain that starts at the apex of the longest "
                "chain between two events at the same place, with the rising half of that chain as its history, and then "
                "follows the longest-chain rule of F5.1c, fall from the apex as phi = -sigma/2 predicts?",
    "method": "3+1 order of F7.1c, w = 4, box t [0, 24], |x| <= 8; launch ring r = 4 at t = 1 (z = 0); longest chain to the "
              "same place 6 later; start at its middle event with the rising half as history; rule of F5.1c (m = H = 8, "
              "cylinder radius 4); radial projection from the apex over 12, fit a/2 tau^2; prediction: fall from rest at "
              "r_apex in phi = -sigma/2; control: longest chains r0 = 4, t0 = 2, dT = 16; flat null",
    "prediction": "a free state continued as a geodesic of the order falls from the apex with a = -dphi/dr(r_apex), "
                  "phi = -sigma/2",
    "criterion": "falls as predicted 1/2 < R -+ 2 SE < 2; falls too strongly R - 2 SE >= 2; does not fall R + 2 SE < 1/2. "
                 "INSTRUMENT: flat |R| + 2 SE <= 1/2 (walk and longest chains); loaded longest chains "
                 "1/2 < R_max -+ 2 SE < 2. PASSED falls as predicted, FAILED does not fall or falls too strongly, "
                 "UNDECIDED otherwise.",
    "null_model": "flat order (sigma0 = 0): no fall",
    "level": "mechanism",
    "meters": ["K.16c"],
    "depends_on": ["F5.1d", "F7.1c", "K.16c"],
    "params": {"T": 24.0, "half": 8.0, "width": 4.0, "sigma0": 0.3, "r_launch": 4.0, "t_launch": 1.0,
               "dT_launch": 6.0, "span": 12.0, "n_ring": 60, "memory": 8, "horizon": 8, "region": 4.0,
               "r0": 4.0, "t0": 2.0, "dT": 16.0, "n_max": 40, "n_seeds": 16, "n_se": 2.0, "boot_B": 4000},
    "parts": {"main": {}},
    "quick": {"n_seeds": 1, "n_ring": 4, "n_max": 4, "boot_B": 200, "half": 6.0, "T": 14.0, "span": 4.0,
              "r_launch": 3.0, "region": 3.0, "r0": 3.0, "dT": 8.0},
    "estimate": "one order took ~15-20 min per core in the container (~100 000 events); 32 orders; memory ~1.5 GB per "
                "process: use --workers 6 on a 16 GB machine, 10 on 32 GB; the walks are cheap (~2-3 s per chain)",
}


# ===================================================================================== geometry and prediction (3+1)

def accel3(r, p):
    """Coordinate acceleration (radial, outward positive) of a slow body: -dphi/dr with phi = -sigma/2."""
    s0, w = float(p["sigma0"]), float(p["width"])
    return -(s0 * r / (2.0 * w * w)) * math.exp(-r * r / (2.0 * w * w))


def fall_curve3(p, ts):
    """Numerical fall from rest at r0 (RK4), r(t) - r0 at the times ts - t0 (F5.1 fall_curve with accel3)."""
    r, v, t, dt = float(p["r0"]), 0.0, 0.0, 0.01
    out = []
    ts = np.asarray(ts, float) - float(p["t0"])
    for tt in ts:
        while t < tt - 1e-12:
            h = min(dt, tt - t)
            k1v, k1x = accel3(r, p), v
            k2v, k2x = accel3(r + 0.5 * h * k1x, p), v + 0.5 * h * k1v
            k3v, k3x = accel3(r + 0.5 * h * k2x, p), v + 0.5 * h * k2v
            k4v, k4x = accel3(r + h * k3x, p), v + h * k3v
            r += h * (k1x + 2 * k2x + 2 * k3x + k4x) / 6.0
            v += h * (k1v + 2 * k2v + 2 * k3v + k4v) / 6.0
            t += h
        out.append(r - float(p["r0"]))
    return np.array(out)


def sag_prediction3(p):
    """Mid-time outward deviation of the proper-time-maximising path between (t0, r0) and (t0 + dT, r0) in
    phi = -sigma/2: the radial trajectory that returns to r0 at dT (shooting), its deviation at dT/2."""
    dT = float(p["dT"])

    def end(v0):
        r, v, t, h = float(p["r0"]), v0, 0.0, 0.01
        mid = None
        while t < dT - 1e-12:
            v += accel3(r, p) * h
            r += v * h
            t += h
            if mid is None and t >= dT / 2:
                mid = r
        return r - float(p["r0"]), mid - float(p["r0"])
    lo, hi = 0.0, 1.0
    for _ in range(60):
        m = 0.5 * (lo + hi)
        if end(m)[0] > 0:
            hi = m
        else:
            lo = m
    return end(0.5 * (lo + hi))[1]


# ===================================================================================== starts (3+1)

def ring_starts3(X, r0, t0, n):
    """The events nearest to n equally spaced points of the ring r = r0, z = 0 at time t0."""
    out = []
    for j in range(n):
        ang = 2.0 * math.pi * j / n
        target = np.array([t0, r0 * math.cos(ang), r0 * math.sin(ang), 0.0])
        out.append(int(np.argmin(((X - target) ** 2).sum(axis=1))))
    return out


def same_place_target(order, s_, dT):
    """The event nearest to the place of s_ dT later, inside the kept-link future of s_; None if there is none."""
    X = order.X
    tq = X[s_, 0] + dT
    target = X[s_].copy()
    target[0] = tq
    fut = order.future([s_], tq + 1.0)
    cand = np.where(fut & (np.abs(X[:, 0] - tq) < 1.0))[0]
    if cand.size == 0:
        return None
    return int(cand[np.argmin(((X[cand] - target) ** 2).sum(axis=1))])


def apex_history3(order, s_, p):
    """The rising half of the longest chain from s_ to the same place dT_launch later (F5.1d, in any dimension)."""
    q = same_place_target(order, s_, float(p["dT_launch"]))
    if q is None:
        return None
    chain = longest_chain(order, s_, q)
    if chain is None or len(chain) < 3:
        return None
    return chain[: len(chain) // 2 + 1]


# ===================================================================================== one sprinkling

def seed_task(p, sigma0, seed):
    rng = np.random.default_rng(seed)
    T, Hh, w = float(p["T"]), float(p["half"]), float(p["width"])
    X = sprinkle_density(sigma0, w, T, Hh, rng)
    order = CapacityOrder3(X, sigma0, w)
    span = float(p["span"])
    tau = np.linspace(0.0, span, 25)
    cols, preds, full, r_apex, v_apex, v2 = [], [], 0, [], [], []
    for s_ in ring_starts3(X, float(p["r_launch"]), float(p["t_launch"]), int(p["n_ring"])):
        hist = apex_history3(order, s_, p)
        if hist is None:
            continue
        ap = hist[-1]
        ta, ra = float(X[ap, 0]), X[ap, 1:]
        path = walk_history(order, hist, ta + span + 1.0, int(p["memory"]), int(p["horizon"]), float(p["region"]))
        P = X[path[len(hist) - 1:]]                         # from the apex on
        dP = np.diff(P, axis=0)
        if len(dP):
            v2.append(float(np.mean((dP[:, 1:] ** 2).sum(axis=1) / dP[:, 0] ** 2)))   # squared speed of each link
        rhat = ra / np.linalg.norm(ra)
        proj = (P[:, 1:] - ra) @ rhat
        y = np.interp(ta + tau, P[:, 0], proj)
        y[ta + tau > P[-1, 0]] = np.nan                      # a chain that ended counts up to its end
        cols.append(y)
        full += int(P[-1, 0] >= ta + span)
        r_apex.append(float(np.linalg.norm(ra)))
        H0 = X[hist]                                         # radial velocity of the rising half near the apex
        k = max(0, len(hist) - 4)
        v_apex.append(float(((H0[-1, 1:] - H0[k, 1:]) @ rhat) / max(H0[-1, 0] - H0[k, 0], 1e-9)))
        preds.append(fall_curve3(dict(p, r0=float(np.linalg.norm(ra)), t0=0.0), tau))
    with np.errstate(all="ignore"):
        ybar = np.nanmean(cols, axis=0) if cols else np.full(len(tau), np.nan)
    ok = np.isfinite(ybar)
    a, v0 = fit_a(tau[ok], ybar[ok], 0.0) if ok.sum() >= 3 else (float("nan"), float("nan"))
    a_pred, _ = fit_a(tau, np.mean(preds, axis=0), 0.0) if preds else (float("nan"), 0.0)
    out = {"N": len(X), "links_kept": order.n_kept / max(order.n_links, 1), "sigma0": sigma0, "a": a, "v0": v0,
           "a_pred": a_pred, "n_chains": len(cols), "n_full": full, "series": ybar.tolist(),
           "r_apex": float(np.mean(r_apex)) if r_apex else float("nan"),
           "v_apex": float(np.mean(v_apex)) if v_apex else float("nan"), "v2": float(np.mean(v2)) if v2 else float("nan")}
    # control: longest chains between two events at the same place (as F5.1, in 3+1)
    devs = []
    for s_ in ring_starts3(X, float(p["r0"]), float(p["t0"]), int(p["n_max"])):
        q = same_place_target(order, s_, float(p["dT"]))
        if q is None:
            continue
        chain = longest_chain(order, s_, q)
        if chain is None:
            continue
        P = X[chain]
        rhat = P[0, 1:] / np.linalg.norm(P[0, 1:])
        proj = (P[:, 1:] - P[0, 1:]) @ rhat
        tm = 0.5 * (P[0, 0] + P[-1, 0])
        rline = np.interp(tm, [P[0, 0], P[-1, 0]], [proj[0], proj[-1]])
        devs.append(float(np.interp(tm, P[:, 0], proj) - rline))
    out["sag"] = float(np.mean(devs)) if devs else float("nan")
    return out


# ===================================================================================== the test

def run(p, ctx):
    t0 = time.time()
    n, B, nse = int(p["n_seeds"]), int(p["boot_B"]), float(p["n_se"])
    sag_pred = sag_prediction3(p)
    jobs = []
    for s0 in (float(p["sigma0"]), 0.0):
        jobs += [(p, s0, derive_seed(ctx.seed, "F5.1e", ctx.part, s0, k)) for k in range(n)]
    rows = parallel_map(seed_task, jobs, ctx.workers, ctx.progress,
                        lambda i, r_, d, nn: f"[F5.1e] {d}/{nn}: sigma0 {r_['sigma0']:g} a {r_['a']:+.4f} (pred {r_['a_pred']:+.4f}) "
                                             f"sag {r_['sag']:+.3f} ({time.time() - t0:.0f}s)")
    res = {}
    for s0 in (float(p["sigma0"]), 0.0):
        rr = [r_ for r_ in rows if r_["sigma0"] == s0 and np.isfinite(r_["a"])]
        bs = lambda *tag: derive_seed(ctx.seed, "F5.1e-boot", ctx.part, s0, *tag)
        R = boot_se([r_["a"] / r_["a_pred"] for r_ in rr], B, bs("R"))[:2]
        Rm = boot_se([r_["sag"] / sag_pred for r_ in rr if np.isfinite(r_["sag"])], B, bs("Rm"))[:2]
        v0 = boot_se([r_["v0"] for r_ in rr], B, bs("v0"))[:2]
        with np.errstate(all="ignore"):
            ser = np.nanmean([r_["series"] for r_ in rr], axis=0).tolist()
        res[s0] = {"R": R, "R_max": Rm, "v0": v0, "series": ser, "N": float(np.mean([r_["N"] for r_ in rr])),
                   "kept": float(np.mean([r_["links_kept"] for r_ in rr])),
                   "a_pred": float(np.mean([r_["a_pred"] for r_ in rr])), "r_apex": float(np.mean([r_["r_apex"] for r_ in rr])),
                   "v_apex": float(np.mean([r_["v_apex"] for r_ in rr])), "v2": float(np.mean([r_["v2"] for r_ in rr])),
                   "chains": float(np.mean([r_["n_chains"] for r_ in rr])), "full": float(np.mean([r_["n_full"] for r_ in rr]))}
    s0 = float(p["sigma0"])
    lines = [f"3+1 order of F7.1c: sigma0 = {s0:g}, w = {p['width']:g}; box t [0, {p['T']:g}], |x|, |y|, |z| <= {p['half']:g} "
             f"(~{res[s0]['N']:.0f} events, links kept {res[s0]['kept']:.2f}); launch ring r = {p['r_launch']:g} at "
             f"t = {p['t_launch']:g}, {p['n_ring']} chains, apex of the longest chain over dT_launch = {p['dT_launch']:g}; "
             f"longest-chain rule m {p['memory']}, H {p['horizon']}, cylinder {p['region']:g}; fit over {p['span']:g} from the apex",
             f"  prediction (phi = -sigma/2, fall from rest at r_apex): a_pred = {res[s0]['a_pred']:+.5f}; control sag "
             f"{sag_pred:+.4f}"]
    for key, lab in ((s0, "loaded"), (0.0, "flat  ")):
        r_ = res[key]
        lines.append(f"  {lab}: R = {r_['R'][0]:+.3f} +- {r_['R'][1]:.3f} ({classify_band(*r_['R'], nse)}); v0 {r_['v0'][0]:+.4f} +- "
                     f"{r_['v0'][1]:.4f}; apex r {r_['r_apex']:.2f}, radial velocity of the rising half {r_['v_apex']:+.3f}; "
                     f"longest chains R_max = {r_['R_max'][0]:+.3f} +- {r_['R_max'][1]:.3f}; chains {r_['chains']:.1f} of "
                     f"{p['n_ring']}, to the end {r_['full']:.1f}; links <v^2> {r_['v2']:.3f}")
        lines.append("          displacement from the apex: " + " ".join(f"{x:+.2f}" for x in r_["series"][::2]))
    Rl, Rml = res[s0]["R"][0], res[s0]["R_max"][0]
    lines.append(f"  reported: R / R_max = {Rl / Rml:+.3f}" if Rml else "  reported: R_max = 0")
    f_ = res[0.0]
    inst = {"flat: no drift (|R| + 2 SE <= 1/2)": abs(f_["R"][0]) + nse * f_["R"][1] <= 0.5,
            "flat longest chains: no sag (|R_max| + 2 SE <= 1/2)": abs(f_["R_max"][0]) + nse * f_["R_max"][1] <= 0.5,
            "loaded longest chains read the geometry (1/2 < R_max -+ 2 SE < 2)":
                res[s0]["R_max"][0] - nse * res[s0]["R_max"][1] > 0.5 and res[s0]["R_max"][0] + nse * res[s0]["R_max"][1] < 2.0}
    lines.append("instrument: " + "; ".join(f"{k}: {'OK' if ok else 'FAILS'}" for k, ok in inst.items()))
    cls = classify_band(*res[s0]["R"], nse)
    status = "INSTRUMENT_FAILED" if not all(inst.values()) else {"falls as predicted": "PASSED", "does not fall": "FAILED",
                                                                   "falls too strongly": "FAILED"}.get(cls, "UNDECIDED")
    data = {"sag_pred": sag_pred, "res": {str(k): v for k, v in res.items()}, "rows": rows}
    return Result(status, lines, data)
