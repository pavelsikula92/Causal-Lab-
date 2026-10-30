"""
F5.1d -- Free fall from the apex: does a longest-chain history that starts at the top of a thrown ball fall?
===========================================================================================================

Why this test exists (30 Oct 2026, after F5.1c FAILED; its verdict and criteria stay, R132)
-----------------------------------------------------------------------------------------
F5.1c (rule: continue the own history as a longest chain of the order) moved the chains inwards (-4.0 over 22; the
geometry -2.7) but not with the Newtonian acceleration: R = 0.337 +- 0.049. Its displacement shows three phases
(descriptive, after the run): up to t ~ 8 about 2 a_pred, while a VIRTUAL history at rest stood for the past of the
chain; t = 8-16 about 0.5 a_pred; from t ~ 16 a slowdown, where the horizon of H links reaches over the top of the box.
Two recognised distortions, both of the set-up rather than of the rule:
  1. the initial condition: an event at the same place earlier is a body HELD at rest, not a free one; its history is
     not a geodesic of the order, and the rule continues it as if it were;
  2. the top of the box: near t = T the rule cannot see H links ahead.

The initial condition from the order alone (input of this test)
--------------------------------------------------------------
The longest chain between two events at the same place, dT_launch apart, is the geodesic of the order between them:
a ball thrown up that comes back (F5.1 shows these chains bulge outwards as the geometry predicts, R_max = 0.63-0.86).
At its middle link the ball is momentarily at rest and in free fall. The chain starts there (the apex = the event
after half of the links), and its history is the rising half of that longest chain: a genuine piece of a geodesic of
the order, with no virtual events and no coordinates (the launch ring is chosen in coordinates, as before).
From the apex on, the rule of F5.1c is applied unchanged (p = m steps back in the history, or the launch event while
the history is shorter; Q = {q : L(v, q) = H}; q* = argmin L(p, q), ties min n(p, q); one step along a longest chain
to q*; rolling horizon). The box is taller (T = 46) so that the horizon fits up to the end of the readout.

Readout and prediction
----------------------
For each chain the displacement from the apex, projected on the outward radial unit vector at the apex, at the times
tau = t - t_apex on a grid 0..span (23 points). The mean over chains is fitted by v0 tau + a tau^2 / 2. The prediction is
the mean over chains of the fall from rest at r_apex in the geometry phi = -sigma/3 (RK4, as F5.1), fitted the same way:
a_pred per sprinkling; R = a / a_pred. The flat order uses the a_pred of the loaded geometry for the same apex radii
(the null: no fall). The control of F5.1 is kept as it was (longest chains between events at r0 = 8, t0 = 7 and t0 + 22;
R_max).

Criteria (stated before any full run)
------------------------------------
Classification by R over sprinklings:
  falls as predicted   1/2 < R - 2 SE  and  R + 2 SE < 2   (the band of the R_max instrument);
  falls too strongly   R - 2 SE >= 2;
  does not fall        R + 2 SE < 1/2;
  inconclusive         otherwise.
Instrument (INSTRUMENT_FAILED if any fails): flat |R| + 2 SE <= 1/2 for the walk and for the longest chains; loaded
longest chains 1/2 < R_max -+ 2 SE < 2.
PASSED if the walk falls as predicted; FAILED if it does not fall or falls too strongly; UNDECIDED otherwise.
Reported: R / R_max, v0 (the apex should be at rest: v0 near 0), the fraction of chains that reach the end, the
displacement series, and the mean squared coordinate speed <v^2> of the links of the walk.
The upper bound R < 2 was added on 30 Oct after the development rows below and before any full run. F5.1-F5.1c had
only the lower bound ("falls"), because none of their rules fell at all. The bound makes the test stricter, not
looser: the development rows lie near R = 2, and "as the geometry predicts" needs a magnitude, not only a direction.

Development rows (10 chains per order, not deciding):
  loaded seed 21: R = +1.67, v0 = +0.010, apex r 8.18; displacement from the apex (tau = 0, 2, ..., 22)
     0.00 -0.01 -0.15 -0.27 -0.45 -0.77 -1.26 -1.68 -2.24 -2.84 -3.49 -4.07;
  loaded seed 22: R = +2.23, v0 = +0.065; 0.00 +0.07 +0.22 +0.09 -0.21 -0.59 -0.99 -1.67 -2.28 -2.99 -3.70 -4.38;
  flat seed 21:   R = -0.21, v0 = -0.120; 0.00 +0.03 -0.22 -0.54 -0.89 -1.21 -1.19 -1.67 -1.78 -1.78 -1.95 -1.79.
  From the apex the loaded chains start at rest and accelerate inwards (the first rule of the series that does), at
  about twice the predicted rate (the geometry gives -2.7 over 22). A chain of links is not a slow body: in the quick
  run the links have <v^2> = 0.37 (loaded) to 0.57 (flat), only the mean velocity is zero. Whether the overshoot
  depends on the link speeds is open; <v^2> is reported for that question and is not a criterion. (For a body with a
  persistent velocity the weak-field geodesic averaged over isotropic directions in 2+1 gives -grad(phi) (1 - <v^2>),
  a reduction; a zigzag of links is not such a body, so no prediction is declared.)
Reading declared now. PASSED: started from a free state of the order, a history continued as a longest chain falls
with the Newtonian acceleration; F5.1c failed through its initial condition and box, not through the rule; free fall
follows from inertia plus the geodesic of the order. FAILED (does not fall): even from a free initial state the
step-by-step continuation of a longest chain does not carry the acceleration that the longest chains between fixed
events show; the rule loses the gravity of the order on the way (a property of the rule, m = H = 8). FAILED (falls too
strongly): the rule falls, but faster than a slow body in the geometry of the order; the magnitude is not Newtonian.
"""

from __future__ import annotations

import time

import numpy as np

from core.order.sprinkle import derive_seed
from core.stats.parallel import parallel_map
from core.stats.stats import boot_se
from core.testing import Result
from phases.F5_gravity.F5_01_free_fall_order import fall_curve, fit_a, longest_chain, ring_starts, sag_prediction
from phases.F5_gravity.F5_01c_free_fall_longest import LocalLinks, _popcount
from phases.F7_light.F7_01a_capacity_cones import CapacityOrder
from phases.F7_light.F7_01b_throughput_cones import sprinkle_saturated

TEST = {
    "id": "F5.1d",
    "phase": "F5",
    "title": "Free fall from the apex: does a longest-chain history that starts at the top of a thrown ball fall?",
    "theory_part": "F5 gravity and equivalence; free fall as inertia plus the geodesic of the order, from a free initial "
                   "state read from the order (the apex of a longest chain)",
    "question": "In the order shaped by a static source (F7.1b), does a chain that starts at the apex of the longest chain "
                "between two events at the same place, with the rising half of that chain as its history, and then follows "
                "the longest-chain rule of F5.1c, fall from the apex as the geometry of the order predicts?",
    "method": "2+1 order of F7.1b in a taller box (T = 46); launch ring r = 8 at t = 1.5; longest chain to the same place "
              "dT_launch = 14 later; start at its middle link with the rising half as history; rule of F5.1c (m = H = 8); "
              "radial projection from the apex over 22, fit a/2 tau^2; prediction: fall from rest at r_apex; control and "
              "flat null as F5.1",
    "prediction": "a free state continued as a geodesic of the order falls from the apex with a = -dphi/dr(r_apex)",
    "criterion": "falls as predicted 1/2 < R -+ 2 SE < 2; falls too strongly R - 2 SE >= 2; does not fall R + 2 SE < 1/2. "
                 "INSTRUMENT: flat |R| + 2 SE <= 1/2 (walk and longest chains); "
                 "loaded longest chains 1/2 < R_max -+ 2 SE < 2. PASSED falls as predicted, FAILED does not fall or falls too strongly, "
                 "UNDECIDED otherwise.",
    "null_model": "flat order (sigma0 = 0): no fall",
    "level": "mechanism",
    "meters": ["K.16b"],
    "depends_on": ["F5.1", "F5.1c", "F7.1b", "F2.5"],
    "params": {"T": 46.0, "half": 18.0, "width": 8.0, "sigma0": 0.45, "r_launch": 8.0, "t_launch": 1.5,
               "dT_launch": 14.0, "span": 22.0, "n_ring": 40, "memory": 8, "horizon": 8, "region": 12.0,
               "r0": 8.0, "t0": 7.0, "dT": 22.0, "n_max": 40, "n_seeds": 32, "n_se": 2.0, "boot_B": 4000},
    "parts": {"main": {}},
    "quick": {"n_seeds": 2, "n_ring": 6, "n_max": 6, "boot_B": 200},
    "estimate": "~2.5-3.5 h with 10 processes (64 orders of ~69 000 events; per step three longest-chain sweeps over a cylinder)",
}


# ===================================================================================== the walk

def walk_history(order, hist, t_end, m, H, region):
    """The longest-chain rule of F5.1c, started with a real history hist (a list of events ending at the start)."""
    X = order.X
    v = hist[-1]
    ll = LocalLinks(order, X[v, 1:], region, float(X[hist[0], 0]) - 1.0, float(X[:, 0].max()) + 1.0)
    pos = ll.pos
    path = list(hist)
    while X[v, 0] < t_end:
        jv = pos[v]
        src = path[len(path) - 1 - m] if len(path) > m else path[0]
        jp = pos[src]
        if jv < 0 or jp < 0:
            break
        Lv = ll.forward(jv)
        Q = np.flatnonzero(Lv == H)                       # 1. the horizon H links ahead
        if Q.size == 0:
            Hd = int(Lv.max())
            if Hd < 2:
                break
            Q = np.flatnonzero(Lv == Hd)
        Lp = ll.forward(jp)                               # 2. the target: no boost relative to the history
        Q = Q[Lp[Q] >= 0]
        if Q.size == 0:
            break
        sc = Lp[Q]
        Q = Q[sc == sc.min()]
        if Q.size > 1:
            fp = np.packbits(Lp >= 0)
            nq = np.array([ll.between(fp, jp, int(q)) for q in Q])
            Q = Q[nq == nq.min()]
        jq = int(Q.min())
        Lb = ll.backward(jq)                              # 3. one step along a longest chain to q*
        kids = ll.c_idx[ll.c_ptr[jv]:ll.c_ptr[jv + 1]]
        kids = kids[Lb[kids] == Lb[jv] - 1]
        if kids.size == 0:
            break
        if kids.size > 1:
            nk = np.array([_popcount(ll.C[jq] & np.packbits(ll.forward(int(k)) >= 0)) for k in kids])
            kids = kids[nk == nk.max()]
        v = int(ll.ids[int(kids.min())])
        path.append(v)
    return path


def classify_band(R, se, nse):
    """The classification of this test (module notes): a direction and a magnitude."""
    lo, hi = R - nse * se, R + nse * se
    if hi < 0.5:
        return "does not fall"
    if lo >= 2.0:
        return "falls too strongly"
    if lo > 0.5 and hi < 2.0:
        return "falls as predicted"
    return "inconclusive"


def launch_ring(X, p, n):
    """The events nearest to n equally spaced points of the launch ring r = r_launch at t = t_launch."""
    return ring_starts(X, {"r0": p["r_launch"], "t0": p["t_launch"]}, n)


def apex_history(order, s_, p):
    """The longest chain from s_ to the event nearest to the same place dT_launch later; its rising half (to the
    middle link) is the history, its middle event the apex. None if the chain does not exist."""
    X = order.X
    tq = X[s_, 0] + float(p["dT_launch"])
    target = np.array([tq, X[s_, 1], X[s_, 2]])
    fut = order.future([s_], tq + 1.0)
    cand = np.where(fut & (np.abs(X[:, 0] - tq) < 1.0))[0]
    if cand.size == 0:
        return None
    q = int(cand[np.argmin(((X[cand] - target) ** 2).sum(axis=1))])
    chain = longest_chain(order, s_, q)
    if chain is None or len(chain) < 3:
        return None
    return chain[: len(chain) // 2 + 1]


# ===================================================================================== one sprinkling

def seed_task(p, sigma0, seed):
    rng = np.random.default_rng(seed)
    T, Hh, w = float(p["T"]), float(p["half"]), float(p["width"])
    X = sprinkle_saturated(sigma0, w, T, Hh, rng)
    order = CapacityOrder(X, sigma0, w)
    span = float(p["span"])
    tau = np.linspace(0.0, span, 23)
    cols, preds, full, r_apex, v_apex, v2 = [], [], 0, [], [], []
    for s_ in launch_ring(X, p, int(p["n_ring"])):
        hist = apex_history(order, s_, p)
        if hist is None:
            continue
        ap = hist[-1]
        ta, ra = float(X[ap, 0]), X[ap, 1:]
        path = walk_history(order, hist, ta + span + 1.0, int(p["memory"]), int(p["horizon"]), float(p["region"]))
        P = X[path[len(hist) - 1:]]                         # from the apex on
        dP = np.diff(P, axis=0)
        if len(dP):
            v2.append(float(np.mean((dP[:, 1] ** 2 + dP[:, 2] ** 2) / dP[:, 0] ** 2)))   # squared speed of each link
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
        preds.append(fall_curve(dict(p, r0=float(np.linalg.norm(ra)), t0=0.0), tau))
    with np.errstate(all="ignore"):
        ybar = np.nanmean(cols, axis=0) if cols else np.full(len(tau), np.nan)
    ok = np.isfinite(ybar)
    a, v0 = fit_a(tau[ok], ybar[ok], 0.0) if ok.sum() >= 3 else (float("nan"), float("nan"))
    a_pred, _ = fit_a(tau, np.mean(preds, axis=0), 0.0) if preds else (float("nan"), 0.0)
    out = {"N": len(X), "sigma0": sigma0, "a": a, "v0": v0, "a_pred": a_pred, "n_chains": len(cols), "n_full": full,
           "series": ybar.tolist(), "r_apex": float(np.mean(r_apex)) if r_apex else float("nan"),
           "v_apex": float(np.mean(v_apex)) if v_apex else float("nan"), "v2": float(np.mean(v2)) if v2 else float("nan")}
    # control: longest chains between two events at the same place (as F5.1)
    dT = float(p["dT"])
    devs = []
    for s_ in ring_starts(X, p, int(p["n_max"])):
        tq = X[s_, 0] + dT
        target = np.array([tq, X[s_, 1], X[s_, 2]])
        fut = order.future([s_], tq + 1.0)
        cand = np.where(fut & (np.abs(X[:, 0] - tq) < 1.0))[0]
        if cand.size == 0:
            continue
        q = int(cand[np.argmin(((X[cand] - target) ** 2).sum(axis=1))])
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
    sag_pred = sag_prediction(p)
    jobs = []
    for s0 in (float(p["sigma0"]), 0.0):
        jobs += [(p, s0, derive_seed(ctx.seed, "F5.1d", ctx.part, s0, k)) for k in range(n)]
    rows = parallel_map(seed_task, jobs, ctx.workers, ctx.progress,
                        lambda i, r_, d, nn: f"[F5.1d] {d}/{nn}: sigma0 {r_['sigma0']:g} a {r_['a']:+.4f} (pred {r_['a_pred']:+.4f}) "
                                             f"sag {r_['sag']:+.3f} ({time.time() - t0:.0f}s)")
    res = {}
    for s0 in (float(p["sigma0"]), 0.0):
        rr = [r_ for r_ in rows if r_["sigma0"] == s0 and np.isfinite(r_["a"])]
        bs = lambda *tag: derive_seed(ctx.seed, "F5.1d-boot", ctx.part, s0, *tag)
        R = boot_se([r_["a"] / r_["a_pred"] for r_ in rr], B, bs("R"))[:2]
        Rm = boot_se([r_["sag"] / sag_pred for r_ in rr if np.isfinite(r_["sag"])], B, bs("Rm"))[:2]
        v0 = boot_se([r_["v0"] for r_ in rr], B, bs("v0"))[:2]
        with np.errstate(all="ignore"):
            ser = np.nanmean([r_["series"] for r_ in rr], axis=0).tolist()
        res[s0] = {"R": R, "R_max": Rm, "v0": v0, "series": ser, "N": float(np.mean([r_["N"] for r_ in rr])),
                   "a_pred": float(np.mean([r_["a_pred"] for r_ in rr])), "r_apex": float(np.mean([r_["r_apex"] for r_ in rr])),
                   "v_apex": float(np.mean([r_["v_apex"] for r_ in rr])), "v2": float(np.mean([r_["v2"] for r_ in rr])),
                   "chains": float(np.mean([r_["n_chains"] for r_ in rr])), "full": float(np.mean([r_["n_full"] for r_ in rr]))}
    s0 = float(p["sigma0"])
    lines = [f"2+1 order of F7.1b: sigma0 = {s0:g}, w = {p['width']:g}; box t [0, {p['T']:g}], |x|, |y| <= {p['half']:g} "
             f"(~{res[s0]['N']:.0f} events); launch ring r = {p['r_launch']:g} at t = {p['t_launch']:g}, {p['n_ring']} chains, "
             f"apex of the longest chain over dT_launch = {p['dT_launch']:g}; longest-chain rule m {p['memory']}, H {p['horizon']}; "
             f"fit over {p['span']:g} from the apex",
             f"  prediction (phi = -sigma/3, fall from rest at r_apex): a_pred = {res[s0]['a_pred']:+.5f}; control sag "
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
