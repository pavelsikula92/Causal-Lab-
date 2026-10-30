"""
F5.1b -- Free fall by straightness: does a flow that keeps its history straight in the order fall?
==================================================================================================

Why this test exists (29 Oct 2026, after F5.1 intrinsic FAILED; the author's picture)
------------------------------------------------------------------------------------
F5.1 (intrinsic part, FAILED; its verdict and criteria stay) showed that the order shaped by a static source
carries gravity (longest chains between fixed events sag as the geometry predicts, R_max = 0.83) but that a chain
which only remembers its ORIGIN and keeps to the shells of longest chains from it does not fall: it drifts away
from the source (R = -0.30). Free fall needs a memory of the VELOCITY: a discrete parallel transport of the
direction, built from the order alone.

The author's picture (29 Oct): a body held at rest in a field is accelerated; an accelerated observer has a horizon
behind it and sees the landscape distorted, as by a Doppler effect, blue ahead and red behind. A freely falling body
has no such distortion: around it the order looks locally balanced (the equivalence principle). In terms of the
order: a history that is locally straight has its surroundings balanced; a bent history has them lopsided.

The rule (input of this test; order and counts only, principle M3)
-----------------------------------------------------------------
A chain carries the last m events of its own history. Let p be the event m steps back and c the event m/2 steps back.
From the current event v it goes to the kept child k that makes p -> c -> k as straight as the order allows, i.e.
maximises the straightness
    s(k) = [ n(p, c)^(1/D) + n(c, k)^(1/D) ] / n(p, k)^(1/D),       D = 3 (2+1),
with n(x, y) the number of events strictly between x and y. For an interval the number of events is proportional to
its proper volume, i.e. to tau^D, so s is the ratio (tau(p, c) + tau(c, k)) / tau(p, k) of the reverse triangle
inequality: s <= 1, with equality when c lies on the geodesic from p to k. Maximising s continues the geodesic through
p and c: a discrete parallel transport of the direction, read by counting. It rejects the long, nearly null links
that a maximisation of an interval alone rewards (F5.1 development), because they bend the history. Ties: the earliest
child. The chain starts with a virtual history at rest: two events of the start's past at the same place, anchor_dt and
anchor_dt / 2 earlier, stand for p and c while the chain has fewer than m steps.
No dead ends (development, 29 Oct, before any full run): a candidate child must itself have a kept child inside the
cylinder (or lie past the end time). In the first development rows about half of the loaded chains stopped at an event
with no kept child (a sink of the loaded order); the flat order has almost none. A flow cannot continue from a sink, so
the one-step look-ahead is part of the rule, read from the order alone.
Development rows (16 chains, not deciding; with the look-ahead): loaded R = -0.01 / -0.12, v0 = -0.06 / -0.15, 15 / 14
of 16 chains reach the end; flat R = +0.26 / -0.30, v0 = +0.04 / -0.03, 16 / 16. The loaded chains move inwards by
-1.3 / -2.9 over 22 (the geometry predicts -2.7), but with a constant velocity after an early kick, not with an
acceleration. Memory m = 4 / 16 / 32 gives the same shape (R = -0.36 / -0.26 / -0.38). Expected outcome by the stated
criterion: FAILED. The criterion is not changed (R132).
Reading of v0 declared now (reported, not a criterion): an inward v0 with R near 0 means that straightness carries a
velocity (inertia) and that the order acts on the history only once, when the virtual history at rest is first
straightened, and not continuously as gravity.

The order, readout, predictions and controls are those of F5.1 (2+1 order of F7.1b, density exp(sigma), cones
exp(-sigma), sigma0 = 0.45, w = 8; a ring of 40 chains at rest at r0 = 8; the radial projection of the displacement
fitted by v0 t + a t^2 / 2 over 22; R = a / a_pred with a_pred from the fall in the geometry phi = -sigma/3; longest
chains between two events at the same place, R_max; the flat order as the null).

Criteria (stated before any full run; as F5.1)
---------------------------------------------
Classification by R over sprinklings: falls R - 2 SE > 1/2; does not fall R + 2 SE < 1/2; inconclusive otherwise.
Instrument (INSTRUMENT_FAILED if any fails): flat |R| + 2 SE <= 1/2 for the walk and for the longest chains; loaded
longest chains 1/2 < R_max -+ 2 SE < 2.
PASSED if the walk falls; FAILED if it does not fall; UNDECIDED otherwise. Reported: R / R_max (the walk against the
longest chains of the same orders), the fraction of chains that reach the end.
Reading declared now. PASSED: matter that keeps its own history straight in the order falls in the gravity of the
order; the equivalence principle follows from straightness read by counting, with no bending rule and no
coordinates. FAILED: straightness at this scale (m steps) does not transmit the gravity of the order to a history.
"""

from __future__ import annotations

import math
import time

import numpy as np

from core.order.sprinkle import derive_seed
from core.stats.parallel import parallel_map
from core.stats.stats import boot_se
from core.testing import Result
from phases.F5_gravity.F5_01_free_fall_order import (classify, fall_curve, fit_a, longest_chain, radial_series,
                                                        ring_starts, sag_prediction)
from phases.F7_light.F7_01a_capacity_cones import CapacityOrder
from phases.F7_light.F7_01b_throughput_cones import sprinkle_saturated

TEST = {
    "id": "F5.1b",
    "phase": "F5",
    "title": "Free fall by straightness: does a flow that keeps its history straight in the order fall?",
    "theory_part": "F5 gravity and equivalence; the equivalence principle as straightness of a history read by counting (M3)",
    "question": "In the order shaped by a static source (F7.1b), does a chain that goes to the child maximising the straightness "
                "(tau(p, c) + tau(c, k)) / tau(p, k) of its own recent history (tau from counts) fall towards the source as the "
                "geometry of the order predicts?",
    "method": "as F5.1 (2+1 order of F7.1b, ring of 40 chains at rest, radial projection, fit a/2 t^2, longest chains as the "
              "control, flat null); the rule: carried history of m steps, p m steps back, c m/2 steps back, child maximising "
              "s = [n(p,c)^(1/3) + n(c,k)^(1/3)] / n(p,k)^(1/3)",
    "prediction": "a history that is straight in the order is a geodesic of the order: it falls with a = -dphi/dr",
    "criterion": "falls R - 2 SE > 1/2, does not fall R + 2 SE < 1/2. INSTRUMENT: flat |R| + 2 SE <= 1/2 (walk and longest chains); "
                 "loaded longest chains 1/2 < R_max -+ 2 SE < 2. PASSED falls, FAILED does not fall, UNDECIDED otherwise.",
    "null_model": "flat order (sigma0 = 0): no fall",
    "level": "mechanism",
    "meters": ["K.16b"],
    "depends_on": ["F5.1", "F7.1b", "F2.5"],
    "params": {"T": 32.0, "half": 18.0, "width": 8.0, "sigma0": 0.45, "r0": 8.0, "t0": 7.0, "span": 22.0,
               "n_ring": 40, "memory": 8, "anchor_dt": 6.0, "region": 12.0, "dT": 22.0, "n_max": 40, "n_seeds": 32,
               "n_se": 2.0, "boot_B": 4000},
    "parts": {"main": {}},
    "quick": {"n_seeds": 2, "n_ring": 12, "n_max": 6, "boot_B": 200},
    "estimate": "~30-60 min with 10 processes (64 orders of ~48 000 events; per chain a bitset closure of a cylinder around it)",
}

D = 3


class LocalOrder:
    """Past closures (bitsets) of the kept-link order restricted to a cylinder around a chain: n(x, y) by counting."""

    def __init__(self, order, centre, radius, t_lo, t_hi):
        X = order.X
        m = (((X[:, 1:] - centre) ** 2).sum(axis=1) <= radius * radius) & (X[:, 0] >= t_lo) & (X[:, 0] <= t_hi)
        ids = np.where(m)[0]
        ids = ids[np.argsort(X[ids, 0], kind="stable")]
        self.pos = np.full(order.N, -1, dtype=np.int64)
        self.pos[ids] = np.arange(ids.size)
        self.ids = ids
        nb = (ids.size + 7) // 8
        C = np.zeros((ids.size, nb), dtype=np.uint8)
        for j, w in enumerate(ids):
            pj = self.pos[order.pidx[order.pptr[w]:order.pptr[w + 1]]]
            pj = pj[pj >= 0]
            if pj.size:
                C[j] = np.bitwise_or.reduce(C[pj], axis=0)
            C[j, j >> 3] |= np.uint8(1 << (7 - (j & 7)))   # closure includes the event itself
        self.C = C

    def between(self, x, y):
        """Number of events strictly between x and y (0 if not related or outside the cylinder)."""
        jx, jy = self.pos[x], self.pos[y]
        if jx < 0 or jy < 0:
            return 0
        past_y = np.flatnonzero(np.unpackbits(self.C[jy])[: self.ids.size])   # events in the past of y (y included)
        if past_y.size == 0:
            return 0
        bit = (self.C[past_y, jx >> 3] >> (7 - (jx & 7))) & 1                 # x in the past of w (x included)
        cnt = int(bit.sum())
        return max(cnt - 2, 0) if bit[np.searchsorted(past_y, jy)] else 0      # minus x and y themselves


def walk_straight(order, v0, hist0, t_end, m, region):
    """The straightness rule (module notes). hist0: the virtual history [p0, c0] at rest before the start."""
    X, cptr, cidx = order.X, order.cptr, order.cidx
    lo = LocalOrder(order, X[v0, 1:], region, float(X[hist0[0], 0]) - 1.0, t_end + 2.0)
    path, v = [v0], v0
    while X[v, 0] < t_end:
        a, b = cptr[v], cptr[v + 1]
        if b == a:
            break
        ks = cidx[a:b]
        ks = ks[lo.pos[ks] >= 0]
        if ks.size == 0:
            break
        L = len(path)
        if L > m:
            p_, c_ = path[L - 1 - m], path[L - 1 - m // 2]
        elif L > m // 2:
            p_, c_ = hist0[0], path[L - 1 - m // 2]
        else:
            p_, c_ = hist0[0], hist0[1]
        npc = lo.between(p_, c_) ** (1.0 / D)
        best, arg = -1.0, []
        for k in ks:
            # A flow cannot enter a dead end: the child must itself have a kept child inside the
            # cylinder, unless it already lies past the end time (one-step look-ahead, order only).
            if X[k, 0] < t_end:
                kc = cidx[cptr[k]:cptr[k + 1]]
                if not np.any(lo.pos[kc] >= 0):
                    continue
            npk = lo.between(p_, int(k))
            if npk <= 0:
                continue
            s = (npc + lo.between(c_, int(k)) ** (1.0 / D)) / npk ** (1.0 / D)
            if s > best + 1e-12:
                best, arg = s, [int(k)]
            elif abs(s - best) <= 1e-12:
                arg.append(int(k))
        if not arg:
            break
        v = min(arg, key=lambda kk: X[kk, 0])
        path.append(v)
    return path


def seed_task(p, sigma0, seed):
    rng = np.random.default_rng(seed)
    T, H, w = float(p["T"]), float(p["half"]), float(p["width"])
    X = sprinkle_saturated(sigma0, w, T, H, rng)
    order = CapacityOrder(X, sigma0, w)
    t0, span = float(p["t0"]), float(p["span"])
    grid = np.linspace(t0, t0 + span, 23)
    starts = ring_starts(X, p, int(p["n_ring"]))
    paths = []
    adt = float(p["anchor_dt"])
    for s_ in starts:
        hist = []
        for dt_ in (adt, adt / 2.0):
            ta = X[s_, 0] - dt_
            past = order._bfs([s_], order.pptr, order.pidx, X[:, 0] >= ta - 1.0)
            cand = np.where(past & (np.abs(X[:, 0] - ta) < 0.75))[0]
            hist.append(int(cand[np.argmin(((X[cand, 1:] - X[s_, 1:]) ** 2).sum(axis=1))]) if cand.size else s_)
        paths.append(walk_straight(order, s_, hist, t0 + span + 1.0, int(p["memory"]), float(p["region"])))
    y, n_ok = radial_series(X, paths, p, grid)
    ok = np.isfinite(y)
    a, v0 = fit_a(grid[ok], y[ok], t0)
    out = {"N": len(X), "sigma0": sigma0, "a": a, "v0": v0, "n_chains": n_ok, "series": y.tolist()}
    dT = float(p["dT"])
    devs = []
    for s_ in starts[: int(p["n_max"])]:
        tq = X[s_, 0] + dT
        target = np.array([tq, X[s_, 1], X[s_, 2]])
        fut = order.future([s_], tq + 1.0)
        cand = np.where(fut & (np.abs(X[:, 0] - tq) < 1.0))[0]
        if cand.size == 0:
            continue
        q = int(cand[np.argmin(((X[cand] - target) ** 2).sum(axis=1))])
        path = longest_chain(order, s_, q)
        if path is None:
            continue
        P = X[path]
        rhat = P[0, 1:] / np.linalg.norm(P[0, 1:])
        proj = (P[:, 1:] - P[0, 1:]) @ rhat
        tm = 0.5 * (P[0, 0] + P[-1, 0])
        rline = np.interp(tm, [P[0, 0], P[-1, 0]], [proj[0], proj[-1]])
        devs.append(float(np.interp(tm, P[:, 0], proj) - rline))
    out["sag"] = float(np.mean(devs)) if devs else float("nan")
    return out


def run(p, ctx):
    t0 = time.time()
    n, B, nse = int(p["n_seeds"]), int(p["boot_B"]), float(p["n_se"])
    grid = np.linspace(float(p["t0"]), float(p["t0"]) + float(p["span"]), 23)
    a_pred, _ = fit_a(grid, fall_curve(p, grid), float(p["t0"]))
    sag_pred = sag_prediction(p)
    jobs = []
    for s0 in (float(p["sigma0"]), 0.0):
        jobs += [(p, s0, derive_seed(ctx.seed, "F5.1b", ctx.part, s0, k)) for k in range(n)]
    rows = parallel_map(seed_task, jobs, ctx.workers, ctx.progress,
                        lambda i, r_, d, nn: f"[F5.1b] {d}/{nn}: sigma0 {r_['sigma0']:g} a {r_['a']:+.4f} sag {r_['sag']:+.3f} "
                                             f"({time.time() - t0:.0f}s)")
    res = {}
    for s0 in (float(p["sigma0"]), 0.0):
        rr = [r_ for r_ in rows if r_["sigma0"] == s0]
        bs = lambda *tag: derive_seed(ctx.seed, "F5.1b-boot", ctx.part, s0, *tag)
        R = boot_se([r_["a"] / a_pred for r_ in rr], B, bs("R"))[:2]
        Rm = boot_se([r_["sag"] / sag_pred for r_ in rr if np.isfinite(r_["sag"])], B, bs("Rm"))[:2]
        v0 = boot_se([r_["v0"] for r_ in rr], B, bs("v0"))[:2]
        with np.errstate(all="ignore"):
            ser = np.nanmean([r_["series"] for r_ in rr], axis=0).tolist()
        res[s0] = {"R": R, "R_max": Rm, "v0": v0, "series": ser, "N": float(np.mean([r_["N"] for r_ in rr])),
                   "chains": float(np.mean([r_["n_chains"] for r_ in rr]))}
    s0 = float(p["sigma0"])
    lines = [f"2+1 order of F7.1b: sigma0 = {s0:g}, w = {p['width']:g}; box t [0, {p['T']:g}], |x|, |y| <= {p['half']:g} "
             f"(~{res[s0]['N']:.0f} events); ring r0 = {p['r0']:g} at t0 = {p['t0']:g}, {p['n_ring']} chains, fit over "
             f"{p['span']:g}; straightness rule, memory {p['memory']} steps",
             f"  prediction (phi = -sigma/3): a_pred = {a_pred:+.5f}; longest-chain sag {sag_pred:+.4f}"]
    for key, lab in ((s0, "loaded"), (0.0, "flat  ")):
        r_ = res[key]
        lines.append(f"  {lab}: R = {r_['R'][0]:+.3f} +- {r_['R'][1]:.3f} ({classify(*r_['R'], nse)}); v0 {r_['v0'][0]:+.4f} +- "
                     f"{r_['v0'][1]:.4f}; longest chains R_max = {r_['R_max'][0]:+.3f} +- {r_['R_max'][1]:.3f}; chains to the "
                     f"end {r_['chains']:.1f} of {p['n_ring']}")
        lines.append("          displacement: " + " ".join(f"{x:+.2f}" for x in r_["series"][::2]))
    Rl, Rml = res[s0]["R"][0], res[s0]["R_max"][0]
    lines.append(f"  reported: R / R_max = {Rl / Rml:+.3f}" if Rml else "  reported: R_max = 0")
    f_ = res[0.0]
    inst = {"flat: no drift (|R| + 2 SE <= 1/2)": abs(f_["R"][0]) + nse * f_["R"][1] <= 0.5,
            "flat longest chains: no sag (|R_max| + 2 SE <= 1/2)": abs(f_["R_max"][0]) + nse * f_["R_max"][1] <= 0.5,
            "loaded longest chains read the geometry (1/2 < R_max -+ 2 SE < 2)":
                res[s0]["R_max"][0] - nse * res[s0]["R_max"][1] > 0.5 and res[s0]["R_max"][0] + nse * res[s0]["R_max"][1] < 2.0}
    lines.append("instrument: " + "; ".join(f"{k}: {'OK' if ok else 'FAILS'}" for k, ok in inst.items()))
    cls = classify(*res[s0]["R"], nse)
    status = "INSTRUMENT_FAILED" if not all(inst.values()) else {"falls": "PASSED", "does not fall": "FAILED"}.get(cls, "UNDECIDED")
    data = {"a_pred": a_pred, "sag_pred": sag_pred, "res": {str(k): v for k, v in res.items()}, "rows": rows}
    return Result(status, lines, data)
