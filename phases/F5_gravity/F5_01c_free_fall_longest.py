"""
F5.1c -- Free fall by longest chains: does a flow that continues its history as a longest chain fall?
====================================================================================================

Why this test exists (29 Oct 2026, after the F5.1b development rows; the author's picture of a flowing substrate)
---------------------------------------------------------------------------------------------------------------
F5.1 showed that the order shaped by a static source carries gravity: the longest chains between two events at the
same place bulge as the geometry predicts (R_max = 0.83). Two step-by-step rules failed to make a history fall:
F5.1 intrinsic (memory of the origin only; FAILED, drifts outwards) and, in development, F5.1b (straightness of the
recent history read from interval VOLUMES; it moves inwards once, early, and then keeps a constant velocity). The
reverse triangle inequality read from volumes over a few links is a second-order effect buried in the counting noise.

The author's picture (29 Oct): the surroundings are dragged, as if the light cones were tilted towards the source
(the "river" picture of a static field). Worked out on paper before this test (recorded in docs/findings.md): being
carried by the river alone does not give Newton (a body released from rest has (x' . grad) v = 0 at the start); the
river picture gives x'' = grad(v^2 / 2), which is the geodesic equation. And the frame in which the surroundings look
balanced is, by the equivalence principle, every freely falling frame whatever its velocity: the surroundings fix
the frame only together with a memory of the velocity. Both routes lead to: a history that continues itself as a
geodesic of the order. The longest chain is the geodesic of an order, and F5.1 shows it is read correctly.

The rule (input of this test; order and chain lengths only, principle M3)
-----------------------------------------------------------------------
L(x, y) is the number of links of the longest chain from x to y (kept links of the order; -1 if y is not in the
future of x). A chain carries its own history; p is the event m steps back (the virtual history at rest, an event of
the start's past at the same place anchor_dt earlier, stands for p while the chain has fewer than m steps).
At the current event v:
  1. the horizon: Q = { q : L(v, q) = H } (events H links ahead; if none, the largest depth H' >= 2 available);
  2. the target: q* in Q with the SHORTEST longest chain from p, min L(p, q); ties by the fewest events n(p, q)
     between p and q, then the earliest. By the reverse triangle inequality L(p, q) >= L(p, v) + L(v, q), with
     equality when v lies on a longest chain from p to q; on the hyperboloid tau(v, q) = H the proper time
     tau(p, q)^2 = tau(p, v)^2 + H^2 + 2 tau(p, v) H cosh(eta) is smallest for eta = 0, i.e. for the q that continues the
     history p -> v without a boost. Fixing L(v, q) = H excludes the long, nearly null links that a plain maximisation
     rewards (F5.1 development).
     Development (29 Oct, before any run): the first version took max L(p, q), the wrong side of the inequality; it
     picks the most boosted q and the chains ran outwards at v0 = +0.2 to +0.5 (loaded) -- recorded, then corrected.
  3. the step: to the kept child k of v that lies on a longest chain from v to q* (L(k, q*) = L(v, q*) - 1); ties by
     the most events n(k, q*), then the earliest child.
The target is chosen again after every step (a rolling horizon); nothing is kept but the chain's own history.
Coordinates enter only the readout and the choice of the start ring.

Development rows (8-16 chains per order, not deciding; after the correction of step 2):
  loaded seeds 22 / 21 / 23: R = +1.20 / +0.27 / -0.12, v0 = +0.04 / -0.19 / -0.28, displacement over 22 -2.2 / -4.4 /
  -5.3 (the geometry predicts -2.7); flat seeds 21 / 22 / 23: R = +0.21 / +0.12 / -0.25, displacement -1.7 / -0.6 / -0.5.
  The loaded chains move inwards, in two of three orders with an early inward kick as in F5.1b; the kick does not change
  with anchor_dt = 3, H = 4 or m = H = 16 (seed 23: R = -0.31 / +0.40 / -0.22). Parameters stay m = H = 8 as first set.
  The control R_max is that of F5.1; F5.1 carried failed it narrowly (0.630 +- 0.074); the band is not changed.

The order, readout, predictions and controls are those of F5.1 (2+1 order of F7.1b, density exp(sigma), cones
exp(-sigma), sigma0 = 0.45, w = 8; a ring of 40 chains at rest at r0 = 8; the radial projection of the displacement
fitted by v0 t + a t^2 / 2 over 22; R = a / a_pred with a_pred from the fall in the geometry phi = -sigma/3; longest
chains between two events at the same place, R_max; the flat order as the null).

Criteria (stated before any full run; as F5.1 and F5.1b)
-------------------------------------------------------
Classification by R over sprinklings: falls R - 2 SE > 1/2; does not fall R + 2 SE < 1/2; inconclusive otherwise.
Instrument (INSTRUMENT_FAILED if any fails): flat |R| + 2 SE <= 1/2 for the walk and for the longest chains; loaded
longest chains 1/2 < R_max -+ 2 SE < 2.
PASSED if the walk falls; FAILED if it does not fall; UNDECIDED otherwise. Reported: R / R_max, v0, the fraction of
chains that reach the end.
Reading declared now. PASSED: a history that continues itself as a longest chain of the order falls with the
Newtonian acceleration; free fall (and with it the equivalence principle) follows from inertia (the carried history)
plus the geodesic of the order, both read by counting. FAILED: the step-by-step continuation of a longest chain with
memory m and horizon H does not transmit the gravity that the longest chains between fixed events carry.
"""

from __future__ import annotations

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
    "id": "F5.1c",
    "phase": "F5",
    "title": "Free fall by longest chains: does a flow that continues its history as a longest chain fall?",
    "theory_part": "F5 gravity and equivalence; free fall as inertia (a carried history) plus the geodesic of the order "
                   "(the longest chain), read by counting (M3)",
    "question": "In the order shaped by a static source (F7.1b), does a chain that at every step aims at the event H links "
                "ahead with the shortest longest chain from its own past (m steps back; no boost), and steps along a longest chain to it, fall "
                "towards the source as the geometry of the order predicts?",
    "method": "as F5.1 (2+1 order of F7.1b, ring of 40 chains at rest, radial projection, fit a/2 t^2, longest chains as the "
              "control, flat null); the rule: p = m steps back, Q = {q : L(v,q) = H}, q* = argmin L(p,q) (ties min n(p,q)), "
              "step to the child on a longest chain to q*; rolling horizon",
    "prediction": "a history continued as a longest chain is a geodesic of the order: it falls with a = -dphi/dr",
    "criterion": "falls R - 2 SE > 1/2, does not fall R + 2 SE < 1/2. INSTRUMENT: flat |R| + 2 SE <= 1/2 (walk and longest chains); "
                 "loaded longest chains 1/2 < R_max -+ 2 SE < 2. PASSED falls, FAILED does not fall, UNDECIDED otherwise.",
    "null_model": "flat order (sigma0 = 0): no fall",
    "level": "mechanism",
    "meters": ["K.16b"],
    "depends_on": ["F5.1", "F5.1b", "F7.1b", "F2.5"],
    "params": {"T": 32.0, "half": 18.0, "width": 8.0, "sigma0": 0.45, "r0": 8.0, "t0": 7.0, "span": 22.0,
               "n_ring": 40, "memory": 8, "horizon": 8, "anchor_dt": 6.0, "region": 12.0, "dT": 22.0, "n_max": 40,
               "n_seeds": 32, "n_se": 2.0, "boot_B": 4000},
    "parts": {"main": {}},
    "quick": {"n_seeds": 2, "n_ring": 12, "n_max": 6, "boot_B": 200},
    "estimate": "~1.5-2 h with 10 processes (64 orders of ~48 000 events; per step three longest-chain sweeps over a cylinder)",
}


# ================================================================================ the order around one chain

def _popcount(a):
    """Number of set bits of a packed uint8 array."""
    return int(np.unpackbits(a).sum())


class LocalLinks:
    """The kept links of the order restricted to a cylinder around a chain.

    Longest chains are computed by a vectorised relaxation over the links (Bellman-Ford with max instead of min; the
    links point forward in time, so it converges after as many sweeps as the longest chain has links). Counts of events
    between two events use past closures stored as bitsets.
    """

    def __init__(self, order, centre, radius, t_lo, t_hi):
        X = order.X
        m = (((X[:, 1:] - centre) ** 2).sum(axis=1) <= radius * radius) & (X[:, 0] >= t_lo) & (X[:, 0] <= t_hi)
        ids = np.where(m)[0]
        ids = ids[np.argsort(X[ids, 0], kind="stable")]
        self.ids, self.n = ids, ids.size
        self.pos = np.full(order.N, -1, dtype=np.int64)
        self.pos[ids] = np.arange(ids.size)
        # links inside the cylinder, in local indices (eu earlier, ew later)
        cnt = order.cptr[ids + 1] - order.cptr[ids]
        eu = np.repeat(np.arange(ids.size), cnt)
        ew = self.pos[np.concatenate([order.cidx[order.cptr[w]:order.cptr[w + 1]] for w in ids])] if ids.size else eu
        ok = ew >= 0
        eu, ew = eu[ok], ew[ok]
        # sorted by the later end (forward sweeps) and by the earlier end (backward sweeps)
        o = np.argsort(ew, kind="stable")
        self.f_u, self.f_w = eu[o], ew[o]
        o = np.argsort(eu, kind="stable")
        self.b_u, self.b_w = eu[o], ew[o]
        # children of each local event (for the step), CSR in local indices
        self.c_ptr = np.concatenate([[0], np.cumsum(np.bincount(self.b_u, minlength=self.n))]).astype(np.int64)
        self.c_idx = self.b_w
        # past closures (event included) as bitsets, for counting
        nb = (self.n + 7) // 8
        C = np.zeros((self.n, nb), dtype=np.uint8)
        p_ptr = np.concatenate([[0], np.cumsum(np.bincount(self.f_w, minlength=self.n))]).astype(np.int64)
        for j in range(self.n):
            pj = self.f_u[p_ptr[j]:p_ptr[j + 1]]
            if pj.size:
                C[j] = np.bitwise_or.reduce(C[pj], axis=0)
            C[j, j >> 3] |= np.uint8(1 << (7 - (j & 7)))
        self.C = C

    @staticmethod
    def _sweep(n, src, srt_from, srt_to, lo_ok):
        """Longest chain lengths from src along links (from -> to), keeping only links with lo_ok."""
        a, b = srt_from[lo_ok], srt_to[lo_ok]
        L = np.full(n, -1, dtype=np.int64)
        L[src] = 0
        if b.size == 0:
            return L
        starts = np.flatnonzero(np.r_[True, b[1:] != b[:-1]])
        tg = b[starts]
        while True:
            la = L[a]
            vals = np.where(la >= 0, la + 1, -1)
            red = np.maximum.reduceat(vals, starts)
            new = np.maximum(L[tg], red)
            if np.array_equal(new, L[tg]):
                return L
            L[tg] = new

    def forward(self, js):
        """L(x, y) for all y in the cylinder (local index js of x); -1 outside its future."""
        return self._sweep(self.n, js, self.f_u, self.f_w, self.f_u >= js)

    def backward(self, jq):
        """L(y, q) for all y in the cylinder (local index jq of q); -1 outside its past."""
        return self._sweep(self.n, jq, self.b_w, self.b_u, self.b_w <= jq)

    def between(self, fut_bits, jx, jy):
        """Events strictly between x and y, given the packed future closure of x (x included)."""
        return max(_popcount(self.C[jy] & fut_bits) - 2, 0)


def walk_longest(order, v0, p0, t_end, m, H, region):
    """The longest-chain rule (module notes). p0: the virtual history at rest (one event of the start's past)."""
    X = order.X
    ll = LocalLinks(order, X[v0, 1:], region, float(X[p0, 0]) - 1.0, float(X[:, 0].max()) + 1.0)
    pos = ll.pos
    path = [v0]
    v = v0
    while X[v, 0] < t_end:
        jv = pos[v]
        Lp_src = path[len(path) - 1 - m] if len(path) > m else p0
        jp = pos[Lp_src]
        if jv < 0 or jp < 0:
            break
        Lv = ll.forward(jv)
        # 1. the horizon H links ahead (or the largest depth >= 2 still inside the box)
        Q = np.flatnonzero(Lv == H)
        if Q.size == 0:
            Hd = int(Lv.max())
            if Hd < 2:
                break
            Q = np.flatnonzero(Lv == Hd)
        # 2. the target: the shortest longest chain from p (no boost), then the fewest events between p and q
        Lp = ll.forward(jp)
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
        # 3. the step: a child of v on a longest chain to q*
        Lb = ll.backward(jq)
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


# ===================================================================================== one sprinkling

def seed_task(p, sigma0, seed):
    rng = np.random.default_rng(seed)
    T, Hh, w = float(p["T"]), float(p["half"]), float(p["width"])
    X = sprinkle_saturated(sigma0, w, T, Hh, rng)
    order = CapacityOrder(X, sigma0, w)
    t0, span = float(p["t0"]), float(p["span"])
    grid = np.linspace(t0, t0 + span, 23)
    starts = ring_starts(X, p, int(p["n_ring"]))
    paths = []
    adt = float(p["anchor_dt"])
    for s_ in starts:
        # the virtual history at rest: the event of the start's past nearest to the same place, anchor_dt earlier
        ta = X[s_, 0] - adt
        past = order._bfs([s_], order.pptr, order.pidx, X[:, 0] >= ta - 1.0)
        cand = np.where(past & (np.abs(X[:, 0] - ta) < 0.75))[0]
        p0 = int(cand[np.argmin(((X[cand, 1:] - X[s_, 1:]) ** 2).sum(axis=1))]) if cand.size else s_
        paths.append(walk_longest(order, s_, p0, t0 + span + 1.0, int(p["memory"]), int(p["horizon"]),
                                  float(p["region"])))
    y, n_ok = radial_series(X, paths, p, grid)
    ok = np.isfinite(y)
    a, v0 = fit_a(grid[ok], y[ok], t0)
    out = {"N": len(X), "sigma0": sigma0, "a": a, "v0": v0, "n_chains": n_ok, "series": y.tolist(),
           "steps": float(np.mean([len(q) for q in paths]))}
    # control: longest chains between two events at the same place, dT apart (as F5.1)
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


# ===================================================================================== the test

def run(p, ctx):
    t0 = time.time()
    n, B, nse = int(p["n_seeds"]), int(p["boot_B"]), float(p["n_se"])
    grid = np.linspace(float(p["t0"]), float(p["t0"]) + float(p["span"]), 23)
    a_pred, _ = fit_a(grid, fall_curve(p, grid), float(p["t0"]))
    sag_pred = sag_prediction(p)
    jobs = []
    for s0 in (float(p["sigma0"]), 0.0):
        jobs += [(p, s0, derive_seed(ctx.seed, "F5.1c", ctx.part, s0, k)) for k in range(n)]
    rows = parallel_map(seed_task, jobs, ctx.workers, ctx.progress,
                        lambda i, r_, d, nn: f"[F5.1c] {d}/{nn}: sigma0 {r_['sigma0']:g} a {r_['a']:+.4f} sag {r_['sag']:+.3f} "
                                             f"({time.time() - t0:.0f}s)")
    res = {}
    for s0 in (float(p["sigma0"]), 0.0):
        rr = [r_ for r_ in rows if r_["sigma0"] == s0]
        bs = lambda *tag: derive_seed(ctx.seed, "F5.1c-boot", ctx.part, s0, *tag)
        R = boot_se([r_["a"] / a_pred for r_ in rr], B, bs("R"))[:2]
        Rm = boot_se([r_["sag"] / sag_pred for r_ in rr if np.isfinite(r_["sag"])], B, bs("Rm"))[:2]
        v0 = boot_se([r_["v0"] for r_ in rr], B, bs("v0"))[:2]
        with np.errstate(all="ignore"):
            ser = np.nanmean([r_["series"] for r_ in rr], axis=0).tolist()
        res[s0] = {"R": R, "R_max": Rm, "v0": v0, "series": ser, "N": float(np.mean([r_["N"] for r_ in rr])),
                   "chains": float(np.mean([r_["n_chains"] for r_ in rr])),
                   "steps": float(np.mean([r_["steps"] for r_ in rr]))}
    s0 = float(p["sigma0"])
    lines = [f"2+1 order of F7.1b: sigma0 = {s0:g}, w = {p['width']:g}; box t [0, {p['T']:g}], |x|, |y| <= {p['half']:g} "
             f"(~{res[s0]['N']:.0f} events); ring r0 = {p['r0']:g} at t0 = {p['t0']:g}, {p['n_ring']} chains, fit over "
             f"{p['span']:g}; longest-chain rule, memory {p['memory']}, horizon {p['horizon']} links",
             f"  prediction (phi = -sigma/3): a_pred = {a_pred:+.5f}; longest-chain sag {sag_pred:+.4f}"]
    for key, lab in ((s0, "loaded"), (0.0, "flat  ")):
        r_ = res[key]
        lines.append(f"  {lab}: R = {r_['R'][0]:+.3f} +- {r_['R'][1]:.3f} ({classify(*r_['R'], nse)}); v0 {r_['v0'][0]:+.4f} +- "
                     f"{r_['v0'][1]:.4f}; longest chains R_max = {r_['R_max'][0]:+.3f} +- {r_['R_max'][1]:.3f}; chains to the "
                     f"end {r_['chains']:.1f} of {p['n_ring']}; {r_['steps']:.1f} links per chain")
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
