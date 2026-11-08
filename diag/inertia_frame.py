"""
Diagnostic: is the computation cylinder the "medium" that brakes the 3+1 history? (development; decides nothing)
===============================================================================================================

Why (8 Nov 2026, after diag/fall_layers in 2+1 and 3+1)
-------------------------------------------------------
The three-layer diagnostic showed that in 3+1 the longest chains read the geometry (tube sag loaded - flat 0.47 +- 0.12
against 0.61) and that the tie-breakers do not decide the fall, but the history continued by the rule of F5.1c/d/e
behaves like a body in a medium at rest in the laboratory: the geometry turns it inwards at the apex (jump -0.17), then
its velocity relaxes within ~4 units of time, to rest in the flat order and to a terminal inward velocity in the loaded
one. Flat throws lose most of their velocity (0.18 -> 0.06).

The candidate for the frame at rest: the rule reads the order only inside a cylinder of LocalLinks that is fixed at the
start (radius 4 in 3+1) and has its axis along the laboratory time. The horizon H = 8 links reaches ~7.4 units of time,
so its candidates (a hyperboloid of events with L(v, q) = H) leave the cylinder; the cut removes more of the candidates
that continue the velocity and pulls the target back towards the axis.

The test: the flat 3+1 order (density 1, cones c = 1, as the flat arm of F5.1e), histories thrown with a velocity v
along +-x, and the registered rule unchanged (memory m = 8, horizon H = 8, the same tie-breakers); only the region in
which it reads the order changes: a cylinder of radius R whose axis moves with a velocity u along the throw,
x_axis(t) = x_start + u (t - t_start). Arms (v, u, R):

    (0, 0, 4)     at rest in a cylinder at rest                 (the domain of F5.1e)
    (0.3, 0, 4)   thrown in a cylinder at rest                  (the domain of F5.1e)
    (0.3, 0.3, 4) thrown in a cylinder that moves with it
    (0, 0.3, 4)   at rest in a cylinder that moves away
    and the same four at R = 6 (wider than the reach of most of the horizon set, narrower than all of it)

Every arm of one start uses the same start event and, for the same v, the same history. The top of the box (T = 28)
leaves >= 10 units of time above the end of every walk (F5.1e: ~7), so only the cylinder cuts the horizon.

Reading declared before the run
-------------------------------
  * The cylinder is the medium: the late velocity of a history goes towards u. A thrown history keeps its velocity in
    the co-moving cylinder and loses it in the cylinder at rest (paired difference co - rest > 0); a history at rest is
    pulled along by the moving cylinder (pull - rest > 0); both effects are smaller at R = 6 than at R = 4.
  * The laboratory itself is the medium (the box, a coordinate-time rule somewhere): co and rest arms the same, the
    pulled history stays at rest.
  * The time top of F5.1e was the cut: here the thrown history keeps its velocity already in the cylinder at rest.
  * Inertia is lost in the rule itself: all thrown histories slow down, whatever u and R.
The diagnostic decides nothing; a rule with a changed region would be a new registered test (F5.1f) with its own
criteria (R132). Also recorded: how often the horizon set, the tie set of min L(p, q) and the chosen target lie within
0.5 of the wall of the cylinder.

Run from the root of the laboratory:
    python -m diag.inertia_frame --seeds 8 --starts 16 --workers 8
    python -m diag.inertia_frame --smoke            (a small box only to check that the code runs)
Output: a summary on the screen and the raw records in diag/out/inertia_frame_{stamp}.json; every sprinkling is also
saved in diag/out/rows/inertia_*.json as soon as it is done (--summarize reads them after a stopped run).
"""

from __future__ import annotations

import argparse
import glob
import json
import math
import os
import time

import numpy as np

from core.order.sprinkle import derive_seed
from core.stats.parallel import parallel_map
from phases.F5_gravity.F5_01c_free_fall_longest import LocalLinks, _popcount
from phases.F5_gravity.F5_01d_free_fall_apex import walk_history
from phases.F7_light.F7_01c_throughput_cones_3d import CapacityOrder3
from diag.fall_layers import chain_between

# the flat order and the rule of F5.1e; a longer box along x (the direction of the throws) and a taller one in time
PARAMS = {"T": 28.0, "half_x": 11.0, "half_yz": 7.5, "width": 4.0, "memory": 8, "horizon": 8, "span": 12.0,
          "t_launch": 1.0, "dT_launch": 6.0, "x_start": 2.5, "jitter": 1.0, "region_hist": 4.0}
SMOKE = {"T": 14.0, "half_x": 6.0, "half_yz": 4.0, "span": 4.0, "dT_launch": 4.0, "x_start": 1.0, "jitter": 0.5,
         "region_hist": 3.0}
ARMS = ((0.0, 0.0, 4.0), (0.3, 0.0, 4.0), (0.3, 0.3, 4.0), (0.0, 0.3, 4.0),
        (0.0, 0.0, 6.0), (0.3, 0.0, 6.0), (0.3, 0.3, 6.0), (0.0, 0.3, 6.0))
SMOKE_R = {4.0: 3.0, 6.0: 3.5}                         # smaller cylinders in the smoke box
WALL = 0.5                                             # "at the wall": within this distance of the cylinder surface


def arm_key(v, u, R):
    return f"v{v:g}_u{u:g}_R{R:g}"


# ============================================================================================== the order

def sprinkle_flat_box(T, hx, hyz, rng):
    """Poisson events of density 1 in t [0, T], |x| <= hx, |y|, |z| <= hyz, sorted by time."""
    n = rng.poisson(T * (2 * hx) * (2 * hyz) ** 2)
    X = np.column_stack([rng.uniform(0, T, n), rng.uniform(-hx, hx, n),
                         rng.uniform(-hyz, hyz, n), rng.uniform(-hyz, hyz, n)])
    return X[np.argsort(X[:, 0], kind="stable")]


class MovingLinks(LocalLinks):
    """LocalLinks of F5.1c on a cylinder whose axis moves with the velocity u: |x - centre - u (t - t_c)| <= radius.

    For u = 0 it is LocalLinks exactly (the same events in the same order, the same links and bitsets); the body below
    is a copy of LocalLinks.__init__ with only the membership changed."""

    def __init__(self, order, centre, t_c, u, radius, t_lo, t_hi):
        X = order.X
        self.centre, self.t_c, self.u, self.radius = np.asarray(centre, float), float(t_c), np.asarray(u, float), \
            float(radius)
        axis = self.centre[None, :] + (X[:, :1] - self.t_c) * self.u[None, :]
        m = (((X[:, 1:] - axis) ** 2).sum(axis=1) <= radius * radius) & (X[:, 0] >= t_lo) & (X[:, 0] <= t_hi)
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
        self._X = X

    def dist(self, local):
        """Distance of local events from the axis at their own time."""
        P = self._X[self.ids[np.atleast_1d(local)]]
        axis = self.centre[None, :] + (P[:, :1] - self.t_c) * self.u[None, :]
        return np.sqrt(((P[:, 1:] - axis) ** 2).sum(axis=1))


# ============================================================================================== the walk

def walk_frame(order, ll, hist, t_end, m, H, rec):
    """walk_history of F5.1d (the registered rule, line by line) on a given region ll; rec collects, per step, the
    fraction of the horizon set Q and of the tie set of min L(p, q) at the wall, and the distance of the chosen
    target from the axis."""
    X = order.X
    pos = ll.pos
    path = list(hist)
    v = hist[-1]
    R = ll.radius
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
        Q1 = Q[sc == sc.min()]
        Q2 = Q1
        if Q2.size > 1:
            fp = np.packbits(Lp >= 0)
            nq = np.array([ll.between(fp, jp, int(q)) for q in Q2])
            Q2 = Q2[nq == nq.min()]
        jq = int(Q2.min())
        dQ, dQ1, dq = ll.dist(Q), ll.dist(Q1), float(ll.dist(jq)[0])
        rec.append({"fQ": float(np.mean(dQ > R - WALL)), "fQ1": float(np.mean(dQ1 > R - WALL)),
                    "q_wall": float(dq > R - WALL), "dq": dq / R, "nQ": int(Q.size), "n1": int(Q1.size)})
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


def _slope(s, y):
    """Least-squares slope of y(s); nan for fewer than 3 points."""
    if s.size < 3 or np.ptp(s) <= 0:
        return float("nan")
    return float(np.polyfit(s, y, 1)[0])


def measure(X, hist, path, d, span):
    """Velocities along the throw direction d: of the history, of the first three steps, of a quadratic fit
    (v0, a) and of the late half of the walk (s in [span/2, span])."""
    H0 = X[hist]
    v_hist = float(((H0[-1, 1:] - H0[0, 1:]) @ d) / max(H0[-1, 0] - H0[0, 0], 1e-9))
    P = X[path[len(hist) - 1:]]
    s = P[:, 0] - P[0, 0]
    y = (P[:, 1:] - P[0, 1:]) @ d
    j = min(len(P) - 1, 3)
    v_early = float(y[j] / s[j]) if j > 0 and s[j] > 0 else float("nan")
    if len(P) >= 4:
        c = np.polyfit(s, y, 2)
        a, v0 = float(2 * c[0]), float(c[1])
    else:
        a, v0 = float("nan"), float("nan")
    late = (s >= 0.5 * span) & (s <= span)
    return {"v_hist": v_hist, "v_early": v_early, "v0": v0, "a": a, "v_late": _slope(s[late], y[late]),
            "s_end": float(s[-1]), "n_steps": int(len(P) - 1)}


# ============================================================================================== one sprinkling

def seed_task(p, seed, n_starts, arms):
    t_start = time.time()
    rng = np.random.default_rng(seed)
    X = sprinkle_flat_box(float(p["T"]), float(p["half_x"]), float(p["half_yz"]), rng)
    order = CapacityOrder3(X, 0.0, float(p["width"]))
    m, H, span = int(p["memory"]), int(p["horizon"]), float(p["span"])
    t_hi = float(X[:, 0].max()) + 1.0                  # as in walk_history: the whole order above the start
    out = {"seed": int(seed), "N": int(len(X)), "walks": [], "reproduces_rule": None}
    for j in range(n_starts):
        sgn = 1.0 if j % 2 == 0 else -1.0              # throws along +x and -x in turn
        d = np.array([sgn, 0.0, 0.0])
        jit = rng.uniform(-float(p["jitter"]), float(p["jitter"]), 2)
        tgt0 = np.array([float(p["t_launch"]) + 0.5 * rng.random(), -sgn * float(p["x_start"]), jit[0], jit[1]])
        s_ = int(np.argmin(((X - tgt0) ** 2).sum(axis=1)))
        hists = {}
        for v in sorted({a[0] for a in arms}):         # one history per velocity, the same for all its arms
            dT = float(p["dT_launch"])
            tgt = X[s_].copy()
            tgt[0] += dT
            tgt[1:] += v * dT * d
            cand = np.where(np.abs(X[:, 0] - tgt[0]) < 1.0)[0]
            q = int(cand[np.argmin(((X[cand] - tgt) ** 2).sum(axis=1))])
            ch = chain_between(order, s_, q, float(p["region_hist"]))
            if ch is not None and len(ch) >= 3:
                hists[v] = ch[: len(ch) // 2 + 1]
        for (v, u, R) in arms:
            if v not in hists:
                continue
            hist = hists[v]
            ve = hist[-1]
            t0 = float(X[ve, 0])
            ll = MovingLinks(order, X[ve, 1:], t0, u * d, R, float(X[hist[0], 0]) - 1.0, t_hi)
            rec = []
            path = walk_frame(order, ll, hist, t0 + span + 1.0, m, H, rec)
            if out["reproduces_rule"] is None and u == 0.0 and v > 0.0:
                out["reproduces_rule"] = bool(walk_history(order, hist, t0 + span + 1.0, m, H, R) == path)
            w = {"j": j, "v": v, "u": u, "R": R, "arm": arm_key(v, u, R)}
            w.update(measure(X, hist, path, d, span))
            for k in ("fQ", "fQ1", "q_wall", "dq", "nQ", "n1"):
                w[k] = float(np.mean([r[k] for r in rec])) if rec else float("nan")
            out["walks"].append(w)
    out["seconds"] = time.time() - t_start
    rd = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "rows")
    os.makedirs(rd, exist_ok=True)
    pre = "smoke_inertia" if p.get("smoke") else "inertia"   # smoke rows never enter --summarize
    with open(os.path.join(rd, f"{pre}_{seed}.json"), "w", encoding="utf-8") as f:
        json.dump(out, f)
    return out


# ============================================================================================== summary

def _cluster(rows, key, arm, other=None):
    """Mean over sprinklings of the per-sprinkling mean of `key` in `arm` (or of the paired difference arm - other,
    walk by walk with the same start), and its SE over sprinklings."""
    per = []
    for r in rows:
        W = {(w["j"], w["arm"]): w for w in r["walks"]}
        vals = []
        for (j, a), w in W.items():
            if a != arm or not np.isfinite(w[key]):
                continue
            if other is None:
                vals.append(w[key])
            elif (j, other) in W and np.isfinite(W[(j, other)][key]):
                vals.append(w[key] - W[(j, other)][key])
        if vals:
            per.append(float(np.mean(vals)))
    per = np.asarray(per, float)
    if per.size == 0:
        return float("nan"), float("nan"), 0
    se = float(per.std(ddof=1) / math.sqrt(per.size)) if per.size > 1 else float("nan")
    return float(per.mean()), se, int(per.size)


def summary(p, rows, arms):
    rep = [r["reproduces_rule"] for r in rows if r["reproduces_rule"] is not None]
    L = [f"Inertia and the frame of the region, flat 3+1, {len(rows)} sprinklings of ~{np.mean([r['N'] for r in rows]):.0f}"
         f" events, rule of F5.1e (m {p['memory']}, H {p['horizon']}), walk span {p['span']:g}",
         f"  the cylinder at rest reproduces the registered rule: {all(rep) if rep else 'not checked'} ({len(rep)} checks)",
         "  arm (v, u, R)      walks  v_hist   v_early  v0       a         v_late           at the wall: Q  ties  q*   d(q*)/R",
         "  (+- = SE over sprinklings; velocities along the throw)"]
    for (v, u, R) in arms:
        k = arm_key(v, u, R)
        n = sum(1 for r in rows for w in r["walks"] if w["arm"] == k)
        c = {x: _cluster(rows, x, k) for x in ("v_hist", "v_early", "v0", "a", "v_late", "fQ", "fQ1", "q_wall", "dq")}
        L.append(f"  ({v:.1f}, {u:.1f}, {R:g})  {n:6d}  {c['v_hist'][0]:+.3f}   {c['v_early'][0]:+.3f}   "
                 f"{c['v0'][0]:+.3f}   {c['a'][0]:+.4f}   {c['v_late'][0]:+.3f} +- {c['v_late'][1]:.3f}   "
                 f"{c['fQ'][0]:.2f}  {c['fQ1'][0]:.2f}  {c['q_wall'][0]:.2f}  {c['dq'][0]:.2f}")
    L.append("  paired differences of v_late (same start, same history where v is the same):")
    for R in sorted({a[2] for a in arms}):
        for name, a1, a2 in (("thrown: co-moving - at rest", (0.3, 0.3, R), (0.3, 0.0, R)),
                             ("at rest: pulled - at rest  ", (0.0, 0.3, R), (0.0, 0.0, R)),
                             ("cylinder at rest: thrown - at rest", (0.3, 0.0, R), (0.0, 0.0, R))):
            if a1 in arms and a2 in arms:
                m_, s_, n_ = _cluster(rows, "v_late", arm_key(*a1), arm_key(*a2))
                L.append(f"    R {R:g}  {name}: {m_:+.3f} +- {s_:.3f}  ({n_} sprinklings)")
    L.append("  inertia would give v_late = v_hist in every arm; a cylinder acting as the medium gives v_late -> u")
    return L


# ============================================================================================== main

def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--seeds", type=int, default=8, help="flat sprinklings")
    ap.add_argument("--starts", type=int, default=16, help="starts per sprinkling (each runs every arm)")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--base-seed", type=int, default=0)
    ap.add_argument("--summarize", action="store_true", help="only summarise diag/out/rows/inertia_*.json")
    ap.add_argument("--smoke", action="store_true", help="a small box only to check that the code runs; its numbers "
                    "mean nothing")
    a = ap.parse_args()
    p = dict(PARAMS)
    arms = ARMS
    if a.smoke:
        p.update(SMOKE, smoke=True)
        arms = tuple((v, u, SMOKE_R[R]) for (v, u, R) in ARMS)
    if a.summarize:
        rd = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "rows")
        rows = []
        for fn in sorted(glob.glob(os.path.join(rd, "inertia_*.json"))):
            with open(fn, encoding="utf-8") as f:
                rows.append(json.load(f))
        print("\n".join(summary(p, rows, arms)))
        return
    jobs = [(p, derive_seed(a.base_seed, "diag-inertia-frame", "smoke" if a.smoke else "main", k), a.starts, arms)
            for k in range(a.seeds)]
    t0 = time.time()
    rows = parallel_map(seed_task, jobs, a.workers, print,
                        lambda i, r, dn, n: f"[diag] {dn}/{n}: N {r['N']}, {len(r['walks'])} walks "
                                            f"({time.time() - t0:.0f}s)")
    lines = summary(p, rows, arms)
    print("\n".join(lines))
    od = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
    os.makedirs(od, exist_ok=True)
    fn = os.path.join(od, f"inertia_frame{'_smoke' if a.smoke else ''}_{time.strftime('%Y%m%dT%H%M%S')}.json")
    with open(fn, "w", encoding="utf-8") as f:
        json.dump({"params": p, "arms": [list(x) for x in arms], "args": vars(a), "summary": lines, "rows": rows}, f)
    print(f"saved: {fn}  ({time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main()
