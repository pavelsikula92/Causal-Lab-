"""
F5.1-diag3: does the longest-chain rule keep its velocity when no region cuts its diamonds? (development; decides
nothing; the author's factorial design of 8 Nov, revised after the geometry of the diamonds)
=====================================================================================================================

Why
---
F5.1-diag2 showed that in flat 3+1 the history continued by the rule of F5.1c/d/e takes the velocity of the cylinder in
which the rule reads the order (co-moving - at rest +0.35 +- 0.04, pulled - at rest +0.28 +- 0.05 at R = 4), and the
chosen target lies at the wall in 3/4 of the steps. Geometry explains why R = 6 did not help: a target of min L(p, q)
has L(p, q) = m + H = 16 links, i.e. tau(p, q) ~ 13, so the diamond I(p, q) has a waist of radius ~6.5; every diamond is
cut at R = 4 and at R = 6, and both minima of the rule (L(p, q), then the number of events n(p, q)) are lowered by the
cut. A guard band (C) or complete diamonds (B) need a region larger than the diamonds; with the bitsets of LocalLinks
that is ~1.5 GB per process. This script computes the same rule without bitsets.

The engine (exactly the registered rule, a different computation)
----------------------------------------------------------------
At every step the order is read in the future cone of the memory event p (in the flat order the relation lies inside
the light cone, so nothing reachable from p is lost), up to a time cap t_v + cap. Longest chains are the vectorised
relaxations of LocalLinks. The number of events between p and each tied target q is counted by propagating one bit per
target backwards along the links inside the future of p (64 targets per machine word); the tie of the step to a child
the same way forwards inside the past of q*. Earliest = the smallest global index (the order is sorted by time), as in
LocalLinks. In the mode "cyl" (cylinder R = 4 at rest, no cap) the engine must give exactly the path of walk_history;
that is checked on every sprinkling.

Arms (flat 3+1, density 1, box t [0, 26], |x| <= 11, |y|, |z| <= 9.5; thrown histories v = 0 and 0.3 along +-x; the same
start, history and sprinkling in all arms of one start):
    A0      the registered region: cylinder R = 4 at rest around the start              (reference, F5.1-diag2)
    Binf    the whole box; a candidate q is admissible only if its diamond I(p, q) lies inside the box
    Bplain  the whole box, every candidate (the box cuts only far, boosted diamonds)
    C2      the whole box; admissible only if the diamond keeps >= 2 from the faces of the box
The diamond of p, q projects onto the spheroid |x - x_p| + |x - x_q| <= t_q - t_p, whose extent along each axis is
computed exactly. The horizon set is limited by the time cap t_v + 9 in all box arms; the share of ties within 0.5 of
the cap is recorded (the cap is a boundary of the laboratory too).

Frozen before the run
---------------------
  primary   the late velocity along the throw (slope over s in [6, 12] after the start of the walk), as in F5.1-diag2;
  secondary per step: the horizon set Q, the admissible part and the share removed by the boundary, the ties of min
            L(p, q) and of the fewest events, the excess L(p, q) - L(p, v) - L(v, q) of the admissible set (min, mean;
            0 = aligned with the history), n(p, q) of the ties
            (min, median), the velocity of the target relative to the history (u_q - u_hist), the margin of the
            diamond of q* to the box, the share of ties at the time cap;
  readings  inertia is restored in Binf if |v_late - v_hist| + 2 SE < 0.1 for both v (SE over sprinklings);
            admissibility matters if |Binf - Bplain| (paired v_late) > 2 SE;
            the boundary still reaches the decision if |C2 - Binf| > 2 SE;
            if Binf does not keep the velocity, the loss of inertia is in the rule itself (a new hypothesis).
The diagnostic decides no registered test (R132); F5.1f would be a new card with its own criteria.

Run from the root of the laboratory:
    python -m diag.inertia_box --seeds 8 --starts 8 --workers 6         (~1 h per sprinkling and ~2 GB per process)
    python -m diag.inertia_box --smoke
Each sprinkling is saved to diag/out/rows/box_*.json when done; --summarize reads them after a stopped run.
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
from phases.F5_gravity.F5_01c_free_fall_longest import LocalLinks
from phases.F5_gravity.F5_01d_free_fall_apex import walk_history
from phases.F7_light.F7_01c_throughput_cones_3d import CapacityOrder3
from diag.fall_layers import chain_between
from diag.inertia_frame import sprinkle_flat_box, measure, _cluster

# the box: the waist of a diamond (~6.5) + the jitter of the starts (0.5) + the guard of C2 (2) fits in |y|, |z| <= 9.5;
# along x the walks move by <= 5; the cap t_v + 9 stays below T = 26 (the walks end at t ~ 17). The order of ~206 000
# events is the cost (building it grows ~N^1.9).
PARAMS = {"T": 26.0, "half_x": 11.0, "half_yz": 9.5, "width": 4.0, "memory": 8, "horizon": 8, "span": 12.0,
          "t_launch": 1.0, "dT_launch": 6.0, "x_start": 2.5, "jitter": 0.5, "region_hist": 4.0, "R_cyl": 4.0,
          "cap": 9.0}
SMOKE = {"T": 14.0, "half_x": 6.0, "half_yz": 5.0, "memory": 4, "horizon": 4, "span": 4.0, "dT_launch": 4.0,
         "x_start": 1.0, "jitter": 0.5, "region_hist": 3.0, "R_cyl": 3.0, "cap": 5.0}
ARMS = (("A0", "cyl", 0.0), ("Binf", "adm", 0.0), ("Bplain", "plain", 0.0), ("C2", "adm", 2.0))
V_THROW = (0.0, 0.3)
_SWEEP = LocalLinks._sweep                               # the relaxation of LocalLinks (static)
_ONE = np.uint64(1)


# ============================================================================================== the window

class Window:
    """The kept links of the order among the events sel (global indices, ascending = time order), in local indices:
    f_u -> f_w sorted by the later end (forward sweeps), b_u -> b_w sorted by the earlier end (backward sweeps), and
    the children of every event (CSR) for the step."""

    def __init__(self, order, sel):
        self.ids, self.n = sel, sel.size
        lo = order.pptr[sel]
        cnt = order.pptr[sel + 1] - lo
        tot = int(cnt.sum())
        child = np.repeat(np.arange(self.n), cnt)
        offs = np.repeat(lo - np.concatenate([[0], np.cumsum(cnt)[:-1]]), cnt) + np.arange(tot)
        par_g = order.pidx[offs]                          # parents of every event of sel
        par = np.minimum(np.searchsorted(sel, par_g), max(self.n - 1, 0))
        ok = sel[par] == par_g                            # keep links with both ends in sel
        self.f_u, self.f_w = par[ok], child[ok]           # already sorted by the child
        o = np.argsort(self.f_u, kind="stable")
        self.perm_b = o                                   # forward order -> backward order of the links
        self.b_u, self.b_w = self.f_u[o], self.f_w[o]
        self.c_ptr = np.concatenate([[0], np.cumsum(np.bincount(self.b_u, minlength=self.n))]).astype(np.int64)
        self.c_idx = self.b_w

    def loc(self, g):
        """Local index of the global event g, -1 if it is not in the window."""
        j = int(np.searchsorted(self.ids, g))
        return j if j < self.n and self.ids[j] == g else -1


def _propagate(n, seeds, frm, to, alive):
    """Bit masks (one bit per seed, uint64 words) propagated along links frm -> to (sorted by to) inside the events
    with alive = True; a node ends with the bits of every seed it reaches (or is reached from) through alive events."""
    W = (len(seeds) + 63) // 64
    M = np.zeros((n, W), dtype=np.uint64)
    for k, s in enumerate(seeds):
        M[s, k >> 6] |= _ONE << np.uint64(k & 63)
    ok = alive[frm] & alive[to]
    a, b = frm[ok], to[ok]
    if b.size == 0:
        return M
    starts = np.flatnonzero(np.r_[True, b[1:] != b[:-1]])
    tg = b[starts]
    while True:
        new = M[tg] | np.bitwise_or.reduceat(M[a], starts, axis=0)
        if np.array_equal(new, M[tg]):
            return M
        M[tg] = new


def _counts(M, nodes, K):
    """Number of nodes carrying each of the K bits."""
    out = np.zeros(K, dtype=np.int64)
    sh = np.arange(64, dtype=np.uint64)
    for w in range(M.shape[1]):
        col = M[nodes, w]
        c = ((col[:, None] >> sh[None, :]) & _ONE).sum(axis=0)
        k = min(64, K - 64 * w)
        out[64 * w:64 * w + k] = c[:k]
    return out


def diamond_margin(xp, tp, Xq, half):
    """Smallest distance of the spatial projection of the diamonds I(p, q) (spheroids with foci x_p, x_q and the sum
    of distances t_q - t_p) from the faces of the box |x_i| <= half_i; negative = the diamond leaves the box."""
    xq, tq = Xq[:, 1:], Xq[:, 0]
    a = 0.5 * (tq - tp)
    dv = xq - xp[None, :]
    c = 0.5 * np.linalg.norm(dv, axis=1)
    f = np.where(c[:, None] > 0, dv / np.maximum(2 * c[:, None], 1e-300), 0.0)
    b2 = np.maximum(a * a - c * c, 0.0)
    ext = np.sqrt(a[:, None] ** 2 * f ** 2 + b2[:, None] * (1.0 - f ** 2))
    cen = 0.5 * (xq + xp[None, :])
    return np.min(half[None, :] - np.abs(cen) - ext, axis=1)


# ============================================================================================== the walk

def walk_box(order, hist, t_end, m, H, mode, guard, p, d, rec):
    """The registered rule of walk_history (F5.1d) computed without bitsets. mode: "cyl" (cylinder R_cyl at rest
    around the start, no cap: the registered region), "plain" (whole box, cap), "adm" (whole box, cap, only targets
    whose diamond keeps >= guard from the faces of the box)."""
    X = order.X
    Tt = X[:, 0]
    half = np.array([p["half_x"], p["half_yz"], p["half_yz"]], float)
    path = list(hist)
    v = hist[-1]
    centre = X[v, 1:].copy()
    R = float(p["R_cyl"])
    while X[v, 0] < t_end:
        src = path[len(path) - 1 - m] if len(path) > m else path[0]
        tp, xp = float(X[src, 0]), X[src, 1:]
        t_hi = float(Tt[-1]) + 1.0 if mode == "cyl" else float(X[v, 0]) + float(p["cap"])
        i0, i1 = int(np.searchsorted(Tt, tp, "left")), int(np.searchsorted(Tt, t_hi, "right"))
        idx = np.arange(i0, i1)
        P = X[idx]
        keep = np.sqrt(((P[:, 1:] - xp) ** 2).sum(axis=1)) <= (P[:, 0] - tp) + 1e-12   # future cone of p
        if mode == "cyl":
            keep &= ((P[:, 1:] - centre) ** 2).sum(axis=1) <= R * R
        sel = idx[keep]
        win = Window(order, sel)
        jv, jp = win.loc(v), win.loc(src)
        if jv < 0 or jp < 0:
            break
        Y = X[win.ids]
        cone_v = np.sqrt(((Y[:, 1:] - X[v, 1:]) ** 2).sum(axis=1)) <= (Y[:, 0] - X[v, 0]) + 1e-12
        ev = cone_v[win.f_u] & cone_v[win.f_w]            # links inside the future cone of v
        Lv = _SWEEP(win.n, jv, win.f_u, win.f_w, ev & (win.f_u >= jv))
        Q = np.flatnonzero(Lv == H)                       # 1. the horizon H links ahead
        if Q.size == 0:
            Hd = int(Lv.max())
            if Hd < 2:
                break
            Q = np.flatnonzero(Lv == Hd)
        Lp = _SWEEP(win.n, jp, win.f_u, win.f_w, win.f_u >= jp)
        Q = Q[Lp[Q] >= 0]
        if Q.size == 0:
            break
        nQ = int(Q.size)
        at_cap = Y[Q, 0] > t_hi - 0.5 if mode != "cyl" else np.zeros(Q.size, bool)
        marg = diamond_margin(xp, tp, Y[Q], half)
        if mode == "adm":                                 # only diamonds inside the box (with the guard band)
            okq = marg >= guard
            Q, marg, at_cap = Q[okq], marg[okq], at_cap[okq]
            if Q.size == 0:
                rec.append({"stop": 1.0})
                break
        sc = Lp[Q]                                        # 2. no boost: the shortest longest chain from p
        exc = sc - Lp[jv] - Lv[Q]                         # excess over the chain through v (0 = aligned with p -> v)
        t1 = sc == sc.min()
        Q1 = Q[t1]
        Q2 = Q1
        n_med, n_min = float("nan"), float("nan")
        if Q1.size > 1:                                   # then the fewest events between p and q
            M = _propagate(win.n, list(Q1), win.b_w, win.b_u, Lp >= 0)
            nq = np.maximum(_counts(M, np.flatnonzero(Lp >= 0), Q1.size) - 2, 0)
            Q2 = Q1[nq == nq.min()]
            n_med, n_min = float(np.median(nq)), float(nq.min())
        jq = int(Q2.min())                                # earliest (global indices are time-ordered)
        Lb = _SWEEP(win.n, jq, win.b_w, win.b_u, ev[win.perm_b] & (win.b_w <= jq))
        kids = win.c_idx[win.c_ptr[jv]:win.c_ptr[jv + 1]]  # 3. one step along a longest chain to q*
        kids = kids[Lb[kids] == Lb[jv] - 1]
        if kids.size == 0:
            break
        if kids.size > 1:                                 # the child with the most events up to q*
            alive = (Lb >= 0) & cone_v
            Mk = _propagate(win.n, list(kids), win.f_u, win.f_w, alive)
            nk = _counts(Mk, np.flatnonzero(alive), kids.size)
            kids = kids[nk == nk.max()]
        xv = X[v]
        u_hist = float(((xv[1:] - xp) @ d) / max(xv[0] - tp, 1e-12))
        yq = Y[jq]
        u_q = float(((yq[1:] - xv[1:]) @ d) / max(yq[0] - xv[0], 1e-12))
        rec.append({"nQ": nQ, "n_adm": int(Q.size), "removed": 1.0 - Q.size / nQ, "n1": int(Q1.size),
                    "n2": int(Q2.size), "dL_min": float(exc.min()), "dL_mean": float(exc.mean()),
                    "n_min": n_min, "n_med": n_med, "du": u_q - u_hist, "marg_q": float(marg[Q == jq][0]),
                    "cap1": float(at_cap[t1].mean()), "stop": 0.0, "n_win": int(win.n)})
        v = int(win.ids[int(kids.min())])
        path.append(v)
    return path


# ============================================================================================== one sprinkling

STEP_KEYS = ("nQ", "n_adm", "removed", "n1", "n2", "dL_min", "dL_mean", "n_min", "n_med", "du", "marg_q", "cap1",
             "stop", "n_win")


def seed_task(p, seed, n_starts):
    t_start = time.time()
    rng = np.random.default_rng(seed)
    X = sprinkle_flat_box(float(p["T"]), float(p["half_x"]), float(p["half_yz"]), rng)
    order = CapacityOrder3(X, 0.0, float(p["width"]))
    m, H, span = int(p["memory"]), int(p["horizon"]), float(p["span"])
    out = {"seed": int(seed), "N": int(len(X)), "walks": [], "checks": []}
    for j in range(n_starts):
        sgn = 1.0 if j % 2 == 0 else -1.0              # throws along +x and -x in turn
        d = np.array([sgn, 0.0, 0.0])
        jit = rng.uniform(-float(p["jitter"]), float(p["jitter"]), 2)
        tgt0 = np.array([float(p["t_launch"]) + 0.5 * rng.random(), -sgn * float(p["x_start"]), jit[0], jit[1]])
        s_ = int(np.argmin(((X - tgt0) ** 2).sum(axis=1)))
        for v in V_THROW:
            dT = float(p["dT_launch"])
            tgt = X[s_].copy()
            tgt[0] += dT
            tgt[1:] += v * dT * d
            cand = np.where(np.abs(X[:, 0] - tgt[0]) < 1.0)[0]
            q = int(cand[np.argmin(((X[cand] - tgt) ** 2).sum(axis=1))])
            ch = chain_between(order, s_, q, float(p["region_hist"]))
            if ch is None or len(ch) < 3:
                continue
            hist = ch[: len(ch) // 2 + 1]
            t0 = float(X[hist[-1], 0])
            for (name, mode, guard) in ARMS:
                rec = []
                path = walk_box(order, hist, t0 + span + 1.0, m, H, mode, guard, p, d, rec)
                if mode == "cyl" and len(out["checks"]) < 2:   # the engine must be the registered rule
                    ref = walk_history(order, hist, t0 + span + 1.0, m, H, float(p["R_cyl"]))
                    out["checks"].append(bool(ref == path))
                w = {"j": j, "v": v, "arm": f"{name}_v{v:g}", "name": name}
                w.update(measure(X, hist, path, d, span))
                for k in STEP_KEYS:
                    vals = [r[k] for r in rec if k in r and np.isfinite(r[k])]
                    w[k] = float(np.mean(vals)) if vals else float("nan")
                out["walks"].append(w)
    out["seconds"] = time.time() - t_start
    rd = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "rows")
    os.makedirs(rd, exist_ok=True)
    pre = "smoke_box" if p.get("smoke") else "box"     # smoke rows never enter --summarize
    with open(os.path.join(rd, f"{pre}_{seed}.json"), "w", encoding="utf-8") as f:
        json.dump(out, f)
    return out


# ============================================================================================== summary

def summary(p, rows):
    chk = [c for r in rows for c in r["checks"]]
    L = [f"F5.1-diag3: inertia without a cutting region, flat 3+1, {len(rows)} sprinklings of "
         f"~{np.mean([r['N'] for r in rows]):.0f} events, rule of F5.1e (m {p['memory']}, H {p['horizon']}), "
         f"cap {p['cap']:g}, span {p['span']:g}",
         f"  the engine reproduces the registered rule (cylinder R {p['R_cyl']:g}): {all(chk) if chk else 'not checked'}"
         f" ({sum(chk)} of {len(chk)})",
         "  arm      v    walks  v_hist  v_late (+- SE)    late-hist (+- SE)   Q     adm   removed ties  few  "
         "dL_min dL_mean n_min  du      margin cap  stop"]
    for (name, _, _) in ARMS:
        for v in V_THROW:
            k = f"{name}_v{v:g}"
            n = sum(1 for r in rows for w in r["walks"] if w["arm"] == k)
            if n == 0:
                continue
            for r in rows:                               # late - history per walk, for the inertia reading
                for w in r["walks"]:
                    w["late_minus_hist"] = w["v_late"] - w["v_hist"]
            c = {x: _cluster(rows, x, k) for x in ("v_hist", "v_late", "late_minus_hist") + STEP_KEYS}
            L.append(f"  {name:7s} {v:.1f} {n:6d}  {c['v_hist'][0]:+.3f}  {c['v_late'][0]:+.3f} +- {c['v_late'][1]:.3f}"
                     f"   {c['late_minus_hist'][0]:+.3f} +- {c['late_minus_hist'][1]:.3f}    {c['nQ'][0]:5.0f} "
                     f"{c['n_adm'][0]:5.0f}  {c['removed'][0]:.2f}   {c['n1'][0]:4.1f} {c['n2'][0]:4.2f} "
                     f"{c['dL_min'][0]:+.2f}  {c['dL_mean'][0]:+.2f}  {c['n_min'][0]:5.0f} {c['du'][0]:+.3f}  "
                     f"{c['marg_q'][0]:+.2f}  {c['cap1'][0]:.2f} {c['stop'][0]:.2f}")
    L.append("  paired differences of v_late (same start and history):")
    for v in V_THROW:
        for a1, a2 in (("Binf", "A0"), ("Binf", "Bplain"), ("C2", "Binf")):
            m_, s_, n_ = _cluster(rows, "v_late", f"{a1}_v{v:g}", f"{a2}_v{v:g}")
            L.append(f"    v {v:.1f}  {a1} - {a2}: {m_:+.3f} +- {s_:.3f}  ({n_} sprinklings)")
    for v in V_THROW:
        m_, s_, _ = _cluster(rows, "late_minus_hist", f"Binf_v{v:g}")
        ok = abs(m_) + 2 * s_ < 0.1 if np.isfinite(s_) else False
        L.append(f"  reading: Binf v {v:.1f} keeps its velocity (|late - hist| + 2 SE < 0.1): {ok}")
    return L


# ============================================================================================== main

def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--seeds", type=int, default=8, help="flat sprinklings")
    ap.add_argument("--starts", type=int, default=8, help="starts per sprinkling (each: v 0 and 0.3, every arm)")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--base-seed", type=int, default=0)
    ap.add_argument("--summarize", action="store_true", help="only summarise diag/out/rows/box_*.json")
    ap.add_argument("--smoke", action="store_true", help="a small box only to check that the code runs; its numbers "
                    "mean nothing")
    a = ap.parse_args()
    p = dict(PARAMS)
    if a.smoke:
        p.update(SMOKE, smoke=True)
    if a.summarize:
        rd = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "rows")
        rows = []
        for fn in sorted(glob.glob(os.path.join(rd, "box_*.json"))):
            with open(fn, encoding="utf-8") as f:
                rows.append(json.load(f))
        print("\n".join(summary(p, rows)))
        return
    jobs = [(p, derive_seed(a.base_seed, "diag-inertia-box", "smoke" if a.smoke else "main", k), a.starts)
            for k in range(a.seeds)]
    t0 = time.time()
    rows = parallel_map(seed_task, jobs, a.workers, print,
                        lambda i, r, dn, n: f"[diag] {dn}/{n}: N {r['N']}, {len(r['walks'])} walks "
                                            f"({time.time() - t0:.0f}s)")
    lines = summary(p, rows)
    print("\n".join(lines))
    od = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
    os.makedirs(od, exist_ok=True)
    fn = os.path.join(od, f"inertia_box{'_smoke' if a.smoke else ''}_{time.strftime('%Y%m%dT%H%M%S')}.json")
    with open(fn, "w", encoding="utf-8") as f:
        json.dump({"params": p, "args": vars(a), "summary": lines, "rows": rows}, f)
    print(f"saved: {fn}  ({time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main()
