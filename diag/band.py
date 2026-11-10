"""
F5.1-diag4: is the set of admissible continuations of a history finite without a wall? (development; decides nothing)
====================================================================================================================

Why (9 Nov 2026, the author's reading of superposition)
------------------------------------------------------
Every rule tried so far chose one continuation of a history, and the most extreme one (fewest or most events); that
choice is a measurement in every step, and an extreme of a large random set lies at its edge, which the edge of the
region then decides (F5.1-diag2, -diag3, the withdrawn draft of F5.1f). The author's reading: only what has happened
is fixed, i.e. the memory event p and the last event v of the history; the next layer stays a superposition of the
admissible continuations, and the flow is shared among them instead of being sent along one line. A rule of this kind
is finite (M5) only if the admissible set is finite: the children c of v with no deficit against the memory,

    L(p, c) = L(p, v) + 1      (the history through v stays a longest chain from p),

or, if no child reaches it, those with the smallest deficit L(p, c) - L(p, v) - 1. In the continuum a child at a
distance d along the light cone of v has a deficit of about d, so the set should be confined to about a link. This
diagnostic measures it.

What it measures
----------------
Flat order (2+1 or 3+1, density 1), histories that are the first half of a longest chain (chain_between, cylinder radius
4 in 2+1 / 3 in 3+1 relative to the launch) from a start to the event nearest to the start shifted by v dT along +-x,
dT later; v = 0 and 0.3; memory m = 8 links (p = the event m steps back in the history). For every history the
admissible set is computed in windows of growing size: the future cone of p, up to t(v) + dt_w, dt_w in DT_W, inside
a box of half size h, two sizes; a child counts only if its diamond I(p, c) lies inside the box (then L is exact).
Reported per window: the children of v in the window, the admissible children, the smallest deficit and how often it is
> 0, the size of the admissible set, its largest time and spatial offset from v, and the largest rapidity of a child
relative to the history (from the velocity of p -> v); and the saturation: how often the admissible set in a window
equals the set in the largest window (the same box), and how often the sets of the two boxes agree.

Run from the root of the laboratory:
    python -m diag.band --dim 2 --seeds 4 --starts 40 --workers 2 --dt 6,8,10,13 --T 40
    python -m diag.band --dim 3 --seeds 8 --starts 20 --workers 8
Output: a summary on the screen and the raw records in diag/out/band_d{dim}_{stamp}.json.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import time

import numpy as np

from core.order.sprinkle import derive_seed
from core.stats.parallel import parallel_map
from phases.F5_gravity.F5_window import Window, diamond_margin, chain_between, _SWEEP
from phases.F7_light.F7_01a_capacity_cones import CapacityOrder
from phases.F7_light.F7_01c_throughput_cones_3d import CapacityOrder3

PARAMS = {
    2: {"T": 34.0, "half": 16.0, "half_small": 12.0, "dT_launch": 22.0, "region_hist": 4.0, "t_launch": 1.0,
        "ring": 2.0, "memory": 8, "dt_w": (1.5, 2.5, 4.0, 6.0, 8.0)},
    3: {"T": 20.0, "half": 8.0, "half_small": 6.0, "dT_launch": 16.0, "region_hist": 3.0, "t_launch": 1.0,
        "ring": 1.5, "memory": 8, "dt_w": (2.0, 3.0, 4.5, 6.0, 8.0)},
}
SMOKE = {2: {"T": 20.0, "half": 8.0, "half_small": 6.0, "dT_launch": 10.0, "region_hist": 3.0, "memory": 4,
             "dt_w": (1.0, 2.0, 3.0)},
         3: {"T": 14.0, "half": 5.0, "half_small": 4.0, "dT_launch": 8.0, "region_hist": 3.0, "memory": 4,
             "dt_w": (1.0, 2.0)}}
V_THROW = (0.0, 0.3)


def sprinkle_flat(dim, T, h, rng):
    """Poisson events of density 1 in t [0, T], |x_i| <= h, sorted by time."""
    n = rng.poisson(T * (2 * h) ** dim)
    X = np.column_stack([rng.uniform(0, T, n)] + [rng.uniform(-h, h, n) for _ in range(dim)])
    return X[np.argsort(X[:, 0], kind="stable")]


def rapidity(u, k):
    """Rapidity of the timelike or null vector k relative to the unit timelike u (Minkowski, signature -+..);
    inf for a null or spacelike k."""
    kk = k[0] ** 2 - float(k[1:] @ k[1:])
    if kk <= 0:
        return float("inf")
    g = (u[0] * k[0] - float(u[1:] @ k[1:])) / math.sqrt(kk)
    return float(math.acosh(max(g, 1.0)))


def admissible(order, p, v, dt_w, half):
    """The admissible children of v in the window (future cone of p, t <= t(v) + dt_w, box |x_i| <= half), counting
    only children whose diamond I(p, c) lies inside the box; returns (children in window, admissible, min deficit,
    set of the admissible global ids with the smallest deficit)."""
    X = order.X
    Tt = X[:, 0]
    tp, xp = float(X[p, 0]), X[p, 1:]
    i0 = int(np.searchsorted(Tt, tp, "left"))
    i1 = int(np.searchsorted(Tt, float(X[v, 0]) + dt_w, "right"))
    idx = np.arange(i0, i1)
    P = X[idx]
    hv = np.full(X.shape[1] - 1, half)
    keep = (np.sqrt(((P[:, 1:] - xp) ** 2).sum(axis=1)) <= (P[:, 0] - tp) + 1e-12) & \
           np.all(np.abs(P[:, 1:]) <= hv[None, :], axis=1)
    win = Window(order, idx[keep])
    jv, jp = win.loc(v), win.loc(p)
    if jv < 0 or jp < 0:
        return 0, 0, None, set()
    Lp = _SWEEP(win.n, jp, win.f_u, win.f_w, win.f_u >= jp)
    kids = win.c_idx[win.c_ptr[jv]:win.c_ptr[jv + 1]]
    nk = int(kids.size)
    if kids.size:
        kids = kids[diamond_margin(xp, tp, X[win.ids[kids]], hv) >= 0.0]
    if kids.size == 0:
        return nk, 0, None, set()
    dfc = Lp[kids] - Lp[jv] - 1
    mn = int(dfc.min())
    return nk, int(kids.size), mn, set(int(g) for g in win.ids[kids[dfc == mn]])


def seed_task(dim, p, seed, n_starts):
    t_start = time.time()
    rng = np.random.default_rng(seed)
    X = sprinkle_flat(dim, float(p["T"]), float(p["half"]), rng)
    order = CapacityOrder(X, 0.0, 4.0) if dim == 2 else CapacityOrder3(X, 0.0, 4.0)
    m = int(p["memory"])
    rows, used = [], set()
    for v in V_THROW:
        phase = rng.uniform(0.0, 2.0 * math.pi)
        for j in range(n_starts):
            sgn = 1.0 if j % 2 == 0 else -1.0
            d = np.zeros(dim)
            d[0] = sgn
            ang = phase + 2.0 * math.pi * (j // 2) / max(n_starts // 2, 1)
            tgt0 = np.zeros(dim + 1)
            tgt0[0] = float(p["t_launch"]) + 0.5 * rng.random()
            tgt0[1] = -sgn * v * float(p["dT_launch"]) / 2.0
            if dim == 2:
                tgt0[2] = float(p["ring"]) * math.sin(ang)
            else:
                tgt0[2], tgt0[3] = float(p["ring"]) * math.cos(ang), float(p["ring"]) * math.sin(ang)
            near = np.argsort(((X - tgt0) ** 2).sum(axis=1))[:64]
            s_ = next(int(i) for i in near if int(i) not in used)
            used.add(s_)
            tgt = X[s_].copy()
            tgt[0] += float(p["dT_launch"])
            tgt[1:] += v * float(p["dT_launch"]) * d
            cand = np.where(np.abs(X[:, 0] - tgt[0]) < 1.0)[0]
            q = int(cand[np.argmin(((X[cand] - tgt) ** 2).sum(axis=1))])
            # the cylinder of chain_between must contain the target of the throw
            ch = chain_between(order, s_, q, max(float(p["region_hist"]), v * float(p["dT_launch"]) + 2.0))
            if ch is None or len(ch) < 4:
                continue
            hist = ch[: len(ch) // 2 + 1]
            vv = hist[-1]
            pp = hist[len(hist) - 1 - m] if len(hist) > m else hist[0]
            k_ph = X[vv] - X[pp]
            u = k_ph / math.sqrt(max(k_ph[0] ** 2 - float(k_ph[1:] @ k_ph[1:]), 1e-12))
            r = {"v": v, "mem": int(len(hist) - 1 - hist.index(pp)), "win": []}
            sets = {}
            for h in (float(p["half_small"]), float(p["half"])):
                for dtw in p["dt_w"]:
                    nk, na, mn, Z = admissible(order, pp, vv, float(dtw), h)
                    sets[(h, dtw)] = Z
                    off = [X[c] - X[vv] for c in Z]
                    r["win"].append({"h": h, "dt": float(dtw), "kids": nk, "adm": na, "min_def": mn,
                                     "size": len(Z),
                                     "dt_max": float(max(o[0] for o in off)) if off else float("nan"),
                                     "dx_max": float(max(np.linalg.norm(o[1:]) for o in off)) if off else float("nan"),
                                     "rap_max": float(max(rapidity(u, o) for o in off)) if off else float("nan")})
            big = (float(p["half"]), float(p["dt_w"][-1]))
            for w in r["win"]:
                w["same_as_largest"] = bool(sets[(w["h"], w["dt"])] == sets[big])
                w["same_other_box"] = bool(sets[(w["h"], w["dt"])] == sets[(
                    float(p["half"]) if w["h"] == float(p["half_small"]) else float(p["half_small"]), w["dt"])])
            rows.append(r)
    return {"seed": int(seed), "N": int(len(X)), "rows": rows, "seconds": time.time() - t_start}


def summary(dim, p, res):
    rows = [r for s in res for r in s["rows"]]
    L = [f"F5.1-diag4: admissible continuations of a history, flat {dim}+1, {len(res)} sprinklings of "
         f"~{np.mean([s['N'] for s in res]):.0f} events, memory {p['memory']} (actual mean "
         f"{np.mean([r['mem'] for r in rows]):.1f}), {len(rows)} histories",
         "  v    box   dt_w   children  admissible  min deficit (>0)   set size  dt_max  dx_max  rapidity_max   "
         "= largest window   = other box"]
    for v in V_THROW:
        R = [r for r in rows if r["v"] == v]
        if not R:
            continue
        for k in range(len(R[0]["win"])):
            W = [r["win"][k] for r in R]
            md = [w["min_def"] for w in W if w["min_def"] is not None]
            rap = np.array([w["rap_max"] for w in W], float)
            L.append(f"  {v:.1f}  {W[0]['h']:5.1f}  {W[0]['dt']:4.1f}   {np.mean([w['kids'] for w in W]):8.1f}  "
                     f"{np.mean([w['adm'] for w in W]):10.1f}  {np.mean(md) if md else float('nan'):+.2f} "
                     f"({np.mean([x > 0 for x in md]) if md else float('nan'):.2f})      "
                     f"{np.mean([w['size'] for w in W]):6.1f}  {np.nanmean([w['dt_max'] for w in W]):6.2f}  "
                     f"{np.nanmean([w['dx_max'] for w in W]):6.2f}  {np.nanmedian(rap):5.2f} (inf {np.mean(~np.isfinite(rap)):.2f})"
                     f"   {np.mean([w['same_as_largest'] for w in W]):.2f}             "
                     f"{np.mean([w['same_other_box'] for w in W]):.2f}")
    return L


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--dim", type=int, choices=(2, 3), required=True)
    ap.add_argument("--seeds", type=int, default=4)
    ap.add_argument("--starts", type=int, default=40, help="starts per sprinkling and velocity")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--base-seed", type=int, default=0)
    ap.add_argument("--smoke", action="store_true", help="a small box only to check that the code runs")
    ap.add_argument("--dt", type=str, default=None, help="comma-separated windows dt_w (default of the dimension)")
    ap.add_argument("--T", type=float, default=None, help="height of the box (default of the dimension)")
    a = ap.parse_args()
    p = dict(PARAMS[a.dim])
    if a.smoke:
        p.update(SMOKE[a.dim])
    if a.dt:
        p["dt_w"] = tuple(float(x) for x in a.dt.split(","))
    if a.T:
        p["T"] = a.T
    jobs = [(a.dim, p, derive_seed(a.base_seed, "diag-band", a.dim, "smoke" if a.smoke else "main", k), a.starts)
            for k in range(a.seeds)]
    t0 = time.time()
    res = parallel_map(seed_task, jobs, a.workers, print,
                       lambda i, r, dn, n: f"[diag] {dn}/{n}: N {r['N']}, {len(r['rows'])} histories "
                                           f"({time.time() - t0:.0f}s)")
    lines = summary(a.dim, p, res)
    print("\n".join(lines))
    od = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
    os.makedirs(od, exist_ok=True)
    fn = os.path.join(od, f"band_d{a.dim}{'_smoke' if a.smoke else ''}_{time.strftime('%Y%m%dT%H%M%S')}.json")
    with open(fn, "w", encoding="utf-8") as f:
        json.dump({"params": {k: list(v) if isinstance(v, tuple) else v for k, v in p.items()}, "args": vars(a),
                   "summary": lines, "res": res}, f)
    print(f"saved: {fn}  ({time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main()
