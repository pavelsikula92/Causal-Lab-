"""
F5.1f-a -- Inertia of a flow shared among its admissible continuations: does the centre of a free flow keep the
velocity of its history in a flat 3+1 order, and is the rule finite without a wall?
==================================================================================================================

Why this test exists (9 Nov 2026, the author's decisions)
--------------------------------------------------------
The fall rule of F5.1c/d/e is a proxy by the axiom audit: it reads the order H = 8 links in the future (M4) and its
result depends on the wall of the region (M5; F5.1-diag2, -diag3, -diag3c). A first draft of this test chose, among
the children of the last event, the one with the most events from the memory; a quick run showed that it depends on
the window (0 of 16 identical paths): maximising proper time towards a free endpoint has no maximum. The author's
reading of superposition (ontology: superposition is the coexistence of admissible continuations, measurement biases
the routing towards one of them) gives the rule of this test: only what has happened is fixed, the memory event p and
the last event v; the next layer is not chosen, the flow is shared among all admissible continuations, and what is read
is the centre of the flow on its realised past (M2: the reading does not act on the dynamics).

The admissible continuations are the children c of v with no deficit against the memory, L(p, c) = L(p, v) + 1 (the
history through v stays a longest chain from p; the smallest deficit if no child reaches it). F5.1-diag4 measured this
set: in 2+1 ~8 children within ~3 units, identical for windows dt_w = 8 to 13 and for two boxes; in 3+1 ~21 (v 0) and
~25 (v 0.3) children, identical in 0.96-0.99 of the histories between dt_w = 6 and 8 in a box |x| <= 8, but a box of 6
was too small for the diamonds (0.12-0.17 identical), so the spatial saturation in 3+1 is tested here (M5).

The rule (fixed before any run; it uses only the order and counts)
------------------------------------------------------------------
A flow with history ... -> p -> ... -> v, p the event m = 8 steps back (the first event while the history is shorter):
  the admissible continuations Z(p, v) = the children c of v with the smallest L(p, c) - L(p, v) - 1;
  the flow is shared equally among Z(p, v); every continuation carries its own history on.
The superposition is computed by K walkers per history: walker k continues to the member of Z(p, v) with the smallest
hash of (k, identity of the event), a choice uniform over Z that depends neither on coordinates nor on the order of
computation. The centre of the flow is the mean over the walkers.

Computation (not part of the rule)
----------------------------------
L(p, .) by the relaxation of LocalLinks inside the future cone of p up to t(v) + dt_w; a child counts only if its
diamond I(p, c) lies inside the box of the region (then L(p, c) is exact). Region B: |x| <= 10, |y|, |z| <= 8.5,
dt_w = 7; region A: |x| <= 8.5, |y|, |z| <= 7, dt_w = 5.5. The admissible sets of one (p, v) are cached.

Set-up
------
Flat 3+1 order of density 1 (CapacityOrder3, sigma0 = 0), box t [0, 27], |x| <= 10, |y|, |z| <= 8.5 (~156 000 events).
Histories: the first half of the longest chain (chain_between) from a start event to the event nearest to the start
shifted by v * 16 along +-x sixteen units later (~9 links, so the memory is full from the first step); 8 distinct starts
per sprinkling and velocity, centred so that history and walk are symmetric in the box, ring of radius 1.5 in (y, z);
v = 0 and 0.3; K = 8 walkers per history; walk 10 units; 8 sprinklings. Velocity of the history v_hist (first to last
event, along the throw), of a walker v_walk (least-squares slope along the throw over s in [0, 10]); Delta_h of a
history = the mean of v_walk over its walkers minus v_hist.
Instrument: the registered proxy (walk_history of F5.1d, cylinder of radius 4) on the same histories, whose braking is
known (F5.1-diag2/3/3b/3c).

Criteria (stated before any run; SE by bootstrap over sprinklings)
-----------------------------------------------------------------
INSTRUMENT_FAILED  if the estimator does not see the braking of the proxy: Delta_h(proxy, 0.3) + 2 SE >= -0.05.
finite (M5)        walker 0 of every history is also computed in region A; its path in A equals its path in B in
                   >= 95 % of the histories (both v together).
kept               |Delta_h(v)| + 2 SE < 0.05 for v = 0 and v = 0.3 (region B).
lost               Delta_h(0.3) + 2 SE < -0.05, or |Delta_h(0)| - 2 SE > 0.05 (a flow at rest starts to move).
PASSED if finite and kept; FAILED if not finite or lost; UNDECIDED otherwise.
Reported: the identical fraction per v; the size of Z per step and its smallest deficit; the spread of the flow (SD over
the walkers of the displacement along the throw and across it at the end, mean over histories); the per-walker SD of
v_walk; the share of walkers that stop.

Reading declared now. PASSED: a flow that keeps its past and shares its future among the admissible continuations
keeps its velocity without any wall; it is the candidate for the fall (F5.1f-b, 3+1, criteria of F5.1e, the same rule
on the loaded order). FAILED, not finite: the admissible set in 3+1 is not confined by min L in space; a covariant
confinement is needed. FAILED, lost: equal sharing does not keep the centre; the weights or the memory are wrong (a new
hypothesis). UNDECIDED: the precision does not decide. The spread is the first measurement of how fast a free flow
without binding disperses.
Limits known in advance: one memory m = 8; equal weights; the launch history comes from chain_between (stepping to the
earliest child); flat order only.
"""

from __future__ import annotations

import math
import time

import numpy as np

from core.order.sprinkle import derive_seed
from core.stats.parallel import parallel_map
from core.stats.stats import boot_se
from core.testing import Result
from phases.F5_gravity.F5_01d_free_fall_apex import walk_history
from phases.F5_gravity.F5_window import Window, diamond_margin, chain_between, id_hash, _SWEEP
from phases.F7_light.F7_01a_capacity_cones import CapacityOrder
from phases.F7_light.F7_01c_throughput_cones_3d import CapacityOrder3

TEST = {
    "id": "F5.1f",
    "phase": "F5",
    "title": "Inertia of a flow shared among its admissible continuations (F5.1f-a): does the centre of a free flow "
             "keep its velocity in a flat 3+1 order, and is the rule finite without a wall?",
    "theory_part": "F5 gravity and equivalence; ontology (superposition = coexisting admissible continuations), Axiom "
                   "of the action (longest chains), Axiom of momentum (memory of the flow), M2, M4, M5",
    "question": "Does a flow whose history is fixed up to its last event v and whose next layer is shared equally among "
                "the children c of v with L(p, c) = L(p, v) + 1 (p = 8 steps back) keep the velocity of its history "
                "in its centre, and does the rule not depend on the region in which the order is read?",
    "method": "flat 3+1 order, density 1, box t [0, 27], |x| <= 10, |y|, |z| <= 8.5; histories: first half of the "
              "longest chain to the start shifted by v * 16 along +-x sixteen units later; 8 distinct starts per "
              "sprinkling and v, v = 0 and 0.3, K = 8 walkers per history (uniform choice by a hash), walk 10 units, "
              "8 sprinklings; region B (|x| <= 10, |y|, |z| <= 8.5, dt_w 7), walker 0 also in region A (|x| <= 8.5, "
              "|y|, |z| <= 7, dt_w 5.5); only children whose diamond lies inside; Delta_h = mean slope of the walkers "
              "over [0, 10] minus the velocity of the history; instrument: the proxy of F5.1d (cylinder 4)",
    "prediction": "a flow that keeps its past and shares its future keeps its velocity (Delta_h = 0) whatever the "
                  "region",
    "criterion": "INSTRUMENT_FAILED if Delta_h(proxy, 0.3) + 2 SE >= -0.05. finite: path of walker 0 in A = in B in "
                 ">= 95 % of histories. kept: |Delta_h(v)| + 2 SE < 0.05 for v = 0 and 0.3. lost: Delta_h(0.3) + 2 SE "
                 "< -0.05 or |Delta_h(0)| - 2 SE > 0.05. PASSED finite and kept; FAILED not finite or lost; UNDECIDED "
                 "otherwise.",
    "null_model": "the proxy of F5.1d in a cylinder at rest, which brakes a thrown history (F5.1-diag2/3/3b/3c)",
    "level": "mechanism",
    "meters": [],
    "depends_on": ["F5.1e", "F5.1-diag4", "F3.6"],
    "params": {"T": 27.0, "half_x": 10.0, "half_yz": 8.5, "half_x_A": 8.5, "half_yz_A": 7.0, "dt_A": 5.5,
               "dt_B": 7.0, "width": 4.0, "memory": 8, "span": 10.0, "t_launch": 1.0, "t_spread": 0.75,
               "dT_launch": 16.0, "ring": 1.5, "region_hist": 3.0, "R_cyl": 4.0, "proxy_H": 8, "n_starts": 8,
               "n_walkers": 8, "n_seeds": 8, "n_se": 2.0, "tol": 0.05, "ident": 0.95, "boot_B": 4000},
    "parts": {"main": {}},
    "quick": {"n_seeds": 2, "n_starts": 2, "n_walkers": 3, "boot_B": 200, "T": 16.0, "half_x": 6.0,
              "half_yz": 5.0, "half_x_A": 5.0, "half_yz_A": 4.0, "dt_A": 2.0, "dt_B": 3.0, "memory": 4,
              "proxy_H": 4, "span": 4.0, "dT_launch": 8.0, "ring": 0.8, "R_cyl": 3.0},
    "estimate": "~3 h with 8 processes (8 orders of ~156 000 events; per order 128 walkers of ~12 steps in region B, "
                "16 in region A and 16 proxy walks); ~2 GB per process",
}

V_THROW = (0.0, 0.3)
STEP_KEYS = ("size", "min_def", "n_kids", "n_adm", "n_win")


# ============================================================================================== the rule

def admissible(order, p, v, half, dt_w):
    """Z(p, v): the children c of v with the smallest L(p, c) - L(p, v) - 1, read in the future cone of p up to
    t(v) + dt_w inside the box |x_i| <= half_i, counting only children whose diamond I(p, c) lies inside the box.
    Returns (global ids of Z, statistics)."""
    X = order.X
    Tt = X[:, 0]
    tp, xp = float(X[p, 0]), X[p, 1:]
    i0 = int(np.searchsorted(Tt, tp, "left"))
    i1 = int(np.searchsorted(Tt, float(X[v, 0]) + dt_w, "right"))
    idx = np.arange(i0, i1)
    P = X[idx]
    keep = (np.sqrt(((P[:, 1:] - xp) ** 2).sum(axis=1)) <= (P[:, 0] - tp) + 1e-12) & \
           np.all(np.abs(P[:, 1:]) <= half[None, :], axis=1)
    win = Window(order, idx[keep])
    jv, jp = win.loc(v), win.loc(p)
    if jv < 0 or jp < 0:
        return np.zeros(0, dtype=np.int64), {"stop": 1.0}
    Lp = _SWEEP(win.n, jp, win.f_u, win.f_w, win.f_u >= jp)
    kids = win.c_idx[win.c_ptr[jv]:win.c_ptr[jv + 1]]
    n_kids = int(kids.size)
    if kids.size:
        kids = kids[diamond_margin(xp, tp, X[win.ids[kids]], half) >= 0.0]
    if kids.size == 0:
        return np.zeros(0, dtype=np.int64), {"stop": 1.0}
    dfc = Lp[kids] - Lp[jv] - 1
    mn = int(dfc.min())
    Z = np.sort(win.ids[kids[dfc == mn]]).astype(np.int64)
    return Z, {"stop": 0.0, "size": int(Z.size), "min_def": mn, "n_kids": n_kids, "n_adm": int(kids.size),
               "n_win": int(win.n)}


def walker(order, hist, t_end, m, half, dt_w, key, cache, rec):
    """One walker of the shared flow: at every step the member of Z(p, v) with the smallest hash of (key, identity)."""
    X = order.X
    path = list(hist)
    v = hist[-1]
    kk = np.uint64(key)
    while X[v, 0] < t_end:
        p = path[len(path) - 1 - m] if len(path) > m else path[0]
        ck = (int(p), int(v))
        if ck not in cache:
            cache[ck] = admissible(order, p, v, half, dt_w)
        Z, st = cache[ck]
        rec.append(st)
        if Z.size == 0:
            break
        v = int(Z[np.argmin(id_hash(Z.astype(np.uint64) ^ kk))])
        path.append(v)
    return path


# ============================================================================================== readouts

def velocities(X, hist, path, d, span):
    """v_hist from the first to the last event of the history; v_walk the least-squares slope along d over s in
    [0, span] from the end of the history; the displacement along d and across it at the end of the walk."""
    H0 = X[hist]
    v_hist = float(((H0[-1, 1:] - H0[0, 1:]) @ d) / max(H0[-1, 0] - H0[0, 0], 1e-9))
    P = X[path[len(hist) - 1:]]
    s = P[:, 0] - P[0, 0]
    rel = P[:, 1:] - P[0, 1:]
    y = rel @ d
    k = s <= span
    v_walk = float(np.polyfit(s[k], y[k], 1)[0]) if k.sum() >= 3 and np.ptp(s[k]) > 0 else float("nan")
    return v_hist, v_walk, float(y[-1]), rel[-1] - y[-1] * d, float(s[-1])


def seed_task(p, seed):
    t_start = time.time()
    rng = np.random.default_rng(seed)
    X = sprinkle(p, rng)
    dim = int(p.get("dim", 3))                          # 3 (the test); 2 only for development rows
    order = (CapacityOrder3 if dim == 3 else CapacityOrder)(X, 0.0, float(p["width"]))
    m, span, K = int(p["memory"]), float(p["span"]), int(p["n_walkers"])
    hB = np.array([p["half_x"]] + [p["half_yz"]] * (dim - 1), float)
    hA = np.array([p["half_x_A"]] + [p["half_yz_A"]] * (dim - 1), float)
    n_st, half = int(p["n_starts"]), max(int(p["n_starts"]) // 2, 1)
    dT = float(p["dT_launch"])
    used, hists = set(), []
    for v in V_THROW:
        phase = rng.uniform(0.0, 2.0 * math.pi)
        x0 = v * (dT / 2.0 + span + 1.0) / 2.0         # history and walk symmetric in the box
        for j in range(n_st):
            sgn = 1.0 if j % 2 == 0 else -1.0
            d = np.zeros(dim)
            d[0] = sgn
            k = j // 2
            ang = phase + 2.0 * math.pi * k / half
            ring = [float(p["ring"]) * math.cos(ang), float(p["ring"]) * math.sin(ang)][:dim - 1] if dim == 3 \
                else [float(p["ring"]) * math.sin(ang)]
            tgt0 = np.array([float(p["t_launch"]) + float(p["t_spread"]) * k / half, -sgn * x0] + ring)
            near = np.argsort(((X - tgt0) ** 2).sum(axis=1))[:64]
            s_ = next(int(i) for i in near if int(i) not in used)
            used.add(s_)
            tgt = X[s_].copy()
            tgt[0] += dT
            tgt[1:] += v * dT * d
            cand = np.where(np.abs(X[:, 0] - tgt[0]) < 1.0)[0]
            q = int(cand[np.argmin(((X[cand] - tgt) ** 2).sum(axis=1))])
            ch = chain_between(order, s_, q, max(float(p["region_hist"]), v * dT + 2.0))
            if ch is None or len(ch) < 4:
                continue
            hist = ch[: len(ch) // 2 + 1]
            t0 = float(X[hist[-1], 0])
            cache_B, cache_A = {}, {}
            h = {"v": v, "j": j, "start": s_, "mem": min(m, len(hist) - 1), "walkers": []}
            for kw in range(K):
                key = int(id_hash(np.array([(seed % 4294967291) * 1000 + 7919 * kw + 13 * j + 100003 * int(v * 10)],
                                             dtype=np.uint64))[0])
                rec = []
                path = walker(order, hist, t0 + span + 1.0, m, hB, float(p["dt_B"]), key, cache_B, rec)
                vh, vw, y_end, perp, s_end = velocities(X, hist, path, d, span)
                wk = {"v_walk": vw, "y_end": y_end, "perp": [float(x) for x in perp], "s_end": s_end,
                      "stopped": float(s_end < span - 1.0)}
                for kk in STEP_KEYS:
                    vals = [r[kk] for r in rec if kk in r]
                    wk[kk] = float(np.mean(vals)) if vals else float("nan")
                if kw == 0:                              # M5: the same walker in region A
                    pathA = walker(order, hist, t0 + span + 1.0, m, hA, float(p["dt_A"]), key, cache_A, [])
                    h["identical"] = bool(pathA == path)
                    h["v_walk_A"] = velocities(X, hist, pathA, d, span)[1]
                h["walkers"].append(wk)
            h["v_hist"] = vh
            h["v_walk"] = float(np.nanmean([w["v_walk"] for w in h["walkers"]]))
            h["delta"] = h["v_walk"] - vh
            ys = np.array([w["y_end"] for w in h["walkers"]])
            pp = np.array([w["perp"] for w in h["walkers"]])
            h["spread_along"] = float(ys.std(ddof=1)) if ys.size > 1 else float("nan")
            h["spread_across"] = float(np.sqrt(((pp - pp.mean(axis=0)) ** 2).sum(axis=1).mean())) \
                if len(pp) > 1 else float("nan")
            proxy = walk_history(order, hist, t0 + span + 1.0, m, int(p["proxy_H"]), float(p["R_cyl"]))
            h["delta_proxy"] = velocities(X, hist, proxy, d, span)[1] - vh
            hists.append(h)
    return {"seed": int(seed), "N": int(len(X)), "hists": hists, "seconds": time.time() - t_start}


def sprinkle(p, rng):
    """Poisson events of density 1 in t [0, T], |x| <= half_x, |y|, |z| <= half_yz, sorted by time."""
    T, hx, hyz = float(p["T"]), float(p["half_x"]), float(p["half_yz"])
    dim = int(p.get("dim", 3))
    n = rng.poisson(T * (2 * hx) * (2 * hyz) ** (dim - 1))
    X = np.column_stack([rng.uniform(0, T, n), rng.uniform(-hx, hx, n)] +
                        [rng.uniform(-hyz, hyz, n) for _ in range(dim - 1)])
    return X[np.argsort(X[:, 0], kind="stable")]


# ============================================================================================== the test

def run(p, ctx):
    t0 = time.time()
    n, B, nse, tol = int(p["n_seeds"]), int(p["boot_B"]), float(p["n_se"]), float(p["tol"])
    jobs = [(p, derive_seed(ctx.seed, "F5.1f", ctx.part, k)) for k in range(n)]
    rows = parallel_map(seed_task, jobs, ctx.workers, ctx.progress,
                        lambda i, r_, dn, nn: f"[F5.1f] {dn}/{nn}: N {r_['N']}, {len(r_['hists'])} histories, "
                                              f"identical {np.mean([h['identical'] for h in r_['hists']]):.2f} "
                                              f"({time.time() - t0:.0f}s)")

    def bs(key, v):
        vals = [float(np.nanmean([h[key] for h in r_["hists"] if h["v"] == v])) for r_ in rows
                if any(h["v"] == v for h in r_["hists"])]
        return boot_se(vals, B, derive_seed(ctx.seed, "F5.1f-boot", ctx.part, key, v))[:2]

    allh = [h for r_ in rows for h in r_["hists"]]
    ident = float(np.mean([h["identical"] for h in allh])) if allh else float("nan")
    res = {}
    for v in V_THROW:
        H_ = [h for h in allh if h["v"] == v]
        W_ = [w for h in H_ for w in h["walkers"]]
        res[v] = {"n": len(H_), "nw": len(W_), "ident": float(np.mean([h["identical"] for h in H_])) if H_ else
                  float("nan"), "d": bs("delta", v), "dP": bs("delta_proxy", v), "v_hist": bs("v_hist", v),
                  "v_walk": bs("v_walk", v), "sd_walk": float(np.nanstd([w["v_walk"] for w in W_], ddof=1))
                  if len(W_) > 1 else float("nan"),
                  "along": float(np.nanmean([h["spread_along"] for h in H_])) if H_ else float("nan"),
                  "across": float(np.nanmean([h["spread_across"] for h in H_])) if H_ else float("nan"),
                  "mem": float(np.mean([h["mem"] for h in H_])) if H_ else float("nan")}
        for kk in STEP_KEYS + ("stopped",):
            res[v][kk] = float(np.nanmean([w[kk] for w in W_])) if W_ else float("nan")
    lines = [f"flat 3+1 order, density 1, box t [0, {p['T']:g}], |x| <= {p['half_x']:g}, |y|, |z| <= {p['half_yz']:g} "
             f"(~{np.mean([r_['N'] for r_ in rows]):.0f} events); the flow is shared among the children c of v with "
             f"min L(p, c) (memory {p['memory']}); {p['n_walkers']} walkers per history; walk {p['span']:g} + 1",
             f"  finite (M5): walker 0 has the same path in region A as in region B in {ident:.3f} of {len(allh)} "
             f"histories (v 0: {res[0.0]['ident']:.3f}, v 0.3: {res[0.3]['ident']:.3f}); criterion >= {p['ident']:g}"]
    for v in V_THROW:
        r_ = res[v]
        lines.append(f"  v {v:.1f} ({r_['n']} histories, {r_['nw']} walkers, memory {r_['mem']:.1f}): v_hist "
                     f"{r_['v_hist'][0]:+.3f}, centre {r_['v_walk'][0]:+.3f}; Delta_h {r_['d'][0]:+.3f} +- "
                     f"{r_['d'][1]:.3f}; proxy {r_['dP'][0]:+.3f} +- {r_['dP'][1]:.3f}")
        lines.append(f"          per step: admissible set {r_['size']:.1f} (deficit {r_['min_def']:+.2f}), children "
                     f"{r_['n_kids']:.1f}, inside {r_['n_adm']:.1f}, window {r_['n_win']:.0f}; per-walker SD of "
                     f"v_walk {r_['sd_walk']:.3f}; spread at the end along {r_['along']:.2f}, across {r_['across']:.2f};"
                     f" stopped {r_['stopped']:.2f}")
    dP = res[0.3]["dP"]
    instrument_ok = np.isfinite(dP[1]) and dP[0] + nse * dP[1] < -tol
    lines.append(f"instrument: the proxy brakes (Delta_h + 2 SE < -{tol:g}): {'OK' if instrument_ok else 'FAILS'}")
    finite = ident >= float(p["ident"])
    d0, d3 = res[0.0]["d"], res[0.3]["d"]
    kept = abs(d0[0]) + nse * d0[1] < tol and abs(d3[0]) + nse * d3[1] < tol
    lost = d3[0] + nse * d3[1] < -tol or abs(d0[0]) - nse * d0[1] > tol
    lines.append(f"  finite: {finite}; kept: {kept}; lost: {lost}")
    if not instrument_ok:
        status = "INSTRUMENT_FAILED"
    elif not finite or lost:
        status = "FAILED"
    elif kept:
        status = "PASSED"
    else:
        status = "UNDECIDED"
    return Result(status, lines, {"res": {str(k): v for k, v in res.items()}, "identical": ident, "rows": rows})
