"""
F5.1-diag3c: does the longest-chain rule keep the velocity of the history it is given, without a cutting region?
(development; decides nothing; after F5.1-diag3b UNDECIDED, 9 Nov)
================================================================================================================

Why
---
F5.1-diag3b (engine 16 of 16, 128 distinct starts per arm) gave UNDECIDED by its frozen criteria: in the whole box
the thrown history deviates from the velocity of its launch chain by Delta = -0.099 +- 0.017 (the cylinder: -0.285 +-
0.011). After the run it showed that the reference was wrong, not the walk: the history handed to the rule (the first
half of the launch chain) moves at v_hist = +0.204 +- 0.017, while the whole chain moves at v_ref = +0.295 +- 0.003. The
launch chain is bent: its first half is slower than the whole chain. A check on a small flat order (40 chains, 9 Nov)
gave the same bend with the cylinder of chain_between at radius 4 and at radius 7 (-0.222 +- 0.037 in both), so the
bend is not the cut; it belongs to chain_between itself (it steps along the tube of maximal chains to the earliest
child, a coordinate-time choice). Relative to the history it was given, the walk keeps its velocity (v_walk - v_hist =
-0.009 +- 0.016), but that comparison was not frozen and the verdict of F5.1-diag3b stays (R132). The history at rest
also drifted towards the middle of the box (a = +0.018 +- 0.005), because its start lay at x = -+2.5.

What changes (fixed before the run)
-----------------------------------
  * the launch chain is built as before (chain_between, radius 4); recorded: the velocity of its first half (v_hist,
    the history handed to the rule), of its second half (v_hist2) and of the whole chain (v_ref);
  * the starts lie where the walk is symmetric in the box: x = -v (span + 1) / 2 along the throw (0 for v = 0, -1.95
    for v = 0.3), ring of radius 1.5 in (y, z), launch times spread over 0.75, 16 distinct starts per sprinkling;
  * the primary deviation is relative to the history handed to the rule: Delta_h = v_walk - v_hist, with v_walk the
    least-squares slope over the whole walk (s in [0, 12]) as in F5.1-diag3b;
  * the tolerance is 0.05 (one sixth of the throw), since the precision is now ~0.02;
  * arms A0 (cylinder R = 4 at rest, instrument) and Binf (whole box, only diamonds inside), v 0 and 0.3, 8 new
    sprinklings; rule, engine, box and order of F5.1-diag3.

Criteria (frozen before the run; SE over sprinklings)
-----------------------------------------------------
  INSTRUMENT_FAILED  if Delta_h(A0, 0.3) + 2 SE >= -0.05 (the estimator does not see the braking of the cylinder);
  KEPT               if |Delta_h(Binf, v)| + 2 SE < 0.05 for both v = 0 and v = 0.3;
  LOST               if Delta_h(Binf, 0.3) + 2 SE < -0.05;
  UNDECIDED          otherwise.
Reported, not evaluated: the straightness of the launch chain (v_hist - v_ref, v_hist2 - v_ref), Delta relative to
v_ref, the acceleration of a quadratic fit, the late slope. KEPT: the rule keeps the velocity it is given when no
region cuts its diamonds, and it is a candidate for F5.1f. LOST: the rule itself loses inertia.

Run from the root of the laboratory:
    python -m diag.inertia_box_c --seeds 8 --starts 16 --workers 8         (~3 h, ~2 GB per process)
    python -m diag.inertia_box_c --smoke
Each sprinkling is saved to diag/out/rows/boxc_*.json when done; --summarize reads them after a stopped run.
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
from phases.F5_gravity.F5_01d_free_fall_apex import walk_history
from phases.F7_light.F7_01c_throughput_cones_3d import CapacityOrder3
from diag.fall_layers import chain_between
from diag.inertia_frame import sprinkle_flat_box, measure
from diag.inertia_box import STEP_KEYS, walk_box
from diag.inertia_box_b import PARAMS as PARAMS_B, SMOKE as SMOKE_B, walk_slope, cluster

PARAMS = dict(PARAMS_B)
SMOKE = dict(SMOKE_B)
ARMS = (("A0", "cyl", 0.0), ("Binf", "adm", 0.0))
V_THROW = (0.0, 0.3)
LIMIT = 0.05                                           # the tolerance of the criteria


def pick_starts(X, p, n_starts, v, rng, used):
    """n_starts distinct start events for the throw velocity v: half for each direction (+-x), centred so that the walk
    is symmetric in the box (x = -v (span + 1) / 2 along the throw), on a ring in (y, z), launch times spread."""
    out = []
    half = n_starts // 2
    phase = rng.uniform(0.0, 2.0 * math.pi)
    x0 = v * (float(p["span"]) + 1.0) / 2.0
    for j in range(n_starts):
        sgn = 1.0 if j % 2 == 0 else -1.0
        k = j // 2
        ang = phase + 2.0 * math.pi * k / max(half, 1)
        tgt0 = np.array([float(p["t_launch"]) + float(p["t_spread"]) * k / max(half, 1), -sgn * x0,
                         float(p["ring"]) * math.cos(ang), float(p["ring"]) * math.sin(ang)])
        order_ = np.argsort(((X - tgt0) ** 2).sum(axis=1))
        s_ = next(int(i) for i in order_[:64] if int(i) not in used)
        used.add(s_)
        out.append((j, sgn, s_))
    return out


def _vel(X, a, b, d):
    """Velocity along d from event a to event b."""
    return float(((X[b, 1:] - X[a, 1:]) @ d) / max(X[b, 0] - X[a, 0], 1e-12))


def seed_task(p, seed, n_starts):
    t_start = time.time()
    rng = np.random.default_rng(seed)
    X = sprinkle_flat_box(float(p["T"]), float(p["half_x"]), float(p["half_yz"]), rng)
    order = CapacityOrder3(X, 0.0, float(p["width"]))
    m, H, span = int(p["memory"]), int(p["horizon"]), float(p["span"])
    out = {"seed": int(seed), "N": int(len(X)), "walks": [], "checks": []}
    used = set()
    for v in V_THROW:
        for (j, sgn, s_) in pick_starts(X, p, n_starts, v, rng, used):
            d = np.array([sgn, 0.0, 0.0])
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
            v_ref, v_h2 = _vel(X, s_, q, d), _vel(X, hist[-1], q, d)
            t0 = float(X[hist[-1], 0])
            for (name, mode, guard) in ARMS:
                rec = []
                path = walk_box(order, hist, t0 + span + 1.0, m, H, mode, guard, p, d, rec)
                if mode == "cyl" and len(out["checks"]) < 2:   # the engine must be the registered rule
                    ref = walk_history(order, hist, t0 + span + 1.0, m, H, float(p["R_cyl"]))
                    out["checks"].append(bool(ref == path))
                w = {"j": j, "v": v, "arm": f"{name}_v{v:g}", "name": name, "start": s_, "v_ref": v_ref,
                     "v_hist2": v_h2, "walk": [int(e) for e in path[len(hist):]]}
                w.update(measure(X, hist, path, d, span))
                w["v_walk"] = walk_slope(X, path, len(hist), d, span)
                w["delta_h"] = w["v_walk"] - w["v_hist"]
                w["delta_ref"] = w["v_walk"] - v_ref
                w["bend1"] = w["v_hist"] - v_ref
                w["bend2"] = v_h2 - v_ref
                for k in STEP_KEYS:
                    vals = [r[k] for r in rec if k in r and np.isfinite(r[k])]
                    w[k] = float(np.mean(vals)) if vals else float("nan")
                out["walks"].append(w)
    for w in out["walks"]:                             # the share of a walk's events visited by another walk
        mine = set(w["walk"])
        others = set()
        for u in out["walks"]:
            if u is not w and u["arm"] == w["arm"]:
                others.update(u["walk"])
        w["shared"] = len(mine & others) / max(len(mine), 1)
    out["seconds"] = time.time() - t_start
    rd = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "rows")
    os.makedirs(rd, exist_ok=True)
    pre = "smoke_boxc" if p.get("smoke") else "boxc"   # smoke rows never enter --summarize
    with open(os.path.join(rd, f"{pre}_{seed}.json"), "w", encoding="utf-8") as f:
        json.dump(out, f)
    return out


def verdict(rows):
    a0, s0, _, _ = cluster(rows, "delta_h", "A0_v0.3")
    if not (np.isfinite(s0) and a0 + 2 * s0 < -LIMIT):
        return "INSTRUMENT_FAILED"
    d0, e0, _, _ = cluster(rows, "delta_h", "Binf_v0")
    d3, e3, _, _ = cluster(rows, "delta_h", "Binf_v0.3")
    if abs(d0) + 2 * e0 < LIMIT and abs(d3) + 2 * e3 < LIMIT:
        return "KEPT"
    if d3 + 2 * e3 < -LIMIT:
        return "LOST"
    return "UNDECIDED"


def summary(p, rows):
    chk = [c for r in rows for c in r["checks"]]
    L = [f"F5.1-diag3c: inertia relative to the given history, flat 3+1, {len(rows)} sprinklings of "
         f"~{np.mean([r['N'] for r in rows]):.0f} events, rule of F5.1e (m {p['memory']}, H {p['horizon']}), "
         f"cap {p['cap']:g}, span {p['span']:g}, launch chain in a cylinder of radius {p['region_hist']:g}",
         f"  the engine reproduces the registered rule (cylinder R {p['R_cyl']:g}): {all(chk) if chk else 'not checked'}"
         f" ({sum(chk)} of {len(chk)})",
         f"  distinct starts: {sum(len({w['start'] for w in r['walks']}) for r in rows)} "
         f"(of {sum(len({(w['v'], w['j']) for w in r['walks']}) for r in rows)} starts)",
         "  arm   v    walks  v_ref   v_hist  v_hist2  v_walk   Delta_h = walk - hist (SE spr., SE walks)  "
         "walk - ref     a        shared"]
    for (name, _, _) in ARMS:
        for v in V_THROW:
            k = f"{name}_v{v:g}"
            c = {x: cluster(rows, x, k) for x in ("v_ref", "v_hist", "v_hist2", "v_walk", "delta_h", "delta_ref", "a",
                                                    "shared")}
            if c["delta_h"][3] == 0:
                continue
            L.append(f"  {name:5s} {v:.1f} {c['delta_h'][3]:6d}  {c['v_ref'][0]:+.3f}  {c['v_hist'][0]:+.3f}  "
                     f"{c['v_hist2'][0]:+.3f}   {c['v_walk'][0]:+.3f}   {c['delta_h'][0]:+.3f} +- {c['delta_h'][1]:.3f}"
                     f" ({c['delta_h'][2]:.3f})            {c['delta_ref'][0]:+.3f} +- {c['delta_ref'][1]:.3f}  "
                     f"{c['a'][0]:+.4f}  {c['shared'][0]:.2f}")
    for v in V_THROW:
        b1, e1, _, _ = cluster(rows, "bend1", f"Binf_v{v:g}")
        b2, e2, _, _ = cluster(rows, "bend2", f"Binf_v{v:g}")
        L.append(f"  launch chain v {v:.1f}: first half - whole {b1:+.3f} +- {e1:.3f}, second half - whole "
                 f"{b2:+.3f} +- {e2:.3f} (0 = straight)")
    L.append("  paired Binf - A0 (Delta_h):")
    for v in V_THROW:
        m_, s_, sw, n_ = cluster(rows, "delta_h", f"Binf_v{v:g}", f"A0_v{v:g}")
        L.append(f"    v {v:.1f}: {m_:+.3f} +- {s_:.3f} (walks {sw:.3f}, {n_})")
    L.append(f"  criteria (frozen): instrument Delta_h(A0, 0.3) + 2 SE < -{LIMIT}; KEPT |Delta_h(Binf, v)| + 2 SE < "
             f"{LIMIT} for both v; LOST Delta_h(Binf, 0.3) + 2 SE < -{LIMIT}")
    L.append(f"  VERDICT: {verdict(rows)}")
    return L


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--seeds", type=int, default=8, help="flat sprinklings")
    ap.add_argument("--starts", type=int, default=16, help="distinct starts per sprinkling and velocity (even)")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--base-seed", type=int, default=0)
    ap.add_argument("--summarize", action="store_true", help="only summarise diag/out/rows/boxc_*.json")
    ap.add_argument("--smoke", action="store_true", help="a small box only to check that the code runs; its numbers "
                    "mean nothing")
    a = ap.parse_args()
    p = dict(PARAMS)
    if a.smoke:
        p.update(SMOKE, smoke=True)
    if a.summarize:
        rd = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "rows")
        rows = []
        for fn in sorted(glob.glob(os.path.join(rd, "boxc_*.json"))):
            with open(fn, encoding="utf-8") as f:
                rows.append(json.load(f))
        print("\n".join(summary(p, rows)))
        return
    jobs = [(p, derive_seed(a.base_seed, "diag-inertia-box-c", "smoke" if a.smoke else "main", k), a.starts)
            for k in range(a.seeds)]
    t0 = time.time()
    rows = parallel_map(seed_task, jobs, a.workers, print,
                        lambda i, r, dn, n: f"[diag] {dn}/{n}: N {r['N']}, {len(r['walks'])} walks "
                                            f"({time.time() - t0:.0f}s)")
    lines = summary(p, rows)
    print("\n".join(lines))
    od = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
    os.makedirs(od, exist_ok=True)
    fn = os.path.join(od, f"inertia_box_c{'_smoke' if a.smoke else ''}_{time.strftime('%Y%m%dT%H%M%S')}.json")
    with open(fn, "w", encoding="utf-8") as f:
        json.dump({"params": p, "args": vars(a), "summary": lines, "rows": rows}, f)
    print(f"saved: {fn}  ({time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main()
