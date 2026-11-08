"""
F5.1-diag3b: does the longest-chain rule keep its velocity without a cutting region? A more precise repetition of
F5.1-diag3 (development; decides nothing; the author's choice of 8 Nov)
================================================================================================================

Why
---
F5.1-diag3 (the engine without bitsets, 16 of 16 checks) showed that without the cylinder the thrown history keeps
most of its velocity (Binf - A0 = +0.24 +- 0.08 at v 0.3), but the frozen criterion |late - history| + 2 SE < 0.1 was
not met: the late slope is noisy (SD ~0.35 per walk), and the starts collapsed: the start was the event nearest to a
point with a jitter of 0.5, so in 22 of 64 walks (Binf, v 0.3) two starts of the same sprinkling were the same event and
gave the same walk. This repetition changes only the estimator, the starts and the arms; the rule, the engine, the box
and the order are those of F5.1-diag3.

What changes (fixed before the run)
-----------------------------------
  * starts: 16 per sprinkling, 8 for each direction of the throw (+-x), on a ring of radius 1.5 in (y, z) around
    x = -+2.5, with launch times spread over [1, 1.75]; every start is a different event (the nearest unused one) and the
    shared events between walks are recorded;
  * reference velocity: of the whole launch chain s -> q (6 units of time), v_ref = (x_q - x_s).d / (t_q - t_s), instead
    of its first half; reported also v_hist of F5.1-diag3;
  * primary estimator: the least-squares slope of the walk along the throw over the whole walk, s in [0, 12], v_walk;
    the deviation Delta = v_walk - v_ref per walk; mean per sprinkling, SE over sprinklings;
  * arms: A0 (the registered cylinder R = 4 at rest: instrument and reference) and Binf (the whole box, only diamonds
    inside it); Bplain and C2 are dropped (F5.1-diag3: Bplain = Binf within 1 SE; C2 measured the asymmetry of the
    guard band); v = 0 and 0.3; 8 sprinklings, new seeds.

Criteria (frozen before the run; SE over sprinklings)
-----------------------------------------------------
  INSTRUMENT_FAILED  if the estimator does not see the known braking in the cylinder: Delta(A0, 0.3) + 2 SE >= -0.1;
  KEPT               if |Delta(Binf, v)| + 2 SE < 0.1 for both v = 0 and v = 0.3;
  LOST               if Delta(Binf, 0.3) + 2 SE < -0.1 (the thrown history loses more than 0.1 without any cylinder);
  UNDECIDED          otherwise.
KEPT means that the loss of inertia of F5.1e / F5.1-diag2 was the computation (the cylinder), and the rule computed
in a large region is a candidate for F5.1f. LOST means that the rule itself loses inertia (a new hypothesis about the
rule). The diagnostic decides no registered test (R132).

Run from the root of the laboratory:
    python -m diag.inertia_box_b --seeds 8 --starts 16 --workers 8         (~3-4 h, ~2 GB per process)
    python -m diag.inertia_box_b --smoke
Each sprinkling is saved to diag/out/rows/boxb_*.json when done; --summarize reads them after a stopped run.
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
from diag.inertia_box import PARAMS as PARAMS3, SMOKE as SMOKE3, STEP_KEYS, walk_box

PARAMS = dict(PARAMS3, ring=1.5, t_spread=0.75)
SMOKE = dict(SMOKE3, ring=0.8, t_spread=0.5)
ARMS = (("A0", "cyl", 0.0), ("Binf", "adm", 0.0))
V_THROW = (0.0, 0.3)
LIMIT = 0.1                                            # the tolerance of the criteria


def walk_slope(X, path, n_hist, d, span):
    """Least-squares slope along d of the walk (from the end of the history) over s in [0, span]."""
    P = X[path[n_hist - 1:]]
    s = P[:, 0] - P[0, 0]
    y = (P[:, 1:] - P[0, 1:]) @ d
    k = s <= span
    if k.sum() < 3 or np.ptp(s[k]) <= 0:
        return float("nan")
    return float(np.polyfit(s[k], y[k], 1)[0])


def pick_starts(X, p, n_starts, rng):
    """n_starts distinct start events: half for each direction (+-x), on a ring in (y, z) around x = -+x_start, with
    launch times spread over [t_launch, t_launch + t_spread]; the nearest unused event to each point."""
    used, out = set(), []
    half = n_starts // 2
    phase = rng.uniform(0.0, 2.0 * math.pi)
    for j in range(n_starts):
        sgn = 1.0 if j % 2 == 0 else -1.0
        k = j // 2
        ang = phase + 2.0 * math.pi * k / max(half, 1)
        tgt0 = np.array([float(p["t_launch"]) + float(p["t_spread"]) * k / max(half, 1), -sgn * float(p["x_start"]),
                         float(p["ring"]) * math.cos(ang), float(p["ring"]) * math.sin(ang)])
        order_ = np.argsort(((X - tgt0) ** 2).sum(axis=1))
        s_ = next(int(i) for i in order_[:64] if int(i) not in used)
        used.add(s_)
        out.append((j, sgn, s_))
    return out


def seed_task(p, seed, n_starts):
    t_start = time.time()
    rng = np.random.default_rng(seed)
    X = sprinkle_flat_box(float(p["T"]), float(p["half_x"]), float(p["half_yz"]), rng)
    order = CapacityOrder3(X, 0.0, float(p["width"]))
    m, H, span = int(p["memory"]), int(p["horizon"]), float(p["span"])
    out = {"seed": int(seed), "N": int(len(X)), "walks": [], "checks": []}
    for (j, sgn, s_) in pick_starts(X, p, n_starts, rng):
        d = np.array([sgn, 0.0, 0.0])
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
            v_ref = float(((X[q, 1:] - X[s_, 1:]) @ d) / (X[q, 0] - X[s_, 0]))
            hist = ch[: len(ch) // 2 + 1]
            t0 = float(X[hist[-1], 0])
            for (name, mode, guard) in ARMS:
                rec = []
                path = walk_box(order, hist, t0 + span + 1.0, m, H, mode, guard, p, d, rec)
                if mode == "cyl" and len(out["checks"]) < 2:   # the engine must be the registered rule
                    ref = walk_history(order, hist, t0 + span + 1.0, m, H, float(p["R_cyl"]))
                    out["checks"].append(bool(ref == path))
                w = {"j": j, "v": v, "arm": f"{name}_v{v:g}", "name": name, "start": s_, "v_ref": v_ref,
                     "walk": [int(e) for e in path[len(hist):]]}
                w.update(measure(X, hist, path, d, span))
                w["v_walk"] = walk_slope(X, path, len(hist), d, span)
                w["delta"] = w["v_walk"] - v_ref
                w["delta_late"] = w["v_late"] - w["v_hist"]
                for k in STEP_KEYS:
                    vals = [r[k] for r in rec if k in r and np.isfinite(r[k])]
                    w[k] = float(np.mean(vals)) if vals else float("nan")
                out["walks"].append(w)
    # shared events: for every walk, the share of its events that another walk of the same arm and v also visits
    for w in out["walks"]:
        mine = set(w["walk"])
        others = set()
        for u in out["walks"]:
            if u is not w and u["arm"] == w["arm"]:
                others.update(u["walk"])
        w["shared"] = len(mine & others) / max(len(mine), 1)
    out["seconds"] = time.time() - t_start
    rd = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "rows")
    os.makedirs(rd, exist_ok=True)
    pre = "smoke_boxb" if p.get("smoke") else "boxb"   # smoke rows never enter --summarize
    with open(os.path.join(rd, f"{pre}_{seed}.json"), "w", encoding="utf-8") as f:
        json.dump(out, f)
    return out


# ============================================================================================== summary

def cluster(rows, key, arm, other=None):
    """Mean over sprinklings of the per-sprinkling mean of key in arm (or of the paired difference arm - other on the
    same start), SE over sprinklings, and the SE over walks (for comparison)."""
    per, allv = [], []
    for r in rows:
        Wd = {(w["j"], w["arm"]): w for w in r["walks"]}
        vals = []
        for (j, a), w in Wd.items():
            if a != arm or not np.isfinite(w[key]):
                continue
            if other is None:
                vals.append(w[key])
            elif (j, other) in Wd and np.isfinite(Wd[(j, other)][key]):
                vals.append(w[key] - Wd[(j, other)][key])
        if vals:
            per.append(float(np.mean(vals)))
            allv.extend(vals)
    per, allv = np.asarray(per, float), np.asarray(allv, float)
    if per.size == 0:
        return float("nan"), float("nan"), float("nan"), 0
    se = float(per.std(ddof=1) / math.sqrt(per.size)) if per.size > 1 else float("nan")
    sw = float(allv.std(ddof=1) / math.sqrt(allv.size)) if allv.size > 1 else float("nan")
    return float(per.mean()), se, sw, int(allv.size)


def verdict(rows):
    a0, s0, _, _ = cluster(rows, "delta", "A0_v0.3")
    if not (np.isfinite(s0) and a0 + 2 * s0 < -LIMIT):
        return "INSTRUMENT_FAILED"
    d0, e0, _, _ = cluster(rows, "delta", "Binf_v0")
    d3, e3, _, _ = cluster(rows, "delta", "Binf_v0.3")
    if abs(d0) + 2 * e0 < LIMIT and abs(d3) + 2 * e3 < LIMIT:
        return "KEPT"
    if d3 + 2 * e3 < -LIMIT:
        return "LOST"
    return "UNDECIDED"


def summary(p, rows):
    chk = [c for r in rows for c in r["checks"]]
    L = [f"F5.1-diag3b: inertia without a cutting region, flat 3+1, {len(rows)} sprinklings of "
         f"~{np.mean([r['N'] for r in rows]):.0f} events, rule of F5.1e (m {p['memory']}, H {p['horizon']}), "
         f"cap {p['cap']:g}, span {p['span']:g}",
         f"  the engine reproduces the registered rule (cylinder R {p['R_cyl']:g}): {all(chk) if chk else 'not checked'}"
         f" ({sum(chk)} of {len(chk)})",
         f"  distinct starts: {sum(len({w['start'] for w in r['walks']}) for r in rows)} "
         f"(of {sum(len({w['j'] for w in r['walks']}) for r in rows)} starts)",
         "  arm   v    walks  v_ref   v_hist  v_walk   Delta = walk - ref (SE spr., SE walks)   late - hist (SE)    "
         "shared  ties  du"]
    for (name, _, _) in ARMS:
        for v in V_THROW:
            k = f"{name}_v{v:g}"
            c = {x: cluster(rows, x, k) for x in ("v_ref", "v_hist", "v_walk", "delta", "delta_late", "shared", "n1",
                                                    "du")}
            if c["delta"][3] == 0:
                continue
            L.append(f"  {name:5s} {v:.1f} {c['delta'][3]:6d}  {c['v_ref'][0]:+.3f}  {c['v_hist'][0]:+.3f}  "
                     f"{c['v_walk'][0]:+.3f}   {c['delta'][0]:+.3f} +- {c['delta'][1]:.3f} ({c['delta'][2]:.3f})"
                     f"            {c['delta_late'][0]:+.3f} +- {c['delta_late'][1]:.3f}   {c['shared'][0]:.2f}   "
                     f"{c['n1'][0]:5.1f} {c['du'][0]:+.3f}")
    L.append("  paired Binf - A0 (Delta):")
    for v in V_THROW:
        m_, s_, sw, n_ = cluster(rows, "delta", f"Binf_v{v:g}", f"A0_v{v:g}")
        L.append(f"    v {v:.1f}: {m_:+.3f} +- {s_:.3f} (walks {sw:.3f}, {n_})")
    L.append(f"  criteria (frozen): instrument Delta(A0, 0.3) + 2 SE < -{LIMIT}; KEPT |Delta(Binf, v)| + 2 SE < {LIMIT} "
             f"for both v; LOST Delta(Binf, 0.3) + 2 SE < -{LIMIT}")
    L.append(f"  VERDICT: {verdict(rows)}")
    return L


# ============================================================================================== main

def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--seeds", type=int, default=8, help="flat sprinklings")
    ap.add_argument("--starts", type=int, default=16, help="distinct starts per sprinkling (even)")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--base-seed", type=int, default=0)
    ap.add_argument("--summarize", action="store_true", help="only summarise diag/out/rows/boxb_*.json")
    ap.add_argument("--smoke", action="store_true", help="a small box only to check that the code runs; its numbers "
                    "mean nothing")
    a = ap.parse_args()
    p = dict(PARAMS)
    if a.smoke:
        p.update(SMOKE, smoke=True)
    if a.summarize:
        rd = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "rows")
        rows = []
        for fn in sorted(glob.glob(os.path.join(rd, "boxb_*.json"))):
            with open(fn, encoding="utf-8") as f:
                rows.append(json.load(f))
        print("\n".join(summary(p, rows)))
        return
    jobs = [(p, derive_seed(a.base_seed, "diag-inertia-box-b", "smoke" if a.smoke else "main", k), a.starts)
            for k in range(a.seeds)]
    t0 = time.time()
    rows = parallel_map(seed_task, jobs, a.workers, print,
                        lambda i, r, dn, n: f"[diag] {dn}/{n}: N {r['N']}, {len(r['walks'])} walks "
                                            f"({time.time() - t0:.0f}s)")
    lines = summary(p, rows)
    print("\n".join(lines))
    od = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
    os.makedirs(od, exist_ok=True)
    fn = os.path.join(od, f"inertia_box_b{'_smoke' if a.smoke else ''}_{time.strftime('%Y%m%dT%H%M%S')}.json")
    with open(fn, "w", encoding="utf-8") as f:
        json.dump({"params": p, "args": vars(a), "summary": lines, "rows": rows}, f)
    print(f"saved: {fn}  ({time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main()
