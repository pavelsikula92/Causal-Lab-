"""
Diagnostic of the longest-chain fall rule in three layers (development; decides nothing, not a registered test)
================================================================================================================

Why (1 Nov 2026, after F5.1d PASSED in 2+1 twice and F5.1e FAILED in 3+1)
--------------------------------------------------------------------------
In 3+1 the longest chains between fixed events read the geometry (R_max = 0.85), yet the history continued by the rule
of F5.1c/F5.1d does not accelerate; in both dimensions it receives a reproducible inward velocity jump at the apex. The
author's plan: separate a failure of reading from a failure of dynamics, and never change the longest-chain rule and the
rule of motion at the same time. Three layers:

  1. geometry   between fixed events at the same place: the tube of all maximal chains (events x with
                L(s, x) + L(x, q) = L(s, q)), its size, its outward sag at mid-time against the geometry and its width,
                at the scale at which the rule decides (m + H = 16 links) and at the scale of the control (dT);
  2. decision   at every step of the rule as it is: the size of the horizon set Q, the ties of min L(p, q) (level 1)
                and of the fewest events n(p, q) (level 2), how often the last tie-breaker (the earliest event, a
                coordinate-time rule) decides, the same for the step to a child; the radial velocity of the history
                (p -> v), of the chosen target (v -> q*) and of the mean of each tie set;
  3. dynamics   the same sprinklings and starts, only the last tie-breaker changed: "earliest" (the rule of F5.1d/e,
                reproduced exactly), "latest", and "random" (a pseudo-random order by the identity of the event, which
                does not depend on coordinates); the fall R = a / a_pred for each, loaded and flat; and in the flat order
                histories thrown with a velocity v (0 and 0.3), which a rule that respects relativity must keep.
                One more arm, "wide", belongs to layer 1: the registered rule unchanged, computed on a wider cylinder
                (x2 in 2+1, x1.5 in 3+1). The rule reads longest chains only inside a cylinder around the start (radius 4 in 3+1,
                while the window p -> q spans ~7 units of time); a truncated order shortens the chains near its edge.

The hypothesis under test (not a conclusion): the last tie-breaker prefers the event at rest in the laboratory, which
acts as a friction towards the frame of the source and can cancel the fall when ties are frequent (more in 3+1).

Run from the root of the laboratory:
    python -m diag.fall_layers --dim 2 --seeds 2 --n-ring 10 --workers 2
    python -m diag.fall_layers --dim 3 --seeds 8 --workers 8
Output: a summary on the screen and the raw records in diag/out/fall_layers_d{dim}_{stamp}.json.
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
from core.stats.stats import boot_se
from phases.F5_gravity.F5_01_free_fall_order import fall_curve, fit_a, sag_prediction
from phases.F5_gravity.F5_01c_free_fall_longest import LocalLinks, _popcount
from phases.F5_gravity.F5_01d_free_fall_apex import walk_history
from phases.F5_gravity.F5_01e_free_fall_apex_3d import (apex_history3, fall_curve3, ring_starts3, sag_prediction3,
                                                        same_place_target)
from phases.F7_light.F7_01a_capacity_cones import CapacityOrder
from phases.F7_light.F7_01b_throughput_cones import sprinkle_saturated
from phases.F7_light.F7_01c_throughput_cones_3d import CapacityOrder3, sprinkle_density

# the parameters of the registered tests (F5.1d for 2+1, F5.1e for 3+1), unchanged
PARAMS = {
    2: {"T": 46.0, "half": 18.0, "width": 8.0, "sigma0": 0.45, "r_launch": 8.0, "t_launch": 1.5, "dT_launch": 14.0,
        "span": 22.0, "memory": 8, "horizon": 8, "region": 12.0, "r0": 8.0, "t0": 7.0, "dT": 22.0, "wide": 2.0},
    3: {"T": 24.0, "half": 8.0, "width": 4.0, "sigma0": 0.3, "r_launch": 4.0, "t_launch": 1.0, "dT_launch": 6.0,
        "span": 12.0, "memory": 8, "horizon": 8, "region": 4.0, "r0": 4.0, "t0": 2.0, "dT": 16.0, "wide": 1.5},
}
VARIANTS = ("earliest", "latest", "random", "wide")
V_THROW = (0.0, 0.3)


# ============================================================================================== the order

def make_order(dim, sigma0, p, rng):
    """The throughput order of F7.1b (2+1) or F7.1c (3+1), exactly as in F5.1d / F5.1e."""
    if dim == 2:
        X = sprinkle_saturated(sigma0, float(p["width"]), float(p["T"]), float(p["half"]), rng)
        return CapacityOrder(X, sigma0, float(p["width"]))
    X = sprinkle_density(sigma0, float(p["width"]), float(p["T"]), float(p["half"]), rng)
    return CapacityOrder3(X, sigma0, float(p["width"]))


def starts(dim, X, r0, t0, n):
    """n events on a ring r = r0 at t0 (plane z = 0 in 3+1)."""
    if dim == 3:
        return ring_starts3(X, r0, t0, n)
    out = []
    for j in range(n):
        ang = 2.0 * math.pi * j / n
        target = np.array([t0, r0 * math.cos(ang), r0 * math.sin(ang)])
        out.append(int(np.argmin(((X - target) ** 2).sum(axis=1))))
    return out


def pred_curve(dim, p, r_apex, tau):
    """Fall from rest at r_apex in the geometry of the loaded order (phi = -sigma/3 in 2+1, -sigma/2 in 3+1)."""
    q = dict(p, r0=r_apex, t0=0.0)
    return fall_curve(q, tau) if dim == 2 else fall_curve3(q, tau)


def pred_sag(dim, p, r0, dT):
    q = dict(p, r0=r0, dT=dT)
    return sag_prediction(q) if dim == 2 else sag_prediction3(q)


def chain_between(order, s_, q, region):
    """The longest chain s_ -> q inside a cylinder (for histories); None if q is not reachable."""
    X = order.X
    ll = LocalLinks(order, X[s_, 1:], region, float(X[s_, 0]) - 1.0, float(X[q, 0]) + 1.0)
    js, jq = ll.pos[s_], ll.pos[q]
    if js < 0 or jq < 0:
        return None
    Lf, Lb = ll.forward(js), ll.backward(jq)
    if Lf[jq] < 2:
        return None
    path, j = [js], js
    while j != jq:                                     # step along the tube, earliest child (as longest_chain)
        kids = ll.c_idx[ll.c_ptr[j]:ll.c_ptr[j + 1]]
        kids = kids[(Lf[kids] == Lf[j] + 1) & (Lb[kids] == Lb[j] - 1)]
        if kids.size == 0:
            return None
        j = int(kids.min())
        path.append(j)
    return [int(ll.ids[k]) for k in path]


# ============================================================================================== layer 2 + 3: the walk

def _key(ll, js, tie, salt):
    """Order of candidates for the last tie-breaker: earliest = smallest local index (time order), latest = largest,
    random = a hash of the event's identity and a salt (independent of coordinates)."""
    js = np.asarray(js, dtype=np.int64)
    if tie == "earliest":
        return js
    if tie == "latest":
        return -js
    g = ll.ids[js].astype(np.uint64)
    h = (g * np.uint64(0x9E3779B97F4A7C15)) ^ np.uint64(salt & 0xFFFFFFFFFFFFFFFF)
    h ^= h >> np.uint64(29)
    h *= np.uint64(0xBF58476D1CE4E5B9)
    h ^= h >> np.uint64(32)
    return h


def walk_variant(order, hist, t_end, m, H, region, tie="earliest", salt=0, rec=None, rdir=None, wide=2.0):
    """walk_history of F5.1d with the last tie-breaker chosen by `tie` (identical to it for "earliest") and an optional
    record of every decision (layer 2). rdir: a fixed unit direction for the velocities (None: outward radial at v)."""
    X = order.X
    if tie == "wide":                                  # the registered rule on a wider cylinder (layer 1 check of the
        tie, region = "earliest", wide * region        # truncation of the order by the cylinder; x2 in 2+1, x1.5 in 3+1
                                                       # where x2 would need ~0.35 GB of bitsets per walk)
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
        Q = np.flatnonzero(Lv == H)
        if Q.size == 0:
            Hd = int(Lv.max())
            if Hd < 2:
                break
            Q = np.flatnonzero(Lv == Hd)
        Lp = ll.forward(jp)
        Q = Q[Lp[Q] >= 0]
        if Q.size == 0:
            break
        nQ = int(Q.size)
        sc = Lp[Q]
        Q1 = Q[sc == sc.min()]
        Q2 = Q1
        if Q1.size > 1:
            fp = np.packbits(Lp >= 0)
            nq = np.array([ll.between(fp, jp, int(q)) for q in Q1])
            Q2 = Q1[nq == nq.min()]
        jq = int(Q2[int(np.argmin(_key(ll, Q2, tie, salt)))])
        Lb = ll.backward(jq)
        kids = ll.c_idx[ll.c_ptr[jv]:ll.c_ptr[jv + 1]]
        kids = kids[Lb[kids] == Lb[jv] - 1]
        if kids.size == 0:
            break
        k1 = int(kids.size)
        if kids.size > 1:
            nk = np.array([_popcount(ll.C[jq] & np.packbits(ll.forward(int(k)) >= 0)) for k in kids])
            kids = kids[nk == nk.max()]
        k2 = int(kids.size)
        jk = int(kids[int(np.argmin(_key(ll, kids, tie, salt + 1)))])
        if rec is not None:
            xv = X[v]
            rh = rdir if rdir is not None else xv[1:] / max(np.linalg.norm(xv[1:]), 1e-12)

            def u(ev):                                 # radial (or rdir) velocity from v to the events ev
                d = X[np.atleast_1d(ev)] - xv
                return (d[:, 1:] @ rh) / np.maximum(d[:, 0], 1e-12)
            xp = X[src]
            u_hist = float(((xv[1:] - xp[1:]) @ rh) / max(xv[0] - xp[0], 1e-12))
            rec.append({"t": float(xv[0]), "nQ": nQ, "n1": int(Q1.size), "n2": int(Q2.size), "k1": k1, "k2": k2,
                        "u_hist": u_hist, "u_q": float(u(ll.ids[jq])[0]),
                        "u_Q1": float(u(ll.ids[Q1]).mean()), "u_Q2": float(u(ll.ids[Q2]).mean()),
                        "dt_q": float(X[ll.ids[jq], 0] - xv[0])})
        v = int(ll.ids[jk])
        path.append(v)
    return path


# ============================================================================================== layer 1: the tube

def tube_stats(order, s_, q, region, rhat):
    """The tube of maximal chains between s_ and q: L, number of tube events, near-maximal counts, and the outward
    offset (mean, sd) of the tube events at mid-time relative to the place of s_."""
    X = order.X
    ll = LocalLinks(order, X[s_, 1:], region, float(X[s_, 0]) - 1.0, float(X[q, 0]) + 1.0)
    js, jq = ll.pos[s_], ll.pos[q]
    if js < 0 or jq < 0:
        return None
    Lf, Lb = ll.forward(js), ll.backward(jq)
    L = int(Lf[jq])
    if L < 2:
        return None
    ok = (Lf >= 0) & (Lb >= 0)
    S = np.where(ok, Lf + Lb, -1)
    tube = np.flatnonzero(S == L)
    P = X[ll.ids[tube]]
    tm = 0.5 * (X[s_, 0] + X[q, 0])
    mid = np.abs(P[:, 0] - tm) < 1.0
    off = (P[mid, 1:] - X[s_, 1:]) @ rhat
    return {"L": L, "n_tube": int(tube.size), "n_m1": int((S == L - 1).sum()), "n_m2": int((S == L - 2).sum()),
            "mid_n": int(mid.sum()), "sag": float(off.mean()) if off.size else float("nan"),
            "width": float(off.std()) if off.size > 1 else float("nan")}


# ============================================================================================== one sprinkling

def seed_task(dim, p, arm, seed, n_ring, n_pairs, n_throw):
    rng = np.random.default_rng(seed)
    sigma0 = float(p["sigma0"]) if arm == "loaded" else 0.0
    order = make_order(dim, sigma0, p, rng)
    X = order.X
    m, H, reg, span = int(p["memory"]), int(p["horizon"]), float(p["region"]), float(p["span"])
    tau = np.linspace(0.0, span, 23)
    out = {"arm": arm, "seed": int(seed), "N": int(len(X)), "fall": {}, "steps": {}, "throw": {}, "tube": {}}
    # ---- layers 2 and 3: the apex histories, every variant on the same history
    hists = []
    for s_ in starts(dim, X, float(p["r_launch"]), float(p["t_launch"]), n_ring):
        h = apex_history3(order, s_, p)                 # works in any dimension (F5.1d apex_history)
        if h is not None:
            hists.append(h)
    out["n_hist"] = len(hists)
    checked = False
    dts = []
    for tie in VARIANTS:
        cols, preds, recs, vjump = [], [], [], []
        for i, hist in enumerate(hists):
            ap = hist[-1]
            ta, ra = float(X[ap, 0]), X[ap, 1:]
            rec = [] if tie == "earliest" else None
            path = walk_variant(order, hist, ta + span + 1.0, m, H, reg, tie, derive_seed(seed, "tie", i), rec,
                                wide=float(p["wide"]))
            if tie == "earliest" and not checked:     # the earliest variant must be the registered rule exactly
                ref = walk_history(order, hist, ta + span + 1.0, m, H, reg)
                out["reproduces_rule"] = bool(ref == path)
                checked = True
            P = X[path[len(hist) - 1:]]
            if len(P) > 1:
                dts.append(float(np.mean(np.diff(P[:, 0]))))
            rhat = ra / np.linalg.norm(ra)
            proj = (P[:, 1:] - ra) @ rhat
            y = np.interp(ta + tau, P[:, 0], proj)
            y[ta + tau > P[-1, 0]] = np.nan
            cols.append(y)
            preds.append(pred_curve(dim, p, float(np.linalg.norm(ra)), tau))
            H0 = X[hist]                                # velocity of the rising half before the apex
            k = max(0, len(hist) - 4)
            vin = float(((H0[-1, 1:] - H0[k, 1:]) @ rhat) / max(H0[-1, 0] - H0[k, 0], 1e-9))
            j = min(len(P) - 1, 3)                      # and of the first three steps after it
            vout = float(((P[j, 1:] - P[0, 1:]) @ rhat) / max(P[j, 0] - P[0, 0], 1e-9)) if j > 0 else float("nan")
            vjump.append((vin, vout))
            if rec is not None:
                for r_ in rec:
                    r_["chain"] = i
                    r_["tau"] = r_["t"] - ta
                recs.extend(rec)
        with np.errstate(all="ignore"):
            ybar = np.nanmean(cols, axis=0) if cols else np.full(len(tau), np.nan)
        ok = np.isfinite(ybar)
        a, v0 = fit_a(tau[ok], ybar[ok], 0.0) if ok.sum() >= 3 else (float("nan"), float("nan"))
        a_pred, _ = fit_a(tau, np.mean(preds, axis=0), 0.0) if preds else (float("nan"), 0.0)
        vj = np.array(vjump, float) if vjump else np.full((1, 2), np.nan)
        out["fall"][tie] = {"a": a, "v0": v0, "a_pred": a_pred, "series": ybar.tolist(),
                            "v_in": float(np.nanmean(vj[:, 0])), "v_out": float(np.nanmean(vj[:, 1]))}
        if recs:
            out["steps"][tie] = recs
    dt_link = float(np.mean(dts)) if dts else float("nan")
    out["dt_link"] = dt_link
    # ---- layer 3 (flat only): histories thrown with a velocity v, which a covariant rule must keep
    if arm == "flat":
        for vth in V_THROW:
            res = {t: [] for t in VARIANTS}
            for j in range(n_throw):
                ang = 2.0 * math.pi * (j + 0.5) / n_throw
                d = np.zeros(dim)
                d[0], d[1] = math.cos(ang), math.sin(ang)
                st = np.zeros(dim)                      # an independent start for every throw: a ring of radius 2
                st[0], st[1] = 2.0 * math.cos(ang + 1.0), 2.0 * math.sin(ang + 1.0)
                s_ = int(np.argmin(((X - np.r_[float(p["t_launch"]), st]) ** 2).sum(axis=1)))
                dT = float(p["dT_launch"])
                tgt = X[s_].copy()
                tgt[0] += dT
                tgt[1:] += vth * dT * d
                cand = np.where(np.abs(X[:, 0] - tgt[0]) < 1.0)[0]
                q = int(cand[np.argmin(((X[cand] - tgt) ** 2).sum(axis=1))])
                ch = chain_between(order, s_, q, reg)
                if ch is None or len(ch) < 3:
                    continue
                hist = ch[: len(ch) // 2 + 1]
                H0 = X[hist]
                vin = float(((H0[-1, 1:] - H0[0, 1:]) @ d) / max(H0[-1, 0] - H0[0, 0], 1e-9))
                t0 = float(H0[-1, 0])
                for tie in VARIANTS:
                    path = walk_variant(order, hist, t0 + span + 1.0, m, H, reg, tie, derive_seed(seed, "thr", j),
                                        wide=float(p["wide"]))
                    P = X[path[len(hist) - 1:]]
                    if len(P) < 3:
                        continue
                    s = P[:, 0] - P[0, 0]
                    y = (P[:, 1:] - P[0, 1:]) @ d
                    c = np.polyfit(s, y, 2)                 # y = a/2 s^2 + v s + y0
                    res[tie].append({"v_in": vin, "v_out": float(c[1]), "a": float(2 * c[0])})
            out["throw"][f"{vth:g}"] = res
    # ---- layer 1: tubes at the scale of the rule (16 links) and of the control (dT)
    for name, dT in (("rule", 2 * H * dt_link if np.isfinite(dt_link) else float("nan")), ("control", float(p["dT"]))):
        if not np.isfinite(dT):
            continue
        rows = []
        for s_ in starts(dim, X, float(p["r0"]), float(p["t0"]), n_pairs):
            q = same_place_target(order, s_, dT)
            if q is None:
                continue
            rh = X[s_, 1:] / np.linalg.norm(X[s_, 1:])
            st = tube_stats(order, s_, q, reg, rh)
            if st is not None:
                rows.append(st)
        out["tube"][name] = {"dT": dT, "pred": pred_sag(dim, p, float(p["r0"]), dT), "rows": rows}
    return out


# ============================================================================================== summary

def _ms(x):
    x = np.asarray([v for v in x if np.isfinite(v)], float)
    if x.size == 0:
        return float("nan"), float("nan")
    return float(x.mean()), (float(x.std(ddof=1) / math.sqrt(x.size)) if x.size > 1 else float("nan"))


def summary(dim, p, rows):
    L = [f"Fall rule diagnostic, {dim}+1, parameters of {'F5.1d' if dim == 2 else 'F5.1e'}; "
         f"{sum(r['arm'] == 'loaded' for r in rows)} loaded + {sum(r['arm'] == 'flat' for r in rows)} flat sprinklings"]
    L.append("  earliest variant reproduces the registered rule: "
             + str(all(r.get("reproduces_rule", True) for r in rows)))
    for arm in ("loaded", "flat"):
        R_ = [r for r in rows if r["arm"] == arm]
        if not R_:
            continue
        L.append(f"[{arm}]  events ~{np.mean([r['N'] for r in R_]):.0f}, histories {np.mean([r['n_hist'] for r in R_]):.1f}"
                 f" per sprinkling, time per link {np.mean([r['dt_link'] for r in R_]):.3f}")
        # layer 1
        for name in ("rule", "control"):
            tb = [r["tube"][name] for r in R_ if name in r["tube"]]
            if not tb:
                continue
            rr = [x for t in tb for x in t["rows"]]
            pred = float(np.mean([t["pred"] for t in tb]))
            sag, se = _ms([x["sag"] for x in rr])
            L.append(f"  L1 tube ({name}, dT {np.mean([t['dT'] for t in tb]):.2f}): L {np.mean([x['L'] for x in rr]):.1f}, "
                     f"tube events {np.mean([x['n_tube'] for x in rr]):.0f} (L-1: {np.mean([x['n_m1'] for x in rr]):.0f}, "
                     f"L-2: {np.mean([x['n_m2'] for x in rr]):.0f}), at mid-time {np.mean([x['mid_n'] for x in rr]):.1f}; "
                     f"sag {sag:+.3f} +- {se:.3f} against geometry {pred:+.3f}; width {np.nanmean([x['width'] for x in rr]):.3f}")
        # layer 2
        st = [s for r in R_ for s in r["steps"].get("earliest", [])]
        if st:
            A = {k: np.array([s[k] for s in st], float) for k in st[0] if k != "chain"}
            tie2 = A["n2"] > 1
            L.append(f"  L2 steps {len(st)}: |Q| {A['nQ'].mean():.1f}; level-1 ties (min L) {A['n1'].mean():.2f} "
                     f"(>1 in {np.mean(A['n1'] > 1):.2f}); level-2 ties (fewest events) >1 in {tie2.mean():.2f} -> the "
                     f"earliest decides; child ties {np.mean(A['k1'] > 1):.2f}, after events {np.mean(A['k2'] > 1):.2f}")
            m1, s1 = _ms(A["u_Q1"] - A["u_hist"])
            mq, sq = _ms(A["u_q"] - A["u_hist"])
            mb, sb = _ms((A["u_q"] - A["u_Q2"])[tie2])
            mr, sr = _ms(np.abs(A["u_q"][tie2]) - np.abs(A["u_Q2"][tie2]))
            L.append(f"      radial velocity, target minus history: chosen {mq:+.4f} +- {sq:.4f}; mean of the min-L set "
                     f"{m1:+.4f} +- {s1:.4f}")
            L.append(f"      last tie-breaker: chosen minus mean of the tie set {mb:+.4f} +- {sb:.4f}; |chosen| - |mean| "
                     f"{mr:+.4f} +- {sr:.4f} (negative = towards rest in the laboratory)")
        # layer 3: fall per variant
        for tie in VARIANTS:
            fr = [r["fall"][tie] for r in R_ if tie in r["fall"]]
            Rv = [f["a"] / f["a_pred"] for f in fr]
            Rm, Rs = _ms(Rv)
            L.append(f"  L3 {tie:8s}: R = a / a_pred {Rm:+.3f} +- {Rs:.3f}; v0 {np.mean([f['v0'] for f in fr]):+.4f}; "
                     f"velocity before / after the apex {np.mean([f['v_in'] for f in fr]):+.4f} / "
                     f"{np.mean([f['v_out'] for f in fr]):+.4f}")
        if arm == "flat":
            for vth in V_THROW:
                for tie in VARIANTS:
                    th = [x for r in R_ for x in r["throw"].get(f"{vth:g}", {}).get(tie, [])]
                    if not th:
                        continue
                    vi, _ = _ms([x["v_in"] for x in th])
                    vo, so = _ms([x["v_out"] for x in th])
                    aa, sa = _ms([x["a"] for x in th])
                    L.append(f"  L3 throw v {vth:.1f} {tie:8s}: history {vi:+.3f} -> walk {vo:+.3f} +- {so:.3f}, "
                             f"acceleration along {aa:+.4f} +- {sa:.4f} ({len(th)} walks)")
    return L


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--dim", type=int, choices=(2, 3), required=True)
    ap.add_argument("--seeds", type=int, default=2, help="sprinklings per arm (loaded and flat)")
    ap.add_argument("--n-ring", type=int, default=None, help="apex histories per sprinkling (default 40 / 60)")
    ap.add_argument("--n-pairs", type=int, default=12, help="pairs per scale for the tubes")
    ap.add_argument("--n-throw", type=int, default=6, help="thrown histories per velocity (flat)")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--base-seed", type=int, default=0)
    ap.add_argument("--smoke", action="store_true", help="a small box only to check that the code runs (3+1: the "
                    "quick parameters of F5.1e); its numbers mean nothing")
    a = ap.parse_args()
    p = dict(PARAMS[a.dim])
    if a.smoke and a.dim == 3:
        p.update({"half": 6.0, "T": 14.0, "span": 4.0, "r_launch": 3.0, "region": 3.0, "r0": 3.0, "dT": 8.0})
    n_ring = a.n_ring or (40 if a.dim == 2 else 60)
    jobs = [(a.dim, p, arm, derive_seed(a.base_seed, "diag-fall-layers", a.dim, arm, k), n_ring, a.n_pairs, a.n_throw)
            for k in range(a.seeds) for arm in ("loaded", "flat")]
    t0 = time.time()
    rows = parallel_map(seed_task, jobs, a.workers, print,
                        lambda i, r, d, n: f"[diag] {d}/{n}: {r['arm']} N {r['N']} ({time.time() - t0:.0f}s)")
    lines = summary(a.dim, p, rows)
    print("\n".join(lines))
    od = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
    os.makedirs(od, exist_ok=True)
    fn = os.path.join(od, f"fall_layers_d{a.dim}_{time.strftime('%Y%m%dT%H%M%S')}.json")
    with open(fn, "w", encoding="utf-8") as f:
        json.dump({"dim": a.dim, "params": p, "args": vars(a), "summary": lines, "rows": rows}, f)
    print(f"saved: {fn}  ({time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main()
