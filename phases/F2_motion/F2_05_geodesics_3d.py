"""
F2.5 -- Geodesic motion in 3+1 and its dependence on N (weighted representation)
================================================================================

Where this sits (phase F2, gate of the roadmap): geodesic behaviour is
verified in 2+1; the predictive value of the action (axiom 6, W = sum Pi)
was shown only in the WEIGHTED representation in 2+1 (F2.4, E46b part B:
SPECIFIC, replicated -- the heaviest chains of sum Pi^2 read 0.0770 against
the geodesic of Pi^4 eta, 0.0776). Open: 3+1, and the dependence on N.
In the CONFORMAL representation the bending was universal (every maximal
chain follows the density, F2.2/F2.3, E46b part A), so it cannot carry the
specificity question; this test uses the weighted representation only.

Question
--------
d4 (the decision): in 3+1, do the heaviest chains of sum Pi (axiom 6)
    follow the geodesic of Pi^2 eta, while the competing functional sum Pi^2
    follows the geodesic of Pi^4 eta, clearly apart? I.e. does the form of
    the action fix the path as the continuum predicts, as in 2+1?
d3 (the 10 % question): in 2+1 the heaviest chains of sum Pi read about
    10 % ABOVE the continuum geodesic, those of sum Pi^2 on it (unexplained,
    docs/findings.md). Is that excess a finite-size effect that shrinks with
    N, or does it persist? Three sizes, many seeds.

Setup (the E46/E46b field and reading, unchanged)
-------------------------------------------------
Field Pi_m(x) = 1 - eps / sqrt(|x - m|^2 + soft^2), eps = 0.1, soft = 0.2,
source at m = -b e or +b e, b = 0.5, e the first spatial axis. End events
p = (0, 0) and q = (T, 0), T = 1. N0 points uniform in the causal interval
I(p, q); weights w = Pi_m (H1) or Pi_m^2 (H2) on the FLAT order (all
points). Reading: the mid-time displacement "away from the source", the
uniform average over all heaviest chains (core.chains.chains.
heaviest_between), mirror-antisymmetrised with common points:
delta = (d(source at -b e) + d(source at +b e)) / 2.

Predictions (continuum, core.chains.geodesics; independent of D because
the geodesic stays on the source axis): P1 = 0.0432 (geodesic of Pi^2 eta,
H1), P2 = 0.0776 (geodesic of Pi^4 eta, H2).

Sizes. 3+1: 32 000 and 64 000. The geometry of a 3+1 sprinkling is
readable from the order only from N* = 32 000 (K.1, history E54,
replicated; the meter itself is not ported), so 3+1 statements start
there. 2+1: 8 000, 16 000, 32 000.

Criteria (stated before the run)
--------------------------------
Statistics over seeds (the unit of independence); SE from B bootstrap
resamples of whole seeds, the same resample for H1 and H2.

d4, at the largest N, the rule of E46b part B (rel = rel2 = 0.25, n_se = 2;
Delta = rel * P1 = 0.0108):
  SPECIFIC   H1 in P1 +- Delta, mean(H2) - 2 SE(H2) > P1 + Delta, and
             H2 in P2 +- rel2 * P2
  UNIVERSAL  |H2 - H1| + 2 SE(H2 - H1) <= Delta
  otherwise  NO_STATEMENT
and the E46 condition that the error of H1 does not grow with N:
  |H1(N_big) - P1| <= |H1(N_small) - P1| + 2 sqrt(SE_big^2 + SE_small^2).
PASSED: SPECIFIC and the error does not grow. FAILED: UNIVERSAL, or
SPECIFIC with a growing error. UNDECIDED: NO_STATEMENT.

d3, the excess e(N) = mean(H1) - P1 (reported also relative to P1),
classified by a rule fixed here:
  FINITE_SIZE  e(8000) - e(32000) > 2 SE(difference)   (it shrinks)
  PERSISTENT   e(32000) - 2 SE > 0 and |e(8000) - e(32000)| <= 2 SE(difference)
  UNRESOLVED   otherwise
Status RECORDED (a diagnostic of a known discrepancy, not a test of the
theory). The E46b-B rule at 32 000 is reported as a replication of F2.4.

Instrument (any failing -> INSTRUMENT_FAILED)
---------------------------------------------
  (a) mirror identity at eps = 0 at the largest size: |d(P) + d(MP)| <=
      1e-12 for H1 and H2 on mirror_seeds seeds (M: x1 -> -x1), and the
      readings non-trivial (sd > 0) -- the K.7 / E46 gate at this size;
  (b) path weights sum to 1 within 1e-6 on every seed;
  (c) the 1+1 check of the heaviest-chain code on the E44b linear profile
      (Omega = 1 + 0.6 x, continuum 0.0584): within 25 % and > 3 SE;
  (d) seeds on which a reading is undefined <= 5 %.

Null model: no separate null; the competing functional H2 is the
alternative the decision is about (as in E46b part B).

Seen before this test was written: the E46b part B results in 2+1 at
N0 = 16 000 (SPECIFIC, replicated; H2 0.0770 vs P2 0.0776; H1 about 10 %
above P1). F2.5 uses new seed streams ("F2.5", "F2.5-mirror", "F2.5-boot").
"""

from __future__ import annotations

import math
import time

import numpy as np

from core.order.sprinkle import Interval, derive_seed
from core.chains.chains import heaviest_between
from core.chains.readings import with_ends
from core.chains.geodesics import radial_pi, radial_prediction, radial_prediction_k, geodesic_1p1
from core.stats.parallel import parallel_map
from core.testing import Result

TEST = {
    "id": "F2.5",
    "phase": "F2",
    "title": "Geodesic motion in 3+1 and the dependence on N (weighted representation)",
    "theory_part": "axiom 6 (action W = Σ Π) and motion in a field: the geodesic of the metric Π²η",
    "question": "Do the heaviest chains Σ Π in 3+1 follow the geodesic of Π²η and Σ Π² the geodesic of Π⁴η, clearly "
                "separated (the form of the action determines the path as in 2+1)? And is the ~10 % excess of Σ Π in "
                "2+1 an effect of finite N?",
    "method": "field of a soft point source (eps 0.1, soft 0.2, b 0.5, T 1); flat sprinkling with weights Π and Π²; "
              "average over all heaviest chains; mirror-antisymmetrised displacement at mid-time; 3+1 at N0 32 000 "
              "and 64 000, 2+1 at 8 000, 16 000, 32 000; bootstrap over seeds",
    "prediction": "H1 → P1 = 0.0432 (geodesic of Π²η), H2 → P2 = 0.0776 (geodesic of Π⁴η); the error of H1 does not "
                  "grow with N",
    "criterion": "INSTRUMENT: mirror identity at eps = 0 at the largest N (<= 1e-12, sd > 0), chain weights 1 ± "
                 "1e-6, 1+1 check (E44b), invalid seeds <= 5 %. d4: rule of E46b part B at the largest N (SPECIFIC / "
                 "UNIVERSAL / NO_STATEMENT, Delta = 0.25 P1, 2 SE) and the error of H1 does not grow with N. PASSED: "
                 "SPECIFIC and the error does not grow. FAILED: UNIVERSAL or a growing error. UNDECIDED: "
                 "NO_STATEMENT. d3: excess e(N) = H1 − P1 classified FINITE_SIZE / PERSISTENT / UNRESOLVED "
                 "(RECORDED).",
    "null_model": "the competing functional Σ Π² (H2) is the alternative being decided (as in E46b part B)",
    "level": "mechanism",
    "meters": ["K.7"],
    "depends_on": ["K.7", "F2.4"],
    "params": {"eps": 0.1, "soft": 0.2, "b": 0.5, "T": 1.0, "rel": 0.25, "rel2": 0.25, "n_se": 2.0, "boot_B": 4000,
               "invalid_max": 0.05, "mirror_tol": 1e-12, "weight_tol": 1e-6, "mirror_seeds": 4,
               "e44b_n": 20000, "e44b_seeds": 16},
    "parts": {"d4": {"D": 4, "sizes": [32000, 64000], "n_seeds": 256},
              "d3": {"D": 3, "sizes": [8000, 16000, 32000], "n_seeds": 768}},
    "quick": {"n_seeds": 6, "boot_B": 200, "mirror_seeds": 2, "e44b_n": 2000, "e44b_seeds": 4},
    "quick_parts": {"d4": {"sizes": [2000, 4000]}, "d3": {"sizes": [1500, 3000, 6000]}},
    "estimate": "d4: ~600 s per seed at 64 000 and ~150 s at 32 000 (4 heaviest chains per seed): ~50 h CPU, ~3–4 h "
                "with 16 processes; d3: ~40 h CPU, ~3 h with 16 processes",
}


def one_seed(D, N0, seed, eps, soft, b, T):
    """Mirror-antisymmetrised H1 and H2 readings of one sprinkling, and the worst path-weight error."""
    rng = np.random.default_rng(seed)
    P = Interval(D, T).propose(rng, N0)
    e = np.zeros(D - 1)
    e[0] = 1.0
    Pi = radial_pi(eps, soft)
    X = with_ends(P, D, T)
    q = len(X) - 1
    d = {"H1": {}, "H2": {}}
    wchk = []
    for sgn in (1, -1):
        m = -sgn * b * e                      # source on the -e side for sgn = +1: "away" = +e
        w = Pi(np.sqrt(((X[:, 1:] - m) ** 2).sum(axis=1)))
        h1 = heaviest_between(X, w, 0, q, sgn * e)
        h2 = heaviest_between(X, w * w, 0, q, sgn * e)
        d["H1"][sgn], d["H2"][sgn] = h1["disp"], h2["disp"]
        wchk += [h1.get("weight_check", 0.0), h2.get("weight_check", 0.0)]
    return {"H1": 0.5 * (d["H1"][1] + d["H1"][-1]), "H2": 0.5 * (d["H2"][1] + d["H2"][-1]),
            "weight_worst": max(abs(c - 1.0) for c in wchk)}


def mirror_task(D, N0, seed, T):
    """eps = 0, one-sided readings on P and on its mirror MP: |d(P) + d(MP)| for H1 and H2, and d(P)."""
    rng = np.random.default_rng(seed)
    P = Interval(D, T).propose(rng, N0)
    PM = P.copy()
    PM[:, 1] = -PM[:, 1]
    e = np.zeros(D - 1)
    e[0] = 1.0
    out = {}
    for name, power in (("H1", 1), ("H2", 2)):
        vals = []
        for pts in (P, PM):
            X = with_ends(pts, D, T)
            w = np.ones(len(X)) ** power             # eps = 0: Pi = 1, so w = 1 for both functionals
            vals.append(heaviest_between(X, w, 0, len(X) - 1, e)["disp"])
        out[name] = (abs(vals[0] + vals[1]), vals[0])
    return out


def e44b_check(seed_base, n, nseed):
    """1+1 linear profile Omega = 1 + 0.6 x between (0.05, 0) and (0.95, 0), weighted heaviest chains, mirrored."""
    T = 0.9
    om = lambda x: 1.0 + 0.6 * x
    pred = geodesic_1p1(lambda x: float(om(x)), T, 0.0, fracs=(0.5,))["disp"][0.5]
    ds = []
    for rep in range(nseed):
        rng = np.random.default_rng(derive_seed(seed_base, "F2.5-1p1", n, rep))
        X = with_ends(Interval(2, T).propose(rng, n), 2, T)
        dd = {}
        for sgn in (1, -1):
            dd[sgn] = heaviest_between(X, om(sgn * X[:, 1]), 0, len(X) - 1, [float(sgn)])["disp"]
        ds.append(0.5 * (dd[1] + dd[-1]))
    ds = np.array(ds)
    return {"prediction": pred, "mean": float(ds.mean()), "se": float(ds.std(ddof=1) / math.sqrt(len(ds)))}


def bootstrap(h1, h2, B, seed):
    """Means and bootstrap SE (whole seeds, same resample) of H1, H2 and H2 - H1."""
    n = h1.size
    idx = np.random.default_rng(seed).integers(0, n, size=(B, n))
    st = {}
    for k, v in (("H1", h1), ("H2", h2), ("H2-H1", h2 - h1)):
        bm = v[idx].mean(axis=1)
        st[k] = {"mean": float(v.mean()), "se": float(bm.std(ddof=1))}
    return st, idx


def decide_B(st, P1, P2, rel, rel2, nse):
    """The E46b part B rule."""
    Delta = rel * P1
    h1_band = abs(st["H1"]["mean"] - P1) <= Delta
    h2_above = st["H2"]["mean"] - nse * st["H2"]["se"] > P1 + Delta
    h2_p2 = abs(st["H2"]["mean"] - P2) <= rel2 * P2
    universal = abs(st["H2-H1"]["mean"]) + nse * st["H2-H1"]["se"] <= Delta
    v = "SPECIFIC" if (h1_band and h2_above and h2_p2) else "UNIVERSAL" if universal else "NO_STATEMENT"
    return {"verdict": v, "Delta": Delta, "H1_in_band": h1_band, "H2_above_band": h2_above, "H2_matches_P2": h2_p2,
            "equivalent": universal}


def run(p, ctx):
    D = int(p["D"])
    sizes = [int(x) for x in p["sizes"]]
    n = int(p["n_seeds"])
    eps, soft, b, T = (float(p[k]) for k in ("eps", "soft", "b", "T"))
    rel, rel2, nse, B = float(p["rel"]), float(p["rel2"]), float(p["n_se"]), int(p["boot_B"])
    P1 = radial_prediction(eps, soft, b, T)["mid"]
    P2 = radial_prediction_k(eps, soft, b, T, 2.0)["mid"]
    t0 = time.time()
    per = {}
    wworst = 0.0
    n_invalid_max = 0
    for N0 in sizes:
        args = [(D, N0, derive_seed(ctx.seed, "F2.5", D, N0, k), eps, soft, b, T) for k in range(n)]
        rows = parallel_map(one_seed, args, ctx.workers, ctx.progress,
                            lambda i, r_, dn, nn: f"[F2.5] {D-1}+1 N0={N0}: seed {dn}/{nn} H1 {r_['H1']:+.4f} "
                                                  f"H2 {r_['H2']:+.4f} ({time.time()-t0:.0f}s)")
        wworst = max([wworst] + [r_["weight_worst"] for r_ in rows])
        valid = [r_ for r_ in rows if r_["H1"] == r_["H1"] and r_["H2"] == r_["H2"]]
        n_invalid_max = max(n_invalid_max, n - len(valid))
        h1 = np.array([r_["H1"] for r_ in valid])
        h2 = np.array([r_["H2"] for r_ in valid])
        st, _ = bootstrap(h1, h2, B, derive_seed(ctx.seed, "F2.5-boot", D, N0))
        per[N0] = {"stats": st, "h1": h1, "h2": h2, "n_valid": len(valid), "decision": decide_B(st, P1, P2, rel, rel2, nse)}
    # instrument
    ms = [derive_seed(ctx.seed, "F2.5-mirror", D, sizes[-1], k) for k in range(int(p["mirror_seeds"]))]
    mres = parallel_map(mirror_task, [(D, sizes[-1], s, T) for s in ms], ctx.workers, ctx.progress,
                        lambda i, r_, dn, nn: f"[F2.5] {D-1}+1 mirror identity {dn}/{nn} ({time.time()-t0:.0f}s)")
    mworst = {k: max(r_[k][0] for r_ in mres) for k in ("H1", "H2")}
    msd = {k: float(np.std([r_[k][1] for r_ in mres], ddof=1)) if len(mres) > 1 else float("nan") for k in ("H1", "H2")}
    e44 = e44b_check(ctx.seed, int(p["e44b_n"]), int(p["e44b_seeds"]))
    e44_ok = e44["mean"] >= 3.0 * e44["se"] and abs(e44["mean"] - e44["prediction"]) <= max(rel * e44["prediction"], 2 * e44["se"])
    inst = {"mirror identity (eps = 0, largest N)": all(v <= float(p["mirror_tol"]) for v in mworst.values())
            and all(v > 0 for v in msd.values()),
            "path weights sum to 1": wworst <= float(p["weight_tol"]),
            "1+1 E44b check": e44_ok,
            "invalid seeds": n_invalid_max <= float(p["invalid_max"]) * n}
    big, small = per[sizes[-1]], per[sizes[0]]
    s_b, s_s = big["stats"]["H1"], small["stats"]["H1"]
    grow_ok = abs(s_b["mean"] - P1) <= abs(s_s["mean"] - P1) + 2.0 * math.sqrt(s_b["se"] ** 2 + s_s["se"] ** 2)
    # excess and its change with N (difference of independent sizes: SEs add in quadrature)
    exc = {N0: (per[N0]["stats"]["H1"]["mean"] - P1, per[N0]["stats"]["H1"]["se"]) for N0 in sizes}
    dex = exc[sizes[0]][0] - exc[sizes[-1]][0]
    sdex = math.sqrt(exc[sizes[0]][1] ** 2 + exc[sizes[-1]][1] ** 2)
    if dex > 2.0 * sdex:
        trend = "FINITE_SIZE"
    elif exc[sizes[-1]][0] - 2.0 * exc[sizes[-1]][1] > 0 and abs(dex) <= 2.0 * sdex:
        trend = "PERSISTENT"
    else:
        trend = "UNRESOLVED"
    dec = big["decision"]["verdict"]
    if not all(inst.values()):
        status = "INSTRUMENT_FAILED"
    elif D == 3:
        status = "RECORDED"
    elif dec == "SPECIFIC" and grow_ok:
        status = "PASSED"
    elif dec == "UNIVERSAL" or (dec == "SPECIFIC" and not grow_ok):
        status = "FAILED"
    else:
        status = "UNDECIDED"
    lines = [f"{D-1}+1 weighted representation: eps = {eps}, soft = {soft}, b = {b}, T = {T}; P1 = {P1:.4f} (Pi^2 eta), "
             f"P2 = {P2:.4f} (Pi^4 eta); Delta = {rel * P1:.4f}; {n} seeds per size"]
    for N0 in sizes:
        st = per[N0]["stats"]
        e_, se_ = exc[N0]
        lines.append(f"  N0={N0:<7d} H1 {st['H1']['mean']:+.4f} +- {st['H1']['se']:.4f}   H2 {st['H2']['mean']:+.4f} +- "
                     f"{st['H2']['se']:.4f}   H2-H1 {st['H2-H1']['mean']:+.4f} +- {st['H2-H1']['se']:.4f}   "
                     f"excess H1-P1 {e_:+.4f} +- {se_:.4f} ({e_ / P1:+.1%})   H2-P2 {st['H2']['mean'] - P2:+.4f}   "
                     f"rule B: {per[N0]['decision']['verdict']} (valid {per[N0]['n_valid']})")
    db = big["decision"]
    lines += [f"  largest N: H1 in P1 +- Delta {'yes' if db['H1_in_band'] else 'NO'}; H2 above the band {'yes' if db['H2_above_band'] else 'NO'}; "
              f"H2 matches P2 {'yes' if db['H2_matches_P2'] else 'NO'}; H2 = H1 within Delta {'yes' if db['equivalent'] else 'no'} -> {dec}",
              f"  error of H1 does not grow with N: {'yes' if grow_ok else 'NO'} (|H1-P1| {abs(s_s['mean'] - P1):.4f} at "
              f"{sizes[0]} -> {abs(s_b['mean'] - P1):.4f} at {sizes[-1]})",
              f"  excess trend: e({sizes[0]}) - e({sizes[-1]}) = {dex:+.4f} +- {sdex:.4f} -> {trend}",
              f"instrument: " + "; ".join(f"{k}: {'OK' if v else 'FAILS'}" for k, v in inst.items())
              + f" (mirror worst H1 {mworst['H1']:.1e}, H2 {mworst['H2']:.1e}; 1+1 {e44['mean']:.4f} +- {e44['se']:.4f} "
              f"vs {e44['prediction']:.4f}; weights within {wworst:.1e})"]
    data = {"P1": P1, "P2": P2, "sizes": sizes,
            "per_size": {str(N0): {"stats": v["stats"], "decision": v["decision"], "n_valid": v["n_valid"],
                                   "H1": v["h1"], "H2": v["h2"]} for N0, v in per.items()},
            "excess": {str(k): v for k, v in exc.items()}, "trend": trend, "error_not_growing": grow_ok,
            "instrument": inst, "mirror": {"worst": mworst, "sd": msd}, "e44b": e44}
    return Result(status, lines, data)
