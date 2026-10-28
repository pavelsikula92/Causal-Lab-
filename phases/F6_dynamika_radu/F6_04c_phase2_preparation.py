"""
F6.4c -- Does the crystalline phase of 2d orders (beta = 3.5) depend on how the chain was prepared?
====================================================================================================

Why this test exists (27 Oct 2026, author's decision after F6.4a)
-----------------------------------------------------------------
F6.4a (calibration against Surya 2012, N = 50, eps = 0.12) FAILED on one number:
phase I and the transition reproduce the paper to the second decimal, phase II
(beta = 3.5) has the action and the ordering fraction within tolerance, but the
height h (links of the longest chain) is 3.70 against the published 4.18
(tolerance 0.3). Phase II accepts 2 % of the moves and its ordering fraction is
not monotone in beta (0.593 / 0.576 / 0.605 at beta = 3.5 / 4.0 / 4.5): signs of
frozen states. New hypothesis after that failure (R132: F6.4a keeps its verdict
and its tolerances; this is a separate test with its own question): the height in
phase II depends on the preparation of the chain (thermalisation length, starting
state, cooling), so a single value from one preparation is not an equilibrium
quantity.

Model (locked before any run; the instrument of F6.4a, imported unchanged)
-------------------------------------------------------------------------
2d orders, N = 50, eps = 0.12, beta = 3.5, Metropolis-Hastings with exchanges of
two elements in U or V (one sweep = N(N-1)/2 moves), action and height from
phases/F6_dynamika_radu/F6_04a_bd_2d_orders.py. Five preparations, 8 independent
chains each (own seeds), then 3000 measured sweeps at beta = 3.5:
  P1 replicate: random start, 500 sweeps at beta = 3.5 (exactly F6.4a);
  P2 long:      random start, 5000 sweeps at beta = 3.5;
  P3 annealed:  random start, beta raised linearly from 2.0 to 3.5 over 5000 sweeps
                (slow cooling through the transition at ~2.4);
  P4 chain:     start from a total order (u = v: one chain of N elements, h = N - 1),
                5000 sweeps at beta = 3.5;
  P5 antichain: start from no relations at all (v = reversed u, h = 0), 5000 sweeps.
P2-P5 are the long preparations (ten times the thermalisation of F6.4a), coming
from the continuum side (P2, P3) and from the two extreme orders (P4, P5).
Observables per chain: mean h, S/hbar, ordering fraction chi over the measurement,
and the mean h of its first and second half (drift).

Criteria (stated before any run)
--------------------------------
Instrument:
  (i)  replicate: P1 reproduces F6.4a at beta = 3.5 (h 3.699, S -41.074):
       |h - 3.699| <= 0.3 and |S + 41.074| <= 1.0;
  (ii) the exact detailed-balance check of F6.4a (N = 4, eps = 0.5, beta = 1, TV <= 0.01).
  INSTRUMENT_FAILED if (i) or (ii) fails.
Preparation dependence (long preparations P2-P5; arm mean and SE over its 8 chains):
  independent  max over pairs |h_a - h_b| <= 0.15 (half the F6.4a tolerance) AND
               in every long arm |second half - first half| <= 0.15 (no drift);
  dependent    max over pairs (|h_a - h_b| - 2 SE of the difference) > 0.3.
Status:
  PASSED     independent and |h_eq - 4.18| <= 0.3 (h_eq = mean of P2-P5): phase II
             is an equilibrium phase and its published height is reproduced once the
             chain is mixed; the F6.4a mismatch was short thermalisation;
  FAILED     independent and |h_eq - 4.18| > 0.3: phase II is reproducible but its
             height differs from the paper (a difference of definition or of the
             paper's preparation; the laboratory's value stands);
  UNDECIDED  dependent (phase II is non-ergodic at this N and these run lengths: its
             height is not a calibration quantity), or neither independent nor
             dependent.
Reading declared now. Neither outcome changes F6.4a. F6.4b needs the continuum phase
and the transition, which F6.4a reproduced descriptively; whether that suffices for
F6.4b is the author's decision after this test, recorded as such.

Naming: the ID F6.4c (not "F6.4a2") because the registry accepts one letter after the number; F6.4b stays reserved
for matter in the action.

Quick run (after everything above was written; 2 chains, 50 / 150 preparation sweeps, 100 measured; decide nothing):
mechanics in order (replicate h 3.65, S -41.24 against F6.4a; exact balance TV 0.005); h of the two chains of P2 3.86
and 4.46, the other preparations 3.5-3.8; the short runs are "neither" (spread 0.53, drift up to 0.29). ~0.07 s per
sweep. No criterion changed.
"""

from __future__ import annotations

import itertools
import math
import time

import numpy as np

from core.order.sprinkle import derive_seed
from core.stats.parallel import parallel_map
from core.testing import Result
from phases.F6_dynamika_radu.F6_04a_bd_2d_orders import action, exact_check, f_weights, height_links, relations

TEST = {
    "id": "F6.4c",
    "phase": "F6",
    "title": "Závisí krystalická fáze 2D uspořádání (β = 3,5) na přípravě řetězce?",
    "theory_part": "fáze C, měřidlo pro F6.4b; nová hypotéza po F6.4a (výška fáze II 3,70 proti 4,18)",
    "question": "Sejde výška h ve fázi II (N = 50, ε = 0,12, β = 3,5) při desetkrát delší termalizaci ke stejné hodnotě "
                "z náhodného startu, pomalého ochlazení z β = 2 a z obou extrémních uspořádání (řetězec, antiřetězec)? "
                "Pokud ano, je to publikovaných 4,18?",
    "method": "měřidlo F6.4a beze změny; pět příprav po 8 řetězcích: P1 opakování F6.4a (500 sweepů), P2 náhodný start "
              "5000 sweepů, P3 ochlazení β 2,0 → 3,5 za 5000 sweepů, P4 start z řetězce, P5 start z antiřetězce (5000 "
              "sweepů); měření 3000 sweepů; h, S, χ a posun h mezi polovinami měření",
    "prediction": "neurčena předem (hypotéza: výška ve fázi II závisí na přípravě)",
    "criterion": "MĚŘIDLO: P1 neodpovídá F6.4a (|h − 3,699| > 0,3 nebo |S + 41,074| > 1) nebo přesná rovnováha N = 4 "
                 "(TV > 0,01). Nezávislé na přípravě: max |h_a − h_b| <= 0,15 mezi P2–P5 a posun polovin <= 0,15; "
                 "závislé: max (|h_a − h_b| − 2 SE) > 0,3. PROŠEL nezávislé a |h_eq − 4,18| <= 0,3; NEPROŠEL nezávislé "
                 "a |h_eq − 4,18| > 0,3; NEROZHODNUTO závislé nebo ani jedno.",
    "null_model": "fáze II je rovnovážná: všechny dlouhé přípravy dají stejnou výšku",
    "level": "calibration",
    "meters": ["F6.4a"],
    "depends_on": ["F6.4a"],
    "params": {"N": 50, "eps": 0.12, "beta": 3.5, "chains": 8, "meas": 3000, "short": 500, "long": 5000,
               "anneal_from": 2.0, "ref_paper_h": 4.18, "tol_paper": 0.3, "ref_f64a": {"h": 3.699, "S": -41.074},
               "tol_rep": {"h": 0.3, "S": 1.0}, "tol_indep": 0.15, "tol_dep": 0.3,
               "exact": {"N": 4, "eps": 0.5, "beta": 1.0, "moves": 300000, "tv_tol": 0.01}},
    "parts": {"main": {}},
    "quick": {"chains": 2, "meas": 100, "short": 50, "long": 150,
              "exact": {"N": 4, "eps": 0.5, "beta": 1.0, "moves": 60000, "tv_tol": 0.02}},
    "estimate": "~30–40 min s 10 procesy (40 řetězců, 35 500 sweepů celkem na řetězec × 8; ~0,07 s na sweep)",
}

PREPS = ("P1", "P2", "P3", "P4", "P5")
PREP_CS = {"P1": "opakování F6.4a (náhodný start, 500)", "P2": "náhodný start, dlouhá termalizace",
           "P3": "ochlazení z β = 2,0", "P4": "start z řetězce", "P5": "start z antiřetězce"}


def start_state(kind: str, N: int, rng: np.random.Generator):
    """Positions (u, v) of the elements in the two total orders for a starting state."""
    if kind == "chain":                                   # u = v: every pair related, one chain of N elements
        u = rng.permutation(N)
        return u, u.copy()
    if kind == "antichain":                               # v reversed against u: no pair related
        u = rng.permutation(N)
        return u, (N - 1) - u
    return rng.permutation(N), rng.permutation(N)


def prep_task(N, eps, beta, prep, n_prep, meas, anneal_from, seed):
    """One chain: a preparation of n_prep sweeps, then meas sweeps at beta; per-chain means and half-means of h."""
    rng = np.random.default_rng(seed)
    fw = f_weights(N, eps)
    kind = {"P4": "chain", "P5": "antichain"}.get(prep, "random")
    u, v = start_state(kind, N, rng)
    R = relations(u, v)
    S = action(R, eps, fw)
    moves = N * (N - 1) // 2
    hs, Ss, chis = [], [], []
    acc = 0
    for sweep in range(n_prep + meas):
        if prep == "P3" and sweep < n_prep:               # linear cooling through the transition
            b = anneal_from + (beta - anneal_from) * sweep / max(n_prep - 1, 1)
        else:
            b = beta
        for _ in range(moves):
            w = u if rng.random() < 0.5 else v
            i, j = rng.choice(N, 2, replace=False)
            w[i], w[j] = w[j], w[i]
            R2 = relations(u, v)
            S2 = action(R2, eps, fw)
            if S2 <= S or rng.random() < math.exp(-b * (S2 - S)):
                S, R = S2, R2
                if sweep >= n_prep:
                    acc += 1
            else:
                w[i], w[j] = w[j], w[i]
        if sweep >= n_prep:
            hs.append(height_links(u, v))
            Ss.append(S)
            chis.append(R.sum() / (N * (N - 1) / 2))
    hs = np.array(hs, float)
    half = len(hs) // 2
    return {"h": float(hs.mean()), "h1": float(hs[:half].mean()), "h2": float(hs[half:].mean()),
            "S": float(np.mean(Ss)), "chi": float(np.mean(chis)), "acc": acc / max(meas * moves, 1)}


def _dispatch(kind, args):
    return prep_task(*args) if kind == "chain" else exact_check(*args)


def run(p, ctx):
    t0 = time.time()
    N, eps, beta = int(p["N"]), float(p["eps"]), float(p["beta"])
    nc, meas = int(p["chains"]), int(p["meas"])
    n_prep = {"P1": int(p["short"]), "P2": int(p["long"]), "P3": int(p["long"]), "P4": int(p["long"]), "P5": int(p["long"])}
    jobs = [("chain", (N, eps, beta, pr, n_prep[pr], meas, float(p["anneal_from"]), derive_seed(ctx.seed, "F6.4c", pr, k)))
            for pr in PREPS for k in range(nc)]
    jobs.append(("exact", (p["exact"], derive_seed(ctx.seed, "F6.4c-exact"))))
    res = parallel_map(_dispatch, jobs, ctx.workers, ctx.progress,
                       lambda i, r_, d, n: f"[F6.4c] task {d}/{n} ({time.time() - t0:.0f}s)")
    exact = res[-1]
    arms = {}
    for a, pr in enumerate(PREPS):
        cs = res[a * nc:(a + 1) * nc]
        st = {}
        for k in ("h", "h1", "h2", "S", "chi", "acc"):
            vals = np.array([c[k] for c in cs])
            st[k] = (float(vals.mean()), float(vals.std(ddof=1) / math.sqrt(nc)) if nc > 1 else float("nan"))
        st["h_chains"] = [c["h"] for c in cs]
        arms[pr] = st

    lines = [f"2d orders, N = {N}, eps = {eps:g}, beta = {beta:g}; {nc} chains per preparation; {meas} measured sweeps; "
             f"h = links of the longest chain",
             "  prep  h (± SE)         h 1st/2nd half    S/hbar (± SE)       chi      acc     h of the chains"]
    for pr in PREPS:
        st = arms[pr]
        lines.append(f"  {pr}   {st['h'][0]:.3f} ± {st['h'][1]:.3f}   {st['h1'][0]:.3f} / {st['h2'][0]:.3f}    "
                     f"{st['S'][0]:+8.3f} ± {st['S'][1]:.3f}   {st['chi'][0]:.4f}  {st['acc'][0]:.4f}   "
                     + " ".join(f"{x:.2f}" for x in st["h_chains"]) + f"   ({PREP_CS[pr]})")

    # instrument
    ref, tr = p["ref_f64a"], p["tol_rep"]
    rep_ok = abs(arms["P1"]["h"][0] - ref["h"]) <= tr["h"] and abs(arms["P1"]["S"][0] - ref["S"]) <= tr["S"]
    ex_ok = exact["tv"] <= float(p["exact"]["tv_tol"])
    lines.append(f"  instrument (i) replicate of F6.4a: h {arms['P1']['h'][0]:.3f} vs {ref['h']} (±{tr['h']}), S "
                 f"{arms['P1']['S'][0]:+.3f} vs {ref['S']} (±{tr['S']}) -> {'OK' if rep_ok else 'FAILS'}")
    lines.append(f"  instrument (ii) exact balance N = {p['exact']['N']}: total variation {exact['tv']:.4f} "
                 f"(tolerance {p['exact']['tv_tol']}) -> {'OK' if ex_ok else 'FAILS'}")

    # preparation dependence among the long preparations
    long_ = ["P2", "P3", "P4", "P5"]
    ti, td = float(p["tol_indep"]), float(p["tol_dep"])
    pair_max, pair_sig = 0.0, -math.inf
    for a, b in itertools.combinations(long_, 2):
        d = abs(arms[a]["h"][0] - arms[b]["h"][0])
        se = math.hypot(arms[a]["h"][1], arms[b]["h"][1])
        pair_max = max(pair_max, d)
        pair_sig = max(pair_sig, d - 2.0 * se)
    drift = max(abs(arms[a]["h2"][0] - arms[a]["h1"][0]) for a in long_)
    independent = pair_max <= ti and drift <= ti
    dependent = pair_sig > td
    h_eq = float(np.mean([arms[a]["h"][0] for a in long_]))
    lines.append(f"  long preparations P2-P5: largest difference of h {pair_max:.3f} (independent if <= {ti}); largest "
                 f"(difference - 2 SE) {pair_sig:+.3f} (dependent if > {td}); largest drift between halves {drift:.3f} "
                 f"(<= {ti}) -> {'INDEPENDENT' if independent else 'DEPENDENT' if dependent else 'neither'}")
    lines.append(f"  h_eq (mean of P2-P5) = {h_eq:.3f}; published {p['ref_paper_h']} (tolerance {p['tol_paper']})")
    if not (rep_ok and ex_ok):
        status = "INSTRUMENT_FAILED"
    elif independent:
        status = "PASSED" if abs(h_eq - float(p["ref_paper_h"])) <= float(p["tol_paper"]) else "FAILED"
    else:
        status = "UNDECIDED"
    data = {"arms": {k: {kk: (list(vv) if isinstance(vv, tuple) else vv) for kk, vv in v.items()} for k, v in arms.items()},
            "pair_max": pair_max, "pair_sig": pair_sig, "drift": drift, "h_eq": h_eq, "independent": independent,
            "dependent": dependent, "exact": exact}
    return Result(status, lines, data)
