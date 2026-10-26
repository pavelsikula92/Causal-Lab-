"""
F6.4a -- Calibration: does the laboratory reproduce the phase structure of 2d orders weighted by the Benincasa-Dowker action?
=============================================================================================================================

Why this test exists (25 Oct 2026, author's decision)
-----------------------------------------------------
The F3 series on a fixed sprinkling has given what it can (F3.6c, F3.11): the
order is fixed and does not respond to matter, so the spatial part of the
metric (gamma) can only be put in. A local growth rule that yields manifold-
like orders is not known (variacni_audit.md, section 5). The route that is
known to work is an ensemble of orders weighted by a global action (phase C,
F6.4). Before any matter is coupled to such an ensemble (F6.4b: phi from the
clocks, psi from the cones, gamma = psi/phi as an output; it needs 3D, since a
2D order has no transverse direction), the laboratory must reproduce the
published 2D result with its own code. This is a calibration of the
instrument, not a test of the theory.

Reference (verified 25 Oct 2026 from the full text): S. Surya, "Evidence for
a continuum limit in causal set dynamics", Class. Quantum Grav. 29 (2012)
132001, arXiv:1110.6244.
  - 2d order: C = U intersect V of two total orders on N elements (e < f iff
    e precedes f in both).
  - Action (2D Benincasa-Dowker, smeared):
        S/hbar = 4 eps (N - 2 eps sum_{n=0}^{N-2} N_n f(n, eps)),
        f(n, eps) = (1 - eps)^n [1 - 2 eps n/(1 - eps) + eps^2 n (n - 1)/(2 (1 - eps)^2)],
    N_n = number of intervals with n elements strictly between (N_0 = links).
  - Weight exp(-beta S/hbar) (analytic continuation); Metropolis-Hastings; move:
    pick U or V at random and exchange the positions of two elements in it;
    one sweep = N(N - 1)/2 moves.
  - N = 50, eps = 0.12: transition "roughly 2.2 < beta_c < 2.5", located by the
    peaks of the fluctuations and of the autocorrelation times.
    Phase I (continuum), beta = 0.1: ordering fraction 0.499, d_MM ~ 2.004,
      height 10.232 +- 0.014, S/hbar = 3.846 +- 0.013.
    Phase II (crystalline, layered), beta = 3.5: ordering fraction 0.579,
      d_MM ~ 1.727, height 4.180 +- 0.018, S/hbar = -41.367 +- 0.054.

Convention fixed before the criteria (development calculation, 25 Oct, 20 000
uniform 2d orders, N = 50, eps = 0.12, no MCMC): S/hbar = 3.976 +- 0.009,
ordering fraction 0.5003, longest chain 11.31 ELEMENTS. The paper's height at
beta = 0.1 (10.23) therefore counts LINKS of the longest chain (elements - 1);
the laboratory reports it the same way. The action formula is consistent with
the paper (beta = 0.1 lowers <S> by about beta var(S) ~ 0.16).

Model and measurement (locked before any run)
---------------------------------------------
N = 50, eps = 0.12. beta grid 0, 0.1, 0.5, 1.0, 1.5, 1.8 ... 2.8 (step 0.1),
3.0, 3.5, 4.0, 4.5; four independent chains per beta (own seeds, random
start), 500 sweeps of thermalisation, 3000 sweeps of measurement, one sample
per sweep. Per chain: mean S/hbar, its variance (fluctuations), ordering
fraction chi = R / (N(N-1)/2), height h (links of the longest chain), number of
links N_0, acceptance, integrated autocorrelation time of S. Per beta: mean
over chains, SE = sd of the chain means / 2. d_MM from chi by inverting
chi = Gamma(d + 1) Gamma(d/2) / (2 Gamma(3 d/2)) (chi = 1/2 for d = 2).

Instrument checks (part of the run)
-----------------------------------
  (i)  exact detailed balance: N = 4, eps = 0.5, beta = 1, all 576 pairs (U, V)
       enumerated; the distribution of S from 300 000 MCMC moves against the exact
       Boltzmann distribution: total variation distance <= 0.01;
  (ii) beta = 0: one chain (100 + 1000 sweeps) against 20 000 directly sampled
       uniform 2d orders: |dS| <= 0.2, |d chi| <= 0.005, |dh| <= 0.15.
INSTRUMENT_FAILED if (i) or (ii) fails.

Criterion (stated before any run; a calibration: FAILED = the instrument is
not calibrated, F6.4b may not use it)
-----------------------------------------------------------------------
PASSED if all of:
  phase I, beta = 0.1:  |chi - 0.499| <= 0.01, |h - 10.232| <= 0.3, |S - 3.846| <= 0.3;
  phase II, beta = 3.5: |chi - 0.579| <= 0.02, |h - 4.180| <= 0.3, |S + 41.367| <= 2.0;
  transition: the beta of the largest fluctuation var(S) on the grid lies in
              [2.0, 2.7] (the published range 2.2-2.5 widened by 0.2, the grid step
              being 0.1).
FAILED otherwise. Tolerances: a few per cent, wider than the paper's quoted
errors, which do not include its autocorrelation or our finite chains; chosen
now, not after the run.
Reported, not decided: the whole scan (S, var S, chi, d_MM, h, links,
acceptance, tau_int) and the interval abundances N_0..N_10 at beta = 0.1 and 3.5.

Quick run (after everything above was written; 2 chains, 50 + 200 sweeps, beta 0.1 / 2.3 / 3.5; decide nothing):
instrument checks OK (exact balance TV 0.008 over 8 action values; beta = 0 chain against direct sampling within
tolerance); beta = 0.1: chi 0.496, h 10.15, S 3.79 (as published); beta = 2.3: var(S) 88 (the fluctuation peak);
beta = 3.5: S -41.5 and chi 0.593 as published but h 3.64 against 4.18, with acceptance 2 % (slow mixing in the
crystalline phase; the full run thermalises ten times longer). ~0.07 s per sweep and chain. No criterion changed.
"""

from __future__ import annotations

import bisect
import itertools
import math
import time

import numpy as np

from core.order.sprinkle import derive_seed
from core.stats.parallel import parallel_map
from core.testing import Result

TEST = {
    "id": "F6.4a",
    "phase": "F6",
    "title": "Kalibrace: reprodukuje laboratoř fázovou strukturu 2D uspořádání s Benincasovou–Dowkerovou akcí?",
    "theory_part": "fáze C (variacni_audit.md, oddíl 5): soubor uspořádání vážený akcí; měřidlo pro F6.4b (hmota v akci, "
                   "γ jako výstup)",
    "question": "Dá vlastní kód laboratoře pro 2D uspořádání (N = 50, ε = 0,12) vážená exp(−βS) publikovanou fázi podobnou "
                "varietě při β = 0,1, krystalickou fázi při β = 3,5 a přechod mezi β ≈ 2,2 a 2,5 (Surya 2012)?",
    "method": "Metropolis–Hastings přes 2D uspořádání (výměna dvou prvků v jednom ze dvou lineárních uspořádání), BD akce "
              "ve 2D s ε = 0,12; mřížka β od 0 do 4,5 (krok 0,1 kolem přechodu), 4 nezávislé řetězce, 500 + 3000 sweepů; "
              "kontroly: přesná rovnováha pro N = 4 a β = 0 proti přímému rovnoměrnému vzorkování",
    "prediction": "fáze I: χ 0,499, h 10,232, S 3,846; fáze II: χ 0,579, h 4,180, S −41,367; přechod 2,2 < β_c < 2,5",
    "criterion": "MĚŘIDLO: přesná rovnováha N = 4 (TV <= 0,01) nebo β = 0 proti přímému vzorkování (|ΔS| <= 0,2, "
                 "|Δχ| <= 0,005, |Δh| <= 0,15). PROŠEL: β = 0,1 |χ − 0,499| <= 0,01, |h − 10,232| <= 0,3, |S − 3,846| <= 0,3; "
                 "β = 3,5 |χ − 0,579| <= 0,02, |h − 4,180| <= 0,3, |S + 41,367| <= 2; maximum var(S) na mřížce v [2,0; 2,7]. "
                 "NEPROŠEL jinak (měřidlo nekalibrované).",
    "null_model": "rovnoměrné 2D uspořádání (β = 0)",
    "level": "calibration",
    "meters": [],
    "depends_on": [],
    "params": {"N": 50, "eps": 0.12,
               "betas": [0.0, 0.1, 0.5, 1.0, 1.5, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 3.0, 3.5, 4.0, 4.5],
               "chains": 4, "therm": 500, "meas": 3000,
               "ref": {"I": {"beta": 0.1, "chi": 0.499, "h": 10.232, "S": 3.846},
                       "II": {"beta": 3.5, "chi": 0.579, "h": 4.180, "S": -41.367}},
               "tol": {"I": {"chi": 0.01, "h": 0.3, "S": 0.3}, "II": {"chi": 0.02, "h": 0.3, "S": 2.0}},
               "beta_c_window": [2.0, 2.7],
               "exact": {"N": 4, "eps": 0.5, "beta": 1.0, "moves": 300000, "tv_tol": 0.01},
               "zero": {"therm": 100, "meas": 1000, "direct": 20000, "tol_S": 0.2, "tol_chi": 0.005, "tol_h": 0.15}},
    "parts": {"main": {}},
    "quick": {"betas": [0.1, 2.3, 3.5], "chains": 2, "therm": 50, "meas": 200,
              "exact": {"N": 4, "eps": 0.5, "beta": 1.0, "moves": 60000, "tv_tol": 0.02},
              "zero": {"therm": 20, "meas": 150, "direct": 2000, "tol_S": 0.5, "tol_chi": 0.01, "tol_h": 0.3}},
    "estimate": "~30–60 min s 10 procesy (80 řetězců × 3500 sweepů; jeden sweep 1225 tahů, ~0,06 s)",
}


# ============================================================================================== the order and its action

def f_weights(N: int, eps: float) -> np.ndarray:
    """f(n, eps) for n = 0 .. N-2 (module notes)."""
    n = np.arange(max(N - 1, 1), dtype=float)
    return (1.0 - eps) ** n * (1.0 - 2.0 * eps * n / (1.0 - eps) + eps ** 2 * n * (n - 1.0) / (2.0 * (1.0 - eps) ** 2))


def relations(u: np.ndarray, v: np.ndarray) -> np.ndarray:
    """R[e, f] = e < f in the 2d order: e precedes f in both total orders (u, v are the positions of the elements)."""
    return (u[:, None] < u[None, :]) & (v[:, None] < v[None, :])


def interval_abundances(R: np.ndarray) -> np.ndarray:
    """N_n: number of related pairs with exactly n elements strictly between them, n = 0 .. N-2."""
    N = R.shape[0]
    Rf = R.astype(np.float32)
    between = (Rf @ Rf)[R]                                     # for each related pair: #k with e < k < f
    return np.bincount(np.rint(between).astype(np.int64), minlength=max(N - 1, 1))[:max(N - 1, 1)]


def action(R: np.ndarray, eps: float, fw: np.ndarray) -> float:
    """2D Benincasa-Dowker smeared action S/hbar (module notes)."""
    N = R.shape[0]
    Nn = interval_abundances(R)
    return 4.0 * eps * (N - 2.0 * eps * float((Nn * fw).sum()))


def height_links(u: np.ndarray, v: np.ndarray) -> int:
    """Links of the longest chain (= its elements - 1): longest increasing run of v in the order of u."""
    seq = v[np.argsort(u)]
    tails = []
    for x in seq:
        k = bisect.bisect_left(tails, x)
        if k == len(tails):
            tails.append(x)
        else:
            tails[k] = x
    return len(tails) - 1


def chi_of_d(d: float) -> float:
    """Ordering fraction of a d-dimensional Minkowski interval (Myrheim-Meyer), chi = 1/2 for d = 2."""
    return math.exp(math.lgamma(d + 1) + math.lgamma(d / 2) - math.log(2.0) - math.lgamma(3 * d / 2))


def d_mm(chi: float) -> float:
    """Myrheim-Meyer dimension by bisection (chi decreases with d)."""
    if not (0.0 < chi < 1.0):
        return float("nan")
    lo, hi = 0.5, 12.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if chi_of_d(mid) > chi:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def tau_int(x: np.ndarray) -> float:
    """Integrated autocorrelation time (sum of the autocorrelation up to its first negative value)."""
    x = np.asarray(x, float) - np.mean(x)
    n = x.size
    var = float(x @ x) / n
    if var <= 0:
        return 0.5
    tau = 0.5
    for lag in range(1, n // 2):
        c = float(x[:-lag] @ x[lag:]) / (n - lag) / var
        if c <= 0:
            break
        tau += c
    return tau


# ============================================================================================== Markov chain

def chain_task(N, eps, beta, therm, meas, seed):
    """One chain from a random 2d order: means, fluctuations and autocorrelation of the observables."""
    rng = np.random.default_rng(seed)
    fw = f_weights(N, eps)
    u, v = rng.permutation(N), rng.permutation(N)
    R = relations(u, v)
    S = action(R, eps, fw)
    moves = N * (N - 1) // 2
    acc = 0
    Ss, chis, hs, links = [], [], [], []
    Nn_sum = np.zeros(max(N - 1, 1))
    for sweep in range(therm + meas):
        for _ in range(moves):
            w = u if rng.random() < 0.5 else v
            a, b = rng.choice(N, 2, replace=False)
            w[a], w[b] = w[b], w[a]                                # exchange the positions of two elements
            R2 = relations(u, v)
            S2 = action(R2, eps, fw)
            if S2 <= S or rng.random() < math.exp(-beta * (S2 - S)):
                S, R = S2, R2
                acc += 1
            else:
                w[a], w[b] = w[b], w[a]                            # rejected: undo
        if sweep >= therm:
            Nn = interval_abundances(R)
            Ss.append(S)
            chis.append(R.sum() / (N * (N - 1) / 2))
            hs.append(height_links(u, v))
            links.append(int(Nn[0]))
            Nn_sum += Nn
    Ss = np.array(Ss)
    return {"S": float(Ss.mean()), "varS": float(Ss.var()), "chi": float(np.mean(chis)), "h": float(np.mean(hs)),
            "links": float(np.mean(links)), "acc": acc / ((therm + meas) * moves), "tau": tau_int(Ss),
            "Nn": (Nn_sum / max(meas, 1))[:11].tolist()}


def exact_check(p, seed):
    """Instrument check (i): MCMC distribution of S for N = 4 against the exact Boltzmann distribution."""
    N, eps, beta, nmov = int(p["N"]), float(p["eps"]), float(p["beta"]), int(p["moves"])
    fw = f_weights(N, eps)
    perms = [np.array(x) for x in itertools.permutations(range(N))]
    vals = {}
    for u in perms:
        for v in perms:
            s = round(action(relations(u, v), eps, fw), 9)
            vals[s] = vals.get(s, 0.0) + math.exp(-beta * s)
    Z = sum(vals.values())
    exact = {k: w / Z for k, w in vals.items()}
    rng = np.random.default_rng(seed)
    u, v = rng.permutation(N), rng.permutation(N)
    S = action(relations(u, v), eps, fw)
    hist = {}
    for _ in range(nmov):
        w = u if rng.random() < 0.5 else v
        a, b = rng.choice(N, 2, replace=False)
        w[a], w[b] = w[b], w[a]
        S2 = action(relations(u, v), eps, fw)
        if S2 <= S or rng.random() < math.exp(-beta * (S2 - S)):
            S = S2
        else:
            w[a], w[b] = w[b], w[a]
        k = round(S, 9)
        hist[k] = hist.get(k, 0) + 1
    tv = 0.5 * sum(abs(hist.get(k, 0) / nmov - exact.get(k, 0.0)) for k in set(hist) | set(exact))
    return {"tv": tv, "n_values": len(exact)}


def direct_task(N, eps, n, seed):
    """Uniform 2d orders sampled directly (independent permutations): means of S, chi, h."""
    rng = np.random.default_rng(seed)
    fw = f_weights(N, eps)
    out = []
    for _ in range(n):
        u, v = rng.permutation(N), rng.permutation(N)
        R = relations(u, v)
        out.append((action(R, eps, fw), R.sum() / (N * (N - 1) / 2), height_links(u, v)))
    a = np.array(out)
    return {"S": float(a[:, 0].mean()), "chi": float(a[:, 1].mean()), "h": float(a[:, 2].mean()),
            "seS": float(a[:, 0].std() / math.sqrt(n))}


# ============================================================================================== run

def run(p, ctx):
    t0 = time.time()
    N, eps = int(p["N"]), float(p["eps"])
    betas = [float(b) for b in p["betas"]]
    nc = int(p["chains"])
    z = p["zero"]
    jobs = [("chain", (N, eps, b, int(p["therm"]), int(p["meas"]), derive_seed(ctx.seed, "F6.4a", b, k)))
            for b in betas for k in range(nc)]
    jobs.append(("chain", (N, eps, 0.0, int(z["therm"]), int(z["meas"]), derive_seed(ctx.seed, "F6.4a-zero"))))
    jobs.append(("direct", (N, eps, int(z["direct"]), derive_seed(ctx.seed, "F6.4a-direct"))))
    jobs.append(("exact", (p["exact"], derive_seed(ctx.seed, "F6.4a-exact"))))
    res = parallel_map(_dispatch, jobs, ctx.workers, ctx.progress,
                       lambda i, r_, d, n: f"[F6.4a] task {d}/{n} ({time.time() - t0:.0f}s)")
    chains = res[:len(betas) * nc]
    zero_chain, direct, exact = res[-3], res[-2], res[-1]

    lines = [f"2d orders, N = {N}, eps = {eps:g}; Metropolis with exchanges in U or V; {nc} chains per beta, "
             f"{p['therm']} + {p['meas']} sweeps (1 sweep = {N * (N - 1) // 2} moves); h = links of the longest chain",
             "    beta      S/hbar        var S     chi      d_MM     h       links   acc    tau_int"]
    scan = {}
    for i, b in enumerate(betas):
        cs = chains[i * nc:(i + 1) * nc]
        m = {k: float(np.mean([c[k] for c in cs])) for k in ("S", "varS", "chi", "h", "links", "acc", "tau")}
        se = {k: float(np.std([c[k] for c in cs], ddof=1) / math.sqrt(nc)) if nc > 1 else float("nan")
              for k in ("S", "chi", "h")}
        m["d_mm"] = d_mm(m["chi"])
        m["Nn"] = np.mean([c["Nn"] for c in cs], axis=0).tolist()
        scan[b] = {"mean": m, "se": se}
        lines.append(f"    {b:4.2f}  {m['S']:+9.3f} ±{se['S']:.3f}  {m['varS']:8.2f}  {m['chi']:.4f}  {m['d_mm']:.3f}  "
                     f"{m['h']:6.3f}  {m['links']:6.1f}  {m['acc']:.3f}  {m['tau']:7.1f}")

    # instrument
    inst_exact = exact["tv"] <= float(p["exact"]["tv_tol"])
    dS, dchi, dh = zero_chain["S"] - direct["S"], zero_chain["chi"] - direct["chi"], zero_chain["h"] - direct["h"]
    inst_zero = abs(dS) <= float(z["tol_S"]) and abs(dchi) <= float(z["tol_chi"]) and abs(dh) <= float(z["tol_h"])
    lines.append(f"  instrument (i) exact balance N = {p['exact']['N']}, beta = {p['exact']['beta']:g}: total variation "
                 f"{exact['tv']:.4f} over {exact['n_values']} action values (tolerance {p['exact']['tv_tol']}) -> "
                 f"{'OK' if inst_exact else 'FAILS'}")
    lines.append(f"  instrument (ii) beta = 0 chain against {z['direct']} direct uniform orders: S {zero_chain['S']:.3f} vs "
                 f"{direct['S']:.3f} (±{direct['seS']:.3f}), chi {zero_chain['chi']:.4f} vs {direct['chi']:.4f}, h "
                 f"{zero_chain['h']:.3f} vs {direct['h']:.3f} -> {'OK' if inst_zero else 'FAILS'}")

    # decision
    ok = {}
    for ph in ("I", "II"):
        ref, tol = p["ref"][ph], p["tol"][ph]
        b = float(ref["beta"])
        if b not in scan:
            ok[ph] = False
            lines.append(f"  phase {ph}: beta = {b:g} not on the grid")
            continue
        m = scan[b]["mean"]
        checks = {k: abs(m[k] - float(ref[k])) <= float(tol[k]) for k in ("chi", "h", "S")}
        ok[ph] = all(checks.values())
        lines.append(f"  phase {ph} (beta = {b:g}): chi {m['chi']:.4f} vs {ref['chi']} (±{tol['chi']}), h {m['h']:.3f} vs "
                     f"{ref['h']} (±{tol['h']}), S {m['S']:+.3f} vs {ref['S']:+} (±{tol['S']}), d_MM {m['d_mm']:.3f} -> "
                     f"{'OK' if ok[ph] else 'NO'}")
        lines.append(f"      interval abundances N_0..N_10: " + " ".join(f"{x:.1f}" for x in m["Nn"]))
    b_peak = max(betas, key=lambda b: scan[b]["mean"]["varS"])
    lo, hi = p["beta_c_window"]
    ok["transition"] = lo <= b_peak <= hi
    lines.append(f"  transition: largest var(S) at beta = {b_peak:g} (window [{lo}, {hi}], published 2.2-2.5) -> "
                 f"{'OK' if ok['transition'] else 'NO'}")
    if not (inst_exact and inst_zero):
        status = "INSTRUMENT_FAILED"
    elif all(ok.values()):
        status = "PASSED"
    else:
        status = "FAILED"
    data = {"scan": {str(b): v for b, v in scan.items()}, "b_peak": b_peak, "exact": exact, "zero_chain": zero_chain,
            "direct": direct, "checks": ok}
    return Result(status, lines, data)


def _dispatch(kind, args):
    """Module-level task router (pickling for the worker processes)."""
    if kind == "chain":
        return chain_task(*args)
    if kind == "direct":
        return direct_task(*args)
    return exact_check(*args)
