"""
F7.1b -- Conserved causal throughput: narrower cones AND more events. Does the discrete order give gamma = 2/(d_s - 1)?
=====================================================================================================================

Why this test exists (27 Oct 2026, author's decisions 12 and the principle of docs/derivations/proc_gamma1.md)
---------------------------------------------------------------------------------------------------------------
F7.1a: a capacity rule that only narrows the cones (keeps the events) gives gamma = 1/d_s (2+1: 0.47 / 0.43). The
principle of conserved causal throughput (the author's "throughput by saturation: the more saturated a level of the
light cone, the more materialised"; Physics of Causality 5.3, 23, 26.1) adds the missing half: where the cones narrow
to c_eff, the level saturates, the coordinate density of events rises to S = 1/c_eff, so the throughput q = S c_eff
stays constant. Continuum consequence (sympy): ln(S c) = 2 phi - (d_s - 1) psi = 0, gamma = 2/(d_s - 1): 1 in 3+1,
2 in 2+1. This test asks whether a DISCRETE order built this way realises the principle: a consistent metric (clocks
and cones agree) with the predicted gamma.

Prediction (stated before any run): gamma = 2 in 2+1 (not 1: the 3+1 value needs the 3+1 version of this test).
Continuum expectations for the readers (2+1, sigma = sigma0 exp(-r^2/128)): density d = sigma, cones s = -sigma,
clocks phi = -sigma/3, psi = -2 sigma/3.

Model (locked before any run)
-----------------------------
2+1 box t in [0, 14], |x|, |y| <= 18. Events: Poisson with coordinate density S(r) = exp(sigma(r)) (the saturated
levels hold more events; the added events are ordinary events, nothing marks them). Links of this sprinkling
(core.order.links with rho_min = 1); the capacity rule of F7.1a: the link v -> k is kept iff its speed |dx|/dt <=
c_eff(r_v) = exp(-sigma(r_v)), in the rest frame of the source; the relation is reachability through kept links;
virtual diamond endpoints as F7.1a. So S(r) c_eff(r) = 1 everywhere: the throughput is conserved by construction,
and the test checks what geometry the discrete order then has.
Readers: K.16b (pooled cone boundary, diamonds of height 3, density), as in F7.1a.

Reader controls in the same run (K.16b generator: exact metrics, the reader blind):
  ctrl_g2  a = 0.1, gamma = 2   (the predicted geometry, outside the range K.16b was calibrated on: -1 .. 1);
  ctrl_g1  a = 0.1, gamma = 1   (general relativity's value, the alternative to be told apart).
sigma0 = 0.3 gives phi(0) = -0.1 if the prediction holds, the amplitude of the controls.
Configurations: flat (sigma0 = 0: uniform events, all links: the null), throughput (sigma0 = 0.3), throughput_strong
(sigma0 = 0.45, reported), ctrl_g2, ctrl_g1; 16 sprinklings each.

Criteria (stated before any run)
--------------------------------
INSTRUMENT_FAILED if
  (i)  flat: some reading of some bin |mean| > 0.01 + 3 SE, or
  (ii) the reader cannot tell gamma = 2 from gamma = 1 on the exact metrics: for some estimate
       |gamma(ctrl_g2) - gamma(ctrl_g1)| - 2 SE <= 0.5.
UNDECIDED (not metric-like) if |gamma_A - gamma_B| - 2 SE > 0.25 for throughput.
Against the exact gamma = 2 metric read by the same reader (its bias cancels):
  PASSED   |gamma(throughput) - gamma(ctrl_g2)| + 2 SE <= 0.25 for both estimates;
  FAILED   |gamma(throughput) - gamma(ctrl_g2)| - 2 SE > 0.25 for both estimates;
  UNDECIDED otherwise.
Reported: the raw gamma against 2 and against the controls, the per-bin profiles against the continuum expectations,
sigma0 = 0.45.
Development (after everything above was written; decide nothing; nothing changed after it):
  quick run (2 small sprinklings): mechanics in order, the statistics too thin for the instrument checks;
  4 full-size sprinklings at sigma0 = 0.3: readings follow the continuum expectations within their scatter
  (d 0.33 / 0.25 / 0.23 / 0.21 against 0.30 / 0.28 / 0.25 / 0.20 in the first four bins; s -0.36 / -0.30 / -0.27 /
  -0.21 against -0.30 / -0.28 / -0.25 / -0.20; phi_clock -0.075 / -0.096 / -0.086 / -0.072 against -0.099 / -0.093 /
  -0.082 / -0.068); gamma per sprinkling A 1.36-1.93, B 1.80-2.36. The central bin is the noisiest (smallest area).
"""

from __future__ import annotations

import math
import time

import numpy as np

from core.order.sprinkle import derive_seed
from core.stats.parallel import parallel_map
from core.testing import Result
from phases.F7_svetlo import F7_01a_capacity_cones as f71a
from phases.K_meridla import K_16_gamma_reader as k16
from phases.K_meridla import K_16b_gamma_reader_pooled as k16b

TEST = {
    "id": "F7.1b",
    "phase": "F7",
    "title": "Zachovaná kauzální propustnost: užší kužele i víc událostí. Dá diskrétní řád γ = 2/(d_s − 1)?",
    "theory_part": "proč γ = 1: princip zachování propustnosti (rozhodnutí 12, docs/derivations/proc_gamma1.md); "
                   "propustnost podle nasycenosti",
    "question": "Když se kužele zúží na c_eff = exp(−σ) a hustota událostí vzroste na S = 1/c_eff (propustnost S·c_eff "
                "zachována), jaké γ přečte slepé měřidlo? Předpověď ve 2+1: γ = 2 (ve 3+1 by to bylo 1).",
    "method": "sprinkling 2+1 se souřadnicovou hustotou exp(σ(r)), σ = σ0·exp(−r²/128); vazby; pravidlo kapacity F7.1a "
              "(vazba ponechána, když rychlost <= exp(−σ) u rodiče); relace = dosažitelnost; čtení K.16b; kontroly měřidla "
              "na přesných metrikách γ = 2 a γ = 1 (a = 0,1); σ0 = 0, 0,3, 0,45; 16 sprinklingů",
    "prediction": "γ = 2 ve 2+1 (ln(S c) = 2φ − (d_s − 1)ψ = 0); čtení: hustota σ, kužele −σ, hodiny φ = −σ/3",
    "criterion": "MĚŘIDLO: plochá konfigurace čte nenulově (|průměr| > 0,01 + 3 SE), nebo měřidlo nerozliší přesné γ = 2 "
                 "od γ = 1 (|γ2 − γ1| − 2 SE <= 0,5). NEROZHODNUTO nemetrické (|γ_A − γ_B| − 2 SE > 0,25). PROŠEL: "
                 "|γ − γ(kontrola γ = 2)| + 2 SE <= 0,25 obě cesty; NEPROŠEL: |γ − γ(kontrola γ = 2)| − 2 SE > 0,25 obě "
                 "cesty; jinak NEROZHODNUTO.",
    "null_model": "σ0 = 0 (rovnoměrné události, všechny vazby)",
    "level": "mechanism",
    "meters": ["K.16b"],
    "depends_on": ["K.16b", "F7.1a"],
    "params": {"T": 14.0, "half": 18.0, "width": 8.0, "bins": [0, 2, 4, 6, 8, 10, 12, 14],
               "configs": [["flat", "throughput", 0.0], ["throughput", "throughput", 0.3],
                           ["throughput_strong", "throughput", 0.45], ["ctrl_g2", "metric", 2.0], ["ctrl_g1", "metric", 1.0]],
               "ctrl_a": 0.1, "n_seeds": 16, "anchors": 100, "diamonds": 200, "T_d": 3.0, "cone_dt": [1.0, 2.0],
               "cone_q": 1.5, "first_dt": 3.0, "tol_gamma": 0.25, "tol_separate": 0.5, "null_abs": 0.01,
               "null_se": 3.0, "n_se": 2.0},
    "parts": {"main": {}},
    "quick": {"n_seeds": 2, "anchors": 30, "diamonds": 40, "half": 12.0, "bins": [0, 2, 4, 6, 8, 10]},
    "estimate": "~5–15 min s 8–10 procesy (vazby 48 sprinklingů po ~18 000–22 000 událostech, 32 kontrolních sprinklingů)",
}


def sprinkle_saturated(sigma0, width, T, H, rng):
    """Poisson events of coordinate density S(r) = exp(sigma0 exp(-r^2 / 2 w^2)) (thinning), sorted by time."""
    smax = math.exp(max(sigma0, 0.0))
    n = rng.poisson(smax * T * (2 * H) ** 2)
    X = np.column_stack([rng.uniform(0, T, n), rng.uniform(-H, H, n), rng.uniform(-H, H, n)])
    r2 = (X[:, 1:] ** 2).sum(axis=1)
    S = np.exp(sigma0 * np.exp(-r2 / (2.0 * width ** 2)))
    X = X[rng.random(n) < S / smax]
    return X[np.argsort(X[:, 0])]


def throughput_task(p, sigma0, seed):
    rng = np.random.default_rng(seed)
    T, H = float(p["T"]), float(p["half"])
    X = sprinkle_saturated(sigma0, p["width"], T, H, rng)
    order = f71a.CapacityOrder(X, sigma0, p["width"])
    edges = np.array(p["bins"], float)
    d = k16.read_density(X, edges, T)
    s = f71a.read_cones(order, edges, T, int(p["anchors"]), p["cone_dt"], float(p["cone_q"]), rng)
    pc = f71a.read_clocks(order, edges, T, int(p["diamonds"]), float(p["T_d"]), float(p["first_dt"]), rng)
    phiA, psiA = (2 * s + d) / 3.0, (s - d) / 3.0
    phiB, psiB = pc, s - pc
    return {"N": len(X), "d": d.tolist(), "s": s.tolist(), "phic": pc.tolist(),
            "gA": k16.gamma_fit(phiA, psiA), "gB": k16.gamma_fit(phiB, psiB),
            "phiA": phiA.tolist(), "psiA": psiA.tolist(), "phiB": phiB.tolist(), "psiB": psiB.tolist()}


def _task(kind, p, val, seed):
    if kind == "throughput":
        return throughput_task(p, float(val), seed)
    return k16b.seed_task(p, float(p["ctrl_a"]), float(val), seed)


def _mean_se(rows, key):
    a = np.array([r_[key] for r_ in rows], float)
    if a.ndim == 1:
        a = a[np.isfinite(a)]
        return float(a.mean()), float(a.std(ddof=1) / math.sqrt(a.size)) if a.size > 1 else float("nan")
    return np.nanmean(a, axis=0), np.nanstd(a, axis=0, ddof=1) / math.sqrt(a.shape[0])


def run(p, ctx):
    t0 = time.time()
    n = int(p["n_seeds"])
    jobs, keys = [], []
    for name, kind, val in p["configs"]:
        for k in range(n):
            jobs.append((kind, p, val, derive_seed(ctx.seed, "F7.1b", name, k)))
            keys.append(name)
    res = parallel_map(_task, jobs, ctx.workers, ctx.progress,
                       lambda i, r_, dd, nn: f"[F7.1b] {keys[i]} sprinkling {dd}/{nn} ({time.time() - t0:.0f}s)")
    edges = np.array(p["bins"], float)
    centres = 0.5 * (edges[:-1] + edges[1:])
    g = np.exp(-centres ** 2 / (2 * p["width"] ** 2))
    lines = [f"2+1, t [0, {p['T']:g}], |x|, |y| <= {p['half']:g}; sigma = sigma0 exp(-r^2/(2*{p['width']:g}^2)); events of "
             f"coordinate density exp(sigma), links kept iff speed <= exp(-sigma) at the parent (throughput S c = 1); readers "
             f"K.16b; controls: exact metrics a = {p['ctrl_a']:g}, gamma = 2 and 1; {n} sprinklings per configuration"]
    st = {}
    for name, kind, val in p["configs"]:
        rows = [r_ for r_, k in zip(res, keys) if k == name]
        st[name] = {k: _mean_se(rows, k) for k in ("d", "s", "phic", "phiA", "psiA", "phiB", "psiB", "gA", "gB", "N")}
        x = st[name]
        lines.append(f"  {name} ({kind} {val:g}): ~{x['N'][0]:.0f} events; gamma A (cones + density) {x['gA'][0]:+.3f} ± "
                     f"{x['gA'][1]:.3f}; gamma B (clocks + cones) {x['gB'][0]:+.3f} ± {x['gB'][1]:.3f}")
        lines.append("        r            " + " ".join(f"{c:7.1f}" for c in centres))
        if kind == "throughput":
            sg = float(val) * g
            lines.append("        d expected   " + " ".join(f"{v:+7.4f}" for v in sg))
            lines.append("        d read       " + " ".join(f"{v:+7.4f}" for v in x["d"][0]))
            lines.append("        s expected   " + " ".join(f"{v:+7.4f}" for v in -sg))
            lines.append("        s read       " + " ".join(f"{v:+7.4f}" for v in x["s"][0]))
            lines.append("        phi expected " + " ".join(f"{v:+7.4f}" for v in -sg / 3))
            lines.append("        phi clock    " + " ".join(f"{v:+7.4f}" for v in x["phic"][0]))
        else:
            lines.append("        phi true     " + " ".join(f"{v:+7.4f}" for v in -float(p["ctrl_a"]) * g))
        lines.append("        phi A / psi A " + " ".join(f"{a:+.3f}/{b:+.3f}" for a, b in zip(x["phiA"][0], x["psiA"][0])))
        lines.append("        phi B / psi B " + " ".join(f"{a:+.3f}/{b:+.3f}" for a, b in zip(x["phiB"][0], x["psiB"][0])))
    nse, tol = float(p["n_se"]), float(p["tol_gamma"])
    # instrument (i): flat reads nothing
    inst_null, worst = True, -math.inf
    for f in ("d", "s", "phic"):
        m, se = st["flat"][f]
        for mb, sb in zip(m, se):
            bad = not (np.isfinite(mb) and np.isfinite(sb)) or abs(mb) > float(p["null_abs"]) + float(p["null_se"]) * sb
            inst_null &= not bad
            if np.isfinite(mb) and np.isfinite(sb):
                worst = max(worst, abs(mb) - float(p["null_se"]) * sb)
    # instrument (ii): the reader tells gamma 2 from gamma 1
    sep_ok = True
    for e in ("gA", "gB"):
        (a2, s2), (a1, s1) = st["ctrl_g2"][e], st["ctrl_g1"][e]
        sep_ok &= abs(a2 - a1) - nse * math.hypot(s2, s1) > float(p["tol_separate"])
    lines.append(f"  instrument (i) flat: every reading |mean| <= {p['null_abs']} + {p['null_se']:g} SE -> "
                 f"{'OK' if inst_null else 'FAILS'} (largest |mean| - {p['null_se']:g} SE = {worst:+.4f})")
    lines.append(f"  instrument (ii) exact metrics: gamma 2 reads A {st['ctrl_g2']['gA'][0]:+.3f}, B {st['ctrl_g2']['gB'][0]:+.3f}; "
                 f"gamma 1 reads A {st['ctrl_g1']['gA'][0]:+.3f}, B {st['ctrl_g1']['gB'][0]:+.3f} -> "
                 f"{'told apart' if sep_ok else 'NOT told apart'} (needs |difference| - 2 SE > {p['tol_separate']})")
    x = st["throughput"]
    (gA, sA), (gB, sB) = x["gA"], x["gB"]
    metric_like = abs(gA - gB) - nse * math.hypot(sA, sB) <= tol
    dev = {}
    for e, (gv, sv) in (("gA", (gA, sA)), ("gB", (gB, sB))):
        c2, cs = st["ctrl_g2"][e]
        dev[e] = (gv - c2, math.hypot(sv, cs))
    passed = all(abs(d_) + nse * s_ <= tol for d_, s_ in dev.values())
    failed = all(abs(d_) - nse * s_ > tol for d_, s_ in dev.values())
    lines.append(f"  throughput (sigma0 = 0.3): |gamma_A - gamma_B| = {abs(gA - gB):.3f} ± {math.hypot(sA, sB):.3f} -> "
                 f"{'metric-like' if metric_like else 'NOT metric-like'}; against the exact gamma = 2 metric: A "
                 f"{dev['gA'][0]:+.3f} ± {dev['gA'][1]:.3f}, B {dev['gB'][0]:+.3f} ± {dev['gB'][1]:.3f}; raw against 2: "
                 f"A {gA - 2:+.3f}, B {gB - 2:+.3f}")
    xs = st["throughput_strong"]
    lines.append(f"  sigma0 = 0.45 (reported): gamma A {xs['gA'][0]:+.3f} ± {xs['gA'][1]:.3f}, B {xs['gB'][0]:+.3f} ± "
                 f"{xs['gB'][1]:.3f}")
    if not (inst_null and sep_ok):
        status = "INSTRUMENT_FAILED"
    elif not metric_like:
        status = "UNDECIDED"
    elif passed:
        status = "PASSED"
    elif failed:
        status = "FAILED"
    else:
        status = "UNDECIDED"
    data = {name: {k: (list(map(float, v[0])), list(map(float, v[1]))) if isinstance(v[0], np.ndarray) else list(v)
                   for k, v in x_.items()} for name, x_ in st.items()}
    return Result(status, lines, data)
