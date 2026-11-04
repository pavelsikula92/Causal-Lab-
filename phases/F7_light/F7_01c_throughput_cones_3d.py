"""
F7.1c -- Conserved causal throughput in 3+1: does the discrete order give gamma = 1?
====================================================================================

Why this test exists (27 Oct 2026, author's decision 12; docs/derivations/why_gamma1.md)
-----------------------------------------------------------------------------------------
The principle of conserved causal throughput (q = S c_eff constant: where the cones narrow, the levels of the light
cone saturate) predicts gamma = 2/(d_s - 1). In 2+1 the discrete order gave ~1.85 in two runs (F7.1b, prediction 2;
alternatives 1/2 and 1 far away). In 3+1, the dimension of our world, the prediction is gamma = 1, general
relativity's value, and the alternatives are -1 (fixed order) and 1/3 (narrowing cones without new events). The
blind reader in 3+1 is K.16c (PASSED: -1.001 / -1.001, 0.325 / 0.249, 0.994 / 1.013).

Prediction (stated before any run): gamma = 1. Continuum expectations for the readers (3+1, sigma = sigma0
exp(-r^2/32)): density d = sigma, cones s = -sigma, phi = psi = -sigma/2. For the cones-only configuration (no new
events): d = 0, gamma = 1/3, phi = -3 sigma/4, psi = -sigma/4.

Model (locked before any run)
-----------------------------
3+1 box t in [0, 9], |x|, |y|, |z| <= 9. Events: Poisson of coordinate density S(r) = exp(sigma(r)) (the saturated
levels hold more events; nothing marks the added ones). Links (core.order.links, rho_min = 1). Capacity rule of
F7.1a in 3+1: the link v -> k is kept iff its speed |dx|/dt <= c_eff(r_v) = exp(-sigma(r_v)), in the rest frame of the
source. The relation is reachability through kept links; virtual diamond endpoints relate directly to the window
within their capacity cone (as F7.1a/b). So S c_eff = 1 everywhere: the throughput is conserved by construction, and
the test reads which geometry the discrete order has.
Readers: K.16c (pooled cone boundary, diamonds of height 3 with the 4-volume law, density per shell).
Reader controls in the same run (K.16c generator, exact metrics, a = 0.1): ctrl_g1 (gamma = 1, the prediction) and
ctrl_g13 (gamma = 1/3, the cones-only alternative).
Configurations: flat (sigma0 = 0: the null), throughput (sigma0 = 0.2: phi(0) = -0.1 if the prediction holds, the
amplitude of the controls), throughput_strong (sigma0 = 0.3, reported), cones_only (sigma0 = 0.2, uniform events:
F7.1a in 3+1, reported, prediction 1/3), ctrl_g1, ctrl_g13; 16 sprinklings each.

Criteria (stated before any run; those of F7.1b)
------------------------------------------------
INSTRUMENT_FAILED if (i) flat: some reading of some shell |mean| > 0.01 + 3 SE, or (ii) the reader cannot tell
gamma = 1 from gamma = 1/3 on the exact metrics: for some estimate |gamma(ctrl_g1) - gamma(ctrl_g13)| - 2 SE <= 1/3.
UNDECIDED (not metric-like) if |gamma_A - gamma_B| - 2 SE > 0.25 for throughput.
Against the exact gamma = 1 metric read by the same reader:
  PASSED   |gamma(throughput) - gamma(ctrl_g1)| + 2 SE <= 0.25 for both estimates;
  FAILED   |gamma(throughput) - gamma(ctrl_g1)| - 2 SE > 0.25 for both estimates;
  UNDECIDED otherwise.
Reported: the raw gamma, the shell profiles against the continuum, sigma0 = 0.3, and the cones-only configuration.

Development (after everything above was written; decide nothing): the quick run (2 sprinklings in a box |x| <= 6) is
mechanically in order (flat reads zero; the cone reading follows -sigma to the third decimal; links kept 37 %); its
controls are too noisy to tell anything (the reader is calibrated on the full box, K.16c). One full-size 3+1 link
computation takes ~2 min for ~55 000 events. The speed of the links is computed in chunks to bound the memory.
No criterion changed.
"""

from __future__ import annotations

import math
import time

import numpy as np

from core.order.links import links_csr
from core.order.sprinkle import derive_seed
from core.stats.parallel import parallel_map
from core.testing import Result
from phases.F7_light.F7_01a_capacity_cones import CapacityOrder
from phases.K_instruments import K_16c_gamma_reader_3d as k16c

TEST = {
    "id": "F7.1c",
    "phase": "F7",
    "title": "Conserved causal throughput in 3+1: does the discrete order give γ = 1?",
    "theory_part": "why γ = 1: the principle of conserved throughput (decision 12) in the dimension of our world",
    "question": "When in 3+1 the cones narrow to c_eff = exp(−σ) and the density of events grows to exp(σ) "
                "(throughput conserved), does the blind instrument K.16c read γ = 1? Alternatives: −1 (fixed order), "
                "1/3 (narrower cones only).",
    "method": "3+1 sprinkling with coordinate density exp(σ(r)), σ = σ0·exp(−r²/32); links; the capacity rule (link "
              "kept when the speed <= exp(−σ) at the parent); relation = reachability; reading K.16c; controls on "
              "exact metrics γ = 1 and 1/3 (a = 0.1); σ0 = 0, 0.2, 0.3 and narrower cones only 0.2; 16 sprinklings",
    "prediction": "γ = 1 in 3+1 (ln(S c) = 2φ − 2ψ = 0); readings: density σ, cones −σ, φ = ψ = −σ/2; narrower cones "
                  "only: γ = 1/3",
    "criterion": "INSTRUMENT: the flat configuration reads non-zero, or the instrument does not separate the exact γ "
                 "= 1 from 1/3 (|difference| − 2 SE <= 1/3). UNDECIDED not metric (|γ_A − γ_B| − 2 SE > 0.25). "
                 "PASSED: |γ − γ(control 1)| + 2 SE <= 0.25 both paths; FAILED: |γ − γ(control 1)| − 2 SE > 0.25 "
                 "both paths; otherwise UNDECIDED.",
    "null_model": "σ0 = 0 (uniform events, all links)",
    "level": "mechanism",
    "meters": ["K.16c"],
    "depends_on": ["K.16c", "F7.1b"],
    "params": {"T": 9.0, "half": 9.0, "width": 4.0, "bins": [0, 1.5, 3.0, 4.5, 6.0, 7.5],
               "configs": [["flat", "throughput", 0.0], ["throughput", "throughput", 0.2],
                           ["throughput_strong", "throughput", 0.3], ["cones_only", "cones", 0.2],
                           ["ctrl_g1", "metric", 1.0], ["ctrl_g13", "metric", 1.0 / 3.0]],
               "ctrl_a": 0.1, "n_seeds": 16, "anchors": 100, "diamonds": 200, "T_d": 3.0, "cone_dt": [1.0, 2.0],
               "cone_q": 1.5, "tol_gamma": 0.25, "tol_separate": 1.0 / 3.0, "null_abs": 0.01, "null_se": 3.0,
               "n_se": 2.0},
    "parts": {"main": {}},
    "quick": {"n_seeds": 2, "anchors": 30, "diamonds": 40, "half": 6.0, "T": 7.0, "bins": [0, 1.5, 3.0, 4.5]},
    "estimate": "~20–40 min with 8–10 processes (3+1 links: ~2 min per sprinkling of ~55 000 events, 64 of them; the "
                "controls are cheap); memory ~0.5 GB per process",
}


class CapacityOrder3(CapacityOrder):
    """CapacityOrder of F7.1a with the speed of a link measured in three spatial dimensions."""

    def __init__(self, X, sigma0, width):
        self.X = X
        self.N = len(X)
        self.sigma0, self.width = float(sigma0), float(width)
        ptr, idx = links_csr(X, 1.0)
        child = np.repeat(np.arange(self.N), np.diff(ptr))
        parent = idx
        keep = np.empty(parent.size, dtype=bool)
        for a in range(0, parent.size, 1_000_000):              # in chunks: a 3+1 sprinkling has ~7e6 links
            pa, ch = parent[a:a + 1_000_000], child[a:a + 1_000_000]
            d = X[ch] - X[pa]
            keep[a:a + 1_000_000] = np.sqrt((d[:, 1:] ** 2).sum(axis=1)) / d[:, 0] <= self.c_at(X[pa, 1:])
        self.n_links, self.n_kept = int(keep.size), int(keep.sum())
        parent, child = parent[keep], child[keep]
        o = np.argsort(parent, kind="stable")
        self.cptr = np.concatenate([[0], np.cumsum(np.bincount(parent, minlength=self.N))]).astype(np.int64)
        self.cidx = child[o]
        o = np.argsort(child, kind="stable")
        self.pptr = np.concatenate([[0], np.cumsum(np.bincount(child, minlength=self.N))]).astype(np.int64)
        self.pidx = parent[o]


def sprinkle_density(sigma0, width, T, H, rng):
    """Poisson events of coordinate density exp(sigma0 exp(-r^2 / 2 w^2)) in 3+1 (thinning), sorted by time."""
    smax = math.exp(max(sigma0, 0.0))
    n = rng.poisson(smax * T * (2 * H) ** 3)
    X = np.column_stack([rng.uniform(0, T, n)] + [rng.uniform(-H, H, n) for _ in range(3)])
    S = np.exp(sigma0 * np.exp(-(X[:, 1:] ** 2).sum(axis=1) / (2.0 * width ** 2)))
    X = X[rng.random(n) < S / smax]
    return X[np.argsort(X[:, 0])]


def order_task(p, sigma0, saturate, seed):
    rng = np.random.default_rng(seed)
    T, H = float(p["T"]), float(p["half"])
    X = sprinkle_density(sigma0 if saturate else 0.0, p["width"], T, H, rng)
    order = CapacityOrder3(X, sigma0, p["width"])
    edges = np.array(p["bins"], float)
    T_d = float(p["T_d"])

    def related_from(i, cand):
        fut = order.future(order.cidx[order.cptr[i]:order.cptr[i + 1]], X[i, 0] + p["cone_dt"][1])
        return fut[cand]

    def diamond_count(pp, qq):
        dxp = np.sqrt(((X[:, 1:] - pp[1:]) ** 2).sum(axis=1))
        c_p = float(order.c_at(pp[None, 1:])[0])
        first = np.flatnonzero((X[:, 0] > pp[0]) & (X[:, 0] < qq[0]) & (dxp <= c_p * (X[:, 0] - pp[0])))
        last = np.flatnonzero((X[:, 0] < qq[0]) & (X[:, 0] > pp[0]) & (dxp <= order.c_at(X[:, 1:]) * (qq[0] - X[:, 0])))
        if first.size == 0 or last.size == 0:
            return 0
        return int((order.future(first, qq[0]) & order.past(last, pp[0])).sum())

    d = k16c.read_density3(X, edges, T)
    s = k16c.read_cones3(X, related_from, edges, T, int(p["anchors"]), p["cone_dt"], float(p["cone_q"]), rng)
    pc = k16c.read_clocks3(X, diamond_count, edges, T, int(p["diamonds"]), T_d, rng)
    out = k16c.estimates3(d, s, pc)
    out["N"] = len(X)
    out["kept"] = order.n_kept / max(order.n_links, 1)
    return out


def _task(kind, p, val, seed):
    if kind == "throughput":
        return order_task(p, float(val), True, seed)
    if kind == "cones":
        return order_task(p, float(val), False, seed)
    out = k16c.seed_task(p, float(p["ctrl_a"]), float(val), seed)
    out["kept"] = float("nan")
    return out


def run(p, ctx):
    t0 = time.time()
    n = int(p["n_seeds"])
    jobs, keys = [], []
    for name, kind, val in p["configs"]:
        for k in range(n):
            jobs.append((kind, p, val, derive_seed(ctx.seed, "F7.1c", name, k)))
            keys.append(name)
    res = parallel_map(_task, jobs, ctx.workers, ctx.progress,
                       lambda i, r_, dd, nn: f"[F7.1c] {keys[i]} sprinkling {dd}/{nn} ({time.time() - t0:.0f}s)")
    edges = np.array(p["bins"], float)
    centres = 0.5 * (edges[:-1] + edges[1:])
    g = np.exp(-centres ** 2 / (2 * p["width"] ** 2))
    lines = [f"3+1, t [0, {p['T']:g}], |x|, |y|, |z| <= {p['half']:g}; sigma = sigma0 exp(-r^2/(2*{p['width']:g}^2)); throughput: "
             f"events of coordinate density exp(sigma) and links kept iff speed <= exp(-sigma) at the parent; cones_only: the "
             f"same links, uniform events; readers K.16c; controls: exact metrics a = {p['ctrl_a']:g}, gamma = 1 and 1/3; "
             f"{n} sprinklings per configuration"]
    st = {}
    for name, kind, val in p["configs"]:
        rows = [r_ for r_, k in zip(res, keys) if k == name]
        st[name] = {k: k16c.mean_se(rows, k) for k in ("d", "s", "phic", "phiA", "psiA", "phiB", "psiB", "gA", "gB", "N")}
        x = st[name]
        kept = float(np.nanmean([r_["kept"] for r_ in rows]))
        lines.append(f"  {name} ({kind} {val:.3g}): ~{x['N'][0]:.0f} events" + (f", links kept {kept:.3f}" if np.isfinite(kept) else "")
                     + f"; gamma A (cones + density) {x['gA'][0]:+.3f} ± {x['gA'][1]:.3f}; gamma B (clocks + cones) "
                     f"{x['gB'][0]:+.3f} ± {x['gB'][1]:.3f}")
        lines.append("        r            " + " ".join(f"{c:7.2f}" for c in centres))
        if kind in ("throughput", "cones"):
            sg = float(val) * g
            d_exp = sg if kind == "throughput" else 0 * sg
            phi_exp = -sg / 2 if kind == "throughput" else -3 * sg / 4
            lines.append("        d expected   " + " ".join(f"{v:+7.4f}" for v in d_exp))
            lines.append("        d read       " + " ".join(f"{v:+7.4f}" for v in x["d"][0]))
            lines.append("        s expected   " + " ".join(f"{v:+7.4f}" for v in -sg))
            lines.append("        s read       " + " ".join(f"{v:+7.4f}" for v in x["s"][0]))
            lines.append("        phi expected " + " ".join(f"{v:+7.4f}" for v in phi_exp))
            lines.append("        phi clock    " + " ".join(f"{v:+7.4f}" for v in x["phic"][0]))
        else:
            lines.append("        phi true     " + " ".join(f"{v:+7.4f}" for v in -float(p["ctrl_a"]) * g))
        lines.append("        phi A / psi A " + " ".join(f"{a:+.3f}/{b:+.3f}" for a, b in zip(x["phiA"][0], x["psiA"][0])))
        lines.append("        phi B / psi B " + " ".join(f"{a:+.3f}/{b:+.3f}" for a, b in zip(x["phiB"][0], x["psiB"][0])))
    nse, tol = float(p["n_se"]), float(p["tol_gamma"])
    inst_null, worst = k16c.null_check(st["flat"], p)
    sep_ok = True
    for e in ("gA", "gB"):
        (a1, s1), (a3, s3) = st["ctrl_g1"][e], st["ctrl_g13"][e]
        sep_ok &= abs(a1 - a3) - nse * math.hypot(s1, s3) > float(p["tol_separate"])
    lines.append(f"  instrument (i) flat: every reading |mean| <= {p['null_abs']} + {p['null_se']:g} SE -> "
                 f"{'OK' if inst_null else 'FAILS'} (largest |mean| - {p['null_se']:g} SE = {worst:+.4f})")
    lines.append(f"  instrument (ii) exact metrics: gamma 1 reads A {st['ctrl_g1']['gA'][0]:+.3f}, B {st['ctrl_g1']['gB'][0]:+.3f}; "
                 f"gamma 1/3 reads A {st['ctrl_g13']['gA'][0]:+.3f}, B {st['ctrl_g13']['gB'][0]:+.3f} -> "
                 f"{'told apart' if sep_ok else 'NOT told apart'} (needs |difference| - 2 SE > 1/3)")
    x = st["throughput"]
    (gA, sA), (gB, sB) = x["gA"], x["gB"]
    metric_like = abs(gA - gB) - nse * math.hypot(sA, sB) <= tol
    dev = {}
    for e, (gv, sv) in (("gA", (gA, sA)), ("gB", (gB, sB))):
        c1, cs = st["ctrl_g1"][e]
        dev[e] = (gv - c1, math.hypot(sv, cs))
    passed = all(abs(d_) + nse * s_ <= tol for d_, s_ in dev.values())
    failed = all(abs(d_) - nse * s_ > tol for d_, s_ in dev.values())
    lines.append(f"  throughput (sigma0 = 0.2): |gamma_A - gamma_B| = {abs(gA - gB):.3f} ± {math.hypot(sA, sB):.3f} -> "
                 f"{'metric-like' if metric_like else 'NOT metric-like'}; against the exact gamma = 1 metric: A "
                 f"{dev['gA'][0]:+.3f} ± {dev['gA'][1]:.3f}, B {dev['gB'][0]:+.3f} ± {dev['gB'][1]:.3f}; raw against 1: "
                 f"A {gA - 1:+.3f}, B {gB - 1:+.3f}")
    for name, label in (("throughput_strong", "sigma0 = 0.3"), ("cones_only", "cones only, sigma0 = 0.2 (prediction 1/3)")):
        y = st[name]
        lines.append(f"  {label} (reported): gamma A {y['gA'][0]:+.3f} ± {y['gA'][1]:.3f}, B {y['gB'][0]:+.3f} ± {y['gB'][1]:.3f}")
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
