"""
K.16b -- The blind gamma reader, second version: pooled cone boundary, smaller diamonds, new sprinklings
========================================================================================================

Why this version exists (27 Oct 2026, after K.16)
-------------------------------------------------
K.16 (calibration) FAILED by a hair: estimate B (clocks + cones) for gamma = +1 gave
1.092 +- 0.080, |dev| + 2 SE = 0.253 against 0.25; every other check passed (null
control OK; gamma = -1: -0.911 / -0.913; gamma = 0: +0.050 / +0.107; gamma = +1 A: +0.905).
A diagnostic on the SAME sprinklings (decides nothing) showed that the shrinking of
|gamma| by ~0.1 is NOT regression dilution (the ratio of sums gives -0.904 and -0.899 for
gamma = -1 as well) but two biases of the readers themselves:
  (1) cones: the per-anchor midpoint between the largest related and the next unrelated
      ratio q = |dx|/dt reads the light speed low by ~0.005 even without a field (null
      s = -0.0041 +- 0.0022 in the first bin). An offset of s enters phi_A = (2 s + d)/3 for
      every bin, and in the outer bins, where |phi| is small, it pulls gamma towards 0;
  (2) clocks: a diamond of height 5 reaches 3.75 from its axis and averages the field
      over that width; in the centre it reads |phi| too small (gamma = 0: -0.087 against
      -0.099), which biases psi_B = s - phi_clock.
K.16 keeps its verdict (R132). This is a revised instrument with new sprinklings (own seed
stream), as K.11 -> K.11b -> K.11c.

Changes against K.16 (and nothing else)
---------------------------------------
  cones   the (q, related) pairs of ALL anchors of a bin are pooled and the boundary is the
          c that minimises the misclassifications #(related, q > c) + #(unrelated, q <= c)
          (100 anchors per bin); with many pairs the gap around the boundary, and the
          midpoint offset with it, shrinks;
  clocks  diamonds of coordinate height T_d = 3 (reach 2.25), 200 per bin;
  stats   32 sprinklings per configuration;
  estimate of gamma: as K.16 (least squares through the origin per sprinkling, mean and SE).
Generator, field, configurations, bins, the density reader and every criterion and
tolerance are those of K.16 (imported unchanged).

Criteria (stated before any run; identical to K.16)
---------------------------------------------------
INSTRUMENT_FAILED if in the null configuration some reading of some bin has
|mean| > 0.01 + 3 SE. PASSED if for gamma_true = -1, 0, +1 (a = 0.1) both estimates satisfy
|gamma_est - gamma_true| + 2 SE <= 0.25; FAILED otherwise. Reported: the per-bin profiles,
the a = 0.03 resolution and the size of the remaining biases (mean minus truth).

Quick run (3 sprinklings; decide nothing): mechanics in order; the null check fails with 3 sprinklings (SE from three
values). Nothing was changed after it.
"""

from __future__ import annotations

import math
import time

import numpy as np

from core.order.sprinkle import derive_seed
from core.stats.parallel import parallel_map
from core.testing import Result
from phases.K_meridla import K_16_gamma_reader as k16

TEST = {
    "id": "K.16b",
    "phase": "K",
    "title": "Slepé čtení γ, druhá verze: sdružená hranice kuželu, menší diamanty, nové sprinklingy",
    "theory_part": k16.TEST["theory_part"],
    "question": "Odstraní sdružená hranice relace (místo středu mezery u každé kotvy) a menší diamanty systematické "
                "posuny K.16 a určí měřidlo γ = −1, 0, +1 v toleranci 0,25 oběma cestami?",
    "method": "generátor, pole, konfigurace, biny a čtení hustoty z K.16 beze změny; kužele: hranice c minimalizující "
              "chybné zařazení všech párů binu (100 kotev); hodiny: diamanty výšky 3, 200 na bin; 32 sprinklingů, nový proud "
              "seedů",
    "prediction": "měřidlo vrátí vložené γ (kalibrace)",
    "criterion": k16.TEST["criterion"],
    "null_model": "a = 0 (bez pole)",
    "level": "calibration",
    "meters": ["K.16"],
    "depends_on": ["K.16"],
    "params": {**k16.TEST["params"], "n_seeds": 32, "anchors": 100, "diamonds": 200, "T_d": 3.0},
    "parts": {"main": {}},
    "quick": {"n_seeds": 3, "anchors": 30, "diamonds": 40},
    "estimate": "~1–2 min s 10 procesy",
}


def read_cones_pooled(X, rel, edges, T, n_anchor, dt_rng, qmax, rng):
    """s_b = ln c_b, c_b the boundary that misclassifies the fewest (q, related) pairs pooled over the bin's anchors."""
    r = np.hypot(X[:, 1], X[:, 2])
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        idx = np.flatnonzero((r >= lo) & (r < hi) & (X[:, 0] >= 1.0) & (X[:, 0] <= T - 1.0 - dt_rng[1]))
        if idx.size == 0:
            out.append(float("nan"))
            continue
        qs, oks = [], []
        for i in rng.choice(idx, min(n_anchor, idx.size), replace=False):
            d = X - X[i]
            dt = d[:, 0]
            dx = np.hypot(d[:, 1], d[:, 2])
            cand = np.flatnonzero((dt >= dt_rng[0]) & (dt <= dt_rng[1]) & (dx <= qmax * dt))
            if cand.size == 0:
                continue
            qs.append(dx[cand] / dt[cand])
            oks.append(rel(np.repeat(X[i][None, :], cand.size, axis=0), X[cand]))
        if not qs:
            out.append(float("nan"))
            continue
        q, ok = np.concatenate(qs), np.concatenate(oks)
        o = np.argsort(q)
        q, ok = q[o], ok[o]
        # errors if the boundary lies just after position k (pairs 0..k are 'inside'):
        #   related beyond k + unrelated up to k
        rel_after = ok[::-1].cumsum()[::-1]                      # related at positions >= k
        unrel_upto = (~ok).cumsum()                              # unrelated at positions <= k
        err = np.append(rel_after[1:], 0) + unrel_upto
        k = int(np.argmin(err))
        c = 0.5 * (q[k] + q[k + 1]) if k + 1 < q.size else q[k]
        out.append(math.log(c))
    return np.array(out)


def seed_task(p, a, gamma, seed):
    rng = np.random.default_rng(seed)
    field = k16.Field(a, gamma, p["width"])
    T, H = float(p["T"]), float(p["half"])
    X = k16.sprinkle(field, T, H, rng)
    rel = k16.make_relation(field)
    edges = np.array(p["bins"], float)
    d = k16.read_density(X, edges, T)
    s = read_cones_pooled(X, rel, edges, T, int(p["anchors"]), p["cone_dt"], float(p["cone_q"]), rng)
    pc = k16.read_clocks(X, rel, edges, T, int(p["diamonds"]), float(p["T_d"]), rng)
    phiA, psiA = (2 * s + d) / 3.0, (s - d) / 3.0
    phiB, psiB = pc, s - pc
    return {"N": len(X), "d": d.tolist(), "s": s.tolist(), "phic": pc.tolist(),
            "gA": k16.gamma_fit(phiA, psiA), "gB": k16.gamma_fit(phiB, psiB),
            "phiA": phiA.tolist(), "psiA": psiA.tolist(), "phiB": phiB.tolist(), "psiB": psiB.tolist()}


def run(p, ctx):
    # the evaluation of K.16 with this module's seed_task and an own seed stream
    t0 = time.time()
    n = int(p["n_seeds"])
    jobs, keys = [], []
    for name, a, g in p["configs"]:
        for k in range(n):
            jobs.append((p, float(a), float(g), derive_seed(ctx.seed, "K.16b", name, k)))
            keys.append(name)
    res = parallel_map(seed_task, jobs, ctx.workers, ctx.progress,
                       lambda i, r_, dd, nn: f"[K.16b] {keys[i]} sprinkling {dd}/{nn} ({time.time() - t0:.0f}s)")
    return k16.evaluate(p, res, keys, n)
