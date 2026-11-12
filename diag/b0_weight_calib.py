"""
B0, development stage 4: calibration of the Green weight of the longest chain (decision 17; decides nothing)
==========================================================================================================
Decision 17 (11 Nov 2026, 18:45):
  - condition E[w(L) | tau] = G(tau) = 1 / (2 pi tau), not E[G | L];
  - scaling w(L; rho) = rho^(1/3) w~(L), s = rho^(1/3) tau;
  - form w~(1) = a, w~(L) = c / (L + beta) for L >= 2, with a > 0, c > 0, beta > -2;
  - two moment conditions on r(s) = E[w~(L) | s] - 1/(2 pi s):  int_0^s_max s r ds = 0,  int_0^s_max s^2 r ds = 0.

Declared here, before the run (to be frozen with the calibration, never changed by source):
  - s_max = 6. In [0, 6] the residual is far from zero pointwise (L = 1, 2, 3 dominate); beyond 6 the scratch residual
    was within 3 %, so the leading error of a smooth source lives in [0, 6].
  - The third parameter: beta minimises the rms of the relative residual E[w~ | s] / E[G | s] - 1 over bins of width
    0.5 in s in [2, 20] (where L is not dominated by links); for each beta, a and c solve the two moment conditions.
  - Calibration sample: 40 000 intervals with s ~ U[0, 6] and 8 000 with s ~ U[2, 20]; an independent sample of the
    same size (another seed) for the out-of-sample check of the moments and of the pointwise residual.
  - Uncertainty of (a, c, beta): bootstrap over the intervals of the calibration sample (200 resamples).
Intervals are sprinkled at unit density in their rest frame (the sprinkling is Lorentz invariant; L depends only on s).
The density scaling and the boost independence are checked separately on pyramid sprinklings with the source DP
(diag/b0_weight_valid.py). L counts links (a link has L = 1).

    python -m diag.b0_weight_calib
"""

from __future__ import annotations

import json
import math
import os
import time

import numpy as np

from diag.b0_chain_calib import _longest

S_MAX = 6.0
TAIL = (2.0, 20.0)
N_MAIN, N_TAIL, N_BOOT = 40000, 8000, 200


def chain_length(rng, s):
    n = rng.poisson(math.pi / 12.0 * s ** 3)
    pts = np.empty((0, 3))
    while pts.shape[0] < n:
        m = 4 * (n - pts.shape[0]) + 10
        tt = rng.uniform(0, s, m)
        xx = rng.uniform(-s / 2, s / 2, m)
        yy = rng.uniform(-s / 2, s / 2, m)
        r = np.hypot(xx, yy)
        ok = (r < tt) & (r < s - tt)
        pts = np.vstack([pts, np.column_stack([tt[ok], xx[ok], yy[ok]])])
    pts = pts[:n]
    tt = np.concatenate([[0.0], pts[:, 0], [s]])
    xx = np.concatenate([[0.0], pts[:, 1], [0.0]])
    yy = np.concatenate([[0.0], pts[:, 2], [0.0]])
    o = np.argsort(tt, kind="stable")
    return int(_longest(tt[o], xx[o], yy[o]))


def sample(seed):
    rng = np.random.default_rng(seed)
    s1 = rng.uniform(0, S_MAX, N_MAIN)
    s2 = rng.uniform(*TAIL, N_TAIL)
    L1 = np.array([chain_length(rng, s) for s in s1])
    L2 = np.array([chain_length(rng, s) for s in s2])
    return s1, L1, s2, L2


def weights(L, a, c, beta):
    return np.where(L == 1, a, c / (L + beta))


def solve_ac(s1, L1, beta):
    e1 = (L1 == 1).astype(float)
    e2 = np.where(L1 >= 2, 1.0 / (L1 + beta), 0.0)
    M = np.array([[S_MAX * np.mean(s1 ** k * e1), S_MAX * np.mean(s1 ** k * e2)] for k in (1, 2)])
    rhs = np.array([S_MAX ** k / (2 * math.pi * k) for k in (1, 2)])
    return np.linalg.solve(M, rhs)


def tail_bins(s1, L1, s2, L2):
    s = np.concatenate([s1[s1 >= TAIL[0]], s2])
    L = np.concatenate([L1[s1 >= TAIL[0]], L2])
    edges = np.arange(TAIL[0], TAIL[1] + 1e-9, 0.5)
    k = np.clip(np.digitize(s, edges) - 1, 0, edges.size - 2)
    return s, L, edges, k


def rel_resid(s, L, edges, k, a, c, beta):
    w = weights(L, a, c, beta)
    out = []
    for b in range(edges.size - 1):
        m = k == b
        g = math.log(edges[b + 1] / edges[b]) / (2 * math.pi * (edges[b + 1] - edges[b]))  # E[G] over the bin
        out.append(np.mean(w[m]) / g - 1.0)
    return np.array(out)


def fit(s1, L1, s2, L2):
    s, L, edges, k = tail_bins(s1, L1, s2, L2)

    def cost(beta):
        a, c = solve_ac(s1, L1, beta)
        if a <= 0 or c <= 0:
            return 1e9
        return float(np.sqrt(np.mean(rel_resid(s, L, edges, k, a, c, beta) ** 2)))

    grid = np.arange(-1.95, 5.0, 0.05)
    cs = [cost(b) for b in grid]
    b0 = grid[int(np.argmin(cs))]
    lo, hi = max(-1.999, b0 - 0.05), b0 + 0.05
    for _ in range(40):                                  # golden section
        m1, m2 = lo + 0.382 * (hi - lo), lo + 0.618 * (hi - lo)
        if cost(m1) < cost(m2):
            hi = m2
        else:
            lo = m1
    beta = 0.5 * (lo + hi)
    a, c = solve_ac(s1, L1, beta)
    return float(a), float(c), float(beta), cost(beta)


def moments(s1, L1, a, c, beta):
    """Out-of-sample moments of r over [0, s_max], relative to the moment of G, with their SE."""
    w = weights(L1, a, c, beta)
    out = []
    for k in (1, 2):
        gk = S_MAX ** k / (2 * math.pi * k)
        v = S_MAX * s1 ** k * w
        out.append(((np.mean(v) - gk) / gk, np.std(v, ddof=1) / math.sqrt(v.size) / gk))
    return out


def pointwise(s1, L1, s2, L2, a, c, beta):
    rows = []
    edges = list(np.arange(0.0, S_MAX + 1e-9, 0.5)) + [7, 8, 10, 12, 15, 20]
    s = np.concatenate([s1, s2])
    L = np.concatenate([L1, L2])
    w = weights(L, a, c, beta)
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (s >= lo) & (s < hi)
        if lo == 0.0:
            g = None                                     # E[G] diverges logarithmically in [0, 0.5]; report w only
        else:
            g = math.log(hi / lo) / (2 * math.pi * (hi - lo))
        # (when the two samples overlap the bin, s is not uniform in it; use E[G | s] per interval instead)
        gs = 1.0 / (2 * math.pi * s[m])
        rows.append((lo, hi, int(m.sum()), float(np.mean(w[m])), float(np.std(w[m], ddof=1) / math.sqrt(m.sum())),
                     float(np.mean(w[m]) / np.mean(gs)) if lo > 0 else float("nan"),
                     float(np.mean(L[m] == 1))))
    return rows


def link_part(s1, L1, a):
    """The links' share of the error: int_0^s_max m(s) [a - G(s)] 1{L = 1} ds, relative to int m G, for m = s, s^2."""
    e1 = (L1 == 1)
    out = []
    for k in (1, 2):
        gk = S_MAX ** k / (2 * math.pi * k)
        v = S_MAX * s1 ** k * e1 * (a - 1.0 / (2 * math.pi * s1))
        out.append((float(np.mean(v) / gk), float(np.std(v, ddof=1) / math.sqrt(v.size) / gk)))
    return out


def main():
    t0 = time.time()
    cal = sample(20261111_17)
    val = sample(20261111_18)
    print(f"samples drawn ({time.time() - t0:.0f}s)")
    a, c, beta, rr = fit(*cal)
    rng = np.random.default_rng(7)
    boot = []
    for _ in range(N_BOOT):
        i1 = rng.integers(0, N_MAIN, N_MAIN)
        i2 = rng.integers(0, N_TAIL, N_TAIL)
        boot.append(fit(cal[0][i1], cal[1][i1], cal[2][i2], cal[3][i2])[:3])
    boot = np.array(boot)
    se = boot.std(axis=0, ddof=1)
    L_ = [f"decision 17: w~(1) = a, w~(L >= 2) = c / (L + beta); s_max {S_MAX:g}; beta from the rms relative residual over "
          f"s in [{TAIL[0]:g}, {TAIL[1]:g}] (bins 0.5)",
          f"a = {a:.4f} +- {se[0]:.4f}; c = {c:.4f} +- {se[1]:.4f}; beta = {beta:.3f} +- {se[2]:.3f} "
          f"(bootstrap {N_BOOT}); rms relative tail residual {rr:.4f} (calibration sample)",
          f"  equivalent alpha of the tail 2 pi c = {2 * math.pi * c:.3f}"]
    for name, smp in (("calibration", cal), ("independent", val)):
        m = moments(smp[0], smp[1], a, c, beta)
        L_.append(f"moments of r over [0, {S_MAX:g}] / moments of G ({name} sample): first {m[0][0]:+.4f} +- {m[0][1]:.4f}, "
                  f"second {m[1][0]:+.4f} +- {m[1][1]:.4f}")
    lp = link_part(val[0], val[1], a)
    L_.append(f"links alone (L = 1, independent sample): their share of the error, relative to the full moment of G: "
              f"first {lp[0][0]:+.4f} +- {lp[0][1]:.4f}, second {lp[1][0]:+.4f} +- {lp[1][1]:.4f} "
              f"(compensated by L >= 2 in the moments)")
    L_.append("pointwise, independent sample: bin of s, n, E[w~], SE, E[w~] / E[G], P(L = 1)")
    for lo, hi, n, mw, sw, ratio, p1 in pointwise(*val, a, c, beta):
        L_.append(f"  [{lo:g}, {hi:g}): n {n}, E[w~] {mw:.4f} +- {sw:.4f}, ratio {ratio:.3f}, P(L=1) {p1:.3f}")
    print("\n".join(L_))
    os.makedirs(os.path.join("diag", "out"), exist_ok=True)
    fn = os.path.join("diag", "out", f"b0_weight_calib_{time.strftime('%Y%m%dT%H%M%S', time.gmtime())}.json")
    with open(fn, "w") as fh:
        json.dump({"a": a, "c": c, "beta": beta, "se": se.tolist(), "s_max": S_MAX, "summary": L_,
                   "boot": boot.tolist()}, fh)
    print(f"raw records: {fn} ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
