"""
C0, development: the normalisation of the smooth kernel S16, calibrated and validated on SEPARATE data (13 Nov; decides
nothing; the author's condition of 13 Nov, 14:36: the constants must not be fitted on the data that judge them).
  calibration:  rho 30 and 60, seeds 100, 101 -> a = 1 / mean over t >= 3 of [c sum w / (t^2/2)]
  validation:   rho 30, 45, 60, seeds 200, 201 (rho 45 not used in the calibration) -> the constant-source ratio with
                a applied, and the linear vacuum and k = 1 waves with a applied (reported: the constant is NOT allowed to
                correct the dispersion; those are judged separately)
Box X 7, T 6 (as the C0 development rows).
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import sys
import time

import numpy as np

from core.stats.parallel import parallel_map, default_workers
from diag.c0_links4 import sprinkle4
from diag.c0_layers4 import layers4, evolve_layers
from diag.c0_dev import kernel


def task(rho, seed, a_norm, X=7.0, T=6.0, waves=False):
    t, x, y, z = sprinkle4(rho, X, T, 7000 + seed)
    coef, lam = kernel("S16")
    coef = coef * a_norm
    a0 = time.time()
    ptr, idx, lay = layers4(t, x, y, z, (24 * lam / (math.pi * rho)) ** 0.25, coef.size - 1)
    c = 1 / (2 * math.pi * math.sqrt(6 * rho))
    inner = (np.abs(x) < X - t - 0.3) & (np.abs(y) < X - t - 0.3) & (np.abs(z) < X - t - 0.3)
    tb = np.arange(0.25, T - 0.4, 0.5)
    tot = np.array([c * np.sum(coef[lay[ptr[v]:ptr[v + 1]]]) for v in range(t.size)])
    const = [float(tot[m].mean() / (0.5 * np.mean(t[m] ** 2))) for m in [inner & (np.abs(t - b) < 0.25) for b in tb]]
    out = {"rho": rho, "seed": seed, "a": a_norm, "t": tb.tolist(), "const": const, "N": int(t.size)}
    if waves:
        nopre = np.zeros(t.size, np.bool_)
        for k in (0.0, 1.0):
            A = 0.05
            om = math.sqrt(k * k + 1)
            ex = A * np.cos(k * x) * np.cos(om * t)
            phi = evolve_layers(ptr, idx, lay, coef, A * np.cos(k * x) * np.cos(k * t), nopre, ex, c, 1.0, True)
            out[f"k{k:g}"] = [[float(((phi - ex) / A)[m].mean()), float(((phi - ex) / A)[m].std())]
                              for m in [inner & (np.abs(t - b) < 0.25) for b in tb]]
    out["sec"] = time.time() - a0
    return out


def cached(*a):
    key = hashlib.sha256(json.dumps([str(q) for q in a] + ["v1"]).encode()).hexdigest()[:20]
    fn = os.path.join("diag", "out", "c0_calib_cache", key + ".json")
    if os.path.exists(fn):
        return json.load(open(fn))
    out = task(*a)
    os.makedirs(os.path.dirname(fn), exist_ok=True)
    json.dump(out, open(fn + ".tmp", "w"))
    os.replace(fn + ".tmp", fn)
    return out


def main():
    w = default_workers(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
    t0 = time.time()
    cal = parallel_map(cached, [(r, s, 1.0, 7.0, 6.0, False) for r in (30.0, 60.0) for s in (100, 101)], w, print,
                       lambda i, r, d, n: f"[calib] {d}/{n} rho {r['rho']} seed {r['seed']} ({time.time() - t0:.0f}s)")
    ratios = [np.mean([v for tt, v in zip(r["t"], r["const"]) if tt >= 3.0]) for r in cal]
    a_norm = float(1.0 / np.mean(ratios))
    print(f"calibration: ratios {np.round(ratios, 4)} -> a = {a_norm:.5f}")
    val = parallel_map(cached, [(r, s, a_norm, 7.0, 6.0, True) for r in (30.0, 45.0, 60.0) for s in (200, 201)], w, print,
                       lambda i, r, d, n: f"[valid] {d}/{n} rho {r['rho']} seed {r['seed']} ({time.time() - t0:.0f}s)")
    for rho in (30.0, 45.0, 60.0):
        R = [r for r in val if r["rho"] == rho]
        late = [i for i, tt in enumerate(R[0]["t"]) if tt >= 3.0]
        cr = [np.mean([r["const"][i] for i in late]) for r in R]
        print(f"validation rho {rho:g}: constant-source ratio {np.mean(cr):.4f} (per seed {np.round(cr, 4)}); "
              + "; ".join(f"k {k}: max|mean| {np.max(np.abs(np.mean([[q[0] for q in r[f'k{k}']] for r in R], 0)[late])):.3f}, "
                          f"sd(t end) {np.mean([r[f'k{k}'][-1][1] for r in R]):.2f}" for k in ("0", "1")))
    fn = os.path.join("diag", "out", f"c0_calib_{time.strftime('%Y%m%dT%H%M%S', time.gmtime())}.json")
    json.dump({"a": a_norm, "cal": cal, "val": val}, open(fn, "w"))
    print("saved", fn)


if __name__ == "__main__":
    main()
