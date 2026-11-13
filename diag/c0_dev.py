"""
C0 development rows: kernels for the 3+1 Green function of the order (decides nothing)
=====================================================================================
Kernels (per pair u < v with n events in its interval; units of c = 1/(2 pi sqrt(6) sqrt(rho)); all intrinsic):
  K1      w(0) = 1                                   (links, Johnston)
  K2      w(0) = 2, w(1) = -2                        (the first moment across the cone removed)
  S<s0>   w(n) = sqrt(pi) 15/(8 sqrt(s0)) [1 - 4 n/s0 + 3 n(n-1)/s0^2], n <= s0
          (a smooth kernel of the expected count s = rho V: zeroth moment of links, first moment zero; falling
           factorials so that E[w(n) | s] is the smooth profile)
Metrics on the inner region of a 3+1 pyramid (X 7, T 6):
  constant source: c sum w / (t^2/2) by time;
  linear field (J = mu phi, mu 1), data A cos(k x), phi_t = 0, k = 0, 1, 2: (phi - exact)/A mean and sd by time.

    python -m diag.c0_dev --workers 2
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import time

import numpy as np

from core.stats.parallel import parallel_map, default_workers
from diag.c0_links4 import sprinkle4
from diag.c0_layers4 import layers4, evolve_layers


def kernel(name):
    if name == "K1":
        return np.array([1.0]), 12.0
    if name == "K2":
        return np.array([2.0, -2.0]), 16.0
    s0 = float(name[1:])
    n = np.arange(int(s0) + 1, dtype=float)
    return math.sqrt(math.pi) * 15 / (8 * math.sqrt(s0)) * (1 - 4 * n / s0 + 3 * n * (n - 1) / s0 ** 2), s0 + 6.0


def task(rho, kname, seed, X=7.0, T=6.0):
    t, x, y, z = sprinkle4(rho, X, T, 1000 + seed)
    coef, lam = kernel(kname)
    a = time.time()
    ptr, idx, lay = layers4(t, x, y, z, (24 * lam / (math.pi * rho)) ** 0.25, coef.size - 1)
    sec = time.time() - a
    c = 1 / (2 * math.pi * math.sqrt(6 * rho))
    inner = (np.abs(x) < X - t - 0.3) & (np.abs(y) < X - t - 0.3) & (np.abs(z) < X - t - 0.3)
    nopre = np.zeros(t.size, np.bool_)
    tb = np.arange(0.25, T - 0.4, 0.5)
    tot = np.zeros(t.size)
    for v in range(t.size):
        tot[v] = c * np.sum(coef[lay[ptr[v]:ptr[v + 1]]])
    const = [float(tot[m].mean() / (0.5 * np.mean(t[m] ** 2))) for m in [inner & (np.abs(t - b) < 0.25) for b in tb]]
    waves = {}
    for k in (0.0, 1.0, 2.0):
        A = 0.05
        om = math.sqrt(k * k + 1)
        phi0 = A * np.cos(k * x) * np.cos(k * t)
        ex = A * np.cos(k * x) * np.cos(om * t)
        phi = evolve_layers(ptr, idx, lay, coef, phi0, nopre, ex, c, 1.0, True)
        rows = []
        for b in tb:
            m = inner & (np.abs(t - b) < 0.25)
            d = (phi[m] - ex[m]) / A
            rows.append([float(d.mean()), float(d.std()), int(m.sum())])
        waves[str(k)] = rows
    return {"rho": rho, "kernel": kname, "seed": seed, "N": int(t.size), "pairs": int(idx.size), "sec": sec,
            "t": tb.tolist(), "const": const, "waves": waves}


def cached(*a):
    key = hashlib.sha256(json.dumps([str(q) for q in a] + ["v1"]).encode()).hexdigest()[:20]
    fn = os.path.join("diag", "out", "c0_dev_cache", key + ".json")
    if os.path.exists(fn):
        return json.load(open(fn))
    out = task(*a)
    os.makedirs(os.path.dirname(fn), exist_ok=True)
    json.dump(out, open(fn + ".tmp", "w"))
    os.replace(fn + ".tmp", fn)
    return out


def summary(res):
    for kn in ("K1", "K2", "S8", "S12", "S16"):
        for rho in sorted({r["rho"] for r in res}):
            R = [r for r in res if r["kernel"] == kn and r["rho"] == rho]
            if not R:
                continue
            tb = R[0]["t"]
            cs = np.mean([r["const"] for r in R], 0)
            print(f"== {kn} rho {rho:g} (n {len(R)}, pairs/event {np.mean([r['pairs'] / r['N'] for r in R]):.0f}, "
                  f"{np.mean([r['sec'] for r in R]):.0f}s): constant " + " ".join(f"{a:.2f}:{b:.3f}" for a, b in zip(tb, cs)))
            for k in ("0.0", "1.0", "2.0"):
                W = np.array([r["waves"][k] for r in R])
                print(f"   k {k}: " + " ".join(f"{a:.2f}:{m:+.3f}({s:.2f})" for a, m, s in
                                             zip(tb, W[:, :, 0].mean(0), W[:, :, 1].mean(0))))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=0)
    ap.add_argument("--seeds", type=int, default=2)
    a = ap.parse_args()
    jobs = [(rho, kn, s) for s in range(a.seeds) for rho in (30.0, 60.0) for kn in ("K1", "K2", "S8", "S12", "S16")]
    t0 = time.time()
    res = parallel_map(cached, jobs, default_workers(a.workers), print,
                       lambda i, r, d, n: f"[C0 dev] {d}/{n}: {r['kernel']} rho {r['rho']} seed {r['seed']} N {r['N']} "
                                          f"{r['sec']:.0f}s ({time.time() - t0:.0f}s)")
    fn = os.path.join("diag", "out", f"c0_dev_{time.strftime('%Y%m%dT%H%M%S', time.gmtime())}.json")
    json.dump(res, open(fn, "w"))
    print("saved", fn)
    summary(res)


if __name__ == "__main__":
    main()
