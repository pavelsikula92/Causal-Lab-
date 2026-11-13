"""
C0, development: long-time stability of the 3+1 kernels K2 and S16 (13 Nov; decides nothing; the card C0 is NOT frozen)
=====================================================================================================================
The author's instructions (13 Nov, 14:36): separate the accumulation of numerical noise from a real instability.
Per kernel and density, on the 3+1 pyramid X 9.5, T 8 (inner region |x|, |y|, |z| < X - t - 0.3):
  wave0    linear rule (J = mu phi), data A cos(0), A 0.05: deviation from A cos(t)
  wave1    linear rule, data A cos(x): deviation from A cos(x) cos(sqrt 2 t)
  probe    linear rule, phi0 = eps xi (independent N(0, 1) per event; no physical data): the response of the discrete
           operator to white noise; a growing mode shows as exponential growth of the response
  nl<A>    the nonlinear rule (J = sin phi), data A cos(0), A 0.05, 0.5, 1.0 (the exact solution: the pendulum
           equation phi'' = -sin phi): does the growth depend on the amplitude?
Per slab of width 0.25: mean and sd of the deviation, its projection on plane waves with |k| = 0, 1, 2, 4 (cos and
sin along x, y, z; least squares), and the correlation of the deviation with that of the previous slab at the nearest
event (time correlation).
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
from diag.c0_dev import kernel


def pendulum(A, T, dt=1e-3):
    n = int(T / dt) + 2
    ph = np.empty(n)
    p, v = A, 0.0
    for i in range(n):
        ph[i] = p
        a = -math.sin(p)
        v += dt * a
        p += dt * v
    return lambda t: np.interp(t, np.arange(n) * dt, ph)


MODES = ("wave0", "wave1", "probe", "nl0.05", "nl0.5", "nl1.0")


def task(rho, kname, seed, X=9.5, T=8.0):
    """All modes on one sprinkling (the pairs are computed once)."""
    t, x, y, z = sprinkle4(rho, X, T, 5000 + seed)
    coef, lam = kernel(kname)
    a0 = time.time()
    ptr, idx, lay = layers4(t, x, y, z, (24 * lam / (math.pi * rho)) ** 0.25, coef.size - 1)
    sec_l = time.time() - a0
    out = {"rho": rho, "kernel": kname, "seed": seed, "N": int(t.size), "pairs": int(idx.size), "sec_layers": sec_l,
           "modes": {}}
    for mode in MODES:
        out["modes"][mode] = _mode(t, x, y, z, ptr, idx, lay, coef, rho, mode, seed, X, T)
    out["sec"] = time.time() - a0
    return out


def _mode(t, x, y, z, ptr, idx, lay, coef, rho, mode, seed, X, T):
    c = 1 / (2 * math.pi * math.sqrt(6 * rho))
    nopre = np.zeros(t.size, np.bool_)
    rng = np.random.default_rng(9000 + seed)
    A = 0.05
    if mode == "wave0":
        phi0, ex, lin = A * np.ones(t.size), A * np.cos(t), True
    elif mode == "wave1":
        phi0, ex, lin = A * np.cos(x) * np.cos(t), A * np.cos(x) * np.cos(math.sqrt(2) * t), True
    elif mode == "probe":
        A = 0.05
        phi0, ex, lin = A * rng.normal(size=t.size), np.zeros(t.size), True
    else:
        A = float(mode[2:])
        P = pendulum(A, T)
        phi0, ex, lin = A * np.ones(t.size), P(t), False
    phi = evolve_layers(ptr, idx, lay, coef, phi0, nopre, ex, c, 1.0, lin)
    dev = (phi - ex) / A
    inner = (np.abs(x) < X - t - 0.3) & (np.abs(y) < X - t - 0.3) & (np.abs(z) < X - t - 0.3)
    rows = []
    prev = None
    for b in np.arange(0.125, T - 0.3, 0.25):
        m = np.nonzero(inner & (np.abs(t - b) < 0.125))[0]
        if m.size < 50:
            break
        d = dev[m]
        cols = [np.ones(m.size)]
        names = ["k0"]
        for kk in (1.0, 2.0, 4.0):
            for q in (x, y, z):
                cols += [np.cos(kk * q[m]), np.sin(kk * q[m])]
                names += [f"k{kk:g}", f"k{kk:g}"]
        M = np.column_stack(cols)
        coefs, *_ = np.linalg.lstsq(M, d, rcond=None)
        power = {}
        for nm, cv in zip(names, coefs):
            power[nm] = power.get(nm, 0.0) + float(cv * cv)
        resid = float(np.std(d - M @ coefs))
        corr = float("nan")
        if prev is not None:
            pm, pd = prev
            # nearest event of the previous slab in space (coarse: sorted by x)
            o = np.argsort(x[pm])
            j = np.clip(np.searchsorted(x[pm][o], x[m]), 0, o.size - 1)
            nb = pm[o][j]
            dd = dev[nb]
            if np.std(dd) > 0 and np.std(d) > 0:
                corr = float(np.corrcoef(d, dd)[0, 1])
        rows.append({"t": float(b), "mean": float(d.mean()), "sd": float(d.std()), "power": power, "resid": resid,
                     "corr": corr, "n": int(m.size)})
        prev = (m, d)
    return rows


def cached(*a):
    key = hashlib.sha256(json.dumps([str(q) for q in a] + ["v2"]).encode()).hexdigest()[:20]
    fn = os.path.join("diag", "out", "c0_stab_cache", key + ".json")
    if os.path.exists(fn):
        return json.load(open(fn))
    out = task(*a)
    os.makedirs(os.path.dirname(fn), exist_ok=True)
    json.dump(out, open(fn + ".tmp", "w"))
    os.replace(fn + ".tmp", fn)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=0)
    ap.add_argument("--seeds", type=int, default=2)
    ap.add_argument("--rhos", type=str, default="30")
    ap.add_argument("--kernels", type=str, default="K2,S16")
    a = ap.parse_args()
    jobs = [(float(r), kn, s) for s in range(a.seeds) for r in a.rhos.split(",") for kn in a.kernels.split(",")]
    t0 = time.time()
    res = parallel_map(cached, jobs, default_workers(a.workers), print,
                       lambda i, r, d, n: f"[C0 stab] {d}/{n}: {r['kernel']} rho {r['rho']} seed {r['seed']} "
                                          f"N {r['N']} {r['sec']:.0f}s ({time.time() - t0:.0f}s)")
    fn = os.path.join("diag", "out", f"c0_stab_{time.strftime('%Y%m%dT%H%M%S', time.gmtime())}.json")
    json.dump(res, open(fn, "w"))
    print("saved", fn)


if __name__ == "__main__":
    main()
