"""
C0, development: the counterterm kernels on the nonlinear response, the higher modes and the long-time noise (13 Nov,
decision 25; decides nothing; card C0 NOT frozen)
=====================================================================================================================
Decision 25: the apex term of a layer kernel is removed by the local counterterm,
    phi(v) - delta J(v) = phi0(v) - c sum_{u < v} w(n) J(u),   delta = (a / (8 sqrt pi)) sum_n w_n H_n   (diag/c0_apex.py),
solved per event (Newton; 1 - delta mu cos phi > 0, unique). The candidates MV2+ct and MV8+ct; references MV2 and MV8
without the counterterm and MVA8, MVA16 (the constraint sum w H_n = 0).

Before the card C0 (decision 25) three things are measured, on 3+1 pyramids X 9.5, T 8 (rho 20 and 30):
  higher modes     linear (J = mu phi, mu 1), data A cos(k x) cos(k t), k = 0, 1, 2, 3, against A cos(k x)
                   cos(sqrt(k^2 + mu) t); the expectation (A3) of the kernel predicts the mean (diag/c0_apex.dyn_pred);
  nonlinear        J = mu sin phi, homogeneous data of amplitude A = 0.5, 1, 2 against the pendulum, and A = 1 cos x
                   against the 1+1 sine-Gordon solution with the same data (a finite-difference solution, periodic in x);
                   the expectation of the kernel without noise feedback is predicted for the homogeneous modes
                   (diag/c0_apex.dyn_pred_nl0); the renormalization of sin by the noise is measured;
  long-time noise  the sd of (phi - exact)/A over events by time, and its growth exponent p (sd ~ t^p over t 1.5 .. 7.5)
                   for every mode and kernel, at two densities.
Every reading is the projection on the mode, cos(k x) (k = 0: the mean), and the sd of the remainder, in slabs
|t - b| < 0.25 of the inner events |x|, |y|, |z| < X - t - 0.3. Coordinates only as a diagnostic.

    python -m diag.c0_ct --predict --rhos 20,30          (the predictions alone; written before the run)
    python -m diag.c0_ct --rhos 20,30 --seeds 0,1 --workers 4
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
from diag.c0_apex import delta, dyn_pred, dyn_pred_nl0, pendulum
from diag.c0_links4 import sprinkle4
from diag.c0_minvar import NMAX_ALL, minvar
from phases.F3_structure import field1d as F

_nb = F._nb
X0, T0 = 9.5, 8.0
MODES = (("lin0", True, 0, 0.05), ("lin1", True, 1, 0.05), ("lin2", True, 2, 0.05), ("lin3", True, 3, 0.05),
         ("nl0_A0.5", False, 0, 0.5), ("nl0_A1", False, 0, 1.0), ("nl0_A2", False, 0, 2.0), ("nl1_A1", False, 1, 1.0))


def rows(rho):
    """(name, w, ct) in the order of the run; ct = delta of (A2) at rho for the +ct rows, 0 otherwise."""
    out = []
    for m in (2, 8):
        w = minvar(m)
        out.append((f"MV{m}+ct", w, delta(w, rho)))
    for m in (2, 8):
        out.append((f"MV{m}", minvar(m), 0.0))
    for m in (8, 16):
        out.append((f"MVA{m}", minvar(m, apex=True), 0.0))
    return out


if _nb is not None:
    @_nb.njit(cache=True)
    def evolve_ct(t, x, y, z, tc, kmax, W, ctd, phi0, lin, c, mu):
        """All (kernel, mode) rows in one pass over the pairs (n <= kmax); rows kernel-major, Mo = len(lin) modes per
        kernel. phi - ctd[kernel] J(phi) = phi0 - c sum w(n) J(u); J = mu phi (lin) or mu sin phi."""
        N = t.size
        K = W.shape[0]
        R = phi0.shape[0]
        Mo = R // K
        phi = np.empty((R, N))
        J = np.empty((R, N))
        s = np.zeros(R)
        tc2 = tc * tc
        cand = np.empty(N, np.int64)
        for j in range(N):
            m = 0
            for k in range(j - 1, -1, -1):
                dt = t[j] - t[k]
                if dt <= 0.0:
                    continue
                q = dt * dt - (x[j] - x[k]) ** 2 - (y[j] - y[k]) ** 2 - (z[j] - z[k]) ** 2
                if q > 0.0 and q <= tc2:
                    cand[m] = k
                    m += 1
            for i in range(R):
                s[i] = 0.0
            for a in range(m):
                u = cand[a]
                nn = 0
                for b in range(a):
                    wv = cand[b]
                    dt = t[wv] - t[u]
                    if dt > 0.0 and dt * dt > (x[wv] - x[u]) ** 2 + (y[wv] - y[u]) ** 2 + (z[wv] - z[u]) ** 2:
                        nn += 1
                        if nn > kmax:
                            break
                if nn <= kmax:
                    for kk in range(K):
                        wg = W[kk, nn]
                        if wg != 0.0:
                            base = kk * Mo
                            for mo in range(Mo):
                                s[base + mo] += wg * J[base + mo, u]
            for i in range(R):
                rhs = phi0[i, j] - c * s[i]
                d = ctd[i // Mo] * mu
                if lin[i % Mo]:
                    p = rhs / (1.0 - d)
                    phi[i, j] = p
                    J[i, j] = mu * p
                else:
                    p = rhs
                    for _ in range(60):
                        dp = (p - d * math.sin(p) - rhs) / (1.0 - d * math.cos(p))
                        p -= dp
                        if abs(dp) < 1e-13:
                            break
                    phi[i, j] = p
                    J[i, j] = mu * math.sin(p)
        return phi


def sg1d(A, k, T, mu=1.0, nx=512, dt=0.002, every=5):
    """phi_tt - phi_xx + mu sin phi = 0 on x in [0, 2 pi / k) periodic, phi(0) = A cos(k x), phi_t(0) = 0 (leapfrog).
    Returns (times, xs, field) sampled every `every` steps."""
    L = 2 * math.pi / k
    xs = np.arange(nx) * L / nx
    dx = L / nx
    p0 = A * np.cos(k * xs)
    lap = lambda p: (np.roll(p, 1) - 2 * p + np.roll(p, -1)) / (dx * dx)
    p1 = p0 + 0.5 * dt * dt * (lap(p0) - mu * np.sin(p0))
    nt = int(round(T / dt)) + 1
    ts, fs = [0.0], [p0.copy()]
    for n in range(1, nt):
        if n % every == 0:
            ts.append(n * dt)
            fs.append(p1.copy())
        p2 = 2 * p1 - p0 + dt * dt * (lap(p1) - mu * np.sin(p1))
        p0, p1 = p1, p2
    return np.array(ts), xs, np.array(fs)


def exact_field(mode, t, x, T, mu=1.0):
    name, lin, k, A = mode
    if lin:
        return A * np.cos(k * x) * np.cos(math.sqrt(k * k + mu) * t)
    if k == 0:
        tg = np.arange(0.0, T + 0.01, 0.001)
        return np.interp(t, tg, pendulum(A, tg, mu))
    ts, xs, f = sg1d(A, k, T + 0.05, mu)
    L = 2 * math.pi / k
    nx = xs.size
    fx = (x % L) / (L / nx)
    i0 = np.floor(fx).astype(int) % nx
    i1 = (i0 + 1) % nx
    ax = fx - np.floor(fx)
    ft = t / (ts[1] - ts[0])
    j0 = np.clip(np.floor(ft).astype(int), 0, ts.size - 2)
    at = ft - j0
    v0 = f[j0, i0] * (1 - ax) + f[j0, i1] * ax
    v1 = f[j0 + 1, i0] * (1 - ax) + f[j0 + 1, i1] * ax
    return v0 * (1 - at) + v1 * at


def task(rho, seed, X=X0, T=T0):
    t, x, y, z = sprinkle4(rho, X, T, 12000 + seed)
    RS = rows(rho)
    W = np.array([r[1] for r in RS])
    ctd = np.array([r[2] for r in RS])
    phi0 = np.array([m[3] * np.cos(m[2] * x) * np.cos(m[2] * t) for _ in RS for m in MODES])
    lin = np.array([m[1] for m in MODES])
    c = 1.0 / (2.0 * math.pi * math.sqrt(6.0 * rho))
    tc = (24.0 * (NMAX_ALL + 6) / (math.pi * rho)) ** 0.25
    a0 = time.time()
    phi = evolve_ct(t, x, y, z, tc, NMAX_ALL, W, ctd, phi0, lin, c, 1.0)
    sec = time.time() - a0
    inner = (np.abs(x) < X - t - 0.3) & (np.abs(y) < X - t - 0.3) & (np.abs(z) < X - t - 0.3)
    tb = np.arange(0.25, T - 0.4, 0.5)
    slabs = [inner & (np.abs(t - b) < 0.25) for b in tb]
    ex = {m[0]: exact_field(m, t, x, T) for m in MODES}
    out = {"rho": rho, "seed": seed, "N": int(t.size), "sec": sec, "t": tb.tolist(), "X": X, "T": T,
           "nslab": [int(m.sum()) for m in slabs], "rows": {}}
    Mo = len(MODES)
    for kk, (name, w, ct) in enumerate(RS):
        R = {"w": w.tolist(), "ct": ct, "modes": {}}
        for mo, m in enumerate(MODES):
            dev = (phi[kk * Mo + mo] - ex[m[0]]) / m[3]
            bas = np.cos(m[2] * x)
            dd = []
            for sl in slabs:
                amp = float(dev[sl] @ bas[sl] / (bas[sl] @ bas[sl]))
                dd.append([amp, float((dev[sl] - amp * bas[sl]).std())])
            R["modes"][m[0]] = dd
        out["rows"][name] = R
    return out


def cached(*a):
    key = hashlib.sha256(json.dumps([str(q) for q in a] + ["ct-v1"]).encode()).hexdigest()[:20]
    fn = os.path.join("diag", "out", "c0_ct_cache", key + ".json")
    if os.path.exists(fn):
        return json.load(open(fn))
    out = task(*a)
    os.makedirs(os.path.dirname(fn), exist_ok=True)
    json.dump(out, open(fn + ".tmp", "w"))
    os.replace(fn + ".tmp", fn)
    return out


_PRED = {}


def prediction(name, w, ct, rho, mode, tb, X=X0):
    """The expectation of the kernel (no noise feedback); None where not computed (nl1)."""
    key = (name, rho, mode[0])
    if key not in _PRED:
        mname, lin, k, A = mode
        if lin:
            _PRED[key] = dyn_pred(w, rho, k, tb, X=X, ct=ct)
        elif k == 0:
            _PRED[key] = dyn_pred_nl0(w, rho, A, tb, X=X, ct=ct)
        else:
            _PRED[key] = None
    return _PRED[key]


def growth(tb, sd, lo=1.5, hi=7.5):
    m = (tb >= lo) & (tb <= hi) & (sd > 0)
    return float(np.polyfit(np.log(tb[m]), np.log(sd[m]), 1)[0])


def summary(res, predict=True):
    rho = res[0]["rho"]
    tb = np.array(res[0]["t"])
    sel = slice(1, None, 2)
    print(f"rho {rho:g}; seeds {sorted(r['seed'] for r in res)}; N {np.mean([r['N'] for r in res]):.0f}; "
          f"{np.mean([r['sec'] for r in res]):.0f} s per sprinkling; events per slab (t = "
          + " ".join(f"{v:.2f}" for v in tb[sel]) + "): " + " ".join(str(v) for v in np.array(res[0]["nslab"])[sel]))
    for name in res[0]["rows"]:
        Rs = [r["rows"][name] for r in res]
        ct = Rs[0]["ct"]
        w = np.array(Rs[0]["w"])
        print(f"== {name}" + (f" (counterterm delta {ct:+.4f})" if ct else ""))
        for m in MODES:
            D = np.array([R["modes"][m[0]] for R in Rs])
            amp = D[:, :, 0].mean(0)
            sd = D[:, :, 1].mean(0)
            P = prediction(name, w, ct, rho, m, tb) if predict else None
            print(f"   {m[0]:9s} amp-dev " + " ".join(f"{v:+.2f}" for v in amp[sel])
                  + ((" (pred " + " ".join(f"{v:+.2f}" for v in np.array(P)[sel]) + ")") if P is not None else "")
                  + "  sd " + " ".join(f"{v:.2f}" for v in sd[sel]) + f"  p {growth(tb, sd):.2f}")


def predict_only(rho):
    tb = np.arange(0.25, T0 - 0.4, 0.5)
    sel = slice(1, None, 2)
    print(f"######## predictions, rho {rho:g} (t = " + " ".join(f"{v:.2f}" for v in tb[sel]) + ")")
    for name, w, ct in rows(rho):
        print(f"== {name}" + (f" (counterterm delta {ct:+.4f})" if ct else ""))
        for m in MODES:
            P = prediction(name, w, ct, rho, m, tb)
            if P is not None:
                print(f"   {m[0]:9s} amp-dev " + " ".join(f"{v:+.3f}" for v in np.array(P)[sel]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=0)
    ap.add_argument("--rhos", type=str, default="20,30")
    ap.add_argument("--seeds", type=str, default="0,1")
    ap.add_argument("--predict", action="store_true")
    ap.add_argument("--no-pred", action="store_true")
    a = ap.parse_args()
    rhos = [float(r) for r in a.rhos.split(",")]
    if a.predict:
        for rho in rhos:
            predict_only(rho)
        return
    jobs = [(rho, int(s)) for rho in sorted(rhos, reverse=True) for s in a.seeds.split(",")]
    t0 = time.time()
    res = parallel_map(cached, jobs, default_workers(a.workers), print,
                       lambda i, r, d, n: f"[C0 ct] {d}/{n}: rho {r['rho']} seed {r['seed']} N {r['N']} "
                                          f"{r['sec']:.0f}s ({time.time() - t0:.0f}s)")
    fn = os.path.join("diag", "out", f"c0_ct_{time.strftime('%Y%m%dT%H%M%S', time.gmtime())}.json")
    json.dump(res, open(fn, "w"))
    print("saved", fn)
    for rho in sorted(rhos):
        R = [r for r in res if r["rho"] == rho]
        print(f"######## rho {rho:g}")
        summary(R, not a.no_pred)


if __name__ == "__main__":
    main()
