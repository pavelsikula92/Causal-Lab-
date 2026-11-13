"""
C0, development: the minimum-variance kernel of the 3+1 Green function (13 Nov; decides nothing; card C0 NOT frozen)
====================================================================================================================
The author's formulation (13 Nov, 17:41). A pair u < v of the order belongs to the class n = the number of events in
its interval (intrinsic; for a pair near the light cone n is Poisson with mean s = rho V, V = pi tau^4/24). Per unit
measure of the cone, in the continuum limit, the expected number of pairs of class n is
    N_n proportional to g_n = Gamma(n + 1/2)/n!,
and the moments of a kernel w across the cone (in sigma = tau^2 ~ s^1/2) are
    zeroth:  sum_n g_n w_n        (normalisation)
    first:   sum_n w_n            (the inward smearing of the cone)
The 3+1 retarded Green function lies on the cone, so the targets are b0 = Gamma(1/2) (the normalisation of the links,
w = (1), which is Johnston's c = 1/(2 pi sqrt 6 sqrt rho)) and b1 = 0. They come from the continuum Green function,
not from any measurement. The kernel of least variance in the diagonal (independent-contribution) model,
    min sum_n g_n w_n^2  subject to  A w = b,   A = [[g_n], [1]],
is w = N^-1 A^T (A N^-1 A^T)^-1 b, i.e. w_n = l0 + l1 / g_n (n = 0 .. n_max). n_max = 1 gives K2 = (2, -2).
The family n_max = 1, 2, 4, 8, 16 (and the references K1 = links, S16 of the C0 development rows).

Measured on 3+1 pyramids X 7, T 6 (rho 30; seeds 0, 1 and, separately, 2, 3):
  normalisation   the constant source: c sum w / (t^2/2) by time;
  eps(r)          the normalisation resolved by the distance r of the pair (coordinates only as a diagnostic): for
                  sample events at t in [3.5, 5.5], (c sum over pairs in the shell [r, r + dr) of w) / (r dr) - 1;
  noise           the variance over events of the constant-source sum at fixed time, against the diagonal model
                  c^2 sum_n w_n^2 E[count_n] (a test of the off-diagonal covariance);
  dynamics        the linear field (mu 1), data A cos(0) and A cos(x): deviation from the exact solution, mean and sd.
The pairs are enumerated once per sprinkling (n <= 16) and every kernel is evaluated on the same pairs.

    python -m diag.c0_minvar --workers 6
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
from phases.F3_structure import field1d as F

_nb = F._nb
NMAX_ALL = 16
FAMILY = (1, 2, 4, 8, 16)


def g(n):
    return math.gamma(n + 0.5) / math.factorial(n)


def minvar(nmax):
    """w_n (n = 0 .. NMAX_ALL, zero beyond nmax) of least sum g w^2 with sum g w = Gamma(1/2), sum w = 0."""
    gn = np.array([g(n) for n in range(nmax + 1)])
    A = np.vstack([gn, np.ones_like(gn)])
    b = np.array([math.gamma(0.5), 0.0])
    Ninv = np.diag(1.0 / gn)
    w = Ninv @ A.T @ np.linalg.solve(A @ Ninv @ A.T, b)
    out = np.zeros(NMAX_ALL + 1)
    out[:nmax + 1] = w
    return out


def kernels():
    K = {"K1": np.r_[1.0, np.zeros(NMAX_ALL)]}
    for m in FAMILY:
        K[f"MV{m}"] = minvar(m)
    n = np.arange(NMAX_ALL + 1, dtype=float)
    K["S16"] = math.sqrt(math.pi) * 15 / (8 * 4.0) * (1 - 4 * n / 16 + 3 * n * (n - 1) / 256)
    return K


if _nb is not None:
    @_nb.njit(cache=True)
    def enumerate_all(t, x, y, z, tc, kmax, W, phi0, sample, rbin, nrb, c, mu):
        """One pass over the pairs (n <= kmax). Returns the per-event class counts (N x (kmax+1)), the r-histogram of
        the sample events (nrb x (kmax+1)) and the fields of all (kernel, mode) combinations: W is (K x (kmax+1)),
        phi0 is (K*Mo x N) with the rows ordered kernel-major; all fields are linear (J = mu phi)."""
        N = t.size
        K = W.shape[0]
        R = phi0.shape[0]
        Mo = R // K
        cnt = np.zeros((N, kmax + 1), np.int32)
        hist = np.zeros((nrb, kmax + 1))
        phi = np.empty((R, N))
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
                    w = cand[b]
                    dt = t[w] - t[u]
                    if dt > 0.0 and dt * dt > (x[w] - x[u]) ** 2 + (y[w] - y[u]) ** 2 + (z[w] - z[u]) ** 2:
                        nn += 1
                        if nn > kmax:
                            break
                if nn <= kmax:
                    cnt[j, nn] += 1
                    if sample[j]:
                        r = math.sqrt((x[j] - x[u]) ** 2 + (y[j] - y[u]) ** 2 + (z[j] - z[u]) ** 2)
                        ib = int(r / rbin)
                        if ib < nrb:
                            hist[ib, nn] += 1.0
                    for kk in range(K):
                        wg = W[kk, nn]
                        if wg != 0.0:
                            for mo in range(Mo):
                                i = kk * Mo + mo
                                s[i] += wg * mu * phi[i, u]
            for i in range(R):
                phi[i, j] = phi0[i, j] - c * s[i]
        return cnt, hist, phi


def task(rho, seed, X=7.0, T=6.0, rbin=0.25, rmax=3.0):
    t, x, y, z = sprinkle4(rho, X, T, 11000 + seed)
    KS = kernels()
    names = list(KS)
    W = np.array([KS[k] for k in names])
    A = 0.05
    modes = (("wave0", A * np.ones(t.size), A * np.cos(t)),
             ("wave1", A * np.cos(x) * np.cos(t), A * np.cos(x) * np.cos(math.sqrt(2) * t)))
    phi0 = np.array([m[1] for _ in names for m in modes])
    inner = (np.abs(x) < X - t - 0.3) & (np.abs(y) < X - t - 0.3) & (np.abs(z) < X - t - 0.3)
    sample = inner & (t >= 3.5) & (t < 5.5)
    c = 1.0 / (2.0 * math.pi * math.sqrt(6.0 * rho))
    tc = (24.0 * (NMAX_ALL + 6) / (math.pi * rho)) ** 0.25
    a0 = time.time()
    cnt, hist, phi = enumerate_all(t, x, y, z, tc, NMAX_ALL, W, phi0, sample, rbin, int(rmax / rbin), c, 1.0)
    sec = time.time() - a0
    tb = np.arange(0.25, T - 0.4, 0.5)
    slabs = [inner & (np.abs(t - b) < 0.25) for b in tb]
    out = {"rho": rho, "seed": seed, "N": int(t.size), "sec": sec, "t": tb.tolist(), "rbin": rbin,
           "nsample": int(sample.sum()), "mean_counts_late": cnt[inner & (t > 3.0)].mean(0).tolist(), "kernels": {}}
    for kk, name in enumerate(names):
        w = W[kk]
        S = c * (cnt @ w)
        norm = [float(S[m].mean() / (0.5 * np.mean(t[m] ** 2))) for m in slabs]
        var_emp = [float(S[m].var()) for m in slabs]
        var_diag = [float(c * c * (cnt[m].mean(0) @ (w * w))) for m in slabs]
        rr = (np.arange(hist.shape[0]) + 0.5) * rbin
        eps = (c * (hist @ w) / sample.sum() / (rr * rbin) - 1.0).tolist()
        dyn = {}
        for mo, (mname, _, ex) in enumerate(modes):
            dev = (phi[kk * len(modes) + mo] - ex) / A
            dyn[mname] = [[float(dev[m].mean()), float(dev[m].std())] for m in slabs]
        out["kernels"][name] = {"w": w.tolist(), "norm": norm, "var_emp": var_emp, "var_diag": var_diag, "eps_r": eps,
                                "dyn": dyn}
    return out


def cached(*a):
    key = hashlib.sha256(json.dumps([str(q) for q in a] + ["v1"]).encode()).hexdigest()[:20]
    fn = os.path.join("diag", "out", "c0_minvar_cache", key + ".json")
    if os.path.exists(fn):
        return json.load(open(fn))
    out = task(*a)
    os.makedirs(os.path.dirname(fn), exist_ok=True)
    json.dump(out, open(fn + ".tmp", "w"))
    os.replace(fn + ".tmp", fn)
    return out


def summary(res):
    seeds = sorted({r["seed"] for r in res})
    print(f"seeds {seeds}; N {np.mean([r['N'] for r in res]):.0f}; {np.mean([r['sec'] for r in res]):.0f} s per sprinkling")
    tb = np.array(res[0]["t"])
    late = tb >= 3.0
    rr = (np.arange(len(res[0]["kernels"]["K1"]["eps_r"])) + 0.5) * res[0]["rbin"]
    for name in res[0]["kernels"]:
        Ks = [r["kernels"][name] for r in res]
        w = np.array(Ks[0]["w"])
        nz = np.nonzero(w)[0]
        norm = np.mean([np.array(k["norm"])[late].mean() for k in Ks])
        ratio = np.mean([np.array(k["var_emp"])[late] / np.array(k["var_diag"])[late] for k in Ks])
        sdn = np.mean([np.sqrt(np.array(k["var_emp"])[late]).mean() for k in Ks])
        eps = np.mean([k["eps_r"] for k in Ks], 0)
        print(f"== {name}: w = [" + ", ".join(f"{w[i]:+.3f}" for i in nz[:6]) + (" ..." if nz.size > 6 else "") + "]")
        print(f"   normalisation (t >= 3) {norm:.4f}; noise of the constant-source sum (sd) {sdn:.3f}; "
              f"var empirical / diagonal model {ratio:.3f}")
        print("   eps(r) r = " + " ".join(f"{a:.2f}:{b:+.3f}" for a, b in zip(rr, eps) if a < 2.6))
        for mname in ("wave0", "wave1"):
            D = np.array([k["dyn"][mname] for k in Ks])
            print(f"   {mname}: mean " + " ".join(f"{v:+.2f}" for v in D[:, :, 0].mean(0)[1::2])
                  + "  sd " + " ".join(f"{v:.2f}" for v in D[:, :, 1].mean(0)[1::2]) + "  (t = "
                  + " ".join(f"{v:.2f}" for v in tb[1::2]) + ")")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=0)
    ap.add_argument("--rho", type=float, default=30.0)
    ap.add_argument("--seeds", type=str, default="0,1,2,3")
    a = ap.parse_args()
    print("kernels (n_max: weights):")
    for name, w in kernels().items():
        nz = np.nonzero(w)[0]
        print(f"  {name}: " + ", ".join(f"{w[i]:+.4f}" for i in nz))
    jobs = [(a.rho, int(s)) for s in a.seeds.split(",")]
    t0 = time.time()
    res = parallel_map(cached, jobs, default_workers(a.workers), print,
                       lambda i, r, d, n: f"[C0 minvar] {d}/{n}: rho {r['rho']} seed {r['seed']} N {r['N']} "
                                          f"{r['sec']:.0f}s ({time.time() - t0:.0f}s)")
    fn = os.path.join("diag", "out", f"c0_minvar_{time.strftime('%Y%m%dT%H%M%S', time.gmtime())}.json")
    json.dump(res, open(fn, "w"))
    print("saved", fn)
    for grp, name in (([0, 1], "derivation seeds 0, 1"), ([2, 3], "validation seeds 2, 3")):
        R = [r for r in res if r["seed"] in grp]
        if R:
            print(f"######## {name}")
            summary(R)


if __name__ == "__main__":
    main()
