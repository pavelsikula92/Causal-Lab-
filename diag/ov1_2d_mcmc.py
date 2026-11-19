"""
OV1 phase control in 2D (decides nothing): an F6.4-type sampler of 2D orders with the LOCAL moves of family I.
=============================================================================================================
The author (19 Nov): the escape from the manifold-like phase is probed by an F6.4-type configuration-space sampler,
not by OV3. Proposals are the adjacent transpositions of family I (docs/derivations/order_variation_OV1.md, 3.2); the
acceptance and sampling rule of F6.4a is kept unchanged:
- a proposal picks U or V with probability 1/2, then a position p uniformly in 0..N-2, and swaps the elements at p
  and p + 1 in that order (F6.4a swaps two elements chosen uniformly; that is the only difference);
- Metropolis: accept if S' <= S or with probability exp(-beta (S' - S));
- one sweep = N (N - 1)/2 proposals; one sample per sweep after the thermalisation.
The proposal is symmetric and the set of adjacent transpositions connects all 2D orders, so the stationary
distribution is the same Boltzmann distribution as in F6.4a; only the mixing differs. This is an exploration of the
configuration space, not a physical time evolution.

Action: the F6.4a convention S = 4 eps (N - 2 eps sum_n N_n f(n, eps)); Delta S of a move is computed locally
(diag/ov1_2d.py, checked against brute force). Observables per sample: S, chi (ordering fraction), d_MM, h (links
of the longest chain), N_0 (links). Per chain: means, var(S), acceptance, tau_int(S) (F6.4a definition).

Instrument checks (as F6.4a): (i) N = 4, eps = 0.5, beta = 1: the distribution of S from 300 000 proposals against
the exact Boltzmann distribution over all 2D orders (total variation <= 0.01); (ii) beta = 0: one chain against
20 000 directly sampled uniform 2D orders (|dS| <= 0.2, |dchi| <= 0.005, |dh| <= 0.15). The local Delta S is also
compared with a full recomputation along a chain.

    python -m diag.ov1_2d_mcmc --check
    python -m diag.ov1_2d_mcmc --scan --N 50 --eps 0.12 --chains 4 --therm 500 --meas 3000 --workers 8
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
import os
import time

import numpy as np

from diag.ov1_2d import fker, n_cut_of, reversed_order
from diag.ov1_kernels import _seed, abundances, rel_count_height, run_sweeps

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "ov1_2d_mcmc")
BETAS_F64A = [0.0, 0.1, 0.5, 1.0, 1.5] + [round(1.8 + 0.1 * k, 1) for k in range(11)] + [3.0, 3.5, 4.0, 4.5]


def _reps(vr):
    N = vr.size
    iv = np.empty(N, np.int64)
    iv[vr] = np.arange(N)
    vrR, ivR = reversed_order(vr)
    tvR, tivR = reversed_order(iv.copy())
    return (vr.copy(), iv.copy(), vrR, ivR), (iv.copy(), vr.copy(), tvR, tivR)


def exact_S(vr, f, ncut, eps):
    N = vr.size
    iv = np.empty(N, np.int64)
    iv[vr] = np.arange(N)
    H = abundances(vr, iv, ncut, np.zeros(N, np.int64), 1)[0]
    return 4 * eps * (N - 2 * eps * float(np.dot(H, f[:H.size]))), int(H[0])


def tau_int(x):
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


def d_mm(chi):
    g = lambda d: math.exp(math.lgamma(d + 1) + math.lgamma(d / 2) - math.log(2.0) - math.lgamma(1.5 * d))
    if not 0 < chi < 1:
        return float("nan")
    lo, hi = 0.5, 12.0
    for _ in range(80):
        m = 0.5 * (lo + hi)
        lo, hi = (m, hi) if g(m) > chi else (lo, m)
    return 0.5 * (lo + hi)


def chain(N, eps, beta, therm, meas, seed, start=None):
    rng = np.random.default_rng(seed)
    _seed(int(rng.integers(2 ** 31)))
    nc = min(n_cut_of(eps), N)
    f = fker(eps, max(nc, N) + 2)
    vr = rng.permutation(N).astype(np.int64) if start is None else np.asarray(start, np.int64).copy()
    A, B = _reps(vr)
    S, _ = exact_S(A[0], f, nc, eps)
    R, _ = rel_count_height(A[0])
    moves = N * (N - 1) // 2
    acc = 0
    Ss, chis, hs, links, drift = [], [], [], [], 0.0
    t0 = time.time()
    for sw in range(therm + meas):
        S, R, a = run_sweeps(A, B, S, R, beta, f, nc, eps, moves)
        acc += a
        if sw >= therm:
            Sx, N0 = exact_S(A[0], f, nc, eps)
            drift = max(drift, abs(Sx - S))
            S = Sx
            Rx, Lc = rel_count_height(A[0])
            assert Rx == R
            Ss.append(S)
            chis.append(R / (N * (N - 1) / 2))
            hs.append(Lc - 1)
            links.append(N0)
    Ss = np.array(Ss)
    return {"N": N, "eps": eps, "beta": beta, "seed": seed, "S": float(Ss.mean()), "varS": float(Ss.var()),
            "chi": float(np.mean(chis)), "d_mm": d_mm(float(np.mean(chis))), "h": float(np.mean(hs)),
            "links": float(np.mean(links)), "acc": acc / ((therm + meas) * moves), "tau": tau_int(Ss),
            "max_dS_drift": drift, "sec": time.time() - t0, "final_vr": A[0].tolist()}


def check():
    ok = True
    # (i) exact Boltzmann distribution of S, N = 4
    N, eps, beta = 4, 0.5, 1.0
    f = fker(eps, N + 2)
    vals = {}
    for p in itertools.permutations(range(N)):
        s = round(exact_S(np.array(p, np.int64), f, N, eps)[0], 9)
        vals[s] = vals.get(s, 0.0) + math.exp(-beta * s)
    Z = sum(vals.values())
    _seed(7)
    A, B = _reps(np.random.default_rng(7).permutation(N).astype(np.int64))
    S, _ = exact_S(A[0], f, N, eps)
    R, _ = rel_count_height(A[0])
    hist = {}
    for _ in range(300000):
        S, R, _a = run_sweeps(A, B, S, R, beta, f, N, eps, 1)
        k = round(S, 9)
        hist[k] = hist.get(k, 0) + 1
    tv = 0.5 * sum(abs(hist.get(k, 0) / 300000 - vals.get(k, 0) / Z) for k in set(hist) | set(vals))
    print(f"(i) N 4, eps 0.5, beta 1: TV {tv:.4f} over {len(vals)} action values (<= 0.01)")
    ok &= tv <= 0.01
    # local Delta S against recomputation along a chain at N 300, beta 2
    N, eps = 300, 0.12
    nc = n_cut_of(eps)
    f = fker(eps, nc + 2)
    _seed(3)
    A, B = _reps(np.random.default_rng(3).permutation(N).astype(np.int64))
    S, _ = exact_S(A[0], f, nc, eps)
    R, _ = rel_count_height(A[0])
    worst = 0.0
    for _ in range(200):
        S, R, _a = run_sweeps(A, B, S, R, 2.0, f, nc, eps, 50)
        Sx, _ = exact_S(A[0], f, nc, eps)
        worst = max(worst, abs(Sx - S))
        assert rel_count_height(A[0])[0] == R and np.all(A[1] == B[0]) and np.all(A[0] == B[1])
    print(f"Delta S accumulated over 10 000 proposals (N 300, beta 2) vs recomputation: max {worst:.1e} (<= 1e-8)")
    ok &= worst <= 1e-8
    # (ii) beta = 0 against direct uniform sampling, N 50, eps 0.12
    N, eps = 50, 0.12
    f = fker(eps, N + 2)
    rng = np.random.default_rng(11)
    D = []
    for _ in range(20000):
        vr = rng.permutation(N).astype(np.int64)
        s, _ = exact_S(vr, f, N, eps)
        r, h = rel_count_height(vr)
        D.append((s, r / (N * (N - 1) / 2), h - 1))
    D = np.array(D).mean(0)
    c = chain(N, eps, 0.0, 100, 1000, 12)
    d = (c["S"] - D[0], c["chi"] - D[1], c["h"] - D[2])
    print(f"(ii) beta 0: chain S {c['S']:.3f} chi {c['chi']:.4f} h {c['h']:.3f}; direct S {D[0]:.3f} chi {D[1]:.4f} "
          f"h {D[2]:.3f}; |dS| {abs(d[0]):.3f} (<= 0.2) |dchi| {abs(d[1]):.4f} (<= 0.005) |dh| {abs(d[2]):.3f} (<= 0.15)")
    ok &= abs(d[0]) <= 0.2 and abs(d[1]) <= 0.005 and abs(d[2]) <= 0.15
    print("checks:", "OK" if ok else "FAILED")


def _job(a):
    N, eps, beta, therm, meas, seed = a
    f = os.path.join(OUT, f"mc_{N}_{eps:g}_{beta:g}_{seed}.json")
    if os.path.exists(f):
        return json.load(open(f))
    r = chain(N, eps, beta, therm, meas, seed)
    json.dump(r, open(f + ".tmp", "w"))
    os.replace(f + ".tmp", f)
    return r


def scan(N, eps, chains, therm, meas, workers, betas):
    os.makedirs(OUT, exist_ok=True)
    jobs = [(N, eps, b, therm, meas, 7919 * c + 104729 * int(round(10 * b)) + N) for b in betas for c in range(chains)]
    from concurrent.futures import ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=workers) as ex:
        res = list(ex.map(_job, jobs))
    print(f"N {N}, eps {eps}, {chains} chains x ({therm} + {meas}) sweeps of {N * (N - 1) // 2} local proposals")
    print(" beta      S (SE)        varS      chi      d_MM     h       links   acc    tau_int  max drift")
    for b in betas:
        r = [x for x in res if x["beta"] == b]
        S = np.array([x["S"] for x in r])
        g = lambda k: np.mean([x[k] for x in r])
        print(f" {b:4.1f} {S.mean():+9.3f} ({S.std(ddof=1) / np.sqrt(len(r)) if len(r) > 1 else 0:.3f}) {g('varS'):9.2f} "
              f"{g('chi'):.4f} {g('d_mm'):.4f} {g('h'):7.3f} {g('links'):7.2f} {g('acc'):.3f} {g('tau'):8.1f} "
              f"{max(x['max_dS_drift'] for x in r):.1e}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--scan", action="store_true")
    ap.add_argument("--N", type=int, default=50)
    ap.add_argument("--eps", type=float, default=0.12)
    ap.add_argument("--chains", type=int, default=4)
    ap.add_argument("--therm", type=int, default=500)
    ap.add_argument("--meas", type=int, default=3000)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--betas", default="")
    a = ap.parse_args()
    if a.check:
        check()
    if a.scan:
        betas = [float(x) for x in a.betas.split(",")] if a.betas else BETAS_F64A
        scan(a.N, a.eps, a.chains, a.therm, a.meas, a.workers, betas)


if __name__ == "__main__":
    main()
