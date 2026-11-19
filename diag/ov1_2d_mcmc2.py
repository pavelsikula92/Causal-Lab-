"""
OV1 phase control P1 (decides nothing): the local-move sampler of diag/ov1_2d_mcmc.py with a dense engine, for the
size scan N = 50, 100, 200, 400 at eps = 0.12 (the author, 19 Nov, 12:34: option 1).
====================================================================================================================
Same chain as diag/ov1_2d_mcmc.chain (proposal: U or V with probability 1/2, an adjacent position uniformly; Metropolis
with exp(-beta S); one sweep = N (N - 1)/2 proposals; one sample per sweep after the thermalisation). The engine keeps
the matrix C[e, f] of the number of elements strictly between every related pair, so Delta S of a proposal costs O(N)
instead of the sparse sweeps of diag/ov1_2d.py (which stay the reference): toggling the link (a, b) changes C only on
the pairs (a, y), y > b and (x, b), x < a (or (b, y), (x, a) when the link b < a is created).

Calibration of the phase bands at each N (the author, 19 Nov, 11:30): uniform 2D orders of the same N drawn directly,
seeds >= 100000, disjoint from the chain seeds; per N 20 000 orders, bands [0.5 %, 99.5 %] of chi, d_MM, h.
Reported per chain: means of S, chi, d_MM, h, N_0, var(S), acceptance, tau_int, and the fraction of samples outside
each band.

beta grid, fixed before the run (the transition moves with N, so the F6.4a grid cannot be reused):
- stage A, the same for every N: 0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.4, 2.0, 2.8, 4.0;
- stage B: 9 equally spaced interior points between the two stage-A neighbours of the stage-A maximum of var(S).
Chains: 4 per beta; sweeps 500 + 3000 (as F6.4a).

    python -m diag.ov1_2d_mcmc2 --check
    python -m diag.ov1_2d_mcmc2 --calib --N 50,100,200,400
    python -m diag.ov1_2d_mcmc2 --stage A --N 50,100,200,400 --workers 8
    python -m diag.ov1_2d_mcmc2 --stage B --N 50,100,200,400 --workers 8
    python -m diag.ov1_2d_mcmc2 --report --N 50,100,200,400
"""
from __future__ import annotations

import argparse
import glob
import json
import math
import os
import time

import numpy as np

from diag.ov1_2d import fker
from diag.ov1_kernels import _F, _init, _rel_count, _seed, height, sweeps
from diag.ov1_2d_mcmc import d_mm, tau_int

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "ov1_2d_mcmc2")
EPS = 0.12
STAGE_A = [0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.4, 2.0, 2.8, 4.0]
CAL0 = 100000


def action_from(pu, pv, C, f, eps):
    s, R, L0 = _F(pu, pv, C, f)
    return 4 * eps * (pu.size - 2 * eps * s), R, L0


def brute_C(pu, pv):
    rel = (pu[:, None] < pu[None, :]) & (pv[:, None] < pv[None, :])
    Rm = rel.astype(float)
    return np.where(rel, np.rint(Rm @ Rm), 0).astype(np.int32), rel


def bands(N):
    f = os.path.join(OUT, f"cal_{N}.json")
    return json.load(open(f)) if os.path.exists(f) else None


def chain(N, beta, therm, meas, seed, eps=EPS):
    rng = np.random.default_rng(seed)
    _seed(int(rng.integers(2 ** 31)))
    f = fker(eps, N + 2)
    pu = rng.permutation(N).astype(np.int64)
    pv = rng.permutation(N).astype(np.int64)
    eu = np.argsort(pu).astype(np.int64)
    ev = np.argsort(pv).astype(np.int64)
    C = _init(pu, pv)
    S, R, _ = action_from(pu, pv, C, f, eps)
    B = bands(N)
    moves = N * (N - 1) // 2
    acc, drift = 0, 0.0
    X = []
    t0 = time.time()
    for sw in range(therm + meas):
        S, R, a = sweeps(pu, pv, eu, ev, C, f, beta, eps, moves, S, R)
        acc += a
        if sw >= therm:
            Sx, Rx, L0 = action_from(pu, pv, C, f, eps)
            drift = max(drift, abs(Sx - S))
            assert Rx == R
            S = Sx
            chi = R / (N * (N - 1) / 2)
            X.append((S, chi, d_mm(chi), height(pu, pv), L0))
    X = np.array(X)
    out = {"N": N, "eps": eps, "beta": beta, "seed": seed, "S": float(X[:, 0].mean()), "varS": float(X[:, 0].var()),
           "chi": float(X[:, 1].mean()), "d_mm": float(X[:, 2].mean()), "h": float(X[:, 3].mean()),
           "links": float(X[:, 4].mean()), "acc": acc / ((therm + meas) * moves), "tau": tau_int(X[:, 0]),
           "tau_chi": tau_int(X[:, 1]), "max_dS_drift": drift, "sec": time.time() - t0}
    if B:
        for k, j in (("chi", 1), ("d_mm", 2), ("h", 3)):
            lo, hi = B[k]
            out[f"out_{k}"] = float(np.mean((X[:, j] < lo) | (X[:, j] > hi)))
    return out


def calibrate(N, n=20000):
    """Uniform 2D orders of N elements (the flat reference at fixed N); chi, d_MM and h only."""
    rng = np.random.default_rng(CAL0 + N)
    rows = []
    for _ in range(n):
        pu = rng.permutation(N).astype(np.int64)
        pv = rng.permutation(N).astype(np.int64)
        chi = _rel_count(pu, pv) / (N * (N - 1) / 2)
        rows.append((chi, d_mm(chi), height(pu, pv)))
    X = np.array(rows)
    out = {"N": N, "n": n, "seed": CAL0 + N}
    for k, j in (("chi", 0), ("d_mm", 1), ("h", 2)):
        lo, hi = np.percentile(X[:, j], [0.5, 99.5])
        out[k] = [float(lo), float(hi)]
        out[k + "_mean"] = float(X[:, j].mean())
    os.makedirs(OUT, exist_ok=True)
    json.dump(out, open(os.path.join(OUT, f"cal_{N}.json"), "w"), indent=1)
    return out


def check():
    ok = True
    rng = np.random.default_rng(5)
    for N, beta in ((12, 1.0), (40, 2.5), (90, 0.8)):
        f = fker(EPS, N + 2)
        pu = rng.permutation(N).astype(np.int64)
        pv = rng.permutation(N).astype(np.int64)
        eu, ev = np.argsort(pu).astype(np.int64), np.argsort(pv).astype(np.int64)
        C = _init(pu, pv)
        S, R, _ = action_from(pu, pv, C, f, EPS)
        _seed(N)
        worst = 0.0
        for _ in range(300):
            S, R, _a = sweeps(pu, pv, eu, ev, C, f, beta, EPS, 20, S, R)
            Cb, rel = brute_C(pu, pv)
            nb_ = Cb[rel]
            Sb = 4 * EPS * (N - 2 * EPS * f[nb_].sum())
            worst = max(worst, abs(Sb - S), float(np.abs(np.where(rel, C, 0) - Cb).max()))
            assert R == rel.sum() and np.all(eu[pu] == np.arange(N)) and np.all(ev[pv] == np.arange(N))
        print(f"N {N} beta {beta}: 6000 proposals, max |S - brute| and C mismatch {worst:.1e}")
        ok &= worst < 1e-9
    # against the sparse engine of diag/ov1_2d_mcmc at N 50: same stationary statistics at beta 0.1 and 2.0
    from diag.ov1_2d_mcmc import chain as chain_sparse
    for beta in (0.1, 2.0):
        a = [chain(50, beta, 200, 1000, 300 + k) for k in range(2)]
        b = [chain_sparse(50, EPS, beta, 200, 1000, 400 + k) for k in range(2)]
        print(f"beta {beta}: dense S {np.mean([x['S'] for x in a]):+.3f} chi {np.mean([x['chi'] for x in a]):.4f} "
              f"h {np.mean([x['h'] for x in a]):.3f} | sparse S {np.mean([x['S'] for x in b]):+.3f} "
              f"chi {np.mean([x['chi'] for x in b]):.4f} h {np.mean([x['h'] for x in b]):.3f}")
    print("checks:", "OK" if ok else "FAILED")


def _job(a):
    N, beta, therm, meas, seed = a
    fn = os.path.join(OUT, f"p1_{N}_{beta:.4f}_{seed}.json")
    if os.path.exists(fn):
        return json.load(open(fn))
    r = chain(N, beta, therm, meas, seed)
    json.dump(r, open(fn + ".tmp", "w"))
    os.replace(fn + ".tmp", fn)
    return r


def _load(N):
    return [json.load(open(f)) for f in glob.glob(os.path.join(OUT, f"p1_{N}_*.json"))]


def _grid_B(N):
    R = [r for r in _load(N) if any(abs(r["beta"] - b) < 1e-9 for b in STAGE_A)]
    v = {b: np.mean([r["varS"] for r in R if abs(r["beta"] - b) < 1e-9]) for b in STAGE_A}
    assert all(not math.isnan(x) for x in v.values()), "stage A incomplete"
    k = int(np.argmax([v[b] for b in STAGE_A]))
    lo, hi = STAGE_A[max(k - 1, 0)], STAGE_A[min(k + 1, len(STAGE_A) - 1)]
    return [round(x, 4) for x in np.linspace(lo, hi, 11)[1:-1] if not any(abs(x - b) < 1e-9 for b in STAGE_A)]


def run_stage(stage, Ns, chains, therm, meas, workers):
    os.makedirs(OUT, exist_ok=True)
    jobs = []
    for N in Ns:
        assert bands(N), f"calibrate N {N} first (--calib)"
        betas = STAGE_A if stage == "A" else _grid_B(N)
        for b in betas:
            for c in range(chains):
                jobs.append((N, b, therm, meas, 7919 * c + 104729 * int(round(1000 * b)) + 31 * N))
    jobs.sort(key=lambda j: -j[0])
    from concurrent.futures import ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for r in ex.map(_job, jobs):
            print(f"N {r['N']} beta {r['beta']:.4f} seed {r['seed']}: S {r['S']:+.3f} chi {r['chi']:.4f} "
                  f"h {r['h']:.2f} tau {r['tau']:.1f} ({r['sec']:.0f} s)", flush=True)


def report(Ns):
    for N in Ns:
        B = bands(N)
        R = _load(N)
        print(f"N {N}, eps {EPS}: calibration (n {B['n']}): chi [{B['chi'][0]:.4f}, {B['chi'][1]:.4f}] "
              f"d_MM [{B['d_mm'][0]:.4f}, {B['d_mm'][1]:.4f}] h [{B['h'][0]:.1f}, {B['h'][1]:.1f}]")
        print("   beta   chains      S (SE)       varS     chi     d_MM      h     links   acc  tau_S tau_chi  "
              "out chi/dMM/h")
        for b in sorted({r["beta"] for r in R}):
            r = [x for x in R if x["beta"] == b]
            S = np.array([x["S"] for x in r])
            g = lambda k: np.mean([x[k] for x in r])
            print(f" {b:7.4f}  {len(r)}  {S.mean():+10.3f} ({S.std(ddof=1) / np.sqrt(len(r)) if len(r) > 1 else 0:6.3f}) "
                  f"{g('varS'):9.2f} {g('chi'):.4f} {g('d_mm'):.4f} {g('h'):7.2f} {g('links'):7.1f} {g('acc'):.3f} "
                  f"{g('tau'):6.1f} {g('tau_chi'):6.1f}  {g('out_chi'):.2f}/{g('out_d_mm'):.2f}/{g('out_h'):.2f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--calib", action="store_true")
    ap.add_argument("--stage", choices=["A", "B"])
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--N", default="50,100,200,400")
    ap.add_argument("--chains", type=int, default=4)
    ap.add_argument("--therm", type=int, default=500)
    ap.add_argument("--meas", type=int, default=3000)
    ap.add_argument("--workers", type=int, default=1)
    a = ap.parse_args()
    Ns = [int(x) for x in a.N.split(",")]
    if a.check:
        check()
    if a.calib:
        for N in Ns:
            B = calibrate(N)
            print(f"N {N}: chi {B['chi_mean']:.4f} [{B['chi'][0]:.4f}, {B['chi'][1]:.4f}]  d_MM [{B['d_mm'][0]:.4f}, "
                  f"{B['d_mm'][1]:.4f}]  h {B['h_mean']:.2f} [{B['h'][0]:.0f}, {B['h'][1]:.0f}]")
    if a.stage:
        run_stage(a.stage, Ns, a.chains, a.therm, a.meas, a.workers)
    if a.report:
        report(Ns)


if __name__ == "__main__":
    main()
