"""
C0, development: the late deviation of the k = 0 rows of F3.16 (14 Nov; decides nothing)
========================================================================================
F3.16 (INSTRUMENT_FAILED) showed about −0.02 from (A3) and from B in the last two slabs of the k = 0 rows. The
predictions P1–P3 of docs/derivations/c0_mean_exactness.md were fixed before this run:
  P1  linear (lin0, lin1): meas − (A3) = 0 within 2 SE in the last two slabs; the mean of the linear rule obeys (A3)
      exactly on a Poisson sprinkling (past-complete domain, disjoint interval, Mecke);
  P2  nonlinear homogeneous: R = meas − B equals (S), the response to the measured direct check D;
  P3  no dependence on the number of events beyond the SE (X 7 and X 8.5 put 1 and about 2.4 times the events into the
      same late slabs).
Configurations: rho 20 and 30, X 7 and 8.5, T 6, 16 fresh sprinklings each (seeds 15000+ for X 7, 16000+ for X 8.5);
kernels MV8+ct and MV2+ct (diag/c0_card.task, unchanged).

    python -m diag.c0_late --workers 4
"""
from __future__ import annotations

import argparse
import json
import math
import os
import time

import numpy as np

from core.stats.parallel import parallel_map, default_workers
from diag import c0_card as C
from diag.c0_apex import F_k, dyn_pred, pendulum

RHOS = (20.0, 30.0)
XS = (7.0, 8.5)
SEED0 = {7.0: 15000, 8.5: 16000}
T = 6.0
N = 16


def job(rho, seed, X):
    return C.cached(rho, seed, X, T, SEED0[X])


def mean_field_D(w, ct, rho, A, tg, s, D, dt, F):
    """B with the measured direct check added to the mean source: E sin phi = exp(-s^2/2) sin f + D(t)  ((S))."""
    q = np.exp(-0.5 * s ** 2)
    f = np.zeros_like(tg)
    J = np.zeros_like(tg)
    for i in range(tg.size):
        acc = 0.5 * F[i] * J[0] if i > 0 else 0.0
        if i > 1:
            acc += np.dot(F[1:i], J[i - 1:0:-1])
        rhs = A - dt * acc
        p = rhs
        for _ in range(50):
            dp = (p - ct * (q[i] * math.sin(p) + D[i]) - rhs) / (1 - ct * q[i] * math.cos(p))
            p -= dp
            if abs(dp) < 1e-14:
                break
        f[i] = p
        J[i] = q[i] * math.sin(p) + D[i]
    return f


def slabavg(dev, tg, tb, X):
    out = []
    for b in tb:
        m = np.abs(tg - b) < 0.25
        wt = np.clip(X - tg[m] - 0.3, 0, None) ** 3
        out.append(float(np.sum(dev[m] * wt) / np.sum(wt)))
    return np.array(out)


def analyse(R, X):
    rho = R[0]["rho"]
    tb = np.array(R[0]["t"])
    n = len(R)
    dt = 0.005
    tg = np.arange(0.0, tb[-1] + 0.3, dt)
    late = slice(len(tb) - 2, len(tb))
    show = slice(len(tb) - 4, len(tb))
    print(f"#### rho {rho:g}, X {X:g}: {n} sprinklings; events per slab (t = "
          + " ".join(f"{v:.2f}" for v in tb[show]) + "): " + " ".join(str(v) for v in np.array(R[0]["nslab"])[show]))
    out = {}
    for name, w, ct in C.rows(rho):
        F0 = F_k(w, rho, 0, tg)
        print(f"== {name}")
        for m in C.MODES:
            mname, lin, k, A = m
            if not (lin and k in (0, 1)) and not (not lin and k == 0):
                continue
            D = np.array([r["k"][name]["modes"][mname] for r in R])
            amp = D[:, :, 0]
            if lin:
                P = np.array(dyn_pred(w, rho, k, tb, X=X, ct=ct))
                res = amp - P
                lab = "meas - (A3)"
                pred = np.zeros_like(P)
            else:
                sd = D[:, :, 1].mean(0)
                s = np.interp(tg, np.r_[0.0, tb], np.r_[0.0, sd * A])
                from diag.c0_ct_renorm import mean_field
                fB = mean_field(w, ct, rho, A, tg, s, dt, F0)
                Dd = np.interp(tg, np.r_[0.0, tb], np.r_[0.0, D[:, :, 4].mean(0)])
                fS = mean_field_D(w, ct, rho, A, tg, s, Dd, dt, F0)
                ex = pendulum(A, tg)
                B = slabavg((fB - ex) / A, tg, tb, X)
                S = slabavg((fS - ex) / A, tg, tb, X)
                res = amp - B
                pred = S - B
                lab = "meas - B"
            mu, se = res.mean(0), res.std(0, ddof=1) / math.sqrt(n)
            lastk = res[:, late].mean(1)
            lm, ls = float(lastk.mean()), float(lastk.std(ddof=1) / math.sqrt(n))
            lp = float(pred[late].mean())
            out[(name, mname)] = (lm, ls, lp)
            print(f"   {mname:9s} {lab:12s} " + " ".join(f"{a:+.4f}({b:.4f})" for a, b in zip(mu[show], se[show]))
                  + f" | last two: {lm:+.4f} +- {ls:.4f}, predicted {lp:+.4f}  ({(lm - lp) / ls if ls > 0 else 0:+.1f} SE)")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=0)
    ap.add_argument("--n", type=int, default=N)
    a = ap.parse_args()
    jobs = [(rho, s, X) for X in sorted(XS, reverse=True) for rho in sorted(RHOS, reverse=True) for s in range(a.n)]
    t0 = time.time()
    res = parallel_map(job, jobs, default_workers(a.workers), print,
                       lambda i, r, d, nn: f"[C0 late] {d}/{nn}: rho {r['rho']:g} X {r['X']:g} seed {r['seed']} "
                                           f"{r['sec']:.0f}s ({time.time() - t0:.0f}s)")
    fn = os.path.join("diag", "out", f"c0_late_{time.strftime('%Y%m%dT%H%M%S', time.gmtime())}.json")
    json.dump(res, open(fn, "w"))
    print("saved", fn)
    allr = {}
    for X in XS:
        for rho in RHOS:
            allr[(rho, X)] = analyse([r for r in res if r["rho"] == rho and r["X"] == X], X)
    print("#### P1 pooled (inverse-variance) over the four configurations, last two slabs, meas - (A3):")
    for name, _ in C.KERNELS:
        for mname in ("lin0", "lin1"):
            v = [allr[c][(name, mname)] for c in allr]
            wts = np.array([1 / q[1] ** 2 for q in v])
            m = float(np.sum(wts * np.array([q[0] for q in v])) / wts.sum())
            print(f"   {name:7s} {mname}: {m:+.4f} +- {1 / math.sqrt(wts.sum()):.4f}")


if __name__ == "__main__":
    main()
