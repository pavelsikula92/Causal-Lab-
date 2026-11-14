"""
C0, card configuration (development rows before freezing; decisions 25 and 26; decides nothing until the card is frozen)
======================================================================================================================
The instrument and every readout of docs/C0_card_draft.md, in the card's configuration: kernels MV8+ct (main) and
MV2+ct (reference), 3+1 pyramids X 7, T 6, densities 30 and 60. The development rows run on seeds 13000+; the card
will run on fresh seeds.

Readouts, per kernel and density:
  eps          max over shells r < 2.5 of |eps(r) - (A1)| (the kernel's sum per unit cone measure, sample events at
               t in [3.5, 5.5));
  norm         the constant-source sum against its expectation at each event's own t, (S / E[S](t_v)) - 1, slabs t >= 3;
  noise S      var(S - E[S](t_v)) against the diagonal model (descriptive);
  modes        the deviation of the mode amplitude from the exact solution (projection on cos(kx), slabs |t - b| < 0.25
               of the inner events); linear k = 0..3, nonlinear homogeneous A = 0.5, 1, 2, nonlinear A = 1 cos x;
               against
                 continuum (A)  the noise-free classical solution (amp-dev itself);
                 (A3)           the finite-density expectation of the kernel, no noise (linear exactly; nonlinear
                                homogeneous by the nonlinear Volterra equation);
                 J_eff (B)      the mean field with mu exp(-s^2/2) sin f, s = the measured sd of the same row pooled
                                over the sprinklings of the density (decision 26);
               residual = measurement - B;
  gauss        for the nonlinear homogeneous rows, per slab: skewness and excess kurtosis of the event noise eta (the
               remainder after the exact solution and the slab's mean deviation), and the direct check
               mean sin(phi) - exp(-s^2/2) mean sin(phi - eta) (zero for zero-mean Gaussian noise independent of phi);
  noise        the sd of the remainder by time, its growth exponent p (sd ~ t^p over t 1.5 .. 5.5), the ratio
               sd(rho 60)/sd(rho 30) per slab.

    python -m diag.c0_card --rhos 30,60 --seeds 0,1,2,3 --workers 4         (development rows)
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
from diag.c0_apex import F_k, delta, dyn_pred, dyn_pred_nl0, eps_pred
from diag.c0_ct import MODES, exact_field
from diag.c0_ct_renorm import mean_field
from diag.c0_apex import pendulum
from diag.c0_links4 import sprinkle4
from diag.c0_minvar import NMAX_ALL, minvar
from phases.F3_structure import field1d as F

_nb = F._nb
X0, T0 = 7.0, 6.0
SEED0 = 13000
KERNELS = (("MV8+ct", 8), ("MV2+ct", 2))
RBIN, NRB = 0.25, 10


def rows(rho):
    out = []
    for name, m in KERNELS:
        w = minvar(m)
        out.append((name, w, delta(w, rho)))
    return out


if _nb is not None:
    @_nb.njit(cache=True)
    def evolve_card(t, x, y, z, tc, kmax, W, ctd, phi0, lin, c, mu, sample, rbin, nrb):
        """As diag/c0_ct.evolve_ct, plus the per-event class counts and the r-histogram of the sample events."""
        N = t.size
        K = W.shape[0]
        R = phi0.shape[0]
        Mo = R // K
        phi = np.empty((R, N))
        J = np.empty((R, N))
        s = np.zeros(R)
        cnt = np.zeros((N, kmax + 1), np.int32)
        hist = np.zeros((nrb, kmax + 1))
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
                    cnt[j, nn] += 1
                    if sample[j]:
                        r = math.sqrt((x[j] - x[u]) ** 2 + (y[j] - y[u]) ** 2 + (z[j] - z[u]) ** 2)
                        ib = int(r / rbin)
                        if ib < nrb:
                            hist[ib, nn] += 1.0
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
        return phi, cnt, hist


def _ES(w, rho, tg):
    F0 = F_k(w, rho, 0, tg)
    dt = tg[1] - tg[0]
    return np.concatenate([[0.0], np.cumsum(0.5 * (F0[1:] + F0[:-1]) * dt)])


def task(rho, seed, X=X0, T=T0):
    t, x, y, z = sprinkle4(rho, X, T, SEED0 + seed)
    RS = rows(rho)
    W = np.array([r[1] for r in RS])
    ctd = np.array([r[2] for r in RS])
    phi0 = np.array([m[3] * np.cos(m[2] * x) * np.cos(m[2] * t) for _ in RS for m in MODES])
    lin = np.array([m[1] for m in MODES])
    c = 1.0 / (2.0 * math.pi * math.sqrt(6.0 * rho))
    tc = (24.0 * (NMAX_ALL + 6) / (math.pi * rho)) ** 0.25
    inner = (np.abs(x) < X - t - 0.3) & (np.abs(y) < X - t - 0.3) & (np.abs(z) < X - t - 0.3)
    sample = inner & (t >= 3.5) & (t < 5.5)
    a0 = time.time()
    phi, cnt, hist = evolve_card(t, x, y, z, tc, NMAX_ALL, W, ctd, phi0, lin, c, 1.0, sample, RBIN, NRB)
    sec = time.time() - a0
    tb = np.arange(0.25, T - 0.4, 0.5)
    slabs = [inner & (np.abs(t - b) < 0.25) for b in tb]
    ex = {m[0]: exact_field(m, t, x, T) for m in MODES}
    tg = np.arange(0.0, T + 0.01, 0.005)
    out = {"rho": rho, "seed": seed, "N": int(t.size), "sec": sec, "t": tb.tolist(), "X": X, "T": T,
           "nslab": [int(m.sum()) for m in slabs], "nsample": int(sample.sum()), "k": {}}
    Mo = len(MODES)
    rr = (np.arange(NRB) + 0.5) * RBIN
    for kk, (name, w, ct) in enumerate(RS):
        S = c * (cnt @ w)
        E = np.interp(t, tg, _ES(w, rho, tg))
        K = {"ct": ct,
             "eps": (c * (hist @ w) / sample.sum() / (rr * RBIN) - 1.0).tolist(),
             "norm": [float(S[m].mean() / E[m].mean() - 1.0) for m in slabs],
             "varS": [float((S[m] - E[m]).var()) for m in slabs],
             "vardiag": [float(c * c * (cnt[m].mean(0) @ (w * w))) for m in slabs],
             "modes": {}}
        for mo, m in enumerate(MODES):
            ph = phi[kk * Mo + mo]
            dev = (ph - ex[m[0]]) / m[3]
            bas = np.cos(m[2] * x)
            dd = []
            for sl in slabs:
                amp = float(dev[sl] @ bas[sl] / (bas[sl] @ bas[sl]))
                eta = (dev[sl] - amp * bas[sl]) * m[3]
                row = [amp, float(eta.std() / m[3])]
                if not m[1] and m[2] == 0:
                    sd = eta.std()
                    z3 = float(np.mean((eta / sd) ** 3)) if sd > 0 else 0.0
                    z4 = float(np.mean((eta / sd) ** 4) - 3.0) if sd > 0 else 0.0
                    direct = float(np.mean(np.sin(ph[sl])) - math.exp(-0.5 * sd * sd) * np.mean(np.sin(ph[sl] - eta)))
                    row += [z3, z4, direct]
                dd.append(row)
            K["modes"][m[0]] = dd
        out["k"][name] = K
    return out


def cached(*a):
    key = hashlib.sha256(json.dumps([str(q) for q in a] + ["card-v1"]).encode()).hexdigest()[:20]
    fn = os.path.join("diag", "out", "c0_card_cache", key + ".json")
    if os.path.exists(fn):
        return json.load(open(fn))
    out = task(*a)
    os.makedirs(os.path.dirname(fn), exist_ok=True)
    json.dump(out, open(fn + ".tmp", "w"))
    os.replace(fn + ".tmp", fn)
    return out


def _slabavg(dev, tg, tb, X):
    out = []
    for b in tb:
        m = np.abs(tg - b) < 0.25
        wt = np.clip(X - tg[m] - 0.3, 0, None) ** 3
        out.append(float(np.sum(dev[m] * wt) / np.sum(wt)))
    return np.array(out)


def growth(tb, sd, lo=1.5, hi=5.5):
    m = (tb >= lo) & (tb <= hi) & (sd > 0)
    return float(np.polyfit(np.log(tb[m]), np.log(sd[m]), 1)[0])


def analyse(R):
    """Every readout of one density, with the SE over sprinklings; returns a dict and prints it."""
    rho = R[0]["rho"]
    X = R[0]["X"]
    tb = np.array(R[0]["t"])
    n = len(R)
    late = tb >= 3.0
    tg = np.arange(0.0, tb[-1] + 0.3, 0.005)
    out = {}
    print(f"######## rho {rho:g}; seeds {sorted(r['seed'] for r in R)}; N {np.mean([r['N'] for r in R]):.0f}; "
          f"{np.mean([r['sec'] for r in R]):.0f} s per sprinkling (t = " + " ".join(f"{v:.2f}" for v in tb[1::2]) + ")")
    for name, w, ct in rows(rho):
        Ks = [r["k"][name] for r in R]
        o = {}
        eps = np.array([k["eps"] for k in Ks])
        ep = np.array(eps_pred(w, rho, RBIN, NRB))
        d = eps.mean(0) - ep
        i = int(np.argmax(np.abs(d)))
        o["eps"] = (float(abs(d[i])), float(eps[:, i].std(ddof=1) / math.sqrt(n)))
        nm = np.array([np.mean(np.array(k["norm"])[late]) for k in Ks])
        o["norm"] = (float(nm.mean()), float(nm.std(ddof=1) / math.sqrt(n)))
        vr = np.mean([np.array(k["varS"])[late].mean() for k in Ks])
        vd = np.mean([np.array(k["vardiag"])[late].mean() for k in Ks])
        print(f"== {name} (counterterm delta {ct:+.4f})")
        print(f"   eps: max |meas - (A1)| {o['eps'][0]:.4f} +- {o['eps'][1]:.4f};  norm - 1: {o['norm'][0]:+.4f} +- "
              f"{o['norm'][1]:.4f};  var(S - E[S]) {vr:.3f} = {vr / vd:.2f} x diagonal")
        o["modes"] = {}
        for m in MODES:
            D = np.array([k["modes"][m[0]] for k in Ks])
            amp = D[:, :, 0]
            mean = amp.mean(0)
            se = amp.std(0, ddof=1) / math.sqrt(n)
            sd = D[:, :, 1].mean(0)
            mname, lin, kk, A = m
            P = None
            B = None
            if lin:
                P = np.array(dyn_pred(w, rho, kk, tb, X=X, ct=ct))
            elif kk == 0:
                P = np.array(dyn_pred_nl0(w, rho, A, tb, X=X, ct=ct))
                s = np.interp(tg, np.r_[0.0, tb], np.r_[0.0, sd * A])
                B = _slabavg((mean_field(w, ct, rho, A, tg, s, 0.005) - pendulum(A, tg)) / A, tg, tb, X)
            mo = {"cont": (float(np.max(np.abs(mean))), float(se[int(np.argmax(np.abs(mean)))])),
                  "sd": sd.tolist(), "p": growth(tb, sd)}
            line = f"   {mname:9s} amp-dev " + " ".join(f"{v:+.2f}" for v in mean[1::2])
            if P is not None:
                dP = mean - P
                j = int(np.argmax(np.abs(dP)))
                mo["A3"] = (float(abs(dP[j])), float(se[j]))
                line += "  (A3) " + " ".join(f"{v:+.2f}" for v in P[1::2])
            if B is not None:
                dB = mean - B
                j = int(np.argmax(np.abs(dB)))
                mo["B"] = (float(abs(dB[j])), float(se[j]))
                mo["renorm"] = float(np.max(np.abs(B - P)))
                line += "  B " + " ".join(f"{v:+.2f}" for v in B[1::2])
            print(line)
            extra = f"      max|A| {mo['cont'][0]:.3f}+-{mo['cont'][1]:.3f}"
            if "A3" in mo:
                extra += f"  max|meas-(A3)| {mo['A3'][0]:.3f}+-{mo['A3'][1]:.3f}"
            if "B" in mo:
                extra += f"  residual max|meas-B| {mo['B'][0]:.3f}+-{mo['B'][1]:.3f}  renorm max|B-(A3)| {mo['renorm']:.3f}"
            extra += "  sd " + " ".join(f"{v:.2f}" for v in sd[1::2]) + f"  p {mo['p']:.2f}"
            print(extra)
            if not lin and kk == 0:
                z3 = D[:, :, 2].mean(0)
                z4 = D[:, :, 3].mean(0)
                dr = D[:, :, 4].mean(0)
                mo["gauss"] = (float(np.max(np.abs(z3[1:]))), float(np.max(np.abs(z4[1:]))), float(np.max(np.abs(dr[1:]))))
                print(f"      gauss: skew " + " ".join(f"{v:+.2f}" for v in z3[1::2]) + "  exc.kurt "
                      + " ".join(f"{v:+.2f}" for v in z4[1::2]) + "  direct " + " ".join(f"{v:+.3f}" for v in dr[1::2]))
            o["modes"][mname] = mo
        out[name] = o
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=0)
    ap.add_argument("--rhos", type=str, default="30,60")
    ap.add_argument("--seeds", type=str, default="0,1,2,3")
    a = ap.parse_args()
    rhos = [float(r) for r in a.rhos.split(",")]
    jobs = [(rho, int(s)) for rho in sorted(rhos, reverse=True) for s in a.seeds.split(",")]
    t0 = time.time()
    res = parallel_map(cached, jobs, default_workers(a.workers), print,
                       lambda i, r, d, n: f"[C0 card dev] {d}/{n}: rho {r['rho']} seed {r['seed']} N {r['N']} "
                                          f"{r['sec']:.0f}s ({time.time() - t0:.0f}s)")
    fn = os.path.join("diag", "out", f"c0_card_dev_{time.strftime('%Y%m%dT%H%M%S', time.gmtime())}.json")
    json.dump(res, open(fn, "w"))
    print("saved", fn)
    A = {}
    for rho in sorted(rhos):
        A[rho] = analyse([r for r in res if r["rho"] == rho])
    if len(rhos) >= 2:
        lo, hi = sorted(rhos)[0], sorted(rhos)[-1]
        print(f"######## density: rho {hi:g} against rho {lo:g}")
        for name, _ in KERNELS:
            for m in MODES:
                a, b = A[lo][name]["modes"][m[0]], A[hi][name]["modes"][m[0]]
                r = np.array(b["sd"]) / np.array(a["sd"])
                print(f"   {name:7s} {m[0]:9s} max|A| {a['cont'][0]:.3f} -> {b['cont'][0]:.3f}   sd ratio (t 2.75, 4.75) "
                      f"{r[5]:.3f} {r[9]:.3f}   (rho^-1/4: {(lo / hi) ** 0.25:.3f})")


if __name__ == "__main__":
    main()
