"""
Direction A, development probe: a causal growth rule of the order with feedback from the field (1+1; decides nothing)
====================================================================================================================

Why (11 Nov 2026; approved by the author at 11:59 as the development step before 2+1, with his corrections)
---------------------------------------------------------------------------------------------------------
The chain phi -> load -> response of the order -> number of events -> fall (decision 15) needs a dynamics of the order.
Candidate growth rule: every candidate event x of a Poisson sprinkling of (proper) density rho, carrying an independent
random mark u_x ~ U(0, 1) attached at the sprinkling, is realised iff u_x < Pi(x)^2, where the clock rate Pi = e^L is
computed from the REALISED causal past of x only, with the normalisation of the retarded Green function of 1+1
(1/2 inside the cone, as for the field):
    L(x) = -(kappa / 2 rho) * sum over realised y < x of q(y),   q = 2 V_A(phi_A) + 2 V_B(phi_B)  (the load, trace of T).
The realised events have coordinate density rho Pi^2, so the sum carries the factor Pi^2 itself and the continuum law is
    (d_t^2 - d_x^2) L = -kappa Pi^2 q,   i.e. box_g L = -kappa q in the metric ds^2 = Pi^2 (dt^2 - dx^2).
In 1+1 every metric is conformally flat: a Poisson process of density rho Pi^2 IS a sprinkling of that metric, with the
causal relation of Minkowski space; the geometry is in the number of events (proper times), not in the cones. (In 2+1
this is no longer enough: the relation itself will have to change.)
Two separate fields phi_A (mu_A) and phi_B (mu_B) live on the realised events, each with the UNCHANGED rule of F3.14a,
    phi(x) = phi_0(x) - (1 / 2 rho) sum over realised y < x of mu sin(phi(y)),
so they interact only through the order (no direct kink-kink force). Kink A at x = -d/2, kink B at +d/2, both at rest.
Expected (1+1, Newton in one dimension): once the response of A reaches B (t > d), B accelerates towards A with
a_B ~ kappa M_A / 2, independent of M_B; A towards B with kappa M_B / 2.

What the probe measures (the author's list):
  A0  causal: the rule reads only the realised past (by construction; checked by I/O: a change of a later mark leaves
      every earlier event unchanged);
  A1  order of computation: two linear extensions (t and t + 0.5 x), the same marks -> the same realised set and fields;
  A2  covariance: reported later (needs a boosted launch; not in this probe);
  A3  density: realised fraction per bin against the mean of Pi^2 over the candidates of the bin; interval statistics:
      for random realised timelike pairs, the realised count inside the interval against rho * (the integral of Pi^2 over
      the coordinate interval) (the proper volume);
  A4  continuum: the coupled equations (d_t^2 - d_x^2) phi_i = -Pi^2 V_i'(phi_i), (d_t^2 - d_x^2) L = -kappa Pi^2 q,
      L(0) = L_t(0) = 0, by finite differences; the centres of A and B and L along x at the end;
  A5  active source: M_A = 4, 8, 16 (mu_A 0.25, 1, 4) with M_B = 8;
  A6  universality: M_B = 4, 8 (mu_B 0.25, 1) with M_A = 8;
  A7  null: kappa = 0 on the same candidates and marks (every candidate realised).
Readout: the kink fit of F3.14a per slab of width 1 for each field; displacement relative to the paired kappa = 0 run.

Run from the root of the laboratory (numba required):
    python -m diag.growth_probe --quick
    python -m diag.growth_probe --seeds 32 --rho 800 --workers 8
"""

from __future__ import annotations

import argparse
import json
import math
import os
import time

import numpy as np

from core.order.sprinkle import derive_seed
from core.stats.parallel import parallel_map, default_workers
from phases.F3_structure import field1d as F

_nb = F._nb

if _nb is not None:
    @_nb.njit(cache=True)
    def _grow(U, V, order, marks, phiA0, phiB0, rho, muA, muB, kappa):
        """Fenwick tree of Fenwick trees over all candidates; only realised events are inserted (3 weight channels)."""
        N = U.size
        ru = np.empty(N, np.int64)
        rv = np.empty(N, np.int64)
        su = np.argsort(U)
        sv = np.argsort(V)
        for k in range(N):
            ru[su[k]] = k + 1
            rv[sv[k]] = k + 1
        cnt = np.zeros(N + 1, np.int64)
        for i in range(N):
            k = ru[i]
            while k <= N:
                cnt[k] += 1
                k += k & (-k)
        off = np.zeros(N + 2, np.int64)
        for k in range(1, N + 1):
            off[k + 1] = off[k] + cnt[k]
        lst = np.empty(off[N + 1], np.int64)
        fill = off.copy()
        for i in range(N):
            k = ru[i]
            while k <= N:
                lst[fill[k]] = rv[i]
                fill[k] += 1
                k += k & (-k)
        for k in range(1, N + 1):
            lst[off[k]:off[k + 1]] = np.sort(lst[off[k]:off[k + 1]])
        size = off[N + 1] + N + 1
        tA = np.zeros(size)
        tB = np.zeros(size)
        tQ = np.zeros(size)
        phiA = np.full(N, np.nan)
        phiB = np.full(N, np.nan)
        Lx = np.empty(N)
        real = np.zeros(N, np.bool_)
        for a in range(N):
            i = order[a]
            sA = 0.0
            sB = 0.0
            sQ = 0.0
            k = ru[i] - 1
            while k > 0:
                lo, hi = off[k], off[k + 1]
                p = np.searchsorted(lst[lo:hi], rv[i])
                base = lo + k
                q = p
                while q > 0:
                    sA += tA[base + q - 1]
                    sB += tB[base + q - 1]
                    sQ += tQ[base + q - 1]
                    q -= q & (-q)
                k -= k & (-k)
            L = -0.5 * kappa * sQ / rho
            Lx[i] = L
            if marks[i] >= math.exp(2.0 * L):
                continue
            real[i] = True
            pa = phiA0[i] - 0.5 * sA / rho
            pb = phiB0[i] - 0.5 * sB / rho
            phiA[i] = pa
            phiB[i] = pb
            wA = muA * math.sin(pa)
            wB = muB * math.sin(pb)
            wQ = 2.0 * muA * (1.0 - math.cos(pa)) + 2.0 * muB * (1.0 - math.cos(pb))
            k = ru[i]
            while k <= N:
                lo, hi = off[k], off[k + 1]
                p = np.searchsorted(lst[lo:hi], rv[i]) + 1
                base = lo + k
                n_k = hi - lo
                q = p
                while q <= n_k:
                    tA[base + q - 1] += wA
                    tB[base + q - 1] += wB
                    tQ[base + q - 1] += wQ
                    q += q & (-q)
                k += k & (-k)
        return phiA, phiB, Lx, real

    @_nb.njit(cache=True)
    def _interval_counts(U, V, real, pairs):
        out = np.zeros(pairs.shape[0], np.int64)
        for k in range(pairs.shape[0]):
            y, x = pairs[k, 0], pairs[k, 1]
            c = 0
            for z in range(U.size):
                if real[z] and U[z] > U[y] and U[z] < U[x] and V[z] > V[y] and V[z] < V[x]:
                    c += 1
            out[k] = c
        return out


def evolve_growth(t, x, marks, muA, muB, xA, xB, rho, kappa, order="t"):
    U, V = (t + x) / math.sqrt(2.0), (t - x) / math.sqrt(2.0)
    od = np.arange(t.size) if order == "t" else np.argsort(t + 0.5 * x, kind="stable")
    phiA0 = F.kink_free(t, x - xA, 0.0, muA)
    phiB0 = F.kink_free(t, x - xB, 0.0, muB)
    return _grow(U, V, od.astype(np.int64), marks, phiA0, phiB0, float(rho), float(muA), float(muB), float(kappa))


def track(t, x, phi, real, mu, x0, T):
    tt, xx, pp = t[real], x[real], phi[real]
    xc = []
    for a in np.arange(0.0, T - 1.0 + 1e-9, 1.0):
        g = xc[-1] if (xc and np.isfinite(xc[-1])) else x0
        r = F.fit_kink(tt, xx, pp, a, a + 1.0, g, mu)
        xc.append(r["xc"] if r else float("nan"))
    return xc


def continuum(muA, muB, xA, xB, kappa, T, dx=0.01, X=60.0):
    xs = np.arange(-X, X + dx / 2, dx)
    dt = 0.5 * dx
    lap = lambda p: np.concatenate([[0.0], (p[2:] - 2 * p[1:-1] + p[:-2]) / dx ** 2, [0.0]])
    a0 = F.kink_exact(0.0, xs - xA, 0.0, muA)
    b0 = F.kink_exact(0.0, xs - xB, 0.0, muB)
    l0 = np.zeros_like(xs)

    def rhs(a, b, l):
        om2 = np.exp(2.0 * l)
        q = 2.0 * muA * (1 - np.cos(a)) + 2.0 * muB * (1 - np.cos(b))
        return lap(a) - om2 * muA * np.sin(a), lap(b) - om2 * muB * np.sin(b), lap(l) - kappa * om2 * q

    ra, rb, rl = rhs(a0, b0, l0)
    pa, pb, pl = a0, b0, l0
    ca, cb, cl = a0 + 0.5 * dt * dt * ra, b0 + 0.5 * dt * dt * rb, l0 + 0.5 * dt * dt * rl
    cenA, cenB, cnt = np.zeros(int(T)), np.zeros(int(T)), np.zeros(int(T))
    nsteps = int(round(T / dt))
    for n in range(1, nsteps):
        tn = n * dt
        s = int(tn)
        if s < int(T) and n % 4 == 0:
            for arr, cen, x0 in ((ca, cenA, xA), (cb, cenB, xB)):
                j = np.nonzero((arr[:-1] - math.pi) * (arr[1:] - math.pi) <= 0)[0]
                if j.size:
                    j0 = j[np.argmin(np.abs(xs[j] - x0))]
                    cen[s] += xs[j0] + (math.pi - arr[j0]) / (arr[j0 + 1] - arr[j0]) * dx
            cnt[s] += 1
        ra, rb, rl = rhs(ca, cb, cl)
        na, nb, nl = 2 * ca - pa + dt * dt * ra, 2 * cb - pb + dt * dt * rb, 2 * cl - pl + dt * dt * rl
        pa, pb, pl, ca, cb, cl = ca, cb, cl, na, nb, nl
    return (cenA / np.maximum(cnt, 1)).tolist(), (cenB / np.maximum(cnt, 1)).tolist(), \
        np.interp(np.linspace(-15, 15, 31), xs, cl).tolist()


def task(rho, seed, muA, muB, kappa, d, T, extra):
    X = T + 16.0
    t, x = F.sprinkle_trapezoid(rho, X, T, seed)
    marks = np.random.default_rng(derive_seed(seed, "marks")).uniform(0.0, 1.0, t.size)
    xA, xB = -0.5 * d, 0.5 * d
    out = {"rho": rho, "seed": seed, "muA": muA, "muB": muB, "kappa": kappa, "N": int(t.size)}
    res = {}
    for kap in (kappa, 0.0):
        phiA, phiB, L, real = evolve_growth(t, x, marks, muA, muB, xA, xB, rho, kap)
        res[kap] = (phiA, phiB, L, real)
        tag = "k" if kap > 0 else "0"
        out[f"xcA_{tag}"] = track(t, x, phiA, real, muA, xA, T)
        out[f"xcB_{tag}"] = track(t, x, phiB, real, muB, xB, T)
    phiA, phiB, L, real = res[kappa]
    out["frac_real"] = float(real.mean())
    # A3: realised fraction against mean Pi^2 per bin (last third of the run, near A and B)
    bins = np.linspace(-12, 12, 13)
    s = t > (2.0 * T / 3.0)
    idx = np.digitize(x[s], bins) - 1
    fr, pr = [], []
    for k in range(len(bins) - 1):
        m = idx == k
        fr.append(float(real[s][m].mean()) if m.any() else float("nan"))
        pr.append(float(np.exp(2 * L[s][m]).mean()) if m.any() else float("nan"))
    out["bin_frac"], out["bin_pi2"] = fr, pr
    # L along x at the end (realised and candidate events in the last slab)
    s2 = t > T - 1.0
    out["L_end"] = [float(L[s2][(np.abs(x[s2] - xg) < 0.5)].mean()) if np.any(np.abs(x[s2] - xg) < 0.5) else float("nan")
                    for xg in np.linspace(-15, 15, 31)]
    if extra:
        # A3 intervals: random realised timelike pairs near A, count of realised events inside / rho / proper volume
        rng = np.random.default_rng(derive_seed(seed, "pairs"))
        ri = np.nonzero(real & (np.abs(x - xA) < 4.0) & (t > T / 2))[0]
        pairs = []
        U, V = (t + x) / math.sqrt(2.0), (t - x) / math.sqrt(2.0)
        while len(pairs) < 60 and ri.size > 10:
            a, b = rng.choice(ri, 2, replace=False)
            if U[a] < U[b] and V[a] < V[b] and (U[b] - U[a]) * (V[b] - V[a]) > 1.0:
                pairs.append((a, b))
        if pairs:
            pairs = np.array(pairs, np.int64)
            n = _interval_counts(U, V, real, pairs)
            ratio, ratio_flat = [], []
            for (a, b), c in zip(pairs, n):
                inside = (U > U[a]) & (U < U[b]) & (V > V[a]) & (V < V[b])
                vol = (U[b] - U[a]) * (V[b] - V[a])                       # coordinate area
                pv = vol * float(np.exp(2 * L[inside]).mean())             # proper area = int Pi^2
                ratio.append(c / (rho * pv))
                ratio_flat.append(c / (rho * vol))
            out["interval_ratio"] = float(np.mean(ratio))
            out["interval_ratio_flat"] = float(np.mean(ratio_flat))
        # A1: a second linear extension, same marks
        pa2, pb2, L2, r2 = evolve_growth(t, x, marks, muA, muB, xA, xB, rho, kappa, order="boost")
        out["A1_same_set"] = bool(np.array_equal(r2, real))
        out["A1_max_diff"] = float(np.nanmax(np.abs(np.nan_to_num(pa2 - phiA)))) if np.array_equal(r2, real) else float("nan")
        # A0: change the marks of the last 10 % of candidates (latest in t): earlier events unchanged
        m2 = marks.copy()
        late = t > np.quantile(t, 0.9)
        m2[late] = np.random.default_rng(1).uniform(0, 1, late.sum())
        pa3, pb3, L3, r3 = evolve_growth(t, x, m2, muA, muB, xA, xB, rho, kappa)
        early = ~late
        out["A0_early_unchanged"] = bool(np.array_equal(r3[early], real[early]) and
                                         np.allclose(np.nan_to_num(pa3[early]), np.nan_to_num(phiA[early])))
    return out


def run_task(*a):
    return task(*a)


def _ms(a):
    a = np.asarray([v for v in a if np.isfinite(v)], float)
    return (float(a.mean()), float(a.std(ddof=1) / math.sqrt(a.size))) if a.size > 1 else (float("nan"), float("nan"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=32)
    ap.add_argument("--rho", type=float, default=800.0)
    ap.add_argument("--kappa", type=float, default=0.003)
    ap.add_argument("--d", type=float, default=6.0)
    ap.add_argument("--T", type=float, default=20.0)
    ap.add_argument("--workers", type=int, default=0)
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--base", type=int, default=20261111)
    a = ap.parse_args()
    if _nb is None:
        raise SystemExit("numba required")
    rho, ns, T = a.rho, a.seeds, a.T
    if a.quick:
        rho, ns, T = 200.0, 4, 14.0
    combos = [(0.25, 1.0), (1.0, 1.0), (4.0, 1.0), (1.0, 0.25)]          # (mu_A, mu_B)
    jobs = [(rho, derive_seed(a.base, "growth", k), mA, mB, a.kappa, a.d, T, k < 2)
            for k in range(ns) for (mA, mB) in combos]
    t0 = time.time()
    res = parallel_map(run_task, jobs, default_workers(a.workers), print,
                       lambda i, r, d, n: f"[growth] {d}/{n}: muA {r['muA']:g} muB {r['muB']:g} N {r['N']} "
                                          f"({time.time() - t0:.0f}s)")
    tm = np.arange(0.5, T, 1.0)
    L = [f"1+1 growth rule, rho {rho:g}, kappa {a.kappa:g}, d {a.d:g} (A at {-a.d / 2:g}, B at {a.d / 2:g}), T {T:g}, "
         f"{ns} sprinklings; displacements relative to the paired kappa = 0 run (same candidates and marks)"]
    for (mA, mB) in combos:
        R = [r for r in res if r["muA"] == mA and r["muB"] == mB]
        MA, MB = 8 * math.sqrt(mA), 8 * math.sqrt(mB)
        cA, cB, Lc = continuum(mA, mB, -a.d / 2, a.d / 2, a.kappa, T)
        c0A, c0B, _ = continuum(mA, mB, -a.d / 2, a.d / 2, 0.0, T)
        dB_c = (cB[-1] - cB[0]) - (c0B[-1] - c0B[0])
        dA_c = (cA[-1] - cA[0]) - (c0A[-1] - c0A[0])
        dB = _ms([(r["xcB_k"][-1] - r["xcB_k"][0]) - (r["xcB_0"][-1] - r["xcB_0"][0]) for r in R])
        dA = _ms([(r["xcA_k"][-1] - r["xcA_k"][0]) - (r["xcA_0"][-1] - r["xcA_0"][0]) for r in R])
        newton_B = -0.5 * a.kappa * MA * 0.5 * max(0.0, tm[-1] - a.d) ** 2
        L.append(f"--- M_A {MA:g}, M_B {MB:g}")
        L.append(f"  B moves towards A by {dB[0]:+.4f} +- {dB[1]:.4f} (continuum {dB_c:+.4f}; 1D Newton estimate "
                 f"{newton_B:+.4f}); ratio to continuum {dB[0] / dB_c:.3f}")
        L.append(f"  A moves by {dA[0]:+.4f} +- {dA[1]:.4f} (continuum {dA_c:+.4f}); towards B is +")
        L.append(f"  realised fraction {np.mean([r['frac_real'] for r in R]):.4f}")
        bf = np.nanmean(np.array([r["bin_frac"] for r in R]), axis=0)
        bp = np.nanmean(np.array([r["bin_pi2"] for r in R]), axis=0)
        L.append("  A3 realised fraction / mean Pi^2 per x bin (t > 2T/3): " + " ".join(f"{f / p:.3f}" for f, p in zip(bf, bp)))
        Le = np.nanmean(np.array([r["L_end"] for r in R]), axis=0)
        L.append("  L at the end, x = -15..15 step 1 (discrete / continuum): " +
                 " ".join(f"{u:+.3f}/{v:+.3f}" for u, v in list(zip(Le, Lc))[::3]))
        ex = [r for r in R if "interval_ratio" in r]
        if ex:
            L.append(f"  A3 intervals near A: count / (rho * proper area) {np.mean([r['interval_ratio'] for r in ex]):.3f}; "
                     f"count / (rho * coordinate area) {np.mean([r['interval_ratio_flat'] for r in ex]):.3f}")
            L.append(f"  A1 same realised set {all(r['A1_same_set'] for r in ex)}, max field difference "
                     f"{max(r['A1_max_diff'] for r in ex):.1e}; A0 earlier events unchanged by later marks "
                     f"{all(r['A0_early_unchanged'] for r in ex)}")
    print("\n".join(L))
    os.makedirs(os.path.join("diag", "out"), exist_ok=True)
    fn = os.path.join("diag", "out", f"growth_probe_{time.strftime('%Y%m%dT%H%M%S', time.gmtime())}.json")
    with open(fn, "w") as fh:
        json.dump({"args": vars(a), "summary": L, "results": res}, fh)
    print(f"raw records: {fn} ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
