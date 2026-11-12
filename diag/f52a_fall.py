"""
F5.2a development rows: free fall of a field defect in a 1+1 order of nonuniform density (development; decides nothing)
=====================================================================================================================

Why (10 Nov 2026; card docs/F5_2a_card_draft.md, revision 1, approved by the author at 12:30: ID F5.2a, profile,
three masses, counting null as an instrument condition; tolerances from these rows)
-------------------------------------------------------------------------------------------------------------------
Every 1+1 metric is conformally flat, ds^2 = Omega^2 (dt^2 - dx^2): the light cones and the causal order are those of
flat spacetime, only the number of events per coordinate area changes (rho Omega^2 for a sprinkling of proper density
rho). The rule of F3.14a, phi(x) = phi_0(x) - (1/2 rho) sum over the causal past of V'(phi), is not changed; on such an
order it is the field equation in the curved spacetime (the massless operator is conformally invariant in 1+1). A kink
that falls falls because of the counting.

Arms per sprinkling (Phi(x) = g L tanh(x / L), Omega = exp(Phi), g = 0.05, L = 10; trapezoid |x| <= 26 - t, T = 10):
  curved   mu = 0.25, 1, 4 (M = 4, 8, 16) on the SAME curved sprinkling (paired);
  flat     the same masses on a flat sprinkling of the same seed (g = 0);
  null     mu = 1 on the curved sprinkling with every term divided by Omega^2(y) (the sum measures coordinate area, not
           the number of events; uses coordinates; a null model);
  up       mu = 1 launched at v = 0.3 (coordinate velocity) towards +x, against the field (reported);
  I4       on the curved order: linear extensions t / t + 0.5 x, Fenwick / brute force, two regions.
Launch: the kink at rest at x0 = 0 with its profile in proper distance s = int Omega dx; phi_0 = (g(x - t) + g(x + t))/2
(d'Alembert with phi_t = 0; for "up" the time derivative of the boosted profile is included numerically).
Readout: the kink fit of F3.14a per slab of width 1 (centres 0.5 ... 9.5), tracking the centre; D = x_c(9.5) - x_c(0.5).
References: the geodesic d2x/dt2 = -Phi'(x)(1 - (dx/dt)^2) from rest at x0, and the continuum field equation
phi_tt - phi_xx + Omega^2 mu sin(phi) = 0 solved by finite differences with the same initial data and the same centre
readout (the extended body, its breathing and radiation included).

Run from the root of the laboratory (numba required):
    python -m diag.f52a_fall --quick
    python -m diag.f52a_fall --seeds 16 --rhos 200,400 --workers 8
Output: summary on screen and raw records in diag/out/f52a_fall_{stamp}.json.
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
G, L, X0 = 0.05, 10.0, 0.0
XB, TB = 26.0, 10.0
MUS = (0.25, 1.0, 4.0)


def Phi(x, g=G):
    return g * L * np.tanh(x / L)


def dPhi(x, g=G):
    return g / np.cosh(x / L) ** 2


def sprinkle_curved(rho, X, T, g, seed):
    """FLAWED (found 10 Nov after the run of F5.2a; kept as the record of these development rows): the points are placed in
    the whole box and those outside the trapezoid dropped, so the density is 0.808 rho exp(2 Phi) for X 26, T 10 while
    the rule normalises with rho. See phases/F5_gravity/field1d_curved.sprinkle_curved2."""
    rng = np.random.default_rng(seed)
    pmax = math.exp(2.0 * g * L) if g > 0 else 1.0
    area = 2.0 * X * T - T * T
    n = int(rng.poisson(rho * pmax * area))
    t = rng.uniform(0.0, T, n)
    x = rng.uniform(-X, X, n)
    keep = (np.abs(x) <= X - t) & (rng.uniform(0.0, 1.0, n) < np.exp(2.0 * Phi(x, g)) / pmax)
    t, x = t[keep], x[keep]
    o = np.argsort(t, kind="stable")
    return t[o], x[o]


_S = {}


def proper(xq, g):
    """s(x) = int_0^x Omega dx' (table)."""
    if g not in _S:
        xs = np.linspace(-80.0, 80.0, 160001)
        om = np.exp(Phi(xs, g))
        s = np.concatenate([[0.0], np.cumsum(0.5 * (om[1:] + om[:-1]) * np.diff(xs))])
        s -= np.interp(0.0, xs, s)
        _S[g] = (xs, s)
    xs, s = _S[g]
    return np.interp(xq, xs, s)


def profile(xq, mu, g, x0=X0):
    return 4.0 * np.arctan(np.exp(math.sqrt(mu) * (proper(xq, g) - proper(np.array([x0]), g)[0])))


def free_wave(t, x, mu, g, v0=0.0):
    """d'Alembert with g(z) = profile(z), h(z) = -v0 profile'(z) (a profile translated with coordinate velocity v0)."""
    gm, gp = profile(x - t, mu, g), profile(x + t, mu, g)
    return 0.5 * (gm + gp) - 0.5 * v0 * (gp - gm)


if _nb is not None:
    @_nb.njit(cache=True)
    def _fenwick_weighted(U, V, order, phi0, cw, rho, mu):
        """F.field1d._fenwick_evolve (conservative, no kick) with every term multiplied by cw[y]."""
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
        tree = np.zeros(off[N + 1] + N + 1, np.float64)
        phi = np.empty(N)
        for a in range(N):
            i = order[a]
            s = 0.0
            k = ru[i] - 1
            while k > 0:
                lo, hi = off[k], off[k + 1]
                p = np.searchsorted(lst[lo:hi], rv[i])
                base = lo + k
                q = p
                while q > 0:
                    s += tree[base + q - 1]
                    q -= q & (-q)
                k -= k & (-k)
            ph = phi0[i] - 0.5 * s / rho
            phi[i] = ph
            w = cw[i] * mu * math.sin(ph)
            k = ru[i]
            while k <= N:
                lo, hi = off[k], off[k + 1]
                p = np.searchsorted(lst[lo:hi], rv[i]) + 1
                base = lo + k
                n_k = hi - lo
                q = p
                while q <= n_k:
                    tree[base + q - 1] += w
                    q += q & (-q)
                k += k & (-k)
        return phi


def track(t, x, phi, mu, x0, T):
    xc = []
    for a in np.arange(0.0, T - 1.0 + 1e-9, 1.0):
        guess = xc[-1] if (xc and np.isfinite(xc[-1])) else x0
        r = F.fit_kink(t, x, phi, a, a + 1.0, guess, mu)
        xc.append(r["xc"] if r else float("nan"))
    return xc


def geodesic(g=G, x0=X0, v0=0.0, T=TB, n=20001):
    tt = np.linspace(0.0, T, n)
    dt = tt[1] - tt[0]
    x, v = x0, v0
    out = np.empty(n)
    for i in range(n):                                   # RK4
        out[i] = x
        f = lambda xx, vv: (vv, -dPhi(xx, g) * (1.0 - vv * vv))
        k1 = f(x, v)
        k2 = f(x + 0.5 * dt * k1[0], v + 0.5 * dt * k1[1])
        k3 = f(x + 0.5 * dt * k2[0], v + 0.5 * dt * k2[1])
        k4 = f(x + dt * k3[0], v + dt * k3[1])
        x += dt / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        v += dt / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
    return tt, out


def continuum(mu, g=G, x0=X0, v0=0.0, T=TB, dx=0.01, X=60.0):
    """phi_tt - phi_xx + Omega^2 mu sin phi = 0, leapfrog; centre per slab = mean over the slab of the pi crossing."""
    xs = np.arange(-X, X + dx / 2, dx)
    dt = 0.5 * dx
    om2 = np.exp(2.0 * Phi(xs, g))
    p0 = profile(xs, mu, g, x0)
    pt = -v0 * np.gradient(p0, dx)
    lap = lambda p: np.concatenate([[0.0], (p[2:] - 2 * p[1:-1] + p[:-2]) / dx ** 2, [0.0]])
    prev = p0.copy()
    cur = p0 + dt * pt + 0.5 * dt * dt * (lap(p0) - om2 * mu * np.sin(p0))
    nsteps = int(round(T / dt))
    cen = np.zeros(int(T))
    cnt = np.zeros(int(T))
    for n_ in range(1, nsteps):
        tnow = n_ * dt
        j = np.nonzero((cur[:-1] - math.pi) * (cur[1:] - math.pi) <= 0)[0]
        if j.size:
            j0 = j[np.argmin(np.abs(xs[j] - (cen[int(tnow) - 1] / max(cnt[int(tnow) - 1], 1) if int(tnow) > 0 else x0)))]
            xcross = xs[j0] + (math.pi - cur[j0]) / (cur[j0 + 1] - cur[j0]) * dx
            if int(tnow) < int(T):
                cen[int(tnow)] += xcross
                cnt[int(tnow)] += 1
        nxt = 2 * cur - prev + dt * dt * (lap(cur) - om2 * mu * np.sin(cur))
        prev, cur = cur, nxt
    return (cen / np.maximum(cnt, 1)).tolist()


def task(rho, seed, arm, quick):
    T = TB
    X = 18.0 if quick else XB
    t0 = time.time()
    kind, mu = arm.split("_")
    mu = float(mu)
    g = 0.0 if kind == "flat" else G
    t, x = sprinkle_curved(rho, X, T, g, seed)
    v0 = 0.3 if kind == "up" else 0.0
    phi0 = free_wave(t, x, mu, g, v0)
    if kind == "null":
        U, V = (t + x) / math.sqrt(2.0), (t - x) / math.sqrt(2.0)
        phi = _fenwick_weighted(U, V, np.arange(t.size, dtype=np.int64), phi0, np.exp(-2.0 * Phi(x, g)),
                                float(rho), mu)
    elif kind == "i4":
        ref = F.evolve(t, x, phi0, rho, mu=mu)
        sc = float(np.max(np.abs(ref)))
        alt = F.evolve(t, x, phi0, rho, mu=mu, order="boost")
        blk = F.evolve(t, x, phi0, rho, mu=mu, engine="block", block=256)
        sub = np.abs(x) <= (X - 2.0) - t
        phs = F.evolve(t[sub], x[sub], phi0[sub], rho, mu=mu)
        return {"rho": rho, "seed": seed, "arm": arm, "N": int(t.size), "order": float(np.max(np.abs(alt - ref)) / sc),
                "engine": float(np.max(np.abs(blk - ref)) / sc), "region": float(np.max(np.abs(phs - ref[sub])) / sc)}
    else:
        phi = F.evolve(t, x, phi0, rho, mu=mu)
    xc = track(t, x, phi, mu, X0, T)
    return {"rho": rho, "seed": seed, "arm": arm, "N": int(t.size), "xc": xc, "sec": time.time() - t0}


def run_task(*a):
    return task(*a)


def _ms(a):
    a = np.array([v for v in a if np.isfinite(v)], float)
    return (float(a.mean()), float(a.std(ddof=1) / math.sqrt(a.size)), float(a.std(ddof=1)), a.size) if a.size > 1 \
        else (float("nan"),) * 3 + (a.size,)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=16)
    ap.add_argument("--rhos", type=str, default="200,400")
    ap.add_argument("--workers", type=int, default=0)
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--base", type=int, default=20261110)
    a = ap.parse_args()
    if _nb is None:
        raise SystemExit("numba required")
    rhos = [float(r) for r in a.rhos.split(",")]
    ns = a.seeds
    if a.quick:
        rhos, ns = [100.0], 2
    arms = [f"curved_{m:g}" for m in MUS] + [f"flat_{m:g}" for m in MUS] + ["null_1", "up_1"]
    jobs = [(rho, derive_seed(a.base, "F5.2a-dev", rho, k), arm, a.quick) for rho in rhos for k in range(ns)
            for arm in arms]
    jobs += [(100.0, derive_seed(a.base, "F5.2a-dev-i4", k), "i4_1", True) for k in range(2)]
    t0 = time.time()
    res = parallel_map(run_task, jobs, default_workers(a.workers), print,
                       lambda i, r, d, n: f"[F5.2a-dev] {d}/{n}: {r['arm']} rho {r['rho']:g} N {r['N']} "
                                          f"({time.time() - t0:.0f}s)")
    tt, geo = geodesic()
    tm = np.arange(0.5, TB, 1.0)
    g_at = np.interp(tm, tt, geo)
    Dgeo = g_at[-1] - g_at[0]
    tt2, geo_up = geodesic(v0=0.3)
    gu = np.interp(tm, tt2, geo_up)
    L_ = [f"geodesic from rest: x(t_m) - x0 = " + " ".join(f"{v:+.3f}" for v in g_at - X0) + f"; D_geo {Dgeo:+.4f}",
          f"geodesic launched at v 0.3: D {gu[-1] - gu[0]:+.4f}"]
    cont = {m: continuum(m) for m in MUS}
    for m in MUS:
        c = np.array(cont[m])
        L_.append(f"continuum field equation mu {m:g}: centres " + " ".join(f"{v:+.3f}" for v in c) +
                  f"; D_cont (slabs 1..9) {c[-1] - c[0]:+.4f}  [slab 0 is t in [0, 1)]")
    i4 = [r for r in res if r["arm"] == "i4_1"]
    L_.append("I4 curved (order / engine / region): " + ", ".join(f"{r['order']:.1e} / {r['engine']:.1e} / "
                                                                    f"{r['region']:.1e}" for r in i4))
    for rho in rhos:
        L_.append(f"--- rho {rho:g}")
        R = {(r["seed"], r["arm"]): r for r in res if r["rho"] == rho and r["arm"] != "i4_1"}
        seeds = sorted({s for s, _ in R})
        D = {arm: [R[(s, arm)]["xc"][-1] - R[(s, arm)]["xc"][0] for s in seeds] for arm in arms}
        for arm in arms:
            m = _ms(D[arm])
            extra = ""
            if arm.startswith("curved"):
                mu = float(arm.split("_")[1])
                c = np.array(cont[mu])
                extra = f"; a = D / D_geo {m[0] / Dgeo:.3f} +- {m[1] / abs(Dgeo):.3f}; D / D_cont " \
                        f"{m[0] / (c[-1] - c[0]):.3f} +- {m[1] / abs(c[-1] - c[0]):.3f}"
            L_.append(f"  {arm:10s}: D {m[0]:+.4f} +- {m[1]:.4f} (sd per sprinkling {m[2]:.3f}, n {m[3]}){extra}")
        for i, j in ((0, 1), (1, 2), (0, 2)):
            ai, aj = f"curved_{MUS[i]:g}", f"curved_{MUS[j]:g}"
            eta = _ms([di / dj - 1.0 for di, dj in zip(D[ai], D[aj])])
            dd = _ms([di - dj for di, dj in zip(D[ai], D[aj])])
            L_.append(f"  WEP M {8 * math.sqrt(MUS[i]):g} vs {8 * math.sqrt(MUS[j]):g}: eta {eta[0]:+.4f} +- {eta[1]:.4f} "
                      f"(sd {eta[2]:.3f}); paired D_i - D_j {dd[0]:+.4f} +- {dd[1]:.4f}")
    print("\n".join(L_))
    os.makedirs(os.path.join("diag", "out"), exist_ok=True)
    fn = os.path.join("diag", "out", f"f52a_fall_{time.strftime('%Y%m%dT%H%M%S', time.gmtime())}.json")
    with open(fn, "w") as fh:
        json.dump({"args": vars(a), "summary": L_, "results": res, "continuum": cont}, fh)
    print(f"raw records: {fn} ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
