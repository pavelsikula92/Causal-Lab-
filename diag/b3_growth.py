"""
B3, development: the causal growth of a 2+1 order whose density and cones follow the potential of matter (B3a)
=============================================================================================================
Decisions 16, 19, 22 and the author's choice of 13 Nov (B3a; the active side; the fall of a test body moved to B4).

Candidates: a Poisson sprinkling of the coordinate pyramid |x| <= X - t, |y| <= Y - t at density rho e^{-Lmin}, each
with an independent mark m. In time order, every candidate v
  1. finds its parents among the REALISED earlier events by the cone rule of B3a (the rest frame of the static source
     wall = the coordinate frame; a declared proxy as F7.1b): u < v iff |x_v - x_u|^2 + |y_v - y_u|^2 < (Pi_u^2 dt)^2
     (the coordinate speed of light Pi^2 of the metric of the theory, read at the parent); links = the maximal ones
     among the candidates of the shell tau^2 = Pi_u^2 dt^2 - d^2/Pi_u^2 <= tau_cut^2 (the proper time at the parent);
  2. gets its longest chains from every earlier realised event through the links (backward dynamic programme);
  3. gets the field and the potential from its realised past with the chain weight of B0 (variant B):
       phi(v) = phi0(v) - (1/rho) sum w(L_uv) mu sin phi(u)
       L(v)   = -(kappa/rho) sum w(L_uv) q(u) - (kappa/2) s(t_v) f(v),   f = phi (2 pi - phi) / 2
     q(u) = s(t_u) [3V + (pi - phi) mu sin(phi) / 2] - s''(t_u) f(u) / 2 [- s'(t_u) d_t f, zero for the static kink]
     (the load = trace of T = 3V + X/2 by the identity of decision 19 with the quadratic f, f'' = -1, that vanishes on
     both vacua: X = (pi - phi) box phi - box f; box f is moved onto the Green function, with the switch s(t) of
     gravity, smoothstep from t1 to t2. The f = phi^2 of decision 19 is 4 pi^2 on the right plateau, and the cancellation
     of its retarded sum amplified the calibration error of the weight: development, 13 Nov, r_L 1.5);
  4. is realised iff m < exp(Lmin - L(v)), i.e. at coordinate density rho Pi^-1 = rho sqrt(-g) of the metric
     -Pi^2 dt^2 + Pi^-2 (dx^2 + dy^2) (more events near matter).
phi0 is the FLAT free wave of the kink data; it is exact because gravity is switched on (t1 = 3) after the fronts of the
free wave have left the wall: where the order is curved, phi0 is constant (diag/b3_continuum.py, variant rule_s: the
fall of a test body within 1-2 % of the exact continuum).
Continuum reference: diag/b3_continuum.solve(variant "exact", the same switch).

    python -m diag.b3_growth --rho 30 --seeds 2 --workers 2
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import time

import numpy as np

from core.order.sprinkle import derive_seed
from core.stats.parallel import parallel_map, default_workers
from phases.F3_structure import field1d as F
from diag import b0_chain_field as B0
from diag.b3_continuum import Bodies, solve as csolve

_nb = F._nb
CALIB = (0.2732723632276477, 0.21313211596642687, 1.2025922117299621)
P0 = {"T": 9.0, "dX": 6.0, "dY": 2.5, "kappa": 0.01, "t1": 3.0, "t2": 4.0, "Lmin": -0.7, "band": 2.0,
      "bt": 0.25, "bx": 0.5, "by": 2.0, "mid": True}


def sprinkle_pyramid(rho, X, Y, T, seed):
    rng = np.random.default_rng(seed)
    n = int(rng.poisson(rho * T * 4 * X * Y))
    t = rng.uniform(0, T, n)
    x = rng.uniform(-X, X, n)
    y = rng.uniform(-Y, Y, n)
    keep = (np.abs(x) <= X - t) & (np.abs(y) <= Y - t)
    t, x, y = t[keep], x[keep], y[keep]
    o = np.argsort(t, kind="stable")
    return t[o], x[o], y[o]


if _nb is not None:
    @_nb.njit(cache=True, inline="always")
    def _ramp(t, t1, t2):
        if t <= t1:
            return 0.0, 0.0, 0.0
        if t >= t2:
            return 1.0, 0.0, 0.0
        h = t2 - t1
        u = (t - t1) / h
        return u * u * (3.0 - 2.0 * u), 6.0 * u * (1.0 - u) / h, (6.0 - 12.0 * u) / (h * h)

    @_nb.njit(cache=True, inline="always")
    def _c2(tu, xu, yu, tv, xv, yv, Piu, Lsum, Lcnt, X, Y, bt, bx, by):
        """Pi^2 at the coordinate midpoint of the pair from the binned potential of the realised events (the cone of a
        long near-null pair is not the cone at its parent when the potential deepens in time); fallback Pi_u^2."""
        it = int(0.5 * (tu + tv) / bt)
        ix = int((0.5 * (xu + xv) + X) / bx)
        iy = int((0.5 * (yu + yv) + Y) / by)
        if it < Lsum.shape[0] and 0 <= ix < Lsum.shape[1] and 0 <= iy < Lsum.shape[2] and Lcnt[it, ix, iy] > 0:
            return math.exp(2.0 * Lsum[it, ix, iy] / Lcnt[it, ix, iy])
        return Piu * Piu

    @_nb.njit(cache=True)
    def grow(t, x, y, mark, phi0, rho, mu, kappa, wL, tc2, Lmin, t1, t2, X, Y, bt, bx, by, mid):
        N = t.size
        nt = int(t[N - 1] / bt) + 2
        Lsum = np.zeros((nt, int(2 * X / bx) + 1, int(2 * Y / by) + 1))
        Lcnt = np.zeros((nt, int(2 * X / bx) + 1, int(2 * Y / by) + 1), np.int64)
        real = np.zeros(N, np.bool_)
        Pi = np.ones(N)
        phi = np.zeros(N)
        Lv = np.zeros(N)
        q = np.zeros(N)
        rl = np.empty(N, np.int64)
        nr = 0
        ls = np.zeros(N, np.int64)
        le = np.zeros(N, np.int64)
        cap = 64 * N
        links = np.empty(cap, np.int64)
        tot = 0
        best = np.zeros(N, np.int64)
        stamp = np.full(N, -1, np.int64)
        cand = np.empty(N, np.int64)
        maxi = np.empty(N, np.int64)
        over = 0
        for v in range(N):
            m = 0
            for a in range(nr - 1, -1, -1):
                u = rl[a]
                dt = t[v] - t[u]
                if dt <= 0.0:
                    continue
                d2 = (x[v] - x[u]) ** 2 + (y[v] - y[u]) ** 2
                if d2 >= dt * dt:
                    continue
                c2 = _c2(t[u], x[u], y[u], t[v], x[v], y[v], Pi[u], Lsum, Lcnt, X, Y, bt, bx, by) if mid else Pi[u] * Pi[u]
                if d2 < c2 * c2 * dt * dt:
                    if c2 * dt * dt - d2 / c2 <= tc2:
                        cand[m] = u
                        m += 1
            nm = 0
            for i in range(m):
                k = cand[i]
                dead = False
                for j in range(nm):
                    mm = maxi[j]
                    dt = t[mm] - t[k]
                    if dt > 0.0:
                        d2 = (x[mm] - x[k]) ** 2 + (y[mm] - y[k]) ** 2
                        c2 = _c2(t[k], x[k], y[k], t[mm], x[mm], y[mm], Pi[k], Lsum, Lcnt, X, Y, bt, bx, by) if mid \
                            else Pi[k] * Pi[k]
                        if d2 < c2 * c2 * dt * dt:
                            dead = True
                            break
                if not dead:
                    maxi[nm] = k
                    nm += 1
            S1 = 0.0
            S2 = 0.0
            for j in range(nm):
                best[maxi[j]] = 1
                stamp[maxi[j]] = v
            for a in range(nr - 1, -1, -1):
                u = rl[a]
                if stamp[u] != v:
                    continue
                b = best[u]
                S1 += wL[b] * mu * math.sin(phi[u])
                S2 += wL[b] * q[u]
                for e in range(ls[u], le[u]):
                    p = links[e]
                    if stamp[p] != v or best[p] < b + 1:
                        best[p] = b + 1
                        stamp[p] = v
            ph = phi0[v] - S1 / rho
            s, s1, s2 = _ramp(t[v], t1, t2)
            fv = 0.5 * ph * (2.0 * math.pi - ph)
            L = -kappa * S2 / rho - 0.5 * kappa * s * fv
            Lv[v] = L
            phi[v] = ph
            if L < Lmin:
                over += 1
            if mark[v] < math.exp(Lmin - L):
                real[v] = True
                Pi[v] = math.exp(L)
                V = mu * (1.0 - math.cos(ph))
                q[v] = s * (3.0 * V + 0.5 * (math.pi - ph) * mu * math.sin(ph)) - 0.5 * s2 * fv
                if tot + nm > cap:
                    new = np.empty(2 * cap, np.int64)
                    new[:tot] = links[:tot]
                    links = new
                    cap *= 2
                ls[v] = tot
                for j in range(nm):
                    links[tot] = maxi[j]
                    tot += 1
                le[v] = tot
                rl[nr] = v
                nr += 1
                it = int(t[v] / bt)
                ix = int((x[v] + X) / bx)
                iy = int((y[v] + Y) / by)
                Lsum[it, ix, iy] += L
                Lcnt[it, ix, iy] += 1
        return real, phi, Lv, Pi, q, ls, le, links[:tot], over

    @_nb.njit(cache=True)
    def future_reach(src, real, t, x, y, ls, le, links, dt_lo, dt_hi, dxmax):
        """For each source event: the largest |dy| among realised events in its causal future (reachability through the
        links) with t - t_src in [dt_lo, dt_hi) and |dx| < dxmax; also the longest chain to the event straight above
        (|dx|, |dy| < dxmax, t - t_src in [dt_lo, dt_hi)) -- returns (reach, chainlen)."""
        N = t.size
        # children lists (forward adjacency) built once
        cnt = np.zeros(N + 1, np.int64)
        for v in range(N):
            if real[v]:
                for e in range(ls[v], le[v]):
                    cnt[links[e] + 1] += 1
        for i in range(N):
            cnt[i + 1] += cnt[i]
        ch = np.empty(cnt[N], np.int64)
        fill = cnt[:-1].copy()
        for v in range(N):
            if real[v]:
                for e in range(ls[v], le[v]):
                    u = links[e]
                    ch[fill[u]] = v
                    fill[u] += 1
        reach = np.zeros(src.size)
        clen = np.full(src.size, -1.0)
        dist = np.full(N, -1, np.int64)
        for s_i in range(src.size):
            s0 = src[s_i]
            dist[:] = -1
            dist[s0] = 0
            best = 0.0
            bl = -1
            for v in range(s0, N):
                if dist[v] < 0:
                    continue
                if t[v] - t[s0] >= dt_hi:
                    continue
                for e in range(cnt[v], cnt[v + 1]):
                    c = ch[e]
                    if dist[c] < dist[v] + 1:
                        dist[c] = dist[v] + 1
                d = t[v] - t[s0]
                if d >= dt_lo and abs(x[v] - x[s0]) < dxmax:
                    if abs(y[v] - y[s0]) > best:
                        best = abs(y[v] - y[s0])
                    if abs(y[v] - y[s0]) < dxmax and dist[v] > bl:
                        bl = dist[v]
            reach[s_i] = best
            clen[s_i] = bl
        return reach, clen


def task(rho, muA, kappa, k, base, P=P0, label="B3dev", n_probe=12):
    T = P["T"]
    X, Y = T + P["dX"], T + P["dY"]
    Lmin = P["Lmin"]
    sd = derive_seed(base, label, rho, k)
    t, x, y = sprinkle_pyramid(rho * math.exp(-Lmin), X, Y, T, sd)
    mark = np.random.default_rng(derive_seed(sd, "marks")).uniform(0, 1, t.size)
    f = lambda z: 4.0 * np.arctan(np.exp(np.clip(math.sqrt(muA) * z, -50, 50)))
    phi0 = 0.5 * (f(x - t) + f(x + t))
    wL = B0.weight_table("B", rho, CALIB)
    tc = 1.15 * (12.0 * 25.0 / (math.pi * rho)) ** (1.0 / 3.0)
    a0 = time.time()
    real, phi, L, Pi, q, ls, le, links, over = grow(t, x, y, mark, phi0, float(rho), float(muA), float(kappa), wL,
                                                     tc * tc, Lmin, P["t1"], P["t2"], X, Y, P["bt"], P["bx"], P["by"],
                                                     P["mid"])
    sec = time.time() - a0
    xg, tg, PH, LL = csolve(Bodies([("kink", 0.0, muA, 0.0)]), kappa, T=T, dx=0.01, X=X + 10, variant="exact",
                            t1=P["t1"], t2=P["t2"])

    def cont(arr, tt, xx):
        ft = np.clip(tt / (tg[1] - tg[0]), 0, tg.size - 1.000001)
        it = ft.astype(int)
        at = ft - it
        fx = np.clip((xx - xg[0]) / (xg[1] - xg[0]), 0, xg.size - 1.000001)
        ix = fx.astype(int)
        ax = fx - ix
        return ((1 - at) * ((1 - ax) * arr[it, ix] + ax * arr[it, ix + 1])
                + at * ((1 - ax) * arr[it + 1, ix] + ax * arr[it + 1, ix + 1]))

    Lc = cont(LL, t, x)
    phic = cont(PH[:, 0], t, x)
    R = real
    band = R & (np.abs(y) < P["band"])
    rows = []
    for a in np.arange(0.0, T - 2.0 + 1e-9, 1.0):
        s = band & (t >= a) & (t < a + 1.0)
        rk = F.fit_kink(t[s], x[s], phi[s], a, a + 1.0, 0.0, 1.0 / math.sqrt(muA))
        rx = F.fit_kink(t[s], x[s], phic[s], a, a + 1.0, 0.0, 1.0 / math.sqrt(muA))
        reg = s & (np.abs(x) < 6.0)
        num = float(np.sum(L[reg] * Lc[reg]))
        den = float(np.sum(Lc[reg] ** 2))
        bins = []
        for xb in np.arange(-6.0, 6.01, 1.0):
            m = reg & (np.abs(x - xb) < 0.5)
            bins.append([float(np.mean(L[m])) if m.any() else float("nan"),
                         float(np.mean(Lc[m])) if m.any() else float("nan"), int(m.sum())])
        rows.append({"t": a + 0.5, "xc": rk["xc"] if rk else float("nan"), "xc_ex": rx["xc"] if rx else float("nan"),
                     "w": rk["w"] if rk else float("nan"), "w_ex": rx["w"] if rx else float("nan"),
                     "LL": num, "LcLc": den, "bins": bins})
    # density: realised against the acceptance probability of every candidate
    acc = np.exp(np.minimum(Lmin - L, 0.0))
    late = t > P["t2"]
    dens = [float(R[late].sum()), float(acc[late].sum())]
    # clocks and cones: sources in the band at t in [4.5, 5), near the wall and far
    rng = np.random.default_rng(derive_seed(sd, "probe"))
    out_probe = {}
    for name, x0 in (("near", 0.0), ("far", 4.5)):
        cand = np.nonzero(R & (np.abs(y) < 1.0) & (np.abs(x - x0) < 0.5) & (t >= 4.5) & (t < 5.0))[0]
        src = np.sort(rng.choice(cand, size=min(n_probe, cand.size), replace=False)) if cand.size else cand
        if src.size:
            reach, clen = future_reach(src.astype(np.int64), R, t, x, y, ls, le, links, 2.5, 3.0, 0.4)
            Pic = np.exp(cont(LL, t[src] + 1.4, x[src]))
            out_probe[name] = {"reach": reach.tolist(), "clen": clen.tolist(), "Pi2_mid": (Pic ** 2).tolist(),
                               "tau_cont": [float(np.mean(np.exp(cont(LL, t[s_] + np.linspace(0, 2.75, 12), np.full(12, x[s_]))))) * 2.75
                                            for s_ in src]}
    return {"rho": rho, "muA": muA, "kappa": kappa, "k": k, "Ncand": int(t.size), "N": int(R.sum()),
            "links": int(links.size), "over": int(over), "sec": sec, "rows": rows, "dens": dens, "probe": out_probe,
            "Lmin_seen": float(L[R].min())}


def cached(rho, muA, kappa, k, base):
    key = hashlib.sha256(json.dumps([rho, muA, kappa, k, base, P0, CALIB, "v3"]).encode()).hexdigest()[:20]
    fn = os.path.join("diag", "out", "b3_dev_cache", key + ".json")
    if os.path.exists(fn):
        return json.load(open(fn))
    out = task(rho, muA, kappa, k, base)
    os.makedirs(os.path.dirname(fn), exist_ok=True)
    json.dump(out, open(fn + ".tmp", "w"))
    os.replace(fn + ".tmp", fn)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rho", type=float, default=30.0)
    ap.add_argument("--seeds", type=int, default=2)
    ap.add_argument("--k0", type=int, default=0)
    ap.add_argument("--mus", type=str, default="1.0,0.25")
    ap.add_argument("--kappas", type=str, default="0.01,0.0")
    ap.add_argument("--workers", type=int, default=0)
    a = ap.parse_args()
    jobs = [(a.rho, float(mu), float(kp), k, 20261113) for k in range(a.k0, a.k0 + a.seeds)
            for mu in a.mus.split(",") for kp in a.kappas.split(",")]
    t0 = time.time()
    res = parallel_map(cached, jobs, default_workers(a.workers), print,
                       lambda i, r, d, n: f"[B3 dev] {d}/{n}: mu {r['muA']} kappa {r['kappa']} k {r['k']} N {r['N']} "
                                          f"(cand {r['Ncand']}) {r['sec']:.0f}s ({time.time() - t0:.0f}s)")
    fn = os.path.join("diag", "out", f"b3_dev_{time.strftime('%Y%m%dT%H%M%S', time.gmtime())}.json")
    json.dump({"args": vars(a), "P": P0, "results": res}, open(fn, "w"))
    print("saved", fn)


if __name__ == "__main__":
    main()
