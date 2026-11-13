"""
B3 diagnostics of the 5-8 % deficit of the potential in F5.2d (13 Nov; decides nothing)
====================================================================================
A copy of the frozen growth kernel of diag/b3_growth.py (left unchanged, md5 6bada302...) with one option: a FIXED
load q = s(t) * 4V of the flat static kink at the event (no field dynamics in the source, no local term), so the
potential sum over the growing order can be compared with the continuum response to the same prescribed load.
Runs: the frozen source at kappa 0.01 and 0.001 (an almost flat order: the weight alone), at rho 20 / 30 / 45 (the
trend with density), and with finer bins of the midpoint potential; the dynamic source at rho 45.

    python -m diag.b3_diag --workers 2
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
from diag import b3_growth as G

_nb = F._nb
CALIB = G.CALIB

if _nb is not None:
    @_nb.njit(cache=True)
    def grow_d(t, x, y, mark, phi0, rho, mu, kappa, wL, tc2, Lmin, t1, t2, X, Y, bt, bx, by, mid, fixed, qfix):
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
                c2 = G._c2(t[u], x[u], y[u], t[v], x[v], y[v], Pi[u], Lsum, Lcnt, X, Y, bt, bx, by) if mid else Pi[u] * Pi[u]
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
                        c2 = G._c2(t[k], x[k], y[k], t[mm], x[mm], y[mm], Pi[k], Lsum, Lcnt, X, Y, bt, bx, by) if mid \
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
            s, s1, s2 = G._ramp(t[v], t1, t2)
            fv = 0.0 if fixed else 0.5 * ph * (2.0 * math.pi - ph)
            L = -kappa * S2 / rho - 0.5 * kappa * s * fv
            Lv[v] = L
            phi[v] = ph
            if L < Lmin:
                over += 1
            if mark[v] < math.exp(Lmin - L):
                real[v] = True
                Pi[v] = math.exp(L)
                V = mu * (1.0 - math.cos(ph))
                q[v] = s * qfix[v] if fixed else s * (3.0 * V + 0.5 * (math.pi - ph) * mu * math.sin(ph)) - 0.5 * s2 * fv
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



def cont_fixed(mu, kappa, T, X, t1, t2, dx=0.01, dt_out=0.05):
    """Continuum response box_g L = kappa s(t) q(x) to the prescribed load q = 4V of the flat static kink."""
    x = np.arange(-X, X + 0.5 * dx, dx)
    n = x.size
    qx = 4.0 * mu * (1.0 - np.cos(4.0 * np.arctan(np.exp(np.clip(math.sqrt(mu) * x, -50, 50)))))
    L = np.zeros(n)
    M = np.zeros(n)
    dt = 0.3 * dx
    every = int(round(dt_out / dt))
    dt = dt_out / every
    nsteps = int(round(T / dt_out)) * every
    out = [L.copy()]
    t = 0.0
    for s_ in range(1, nsteps + 1):
        sw = _ramp_py(t, t1, t2)[0]
        Pi = np.exp(L)
        Ph = 0.5 * (Pi[1:] + Pi[:-1])
        lap = np.zeros(n)
        lap[1:-1] = (Ph[1:] * (L[2:] - L[1:-1]) - Ph[:-1] * (L[1:-1] - L[:-2])) / (dx * dx)
        M += dt * (lap - kappa * sw * qx / Pi)
        L += dt * Pi ** 3 * M
        L[0], L[-1] = L[1], L[-2]
        t += dt
        if s_ % every == 0:
            out.append(L.copy())
    return x, np.arange(len(out)) * dt_out, np.array(out)


def _ramp_py(t, t1, t2):
    from diag.b3_continuum import ramp
    return ramp(t, t1, t2)


def task(rho, mu, kappa, k, base, fixed, bt, bx, n_probe=0):
    P = dict(G.P0, bt=bt, bx=bx)
    T = P["T"]
    X, Y = T + P["dX"], T + P["dY"]
    Lmin = P["Lmin"]
    sd = derive_seed(base, "B3diag", rho, k)
    t, x, y = G.sprinkle_pyramid(rho * math.exp(-Lmin), X, Y, T, sd)
    mark = np.random.default_rng(derive_seed(sd, "marks")).uniform(0, 1, t.size)
    f = lambda z: 4.0 * np.arctan(np.exp(np.clip(math.sqrt(mu) * z, -50, 50)))
    phi0 = 0.5 * (f(x - t) + f(x + t))
    qfix = 4.0 * mu * (1.0 - np.cos(f(x)))
    wL = B0.weight_table("B", rho, CALIB)
    tc = 1.15 * (12.0 * 25.0 / (math.pi * rho)) ** (1.0 / 3.0)
    a0 = time.time()
    real, phi, L, Pi, q, ls, le, links, over = grow_d(t, x, y, mark, phi0, float(rho), float(mu), float(kappa), wL,
                                                       tc * tc, Lmin, P["t1"], P["t2"], X, Y, bt, bx, P["by"], True,
                                                       bool(fixed), qfix)
    sec = time.time() - a0
    if fixed:
        xg, tg, LL = cont_fixed(mu, kappa, T, X + 10, P["t1"], P["t2"])
    else:
        from diag.b3_continuum import Bodies, solve as csolve
        xg, tg, PH, LL = csolve(Bodies([("kink", 0.0, mu, 0.0)]), kappa, T=T, dx=0.01, X=X + 10, variant="exact",
                                t1=P["t1"], t2=P["t2"])
    ft = np.clip(t / (tg[1] - tg[0]), 0, tg.size - 1.000001)
    it = ft.astype(int)
    at = ft - it
    fx = np.clip((x - xg[0]) / (xg[1] - xg[0]), 0, xg.size - 1.000001)
    ix = fx.astype(int)
    ax = fx - ix
    Lc = ((1 - at) * ((1 - ax) * LL[it, ix] + ax * LL[it, ix + 1]) + at * ((1 - ax) * LL[it + 1, ix] + ax * LL[it + 1, ix + 1]))
    band = real & (np.abs(y) < 2.0) & (np.abs(x) < 6.0)
    rows = []
    for a in np.arange(3.0, 8.0, 0.5):
        s = band & (t >= a) & (t < a + 0.5)
        rows.append({"t": a + 0.25, "LL": float(np.sum(L[s] * Lc[s])), "LcLc": float(np.sum(Lc[s] ** 2)),
                     "near": [float(np.mean(L[s & (np.abs(x) < 1)])), float(np.mean(Lc[s & (np.abs(x) < 1)]))],
                     "far": [float(np.mean(L[s & (np.abs(x) > 4)])), float(np.mean(Lc[s & (np.abs(x) > 4)]))]})
    return {"rho": rho, "mu": mu, "kappa": kappa, "k": k, "fixed": bool(fixed), "bt": bt, "bx": bx, "N": int(real.sum()),
            "sec": sec, "rows": rows}


def cached(*a):
    key = hashlib.sha256(json.dumps([str(q) for q in a] + ["v1"]).encode()).hexdigest()[:20]
    fn = os.path.join("diag", "out", "b3_diag_cache", key + ".json")
    if os.path.exists(fn):
        return json.load(open(fn))
    out = task(*a)
    os.makedirs(os.path.dirname(fn), exist_ok=True)
    json.dump(out, open(fn + ".tmp", "w"))
    os.replace(fn + ".tmp", fn)
    return out


def summary(res):
    groups = {}
    for r in res:
        groups.setdefault((r["fixed"], r["rho"], r["kappa"], r["bt"], r["bx"], r["mu"]), []).append(r)
    for g, R in sorted(groups.items(), key=lambda z: str(z[0])):
        rl = []
        for r in R:
            late = [w for w in r["rows"] if 5.0 <= w["t"] < 8.0]
            rl.append(sum(w["LL"] for w in late) / sum(w["LcLc"] for w in late))
        rl = np.array(rl)
        per = np.mean([[w["LL"] / w["LcLc"] if w["LcLc"] > 0 else np.nan for w in r["rows"]] for r in R], 0)
        print(f"fixed {g[0]!s:5} rho {g[1]:g} kappa {g[2]:g} bins {g[3]:g}/{g[4]:g} mu {g[5]:g}: r_L(5-8) {rl.mean():.4f} +- "
              f"{rl.std(ddof=1) / math.sqrt(rl.size) if rl.size > 1 else float('nan'):.4f} (n {rl.size}); per half-slab t 3.25.. "
              + " ".join(f"{v:.3f}" for v in per))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=0)
    ap.add_argument("--seeds", type=int, default=4)
    a = ap.parse_args()
    jobs = []
    for k in range(a.seeds):
        jobs += [(30.0, 1.0, 0.01, k, 7, True, 0.25, 0.5), (30.0, 1.0, 0.001, k, 7, True, 0.25, 0.5),
                 (30.0, 1.0, 0.01, k, 7, True, 0.125, 0.25), (20.0, 1.0, 0.01, k, 7, True, 0.25, 0.5),
                 (45.0, 1.0, 0.01, k, 7, True, 0.25, 0.5), (30.0, 1.0, 0.01, k, 7, False, 0.25, 0.5),
                 (45.0, 1.0, 0.01, k, 7, False, 0.25, 0.5)]
    t0 = time.time()
    res = parallel_map(cached, jobs, default_workers(a.workers), print,
                       lambda i, r, d, n: f"[B3 diag] {d}/{n}: fixed {r['fixed']} rho {r['rho']} kappa {r['kappa']} "
                                          f"bins {r['bt']}/{r['bx']} k {r['k']} {r['sec']:.0f}s ({time.time() - t0:.0f}s)")
    fn = os.path.join("diag", "out", f"b3_diag_{time.strftime('%Y%m%dT%H%M%S', time.gmtime())}.json")
    json.dump(res, open(fn, "w"))
    print("saved", fn)
    summary(res)


if __name__ == "__main__":
    main()
