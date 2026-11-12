"""
D2 / B-b, development stage 2: the kinetic part of the load without a gradient (decision 19; decides nothing)
=============================================================================================================
Signature (+, -, -), box = d_t^2 - grad^2, field rule phi = phi_0 - G*J with box phi = -J, J = V'(phi).
The load is the trace T = -1/2 (d phi)^2 + 3V = -K/4 + 3V with K = 2 (d phi)^2. In the growth rule the load enters only
through its Green integral (L = -kappa G*T). From box(phi^2) = 2 phi box phi + K = -2 phi J + K:

    G*K = phi^2 - chi_0 + 2 G*(phi J),

chi_0 the free (massless) wave of the initial data of phi^2 (as phi_0 is that of phi). So the Green integral of the
kinetic term needs no gradient and no frame: phi^2 at the event, the free wave of the data, and one more Green sum with
the same weights as the field. Then G*T = -G*K / 4 + 3 G*V.

Checked here on the intrinsic field (variant B of diag/b0_chain_field.py, the table of decision 17), at targets near the
centre line of the configuration, against
  (i)  the same discrete Green sum of the EXACT kinetic term evaluated at the events (same weights, same sprinkling):
       isolates the identity from the weights;
  (ii) the continuum Green integral of the exact kinetic term (and of V, and of T).
Arms: the wall at rest, the wall at v = 0.3 (sine-Gordon, exact), a plane wave k = (1, 0) of amplitude 0.1 (sine-Gordon;
the linear wave is exact to ~a^2/8 ~ 0.1 %).

    python -m diag.d2_load --rhos 25,50 --seeds 4 --workers 2
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
from diag.b0_green2p1 import sprinkle_pyramid, modes_free_2d
from diag.b0_chain_alg import links3
from diag.b0_chain_field import weight_table, newest_calib

_nb = F._nb
AMP = 0.1

if _nb is not None:
    @_nb.njit(cache=True)
    def evolve_load(ptr, idx, phi0, rho, mu, wL, Kex):
        """The field rule (J = mu sin phi) with three more Green sums on the same chains:
        S_pJ = G*(phi J), S_V = G*(V), S_K = G*(K_exact); V = mu (1 - cos phi)."""
        N = phi0.size
        acc = np.zeros(N)
        SpJ = np.zeros(N)
        SV = np.zeros(N)
        SK = np.zeros(N)
        phi = np.empty(N)
        L = np.zeros(N, np.int64)
        stamp = np.full(N, -1, np.int64)
        for v in range(N):
            ph = phi0[v] - acc[v]
            phi[v] = ph
            J = mu * math.sin(ph)
            c1 = J / rho
            c2 = ph * J / rho
            c3 = mu * (1.0 - math.cos(ph)) / rho
            c4 = Kex[v] / rho
            L[v] = 0
            stamp[v] = v
            for u in range(v + 1, N):
                b = -1
                for q in range(ptr[u], ptr[u + 1]):
                    p = idx[q]
                    if p >= v and stamp[p] == v and L[p] + 1 > b:
                        b = L[p] + 1
                if b > 0:
                    L[u] = b
                    stamp[u] = v
                    w = wL[b]
                    acc[u] += w * c1
                    SpJ[u] += w * c2
                    SV[u] += w * c3
                    SK[u] += w * c4
        return phi, SpJ, SV, SK


_GT, _GW = np.polynomial.legendre.leggauss(200)


def green_moving_sech2(t, x, v):
    """Continuum G*f for f = sech^2(g (x - v t)) from t = 0 (y-independent: 1+1 reduction, G -> 1/2 theta)."""
    g = 1.0 / math.sqrt(1 - v * v)
    out = np.empty(t.size)
    for i in range(t.size):
        tp = 0.5 * t[i] * (_GT + 1.0)
        a = np.tanh(g * (x[i] + (t[i] - tp) - v * tp)) - np.tanh(g * (x[i] - (t[i] - tp) - v * tp))
        out[i] = 0.5 * 0.5 * t[i] * np.sum(_GW * a) / g
    return out


def chi0_wall(t, x, v):
    """Free massless evolution of the data of phi^2 for the kink launched at v (d'Alembert in x)."""
    g = 1.0 / math.sqrt(1 - v * v)
    X = np.linspace(-40, 40, 160001)
    ph = 4 * np.arctan(np.exp(g * X))
    f = ph ** 2
    gdat = 2 * ph * (-v * g * 2.0 / np.cosh(g * X))
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (gdat[1:] + gdat[:-1]) * np.diff(X))])
    fi = lambda z: np.interp(z, X, f)
    ci = lambda z: np.interp(z, X, cum)
    return 0.5 * (fi(x + t) + fi(x - t)) + 0.5 * (ci(x + t) - ci(x - t))


def wave_parts(t, x, y, kx, ky, th=0.0):
    k = math.hypot(kx, ky)
    om = math.sqrt(k * k + 1.0)
    p0 = kx * x + ky * y + th
    p = p0 - om * t
    a = AMP
    chi0 = a * a / 2 + a * a / 2 * np.cos(2 * p0) * np.cos(2 * k * t) + a * a * om / (2 * k) * np.sin(2 * p0) * np.sin(2 * k * t)
    K = a * a * (1 - np.cos(2 * p))                      # 2 (d phi)^2 for the linear wave, mu = 1
    # G*cos(2p): u_p = -cos(2p)/4 minus the free wave of its data
    Om = 2 * om
    hfree = -np.cos(2 * p0) * np.cos(2 * k * t) / 4 - (Om / 4) * np.sin(2 * p0) * np.sin(2 * k * t) / (2 * k)
    Gcos = -np.cos(2 * p) / 4 - hfree
    GK = a * a * t * t / 2 - a * a * Gcos
    # V = 1 - cos(phi) ~ phi^2 / 2 = a^2/4 (1 + cos 2p): G*V
    GV = a * a / 4 * (t * t / 2 + Gcos)
    return chi0, K, GK, GV


def task(rho, seed, arm, T, calib):
    X, Y = T + 8.0, T + 3.0
    t, x, y = sprinkle_pyramid(rho, X, Y, T, seed)
    if arm in ("rest", "v0.3"):
        v = 0.0 if arm == "rest" else 0.3
        g = 1.0 / math.sqrt(1 - v * v)
        phi0 = F.kink_free(t, x, v)
        z = g * (x - v * t)
        Kex = -8.0 / np.cosh(z) ** 2
    else:
        kx, ky = 1.0, 0.0
        phi0 = modes_free_2d(t, x, y, [(AMP, kx, ky, 0.0)])
        _, Kex, _, _ = wave_parts(t, x, y, kx, ky)
    t0 = time.time()
    ptr, idx = links3(t, x, y, (12.0 * 25.0 / (math.pi * rho)) ** (1.0 / 3.0))
    phi, SpJ, SV, SK = evolve_load(ptr, idx, phi0, float(rho), 1.0, weight_table("B", rho, calib), Kex)
    sec = time.time() - t0
    rng = np.random.default_rng(derive_seed(seed, "targets", arm))
    if arm in ("rest", "v0.3"):
        cand = np.nonzero((t > T - 1.5) & (t < T - 0.5) & (np.abs(x - v * t) < 3.0) & (np.abs(y) < 2.0))[0]
    else:
        cand = np.nonzero((t > T - 1.5) & (t < T - 0.5) & (np.abs(x) < 4.0) & (np.abs(y) < 2.0))[0]
    tg = np.sort(rng.choice(cand, min(400, cand.size), replace=False))
    tt, xx, yy = t[tg], x[tg], y[tg]
    if arm in ("rest", "v0.3"):
        chi0 = chi0_wall(tt, xx, v)
        Gs = green_moving_sech2(tt, xx, v)
        GKc, GVc = -8.0 * Gs, 2.0 * Gs
    else:
        chi0, _, GKc, GVc = wave_parts(tt, xx, yy, 1.0, 0.0)
    GKid = phi[tg] ** 2 - chi0 + 2.0 * SpJ[tg]
    return {"rho": rho, "seed": seed, "arm": arm, "N": int(t.size), "sec": sec, "x": (xx - (v * tt if arm != "k1_0" else 0)).tolist(),
            "GKid": GKid.tolist(), "GKdir": SK[tg].tolist(), "GKc": GKc.tolist(), "GV": SV[tg].tolist(), "GVc": GVc.tolist()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=4)
    ap.add_argument("--rhos", type=str, default="25,50")
    ap.add_argument("--T", type=float, default=4.0)
    ap.add_argument("--workers", type=int, default=0)
    a = ap.parse_args()
    if _nb is None:
        raise SystemExit("numba required")
    ca, cc, cb, fn_cal = newest_calib()
    rhos = [float(r) for r in a.rhos.split(",")]
    arms = ["rest", "v0.3", "k1_0"]
    jobs = [(rho, derive_seed(20261111, "D2load", rho, k), arm, a.T, (ca, cc, cb)) for rho in rhos for k in range(a.seeds)
            for arm in arms]
    jobs.sort(key=lambda j: -j[0])
    t0 = time.time()
    res = parallel_map(task, jobs, default_workers(a.workers), print,
                       lambda i, r, d, n: f"[D2 load] {d}/{n}: {r['arm']} rho {r['rho']:g} N {r['N']} {r['sec']:.0f}s "
                                          f"({time.time() - t0:.0f}s)")
    L_ = [f"G*K by the identity phi^2 - chi_0 + 2 G*(phi J) against (i) the same weights on the exact K and (ii) the "
          f"continuum; intrinsic weight B ({fn_cal}); T {a.T:g}; ratio of means over targets (SE over sprinklings), "
          f"rms per target relative to the rms of the continuum value"]
    for rho in rhos:
        L_.append(f"=== rho {rho:g}")
        for arm in arms:
            R = [r for r in res if r["rho"] == rho and r["arm"] == arm]

            def rat(num, den):
                q = [np.sum(r[num]) / np.sum(r[den]) for r in R]
                return np.mean(q), np.std(q, ddof=1) / math.sqrt(len(q))

            def rms(num, den):
                d = np.concatenate([np.array(r[num]) - np.array(r[den]) for r in R])
                s = np.sqrt(np.mean(np.concatenate([np.array(r[den]) for r in R]) ** 2))
                return np.sqrt(np.mean(d ** 2)) / s

            parts = []
            for num, den, lab in (("GKid", "GKc", "identity / continuum"), ("GKdir", "GKc", "direct sum / continuum"),
                                  ("GKid", "GKdir", "identity / direct sum"), ("GV", "GVc", "G*V / continuum")):
                m, s = rat(num, den)
                parts.append(f"{lab} {m:.4f} +- {s:.4f} (rms {rms(num, den):.3f})")
            GT = [(-0.25 * np.array(r["GKid"]) + 3 * np.array(r["GV"])) for r in R]
            GTc = [(-0.25 * np.array(r["GKc"]) + 3 * np.array(r["GVc"])) for r in R]
            q = [np.sum(a_) / np.sum(b_) for a_, b_ in zip(GT, GTc)]
            parts.append(f"G*T (identity route) / continuum {np.mean(q):.4f} +- {np.std(q, ddof=1) / math.sqrt(len(q)):.4f}")
            L_.append(f"  {arm}: " + "; ".join(parts) + f"  [s/run {np.mean([r['sec'] for r in R]):.0f}]")
    print("\n".join(L_))
    os.makedirs(os.path.join("diag", "out"), exist_ok=True)
    fn = os.path.join("diag", "out", f"d2_load_{time.strftime('%Y%m%dT%H%M%S', time.gmtime())}.json")
    with open(fn, "w") as fh:
        json.dump({"args": vars(a), "calib": fn_cal, "summary": L_, "results": res}, fh)
    print(f"raw records: {fn} ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
