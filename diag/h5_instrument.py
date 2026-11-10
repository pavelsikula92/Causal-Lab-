"""
H5-0 development rows: the instrument of H5 in 1+1 (development; decides nothing)
================================================================================

Why (10 Nov 2026; card docs/H5_card_draft.md, revision 3, decisions D1-D4, D2')
-----------------------------------------------------------------------------
H5 tests whether a topological defect of a field of the order has inertia under a conservative dynamics. In 1+1 the
dynamics is the retarded Green function of the order (D2'): the continuum retarded Green function of d_t^2 - d_x^2 is
1/2 inside the past light cone, on a causal set 1/2 times the causal relation. For the equation

    phi_tt - phi_xx + mu sin(phi) = -f chi(t)        (chi = 1 for t1 <= t < t2, else 0)

the field at an event x is

    phi(x) = phi_0(x) - (1 / 2 rho) * sum over events y < x (the order), t_y >= 0, of [mu sin(phi(y)) + f chi(t_y)],

phi_0 = the free wave of the initial data on the slice t = 0 (d'Alembert's formula, exact continuum functions at the
coordinates of the event: the launch, the only place where coordinates enter). The events are computed in any linear
extension of the order. Only events in the domain of dependence of the slice inside |x| <= X are sprinkled (the
trapezoid |x| <= X - t), so every computed value is exact for its order: no wall enters the dynamics.

The card requires the instrument first (author's condition C1). This diagnostic measures what the frozen criteria of
H5-0 will be set from:
  I1  linear stability of the vacuum: small random smooth modes (amplitude 0.01, k in (0.25, 4]); rms of phi against
      the exact linear solution, per slab;
  I2  dispersion: plane waves k in {0.5, 1, 2, 4}; omega from the phase of the projection on cos kx, sin kx in thin
      slabs; against sqrt(k^2 + mu) and against the massless k (the fraction of the mass shift reproduced);
  I3  convergence: the rest kink, rms error against the exact kink, rho in RHOS;
  I4  order and layers: two linear extensions (t and t + 0.5 x), two engines (the Fenwick engine and blocks of the
      brute-force sum, when numba is installed; else two block sizes), two regions (X and a smaller trapezoid inside it,
      the same events): the largest relative difference;
and, as development rows for H5-A to H5-C (decide nothing): the drift of the rest kink (the open point of probe P-b),
its width and energy (with the same readout applied to the exact field at the same events, to separate the readout from
the dynamics), the free motion at v = 0.3, 0.6, and the kick f = 0.05 for 1 <= t < 3 paired with the rest kink on the
same sprinkling (M_in = 2 pi Q f tau / Delta(gamma v)).

Readouts (coordinates, declared): slabs of width 1 (0.25 for I2); the kink fit
phi = 2 pi n0 + 4 arctan exp((x - x_c - v_s (t - t_m)) / w) on the events with |x - guess| < 4 in the slab (centre,
local velocity, width); Q from the plateaus 5 to 7 units from the centre; the energy
E = integral of [phi_t^2 / 2 + phi_x^2 / 2 + mu (1 - cos phi)] dx over |x - x_c| < 5 at the middle of the slab, with
phi, phi_t, phi_x from local linear regressions on the events of cells 0.5 x 0.5 (step 0.25).

Run from the root of the laboratory (numba strongly recommended: pip install numba):
    python -m diag.h5_instrument --quick                      # smoke test, a minute
    python -m diag.h5_instrument --seeds 16 --workers 8       # development rows
Output: a summary on the screen and the raw records in diag/out/h5_instrument_{stamp}.json.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import time

import numpy as np
from scipy.optimize import least_squares

from core.order.sprinkle import derive_seed
from core.stats.parallel import parallel_map

try:
    import numba as _nb
except ImportError:                     # the block engine is used instead (slower, the same sums)
    _nb = None

MU = 1.0
RHOS = (100.0, 200.0, 400.0)
X_KINK, T_KINK = 26.0, 10.0
X_WAVE, T_WAVE = 20.0, 10.0
VS = (0.3, 0.6)
KICK = dict(f=0.05, t1=1.0, t2=3.0)
I2_K = (0.5, 1.0, 2.0, 4.0)


# ============================================================================================== sprinkling

def sprinkle_trapezoid(rho, X, T, seed):
    """Poisson sprinkling of density rho in the domain of dependence {0 <= t <= T, |x| <= X - t}, sorted by t."""
    rng = np.random.default_rng(seed)
    area = 2.0 * X * T - T * T
    n = int(rng.poisson(rho * area))
    t = np.empty(0)
    x = np.empty(0)
    while t.size < n:
        m = int(1.3 * (n - t.size)) + 100
        tt = rng.uniform(0.0, T, m)
        xx = rng.uniform(-X, X, m)
        ok = np.abs(xx) <= X - tt
        t = np.concatenate([t, tt[ok]])
        x = np.concatenate([x, xx[ok]])
    t, x = t[:n], x[:n]
    o = np.argsort(t, kind="stable")
    return t[o], x[o]


# ============================================================================================== engines

if _nb is not None:
    @_nb.njit(cache=True)
    def _fenwick_evolve(U, V, tt, order, phi0, rho, mu, f, t1, t2):
        """Process the events in the given linear extension; the sum over the past is a 2D dominance sum in the
        null coordinates (u, v), kept in a Fenwick tree of Fenwick trees (offline, compressed)."""
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
        tree = np.zeros(off[N + 1] + N + 1, np.float64)      # inner Fenwick k occupies [off[k] + k, off[k+1] + k]
        phi = np.empty(N)
        for a in range(N):
            i = order[a]
            s = 0.0
            k = ru[i] - 1
            while k > 0:
                lo, hi = off[k], off[k + 1]
                p = np.searchsorted(lst[lo:hi], rv[i])            # elements with v-rank < rv[i]
                base = lo + k
                q = p
                while q > 0:
                    s += tree[base + q - 1]
                    q -= q & (-q)
                k -= k & (-k)
            ph = phi0[i] - 0.5 * s / rho
            phi[i] = ph
            w = mu * math.sin(ph) + (f if (tt[i] >= t1 and tt[i] < t2) else 0.0)
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


def _block_evolve(U, V, tt, order, phi0, rho, mu, f, t1, t2, block):
    """The same sums by brute force: the events in the order given, in blocks; the past outside the block by one
    vectorised product, the past inside the block event by event."""
    N = U.size
    phi = np.empty(N)
    w = np.zeros(N)
    done = np.zeros(0, np.int64)
    for b0 in range(0, N, block):
        idx = order[b0:b0 + block]
        if done.size:
            R = (U[done][None, :] < U[idx][:, None]) & (V[done][None, :] < V[idx][:, None])
            S = R.astype(np.float64) @ w[done]
        else:
            S = np.zeros(idx.size)
        Rb = (U[idx][None, :] < U[idx][:, None]) & (V[idx][None, :] < V[idx][:, None])
        for a, i in enumerate(idx):
            s = S[a] + (w[idx[:a]][Rb[a, :a]].sum() if a else 0.0)
            ph = phi0[i] - 0.5 * s / rho
            phi[i] = ph
            w[i] = mu * math.sin(ph) + (f if t1 <= tt[i] < t2 else 0.0)
        done = np.concatenate([done, idx])
    return phi


def evolve(t, x, phi0, rho, mu=MU, f=0.0, t1=0.0, t2=0.0, order="t", engine="auto", block=256):
    U = (t + x) / math.sqrt(2.0)
    V = (t - x) / math.sqrt(2.0)
    if order == "t":
        od = np.arange(t.size)
    elif order == "boost":                      # t + 0.5 x is a time function: a linear extension of the order
        od = np.argsort(t + 0.5 * x, kind="stable")
    else:
        raise ValueError(order)
    if engine == "auto":
        engine = "fenwick" if _nb is not None else "block"
    if engine == "fenwick":
        return _fenwick_evolve(U, V, t, od.astype(np.int64), phi0, float(rho), float(mu), float(f), float(t1), float(t2))
    return _block_evolve(U, V, t, od, phi0, float(rho), float(mu), float(f), float(t1), float(t2), int(block))


# ============================================================================================== initial data

def kink_exact(t, x, v, mu=MU):
    g = 1.0 / math.sqrt(1.0 - v * v)
    return 4.0 * np.arctan(np.exp(math.sqrt(mu) * g * (x - v * t)))


def kink_free(t, x, v, mu=MU):
    """d'Alembert: (g(x - t) + g(x + t)) / 2 + (1/2) int_{x-t}^{x+t} h, g = kink(0, .), h = phi_t(0, .) = -v g'."""
    gm, gp = kink_exact(0.0, x - t, v, mu), kink_exact(0.0, x + t, v, mu)
    return 0.5 * (gm + gp) - 0.5 * v * (gp - gm)


def modes_exact(t, x, modes, mu=MU):
    out = np.zeros_like(x)
    for a, k, th in modes:
        out += a * np.cos(k * x - math.sqrt(k * k + mu) * t + th)
    return out


def modes_free(t, x, modes, mu=MU):
    """Free wave of the initial data of a cos(kx - omega t + th): g = a cos(kx + th), h = a omega sin(kx + th)."""
    out = np.zeros_like(x)
    for a, k, th in modes:
        om = math.sqrt(k * k + mu)
        out += a * (np.cos(k * x + th) * np.cos(k * t) + (om / k) * np.sin(k * x + th) * np.sin(k * t))
    return out


# ============================================================================================== readouts

def _kink_model(p, tt, xx, tm, n0):
    xc, vs, lw = p
    return 2.0 * math.pi * n0 + 4.0 * np.arctan(np.exp(np.clip((xx - xc - vs * (tt - tm)) / math.exp(lw), -50, 50)))


def fit_kink(t, x, phi, a, b, guess, mu=MU):
    """Centre, local velocity and width of a kink in the slab a <= t < b; Q and the plateau levels."""
    tm = 0.5 * (a + b)
    s = (t >= a) & (t < b)
    near = s & (np.abs(x - guess) < 4.0)
    if near.sum() < 20:
        return None
    lo = s & (x < guess - 5.0) & (x > guess - 7.0)
    hi = s & (x > guess + 5.0) & (x < guess + 7.0)
    if lo.sum() < 5 or hi.sum() < 5:
        return None
    plo, phi_hi = float(np.mean(phi[lo])), float(np.mean(phi[hi]))
    n0 = round(plo / (2.0 * math.pi))
    r = least_squares(lambda p: _kink_model(p, t[near], x[near], tm, n0) - phi[near],
                      x0=[guess, 0.0, math.log(1.0 / math.sqrt(mu))])
    xc, vs, lw = r.x
    res = float(np.sqrt(np.mean(r.fun ** 2)))
    return {"xc": float(xc), "vs": float(vs), "w": float(math.exp(lw)), "res": res,
            "Q": (phi_hi - plo) / (2.0 * math.pi), "plo": plo, "phi": phi_hi}


def energy(t, x, phi, tm, xc, mu=MU, half=5.0, cell=0.5, step=0.25):
    """E = int [phi_t^2/2 + phi_x^2/2 + mu (1 - cos phi)] dx at time tm over |x - xc| < half (local regressions)."""
    s = np.abs(t - tm) < 0.5 * cell
    ts, xs, ps = t[s], x[s], phi[s]
    o = np.argsort(xs)
    ts, xs, ps = ts[o], xs[o], ps[o]
    E = 0.0
    used = 0
    for xg in np.arange(xc - half, xc + half + 1e-9, step):
        i0, i1 = np.searchsorted(xs, [xg - 0.5 * cell, xg + 0.5 * cell])
        if i1 - i0 < 6:
            continue
        A = np.column_stack([np.ones(i1 - i0), ts[i0:i1] - tm, xs[i0:i1] - xg])
        c, *_ = np.linalg.lstsq(A, ps[i0:i1], rcond=None)
        E += (0.5 * c[1] ** 2 + 0.5 * c[2] ** 2 + mu * (1.0 - math.cos(c[0]))) * step
        used += 1
    n_all = int(round(2 * half / step)) + 1
    return E * n_all / used if used else float("nan")


def project(t, x, phi, k, a, b):
    """Least squares phi = A cos kx + B sin kx (+ C) on the slab; returns the phase atan2(B, A) and the amplitude."""
    s = (t >= a) & (t < b)
    M = np.column_stack([np.cos(k * x[s]), np.sin(k * x[s]), np.ones(s.sum())])
    c, *_ = np.linalg.lstsq(M, phi[s], rcond=None)
    return math.atan2(c[1], c[0]), math.hypot(c[0], c[1])


# ============================================================================================== tasks

def kink_track(t, x, phi, v, T, mu=MU, exact=None):
    rows = []
    g = 1.0 / math.sqrt(1.0 - v * v)
    for a in np.arange(0.0, T - 1.0 + 1e-9, 1.0):
        guess = v * (a + 0.5)
        r = fit_kink(t, x, phi, a, a + 1.0, guess, mu)
        if r is None:
            rows.append(None)
            continue
        r["t"] = float(a + 0.5)
        r["E"] = energy(t, x, phi, a + 0.5, r["xc"], mu)
        if exact is not None:
            r["E_ex"] = energy(t, x, exact, a + 0.5, r["xc"], mu)
            rx = fit_kink(t, x, exact, a, a + 1.0, guess, mu)
            r["xc_ex"] = rx["xc"] if rx else float("nan")
            r["w_ex"] = rx["w"] if rx else float("nan")
            s = (t >= a) & (t < a + 1.0) & (np.abs(x - guess) < 6.0)
            r["rms"] = float(np.sqrt(np.mean((phi[s] - exact[s]) ** 2)))
        rows.append(r)
    return rows


def task_kink(rho, seed, arm, quick):
    X, T = (20.0, 8.0) if quick else (X_KINK, T_KINK)
    t, x = sprinkle_trapezoid(rho, X, T, seed)
    t0 = time.time()
    if arm == "rest" or arm == "kick":
        v = 0.0
    else:
        v = float(arm[1:])
    phi0 = kink_free(t, x, v)
    kw = dict(KICK) if arm == "kick" else {}
    phi = evolve(t, x, phi0, rho, **kw)
    ex = kink_exact(t, x, v) if arm != "kick" else None
    rows = kink_track(t, x, phi, v, T, exact=ex)
    return {"kind": "kink", "rho": rho, "seed": seed, "arm": arm, "v": v, "N": int(t.size),
            "sec": time.time() - t0, "rows": rows}


def task_wave(rho, seed, what, quick):
    X, T = (12.0, 5.0) if quick else (X_WAVE, T_WAVE)
    t, x = sprinkle_trapezoid(rho, X, T, seed)
    t0 = time.time()
    rng = np.random.default_rng(derive_seed(seed, "modes"))
    out = {"kind": "wave", "rho": rho, "seed": seed, "what": what, "N": int(t.size)}
    if what == "vacuum":
        ks = rng.uniform(0.25, 4.0, 8)
        modes = [(0.01 / math.sqrt(4.0), float(k), float(rng.uniform(0, 2 * math.pi))) for k in ks]
        ex = modes_exact(t, x, modes)
        phi = evolve(t, x, modes_free(t, x, modes), rho)
        rows = []
        for a in np.arange(0.0, T - 1.0 + 1e-9, 1.0):
            s = (t >= a) & (t < a + 1.0)
            rows.append({"t": float(a + 0.5), "rms_phi": float(np.sqrt(np.mean(phi[s] ** 2))),
                         "rms_ex": float(np.sqrt(np.mean(ex[s] ** 2))),
                         "rms_err": float(np.sqrt(np.mean((phi[s] - ex[s]) ** 2)))})
        out["rows"] = rows
    else:
        k = float(what[1:])
        modes = [(0.01, k, 0.0)]
        phi = evolve(t, x, modes_free(t, x, modes), rho)
        ph, amp, tm = [], [], []
        for a in np.arange(0.0, T - 0.25 + 1e-9, 0.25):
            p, A = project(t, x, phi, k, a, a + 0.25)
            ph.append(p)
            amp.append(A)
            tm.append(a + 0.125)
        # phi = a cos(kx - omega t) -> projection phase = +omega t
        phu = np.unwrap(np.array(ph))
        om = float(np.polyfit(np.array(tm), phu, 1)[0])
        out.update({"k": k, "omega": om, "omega_cont": math.sqrt(k * k + MU), "amp": amp, "t": tm})
    out["sec"] = time.time() - t0
    return out


def task_i4(rho, seed, quick):
    X, T = (8.0, 4.0) if quick else (12.0, 6.0)
    t, x = sprinkle_trapezoid(rho, X, T, seed)
    phi0 = kink_free(t, x, 0.3)
    base = evolve(t, x, phi0, rho, **KICK)
    res = {"kind": "i4", "rho": rho, "seed": seed, "N": int(t.size)}
    scale = float(np.max(np.abs(base)))
    alt = evolve(t, x, phi0, rho, order="boost", **KICK)
    res["order"] = float(np.max(np.abs(alt - base)) / scale)
    if _nb is not None:
        blk = evolve(t, x, phi0, rho, engine="block", block=256, **KICK)
        res["engine"] = float(np.max(np.abs(blk - base)) / scale)
        res["engines"] = "fenwick vs block"
    else:
        blk = evolve(t, x, phi0, rho, engine="block", block=64, **KICK)
        res["engine"] = float(np.max(np.abs(blk - base)) / scale)
        res["engines"] = "block 256 vs block 64"
    sub = np.abs(x) <= (X - 2.0) - t                 # the smaller trapezoid: the same events
    phs = evolve(t[sub], x[sub], phi0[sub], rho, **KICK)
    res["region"] = float(np.max(np.abs(phs - base[sub])) / scale)
    return res


def run_task(kind, *a):
    return {"kink": task_kink, "wave": task_wave, "i4": task_i4}[kind](*a)


# ============================================================================================== summary

def _ms(vals):
    v = np.array([x for x in vals if x is not None and np.isfinite(x)], float)
    if v.size == 0:
        return float("nan"), float("nan"), 0
    return float(v.mean()), float(v.std(ddof=1) / math.sqrt(v.size)) if v.size > 1 else float("nan"), int(v.size)


def _slope(rows, key="xc", tmin=0.0):
    pts = [(r["t"], r[key]) for r in rows if r is not None and r["t"] >= tmin]
    if len(pts) < 3:
        return float("nan")
    a = np.array(pts)
    return float(np.polyfit(a[:, 0], a[:, 1], 1)[0])


def summarise(res, rhos, quick):
    L = []
    L.append(f"engine: {'fenwick (numba)' if _nb is not None else 'block (numpy; install numba for speed)'}")
    i4 = [r for r in res if r["kind"] == "i4"]
    if i4:
        L.append("I4 order / engine / region (largest relative difference): " + ", ".join(
            f"{r['order']:.1e} / {r['engine']:.1e} / {r['region']:.1e}" for r in i4) + f"  [{i4[0]['engines']}]")
    for rho in rhos:
        L.append(f"--- rho {rho:g}")
        vac = [r for r in res if r["kind"] == "wave" and r["what"] == "vacuum" and r["rho"] == rho]
        if vac:
            last = [r["rows"][-1] for r in vac]
            first = [r["rows"][0] for r in vac]
            g = _ms([b["rms_phi"] / a["rms_phi"] for a, b in zip(first, last)])
            e = _ms([b["rms_err"] / b["rms_ex"] for b in last])
            L.append(f"  I1 vacuum: rms last / first {g[0]:.3f} +- {g[1]:.3f} (exact "
                     f"{np.mean([b['rms_ex'] / a['rms_ex'] for a, b in zip(first, last)]):.3f}); rel. error last slab "
                     f"{e[0]:.3f} +- {e[1]:.3f}")
        for k in I2_K:
            w = [r for r in res if r["kind"] == "wave" and r["what"] == f"k{k:g}" and r["rho"] == rho]
            if w:
                om = _ms([r["omega"] for r in w])
                oc = math.sqrt(k * k + MU)
                L.append(f"  I2 k {k:g}: omega {om[0]:.4f} +- {om[1]:.4f} (continuum {oc:.4f}, massless {k:.4f}; "
                         f"ratio {om[0] / oc:.4f}; mass shift reproduced {(om[0] - k) / (oc - k):.3f})")
        for arm in ["rest", "v0.3", "v0.6", "kick"]:
            ks = [r for r in res if r["kind"] == "kink" and r["arm"] == arm and r["rho"] == rho]
            if not ks:
                continue
            v = ks[0]["v"]
            g = 1.0 / math.sqrt(1.0 - v * v)
            ok = [all(rr is not None for rr in r["rows"]) for r in ks]
            Q = _ms([np.mean([rr["Q"] for rr in r["rows"] if rr]) for r in ks])
            sl = _ms([_slope(r["rows"]) for r in ks])
            if arm == "kick":
                rest = {r["seed"]: r for r in res if r["kind"] == "kink" and r["arm"] == "rest" and r["rho"] == rho}
                tmin = KICK["t2"] + 1.0
                dv = []
                for r in ks:
                    if r["seed"] in rest:
                        vk = _slope(r["rows"], tmin=tmin)
                        vr = _slope(rest[r["seed"]]["rows"], tmin=tmin)
                        dv.append(vk / math.sqrt(1 - vk * vk) - vr / math.sqrt(1 - vr * vr))
                d = _ms(dv)
                F = 2 * math.pi * KICK["f"] * (KICK["t2"] - KICK["t1"])
                L.append(f"  kick: Q {Q[0]:.3f}; Delta(gamma v) {d[0]:.4f} +- {d[1]:.4f} (continuum, M = 8: {F / 8:.4f}); "
                         f"M_in = 2 pi Q f tau / Delta(gamma v) = {F * Q[0] / d[0]:.2f}")
                continue
            dx = _ms([(r["rows"][-1]["xc"] - r["rows"][0]["xc"]) if (r["rows"][-1] and r["rows"][0]) else None
                      for r in ks])
            dxe = _ms([(r["rows"][-1]["xc_ex"] - r["rows"][0]["xc_ex"]) if (r["rows"][-1] and r["rows"][0]) else None
                       for r in ks])
            w = _ms([np.mean([rr["w"] for rr in r["rows"] if rr]) for r in ks])
            E = _ms([np.mean([rr["E"] for rr in r["rows"] if rr]) for r in ks])
            Ee = _ms([np.mean([rr["E_ex"] for rr in r["rows"] if rr]) for r in ks])
            rms = [np.mean([rr["rms"] for rr in r["rows"] if rr]) for r in ks]
            rms_last = [r["rows"][-1]["rms"] for r in ks if r["rows"][-1]]
            res_ = _ms([np.mean([rr["res"] for rr in r["rows"] if rr]) for r in ks])
            L.append(f"  {arm:5s}: all slabs read {sum(ok)}/{len(ok)}; Q {Q[0]:.3f}; v {sl[0]:+.4f} +- {sl[1]:.4f} "
                     f"(launch {v:g}); centre shift over the run {dx[0]:+.3f} +- {dx[1]:.3f} (readout on exact "
                     f"{dxe[0]:+.3f}); w*gamma {w[0] * g:.3f} +- {w[1] * g:.3f} (1); E/gamma {E[0] / g:.3f} +- "
                     f"{E[1] / g:.3f} (readout on exact {Ee[0] / g:.3f}; continuum 8); rms err mean {np.mean(rms):.3f}, "
                     f"last slab {np.mean(rms_last):.3f}; fit residual {res_[0]:.3f}")
    return L


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=16)
    ap.add_argument("--workers", type=int, default=0)
    ap.add_argument("--rhos", type=str, default=",".join(f"{r:g}" for r in RHOS))
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--base", type=int, default=20261110)
    a = ap.parse_args()
    from core.stats.parallel import default_workers
    W = default_workers(a.workers)
    rhos = [float(r) for r in a.rhos.split(",")]
    nseed = a.seeds
    if a.quick:
        rhos, nseed = [100.0], 2
    jobs = []
    for rho in rhos:
        for k in range(nseed):
            sd = derive_seed(a.base, "H5-0", rho, k)
            for arm in ["rest", "v0.3", "v0.6", "kick"]:
                jobs.append(("kink", rho, sd, arm, a.quick))
            jobs.append(("wave", rho, sd, "vacuum", a.quick))
            for kk in I2_K:
                jobs.append(("wave", rho, sd, f"k{kk:g}", a.quick))
    for k in range(2 if a.quick else 4):
        jobs.append(("i4", 100.0, derive_seed(a.base, "H5-0-i4", k), a.quick))
    jobs.sort(key=lambda j: -j[1] if j[0] != "i4" else 0)       # the largest first
    t0 = time.time()
    res = parallel_map(run_task, jobs, W, print,
                       lambda i, r, d, n: f"[H5-0] {d}/{n}: {r['kind']} rho {r['rho']:g} N {r['N']} "
                                          f"({time.time() - t0:.0f}s)")
    lines = summarise(res, rhos, a.quick)
    print("\n".join(lines))
    os.makedirs(os.path.join("diag", "out"), exist_ok=True)
    stamp = time.strftime("%Y%m%dT%H%M%S", time.gmtime())
    fn = os.path.join("diag", "out", f"h5_instrument_{stamp}.json")
    with open(fn, "w") as fh:
        json.dump({"args": vars(a), "numba": _nb is not None, "summary": lines, "results": res}, fh)
    print(f"raw records: {fn}  ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
