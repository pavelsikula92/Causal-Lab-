"""
phases.F3_structure.field1d -- a field of the causal order in 1+1 (helpers of F3.14a; no TEST in this file)
==========================================================================================================

Ported unchanged from diag/h5_instrument.py (the development rows of H5-0, 10 Nov 2026), which keeps its own copy as
the record of those rows.

Dynamics (card docs/H5_card_draft.md, decision D2'): the retarded Green function of d_t^2 - d_x^2 in 1+1 is 1/2 inside
the past light cone, on a causal set 1/2 times the causal relation. For phi_tt - phi_xx + mu sin(phi) = -f chi(t),

    phi(x) = phi_0(x) - (1 / 2 rho) * sum over events y < x (the order), t_y >= 0, of [mu sin(phi(y)) + f chi(t_y)],

phi_0 the free wave of the initial data on t = 0 (d'Alembert's formula; the launch). The sum over the causal past is a
dominance sum in the null coordinates, kept in a Fenwick tree of Fenwick trees (numba); a brute-force block engine
computes the same sums (I4). The telegraph null (phi_tt + Gamma phi_t - phi_xx + V' = -f chi) is evolved for
psi = exp(Gamma t / 2) phi with the same sum: psi(x) = psi_0(x) - (1/2 rho) sum [-(Gamma^2/4) psi(y) + exp(Gamma t_y / 2)
(mu sin phi(y) + f chi)], which uses the coordinate time of a medium frame (declared; a null model). The causality-blind
control sums over m_x events drawn uniformly with replacement from all earlier events, m_x the size of the causal past.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.optimize import least_squares

try:
    import numba as _nb
except ImportError:
    _nb = None

MU = 1.0


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
    def _fenwick_evolve(U, V, tt, order, phi0, rho, mu, f, t1, t2, gam, unit):
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
            if unit:                                            # count of the causal past
                phi[i] = s
                w = 1.0
            elif gam > 0.0:                                     # telegraph: phi = exp(-gam t / 2) psi
                psi = phi0[i] - 0.5 * s / rho
                eg = math.exp(0.5 * gam * tt[i])
                ph = psi / eg
                phi[i] = ph
                w = -0.25 * gam * gam * psi + eg * (mu * math.sin(ph) + (f if (tt[i] >= t1 and tt[i] < t2) else 0.0))
            else:
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


if _nb is not None:
    @_nb.njit(cache=True)
    def _blind_evolve(tt, phi0, m, rho, mu, f, t1, t2, seed):
        """Control: each event sums over m[i] events drawn uniformly (with replacement) from all earlier events (the
        same time range 0 <= t_y < t_x as its causal past, the same number of terms, the same normalisation)."""
        np.random.seed(seed)
        N = tt.size
        phi = np.empty(N)
        w = np.zeros(N)
        for i in range(N):
            s = 0.0
            if i > 0:
                for k in range(int(m[i])):
                    s += w[np.random.randint(0, i)]
            ph = phi0[i] - 0.5 * s / rho
            phi[i] = ph
            w[i] = mu * math.sin(ph) + (f if (tt[i] >= t1 and tt[i] < t2) else 0.0)
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


def evolve(t, x, phi0, rho, mu=MU, f=0.0, t1=0.0, t2=0.0, order="t", engine="auto", block=256, gam=0.0,
           unit=False):
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
        return _fenwick_evolve(U, V, t, od.astype(np.int64), phi0, float(rho), float(mu), float(f), float(t1), float(t2),
                               float(gam), bool(unit))
    if gam > 0.0 or unit:
        raise RuntimeError("the telegraph arm and the past counts need numba")
    return _block_evolve(U, V, t, od, phi0, float(rho), float(mu), float(f), float(t1), float(t2), int(block))


# ============================================================================================== initial data

def kink_exact(t, x, v, mu=MU):
    g = 1.0 / math.sqrt(1.0 - v * v)
    return 4.0 * np.arctan(np.exp(math.sqrt(mu) * g * (x - v * t)))


def kink_free(t, x, v, mu=MU):
    """d'Alembert: (g(x - t) + g(x + t)) / 2 + (1/2) int_{x-t}^{x+t} h, g = kink(0, .), h = phi_t(0, .) = -v g'."""
    gm, gp = kink_exact(0.0, x - t, v, mu), kink_exact(0.0, x + t, v, mu)
    return 0.5 * (gm + gp) - 0.5 * v * (gp - gm)


_GRID = None


def _kink_integral(s, v, mu=MU):
    """int_0^s g(z) dz for g = kink profile at t = 0 (cumulative trapezoid on a fine grid, error < 1e-7)."""
    global _GRID
    key = (v, mu)
    if _GRID is None or _GRID[0] != key:
        z = np.linspace(-80.0, 80.0, 320001)
        g = kink_exact(0.0, z, v, mu)
        c = np.concatenate([[0.0], np.cumsum(0.5 * (g[1:] + g[:-1]) * np.diff(z))])
        c -= np.interp(0.0, z, c)
        _GRID = (key, z, c)
    _, z, c = _GRID
    return np.interp(s, z, c)


def kink_free_telegraph(t, x, v, gam, mu=MU):
    """Free wave of psi = exp(gam t / 2) phi: psi(0) = g, psi_t(0) = h + gam g / 2, h = -v g'."""
    base = kink_free(t, x, v, mu)
    return base + 0.25 * gam * (_kink_integral(x + t, v, mu) - _kink_integral(x - t, v, mu))


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
                      x0=[guess, 0.0, math.log(1.0 / math.sqrt(max(mu, 1.0)))])
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


