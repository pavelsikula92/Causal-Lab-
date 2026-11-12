"""
phases.F5_gravity.field1d_curved -- the field of F3.14a on a 1+1 order of nonuniform density (helpers of F5.2a; no TEST)
=====================================================================================================================

Ported unchanged from diag/f52a_fall.py (development rows of F5.2a, 10 Nov 2026), which keeps its own copy as the
record of those rows. Every 1+1 metric is conformally flat, ds^2 = Omega^2 (dt^2 - dx^2): the causal order is that of
flat spacetime and only the number of events per coordinate area changes (rho Omega^2). The rule of F3.14a
(phases/F3_structure/field1d.py) is used unchanged; here: the sprinkling of density rho exp(2 Phi(x)), Phi = g L
tanh(x / L); the launch profile in proper distance; the weighted engine of the counting null (every term divided by
Omega^2(y), coordinates, a null model); the kink tracking; the geodesic; the continuum field equation
phi_tt - phi_xx + Omega^2 mu sin phi = 0 by finite differences with the same initial data and centre readout.
"""

from __future__ import annotations

import math

import numpy as np

from phases.F3_structure import field1d as F

_nb = F._nb
G, L, X0 = 0.05, 10.0, 0.0
XB, TB = 26.0, 10.0


def Phi(x, g=G):
    return g * L * np.tanh(x / L)


def dPhi(x, g=G):
    return g / np.cosh(x / L) ** 2


def sprinkle_curved(rho, X, T, g, seed):
    """FLAWED (found 10 Nov after the run of F5.2a; kept unchanged as the record of that run): the number of points is
    drawn for the area of the trapezoid, but the points are placed in the whole box [0, T] x [-X, X] and those outside
    the trapezoid are dropped, so the density is rho (2XT - T^2) / (2XT) exp(2 Phi), i.e. 0.808 rho exp(2 Phi) for
    X 26, T 10, while the rule normalises with rho: an effective coupling 0.808 mu. Use sprinkle_curved2."""
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


def sprinkle_curved2(rho, X, T, g, seed):
    """Poisson process of coordinate density rho exp(2 Phi(x)) in the trapezoid |x| <= X - t, sorted by t (F5.2a2):
    a homogeneous process of density rho max(exp(2 Phi)) placed uniformly IN the trapezoid (as field1d.sprinkle_trapezoid),
    then thinned with probability exp(2 Phi(x)) / max."""
    rng = np.random.default_rng(seed)
    pmax = math.exp(2.0 * g * L) if g > 0 else 1.0
    n = int(rng.poisson(rho * pmax * (2.0 * X * T - T * T)))
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
    keep = rng.uniform(0.0, 1.0, n) < np.exp(2.0 * Phi(x, g)) / pmax
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


