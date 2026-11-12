"""
phases.F5_gravity.growth1d -- the causal growth rule of the order with feedback from the field (helpers of F5.2b; no TEST)
=======================================================================================================================

Ported from diag/growth_probe.py (development rows of direction A, 11 Nov 2026), with initial velocities added for
the boosted launch (A2). A candidate event x of a Poisson sprinkling of (proper) density rho, carrying an independent
mark u_x fixed at the sprinkling, is realised iff u_x < exp(2 L(x)),
    L(x) = -(kappa / 2 rho) * sum over realised y < x of [2 V_A(phi_A(y)) + 2 V_B(phi_B(y))],
the normalisation of the retarded Green function of 1+1 (1/2 inside the cone). Two separate fields phi_A, phi_B evolve
on the realised events with the unchanged rule of F3.14a. Continuum law: (d_t^2 - d_x^2) L = -kappa Pi^2 q, i.e.
box_g L = -kappa q in ds^2 = Pi^2 (dt^2 - dx^2); (d_t^2 - d_x^2) phi_i = -Pi^2 V_i'(phi_i).
"""

from __future__ import annotations

import math

import numpy as np

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


def evolve_growth(t, x, marks, muA, muB, xA, xB, rho, kappa, order="t", v=0.0):
    U, V = (t + x) / math.sqrt(2.0), (t - x) / math.sqrt(2.0)
    od = np.arange(t.size) if order == "t" else np.argsort(t + 0.5 * x, kind="stable")
    phiA0 = F.kink_free(t, x - xA, v, muA)
    phiB0 = F.kink_free(t, x - xB, v, muB)
    return _grow(U, V, od.astype(np.int64), marks, phiA0, phiB0, float(rho), float(muA), float(muB), float(kappa))


def track(t, x, phi, real, mu, x0, T):
    tt, xx, pp = t[real], x[real], phi[real]
    xc = []
    for a in np.arange(0.0, T - 1.0 + 1e-9, 1.0):
        g = xc[-1] if (xc and np.isfinite(xc[-1])) else x0
        r = F.fit_kink(tt, xx, pp, a, a + 1.0, g, mu)
        xc.append(r["xc"] if r else float("nan"))
    return xc


def continuum(muA, muB, xA, xB, kappa, T, dx=0.01, X=60.0, v=0.0):
    xs = np.arange(-X, X + dx / 2, dx)
    dt = 0.5 * dx
    lap = lambda p: np.concatenate([[0.0], (p[2:] - 2 * p[1:-1] + p[:-2]) / dx ** 2, [0.0]])
    a0 = F.kink_exact(0.0, xs - xA, v, muA)
    b0 = F.kink_exact(0.0, xs - xB, v, muB)
    g = 1.0 / math.sqrt(1.0 - v * v)
    at = -v * g * 2.0 * math.sqrt(muA) / np.cosh(math.sqrt(muA) * g * (xs - xA))
    bt = -v * g * 2.0 * math.sqrt(muB) / np.cosh(math.sqrt(muB) * g * (xs - xB))
    l0 = np.zeros_like(xs)

    def rhs(a, b, l):
        om2 = np.exp(2.0 * l)
        q = 2.0 * muA * (1 - np.cos(a)) + 2.0 * muB * (1 - np.cos(b))
        return lap(a) - om2 * muA * np.sin(a), lap(b) - om2 * muB * np.sin(b), lap(l) - kappa * om2 * q

    ra, rb, rl = rhs(a0, b0, l0)
    pa, pb, pl = a0, b0, l0
    ca, cb, cl = a0 + dt * at + 0.5 * dt * dt * ra, b0 + dt * bt + 0.5 * dt * dt * rb, l0 + 0.5 * dt * dt * rl
    cenA, cenB, cnt = np.zeros(int(T)), np.zeros(int(T)), np.zeros(int(T))
    nsteps = int(round(T / dt))
    for n in range(1, nsteps):
        tn = n * dt
        s = int(tn)
        if s < int(T) and n % 4 == 0:
            for arr, cen, x0 in ((ca, cenA, xA + v * tn), (cb, cenB, xB + v * tn)):
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


