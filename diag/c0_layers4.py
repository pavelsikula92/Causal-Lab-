"""
C0, development: layered link kernels for the 3+1 Green function (decides nothing)
=================================================================================
The link kernel c L0 (Johnston) has the right zeroth moment across the light cone but a first moment (its cone is
smeared inward by sigma = tau^2 ~ rho^-1/2), which shifts the phase of a propagating source and does not vanish with the
depth of the past (development 13 Nov: the vacuum of the linear rule drifts by 0.4 A within t 3; a longer prepared past
made it worse). Layers L_k = pairs u < v with exactly k events in the interval (intrinsic: counts only) have the
profile (rho V)^k e^{-rho V}/k!; per unit of the cone their zeroth moment is ~ Gamma(k + 1/2)/k!, their first moment in
sigma is the same for all k, the second ~ Gamma(k + 3/2)/k!. Hence (in units of c = 1/(2 pi sqrt(6) sqrt(rho))):
    K1 = L0                    (zeroth moment)
    K2 = 2 (L0 - L1)           (zeroth, first moment zero)
    K3 = 3 L0 - 7 L1 + 4 L2    (zeroth, first and second moments zero)
"""
from __future__ import annotations

import math

import numpy as np

from phases.F3_structure import field1d as F

_nb = F._nb

if _nb is not None:
    @_nb.njit(cache=True)
    def layers4(t, x, y, z, tc, kmax):
        """For every event: its past pairs with exactly k = 0..kmax events in the interval (CSR per layer, packed:
        idx and lay arrays)."""
        N = t.size
        ptr = np.zeros(N + 1, np.int64)
        cap = 512 * N
        idx = np.empty(cap, np.int64)
        lay = np.empty(cap, np.int8)
        tot = 0
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
            # cand is in decreasing time; elements of the interval of (u, j) are candidates later than u
            for a in range(m):
                u = cand[a]
                cnt = 0
                for b in range(a):
                    w = cand[b]
                    dt = t[w] - t[u]
                    if dt > 0.0 and dt * dt > (x[w] - x[u]) ** 2 + (y[w] - y[u]) ** 2 + (z[w] - z[u]) ** 2:
                        cnt += 1
                        if cnt > kmax:
                            break
                if cnt <= kmax:
                    if tot >= cap:
                        ni = np.empty(2 * cap, np.int64)
                        ni[:tot] = idx[:tot]
                        nl = np.empty(2 * cap, np.int8)
                        nl[:tot] = lay[:tot]
                        idx, lay, cap = ni, nl, 2 * cap
                    idx[tot] = u
                    lay[tot] = cnt
                    tot += 1
            ptr[j + 1] = tot
        return ptr, idx[:tot], lay[:tot]

    @_nb.njit(cache=True)
    def evolve_layers(ptr, idx, lay, coef, phi0, pre, phi_pre, c, mu, linear):
        N = phi0.size
        phi = np.empty(N)
        for v in range(N):
            if pre[v]:
                phi[v] = phi_pre[v]
                continue
            s = 0.0
            for e in range(ptr[v], ptr[v + 1]):
                a = coef[lay[e]]
                if a != 0.0:
                    u = idx[e]
                    s += a * (mu * phi[u] if linear else mu * math.sin(phi[u]))
            phi[v] = phi0[v] - c * s
        return phi


if _nb is not None:
    @_nb.njit(cache=True)
    def layers_evolve_multi(t, x, y, z, tc, coef, phi0, lin, c, mu):
        """Fused: for every event in time order, find its past pairs with n <= kmax events in the interval
        (kmax = coef.size - 1), and update all M fields at once (phi0: M x N; lin[m]: linear J = mu phi, else
        mu sin phi). No pairs are stored (memory O(M N)). Returns (phi, number of pairs)."""
        M, N = phi0.shape
        kmax = coef.size - 1
        phi = np.empty((M, N))
        tc2 = tc * tc
        cand = np.empty(N, np.int64)
        s = np.zeros(M)
        npairs = 0
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
            for i in range(M):
                s[i] = 0.0
            for a in range(m):
                u = cand[a]
                cnt = 0
                for b in range(a):
                    w = cand[b]
                    dt = t[w] - t[u]
                    if dt > 0.0 and dt * dt > (x[w] - x[u]) ** 2 + (y[w] - y[u]) ** 2 + (z[w] - z[u]) ** 2:
                        cnt += 1
                        if cnt > kmax:
                            break
                if cnt <= kmax:
                    npairs += 1
                    wgt = coef[cnt]
                    if wgt != 0.0:
                        for i in range(M):
                            s[i] += wgt * (mu * phi[i, u] if lin[i] else mu * math.sin(phi[i, u]))
            for i in range(M):
                phi[i, j] = phi0[i, j] - c * s[i]
        return phi, npairs
