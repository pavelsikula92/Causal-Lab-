"""
OV2 sec. 14, T3 pilot: numba kernels for single sprinklings (development; imported, never run as __main__, so the
numba cache works under Windows spawn).

A 2D sprinkling in (t, z), sorted by t; x < y iff t_y - t_x > |z_y - z_x|. The pairwise matter action of T1:
    S_m = sum_{x<y} f(n_xy, eps) L(phi_x, phi_y) + sum_x U(phi_x),   L = -eps^2 (phi_x - phi_y)^2,   U = -V(phi)/rho,
f(n, eps) = (1-eps)^n [1 - 2 eps n/(1-eps) + eps^2 n (n-1)/(2 (1-eps)^2)] (F6.4a), V = m^2 (1 - cos phi).
For a window chi(t) and a direction Phi_t (values on the events) the kernel returns, per sprinkling,
    B   = sum_x chi_x E_x Phi_t(x),             E_x = dS_m/dphi_x (all partners)
    P   = sum_x sum_y chi_x (1 - chi_y) w dL/dphi_x Phi_t(x)       (the pair term of R)
    Q2, Q1: the parts of R = P - sum_pairs (chi'_x chi_y + chi_x chi'_y) w L - sum chi'_x U_x at the ramps t2 and t1
    E (array), the residual per event.
"""
import math

import numba as nb
import numpy as np


@nb.njit(cache=True, inline="always")
def fker(n, eps):
    a = 1.0 - eps
    return a ** n * (1.0 - 2.0 * eps * n / a + eps * eps * n * (n - 1) / (2.0 * a * a))


@nb.njit(cache=True, inline="always")
def chi_of(t, t1, t2, d):
    if t <= t1 - d or t >= t2 + d:
        return 0.0, 0.0
    if t < t1 + d:
        x = (t - (t1 - d)) / (2.0 * d)
        return math.sin(0.5 * math.pi * x) ** 2, math.pi / (4.0 * d) * math.sin(math.pi * x)
    if t > t2 - d:
        x = (t - (t2 - d)) / (2.0 * d)
        return math.cos(0.5 * math.pi * x) ** 2, -math.pi / (4.0 * d) * math.sin(math.pi * x)
    return 1.0, 0.0


@nb.njit(cache=True)
def past_bits(t, z):
    """P[y, word] bitset of the past of y (indices sorted by t)."""
    n = t.size
    nw = (n + 63) // 64
    P = np.zeros((n, nw), dtype=np.uint64)
    for y in range(n):
        for x in range(y):
            if t[y] - t[x] > abs(z[y] - z[x]):
                P[y, x >> 6] |= np.uint64(1) << np.uint64(x & 63)
    return P


@nb.njit(cache=True, inline="always")
def _popcount(x):
    x = x - ((x >> np.uint64(1)) & np.uint64(0x5555555555555555))
    x = (x & np.uint64(0x3333333333333333)) + ((x >> np.uint64(2)) & np.uint64(0x3333333333333333))
    x = (x + (x >> np.uint64(4))) & np.uint64(0x0F0F0F0F0F0F0F0F)
    return int((x * np.uint64(0x0101010101010101)) >> np.uint64(56))


@nb.njit(cache=True)
def interval_count(P, x, y):
    """n_xy = |{w : x < w < y}| = popcount(future(x) & past(y)); future(x) & past(y) = {w in past(y) : x in past(w)}."""
    c = 0
    for w in range(x + 1, y):
        if (P[y, w >> 6] >> np.uint64(w & 63)) & np.uint64(1):
            if (P[w, x >> 6] >> np.uint64(x & 63)) & np.uint64(1):
                c += 1
    return c


@nb.njit(cache=True)
def interval_counts_fast(P, F, x, y):
    """the same with the future bitset F of x: popcount(F[x] & P[y]) over the words between x and y."""
    c = 0
    for k in range((x + 1) >> 6, (y >> 6) + 1):
        c += _popcount(F[x, k] & P[y, k])
    return c


@nb.njit(cache=True)
def future_bits(P):
    n = P.shape[0]
    F = np.zeros_like(P)
    for y in range(n):
        for k in range(P.shape[1]):
            wd = P[y, k]
            while wd:
                b = wd & (~wd + np.uint64(1))
                x = (k << 6) + _popcount(b - np.uint64(1))
                F[x, y >> 6] |= np.uint64(1) << np.uint64(y & 63)
                wd ^= b
    return F


@nb.njit(cache=True)
def noether(t, z, phi, phit, rho, eps, m, t1, t2, d, P, F):
    n = t.size
    mid = 0.5 * (t1 + t2)
    c = eps * eps
    E = np.zeros(n)
    chi = np.zeros(n)
    dchi = np.zeros(n)
    for i in range(n):
        chi[i], dchi[i] = chi_of(t[i], t1, t2, d)
    Pp = 0.0
    Q2 = 0.0
    Q1 = 0.0
    S = 0.0
    for b in range(n):
        for k in range(P.shape[1]):
            wd = P[b, k]
            while wd:
                bit = wd & (~wd + np.uint64(1))
                a = (k << 6) + _popcount(bit - np.uint64(1))
                wd ^= bit
                nab = interval_counts_fast(P, F, a, b)
                w = fker(nab, eps)
                D = phi[a] - phi[b]
                L = -c * D * D
                dLa = -2.0 * c * D
                dLb = 2.0 * c * D
                S += w * L
                E[a] += w * dLa
                E[b] += w * dLb
                t1v = chi[a] * (1.0 - chi[b]) * w * dLa * phit[a]
                t2v = chi[b] * (1.0 - chi[a]) * w * dLb * phit[b]
                Pp += t1v + t2v
                t3a = -dchi[a] * chi[b] * w * L
                t3b = -chi[a] * dchi[b] * w * L
                if t[b] > mid:
                    Q2 += t1v + t3b
                else:
                    Q1 += t1v + t3b
                if t[a] > mid:
                    Q2 += t2v + t3a
                else:
                    Q1 += t2v + t3a
    B = 0.0
    for x in range(n):
        V = m * m * (1.0 - math.cos(phi[x]))
        Up = -m * m * math.sin(phi[x]) / rho
        U = -V / rho
        S += U
        E[x] += Up
        B += chi[x] * E[x] * phit[x]
        if t[x] > mid:
            Q2 += -dchi[x] * U
        else:
            Q1 += -dchi[x] * U
    return B, Pp, Q2, Q1, S, E


# ================================================================================ fast path (2D only)
# In 1+1, x < y iff du > 0 and dv > 0 (u = t - z, v = t + z), and n_xy is the number of events in the open
# rectangle (u_x, u_y) x (v_x, v_y). For each a the future events are visited in increasing u; the v-ranks of the
# future events already visited are kept in a Fenwick tree, so n_ab = #(visited with v < v_b). Cost O(#pairs log N),
# no bitsets; all arms (rows of PHI, PHIT) are accumulated in the same pass over the pairs.

@nb.njit(cache=True)
def noether_multi(t, z, PHI, PHIT, rho, eps, m, t1, t2, d):
    """As noether() for several arms at once: PHI, PHIT have shape (n_arms, N), events in any order.
    Returns B, P, Q2, Q1, S (arrays of length n_arms)."""
    n = t.size
    na = PHI.shape[0]
    mid = 0.5 * (t1 + t2)
    c = eps * eps
    u = t - z
    v = t + z
    ou = np.argsort(u)
    vr = np.empty(n, dtype=np.int64)                 # 1-based ranks of v
    ov = np.argsort(v)
    for r in range(n):
        vr[ov[r]] = r + 1
    chi = np.zeros(n)
    dchi = np.zeros(n)
    for i in range(n):
        chi[i], dchi[i] = chi_of(t[i], t1, t2, d)
    E = np.zeros((na, n))
    B = np.zeros(na)
    Pp = np.zeros(na)
    Q2 = np.zeros(na)
    Q1 = np.zeros(na)
    S = np.zeros(na)
    fen = np.zeros(n + 1, dtype=np.int64)
    ins = np.empty(n, dtype=np.int64)
    for ia in range(n):
        a = ou[ia]
        va = v[a]
        nins = 0
        for ib in range(ia + 1, n):
            b = ou[ib]
            if v[b] <= va:
                continue
            # n_ab = number of inserted (future of a, u < u_b) with v < v_b
            cnt = 0
            r = vr[b] - 1
            while r > 0:
                cnt += fen[r]
                r -= r & (-r)
            w = fker(cnt, eps)
            ta = t[a]
            tb = t[b]
            for k in range(na):
                D = PHI[k, a] - PHI[k, b]
                L = -c * D * D
                dLa = -2.0 * c * D
                dLb = 2.0 * c * D
                S[k] += w * L
                E[k, a] += w * dLa
                E[k, b] += w * dLb
                t1v = chi[a] * (1.0 - chi[b]) * w * dLa * PHIT[k, a]
                t2v = chi[b] * (1.0 - chi[a]) * w * dLb * PHIT[k, b]
                Pp[k] += t1v + t2v
                t3a = -dchi[a] * chi[b] * w * L
                t3b = -chi[a] * dchi[b] * w * L
                if tb > mid:
                    Q2[k] += t1v + t3b
                else:
                    Q1[k] += t1v + t3b
                if ta > mid:
                    Q2[k] += t2v + t3a
                else:
                    Q1[k] += t2v + t3a
            # insert b
            r = vr[b]
            while r <= n:
                fen[r] += 1
                r += r & (-r)
            ins[nins] = vr[b]
            nins += 1
        # undo the insertions of this a
        for j in range(nins):
            r = ins[j]
            while r <= n:
                fen[r] -= 1
                r += r & (-r)
    for k in range(na):
        for x in range(n):
            ph = PHI[k, x]
            U = -m * m * (1.0 - math.cos(ph)) / rho
            S[k] += U
            E[k, x] += -m * m * math.sin(ph) / rho
            B[k] += chi[x] * E[k, x] * PHIT[k, x]
            if t[x] > mid:
                Q2[k] += -dchi[x] * U
            else:
                Q1[k] += -dchi[x] * U
    return B, Pp, Q2, Q1, S
