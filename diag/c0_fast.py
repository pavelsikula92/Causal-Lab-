"""
C0: the same rule as diag/c0_card.evolve_card, laid out for the memory (14 Nov; a computational change only)
===========================================================================================================
Identical arithmetic in identical order; two changes of layout:
  * the candidates' coordinates are gathered into contiguous arrays before the interval counts;
  * the source J and the field are stored event-major, (N, R), so the R rows of one past event are adjacent.
F3.16 showed that one sprinkling at rho 30 took 320 s at X 7 and 4 040 s at X 8.5 (2.2x the events): the run was
bound by cache misses (J[r, u] for 16 rows at random u), not by the arithmetic. Checked bit-identical against
evolve_card (check()).
"""
from __future__ import annotations

import math
import time

import numpy as np

from phases.F3_structure import field1d as F

_nb = F._nb

if _nb is not None:
    @_nb.njit(cache=True)
    def evolve_fast(t, x, y, z, tc, kmax, W, ctd, phi0, lin, c, mu, sample, rbin, nrb):
        N = t.size
        K = W.shape[0]
        R = phi0.shape[0]
        Mo = R // K
        P0 = np.ascontiguousarray(phi0.T)
        phi = np.empty((N, R))
        J = np.empty((N, R))
        s = np.zeros(R)
        cnt = np.zeros((N, kmax + 1), np.int32)
        hist = np.zeros((nrb, kmax + 1))
        tc2 = tc * tc
        cand = np.empty(N, np.int64)
        ct_ = np.empty(N)
        cx = np.empty(N)
        cy = np.empty(N)
        cz = np.empty(N)
        for j in range(N):
            m = 0
            tj, xj, yj, zj = t[j], x[j], y[j], z[j]
            for k in range(j - 1, -1, -1):
                dt = tj - t[k]
                if dt <= 0.0:
                    continue
                q = dt * dt - (xj - x[k]) ** 2 - (yj - y[k]) ** 2 - (zj - z[k]) ** 2
                if q > 0.0 and q <= tc2:
                    cand[m] = k
                    ct_[m] = t[k]
                    cx[m] = x[k]
                    cy[m] = y[k]
                    cz[m] = z[k]
                    m += 1
            for i in range(R):
                s[i] = 0.0
            for a in range(m):
                u = cand[a]
                tu, xu, yu, zu = ct_[a], cx[a], cy[a], cz[a]
                nn = 0
                for b in range(a):
                    dt = ct_[b] - tu
                    if dt > 0.0 and dt * dt > (cx[b] - xu) ** 2 + (cy[b] - yu) ** 2 + (cz[b] - zu) ** 2:
                        nn += 1
                        if nn > kmax:
                            break
                if nn <= kmax:
                    cnt[j, nn] += 1
                    if sample[j]:
                        r = math.sqrt((xj - xu) ** 2 + (yj - yu) ** 2 + (zj - zu) ** 2)
                        ib = int(r / rbin)
                        if ib < nrb:
                            hist[ib, nn] += 1.0
                    for kk in range(K):
                        wg = W[kk, nn]
                        if wg != 0.0:
                            base = kk * Mo
                            for mo in range(Mo):
                                s[base + mo] += wg * J[u, base + mo]
            for i in range(R):
                rhs = P0[j, i] - c * s[i]
                d = ctd[i // Mo] * mu
                if lin[i % Mo]:
                    p = rhs / (1.0 - d)
                    phi[j, i] = p
                    J[j, i] = mu * p
                else:
                    p = rhs
                    for _ in range(60):
                        dp = (p - d * math.sin(p) - rhs) / (1.0 - d * math.cos(p))
                        p -= dp
                        if abs(dp) < 1e-13:
                            break
                    phi[j, i] = p
                    J[j, i] = mu * math.sin(p)
        return np.ascontiguousarray(phi.T), cnt, hist


def check(rho=20.0, X=5.0, T=4.0, seed=3):
    from diag import c0_card as C
    from diag.c0_links4 import sprinkle4
    from diag.c0_minvar import NMAX_ALL
    t, x, y, z = sprinkle4(rho, X, T, seed)
    RS = C.rows(rho)
    W = np.array([r[1] for r in RS])
    ctd = np.array([r[2] for r in RS])
    phi0 = np.array([m[3] * np.cos(m[2] * x) * np.cos(m[2] * t) for _ in RS for m in C.MODES])
    lin = np.array([m[1] for m in C.MODES])
    c = 1.0 / (2.0 * math.pi * math.sqrt(6.0 * rho))
    tc = (24.0 * (NMAX_ALL + 6) / (math.pi * rho)) ** 0.25
    sample = (t > 2.0)
    out = []
    for f in (C.evolve_card, evolve_fast):
        f(t[:50], x[:50], y[:50], z[:50], tc, NMAX_ALL, W, ctd, phi0[:, :50], lin, c, 1.0, sample[:50], 0.25, 10)
        a = time.time()
        r = f(t, x, y, z, tc, NMAX_ALL, W, ctd, phi0, lin, c, 1.0, sample, 0.25, 10)
        out.append((r, time.time() - a))
    (p1, c1, h1), s1 = out[0]
    (p2, c2, h2), s2 = out[1]
    print(f"N {t.size}: evolve_card {s1:.1f}s, evolve_fast {s2:.1f}s; identical: phi {np.array_equal(p1, p2)}, "
          f"cnt {np.array_equal(c1, c2)}, hist {np.array_equal(h1, h2)}")


if __name__ == "__main__":
    check()
