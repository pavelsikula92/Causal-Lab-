"""
B0, development stage 3: longest chains from sources by dynamic programming over the Hasse diagram (decides nothing)
=====================================================================================================================
The author's plan (11 Nov 2026): do not compute the matrix of longest chains for all pairs; for every SOURCE event s
compute L_s(v) for its future events by dynamic programming over the links, in topological order:
    L_s(s) = 0,   L_s(v) = 1 + max over links u -> v with L_s(u) defined of L_s(u),
so a direct link has L = 1. Cost ~ S * E_H (S sources, E_H links), memory O(N) per source. The cost of the Hasse
diagram itself is measured separately.

This stage:
  1. builds the links of a 2+1 sprinkling with a numba version of the shell method of core.order.links (candidates in the
     near-cone shell rho V(k, j) <= lam, links = maximal candidates), and checks it against core.order.links.links_csr on
     a small box;
  2. computes L_s(v) from sources by the DP, and checks it against an independent longest chain inside the interval
     (dynamic programming over the events of the interval with the relation as the only input) on random pairs;
  3. measures the costs: links, and DP per source, against N.
Coordinates enter only to build the relation of the flat sprinkling (as everywhere in the laboratory) and to find the
pairs of the check.

    python -m diag.b0_chain_alg
"""

from __future__ import annotations

import math
import time

import numpy as np

from core.order.links import links_csr
from phases.F3_structure import field1d as F
from diag.b0_green2p1 import sprinkle_pyramid

_nb = F._nb

if _nb is not None:
    @_nb.njit(cache=True)
    def links3(t, x, y, tau_cut):
        """Past links (parents) of every event of a time-sorted 2+1 sprinkling: maximal elements of the near-cone shell
        tau(k, j) <= tau_cut. Returns CSR (ptr, idx), parents sorted by index."""
        N = t.size
        ptr = np.zeros(N + 1, np.int64)
        cap = 64 * N
        idx = np.empty(cap, np.int64)
        tot = 0
        tc2 = tau_cut * tau_cut
        cand = np.empty(N, np.int64)
        maxi = np.empty(N, np.int64)
        for j in range(N):
            m = 0
            for k in range(j - 1, -1, -1):                     # decreasing time
                dt = t[j] - t[k]
                if dt <= 0.0:
                    continue
                q = dt * dt - (x[j] - x[k]) ** 2 - (y[j] - y[k]) ** 2
                if q > 0.0 and q <= tc2:
                    cand[m] = k
                    m += 1
            nm = 0
            for a in range(m):                                   # candidates already in decreasing time
                k = cand[a]
                dead = False
                for b in range(nm):
                    mm = maxi[b]
                    dt = t[mm] - t[k]
                    if dt > 0.0 and dt * dt > (x[mm] - x[k]) ** 2 + (y[mm] - y[k]) ** 2:
                        dead = True
                        break
                if not dead:
                    maxi[nm] = k
                    nm += 1
            if tot + nm > cap:
                new = np.empty(2 * cap, np.int64)
                new[:tot] = idx[:tot]
                idx = new
                cap = 2 * cap
            par = np.sort(maxi[:nm])
            idx[tot:tot + nm] = par
            tot += nm
            ptr[j + 1] = tot
        return ptr, idx[:tot]

    @_nb.njit(cache=True)
    def chain_from(s, ptr, idx, N):
        """L_s(v) for all v (-1 where v is not in the future of s): DP in index (time) order over the parents."""
        L = np.full(N, -1, np.int64)
        L[s] = 0
        for v in range(s + 1, N):
            b = -1
            for a in range(ptr[v], ptr[v + 1]):
                u = idx[a]
                if u >= s and L[u] >= 0 and L[u] + 1 > b:
                    b = L[u] + 1
            L[v] = b
        return L

    @_nb.njit(cache=True)
    def chain_interval(t, x, y, s, v):
        """Independent check: longest chain from s to v among the events of the interval, relation only."""
        m = 0
        sel = np.empty(t.size, np.int64)
        for z in range(s, v + 1):
            if z == s or z == v:
                sel[m] = z
                m += 1
                continue
            d1 = t[z] - t[s]
            d2 = t[v] - t[z]
            if d1 > 0 and d2 > 0 and d1 * d1 > (x[z] - x[s]) ** 2 + (y[z] - y[s]) ** 2 and \
                    d2 * d2 > (x[v] - x[z]) ** 2 + (y[v] - y[z]) ** 2:
                sel[m] = z
                m += 1
        best = np.full(m, -1, np.int64)
        best[0] = 0
        for i in range(1, m):
            zi = sel[i]
            b = -1
            for j in range(i):
                if best[j] >= 0:
                    zj = sel[j]
                    dt = t[zi] - t[zj]
                    if dt > 0 and dt * dt > (x[zi] - x[zj]) ** 2 + (y[zi] - y[zj]) ** 2 and best[j] + 1 > b:
                        b = best[j] + 1
            best[i] = b
        return best[m - 1]


def main():
    if _nb is None:
        raise SystemExit("numba required")
    lam = 25.0
    print("1. links: numba shell method against core.order.links.links_csr (small box)")
    t, x, y = sprinkle_pyramid(50.0, 5.0, 5.0, 3.0, 7)
    tau_cut = (12.0 * lam / (math.pi * 50.0)) ** (1.0 / 3.0)
    p1, i1 = links3(t, x, y, tau_cut)
    X = np.column_stack([t, x, y])
    p2, i2 = links_csr(X, 50.0, lam=lam)
    same = np.array_equal(p1, p2) and np.array_equal(i1, i2)
    print(f"   N {t.size}, links {i1.size} ({i1.size / t.size:.1f} per event); identical to links_csr: {same}")
    print("2. L_s(v) by DP from sources against the independent interval chain (random pairs)")
    rng = np.random.default_rng(3)
    ok, tot = 0, 0
    for _ in range(40):
        s = int(rng.integers(0, t.size // 3))
        L = chain_from(s, p1, i1, t.size)
        fut = np.nonzero(L > 0)[0]
        for v in rng.choice(fut, min(10, fut.size), replace=False):
            tot += 1
            ok += int(L[v] == chain_interval(t, x, y, s, int(v)))
    print(f"   {ok} / {tot} pairs equal")
    print("3. costs against N (rho 100, pyramid T, X = T + 8, Y = T + 3)")
    for T in (3.0, 4.0, 5.0, 6.0):
        t, x, y = sprinkle_pyramid(100.0, T + 8.0, T + 3.0, T, 11)
        tc = (12.0 * lam / (math.pi * 100.0)) ** (1.0 / 3.0)
        a = time.time()
        p, i = links3(t, x, y, tc)
        tl = time.time() - a
        a = time.time()
        srcs = rng.integers(0, t.size // 2, 20)
        for s in srcs:
            chain_from(int(s), p, i, t.size)
        td = (time.time() - a) / 20
        print(f"   T {T:g}: N {t.size}, links {i.size} ({i.size / t.size:.1f} per event), build {tl:.1f} s, DP per "
              f"source {td * 1e3:.1f} ms -> all N sources {td * t.size / 60:.1f} min, sources in |x| < 4 only "
              f"{td * np.sum(np.abs(x) < 4) / 60:.1f} min")


if __name__ == "__main__":
    main()
