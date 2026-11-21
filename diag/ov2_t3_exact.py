"""
OV2 sec. 14, T3 pilot stage 1' (development): exact synthetic tests of the kernels of diag/ov2_t3_kernels.py against
brute force on hand-built posets and a small random sprinkling (decision 41: the weights f(n)(dphi)^2 are checked
here, not by F-b).

Checked, each to rounding:
  1. the relation and the interval counts n_xy (hand-built: a diamond and a chain with n known by hand);
  2. S_m and E_x = dS_m/dphi_x (central finite differences of a brute-force S_m);
  3. the chain rule on one sprinkling: B - P = sum_x Phi_t(x) dS_chi/dphi_x (finite differences of S_chi);
  4. R = Q2 + Q1 = P - sum_pairs (chi'_a chi_b + chi_a chi'_b) w L - sum_x chi'_x U_x (brute force);
  5. the fast 2D path noether_multi (Fenwick interval counts) equals the bitset kernel.

    python -m diag.ov2_t3_exact
"""
import math

import numpy as np

from diag.ov2_t3_kernels import (chi_of, fker, future_bits, interval_count, interval_counts_fast, noether,
                                 noether_multi, past_bits)


def brute(t, z, phi, phit, rho, eps, m, t1, t2, d):
    n = t.size
    rel = [[t[y] - t[x] > abs(z[y] - z[x]) for y in range(n)] for x in range(n)]
    nn = {}
    for x in range(n):
        for y in range(n):
            if rel[x][y]:
                nn[(x, y)] = sum(1 for w in range(n) if rel[x][w] and rel[w][y])
    ch = [chi_of(float(tt), t1, t2, d) for tt in t]

    def S_of(ph, windowed):
        s = 0.0
        for (x, y), k in nn.items():
            wgt = ch[x][0] * ch[y][0] if windowed else 1.0
            s += wgt * float(fker(k, eps)) * (-(eps ** 2) * (ph[x] - ph[y]) ** 2)
        for x in range(n):
            wgt = ch[x][0] if windowed else 1.0
            s += wgt * (-(m * m * (1 - math.cos(ph[x]))) / rho)
        return s

    h = 1e-5
    E = np.zeros(n)
    dSchi = np.zeros(n)
    for x in range(n):
        p1 = phi.copy(); p1[x] += h
        p0 = phi.copy(); p0[x] -= h
        E[x] = (S_of(p1, False) - S_of(p0, False)) / (2 * h)
        dSchi[x] = (S_of(p1, True) - S_of(p0, True)) / (2 * h)
    chain = float(sum(phit[x] * dSchi[x] for x in range(n)))      # dS_chi/dphi_x already contains chi_x
    Rb = 0.0
    for (x, y), k in nn.items():
        wv = float(fker(k, eps))
        D = phi[x] - phi[y]
        L = -(eps ** 2) * D * D
        Rb += ch[x][0] * (1 - ch[y][0]) * wv * (-2 * eps ** 2 * D) * phit[x]
        Rb += ch[y][0] * (1 - ch[x][0]) * wv * (2 * eps ** 2 * D) * phit[y]
        Rb -= (ch[x][1] * ch[y][0] + ch[x][0] * ch[y][1]) * wv * L
    for x in range(n):
        Rb -= ch[x][1] * (-(m * m * (1 - math.cos(phi[x]))) / rho)
    return nn, S_of(phi, False), E, chain, Rb


def check(name, t, z, phi, phit, rho=100.0, eps=0.12, m=1.0, t1=1.0, t2=2.0, d=0.4, expect_n=None):
    o = np.argsort(t, kind="stable")
    t, z, phi, phit = t[o], z[o], phi[o], phit[o]
    P = past_bits(t, z)
    F = future_bits(P)
    B, Pp, Q2, Q1, S, E = noether(t, z, phi, phit, rho, eps, m, t1, t2, d, P, F)
    nn, Sb, Eb, chain, Rb = brute(t, z, phi, phit, rho, eps, m, t1, t2, d)
    ok_n = all(interval_count(P, x, y) == k and interval_counts_fast(P, F, x, y) == k for (x, y), k in nn.items())
    if expect_n is not None:
        inv = {int(i): j for j, i in enumerate(o)}
        ok_n &= all(nn.get((inv[a], inv[b]), -1) == k for (a, b), k in expect_n.items())
    errS = abs(S - Sb)
    errE = float(np.max(np.abs(E - Eb)))
    errC = abs((B - Pp) - chain)
    errR = abs((Q2 + Q1) - Rb)
    Bm, Pm, Q2m, Q1m, Sm = noether_multi(t, z, np.vstack([phi, phi]), np.vstack([phit, phit]), rho, eps, m, t1, t2, d)
    errM = max(abs(Bm[1] - B), abs(Pm[1] - Pp), abs(Q2m[1] - Q2), abs(Q1m[1] - Q1), abs(Sm[1] - S))
    scale = max(1.0, abs(Sb))
    ok = (ok_n and errS < 1e-12 * scale and errE < 1e-8 and errC < 1e-8 and errR < 1e-12 * max(1.0, abs(Rb))
          and errM < 1e-12 * scale)
    print(f"{name:28s} n_xy {'OK' if ok_n else 'FAIL'}; |dS| {errS:.1e}; max|dE| {errE:.1e}; |chain| {errC:.1e}; "
          f"|dR| {errR:.1e}; fast path {errM:.1e}  -> {'OK' if ok else 'FAIL'}", flush=True)
    return ok


def main():
    ok = True
    # a diamond: a < b, a < c, b < d, c < d, b || c; n_ad = 2, the links n = 0
    t = np.array([0.0, 1.0, 1.0, 2.0]); z = np.array([0.0, -0.5, 0.5, 0.0])
    phi = np.array([0.1, 0.7, 2.0, 3.0]); pt = np.array([0.3, -0.2, 0.5, 0.1])
    ok &= check("diamond (n_ad = 2)", t, z, phi, pt, t1=0.6, t2=1.4, d=0.5,
                expect_n={(0, 3): 2, (0, 1): 0, (0, 2): 0, (1, 3): 0, (2, 3): 0})
    # a chain of 5 along t: n(i, j) = j - i - 1
    t = np.arange(5.0) * 0.7; z = np.zeros(5)
    phi = np.array([0.0, 0.5, 1.7, 2.9, 3.1]); pt = np.array([0.1, 0.2, -0.3, 0.4, 0.0])
    ok &= check("chain of 5 (n = j - i - 1)", t, z, phi, pt, t1=0.8, t2=2.0, d=0.5,
                expect_n={(i, j): j - i - 1 for i in range(5) for j in range(i + 1, 5)})
    # random small sprinklings with a kink-like field
    rng = np.random.default_rng(7)
    for nrun, npts in enumerate((40, 70)):
        t = rng.uniform(0, 3, npts); z = rng.uniform(-2, 2, npts)
        phi = 4 * np.arctan(np.exp(1.3 * (z - 0.3 * t))) + 0.05 * rng.normal(size=npts)
        pt = rng.normal(size=npts)
        ok &= check(f"random sprinkling N {npts}", t, z, phi, pt, t1=1.0, t2=2.0, d=0.4)
    print("ALL OK" if ok else "SOME FAILED")


if __name__ == "__main__":
    main()
