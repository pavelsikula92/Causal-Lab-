"""
C0, development: the retarded Green function of the 3+1 order from its links (decides nothing)
=============================================================================================
In 3+1 the retarded Green function of the wave operator lives on the light cone, G = delta(t - r)/(4 pi r). On a
sprinkling of density rho the expected sum over the past links of v of a smooth source,
    c sum_{u link v} J(u),   c = 1 / (2 pi sqrt(6) sqrt(rho)),
equals int G J d^4u (Johnston; the link probability exp(-rho V), V = pi tau^4/24, integrated across the cone gives
sqrt(6)/(2 sqrt(rho)) per unit of the cone's measure). The field rule of the laboratory becomes
    phi(v) = phi0(v) - c sum_{u link v} J(u),   J = mu sin phi   (the mass fills the inside of the cone itself).
The weight reads only the order (whether u -- v is a link); no chains, no coordinates.

Checks (development):
  1. calibration: c sum_links J against the continuum retarded integral for a smooth Gaussian source;
  2. vacuum and dispersion: the linear rule (J = mu phi) with plane-wave data A cos(k x), phi_t = 0, against
     A cos(k x) cos(omega t), omega^2 = k^2 + mu; growth of the vacuum (k = 0);
  3. density dependence and cost.

    python -m diag.c0_links4 --check
"""
from __future__ import annotations

import math
import time

import numpy as np

from phases.F3_structure import field1d as F

_nb = F._nb


def sprinkle4(rho, X, T, seed, Y=None):
    """Poisson sprinkling of the 3+1 pyramid |x| <= X - t, |y|, |z| <= Y - t (Y = X by default), time sorted."""
    Y = X if Y is None else Y
    rng = np.random.default_rng(seed)
    n = int(rng.poisson(rho * T * 2 * X * (2 * Y) ** 2))
    t = rng.uniform(0, T, n)
    x = rng.uniform(-X, X, n)
    y = rng.uniform(-Y, Y, n)
    z = rng.uniform(-Y, Y, n)
    keep = (np.abs(x) <= X - t) & (np.abs(y) <= Y - t) & (np.abs(z) <= Y - t)
    t, x, y, z = t[keep], x[keep], y[keep], z[keep]
    o = np.argsort(t, kind="stable")
    return t[o], x[o], y[o], z[o]


def tau_cut4(rho, lam=12.0):
    return (24.0 * lam / (math.pi * rho)) ** 0.25


if _nb is not None:
    @_nb.njit(cache=True)
    def links4(t, x, y, z, tc):
        """Past links of every event of a time-sorted 3+1 sprinkling (shell rho V <= lam, exact maximality). CSR."""
        N = t.size
        ptr = np.zeros(N + 1, np.int64)
        cap = 256 * N
        idx = np.empty(cap, np.int64)
        tot = 0
        tc2 = tc * tc
        cand = np.empty(N, np.int64)
        maxi = np.empty(N, np.int64)
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
            nm = 0
            for a in range(m):
                k = cand[a]
                dead = False
                for b in range(nm):
                    mm = maxi[b]
                    dt = t[mm] - t[k]
                    if dt > 0.0 and dt * dt > (x[mm] - x[k]) ** 2 + (y[mm] - y[k]) ** 2 + (z[mm] - z[k]) ** 2:
                        dead = True
                        break
                if not dead:
                    maxi[nm] = k
                    nm += 1
            if tot + nm > cap:
                new = np.empty(2 * cap, np.int64)
                new[:tot] = idx[:tot]
                idx = new
                cap *= 2
            idx[tot:tot + nm] = np.sort(maxi[:nm])
            tot += nm
            ptr[j + 1] = tot
        return ptr, idx[:tot]

    @_nb.njit(cache=True)
    def evolve_links(ptr, idx, phi0, c, mu, linear):
        """phi(v) = phi0(v) - c sum_{u link v} J(u), J = mu phi (linear) or mu sin phi."""
        N = phi0.size
        phi = np.empty(N)
        for v in range(N):
            s = 0.0
            for e in range(ptr[v], ptr[v + 1]):
                u = idx[e]
                s += mu * phi[u] if linear else mu * math.sin(phi[u])
            phi[v] = phi0[v] - c * s
        return phi


def retarded_integral(J, tv, xv, yv, zv, nr=200, ndir=4000, seed=0):
    """(1/4 pi) int J(tv - r, xv + r n) r dr dOmega over 0 < r < tv (Monte Carlo directions, trapezoid in r)."""
    rng = np.random.default_rng(seed)
    n = rng.normal(size=(ndir, 3))
    n /= np.linalg.norm(n, axis=1)[:, None]
    rs = np.linspace(0, tv, nr)
    vals = np.array([np.mean(J(tv - r, xv + r * n[:, 0], yv + r * n[:, 1], zv + r * n[:, 2])) * r for r in rs])
    return float(np.trapezoid(vals, rs))


def check(rho_list=(30.0, 60.0), X=5.0, T=3.0):
    for rho in rho_list:
        a = time.time()
        t, x, y, z = sprinkle4(rho, X, T, 11)
        tc = tau_cut4(rho)
        ptr, idx = links4(t, x, y, z, tc)
        sl = time.time() - a
        c = 1.0 / (2.0 * math.pi * math.sqrt(6.0) * math.sqrt(rho))
        print(f"rho {rho:g}: N {t.size}, links per event {idx.size / t.size:.0f} (last 10 %: "
              f"{np.mean(np.diff(ptr)[-t.size // 10:]):.0f}), links {sl:.0f}s")
        # 1. calibration with a Gaussian source centred at (t 0.8, 0, 0, 0), width 0.5
        s0 = 0.5
        J = lambda tt, xx, yy, zz: np.exp(-((tt - 0.8) ** 2 + xx ** 2 + yy ** 2 + zz ** 2) / (2 * s0 * s0))
        Jv = J(t, x, y, z)
        S = np.array([c * Jv[idx[ptr[v]:ptr[v + 1]]].sum() for v in range(t.size)])
        rows = []
        for tv, rv in ((2.2, 0.0), (2.6, 0.0), (2.6, 1.0), (2.8, 1.8)):
            m = (np.abs(t - tv) < 0.1) & (np.abs(np.sqrt(x * x + y * y + z * z) - rv) < 0.25)
            if m.sum() < 5:
                continue
            cont = np.mean([retarded_integral(J, t[i], x[i], y[i], z[i], ndir=800, nr=80) for i in
                            np.nonzero(m)[0][:12]])
            rows.append((tv, rv, float(S[m].mean()), float(S[m].std() / math.sqrt(m.sum())), cont, int(m.sum())))
        for r in rows:
            print(f"   calibration at t {r[0]}, r {r[1]}: discrete {r[2]:.4f} +- {r[3]:.4f} (n {r[5]}), continuum "
                  f"{r[4]:.4f}, ratio {r[2] / r[4]:.3f}")
        # 2. vacuum (k = 0) and dispersion
        mu = 1.0
        inner = (np.abs(x) < X - t - 0.5) & (np.abs(y) < X - t - 0.5) & (np.abs(z) < X - t - 0.5) & (t > 0.3)
        for k in (0.0, 1.0, 2.0):
            A = 0.05
            phi0 = A * np.cos(k * x) * np.cos(k * t)
            phi = evolve_links(ptr, idx, phi0, c, mu, True)
            om_c = math.sqrt(k * k + mu)
            best = None
            for om in np.linspace(0.6 * om_c, 1.4 * om_c, 801):
                basis = np.cos(k * x[inner]) * np.cos(om * t[inner])
                amp = float(basis @ phi[inner] / (basis @ basis))
                res = float(np.sum((phi[inner] - amp * basis) ** 2))
                if best is None or res < best[0]:
                    best = (res, om, amp)
            noise = math.sqrt(best[0] / inner.sum())
            print(f"   k {k}: omega/omega_c - 1 {best[1] / om_c - 1:+.4f}, amplitude/A {best[2] / A:.4f}, residual rms "
                  f"{noise / A:.3f} A (inner events {inner.sum()})")


if __name__ == "__main__":
    import sys
    if "--check" in sys.argv:
        check()
