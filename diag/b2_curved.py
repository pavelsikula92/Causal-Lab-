"""
B2, development: a prescribed curved 2+1 order of the theory and the continuum reference (decides nothing)
=========================================================================================================
Decision 16 (B2: a wall in a prescribed curved order of the theory). The metric of the theory (H1, throughput,
decision 12) is
    ds^2 = -Pi^2 dt^2 + Pi^-2 (dx^2 + dy^2),
narrower cones (coordinate light speed Pi^2) and more events (sqrt(-g) = Pi^-1 per coordinate volume). The profile is
chosen so that the causal relation has a closed form:
    Pi^2 = (x + ell) / ell = H / ell,   H = x + ell,
a uniform-like field (g_tt linear in x). The metric is conformal to -dt^2 + h with the optical metric
h = ell^2 (dH^2 + dy^2) / H^2, the hyperbolic plane of radius ell, so
    u < v  <=>  t_v - t_u > d_h(u, v),   cosh(d_h / ell) = 1 + ((H_u - H_v)^2 + (y_u - y_v)^2) / (2 H_u H_v).
ell = inf is the flat order of B1.

The sprinkling is Poisson with density rho per proper volume (coordinate density rho Pi^-1) in the domain of dependence
D of the initial rectangle -xl <= x <= xr, |y| <= Y: t <= the optical distance to the complement, which is
ell ln(H2/H) and ell ln(H/H1) (the horocycles x = xr, -xl) and ell asinh((Y -+ y)/H) (the geodesics y = +-Y).

Links: the shell method of diag/b0_chain_alg.links3 with the exact relation above; the shell (rho V <= 25) is applied to
the approximate proper time tau^2 ~ Pi(u) Pi(v) (dt^2 - d_h^2) with a safety factor (only the candidate set is
approximate; maximality uses the exact relation). Checked against the brute-force transitive reduction (check()).

Continuum: for data independent of y, box_g phi = mu sin phi reduces to
    phi_tt = Pi^3 d_x(Pi phi_x) - Pi^2 mu sin phi,
solved by a second-order leapfrog in flux form. Initial data: the kink in proper distance s(x) = int_x0^x Pi^-1 dx'
= 2 sqrt(ell) (sqrt(x + ell) - sqrt(x0 + ell)), launched at proper velocity v in the local frame:
    phi(0) = 4 atan(exp(sqrt(mu) gamma s)),  phi_t(0) = -Pi sqrt(mu) gamma v 2 / cosh(sqrt(mu) gamma s).
The same solver with mu = 0 gives the free wave phi_0 of the discrete rule (the initial data are the only place where
coordinates enter the dynamics, as in B1).

    python -m diag.b2_curved --check
"""

from __future__ import annotations

import math

import numpy as np

from phases.F3_structure import field1d as F

_nb = F._nb
INF = float("inf")


def _H(x, ell):
    return x + ell


if _nb is not None:
    @_nb.njit(cache=True, inline="always")
    def _rel(t1, x1, y1, t2, x2, y2, ell):
        """1 < 2 in the order of the metric with Pi^2 = (x + ell)/ell (ell <= 0: flat)."""
        dt = t2 - t1
        if dt <= 0.0:
            return False
        if ell <= 0.0:
            return dt * dt > (x2 - x1) ** 2 + (y2 - y1) ** 2
        H1 = x1 + ell
        H2 = x2 + ell
        return math.cosh(dt / ell) - 1.0 > ((x2 - x1) ** 2 + (y2 - y1) ** 2) / (2.0 * H1 * H2)

    @_nb.njit(cache=True)
    def links_curved(t, x, y, ell, tau_cut):
        """Past links of every event of a time-sorted sprinkling in the curved order (ell <= 0: flat). CSR (ptr, idx)."""
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
            for k in range(j - 1, -1, -1):
                dt = t[j] - t[k]
                if dt <= 0.0:
                    continue
                if ell <= 0.0:
                    q = dt * dt - (x[j] - x[k]) ** 2 - (y[j] - y[k]) ** 2
                    if q > 0.0 and q <= tc2:
                        cand[m] = k
                        m += 1
                else:
                    H1 = x[k] + ell
                    H2 = x[j] + ell
                    r = ((x[j] - x[k]) ** 2 + (y[j] - y[k]) ** 2) / (2.0 * H1 * H2)
                    ch = math.cosh(dt / ell) - 1.0
                    if ch > r:
                        d = ell * math.acosh(1.0 + r)
                        q = math.sqrt(H1 * H2) / ell * (dt * dt - d * d)
                        if q <= tc2:
                            cand[m] = k
                            m += 1
            nm = 0
            for a in range(m):
                k = cand[a]
                dead = False
                for b in range(nm):
                    mm = maxi[b]
                    if _rel(t[k], x[k], y[k], t[mm], x[mm], y[mm], ell):
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
    def links_brute(t, x, y, ell):
        """Transitive reduction of the exact relation (O(N^3); small N only)."""
        N = t.size
        R = np.zeros((N, N), np.bool_)
        for i in range(N):
            for j in range(i + 1, N):
                R[i, j] = _rel(t[i], x[i], y[i], t[j], x[j], y[j], ell)
        cnt = 0
        out = np.zeros((N, N), np.bool_)
        for j in range(N):
            for i in range(j):
                if R[i, j]:
                    ok = True
                    for z in range(i + 1, j):
                        if R[i, z] and R[z, j]:
                            ok = False
                            break
                    if ok:
                        out[i, j] = True
                        cnt += 1
        return out, cnt


def tau_cut(rho, safety=1.15):
    return safety * (12.0 * 25.0 / (math.pi * rho)) ** (1.0 / 3.0)


def dist_boundary(x, y, ell, xl, xr, Y):
    """Optical distance from (x, y) to the complement of the initial rectangle."""
    if not np.isfinite(ell):
        return np.minimum.reduce([xr - x, x + xl, Y - y, Y + y])
    H = x + ell
    H1, H2 = ell - xl, ell + xr
    with np.errstate(invalid="ignore", divide="ignore"):
        d = np.minimum.reduce([ell * np.log(H2 / H), ell * np.log(H / H1), ell * np.arcsinh((Y - y) / H),
                               ell * np.arcsinh((Y + y) / H)])
    return np.where((H > H1) & (H < H2) & (np.abs(y) < Y), d, -1.0)


def inv_pi(x, ell):
    return np.ones_like(x) if not np.isfinite(ell) else np.sqrt(ell / (x + ell))


def sprinkle_curved(rho, ell, xl, xr, Y, T, seed):
    """Poisson sprinkling of proper density rho in the domain of dependence (time sorted)."""
    rng = np.random.default_rng(seed)
    wmax = 1.0 if not np.isfinite(ell) else math.sqrt(ell / (ell - xl))
    box = T * (xl + xr) * 2 * Y
    n = int(rng.poisson(rho * wmax * box))
    t = rng.uniform(0, T, n)
    x = rng.uniform(-xl, xr, n)
    y = rng.uniform(-Y, Y, n)
    keep = (rng.uniform(0, 1, n) * wmax < inv_pi(x, ell)) & (t <= dist_boundary(x, y, ell, xl, xr, Y))
    t, x, y = t[keep], x[keep], y[keep]
    o = np.argsort(t, kind="stable")
    return t[o], x[o], y[o]


def proper_volume(ell, xl, xr, Y, T, n=1200):
    xs = np.linspace(-xl, xr, n + 1)
    ys = np.linspace(-Y, Y, n + 1)
    xm, ym = 0.5 * (xs[1:] + xs[:-1]), 0.5 * (ys[1:] + ys[:-1])
    X_, Y_ = np.meshgrid(xm, ym, indexing="ij")
    h = np.clip(dist_boundary(X_, Y_, ell, xl, xr, Y), 0.0, T)
    return float(np.sum(h * inv_pi(X_, ell)) * (xs[1] - xs[0]) * (ys[1] - ys[0]))


# ------------------------------------------------------------------------------------------------ continuum

if _nb is not None:
    @_nb.njit(cache=True)
    def _leap(phi0, phit0, P, mu, dx, dt, nsteps, every, i0, i1):
        n = phi0.size
        Ph = 0.5 * (P[1:] + P[:-1])
        P3 = P ** 3
        P2 = P ** 2
        nout = nsteps // every + 1
        out = np.empty((nout, i1 - i0))
        prev = phi0.copy()
        rhs = np.zeros(n)
        for i in range(1, n - 1):
            rhs[i] = P3[i] * (Ph[i] * (phi0[i + 1] - phi0[i]) - Ph[i - 1] * (phi0[i] - phi0[i - 1])) / (dx * dx) \
                - P2[i] * mu * math.sin(phi0[i])
        cur = phi0 + dt * phit0 + 0.5 * dt * dt * rhs
        cur[0] = cur[1]
        cur[n - 1] = cur[n - 2]
        out[0] = phi0[i0:i1]
        k = 1
        for s in range(1, nsteps + 1):
            if s % every == 0:
                out[k] = cur[i0:i1]
                k += 1
            nxt = np.empty(n)
            for i in range(1, n - 1):
                lap = P3[i] * (Ph[i] * (cur[i + 1] - cur[i]) - Ph[i - 1] * (cur[i] - cur[i - 1])) / (dx * dx)
                nxt[i] = 2.0 * cur[i] - prev[i] + dt * dt * (lap - P2[i] * mu * math.sin(cur[i]))
            nxt[0] = nxt[1]
            nxt[n - 1] = nxt[n - 2]
            prev = cur
            cur = nxt
        return out


class Continuum:
    """The y-independent continuum solution on a grid; .at(t, x) interpolates, .centre(t) the zero of phi - pi."""

    def __init__(self, ell, mu, v, x0, T, xlo, xhi, dx=0.005, dt_out=0.01, pad_lo=None, pad_hi=60.0, omega=None,
                 mu_dyn=None):
        """omega: a breather wall of frequency omega at maximal amplitude (phi_t = 0) instead of the kink; mu_dyn: the
        mu of the dynamics (0 for the free wave of the same initial data)."""
        self.ell, self.mu, self.v, self.x0, self.omega = ell, mu, v, x0, omega
        lo = (-ell + 1.0) if np.isfinite(ell) else xlo - T - 2.0
        if pad_lo is not None:
            lo = max(lo, xlo - pad_lo)
        hi = xhi + pad_hi
        xg = np.arange(lo, hi + 0.5 * dx, dx)
        P = 1.0 / inv_pi(xg, ell)
        g = 1.0 / math.sqrt(1.0 - v * v)
        if np.isfinite(ell):
            s = 2.0 * math.sqrt(ell) * (np.sqrt(xg + ell) - math.sqrt(x0 + ell))
        else:
            s = xg - x0
        xi = math.sqrt(mu if mu > 0 else 1.0) * g * s
        self.mu_shape = mu if mu > 0 else 1.0
        if omega is None:
            phi0 = 4.0 * np.arctan(np.exp(np.clip(xi, -50, 50)))
            phit0 = -P * math.sqrt(self.mu_shape) * g * v * 2.0 / np.cosh(np.clip(xi, -50, 50))
        else:
            k = math.sqrt(self.mu_shape * (1.0 - omega * omega))
            phi0 = 4.0 * np.arctan(k / (math.sqrt(self.mu_shape) * omega) / np.cosh(np.clip(k * s, -50, 50)))
            phit0 = np.zeros_like(phi0)
        if mu_dyn is not None:
            mu = mu_dyn
        cmax = float(np.max(P * P))
        dt = 0.4 * dx / cmax
        every = max(1, int(round(dt_out / dt)))
        dt = dt_out / every
        nsteps = int(math.ceil(T / dt_out)) * every
        i0 = int(np.searchsorted(xg, xlo - 1.0))
        i1 = int(np.searchsorted(xg, xhi + 1.0))
        self.grid = _leap(phi0, phit0, P, float(mu), dx, dt, nsteps, every, i0, i1)
        self.xg = xg[i0:i1]
        self.tg = np.arange(self.grid.shape[0]) * dt_out
        self.dx, self.dt_out = dx, dt_out

    def at(self, t, x):
        ft = np.clip(t / self.dt_out, 0, self.tg.size - 1.000001)
        it = ft.astype(int)
        at = ft - it
        fx = np.clip((x - self.xg[0]) / self.dx, 0, self.xg.size - 1.000001)
        ix = fx.astype(int)
        ax = fx - ix
        G = self.grid
        return ((1 - at) * ((1 - ax) * G[it, ix] + ax * G[it, ix + 1])
                + at * ((1 - ax) * G[it + 1, ix] + ax * G[it + 1, ix + 1]))

    def centre(self, tm):
        i = int(round(tm / self.dt_out))
        row = self.grid[i] - math.pi
        j = np.nonzero(np.diff(np.sign(row)))[0]
        if j.size == 0:
            return float("nan")
        j = j[np.argmin(np.abs(self.xg[j] - self.x0))]
        return float(self.xg[j] - row[j] * self.dx / (row[j + 1] - row[j]))


def check():
    import time
    print("1. links (shell, exact relation) against the brute-force transitive reduction")
    for ell in (24.0, INF):
        t, x, y = sprinkle_curved(30.0, ell, 3.0, 3.0, 3.0, 2.5, 5)
        a = time.time()
        p, i = links_curved(t, x, y, ell if np.isfinite(ell) else -1.0, tau_cut(30.0))
        tl = time.time() - a
        B, cnt = links_brute(t, x, y, ell if np.isfinite(ell) else -1.0)
        same = cnt == i.size and all(B[k, j] for j in range(t.size) for k in i[p[j]:p[j + 1]])
        print(f"   ell {ell}: N {t.size}, links {i.size}, brute {cnt}; identical: {same} ({tl:.1f}s)")
    print("2. sprinkling density: count / (rho * proper volume)")
    for ell in (24.0, INF):
        ns = [sprinkle_curved(30.0, ell, 13.0, 22.0, 12.5, 8.0, s)[0].size for s in range(4)]
        V = proper_volume(24.0 if np.isfinite(ell) else ell, 13.0, 22.0, 12.5, 8.0)
        print(f"   ell {ell}: N {np.mean(ns):.0f}, ratio {np.mean(ns) / (30 * V):.4f} (Poisson sd {1 / math.sqrt(30 * V * 4):.4f})")
    print("3. continuum: convergence (dx 0.01 against 0.005) and the fall of the centre")
    for ell in (24.0, INF):
        a = time.time()
        C1 = Continuum(ell, 1.0, 0.0, 0.0, 8.0, -13.0, 22.0, dx=0.01)
        C2 = Continuum(ell, 1.0, 0.0, 0.0, 8.0, -13.0, 22.0, dx=0.005)
        tt = np.random.default_rng(1).uniform(0, 8, 5000)
        xx = np.random.default_rng(2).uniform(-8, 9, 5000)
        d = np.max(np.abs(C1.at(tt, xx) - C2.at(tt, xx)))
        cs = [C2.centre(tm) for tm in (0.5, 2.5, 4.5, 6.5)]
        print(f"   ell {ell}: max |phi(dx .01) - phi(dx .005)| {d:.1e}; centre at t 0.5 2.5 4.5 6.5: "
              + " ".join(f"{c:+.4f}" for c in cs) + f"; geodesic -t^2/(4 ell) at 6.5: "
              f"{-6.5 ** 2 / (4 * ell) if np.isfinite(ell) else 0:+.4f} ({time.time() - a:.0f}s)")


if __name__ == "__main__":
    import sys
    if "--check" in sys.argv:
        check()
