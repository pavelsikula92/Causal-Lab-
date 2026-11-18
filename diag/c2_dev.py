"""
C2, development: walls of the C0 field in a prescribed curved 3+1 order of the theory (decides nothing)
=====================================================================================================
Direction C (the author, 18 Nov): the motion of a body must follow from the light cones and the substrate alone. The
body is a wall of the field that the order itself carries (C0: the retarded Green function of the order, MV8 + the
local counterterm, F3.16b; on a curved order F5.3), no rule of motion is put in.

Order. The metric of the theory, ds^2 = -Pi^2 dt^2 + Pi^-2 (dx^2 + dy^2 + dz^2), Pi^2 = (x + ell)/ell = H/ell. Its
optical metric ell^2 (dH^2 + dy^2 + dz^2)/H^2 is hyperbolic 3-space of radius ell, so
    u < v  <=>  cosh((t_v - t_u)/ell) - 1 > ((H_u - H_v)^2 + dy^2 + dz^2) / (2 H_u H_v).
The wall is uniform in y and z, and translations in y and z are isometries, so y and z are periodic with period L
(minimal image; the relation of the quotient order). L > 2 T is required: then no pair has two related lifts and the
quotient order is, within the run, the plane order with periodic sources (L = 3 at T = 5 gave a field noise growing to
0.86 by t 4, an artefact of the wrapped cones). Coordinate density rho Pi^-2 = rho ell/H; the domain of dependence
of the initial interval -xl <= x <= xr: t <= ell ln(H2/H), ell ln(H/H1). ell = inf: the flat order (B1/C0 limit).

Rule. As diag/c0_fast.evolve_fast with the relation above:
    phi(v) - delta J(v) = phi0(v) - c sum_{u < v, n(u, v) <= 8} w_n J(u),   J = mu sin phi,
MV8 weights (diag/c0_minvar.minvar(8)), c = 1/(2 pi sqrt(6 rho)), delta of diag/c0_apex.delta. The candidate set is cut
on the approximate proper time tau^2 ~ Pi_u Pi_v (dt^2 - d_h^2) <= tc^2 with rho V(tc) = 8 + 6 (the margin of C0) and a
safety factor; the interval counts use the exact relation.

Continuum (y- and z-independent): box_g phi = mu sin phi reduces to
    phi_tt = Pi^4 phi_xx - Pi^2 mu sin phi .
The static energy of a wall is int [phi_x^2/2 + Pi^-2 V] dx = int Pi^-1 [phi_s^2/2 + V] ds: a wall whose tension equals
its energy (the kink) is pushed towards larger Pi with the acceleration of a test particle reversed; the passive weight
of a wall in 3+1 is E - 2 sigma (docs/derivations/walls_in_curved_order.md, section 3).

    python -m diag.c2_dev --ell inf --rho 30 --seeds 0,1 --arms rest       (flat torus, kink at rest)
"""
from __future__ import annotations

import argparse
import json
import math
import os
import time

import numpy as np

from diag.c0_apex import delta
from diag.c0_minvar import minvar
from phases.F3_structure import field1d as F

_nb = F._nb
INF = float("inf")
KMAX = 8
OUT = os.path.join(os.path.dirname(__file__), "out", "c2_dev")


# ------------------------------------------------------------------------------------------------ order

def inv_pi2(x, ell):
    return np.ones_like(x) if not np.isfinite(ell) else ell / (x + ell)


def dist_x(x, ell, xl, xr):
    """Optical distance from x to the complement of the initial interval (y, z periodic)."""
    if not np.isfinite(ell):
        return np.minimum(xr - x, x + xl)
    H = x + ell
    return np.minimum(ell * np.log((ell + xr) / H), ell * np.log(H / (ell - xl)))


def sprinkle(rho, ell, xl, xr, L, T, seed):
    rng = np.random.default_rng(seed)
    wmax = 1.0 if not np.isfinite(ell) else ell / (ell - xl)
    n = int(rng.poisson(rho * wmax * T * (xl + xr) * L * L))
    t = rng.uniform(0, T, n)
    x = rng.uniform(-xl, xr, n)
    y = rng.uniform(0, L, n)
    z = rng.uniform(0, L, n)
    keep = (rng.uniform(0, 1, n) * wmax < inv_pi2(x, ell)) & (t <= dist_x(x, ell, xl, xr))
    t, x, y, z = t[keep], x[keep], y[keep], z[keep]
    o = np.argsort(t, kind="stable")
    return t[o], x[o], y[o], z[o]


def proper_volume(ell, xl, xr, L, T, n=20000):
    xs = np.linspace(-xl, xr, n + 1)
    xm = 0.5 * (xs[1:] + xs[:-1])
    h = np.clip(dist_x(xm, ell, xl, xr), 0.0, T)
    return float(np.sum(h * inv_pi2(xm, ell)) * (xs[1] - xs[0]) * L * L)


if _nb is not None:
    @_nb.njit(cache=True, inline="always")
    def _mi(d, L):
        if d > 0.5 * L:
            return d - L
        if d < -0.5 * L:
            return d + L
        return d

    @_nb.njit(cache=True, inline="always")
    def _prec(t1, x1, y1, z1, t2, x2, y2, z2, ell, L):
        """1 < 2; returns the approximate tau^2 (> 0) or -1."""
        dt = t2 - t1
        if dt <= 0.0:
            return -1.0
        dy = _mi(y2 - y1, L)
        dz = _mi(z2 - z1, L)
        D2 = (x2 - x1) ** 2 + dy * dy + dz * dz
        if ell <= 0.0:
            q = dt * dt - D2
            return q if q > 0.0 else -1.0
        H1 = x1 + ell
        H2 = x2 + ell
        a = D2 / (2.0 * H1 * H2)
        if math.cosh(dt / ell) - 1.0 <= a:
            return -1.0
        dh = ell * math.acosh(1.0 + a)
        return math.sqrt(H1 * H2) / ell * (dt * dt - dh * dh)

    @_nb.njit(cache=True)
    def evolve(t, x, y, z, ell, L, tc2, W, ctd, phi0, c, mu):
        N = t.size
        R = phi0.shape[0]
        P0 = np.ascontiguousarray(phi0.T)
        phi = np.empty((N, R))
        J = np.empty((N, R))
        s = np.zeros(R)
        cand = np.empty(N, np.int64)
        cnt = np.zeros(KMAX + 2, np.int64)
        for j in range(N):
            m = 0
            for k in range(j - 1, -1, -1):
                q = _prec(t[k], x[k], y[k], z[k], t[j], x[j], y[j], z[j], ell, L)
                if q > 0.0 and q <= tc2:
                    cand[m] = k
                    m += 1
            for i in range(R):
                s[i] = 0.0
            for a in range(m):
                u = cand[a]
                nn = 0
                for b in range(a):
                    wv = cand[b]
                    if _prec(t[u], x[u], y[u], z[u], t[wv], x[wv], y[wv], z[wv], ell, L) > 0.0:
                        nn += 1
                        if nn > KMAX:
                            break
                if nn <= KMAX:
                    cnt[nn] += 1
                    wg = W[nn]
                    if wg != 0.0:
                        for i in range(R):
                            s[i] += wg * J[u, i]
                else:
                    cnt[KMAX + 1] += 1
            for i in range(R):
                rhs = P0[j, i] - c * s[i]
                d = ctd * mu
                p = rhs
                for _ in range(60):
                    dp = (p - d * math.sin(p) - rhs) / (1.0 - d * math.cos(p))
                    p -= dp
                    if abs(dp) < 1e-13:
                        break
                phi[j, i] = p
                J[j, i] = mu * math.sin(p)
        return np.ascontiguousarray(phi.T), cnt


def tc2_of(rho, safety=1.15):
    tc = (24.0 * (KMAX + 6) / (math.pi * rho)) ** 0.25
    return (safety * tc) ** 2


# ------------------------------------------------------------------------------------------------ continuum

if _nb is not None:
    @_nb.njit(cache=True)
    def _leap(phi0, phit0, P, mu, dx, dt, nsteps, every, i0, i1):
        n = phi0.size
        P4 = P ** 4
        P2 = P ** 2
        nout = nsteps // every + 1
        out = np.empty((nout, i1 - i0))
        prev = phi0.copy()
        rhs = np.zeros(n)
        for i in range(1, n - 1):
            rhs[i] = P4[i] * (phi0[i + 1] - 2.0 * phi0[i] + phi0[i - 1]) / (dx * dx) - P2[i] * mu * math.sin(phi0[i])
        cur = phi0 + dt * phit0 + 0.5 * dt * dt * rhs
        cur[0] = cur[1]
        cur[n - 1] = cur[n - 2]
        out[0] = phi0[i0:i1]
        k = 1
        for st in range(1, nsteps + 1):
            if st % every == 0:
                out[k] = cur[i0:i1]
                k += 1
            nxt = np.empty(n)
            for i in range(1, n - 1):
                lap = P4[i] * (cur[i + 1] - 2.0 * cur[i] + cur[i - 1]) / (dx * dx)
                nxt[i] = 2.0 * cur[i] - prev[i] + dt * dt * (lap - P2[i] * mu * math.sin(cur[i]))
            nxt[0] = nxt[1]
            nxt[n - 1] = nxt[n - 2]
            prev = cur
            cur = nxt
        return out


class Continuum:
    """y-, z-independent continuum solution; .at(t, x) interpolates, .centre(t) the zero of phi - pi near x0."""

    def __init__(self, ell, mu, v, x0, T, xlo, xhi, dx=0.005, dt_out=0.01, pad=40.0, omega=None, mu_dyn=None):
        self.ell, self.mu, self.v, self.x0, self.omega = ell, mu, v, x0, omega
        lo = (-ell + 1.0) if np.isfinite(ell) else xlo - pad
        lo = max(lo, xlo - pad)
        hi = xhi + pad
        xg = np.arange(lo, hi + 0.5 * dx, dx)
        P = 1.0 / np.sqrt(inv_pi2(xg, ell))
        g = 1.0 / math.sqrt(1.0 - v * v)
        s = 2.0 * math.sqrt(ell) * (np.sqrt(xg + ell) - math.sqrt(x0 + ell)) if np.isfinite(ell) else xg - x0
        ms = mu if mu > 0 else 1.0
        if omega is None:
            xi = math.sqrt(ms) * g * s
            phi0 = 4.0 * np.arctan(np.exp(np.clip(xi, -50, 50)))
            phit0 = -P * math.sqrt(ms) * g * v * 2.0 / np.cosh(np.clip(xi, -50, 50))
        else:
            k = math.sqrt(ms * (1.0 - omega * omega))
            phi0 = 4.0 * np.arctan(k / (math.sqrt(ms) * omega) / np.cosh(np.clip(k * s, -50, 50)))
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

    def centre(self, tm, level=math.pi):
        i = int(round(tm / self.dt_out))
        row = self.grid[i] - level
        j = np.nonzero(np.diff(np.sign(row)))[0]
        if j.size == 0:
            return float("nan")
        j = j[np.argmin(np.abs(self.xg[j] - self.x0))]
        return float(self.xg[j] - row[j] * self.dx / (row[j + 1] - row[j]))


# ------------------------------------------------------------------------------------------------ development task

ARMS = {"rest": dict(v=0.0, omega=None), "launch": dict(v=0.3, omega=None), "breather": dict(v=0.0, omega=0.7)}


def task(rho, ell, seed, arms, xl, xr, L, T, mu=1.0, base=21000):
    t, x, y, z = sprinkle(rho, ell, xl, xr, L, T, base + seed)
    V = proper_volume(ell, xl, xr, L, T)
    conts, exs, P0 = {}, {}, []
    for a in arms:
        free = Continuum(ell, mu, ARMS[a]["v"], 0.0, T, -xl, xr, omega=ARMS[a]["omega"], mu_dyn=0.0)
        full = Continuum(ell, mu, ARMS[a]["v"], 0.0, T, -xl, xr, omega=ARMS[a]["omega"])
        P0.append(free.at(t, x))
        conts[a] = full
    phi0 = np.array(P0)
    w = np.asarray(minvar(KMAX), float)
    assert np.all(w[KMAX + 1:] == 0.0)
    W = np.ascontiguousarray(w[:KMAX + 1])
    ctd = delta(w, rho)
    c = 1.0 / (2.0 * math.pi * math.sqrt(6.0 * rho))
    a0 = time.time()
    phi, cnt = evolve(t, x, y, z, ell if np.isfinite(ell) else -1.0, L, tc2_of(rho), W, ctd, phi0, c, mu)
    sec = time.time() - a0
    out = {"rho": rho, "ell": ell, "seed": seed, "N": int(t.size), "Nexp": rho * V, "sec": sec, "cnt": cnt.tolist(),
           "xl": xl, "xr": xr, "L": L, "T": T, "arms": {}}
    tb = np.arange(0.5, T - 0.9, 0.5)
    for i, a in enumerate(arms):
        ex = conts[a].at(t, x)
        rows = []
        for b in tb:
            g = round(conts[a].centre(b), 3) + 0.0   # a guess of ~1e-17 freezes least_squares (relative step)
            r = F.fit_kink(t, x, phi[i], b - 0.25, b + 0.25, g) if ARMS[a]["omega"] is None else None
            r_ex = F.fit_kink(t, x, ex, b - 0.25, b + 0.25, g) if ARMS[a]["omega"] is None else None
            sl = (np.abs(t - b) < 0.25) & (np.abs(x - conts[a].centre(b) if ARMS[a]["omega"] is None else x) < 6)
            dev = phi[i][sl] - ex[sl]
            rows.append({"t": float(b), "fit": r, "fit_ex": r_ex, "c_ex": conts[a].centre(b),
                         "dev_mean": float(dev.mean()) if dev.size else None,
                         "dev_sd": float(dev.std()) if dev.size else None, "n": int(sl.sum())})
        out["arms"][a] = rows
        np.save(os.path.join(OUT, f"phi_{a}_{rho:g}_{ell:g}_{seed}.npy"), np.vstack([t, x, phi[i], ex]).astype(np.float32))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rho", type=float, default=30.0)
    ap.add_argument("--ell", type=float, default=INF)
    ap.add_argument("--seeds", default="0")
    ap.add_argument("--arms", default="rest")
    ap.add_argument("--xl", type=float, default=13.0)
    ap.add_argument("--xr", type=float, default=22.0)
    ap.add_argument("--L", type=float, default=4.0)
    ap.add_argument("--T", type=float, default=8.0)
    ap.add_argument("--workers", type=int, default=1)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    arms = a.arms.split(",")
    for sd in [int(s) for s in a.seeds.split(",")]:
        r = task(a.rho, a.ell, sd, arms, a.xl, a.xr, a.L, a.T)
        f = os.path.join(OUT, f"c2_{a.rho:g}_{a.ell:g}_{a.L:g}_{a.T:g}_{sd}_{'-'.join(arms)}.json")
        json.dump(r, open(f, "w"), indent=1)
        print(f"seed {sd}: N {r['N']} (exp {r['Nexp']:.0f}) {r['sec']:.0f} s; cnt {r['cnt']}", flush=True)
        for arm, rows in r["arms"].items():
            for row in rows:
                f1, f2 = row["fit"], row["fit_ex"]
                print(f"  {arm} t {row['t']:.2f} xc {f1['xc'] if f1 else float('nan'):+.4f} "
                      f"xc_ex {f2['xc'] if f2 else float('nan'):+.4f} c_ex {row['c_ex']:+.4f} "
                      f"w {f1['w'] if f1 else float('nan'):.3f} dev {row['dev_mean']:+.4f} sd {row['dev_sd']:.4f} n {row['n']}",
                      flush=True)


if __name__ == "__main__":
    main()
