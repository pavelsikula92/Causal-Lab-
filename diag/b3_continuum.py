"""
B3, development: the continuum of the growth with a changing relation (2+1, planar; decides nothing)
==================================================================================================
Metric of the theory ds^2 = -Pi^2 dt^2 + Pi^-2 (dx^2 + dy^2), L = ln Pi, everything independent of y.
Fields (separate, interacting only through the order): box_g phi_i = mu_i sin phi_i,
    d_t (Pi^-3 d_t phi) = d_x (Pi d_x phi) - Pi^-1 mu sin phi.
Response (R2; the rule L = -(kappa/rho) sum w q of the order): box_g L = kappa q,
load q = trace of T = 3V + X/2, X = g^ab d_a phi d_b phi = -Pi^-2 phi_t^2 + Pi^2 phi_x^2, summed over the fields
(static kink: q = 4V, integral 2E; dust: q ~ energy).

Two variants:
  exact  the PDEs above;
  rule   what the discrete rule represents when it uses FLAT free waves (the free wave of the curved metric is not
         known before the order grows): phi = phi0_flat - int G_g mu sin phi, and the load by the identity of
         decision 19, L = -kappa int G_g (3V - phi mu sin(phi)/2) + (kappa/4)(phi^2 - h_flat), h_flat the flat free
         wave of the data of phi^2. Solved as psi = phi - phi0_flat with box_g psi = mu sin(psi + phi0_flat), zero data;
         Lsum with box_g Lsum = kappa (3V - phi mu sin phi / 2); L = Lsum + (kappa/4) sum (phi^2 - h_flat).
"""
from __future__ import annotations

import math

import numpy as np


def kink(x, x0, mu, v=0.0):
    g = 1.0 / math.sqrt(1 - v * v)
    xi = np.clip(math.sqrt(mu) * g * (x - x0), -50, 50)
    return 4.0 * np.arctan(np.exp(xi)), -math.sqrt(mu) * g * v * 2.0 / np.cosh(xi)


def breather(x, x0, om, mu=1.0):
    k = math.sqrt(mu * (1 - om * om))
    return 4.0 * np.arctan(k / (math.sqrt(mu) * om) / np.cosh(np.clip(k * (x - x0), -50, 50))), np.zeros_like(x)


def dalembert(x, f0, t):
    """Flat free wave of data (f0, 0) on the grid at time t (linear interpolation, constant extension)."""
    return 0.5 * (np.interp(x - t, x, f0) + np.interp(x + t, x, f0))


class Bodies:
    """Initial data: a list of (kind, x0, mu, extra) with kind 'kink' (extra = v) or 'br' (extra = omega)."""

    def __init__(self, specs):
        self.specs = specs

    def data(self, x):
        out = []
        for kind, x0, mu, e in self.specs:
            f, ft = kink(x, x0, mu, e) if kind == "kink" else breather(x, x0, e, mu)
            out.append((f, ft, mu))
        return out


def ramp(t, t1, t2):
    """s(t): 0 before t1, 1 after t2, smoothstep between; returns s, s', s''."""
    if t2 <= t1:
        return 1.0, 0.0, 0.0
    u = min(max((t - t1) / (t2 - t1), 0.0), 1.0)
    h = t2 - t1
    return u * u * (3 - 2 * u), 6 * u * (1 - u) / h, (6 - 12 * u) / h ** 2 if 0 < u < 1 else 0.0


def solve(bodies, kappa, T=8.0, dx=0.01, X=40.0, dt_out=0.05, variant="exact", t1=0.0, t2=0.0):
    x = np.arange(-X, X + 0.5 * dx, dx)
    n = x.size
    D = bodies.data(x)
    nf = len(D)
    mus = [d[2] for d in D]
    f0 = [d[0].copy() for d in D]
    if variant == "exact":
        ph = [d[0].copy() for d in D]
        P = [d[1].copy() for d in D]           # Pi^-3 phi_t, Pi = 1 at t = 0
    else:
        if any(np.any(d[1] != 0) for d in D):
            raise ValueError("rule variants: data with phi_t = 0 only")
        ph = [np.zeros(n) for _ in D]           # psi
        P = [np.zeros(n) for _ in D]
        h0 = [d[0] ** 2 for d in D]
    L = np.zeros(n)
    M = np.zeros(n)
    dt = 0.3 * dx
    every = int(round(dt_out / dt))
    dt = dt_out / every
    nsteps = int(round(T / dt_out)) * every
    snaps = {"phi": [], "L": [], "t": []}

    def full(i, t):
        return ph[i] if variant == "exact" else ph[i] + dalembert(x, f0[i], t)

    def lap(f, Ph):
        out = np.zeros(n)
        out[1:-1] = (Ph[1:] * (f[2:] - f[1:-1]) - Ph[:-1] * (f[1:-1] - f[:-2])) / (dx * dx)
        return out

    t = 0.0
    Ltot = L
    for s in range(nsteps + 1):
        phis = [full(i, t) for i in range(nf)]
        sw, sw1, sw2 = ramp(t, t1, t2)
        if variant == "exact":
            Ltot = L
        elif variant == "rule":
            Ltot = L + 0.25 * kappa * sum(phis[i] ** 2 - dalembert(x, h0[i], t) for i in range(nf))
        else:                                   # rule_s: switched load, local term s f
            Ltot = L + 0.25 * kappa * sw * sum(phis[i] ** 2 for i in range(nf))
        if s % every == 0:
            snaps["phi"].append(np.array(phis))
            snaps["L"].append(Ltot.copy())
            snaps["t"].append(t)
        if s == nsteps:
            break
        Pi = np.exp(Ltot)
        Ph = 0.5 * (Pi[1:] + Pi[:-1])
        q = np.zeros(n)
        for i in range(nf):
            mu = mus[i]
            p = phis[i]
            if variant == "exact":
                pt = Pi ** 3 * P[i]
                px = np.gradient(p, dx)
                q += sw * (3 * mu * (1 - np.cos(p)) + 0.5 * (-pt * pt / Pi ** 2 + Pi ** 2 * px * px))
            elif variant == "rule":
                q += 3 * mu * (1 - np.cos(p)) - 0.5 * p * mu * np.sin(p)
            else:
                pt = Pi ** 3 * P[i] + (np.gradient(dalembert(x, f0[i], t + 1e-3) - dalembert(x, f0[i], t - 1e-3), 1.0) / 2e-3 if False else 0.0)
                ftot_t = 2 * p * (Pi ** 3 * P[i] + 0.5 * (dalembert(x, f0[i], t + 1e-3) - dalembert(x, f0[i], t - 1e-3)) / 1e-3)
                q += sw * (3 * mu * (1 - np.cos(p)) - 0.5 * p * mu * np.sin(p)) + 0.25 * (2 * sw1 * ftot_t + sw2 * p * p)
        for i in range(nf):
            P[i] += dt * (lap(ph[i], Ph) - mus[i] * np.sin(phis[i]) / Pi)
        M += dt * (lap(L, Ph) - kappa * q / Pi)
        for i in range(nf):
            ph[i] += dt * Pi ** 3 * P[i]
            ph[i][0], ph[i][-1] = ph[i][1], ph[i][-2]
        L += dt * Pi ** 3 * M
        L[0], L[-1] = L[1], L[-2]
        t += dt
    return x, np.array(snaps["t"]), np.array(snaps["phi"]), np.array(snaps["L"])


def bfit(x, row_t, rows, t_lo, t_hi, guess):
    """Profile-fit centre of a breather over a time slab, on the grid (as F5.2c fit_breather)."""
    from phases.F5_gravity.F5_02c_walls_curved_2p1 import fit_breather
    m = (row_t >= t_lo) & (row_t < t_hi)
    sel = np.abs(x - guess) < 5.0
    xs = x[sel][::5]
    tt = np.repeat(row_t[m], xs.size)
    xx = np.tile(xs, int(m.sum()))
    pp = rows[m][:, sel][:, ::5].ravel()
    return fit_breather(tt, xx, pp, t_lo, t_hi, guess)[0]


def kcentre(x, row, guess):
    r = row - math.pi
    j = np.nonzero(np.diff(np.sign(r)))[0]
    j = j[np.argmin(np.abs(x[j] - guess))]
    return float(x[j] - r[j] * (x[1] - x[0]) / (r[j + 1] - r[j]))


if __name__ == "__main__":
    import sys
    import time
    d = float(sys.argv[1]) if len(sys.argv) > 1 else 4.0
    for kap in (0.005, 0.01):
        res = {}
        for muA in (1.0, 0.25):
            for var in ("exact", "rule"):
                a = time.time()
                B = Bodies([("kink", -d / 2, muA, 0.0), ("br", d / 2, 1.0, 0.7)])
                x, tg, PH, LL = solve(B, kap, variant=var)
                B0 = solve(B, 0.0, variant=var)
                xb = [bfit(x, tg, PH[:, 1], a_, a_ + 1, d / 2) for a_ in range(7)]
                xb0 = [bfit(x, B0[1], B0[2][:, 1], a_, a_ + 1, d / 2) for a_ in range(7)]
                kc = [kcentre(x, PH[int(round((a_ + 0.5) / 0.05)), 0], -d / 2) for a_ in range(7)]
                fall = np.array(xb) - np.array(xb0)
                res[(muA, var)] = fall
                i75 = int(round(7.5 / 0.05))
                print(f"d {d} kappa {kap} muA {muA} {var:5s}: L(xA) {LL[i75][np.argmin(abs(x + d / 2))]:+.3f} "
                      f"L(xB) {LL[i75][np.argmin(abs(x - d / 2))]:+.3f}; breather fall (fit centre, kappa - 0) "
                      + " ".join(f"{v:+.3f}" for v in fall) + "; kink centre " + " ".join(f"{v:+.3f}" for v in kc)
                      + f" ({time.time() - a:.0f}s)")
        for var in ("exact", "rule"):
            print(f"   ratio of falls muA 1 / 0.25 ({var}): {res[(1.0, var)][-1] / res[(0.25, var)][-1]:.3f}; "
                  f"rule/exact muA 1: {res[(1.0, 'rule')][-1] / res[(1.0, 'exact')][-1]:.3f}")
