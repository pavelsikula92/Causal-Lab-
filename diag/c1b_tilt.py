"""
C1b, development: the theory metric with tilted cones, a uniform ball (14 Nov; decides nothing)
================================================================================================
g = -Pi^2 dt^2 + Pi^-2 dx^2 with N = Pi^-2 = 1 + 2 a g(r); g = (3 - r^2/R^2)/2 inside the ball, R/r outside (first-order
vacuum outside). Relation: u < v iff t_v - t_u >= T(x_u, x_v), the travel time of the true (bent) ray: the straight-ray
time |x_v - x_u| + 2 a Int g dl (closed form) minus the bending gain, tabulated from a path minimisation (gain_table). The density of events per coordinate volume is rho N. See docs/derivations/c1b_volume_first_order.md.

Exact mean (any relation that respects t):  f(v)(1 - delta mu) = A - c mu Int d^4u rho N(u) K(rho V(u, v)) f(u), with
    V(u, v) = Int d^3w N(w) (dt - tau(u, w) - tau(w, v))_+      (prolate coordinates around the pair).
With s' = dt^2 - tau(u, v)^2 the measure is rho N d^4u = rho pi l^2 / tau ds' dcos(theta) d(dt).
mode 0: the exact mean; mode 1: the same with V -> (pi/24) s'^2 / Nbar^2 (the flat formula in the geodesic proper time,
Nbar = tau / l; what a local kernel in the proper time would give); a = 0: flat.

    python -m diag.c1b_tilt --predict --rho 60 --a 0.05
    python -m diag.c1b_tilt --rho 60 --a 0.05 --seeds 0-7 --workers 8
"""
from __future__ import annotations

import argparse
import contextlib
import hashlib
import io
import json
import math
import os
import time

import numpy as np

from core.stats.parallel import parallel_map, default_workers
from diag.c0_apex import delta
from diag.c0_fast import evolve_fast
from diag.c0_links4 import sprinkle4
from diag.c0_minvar import NMAX_ALL, minvar
from diag.c1_sphere import solve_r
from phases.F3_structure import field1d as F

_nb = F._nb
R_B = 1.5
MU, AMP = 1.0, 0.05
LAM_CUT = 40.0                 # K(40) < 1e-9 for n <= 8: the cut never binds
TB_W, RB_W = 0.5, 0.5

if _nb is not None:
    @_nb.njit(cache=True)
    def _F(x, b2, R):
        ax = abs(x)
        h = math.sqrt(R * R - b2) if b2 < R * R else 0.0
        if ax <= h:
            val = 0.5 * (3.0 - b2 / (R * R)) * ax - ax ** 3 / (6.0 * R * R)
        else:
            fin = 0.5 * (3.0 - b2 / (R * R)) * h - h ** 3 / (6.0 * R * R)
            val = fin + R * math.log((ax + math.sqrt(ax * ax + b2)) / (h + math.sqrt(h * h + b2)))
        return val if x >= 0.0 else -val

    @_nb.njit(cache=True)
    def gseg(px, py, pz, qx, qy, qz, R):
        """Int g dl along the straight segment p -> q."""
        dx, dy, dz = qx - px, qy - py, qz - pz
        L = math.sqrt(dx * dx + dy * dy + dz * dz)
        if L < 1e-14:
            return 0.0
        s0 = -(px * dx + py * dy + pz * dz) / L
        b2 = px * px + py * py + pz * pz - s0 * s0
        if b2 < 0.0:
            b2 = 0.0
        return _F(L - s0, b2, R) - _F(-s0, b2, R)

    @_nb.njit(cache=True)
    def gpt(r, R):
        return 0.5 * (3.0 - (r / R) ** 2) if r < R else R / r

    @_nb.njit(cache=True)
    def tau(px, py, pz, qx, qy, qz, a, R):
        dx, dy, dz = qx - px, qy - py, qz - pz
        return math.sqrt(dx * dx + dy * dy + dz * dz) + 2.0 * a * gseg(px, py, pz, qx, qy, qz, R)

    @_nb.njit(cache=True)
    def _gain_one(r1, r2, gam, a, R, K):
        """tau_straight - (true Fermat time) for p = (0, 0, r1), q = r2 (sin gam, 0, cos gam): the path is a polyline with
        K interior nodes displaced in the plane; preconditioned descent (the Hessian of the length term)."""
        px, pz = 0.0, r1
        qx, qz = r2 * math.sin(gam), r2 * math.cos(gam)
        L = math.sqrt((qx - px) ** 2 + (qz - pz) ** 2)
        if L < 1e-9 or a == 0.0:
            return 0.0
        ex, ez = -(qz - pz) / L, (qx - px) / L
        h = L / (K + 1)
        y = np.zeros(K + 2)
        g = np.zeros(K)
        cp = np.zeros(K)
        dp = np.zeros(K)
        e = 1e-6
        for itr in range(200):
            for i in range(1, K + 1):
                bx, bz = px + (qx - px) * i / (K + 1), pz + (qz - pz) * i / (K + 1)
                lx, lz = px + (qx - px) * (i - 1) / (K + 1) + y[i - 1] * ex, pz + (qz - pz) * (i - 1) / (K + 1) + y[i - 1] * ez
                rx, rz = px + (qx - px) * (i + 1) / (K + 1) + y[i + 1] * ex, pz + (qz - pz) * (i + 1) / (K + 1) + y[i + 1] * ez
                yp, ym = y[i] + e, y[i] - e
                tp = tau(lx, 0.0, lz, bx + yp * ex, 0.0, bz + yp * ez, a, R) + tau(bx + yp * ex, 0.0, bz + yp * ez, rx, 0.0, rz, a, R)
                tm = tau(lx, 0.0, lz, bx + ym * ex, 0.0, bz + ym * ez, a, R) + tau(bx + ym * ex, 0.0, bz + ym * ez, rx, 0.0, rz, a, R)
                g[i - 1] = (tp - tm) / (2.0 * e)
            # solve (1/h) tridiag(-1, 2, -1) d = g
            bb = 2.0 / h
            cc = -1.0 / h
            cp[0] = cc / bb
            dp[0] = g[0] / bb
            for i in range(1, K):
                m = bb - cc * cp[i - 1]
                cp[i] = cc / m
                dp[i] = (g[i] - cc * dp[i - 1]) / m
            dmax = 0.0
            d = dp[K - 1]
            y[K] -= d
            dmax = abs(d)
            for i in range(K - 2, -1, -1):
                d = dp[i] - cp[i] * d
                y[i + 1] -= d
                if abs(d) > dmax:
                    dmax = abs(d)
            if dmax < 1e-10:
                break
        tt = 0.0
        for i in range(K + 1):
            lx, lz = px + (qx - px) * i / (K + 1) + y[i] * ex, pz + (qz - pz) * i / (K + 1) + y[i] * ez
            rx, rz = px + (qx - px) * (i + 1) / (K + 1) + y[i + 1] * ex, pz + (qz - pz) * (i + 1) / (K + 1) + y[i + 1] * ez
            tt += tau(lx, 0.0, lz, rx, 0.0, rz, a, R)
        gain = tau(px, 0.0, pz, qx, 0.0, qz, a, R) - tt
        return gain if gain > 0.0 else 0.0

    @_nb.njit(cache=True, parallel=True)
    def gain_table(nr, hr, ng, a, R, K):
        G = np.zeros((nr, nr, ng))
        hg = math.pi / (ng - 1)
        for i in _nb.prange(nr):
            for j in range(i, nr):
                for k in range(ng):
                    G[i, j, k] = _gain_one(i * hr, j * hr, k * hg, a, R, K)
        for i in range(nr):
            for j in range(i):
                for k in range(ng):
                    G[i, j, k] = G[j, i, k]
        return G

    @_nb.njit(cache=True)
    def tauT(px, py, pz, qx, qy, qz, a, R, G, hr):
        """The travel time of the true (bent) ray: the straight-ray time minus the tabulated gain."""
        ts = tau(px, py, pz, qx, qy, qz, a, R)
        r1 = math.sqrt(px * px + py * py + pz * pz)
        r2 = math.sqrt(qx * qx + qy * qy + qz * qz)
        nr = G.shape[0]
        ng = G.shape[2]
        f1, f2 = r1 / hr, r2 / hr
        i, j = int(f1), int(f2)
        if i >= nr - 1 or j >= nr - 1:
            return ts
        if r1 * r2 > 1e-12:
            cg = (px * qx + py * qy + pz * qz) / (r1 * r2)
            cg = 1.0 if cg > 1.0 else (-1.0 if cg < -1.0 else cg)
            fg = math.acos(cg) / (math.pi / (ng - 1))
        else:
            fg = 0.0
        k = int(fg)
        if k >= ng - 1:
            k = ng - 2
        a1, a2, a3 = f1 - i, f2 - j, fg - k
        g = 0.0
        for di in range(2):
            wi = a1 if di else 1.0 - a1
            for dj in range(2):
                wj = a2 if dj else 1.0 - a2
                g += wi * wj * ((1.0 - a3) * G[i + di, j + dj, k] + a3 * G[i + di, j + dj, k + 1])
        return ts - g

    @_nb.njit(cache=True)
    def volume(ux, uz, vz, dt, a, R, gx, gw, ge, gwe, nphi, G, hr):
        """V(u, v) for u = (ux, 0, uz), v = (0, 0, vz), t_v - t_u = dt >= tau(u, v)."""
        ex, ez = -ux, vz - uz
        l = math.sqrt(ex * ex + ez * ez)
        if l > 1e-12:
            ex /= l
            ez /= l
        else:
            ex, ez, l = 0.0, 1.0, 0.0
        mx, mz = 0.5 * ux, 0.5 * (uz + vz)
        # a1 = (ez, 0, -ex) in the plane, a2 = (0, 1, 0); the integrand is even in phi -> half the circle
        tot = 0.0
        for ie in range(ge.size):
            eta = ge[ie]
            for ip in range(nphi):
                ph = math.pi * (ip + 0.5) / nphi
                cp, sp = math.cos(ph), math.sin(ph)
                # h(s) = dt - tau(u, w) - tau(w, v), decreasing; root in (l, dt]
                s_lo, s_hi = l, dt
                h_lo, h_hi = 1.0, -1.0
                first = True
                s_star = dt
                for itr in range(8):
                    if first:
                        s = dt
                    else:
                        s = s_lo + (s_hi - s_lo) * h_lo / (h_lo - h_hi)
                    zeta = 0.5 * s * eta
                    rp = 0.5 * math.sqrt(max((s * s - l * l) * (1.0 - eta * eta), 0.0))
                    wx = mx + zeta * ex + rp * cp * ez
                    wy = rp * sp
                    wz = mz + zeta * ez - rp * cp * ex
                    h = dt - tauT(ux, 0.0, uz, wx, wy, wz, a, R, G, hr) - tauT(wx, wy, wz, 0.0, 0.0, vz, a, R, G, hr)
                    if first:
                        first = False
                        if h >= 0.0:
                            s_star = dt
                            break
                        h_hi = h
                        # h at s = l: w on the segment
                        h_lo = dt - tauT(ux, 0.0, uz, 0.0, 0.0, vz, a, R, G, hr)
                        if h_lo <= 0.0:
                            s_star = l
                            break
                        continue
                    s_star = s
                    if abs(h) < 1e-13 * (1.0 + dt):
                        break
                    if h > 0.0:
                        s_lo, h_lo = s, h
                    else:
                        s_hi, h_hi = s, h
                if s_star <= l:
                    continue
                acc = 0.0
                for k in range(gx.size):
                    s = l + (s_star - l) * 0.5 * (gx[k] + 1.0)
                    zeta = 0.5 * s * eta
                    rp = 0.5 * math.sqrt(max((s * s - l * l) * (1.0 - eta * eta), 0.0))
                    wx = mx + zeta * ex + rp * cp * ez
                    wy = rp * sp
                    wz = mz + zeta * ez - rp * cp * ex
                    h = dt - tauT(ux, 0.0, uz, wx, wy, wz, a, R, G, hr) - tauT(wx, wy, wz, 0.0, 0.0, vz, a, R, G, hr)
                    if h > 0.0:
                        nw = 1.0 + 2.0 * a * gpt(math.sqrt(wx * wx + wy * wy + wz * wz), R)
                        acc += gw[k] * nw * h * (s * s - l * l * eta * eta) / 8.0
                tot += acc * 0.5 * (s_star - l) * gwe[ie] * (2.0 * math.pi / nphi)
        return tot

    @_nb.njit(cache=True, parallel=True)
    def kernel_table(rs, dts, a, R, rho, tc2, w, c, mu, gs, gws, gc, gwc, gx, gw, ge, gwe, nphi, mode, G, hr):
        nr = rs.size
        nt = dts.size
        DR = rs[1] - rs[0]
        out = np.zeros((nr, nt, nr))
        nn = w.size
        for ir in _nb.prange(nr):
            rv = rs[ir]
            for it in range(1, nt):
                dt = dts[it]
                smax = min(dt * dt, tc2)
                for p in range(gs.size):
                    sig = 0.5 * smax * (gs[p] + 1.0)
                    wsig = 0.5 * smax * gws[p]
                    target = math.sqrt(dt * dt - sig)
                    for q in range(gc.size):
                        ct = gc[q]
                        st = math.sqrt(max(0.0, 1.0 - ct * ct))
                        l = target / (1.0 + 2.0 * a)
                        for _ in range(30):
                            ux, uz = l * st, rv + l * ct
                            f = tauT(ux, 0.0, uz, 0.0, 0.0, rv, a, R, G, hr) - target
                            nu_ = 1.0 + 2.0 * a * gpt(math.sqrt(ux * ux + uz * uz), R)
                            dl = f / nu_
                            l -= dl
                            if abs(dl) < 1e-13:
                                break
                        if l < 0.0:
                            l = 0.0
                        ux, uz = l * st, rv + l * ct
                        ru = math.sqrt(ux * ux + uz * uz)
                        if a == 0.0:
                            V = math.pi / 24.0 * sig * sig
                        elif mode == 1:
                            nbar = target / l if l > 1e-9 else 1.0 + 2.0 * a * gpt(rv, R)
                            V = math.pi / 24.0 * sig * sig / (nbar * nbar)
                        else:
                            V = volume(ux, uz, rv, dt, a, R, gx, gw, ge, gwe, nphi, G, hr)
                        lam = rho * V
                        pk = math.exp(-lam)
                        kk = 0.0
                        for n in range(nn):
                            if n > 0:
                                pk *= lam / n
                            kk += w[n] * pk
                        if target < 1e-12:
                            continue
                        val = c * mu * rho * math.pi * l * l / target * kk * wsig * gwc[q]
                        fj = ru / DR
                        j0 = int(math.floor(fj))
                        fr = fj - j0
                        if j0 < nr:
                            out[ir, it, j0] += val * (1.0 - fr)
                        if j0 + 1 < nr:
                            out[ir, it, j0 + 1] += val * fr
        return out

    @_nb.njit(cache=True)
    def evolve_tilt(t, x, y, z, a, R, tc2, kmax, w, ctd, phi0, c, mu, G, hr):
        """The C0 rule (one kernel, linear) on the straight-ray relation; events sorted by t."""
        N = t.size
        phi = np.empty(N)
        J = np.empty(N)
        cand = np.empty(N, np.int64)
        ct_ = np.empty(N)
        cx = np.empty(N)
        cy = np.empty(N)
        cz = np.empty(N)
        amax = (1.0 + 3.0 * a) ** 2 - 1.0
        for j in range(N):
            m = 0
            tj, xj, yj, zj = t[j], x[j], y[j], z[j]
            for k in range(j - 1, -1, -1):
                dt = tj - t[k]
                if dt <= 0.0:
                    continue
                l2 = (xj - x[k]) ** 2 + (yj - y[k]) ** 2 + (zj - z[k]) ** 2
                q = dt * dt - l2
                if q > 0.0 and q <= tc2 + amax * l2:
                    tu = tauT(x[k], y[k], z[k], xj, yj, zj, a, R, G, hr)
                    q2 = dt * dt - tu * tu
                    if dt >= tu and q2 <= tc2:
                        cand[m] = k
                        ct_[m] = t[k]
                        cx[m] = x[k]
                        cy[m] = y[k]
                        cz[m] = z[k]
                        m += 1
            s = 0.0
            for i in range(m):
                u = cand[i]
                tu, xu, yu, zu = ct_[i], cx[i], cy[i], cz[i]
                nn = 0
                for b in range(i):
                    dt = ct_[b] - tu
                    if dt > 0.0 and dt * dt > (cx[b] - xu) ** 2 + (cy[b] - yu) ** 2 + (cz[b] - zu) ** 2:
                        if dt >= tauT(xu, yu, zu, cx[b], cy[b], cz[b], a, R, G, hr):
                            nn += 1
                            if nn > kmax:
                                break
                if nn <= kmax:
                    s += w[nn] * J[u]
            p = (phi0[j] - c * s) / (1.0 - ctd * mu)
            phi[j] = p
            J[j] = mu * p
        return phi


def gfun(r):
    r = np.asarray(r, float)
    return np.where(r < R_B, 0.5 * (3 - (r / R_B) ** 2), R_B / np.maximum(r, 1e-12))


_GT = {}
HR_G, NG_G, K_G = 0.1, 65, 40


def gains(a):
    """The table of the bending gain (cached in memory and on disk)."""
    a = float(a)
    if a not in _GT:
        fn = os.path.join("diag", "out", f"c1b_gain_a{a}_R{R_B}_{HR_G}_{NG_G}_{K_G}.npy")
        if os.path.exists(fn):
            _GT[a] = np.load(fn)
        else:
            _GT[a] = gain_table(int(round(12.4 / HR_G)) + 1, HR_G, NG_G, a, R_B, K_G)
            os.makedirs(os.path.dirname(fn), exist_ok=True)
            np.save(fn, _GT[a])
    return _GT[a]


def mean_field(rho, a, mode=0, T=6.0, Rg=10.6, DR=0.1, DT=0.05, ns=24, nc=24, nx=6, ne=8, nphi=6):
    w = minvar(8)
    ctd = delta(w, rho)
    rs = np.arange(0.0, Rg + 1e-9, DR)
    ts = np.arange(0.0, T + 1e-9, DT)
    tc2 = math.sqrt(24.0 * LAM_CUT / (math.pi * rho))
    c = 1.0 / (2.0 * math.pi * math.sqrt(6.0 * rho))
    gs, gws = np.polynomial.legendre.leggauss(ns)
    gc, gwc = np.polynomial.legendre.leggauss(nc)
    gx, gw = np.polynomial.legendre.leggauss(nx)
    ge, gwe = np.polynomial.legendre.leggauss(ne)
    Wt = kernel_table(rs, ts, float(a), R_B, rho, tc2, w, c, MU, gs, gws, gc, gwc, gx, gw, ge, gwe, nphi, mode, gains(a), HR_G)
    return ts, rs, solve_r(Wt, rs.size, ts.size, DT, 1.0, ctd, MU)


def predict(rho, a):
    t0 = time.time()
    ts, rs, f0 = mean_field(rho, a, 0)
    _, _, f1 = mean_field(rho, a, 1)
    _, _, ff = mean_field(rho, 0.0, 0)
    print(f"theory metric, uniform ball R {R_B}, a {a}; rho {rho:g}; ({time.time() - t0:.0f}s)")
    rq = [0.0, 1.0, 2.0, 3.0, 4.5]
    for lab, g in (("exact - flat (D)", f0 - ff), ("exact - proper-time kernel", f0 - f1)):
        print(f"  {lab} by t (rows) and r = " + " ".join(f"{r:.1f}" for r in rq))
        for tq in (1.0, 2.0, 3.0, 4.0, 5.0, 6.0):
            i = int(round(tq / (ts[1] - ts[0])))
            print(f"   t {tq:.1f}: " + " ".join(f"{g[i, int(round(r / (rs[1] - rs[0])))]:+.4f}" for r in rq))
    return ts, rs, f0, f1, ff


def task(rho, a, seed, X=7.0, T=6.0):
    nmax = 1.0 + 3.0 * a
    t, x, y, z = sprinkle4(rho * nmax, X, T, 81000 + seed)
    r = np.sqrt(x * x + y * y + z * z)
    U = np.random.default_rng(91000 + seed).uniform(size=t.size)
    w = minvar(8)
    ctd = delta(w, rho)
    c = 1.0 / (2.0 * math.pi * math.sqrt(6.0 * rho))
    tc2 = math.sqrt(24.0 * LAM_CUT / (math.pi * rho))
    out = {"rho": rho, "a": a, "seed": seed, "X": X, "T": T, "R": R_B, "lam_cut": LAM_CUT}
    for lab, m in (("flat", U < 1.0 / nmax), ("curved", U < (1.0 + 2.0 * a * gfun(r)) / nmax)):
        tt, xx, yy, zz = t[m], x[m], y[m], z[m]
        a0 = time.time()
        if lab == "flat":
            phi, _, _ = evolve_fast(tt, xx, yy, zz, math.sqrt(tc2), NMAX_ALL, w[None, :], np.array([ctd]),
                                    AMP * np.ones((1, tt.size)), np.array([True]), c, MU, np.zeros(tt.size, bool), 0.25, 1)
            phi = phi[0]
        else:
            phi = evolve_tilt(tt, xx, yy, zz, float(a), R_B, tc2, NMAX_ALL, w, ctd, AMP * np.ones(tt.size), c, MU, gains(a), HR_G)
        inner = (np.abs(xx) < X - tt - 0.3) & (np.abs(yy) < X - tt - 0.3) & (np.abs(zz) < X - tt - 0.3)
        rr = np.sqrt(xx * xx + yy * yy + zz * zz)
        out[lab] = {"sec": time.time() - a0, "n": int(tt.size),
                    "ev": np.c_[tt[inner], rr[inner], phi[inner] / AMP].astype(np.float32).tolist()}
    return out


def cached(rho, a, seed, X=7.0, T=6.0):
    key = hashlib.sha256(json.dumps([str(rho), str(a), seed, X, T, R_B, LAM_CUT, "c1b-v2-bent"]).encode()).hexdigest()[:20]
    fn = os.path.join("diag", "out", "c1b_tilt_cache", key + ".json")
    if os.path.exists(fn):
        return json.load(open(fn))
    o = task(rho, a, seed, X, T)
    os.makedirs(os.path.dirname(fn), exist_ok=True)
    json.dump(o, open(fn + ".tmp", "w"))
    os.replace(fn + ".tmp", fn)
    return o


def compare(R):
    rho, a, T = R[0]["rho"], R[0]["a"], R[0]["T"]
    ts, rs, f0 = mean_field(rho, a, 0, T=T)
    _, _, f1 = mean_field(rho, a, 1, T=T)
    _, _, ff = mean_field(rho, 0.0, 0, T=T)
    nt, nr = int(T / TB_W) + 1, 9

    def interp(f, ev):
        it = np.clip(ev[:, 0] / (ts[1] - ts[0]), 0, ts.size - 1.0001)
        ir = np.clip(ev[:, 1] / (rs[1] - rs[0]), 0, rs.size - 1.0001)
        i0, j0 = it.astype(int), ir.astype(int)
        a_, b_ = it - i0, ir - j0
        return ((f[i0, j0] * (1 - b_) + f[i0, j0 + 1] * b_) * (1 - a_)
                + (f[i0 + 1, j0] * (1 - b_) + f[i0 + 1, j0 + 1] * b_) * a_)

    Dm, P0, P1 = [], [], []
    for r in R:
        res = {}
        for lab in ("curved", "flat"):
            ev = np.array(r[lab]["ev"])
            ib = (ev[:, 0] / TB_W).astype(int)
            jb = (ev[:, 1] / RB_W).astype(int)
            m = jb < nr
            Nn = np.zeros((nt, nr))
            S = np.zeros((nt, nr))
            np.add.at(Nn, (ib[m], jb[m]), 1)
            np.add.at(S, (ib[m], jb[m]), ev[m, 2])
            p = {}
            for k, f in (("0", f0), ("1", f1), ("f", ff)):
                Pp = np.zeros((nt, nr))
                np.add.at(Pp, (ib[m], jb[m]), interp(f, ev[m]))
                p[k] = Pp
            res[lab] = (S, Nn, p)
        Sc, Nc, pc = res["curved"]
        Sf, Nf, pf = res["flat"]
        with np.errstate(invalid="ignore", divide="ignore"):
            Dm.append(Sc / Nc - Sf / Nf)
            P0.append(pc["0"] / Nc - pf["f"] / Nf)
            P1.append(pc["1"] / Nc - pf["f"] / Nf)
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        Dm, P0, P1 = np.array(Dm), np.nanmean(np.array(P0), 0), np.nanmean(np.array(P1), 0)
        n = len(R)
        mean = np.nanmean(Dm, 0)
        se = np.nanstd(Dm, 0, ddof=1) / math.sqrt(n)
    print(f"######## theory metric, ball R {R_B}, a {a}; rho {rho:g}; {n} pairs; D = <phi>_curved - <phi>_flat by (t, r)")
    for i in range(nt):
        tc_ = (i + 0.5) * TB_W
        if tc_ < 1.0:
            continue
        cells = [f"r {(j + 0.5) * RB_W:.2f}: {mean[i, j]:+.3f}+-{se[i, j]:.3f} (ex {P0[i, j]:+.3f}, pt {P1[i, j]:+.3f})"
                 for j in range(nr) if np.isfinite(mean[i, j]) and np.isfinite(se[i, j])]
        print(f"  t {tc_:.2f}: " + "; ".join(cells))
    sel = np.isfinite(mean) & np.isfinite(se) & (se > 0)
    sel[: int(1.0 / TB_W)] = False
    for lab, cols in (("inside the ball (r < 1.5)", slice(0, 3)), ("vacuum (r > 1.5)", slice(3, nr))):
        s2 = sel.copy()
        s2[:, :] = False
        s2[:, cols] = sel[:, cols]
        if s2.sum() == 0:
            continue
        z0 = np.abs(mean - P0)[s2] / np.hypot(se[s2], 0.005)
        z1 = np.abs(mean - P1)[s2] / np.hypot(se[s2], 0.005)
        print(f"  {lab}: {s2.sum()} bins (t >= 1): chi2/bin exact {np.mean(z0 ** 2):.2f} (max z {z0.max():.2f}); "
              f"proper-time kernel {np.mean(z1 ** 2):.2f} (max z {z1.max():.2f})")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--predict", action="store_true")
    ap.add_argument("--rho", type=float, default=60.0)
    ap.add_argument("--a", type=float, default=0.05)
    ap.add_argument("--seeds", type=str, default="0-7")
    ap.add_argument("--workers", type=int, default=0)
    ap.add_argument("--X", type=float, default=7.0)
    ap.add_argument("--T", type=float, default=6.0)
    a = ap.parse_args()
    gains(a.a)                     # the table is computed once here; the workers load it from the disk
    gains(0.0)
    if a.predict:
        predict(a.rho, a.a)
        return
    lo, hi = (int(v) for v in a.seeds.split("-"))
    t0 = time.time()
    res = parallel_map(cached, [(a.rho, a.a, s, a.X, a.T) for s in range(lo, hi + 1)], default_workers(a.workers), print,
                       lambda i, r, d, n: f"[C1b] {d}/{n}: rho {r['rho']:g} a {r['a']} seed {r['seed']} "
                                          f"{r['flat']['sec'] + r['curved']['sec']:.0f}s ({time.time() - t0:.0f}s)")
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        compare(res)
    print(buf.getvalue())
    fn = os.path.join("diag", "out", f"c1b_tilt_rho{a.rho:g}_a{a.a}_{lo}-{hi}.txt")
    open(fn, "w").write(buf.getvalue())
    print("written", fn)


if __name__ == "__main__":
    main()
