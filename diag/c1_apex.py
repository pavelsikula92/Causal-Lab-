"""
C1, development: the exact mean of the rule on the conformally flat planar well (the (A3) of C1; 14 Nov; decides nothing)
=======================================================================================================================
diag/c1_conf.py sprinkles g = Omega(x)^2 eta (flat cones, coordinate intensity lambda = rho Omega^4). As in C0
(docs/derivations/c0_mean_exactness.md), the pyramid is past-complete, the interval I(u, v) lies outside the past of u,
and Mecke's formula holds for an inhomogeneous Poisson process. Hence the mean field f(t, x) = E phi obeys exactly

    f(t, x) (1 - delta mu) = A - c mu Int d^4u  lambda(u)  K(Lambda(u, v))  f(t_u, x_u),
    K(L) = sum_n w_n Pois(n; L),   Lambda(u, v) = Int_{I(u,v)} lambda = rho V_flat(u, v) <Omega^4>_I(u,v),

with the pairs of the enumeration (flat sigma <= tc^2). Coordinates of u relative to v: dt > 0, dx (along the well),
the transverse radius; with sigma = dt^2 - dx^2 - r_perp^2 the measure is d^4u = pi d(dt) d(dx) d(sigma).

The mean of Omega^4 over the diamond. The x coordinate of a point of the diamond is the mean x_mid plus a linear
functional of a uniform point of the unit diamond: X = dx u + sqrt(sigma + dx^2) R(u) z, with u the rest-frame time
in [-1/2, 1/2], R = 1/2 - |u|, z in [-1, 1]. Its joint density is proportional to R^3 (1 - z^2). This is exact for
any pair (the x axis seen from the pair's rest frame: alpha = dx/tau, beta = sqrt(1 + alpha^2)) and is integrated by a
product Gauss quadrature (Gauss-Legendre in u on each half, Gauss-Gegenbauer in z).

Numerics: v on an x grid (step DX); u at fine Gauss nodes in dx within |dx| < dt, interpolated linearly onto the x
grid; sigma by Gauss-Legendre on [0, min(tc^2, dt^2 - dx^2)]; trapezoid in dt. The same grid computes the flat case
(a = 0), so their difference D is free of most of the discretization error. Checked against the flat Volterra of
diag/c0_apex (check_flat).
"""
from __future__ import annotations

import math

import numpy as np

from phases.F3_structure import field1d as F

_nb = F._nb


def diamond_nodes(nu=6, nz=5):
    """(u, R(u) z, weight) for the unit diamond, density proportional to R^3 (1 - z^2)."""
    gu, wu = np.polynomial.legendre.leggauss(nu)
    us, ws = [], []
    for sgn in (-1.0, 1.0):
        u = sgn * (0.25 + 0.25 * gu)            # map [-1, 1] to [0, 1/2] (or [-1/2, 0])
        w = 0.25 * wu * (0.5 - np.abs(u)) ** 3
        us.append(u)
        ws.append(w)
    us, ws = np.concatenate(us), np.concatenate(ws)
    # Gauss-Gegenbauer (weight 1 - z^2): Gauss-Jacobi alpha = beta = 1
    from scipy.special import roots_jacobi
    gz, wz = roots_jacobi(nz, 1.0, 1.0)
    U, Z, W = [], [], []
    for u, w in zip(us, ws):
        R = 0.5 - abs(u)
        for z, v in zip(gz, wz):
            U.append(u)
            Z.append(R * z)
            W.append(w * v)
    W = np.array(W)
    return np.array(U), np.array(Z), W / W.sum()


if _nb is not None:
    @_nb.njit(cache=True)
    def _om4(x, a, L):
        o = 1.0 - a * math.exp(-x * x / (2.0 * L * L))
        return o * o * o * o

    @_nb.njit(cache=True)
    def kernel_table(xs, dts, J, DX, a, L, rho, tc2, w, c, mu, nU, nZ, nW, gs, gws, gx, gwx, mode=0):
        """W[ix, it, J + j]: the weight of f(t - dt_it, x_ix + j DX) in the mean of the sum at (t, x_ix).
        mode 0: the exact mean of the rule; mode 1: the conformally coupled continuum Green function
        Omega_u^-1 Omega_v^-1 G_eta smeared by the flat kernel (Lambda flat, weight rho Omega_u^3 / Omega_v)."""
        nx = xs.size
        nt = dts.size
        out = np.zeros((nx, nt, 2 * J + 1))
        nn = w.size
        for ix in range(nx):
            xv = xs[ix]
            for it in range(1, nt):
                dt = dts[it]
                for q in range(gx.size):
                    ddx = dt * gx[q]                         # fine node in (-dt, dt)
                    wq = dt * gwx[q]
                    smax = dt * dt - ddx * ddx
                    if smax > tc2:
                        smax = tc2
                    if smax <= 0.0:
                        continue
                    xu = xv + ddx
                    lam_u = rho * _om4(xu, a, L)
                    if mode == 1:
                        ou = 1.0 - a * math.exp(-xu * xu / (2.0 * L * L))
                        ov = 1.0 - a * math.exp(-xv * xv / (2.0 * L * L))
                        lam_u = lam_u / (ou * ov)
                    acc = 0.0
                    for p in range(gs.size):
                        sig = 0.5 * smax * (gs[p] + 1.0)
                        wsig = 0.5 * smax * gws[p]
                        xm = xv + 0.5 * ddx
                        sp = math.sqrt(sig + ddx * ddx)
                        e = 0.0
                        if mode == 1:
                            e = 1.0
                        else:
                            for k in range(nW.size):
                                e += nW[k] * _om4(xm + ddx * nU[k] + sp * nZ[k], a, L)
                        lam = rho * math.pi / 24.0 * sig * sig * e
                        # K = sum_n w_n Pois(n; lam)
                        pk = math.exp(-lam)
                        kk = 0.0
                        for n in range(nn):
                            if n > 0:
                                pk *= lam / n
                            kk += w[n] * pk
                        acc += wsig * kk
                    val = c * mu * math.pi * lam_u * acc * wq
                    # spread onto the grid of x_u = xv + ddx
                    fj = ddx / DX
                    j0 = int(math.floor(fj))
                    fr = fj - j0
                    if -J <= j0 <= J:
                        out[ix, it, J + j0] += val * (1.0 - fr)
                    if -J <= j0 + 1 <= J:
                        out[ix, it, J + j0 + 1] += val * fr
        return out

    @_nb.njit(cache=True)
    def solve(Wt, nx, nt, J, dtau, A, ct, mu, far):
        """f[it, ix]; trapezoid in dt (weight 1/2 at the last term dt = t, i.e. u at t = 0); far[it] for x off the grid."""
        f = np.zeros((nt, nx))
        for ix in range(nx):
            f[0, ix] = A / (1.0 - ct * mu)
        for it in range(1, nt):
            for ix in range(nx):
                s = 0.0
                for k in range(1, it + 1):
                    wk = 0.5 if k == it else 1.0
                    src = it - k
                    for j in range(-J, J + 1):
                        wv = Wt[ix, k, J + j]
                        if wv == 0.0:
                            continue
                        iy = ix + j
                        if 0 <= iy < nx:
                            fv = f[src, iy]
                        else:
                            fv = far[src]
                        s += wk * wv * fv
                f[it, ix] = (A - dtau * s) / (1.0 - ct * mu)
        return f


def mean_field(rho, a, L, w, ct, T=6.0, Xg=9.0, DX=0.1, DT=0.05, ns=10, nfine=48, mu=1.0, A=1.0, lam_cut=22.0, mode=0):
    """f(t, x)/A on the grid (ts, xs) for the well (a, L); a = 0 gives the flat field on the same grid."""
    xs = np.arange(-Xg, Xg + 1e-9, DX)
    ts = np.arange(0.0, T + 1e-9, DT)
    J = int(math.ceil(T / DX)) + 1
    tc2 = math.sqrt(24.0 * lam_cut / (math.pi * rho))
    c = 1.0 / (2.0 * math.pi * math.sqrt(6.0 * rho))
    U, Z, Wd = diamond_nodes()
    gs, gws = np.polynomial.legendre.leggauss(ns)
    gx, gwx = np.polynomial.legendre.leggauss(nfine)
    W = kernel_table(xs, ts, J, DX, a, L, rho, tc2, np.asarray(w, float), c, mu, U, Z, Wd, gs, gws, gx, gwx, mode)
    # far field: the flat field on the same discretization (x off the grid is flat for t <= T, see the module text)
    if a != 0.0:
        Wf = kernel_table(np.array([0.0]), ts, J, DX, 0.0, L, rho, tc2, np.asarray(w, float), c, mu, U, Z, Wd, gs,
                          gws, gx, gwx)
        # the flat field is homogeneous: with one grid point every off-grid neighbour has the same value
        far = _flat_selfconsistent(Wf, ts.size, J, DT, A, ct, mu)
    else:
        far = _flat_selfconsistent(W[:1], ts.size, J, DT, A, ct, mu)
    f = solve(W, xs.size, ts.size, J, DT, A, ct, mu, far)
    return ts, xs, f / A


def _flat_selfconsistent(W1, nt, J, DT, A, ct, mu):
    """The homogeneous flat field from a one-point kernel table: all offsets see the same value."""
    tot = W1[0].sum(axis=1)                          # [nt]
    f = np.zeros(nt)
    f[0] = A / (1.0 - ct * mu)
    for it in range(1, nt):
        s = 0.0
        for k in range(1, it + 1):
            wk = 0.5 if k == it else 1.0
            s += wk * tot[k] * f[it - k]
        f[it] = (A - DT * s) / (1.0 - ct * mu)
    return f


def check_flat(rho=30.0, **kw):
    """The flat case of this solver against the 1D flat Volterra of diag/c0_apex (same kernel MV8+ct)."""
    from diag.c0_apex import F_k, volterra, delta
    from diag.c0_minvar import minvar
    w = minvar(8)
    ct = delta(w, rho)
    ts, xs, f = mean_field(rho, 0.0, 1.5, w, ct, **kw)
    dt = 0.002
    tg = np.arange(0, ts[-1] + dt, dt)
    F0 = F_k(w, rho, 0, tg)
    g = volterra(F0, np.ones_like(tg), dt, 1.0 / (1.0 - ct))
    ref = np.interp(ts, tg, g)
    mid = xs.size // 2
    print(f"flat check rho {rho:g}: max |f_grid - f_1D| over t {np.max(np.abs(f[:, mid] - ref)):.4f}; at t "
          + " ".join(f"{t:.1f}:{f[i, mid] - ref[i]:+.4f}" for i, t in enumerate(ts) if abs(t * 2 - round(t * 2)) < 1e-9 and t > 0))
    return ts, f[:, mid], ref


def grid_at(ts, xs, f, ev):
    ev = np.asarray(ev, float)
    tt, xx = ev[:, 0], ev[:, 1]
    ft = np.clip(tt / (ts[1] - ts[0]), 0, ts.size - 1.000001)
    j0 = np.floor(ft).astype(int)
    at = ft - j0
    fx = np.clip((xx - xs[0]) / (xs[1] - xs[0]), 0, xs.size - 1.000001)
    i0 = np.floor(fx).astype(int)
    ax = fx - i0
    v0 = f[j0, i0] * (1 - ax) + f[j0, i0 + 1] * ax
    v1 = f[j0 + 1, i0] * (1 - ax) + f[j0 + 1, i0 + 1] * ax
    return v0 * (1 - at) + v1 * at


def compare(R, verbose=True):
    """The measured paired difference D against the exact mean (this module) and the continuum xi scan."""
    import math as m_
    from diag import c1_conf as C
    from diag.c0_apex import delta
    from diag.c0_minvar import minvar
    rho = R[0]["rho"]
    X, T = R[0]["X"], R[0]["T"]
    w = minvar(8)
    ct = delta(w, rho)
    ts, xs, fc = mean_field(rho, C.A_WELL, C.L_WELL, w, ct, T=T)
    _, _, ff = mean_field(rho, 0.0, C.L_WELL, w, ct, T=T, Xg=1.0)
    nt, nx = np.array(R[0]["flat"]["S"]).shape
    Dm, Dp = [], []
    for r in R:
        Sf, Nf = np.array(r["flat"]["S"]), np.array(r["flat"]["N"])
        Sc, Nc = np.array(r["curved"]["S"]), np.array(r["curved"]["N"])
        evf, evc = np.array(r["flat"]["ev"]), np.array(r["curved"]["ev"])
        pc = np.zeros((nt, nx))
        pf = np.zeros((nt, nx))
        bc = (np.floor(evc[:, 0] / C.TB_W).astype(int), np.floor((evc[:, 1] + X) / C.XB_W).astype(int))
        bf = (np.floor(evf[:, 0] / C.TB_W).astype(int), np.floor((evf[:, 1] + X) / C.XB_W).astype(int))
        np.add.at(pc, bc, grid_at(ts, xs, fc, evc))
        np.add.at(pf, bf, np.interp(evf[:, 0], ts, ff[:, ff.shape[1] // 2]))
        with np.errstate(invalid="ignore", divide="ignore"):
            Dm.append(Sc / Nc - Sf / Nf)
            Dp.append(pc / Nc - pf / Nf)
    Dm, Dp = np.array(Dm), np.array(Dp)
    n = len(R)
    mean = np.nanmean(Dm, 0)
    se = np.nanstd(Dm, 0, ddof=1) / m_.sqrt(n)
    pred = np.nanmean(Dp, 0)
    xc = (np.arange(nx) + 0.5) * C.XB_W - X
    tcs = (np.arange(nt) + 0.5) * C.TB_W
    out = {}
    if verbose:
        print(f"######## rho {rho:g}, {n} pairs: D at the centre |x| < 0.5, measured / exact mean (A3 of C1)")
    for i in range(nt):
        cen = [j for j in range(nx) if abs(xc[j]) < 0.5 and np.isfinite(mean[i, j])]
        if cen and tcs[i] >= 0.5 and verbose:
            e = m_.sqrt(np.sum(se[i, cen] ** 2)) / len(cen)
            print(f"  t {tcs[i]:.2f}: {np.mean(mean[i, cen]):+.3f} +- {e:.3f}   exact mean {np.mean(pred[i, cen]):+.3f}")
    for lo, hi in ((1.0, 3.5), (3.5, 6.0)):
        c2, nn, mx = 0.0, 0, 0.0
        for i in range(nt):
            if not (lo <= tcs[i] < hi):
                continue
            for j in range(nx):
                if abs(xc[j]) <= 2.5 and np.isfinite(mean[i, j]) and se[i, j] > 0:
                    c2 += (mean[i, j] - pred[i, j]) ** 2 / (se[i, j] ** 2 + 0.005 ** 2)
                    nn += 1
                    mx = max(mx, abs(mean[i, j] - pred[i, j]))
        out[(lo, hi)] = (c2 / nn, mx, nn)
        if verbose:
            print(f"  t {lo}-{hi}: chi2/bin {c2 / nn:.2f} ({nn} bins, SE floor 0.005), max |meas - exact mean| {mx:.3f}")
    return out
