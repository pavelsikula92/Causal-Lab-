"""
F6.8f-dev, step 1b at finite kappa and q: the exact interval volume in the real medium (development, decides nothing)
====================================================================================================================

docs/derivations/f68f_dev.md, step 1b. Optical sector n = exp(b cos(k.w + theta)), 3+1, L = 1, k = kappa,
b = eps / kappa^2, cos(k, pair) = CO = p / kappa, sigma = q / kappa. No expansion in b, q or 1/kappa.

The medium is stratified along k, so the optical distance between two points depends only on their separation
along k (dz), on the horizontal distance (rho) and on the phase at the start:  T = Tau(phase, dz, rho).
For the two families that matter (rays from x, phase theta, and rays from y, phase kappa CO + theta) the rays are
integrated as graphs z(rho) with the exact Lagrangian n(z) sqrt(1 + z'^2) (RK4, the excess tau = T - rho carried
separately so that no O(1) number is subtracted), for a dense set of initial slopes; at every rho node the first
arrival Tau(dz, rho) is the lower envelope over the family (cubic Hermite in dz with dTau/dz = n z' / sqrt(1 + z'^2)).
The envelope is tabulated on a (rho, dz) grid around the rays and interpolated (7-point Lagrange in rho, cubic
Hermite in dz).

For w = (l, X, Y) (along the pair, in the plane of k, out of it):
    from x: dz = CO l + sa X,        rho^2 = l^2 + X^2 + Y^2 - dz^2,
    from y: dz = CO (l - 1) + sa X,  rho^2 = (1 - l)^2 + X^2 + Y^2 - dz^2,   sa = sqrt(1 - CO^2),
    h = dt - Tau_x - Tau_y,  dt = sqrt(T_xy^2 + sigma^2) (fixed sigma), T_xy = Tau_x(CO, sa),
    V = int dl dX dY n(w)^3 h_+,  A = V / (pi sigma^4 / 24).
Y enters only through rho and n(w) does not depend on Y: for every (l, X) the root Y* of h = 0 is found and
int dY h is done by Gauss on [0, Y*]; X on a fine grid over the support; l by Gauss.

Run: python diag/f68f_1b_kappa.py eps p n_theta [levels ...]   (output diag/out/f68f_1b_kappa.txt, appended)
"""
import math
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import f68f_1b as B  # noqa: E402

N_RK = 4000
N_U = 2001
N_L = 24
N_X = 401
N_Y = 8
N_DZ = 41


def family(bb, kap, ph0, umax, rho_max):
    """Rays z(rho) from the origin, z'(0) = u in [-umax, umax], in n = exp(bb cos(kap z + ph0)).
    Returns rho grid and arrays (N_RK + 1, N_U) of z, P = n z'/sqrt(1 + z'^2) (= dTau/dz), tau = T - rho, E = dz/du."""
    u = np.linspace(-umax, umax, N_U)
    y = np.zeros((5, N_U))     # z, z', tau, E, E'
    y[1] = u
    y[4] = 1.0
    h = rho_max / N_RK

    def f(y):
        a = kap * y[0] + ph0
        lnn_z = -bb * kap * np.sin(a)
        lnn_zz = -bb * kap * kap * np.cos(a)
        n = np.exp(bb * np.cos(a))
        g = 1.0 + y[1] ** 2
        sq = np.sqrt(g)
        # z'' = (1 + z'^2) dln n / dz ;  tau' = n sqrt(1 + z'^2) - 1 = (n - 1) sq + (sq - 1)
        dtau = (n - 1.0) * sq + y[1] ** 2 / (sq + 1.0)
        dE2 = g * lnn_zz * y[3] + 2.0 * y[1] * lnn_z * y[4]
        return np.stack([y[1], g * lnn_z, dtau, y[4], dE2])

    out = np.empty((N_RK + 1, 5, N_U))
    out[0] = y
    for i in range(N_RK):
        k1 = f(y)
        k2 = f(y + h / 2 * k1)
        k3 = f(y + h / 2 * k2)
        k4 = f(y + h * k3)
        y = y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        out[i + 1] = y
    rho = np.linspace(0.0, rho_max, N_RK + 1)
    z, zp, tau, E = out[:, 0], out[:, 1], out[:, 2], out[:, 3]
    n = np.exp(bb * np.cos(kap * z + ph0))
    P = n * zp / np.sqrt(1.0 + zp ** 2)
    return rho, z, P, tau, E


def env_1d(z, tau, P, E, x):
    """Lower envelope over the family at one rho: tau(x) and dtau/dz(x)."""
    G = np.full(x.shape, np.inf)
    D = np.zeros(x.shape)
    sgn = np.sign(E)
    cuts = np.nonzero(sgn[1:] != sgn[:-1])[0] + 1
    bounds = np.r_[0, cuts, len(z)]
    for a, b in zip(bounds[:-1], bounds[1:]):
        if b - a < 2:
            continue
        zs, ts, ps = z[a:b], tau[a:b], P[a:b]
        if zs[0] > zs[-1]:
            zs, ts, ps = zs[::-1], ts[::-1], ps[::-1]
        m = (x >= zs[0]) & (x <= zs[-1])
        if not np.any(m):
            continue
        xm = x[m]
        j = np.clip(np.searchsorted(zs, xm) - 1, 0, len(zs) - 2)
        x0, hh = zs[j], zs[j + 1] - zs[j]
        t = (xm - x0) / hh
        h00, h10, h01, h11 = 2 * t ** 3 - 3 * t ** 2 + 1, t ** 3 - 2 * t ** 2 + t, -2 * t ** 3 + 3 * t ** 2, t ** 3 - t ** 2
        val = h00 * ts[j] + h10 * hh * ps[j] + h01 * ts[j + 1] + h11 * hh * ps[j + 1]
        dval = ((6 * t ** 2 - 6 * t) * ts[j] + (3 * t ** 2 - 4 * t + 1) * hh * ps[j] + (-6 * t ** 2 + 6 * t) * ts[j + 1]
                + (3 * t ** 2 - 2 * t) * hh * ps[j + 1]) / hh
        better = val < G[m]
        idx = np.nonzero(m)[0]
        G[idx[better]] = val[better]
        D[idx[better]] = dval[better]
    return G, D


class Table:
    """tau(dz, rho) of the first arrival, tabulated at all rho nodes on a dz window [c(rho) - W(rho), c(rho) + W(rho)]."""

    def __init__(self, fam, center, width):
        rho, z, P, tau, E = fam
        self.rho = rho
        self.h = rho[1] - rho[0]
        self.c = center
        self.w = width
        t = np.linspace(-1.0, 1.0, N_DZ)
        self.t = t
        self.G = np.empty((len(rho), N_DZ))
        self.D = np.empty((len(rho), N_DZ))
        for i in range(len(rho)):
            x = center[i] + width[i] * t
            if rho[i] == 0.0:
                self.G[i] = np.inf
                self.D[i] = 0.0
                continue
            self.G[i], self.D[i] = env_1d(z[i], tau[i], P[i], E[i], x)

    def __call__(self, dz, rho):
        """Interpolated tau at (dz, rho) (arrays)."""
        fi = rho / self.h
        i0 = np.clip(np.floor(fi).astype(int) - 3, 1, len(self.rho) - 7)
        out = np.zeros(dz.shape)
        wsum = np.zeros(dz.shape)
        for a in range(7):
            ia = i0 + a
            w = np.ones(dz.shape)
            for b_ in range(7):
                if b_ != a:
                    w *= (fi - (i0 + b_)) / (a - b_)
            tt = (dz - self.c[ia]) / self.w[ia]
            jf = (tt + 1.0) / 2.0 * (N_DZ - 1)
            j = np.clip(np.floor(jf).astype(int), 0, N_DZ - 2)
            u = jf - j
            hh = 2.0 * self.w[ia] / (N_DZ - 1)
            g0, g1 = self.G[ia, j], self.G[ia, j + 1]
            d0, d1 = self.D[ia, j], self.D[ia, j + 1]
            val = ((2 * u ** 3 - 3 * u ** 2 + 1) * g0 + (u ** 3 - 2 * u ** 2 + u) * hh * d0 + (-2 * u ** 3 + 3 * u ** 2) * g1
                   + (u ** 3 - u ** 2) * hh * d1)
            val = np.where((tt < -1) | (tt > 1), np.inf, val)
            out += w * val
            wsum += w
        return out


def phase_A(eps, p, th, kap, q):
    bb = eps / kap ** 2
    CO = p / kap
    sa = math.sqrt(1.0 - CO * CO)
    sig = q / kap
    # the scaling-limit rays x -> y (all within q^2/2 of the minimum) place the windows
    r = B.branches(eps, p, np.array([th]), 1000)[0]
    vv = B.newton(r, eps, p, np.full(r.size, th), 4000)
    vv = np.unique(np.round(vv, 10))
    _, _, S = B.end_state(vv, eps, p, np.full(vv.size, th), 4000)
    near = vv[S - S.min() < q * q / 2 + 0.05]
    vmax = float(np.max(np.abs(near))) + 2.0 * math.sqrt(eps) + 4.0
    rho_max = 1.0 + 0.01
    fx = family(bb, kap, th, vmax / kap, rho_max)
    fy = family(bb, kap, kap * CO + th, vmax / kap, rho_max)
    paths = B.integrate(near, eps, p, np.full(near.size, th), 2000)[:, 0, :] / kap     # in-plane offsets X(l)
    lg = np.linspace(0, 1, 2001)
    rho = fx[0]
    Xr_lo = np.array([np.interp(np.clip(rho, 0, 1), lg, paths.min(1))]).ravel()
    Xr_hi = np.array([np.interp(np.clip(rho, 0, 1), lg, paths.max(1))]).ravel()
    halfw = 4.0 * sig * np.sqrt(np.clip(rho * (1.0 - rho), 0, None)) + 2.0 * sig + 1e-12
    # dz windows: from x, dz = CO l + sa X with l ~ rho; from y, dz = CO (l - 1) + sa X with 1 - l ~ rho
    cx = CO * rho + sa * 0.5 * (Xr_lo + Xr_hi)
    wx = sa * 0.5 * (Xr_hi - Xr_lo) + halfw
    Xr_lo_y = np.interp(np.clip(1.0 - rho, 0, 1), lg, paths.min(1))
    Xr_hi_y = np.interp(np.clip(1.0 - rho, 0, 1), lg, paths.max(1))
    cy = -CO * rho + sa * 0.5 * (Xr_lo_y + Xr_hi_y)
    wy = sa * 0.5 * (Xr_hi_y - Xr_lo_y) + halfw
    Tx = Table(fx, cx, wx)
    Ty = Table(fy, cy, wy)
    Txy = sa + float(Tx(np.array([CO]), np.array([sa]))[0])
    dt = math.sqrt(Txy * Txy + sig * sig)
    xl, wl = np.polynomial.legendre.leggauss(N_L)
    ls = 0.5 * (xl + 1.0)
    wl = 0.5 * wl
    yg, wyg = np.polynomial.legendre.leggauss(N_Y)
    yg, wyg = 0.5 * (yg + 1.0), 0.5 * wyg
    V = 0.0
    for l, w_l in zip(ls, wl):
        Xc_lo = np.interp(l, lg, paths.min(1))
        Xc_hi = np.interp(l, lg, paths.max(1))
        span = 3.0 * sig * math.sqrt(l * (1 - l)) + 1e-12
        X = np.linspace(Xc_lo - span, Xc_hi + span, N_X)

        def hfun(X, Y2):
            dzx = CO * l + sa * X
            rx = np.sqrt(np.maximum(l * l + X * X + Y2 - dzx ** 2, 0.0))
            dzy = CO * (l - 1.0) + sa * X
            ry = np.sqrt(np.maximum((1 - l) ** 2 + X * X + Y2 - dzy ** 2, 0.0))
            return (dt - rx - ry) - Tx(dzx, rx) - Ty(dzy, ry)

        h0 = hfun(X, np.zeros_like(X))
        if not (h0[0] < 0 and h0[-1] < 0):
            raise RuntimeError("X window")
        ok = h0 > 0
        if not np.any(ok):
            continue
        Xo = X[ok]
        # Y*^2 by Newton-secant on Y^2 (h is close to linear in Y^2)
        a = np.zeros_like(Xo)
        fa = h0[ok]
        b2 = np.full_like(Xo, 4.0 * sig * sig)
        fb = hfun(Xo, b2)
        for _ in range(60):
            grow = fb > 0
            if not np.any(grow):
                break
            b2 = np.where(grow, 2 * b2, b2)
            fb = np.where(grow, hfun(Xo, b2), fb)
        for _ in range(60):
            m = b2 - fb * (b2 - a) / (fb - fa)
            m = np.where((m > a) & (m < b2), m, 0.5 * (a + b2))
            fm = hfun(Xo, m)
            left = fm > 0
            a, fa = np.where(left, m, a), np.where(left, fm, fa)
            b2, fb = np.where(left, b2, m), np.where(left, fb, fm)
            if np.max(b2 - a) < 1e-15 * np.max(b2):
                break
        Ys = np.sqrt(0.5 * (a + b2))
        Iy = np.zeros_like(Xo)
        for yy, wy_ in zip(yg, wyg):
            Iy += wy_ * np.maximum(hfun(Xo, (yy * Ys) ** 2), 0.0)
        Iy *= 2.0 * Ys
        nn = np.exp(3 * bb * np.cos(kap * (CO * l + sa * Xo) + th))
        full = np.zeros_like(X)
        full[ok] = nn * Iy
        dX = X[1] - X[0]
        V += w_l * dX * full.sum()
    return V / (math.pi * sig ** 4 / 24.0)


def main():
    eps, p, n_th = float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3])
    levels = [float(x) for x in sys.argv[4:]] or [10.0, 20.0]
    q = 0.25
    t0 = time.time()
    ths = 2 * math.pi * np.arange(n_th) / n_th
    rows = []
    for lev in levels:
        kap = lev * math.sqrt(max(eps, 1.0))
        A = np.array([phase_A(eps, p, th, kap, q) for th in ths])
        rows.append((lev, kap, A))
        print(f"eps {eps} p {p} kappa {kap:.3f} (b {eps / kap ** 2:.2e}) q {q}: mean A {A.mean():.7f}", flush=True)
    lines = [f"eps {eps:g} p {p:g} q {q:g}, {n_th} phases:"]
    for lev, kap, A in rows:
        lines.append(f"  kappa = {lev:g} sqrt(eps) = {kap:.3f}, b = {eps / kap ** 2:.2e}: mean A {A.mean():.8f}, range "
                     f"{A.min():.6f} - {A.max():.6f}")
    if len(rows) > 1:
        d = rows[-1][2] - rows[-2][2]
        lines.append(f"  difference of the last two levels: mean {d.mean():+.2e}, largest per phase {np.max(np.abs(d)):.2e}")
    lines.append(f"  cost {time.time() - t0:.0f} s")
    txt = "\n".join(lines)
    print(txt)
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "f68f_1b_kappa.txt"), "a") as fh:
        fh.write(txt + "\n")


if __name__ == "__main__":
    main()
