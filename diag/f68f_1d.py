"""
F6.8f-dev, step 1d (development, decides nothing): the needle volume at finite kappa and b in the real medium
===========================================================================================================

docs/derivations/f68f_dev.md, step 1d (the author's 1d: an independent brute-force comparison at moderate eps and
finite L / sigma). Optical sector n = exp(b cos(k.w + theta)), 3+1, units L = 1 (the pair x = 0, y = (dt, 1, 0, 0)),
k = kappa (CO, sa, 0), CO = p / kappa, b = eps / kappa^2, sigma = q / kappa. No expansion in b, q or 1 / kappa.

The medium is stratified along k, so the optical distance between two points depends only on their separation along k
(dz), on the horizontal distance (rho) and on the phase at the start. For the families from x (phase theta) and from y
(phase kappa CO + theta) the rays are integrated as graphs z(rho) with the exact Lagrangian n(z) sqrt(1 + z'^2) (RK4;
the excess tau = T - rho carried separately), for a dense set of initial slopes, with the slope derivatives; at every
rho node the first arrival tau(dz) is the lower envelope over the family (quintic Hermite in z: tau, dtau/dz = P,
d2tau/dz2 = (dP/du) / (dz/du)), tabulated on a window around the ray and interpolated (7-point Lagrange in rho, cubic
Hermite in dz on a fine table). The window at every rho node is the union of a nominal band around the ray and the
exact dz range of the points queried there (the X grids of the l nodes with Y^2 <= 8 sigma^2 l (1 - l), and T_xy):
with rho_x = |sa l - CO X| at Y = 0, a node rho is reached from l = rho / sa, and near the ends from the X span itself,
so a band built on l = rho does not cover large CO = p / kappa.

For w = (l, X, Y) (along the pair, in the plane of k, out of it):
    from x: dz = CO l + sa X,        rho^2 = l^2 + X^2 + Y^2 - dz^2,
    from y: dz = CO (l - 1) + sa X,  rho^2 = (1 - l)^2 + X^2 + Y^2 - dz^2,   sa = sqrt(1 - CO^2),
    h = dt - T_x - T_y,  dt = sqrt(T_xy^2 + sigma^2) (fixed sigma), T_xy = T_x(CO, sa),
    V = int dl dX dY n(w)^3 h_+,  A = V / (pi sigma^4 / 24).
For every (l, X) the root Y* of h = 0 is found (h is close to linear in Y^2) and int dY is done by Gauss on [0, Y*];
X on a fine grid over the support; l by Gauss.

Comparison: the scaling limit (kappa -> infinity) of A at the same (eps, p, theta, q) is the q^2 series of
diag.f68f_1c; the finite-kappa corrections are O(b) = O(eps / kappa^2), so two levels kappa = 10 sqrt(eps),
20 sqrt(eps) (b = 0.01, 0.0025) are extrapolated (A_inf = (4 A_20 - A_10) / 3).

Run from the repository root (one job = one (eps, p, theta, kappa level, q), cached in diag/out/f68f_1d_cache/):
    python -m diag.f68f_1d --points 2:1:1.0,4:2:2.0 --levels 10,20 --q 0.25 --workers 2
Output diag/out/f68f_1d_<UTC time>.json and .txt.
"""
from __future__ import annotations

import argparse
import math
import os
import time

import numpy as np

from core.stats.parallel import default_workers, parallel_map
from diag import f68f_1b as B
from diag import f68f_common as C

N_RK = 4000
N_U = 2001
N_L = 24
N_X = 401
N_Y = 8
N_DZ = 81


def family(bb, kap, ph0, umax, rho_max):
    """Rays z(rho) from the origin with z'(0) = u in [-umax, umax], n = exp(bb cos(kap z + ph0)). Returns the rho grid
    and arrays (N_RK + 1, N_U) of z, tau = T - rho, P = dT/dz at fixed rho, Q = d2T/dz2 along the family, E = dz/du."""
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
        dtau = (n - 1.0) * sq + y[1] ** 2 / (sq + 1.0)        # n sqrt(1 + z'^2) - 1 without cancellation
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
    z, zp, tau, E, dE = out[:, 0], out[:, 1], out[:, 2], out[:, 3], out[:, 4]
    a = kap * z + ph0
    n = np.exp(bb * np.cos(a))
    n_z = -bb * kap * np.sin(a) * n
    g = 1.0 + zp ** 2
    P = n * zp / np.sqrt(g)
    with np.errstate(divide="ignore", invalid="ignore"):
        Q = n_z * zp / np.sqrt(g) + n * (dE / E) / g ** 1.5
    return rho, z, tau, P, Q, E


def env_1d(z, tau, P, Q, E, x):
    """Lower envelope over the family at one rho: tau(x) and dtau/dz(x) (quintic Hermite on every monotone piece)."""
    G = np.full(x.shape, np.inf)
    D = np.zeros(x.shape)
    sgn = np.sign(E)
    cuts = np.nonzero(sgn[1:] != sgn[:-1])[0] + 1
    bounds = np.r_[0, cuts, len(z)]
    for a, b in zip(bounds[:-1], bounds[1:]):
        if b - a < 2:
            continue
        zs, ts, ps, qs = z[a:b], tau[a:b], P[a:b], Q[a:b]
        if zs[0] > zs[-1]:
            zs, ts, ps, qs = zs[::-1], ts[::-1], ps[::-1], qs[::-1]
        m = (x >= zs[0]) & (x <= zs[-1])
        if not np.any(m):
            continue
        xm = x[m]
        j = np.clip(np.searchsorted(zs, xm) - 1, 0, len(zs) - 2)
        x0, hh = zs[j], zs[j + 1] - zs[j]
        t = (xm - x0) / hh
        t2, t3, t4, t5 = t * t, t ** 3, t ** 4, t ** 5
        h0, h1 = 1 - 10 * t3 + 15 * t4 - 6 * t5, t - 6 * t3 + 8 * t4 - 3 * t5
        h2 = 0.5 * t2 - 1.5 * t3 + 1.5 * t4 - 0.5 * t5
        h3, h4, h5 = 10 * t3 - 15 * t4 + 6 * t5, -4 * t3 + 7 * t4 - 3 * t5, 0.5 * t3 - t4 + 0.5 * t5
        d0 = (-30 * t2 + 60 * t3 - 30 * t4)
        d1 = 1 - 18 * t2 + 32 * t3 - 15 * t4
        d2 = t - 4.5 * t2 + 6 * t3 - 2.5 * t4
        d3 = 30 * t2 - 60 * t3 + 30 * t4
        d4 = -12 * t2 + 28 * t3 - 15 * t4
        d5 = 1.5 * t2 - 4 * t3 + 2.5 * t4
        val = (h0 * ts[j] + h1 * hh * ps[j] + h2 * hh * hh * qs[j] + h3 * ts[j + 1] + h4 * hh * ps[j + 1]
               + h5 * hh * hh * qs[j + 1])
        dval = (d0 * ts[j] + d1 * hh * ps[j] + d2 * hh * hh * qs[j] + d3 * ts[j + 1] + d4 * hh * ps[j + 1]
                + d5 * hh * hh * qs[j + 1]) / hh
        better = val < G[m]
        idx = np.nonzero(m)[0]
        G[idx[better]] = val[better]
        D[idx[better]] = dval[better]
    return G, D


class Table:
    """tau(dz, rho) of the first arrival at every rho node, on a dz window [c - w, c + w]."""

    def __init__(self, fam, center, width):
        rho, z, tau, P, Q, E = fam
        self.n = len(rho)
        self.h = rho[1] - rho[0]
        self.c, self.w = center, width
        t = np.linspace(-1.0, 1.0, N_DZ)
        self.G = np.full((self.n, N_DZ), np.inf)
        self.D = np.zeros((self.n, N_DZ))
        for i in range(1, self.n):
            self.G[i], self.D[i] = env_1d(z[i], tau[i], P[i], Q[i], E[i], center[i] + width[i] * t)

    def __call__(self, dz, rho):
        fi = rho / self.h
        i0 = np.clip(np.floor(fi).astype(int) - 3, 1, self.n - 7)
        out = np.zeros(dz.shape)
        for a in range(7):
            ia = i0 + a
            w = np.ones(dz.shape)
            for b_ in range(7):
                if b_ != a:
                    w *= (fi - (i0 + b_)) / (a - b_)
            tt = (dz - self.c[ia]) / self.w[ia]
            jf = (tt + 1.0) / 2.0 * (N_DZ - 1)
            j = np.clip(np.floor(jf).astype(int), 0, N_DZ - 2)
            uu = jf - j
            hh = 2.0 * self.w[ia] / (N_DZ - 1)
            g0, g1 = self.G[ia, j], self.G[ia, j + 1]
            e0, e1 = self.D[ia, j], self.D[ia, j + 1]
            val = ((2 * uu ** 3 - 3 * uu ** 2 + 1) * g0 + (uu ** 3 - 2 * uu ** 2 + uu) * hh * e0
                   + (-2 * uu ** 3 + 3 * uu ** 2) * g1 + (uu ** 3 - uu ** 2) * hh * e1)
            out += w * np.where((tt < -1) | (tt > 1), np.inf, val)
        return out


Y2_FAC = 8.0          # Y^2 range of the windows in units of sigma^2 l (1 - l) (the root bracket doubles from 1)
COVER_MARGIN = 0.05   # relative margin of the query-based windows


def _queries(ls, path, lg, sig, CO, sa, from_y):
    """(dz array, rho_lo array, rho_hi array) per l node: the X grid of volume() at Y = 0 and the Y range."""
    out = []
    for l in ls:
        Xc = float(np.interp(l, lg, path))
        span = 3.0 * sig * math.sqrt(l * (1 - l))
        X = np.linspace(Xc - span, Xc + span, N_X)
        if from_y:
            dz, r0 = CO * (l - 1.0) + sa * X, np.abs(sa * (1.0 - l) + CO * X)
        else:
            dz, r0 = CO * l + sa * X, np.abs(sa * l - CO * X)
        out.append((dz, r0, np.sqrt(r0 * r0 + Y2_FAC * sig * sig * l * (1 - l))))
    return out


def _cover(c0, w0, nodes, rho_max, queries):
    """Union of the window [c0 - w0, c0 + w0] with the dz range of the queries at every node of their 7-point
    Lagrange stencils in rho (Table.__call__), widened by COVER_MARGIN."""
    n = len(nodes)
    h = rho_max / (n - 1)
    lo, hi = c0 - w0, c0 + w0
    for dz, r0, r1 in queries:
        dz, r0, r1 = np.atleast_1d(dz), np.atleast_1d(r0), np.atleast_1d(r1)
        i_a = np.clip(np.floor(r0 / h).astype(int) - 3, 1, n - 7)
        i_b = np.clip(np.floor(r1 / h).astype(int) - 3, 1, n - 7) + 6
        for a, b, z in zip(i_a, i_b, dz):
            lo[a:b + 1] = np.minimum(lo[a:b + 1], z)
            hi[a:b + 1] = np.maximum(hi[a:b + 1], z)
    c, w = 0.5 * (lo + hi), 0.5 * (hi - lo)
    return c, w * (1.0 + COVER_MARGIN) + 1e-14


def volume(eps, p, th, kap, q, slices=False):
    """A = V / V0 of the needle at finite kappa (fixed sigma), and diagnostics."""
    bb = eps / kap ** 2
    CO = p / kap
    sa = math.sqrt(1.0 - CO * CO)
    sig = q / kap
    r = B.branches(eps, p, np.array([th]), 1000)[0]
    vv = B.newton(r, eps, p, np.full(r.size, th), 4000)
    vv = np.unique(np.round(vv, 10))
    _, _, S = B.end_state(vv, eps, p, np.full(vv.size, th), 4000)
    v = vv[np.argmin(S)]
    path = B.integrate(np.array([v]), eps, p, np.array([th]), 2000)[:, 0, 0] / kap      # in-plane offset X(l) ~ psi / kappa
    lg = np.linspace(0, 1, 2001)
    # windows: half-width ~ 4 sigma sqrt(l (1 - l)) + a margin; the families must cover them at every node
    rho_max = 1.01
    nodes = np.linspace(0, rho_max, N_RK + 1)
    Wrho = 6.0 * sig * np.sqrt(np.clip(nodes * (1 - nodes), 1e-6, None)) + 0.5 * sig * np.clip(nodes, 0, 1)
    Xc_x = np.interp(np.clip(nodes, 0, 1), lg, path)
    Xc_y = np.interp(np.clip(1 - nodes, 0, 1), lg, path)
    xl, wl = np.polynomial.legendre.leggauss(N_L)
    ls, wl = 0.5 * (xl + 1.0), 0.5 * wl
    # the windows must contain every (dz, rho) that is queried. A query point w = (l, X, Y) has
    # rho_x^2 = (sa l - CO X)^2 + Y^2 and rho_y^2 = (sa (1 - l) + CO X)^2 + Y^2, so a table node rho is reached from
    # l = rho / sa and, near the ends, from the X span itself; the windows are widened to the union of the nominal
    # window above and the exact query range (X grids of the l nodes, Y^2 <= Y2_FAC sigma^2 l (1 - l), T_xy)
    cx, wx = _cover(CO * nodes + sa * Xc_x, sa * Wrho + 1e-14, nodes, rho_max,
                    _queries(ls, path, lg, sig, CO, sa, False) + [(CO, sa, sa)])
    cy, wy = _cover(-CO * nodes + sa * Xc_y, sa * Wrho + 1e-14, nodes, rho_max,
                    _queries(ls, path, lg, sig, CO, sa, True))
    reach = max(np.max((np.abs(cx[1:]) + wx[1:]) / nodes[1:]), np.max((np.abs(cy[1:]) + wy[1:]) / nodes[1:]))
    umax = reach + 2.0 * math.sqrt(eps) / kap + 4.0 / kap
    fx = family(bb, kap, th, umax, rho_max)
    fy = family(bb, kap, kap * CO + th, umax, rho_max)
    Tx = Table(fx, cx, wx)
    Ty = Table(fy, cy, wy)
    Txy = sa + float(Tx(np.array([CO]), np.array([sa]))[0])
    dt = math.sqrt(Txy * Txy + sig * sig)
    yg, wyg = np.polynomial.legendre.leggauss(N_Y)
    yg, wyg = 0.5 * (yg + 1.0), 0.5 * wyg
    V = 0.0
    sl = []
    for l, w_l in zip(ls, wl):
        Xc = float(np.interp(l, lg, path))
        span = 3.0 * sig * math.sqrt(l * (1 - l))
        X = np.linspace(Xc - span, Xc + span, N_X)

        def hfun(X, Y2):
            dzx = CO * l + sa * X
            rx = np.sqrt(np.maximum(l * l + X * X + Y2 - dzx ** 2, 0.0))
            dzy = CO * (l - 1.0) + sa * X
            ry = np.sqrt(np.maximum((1 - l) ** 2 + X * X + Y2 - dzy ** 2, 0.0))
            return (dt - rx - ry) - Tx(dzx, rx) - Ty(dzy, ry)

        h0 = hfun(X, np.zeros_like(X))
        if not (np.isfinite(h0[0]) and np.isfinite(h0[-1])):
            raise RuntimeError("X window outside the tables")
        if not (h0[0] < 0 and h0[-1] < 0):
            raise RuntimeError("X window too narrow")
        ok = h0 > 0
        if not np.any(ok):
            sl.append((l, 0.0))
            continue
        Xo = X[ok]
        a = np.zeros_like(Xo)
        fa = h0[ok]
        b2 = np.full_like(Xo, sig * sig * l * (1 - l))
        fb = hfun(Xo, b2)
        for _ in range(60):
            grow = fb > 0
            if not np.any(grow):
                break
            b2 = np.where(grow, 2 * b2, b2)
            fb = np.where(grow, hfun(Xo, b2), fb)
        if not np.all(np.isfinite(fb)):
            raise RuntimeError("Y root outside the tables")
        # Illinois regula falsi: plain regula falsi stalls with one end fixed (a = 0), which gave Ys^2 = root / 2
        side = np.zeros(Xo.shape, dtype=int)
        for _ in range(200):
            m = b2 - fb * (b2 - a) / (fb - fa)
            m = np.where((m > a) & (m < b2), m, 0.5 * (a + b2))
            fm = hfun(Xo, m)
            left = fm > 0
            fb = np.where(left & (side == 1), 0.5 * fb, fb)
            fa = np.where(~left & (side == -1), 0.5 * fa, fa)
            a, fa = np.where(left, m, a), np.where(left, fm, fa)
            b2, fb = np.where(left, b2, m), np.where(left, fb, fm)
            side = np.where(left, 1, -1)
            if np.all(b2 - a <= 1e-13 * b2):
                break
        if not np.all(b2 - a <= 1e-10 * b2):
            raise RuntimeError("Y root not converged")
        Ys = np.sqrt(0.5 * (a + b2))
        Iy = np.zeros_like(Xo)
        for yy, wy_ in zip(yg, wyg):
            Iy += wy_ * np.maximum(hfun(Xo, (yy * Ys) ** 2), 0.0)
        Iy *= 2.0 * Ys
        nn = np.exp(3 * bb * np.cos(kap * (CO * l + sa * Xo) + th))
        full = np.zeros_like(X)
        full[ok] = nn * Iy
        slice_ = (X[1] - X[0]) * (full.sum() - 0.5 * (full[0] + full[-1]))
        sl.append((l, slice_))
        V += w_l * slice_
    if slices:
        return V / (math.pi * sig ** 4 / 24.0), [(l, v / (math.pi * sig ** 4 * l * (1 - l) / 4)) for l, v in sl]
    return V / (math.pi * sig ** 4 / 24.0)


NAME = "f68f_1d"
CODE = C.code_md5([os.path.abspath(__file__), B.__file__])


def job(eps, p, th, kap, q):
    return C.cached(NAME, CODE, volume, (eps, p, th, kap, q))


def extrapolate(levs, A):
    """Richardson in 1/kappa^2 for levels doubling in kappa: two levels (4 A2 - A1) / 3; three levels
    (16 R2 - R1) / 15 with R_i the two-level values (removes the O(1/kappa^2) and O(1/kappa^4) terms)."""
    if len(levs) == 2 and abs(levs[1] - 2 * levs[0]) < 1e-12:
        return (4 * A[1] - A[0]) / 3
    if len(levs) == 3 and abs(levs[1] - 2 * levs[0]) < 1e-12 and abs(levs[2] - 2 * levs[1]) < 1e-12:
        r1, r2 = (4 * A[1] - A[0]) / 3, (4 * A[2] - A[1]) / 3
        return (16 * r2 - r1) / 15
    return float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--points", type=str, default="2:1:1.0", help="eps:p:theta,...")
    ap.add_argument("--levels", type=str, default="10,20", help="kappa = level sqrt(eps)")
    ap.add_argument("--q", type=float, default=0.25)
    ap.add_argument("--workers", type=int, default=0)
    a = ap.parse_args()
    t0 = time.time()
    pts = [tuple(float(v) for v in it.split(":")) for it in a.points.split(",")]
    levs = C.floats(a.levels)
    jobs = [(e, p, th, L * math.sqrt(e), a.q) for (e, p, th) in pts for L in levs]
    res = parallel_map(job, jobs, default_workers(a.workers), print,
                       lambda i, r, d, n: f"[f68f_1d] {d}/{n}: eps {jobs[i][0]:g} p {jobs[i][1]:g} theta {jobs[i][2]:g} "
                                          f"kappa {jobs[i][3]:.2f}: A {r:.10f} ({time.time() - t0:.0f}s)")
    from diag import f68f_1c as S
    lines = [f"F6.8f-dev step 1d: needle volume at finite kappa and b (real medium) against the scaling-limit series, "
             f"q = {a.q:g}; development, decides nothing", ""]
    data = []
    for k, (e, p, th) in enumerate(pts):
        Ak = res[k * len(levs):(k + 1) * len(levs)]
        As = S.series_coefficients(e, p, [th])["A"][0]
        ser = As[0] + As[1] * a.q ** 2 + As[2] * a.q ** 4 + As[3] * a.q ** 6
        row = f"  eps {e:g} p {p:g} theta {th:g}: " + "  ".join(
            f"kappa {L * math.sqrt(e):.2f} (b {e / (L * math.sqrt(e)) ** 2:.4f}): {A:.10f}" for L, A in zip(levs, Ak))
        lines.append(row)
        ext = extrapolate(levs, Ak)
        lines.append(f"      series (scaling limit) {ser:.10f};  finite kappa - series: "
                     + "  ".join(f"{A - ser:+.3e}" for A in Ak)
                     + f";  extrapolated (b -> 0) - series: {ext - ser:+.3e}")
        data.append({"eps": e, "p": p, "theta": th, "levels": levs, "A": Ak, "series": ser, "extrapolated": ext})
    lines.append(f"cost {time.time() - t0:.0f} s")
    stem = C.save(NAME, vars(a), CODE, data, lines, t0)
    print("\n".join(lines))
    print("saved", stem + ".json")


if __name__ == "__main__":
    main()
