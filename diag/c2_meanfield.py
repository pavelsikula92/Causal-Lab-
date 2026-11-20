"""C2 diagnostic (decides nothing): the mean-field (C1 exact-mean) prediction for a wall on the curved order.

Mecke: for a y,z-uniform wall the mean field f(t, x) obeys, with J = mu sin f replaced by its mean-field value,
    f(v) - delta mu sin f(v) = phi0(v) - c mu Int d^4u lambda(u) K(Lambda(u, v)) sin f(u),   t_u >= 0,
lambda = rho Pi^-2 (coordinate intensity), K(L) = sum_n w_n Pois(n; L), Lambda = rho * proper volume of I(u, v).
Geometry (optical H^3 of radius ell, H = x + ell): d0 = ell |ln(H_u/H_v)| at zero transverse offset; with the transverse
radius y, cosh(d/ell) = 1 + (dH^2 + y^2)/(2 H_u H_v), sigma = dt^2 - d^2.  2 pi y dy = pi H_u H_v sinh(d/ell)/(ell d) dsigma,
so the weight of (dt, x_u, sigma) is  c mu pi rho Pi_v^2 [sinh(d/ell)/(d/ell)] K(Lambda) dsigma dx_u dt.
Proper volume of the diamond: (pi/24) sigma^2 <Pi^4> with <Pi^4> ~ Pi_u^2 Pi_v^2 (exact for ln Pi linear along the
segment; the H^3 curvature correction to the diamond volume is O(sigma/ell^2)). ell = inf: the flat kernel (C0).
"""
import math, sys, time
import numpy as np
import numba as nb
from diag.c0_minvar import minvar
from diag.c0_apex import delta
import diag.c2_dev as C


@nb.njit(cache=True)
def table(xs, dts, J, DX, ell, rho, lam_cut, w, c, mu, gs, gws, gx, gwx):
    nx, nt, nn = xs.size, dts.size, w.size
    out = np.zeros((nx, nt, 2 * J + 1))
    flat = ell <= 0.0
    for ix in range(nx):
        xv = xs[ix]
        Hv = xv + ell
        P2v = 1.0 if flat else Hv / ell
        for it in range(1, nt):
            dt = dts[it]
            # x_u range: flat |dx| < dt; curved H_v e^{-dt/ell} < H_u < H_v e^{dt/ell}
            if flat:
                lo, hi = xv - dt, xv + dt
            else:
                lo, hi = Hv * math.exp(-dt / ell) - ell, Hv * math.exp(dt / ell) - ell
            for q in range(gx.size):
                xu = 0.5 * (lo + hi) + 0.5 * (hi - lo) * gx[q]
                wq = 0.5 * (hi - lo) * gwx[q]
                if flat:
                    d0 = abs(xu - xv)
                    P2u = 1.0
                else:
                    Hu = xu + ell
                    d0 = ell * abs(math.log(Hu / Hv))
                    P2u = Hu / ell
                smax = dt * dt - d0 * d0
                if smax <= 0.0:
                    continue
                scut = math.sqrt(24.0 * lam_cut / (math.pi * rho * P2u * P2v))
                if smax > scut:
                    smax = scut
                acc = 0.0
                for p in range(gs.size):
                    sig = 0.5 * smax * (gs[p] + 1.0)
                    wsig = 0.5 * smax * gws[p]
                    jac = 1.0
                    if not flat:
                        d = math.sqrt(dt * dt - sig)
                        a = d / ell
                        jac = math.sinh(a) / a if a > 1e-8 else 1.0
                    lam = rho * math.pi / 24.0 * sig * sig * P2u * P2v
                    pk = math.exp(-lam)
                    kk = 0.0
                    for n in range(nn):
                        if n > 0:
                            pk *= lam / n
                        kk += w[n] * pk
                    acc += wsig * jac * kk
                val = c * mu * math.pi * rho * P2v * acc * wq
                fj = (xu - xv) / DX
                j0 = int(math.floor(fj))
                fr = fj - j0
                if -J <= j0 <= J:
                    out[ix, it, J + j0] += val * (1.0 - fr)
                if -J <= j0 + 1 <= J:
                    out[ix, it, J + j0 + 1] += val * fr
    return out


@nb.njit(cache=True)
def solve(Wt, phi0, nx, nt, J, dtau, ct, mu, fl, fr):
    """f[it, ix]; trapezoid in dt; off-grid x take the plateau values fl (left) / fr (right)."""
    f = np.zeros((nt, nx))
    for it in range(nt):
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
                    if iy < 0:
                        fv = fl
                    elif iy >= nx:
                        fv = fr
                    else:
                        fv = f[src, iy]
                    s += wk * wv * math.sin(fv)
            rhs = phi0[it, ix] - dtau * s
            g = phi0[it, ix] if it == 0 else f[it - 1, ix]
            for _ in range(50):
                F = g - ct * mu * math.sin(g) - rhs
                dF = 1.0 - ct * mu * math.cos(g)
                dg = F / dF
                g -= dg
                if abs(dg) < 1e-13:
                    break
            f[it, ix] = g
    return f


def centre(xs, row, x0):
    r = row - math.pi
    j = np.nonzero(np.diff(np.sign(r)))[0]
    j = j[np.argmin(np.abs(xs[j] - x0))]
    return float(xs[j] - r[j] * (xs[j + 1] - xs[j]) / (r[j + 1] - r[j]))


def run(rho, ell, v, T=1.5, X=(-7.0, 9.0), DX=0.05, DT=0.025, ns=12, nfine=64, lam_cut=22.0, mu=1.0):
    xs = np.arange(X[0], X[1] + 1e-9, DX)
    ts = np.arange(0.0, T + 1e-9, DT)
    J = int(math.ceil((ts[-1] * (math.exp(ts[-1] / ell) if np.isfinite(ell) else 1.0) + 1.0) / DX)) + 1
    w = np.asarray(minvar(8), float)
    ct = delta(w, rho)
    c = 1.0 / (2.0 * math.pi * math.sqrt(6.0 * rho))
    gs, gws = np.polynomial.legendre.leggauss(ns)
    gx, gwx = np.polynomial.legendre.leggauss(nfine)
    a0 = time.time()
    W = table(xs, ts, J, DX, ell if np.isfinite(ell) else -1.0, rho, lam_cut, w, c, mu, gs, gws, gx, gwx)
    free = C.Continuum(ell, mu, v, 0.0, T + 0.1, X[0], X[1], mu_dyn=0.0)
    full = C.Continuum(ell, mu, v, 0.0, T + 0.1, X[0], X[1])
    TT, XX = np.meshgrid(ts, xs, indexing="ij")
    phi0 = free.at(TT.ravel(), XX.ravel()).reshape(TT.shape)
    f = solve(W, phi0, xs.size, ts.size, J, DT, ct, mu, 0.0, 2 * math.pi)
    rows = []
    for tb in (0.5, 1.0, 1.5):
        if tb > T + 1e-9:
            continue
        i = int(round(tb / DT))
        xm = centre(xs, f[i], full.centre(tb))
        xcont = full.centre(tb)
        rows.append((tb, xm, xcont))
    return rows, time.time() - a0


if __name__ == "__main__":
    rho = float(sys.argv[1]); ell = float(sys.argv[2]); v = float(sys.argv[3])
    DX = float(sys.argv[4]) if len(sys.argv) > 4 else 0.05
    DT = float(sys.argv[5]) if len(sys.argv) > 5 else 0.025
    rows, sec = run(rho, ell, v, DX=DX, DT=DT)
    print(f"rho {rho:g} ell {ell:g} v {v:+g} DX {DX} DT {DT} ({sec:.0f} s)")
    for tb, xm, xc in rows:
        print(f"  t {tb:.2f} mf {xm:+.5f} cont {xc:+.5f} diff {xm - xc:+.5f} ratio {xm / xc if abs(xc) > 1e-9 else float('nan'):.4f}")


def fit_readout(rho, ell, v, npy, T=1.8, **kw):
    """The sims' readout (fit_kink in slabs at the events of a sprinkling) applied to the mean field and the exact field."""
    from phases.F3_structure import field1d as F
    X = (-7.5, 9.5)
    DX, DT = kw.get("DX", 0.05), kw.get("DT", 0.025)
    xs = np.arange(X[0], X[1] + 1e-9, DX)
    ts = np.arange(0.0, T + 1e-9, DT)
    J = int(math.ceil((ts[-1] * (math.exp(ts[-1] / ell) if np.isfinite(ell) else 1.0) + 1.0) / DX)) + 1
    w = np.asarray(minvar(8), float)
    ct = delta(w, rho)
    c = 1.0 / (2.0 * math.pi * math.sqrt(6.0 * rho))
    gs, gws = np.polynomial.legendre.leggauss(kw.get("ns", 12))
    gx, gwx = np.polynomial.legendre.leggauss(kw.get("nfine", 64))
    W = table(xs, ts, J, DX, ell if np.isfinite(ell) else -1.0, rho, kw.get("lam_cut", 22.0), w, c, 1.0, gs, gws, gx, gwx)
    free = C.Continuum(ell, 1.0, v, 0.0, T + 0.1, X[0], X[1], mu_dyn=0.0)
    full = C.Continuum(ell, 1.0, v, 0.0, T + 0.1, X[0], X[1])
    TT, XX = np.meshgrid(ts, xs, indexing="ij")
    phi0 = free.at(TT.ravel(), XX.ravel()).reshape(TT.shape)
    f = solve(W, phi0, xs.size, ts.size, J, DT, ct, 1.0, 0.0, 2 * math.pi)
    out = []
    for fn in (npy if isinstance(npy, (list, tuple)) else [npy]):
        out.append(_readout(f, full, fn, T, X, DX, DT, xs, ts, kw.get('bs', (0.5, 1.0, 1.5))))
    return out


def _readout(f, full, npy, T, X, DX, DT, xs, ts, bs=(0.5, 1.0, 1.5)):
    from phases.F3_structure import field1d as F
    ev = np.load(npy)
    t, x = ev[0].astype(float), ev[1].astype(float)
    m = (t <= T) & (x > X[0]) & (x < X[1] - 1e-6)
    t, x = t[m], x[m]
    ft = t / DT; it = np.minimum(ft.astype(int), ts.size - 2); at = ft - it
    fx = (x - X[0]) / DX; ix = np.minimum(fx.astype(int), xs.size - 2); ax = fx - ix
    fm = ((1 - at) * ((1 - ax) * f[it, ix] + ax * f[it, ix + 1]) + at * ((1 - ax) * f[it + 1, ix] + ax * f[it + 1, ix + 1]))
    ex = full.at(t, x)
    out = []
    for b in bs:
        g = round(full.centre(b), 3) + 0.0
        r1 = F.fit_kink(t, x, fm, b - 0.25, b + 0.25, g)
        r2 = F.fit_kink(t, x, ex, b - 0.25, b + 0.25, g)
        out.append((b, r1["xc"], r2["xc"], r1["w"], r2["w"]))
    return out
