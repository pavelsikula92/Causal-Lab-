"""Numba kernels of OV2 in 2D (diag/ov2_2d.py); a module never run as __main__ (see diag/ov1_kernels.py)."""
import math

import numba as nb
import numpy as np

from diag.ov1_kernels import small_future


@nb.njit(cache=True)
def matter_sums(vr, iv, ncut, f, phi):
    """Over all related pairs y < x (x later) with n(x, y) <= ncut, for the elements in u-rank order:
    D = sum f(n) (phi_x - phi_y)^2           (pairwise, Dirichlet form)
    P = sum f(n) phi_x phi_y                  (the cross term of the plain form)
    R = sum f(n) (phi_y - phi_x) phi_x        (row-difference form, x the later element)
    F = per element: sum over its past of f(n) (for B1)."""
    N = vr.size
    tree = np.zeros(N + 1, np.int64)
    undo = np.empty(N, np.int64)
    oy = np.empty(N, np.int64)
    on = np.empty(N, np.int64)
    D = 0.0
    P = 0.0
    R = 0.0
    Fpast = np.zeros(N)
    for x in range(N):                       # x earlier, y later
        m = small_future(x, vr, iv, ncut, tree, undo, oy, on)
        for k in range(m):
            y = oy[k]
            w = f[on[k]]
            d = phi[y] - phi[x]
            D += w * d * d
            P += w * phi[x] * phi[y]
            R += w * (phi[x] - phi[y]) * phi[y]
            Fpast[y] += w
    return D, P, R, Fpast


# ------------------------------------------------------------------------------------------------- quadrature (M-b, 2D)
@nb.njit(cache=True, inline="always")
def _sig(u, v, a, s):
    return a * np.exp(-((u - 0.5) ** 2 + (v - 0.5) ** 2) / (2.0 * s * s))


@nb.njit(cache=True, inline="always")
def _phi(u, v, m):
    return 4.0 * np.arctan(np.exp(m * 0.5 * (v - u)))


@nb.njit(cache=True, inline="always")
def _rect(C, M, u1, v1, u2, v2):
    """int_{u1}^{u2} int_{v1}^{v2} e^{2 sigma} du dv from the cumulative table C (bilinear)."""
    def cum(u, v):
        x = u * M
        y = v * M
        i = min(int(x), M - 1)
        j = min(int(y), M - 1)
        fx = x - i
        fy = y - j
        return ((1 - fx) * (1 - fy) * C[i, j] + fx * (1 - fy) * C[i + 1, j] + (1 - fx) * fy * C[i, j + 1]
                + fx * fy * C[i + 1, j + 1])
    return cum(u2, v2) - cum(u1, v2) - cum(u2, v1) + cum(u1, v1)


@nb.njit(cache=True)
def kin_quad(n_exp, eps, m, a, s, C, M, xo, wo, xt, wt, xw, ww, pmin, wcap):
    """<sum_{pairs} f(n) (phi_x - phi_y)^2> = int int lam_x lam_y g(eps Lam) (dphi)^2 over x < y (Mecke), with
    lam = n_exp e^{2 sigma} per du dv, Lam = n_exp int_I e^{2 sigma}, g(z) = e^{-z} (1 - 2 z + z^2 / 2).
    Outer Gauss-Legendre (xo, wo) on [0, 1]^2; inner: t = ln p (xt, wt on [0, 1] mapped), w = eps n_exp p q
    (xw, ww on [0, 1] mapped to [0, wmax])."""
    tot = 0.0
    no = xo.size
    for i in range(no):
        u1 = xo[i]
        for j in range(no):
            v1 = xo[j]
            lx = n_exp * np.exp(2.0 * _sig(u1, v1, a, s))
            px = _phi(u1, v1, m)
            tlo = np.log(pmin)
            thi = np.log(1.0 - u1)
            if thi <= tlo:
                continue
            acc = 0.0
            for k in range(xt.size):
                t = tlo + (thi - tlo) * xt[k]
                p = np.exp(t)
                dp = p * (thi - tlo) * wt[k]
                wmax = min(wcap, eps * n_exp * p * (1.0 - v1))
                if wmax <= 0.0:
                    continue
                for l in range(xw.size):
                    w = wmax * xw[l]
                    q = w / (eps * n_exp * p)
                    dq = wmax * ww[l] / (eps * n_exp * p)
                    u2 = u1 + p
                    v2 = v1 + q
                    lam = n_exp * _rect(C, M, u1, v1, u2, v2)
                    z = eps * lam
                    g = np.exp(-z) * (1.0 - 2.0 * z + 0.5 * z * z)
                    ly = n_exp * np.exp(2.0 * _sig(u2, v2, a, s))
                    d = _phi(u2, v2, m) - px
                    acc += ly * g * d * d * dp * dq
            tot += lx * acc * wo[i] * wo[j]
    return tot


@nb.njit(cache=True, inline="always")
def _gint(a1, a2, s):
    """int_{a1}^{a2} exp(-(x - 1/2)^2 / (2 s^2)) dx, exact (erf), with the midpoint rule for tiny intervals."""
    d = a2 - a1
    if d < 1e-7 * s:
        xm = 0.5 * (a1 + a2) - 0.5
        return d * np.exp(-xm * xm / (2.0 * s * s))
    c = s * np.sqrt(np.pi / 2.0)
    r = 1.0 / (s * np.sqrt(2.0))
    return c * (math.erf((a2 - 0.5) * r) - math.erf((a1 - 0.5) * r))


@nb.njit(cache=True)
def kin_quad_d1(n_exp, eps, m, s, xo, wo, xt, wt, xw, ww, pmin, wcap):
    """d/da at a = 0 of <sum f(n) (dphi)^2> for sigma = a sh, sh = exp(-r^2 / 2 s^2) (flat order, analytic
    derivative under the integral): int int N^2 [2 (sh_x + sh_y) g(z) + g'(z) eps N 2 int_I sh] (dphi)^2,
    z = eps N p q, g'(z) = e^{-z} (-3 + 3 z - z^2 / 2)."""
    tot = 0.0
    no = xo.size
    for i in range(no):
        u1 = xo[i]
        for j in range(no):
            v1 = xo[j]
            shx = np.exp(-((u1 - 0.5) ** 2 + (v1 - 0.5) ** 2) / (2.0 * s * s))
            px = _phi(u1, v1, m)
            tlo = np.log(pmin)
            thi = np.log(1.0 - u1)
            if thi <= tlo:
                continue
            acc = 0.0
            for k in range(xt.size):
                t = tlo + (thi - tlo) * xt[k]
                p = np.exp(t)
                dp = p * (thi - tlo) * wt[k]
                wmax = min(wcap, eps * n_exp * p * (1.0 - v1))
                if wmax <= 0.0:
                    continue
                Iu = _gint(u1, u1 + p, s)
                for l in range(xw.size):
                    w = wmax * xw[l]
                    q = w / (eps * n_exp * p)
                    dq = wmax * ww[l] / (eps * n_exp * p)
                    u2 = u1 + p
                    v2 = v1 + q
                    z = w
                    ez = np.exp(-z)
                    g = ez * (1.0 - 2.0 * z + 0.5 * z * z)
                    gp = ez * (-3.0 + 3.0 * z - 0.5 * z * z)
                    shy = np.exp(-((u2 - 0.5) ** 2 + (v2 - 0.5) ** 2) / (2.0 * s * s))
                    dL = 2.0 * n_exp * Iu * _gint(v1, v2, s)
                    d = _phi(u2, v2, m) - px
                    acc += (2.0 * (shx + shy) * g + gp * eps * dL) * d * d * dp * dq
            tot += n_exp * n_exp * acc * wo[i] * wo[j]
    return tot
