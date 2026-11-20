"""Numba kernels of OV2 T2 in 3+1 (diag/ov2_t2.py); a module never run as __main__ (see diag/ov1_kernels.py).

Units eps rho = 1 (l_k = 1, as F6.8). A static wall phi(z) = 4 arctan(e^{m z}) on a flat 3+1 order of age T (the past
of every x is cut at the time T before x). For x at height z and y in its past at spatial offset rho_s n (n . z = mu)
and proper time tau: w = eps rho V = kappa tau^4 (kappa = pi / 24), s = sqrt(kappa) tau^2 = sqrt(w),
dV = rho_s^2 drho_s 2 pi dmu ds / (2 sqrt(kappa) dt), dt = sqrt(rho_s^2 + tau^2), s <= sqrt(kappa) (T^2 - rho_s^2).
"""
import math

import numba as nb
import numpy as np

KAPPA = math.pi / 24.0
SQK = math.sqrt(KAPPA)


@nb.njit(cache=True, inline="always")
def gk(w):
    return math.exp(-w) * (-1.0 + 9.0 * w - 8.0 * w * w + (4.0 / 3.0) * w * w * w)


@nb.njit(cache=True, inline="always")
def gkp(w):
    return math.exp(-w) * (10.0 - 25.0 * w + 12.0 * w * w - (4.0 / 3.0) * w * w * w)


@nb.njit(cache=True, inline="always")
def phi_w(z, m):
    return 4.0 * math.atan(math.exp(m * z))


@nb.njit(cache=True)
def _Ar(r, T, sedges, g8x, g8w, mode):
    """A(r) = int_0^smax K(s^2) / (2 sqrt(kappa) dt) ds, smax = sqrt(kappa) (T^2 - r^2), on the panels sedges
    clipped exactly at smax (8-point Gauss-Legendre per panel)."""
    smax = SQK * (T * T - r * r)
    A = 0.0
    if smax <= 0.0:
        return 0.0
    for k in range(sedges.size - 1):
        a = sedges[k]
        if a >= smax:
            break
        b = min(sedges[k + 1], smax)
        h = 0.5 * (b - a)
        c = 0.5 * (b + a)
        for q in range(8):
            s = c + h * g8x[q]
            w = s * s
            dt = math.sqrt(r * r + s / SQK)
            if mode == 0:
                K = gk(w)
            else:
                K = 8.0 * gk(w) + 4.0 * w * gkp(w)
            A += h * g8w[q] * K / (2.0 * SQK * dt)
        if sedges[k + 1] >= smax:
            break
    # beyond the last edge (s > sqrt(80)) the kernel is below e^-80: omitted
    return A


@nb.njit(cache=True)
def _H(z, r, m, kb, g8x, g8w):
    """H = int_{-1}^{1} (phi(z + r mu) - phi(z))^2 dmu = (1/r) int_{z-r}^{z+r} (phi(zeta) - phi(z))^2 dzeta,
    composite Gauss-Legendre with breakpoints kb (around the kink) inside [z - r, z + r]."""
    lo = z - r
    hi = z + r
    pz = phi_w(z, m)
    pts = np.empty(kb.size + 2)
    n = 0
    pts[n] = lo
    n += 1
    for i in range(kb.size):
        if lo < kb[i] < hi:
            pts[n] = kb[i]
            n += 1
    pts[n] = hi
    n += 1
    H = 0.0
    for i in range(n - 1):
        a = pts[i]
        b = pts[i + 1]
        h = 0.5 * (b - a)
        c = 0.5 * (a + b)
        for q in range(8):
            d = phi_w(c + h * g8x[q], m) - pz
            H += h * g8w[q] * d * d
    return H / r


@nb.njit(cache=True)
def wall_kin(m, T, zs, wz, rs, wr, sedges, kb, g8x, g8w, mode):
    """sum_z wz int dV(y) K (phi(y) - phi(x))^2 (without the prefactor a_4/4); K by mode:
    0: gk(w), the flat kinetic action; 1: 8 gk + 4 w gk', d/db at b = 0 for a uniform conformal factor e^b."""
    tot = 0.0
    for iz in range(zs.size):
        z = zs[iz]
        acc = 0.0
        for ir in range(rs.size):
            r = rs[ir]
            A = _Ar(r, T, sedges, g8x, g8w, mode)
            if A == 0.0:
                continue
            acc += wr[ir] * r * r * 2.0 * math.pi * A * _H(z, r, m, kb, g8x, g8w)
        tot += wz[iz] * acc
    return tot


# ------------------------------------------------------------------------------------- stage 2: plane-wave conformal
@nb.njit(cache=True, inline="always")
def _sinc(x):
    if abs(x) < 1e-4:
        return 1.0 - x * x / 6.0
    return math.sin(x) / x


@nb.njit(cache=True, inline="always")
def _dsinc(x):
    if abs(x) < 1e-3:
        return -x / 3.0 + x * x * x / 30.0
    return (x * math.cos(x) - math.sin(x)) / (x * x)


@nb.njit(cache=True)
def dft4(om, q, h, g16x, g16w):
    """<cos(om t + q.x)> over the 3+1 causal diamond |x| <= h - |t| in its rest frame (port of diamond_ft4 of
    phases/F6_order_dynamics/F6_08a_bd_action_conformal.py, checked there)."""
    V = 2.0 * math.pi * h ** 4 / 3.0
    if q * h < 0.05:
        tot = 0.0
        for i in range(g16x.size):
            u = 0.5 * h * (g16x[i] + 1.0)
            qu = q * u
            B = (4.0 * math.pi / 3.0) * u ** 3 * (1.0 - qu * qu / 10.0 + qu ** 4 / 280.0)
            tot += 0.5 * h * g16w[i] * math.cos(om * (h - u)) * B
        return 2.0 / V * tot
    tot = 0.0
    for sgn in (-1.0, 1.0):
        alpha = q + sgn * om
        ph = -sgn * om * h
        x = 0.5 * alpha * h
        A = h * math.sin(ph + x) * _sinc(x)
        dA = h * (0.5 * h) * (math.cos(ph + x) * _sinc(x) + math.sin(ph + x) * _dsinc(x))
        tot += 0.5 * (A - q * dA)
    return (2.0 / V) * 4.0 * math.pi * tot / (q * q * q)


@nb.njit(cache=True)
def wall_kin_k(m, kw, T, zs, wz, rs, wr, sedges, kb, g8x, g8w, g16x, g16w):
    """d/db at b = 0 of the kinetic sum for Omega = e^{b cos(kw z)}: int dV [4 (c_x + c_y) g(w)
    + 4 w g'(w) <cos(kw z)>_I / 1] (dphi)^2, with <cos>_I = cos(kw z_c) F the diamond average (boosted to the frame of
    the interval: om' = kw dz / tau, q' = kw sqrt(1 + (dz / tau)^2), h = tau / 2)."""
    tot = 0.0
    for iz in range(zs.size):
        z = zs[iz]
        pz = phi_w(z, m)
        cx = math.cos(kw * z)
        acc = 0.0
        for ir in range(rs.size):
            r = rs[ir]
            smax = SQK * (T * T - r * r)
            if smax <= 0.0:
                continue
            lo = z - r
            hi = z + r
            pts = np.empty(kb.size + 2)
            n = 0
            pts[n] = lo
            n += 1
            for i in range(kb.size):
                if lo < kb[i] < hi:
                    pts[n] = kb[i]
                    n += 1
            pts[n] = hi
            n += 1
            inner = 0.0
            for ip in range(n - 1):
                a_ = pts[ip]
                b_ = pts[ip + 1]
                hh = 0.5 * (b_ - a_)
                cc = 0.5 * (a_ + b_)
                for qz in range(8):
                    zy = cc + hh * g8x[qz]
                    wq = hh * g8w[qz] / r                 # d mu = d z_y / r
                    d = phi_w(zy, m) - pz
                    if d == 0.0:
                        continue
                    dz = zy - z
                    cy = math.cos(kw * zy)
                    czc = math.cos(kw * 0.5 * (z + zy))
                    sacc = 0.0
                    for k in range(sedges.size - 1):
                        sa = sedges[k]
                        if sa >= smax:
                            break
                        sb = min(sedges[k + 1], smax)
                        sh = 0.5 * (sb - sa)
                        sc = 0.5 * (sb + sa)
                        for qs in range(8):
                            s = sc + sh * g8x[qs]
                            w = s * s
                            tau2 = s / SQK
                            tau = math.sqrt(tau2)
                            dt = math.sqrt(r * r + tau2)
                            U = dz / tau
                            F = dft4(kw * U, kw * math.sqrt(1.0 + U * U), 0.5 * tau, g16x, g16w)
                            K = 4.0 * (cx + cy) * gk(w) + 4.0 * w * gkp(w) * czc * F
                            sacc += sh * g8w[qs] * K / (2.0 * SQK * dt)
                        if sedges[k + 1] >= smax:
                            break
                    inner += wq * d * d * sacc
            acc += wr[ir] * r * r * 2.0 * math.pi * inner
        tot += wz[iz] * acc
    return tot
