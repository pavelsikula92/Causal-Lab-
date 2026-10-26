"""
core.operators.bd_expect_h1
===========================

The EXPECTED smeared Benincasa-Dowker operator of a causal set sprinkled
into the static geometry of hypothesis H1 (anisotropic, cones of speed
Pi^(1 + gamma)), to FIRST ORDER in the field. Instrument of K.13a.

Definitions (stated before any computation, 22 Oct 2026)
--------------------------------------------------------
Geometry (n = d_s spatial dimensions, D = n + 1, signature - + + +):

    ds^2 = -Pi^2 dt^2 + Pi^(-2 gamma) dx^2,     Pi = 1 + e u(|x|),

gamma = +1 is H1, gamma = -1 the conformal geometry Pi^2 eta of K.11.
The family is linear in gamma at first order, so one computation covers
every gamma (below).

Causal order. y precedes x iff a future-directed causal curve joins them.
The metric is static, so this holds iff

    t_x - t_y >= T(y_s, x_s),   T = inf over spatial paths of int Pi^-(1 + gamma) dl

(null curves move with coordinate speed dl/dt = Pi / Pi^-gamma =
Pi^(1 + gamma); gamma = 1 gives c_eff = Pi^2, the H1 law of F6.3b).
Fermat: a change of the path changes T only at second order, so

    T = L - (1 + gamma) e U(y_s, x_s) + O(e^2),   U = integral of u along the straight segment,

L = |x_s - y_s|. The test K.13a checks this against a numerically
optimised path (part geom).

Sprinkling. Poisson with density rho per unit PROPER volume, i.e.
rho sqrt(-g) per coordinate volume, sqrt(-g) = Pi^(1 - gamma n)
= 1 + (1 - gamma n) e u. The number of elements in the interval [y, x] is
Poisson with mean rho V(x, y), V the proper volume of the interval. For a
static metric the time extent of the interval above a spatial point z is
Delta t - T(y_s, z) - T(z, x_s) when positive, hence

    V(x, y) = int d^n z sqrt(-g)(z) [Delta t - T(y_s, z) - T(z, x_s)]_+ .

(With Pi = 1 this is the flat diamond projected onto the prolate spheroid
d1 + d2 <= Delta t, as in core.operators.bd_expect_curved.)

The operator: the same smeared BD operator as K.4/K.5/K.11 (layers C_i of
core.operators.bd.BD_COEFFS, d = 4: -1, 9, -16, 8; smearing eps;
normalisation a_d = BD_NORM, d = 4: 4/sqrt 6; fundamental density rho in
the prefactor). Its expectation for a static field w:

    E[B w](x) = -a_d (eps rho)^(2/d) [ w(x) + rho eps
                int_{J^-(x)} w(y) gk(eps rho V(x, y)) sqrt(-g)(y) d^D y ],
    gk(w) = sum_i C_i w^i / i! e^-w .

Coordinates that keep the past cone fixed. For y in the past of x take
the spatial offset rho_s = y_s - x_s and s = l2 (Delta t^2 - T^2),
l2 = (eps rho kappa_D)^(2/D), kappa_D = interval_volume(D, 1). The past
cone is s >= 0 for every geometry, so its deformation enters only through
Delta t = sqrt(T^2 + s / l2) and d(Delta t) = ds / (2 l2 Delta t): no
boundary term. With T = L this is exactly the flat parametrisation of
core.operators.bd_expect (same nodes, same panels).

First order. With e factored out and Delta t_0, V_0 the flat values,
dE[w](x) = -a_d (eps rho)^(2/d) rho eps dI[w], and

    dI[w] = (1 - gamma n) P[w] + (1 + gamma) Q[w],

    P[w] = int m0 w(y) [ gk(w0) u(y_s) + gk'(w0) eps rho A1 ]
    Q[w] = int m0 w(y) [ gk(w0) L U_yx / Delta t_0^2 + gk'(w0) eps rho A2 ]

    A1 = int_spheroid u(z) l0(z) d^n z                      (density inside the interval)
    A2 = int_spheroid [U(y_s, z) + U(z, x_s) - (L / Delta t_0) U(y_s, x_s)] d^n z
                                                            (cone shift inside the interval)

m0 is the flat measure of the node, w0 = eps rho V_0, l0 = Delta t_0 - d1 - d2.
Origin of the terms: sqrt(-g)(y) gives gk u(y); sqrt(-g) inside the interval
gives A1; T(y, x) in Delta t gives the measure term L U / Delta t^2 and,
through V_0(Delta t), the last term of A2 (dV_0 / dDelta t = int_spheroid d^n z);
T(y, z) + T(z, x) inside the interval gives the first two terms of A2. The
combination in A2 is kept in one integrand because for long, nearly null
intervals the two parts almost cancel.

P is the first-order response to the density alone with flat cones (the
conformal geometry with Omega^D = sqrt(-g), the K.11 instrument), Q the
response to the change of cone speed. Their local limits (sympy, C7 of
docs/derivations/local_form_symbolic.py):

    P -> (-2 u lap w + (n-1) grad u.grad w + n w lap u) / D
    Q -> ( 2 u lap w +     2 grad u.grad w +   w lap u) / D

and (1 - gamma n) P + (1 + gamma) Q -> d/de (box_g - R/2) w (C6).

Constant part. Near the apex every term above is of size u(x), and the
terms cancel to the response, which is smaller by (l_k / sigma)^2. The
field is therefore split as u = u(x) + (u - u(x)): the deviation goes
through the quadrature, the constant part is added exactly. For constant
u, A1 = u V_0, U(y, x) = u L and A2 = n u V_0 (A2 = int (d1 + d2 - L^2/dt)
d^n z = (dt^2 - L^2)/dt dV_0/d(dt) - V_0 = (D - 1) V_0), so the constant
part of the integrands is u(x) [gk + gk' w0] for P and
u(x) [gk L^2/dt^2 + gk' n w0] for Q, both closed form in s. Long, nearly
null intervals need two more cancellation-free forms (dt - L, and
S = U(y,z) + U(z,x) - U(y,x)): see _spheroid_A and
GaussianField.broken_excess.

Field u: a Gaussian u(r) = exp(-r^2 / 2 sigma^2) (segment integrals in
closed form with erf). Region: |y_s| <= R, rows at t_x large enough that
the bottom t = 0 never enters (checked).

Functions
---------
    GaussianField(sigma)                 u, its segment integral and derivatives
    first_order_row(D, rho, eps, t_x, r_x, R, u, w, w_x, quad)
            P[1], Q[1], P[w - w(x)], Q[w - w(x)] and the flat I_one, I_diff at one row
    response(parts, gamma, D, rho, eps, w_x)
            dE[w] for a given gamma from the parts
    continuum_parts(n, r, u, w)          the local limits of P[w], Q[w] (w a GaussianField or None for w = 1)
"""

from __future__ import annotations

import math
from typing import Dict, Optional, Sequence

import numpy as np
from scipy.special import erf

from .bd import BD_COEFFS, BD_NORM
from .bd_expect import _directions, _s_edges, _wall_distance, g_kernel
from ..order.sprinkle import interval_volume


# =============================================================================
# the field
# =============================================================================

_GL8 = np.polynomial.legendre.leggauss(8)


def _erfdiff(a: np.ndarray, h: np.ndarray) -> np.ndarray:
    """erf(a + h) - erf(a) without cancellation: 8-point Gauss-Legendre of (2/sqrt pi) e^-t^2 for |h| < 0.5."""
    a = np.asarray(a, float)
    h = np.asarray(h, float)
    small = np.abs(h) < 0.5
    hs = np.where(small, h, 0.0)
    acc = np.zeros(np.broadcast(a, h).shape)
    for xi, wi in zip(*_GL8):
        t = a + 0.5 * hs * (1.0 + xi)
        acc = acc + wi * np.exp(-t * t)
    quad = acc * 0.5 * hs * (2.0 / math.sqrt(math.pi))
    return np.where(small, quad, erf(a + h) - erf(a))


class GaussianField:
    """u(r) = exp(-r^2 / 2 sigma^2) in n spatial dimensions: values, segment integrals and the local derivatives."""

    def __init__(self, sigma: float):
        self.sigma = float(sigma)

    def value(self, r2: np.ndarray) -> np.ndarray:
        """u at squared radius r2."""
        return np.exp(-np.asarray(r2, float) / (2.0 * self.sigma ** 2))

    def segment(self, a: Sequence[np.ndarray], b: Sequence[np.ndarray]) -> np.ndarray:
        """
        int_0^1 u(a + t (b - a)) |b - a| dt for vectors given by components (a[0], a[1], ...),
        each component an array (broadcast together). Closed form: with d = (b - a) / L,
        m^2 = |a|^2 - (a.d)^2 (squared distance of the line from the origin),
        U = sigma sqrt(pi/2) exp(-m^2 / 2 sigma^2) [erf(b.d / sqrt2 sigma) - erf(a.d / sqrt2 sigma)].
        """
        s = self.sigma
        diff = [bb - aa for aa, bb in zip(a, b)]
        L2 = sum(dd * dd for dd in diff)
        L = np.sqrt(L2)
        ok = L > 1e-14
        Ls = np.where(ok, L, 1.0)
        ad = sum(aa * dd for aa, dd in zip(a, diff)) / Ls
        bd = sum(bb * dd for bb, dd in zip(b, diff)) / Ls
        a2 = sum(aa * aa for aa in a)
        m2 = np.maximum(a2 - ad * ad, 0.0)
        k = 1.0 / (math.sqrt(2.0) * s)
        U = s * math.sqrt(math.pi / 2.0) * np.exp(-m2 / (2.0 * s * s)) * (erf(bd * k) - erf(ad * k))
        return np.where(ok, U, 0.0)

    def broken_excess(self, y: Sequence[np.ndarray], z: Sequence[np.ndarray], x: Sequence[np.ndarray]) -> np.ndarray:
        """
        S = U(y, z) + U(z, x) - U(y, x) without cancellation (needed for long, thin intervals, where the
        three integrals are ~ sigma and S is many orders smaller). With unit vectors d0 (y -> x), d1 (y -> z),
        d2 (z -> x), k = 1 / (sqrt2 sigma), K = sigma sqrt(pi/2), e_i = exp(-m_i^2 / 2 sigma^2):
            S = K [ (e1 - e0) F1 + e0 (F1 - G1) ] + K [ (e2 - e0) F2 + e0 (F2 - G2) ],
            F1 = erf(k z.d1) - erf(k y.d1),  G1 = erf(k z.d0) - erf(k y.d0)   (U(y, x) split at the foot of z),
            F2 = erf(k x.d2) - erf(k z.d2),  G2 = erf(k x.d0) - erf(k z.d0),
        where every difference of nearby quantities is formed from small vectors (d1 - d0, d2 - d0):
            m1^2 - m0^2 = (y.(d0 - d1)) (y.(d0 + d1)),  m2^2 - m0^2 = (x.(d0 - d2)) (x.(d0 + d2)),
            F1 - G1 = erfdiff(k z.d0, k z.(d1 - d0)) - erfdiff(k y.d0, k y.(d1 - d0)),  similarly F2 - G2.
        """
        s = self.sigma
        k = 1.0 / (math.sqrt(2.0) * s)
        K = s * math.sqrt(math.pi / 2.0)

        def unit(a, b):
            d = [bb - aa for aa, bb in zip(a, b)]
            Ln = np.sqrt(sum(dd * dd for dd in d))
            Ls = np.where(Ln > 0, Ln, 1.0)
            return [dd / Ls for dd in d], Ln

        def dot(a, b):
            return sum(aa * bb for aa, bb in zip(a, b))

        d0, L0 = unit(y, x)
        d1, L1 = unit(y, z)
        d2, L2 = unit(z, x)
        d10 = [p - q for p, q in zip(d1, d0)]
        d20 = [p - q for p, q in zip(d2, d0)]
        y2, x2 = dot(y, y), dot(x, x)
        yd0, xd0, zd0 = dot(y, d0), dot(x, d0), dot(z, d0)
        m0 = np.maximum(y2 - yd0 * yd0, 0.0)
        e0 = np.exp(-m0 / (2.0 * s * s))
        # first piece (y -> z against y -> foot of z)
        yd10, zd10 = dot(y, d10), dot(z, d10)
        dm1 = -yd10 * (2.0 * yd0 + yd10)                       # m1^2 - m0^2 = (y.d0)^2 - (y.d1)^2
        F1 = erf(k * (zd0 + zd10)) - erf(k * (yd0 + yd10))
        FG1 = _erfdiff(k * zd0, k * zd10) - _erfdiff(k * yd0, k * yd10)
        piece1 = e0 * np.expm1(-dm1 / (2.0 * s * s)) * F1 + e0 * FG1
        # second piece (z -> x against foot of z -> x)
        xd20, zd20 = dot(x, d20), dot(z, d20)
        dm2 = -xd20 * (2.0 * xd0 + xd20)                       # m2^2 - m0^2 = (x.d0)^2 - (x.d2)^2
        F2 = erf(k * (xd0 + xd20)) - erf(k * (zd0 + zd20))
        FG2 = _erfdiff(k * xd0, k * xd20) - _erfdiff(k * zd0, k * zd20)
        piece2 = e0 * np.expm1(-dm2 / (2.0 * s * s)) * F2 + e0 * FG2
        ok = (L0 > 0) & (L1 > 0) & (L2 > 0)
        return np.where(ok, K * (piece1 + piece2), 0.0)

    def radial(self, r: np.ndarray, n: int):
        """(u, u'(r), lap u) at radius r in n spatial dimensions."""
        r = np.asarray(r, float)
        s2 = self.sigma ** 2
        u = np.exp(-r * r / (2.0 * s2))
        return u, -u * r / s2, u * (r * r / (s2 * s2) - n / s2)


def continuum_parts(n: int, r: np.ndarray, u: GaussianField, w: Optional[GaussianField]) -> Dict[str, np.ndarray]:
    """
    Local limits of P[w] and Q[w] (module notes; w = None means w = 1) at radii r,
    and of the total d/de (box_g - R/2) w for a given gamma via total(gamma).
    """
    D = n + 1
    uu, du, lu = u.radial(r, n)
    if w is None:
        ww, dw, lw = np.ones_like(uu), np.zeros_like(uu), np.zeros_like(uu)
    else:
        ww, dw, lw = w.radial(r, n)
    gug = du * dw
    P = (-2.0 * uu * lw + (n - 1) * gug + n * ww * lu) / D
    Q = (2.0 * uu * lw + 2.0 * gug + ww * lu) / D
    return {"P": P, "Q": Q, "box_only": lambda g: 2.0 * g * uu * lw + (1.0 - g * (n - 2)) * gug,
            "total": lambda g: (1.0 - g * n) * P + (1.0 + g) * Q}


# =============================================================================
# kernel
# =============================================================================

def g_kernel_prime(w: np.ndarray, d: int) -> np.ndarray:
    """d/dw of gk(w) = sum_i C_i w^i/i! e^-w: sum_i C_i (w^(i-1)/(i-1)! - w^i/i!) e^-w."""
    s = np.zeros_like(w)
    for i, c in enumerate(BD_COEFFS[d]):
        term = -w ** i / math.factorial(i)
        if i > 0:
            term = term + w ** (i - 1) / math.factorial(i - 1)
        s = s + c * term
    return s * np.exp(-w)


# =============================================================================
# the spheroid integrals A1, A2 over flat diamonds
# =============================================================================

def _spheroid_A(u: GaussianField, n: int, r_x: float, CO: float, L: np.ndarray, dt: np.ndarray,
                excess: np.ndarray, n_mu: int, n_nu: int, n_phi: int, u_x: float = 0.0):
    """
    A1 = int u(z) l0(z) d^n z and A2 = int [U(y,z) + U(z,x) - (L/dt) U(y,x)] d^n z over the flat diamond
    between y = x + L nhat and x (time difference dt), nhat = (CO, SN, 0), for the DEVIATION field
    u - u_x (u_x = u at x; first_order_row adds the constant part exactly, see there). Same prolate
    spheroidal coordinates and quadrature as core.operators.bd_expect_curved.diamond_phi_spheroid
    (foci y and x): sigma = (d1 + d2)/L, te = (d1 - d2)/L, d1 = c (sigma + te), d2 = c (sigma - te).
    Shapes: L (P, 1), dt (P, Q) -> A1, A2, U(y, x) of the deviation (P, Q).

    Precision. excess = dt - L is passed in, computed without cancellation ((s/l2) / (dt + L)). Long,
    nearly null intervals have dt - L many orders below L; sigma_max - 1 = excess / L, the time extent
    l0 = excess (1 - x) (n = 3) and the constant-field part d1 + d2 - (L/dt) L = L (sigma - 1) + L excess/dt
    are formed from it, never as differences of numbers of size L (a development run, 22 Oct, showed
    garbage needle averages for L ~ 500 when sigma_max - 1 was computed as dt/L - 1).
    """
    SN = math.sqrt(max(1.0 - CO * CO, 0.0))
    L = np.broadcast_to(L, dt.shape)
    excess = np.broadcast_to(excess, dt.shape)
    c = 0.5 * L
    sm1 = excess / np.maximum(L, 1e-300)                 # sigma_max - 1
    # end points as components (x on axis 1, y = x + L nhat)
    xv = (np.full(dt.shape, r_x), np.zeros(dt.shape)) + ((np.zeros(dt.shape),) if n == 3 else ())
    yv = (r_x + L * CO, L * SN) + ((np.zeros(dt.shape),) if n == 3 else ())
    U_full = u.segment(yv, xv)
    ratio = L / dt
    one_m_ratio = excess / dt                            # 1 - L/dt
    xg, wg = np.polynomial.legendre.leggauss(n_mu)
    xg, wg = 0.5 * (xg + 1.0), 0.5 * wg
    A1 = np.zeros(dt.shape)
    A2 = np.zeros(dt.shape)

    def point(alpha, p1, p2):
        """z = x + alpha nhat + p1 e1perp + p2 e2perp, e1perp = (SN, -CO, 0), e2perp = (0, 0, 1)."""
        z = (r_x + alpha * CO + p1 * SN, alpha * SN - p1 * CO)
        return z + ((p2,) if n == 3 else ())

    def a2_integrand(z, sig_m1):
        """U(y,z) + U(z,x) - (L/dt) U(y,x) for u - u_x = S + (1 - L/dt) U(y,x) - u_x L (sigma - 1 + 1 - L/dt),
        S = U(y,z) + U(z,x) - U(y,x) from GaussianField.broken_excess (all terms without cancellation)"""
        return u.broken_excess(yv, z, xv) + one_m_ratio * U_full - u_x * L * (sig_m1 + one_m_ratio)

    if n == 3:
        tg, wt = np.polynomial.legendre.leggauss(n_nu)
        phis = 2.0 * np.pi * (np.arange(n_phi) + 0.5) / n_phi
        for x, wx in zip(xg, wg):
            sig_m1 = sm1 * x                             # sigma - 1
            sig = 1.0 + sig_m1
            jac = sm1 * wx
            ell = excess * (1.0 - x)                     # dt - L sigma
            for te, wte in zip(tg, wt):
                xi = c * (1.0 + sig * te)                # distance from y along the axis towards x
                rp = c * np.sqrt(np.maximum(sig_m1 * (sig + 1.0) * (1.0 - te * te), 0.0))
                wvol = jac * wte * c ** 3 * ((1.0 - te * te) + sig_m1 * (sig + 1.0)) * (2.0 * np.pi / n_phi)
                alpha = L - xi                           # offset from x towards y
                for ph in phis:
                    z = point(alpha, rp * math.cos(ph), rp * math.sin(ph))
                    r2 = sum(zz * zz for zz in z)
                    A1 += wvol * ell * (u.value(r2) - u_x)
                    A2 += wvol * a2_integrand(z, sig_m1)
    elif n == 2:
        nus = 2.0 * np.pi * (np.arange(n_nu) + 0.5) / n_nu
        for x, wx in zip(xg, wg):
            ch_m1 = sm1 * x * x                          # cosh(mu) - 1
            ch = 1.0 + ch_m1
            sh = np.sqrt(np.maximum(ch_m1 * (ch + 1.0), 0.0))
            jac = wx * 2.0 * np.sqrt(sm1) / np.sqrt(ch + 1.0)
            ell = excess * (1.0 - x * x)                 # dt - L cosh(mu)
            for nu in nus:
                xi = c * (1.0 + ch * math.cos(nu))
                rp = c * sh * math.sin(nu)               # signed: both sides of the axis
                wvol = jac * c * c * (math.sin(nu) ** 2 + ch_m1 * (ch + 1.0)) * (2.0 * np.pi / n_nu)
                z = point(L - xi, rp, None)
                r2 = sum(zz * zz for zz in z)
                A1 += wvol * ell * (u.value(r2) - u_x)
                A2 += wvol * a2_integrand(z, ch_m1)
    else:
        raise ValueError("D must be 3 or 4")
    return A1, A2, U_full - u_x * L


# =============================================================================
# one row of the first-order operator
# =============================================================================

def first_order_row(D: int, rho: float, eps: float, t_x: float, r_x: float, R: float, u: GaussianField,
                    w: Optional[GaussianField], quad: Sequence[int]) -> Dict[str, float]:
    """
    Parts of the first-order operator at x = (t_x, r_x on axis 1), quad = (n_rho, n_s, n_ang, n_mu, n_nu, n_phi):
        P1, Q1   = P[1], Q[1]
        Pw, Qw   = P[w - w(x)], Q[w - w(x)]   (0 if w is None)
        I_one, I_diff = the flat integrals of core.operators.bd_expect (Pi = 1 check)
    Nodes in (rho_s, s) exactly as core.operators.bd_expect._nodes.
    """
    n_rho, n_s, n_ang, n_mu, n_nu, n_phi = (int(q) for q in quad)
    n = D - 1
    lam = eps * rho * interval_volume(D, 1.0)
    l2 = lam ** (2.0 / D)
    if r_x == 0.0:
        # the row on the axis of the radial fields: every direction gives the same integrand (exact symmetry),
        # one direction with the full solid angle replaces the angular quadrature
        cos, wa = np.array([1.0]), np.array([float(_directions(n, n_ang)[1].sum())])
    else:
        cos, wa = _directions(n, n_ang)
    z, wz = np.polynomial.legendre.leggauss(n_rho)
    z, wz = 0.5 * (z + 1.0), 0.5 * wz
    uu, wu = np.polynomial.legendre.leggauss(n_s)
    uu, wu = 0.5 * (uu + 1.0), 0.5 * wu
    edges = _s_edges(D)
    w_x = 1.0 if w is None else float(w.value(np.array(r_x * r_x)))
    u_x = float(u.value(np.array(r_x * r_x)))
    out = {k: 0.0 for k in ("P1", "Q1", "Pw", "Qw", "I_one", "I_diff")}
    for a in range(cos.size):
        CO = float(cos[a])
        SN = math.sqrt(max(1.0 - CO * CO, 0.0))
        rho_lim = min(t_x, float(_wall_distance(r_x, np.array([CO]), R)[0]))
        rs = rho_lim * z ** 2
        wr = wz * 2.0 * rho_lim * z
        s_max = l2 * np.maximum(t_x * t_x - rs * rs, 0.0)
        if np.any(s_max < edges[-1]):
            raise ValueError("the bottom t = 0 enters the past of the row: raise t_x")
        S, WS = [], []
        for lo, hi in zip(edges[:-1], edges[1:]):
            S.append(lo + (hi - lo) * uu[None, :] + 0.0 * rs[:, None])
            WS.append(np.broadcast_to((hi - lo) * wu, (rs.size, n_s)))
        s = np.concatenate(S, axis=1)
        ws = np.concatenate(WS, axis=1)
        RS = rs[:, None]
        dt = np.sqrt(RS * RS + s / l2)
        radial = RS ** (n - 1) if n > 1 else np.ones_like(RS)
        m0 = wa[a] * wr[:, None] * ws * radial / (l2 * 2.0 * dt)
        w0 = s ** (D / 2.0)                                      # eps rho V_0
        gk, gkp = g_kernel(w0, D), g_kernel_prime(w0, D)
        r2_y = r_x * r_x + RS * RS + 2.0 * r_x * RS * CO
        u_y = u.value(r2_y)
        excess = (s / l2) / (dt + RS)                            # dt - L without cancellation
        A1, A2, U_yx = _spheroid_A(u, n, r_x, CO, RS, dt, excess, n_mu, n_nu, n_phi, u_x)
        # deviation part (u - u_x) by quadrature + constant part u_x exactly (module notes, "Constant part")
        p_int = gk * (u_y - u_x) + gkp * eps * rho * A1 + u_x * (gk + gkp * w0)
        q_int = gk * RS * U_yx / (dt * dt) + gkp * eps * rho * A2 + u_x * (gk * RS * RS / (dt * dt) + gkp * n * w0)
        wdiff = np.zeros_like(dt) if w is None else (w.value(r2_y) - w_x) * np.ones_like(dt)
        out["P1"] += float((m0 * p_int).sum())
        out["Q1"] += float((m0 * q_int).sum())
        out["Pw"] += float((m0 * wdiff * p_int).sum())
        out["Qw"] += float((m0 * wdiff * q_int).sum())
        out["I_one"] += rho * eps * float((m0 * gk).sum())
        out["I_diff"] += rho * eps * float((m0 * wdiff * gk).sum())
    out["w_x"] = w_x
    return out


def response(parts: Dict[str, float], gamma: float, D: int, rho: float, eps: float, field: str = "w") -> Dict[str, float]:
    """
    dE (per unit e) at the row for a given gamma: field 'w' (the test field) or 'one' (w = 1).
    Returns the total and its density part (1 - gamma n) P and cone part (1 + gamma) Q.
    """
    n = D - 1
    pref = -BD_NORM[D] * (eps * rho) ** (2.0 / D) * rho * eps
    if field == "one":
        P, Q = parts["P1"], parts["Q1"]
    else:
        P = parts["w_x"] * parts["P1"] + parts["Pw"]
        Q = parts["w_x"] * parts["Q1"] + parts["Qw"]
    dens = pref * (1.0 - gamma * n) * P
    cone = pref * (1.0 + gamma) * Q
    return {"total": dens + cone, "density": dens, "cone": cone, "P": pref * P, "Q": pref * Q}


def flat_value(parts: Dict[str, float], D: int, rho: float, eps: float) -> float:
    """E[B w](x) of the flat operator from the same nodes (Pi = 1 check against bd_expect)."""
    pref = -BD_NORM[D] * (eps * rho) ** (2.0 / D)
    return pref * (parts["w_x"] * (1.0 + parts["I_one"]) + parts["I_diff"])
