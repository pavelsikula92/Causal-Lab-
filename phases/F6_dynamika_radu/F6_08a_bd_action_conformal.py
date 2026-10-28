"""
F6.8a -- Static conformal quadratic form of the expected Benincasa-Dowker action in 3+1
=======================================================================================

Why this test exists (28 Oct 2026, author's decision: route F6.8, values must come out of the
physics, not be put in)
-------------------------------------------------------------------------------------------
docs/derivations/pravidlo_z_akce.md derives the substrate rule of F7.4a from two inputs: the
proper time of a flow counted by the order (A6) and an action of the order whose static quadratic
form fixes gamma = -c_PhiPsi / (2 c_PsiPsi). For the Einstein-Hilbert (EH) action gamma = 1 in
3+1. The candidate action of the order is the smeared Benincasa-Dowker (BD) action
S = -(2/rho) sum_x (B_eps 1)(x), whose expectation is known to tend to the EH action in the
continuum limit. Whether the DISCRETE action really carries the EH stiffness at second order in
a static perturbation, with no EH number put into the computation, is what F6.8 measures.

F6.8a is the first stage: the conformal direction g = e^(2 a f) eta (Phi = a f, Psi = -a f).
A conformal factor does not change the causal order, so the expected action is an exact
quadrature (core.operators.bd_expect_curved). The non-conformal sectors (F6.8b) move the light
cones and need more geometry.

What is computed (nothing from EH enters)
-----------------------------------------
Expected action per unit time and unit volume, static field, in the bulk of a region of age T
(the past of every point x is cut at the time T before x):

    E[S] = -2 int sqrt(-g) E[B_g 1](x) d^3x dt,
    E[B_g 1](x) = -a_4 (eps rho)^(1/2) [ 1 + rho eps int_{J^-(x)} sqrt(-g(y)) gk(eps rho V_g(x, y)) d^4y ],

gk(w) = sum_i C_i w^i / i! e^-w with the BD coefficients C = (-1, 9, -16, 8), V_g = the proper
volume of the interval (flat diamond times the diamond average of e^(4 a f)). The second
derivative in a at a = 0 is taken ANALYTICALLY under the integral (P = e^(4af), <.> = diamond
average, g1 = w gk', g2 = w^2 gk''):

    d^2/da^2 [P_x P_y gk(w <P>)] = 16 [ (f_x + f_y)^2 gk + 2 (f_x + f_y) <f> g1 + <f>^2 g2 + <f^2> g1 ]

The field is a static plane wave f = cos(k z1), averaged over its phase at x (exactly, in closed
form), so no wall and no slowly converging per-point tails enter (a localized Gaussian has per-
point tails that make the radial sum converge slowly: development, below). The diamond averages
of a plane wave are the Fourier transform of the 3+1 causal diamond, used in closed form
(diamond_ft4, exact; checked against the spheroid quadrature of the core). Units: eps rho = 1,
so l_k = (eps rho)^(-1/4) = 1; the expectation depends on eps and rho only through eps rho.

The measured quantity is the ratio to the EH value of the same second-order term,

    K_hat(k, T) = [d^2/da^2 of the action density] / [EH: 2 (D-1)(D-2) <|grad f|^2> = 6 k^2],

and the conformal stiffness of the BD action is c_conf = 6 K_hat. The EH number 6 enters ONLY as
this normalisation of the comparison. The code computes the action from the kernel and the
geometry alone; the null model and the instruments below check that.

Development (27-28 Oct 2026; numbers on record, not used for the criterion)
---------------------------------------------------------------------------
- Algebra: the analytic second derivative agrees with finite differences of the core operator
  (curved_radial_rows, Gaussian field sigma = 3, a = 0.02) to 3e-4 relative (the O(a^2) error of
  the difference).
- A localized Gaussian has per-point tails of the second-order density: the radial sum moved from
  K = -1.14 to -1.31 between r <= 5 sigma and 7 sigma (2+1, sigma = 3). Hence the plane wave.
- Quadrature: the past of x must be integrated in three pieces (the apex region, the light-cone
  region with panels ~1/k long, and the layer at the past boundary in the variable s itself);
  without the layer the null J0 failed by O(1) at T >= 40. With it J0 ~ 1e-9 at every T.
- Result of the development scan (level 1, checked at level 2 for T <= 40):
    K_hat(k, T) = K_loc(k) - c_T k^3 T, c_T = 0.10 (linear in T to 0.5 %, k^3 to 2 %),
    K_loc: k = 0.025 0.990; 0.05 0.968; 0.1 0.899; 0.2 0.703; 0.3 0.470; 0.5 -0.015,
    1 - K_loc ~ k^2 (3.9 ln(1/k) + 1.1); fit K0 = 1.001.
  The age term was predicted before it was seen in the numbers (thin long intervals running
  perpendicular to k stay in phase with the field over their whole length; the BD moment
  conditions kill their contribution up to tau^4 but not tau^6, giving ~ k^5 per unit length of
  history, i.e. ~ k^3 T relative to k^2). Its asymptotic coefficient follows from a 2D integral
  over (phase, tau) of the exact diamond transform (age_coefficient_asymptotic): 0.102 for
  k -> 0. The full quadrature gives 0.0986-0.1004.

Criteria (stated before the full run)
-------------------------------------
Rows: for every k of the grid three ages T = T_a, 2 T_a, 4 T_a with T_a = max(T_min, n_kT / k)
(k T >= n_kT: the age law holds once the past holds several wavelengths). Level 1 quadrature for
all rows, level 2 for the middle age of three k.

Instrument (INSTRUMENT_FAILED if any fails):
  I1 closed-form diamond transform against the core spheroid quadrature
     (diamond_phi_spheroid at r_x = 1e9 with Pi^4 = 1 + cos(q (r - r_x))/2, so that the radial
     field is a plane wave over the diamond; at r_x = 1e7 the curvature of the reference sphere,
     o_perp^2 / 2 r_x, limits the comparison to 7e-6 for a 40 l_k sliver): max |diff| <= 1e-6;
  I2 analytic second derivative against finite differences of the core curved operator
     (Gaussian sigma = 3, rows r = 0 and 2, a = 0.02): relative difference <= 2e-3;
  I3 phase average: (J(0) + J(pi/2k))/2 from the pointwise formula equals the averaged J
     (k = 0.2, T = 20): relative difference <= 1e-8;
  I4 null (constant conformal factor, k = 0, same quadrature): |J0| <= 1e-3 |J| in every row;
  I5 quadrature: |K_hat(level 2) - K_hat(level 1)| <= 2e-3.
Reading:
  per k: K_hat(k, T) = K_loc(k) - s(k) T fitted over its three ages (linearity residual reported);
  K_loc(k) = K0 - k^2 (c1 ln(1/k) + c0) fitted over k <= fit_kmax (rms residual reported).
Status:
  PASSED     |K0 - 1| <= 0.01 and rms residual <= 2e-3: the EH conformal stiffness 6 comes out
             of the discrete BD action at long wavelengths;
  FAILED     |K0 - 1| > 0.03;
  UNDECIDED  otherwise (incl. a poor fit).
Reported, not judged: c1, c0; c_T(k) = s(k)/k^3 against the asymptotic formula; the wavenumber
where K_loc changes sign; c_conf = 6 K0.

What a pass means and what it does not
--------------------------------------
A pass says that the expected discrete BD action, in its local part, has the EH conformal
stiffness, with a non-analytic finite-l_k correction ~ k^2 ln k. With the diffeomorphism
invariance of the expected action (exact: it is defined on the order), a LOCAL two-derivative
limit can only be alpha R + beta: beta = 0 (the null I4), alpha = 1 (this test), so every static
sector is EH and gamma = 1. That step assumes locality in the other sectors; F6.8b (the mixed
sector Phi x (Phi + Psi), which needs only first-order light cones) tests it independently.
The age term is NOT local: it grows with the duration of the static field. It is a property of
the smeared BD action at finite l_k, of relative size 0.1 (k l_k)^3 (T / l_k).
"""

from __future__ import annotations

import math
import time

import numpy as np

from core.operators.bd import BD_COEFFS, BD_NORM
from core.operators.bd_expect import _directions, _wall_distance
from core.operators.bd_expect_curved import _s_edges_extended, curved_radial_rows, diamond_phi_spheroid
from core.order.sprinkle import interval_volume
from core.stats.parallel import parallel_map
from core.testing import Result

TEST = {
    "id": "F6.8a",
    "phase": "F6",
    "title": "Static conformal quadratic form of the expected BD action in 3+1 (does the EH stiffness 6 emerge?)",
    "theory_part": "derivation of the substrate rule from the order's action (pravidlo_z_akce.md, step 4); "
                   "first stage of F6.8",
    "question": "Does the expected smeared Benincasa-Dowker action of a causal set sprinkled into "
                "g = e^(2af) eta (static f), at second order in a, carry the Einstein-Hilbert conformal stiffness "
                "(D-1)(D-2) = 6 at long wavelengths, with nothing from EH put into the computation?",
    "method": "deterministic quadrature in 3+1, eps rho = 1 (l_k = 1); E[S] = -2 int sqrt(-g) E[B_g 1]; second "
              "derivative in a at a = 0 taken analytically from the kernel gk and the proper interval volumes; "
              "static plane wave f = cos(k z1), phase-averaged; exact closed-form Fourier transform of the causal "
              "diamond; past of x cut at age T; per k three ages with k T >= 8; K_hat = ratio to the EH value",
    "prediction": "EH: K_hat -> 1 as k l_k -> 0 (c_conf = 6). Finite-l_k corrections and the age dependence are "
                  "reported (development: 1 - K_loc ~ k^2 (3.9 ln(1/k) + 1.1), age term -0.10 k^3 T)",
    "criterion": "INSTRUMENT: I1 closed-form diamond transform vs core spheroid quadrature (<= 1e-6); I2 analytic "
                 "second derivative vs finite differences of the core curved operator (<= 2e-3 rel.); I3 phase "
                 "average (<= 1e-8 rel.); I4 null k = 0 (|J0| <= 1e-3 |J|); I5 level 2 vs 1 (|dK| <= 2e-3). "
                 "Fit K_hat = K_loc(k) - s(k) T per k, K_loc = K0 - k^2 (c1 ln(1/k) + c0) for k <= 0.1. "
                 "PASSED |K0 - 1| <= 0.01 and rms <= 2e-3; FAILED |K0 - 1| > 0.03; UNDECIDED otherwise.",
    "null_model": "constant conformal factor (k = 0): the bulk second-order response vanishes (J0)",
    "level": "chain",
    "meters": [],
    "depends_on": ["K.11", "K.11c"],
    "params": {"ks": [0.025, 0.035, 0.05, 0.07, 0.1, 0.14, 0.2, 0.3, 0.5], "T_min": 40.0, "n_kT": 8.0,
               "ages": [1, 2, 4], "fit_kmax": 0.1, "check_ks": [0.025, 0.1, 0.3],
               "tol_ft": 1e-6, "tol_fd": 2e-3, "tol_avg": 1e-8, "tol_null": 1e-3, "tol_level": 2e-3,
               "tol_pass": 0.01, "tol_fail": 0.03, "tol_rms": 2e-3},
    "parts": {"main": {}},
    "quick": {"ks": [0.05, 0.1, 0.2], "T_min": 20.0, "n_kT": 2.0, "ages": [1, 2], "check_ks": [0.2]},
    "estimate": "~30 CPU-min, ~5-10 min with 8 processes (27 rows of level 1 at 10-90 s, 3 rows of level 2 at "
                "2-6 min, instruments and the asymptotic age coefficients ~3 min); memory < 0.5 GB per process",
}

D = 4                                    # spacetime dimension (the closed-form diamond transform is 3+1)
_L2 = interval_volume(D, 1.0) ** (2.0 / D)   # s = l2 tau^2 with eps rho = 1 (w = s^(D/2) = eps rho V_flat)
_A4 = BD_NORM[D]                         # a_4 of the operator normalisation
_EH = (D - 1) * (D - 2)                  # 6: used ONLY to normalise the reported ratio K_hat


# ============================================================================================ kernel

def kernel_derivs(w):
    """gk(w), w gk'(w), w^2 gk''(w) for gk(w) = e^-w sum_i C_i w^i / i! (BD coefficients of 3+1)."""
    p = np.polynomial.Polynomial([c / math.factorial(i) for i, c in enumerate(BD_COEFFS[D])])
    p1, p2 = p.deriv(1), p.deriv(2)
    e = np.exp(-w)
    return e * p(w), w * e * (p1(w) - p(w)), w * w * e * (p2(w) - 2.0 * p1(w) + p(w))


def _gl(a, b, n):
    """Gauss-Legendre nodes and weights on [a, b]."""
    x, w = np.polynomial.legendre.leggauss(int(n))
    return 0.5 * (b - a) * (x + 1.0) + a, 0.5 * (b - a) * w


# ================================================================= Fourier transform of the diamond

def _sinc(x):
    """sin(x)/x, stable at x -> 0."""
    out = np.ones_like(x)
    big = np.abs(x) > 1e-4
    out[big] = np.sin(x[big]) / x[big]
    out[~big] = 1.0 - x[~big] ** 2 / 6.0
    return out


def _dsinc(x):
    """d/dx (sin(x)/x), stable at x -> 0."""
    out = np.empty_like(x)
    big = np.abs(x) > 1e-3
    xb = x[big]
    out[big] = (xb * np.cos(xb) - np.sin(xb)) / (xb * xb)
    xs = x[~big]
    out[~big] = -xs / 3.0 + xs ** 3 / 30.0
    return out


def diamond_ft4(omega, q, h):
    """
    <exp(i (omega t + q.x))> over the 3+1 causal diamond |x| <= h - |t| in its rest frame (centred),
    q = |q|. With B(q, u) = 4 pi (sin qu - qu cos qu)/q^3 (Fourier transform of a ball) and
    V = 2 pi h^4 / 3:  F = (2/V) int_0^h cos(omega (h - u)) B(q, u) du, in closed form:
    cos(omega(h-u)) sin(qu) and u cos(omega(h-u)) cos(qu) split into sin/cos of alpha u + phi with
    (alpha, phi) = (q - omega, omega h), (q + omega, -omega h); int_0^h sin(alpha u + phi) du =
    h sin(phi + alpha h/2) sinc(alpha h/2), and int_0^h u cos(...) du is its alpha-derivative.
    Near-null momenta (q - omega -> 0) stay stable through sinc. q h < 0.05 uses the series of B.
    """
    omega, q, h = np.broadcast_arrays(np.asarray(omega, float), np.asarray(q, float), np.asarray(h, float))
    tot = np.zeros(omega.shape)
    for alpha, phi in ((q - omega, omega * h), (q + omega, -omega * h)):
        x = 0.5 * alpha * h
        A = h * np.sin(phi + x) * _sinc(x)
        dA = h * (0.5 * h) * (np.cos(phi + x) * _sinc(x) + np.sin(phi + x) * _dsinc(x))
        tot = tot + 0.5 * (A - q * dA)
    V = 2.0 * math.pi * h ** 4 / 3.0
    small = q * h < 0.05
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        F = (2.0 / V) * 4.0 * math.pi * tot / np.where(small, 1.0, q) ** 3
    if small.any():
        uu, wu = np.polynomial.legendre.leggauss(40)
        for i in np.flatnonzero(small.ravel()):
            hh, om, qq = h.flat[i], omega.flat[i], q.flat[i]
            u = 0.5 * hh * (uu + 1.0)
            qu = qq * u
            B = (4.0 * math.pi / 3.0) * u ** 3 * (1.0 - qu ** 2 / 10.0 + qu ** 4 / 280.0)
            F.flat[i] = (2.0 / V.flat[i]) * float((0.5 * hh * wu * np.cos(om * (hh - u)) * B).sum())
    return F


def mean_exp_diamond(k, RS, CO, dt):
    """
    <exp(i k o1)> over the lab-frame diamond from x to y (spatial offset RS along n, cos(n, e1) = CO,
    time difference dt), o1 = e1 component of the offset from x. The diamond centre is at RS/2 along
    n; in the rest frame (gamma v = RS/tau) the static wave vector becomes omega' = (RS/tau) k CO,
    |q'| = sqrt(k^2 + omega'^2), and the average is exp(i k CO RS/2) F(omega', |q'|, tau/2).
    """
    tau = np.maximum(np.sqrt(np.maximum(dt * dt - RS * RS, 0.0)), 1e-6)
    kpar = k * CO
    omega = (RS / tau) * kpar
    return np.exp(1j * kpar * RS / 2.0) * diamond_ft4(omega, np.sqrt(k * k + omega * omega), tau / 2.0)


# ======================================================================== quadrature of the past

def s_nodes(s_max, n_s):
    """Panels in s = l2 tau^2 (the flat panels of the core) cut at s_max per rs node: (s, ws), shape (P, Q)."""
    edges = _s_edges_extended(D, 1.0)
    u, wu = _gl(0.0, 1.0, n_s)
    S, WS = [], []
    for lo, hi_e in zip(edges[:-1], edges[1:]):
        hi = np.minimum(hi_e, s_max)[:, None]
        span = np.maximum(hi - lo, 0.0)
        S.append(lo + span * u)
        WS.append(span * wu)
    return np.concatenate(S, axis=1), np.concatenate(WS, axis=1)


def rs_nodes(T, k, mult=1.0, s_layer=16.0):
    """
    Spatial distances rs from x to y in [0, T] (past cut at age T) and weights, in three pieces:
    the apex region [0, 4] (rs = 4 z^2, Gauss in z); the light-cone region up to the layer, in
    panels of length <= min(1/k, 8) (8 Gauss nodes each, resolves exp(i k rs cos)); and the layer
    at the past boundary, where s_max = l2 (T^2 - rs^2) cuts the kernel range: there the variable
    is S = s_max itself (rs = sqrt(T^2 - S/l2)), in the panels of s up to s_layer (w = 256).
    Without the layer the moment cancellation near the past boundary is not resolved (development).
    """
    R0 = min(4.0, T)
    out, wts = [], []
    z, wz = _gl(0.0, 1.0, int(48 * mult))
    out.append(R0 * z ** 2)
    wts.append(2.0 * R0 * z * wz)
    s_lay = min(s_layer, _L2 * (T * T - R0 * R0))
    r_lay = math.sqrt(T * T - s_lay / _L2)
    if r_lay > R0:
        L = r_lay - R0
        npan = max(4, int(math.ceil(L * max(k, 0.05) * mult)), int(math.ceil(L / 8.0)))
        pe = np.linspace(R0, r_lay, npan + 1)
        for a, b in zip(pe[:-1], pe[1:]):
            x, w = _gl(a, b, int(8 * mult))
            out.append(x)
            wts.append(w)
    e = _s_edges_extended(D, 1.0)
    e = np.r_[0.0, e[(e > 1e-3) & (e < s_lay)], s_lay]
    for a, b in zip(e[:-1], e[1:]):
        S, w = _gl(a, b, int(8 * mult))
        r = np.sqrt(T * T - S / _L2)
        out.append(r)
        wts.append(w / (2.0 * _L2 * r))
    return np.concatenate(out), np.concatenate(wts)


def angle_nodes(k, T, mult=1.0):
    """cos(n, e1) in [-1, 1]: Gauss panels of 8 nodes, ~1.5 mult nodes per radian of the phase k T cos."""
    n_pan = max(4, int(math.ceil(2.0 * k * T * 1.5 * mult / 8.0)))
    pe = np.linspace(-1.0, 1.0, n_pan + 1)
    cs, wc = [], []
    for a, b in zip(pe[:-1], pe[1:]):
        x, w = _gl(a, b, 8)
        cs.append(x)
        wc.append(w)
    return np.concatenate(cs), 2.0 * math.pi * np.concatenate(wc)


def _past(T, k, mult):
    """Shared nodes of the past of x: RS (P,1), s, dt, the measure without the direction weight, kernels."""
    rs, wr = rs_nodes(T, k, mult)
    s, ws = s_nodes(_L2 * np.maximum(T * T - rs * rs, 0.0), int(16 * mult))
    RS = rs[:, None]
    dt = np.sqrt(RS * RS + s / _L2)
    base = wr[:, None] * ws * RS ** 2 / (_L2 * 2.0 * dt)        # d^3y dt_y = rs^2 drs dOmega ds / (2 l2 dt)
    g0, g1, g2 = kernel_derivs(s ** (D / 2.0))
    return RS, dt, base, g0, g1, g2


# ================================================================================ the rows

def plane_wave_row(k, T, mult=1.0):
    """
    Phase-averaged second-order density for f = cos(k z1) at age T. Averages over the phase k xi of x:
    <f_x^2> = 1/2, <(f_x + f_y)^2> = 1 + cos(k delta), 2 <(f_x + f_y) <f>> = Re[(1 + e^{ik delta}) conj M_k],
    <<f>^2> = |M_k|^2 / 2, <<f^2>> = 1/2 (delta = y1 - x1, M_k the diamond average of e^{ik o1}).
    J = 1/2 + rho eps int [...] and K_hat = a_4 16 J / (6 k^2 / 2). J0: the same with k = 0 (null).
    """
    t0 = time.time()
    RS, dt, base, g0, g1, g2 = _past(T, k, mult)
    cos, wa = angle_nodes(k, T, mult)
    tot = 0.0
    for a in range(cos.size):
        CO = float(cos[a])
        Mk = mean_exp_diamond(k, RS, CO, dt)
        ph = np.exp(1j * k * RS * CO)
        integ = (1.0 + ph.real) * g0 + np.real((1.0 + ph) * np.conj(Mk)) * g1 + 0.5 * np.abs(Mk) ** 2 * g2 + 0.5 * g1
        tot += wa[a] * float((base * integ).sum())
    J = 0.5 + tot
    J0 = 0.5 + 4.0 * math.pi * float((base * (2.0 * g0 + 2.5 * g1 + 0.5 * g2)).sum())
    Khat = _A4 * D * D * J / (_EH * 0.5 * k * k)
    return {"k": k, "T": T, "mult": mult, "Khat": Khat, "J": J, "J0": J0, "n_ang": int(cos.size),
            "n_rs": int(RS.shape[0]), "time": time.time() - t0}


def plane_wave_point(k, T, xi, mult=1.0):
    """
    Pointwise (not phase-averaged) J at x1 = xi, for the instrument I3: f_x = cos(k xi),
    f_y = cos(k (xi + delta)), <f> = Re(e^{ik xi} M_k), <f^2> = 1/2 + Re(e^{2ik xi} M_2k)/2.
    """
    RS, dt, base, g0, g1, g2 = _past(T, k, mult)
    cos, wa = angle_nodes(k, T, mult)
    fx = math.cos(k * xi)
    tot = 0.0
    for a in range(cos.size):
        CO = float(cos[a])
        fy = np.cos(k * (xi + RS * CO))
        m1 = np.real(np.exp(1j * k * xi) * mean_exp_diamond(k, RS, CO, dt))
        m2 = 0.5 + 0.5 * np.real(np.exp(2j * k * xi) * mean_exp_diamond(2.0 * k, RS, CO, dt))
        fs = fx + fy
        tot += wa[a] * float((base * (fs * fs * g0 + 2.0 * fs * m1 * g1 + m1 * m1 * g2 + m2 * g1)).sum())
    return fx * fx + tot


def age_coefficient_asymptotic(k, phimax=400.0, dphi=0.05, n_s=16):
    """
    Independent prediction of the age term (development). For a long interval of length L >> l_k the
    integrand depends only on the phase phi = k L cos(theta) and on tau (the diamond transform is a
    function of phi and k tau); d^3y = (2 pi / k) L dL dphi, dt = tau dtau / L, so every unit of L
    adds the same amount:
        dJ/dT = (2 pi / k) int dphi int tau dtau [2 cos(phi/2) F g1 + F^2 g2 / 2]
    (the tau-independent parts vanish by the BD moment conditions). Returns c_T = -(dK_hat/dT) / k^3.
    """
    e = _s_edges_extended(D, 1.0)
    e = e[e <= 12.0]
    S, W = [], []
    for a, b in zip(e[:-1], e[1:]):
        x, w = _gl(a, b, n_s)
        S.append(x)
        W.append(w)
    s, ws = np.concatenate(S), np.concatenate(W)
    tau = np.sqrt(s / _L2)
    _, g1, g2 = kernel_derivs(s ** (D / 2.0))
    pe = np.linspace(-phimax, phimax, int(2 * phimax / (8 * dphi)) + 1)
    tot = 0.0
    for a, b in zip(pe[:-1], pe[1:]):
        ph, wp = _gl(a, b, 8)
        P, TT = ph[:, None], tau[None, :]
        F = diamond_ft4(P / TT, np.sqrt(k * k + (P / TT) ** 2), TT / 2.0)
        G = 2.0 * np.cos(P / 2.0) * F * g1[None, :] + 0.5 * F * F * g2[None, :]
        tot += float((wp[:, None] * ws[None, :] / (2.0 * _L2) * G).sum())
    dJ = (2.0 * math.pi / k) * tot
    return -(_A4 * D * D * dJ / (_EH * 0.5 * k * k)) / k ** 3


# ============================================================================ instruments

def check_ft(n_mu=12, n_nu=96, n_phi=12):
    """I1: the closed-form transform against the core spheroid quadrature of the diamond (r_x = 1e9)."""
    R0, q = 1e9, 0.7
    RS = np.array([[0.3], [2.0], [10.0], [40.0]])
    dt = np.sqrt(RS ** 2 + np.array([[1.0], [0.5], [2.0], [0.7]]) ** 2)
    tau = np.sqrt(dt * dt - RS * RS)
    worst = 0.0
    for CO in (1.0, 0.6, 0.05, 0.0, -0.4):
        core = diamond_phi_spheroid(lambda r: (1.0 + 0.5 * np.cos(q * (r - R0))) ** (1.0 / D), D, R0, RS, CO, dt,
                                    tau, n_mu, n_nu, n_phi)
        mine = 1.0 + 0.5 * mean_exp_diamond(q, RS, CO, dt).real
        worst = max(worst, float(np.max(np.abs(core - mine))))
    return worst


def _gauss_row_alg(r_x, sigma, t_x, R, quad):
    """Second-order J for a radial Gaussian f, with the core geometry (spheroid diamond averages)."""
    n_rho, n_s, n_ang, n_mu, n_nu, n_phi = quad
    n = D - 1
    f = lambda r: np.exp(-0.5 * (r / sigma) ** 2)
    cos, wa = _directions(n, n_ang)
    z, wz = _gl(0.0, 1.0, n_rho)
    fx = float(f(r_x))
    tot = 0.0
    for a in range(cos.size):
        CO = cos[a]
        rho_lim = min(t_x, float(_wall_distance(r_x, np.array([CO]), R)[0]))
        rs = rho_lim * z ** 2
        wr = wz * 2.0 * rho_lim * z
        s, ws = s_nodes(_L2 * np.maximum(t_x * t_x - rs * rs, 0.0), n_s)
        RS = rs[:, None]
        dt = np.sqrt(RS * RS + s / _L2)
        tau = np.sqrt(s / _L2)
        r_y = np.sqrt(np.maximum(r_x * r_x + RS * RS + 2.0 * r_x * RS * CO, 0.0))
        meas = wa[a] * wr[:, None] * ws * RS ** (n - 1) / (_L2 * 2.0 * dt)
        m1 = diamond_phi_spheroid(lambda r: f(r) ** (1.0 / D), D, r_x, RS, CO, dt, tau, n_mu, n_nu, n_phi)
        m2 = diamond_phi_spheroid(lambda r: f(r) ** (2.0 / D), D, r_x, RS, CO, dt, tau, n_mu, n_nu, n_phi)
        g0, g1, g2 = kernel_derivs(s ** (D / 2.0))
        fs = fx + f(r_y)
        tot += float((meas * (fs * fs * g0 + 2.0 * fs * m1 * g1 + m1 * m1 * g2 + m2 * g1)).sum())
    return fx * fx + tot


def _gauss_row_fd(r_x, sigma, t_x, R, quad, a):
    """The same J from finite differences of the core curved operator E[B_g 1] (geometry e^(a f), spline)."""
    rho, eps = 10.0, 0.1                       # eps rho = 1
    rg = np.linspace(0.0, 1.2 * R + t_x, 6001)
    vals = []
    for aa in (-a, 0.0, a):
        gv = np.exp(aa * np.exp(-0.5 * (rg / sigma) ** 2))
        A = curved_radial_rows(D, rho, eps, t_x, [r_x], R, np.array([0.0, 2.0 * R + t_x]), np.ones(2),
                               *[int(q) for q in quad], geom="spline", geom_nodes=rg, geom_values=gv)
        vals.append(math.exp(D * aa * math.exp(-0.5 * (r_x / sigma) ** 2)) * A[0].sum())
    pref = -_A4 * (eps * rho) ** 0.5
    return (vals[0] + vals[2] - 2.0 * vals[1]) / (a * a) / (pref * D * D)


def check_fd(r_x, sigma=3.0, a=0.02, quad=(160, 8, 16, 4, 12, 4)):
    """I2: analytic second derivative against finite differences of the core operator (relative difference)."""
    R = 6.0 * sigma
    t_x = 3.0 * R
    ja = _gauss_row_alg(r_x, sigma, t_x, R, quad)
    jf = _gauss_row_fd(r_x, sigma, t_x, R, quad, a)
    return {"r": r_x, "alg": ja, "fd": jf, "rel": abs(ja - jf) / abs(jf)}


def check_avg(k=0.2, T=20.0):
    """I3: (J(0) + J(pi/2k))/2 from the pointwise formula against the phase-averaged row."""
    j_avg = plane_wave_row(k, T)["J"]
    j_pts = 0.5 * (plane_wave_point(k, T, 0.0) + plane_wave_point(k, T, 0.5 * math.pi / k))
    return {"avg": j_avg, "pts": j_pts, "rel": abs(j_avg - j_pts) / abs(j_avg)}


def _dispatch(kind, args):
    """Module-level task for the worker processes."""
    if kind == "row":
        return plane_wave_row(*args)
    if kind == "ft":
        return check_ft()
    if kind == "fd":
        return check_fd(*args)
    if kind == "avg":
        return check_avg()
    if kind == "asym":
        return {"k": args[0], "c_T": age_coefficient_asymptotic(args[0])}
    raise ValueError(kind)


# ================================================================================== run

def _ages(k, p):
    """The ages of the rows of one k: T_a, 2 T_a, ... with T_a = max(T_min, n_kT / k)."""
    Ta = max(float(p["T_min"]), float(p["n_kT"]) / k)
    return [Ta * float(m) for m in p["ages"]]


def run(p, ctx):
    t0 = time.time()
    ks = [float(k) for k in p["ks"]]
    jobs = []
    for k in ks:
        for T in _ages(k, p):
            jobs.append(("row", (k, T, 1.0)))
    n_l1 = len(jobs)
    checks = [float(k) for k in p["check_ks"] if float(k) in ks]
    for k in checks:
        ages = _ages(k, p)
        jobs.append(("row", (k, ages[len(ages) // 2], 2.0)))
    jobs += [("ft", ()), ("fd", (0.0,)), ("fd", (2.0,)), ("avg", ())]
    jobs += [("asym", (k,)) for k in ks]
    res = parallel_map(_dispatch, jobs, ctx.workers, ctx.progress,
                       lambda i, r_, d, n: f"[F6.8a] task {d}/{n} ({time.time() - t0:.0f}s)")
    rows = res[:n_l1]
    lev2 = res[n_l1:n_l1 + len(checks)]
    ft_err, fd0, fd2, avg = res[n_l1 + len(checks):n_l1 + len(checks) + 4]
    asym = {r_["k"]: r_["c_T"] for r_ in res[n_l1 + len(checks) + 4:]}

    # ---- instruments
    lines = ["3+1, eps rho = 1 (l_k = 1); f = cos(k z1), phase-averaged; K_hat = second-order BD action density / "
             "EH value 6 k^2; past of x cut at age T",
             f"  I1 diamond transform vs core spheroid quadrature: max |diff| = {ft_err:.2e} (tol {p['tol_ft']:g})"]
    ok1 = ft_err <= float(p["tol_ft"])
    ok2 = True
    for c in (fd0, fd2):
        ok2 &= c["rel"] <= float(p["tol_fd"])
        lines.append(f"  I2 Gaussian row r = {c['r']:g}: analytic J = {c['alg']:+.6e}, finite differences of the core "
                     f"operator {c['fd']:+.6e}, rel. {c['rel']:.1e} (tol {p['tol_fd']:g})")
    ok3 = avg["rel"] <= float(p["tol_avg"])
    lines.append(f"  I3 phase average: averaged J = {avg['avg']:+.9e}, (J(0) + J(pi/2k))/2 = {avg['pts']:+.9e}, "
                 f"rel. {avg['rel']:.1e} (tol {p['tol_avg']:g})")
    null_worst = max(abs(r_["J0"]) / abs(r_["J"]) for r_ in rows + lev2)
    ok4 = null_worst <= float(p["tol_null"])
    lines.append(f"  I4 null k = 0: max |J0/J| over the rows = {null_worst:.1e} (tol {p['tol_null']:g})")
    ok5 = True
    for r2 in lev2:
        r1 = next(r_ for r_ in rows if r_["k"] == r2["k"] and r_["T"] == r2["T"])
        dK = abs(r2["Khat"] - r1["Khat"])
        ok5 &= dK <= float(p["tol_level"])
        lines.append(f"  I5 k = {r2['k']:g}, T = {r2['T']:g}: K_hat level 1 {r1['Khat']:+.6f}, level 2 "
                     f"{r2['Khat']:+.6f}, |diff| = {dK:.1e} (tol {p['tol_level']:g})")
    inst_ok = ok1 and ok2 and ok3 and ok4 and ok5

    # ---- per k: K_hat(T) = K_loc - s T
    lines.append("  k        ages T                 K_hat(T)                                 K_loc      c_T = s/k^3"
                 "   c_T asympt.  lin. resid.")
    Kloc, data_k = {}, {}
    for k in ks:
        rk = sorted([r_ for r_ in rows if r_["k"] == k], key=lambda r_: r_["T"])
        T = np.array([r_["T"] for r_ in rk])
        K = np.array([r_["Khat"] for r_ in rk])
        if len(T) >= 2:
            A = np.c_[np.ones_like(T), -T]
            (kl, sl), *_ = np.linalg.lstsq(A, K, rcond=None)
            resid = float(np.max(np.abs(K - A @ np.array([kl, sl])))) if len(T) > 2 else float("nan")
        else:
            kl, sl, resid = float(K[0]), float("nan"), float("nan")
        Kloc[k] = float(kl)
        cT = float(sl) / k ** 3
        data_k[str(k)] = {"T": T.tolist(), "Khat": K.tolist(), "K_loc": float(kl), "slope": float(sl),
                          "c_T": cT, "c_T_asym": asym[k], "lin_resid": resid}
        lines.append(f"  {k:<8g} {', '.join(f'{x:g}' for x in T):22s} {', '.join(f'{x:+.5f}' for x in K):40s} "
                     f"{kl:+.5f}   {cT:+.4f}      {asym[k]:+.4f}      {resid:.1e}")

    # ---- K_loc(k) -> K0
    kf = np.array([k for k in ks if k <= float(p["fit_kmax"]) + 1e-12])
    status = "UNDECIDED"
    fit = {}
    if kf.size >= 4:
        y = np.array([Kloc[k] for k in kf])
        A = np.c_[np.ones_like(kf), -kf ** 2 * np.log(1.0 / kf), -kf ** 2]
        c, *_ = np.linalg.lstsq(A, y, rcond=None)
        rms = float(np.sqrt(np.mean((y - A @ c) ** 2)))
        K0, c1, c0 = (float(x) for x in c)
        fit = {"K0": K0, "c1": c1, "c0": c0, "rms": rms, "k": kf.tolist()}
        lines.append(f"  fit over k <= {p['fit_kmax']:g}: K_loc = K0 - k^2 (c1 ln(1/k) + c0): K0 = {K0:.5f}, "
                     f"c1 = {c1:.3f}, c0 = {c0:.3f}, rms {rms:.1e}; conformal stiffness c_conf = 6 K0 = {6 * K0:.4f} "
                     f"(EH 6)")
        if abs(K0 - 1.0) <= float(p["tol_pass"]) and rms <= float(p["tol_rms"]):
            status = "PASSED"
        elif abs(K0 - 1.0) > float(p["tol_fail"]):
            status = "FAILED"
    else:
        lines.append(f"  too few k <= {p['fit_kmax']:g} for the fit ({kf.size}); K0 not determined")
    # sign change of K_loc (descriptive)
    kk = np.array(ks)
    kv = np.array([Kloc[k] for k in ks])
    sc = [float(kk[i] + (kk[i + 1] - kk[i]) * kv[i] / (kv[i] - kv[i + 1])) for i in range(len(ks) - 1)
          if kv[i] > 0 >= kv[i + 1]]
    if sc:
        lines.append(f"  K_loc changes sign at k ~ {sc[0]:.3f} / l_k (wavelength ~ {2 * math.pi / sc[0]:.1f} l_k)")
    lines.append("  age term (not local): K_hat(k, T) = K_loc(k) - c_T k^3 T; relative size 0.1 (k l_k)^3 (T / l_k)")
    if not inst_ok:
        status = "INSTRUMENT_FAILED"
    data = {"rows": rows, "level2": lev2, "ft_err": ft_err, "fd": [fd0, fd2], "avg": avg, "per_k": data_k, "fit": fit,
            "c_T_asym": {str(k): v for k, v in asym.items()}}
    return Result(status, lines, data)
