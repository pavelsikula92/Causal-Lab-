"""
F6.8c -- Optical (nu^2, spin-2) static sector of the expected BD action in 3+1
==============================================================================

Why this test exists (29 Oct 2026, after F6.8b2 PASSED; docs/derivations/path_to_predictions.md)
---------------------------------------------------------------------------------------------
Write a static metric as g = Omega^2 [-dt^2 + n^2 dx^2] with omega = ln Omega (conformal factor: only the
density of the sprinkling) and nu = ln n (optical index: the light cones, i.e. the ORDER). The
Einstein-Hilbert (EH) static quadratic form is, per k^2/2 for plane waves omega = a cos kz, nu = b cos kz,

    A a^2 + B a b + C b^2   with   A = 6, B = 8, C = 2.

Matter couples to omega only. Eliminating nu (nu = -B omega / 2C) gives, WITHOUT assuming covariance,

    gamma = B / (2 C) - 1,        Phi / Phi_GR = -2 / (A - B^2 / 4C).

F6.8a computed A = 6 K_hat (K0 = 1.00005), F6.8b2 computed B = 8 K_hat_ab (K0_ab = 0.9979). The last
coefficient C is the nu^2 sector: the light cones to SECOND order (bending of rays). With it the
field of a point mass near l_k follows (F8.P1), and gamma = 1 at large scales follows from the
action of the order if C -> 2. This test computes C.

What is computed (nothing from EH or from covariance enters)
------------------------------------------------------------
d^2/db^2 at b = 0 of the expected action per unit time and volume, E[S] = -2 <sqrt(-g) E[B_g 1]>, for
g = -dt^2 + e^(2bh) dx^2, h = cos(k z1), phase-averaged; eps rho = 1 (l_k = 1); past of x cut at age T
(dt <= T, as F6.8a and F6.8b2). With P = sqrt(-g) = e^(3bh):

  * the order: y precedes x iff dt >= T(x, y), the optical distance, now to second order. Fermat: the
    ray from y to x minimises int e^(bh) dl with fixed ends; with T = L + b T1 + b^2 T2,
        T1 = int h dl,   T2 = 1/2 int h^2 dl - 1/2 int int grad_perp h(l) G(l, l') grad_perp h(l') dl dl',
    G(u, u') = u_< (d - u_>) / d the Green function of -d^2/dl^2 with fixed ends (the bending of the ray
    lowers the optical distance). Phase-averaged, for a segment of length d with phase advance p = k.e d:
        <T2> = d/4 - 1/4 (k^2 d^2 - p^2) d J(p),   J(x) = (x^2 + 2 cos x - 2) / x^4  (J(0) = 1/12);
  * the past of x is parametrised by s = l2 (dt^2 - T(x,y)^2) (fixed cone): dt(b) = sqrt(T(x,y)^2 + s/l2),
    so dt' = L T1/dt, dt'' = (2 L T2 + T1^2 tau^2/dt^2)/dt; the measure m = 1/(2 l2 dt) gives
    (ln m)' = -L T1/dt^2, (ln m)'' = -2 L T2/dt^2 - T1^2 (dt^2 - 2 L^2)/dt^4;
  * the interval volume V = int d^3z P(z) [dt - T(y,z) - T(z,x)]_+ :
        V'' = int [9 h^2 l0 + 6 h l' + l''] theta(l0) + int delta(l0) l'^2,
    l' = dt' - T1_yz - T1_zx (as F6.8b2), l'' = dt'' - 2 T2_yz - 2 T2_zx; the delta term is the surface of
    the spatial projection of the diamond (a prolate spheroid). The differences of near-equal path
    lengths (l', and the h^2 part of l'') are formed in a form that does not cancel O(L) numbers:
    l' = E_xy (L^2/dt - d1 - d2) + d1 (E_xy - E_1) + d2 (E_xy - E_2) with E the segment averages of
    e^(ik.r), and L^2/dt - d1 - d2 = -L [(sigma - 1) + (dt - L)/dt] from tau^2 directly;
  * d^2/db^2 of P_x P_y gk(V) m expanded analytically with g1 = w gk', g2 = w^2 gk'';
  * the age cut s <= s_max = l2 (T^2 - T(x,y)^2) depends on b: boundary terms
        2 F_b(s_max) s_max' + F_s(s_max) s_max'^2 + F(s_max) s_max'',
    integrated over the layer directly in S = s_max on their own panels (as F6.8b2).
The EH number 2 enters ONLY as the normalisation of the reported ratio K_hat_bb = (d^2 E[S]/db^2)/(2 k^2),
so that C = 2 K_hat_bb. The conformal K_hat of F6.8a is computed on the same nodes in the same pass.

Development (29 Oct 2026; numbers on record, not used for the criterion)
------------------------------------------------------------------------
- <T2> closed form against numerically minimised rays (k = 0.9, d = 2-5): 4e-7 to 8e-5.
- d2V/db2 of the test diamond (k = 0.7, CO = 0.4, L = 2, tau = 1.2): -0.0159484 (prolate, delta term as a
  surface) against -0.0159900 / -0.0159590 (numerical roots + finite differences, 24 / 48 nodes): converging.
- The edge term (pairs leaving the past where S(b) < 0) was missing in the first version: the boundary
  was ~2 pi T^4. With it, combined by parts, the boundary falls as T^-4 (k = 0: -4.0e-5 / -2.5e-6 / -1.6e-7 /
  -9.7e-9 / -3.9e-10 per unit solid angle at T = 20 / 40 / 80 / 160 / 320; 32 and 64 nodes agree).
- Null k = 0: J_bb = +8.2e-5 at T = 20, ~ +5e-6 at T = 40 (bulk and boundary of opposite sign).
- First row: k = 0.2, T = 20 (boundary regime): K_hat_bb = 0.5500, K_hat = 0.6603.

Reading
-------
Per k, K_hat_bb(k, T) = K_loc_bb(k) - s_bb(k) T over three ages (the grid of F6.8b2, so that gamma(k, T)
can be formed row by row with the recorded K_hat_ab of F6.8b2).

Criteria (stated before the full run)
-------------------------------------
Rows: for every k of the grid three ages T = T_a, 1.5 T_a, 2 T_a with T_a = max(T_min, n_kT / k)
(k = 0.05: T = 160, 240, 320), as F6.8b2.
Instrument (INSTRUMENT_FAILED if any fails):
  I1 geometry:
     (a) the closed form <T2> against the optical length of a ray found by direct numerical
         minimisation (polyline, finite differences in b), three segments: relative difference <= 2e-3;
     (b) the phase-averaged d2V/db2 of the test diamond (prolate quadrature, delta term as a surface)
         against an independent computation: V(b) per phase, integrated on the rays sigma of prolate
         coordinates up to the root of l(b) = 0 found numerically on every ray, the bending of every
         segment by numerical quadrature of the Green function, second finite difference in b (no delta
         term, no closed forms of the optics): relative difference <= 2e-3;
     (c) the uniform null: a constant index (k = 0) is a rescaling of space at fixed optical proper time,
         so dV/db = d2V/db2 = 0 exactly; four diamonds from apex to sliver: |.| <= 1e-9;
  I2 the conformal K_hat on this grid against F6.8a's own row (k = 0.2, T = 40): |diff| <= 1e-3;
  I3 null k = 0 at T = 40, 160, 320: the null is the finite-age residual of flat space (a constant index
     rescales space only, the cut dt <= T stays), and falls as T^-4; each |J_bb(k = 0, T)| is expressed in
     K_hat_bb units at the smallest k whose rows reach down to that age (k = 0.2 at T = 40, 0.05 at 160 and
     320; a_4 |J| / k^2) and must be <= 2e-3;
  I4 quadrature: |K_hat_bb(level 2) - K_hat_bb(level 1)| <= 2e-3 on the check rows.
Status (fit K_loc_bb = K0_bb - k^2 (c1 ln(1/k) + c0) over k <= 0.14):
  PASSED     |K0_bb - 1| <= 0.01: the third coefficient of the EH form emerges, C(k -> 0) = 2; with
             F6.8a and F6.8b2, gamma = B/(2C) - 1 -> 1 at large scales follows from the action of the order;
  FAILED     |K0_bb - 1| > 0.03: C(0) differs from 2, gamma = 1 does NOT follow from the action of the order;
  UNDECIDED  otherwise.
Reported, not judged: per k and T, gamma(k, T) = 2 K_hat_ab / K_hat_bb - 1 and Phi/Phi_GR =
-2 / (6 K_hat - 8 K_hat_ab^2 / K_hat_bb) with K_hat_ab from the record of F6.8b2 (same k, same T);
the spin-2 form factor that covariance would assign, phi2 = 2 K_loc_bb - 2 - (8/3)(K_loc - 1); the age
coefficient c_T,bb = s_bb / k^3 against c_T.
"""

from __future__ import annotations

import math
import time

import numpy as np

from core.stats.parallel import parallel_map
from core.testing import Result
from phases.F6_order_dynamics.F6_08a_bd_action_conformal import (_A4, _L2, _dsinc, _gl, _sinc, angle_nodes,
                                                                  kernel_derivs, mean_exp_diamond, plane_wave_row,
                                                                  rs_nodes)
from phases.F6_order_dynamics.F6_08b2_bd_action_mixed_rev import _ediv, _s_nodes

TEST = {
    "id": "F6.8c",
    "phase": "F6",
    "title": "Optical (nu^2, spin-2) static sector of the expected BD action: does C -> 2 (EH) emerge?",
    "theory_part": "derivation of the substrate rule from the order's action (rule_from_action.md, path_to_predictions.md); "
                   "third stage of F6.8; input of F8.P1",
    "question": "Does the nu^2 static coefficient C of the expected discrete BD action (light cones to second order, "
                "bending of rays) tend to the Einstein-Hilbert value 2 at long wavelength, so that gamma = B/(2C) - 1 = 1 "
                "follows from the action of the order, with nothing from EH put into the computation?",
    "method": "deterministic quadrature in 3+1, eps rho = 1 (l_k = 1); g = -dt^2 + e^(2bh)dx^2, h = cos(k z1), "
              "phase-averaged; order through the optical distance to second order (Fermat: bent ray, Green function "
              "with fixed ends, closed form); interval volume with the delta (surface) term; d^2/db^2 analytic; past of "
              "x cut at age T with its boundary terms; the grid of F6.8b2 (k T = 8, 12, 16); K_hat_bb = ratio to the "
              "EH value 2 (2 k^2); the conformal K_hat (F6.8a) on the same nodes",
    "prediction": "EH: K_hat_bb -> 1 as k -> 0 (C -> 2), hence gamma -> 1",
    "criterion": "INSTRUMENT: I1 (a) closed-form <T2> vs numerical Fermat ray (<= 2e-3 rel.), (b) d2V/db2 vs independent "
                 "quadrature with numerical roots, numerical Green-function bending and finite differences in b (<= 2e-3 rel.), (c) uniform "
                 "null (<= 1e-9); I2 conformal K_hat vs F6.8a row (<= 1e-3); I3 null k = 0 in K_hat_bb units at the "
                 "smallest k whose rows reach that age (a_4 |J0| / k^2 <= 2e-3); I4 level 2 vs 1 (<= 2e-3). Fit K_loc_bb = K0_bb - k^2 (c1 ln(1/k) "
                 "+ c0) over k <= 0.14: PASSED |K0_bb - 1| <= 0.01; FAILED > 0.03; UNDECIDED otherwise.",
    "null_model": "constant index (k = 0): flat space rescaled at fixed optical proper time, the response vanishes",
    "level": "chain",
    "meters": [],
    "depends_on": ["F6.8a", "F6.8b2"],
    "params": {"ks": [0.05, 0.07, 0.1, 0.14, 0.2, 0.3], "T_min": 40.0, "n_kT": 8.0, "ages": [1, 1.5, 2],
               "fit_kmax": 0.14, "check": [[0.05, 0], [0.1, 0], [0.2, 0]],
               "null_ages": [40.0, 160.0, 320.0], "quad": [4, 16, 0.8, 8, 16], "bnd_nodes": 32,
               "tol_t2": 2e-3, "tol_geom": 2e-3, "tol_unif": 1e-9, "tol_conf": 1e-3, "tol_null": 2e-3,
               "tol_level": 2e-3, "tol_pass": 0.01, "tol_fail": 0.03},
    "parts": {"main": {}},
    "quick": {"ks": [0.1, 0.2], "T_min": 40.0, "n_kT": 1.0, "ages": [1, 2], "check": [[0.2, 0]],
              "null_ages": [40.0], "fit_kmax": 0.2},
    "estimate": "overnight like F6.8b2 (~1-2 h wall with 10 processes; ~1.5x the cost of F6.8b2 per row); "
                "memory < 1 GB per process",
}

# K_hat_ab(k, T) of F6.8b2 (full run 28-29 Oct 2026, results/F6.8b2), the same grid: used ONLY for the
# reported gamma(k, T) and Phi/Phi_GR, never for the status.
_B2_ROWS = {
    (0.05, 160.0): 0.965331271791107, (0.05, 240.0): 0.9633688658929301, (0.05, 320.0): 0.9622327688898147,
    (0.07, 114.29): 0.9385975093690221, (0.07, 171.43): 0.9349952634530874, (0.07, 228.57): 0.9308905434681602,
    (0.1, 80.0): 0.8887727875872902, (0.1, 120.0): 0.8813897657841926, (0.1, 160.0): 0.8731525949103968,
    (0.14, 57.14): 0.8078745091318962, (0.14, 85.71): 0.7934174001307007, (0.14, 114.29): 0.7772865867498508,
    (0.2, 40.0): 0.6632305031666158, (0.2, 60.0): 0.6337918228408992, (0.2, 80.0): 0.6009577355284168,
    (0.3, 40.0): 0.31477003007099524, (0.3, 60.0): 0.19688086163304738, (0.3, 80.0): 0.028657065588118365,
}
_B2_LOC = {0.05: 0.9682920565432223, 0.07: 0.9463882209480492, 0.1: 0.9045353384426327,
           0.14: 0.8387413822438841, 0.2: 0.726069171969276, 0.3: 0.6092720991550357}


def _b2(k, T):
    return _B2_ROWS.get((round(k, 4), round(T, 2)))


# ====================================================================== second-order optics

def jbend(x):
    """J(x) = (x^2 + 2 cos x - 2) / x^4 = int int cos(x (u - u')) G(u, u') du du' on [0, 1] (J(0) = 1/12)."""
    x = np.asarray(x, float)
    out = np.empty_like(x)
    sm = np.abs(x) < 0.5
    xs = x[sm] ** 2
    out[sm] = 1.0 / 12.0 - xs / 360.0 + xs ** 2 / 20160.0 - xs ** 3 / 1814400.0 + xs ** 4 / 239500800.0
    xb = x[~sm]
    out[~sm] = (xb * xb + 2.0 * np.cos(xb) - 2.0) / xb ** 4
    return out


def bend(k, d, p):
    """Phase-averaged bending part of <T2> for a segment of length d with phase advance p (<= 0):
    -1/4 (k^2 d^2 - p^2) d J(p)."""
    return -0.25 * np.maximum(k * k * d * d - p * p, 0.0) * d * jbend(p)


def t2_mean(k, d, p):
    """<T2> = d/4 + bend (phase-averaged second-order optical length of a segment)."""
    return 0.25 * d + bend(k, d, p)


# ====================================================================== geometry of the diamond

def spheroid_pieces2(k, CO, RS, dt, tau2, n_sig, n_te, n_phi):
    """
    Uniform averages over the spatial projection of the diamond x -> y (prolate spheroid d1 + d2 <= dt,
    foci y and x, L = RS, cos(n, e1) = CO); amplitudes relative to x (o = z - x):
        Dm = mean[l'],  Qm = mean[e^{ik.o} conj(l')],  Lm = mean[l''_z]  (l'' without the z-independent
        T1^2 part of dt''),  Sf = (1/2) int delta(l0) |l'|^2 d^3z / volume.
    l' = (L/dt) D_xy - T1_yz - T1_zx in a form free of cancellation of O(L) numbers (see the header).
    Prolate coordinates: sigma in [1, dt/L], te in [-1, 1], phi periodic; d1 = c (sigma + te) to y,
    d2 = c (sigma - te) to x, axial distance from y xi = c (1 + sigma te), weight c^3 (sigma^2 - te^2).
    """
    L = RS[..., None, None]
    dtt = dt[..., None, None]
    t2 = tau2[..., None, None]
    c = 0.5 * L
    Lsafe = np.maximum(L, 1e-300)
    em1 = t2 / (Lsafe * (dtt + L))                     # sigma_max - 1 = (dt - L)/L, from tau^2
    rdt = t2 / (dtt * (dtt + L))                       # (dt - L)/dt
    kpar = k * CO
    kperp = k * math.sqrt(max(1.0 - CO * CO, 0.0))
    Exy = (np.exp(1j * kpar * RS / 2.0) * np.sinc(kpar * RS / 2.0 / math.pi))[..., None, None]
    Dxy = RS * Exy[..., 0, 0]
    Bxy = bend(k, L, kpar * L)
    xs, ws = _gl(0.0, 1.0, n_sig)
    ts, wt = _gl(-1.0, 1.0, n_te)
    te = ts[:, None]
    wte = wt[:, None]
    cph = np.cos(2.0 * math.pi * (np.arange(n_phi) + 0.5) / n_phi)[None, :]
    eL = np.exp(1j * kpar * L)
    num_D = np.zeros(RS.shape, complex)
    num_Q = np.zeros(RS.shape, complex)
    num_L = np.zeros(RS.shape)
    den = np.zeros(RS.shape)

    def at(sm1):
        """Pieces at sigma = 1 + sm1 (sm1 >= 0, given directly to keep the sliver digits)."""
        sig = 1.0 + sm1
        xi = c * (1.0 + sig * te)
        d1 = c * (sig + te)
        d2 = c * (sig - te)
        rp = c * np.sqrt(np.maximum(sm1 * (2.0 + sm1) * (1.0 - te * te), 0.0))
        ko = kpar * (L - xi) + kperp * rp * cph            # k . (z - x)
        kzy = ko - kpar * L                                 # k . (z - y)
        eko = np.exp(1j * ko)
        E1 = _ediv(eko - eL, kzy, eL)
        E2 = _ediv(eko - 1.0, ko, 1.0)
        q = -L * (sm1 + rdt)                                # L^2/dt - d1 - d2
        dlt = Exy * q + d1 * (Exy - E1) + d2 * (Exy - E2)
        return sig, d1, d2, ko, kzy, eko, q, dlt

    for x, wx in zip(xs, ws):
        sm1 = em1 * x
        sig, d1, d2, ko, kzy, eko, q, dlt = at(sm1)
        wb = (em1 * wx) * wte * c ** 3 * (sig * sig - te * te) / n_phi
        l2z = 0.5 * q + (2.0 * L / dtt) * Bxy - 2.0 * bend(k, d1, kzy) - 2.0 * bend(k, d2, ko)
        wbb = np.broadcast_to(wb, dlt.shape)
        num_D += (wbb * dlt).sum(axis=(-2, -1))
        num_Q += (wbb * eko * np.conj(dlt)).sum(axis=(-2, -1))
        num_L += (wbb * l2z).sum(axis=(-2, -1))
        den += wbb.sum(axis=(-2, -1))
    # surface sigma = dt/L: delta(dt - L sigma) gives 1/L; there q = L^2/dt - dt = -tau^2/dt exactly
    sig, d1, d2, ko, kzy, eko, _, _ = at(em1)
    Exy_b = np.broadcast_to(Exy, eko.shape)
    E1 = _ediv(eko - eL, kzy, eL)
    E2 = _ediv(eko - 1.0, ko, 1.0)
    dlt_s = Exy_b * (-t2 / dtt) + d1 * (Exy_b - E1) + d2 * (Exy_b - E2)
    ws_ = np.broadcast_to(wte * c ** 3 * (sig * sig - te * te) / (n_phi * Lsafe), dlt_s.shape)
    surf = (ws_ * 0.5 * np.abs(dlt_s) ** 2).sum(axis=(-2, -1))
    ok = den > 0
    dd = np.where(ok, den, 1.0)
    z = lambda a: np.where(ok, a / dd, 0.0)
    return z(num_D), z(num_Q), z(num_L), z(surf), Dxy, Bxy[..., 0, 0]


def second_order_volume(k, CO, RS, dt, tau2, n_sig, n_te, n_phi):
    """beta = <V'>/V0 amplitude and gamma = <V''>/V0 (phase-averaged) of the interval x -> y at fixed s."""
    M = mean_exp_diamond(k, RS, CO, dt) if k > 0 else np.ones(RS.shape, complex)
    Dm, Qm, Lm, Sf, Dxy, Bxy = spheroid_pieces2(k, CO, RS, dt, tau2, n_sig, n_te, n_phi)
    AV = np.where(tau2 > 0, 4.0 * dt / np.maximum(tau2, 1e-300), 0.0)    # A / V0 = 4 dt / tau^2
    beta = 3.0 * M + AV * Dm
    dD = 0.5 * np.abs(Dxy) ** 2 * tau2 / dt ** 3                          # T1^2 tau^2/dt^3 part of dt''
    gam = 4.5 + 3.0 * AV * np.real(Qm) + AV * (Lm + dD) + AV * Sf
    return M, beta, gam, Dxy, Bxy


# ======================================================================== quadrature of the past

def _pair_terms(k, CO, RS, s, dt, quad):
    """Phase-averaged nu^2 integrand and the conformal (F6.8a) integrand on the same nodes."""
    n_sig, te0, te_per_rad, n_phi = quad
    tau2 = s / _L2
    g0, g1, g2 = kernel_derivs(s * s)
    ph = np.exp(1j * k * CO * RS)
    one = 1.0 + ph
    n_te = int(te0 + math.ceil(te_per_rad * abs(k * CO) * (float(np.max(RS)) if RS.size else 0.0)))
    M, beta, gam, Dxy, Bxy = second_order_volume(k, CO, RS, dt, tau2, n_sig, n_te, n_phi)
    mu = -RS * Dxy / (dt * dt)                                            # (ln m)' amplitude
    t2xy = 0.25 * RS + Bxy
    mup = -2.0 * RS * t2xy / dt ** 2 - 0.5 * np.abs(Dxy) ** 2 * (dt * dt - 2.0 * RS * RS) / dt ** 4
    X1 = 3.0 * one + mu                                                   # (ln X)' amplitude, X = P_x P_y m
    bb = ((0.5 * np.abs(X1) ** 2 + mup) * g0 + np.real(X1 * np.conj(beta)) * g1
          + gam * g1 + 0.5 * np.abs(beta) ** 2 * g2)
    conf = (1.0 + ph.real) * g0 + np.real(one * np.conj(M)) * g1 + 0.5 * np.abs(M) ** 2 * g2 + 0.5 * g1
    return bb, conf, n_te


def _boundary_term(k, CO, T, quad, n_per=32, s_top=16.0):
    """
    Boundary terms of the age cut for one direction (without the direction weight). The past of x is cut at
    dt = T, i.e. s <= S = l2 (T^2 - T_xy^2), and T_xy depends on b:
        d2/db2 int_0^S F ds = int F_bb + 2 F_b(S) S' + F_s(S) S'^2 + F(S) S'',
    S' = -2 l2 L T1_xy, S'' = -2 l2 (T1_xy^2 + 2 L T2_xy). Integrated over the layer directly in S
    (rs = sqrt(T^2 - S/l2), d rs = dS / (2 l2 rs)) on the s panels up to s_top with n_per Gauss nodes each.
    Where S(b) crosses zero the pair leaves the past: the kink of int_0^max(S, 0) gives the edge term
    F(0) S'^2 delta(S). The edge term and the part of F_s S'^2 from dgk/ds are each ~ T^4 and cancel; they
    are combined analytically (by parts) into -(1/(2 l2 T)) int gk(S^2) [W'(S) - W'(0)] dS with
    W = l2 rs^3 <|T1_xy|^2> = l2 rs^5 sinc^2(k_par rs / 2) (int_0^inf gk(S^2) dS = 0, the kernel moment).
    """
    from core.operators.bd_expect_curved import _s_edges_extended
    n_sig, te0, te_per_rad, n_phi = quad
    e = _s_edges_extended(4, 1.0)
    e = np.r_[0.0, e[(e > 1e-3) & (e < s_top)], s_top]
    tot = 0.0
    kpar = k * CO
    dW0 = float(_dW(kpar, np.array(T)))
    for lo, hi in zip(e[:-1], e[1:]):
        S, wS = _gl(lo, hi, n_per)
        r = np.sqrt(np.maximum(T * T - S / _L2, 1e-300))[:, None]
        wr = (wS / (2.0 * _L2 * r[:, 0]))[:, None]
        tau2 = (S / _L2)[:, None]
        dtT = np.full_like(r, T)
        g0, g1, _ = kernel_derivs((S * S)[:, None])
        n_te = int(te0 + math.ceil(te_per_rad * abs(k * CO) * float(np.max(r))))
        M, beta, _, Dxy, Bxy = second_order_volume(k, CO, r, dtT, tau2, n_sig, n_te, n_phi)
        one = 1.0 + np.exp(1j * k * CO * r)
        mu = -r * Dxy / (T * T)
        Fb = ((3.0 * one + mu) * g0 + g1 * beta) / (2.0 * _L2 * T)
        F0 = g0 / (2.0 * _L2 * T)
        Fs_m = -g0 / (4.0 * _L2 * _L2 * T ** 3)                         # d/ds of m = 1/(2 l2 dt)
        Sp = -2.0 * _L2 * r * Dxy
        Spp = -2.0 * _L2 * (0.5 * np.abs(Dxy) ** 2 + 2.0 * r * (0.25 * r + Bxy))
        dens = np.real(Fb * np.conj(Sp)) + Fs_m * 0.5 * np.abs(Sp) ** 2 + F0 * Spp
        tot += float((wr * r ** 2 * dens).sum())
        # d gk/ds part of F_s S'^2 together with the edge term, by parts
        tot += float((-(wS[:, None] / (2.0 * _L2 * T)) * g0 * (_dW(kpar, r) - dW0)).sum())
    return tot


def _dW(kp, r):
    """dW/dS for W = l2 r^5 sinc^2(kp r/2), r = sqrt(T^2 - S/l2): -(5 r^3 si^2 + 2 r^4 si si')/2."""
    x = kp * r / 2.0
    si = _sinc(np.asarray(x, float))
    dsi = 0.5 * kp * _dsinc(np.asarray(x, float))
    return -0.5 * (5.0 * r ** 3 * si * si + 2.0 * r ** 4 * si * dsi)


def optical_row(k, T, mult=1.0, quad=(4, 16, 0.8, 8, 16), bnd_nodes=32):
    """
    Phase-averaged nu^2 second-order density J_bb at age T (and the conformal J of F6.8a on the same nodes).
    J_bb = 9/2 (from P_x) + int d^3y ds/(2 l2 dt) [bb] + boundary terms of the age cut.
    K_hat_bb = a_4 J_bb / k^2 (EH 2 k^2).
    """
    t0 = time.time()
    q = (int(quad[0]), int(quad[1] * mult), float(quad[2]) * mult, int(quad[3]))
    rs, wr = rs_nodes(T, k, mult)
    cos, wa = angle_nodes(k, T, mult)
    n_s = int((quad[4] if len(quad) > 4 else 16) * mult)
    inner = rs < 4.0
    Jbb, Jc, bnd, nte_max = 4.5, 0.5, 0.0, 0
    smx_all = _L2 * np.maximum(T * T - rs * rs, 0.0)
    for a in range(cos.size):
        CO = float(cos[a])
        for sel, far in ((inner, False), (~inner, True)):
            r = rs[sel]
            if r.size == 0:
                continue
            w_r = wr[sel]
            smx = smx_all[sel]
            order = np.argsort(r)
            for chunk in np.array_split(order, max(1, int(math.ceil(r.size / 24)))):
                s, ws = _s_nodes(smx[chunk], n_s, far)
                RS = r[chunk][:, None] * np.ones_like(s)
                dt = np.sqrt(RS * RS + s / _L2)
                base = w_r[chunk][:, None] * ws * RS ** 2 / (_L2 * 2.0 * dt)
                bb, conf, nte = _pair_terms(k, CO, RS, s, dt, q)
                nte_max = max(nte_max, nte)
                Jbb += wa[a] * float((base * bb).sum())
                Jc += wa[a] * float((base * conf).sum())
        bnd += wa[a] * _boundary_term(k, CO, T, q, bnd_nodes)
    Jbb += bnd
    out = {"k": k, "T": T, "mult": mult, "J_bb": Jbb, "J_c": Jc, "bnd": bnd, "n_ang": int(cos.size),
           "n_rs": int(rs.size), "n_te_max": nte_max, "time": time.time() - t0}
    if k > 0:
        out["Kbb"] = _A4 * Jbb / (k * k)
        out["Kc"] = _A4 * 16.0 * Jc / (3.0 * k * k)
    return out


# ============================================================================ instruments

def _fermat_T(A, B, kvec, th, b, n=200):
    """Optical length of the ray A -> B minimising int e^(b cos(k.r + th)) dl (polyline, trapezoid in n)."""
    from scipy.optimize import minimize
    d = B - A
    Ln = float(np.linalg.norm(d))
    e = d / Ln
    tmp = np.array([1.0, 0.0, 0.0]) if abs(e[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    p1 = np.cross(e, tmp)
    p1 /= np.linalg.norm(p1)
    p2 = np.cross(e, p1)
    u = np.linspace(0.0, 1.0, n + 1)
    base = A[None, :] + u[:, None] * d[None, :]

    def path(v):
        w = np.zeros((n + 1, 2))
        w[1:-1] = v.reshape(n - 1, 2)
        return base + w[:, :1] * p1[None, :] + w[:, 1:] * p2[None, :]

    def f(v):
        P = path(v)
        ph = P @ kvec + th
        nn = np.exp(b * np.cos(ph))
        seg = P[1:] - P[:-1]
        ln = np.linalg.norm(seg, axis=1)
        val = float((ln * 0.5 * (nn[1:] + nn[:-1])).sum())
        # gradient
        gn = -b * np.sin(ph)[:, None] * nn[:, None] * kvec[None, :]      # d n / d P
        g = np.zeros_like(P)
        ue = seg / ln[:, None]
        av = 0.5 * (nn[1:] + nn[:-1])
        g[1:] += ue * av[:, None]
        g[:-1] -= ue * av[:, None]
        g[1:] += 0.5 * ln[:, None] * gn[1:]
        g[:-1] += 0.5 * ln[:, None] * gn[:-1]
        gv = np.stack([g[1:-1] @ p1, g[1:-1] @ p2], axis=1)
        return val, gv.ravel()

    r = minimize(f, np.zeros(2 * (n - 1)), jac=True, method="L-BFGS-B",
                 options={"ftol": 1e-15, "gtol": 1e-13, "maxiter": 20000})
    return r.fun


def check_t2():
    """I1 (a): closed-form <T2> against numerically minimised rays (second difference in b, phases 0, pi/2)."""
    k = 0.9
    out = []
    for d, ang in ((3.0, 0.4), (5.0, 1.1), (2.0, 1.5)):
        kvec = k * np.array([math.cos(ang), math.sin(ang), 0.0])
        A = np.array([0.3, -0.2, 0.1])
        B = A + d * np.array([1.0, 0.0, 0.0])
        tt = []
        for th in (0.0, 0.5 * math.pi):
            T0 = _fermat_T(A, B, kvec, th, 0.0)
            sd = {}
            for eps in (0.01, 0.005):                                           # T2 = T''/2 for this phase
                Tp, Tm = (_fermat_T(A, B, kvec, th, bb) for bb in (eps, -eps))
                sd[eps] = (Tp - 2.0 * T0 + Tm) / (eps * eps) / 2.0
            tt.append((4.0 * sd[0.005] - sd[0.01]) / 3.0)                       # Richardson: removes b^2 T4
        num = 0.5 * (tt[0] + tt[1])
        p = k * math.cos(ang) * d
        ana = float(t2_mean(k, np.array(d), np.array(p)))
        out.append({"d": d, "angle": ang, "num": num, "ana": ana, "rel": abs(num - ana) / abs(ana)})
    return out


def _seg_T12(Ps, Pe, kvec, th, n=96):
    """Per-phase T1 and T2 of the straight segments Ps -> Pe (arrays (..., 3)); T2 bending by numerical
    quadrature of the Green function (cumulative trapezoid on n intervals)."""
    d = Pe - Ps
    Ln = np.linalg.norm(d, axis=-1)
    u = np.linspace(0.0, 1.0, n + 1)
    P = Ps[..., None, :] + u[:, None] * d[..., None, :]
    ph = P @ kvec + th
    h = np.cos(ph)
    du = (Ln / n)[..., None]
    trap = lambda f: ((f[..., 1:] + f[..., :-1]) * 0.5 * du).sum(-1)
    T1 = trap(h)
    e = d / np.maximum(Ln, 1e-300)[..., None]
    kpe = (e @ kvec)
    kp2 = np.maximum(kvec @ kvec - kpe ** 2, 0.0)
    gsq = np.sin(ph)                                  # grad_perp h = -k_perp sin(ph) (one direction)
    l = u[None, :] * Ln[..., None] if Ln.ndim else u * Ln
    # xi(l) = ((d - l)/d) int_0^l l' g dl' + (l/d) int_l^d (d - l') g dl'
    def cum(f):
        c = np.zeros_like(f)
        c[..., 1:] = np.cumsum((f[..., 1:] + f[..., :-1]) * 0.5 * du, axis=-1)
        return c
    Lx = Ln[..., None]
    a1 = cum(l * gsq)
    b1 = cum((Lx - l) * gsq)
    xi = ((Lx - l) / np.maximum(Lx, 1e-300)) * a1 + (l / np.maximum(Lx, 1e-300)) * (b1[..., -1:] - b1)
    bendint = trap(gsq * xi)
    T2 = 0.5 * trap(h * h) - 0.5 * kp2 * bendint
    return T1, T2


def check_geometry(nte=48, nph=48, nsig=24, nseg=192):
    """I1 (b): phase-averaged d2V/db2 of a diamond (prolate formula with the delta term as a surface) against
    an independent computation: V(b) integrated on rays sigma at fixed (te, phi) of prolate coordinates up to
    the root of l(b) = 0 found numerically on every ray, per phase (0, pi/2), with the bending of every segment
    by numerical quadrature of the Green function, and a second finite difference in b; I1 (c): the uniform null."""
    k, CO, L, tau = 0.7, 0.4, 2.0, 1.2
    dt0 = math.sqrt(L * L + tau * tau)
    RSa, dta, t2a = np.array([[L]]), np.array([[dt0]]), np.array([[tau * tau]])
    _, beta_a, gam_a, _, _ = second_order_volume(k, CO, RSa, dta, t2a, 16, 128, 32)
    gam_a = float(gam_a[0, 0])
    V0 = math.pi * tau ** 4 / 24.0
    nvec = np.array([CO, math.sqrt(1 - CO * CO), 0.0])
    kvec = np.array([k, 0.0, 0.0])
    y = L * nvec
    x0 = np.zeros(3)
    c = 0.5 * L
    e2 = np.array([-nvec[1], nvec[0], 0.0])
    e3 = np.cross(nvec, e2)
    tes, wts = _gl(-1.0, 1.0, nte)
    phs = 2.0 * math.pi * (np.arange(nph) + 0.5) / nph
    TE, PH = np.meshgrid(tes, phs, indexing="ij")
    W2 = wts[:, None] * np.ones((1, nph)) * (2.0 * math.pi / nph)
    xs_, wxs = _gl(0.0, 1.0, nsig)

    def point(sig, te_, ph_):
        xi = c * (1.0 + sig * te_)                                   # axial distance from y towards x
        rp = c * np.sqrt(np.maximum((sig * sig - 1.0) * (1.0 - te_ * te_), 0.0))
        return (y - xi[..., None] * nvec + rp[..., None] * (np.cos(ph_)[..., None] * e2
                                                            + np.sin(ph_)[..., None] * e3))

    res = []
    for th in (0.0, 0.5 * math.pi):
        T1xy, T2xy = _seg_T12(y, x0, kvec, th, 4 * nseg)

        def ell(Pz, b):
            T1a, T2a = _seg_T12(np.broadcast_to(y, Pz.shape), Pz, kvec, th, nseg)
            T1b, T2b = _seg_T12(Pz, np.broadcast_to(x0, Pz.shape), kvec, th, nseg)
            Txy = L + b * T1xy + b * b * T2xy
            dtb = math.sqrt(Txy * Txy + tau * tau)
            return dtb - (np.linalg.norm(Pz - y, axis=-1) + b * T1a + b * b * T2a) \
                       - (np.linalg.norm(Pz, axis=-1) + b * T1b + b * b * T2b)

        Vb = []
        for b in (-1e-3, 0.0, 1e-3):
            s0 = dt0 / L
            Sa, Sb = np.full(TE.shape, s0 * 0.999), np.full(TE.shape, s0 * 1.001)
            fa, fb = ell(point(Sa, TE, PH), b), ell(point(Sb, TE, PH), b)
            for _ in range(40):                                        # secant for the root of l(b) on the ray
                den_ = fb - fa
                Sc = Sb - np.where(den_ != 0, fb * (Sb - Sa) / np.where(den_ != 0, den_, 1.0), 0.0)
                fc = ell(point(Sc, TE, PH), b)
                Sa, fa, Sb, fb = Sb, fb, Sc, fc
                if np.max(np.abs(fc)) < 1e-15:
                    break
            sig = 1.0 + (Sb[..., None] - 1.0) * xs_                   # (nte, nph, nsig)
            Pz = point(sig, TE[..., None], PH[..., None])
            hz = np.cos(Pz @ kvec + th)
            f = np.exp(3.0 * b * hz) * ell(Pz, b) * c ** 3 * (sig * sig - TE[..., None] ** 2)
            Vb.append(float((W2[..., None] * (Sb[..., None] - 1.0) * wxs * f).sum()))
        res.append((Vb[0] - 2.0 * Vb[1] + Vb[2]) / (1e-3) ** 2 / V0)
    gam_s = 0.5 * (res[0] + res[1])
    # uniform null: four diamonds from near-apex to sliver
    RSu = np.array([[0.01], [0.5], [3.0], [30.0]])
    tu = np.array([[1.0], [2.0], [0.7], [1.5]])
    dtu = np.sqrt(RSu ** 2 + tu ** 2)
    _, bu, gu, _, _ = second_order_volume(0.0, 0.3, RSu, dtu, tu ** 2, 6, 48, 8)
    unif = float(max(np.max(np.abs(bu)), np.max(np.abs(gu))))
    return {"gam": gam_a, "gam_sph": gam_s, "rel_gam": abs(gam_a - gam_s) / abs(gam_s), "unif": unif}


def check_conformal(k=0.2, T=40.0):
    """I2: the conformal K_hat on this test's nodes against F6.8a's own row."""
    return {"k": k, "T": T, "Kc_8a": plane_wave_row(k, T)["Khat"]}


def _dispatch(kind, args):
    if kind == "row":
        return optical_row(*args)
    if kind == "t2":
        return check_t2()
    if kind == "geom":
        return check_geometry()
    if kind == "conf":
        return check_conformal()
    raise ValueError(kind)


# ================================================================================== run

def _ages(k, p):
    Ta = max(float(p["T_min"]), float(p["n_kT"]) / k)
    return [Ta * float(m) for m in p["ages"]]


def _linfit(T, K):
    if len(T) >= 2:
        A = np.c_[np.ones_like(T), -T]
        (kl, sl), *_ = np.linalg.lstsq(A, K, rcond=None)
        res = float(np.max(np.abs(K - A @ np.array([kl, sl])))) if len(T) > 2 else float("nan")
        return float(kl), float(sl), res
    return float(K[0]), float("nan"), float("nan")


def run(p, ctx):
    t0 = time.time()
    ks = [float(k) for k in p["ks"]]
    quad = tuple(p["quad"])
    bn = int(p["bnd_nodes"])
    jobs = [("row", (k, T, 1.0, quad, bn)) for k in ks for T in _ages(k, p)]
    n1 = len(jobs)
    checks = [(float(c[0]), int(c[1])) for c in p["check"] if float(c[0]) in ks]
    jobs += [("row", (k, _ages(k, p)[i], 2.0, quad, 2 * bn)) for k, i in checks]
    jobs += [("row", (0.0, float(T), 1.0, quad, bn)) for T in p["null_ages"]]
    jobs += [("t2", ()), ("geom", ()), ("conf", ())]
    res = parallel_map(_dispatch, jobs, ctx.workers, ctx.progress,
                       lambda i, r_, d, n: f"[F6.8c] task {d}/{n} ({time.time() - t0:.0f}s)")
    rows = res[:n1]
    lev2 = res[n1:n1 + len(checks)]
    nulls = res[n1 + len(checks):n1 + len(checks) + len(p["null_ages"])]
    t2c, geom, conf = res[-3], res[-2], res[-1]

    lines = ["3+1, eps rho = 1 (l_k = 1); g = -dt^2 + e^(2bh)dx^2, h = cos(k z1), phase-averaged; light cones to second "
             "order (bent rays); K_hat_bb = d2E[S]/db2 / EH 2k^2 (C = 2 K_hat_bb); K_hat = conformal (F6.8a) on the same nodes"]
    ok1a = all(c_["rel"] <= float(p["tol_t2"]) for c_ in t2c)
    lines.append("  I1a <T2> closed form vs numerical Fermat ray: " + "; ".join(
        f"d = {c_['d']:g}, angle {c_['angle']:g}: {c_['ana']:+.6f} vs {c_['num']:+.6f} (rel {c_['rel']:.1e})" for c_ in t2c)
        + f" (tol {p['tol_t2']:g})")
    ok1b = geom["rel_gam"] <= float(p["tol_geom"]) and geom["unif"] <= float(p["tol_unif"])
    lines.append(f"  I1b d2V/db2 {geom['gam']:+.7f} vs numerical roots + finite differences {geom['gam_sph']:+.7f} "
                 f"(rel {geom['rel_gam']:.1e}); uniform null {geom['unif']:.1e} (tol {p['tol_geom']:g}, {p['tol_unif']:g})")
    rc = next((r_ for r_ in rows if r_["k"] == conf["k"] and r_["T"] == conf["T"]), None)
    if rc is None:
        rc = optical_row(conf["k"], conf["T"], 1.0, quad, bn)
    dconf = abs(rc["Kc"] - conf["Kc_8a"])
    ok2 = dconf <= float(p["tol_conf"])
    lines.append(f"  I2 conformal K_hat on this grid {rc['Kc']:+.6f} vs F6.8a row {conf['Kc_8a']:+.6f} "
                 f"(k = {conf['k']:g}, T = {conf['T']:g}): |diff| {dconf:.1e} (tol {p['tol_conf']:g})")
    def k_ref(Tn):
        # the smallest k whose rows reach down to the age Tn (its youngest row is not older than Tn)
        cand = [k for k in ks if min(_ages(k, p)) <= Tn * (1 + 1e-9)]
        return min(cand) if cand else max(ks)
    nunits = [(n_["T"], k_ref(n_["T"]), abs(n_["J_bb"]) * _A4 / k_ref(n_["T"]) ** 2) for n_ in nulls]
    nworst = max(u_[2] for u_ in nunits)
    ok3 = nworst <= float(p["tol_null"])
    lines.append(f"  I3 null k = 0: " + "; ".join(
        f"T = {n_['T']:g}: J_bb {n_['J_bb']:+.1e} (boundary part {n_['bnd']:+.1e}), {u_[2]:.1e} in K_hat_bb units at k = {u_[1]:g}"
        for n_, u_ in zip(nulls, nunits)) + f" (tol {p['tol_null']:g})")
    ok4 = True
    for r2 in lev2:
        r1 = next(r_ for r_ in rows if r_["k"] == r2["k"] and r_["T"] == r2["T"])
        d = abs(r2["Kbb"] - r1["Kbb"])
        ok4 &= d <= float(p["tol_level"])
        lines.append(f"  I4 k = {r2['k']:g}, T = {r2['T']:g}: K_hat_bb level 1 {r1['Kbb']:+.6f}, level 2 {r2['Kbb']:+.6f}, "
                     f"|diff| {d:.1e} (tol {p['tol_level']:g})")
    inst_ok = ok1a and ok1b and ok2 and ok3 and ok4

    lines.append("  k       ages T                   K_hat_bb(T) | K_hat(T)                              K_loc_bb  K_loc     "
                 "c_T,bb   c_T      phi2")
    per = {}
    for k in ks:
        rk = sorted([r_ for r_ in rows if r_["k"] == k], key=lambda r_: r_["T"])
        T = np.array([r_["T"] for r_ in rk])
        Kbb = np.array([r_["Kbb"] for r_ in rk])
        Kc = np.array([r_["Kc"] for r_ in rk])
        lb, sb, rb = _linfit(T, Kbb)
        lc, sc, rc_ = _linfit(T, Kc)
        phi2 = 2.0 * lb - 2.0 - (8.0 / 3.0) * (lc - 1.0)
        per[k] = {"T": T.tolist(), "Kbb": Kbb.tolist(), "Kc": Kc.tolist(), "Kloc_bb": lb, "Kloc": lc,
                  "cT_bb": sb / k ** 3, "cT": sc / k ** 3, "phi2": phi2, "res_bb": rb, "res_c": rc_}
        lines.append(f"  {k:<7g} {', '.join(f'{x:g}' for x in T):24s} "
                     f"{' '.join(f'{a:+.4f}|{b:+.4f}' for a, b in zip(Kbb, Kc)):40s} {lb:+.5f}  {lc:+.5f}  "
                     f"{sb / k ** 3:+.4f}  {sc / k ** 3:+.4f}  {phi2:+.5f}")
    small = [k for k in ks if k <= float(p["fit_kmax"]) + 1e-12]
    fit = {}
    if len(small) >= 4:
        kf = np.array(small)
        yv = np.array([per[k]["Kloc_bb"] for k in small])
        A = np.c_[np.ones_like(kf), -kf ** 2 * np.log(1.0 / kf), -kf ** 2]
        c, *_ = np.linalg.lstsq(A, yv, rcond=None)
        rms = float(np.sqrt(np.mean((yv - A @ c) ** 2)))
        fit = {"K0_bb": float(c[0]), "c1": float(c[1]), "c0": float(c[2]), "rms": rms}
        lines.append(f"  fit K_loc_bb = K0 - k^2 (c1 ln(1/k) + c0) over k <= {p['fit_kmax']:g}: K0_bb = {c[0]:.5f}, "
                     f"c1 = {c[1]:.3f}, c0 = {c[2]:.3f}, rms {rms:.1e}; C(k -> 0) = 2 K0_bb = {2 * c[0]:.4f} (EH 2)")
    # reported, not judged: gamma and Phi with K_hat_ab of F6.8b2 (same k, same T)
    gam_rows = []
    for r_ in rows:
        kab = _b2(r_["k"], r_["T"])
        if kab is None or r_["Kbb"] == 0:
            continue
        g_ = 2.0 * kab / r_["Kbb"] - 1.0
        den_ = 6.0 * r_["Kc"] - 8.0 * kab * kab / r_["Kbb"]
        gam_rows.append({"k": r_["k"], "T": r_["T"], "gamma": g_, "Phi_ratio": -2.0 / den_ if den_ != 0 else float("nan")})
    if gam_rows:
        lines.append("  reported (K_hat_ab of F6.8b2, same k and T): gamma(k, T) = 2 K_ab/K_bb - 1, Phi/Phi_GR = "
                     "-2/(6 K - 8 K_ab^2/K_bb):")
        for g_ in gam_rows:
            lines.append(f"    k = {g_['k']:g}, T = {g_['T']:g}: gamma {g_['gamma']:+.5f}, Phi/Phi_GR {g_['Phi_ratio']:+.5f}")
    loc = []
    for k in ks:
        if k in _B2_LOC and per[k]["Kloc_bb"] != 0:
            loc.append(f"k = {k:g}: gamma_loc {2.0 * _B2_LOC[k] / per[k]['Kloc_bb'] - 1.0:+.5f}")
    if loc:
        lines.append("  reported, local parts: " + "; ".join(loc))

    if not inst_ok:
        status = "INSTRUMENT_FAILED"
    elif not fit:
        status = "UNDECIDED"
    elif abs(fit["K0_bb"] - 1.0) <= float(p["tol_pass"]):
        status = "PASSED"
    elif abs(fit["K0_bb"] - 1.0) > float(p["tol_fail"]):
        status = "FAILED"
    else:
        status = "UNDECIDED"
    data = {"rows": rows, "level2": lev2, "nulls": nulls, "t2": t2c, "geom": geom, "conf": conf,
            "per_k": {str(k): v for k, v in per.items()}, "fit": fit, "gamma_rows": gam_rows}
    return Result(status, lines, data)
