"""
F6.8b -- Mixed static sector (conformal x optical) of the expected Benincasa-Dowker action in 3+1
==================================================================================================

Why this test exists (28 Oct 2026, after F6.8a PASSED; docs/derivations/cesta_k_predikcim.md)
---------------------------------------------------------------------------------------------
Write a static metric as g = Omega^2 [-dt^2 + n^2 dx^2] with omega = ln Omega (the conformal
factor: changes only the density of the sprinkling) and nu = ln n (the optical index: changes the
light cones, i.e. the ORDER). The Einstein-Hilbert (EH) static quadratic form is, per k^2/2 for
plane waves omega = a cos kz, nu = b cos kz (sympy, cesta_k_predikcim.md section 2),

    6 a^2 + 8 a b + 2 b^2.

F6.8a measured the conformal coefficient of the discrete BD action: 6 K_hat(k, T), K0 = 1.00005.
If the expected BD action is COVARIANT at quadratic order (any covariant quadratic action is
R + R F0(box) R + C F2(box) C, and the Weyl tensor of these waves comes from nu alone), then

    A = 6 + 9 phi0,   B = 8 + 12 phi0,   C = 2 + 4 phi0 + phi2,

so the mixed coefficient is FIXED by the conformal one with no free parameter:

    B / 8 = A / 6 = K_hat(k)          (covariance prediction of this test)

This is a sharp test that the discrete action behaves covariantly, and, by the chain of
pravidlo_z_akce.md, the step that turns F6.8a into gamma = 1 without assuming locality in the
non-conformal sectors. F6.8b also builds the geometry of the light cones to first order, which F6.8c
(the nu^2 sector, the spin-2 form factor phi2 that the predictions need) extends to second order.

What is computed (nothing from EH or from covariance enters)
------------------------------------------------------------
The mixed derivative d^2/(da db) at a = b = 0 of the expected action per unit time and volume,
E[S] = -2 <sqrt(-g) E[B_g 1]>, for g = e^(2af) (-dt^2 + e^(2bh) dx^2), f = h = cos(k z1), phase-
averaged; eps rho = 1 (l_k = 1); past of x cut at age T. With P = sqrt(-g) = e^(4af + 3bh):

  * the order: y precedes x iff dt >= T(x, y), the optical distance; to first order in b it is the
    straight segment, T = L + b int_seg h dl (Fermat: the path change is second order). The past of
    x is parametrised by s = l2 (dt^2 - T(x,y)^2) (fixed cone), so the measure dt = ds / (2 l2 dt)
    and the age cut s_max = l2 (T^2 - T_xy^2) depend on b (the latter gives a boundary term);
  * the interval volume V = int d^3z P(z) [dt - T(y,z) - T(z,x)]_+ (static: the time extent above z);
    dV/db = int [3 h l0 + theta(l0) ((L/dt) delta_xy - delta_yz - delta_zx)], d2V/da db similarly,
    where delta_.. are the optical delays along the straight segments (closed form for a plane wave)
    and theta(l0) is the spatial projection of the diamond, a prolate spheroid, integrated in
    prolate spheroidal coordinates (the pointwise difference (L/dt) delta_xy - broken path is formed
    before averaging, so the long thin slivers lose no digits);
  * d^2/da db of P_x P_y gk(V) m(b) expanded analytically with g1 = w gk', g2 = w^2 gk''.
The EH number 8 enters ONLY as the normalisation of the reported ratio K_hat_ab = (d^2 E[S]/da db)/(4 k^2).
The conformal ratio K_hat of F6.8a is computed on the same nodes in the same pass (its own
integrand), so the covariance comparison K_hat_ab - K_hat shares the quadrature.

Development (28 Oct 2026; numbers on record, not used for the criterion)
------------------------------------------------------------------------
- Uniform h (k = 0): dV/db = 0 to 1e-15 (the proper volume of an interval is invariant under a
  constant rescaling of space at fixed optical proper time) -- the geometry is right.
- k = 0.7, CO = 0.4, L = 2, tau = 1.2: dV/db / V0 = 0.0573046 (prolate quadrature, converged to
  1e-10) against 0.0574 (Cartesian grid n = 400, O(h)) and 0.057295 (spherical coordinates around
  the centre, 96 x 64 x 64, slowly converging because of the kinks of d1 + d2 at the foci); the
  phase-averaged d2V/da db: -0.073771 against -0.07351 and -0.073793.
- Null k = 0 (constant omega and nu: flat space rescaled): J_ab = -7e-7 against 0.032 at k = 0.2.
- The null is exact only in the bulk: a constant rescaling changes the proper age e^a T of the cut, so
  at finite T the null equals the T-derivative of the flat finite-age residual (J_ab(k = 0) = -5.9e-5
  at T = 10, -7e-7 at T = 20); the full run checks it at T >= 40.
- First rows: k = 0.2, T = 20: K_hat_ab = 0.65992, K_hat = 0.66027 (F6.8a 0.66028);
  k = 0.1, T = 20 (k T = 2, boundary regime): 0.8934 against 0.8874.

Criteria (stated before the full run)
-------------------------------------
Rows: for every k of the grid three ages T = T_a, 2 T_a, 4 T_a with T_a = max(T_min, n_kT / k).
Instrument (INSTRUMENT_FAILED if any fails):
  I1 geometry: dV/db and the phase-averaged d2V/da db of the test diamond against an independent
     quadrature in spherical coordinates around the diamond centre: relative difference <= 2e-3;
     and the uniform null |dV/db (k = 0)| / V0 <= 1e-9 for four diamonds from apex to sliver;
  I2 the conformal K_hat on the reduced far s-grid of this test against F6.8a's own row
     (k = 0.2, T = 40): |diff| <= 1e-3;
  I3 null k = 0 at three ages: |J_ab(k = 0, T)| <= 1e-3 min |J_ab| over the rows;
  I4 quadrature: |K_hat_ab(level 2) - K_hat_ab(level 1)| <= 2e-3 on the check rows.
Reading: per k, K_hat_ab(k, T) = K_loc_ab(k) - s_ab(k) T and K_hat(k, T) = K_loc(k) - s(k) T fitted
over the three ages; Delta(k) = K_loc_ab(k) - K_loc(k).
Status:
  PASSED     max over k <= 0.1 of |Delta(k)| <= 0.01: the mixed coefficient follows the conformal one,
             as covariance requires (8 K_loc in place of 8);
  FAILED     max over k <= 0.1 of |Delta(k)| > 0.03: the discrete action is not covariant at
             quadratic order in the static sector (then gamma = 1 does NOT follow from F6.8a);
  UNDECIDED  otherwise.
Reported, not judged: Delta(k) for the larger k; the age coefficients c_T,ab = s_ab/k^3 against
c_T = s/k^3 (whether the non-local age term is also covariant); K0_ab from the fit
K_loc_ab = K0_ab - k^2 (c1 ln(1/k) + c0) over k <= 0.1 and c_mix = 8 K0_ab (EH 8).
"""

from __future__ import annotations

import math
import time

import numpy as np

from core.operators.bd_expect_curved import _s_edges_extended
from core.stats.parallel import parallel_map
from core.testing import Result
from phases.F6_dynamika_radu.F6_08a_bd_action_conformal import (_A4, _L2, _gl, angle_nodes, kernel_derivs,
                                                                  mean_exp_diamond, plane_wave_row, rs_nodes)

TEST = {
    "id": "F6.8b",
    "phase": "F6",
    "title": "Mixed static sector (conformal x optical) of the expected BD action: is B/8 = K_hat (covariance)?",
    "theory_part": "derivation of the substrate rule from the order's action (pravidlo_z_akce.md, cesta_k_predikcim.md); "
                   "second stage of F6.8",
    "question": "Does the mixed (omega x nu) static coefficient of the expected discrete BD action follow the conformal "
                "one as covariance requires (B/8 = A/6 = K_hat at every k), i.e. is the discrete action covariant at "
                "quadratic order, with nothing from EH or covariance put into the computation?",
    "method": "deterministic quadrature in 3+1, eps rho = 1 (l_k = 1); g = e^(2af)(-dt^2 + e^(2bh)dx^2), f = h = "
              "cos(k z1), phase-averaged; order through the optical distance to first order (straight segments); "
              "interval volume through the spatial projection of the diamond with the optical delays; d^2/da db "
              "analytic; past of x cut at age T; per k three ages with k T >= 8; K_hat_ab = ratio to the EH value "
              "8 (4 k^2); the conformal K_hat (F6.8a) on the same nodes",
    "prediction": "covariance: K_hat_ab(k) = K_hat(k) (local parts); EH: K_hat_ab -> 1 as k -> 0",
    "criterion": "INSTRUMENT: I1 geometry vs independent spherical quadrature (<= 2e-3 rel.) and uniform null "
                 "(<= 1e-9); I2 conformal K_hat vs F6.8a row (<= 1e-3); I3 null k = 0 (<= 1e-3); I4 level 2 vs 1 "
                 "(<= 2e-3). Delta(k) = K_loc_ab - K_loc from linear fits in T. PASSED max_{k<=0.1} |Delta| <= 0.01; "
                 "FAILED > 0.03; UNDECIDED otherwise.",
    "null_model": "constant conformal factor and constant index (k = 0): flat space rescaled, the response vanishes",
    "level": "chain",
    "meters": [],
    "depends_on": ["F6.8a"],
    "params": {"ks": [0.035, 0.05, 0.07, 0.1, 0.14, 0.2, 0.3], "T_min": 40.0, "n_kT": 8.0, "ages": [1, 2, 4],
               "fit_kmax": 0.1, "check": [[0.05, 0], [0.1, 0], [0.2, 0]], "null_ages": [40.0, 160.0, 640.0],
               "quad": [4, 16, 0.8, 8],
               "tol_geom": 2e-3, "tol_unif": 1e-9, "tol_conf": 1e-3, "tol_null": 1e-3, "tol_level": 2e-3,
               "tol_pass": 0.01, "tol_fail": 0.03},
    "parts": {"main": {}},
    "quick": {"ks": [0.1, 0.2], "T_min": 20.0, "n_kT": 1.0, "ages": [1, 2], "check": [[0.2, 0]],
              "null_ages": [20.0]},
    "estimate": "overnight: ~5-10 CPU-hours (21 rows of level 1 at 2-30 min, 3 rows of level 2 at ~1-2 h, 3 null rows, "
                "instruments ~5 min); with 12 processes ~2-3 h wall; memory < 1 GB per process",
}


# ======================================================================== geometry of the delays

def _ediv(num, p, base):
    """num / (i p) where num = base (e^{ip} - 1); the series base (1 + i p/2 - p^2/6) where |p| < 1e-4."""
    num, p = np.broadcast_arrays(num, p)
    small = np.abs(p) < 1e-4
    out = np.asarray(base * (1.0 + 0.5j * p - p * p / 6.0), complex) * np.ones(p.shape)
    np.divide(num, 1j * p, out=out, where=~small)
    return out


def spheroid_pieces(k, CO, RS, dt, n_sig, n_te, n_phi):
    """
    Uniform averages over the spatial projection of the diamond x -> y (prolate spheroid d1 + d2 <= dt;
    foci y and x, L = RS, cos(n, e1) = CO), amplitudes of e^{ik.o} relative to x (o = z - x):
        Dm = mean[(L/dt) D_xy - broken(z)],     Qm = mean[e^{ik.o} conj((L/dt) D_xy - broken(z))]
    broken(z) = delay along y -> z plus along z -> x, D_xy = delay along x -> y (closed forms:
    int_0^d e^{i(phi0 + p l/d)} dl = d e^{i phi0} (e^{ip} - 1)/(ip)). Prolate coordinates as in the core:
    sigma in [1, dt/L] (Gauss), te in [-1, 1] (Gauss), phi (periodic trapezoid); d1 = c (sigma + te) to y,
    d2 = c (sigma - te) to x, axial distance from y xi = c (1 + sigma te), r_perp = c sqrt((s^2-1)(1-te^2)),
    weight c^3 (sigma^2 - te^2). Vectorised over te and phi; loop over sigma.
    """
    L = RS[..., None, None]
    dtt = dt[..., None, None]
    c = 0.5 * L
    smax = np.maximum(dtt / np.maximum(L, 1e-300), 1.0)
    kpar = k * CO
    kperp = k * math.sqrt(max(1.0 - CO * CO, 0.0))
    Dxy = RS * np.exp(1j * kpar * RS / 2.0) * np.sinc(kpar * RS / 2.0 / math.pi)
    LD = (RS / dt * Dxy)[..., None, None]
    xs, ws = _gl(0.0, 1.0, n_sig)
    ts, wt = _gl(-1.0, 1.0, n_te)
    te = ts[:, None]
    wte = wt[:, None]
    cph = np.cos(2.0 * math.pi * (np.arange(n_phi) + 0.5) / n_phi)[None, :]
    eL = np.exp(1j * kpar * L)
    num_D = np.zeros(RS.shape, complex)
    num_Q = np.zeros(RS.shape, complex)
    den = np.zeros(RS.shape)
    for x, wx in zip(xs, ws):
        sig = 1.0 + (smax - 1.0) * x
        jac = (smax - 1.0) * wx
        xi = c * (1.0 + sig * te)
        d1 = c * (sig + te)
        d2 = c * (sig - te)
        rp = c * np.sqrt(np.maximum((sig * sig - 1.0) * (1.0 - te * te), 0.0))
        wb = jac * wte * c ** 3 * (sig * sig - te * te) / n_phi
        ko = kpar * (L - xi) + kperp * rp * cph          # k . (z - x)
        kzy = ko - kpar * L                               # k . (z - y)
        eko = np.exp(1j * ko)
        broken = d1 * _ediv(eko - eL, kzy, eL) + d2 * _ediv(eko - 1.0, ko, 1.0)
        dlt = LD - broken
        wbb = np.broadcast_to(wb, dlt.shape)
        num_D += (wbb * dlt).sum(axis=(-2, -1))
        num_Q += (wbb * eko * np.conj(dlt)).sum(axis=(-2, -1))
        den += wbb.sum(axis=(-2, -1))
    ok = den > 0
    dd = np.where(ok, den, 1.0)
    return np.where(ok, num_D / dd, 0.0), np.where(ok, num_Q / dd, 0.0), Dxy


# ======================================================================== quadrature of the past

def _s_nodes(s_max, n_s, far):
    """Panels in s cut at s_max per rs node. far = True drops the apex panels (s < 0.2) and stops at s = 12."""
    e = _s_edges_extended(4, 1.0)
    if far:
        e = np.r_[0.0, e[(e > 0.2) & (e < 12.0)], 12.0]
    u, wu = _gl(0.0, 1.0, n_s)
    S, WS = [], []
    for lo, hi_e in zip(e[:-1], e[1:]):
        hi = np.minimum(hi_e, s_max)[:, None]
        span = np.maximum(hi - lo, 0.0)
        S.append(lo + span * u)
        WS.append(span * wu)
    return np.concatenate(S, axis=1), np.concatenate(WS, axis=1)


def _pair_terms(k, CO, RS, s, dt, quad):
    """Phase-averaged mixed integrand and the conformal (F6.8a) integrand on the same nodes."""
    n_sig, te0, te_per_rad, n_phi = quad
    tau2 = s / _L2
    g0, g1, g2 = kernel_derivs(s * s)
    ph = np.exp(1j * k * CO * RS)                  # e^{ik delta}, delta = y1 - x1
    one = 1.0 + ph
    M = mean_exp_diamond(k, RS, CO, dt) if k > 0 else np.ones(RS.shape, complex)
    n_te = int(te0 + math.ceil(te_per_rad * abs(k * CO) * (float(np.max(RS)) if RS.size else 0.0)))
    Dm, Qm, Dxy = spheroid_pieces(k, CO, RS, dt, n_sig, n_te, n_phi)
    AV = np.where(tau2 > 0, 4.0 * dt / np.maximum(tau2, 1e-300), 0.0)    # A / V0 = 4 dt / tau^2
    alpha = 4.0 * M                                    # dV/da / V0
    beta = 3.0 * M + AV * Dm                           # dV/db / V0 (h amplitude)
    gam = 6.0 + 2.0 * AV * np.real(Qm)                 # phase-averaged d2V/da db / V0
    mu = -RS * Dxy / (dt * dt)                         # d ln m / db, m = 1/(2 l2 dt)
    mixed = (12.0 * (1.0 + ph.real) * g0                                  # E_ab g0
             + 2.0 * g1 * np.real(one * np.conj(beta))                    # E_a G_b
             + 1.5 * g1 * np.real(alpha * np.conj(one))                   # E_b G_a
             + 0.5 * g2 * np.real(alpha * np.conj(beta)) + g1 * gam       # G_ab
             + 0.5 * np.real((4.0 * one * g0 + g1 * alpha) * np.conj(mu)))  # m_b (E_a g0 + G_a)
    conf = (1.0 + ph.real) * g0 + np.real(one * np.conj(M)) * g1 + 0.5 * np.abs(M) ** 2 * g2 + 0.5 * g1
    return mixed, conf, n_te


def mixed_row(k, T, mult=1.0, quad=(4, 16, 0.8, 8)):
    """
    Phase-averaged mixed second-order density J_ab at age T (and the conformal J of F6.8a on the same
    nodes). J_ab = 6 (from P_x) + int d^3y ds/(2 l2 dt) [mixed] + boundary term of the age cut
    (d s_max/db = -2 l2 L delta_xy times dF/da at s_max). K_hat_ab = a_4 J_ab / (2 k^2).
    """
    t0 = time.time()
    q = (int(quad[0]), int(quad[1] * mult), float(quad[2]) * mult, int(quad[3]))
    rs, wr = rs_nodes(T, k, mult)
    cos, wa = angle_nodes(k, T, mult)
    n_s = int((quad[4] if len(quad) > 4 else 16) * mult)       # Gauss nodes per s panel (F6.8a: 16)
    inner = rs < 4.0
    Jab, Jc, bnd, nte_max = 6.0, 0.5, 0.0, 0
    smx_all = _L2 * np.maximum(T * T - rs * rs, 0.0)
    lay = smx_all < 30.0
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
                mixed, conf, nte = _pair_terms(k, CO, RS, s, dt, q)
                nte_max = max(nte_max, nte)
                Jab += wa[a] * float((base * mixed).sum())
                Jc += wa[a] * float((base * conf).sum())
        if lay.any():
            r = rs[lay][:, None]
            sm = smx_all[lay][:, None]
            g0, g1, _ = kernel_derivs(sm * sm)
            M = mean_exp_diamond(k, r, CO, np.full_like(r, T)) if k > 0 else np.ones(r.shape, complex)
            one = 1.0 + np.exp(1j * k * CO * r)
            dFa = (4.0 * one * g0 + g1 * 4.0 * M) / (2.0 * _L2 * T)
            Dxy = r * np.exp(1j * k * CO * r / 2.0) * np.sinc(k * CO * r / 2.0 / math.pi)
            dsb = -2.0 * _L2 * r * Dxy
            bnd += wa[a] * float((wr[lay][:, None] * r ** 2 * 0.5 * np.real(dFa * np.conj(dsb))).sum())
    Jab += bnd
    out = {"k": k, "T": T, "mult": mult, "J_ab": Jab, "J_c": Jc, "bnd": bnd, "n_ang": int(cos.size),
           "n_rs": int(rs.size), "n_te_max": nte_max, "time": time.time() - t0}
    if k > 0:
        out["Kab"] = _A4 * Jab / (2.0 * k * k)
        out["Kc"] = _A4 * 16.0 * Jc / (3.0 * k * k)
    return out


# ============================================================================ instruments

def check_geometry():
    """I1: the prolate quadrature against spherical coordinates around the diamond centre; the uniform null."""
    k, CO, L, tau = 0.7, 0.4, 2.0, 1.2
    dt = math.sqrt(L * L + tau * tau)
    RSa = np.array([[L]])
    dta = np.array([[dt]])
    M = mean_exp_diamond(k, RSa, CO, dta)
    AV = 4.0 * dt / tau ** 2
    Dm, Qm, _ = spheroid_pieces(k, CO, RSa, dta, 16, 128, 32)
    beta_a = float((3.0 * M + AV * Dm).real[0, 0])
    gam_a = float((6.0 + 2.0 * AV * np.real(Qm))[0, 0])
    # independent: spherical coordinates around the centre, exact ellipsoid boundary, straight-segment delays
    V0 = math.pi * tau ** 4 / 24.0
    nvec = np.array([CO, math.sqrt(1 - CO * CO), 0.0])
    y = L * nvec
    cen = 0.5 * L * nvec
    e2 = np.array([-nvec[1], nvec[0], 0.0])
    e3 = np.cross(nvec, e2)
    a_ax, b_ax = dt / 2.0, tau / 2.0

    def li(A, B, th):
        d = B - A
        Ln = np.linalg.norm(d, axis=-1)
        kp = k * d[..., 0]
        return Ln, Ln * np.real(np.exp(1j * (k * A[..., 0] + th)) * np.exp(1j * kp / 2) * np.sinc(kp / 2 / np.pi))

    nr, nc, nph = 128, 96, 96
    cs, wcs = _gl(-1.0, 1.0, nc)
    phs = 2.0 * math.pi * (np.arange(nph) + 0.5) / nph
    res = {}
    for th in (0.0, 0.5 * math.pi):
        _, dxy = li(np.zeros(3), y, th)
        sb = sab = 0.0
        for ci, wci in zip(cs, wcs):
            si = math.sqrt(1 - ci * ci)
            for p in phs:
                m = ci * nvec + si * math.cos(p) * e2 + si * math.sin(p) * e3
                Rm = 1.0 / math.sqrt((ci / a_ax) ** 2 + (si / b_ax) ** 2)
                r, wr = _gl(0.0, Rm, nr)
                P = cen[None, :] + r[:, None] * m[None, :]
                L1, d1 = li(y[None, :], P, th)
                L2, d2 = li(P, np.zeros((1, 3)), th)
                ell0 = dt - L1 - L2
                ellb = (L / dt) * dxy - d1 - d2
                h = np.cos(k * P[:, 0] + th)
                W = wci * (2 * math.pi / nph) * wr * r * r
                sb += float((W * (3 * h * ell0 + ellb)).sum())
                sab += float((W * (12 * h * h * ell0 + 4 * h * ellb)).sum())
        res[th] = (sb / V0, sab / V0)
    beta_s = res[0.0][0]
    gam_s = 0.5 * (res[0.0][1] + res[0.5 * math.pi][1])
    # uniform null: four diamonds from near-apex to sliver
    RSu = np.array([[0.01], [0.5], [3.0], [30.0]])
    tu = np.array([[1.0], [2.0], [0.7], [1.5]])
    dtu = np.sqrt(RSu ** 2 + tu ** 2)
    Dmu, _, _ = spheroid_pieces(0.0, 0.3, RSu, dtu, 6, 48, 8)
    unif = float(np.max(np.abs(3.0 + 4.0 * dtu / tu ** 2 * Dmu)))
    return {"beta": beta_a, "beta_sph": beta_s, "gam": gam_a, "gam_sph": gam_s,
            "rel_beta": abs(beta_a - beta_s) / abs(beta_s), "rel_gam": abs(gam_a - gam_s) / abs(gam_s), "unif": unif}


def check_conformal(k=0.2, T=40.0):
    """I2: the conformal K_hat on this test's nodes against F6.8a's own row."""
    return {"k": k, "T": T, "Kc_8a": plane_wave_row(k, T)["Khat"]}


def _dispatch(kind, args):
    if kind == "row":
        return mixed_row(*args)
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
    jobs = [("row", (k, T, 1.0, quad)) for k in ks for T in _ages(k, p)]
    n1 = len(jobs)
    checks = [(float(c[0]), int(c[1])) for c in p["check"] if float(c[0]) in ks]
    jobs += [("row", (k, _ages(k, p)[i], 2.0, quad)) for k, i in checks]
    jobs += [("row", (0.0, float(T), 1.0, quad)) for T in p["null_ages"]]
    jobs += [("geom", ()), ("conf", ())]
    res = parallel_map(_dispatch, jobs, ctx.workers, ctx.progress,
                       lambda i, r_, d, n: f"[F6.8b] task {d}/{n} ({time.time() - t0:.0f}s)")
    rows = res[:n1]
    lev2 = res[n1:n1 + len(checks)]
    nulls = res[n1 + len(checks):n1 + len(checks) + len(p["null_ages"])]
    geom, conf = res[-2], res[-1]

    lines = ["3+1, eps rho = 1 (l_k = 1); g = e^(2af)(-dt^2 + e^(2bh)dx^2), f = h = cos(k z1), phase-averaged; "
             "K_hat_ab = d2E[S]/da db / EH 4k^2; K_hat = conformal (F6.8a) on the same nodes"]
    ok1 = (geom["rel_beta"] <= float(p["tol_geom"]) and geom["rel_gam"] <= float(p["tol_geom"])
           and geom["unif"] <= float(p["tol_unif"]))
    lines.append(f"  I1 geometry: dV/db {geom['beta']:+.7f} vs spherical {geom['beta_sph']:+.7f} (rel {geom['rel_beta']:.1e}); "
                 f"d2V/da db {geom['gam']:+.7f} vs {geom['gam_sph']:+.7f} (rel {geom['rel_gam']:.1e}); uniform null "
                 f"{geom['unif']:.1e} (tol {p['tol_geom']:g}, {p['tol_unif']:g})")
    rc = next((r_ for r_ in rows if r_["k"] == conf["k"] and r_["T"] == conf["T"]), None)
    if rc is None:
        rc = mixed_row(conf["k"], conf["T"], 1.0, quad)
    dconf = abs(rc["Kc"] - conf["Kc_8a"])
    ok2 = dconf <= float(p["tol_conf"])
    lines.append(f"  I2 conformal K_hat on this grid {rc['Kc']:+.6f} vs F6.8a row {conf['Kc_8a']:+.6f} "
                 f"(k = {conf['k']:g}, T = {conf['T']:g}): |diff| {dconf:.1e} (tol {p['tol_conf']:g})")
    jmin = min(abs(r_["J_ab"]) for r_ in rows)
    nworst = max(abs(n_["J_ab"]) for n_ in nulls) / jmin
    ok3 = nworst <= float(p["tol_null"])
    lines.append(f"  I3 null k = 0 at T = {', '.join(f'{n_['T']:g}' for n_ in nulls)}: J_ab = "
                 f"{', '.join(f'{n_['J_ab']:+.1e}' for n_ in nulls)}; max |J0| / min |J_ab| = {nworst:.1e} "
                 f"(tol {p['tol_null']:g})")
    ok4 = True
    for r2 in lev2:
        r1 = next(r_ for r_ in rows if r_["k"] == r2["k"] and r_["T"] == r2["T"])
        d = abs(r2["Kab"] - r1["Kab"])
        ok4 &= d <= float(p["tol_level"])
        lines.append(f"  I4 k = {r2['k']:g}, T = {r2['T']:g}: K_hat_ab level 1 {r1['Kab']:+.6f}, level 2 {r2['Kab']:+.6f}, "
                     f"|diff| {d:.1e} (tol {p['tol_level']:g})")
    inst_ok = ok1 and ok2 and ok3 and ok4

    lines.append("  k       ages T                   K_hat_ab(T) | K_hat(T)                              K_loc_ab  K_loc     "
                 "Delta     c_T,ab   c_T")
    per = {}
    for k in ks:
        rk = sorted([r_ for r_ in rows if r_["k"] == k], key=lambda r_: r_["T"])
        T = np.array([r_["T"] for r_ in rk])
        Kab = np.array([r_["Kab"] for r_ in rk])
        Kc = np.array([r_["Kc"] for r_ in rk])
        la, sa, ra = _linfit(T, Kab)
        lc, sc, rc_ = _linfit(T, Kc)
        per[k] = {"T": T.tolist(), "Kab": Kab.tolist(), "Kc": Kc.tolist(), "Kloc_ab": la, "Kloc": lc, "Delta": la - lc,
                  "cT_ab": sa / k ** 3, "cT": sc / k ** 3, "res_ab": ra, "res_c": rc_}
        lines.append(f"  {k:<7g} {', '.join(f'{x:g}' for x in T):24s} "
                     f"{' '.join(f'{a:+.4f}|{b:+.4f}' for a, b in zip(Kab, Kc)):40s} {la:+.5f}  {lc:+.5f}  "
                     f"{la - lc:+.5f}  {sa / k ** 3:+.4f}  {sc / k ** 3:+.4f}")
    small = [k for k in ks if k <= float(p["fit_kmax"]) + 1e-12]
    dmax = max(abs(per[k]["Delta"]) for k in small) if small else float("nan")
    lines.append(f"  max over k <= {p['fit_kmax']:g} of |Delta| = {dmax:.2e}")
    fit = {}
    if len(small) >= 4:
        kf = np.array(small)
        y = np.array([per[k]["Kloc_ab"] for k in small])
        A = np.c_[np.ones_like(kf), -kf ** 2 * np.log(1.0 / kf), -kf ** 2]
        c, *_ = np.linalg.lstsq(A, y, rcond=None)
        rms = float(np.sqrt(np.mean((y - A @ c) ** 2)))
        fit = {"K0_ab": float(c[0]), "c1": float(c[1]), "c0": float(c[2]), "rms": rms}
        lines.append(f"  fit K_loc_ab = K0 - k^2 (c1 ln(1/k) + c0) over k <= {p['fit_kmax']:g}: K0_ab = {c[0]:.5f}, "
                     f"c1 = {c[1]:.3f}, c0 = {c[2]:.3f}, rms {rms:.1e}; mixed coefficient c_mix = 8 K0_ab = {8 * c[0]:.4f} (EH 8)")
    if not inst_ok:
        status = "INSTRUMENT_FAILED"
    elif not small or not math.isfinite(dmax):
        status = "UNDECIDED"
    elif dmax <= float(p["tol_pass"]):
        status = "PASSED"
    elif dmax > float(p["tol_fail"]):
        status = "FAILED"
    else:
        status = "UNDECIDED"
    data = {"rows": rows, "level2": lev2, "nulls": nulls, "geom": geom, "conf": conf,
            "per_k": {str(k): v for k, v in per.items()}, "fit": fit}
    return Result(status, lines, data)
