"""
F6.10b -- The response of the pairwise matter action of a wall to the static metric in 3+1 (OV2 T2)
==================================================================================================

DRAFT, not frozen (21 Nov 2026, prepared by Claude at the author's request; card docs/OV2_coupling_cards_draft.md,
T2). The criterion below is a proposal until the author freezes the card.

Why this test exists
--------------------
OV2 (docs/derivations/order_variation_OV2.md, decision 39): matter couples to the order through a matter action that
depends on the order only, the pairwise (Dirichlet) kinetic action plus the potential as a measure. It acts on the
order as matter must only if its first variation with the metric is (1/2) int T^{ab} delta g_ab sqrt|g| at fixed matter.
In 3+1 the kinetic part is not conformally invariant, so it is tested in both static sectors of family M:
  omega (conformal, g = e^{2 b cos kz} eta): continuum -int cos(kz) phi'^2 dz = -4 pi k / sinh(pi k / 2m)
      (the density and every interval volume change; the order does not);
  nu (the cones, g = -dt^2 + e^{2 b cos kz} dx^2): continuum -(1/2) int cos(kz) phi'^2 dz = -2 pi k / sinh(pi k / 2m)
      (the order itself changes: first-order cones and interval volumes of F6.8b2).
The potential parts are the measure and exact; they are not judged.

What is computed (ensemble level, deterministic; no sprinkling)
---------------------------------------------------------------
Units eps rho = 1 (l_k = 1, as F6.8). A static sine-Gordon wall phi = 4 arctan(e^{m z}) on a flat 3+1 order of age T.
By the Mecke formula the expected kinetic action per unit time and area is
    s_kin = (a_4/4) int dz int_{J^-(x), age T} g(eps rho V) (phi_y - phi_x)^2 dV,   g(w) = e^{-w}(-1 + 9w - 8w^2 + 4w^3/3),
and its first derivative in b is computed analytically under the integral:
  omega: the diamond average of cos kz boosted into the rest frame of each interval (port of F6.8a diamond_ft4);
  nu:    the first-order cones (optical length of the straight segment), the past parametrised by s with the measure
         1/(2 l2 dt) and the age cut depending on b, dV/db from F6.8b2 spheroid_pieces, and the boundary term of the
         age cut in S = s_max.
The quadrature code is the development code of diag/ov2_kernels4.py (md5 1b618890), diag/ov2_t2.py (ed9d0f96) and
diag/ov2_t2nu.py (2a2e2144), copied here verbatim (one rename: the numpy wall profile of ov2_t2nu is phi_np) so the
frozen test does not depend on files that may change.

Parameters (fresh; development used width 1/m = 10 and k = 0.02, 0.05, 0.1, 0.2): width 1/m = 8 l_k (m = 0.125),
k l_k = 0, 0.03, 0.06, 0.12; ages T = 120 (judged) and 80 (I1); grid multiplier 1 (judged) and 1.5 (I2).

Criterion (proposed)
--------------------
ratio = d s_kin/db / continuum, per sector and k, at T = 120 and grid 1. Delta = |ratio(T 120) - ratio(T 80)| +
|ratio(grid 1.5) - ratio(grid 1)|.
  PASSED if in both sectors |ratio(0) - 1| + Delta <= 0.02, the ratios do not increase with k
         (ratio(0) >= ratio(0.03) >= ratio(0.06) >= ratio(0.12)), and |ratio(k) - 1| + Delta <= 0.15 at k <= 0.06;
  FAILED if in either sector |ratio(0) - 1| - Delta > 0.05 or |ratio(k) - 1| - Delta > 0.25 at some k <= 0.06;
  UNDECIDED otherwise.
Instrument (INSTRUMENT_FAILED if any fails):
  I1 age: every row at T = 80 within 1 % of T = 120;
  I2 grid: every row at grid 1.5 within 1 % of grid 1;
  I3 the omega row at k = 0 equals the uniform conformal derivative of the independent stage-1 kernel (wall_kin mode 1)
     to 1e-5 relative;
  I4 the geometry of F6.8b2 (its check_geometry: dV/db and d2V/da db against an independent spherical quadrature,
     <= 2e-3 relative; the uniform null <= 1e-9).
Reported, not judged: the form factors against K_loc(k) of F6.8a (omega) and the development rows (OV2 sec. 12).
Development rows at this test's fresh width (21 Nov, while timing the draft): k 0.06, T 80: omega 0.954, nu 0.949.
"""

from __future__ import annotations

import math
import time

import numba as nb
import numpy as np

from core.operators.bd_expect_curved import _s_edges_extended
from core.stats.parallel import parallel_map
from core.testing import Result
from phases.F6_order_dynamics.F6_08a_bd_action_conformal import (_A4, _L2, _gl, angle_nodes, kernel_derivs,
                                                                  mean_exp_diamond, rs_nodes, s_nodes)
from phases.F6_order_dynamics.F6_08b2_bd_action_mixed_rev import check_geometry, spheroid_pieces

TEST = {
    "id": "F6.10b",
    "phase": "F6",
    "title": "Response of the pairwise matter action of a wall to the static 3+1 metric, omega and nu sectors "
             "(OV2 T2; draft)",
    "theory_part": "OV2 (docs/derivations/order_variation_OV2.md): the coupling of matter to the order; T_ab from "
                   "the first variation of the matter action at fixed matter",
    "question": "Does the first-order response of the expected pairwise kinetic action of a static wall to a conformal "
                "(omega) and a cone (nu) perturbation cos(kz) equal the continuum (1/2) int T^ab delta g_ab times a "
                "form factor that tends to 1 as k l_k -> 0?",
    "method": "Mecke formula with the exact kernel average g(eps rho V); analytic first derivative; omega by the "
              "boosted diamond average, nu by the first-order cones and prolate volumes of F6.8b2; eps rho = 1, "
              "width 8 l_k, k 0/0.03/0.06/0.12, T 120 and 80, grid 1 and 1.5",
    "prediction": "ratios -> 1 at k = 0 (development: omega 0.998, nu 1.0045 at width 10), falling with k "
                  "(development k 0.05: 0.965 and 0.961)",
    "criterion": "PROPOSED: PASSED if in both sectors |ratio(0) - 1| + Delta <= 0.02, ratios non-increasing in k and "
                 "|ratio(k) - 1| + Delta <= 0.15 at k <= 0.06; FAILED if |ratio(0) - 1| - Delta > 0.05 or "
                 "|ratio(k) - 1| - Delta > 0.25 at k <= 0.06 in either sector; else UNDECIDED; INSTRUMENT_FAILED if "
                 "I1 (age 1 %), I2 (grid 1 %), I3 (omega k = 0 vs stage-1 kernel 1e-5) or I4 (F6.8b2 geometry) fails",
    "null_model": "the potential part (the measure), exact by construction",
    "level": "mechanism",
    "meters": [],
    "depends_on": ["F6.8a", "F6.8b2", "F6.10a"],
    "params": {"m": 0.125, "ks": [0.0, 0.03, 0.06, 0.12], "T": 120.0, "T_age": 80.0, "mult": 1.0, "mult_grid": 1.5,
               "tol0": 0.02, "tolk": 0.15, "kmax_judge": 0.06, "fail0": 0.05, "failk": 0.25, "tol_age": 0.01,
               "tol_grid": 0.01, "tol_i3": 1e-5, "tol_geom": 2e-3, "tol_unif": 1e-9},
    "parts": {"main": {}},
    "quick": {"ks": [0.0, 0.06], "T": 40.0, "T_age": 30.0, "mult_grid": 1.25},
    "estimate": "~25 CPU-min per sector and age at grid 1, ~2-3x that at grid 1.5; ~3 CPU-h total, ~20 min wall "
                "with 12 workers",
}


# ================================================================================== kernels (diag/ov2_kernels4.py, verbatim)

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


# ================================================================================== omega sector (diag/ov2_t2.py, verbatim; run -> run_stage1)

A4 = 4.0 / math.sqrt(6.0)


def panels(edges, n=8):
    x, w = np.polynomial.legendre.leggauss(n)
    xs, ws = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        xs.append(0.5 * (b + a) + 0.5 * (b - a) * x)
        ws.append(0.5 * (b - a) * w)
    return np.concatenate(xs), np.concatenate(ws)


def s_edges():
    w_edges = np.array([0.05, 0.25, 1.0, 2.5, 5.0, 10.0, 20.0, 40.0, 80.0])
    geo = [10.0 ** k for k in range(-14, 0) if 10.0 ** k < math.sqrt(0.05)]
    return np.array([0.0] + geo + list(np.sqrt(w_edges)))


def run_stage1(m, T, nr=60, nz=40):
    kb = np.array(sorted({k / m for k in range(-10, 11)}))
    zedges = np.unique(np.concatenate([kb, -np.geomspace(10 / m, T, nz // 2), np.geomspace(10 / m, T, nz // 2)]))
    zedges = zedges[(zedges >= -T) & (zedges <= T)]
    zs, wz = panels(zedges)
    # near r = T the age cut s <= sqrt(kappa) (T^2 - r^2) clips the kernel: panels in r that follow the s edges
    # there (r_k = sqrt(T^2 - sigma_k / sqrt(kappa))), so the clipped region is resolved as the s-integral is
    sig = np.concatenate([s_edges()[1:], np.linspace(math.sqrt(80.0), 40.0, 40)[1:]])
    r_age = np.sqrt(np.clip(T * T - sig / math.sqrt(math.pi / 24.0), 0.0, None))
    r_age = r_age[r_age > 0.5 * T]
    redges = np.unique(np.concatenate([[0.0], np.geomspace(1e-4, T, nr), r_age, [T]]))
    rs, wr = panels(redges)
    g8x, g8w = np.polynomial.legendre.leggauss(8)
    se = s_edges()
    out = {}
    for mode in (0, 1):
        out[mode] = A4 / 4.0 * wall_kin(m, T, zs, wz, rs, wr, se, kb, g8x, g8w, mode)
    return out


def run_k(m, T, kw, nr=60, nz=40):
    """d s_kin/db for Omega = e^{b cos(kw z)}; continuum -int cos(kw z) phi'^2 dz = -4 pi kw / sinh(pi kw / (2 m))."""
    kb = np.array(sorted({k / m for k in range(-10, 11)}))
    zedges = np.unique(np.concatenate([kb, -np.geomspace(10 / m, T, nz // 2), np.geomspace(10 / m, T, nz // 2)]))
    zedges = zedges[(zedges >= -T) & (zedges <= T)]
    zs, wz = panels(zedges)
    sig = np.concatenate([s_edges()[1:], np.linspace(math.sqrt(80.0), 40.0, 40)[1:]])
    r_age = np.sqrt(np.clip(T * T - sig / math.sqrt(math.pi / 24.0), 0.0, None))
    r_age = r_age[r_age > 0.5 * T]
    redges = np.unique(np.concatenate([[0.0], np.geomspace(1e-4, T, nr), r_age, [T]]))
    rs, wr = panels(redges)
    g8x, g8w = np.polynomial.legendre.leggauss(8)
    g16x, g16w = np.polynomial.legendre.leggauss(16)
    d = A4 / 4.0 * wall_kin_k(m, kw, T, zs, wz, rs, wr, s_edges(), kb, g8x, g8w, g16x, g16w)
    cont = -4 * math.pi * kw / math.sinh(math.pi * kw / (2 * m)) if kw > 0 else -8 * m
    return d, cont


# ================================================================================== nu sector (diag/ov2_t2nu.py, verbatim; phi_w -> phi_np)

def phi_np(z, m):
    return 4.0 * np.arctan(np.exp(m * z))


def Phi_of(k, m, Delta):
    """Phi(Delta) = int dz e^{ikz} (phi(z + Delta) - phi(z))^2, composite Gauss-Legendre over the support."""
    Delta = np.asarray(Delta, float)
    out = np.zeros(Delta.shape, complex)
    x8, w8 = np.polynomial.legendre.leggauss(16)
    for idx, D in np.ndenumerate(Delta):
        lo = min(0.0, -D) - 30.0 / m
        hi = max(0.0, -D) + 30.0 / m
        npan = int(math.ceil((hi - lo) / min(0.5 / m, (1.0 / k) if k > 0 else 1e9)))
        e = np.linspace(lo, hi, npan + 1)
        z = (0.5 * (e[1:] + e[:-1])[:, None] + 0.5 * (e[1:] - e[:-1])[:, None] * x8[None, :]).ravel()
        w = (0.5 * (e[1:] - e[:-1])[:, None] * w8[None, :]).ravel()
        d = phi_np(z + D, m) - phi_np(z, m)
        out[idx] = np.sum(w * np.exp(1j * k * z) * d * d)
    return out


def nu_response(k, m, T, mult=1.0, quad=(4, 16, 0.8, 8)):
    t0 = time.time()
    rs, wr = rs_nodes(T, k, mult)
    cos, wa = angle_nodes(k, T, mult)
    s, ws = s_nodes(_L2 * np.maximum(T * T - rs * rs, 0.0), int(16 * mult))
    RS = rs[:, None] * np.ones_like(s)
    dt = np.sqrt(RS * RS + s / _L2)
    base = wr[:, None] * ws * RS ** 2 / (_L2 * 2.0 * dt)
    g0, g1, _ = kernel_derivs(s * s)
    tau2 = s / _L2
    AV = np.where(tau2 > 0, 4.0 * dt / np.maximum(tau2, 1e-300), 0.0)
    n_sig, te0, te_per_rad, n_phi = quad
    bulk = 0.0
    for a in range(cos.size):
        CO = float(cos[a])
        PhiR = Phi_of(k, m, rs * CO)[:, None]
        ph = np.exp(1j * k * CO * RS)
        M = mean_exp_diamond(k, RS, CO, dt) if k > 0 else np.ones(RS.shape, complex)
        n_te = int(te0 + math.ceil(te_per_rad * abs(k * CO) * float(np.max(rs))))
        acc = 0.0
        for chunk in np.array_split(np.arange(rs.size), max(1, rs.size // 24)):
            Dm, _, Dxy = spheroid_pieces(k, CO, RS[chunk], dt[chunk], n_sig, n_te, n_phi)
            beta = 3.0 * M[chunk] + AV[chunk] * Dm
            mu = -RS[chunk] * Dxy / (dt[chunk] ** 2)
            A = 3.0 * (1.0 + ph[chunk]) * g0[chunk] + g1[chunk] * beta + g0[chunk] * mu
            acc += float((base[chunk] * np.real(A * PhiR[chunk])).sum())
        bulk += wa[a] * acc
    # boundary term of the age cut, in S = s_max (as F6.8b2 _boundary_term, first order)
    e = _s_edges_extended(4, 1.0)
    e = np.r_[0.0, e[(e > 1e-3) & (e < 16.0)], 16.0]
    bnd = 0.0
    for a in range(cos.size):
        CO = float(cos[a])
        for lo, hi in zip(e[:-1], e[1:]):
            S, wS = _gl(lo, hi, 32)
            r = np.sqrt(np.maximum(T * T - S / _L2, 1e-300))
            wrr = wS / (2.0 * _L2 * r)
            gS, _, _ = kernel_derivs(S * S)
            Dxy = r * np.exp(1j * k * CO * r / 2.0) * np.sinc(k * CO * r / 2.0 / math.pi)
            dsb = -2.0 * _L2 * r * Dxy
            bnd += wa[a] * float((wrr * r ** 2 * gS / (2.0 * _L2 * T) * np.real(dsb * Phi_of(k, m, r * CO))).sum())
    d = _A4 / 4.0 * (bulk + bnd)
    cont = -2.0 * math.pi * k / math.sinh(math.pi * k / (2.0 * m)) if k > 0 else -4.0 * m
    return d, cont, _A4 / 4.0 * bnd, time.time() - t0


# ================================================================================== run

def _job(kind, args):
    if kind == "omega":
        m, T, k, mult = args
        d, c = run_k(m, T, k, nr=int(round(60 * mult)), nz=int(round(40 * mult)))
        return {"d": d, "c": c}
    if kind == "nu":
        m, T, k, mult = args
        d, c, b, _ = nu_response(k, m, T, mult)
        return {"d": d, "c": c, "bnd": b}
    if kind == "stage1":
        m, T = args
        return {"d": run_stage1(m, T)[1]}
    if kind == "geom":
        return check_geometry()
    raise ValueError(kind)


def run(p, ctx):
    t0 = time.time()
    m = float(p["m"])
    ks = [float(k) for k in p["ks"]]
    T, Ta = float(p["T"]), float(p["T_age"])
    mu, mg = float(p["mult"]), float(p["mult_grid"])
    keys = []
    for sec in ("omega", "nu"):
        for k in ks:
            for tag, TT, mm in (("main", T, mu), ("age", Ta, mu), ("grid", T, mg)):
                keys.append((sec, k, tag))
    jobs = [(sec, (m, {"main": T, "age": Ta, "grid": T}[tag], k, {"main": mu, "age": mu, "grid": mg}[tag]))
            for sec, k, tag in keys]
    jobs += [("stage1", (m, T)), ("geom", ())]
    res = parallel_map(_job, jobs, ctx.workers, ctx.progress,
                       lambda i, r_, d, n: f"[F6.10b] task {d}/{n} ({time.time() - t0:.0f}s)")
    R = {key: r_ for key, r_ in zip(keys, res[:len(keys)])}
    st1, geom = res[len(keys)], res[len(keys) + 1]

    lines = [f"3+1, eps rho = 1; static wall width 1/m = {1 / m:g} l_k; ratio = d s_kin/db / continuum; judged at "
             f"T = {T:g}, grid {mu:g}; I1 age T = {Ta:g}, I2 grid {mg:g}"]
    i1 = i2 = True
    verdicts = {}
    table = {}
    for sec in ("omega", "nu"):
        rat = {}
        dl = {}
        for k in ks:
            r0 = R[(sec, k, "main")]["d"] / R[(sec, k, "main")]["c"]
            ra = R[(sec, k, "age")]["d"] / R[(sec, k, "age")]["c"]
            rg = R[(sec, k, "grid")]["d"] / R[(sec, k, "grid")]["c"]
            da, dg = abs(r0 - ra), abs(r0 - rg)
            i1 &= da <= float(p["tol_age"]) * abs(r0)
            i2 &= dg <= float(p["tol_grid"]) * abs(r0)
            rat[k], dl[k] = r0, da + dg
            extra = f"; boundary part {R[(sec, k, 'main')]['bnd']:+.1e}" if sec == "nu" else ""
            lines.append(f"  {sec:5s} k {k:<5g}: d s_kin/db {R[(sec, k, 'main')]['d']:+.5f} (continuum "
                         f"{R[(sec, k, 'main')]['c']:+.5f}) ratio {r0:.4f}; T {Ta:g}: {ra:.4f}; grid {mg:g}: {rg:.4f}; "
                         f"Delta {da + dg:.4f}{extra}")
        kj = [k for k in ks if 0 < k <= float(p["kmax_judge"])]
        mono = all(rat[a] >= rat[b] for a, b in zip(ks[:-1], ks[1:]))
        ok0 = abs(rat[0.0] - 1) + dl[0.0] <= float(p["tol0"])
        okk = all(abs(rat[k] - 1) + dl[k] <= float(p["tolk"]) for k in kj)
        bad = (abs(rat[0.0] - 1) - dl[0.0] > float(p["fail0"]) or
               any(abs(rat[k] - 1) - dl[k] > float(p["failk"]) for k in kj))
        v = "PASSED" if (ok0 and okk and mono) else "FAILED" if bad else "UNDECIDED"
        verdicts[sec] = v
        table[sec] = {str(k): {"ratio": rat[k], "delta": dl[k]} for k in ks}
        lines.append(f"  {sec}: k = 0 |ratio - 1| + Delta = {abs(rat[0.0] - 1) + dl[0.0]:.4f} (tol {p['tol0']}); "
                     f"non-increasing in k: {mono}; k <= {p['kmax_judge']}: "
                     + ", ".join(f"{abs(rat[k] - 1) + dl[k]:.4f}" for k in kj) + f" (tol {p['tolk']}) -> {v}")
    w0 = R[("omega", 0.0, "main")]["d"] if 0.0 in ks else float("nan")
    i3 = abs(w0 - st1["d"]) <= float(p["tol_i3"]) * abs(st1["d"])
    i4 = (geom["rel_beta"] <= float(p["tol_geom"]) and geom["rel_gam"] <= float(p["tol_geom"])
          and geom["unif"] <= float(p["tol_unif"]))
    inst = {"I1": bool(i1), "I2": bool(i2), "I3": bool(i3), "I4": bool(i4)}
    lines.append(f"  I1 age (1 %) -> {'OK' if i1 else 'FAILS'}; I2 grid (1 %) -> {'OK' if i2 else 'FAILS'}; I3 omega k = 0 "
                 f"{w0:+.7f} vs stage-1 kernel {st1['d']:+.7f} -> {'OK' if i3 else 'FAILS'}; I4 geometry rel "
                 f"{geom['rel_beta']:.1e} / {geom['rel_gam']:.1e}, uniform null {geom['unif']:.1e} -> "
                 f"{'OK' if i4 else 'FAILS'}")
    if not all(inst.values()):
        status = "INSTRUMENT_FAILED"
    elif all(v == "PASSED" for v in verdicts.values()):
        status = "PASSED"
    elif any(v == "FAILED" for v in verdicts.values()):
        status = "FAILED"
    else:
        status = "UNDECIDED"
    lines.append(f"cost: {time.time() - t0:.0f} s")
    return Result(status, lines, {"ratios": table, "verdicts": verdicts, "instrument": inst,
                                  "stage1": st1["d"], "geometry": geom})
