"""
OV2 T2 in 3+1, stage 3: the cone sector nu (development, decides nothing).
==========================================================================
g = -dt^2 + e^{2 b h} dx^2, h = cos(k z) (the optical index of F6.8b2; the order changes). A static wall
phi(z) = 4 arctan(e^{m z}) at fixed coordinates; the pairwise kinetic action of OV2 sec. 3.1 by the Mecke formula,
eps rho = 1, age T, per unit time and area:
    s_kin = (a_4/4) int dz_x int_{y < x} P_x P_y g(eps rho V) (phi_y - phi_x)^2 d^4y,   P = sqrt(-g) = e^{3 b h}.
First order in b, with the first-order cones of F6.8b2 (the order: y < x iff dt >= T_xy = L + b delta_xy, the
optical length of the straight segment; the past parametrised by s = l2 (dt^2 - T_xy^2), so the measure
1/(2 l2 dt) and the age cut s_max = l2 (T^2 - T_xy^2) depend on b; the interval volume V by the prolate quadrature
spheroid_pieces of F6.8b2). Amplitudes relative to e^{i k z_x} (complex, as in F6.8b2) are combined with the wall
through Phi(Delta) = int dz e^{i k z} (phi(z + Delta) - phi(z))^2, Delta = z_y - z_x = RS cos:
    bulk  = int dOmega RS^2 dRS ds/(2 l2 dt) Re[(3 (1 + e^{ik Delta}) g + g1 beta + g mu) Phi(Delta)],
            beta = dV/db / V0 = 3 M + (A/V0) Dm, mu = d ln m/db = -RS D_xy / dt^2, g1 = w g'(w);
    bound = int dOmega r^2 dr [g(S^2) / (2 l2 T)] Re[(-2 l2 r D_xy) Phi(r cos)]   at s = s_max (the age cut).
The continuum kinetic part of (1/2) int T^{ab} delta g_ab for delta g_ij = 2 b h delta_ij is -(1/2) int h phi'^2 dz
= -2 pi k / sinh(pi k / 2m) (k -> 0: -4 m); the potential part -3 int h V is the measure, exact.

    python -m diag.ov2_t2nu --m 0.1 --T 80 --k 0,0.05,0.1
"""
import argparse
import math
import time

import numpy as np

from phases.F6_order_dynamics.F6_08a_bd_action_conformal import (_A4, _L2, _gl, angle_nodes, kernel_derivs,
                                                                   mean_exp_diamond, rs_nodes, s_nodes)
from phases.F6_order_dynamics.F6_08b2_bd_action_mixed_rev import spheroid_pieces
from core.operators.bd_expect_curved import _s_edges_extended


def phi_w(z, m):
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
        d = phi_w(z + D, m) - phi_w(z, m)
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--m", type=float, default=0.1)
    ap.add_argument("--T", default="80")
    ap.add_argument("--k", default="0")
    ap.add_argument("--mult", default="1")
    q = ap.parse_args()
    for T in [float(x) for x in q.T.split(",")]:
        for k in [float(x) for x in q.k.split(",")]:
            for mult in [float(x) for x in q.mult.split(",")]:
                d, c, b, sec = nu_response(k, q.m, T, mult)
                print(f"m {q.m} T {T:g} k {k:g} mult {mult:g}: d s_kin/db (nu) {d:+.5f} (continuum {c:+.5f}, ratio "
                      f"{d / c:.4f}; boundary part {b:+.2e})  ({sec:.0f} s)", flush=True)


if __name__ == "__main__":
    main()
