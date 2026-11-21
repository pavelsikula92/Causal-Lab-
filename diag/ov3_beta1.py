"""
OV3 section 10, stage beta1 (development, decides nothing; author_decisions 48): the conformal sector of the expected
BD action at CUBIC order, for the two-mode static field u(z) = A cos kz + C cos 2kz, at fixed age n = kT.
====================================================================================================================

The conformal metric g = e^{2u} eta does not change the causal order; only sqrt(-g) = e^{4u} = e^h changes. As in
F6.8a (its quadrature, kernel and closed-form diamond transform are imported unchanged), the action density is

    D(u) = a_4 [ e^{h_x} + int_past base(y) e^{h_x} e^{h_y} g_k(w <e^h>_diamond) ],     eps rho = 1,

with w the flat interval volume and <.> the average over the causal diamond from y to x. Expanded to third order in h
(g_n = w^n g_k^(n); X = h_x + h_y; m_n = <h^n>):

    order 2:  X^2/2 g0 + X m1 g1 + m2/2 g1 + m1^2/2 g2                    (+ h_x^2/2 from the first term)
    order 3:  X^3/6 g0 + X^2/2 m1 g1 + X (m2/2 g1 + m1^2/2 g2) + m3/6 g1 + m1 m2/2 g2 + m1^3/6 g3   (+ h_x^3/6)

The diamond averages of powers of u are exact: u(theta)^n is a trigonometric polynomial in theta = k (xi + o1); its
Fourier coefficients (FFT on 32 points, exact for degree <= 6) times the diamond transform M_{mk} = <e^{i m k o1}>.
The phase xi of x is averaged by 16 equally spaced samples over one period of k (exact for harmonics below 16 k).

EH reference, the normalisation of F6.8a (density 3 e^{2u} (grad u)^2; K_hat = 1 for EH):
- quadratic: 3 <u'^2> = (3/2) A^2 k^2 + 6 C^2 k^2;
- cubic: 6 <u u'^2>, whose A^2 C coefficient is 6 (3/4) k^2 = 4.5 k^2.

Reported:
- K2_A = the A^2 coefficient / (1.5 k^2) (equals F6.8a's K_hat at the same k, T: a check of the code);
- K3_c = the A^2 C coefficient / (4.5 k^2) (EH 1);
- the A^3 and C^3 coefficients (zero by the phase average; reported as a check).

    python -m diag.ov3_beta1 --ks 0.06,0.12 --n 4 [--mult 1] [--check]
"""
from __future__ import annotations

import argparse
import math
import time

import numpy as np

from phases.F6_order_dynamics.F6_08a_bd_action_conformal import (_A4, _L2, D, angle_nodes, kernel_derivs,
                                                                   mean_exp_diamond, plane_wave_row, rs_nodes,
                                                                   s_nodes)

NTH, NXI, MMAX = 32, 16, 6


def _g3(w):
    """w^3 g_k'''(w) for g_k = e^-w p(w), p the BD polynomial of 3+1."""
    from core.operators.bd import BD_COEFFS
    p = np.polynomial.Polynomial([c / math.factorial(i) for i, c in enumerate(BD_COEFFS[D])])
    p1, p2, p3 = p.deriv(1), p.deriv(2), p.deriv(3)
    e = np.exp(-w)
    return w ** 3 * e * (p3(w) - 3.0 * p2(w) + 3.0 * p1(w) - p(w))


def _past(T, kres, mult):
    rs, wr = rs_nodes(T, kres, mult)
    s, ws = s_nodes(_L2 * np.maximum(T * T - rs * rs, 0.0), int(16 * mult))
    RS = rs[:, None]
    dt = np.sqrt(RS * RS + s / _L2)
    base = wr[:, None] * ws * RS ** 2 / (_L2 * 2.0 * dt)
    w = s ** (D / 2.0)
    g0, g1, g2 = kernel_derivs(w)
    return RS, dt, base, g0, g1, g2, _g3(w)


def _coeffs(A, C):
    """Fourier coefficients c[n][m] (m = -MMAX..MMAX) of u(theta)^n, n = 1, 2, 3, u = A cos theta + C cos 2 theta."""
    th = 2.0 * math.pi * np.arange(NTH) / NTH
    u = A * np.cos(th) + C * np.cos(2.0 * th)
    out = {}
    for n in (1, 2, 3):
        F = np.fft.fft(u ** n) / NTH
        out[n] = {m: F[m % NTH] for m in range(-MMAX, MMAX + 1)}
    return out


def row(k, T, mult=1.0, amps=((1.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, -1.0))):
    """Second- and third-order parts of the phase-averaged action density for each (A, C) in amps."""
    t0 = time.time()
    RS, dt, base, g0, g1, g2, g3 = _past(T, 3.0 * k, mult)
    cos, wa = angle_nodes(3.0 * k, T, mult)
    xis = 2.0 * math.pi / k * np.arange(NXI) / NXI
    cf = {a: _coeffs(*a) for a in amps}
    o2 = {a: 0.0 for a in amps}
    o3 = {a: 0.0 for a in amps}
    for i in range(cos.size):
        CO = float(cos[i])
        M = {m: mean_exp_diamond(m * k, RS, CO, dt) for m in range(1, MMAX + 1)}
        M[0] = np.ones_like(RS * dt)
        for m in range(1, MMAX + 1):
            M[-m] = np.conj(M[m])
        delta = RS * CO
        for a in amps:
            A, C = a
            acc2 = np.zeros_like(base)
            acc3 = np.zeros_like(base)
            for xi in xis:
                hx = 4.0 * (A * math.cos(k * xi) + C * math.cos(2.0 * k * xi))
                hy = 4.0 * (A * np.cos(k * (xi + delta)) + C * np.cos(2.0 * k * (xi + delta)))
                mm = {}
                for n in (1, 2, 3):
                    z = sum(cf[a][n][m] * np.exp(1j * m * k * xi) * M[m] for m in range(-MMAX, MMAX + 1))
                    mm[n] = 4.0 ** n * np.real(z)
                X = hx + hy
                m1, m2, m3 = mm[1], mm[2], mm[3]
                acc2 += X * X / 2 * g0 + X * m1 * g1 + m2 / 2 * g1 + m1 * m1 / 2 * g2
                acc3 += (X ** 3 / 6 * g0 + X * X / 2 * m1 * g1 + X * (m2 / 2 * g1 + m1 * m1 / 2 * g2) + m3 / 6 * g1
                         + m1 * m2 / 2 * g2 + m1 ** 3 / 6 * g3)
            o2[a] += wa[i] * float((base * acc2).sum()) / NXI
            o3[a] += wa[i] * float((base * acc3).sum()) / NXI
    for a in amps:
        A, C = a
        hx = 4.0 * (A * np.cos(k * xis) + C * np.cos(2.0 * k * xis))
        o2[a] = _A4 * (o2[a] + float(np.mean(hx ** 2)) / 2.0)
        o3[a] = _A4 * (o3[a] + float(np.mean(hx ** 3)) / 6.0)
    c30, c03 = o3[(1.0, 0.0)], o3[(0.0, 1.0)]
    c21 = 0.5 * (o3[(1.0, 1.0)] - o3[(1.0, -1.0)]) - c03
    c12 = 0.5 * (o3[(1.0, 1.0)] + o3[(1.0, -1.0)]) - c30
    return {"k": k, "T": T, "mult": mult, "K2_A": o2[(1.0, 0.0)] / (1.5 * k * k), "K2_C": o2[(0.0, 1.0)] / (6.0 * k * k),
            "K3_c": c21 / (4.5 * k * k), "c30": c30, "c03": c03, "c12": c12, "c21": c21, "n_ang": int(cos.size),
            "n_rs": int(RS.shape[0]), "time": time.time() - t0}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ks", default="0.06,0.12")
    ap.add_argument("--n", type=float, default=4.0)
    ap.add_argument("--mult", type=float, default=1.0)
    ap.add_argument("--check", action="store_true", help="also F6.8a's plane_wave_row at the same (k, T)")
    q = ap.parse_args()
    for k in [float(x) for x in q.ks.split(",")]:
        T = q.n / k
        r = row(k, T, q.mult)
        line = (f"k {k:g} T {T:.2f} mult {q.mult:g}: K2_A {r['K2_A']:.6f}  K2_C {r['K2_C']:.6f}  K3_c {r['K3_c']:.6f}  "
                f"c30 {r['c30']:+.2e} c03 {r['c03']:+.2e} c12/k2 {r['c12'] / k / k:+.4f}  ({r['n_ang']} ang, "
                f"{r['n_rs']} rs, {r['time']:.0f} s)")
        if q.check:
            p = plane_wave_row(k, T, q.mult)
            line += f"  | F6.8a K_hat {p['Khat']:.6f}"
        print(line, flush=True)


if __name__ == "__main__":
    main()
