"""
OV2 sec. 14, T3 pilot stage 0 (development, decides nothing): the windowed Noether identity at the ensemble level (L1).
=====================================================================================================================
2D light-cone square [0, 1]^2 (ds^2 = -du dv, proper density 2N), the pairwise action of T1 (F6.10a):
S_m = -eps^2 sum_{x<y} f(n)(phi_x - phi_y)^2 - (1/rho) sum V(phi), V = m^2 (1 - cos phi), <f> = g(eps N p q).
Field: a moving sine-Gordon kink Phi = 4 arctan(exp(m gamma (z - v t))), t = (u + v)/2, z = (v - u)/2 (a fixed function).
Window chi(t): cos^2 ramps of half-width d at t1 (up) and t2 (down).
Computed (expectations, Mecke formula, deterministic quadrature):
  B  = sum_x chi_x E_x Phi_t(x)                                       (bulk residual along xi = d_t)
  R  = sum_x sum_y chi_x (1 - chi_y) w dL/dphi_x Phi_t(x) - sum_pairs (chi'_x chi_y + chi_x chi'_y) w L - sum chi'_x U_x
split into the ramp parts Q2 (near t2) and Q1 (near t1): R = Q2 + Q1 with Q2 = Q(t2), Q1 = -Q(t1).
Expected: <B> = <R> (the identity, exact in expectation up to the edges of the square); -Q2 = H(t2) and Q1 = H(t1),
H = int T_tt dz on the slice (-> 8 m gamma).

    python -m diag.ov2_t3_l1 --m 20 --v 0,0.3 --N 66560,266240
"""
import argparse
import math
import time

import numba as nb
import numpy as np


@nb.njit(cache=True, inline="always")
def _kink(u, v, m, gam, vel):
    t = 0.5 * (u + v)
    z = 0.5 * (v - u)
    X = m * gam * (z - vel * t)
    ph = 4.0 * math.atan(math.exp(X)) if X < 30.0 else 2.0 * math.pi - 4.0 * math.atan(math.exp(-X))
    pt = -2.0 * m * gam * vel / math.cosh(X) if abs(X) < 300.0 else 0.0
    return ph, pt


@nb.njit(cache=True, inline="always")
def _chi(t, t1, t2, d):
    """window value and derivative (cos^2 ramps)."""
    if t <= t1 - d or t >= t2 + d:
        return 0.0, 0.0
    if t < t1 + d:
        x = (t - (t1 - d)) / (2.0 * d)
        return math.sin(0.5 * math.pi * x) ** 2, math.pi / (2.0 * d) * math.sin(math.pi * x) * 0.5
    if t > t2 - d:
        x = (t - (t2 - d)) / (2.0 * d)
        return math.cos(0.5 * math.pi * x) ** 2, -math.pi / (2.0 * d) * math.sin(math.pi * x) * 0.5
    return 1.0, 0.0


@nb.njit(cache=True)
def _pairs(n_exp, eps, m, gam, vel, t1, t2, d, xo, wo, xt, wt, xw, ww, pmin, wcap):
    """returns B_pair, Q2_pair, Q1_pair (pair parts; times N^2, without the potential terms)."""
    mid = 0.5 * (t1 + t2)
    B = 0.0
    Q2 = 0.0
    Q1 = 0.0
    c = eps * eps
    for i in range(xo.size):
        u1 = xo[i]
        for j in range(xo.size):
            v1 = xo[j]
            ta = 0.5 * (u1 + v1)
            ca, dca = _chi(ta, t1, t2, d)
            pa, pta = _kink(u1, v1, m, gam, vel)
            tlo = math.log(pmin)
            thi = math.log(1.0 - u1)
            if thi <= tlo:
                continue
            aB = 0.0
            a2 = 0.0
            a1 = 0.0
            for k in range(xt.size):
                tt = tlo + (thi - tlo) * xt[k]
                p = math.exp(tt)
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
                    tb = 0.5 * (u2 + v2)
                    cb, dcb = _chi(tb, t1, t2, d)
                    if ca == 0.0 and cb == 0.0 and dca == 0.0 and dcb == 0.0:
                        continue
                    pb, ptb = _kink(u2, v2, m, gam, vel)
                    g = math.exp(-w) * (1.0 - 2.0 * w + 0.5 * w * w)
                    D = pa - pb
                    L = -c * D * D
                    dLa = -2.0 * c * D        # dL/dphi_a
                    dLb = 2.0 * c * D         # dL/dphi_b
                    meas = g * dp * dq
                    # bulk residual (pair part of E): chi_a dL_a Phi_t(a) + chi_b dL_b Phi_t(b)
                    aB += meas * (ca * dLa * pta + cb * dLb * ptb)
                    # ramp terms
                    t1v = meas * ca * (1.0 - cb) * dLa * pta
                    t2v = meas * cb * (1.0 - ca) * dLb * ptb
                    t3a = -meas * dca * cb * L
                    t3b = -meas * ca * dcb * L
                    if tb > mid:
                        a2 += t1v + t3b
                    else:
                        a1 += t1v + t3b
                    if ta > mid:
                        a2 += t2v + t3a
                    else:
                        a1 += t2v + t3a
            W = n_exp * n_exp * wo[i] * wo[j]
            B += W * aB
            Q2 += W * a2
            Q1 += W * a1
    return B, Q2, Q1


def _gl01(n):
    x, w = np.polynomial.legendre.leggauss(int(n))
    return 0.5 * (x + 1), 0.5 * w


def _composite(n_pan, k=8, lo=0.0, hi=1.0):
    x8, w8 = np.polynomial.legendre.leggauss(k)
    e = np.linspace(lo, hi, n_pan + 1)
    xs = (0.5 * (e[1:] + e[:-1])[:, None] + 0.5 * (e[1:] - e[:-1])[:, None] * x8[None, :]).ravel()
    ws = (0.5 * (e[1:] - e[:-1])[:, None] * w8[None, :]).ravel()
    return xs, ws


def _potential(m, gam, vel, t1, t2, d, n=2000):
    """potential parts (exact integrals on the square): B_U = -1/2 int chi V' Phi_t, Q_U = 1/2 int chi' V per ramp;
    and the continuum energy H on the slices t1, t2."""
    g = (np.arange(n) + 0.5) / n
    U, Vv = np.meshgrid(g, g, indexing="ij")
    t = 0.5 * (U + Vv)
    z = 0.5 * (Vv - U)
    X = m * gam * (z - vel * t)
    ph = 4 * np.arctan(np.exp(np.clip(X, -700, 700)))
    pt = -2 * m * gam * vel / np.cosh(np.clip(X, -700, 700))
    pz = 2 * m * gam / np.cosh(np.clip(X, -700, 700))
    Vp = m * m * (1 - np.cos(ph))
    dV = m * m * np.sin(ph)
    ch = np.zeros_like(t)
    dch = np.zeros_like(t)
    for idx in np.ndindex(t.shape):
        ch[idx], dch[idx] = _chi(float(t[idx]), t1, t2, d)
    da = 1.0 / n ** 2
    BU = float(-0.5 * np.sum(ch * dV * pt) * da)
    mid = 0.5 * (t1 + t2)
    Q2U = float(0.5 * np.sum(np.where(t > mid, dch * Vp, 0.0)) * da)
    Q1U = float(0.5 * np.sum(np.where(t <= mid, dch * Vp, 0.0)) * da)
    H = []
    for tc in (t1, t2):
        hw = 0.5 - abs(tc - 0.5)
        zz = np.linspace(-hw, hw, 20001)
        XX = m * gam * (zz - vel * tc)
        ph_t = -2 * m * gam * vel / np.cosh(XX)
        ph_z = 2 * m * gam / np.cosh(XX)
        Vz = m * m * (1 - np.cos(4 * np.arctan(np.exp(XX))))
        H.append(float(np.trapezoid(0.5 * ph_t ** 2 + 0.5 * ph_z ** 2 + Vz, zz)))
    return BU, Q2U, Q1U, H


def run(m, vel, n_exp, eps=0.12, t1=0.4, t2=0.6, d=0.04, n_pan=48, nt=96, nw=64):
    t0 = time.time()
    gam = 1.0 / math.sqrt(1 - vel * vel)
    xo, wo = _composite(n_pan)
    xt, wt = _gl01(nt)
    xw, ww = _gl01(nw)
    Bp, Q2p, Q1p = _pairs(float(n_exp), eps, m, gam, vel, t1, t2, d, xo, wo, xt, wt, xw, ww, 1e-10 / (eps * n_exp), 60.0)
    BU, Q2U, Q1U, H = _potential(m, gam, vel, t1, t2, d)
    B, Q2, Q1 = Bp + BU, Q2p + Q2U, Q1p + Q1U
    return {"B": B, "Q2": Q2, "Q1": Q1, "R": Q2 + Q1, "H1": H[0], "H2": H[1], "H0": 8 * m * gam,
            "sec": time.time() - t0}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--m", type=float, default=20.0)
    ap.add_argument("--v", default="0,0.3")
    ap.add_argument("--N", default="66560")
    ap.add_argument("--pan", default="48")
    q = ap.parse_args()
    for N in [int(x) for x in q.N.split(",")]:
        for vel in [float(x) for x in q.v.split(",")]:
            for pan in [int(x) for x in q.pan.split(",")]:
                r = run(q.m, vel, N, n_pan=pan)
                print(f"m {q.m:g} v {vel:g} N {N} pan {pan}: -Q2 {-r['Q2']:+.4f}  Q1 {r['Q1']:+.4f}  (H slice "
                      f"{r['H2']:.4f} / {r['H1']:.4f}, 8 m gamma {r['H0']:.4f});  <B> {r['B']:+.5f}  <R> {r['R']:+.5f}  "
                      f"({r['sec']:.0f} s)", flush=True)


if __name__ == "__main__":
    main()
