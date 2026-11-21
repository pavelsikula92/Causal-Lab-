"""
OV2 T1 at the ensemble level by quadrature (family M-b in 2D; decides nothing).
=============================================================================
For a Poisson sprinkling, E[f(n, eps)] with n ~ Poisson(Lam) is exactly g(eps Lam), g(z) = e^{-z} (1 - 2z + z^2/2)
(the smeared kernel), so by the Mecke formula
    <S_kin>[sigma] = -eps^2 int int_{x < y} lam_x lam_y g(eps Lam(x, y)) (phi_x - phi_y)^2,
    <S_pot>[sigma] = -(1/2) int e^{2 sigma} V du dv,
with lam = N_exp e^{2 sigma} per du dv and Lam the expected count of the interval. Same model as diag/ov2_2d.py
(the kink m = 12, the bump a exp(-r^2 / 2 s^2), s = 0.12, the light-cone square). Convergence is checked by
refining every grid; the flat value is compared with the sprinkled means (an instrument check of both).

    python -m diag.ov2_quad2d --nexp 4160,16640,66560 --a 0.2
"""
import argparse
import time

import numpy as np

from diag.ov2_kernels import kin_quad, kin_quad_d1


def gl01(n):
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5 * (x + 1), 0.5 * w


def table(a, s, M):
    g = np.linspace(0, 1, M + 1)
    U, V = np.meshgrid(g, g, indexing="ij")
    E = np.exp(2 * a * np.exp(-((U - 0.5) ** 2 + (V - 0.5) ** 2) / (2 * s * s)))
    h = 1.0 / M
    Ec = 0.25 * (E[:-1, :-1] + E[1:, :-1] + E[:-1, 1:] + E[1:, 1:]) * h * h
    C = np.zeros((M + 1, M + 1))
    C[1:, 1:] = np.cumsum(np.cumsum(Ec, 0), 1)
    return C


def outer_nodes(no, n_exp, eps):
    """Nodes on [0, 1] clustered at the future edge (u -> 1), where the inner range is cut by the boundary:
    1 - u = exp(-tau), tau in [0, ln(10^3 eps n_exp)], composite Gauss-Legendre (8 nodes per panel)."""
    tmax = np.log(1e3 * eps * n_exp)
    k = max(1, no // 8)
    x8, w8 = np.polynomial.legendre.leggauss(8)
    xs, ws = [], []
    for i in range(k):
        a_, b_ = tmax * i / k, tmax * (i + 1) / k
        t = a_ + (b_ - a_) * 0.5 * (x8 + 1)
        wt = (b_ - a_) * 0.5 * w8
        xs.append(1 - np.exp(-t))
        ws.append(np.exp(-t) * wt)
    return np.concatenate(xs), np.concatenate(ws)


def s_kin(n_exp, eps, m, a, s, res, M=2048, clustered=True):
    no, nt, nw = res
    xo, wo = outer_nodes(no, n_exp, eps) if clustered else gl01(no)
    xt, wt = gl01(nt)
    xw, ww = gl01(nw)
    C = table(a, s, M)
    # pmin far below the scale 1/(eps n_exp) of the near-null pairs: below it g = 1 and the omitted part is
    # eps^2 n_exp^2 pmin int (1 - v) dphi^2, < 1e-6 here (a fixed pmin = 1e-7 was the error of the first version)
    return -eps ** 2 * kin_quad(float(n_exp), eps, m, a, s, C, M, xo, wo, xt, wt, xw, ww, 1e-10 / (eps * n_exp), 60.0)


def s_kin_d1(n_exp, eps, m, s, res):
    """d<S_kin>/da at a = 0 (analytic derivative; exact rectangle integrals of the Gaussian)."""
    no, nt, nw = res
    xo, wo = outer_nodes(no, n_exp, eps)
    xt, wt = gl01(nt)
    xw, ww = gl01(nw)
    return -eps ** 2 * kin_quad_d1(float(n_exp), eps, m, s, xo, wo, xt, wt, xw, ww, 1e-10 / (eps * n_exp), 60.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nexp", default="4160")
    ap.add_argument("--eps", type=float, default=0.12)
    ap.add_argument("--m", type=float, default=12.0)
    ap.add_argument("--a", type=float, default=0.2)
    ap.add_argument("--s", type=float, default=0.12)
    ap.add_argument("--res", default="64,48,48;96,64,64")
    ap.add_argument("--d1", action="store_true", help="the first-order derivative only")
    q = ap.parse_args()
    if q.d1:
        for ne in [int(x) for x in q.nexp.split(",")]:
            for res in q.res.split(";"):
                r = tuple(int(x) for x in res.split(","))
                t0 = time.time()
                d = s_kin_d1(ne, q.eps, q.m, q.s, r)
                print(f"N_exp {ne} res {r}: d<S_kin>/da = {d:+.5f}; first order at a = {q.a}: {q.a * d:+.5f}  "
                      f"({time.time() - t0:.0f} s)", flush=True)
        return
    for ne in [int(x) for x in q.nexp.split(",")]:
        for res in q.res.split(";"):
            r = tuple(int(x) for x in res.split(","))
            t0 = time.time()
            S0 = s_kin(ne, q.eps, q.m, 0.0, q.s, r)
            Sp = s_kin(ne, q.eps, q.m, q.a, q.s, r)
            Sm = s_kin(ne, q.eps, q.m, -q.a, q.s, r)
            print(f"N_exp {ne} res {r}: <S_kin> flat {S0:+.4f}  d(+a) {Sp - S0:+.4f}  d(-a) {Sm - S0:+.4f}  "
                  f"first order {(Sp - Sm) / 2:+.4f}  second {Sp + Sm - 2 * S0:+.4f}  ({time.time() - t0:.0f} s)", flush=True)


if __name__ == "__main__":
    main()
