"""
OV2 T2 in 3+1, stage 1: the uniform conformal sector (development, decides nothing).
==================================================================================
A static wall phi(z) = 4 arctan(e^{m z}), V = m^2 (1 - cos phi), on a flat 3+1 order of age T, units eps rho = 1
(l_k = 1). The pairwise (Dirichlet) kinetic action of OV2 sec. 3.1, by the Mecke formula, per unit time and area:
    s_kin = (a_4 / 4) int dz int_{J^-(x), age T} g(eps rho V) (phi(y) - phi(x))^2 dV,   a_4 = 4/sqrt(6),
(it tends to -(1/2) int (d phi)^2 = -(1/2) int phi'^2 dz = -4 m). A uniform conformal factor Omega = e^b multiplies the
density by Omega^4 and every interval volume by Omega^4; d s_kin/db at b = 0 must tend to the kinetic part of
int T^a_a: -int phi'^2 dz = -8 m (the potential part, -4 int V = -16 m, is the measure, exact as in 2D).

    python -m diag.ov2_t2 --m 0.1 --T 20,40,80
"""
import argparse
import math
import time

import numpy as np

from diag.ov2_kernels4 import wall_kin

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


def run(m, T, nr=60, nz=40):
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--m", default="0.1")
    ap.add_argument("--T", default="20,40")
    ap.add_argument("--nr", default="60")
    q = ap.parse_args()
    for m in [float(x) for x in q.m.split(",")]:
        for T in [float(x) for x in q.T.split(",")]:
            for nr in [int(x) for x in q.nr.split(",")]:
                t0 = time.time()
                o = run(m, T, nr=nr)
                print(f"m {m} (width {1 / m:g} l_k) T {T:g} nr {nr}: s_kin {o[0]:+.5f} (continuum {-4 * m:+.5f}, ratio "
                      f"{o[0] / (-4 * m):.4f});  d s_kin/db {o[1]:+.5f} (continuum {-8 * m:+.5f}, ratio "
                      f"{o[1] / (-8 * m):.4f})  ({time.time() - t0:.0f} s)", flush=True)


if __name__ == "__main__":
    main()
