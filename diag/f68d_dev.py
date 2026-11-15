"""
F6.8d development (15 Nov 2026; decides nothing): where along the pair length the mixed sector departs from the
conformal one. Angle- and s-integrated bulk densities of the mixed (F6.8b2) and conformal (F6.8a) second-order action
as functions of the spatial pair length r (full s range: no age cut), and the cumulative difference in K_hat units

    D(R) = [a4 / (2 k^2)] (6 + int_0^R rho_ab dr) - [a4 16 / (3 k^2)] (1/2 + int_0^R rho_c dr),

which tends to K_hat_ab - K_hat at age ~R (without the layer and boundary terms of the cut).

    python -m diag.f68d_dev --k 0.1 --R 120
"""
from __future__ import annotations

import argparse
import math
import time

import numpy as np

from phases.F6_order_dynamics.F6_08a_bd_action_conformal import _A4, _L2, _gl
from phases.F6_order_dynamics.F6_08b2_bd_action_mixed_rev import _pair_terms, _s_nodes


def r_nodes(k, R, mult=1.0):
    out, wts = [], []
    z, wz = _gl(0.0, 1.0, int(48 * mult))
    R0 = min(4.0, R)
    out.append(R0 * z ** 2)
    wts.append(2.0 * R0 * z * wz)
    if R > R0:
        npan = max(4, int(math.ceil((R - R0) * max(k, 0.05) * mult)), int(math.ceil((R - R0) / 8.0)))
        pe = np.linspace(R0, R, npan + 1)
        for a, b in zip(pe[:-1], pe[1:]):
            x, w = _gl(a, b, int(8 * mult))
            out.append(x)
            wts.append(w)
    return np.concatenate(out), np.concatenate(wts)


def shell(k, r, mult=1.0, quad=(4, 16, 0.8, 8, 16)):
    """rho_ab(r), rho_c(r): angle- and s-integrated densities at spatial distance r (full s range)."""
    q = (int(quad[0]), int(quad[1] * mult), float(quad[2]) * mult, int(quad[3]))
    n_s = int(quad[4] * mult)
    npan = max(4, int(math.ceil(2.0 * k * r * 1.5 * mult / 8.0)))
    pe = np.linspace(-1.0, 1.0, npan + 1)
    ab = c = 0.0
    far = r >= 4.0
    for a, b in zip(pe[:-1], pe[1:]):
        cs, wc = _gl(a, b, 8)
        for CO, w in zip(cs, wc):
            s, ws = _s_nodes(np.array([1e12]), n_s, far)
            RS = np.full_like(s, r)
            dt = np.sqrt(RS * RS + s / _L2)
            base = ws * RS ** 2 / (_L2 * 2.0 * dt)
            mixed, conf, _ = _pair_terms(k, float(CO), RS, s, dt, q)
            ab += 2.0 * math.pi * w * float((base * mixed).sum())
            c += 2.0 * math.pi * w * float((base * conf).sum())
    return ab, c


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=float, default=0.1)
    ap.add_argument("--R", type=float, default=120.0)
    ap.add_argument("--mult", type=float, default=1.0)
    a = ap.parse_args()
    k = a.k
    rs, wr = r_nodes(k, a.R, a.mult)
    t0 = time.time()
    Iab, Ic = 6.0, 0.5
    print(f"k {k:g}; D(R) = K_hat_ab - K_hat accumulated over pair lengths up to R (no age cut)")
    marks = sorted(set([1, 2, 4, 8, 16, 32, 64, 128, 256, 512]) & set(range(int(a.R) + 1)))
    j = 0
    for r, w in zip(rs, wr):
        rab, rc = shell(k, float(r), a.mult)
        Iab += w * rab
        Ic += w * rc
        while j < len(marks) and r >= marks[j]:
            D = _A4 * Iab / (2 * k * k) - _A4 * 16.0 * Ic / (3 * k * k)
            print(f"  R ~ {r:7.2f}: K_ab {_A4 * Iab / (2 * k * k):+.6f}  K_c {_A4 * 16 * Ic / (3 * k * k):+.6f}  "
                  f"D {D:+.6f}  ({time.time() - t0:.0f}s)", flush=True)
            j += 1
    D = _A4 * Iab / (2 * k * k) - _A4 * 16.0 * Ic / (3 * k * k)
    print(f"  R = {a.R:g}: K_ab {_A4 * Iab / (2 * k * k):+.6f}  K_c {_A4 * 16 * Ic / (3 * k * k):+.6f}  D {D:+.6f}")


if __name__ == "__main__":
    main()
