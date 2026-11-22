"""
F6.8i-dev, D1 check 1 (development, decides nothing): the two-well needle volume near a cut locus against the exact
needle. docs/derivations/f68i_dev.md section 5.

Near a cut locus theta_c the first-arrival rays a and b tie. At finite q both wells count while d = |S_b - S_a| < c:

    A(theta; q) = sum_j A_j^a (2c)^j + [c > d] sum_j A_j^b 2^j (c - d)^(j + 2) / c^2,      c = q^2 / 2,

with A_j^r the q^2-series coefficients of ray r (diag/f68f_1c, evaluated around ray r itself). The check compares this
with the exact needle (diag/f68f_1b_exact: lower envelopes over all rays, no ray selected) at phases across the window.

    python -m diag.f68i_d1 --eps 12 --p 1 --thc 5.783185 --offsets -0.02,-0.005,0,0.005,0.02 --workers 2
"""
from __future__ import annotations

import argparse
import math
import time
from multiprocessing import Pool

import numpy as np

from diag import f68f_1b as B
from diag import f68f_1b_cut as CUT
from diag import f68f_1b_exact as EX
from diag import f68f_1c as SER


def series_for_ray(eps, p, th, v):
    """A_0..A_3 of the needle around the ray with start slope v (not necessarily the first arrival)."""
    saved = SER.rays

    def fixed(eps_, p_, ths):
        ths = np.atleast_1d(np.asarray(ths, float))
        vv = np.full(ths.size, float(v))
        out = B.integrate(vv, eps_, p_, ths, 4000)
        return vv, out[-1, 1], np.ones(ths.size, int)

    SER.rays = fixed
    try:
        return SER.series_coefficients(eps, p, np.array([th]))["A"][0]
    finally:
        SER.rays = saved


def two_well(Aa, Ab, d, q):
    c = q * q / 2
    A = sum(Aa[j] * (2 * c) ** j for j in range(4))
    if c > d:
        A += sum(Ab[j] * 2 ** j * (c - d) ** (j + 2) / c ** 2 for j in range(4))
    return A


def one_phase(args):
    eps, p, th = args
    vv, S = CUT.rays(eps, p, th)
    o = np.argsort(S)
    va, vb = vv[o[0]], vv[o[1]]
    d = float(S[o[1]] - S[o[0]])
    Aa = series_for_ray(eps, p, th, va)
    Ab = series_for_ray(eps, p, th, vb)
    ex = EX.phase_value(eps, p, th)
    return {"th": th, "d": d, "Aa": list(map(float, Aa)), "Ab": list(map(float, Ab)), "exact": ex["A"],
            "minPhi": ex["minPhi"], "nray": int(vv.size)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--eps", type=float, default=12.0)
    ap.add_argument("--p", type=float, default=1.0)
    ap.add_argument("--thc", type=float, default=5.783185)
    ap.add_argument("--offsets", type=str, default="-0.02,-0.005,0,0.005,0.02")
    ap.add_argument("--workers", type=int, default=2)
    q = ap.parse_args()
    t0 = time.time()
    ths = [q.thc + float(x) for x in q.offsets.split(",")]
    with Pool(q.workers) as pool:
        rows = pool.map(one_phase, [(q.eps, q.p, th) for th in ths], chunksize=1)
    print(f"eps {q.eps:g} p {q.p:g} theta_c {q.thc}: two-well series against the exact needle")
    for r in rows:
        line = f"  th {r['th']:.6f} d {r['d']:.3e} rays {r['nray']}:"
        for qq, ex in zip(EX.QS, r["exact"]):
            tw = two_well(r["Aa"], r["Ab"], r["d"], qq)
            one = two_well(r["Aa"], r["Ab"], 1e9, qq)
            line += f"  q {qq:g}: exact {ex:.7f} two-well {tw:.7f} (diff {tw - ex:+.1e}; one ray {one - ex:+.1e})"
        print(line + f"  minPhi {r['minPhi']:.1e}")
    print(f"({time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main()
