"""
F6.8f-dev step 0' (development, decides nothing): the age slope of the F6.8c2 code as a derivative, not a secant
================================================================================================================

optical_row of F6.8c (the code of the F6.8c2 card: the same quad and boundary nodes) at k and the ages T given;
J_bb, its boundary part, the interior (J_bb - bnd) and, for every two consecutive ages, the central-difference slope
-(a4 / k^4) dJ/dn (units of k^2 per unit n) of the total, the interior and the boundary. To compare with the bulk +
layer of diag.f68f_layer at the same n. Choose ages at which the node counts of the angle and r panels are equal
(16 Nov: T = 197.5, 199.9 at k = 0.12, both 72 angles and 344 r nodes).

Run from the repository root (one job = one age, cached in diag/out/f68f_c2fd_cache/):
    python -m diag.f68f_c2fd --k 0.12 --T 197.5,199.9 --workers 2
    python -m diag.f68f_c2fd --k 0.12 --T 197.5,199.9 --mult 2 --workers 2     the same at the double quadrature level
Output diag/out/f68f_c2fd_<UTC time>.json and .txt. The run of 16 Nov before this form is diag/out/f68f_c2fd.txt.
"""
from __future__ import annotations

import argparse
import os
import time

from core.stats.parallel import default_workers, parallel_map
from diag import f68f_common as Cm
from phases.F6_order_dynamics.F6_08a_bd_action_conformal import _A4
from phases.F6_order_dynamics.F6_08c_bd_action_optical import optical_row
import phases.F6_order_dynamics.F6_08c_bd_action_optical as F68c

NAME = "f68f_c2fd"
CODE = Cm.code_md5([os.path.abspath(__file__), F68c.__file__])


def row(k, T, mult):
    q = (4, 16, 0.8, 8, 16)
    return optical_row(k, T, mult, q, int(32 * mult))


def job(*args):
    return Cm.cached(NAME, CODE, row, args)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=float, default=0.12)
    ap.add_argument("--T", type=str, default="197.5,199.9")
    ap.add_argument("--mult", type=float, default=1.0)
    ap.add_argument("--workers", type=int, default=0)
    a = ap.parse_args()
    t0 = time.time()
    k, Ts = a.k, Cm.floats(a.T)
    jobs = [(k, T, a.mult) for T in Ts]
    r = parallel_map(job, jobs, default_workers(a.workers), print,
                     lambda i, x, d, n: f"[f68f_c2fd] {d}/{n}: T {jobs[i][1]:g} ({time.time() - t0:.0f}s)")
    lines = [f"F6.8f-dev: F6.8c2 code (optical_row of F6.8c), k = {k:g}, mult {a.mult:g}; development", ""]
    for x in r:
        lines.append(f"  T {x['T']:.2f} (n {x['T'] * k:.3f}): J_bb {x['J_bb']:.12e}  bnd {x['bnd']:+.6e}  interior "
                     f"{x['J_bb'] - x['bnd']:.12e}  Kbb {x['Kbb']:.9f}  n_ang {x['n_ang']} n_rs {x['n_rs']} "
                     f"n_te_max {x['n_te_max']}  {x['time']:.0f} s")
    lines.append("")
    slopes = []
    for x0, x1 in zip(r[:-1], r[1:]):
        dn = (x1["T"] - x0["T"]) * k
        tot = -_A4 * (x1["J_bb"] - x0["J_bb"]) / dn / k ** 4
        inn = -_A4 * ((x1["J_bb"] - x1["bnd"]) - (x0["J_bb"] - x0["bnd"])) / dn / k ** 4
        bnd = -_A4 * (x1["bnd"] - x0["bnd"]) / dn / k ** 4
        nm = 0.5 * (x0["T"] + x1["T"]) * k
        slopes.append({"n": nm, "total": tot, "interior": inn, "boundary": bnd})
        lines.append(f"  central-difference slope at n = {nm:.3f}: total {tot:.6f}, interior {inn:.6f}, boundary "
                     f"{bnd:+.6f} (units of k^2 per unit n)")
    lines.append(f"  cost {time.time() - t0:.0f} s")
    stem = Cm.save(NAME, vars(a), CODE, {"rows": r, "slopes": slopes}, lines, t0)
    print("\n".join(lines))
    print("saved", stem + ".json")


if __name__ == "__main__":
    main()
