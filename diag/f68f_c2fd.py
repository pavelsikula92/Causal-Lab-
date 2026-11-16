"""
F6.8f-dev step 0' (development, decides nothing): the age slope of the F6.8c2 code as a derivative, not a secant.
optical_row of F6.8c (the code of F6.8c2, same quad and boundary nodes) at k = 0.12 and T = 197.5, 199.9 (the node
counts of the angle and r panels are the same at both ages). Reports J_bb, the boundary part, the interior
(J_bb - bnd) and the central-difference slope -(a4 / k^2) dJ/dn / k^2 at n = 23.844, to compare with the bulk + layer
of diag/f68f_layer.py at the same n.
Run: python diag/f68f_c2fd.py   (output diag/out/f68f_c2fd.txt)
"""
import os
import sys
import time
from multiprocessing import Pool

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from phases.F6_order_dynamics.F6_08c_bd_action_optical import optical_row  # noqa: E402
from phases.F6_order_dynamics.F6_08a_bd_action_conformal import _A4  # noqa: E402

K = 0.12
TS = [197.5, 199.9]


def row(T):
    return optical_row(K, T, 1.0, (4, 16, 0.8, 8, 16), 32)


def main():
    t0 = time.time()
    with Pool(2) as pool:
        r = pool.map(row, TS)
    dn = (TS[1] - TS[0]) * K
    lines = ["F6.8f-dev: F6.8c2 code (optical_row of F6.8c), k = 0.12, mult 1; development", ""]
    for x in r:
        lines.append(f"  T {x['T']:.2f} (n {x['T'] * K:.3f}): J_bb {x['J_bb']:.12e}  bnd {x['bnd']:+.6e}  interior "
                     f"{x['J_bb'] - x['bnd']:.12e}  Kbb {x['Kbb']:.9f}  n_ang {x['n_ang']} n_rs {x['n_rs']} "
                     f"n_te_max {x['n_te_max']}  {x['time']:.0f} s")
    sl_tot = -_A4 * (r[1]['J_bb'] - r[0]['J_bb']) / dn / K ** 4
    sl_int = -_A4 * ((r[1]['J_bb'] - r[1]['bnd']) - (r[0]['J_bb'] - r[0]['bnd'])) / dn / K ** 4
    sl_bnd = -_A4 * (r[1]['bnd'] - r[0]['bnd']) / dn / K ** 4
    lines.append("")
    lines.append(f"  central-difference slope at n = {0.5 * (TS[0] + TS[1]) * K:.3f}: total {sl_tot:.6f}, interior "
                 f"{sl_int:.6f}, boundary {sl_bnd:+.6f} (units of k^2 per unit n)")
    lines.append(f"  cost {time.time() - t0:.0f} s")
    txt = "\n".join(lines)
    print(txt)
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "f68f_c2fd.txt"), "w") as fh:
        fh.write(txt + "\n")


if __name__ == "__main__":
    main()
