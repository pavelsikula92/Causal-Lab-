"""
F6.8e development (15 Nov 2026; decides nothing): the bulk pair density of the optical sector at long separations
with the stable instrument of phases/F6_order_dynamics/F6_08e_bd_action_optical_asymptotic.py.

    python -m diag.f68e_dev --k 0.12 --r 50,100,200,400,800,1600,3200,6400 --p 0
    python -m diag.f68e_dev --k 0.12 --r 800 --p 0,2,4,6,10,15,20,30
"""
from __future__ import annotations

import argparse
import time

from phases.F6_order_dynamics.F6_08e_bd_action_optical_asymptotic import pair_density


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=float, default=0.12)
    ap.add_argument("--r", type=str, default="400,800")
    ap.add_argument("--p", type=str, default="0")
    ap.add_argument("--level", type=int, default=1)
    a = ap.parse_args()
    lev = ((16, (24, 16, 16)), (24, (32, 24, 16)))[a.level - 1]
    for r in (float(v) for v in a.r.split(",")):
        for p in (float(v) for v in a.p.split(",")):
            t0 = time.time()
            b, c = pair_density(a.k, p / (a.k * r), r, lev[0], lev[1])
            print(f"k {a.k:g} r {r:7g} p {p:5.1f}: optical {b:+.6e}  conformal {c:+.6e}  ({time.time() - t0:.0f}s)",
                  flush=True)


if __name__ == "__main__":
    main()
