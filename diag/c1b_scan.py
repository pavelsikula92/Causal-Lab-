"""
C1b: the exact mean and the proper-time kernel for several source strengths a (15 Nov; the author's decision of 08:04:
theory first, a new run only if the predicted difference exceeds the noise). Writes diag/out/c1b_scan_rho60.npz
(f0[a], f1[a], flat) and diag/out/c1b_scan_rho60.txt.

    python -m diag.c1b_scan
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np

from diag import c1b_tilt as B

RHO = 60.0
AS = (0.025, 0.05, 0.10, 0.15)


def main(conv=True):
    out = os.path.join("diag", "out", "c1b_scan_rho60.npz")
    t0 = time.time()
    ts, rs, ff = B.mean_field(RHO, 0.0, 0)
    res = {"ts": ts, "rs": rs, "flat": ff}
    for a in AS:
        _, _, f0 = B.mean_field(RHO, a, 0)
        _, _, f1 = B.mean_field(RHO, a, 1)
        res[f"f0_{a}"], res[f"f1_{a}"] = f0, f1
        np.savez(out, **res)
        print(f"a {a}: done ({time.time() - t0:.0f}s)", flush=True)
    if conv:
        for lab, kw in (("vol 10/12/10", dict(nx=10, ne=12, nphi=10)), ("ns 36 nc 36", dict(ns=36, nc=36))):
            _, _, f0 = B.mean_field(RHO, 0.15, 0, **kw)
            _, _, f1 = B.mean_field(RHO, 0.15, 1, **kw)
            res[f"f0_0.15_{lab}"], res[f"f1_0.15_{lab}"] = f0, f1
            np.savez(out, **res)
            print(f"a 0.15 {lab}: done ({time.time() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main("--noconv" not in sys.argv)
