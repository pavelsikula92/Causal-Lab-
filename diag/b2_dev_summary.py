"""Summary of the B2 development rows (diag/out/b2_dev_*.json, newest or given). Decides nothing."""
from __future__ import annotations

import glob
import json
import math
import sys

import numpy as np


def col(rows, k):
    return np.array([r.get(k, float("nan")) for r in rows], float)


def ms(a):
    a = np.asarray(a, float)
    a = a[np.isfinite(a)]
    return (float(a.mean()), float(a.std(ddof=1) / math.sqrt(a.size)), float(a.std(ddof=1)), a.size) if a.size > 1 \
        else (float("nan"),) * 3 + (a.size,)


def f(m):
    return f"{m[0]:+.4f} +- {m[1]:.4f} (sd {m[2]:.4f}, n {m[3]})"


def main():
    fn = sys.argv[1] if len(sys.argv) > 1 else sorted(glob.glob("diag/out/b2_dev_2*.json"))[-1]
    d = json.load(open(fn))
    res = d["results"]
    print(fn)
    for arm, ell in sorted({(r["arm"], str(r["ell"])) for r in res}):
        R = [r for r in res if r["arm"] == arm and str(r["ell"]) == ell]
        print(f"== {arm} ell {ell}: N {np.mean([r['N'] for r in R]):.0f}, {np.mean([r['sec'] for r in R]):.0f} s")
        if arm.startswith("kink"):
            tm = col(R[0]["rows"], "t")
            D = [col(r["rows"], "xc")[-1] - col(r["rows"], "xc")[0] for r in R]
            Dx = [col(r["rows"], "xc_ex")[-1] - col(r["rows"], "xc_ex")[0] for r in R]
            sl = [np.polyfit(tm, col(r["rows"], "xc"), 1)[0] for r in R]
            slx = [np.polyfit(tm, col(r["rows"], "xc_ex"), 1)[0] for r in R]
            print(f"   displacement last - first: {f(ms(D))}; exact {np.mean(Dx):+.4f}; difference "
                  f"{f(ms(np.array(D) - np.array(Dx)))}")
            print(f"   slope - exact slope {f(ms(np.array(sl) - np.array(slx)))} (exact slope {np.mean(slx):+.4f})")
            print(f"   mean (xc - xc_ex) over slabs {f(ms([np.nanmean(col(r['rows'], 'xc') - col(r['rows'], 'xc_ex')) for r in R]))}")
            print(f"   w/w_ex - 1 {f(ms([np.nanmean(col(r['rows'], 'w')) / np.nanmean(col(r['rows'], 'w_ex')) - 1 for r in R]))}"
                  f"; Q - 1 {f(ms([np.nanmean(col(r['rows'], 'Q')) - 1 for r in R]))}")
            print("   per slab xc - xc_ex: " + " ".join(f"{v:+.3f}" for v in np.mean([col(r['rows'], 'xc') - col(r['rows'], 'xc_ex') for r in R], axis=0)))
            if ell != "inf":
                part = -(tm[-1] ** 2 - tm[0] ** 2) / (4 * float(ell))
                print(f"   (a test particle at rest would move {part:+.4f} between the first and last read slab)")
        else:
            rB, dd = [], []
            for r in R:
                c, cx = col(r["rows"], "c"), col(r["rows"], "c_ex")
                dc, dx = c - c[0], cx - cx[0]
                rB.append(np.nansum(dc[1:] * dx[1:]) / np.nansum(dx[1:] ** 2))
                dd.append(dc[-1] - dx[-1])
            print(f"   centroid per slab (mean): " + " ".join(f"{v:+.3f}" for v in np.nanmean([col(r['rows'], 'c') for r in R], axis=0)))
            print(f"   exact per slab (mean):    " + " ".join(f"{v:+.3f}" for v in np.nanmean([col(r['rows'], 'c_ex') for r in R], axis=0)))
            print(f"   r_B (projection on the exact fall) {f(ms(rB))}; last - exact last {f(ms(dd))}")
            print(f"   E / E_ex per slab: " + " ".join(f"{v:.3f}" for v in np.nanmean([col(r['rows'], 'E') / col(r['rows'], 'E_ex') for r in R], axis=0)))
            print(f"   amplitude per slab: " + " ".join(f"{v:.2f}" for v in np.nanmean([col(r['rows'], 'amp') for r in R], axis=0)))


if __name__ == "__main__":
    main()
