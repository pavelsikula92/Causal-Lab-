"""B2 development: readouts of the breather centre on the saved band fields (decides nothing).

Compares, per sprinkling, three estimators of the centre per slab on the discrete field and on the exact field at the
same events: the energy centroid over |x| < 7 and |x| < 4, and a fit of the breather profile
phi = 4 atan(A(t) / cosh(k (x - xc))), A(t) quadratic in t - tm, over |x| < 5 of the slab.
"""
from __future__ import annotations

import glob
import math
import sys

import numpy as np
from scipy.optimize import least_squares

from diag.b2_dev import e_centroid


def fit_breather(t, x, phi, a, b, guess=0.0, k0=0.714):
    s = (t >= a) & (t < b) & (np.abs(x - guess) < 5.0)
    tt, xx, pp = t[s] - 0.5 * (a + b), x[s], phi[s]
    A0 = float(np.tan(np.clip(np.mean(pp[np.abs(xx - guess) < 0.5]), -6, 6) / 4.0))

    def f(p):
        xc, lk, a0, a1, a2 = p
        return 4.0 * np.arctan((a0 + a1 * tt + a2 * tt * tt) / np.cosh(math.exp(lk) * (xx - xc))) - pp

    r = least_squares(f, [guess, math.log(k0), A0, 0.0, 0.0])
    return float(r.x[0]), float(math.exp(r.x[1])), float(r.x[2])


def main():
    pat = sys.argv[1] if len(sys.argv) > 1 else "diag/out/b2_dev_fields/r30_*br0.7.npz"
    out = {}
    for fn in sorted(glob.glob(pat)):
        d = np.load(fn)
        t, x, phi, ex = d["t"], d["x"], d["phi"], d["exact"]
        rows = []
        for a in np.arange(0.0, 6.0 + 1e-9, 1.0):
            tm = a + 0.5
            c7, _ = e_centroid(t, x, phi, tm)
            c7x, _ = e_centroid(t, x, ex, tm)
            c4, _ = e_centroid(t, x, phi, tm, half=4.0)
            c4x, _ = e_centroid(t, x, ex, tm, half=4.0)
            fb = fit_breather(t, x, phi, a, a + 1.0)
            fbx = fit_breather(t, x, ex, a, a + 1.0)
            rows.append([c7, c7x, c4, c4x, fb[0], fbx[0], fb[2], fbx[2]])
        out[fn] = np.array(rows)
    for name, i in (("energy centroid |x|<7", 0), ("energy centroid |x|<4", 2), ("breather fit", 4)):
        print(f"== {name}")
        for kind in ("24.0", "inf"):
            M = [v for f, v in out.items() if f"_{kind}_" in f]
            if not M:
                continue
            rB, last = [], []
            for v in M:
                d, dx = v[1:, i] - v[0, i], v[1:, i + 1] - v[0, i + 1]
                rB.append(np.sum(d * dx) / np.sum(dx * dx))
                last.append(d[-1] - dx[-1])
            print(f"   ell {kind}: n {len(M)}; mean per slab " + " ".join(f"{q:+.3f}" for q in np.mean([v[:, i] for v in M], 0))
                  + "; exact " + " ".join(f"{q:+.3f}" for q in np.mean([v[:, i + 1] for v in M], 0)))
            rB = np.array(rB)
            last = np.array(last)
            print(f"      r_B {rB.mean():+.3f} +- {rB.std(ddof=1) / math.sqrt(rB.size):.3f} (sd {rB.std(ddof=1):.3f}); "
                  f"d(last) - d_ex(last) {last.mean():+.3f} +- {last.std(ddof=1) / math.sqrt(last.size):.3f} (sd "
                  f"{last.std(ddof=1):.3f})" if rB.size > 1 else f"      r_B {rB}")
    np.save("diag/out/b2_breather_readout.npy", np.array([[f] for f in out]), allow_pickle=True)


if __name__ == "__main__":
    main()
