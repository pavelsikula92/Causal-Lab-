"""
B0: is the larger rms of the wall under the intrinsic weight noise or a distortion of the profile? (development)
=================================================================================================================
Wall at rest, T 6, the last slab t in [5, 6), |y| < 2; bins of x of width 0.25 over |x| < 4. Per sprinkling: the mean of
phi - phi_exact per bin. Over sprinklings: the mean (bias = distortion of the profile) and the scatter of the bin means;
and the per-event rms inside a bin (noise of single events). Variants coord and B on the same sprinklings.

    python -m diag.b0_wall_profile --rho 25 --seeds 4 --workers 2
"""
from __future__ import annotations
import argparse, math, time
import numpy as np
from core.order.sprinkle import derive_seed
from core.stats.parallel import parallel_map, default_workers
from phases.F3_structure import field1d as F
from diag.b0_green2p1 import _evolve3, sprinkle_pyramid
from diag.b0_chain_alg import links3
from diag.b0_chain_field import evolve_chain, weight_table, newest_calib

EDGES = np.arange(-4.0, 4.0 + 1e-9, 0.25)


def task(rho, seed, variant, T, calib):
    t, x, y = sprinkle_pyramid(rho, T + 8.0, T + 3.0, T, seed)
    phi0 = F.kink_free(t, x, 0.0)
    if variant == "coord":
        phi = _evolve3(t, x, y, phi0, float(rho), 1.0, 0.68 * rho ** (-1 / 3))
    else:
        ptr, idx = links3(t, x, y, (12.0 * 25.0 / (math.pi * rho)) ** (1 / 3))
        phi, _, _ = evolve_chain(ptr, idx, phi0, float(rho), 1.0, weight_table("B", rho, calib), 0.0, np.zeros(1), True)
    e = phi - F.kink_exact(t, x, 0.0)
    s = (t >= T - 1) & (np.abs(y) < 2) & (np.abs(x) < 4)
    k = np.digitize(x[s], EDGES) - 1
    mean = np.array([np.mean(e[s][k == b]) for b in range(EDGES.size - 1)])
    rms = np.array([np.sqrt(np.mean(e[s][k == b] ** 2)) for b in range(EDGES.size - 1)])
    n = np.array([np.sum(k == b) for b in range(EDGES.size - 1)])
    return {"variant": variant, "seed": seed, "mean": mean, "rms": rms, "n": n}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rho", type=float, default=25.0)
    ap.add_argument("--seeds", type=int, default=4)
    ap.add_argument("--T", type=float, default=6.0)
    ap.add_argument("--workers", type=int, default=0)
    a = ap.parse_args()
    calib = newest_calib()[:3]
    jobs = [(a.rho, derive_seed(20261111, "B0", a.rho, k), v, a.T, calib) for k in range(a.seeds) for v in ("B", "coord")]
    t0 = time.time()
    res = parallel_map(task, jobs, default_workers(a.workers), print, lambda i, r, d, n: f"[profile] {d}/{n} ({time.time() - t0:.0f}s)")
    xc = 0.5 * (EDGES[1:] + EDGES[:-1])
    print(f"wall at rest, rho {a.rho:g}, last slab; per variant: |bias| (mean over sprinklings of the bin mean of phi - exact), "
          f"scatter of the bin means between sprinklings, per-event rms; summarised over the core |x| < 2 and the tails")
    for v in ("coord", "B"):
        R = [r for r in res if r["variant"] == v]
        M = np.array([r["mean"] for r in R])
        bias = M.mean(axis=0)
        scat = M.std(axis=0, ddof=1)
        rms = np.sqrt(np.mean([r["rms"] ** 2 for r in R], axis=0))
        for lab, m in (("core |x|<2", np.abs(xc) < 2), ("tails 2<|x|<4", np.abs(xc) >= 2)):
            print(f"  {v:5s} {lab}: rms of the bias {np.sqrt(np.mean(bias[m] ** 2)):.4f} (SE ~{np.sqrt(np.mean(scat[m] ** 2)) / math.sqrt(len(R)):.4f}); "
                  f"scatter of bin means {np.sqrt(np.mean(scat[m] ** 2)):.4f}; per-event rms {np.sqrt(np.mean(rms[m] ** 2)):.4f}")
        print(f"  {v:5s} bias by x: " + " ".join(f"{b:+.3f}" for b in bias[::2]))
    print(f"({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
