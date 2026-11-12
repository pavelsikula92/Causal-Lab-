"""
B0: long test of the stability of the vacuum in 2+1 (coordinate tau, development; decides nothing)
==================================================================================================
The author (11 Nov): the ~10 % excess of the amplitude at the end of B0 stage 1 does not decide stability; measure the
time course A(t) and the growth exponent lambda_eff(t) = d/dt ln(A / A_exact) over a longer interval and for several
amplitudes, to tell exponential growth (instability) from a constant relative bias, an initial noise that settles, or a
drift. Field rule and Green function as diag/b0_green2p1.py (coordinate tau, cap 0.68 rho^-1/3); 8 random plane waves
(|k| in [0.3, 3], random directions) of total rms amplitude A0 in {0.003, 0.01, 0.03}; pyramid |x|, |y| <= T + 4 - t;
readout on |x|, |y| < 4 per slab of width 1: rms of phi, rms of the exact linear solution, rms of the difference.
    python -m diag.b0_vacuum_long --seeds 4 --T 12 --rho 50 --workers 2
"""
from __future__ import annotations
import argparse, json, math, os, time
import numpy as np
from core.order.sprinkle import derive_seed
from core.stats.parallel import parallel_map, default_workers
from diag.b0_green2p1 import _evolve3, sprinkle_pyramid, modes_free_2d


def task(rho, A0, seed, T):
    X = T + 4.0
    t, x, y = sprinkle_pyramid(rho, X, X, T, seed)
    rng = np.random.default_rng(derive_seed(seed, "modes"))
    modes = []
    for _ in range(8):
        k = rng.uniform(0.3, 3.0); th = rng.uniform(0, 2 * math.pi)
        modes.append((A0 / 2.0, k * math.cos(th), k * math.sin(th), rng.uniform(0, 2 * math.pi)))
    phi = _evolve3(t, x, y, modes_free_2d(t, x, y, modes), float(rho), 1.0, 0.68 * rho ** (-1 / 3))
    ex = np.zeros_like(x)
    for a, kx, ky, th in modes:
        ex += a * np.cos(kx * x + ky * y - math.sqrt(kx * kx + ky * ky + 1.0) * t + th)
    c = (np.abs(x) < 4.0) & (np.abs(y) < 4.0)
    A, Ae, E = [], [], []
    for s in np.arange(0.0, T - 1.0 + 1e-9, 1.0):
        m = c & (t >= s) & (t < s + 1.0)
        A.append(float(np.sqrt(np.mean(phi[m] ** 2)))); Ae.append(float(np.sqrt(np.mean(ex[m] ** 2))))
        E.append(float(np.sqrt(np.mean((phi[m] - ex[m]) ** 2))))
    return {"rho": rho, "A0": A0, "seed": seed, "N": int(t.size), "A": A, "Ae": Ae, "E": E}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=4); ap.add_argument("--T", type=float, default=12.0)
    ap.add_argument("--rho", type=float, default=50.0); ap.add_argument("--workers", type=int, default=0)
    ap.add_argument("--amps", type=str, default="0.003,0.01,0.03")
    a = ap.parse_args()
    amps = [float(v) for v in a.amps.split(",")]
    jobs = [(a.rho, A0, derive_seed(20261111, "vac", A0, k), a.T) for A0 in amps for k in range(a.seeds)]
    t0 = time.time()
    res = parallel_map(task, jobs, default_workers(a.workers), print,
                       lambda i, r, d, n: f"[vac] {d}/{n}: A0 {r['A0']} N {r['N']} ({time.time() - t0:.0f}s)")
    tm = np.arange(0.5, a.T, 1.0)
    L = [f"2+1 vacuum, rho {a.rho:g}, T {a.T:g}, {a.seeds} sprinklings per amplitude; ratio A/A_exact per slab and "
         f"lambda_eff = slope of ln(A/A_exact) over the second half"]
    for A0 in amps:
        R = [r for r in res if r["A0"] == A0]
        ratio = np.mean([np.array(r["A"]) / np.array(r["Ae"]) for r in R], axis=0)
        err = np.mean([np.array(r["E"]) / np.array(r["Ae"]) for r in R], axis=0)
        h = tm >= a.T / 2
        lam = [np.polyfit(tm[h], np.log(np.array(r["A"])[h] / np.array(r["Ae"])[h]), 1)[0] for r in R]
        lam1 = [np.polyfit(tm[~h], np.log(np.array(r["A"])[~h] / np.array(r["Ae"])[~h]), 1)[0] for r in R]
        L.append(f"--- A0 {A0}: A/A_exact " + " ".join(f"{v:.3f}" for v in ratio))
        L.append(f"  relative error " + " ".join(f"{v:.3f}" for v in err))
        L.append(f"  lambda_eff first half {np.mean(lam1):+.4f} +- {np.std(lam1, ddof=1) / math.sqrt(len(lam1)):.4f}; "
                 f"second half {np.mean(lam):+.4f} +- {np.std(lam, ddof=1) / math.sqrt(len(lam)):.4f} (per unit time)")
    print("\n".join(L))
    os.makedirs(os.path.join("diag", "out"), exist_ok=True)
    fn = os.path.join("diag", "out", f"b0_vacuum_long_{time.strftime('%Y%m%dT%H%M%S', time.gmtime())}.json")
    json.dump({"args": vars(a), "summary": L, "results": res}, open(fn, "w"))
    print(f"raw records: {fn} ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
