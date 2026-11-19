"""
OV1, 2D development (decides nothing): flat vacuum, balance (S-b) and phase nulls of the smeared BD action on causal
diamonds (docs/derivations/order_variation_OV1.md).
=====================================================================================================================
Order. A flat 2D causal diamond is the square [0, 1]^2 in light-cone coordinates (u, v); a Poisson sprinkling of
mean N_exp, relation x < y iff u_x < u_y and v_x < v_y (a 2D order, the configuration space of F6.4a/c). Units: the
expected number of elements of an interval of proper time tau is N_exp tau^2 (tau^2 = du dv of the rectangle).

Action (F6.4a convention, Surya 2012):  S = 4 eps (N - 2 eps sum_n N_n f(n, eps)),
f(n, eps) = (1 - eps)^n [1 - 2 eps n/(1 - eps) + eps^2 n (n - 1)/(2 (1 - eps)^2)], N_n = intervals with n elements
strictly inside. Pairs with |f(n)| < 1e-14 (n > n_cut) are not enumerated; their exact count is never needed.

Kernel range (frozen definitions of OV1): n_loc(eps) = the smallest n with sum_{m <= n} |f(m)| >= 0.99 sum |f|;
w = the proper time of an interval of n_loc expected elements, w^2 N_exp = n_loc. Nested diamonds with a common
future tip: the top-corner squares of side s = 1, 1/2, 1/4 of a diamond with T/w = 32, so T/w = 32, 16, 8.

Family I moves (2D): adjacent transpositions in U or in V. Each toggles exactly one relation, between two elements
whose interval is empty (a link created or removed); N and the class of 2D orders are kept; each move is an
involution. The symmetric pair of an element a in order k: H+ = swap a with its successor, H- = swap a with its
predecessor (a rank shift of a by +-1). Centred and second differences (S-b):
    D S = (S[H+] - S[H-]) / 2,     D2 S = S[H+] + S[H-] - 2 S[H].
The first and last element of an order have no symmetric pair; they are excluded and counted (one-sided).

Outputs per eps and level (mean +- SE over sprinklings):
- vacuum: <S>, N, the paired differences S(16) - S(8), S(32) - S(16), and the per-sprinkling least-squares slope c of
  S against N over the three levels (a volume term; the continuum mean of a flat interval has none);
- balance (on the level's own order): D S over elements in the interior (collar w), mean, SE, its trend with the time
  coordinate, the profile over the whole diamond; D2 S likewise (reported);
- phase nulls: ordering fraction chi, d_MM, height h (links of the longest chain), N_n (n <= 3 n_loc).

    python -m diag.ov1_2d --check                       (brute-force checks of S and of the local Delta S)
    python -m diag.ov1_2d --eps 0.12 --seeds 0-31 --levels 8,16,32 --drift 8,16 --workers 8
"""
from __future__ import annotations

import argparse
import json
import math
import os
import time

import numpy as np

from diag.ov1_kernels import (_add, _dF, _past_lists, _pre, abundances, balance, rel_count_height,  # noqa: F401
                               small_future)

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "ov1_2d")


# ----------------------------------------------------------------------------------------------- kernel
def fker(eps, nmax):
    n = np.arange(nmax + 2, dtype=float)
    return (1 - eps) ** n * (1 - 2 * eps * n / (1 - eps) + eps ** 2 * n * (n - 1) / (2 * (1 - eps) ** 2))


def n_loc_of(eps):
    a = np.abs(fker(eps, 200000))
    return int(np.argmax(np.cumsum(a) >= 0.99 * a.sum()))


def n_cut_of(eps, tol=1e-14):
    f = fker(eps, 200000)
    return int(np.nonzero(np.abs(f) > tol)[0].max())


# ----------------------------------------------------------------------------------------------- order
def sprinkle(n_exp, seed):
    rng = np.random.default_rng(seed)
    n = rng.poisson(n_exp)
    u, v = rng.uniform(0, 1, n), rng.uniform(0, 1, n)
    o = np.argsort(u)
    return u[o], v[o]


def ranks(u, v):
    """Order by u (index = u rank); vr[i] = v rank of element i; iv = inverse of vr."""
    vr = np.empty(u.size, np.int64)
    vr[np.argsort(v)] = np.arange(u.size)
    iv = np.empty(u.size, np.int64)
    iv[vr] = np.arange(u.size)
    return vr, iv


def reversed_order(vr):
    N = vr.size
    vrR = (N - 1 - vr)[::-1].copy()
    ivR = np.empty(N, np.int64)
    ivR[vrR] = np.arange(N)
    return vrR, ivR


def action(N, hist_row, f, eps):
    return 4 * eps * (N - 2 * eps * float(np.dot(hist_row, f[:hist_row.size])))


def d_mm(chi):
    g = lambda d: math.gamma(d + 1) * math.gamma(d / 2) / (2 * math.gamma(1.5 * d))
    a, b = 1.0, 6.0
    for _ in range(80):
        c = 0.5 * (a + b)
        a, b = (c, b) if g(c) > chi else (a, c)
    return 0.5 * (a + b)


# ----------------------------------------------------------------------------------------------- one sprinkling
def task(eps, seed, levels, drift_levels, n_bins=8):
    nl = n_loc_of(eps)
    nc = n_cut_of(eps)
    f = fker(eps, nc + 2)
    Tw_max = max(levels)
    n_exp = Tw_max ** 2 * nl
    wc = 1.0 / Tw_max                                # w in coordinates of the largest diamond
    t0 = time.time()
    u, v = sprinkle(n_exp, 1_000_003 * seed + int(round(eps * 1e6)))
    sides = [L / Tw_max for L in sorted(levels, reverse=True)]       # 1, 1/2, 1/4
    lev = np.zeros(u.size, np.int64)
    for k, s in enumerate(sides):
        lev[(u >= 1 - s) & (v >= 1 - s)] = k
    vr, iv = ranks(u, v)
    H = abundances(vr, iv, nc, lev, len(sides))
    out = {"eps": eps, "seed": seed, "n_loc": nl, "n_cut": nc, "w_coord": wc, "levels": {}}
    for k, s in enumerate(sides):
        Tw = int(round(s * Tw_max))
        m = (u >= 1 - s) & (v >= 1 - s)
        N = int(m.sum())
        us, vs = u[m], v[m]
        vrs, ivs = ranks(us, vs)
        R, hgt = rel_count_height(vrs)
        chi = R / (N * (N - 1) / 2)
        row = {"N": N, "S": action(N, H[k], f, eps), "chi": chi, "d_mm": d_mm(chi), "h": int(hgt - 1),
               "Nn": H[k][:3 * nl + 1].tolist()}
        if Tw in drift_levels:
            row["balance"] = _balance_stats(us, vs, 1 - s, s, wc, f, nc, eps, n_bins)
        out["levels"][str(Tw)] = row
    out["sec"] = time.time() - t0
    return out


def _balance_stats(u, v, lo, s, wc, f, nc, eps, n_bins):
    """D S and D2 S for rank shifts in U and in V of every element with both neighbours."""
    res = {}
    for name, (a, b) in (("U", (u, v)), ("V", (v, u))):
        o = np.argsort(a)
        a, b = a[o], b[o]
        vr, iv = ranks(a, b)
        vrR, ivR = reversed_order(vr)
        dS = balance(vr, iv, vrR, ivR, f, nc, eps)
        Dp, Dm = dS[1:], dS[:-1]            # element i = 1..N-2: H+ = swap (i, i+1), H- = swap (i-1, i)
        D = 0.5 * (Dp - Dm)
        D2 = Dp + Dm
        aa, bb = a[1:-1], b[1:-1]
        inner = (aa > lo + wc) & (aa < 1 - wc) & (bb > lo + wc) & (bb < 1 - wc)
        tt = ((aa - lo) + (bb - lo)) / (2 * s)                   # time from the past tip, 0..1
        prof = []
        for j in range(n_bins):
            sel = (tt >= j / n_bins) & (tt < (j + 1) / n_bins)
            prof.append([float(D[sel].mean()) if sel.any() else None, float(D2[sel].mean()) if sel.any() else None,
                         int(sel.sum())])
        x = tt[inner]
        res[name] = {"n_inner": int(inner.sum()), "one_sided": 2,
                     "D_mean": float(D[inner].mean()), "D_sd": float(D[inner].std()),
                     "D_tslope": float(np.polyfit(x, D[inner], 1)[0]) if inner.sum() > 2 else None,
                     "D2_mean": float(D2[inner].mean()), "D2_sd": float(D2[inner].std()),
                     "D_all_mean": float(D.mean()), "profile": prof}
    return res


# ----------------------------------------------------------------------------------------------- checks
def brute_S(u, v, eps):
    N = u.size
    rel = (u[:, None] < u[None, :]) & (v[:, None] < v[None, :])
    R = rel.astype(np.float64)
    cnt = R @ R                                                  # cnt[x, y] = # z with x < z < y (exact in float)
    n = np.rint(cnt[rel]).astype(np.int64)
    f = fker(eps, N + 2)
    return 4 * eps * (N - 2 * eps * f[n].sum())


def check():
    ok = True
    for eps, n_exp, seed in ((0.12, 150, 1), (0.24, 220, 2), (0.5, 120, 3), (0.06, 260, 4)):
        nc = n_cut_of(eps)
        f = fker(eps, nc + 2)
        u, v = sprinkle(n_exp, seed)
        vr, iv = ranks(u, v)
        H = abundances(vr, iv, nc, np.zeros(u.size, np.int64), 1)
        S1, S0 = action(u.size, H[0], f, eps), brute_S(u, v, eps)
        vrR, ivR = reversed_order(vr)
        dS = balance(vr, iv, vrR, ivR, f, nc, eps)
        worst = 0.0
        for p in range(u.size - 1):
            uu = u.copy()
            uu[p], uu[p + 1] = u[p + 1], u[p]                   # swap U ranks of p and p+1
            worst = max(worst, abs((brute_S(uu, v, eps) - S0) - dS[p]))
        print(f"eps {eps} N {u.size}: S sweep {S1:+.10f} brute {S0:+.10f} (diff {S1 - S0:+.1e}); "
              f"max |Delta S local - brute| over {u.size - 1} swaps {worst:.1e}")
        ok &= abs(S1 - S0) < 1e-9 and worst < 1e-9
    # larger orders: the stopping rule of the sweeps and the cut n_cut are active (n_cut < N)
    rng = np.random.default_rng(9)
    for eps, n_exp, seed in ((0.5, 2500, 5), (0.24, 1800, 6)):
        nc = n_cut_of(eps)
        f = fker(eps, nc + 2)
        u, v = sprinkle(n_exp, seed)
        vr, iv = ranks(u, v)
        H = abundances(vr, iv, nc, np.zeros(u.size, np.int64), 1)
        S1, S0 = action(u.size, H[0], f, eps), brute_S(u, v, eps)
        vrR, ivR = reversed_order(vr)
        dS = balance(vr, iv, vrR, ivR, f, nc, eps)
        worst = 0.0
        for p in rng.choice(u.size - 1, 12, replace=False):
            uu = u.copy()
            uu[p], uu[p + 1] = u[p + 1], u[p]
            worst = max(worst, abs((brute_S(uu, v, eps) - S0) - dS[p]))
        print(f"eps {eps} N {u.size} (n_cut {nc}): S sweep {S1:+.10f} brute {S0:+.10f} (diff {S1 - S0:+.1e}); "
              f"max |Delta S local - brute| over 12 swaps {worst:.1e}")
        ok &= abs(S1 - S0) < 1e-8 and worst < 1e-9
    print("checks:", "OK" if ok else "FAILED")
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--eps", type=float, default=0.12)
    ap.add_argument("--seeds", default="0-3")
    ap.add_argument("--levels", default="8,16,32")
    ap.add_argument("--drift", default="8")
    ap.add_argument("--workers", type=int, default=1)
    a = ap.parse_args()
    if a.check:
        check()
        return
    lo, _, hi = a.seeds.partition("-")
    seeds = list(range(int(lo), int(hi or lo) + 1))
    levels = [int(x) for x in a.levels.split(",")]
    drift = [int(x) for x in a.drift.split(",") if x]
    os.makedirs(OUT, exist_ok=True)
    jobs = [s for s in seeds if not os.path.exists(os.path.join(OUT, f"ov1_{a.eps:g}_{max(levels)}_{s}.json"))]
    from concurrent.futures import ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        for r in ex.map(_run, [(a.eps, s, levels, drift) for s in jobs]):
            print(f"eps {r['eps']} seed {r['seed']}: " + "  ".join(
                f"T/w {k}: N {x['N']} S {x['S']:+.3f}" for k, x in r["levels"].items()) + f"  ({r['sec']:.0f} s)",
                flush=True)


def _run(args):
    eps, s, levels, drift = args
    r = task(eps, s, levels, drift)
    f = os.path.join(OUT, f"ov1_{eps:g}_{max(levels)}_{s}.json")
    json.dump(r, open(f + ".tmp", "w"))
    os.replace(f + ".tmp", f)
    return r


if __name__ == "__main__":
    main()
