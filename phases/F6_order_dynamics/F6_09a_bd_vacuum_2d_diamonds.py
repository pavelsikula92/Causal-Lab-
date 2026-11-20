"""
F6.9a -- No volume term in the smeared BD action of flat 2D causal diamonds (the vacuum control of OV1)
======================================================================================================

Why this test exists (20 Nov 2026, the author's decision; card docs/OV1_vacuum_2d_card_draft.md, decision 38)
-------------------------------------------------------------------------------------------------------------
OV1 (docs/derivations/order_variation_OV1.md) specifies the admissible variations of the order and the variational
rule of the BD action, without matter. Its vacuum control (section 8) is blocking: Lorentz invariance allows the
linear response delta<S> = c(T) int tr h dV, and the flat order is stationary only if c(T) = 0. In 2D the readout is
the volume term of the ensemble mean over nested flat diamonds with a common future tip, T/w = 8, 16, 32 (volumes
1 : 4 : 16), at eps = 0.06, 0.12, 0.24 (the author, 19 Nov). The continuum expects a constant, the joint term of the
two light-cone boundaries (Buck, Dowker, Jubb, Surya, arXiv:1502.05388); the constant is reported, not judged.

Model (as diag/ov1_2d.py, development and validation 19-20 Nov)
-----------------------------------------------------------------
- Order: a Poisson sprinkling of the light-cone square [0, 1]^2 of mean N_exp = 1024 n_loc(eps), x < y iff u_x < u_y
  and v_x < v_y. n_loc(eps) = the smallest n with sum_{m <= n} |f(m)| >= 0.99 sum |f|; w^2 N_exp = n_loc (OV1 sec. 4).
  Nested diamonds = the top-corner squares of side 1, 1/2, 1/4 (T/w = 32, 16, 8).
- Action (F6.4a convention): S = 4 eps (N - 2 eps sum_n N_n f(n, eps)), N_n counted exactly for n <= n_cut
  (|f| >= 1e-14), by sweeps with an exact stopping rule.
- Readout per sprinkling: c = the least-squares slope of S against N over the three levels; per eps: c/(4 eps), the
  volume term per element in units of the bare term, mean +- SE over sprinklings.

Criterion (frozen before the run; decision 38)
----------------------------------------------
Per eps: PASSED if |c/(4 eps)| + 2 SE <= 1e-3; FAILED if |c/(4 eps)| - 2 SE > 1e-3; UNDECIDED otherwise.
The test is PASSED if all three eps pass, FAILED if any fails, UNDECIDED otherwise. Sprinklings: 5400 / 1200 / 280 at
eps 0.24 / 0.12 / 0.06 (2 SE <= tol/3 by the per-sprinkling sd of the validation rows); seeds
derive_seed(base, "F6.9a", eps, k), independent of the development, validation and calibration seeds.
Instrument (INSTRUMENT_FAILED if any fails):
  I1 S and every Delta S of an adjacent transposition against brute force (orders of N = 100-300 and two with
     n_cut < N): |difference| <= 1e-8;
  I2 the mean number of elements per level within 3 SE of N_exp s^2;
  I3 the phase diagnostics (chi, d_MM, h) inside the calibration bands of OV1 sec. 13 (seeds >= 100000, frozen here):
     at most 5 % of the sprinklings outside at any eps, level and diagnostic (nominal 1 %; with the 256-512
     calibration seeds per eps the realised rate in the validation rows was 0.4-2.0 %, and 5 % keeps the false alarm
     below 1e-3 per item at 280 sprinklings, while an escape from the phase moves every sprinkling out);
  I4 time-reversal null: the mean centred difference <D S> of the interior elements at T/w = 8 (all sprinklings)
     |x| + 2 SE <= 1e-3 in both orders (guaranteed by the symmetry of the flat diamond: a code check, not physics;
     with all sprinklings 2 SE <= 1.4e-4 at every eps by the validation rows).
Reported, not judged: <S> per level, the paired differences, the constant <S> at T/w = 8, <D2 S>.

Reading declared now: PASSED = the BD action of the flat 2D order carries no volume term at the three eps; OV1's
vacuum control holds in 2D. FAILED at some eps = the construction fails the vacuum control there; the term is
explained first and not subtracted. A PASS says nothing about the stability of the flat order against changes of the
order (the phase sampler P1: it is not stable above beta_e(N)), nor about 2+1 or 3+1.
"""

from __future__ import annotations

import math
import time

import numba as nb
import numpy as np

from core.order.sprinkle import derive_seed
from core.stats.parallel import parallel_map
from core.testing import Result

BANDS = {  # calibration bands [0.5 %, 99.5 %] (OV1 sec. 13; diag/out/ov1_2d/summary.txt md5 91afb053), per eps, T/w
    0.06: {8: {"chi": (0.49100, 0.50931), "d_mm": (1.97536, 2.02419), "h": (167.55, 188.0)},
           16: {"chi": (0.49614, 0.50479), "d_mm": (1.98729, 2.01032), "h": (348.0, 375.0)},
           32: {"chi": (0.49769, 0.50207), "d_mm": (1.99449, 2.00617), "h": (715.275, 749.625)}},
    0.12: {8: {"chi": (0.48655, 0.51190), "d_mm": (1.96859, 2.03630), "h": (113.0, 130.0)},
           16: {"chi": (0.49380, 0.50630), "d_mm": (1.98328, 2.01662), "h": (238.555, 263.0)},
           32: {"chi": (0.49671, 0.50312), "d_mm": (1.99170, 2.00880), "h": (490.555, 519.0)}},
    0.24: {8: {"chi": (0.47908, 0.52087), "d_mm": (1.94531, 2.05681), "h": (74.0, 89.89)},
           16: {"chi": (0.48990, 0.50902), "d_mm": (1.97613, 2.02716), "h": (158.0, 180.0)},
           32: {"chi": (0.49567, 0.50468), "d_mm": (1.98756, 2.01160), "h": (328.555, 354.445)}},
}

TEST = {
    "id": "F6.9a",
    "phase": "F6",
    "title": "No volume term in the smeared BD action of flat 2D causal diamonds (OV1 vacuum control)",
    "theory_part": "OV1 (docs/derivations/order_variation_OV1.md): the variational rule of the order without matter; "
                   "the blocking vacuum control c(T) = 0",
    "question": "Is the volume term of the mean smeared BD action of flat 2D causal diamonds (nested, T/w = 8, 16, "
                "32) zero to 0.1 % per element at eps = 0.06, 0.12 and 0.24?",
    "method": "Poisson sprinklings of the light-cone square, three nested top-corner diamonds; exact interval "
              "abundances (n <= n_cut); S = 4 eps (N - 2 eps sum N_n f); per sprinkling the slope of S against N; "
              "5400 / 1200 / 280 sprinklings; instruments: brute force, Poisson counts, calibrated phase bands, "
              "time-reversal null",
    "prediction": "c = 0 (continuum: the mean action of a flat interval tends to its joint term, a constant)",
    "criterion": "per eps PASSED if |c/(4 eps)| + 2 SE <= 1e-3, FAILED if |c/(4 eps)| - 2 SE > 1e-3, else "
                 "UNDECIDED; all three PASSED -> PASSED, any FAILED -> FAILED. INSTRUMENT_FAILED if I1 (brute force "
                 "1e-8), I2 (N within 3 SE), I3 (<= 5 % outside the calibration bands) or I4 (time-reversal null "
                 "|<D S>| + 2 SE <= 1e-3, all sprinklings) fails",
    "null_model": "time reversal of the flat diamond (I4); flat calibration bands (I3)",
    "level": "mechanism",
    "meters": ["F6.4a"],
    "depends_on": [],
    "params": {"eps": [0.24, 0.12, 0.06], "n": {"0.24": 5400, "0.12": 1200, "0.06": 280}, "levels": [8, 16, 32],
               "tol": 1e-3, "n_se": 2.0, "n_balance": 100000, "tol_balance": 1e-3, "band_frac": 0.05,
               "brute_tol": 1e-8},
    "parts": {"main": {}},
    "quick": {"n": {"0.24": 6, "0.12": 4, "0.06": 2}},
    "estimate": "~40 core-hours (about 5 h on 8 workers); eps 0.06 about 180 s, 0.12 35 s, 0.24 9 s per sprinkling",
}


# ============================================================================================== kernels
# copied unchanged from diag/ov1_kernels.py (checked against brute force there and in I1 below)

@nb.njit(cache=True, inline="always")
def _add(tree, i, d):
    i += 1
    while i < tree.size:
        tree[i] += d
        i += i & (-i)


@nb.njit(cache=True, inline="always")
def _pre(tree, i):          # sum over [0, i]
    s = 0
    i += 1
    while i > 0:
        s += tree[i]
        i -= i & (-i)
    return s


@nb.njit(cache=True)
def small_future(x, vr, iv, ncut, tree, undo, oy, on):
    """All y > x with n(x, y) <= ncut; exact. Writes (y, n) into oy, on; returns the count. tree must be zero."""
    N = vr.size
    vx = vr[x]
    m = 0
    # case A: dv >= du (in ranks), sweep in u
    nu = 0
    cs = 0
    d = 1
    while True:
        i = x + d
        if i >= N:
            break
        if d >= 2:                       # grow the square (x, x+d) x (vx, vx+d)
            j = x + d - 1
            if vx < vr[j] < vx + d:
                cs += 1
            r = vx + d - 1
            if r > vx and r < N:
                jj = iv[r]
                if x < jj < x + d - 1:
                    cs += 1
        if cs > ncut:
            break
        if vr[i] > vx and vr[i] - vx >= d:
            n = _pre(tree, vr[i] - 1) - _pre(tree, vx)
            if n <= ncut:
                oy[m] = i
                on[m] = n
                m += 1
        if vr[i] > vx:
            _add(tree, vr[i], 1)
            undo[nu] = vr[i]
            nu += 1
        d += 1
    for k in range(nu):
        _add(tree, undo[k], -1)
    # case B: du > dv, sweep in v
    nu = 0
    cs = 0
    d = 1
    while True:
        r = vx + d
        if r >= N:
            break
        if d >= 2:
            j = x + d - 1
            if j < N and vx < vr[j] < vx + d - 1:
                cs += 1
            jj = iv[vx + d - 1]
            if x < jj < x + d:
                cs += 1
        if cs > ncut:
            break
        j = iv[r]
        if j > x and j - x > d:
            n = _pre(tree, j - 1) - _pre(tree, x)
            if n <= ncut:
                oy[m] = j
                on[m] = n
                m += 1
        if j > x:
            _add(tree, j, 1)
            undo[nu] = j
            nu += 1
        d += 1
    for k in range(nu):
        _add(tree, undo[k], -1)
    return m


@nb.njit(cache=True)
def abundances(vr, iv, ncut, lev_of, nlev):
    """N_n (n <= ncut) per nested level; a pair belongs to every level that contains its lower element."""
    N = vr.size
    H = np.zeros((nlev, ncut + 1), np.int64)
    tree = np.zeros(N + 1, np.int64)
    undo = np.empty(N, np.int64)
    oy = np.empty(N, np.int64)
    on = np.empty(N, np.int64)
    for x in range(N):
        m = small_future(x, vr, iv, ncut, tree, undo, oy, on)
        for k in range(m):
            for L in range(lev_of[x] + 1):
                H[L, on[k]] += 1
    return H


@nb.njit(cache=True)
def rel_count_height(vr):
    """Number of relations and the longest chain (elements)."""
    N = vr.size
    tree = np.zeros(N + 1, np.int64)
    R = 0
    for i in range(N):
        R += _pre(tree, vr[i] - 1) if vr[i] > 0 else 0
        _add(tree, vr[i], 1)
    tails = np.empty(N, np.int64)
    L = 0
    for i in range(N):
        a, b = 0, L
        while a < b:
            c = (a + b) // 2
            if tails[c] < vr[i]:
                a = c + 1
            else:
                b = c
        tails[a] = vr[i]
        if a == L:
            L += 1
    return R, L


@nb.njit(cache=True)
def _past_lists(x, vrR, ivR, ncut, tree, undo, oy, on):
    N = vrR.size
    m = small_future(N - 1 - x, vrR, ivR, ncut, tree, undo, oy, on)
    for k in range(m):
        oy[k] = N - 1 - oy[k]
    return m


@nb.njit(cache=True)
def _dF(p, q, vr, f, fy, fn, mf, py, pn, mp, gy, gn, mg, qy, qn, mq):
    """Change of sum_pairs f(n) when the U-adjacent pair (p, q = p+1) is swapped.
    (fy, fn, mf) future of p, (py, pn, mp) past of p, (gy, gn, mg) future of q, (qy, qn, mq) past of q."""
    nf = f.size - 2
    if vr[q] > vr[p]:                    # link p < q removed
        dF = -f[0]
        for k in range(mf):
            y = fy[k]
            if y != q and vr[y] > vr[q]:
                n = fn[k]
                dF += f[n - 1] - f[n]
        for k in range(mq):
            x = qy[k]
            if x != p and vr[x] < vr[p]:
                n = qn[k]
                dF += f[n - 1] - f[n]
    else:                                # link q < p added
        dF = f[0]
        for k in range(mg):
            y = gy[k]
            if vr[y] > vr[p]:
                n = gn[k]
                if n + 1 <= nf:
                    dF += f[n + 1] - f[n]
        for k in range(mp):
            x = py[k]
            if vr[x] < vr[q]:
                n = pn[k]
                if n + 1 <= nf:
                    dF += f[n + 1] - f[n]
    return dF


@nb.njit(cache=True)
def balance(vr, iv, vrR, ivR, f, ncut, eps):
    """Delta S of the swap of every U-adjacent pair (p, p+1); returns dS[p] for p = 0..N-2."""
    N = vr.size
    tree = np.zeros(N + 1, np.int64)
    undo = np.empty(N, np.int64)
    A = [np.empty(N, np.int64) for _ in range(8)]
    fy, fn, py, pn, gy, gn, qy, qn = A[0], A[1], A[2], A[3], A[4], A[5], A[6], A[7]
    dS = np.empty(N - 1)
    mf = small_future(0, vr, iv, ncut, tree, undo, fy, fn)
    mp = _past_lists(0, vrR, ivR, ncut, tree, undo, py, pn)
    for p in range(N - 1):
        q = p + 1
        mg = small_future(q, vr, iv, ncut, tree, undo, gy, gn)
        mq = _past_lists(q, vrR, ivR, ncut, tree, undo, qy, qn)
        dS[p] = -8.0 * eps * eps * _dF(p, q, vr, f, fy, fn, mf, py, pn, mp, gy, gn, mg, qy, qn, mq)
        fy, gy = gy, fy
        fn, gn = gn, fn
        py, qy = qy, py
        pn, qn = qn, pn
        mf, mp = mg, mq
    return dS


# ============================================================================================== order, action

def fker(eps, nmax):
    n = np.arange(nmax + 2, dtype=float)
    return (1 - eps) ** n * (1 - 2 * eps * n / (1 - eps) + eps ** 2 * n * (n - 1) / (2 * (1 - eps) ** 2))


def n_loc_of(eps):
    a = np.abs(fker(eps, 200000))
    return int(np.argmax(np.cumsum(a) >= 0.99 * a.sum()))


def n_cut_of(eps, tol=1e-14):
    f = fker(eps, 200000)
    return int(np.nonzero(np.abs(f) > tol)[0].max())


def ranks(u, v):
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


def d_mm(chi):
    g = lambda d: math.exp(math.lgamma(d + 1) + math.lgamma(d / 2) - math.log(2.0) - math.lgamma(1.5 * d))
    if not 0 < chi < 1:
        return float("nan")
    lo, hi = 0.5, 12.0
    for _ in range(80):
        m = 0.5 * (lo + hi)
        lo, hi = (m, hi) if g(m) > chi else (lo, m)
    return 0.5 * (lo + hi)


def sprinkle(n_exp, seed):
    rng = np.random.default_rng(seed)
    n = rng.poisson(n_exp)
    u, v = rng.uniform(0, 1, n), rng.uniform(0, 1, n)
    o = np.argsort(u)
    return u[o], v[o]


def _balance_inner(u, v, lo, wc, f, nc, eps):
    """<D S> and <D2 S> over the interior elements (collar w), for rank shifts in U and in V."""
    out = {}
    for name, (a, b) in (("U", (u, v)), ("V", (v, u))):
        o = np.argsort(a)
        a, b = a[o], b[o]
        vr, iv = ranks(a, b)
        vrR, ivR = reversed_order(vr)
        dS = balance(vr, iv, vrR, ivR, f, nc, eps)
        D = 0.5 * (dS[1:] - dS[:-1])
        D2 = dS[1:] + dS[:-1]
        aa, bb = a[1:-1], b[1:-1]
        inner = (aa > lo + wc) & (aa < 1 - wc) & (bb > lo + wc) & (bb < 1 - wc)
        out[name] = (float(D[inner].mean()), float(D2[inner].mean()))
    return out


def sprinkling_task(eps, seed, levels, with_balance):
    nl, nc = n_loc_of(eps), n_cut_of(eps)
    f = fker(eps, nc + 2)
    Tmax = max(levels)
    n_exp = Tmax ** 2 * nl
    u, v = sprinkle(n_exp, seed)
    order = sorted(levels, reverse=True)
    sides = [L / Tmax for L in order]
    lev = np.zeros(u.size, np.int64)
    for k, s in enumerate(sides):
        lev[(u >= 1 - s) & (v >= 1 - s)] = k
    vr, iv = ranks(u, v)
    H = abundances(vr, iv, nc, lev, len(sides))
    rows = {}
    for k, (L, s) in enumerate(zip(order, sides)):
        m = (u >= 1 - s) & (v >= 1 - s)
        N = int(m.sum())
        vrs, _ = ranks(u[m], v[m])
        R, hgt = rel_count_height(vrs)
        chi = R / (N * (N - 1) / 2)
        S = 4 * eps * (N - 2 * eps * float(np.dot(H[k], f[:H[k].size])))
        row = {"N": N, "S": S, "chi": chi, "d_mm": d_mm(chi), "h": int(hgt - 1)}
        if with_balance and L == min(levels):
            row["bal"] = _balance_inner(u[m], v[m], 1 - s, 1.0 / Tmax, f, nc, eps)
        rows[L] = row
    return {"rows": rows, "n_exp": n_exp, "sides": dict(zip(order, sides))}


def brute_task(eps, n_exp, seed, n_swaps):
    """I1: S and Delta S of adjacent transpositions against brute force (matrix product of the relation)."""
    nc = n_cut_of(eps)
    f = fker(eps, max(nc, 4000) + 2)

    def brute(u, v):
        rel = (u[:, None] < u[None, :]) & (v[:, None] < v[None, :])
        Rm = rel.astype(np.float64)
        n = np.rint((Rm @ Rm)[rel]).astype(np.int64)
        return 4 * eps * (u.size - 2 * eps * f[n].sum())

    u, v = sprinkle(n_exp, seed)
    vr, iv = ranks(u, v)
    H = abundances(vr, iv, nc, np.zeros(u.size, np.int64), 1)[0]
    S1 = 4 * eps * (u.size - 2 * eps * float(np.dot(H, f[:H.size])))
    S0 = brute(u, v)
    vrR, ivR = reversed_order(vr)
    dS = balance(vr, iv, vrR, ivR, f, nc, eps)
    rng = np.random.default_rng(seed + 1)
    ps = rng.choice(u.size - 1, min(n_swaps, u.size - 1), replace=False)
    worst = 0.0
    for p in ps:
        uu = u.copy()
        uu[p], uu[p + 1] = u[p + 1], u[p]
        worst = max(worst, abs((brute(uu, v) - S0) - dS[p]))
    return {"eps": eps, "N": int(u.size), "n_cut": nc, "dS_tot": abs(S1 - S0), "dS_swap": worst}


def _dispatch(kind, args):
    if kind == "spr":
        return sprinkling_task(*args)
    return brute_task(*args)


# ============================================================================================== run

def run(p, ctx):
    t0 = time.time()
    eps_list = [float(e) for e in p["eps"]]
    levels = [int(x) for x in p["levels"]]
    nmap = {float(k): int(v) for k, v in p["n"].items()}
    nb_ = int(p["n_balance"])
    jobs = []
    for eps in eps_list:
        for k in range(nmap[eps]):
            jobs.append(("spr", (eps, derive_seed(ctx.seed, "F6.9a", eps, k), levels, k < nb_)))
    brute_cases = [(0.12, 150, 1, 149), (0.24, 220, 2, 219), (0.06, 260, 4, 259), (0.5, 2500, 5, 12), (0.24, 1800, 6, 12)]
    for c in brute_cases:
        jobs.append(("brute", (c[0], c[1], derive_seed(ctx.seed, "F6.9a-brute", *c), c[3])))
    jobs.sort(key=lambda j: -(j[1][0] < 0.1))            # the slow eps first
    res = parallel_map(_dispatch, jobs, ctx.workers, ctx.progress,
                       lambda i, r_, d, n: f"[F6.9a] task {d}/{n} ({time.time() - t0:.0f}s)")
    spr = {eps: [r for j, r in zip(jobs, res) if j[0] == "spr" and j[1][0] == eps] for eps in eps_list}
    brute = [r for j, r in zip(jobs, res) if j[0] == "brute"]

    tol, nse = float(p["tol"]), float(p["n_se"])
    lines = [f"2D causal diamonds, nested T/w = {levels} (common future tip); S = 4 eps (N - 2 eps sum N_n f); per "
             f"sprinkling c = slope of S against N; c/(4 eps) per element; mean +- SE over sprinklings"]
    verdicts, data = {}, {"eps": {}}
    inst = {}
    # I1
    w1 = max(max(b["dS_tot"], b["dS_swap"]) for b in brute)
    inst["I1"] = w1 <= float(p["brute_tol"])
    lines.append("  I1 brute force: " + "; ".join(f"eps {b['eps']} N {b['N']} (n_cut {b['n_cut']}) S {b['dS_tot']:.1e} "
                                                  f"Delta S {b['dS_swap']:.1e}" for b in brute)
                 + f" -> max {w1:.1e} ({'OK' if inst['I1'] else 'FAILS'})")
    i2, i3, i4 = True, True, True
    for eps in eps_list:
        R = spr[eps]
        n = len(R)
        sides = R[0]["sides"]
        n_exp = R[0]["n_exp"]
        c = []
        for r in R:
            Ns = np.array([r["rows"][L]["N"] for L in levels], float)
            Ss = np.array([r["rows"][L]["S"] for L in levels])
            c.append(np.polyfit(Ns, Ss, 1)[0] / (4 * eps))
        c = np.array(c)
        m, se = float(c.mean()), float(c.std(ddof=1) / math.sqrt(n)) if n > 1 else float("nan")
        if abs(m) + nse * se <= tol:
            v = "PASSED"
        elif abs(m) - nse * se > tol:
            v = "FAILED"
        else:
            v = "UNDECIDED"
        verdicts[eps] = v
        lines.append(f"  eps {eps:g} (n_loc {n_loc_of(eps)}, {n} sprinklings): c/(4 eps) = {m:+.3e} +- {se:.3e}; "
                     f"|x| + {nse:g} SE = {abs(m) + nse * se:.2e} against {tol:g} -> {v}")
        ed = {"c4eps": m, "se": se, "verdict": v, "levels": {}}
        for L in levels:
            S = np.array([r["rows"][L]["S"] for r in R])
            N = np.array([r["rows"][L]["N"] for r in R], float)
            NE = n_exp * sides[L] ** 2
            ok2 = abs(N.mean() - NE) <= 3 * math.sqrt(NE / n)
            i2 &= ok2
            outs = {}
            for q in ("chi", "d_mm", "h"):
                lo_, hi_ = BANDS[eps][L][q]
                x = np.array([r["rows"][L][q] for r in R], float)
                outs[q] = float(np.mean((x < lo_) | (x > hi_)))
            ok3 = max(outs.values()) <= float(p["band_frac"])
            i3 &= ok3
            seS = S.std(ddof=1) / math.sqrt(n) if n > 1 else float("nan")
            lines.append(f"    T/w {L:>2}: N {N.mean():9.1f} (expected {NE:9.1f}; I2 {'OK' if ok2 else 'FAILS'})  "
                         f"<S> {S.mean():+8.3f} +- {seS:.3f}  outside bands chi/dMM/h {outs['chi']:.3f}/"
                         f"{outs['d_mm']:.3f}/{outs['h']:.3f} (I3 {'OK' if ok3 else 'FAILS'})")
            ed["levels"][L] = {"N": float(N.mean()), "S": float(S.mean()), "seS": float(seS), "out": outs}
        for a, b in zip(levels[:-1], levels[1:]):
            d = np.array([r["rows"][b]["S"] - r["rows"][a]["S"] for r in R])
            lines.append(f"    S({b}) - S({a}) = {d.mean():+.3f} +- {d.std(ddof=1) / math.sqrt(n) if n > 1 else float('nan'):.3f}")
        B = [r["rows"][min(levels)]["bal"] for r in R if "bal" in r["rows"][min(levels)]]
        for k in ("U", "V"):
            D = np.array([b[k][0] for b in B])
            D2 = np.array([b[k][1] for b in B])
            seD = D.std(ddof=1) / math.sqrt(D.size) if D.size > 1 else float("nan")
            ok4 = abs(D.mean()) + nse * seD <= float(p["tol_balance"]) if D.size > 1 else False
            i4 &= ok4
            lines.append(f"    I4 T/w {min(levels)} {k} ({D.size} sprinklings): <D S> {D.mean():+.2e} +- {seD:.2e} "
                         f"({'OK' if ok4 else 'FAILS'}); reported <D2 S> {D2.mean():+.2e}")
        data["eps"][str(eps)] = ed
    inst.update(I2=i2, I3=i3, I4=i4)
    lines.append("  instrument: " + "; ".join(f"{k} {'OK' if v else 'FAILS'}" for k, v in inst.items()))
    if not all(inst.values()):
        status = "INSTRUMENT_FAILED"
    elif all(v == "PASSED" for v in verdicts.values()):
        status = "PASSED"
    elif any(v == "FAILED" for v in verdicts.values()):
        status = "FAILED"
    else:
        status = "UNDECIDED"
    lines.append(f"cost: {time.time() - t0:.0f} s")
    data.update(instrument=inst, brute=brute)
    return Result(status, lines, data)
