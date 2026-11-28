"""
Q.1 -- Johnston's sum-over-chains propagator on a single 1+1 sprinkling: self-averaging at the rate rho^(-1/2)
==============================================================================================================

Why this test exists (28 Nov 2026; card docs/quantum/Q1_card_draft.md draft 1, frozen by the author, decision 50)
---------------------------------------------------------------------------------------------------------------
The first test of the quantum section (decision 49). A particle on a causal set propagates from x to y as a sum over
chains with an amplitude a per relation and b per intermediate element (S. Johnston, arXiv:0806.3083, section 3):
K = a C (1 - a b C)^(-1), a = 1/2, b = -m^2/rho in 1+1. The expected number of chains with k intermediate elements is
(rho V)^k / (k!)^2, V = tau^2/2, so E K(x, y) = (1/2) J0(m tau) exactly at every density (Johnston eq. 3.23). The
question is about ONE order: does a single sprinkled order carry the propagator, with its error falling as
rho^(-1/2), or does the propagator appear only in the average? Q.1 tests the numerical realisation of a propagator on
a causal order, not quantum mechanics as a whole; interference against the non-negative flow is Q.2.

Model (card section 1)
----------------------
- Fixed region: the light-cone square u, v in [0, 1] (u = t + x, v = t - x; dt dx = du dv / 2); fixed elements
  x = (0, 0), y = (1, 1): tau_xy = 1, V = 1/2.
- Fixed mass m in {2, 5, 10}; reference (1/2) J0(m).
- The only varied quantity: rho = 2 N_exp, N_exp in {1000, 4000, 16000} (the expected number of interior elements).
- Each realisation: an independent Poisson sprinkling, N ~ Poisson(N_exp), uniform in (u, v), no conditioning on N;
  z < w iff u_z < u_w and v_z < v_w; K(x, z) for all z by the forward recursion
  k_j = a (C_xj + b sum_{i<j} k_i C_ij) in the order of u; K(x, y) the last entry.
- R = 2000 realisations per (m, N_exp); seeds derive_seed(base, "Q.1", m, N_exp, r).

Metrics (card section 2)
------------------------
e_r = K_r(x, y) - (1/2) J0(m); E_r = rms over the elements z of r with tau_xz >= 1/2 of K_r(x, z) - (1/2) J0(m tau_xz).
M1 = sd(e_r) (ddof 1); M2 = median(E_r). Slope s = least-squares slope of ln M against ln N_exp over the three
densities; bootstrap (R realisations resampled with replacement independently per density, B = 2000; seeds
derive_seed(base, "Q.1-bootstrap", m, metric)) gives the median and the 2.5 - 97.5 % interval [lo, hi].
Reported, not judged: median |e_r|, q90 |e_r| and their slopes; the distributions of e_r and E_r (in the data).

Criterion (frozen; card section 3)
----------------------------------
Per pair (m, metric): in band if -0.65 <= lo and hi <= -0.40; fail if lo > -0.25; deviation if hi < -0.65;
inconclusive otherwise. FAILED if any pair fails; PASSED if all six pairs are in band; otherwise UNDECIDED with the
sub-label "deviation" (any deviation, no fail) or "inconclusive". The bands are never widened after a run.

Instruments (card section 4; INSTRUMENT_FAILED if any fails)
------------------------------------------------------------
I1 the recursion against the dense solve a C (1 - a b C)^(-1) on 20 orders of N = 50 - 300: max |difference| <= 1e-12.
I2 the exact mean: at each of the nine points |mean e_r| <= 3.5 SE (SE = M1 / sqrt(R)); an implementation check.
I3 the null with non-negative weights, b = +m^2/rho, m = 5, N_exp = 4000, R = 2000: |mean - (1/2) I0(5)| <= 3.5 SE.

Reading declared now (card section 5)
-------------------------------------
PASSED: one sprinkled 1+1 order carries Johnston's propagator, its error falling as rho^(-1/2) at fixed m tau both at
a fixed pair and over its own elements; a calibration of the first instrument of the quantum section, no physical
conclusion of its own; nothing about 3+1 (Q.1b), interference (Q.2), the Born rule (Q.3) or m^2/rho -> 1.
FAILED: the propagator appears only in the average; Q.2 - Q.4 would need an ensemble or another instrument.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.special import i0, j0

from core.order.sprinkle import derive_seed
from core.stats.parallel import parallel_map
from core.testing import Result

TEST = {
    "id": "Q.1",
    "phase": "Q",
    "title": "Johnston's sum-over-chains propagator on a single 1+1 sprinkling: self-averaging at rho^(-1/2)",
    "theory_part": "the quantum section (decision 49): amplitudes over chains of the order instead of non-negative "
                   "flows; the first quantum instrument",
    "question": "Does a single sprinkled 1+1 order carry Johnston's propagator (1/2) J0(m tau), with its error "
                "falling at the rate rho^(-1/2), or does the propagator appear only in the average over orders?",
    "method": "fixed light-cone square, x = (0,0), y = (1,1), m = 2, 5, 10; Poisson sprinklings of N_exp = 1000, "
              "4000, 16000; K = aC(1 - abC)^(-1) by forward recursion, a = 1/2, b = -m^2/rho; 2000 realisations per "
              "point; M1 = sd of the error at (x, y), M2 = median over realisations of the rms error over the "
              "elements with tau >= 1/2; slopes against ln N_exp with bootstrap intervals",
    "prediction": "slopes -1/2 (Poisson self-averaging; development: -0.50 ... -0.53)",
    "criterion": "per (m, metric) in band if -0.65 <= lo and hi <= -0.40 (bootstrap 95 % of the slope), fail if "
                 "lo > -0.25, deviation if hi < -0.65, inconclusive otherwise; FAILED if any fail, PASSED if all six "
                 "in band, else UNDECIDED (deviation / inconclusive). INSTRUMENT_FAILED if I1 (recursion vs dense "
                 "1e-12), I2 (exact mean within 3.5 SE at all nine points) or I3 (non-negative null mean (1/2) I0(5) "
                 "within 3.5 SE) fails",
    "null_model": "the same sum with non-negative weights, b = +m^2/rho: mean (1/2) I0(m tau), no zeros (I3)",
    "level": "calibration",
    "meters": [],
    "depends_on": [],
    "params": {"ms": [2.0, 5.0, 10.0], "n_exp": [1000, 4000, 16000], "reps": 2000, "chunk": 100, "boot": 2000,
               "tau_min": 0.5, "band": [-0.65, -0.40], "fail_above": -0.25, "n_se": 3.5,
               "null_m": 5.0, "null_n": 4000, "i1_orders": 20, "i1_tol": 1e-12},
    "parts": {"main": {}},
    "quick": {"reps": 40, "chunk": 20, "boot": 200, "i1_orders": 4},
    "estimate": "about 10 min on 6 workers (dominated by N_exp = 16000: ~0.3 s per realisation)",
}

A = 0.5


def _k_row(u, v, b):
    """K(x, z) for the source x = (0, 0) below every element; u sorted increasingly (ties have probability 0)."""
    n = u.size
    k = np.empty(n)
    for j in range(n):
        rel = v[:j] < v[j]
        k[j] = A * (1.0 + b * float(k[:j][rel].sum()))
    return k


def _realise(m, n_exp, seed, tau_min, sign=-1.0):
    rng = np.random.default_rng(seed)
    n = rng.poisson(n_exp)
    u = rng.random(n); v = rng.random(n)
    o = np.argsort(u)
    u = np.append(u[o], 1.0); v = np.append(v[o], 1.0)          # y = (1, 1) is the last element
    b = sign * m * m / (2.0 * n_exp)
    k = _k_row(u, v, b)
    if sign < 0:
        tau = np.sqrt(u[:-1] * v[:-1])
        far = tau >= tau_min
        err = k[:-1][far] - 0.5 * j0(m * tau[far])
        E = float(np.sqrt(np.mean(err ** 2))) if far.any() else float("nan")
    else:
        E = float("nan")
    return float(k[-1]), E


def _chunk(m, n_exp, r0, r1, base, tau_min, sign):
    tag = "Q.1" if sign < 0 else "Q.1-null"
    out = [_realise(m, n_exp, derive_seed(base, tag, m, n_exp, r), tau_min, sign) for r in range(r0, r1)]
    return {"m": m, "n_exp": n_exp, "r0": r0, "sign": sign, "K": [o[0] for o in out], "E": [o[1] for o in out]}


def _i1(base, n_orders, tol):
    """The recursion against the dense solve on small orders (random m, N)."""
    worst = 0.0
    for t in range(n_orders):
        rng = np.random.default_rng(derive_seed(base, "Q.1-I1", t))
        n = int(rng.integers(50, 301)); m = float(rng.uniform(1.0, 10.0)); rho = 2.0 * n
        u = rng.random(n); v = rng.random(n)
        o = np.argsort(u); u = np.append(0.0, u[o]); v = np.append(0.0, v[o])   # x = (0, 0) first
        b = -m * m / rho
        C = ((u[:, None] < u[None, :]) & (v[:, None] < v[None, :])).astype(float)
        Kd = A * C @ np.linalg.inv(np.eye(n + 1) - A * b * C)
        kr = _k_row(u[1:], v[1:], b)
        worst = max(worst, float(np.max(np.abs(Kd[0, 1:] - kr))))
    return {"max_diff": worst, "ok": worst <= tol}


def _dispatch(kind, args):
    return _chunk(*args) if kind == "chunk" else _i1(*args)


def _slope(x, y):
    return float(np.polyfit(x, y, 1)[0])


def run(p, ctx):
    ms, ns, R, ch = p["ms"], p["n_exp"], int(p["reps"]), int(p["chunk"])
    jobs = [("i1", (ctx.seed, int(p["i1_orders"]), p["i1_tol"]))]
    for m in ms:
        for n in ns:
            for r0 in range(0, R, ch):
                jobs.append(("chunk", (m, n, r0, min(R, r0 + ch), ctx.seed, p["tau_min"], -1.0)))
    for r0 in range(0, R, ch):
        jobs.append(("chunk", (p["null_m"], p["null_n"], r0, min(R, r0 + ch), ctx.seed, p["tau_min"], +1.0)))
    res = parallel_map(_dispatch, jobs, ctx.workers, ctx.progress,
                       lambda i, r, d, n: f"  {d}/{n} tasks")
    i1 = res[0]
    pts = {}
    null = []
    for r in res[1:]:
        if r["sign"] < 0:
            pts.setdefault((r["m"], r["n_exp"]), []).append(r)
        else:
            null.append(r)
    L = ["Q.1: Johnston's sum-over-chains propagator on one 1+1 sprinkling (card docs/quantum/Q1_card_draft.md, "
         "decision 50)", f"  region u, v in [0, 1], x = (0, 0), y = (1, 1), tau = 1; R = {R} realisations per point"]
    data = {"points": [], "slopes": [], "i1": i1}
    i2_ok = True
    stats = {}
    for m in ms:
        ex = 0.5 * float(j0(m))
        for n in ns:
            rows = sorted(pts[(m, n)], key=lambda r: r["r0"])
            K = np.concatenate([r["K"] for r in rows]); E = np.concatenate([r["E"] for r in rows])
            e = K - ex
            M1 = float(np.std(e, ddof=1)); M2 = float(np.median(E)); se = M1 / math.sqrt(K.size)
            z = float(np.mean(e)) / se
            i2_ok &= abs(z) <= p["n_se"]
            stats[(m, n)] = (e, E)
            data["points"].append({"m": m, "n_exp": n, "exact": ex, "mean_e": float(np.mean(e)), "se": se, "z": z,
                                   "M1": M1, "M2": M2, "median_abs_e": float(np.median(np.abs(e))),
                                   "q90_abs_e": float(np.quantile(np.abs(e), 0.9)), "e": e.tolist(),
                                   "E": E.tolist()})
            L.append(f"  m {m:g} N_exp {n:6d} (m^2/rho {m*m/(2*n):.4f}): mean e {np.mean(e):+.6f} ({z:+.2f} SE)  "
                     f"M1 {M1:.6f} ({M1/abs(ex)*100:.2f} % of |K|)  M2 {M2:.6f}  median |e| "
                     f"{np.median(np.abs(e)):.6f}  q90 |e| {np.quantile(np.abs(e), 0.9):.6f}")
    x = np.log(np.array(ns, float))
    funcs = {"M1": lambda e, E: np.std(e, ddof=1), "M2": lambda e, E: np.median(E),
             "median |e|": lambda e, E: np.median(np.abs(e)), "q90 |e|": lambda e, E: np.quantile(np.abs(e), 0.9)}
    bands = []
    lo_b, hi_b = p["band"]
    for m in ms:
        for name, f in funcs.items():
            y = [math.log(f(*stats[(m, n)])) for n in ns]
            s0 = _slope(x, y)
            rng = np.random.default_rng(derive_seed(ctx.seed, "Q.1-bootstrap", m, name))
            bs = []
            for _ in range(int(p["boot"])):
                yb = []
                for n in ns:
                    e, E = stats[(m, n)]
                    i = rng.integers(0, e.size, e.size)
                    yb.append(math.log(f(e[i], E[i])))
                bs.append(_slope(x, yb))
            lo, md, hi = (float(q) for q in np.quantile(bs, [0.025, 0.5, 0.975]))
            judged = name in ("M1", "M2")
            if judged:
                band = ("in band" if (lo >= lo_b and hi <= hi_b) else "fail" if lo > p["fail_above"]
                        else "deviation" if hi < lo_b else "inconclusive")
                bands.append(band)
            data["slopes"].append({"m": m, "metric": name, "slope": s0, "boot_median": md, "lo": lo, "hi": hi,
                                   "judged": judged, "band": band if judged else None})
            L.append(f"  m {m:g} {name:10s}: slope {s0:+.4f}  bootstrap median {md:+.4f}  95 % [{lo:+.4f}, {hi:+.4f}]"
                     + (f"  -> {band}" if judged else "  (reported)"))
    # I3: the non-negative null
    Kn = np.concatenate([r["K"] for r in sorted(null, key=lambda r: r["r0"])])
    exn = 0.5 * float(i0(p["null_m"]))
    sen = float(np.std(Kn, ddof=1)) / math.sqrt(Kn.size)
    zn = (float(np.mean(Kn)) - exn) / sen
    i3_ok = abs(zn) <= p["n_se"]
    data["i3"] = {"mean": float(np.mean(Kn)), "exact": exn, "se": sen, "z": zn, "sd": float(np.std(Kn, ddof=1)),
                  "ok": bool(i3_ok)}
    L.append(f"  I1 recursion vs dense solve: max |diff| {i1['max_diff']:.2e} (tol {p['i1_tol']:g}) "
             f"{'OK' if i1['ok'] else 'FAILED'}")
    L.append(f"  I2 exact mean: max |z| {max(abs(q['z']) for q in data['points']):.2f} (tol {p['n_se']}) "
             f"{'OK' if i2_ok else 'FAILED'}")
    L.append(f"  I3 non-negative null (b = +m^2/rho, m {p['null_m']:g}, N_exp {p['null_n']}): mean {np.mean(Kn):.4f} "
             f"against (1/2) I0 = {exn:.4f} ({zn:+.2f} SE), sd {data['i3']['sd']:.4f} {'OK' if i3_ok else 'FAILED'}")
    if not (i1["ok"] and i2_ok and i3_ok):
        status, sub = "INSTRUMENT_FAILED", ""
    elif "fail" in bands:
        status, sub = "FAILED", ""
    elif all(b == "in band" for b in bands):
        status, sub = "PASSED", ""
    else:
        status, sub = "UNDECIDED", " (deviation)" if "deviation" in bands else " (inconclusive)"
    data["sub_label"] = sub.strip(" ()") or None
    L.append(f"  verdict: {status}{sub}" + ("  [quick run: decides nothing]" if ctx.quick else ""))
    return Result(status, L, data)
