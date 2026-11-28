"""
Q.1c -- the vacuum of one 1+1 sprinkled order: Johnston's Feynman propagator against the continuum SJ state
=========================================================================================================

Why this test exists (28 Nov 2026; card docs/quantum/Q1c_card_draft.md draft 1, NOT YET FROZEN)
---------------------------------------------------------------------------------------------------------------
Q.1 showed that one sprinkled order carries the propagation of a free scalar (the retarded propagator). Q.1c asks
whether one order also carries its vacuum: the imaginary part of Johnston's Feynman propagator (arXiv:0909.0944),
built from the order alone, against the continuum Sorkin-Johnston state of the same diamond. The author's route
"C, then A" (28 Nov): the complex propagator is calibrated before the interference test.

Construction (only from the matrices of the order; card "Construction")
-----------------------------------------------------------------------
K_R = aC(1 - abC)^(-1), a = 1/2, b = -m^2/rho; K_A = K_R^T; Delta = K_R - K_A; Q = the positive part of i Delta over
the full spectrum (no eigenvalue dropped); K_F = K_R + iQ; Im K_F = Re Q = |Delta|/2 with |Delta| = (Delta^T Delta)^(1/2),
computed as V S V^T from the singular value decomposition of Delta (the identity is instrument I1).

Reference (card "Reference")
----------------------------
The continuum SJ state of the light-cone square, computed in this test independently of any causal set: the exact
continuum Pauli-Jordan kernel ((1/2) J0(m tau) inside the light cone, half the jump on it) on the nested product
trapezoid grids n = 24, 48, 96, its positive part, extrapolated in n with the order estimated from the three grids;
uncertainty = the last extrapolation step. The Minkowski vacuum (-Y0/4 timelike, K0/(2 pi) spacelike) is reported,
not judged.

Model (card "Model")
--------------------
Light-cone square u, v in [0, 1]; m in {5, 10, 20}; N_exp in {1000, 2000, 4000, 8000} (the only scaled quantity);
12 probe pairs at u, v = 1/2 + k/24 (4 timelike du = dv, 4 timelike du = 2 dv, 4 spacelike), added as fixed elements,
in 8 groups whose points share no u and no v (general position) and lie outside the causal intervals of the timelike
pairs of their group, one group per sprinkling; R = 48 realisations per
(m, N_exp, group); seeds derive_seed(base, "Q.1c", m, N_exp, group, r).

Metrics and criterion (card "Metrics", "Criterion"; frozen)
----------------------------------------------------------
d_rj = (Re Q_r(p_j, q_j) - W_SJ(p_j, q_j)) / (1/4). Bias b_j = mean_r d_rj, SE_j^2 = var_r d_rj / R + unc_j^2;
chi^2 = sum_j (b_j / SE_j)^2 over the 12 pairs at N_exp = 8000. Fluctuation S = median over pairs of sd_r d_rj; slope s
of ln S against ln N_exp over the four densities, bootstrap 95 % [lo, hi] (realisations resampled within each pair and
density, B = 2000). Per m: bias consistent if p >= 0.01, deviant if p < 0.001, unclear otherwise; fluctuation falls if
hi < 0, grows if lo > 0, unclear otherwise. FAILED if any bias deviant or any fluctuation grows; PASSED if all biases
consistent and all fluctuations fall; else UNDECIDED (inconclusive). INSTRUMENT_FAILED if I1-I4 fails:
I1 |Delta|/2 against the explicit complex eigendecomposition of i Delta on 10 orders of N = 50-300, max |diff| <= 1e-10,
   and Im Q = Delta/2 <= 1e-12;
I2 at every timelike probe pair (every m and N_exp) the mean of K_R(p, q) equals (1/2) J0(m tau) within 3.5 SE;
I3 no two inserted points share u or v;
I4 the same construction with b = +m^2/rho at m = 10, N_exp = 2000 (5 groups x R) differs from W_SJ: chi^2 p < 1e-6.

Reading declared now (card)
---------------------------
PASSED: one sprinkled 1+1 order carries the vacuum of the same region (the SJ state built from the order alone); no
systematic difference from the continuum SJ state remains at N_exp = 8000 within the statistics and the reference, and
the error of one order falls with the density. Nothing about 3+1, interference or the Born rule.
FAILED: a systematic difference remains, or one order does not settle.

Crash resilience: every task is cached in results/Q_quantum/Q.1c/cache, keyed by its inputs and the hash of this file,
so an interrupted run continues where it stopped (quick runs are not cached). Memory: about 3 GB per worker at
N_exp = 8000, below 0.5 GB at N_exp <= 4000.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import j0, k0, y0
from scipy.stats import chi2 as chi2_dist

from core import provenance
from core.order.sprinkle import derive_seed
from core.stats.parallel import parallel_map
from core.testing import Result

TEST = {
    "id": "Q.1c",
    "phase": "Q",
    "title": "The vacuum of one 1+1 sprinkled order: Johnston's Feynman propagator against the continuum SJ state",
    "theory_part": "the quantum section (decision 49): the ground state of a free field from the order alone (the "
                   "Sorkin-Johnston state); the second quantum instrument, before interference (route C then A)",
    "question": "Does one sprinkled 1+1 order carry the vacuum: is Im K_F = Re Q of Johnston's Feynman propagator "
                "(built from the order alone, full spectrum) equal to the continuum SJ state of the same diamond, with "
                "no systematic difference at the highest density and a falling error of one order?",
    "method": "light-cone square, m = 5, 10, 20, N_exp = 1000, 2000, 4000, 8000; 12 probe pairs in general position (8 "
              "groups, one per sprinkling), R = 48 per (m, N_exp, group); Re Q = |Delta|/2; reference: the continuum "
              "SJ state by trapezoid quadrature of the exact kernel on n = 24, 48, 96, extrapolated",
    "prediction": "no systematic difference at N_exp = 8000 (development: chi^2/12 = 0.94 - 1.47 at 4000); a falling "
                  "error of one order (development: slopes about -0.33 ... -0.48)",
    "criterion": "per m: bias consistent if chi^2 (12 pairs, N_exp 8000) p >= 0.01, deviant if p < 0.001; fluctuation "
                 "falls if the bootstrap 95 % interval of the slope of ln S against ln N_exp has hi < 0, grows if "
                 "lo > 0. FAILED if any deviant or growing; PASSED if all consistent and falling; else UNDECIDED "
                 "(inconclusive). INSTRUMENT_FAILED if I1 (identity 1e-10, Im Q 1e-12), I2 (K_R mean at timelike "
                 "pairs within 3.5 SE), I3 (general position) or I4 (non-negative null must differ, p < 1e-6) fails",
    "null_model": "the same construction with non-negative weights b = +m^2/rho (I4): its 'vacuum' must differ",
    "level": "calibration",
    "meters": ["Q.1"],
    "depends_on": ["Q.1"],
    "params": {"ms": [5.0, 10.0, 20.0], "n_exp": [1000, 2000, 4000, 8000], "reps": 48, "chunk": {"1000": 12,
               "2000": 6, "4000": 2, "8000": 1}, "grids": [24, 48, 96], "boot": 2000, "p_ok": 0.01, "p_bad": 0.001,
               "n_se": 3.5, "null_m": 10.0, "null_n": 2000, "null_p": 1e-6, "i1_orders": 10, "i1_tol": 1e-10,
               "i1_tol_im": 1e-12, "cache": True},
    "parts": {"main": {}},
    "quick": {"n_exp": [400, 800], "reps": 6, "chunk": {"400": 3, "800": 3}, "grids": [24, 48], "boot": 200,
              "null_n": 400, "i1_orders": 3, "cache": False},
    "estimate": "about 10 h on 10 workers (N_exp = 8000: about 5 min per realisation, about 3 GB per worker); "
                "the reference about 10 min; an interrupted run continues from its cache",
}

A = 0.5
BASE = 24
_FILE = Path(__file__).resolve()
_CACHE = _FILE.parents[2] / "results" / "Q_quantum" / "Q.1c" / "cache"

# the 12 probe pairs (u, v offsets in units of 1/BASE from the centre): (kind, k) -> (p, q)
def _pair(kind, k):
    if kind == "time":
        return (-k, -k), (k, k)
    if kind == "space":
        return (-k, k), (k, -k)
    return (-2 * k, -k), (0, 0)                     # "skew": timelike, du = 2 dv


PAIRS = [(kd, k) for k in (1, 2, 3, 4) for kd in ("time", "space", "skew")]
# 8 groups: in each, all u and all v of the inserted points are distinct (general position), and no inserted point
# lies inside or on the causal interval of a timelike pair of its group (so that E K_R(p, q) = (1/2) J0(m tau), I2)
GROUPS = [[("time", 1), ("space", 2)], [("time", 2), ("space", 3)], [("time", 3), ("space", 4)], [("time", 4)],
          [("skew", 1)], [("skew", 2), ("space", 1)], [("skew", 3)], [("skew", 4)]]
GROUP_PAIRS = GROUPS


# =============================================================================
# cache
# =============================================================================

def _cached(tag, args, fn, use_cache):
    key = hashlib.sha256(json.dumps({"tag": tag, "args": args, "file": provenance.sha256_file(_FILE)},
                                    sort_keys=True).encode()).hexdigest()[:32]
    path = _CACHE / f"{key}.json"
    if use_cache and path.exists():
        try:
            return json.loads(path.read_text())
        except (OSError, ValueError):
            pass
    out = fn(*args)
    if use_cache:
        try:
            _CACHE.mkdir(parents=True, exist_ok=True)
            tmp = path.with_suffix(".tmp")
            tmp.write_text(json.dumps(out))
            tmp.replace(path)
        except OSError:
            pass
    return out


# =============================================================================
# the continuum SJ reference
# =============================================================================

def _kernel(du, dv, m):
    g = 0.5 * j0(m * np.sqrt(np.abs(du * dv)))
    fwd = (du > 0) & (dv > 0); bwd = (du < 0) & (dv < 0)
    nf = ((du == 0) & (dv > 0)) | ((dv == 0) & (du > 0))
    nb = ((du == 0) & (dv < 0)) | ((dv == 0) & (du < 0))
    return np.where(fwd, g, 0.0) - np.where(bwd, g, 0.0) + 0.5 * (np.where(nf, g, 0.0) - np.where(nb, g, 0.0))


def _sj_grid(m, n):
    x = np.arange(n + 1) / n
    w1 = np.full(n + 1, 1.0 / n); w1[0] = w1[-1] = 0.5 / n
    U, V = np.meshgrid(x, x, indexing="ij"); U = U.ravel(); V = V.ravel()
    sw = np.sqrt(np.outer(w1, w1).ravel() / 2.0)
    B = _kernel(U[None, :] - U[:, None], V[None, :] - V[:, None], m)
    B *= sw[:, None]; B *= sw[None, :]
    lam, E = np.linalg.eigh(B.T @ B)
    del B
    absB = (E * np.sqrt(np.clip(lam, 0.0, None))) @ E.T
    del E
    c = n // 2; f = n // BASE
    out = []
    for kd, k in PAIRS:
        p, q = _pair(kd, k)
        ip = (c + p[0] * f) * (n + 1) + (c + p[1] * f); iq = (c + q[0] * f) * (n + 1) + (c + q[1] * f)
        out.append(float(0.5 * absB[ip, iq] / (sw[ip] * sw[iq])))
    return out


def reference_task(m, grids):
    vals = [_sj_grid(m, n) for n in grids]
    res = []
    for j, (kd, k) in enumerate(PAIRS):
        v = [vv[j] for vv in vals]
        if len(v) >= 3:
            d1, d2 = v[-2] - v[-3], v[-1] - v[-2]
            r = d2 / d1 if d1 != 0 and 0 < d2 / d1 < 1 else 0.5
        else:
            d2, r = v[-1] - v[-2], 0.5
        lim = v[-1] + d2 * r / (1 - r)
        p, q = _pair(kd, k)
        du = (q[0] - p[0]) / BASE; dv = (q[1] - p[1]) / BASE; s = math.sqrt(abs(du * dv))
        mink = -0.25 * float(y0(m * s)) if kd != "space" else float(k0(m * s)) / (2 * math.pi)
        res.append({"W": lim, "unc": abs(lim - v[-1]), "grid_values": v, "ms": m * s, "mink": mink})
    return res


# =============================================================================
# one sprinkling with one probe group
# =============================================================================

def _realise(m, n_exp, seed, pairs, sign):
    pts = sorted({p for kd, k in pairs for p in _pair(kd, k)})
    us = {p[0] for p in pts}; vs = {p[1] for p in pts}
    general = len(us) == len(pts) and len(vs) == len(pts)            # I3
    rng = np.random.default_rng(seed)
    n = rng.poisson(n_exp)
    u = np.concatenate([rng.random(n), [0.5 + p[0] / BASE for p in pts]])
    v = np.concatenate([rng.random(n), [0.5 + p[1] / BASE for p in pts]])
    o = np.argsort(u, kind="stable"); u = u[o]; v = v[o]
    inv = np.empty_like(o); inv[o] = np.arange(o.size)
    pos = {p: int(inv[n + i]) for i, p in enumerate(pts)}
    N = u.size
    b = sign * m * m / (2.0 * n_exp)
    C = ((u[:, None] < u[None, :]) & (v[:, None] < v[None, :])).astype(float)
    M = np.eye(N) - A * b * C
    KR = A * np.linalg.solve(M.T, C.T).T                             # K_R = a C (1 - abC)^(-1)
    del M, C
    reQ, kr = [], []
    D = KR - KR.T
    for kd, k in pairs:
        p, q = _pair(kd, k)
        kr.append(float(KR[pos[p], pos[q]]))
    del KR
    _, sv, Vt = np.linalg.svd(D)                                     # |Delta| = V S V^T (no squaring)
    del D, _
    absD = (Vt.T * sv) @ Vt
    del Vt
    for kd, k in pairs:
        p, q = _pair(kd, k)
        reQ.append(float(0.5 * absD[pos[p], pos[q]]))
    return {"reQ": reQ, "kr": kr, "n": int(n), "general": bool(general)}


def chunk_task(m, n_exp, g, r0, r1, base, sign):
    tag = "Q.1c" if sign < 0 else "Q.1c-null"
    pairs = GROUP_PAIRS[g]
    rows = [_realise(m, n_exp, derive_seed(base, tag, m, n_exp, g, r), pairs, sign) for r in range(r0, r1)]
    return {"reQ": [r["reQ"] for r in rows], "kr": [r["kr"] for r in rows], "general": all(r["general"] for r in rows)}


def i1_task(base, n_orders):
    """|Delta|/2 against the explicit positive part of i Delta; Im Q against Delta/2."""
    worst_re = worst_im = 0.0
    for t in range(n_orders):
        rng = np.random.default_rng(derive_seed(base, "Q.1c-I1", t))
        n = int(rng.integers(50, 301)); m = float(rng.uniform(2.0, 20.0)); rho = 2.0 * n
        u = rng.random(n); v = rng.random(n)
        C = ((u[:, None] < u[None, :]) & (v[:, None] < v[None, :])).astype(float)
        KR = A * C @ np.linalg.inv(np.eye(n) - A * (-m * m / rho) * C)
        D = KR - KR.T
        lam, U = np.linalg.eigh(1j * D)
        P = lam > 0
        Q = (U[:, P] * lam[P]) @ U[:, P].conj().T
        _, sv, Vt = np.linalg.svd(D)
        absD = (Vt.T * sv) @ Vt
        worst_re = max(worst_re, float(np.max(np.abs(Q.real - 0.5 * absD))))
        worst_im = max(worst_im, float(np.max(np.abs(Q.imag - 0.5 * D))))
    return {"re": worst_re, "im": worst_im}


def _dispatch(kind, args, use_cache):
    if kind == "ref":
        return _cached("ref", list(args), reference_task, use_cache)
    if kind == "chunk":
        return _cached("chunk", list(args), chunk_task, use_cache)
    return i1_task(*args)


# =============================================================================
# run
# =============================================================================

def run(p, ctx):
    ms, ns, R = p["ms"], p["n_exp"], int(p["reps"])
    use_cache = bool(p["cache"]) and not ctx.quick
    jobs = [("i1", (ctx.seed, int(p["i1_orders"])), False)]
    for m in ms:
        jobs.append(("ref", (m, list(p["grids"])), use_cache))
    for m in ms:
        for n in ns:
            ch = int(p["chunk"][str(n)])
            for g in range(len(GROUPS)):
                for r0 in range(0, R, ch):
                    jobs.append(("chunk", (m, n, g, r0, min(R, r0 + ch), ctx.seed, -1.0), use_cache))
    for g in range(len(GROUPS)):
        step = int(p["chunk"][str(p["null_n"])]) if str(p["null_n"]) in p["chunk"] else 6
        for r0 in range(0, R, step):
            jobs.append(("chunk", (p["null_m"], p["null_n"], g, r0, min(R, r0 + step), ctx.seed, +1.0), use_cache))
    res = parallel_map(_dispatch, jobs, ctx.workers, ctx.progress, lambda i, r, d, n: f"  {d}/{n} tasks")
    i1 = res[0]
    refs = {m: res[1 + i] for i, m in enumerate(ms)}
    k = 1 + len(ms)
    cells = {}
    for m in ms:
        for n in ns:
            ch = int(p["chunk"][str(n)])
            for g in range(len(GROUPS)):
                for r0 in range(0, R, ch):
                    cells.setdefault((m, n, g), []).append(res[k]); k += 1
    null = {}
    for g in range(len(GROUPS)):
        step = int(p["chunk"][str(p["null_n"])]) if str(p["null_n"]) in p["chunk"] else 6
        for r0 in range(0, R, step):
            null.setdefault(g, []).append(res[k]); k += 1

    pidx = {pk: j for j, pk in enumerate(PAIRS)}

    def gather(cellmap, m, n=None):
        """d[j] (array over realisations) and K_R[j] per pair index j."""
        d, kr, gen = {}, {}, True
        for g, plist in enumerate(GROUP_PAIRS):
            chunks = cellmap[(m, n, g)] if n is not None else cellmap[g]
            rq = [row for c in chunks for row in c["reQ"]]; kk = [row for c in chunks for row in c["kr"]]
            gen &= all(c["general"] for c in chunks)
            for i, pk in enumerate(plist):
                j = pidx[pk]
                d[j] = (np.array([row[i] for row in rq]) - refs[m][j]["W"]) / 0.25
                kr[j] = np.array([row[i] for row in kk])
        return d, kr, gen

    L = ["Q.1c: the vacuum of one 1+1 sprinkled order (card docs/quantum/Q1c_card_draft.md)",
         f"  light-cone square; 12 probe pairs in 8 groups (general position); R = {R} per (m, N_exp, group); "
         f"d = (Re Q - W_SJ)/(1/4)"]
    data = {"i1": i1, "ref": {str(m): refs[m] for m in ms}, "cells": [], "per_m": []}
    i2_ok = True; i3_ok = True; i2_worst = 0.0
    verdicts = []
    xN = np.log(np.array(ns, float))
    for m in ms:
        L.append(f"  m {m:g}: reference W_SJ (grids {p['grids']}), max uncertainty "
                 f"{max(r['unc'] for r in refs[m]) / 0.25:.4f}; offset against Minkowski (reported) "
                 + ", ".join(f"{(r['W'] - r['mink']) / 0.25:+.4f}" for r in refs[m]))
        Ds = {}
        for n in ns:
            d, kr, gen = gather(cells, m, n)
            i3_ok &= gen
            Ds[n] = d
            b = np.array([d[j].mean() for j in range(len(PAIRS))])
            sd = np.array([d[j].std(ddof=1) for j in range(len(PAIRS))])
            unc = np.array([refs[m][j]["unc"] / 0.25 for j in range(len(PAIRS))])
            se = np.sqrt(sd ** 2 / R + unc ** 2)
            chi2 = float(np.sum((b / se) ** 2)); pval = float(chi2_dist.sf(chi2, len(PAIRS)))
            for j, (kd, kk) in enumerate(PAIRS):
                if kd == "space":
                    continue
                pp, qq = _pair(kd, kk)
                tau = math.sqrt(((qq[0] - pp[0]) / BASE) * ((qq[1] - pp[1]) / BASE))
                z = (kr[j].mean() - 0.5 * float(j0(m * tau))) / (kr[j].std(ddof=1) / math.sqrt(kr[j].size))
                i2_worst = max(i2_worst, abs(z)); i2_ok &= abs(z) <= p["n_se"]
            data["cells"].append({"m": m, "n_exp": n, "bias": b.tolist(), "sd": sd.tolist(), "se": se.tolist(),
                                  "chi2": chi2, "p": pval, "d": {str(j): d[j].tolist() for j in d}})
            L.append(f"    N_exp {n:5d} (m^2/rho {m*m/(2*n):.4f}): rms bias {np.sqrt(np.mean(b**2)):.4f}  chi2/12 "
                     f"{chi2/12:.2f} (p {pval:.3g})  median sd of one order {np.median(sd):.4f}")
        # fluctuation slope
        S0 = [math.log(np.median([Ds[n][j].std(ddof=1) for j in range(len(PAIRS))])) for n in ns]
        s0 = float(np.polyfit(xN, S0, 1)[0])
        rng = np.random.default_rng(derive_seed(ctx.seed, "Q.1c-bootstrap", m))
        bs = []
        for _ in range(int(p["boot"])):
            yb = []
            for n in ns:
                sds = []
                for j in range(len(PAIRS)):
                    x = Ds[n][j]
                    sds.append(x[rng.integers(0, x.size, x.size)].std(ddof=1))
                yb.append(math.log(np.median(sds)))
            bs.append(float(np.polyfit(xN, yb, 1)[0]))
        lo, md, hi = (float(q) for q in np.quantile(bs, [0.025, 0.5, 0.975]))
        last = data["cells"][-1]
        bias_band = "consistent" if last["p"] >= p["p_ok"] else "deviant" if last["p"] < p["p_bad"] else "unclear"
        fl_band = "falls" if hi < 0 else "grows" if lo > 0 else "unclear"
        verdicts.append((bias_band, fl_band))
        rms = [float(np.sqrt(np.mean(np.array(c["bias"]) ** 2))) for c in data["cells"] if c["m"] == m]
        data["per_m"].append({"m": m, "bias_band": bias_band, "p_last": last["p"], "fluct_slope": s0,
                              "boot_median": md, "lo": lo, "hi": hi, "fluct_band": fl_band,
                              "rms_bias_slope": float(np.polyfit(xN, np.log(rms), 1)[0])})
        L.append(f"    bias at N_exp {ns[-1]}: p {last['p']:.3g} -> {bias_band}; fluctuation slope {s0:+.3f} "
                 f"(bootstrap {md:+.3f}, 95 % [{lo:+.3f}, {hi:+.3f}]) -> {fl_band}; rms-bias slope (reported) "
                 f"{data['per_m'][-1]['rms_bias_slope']:+.3f}")
    # I4: the non-negative null
    dn, _, gen = gather(null, p["null_m"])
    i3_ok &= gen
    bn = np.array([dn[j].mean() for j in range(len(PAIRS))]); sdn = np.array([dn[j].std(ddof=1) for j in range(len(PAIRS))])
    uncn = np.array([refs[p["null_m"]][j]["unc"] / 0.25 for j in range(len(PAIRS))])
    chi2n = float(np.sum((bn / np.sqrt(sdn ** 2 / R + uncn ** 2)) ** 2)); pn = float(chi2_dist.sf(chi2n, len(PAIRS)))
    i4_ok = pn < p["null_p"]
    i1_ok = i1["re"] <= p["i1_tol"] and i1["im"] <= p["i1_tol_im"]
    data.update({"i2_max_z": i2_worst, "i3": bool(i3_ok), "i4": {"chi2": chi2n, "p": pn, "rms_bias": float(np.sqrt(np.mean(bn ** 2)))}})
    L.append(f"  I1 identity: |Re Q - |Delta|/2| {i1['re']:.1e} (tol {p['i1_tol']:g}), |Im Q - Delta/2| {i1['im']:.1e} "
             f"(tol {p['i1_tol_im']:g}) {'OK' if i1_ok else 'FAILED'}")
    L.append(f"  I2 K_R mean at the timelike pairs: max |z| {i2_worst:.2f} (tol {p['n_se']}) {'OK' if i2_ok else 'FAILED'}")
    L.append(f"  I3 general position: {'OK' if i3_ok else 'FAILED'}")
    L.append(f"  I4 non-negative null (b = +m^2/rho, m {p['null_m']:g}, N_exp {p['null_n']}): rms bias "
             f"{data['i4']['rms_bias']:.3f}, chi2/12 {chi2n/12:.1f}, p {pn:.2g} (must be < {p['null_p']:g}) "
             f"{'OK' if i4_ok else 'FAILED'}")
    if not (i1_ok and i2_ok and i3_ok and i4_ok):
        status, sub = "INSTRUMENT_FAILED", ""
    elif any(b == "deviant" or f == "grows" for b, f in verdicts):
        status, sub = "FAILED", ""
    elif all(b == "consistent" and f == "falls" for b, f in verdicts):
        status, sub = "PASSED", ""
    else:
        status, sub = "UNDECIDED", " (inconclusive)"
    data["sub_label"] = sub.strip(" ()") or None
    L.append(f"  verdict: {status}{sub}" + ("  [quick run: decides nothing]" if ctx.quick else ""))
    return Result(status, L, data)
