"""
C0, revision b of the card (F3.16b): development rows, per-shell statistics and the tolerances (14 Nov; the
author's plan of 08:52 after F3.16 INSTRUMENT_FAILED; decides nothing until the card is frozen)
===============================================================================================================
Changes against card revision 1 (F3.16), each a design error found in F3.16 (docs/findings.md):
  1. eps per shell: each shell r < 2.5 on its own, |eps - (A1)| + 2 SE < its own tolerance (F3.16 took the largest
     of 10 shells with the SE of the chosen shell, against a tolerance derived from the SE of another shell);
  2. the instrument criterion against (A3) on the linear rows k 0..3 only (the nonlinear A 0.5 against the
     noise-free (A3) contains the renormalization, which is physics); the nonlinear rows against (A3) are reported;
  3. rho 60 with 12 sprinklings (F3.16: 4);
  4. tolerances from development rows with the card's own numbers of sprinklings (16 at rho 30, 12 at rho 60;
     F3.16's were derived from 4), on seeds never used for a verdict (19000+); the card runs on fresh seeds.
Model, kernels, readouts, the predictions (A1), (A3), B and the rule of 10 Nov unchanged. The rule runs through
diag/c0_fast.evolve_fast, bit-identical to diag/c0_card.evolve_card and about 2x faster.

    python -m diag.c0_card_b --workers 4          (development rows: rho 60 x 12, rho 30 x 16; then the tolerances)
"""
from __future__ import annotations

import argparse
import json
import math
import os
import time

import numpy as np

from core.stats.parallel import parallel_map, default_workers
from diag.c0_card import (KERNELS, MODES, NRB, RBIN, T0, X0, _ES, analyse, card_stats, rows, sd_ratio)
from diag.c0_ct import exact_field
from diag.c0_fast import evolve_fast
from diag.c0_links4 import sprinkle4
from diag.c0_minvar import NMAX_ALL
from diag.c0_apex import eps_pred

SEED0 = 19000
N_CARD = {30.0: 16, 60.0: 12}
MAIN = "MV8+ct"


def task(rho, seed, X=X0, T=T0, seed0=SEED0):
    """diag/c0_card.task with the bit-identical evolve_fast (diag/c0_fast.py)."""
    t, x, y, z = sprinkle4(rho, X, T, seed0 + seed)
    RS = rows(rho)
    W = np.array([r[1] for r in RS])
    ctd = np.array([r[2] for r in RS])
    phi0 = np.array([m[3] * np.cos(m[2] * x) * np.cos(m[2] * t) for _ in RS for m in MODES])
    lin = np.array([m[1] for m in MODES])
    c = 1.0 / (2.0 * math.pi * math.sqrt(6.0 * rho))
    tc = (24.0 * (NMAX_ALL + 6) / (math.pi * rho)) ** 0.25
    inner = (np.abs(x) < X - t - 0.3) & (np.abs(y) < X - t - 0.3) & (np.abs(z) < X - t - 0.3)
    sample = inner & (t >= 3.5) & (t < 5.5)
    a0 = time.time()
    phi, cnt, hist = evolve_fast(t, x, y, z, tc, NMAX_ALL, W, ctd, phi0, lin, c, 1.0, sample, RBIN, NRB)
    sec = time.time() - a0
    tb = np.arange(0.25, T - 0.4, 0.5)
    slabs = [inner & (np.abs(t - b) < 0.25) for b in tb]
    ex = {m[0]: exact_field(m, t, x, T) for m in MODES}
    tg = np.arange(0.0, T + 0.01, 0.005)
    out = {"rho": rho, "seed": seed, "N": int(t.size), "sec": sec, "t": tb.tolist(), "X": X, "T": T,
           "nslab": [int(m.sum()) for m in slabs], "nsample": int(sample.sum()), "k": {}}
    Mo = len(MODES)
    rr = (np.arange(NRB) + 0.5) * RBIN
    for kk, (name, w, ct) in enumerate(RS):
        S = c * (cnt @ w)
        E = np.interp(t, tg, _ES(w, rho, tg))
        K = {"ct": ct,
             "eps": (c * (hist @ w) / sample.sum() / (rr * RBIN) - 1.0).tolist(),
             "norm": [float(S[m].mean() / E[m].mean() - 1.0) for m in slabs],
             "varS": [float((S[m] - E[m]).var()) for m in slabs],
             "vardiag": [float(c * c * (cnt[m].mean(0) @ (w * w))) for m in slabs],
             "modes": {}}
        for mo, m in enumerate(MODES):
            ph = phi[kk * Mo + mo]
            dev = (ph - ex[m[0]]) / m[3]
            bas = np.cos(m[2] * x)
            dd = []
            for sl in slabs:
                amp = float(dev[sl] @ bas[sl] / (bas[sl] @ bas[sl]))
                eta = (dev[sl] - amp * bas[sl]) * m[3]
                row = [amp, float(eta.std() / m[3])]
                if not m[1] and m[2] == 0:
                    sd = eta.std()
                    z3 = float(np.mean((eta / sd) ** 3)) if sd > 0 else 0.0
                    z4 = float(np.mean((eta / sd) ** 4) - 3.0) if sd > 0 else 0.0
                    direct = float(np.mean(np.sin(ph[sl])) - math.exp(-0.5 * sd * sd) * np.mean(np.sin(ph[sl] - eta)))
                    row += [z3, z4, direct]
                dd.append(row)
            K["modes"][m[0]] = dd
        out["k"][name] = K
    return out


def cached(*a):
    import hashlib
    key = hashlib.sha256(json.dumps([str(q) for q in a] + ["card-b-v1"]).encode()).hexdigest()[:20]
    fn = os.path.join("diag", "out", "c0_card_b_cache", key + ".json")
    if os.path.exists(fn):
        return json.load(open(fn))
    out = task(*a)
    os.makedirs(os.path.dirname(fn), exist_ok=True)
    json.dump(out, open(fn + ".tmp", "w"))
    os.replace(fn + ".tmp", fn)
    return out


def shells(R, name):
    """Per shell r < 2.5: eps - (A1), mean and SE over sprinklings."""
    rho = R[0]["rho"]
    w = [q for q in rows(rho) if q[0] == name][0][1]
    ep = np.array(eps_pred(w, rho, RBIN, NRB))
    E = np.array([r["k"][name]["eps"] for r in R]) - ep
    return E.mean(0), E.std(0, ddof=1) / math.sqrt(len(R))


def ceil_to(v):
    if not math.isfinite(v):
        return float('nan')
    step = 0.001 if v < 0.02 else (0.005 if v < 1 else 0.05)
    return round(math.ceil(v / step - 1e-9) * step, 6)


def derive(res):
    """Tolerances of revision b by the rule of 10 Nov: ceil(1.25 x (|value| + 2 SE)), the development at the card's n."""
    rhos = (30.0, 60.0)
    A = {rho: card_stats([r for r in res if r["rho"] == rho]) for rho in rhos}
    nd = {rho: len([r for r in res if r["rho"] == rho]) for rho in rhos}
    tol, table = {}, []

    def put(key, label, rho, items):
        f = math.sqrt(nd[rho] / N_CARD[rho])
        best = max(items, key=lambda q: q[1] + 2 * q[2] * f)
        t = ceil_to(1.25 * (best[1] + 2 * best[2] * f))
        tol[key] = t
        table.append((label, f"{rho:g}", best[0], best[1], 2 * best[2] * f, t))

    for rho in rhos:
        g = f"{rho:g}"
        o = A[rho][MAIN]
        ms = o["modes"]
        mean, se = shells([r for r in res if r["rho"] == rho], MAIN)
        for i in range(NRB):
            put(f"eps_{g}_{i}", f"|eps - (A1)|, shell {i} (r {(i + 0.5) * RBIN:.3f})", rho, [(f"shell {i}", abs(mean[i]), se[i])])
        put(f"norm_{g}", "|normalisation - 1|", rho, [("t >= 3", abs(o["norm"][0]), o["norm"][1])])
        put(f"A3_{g}", "max |amp-dev - (A3)|, linear k 0-3", rho, [(k, *ms[k]["A3"]) for k in ("lin0", "lin1", "lin2", "lin3")])
        for a in ("0.5", "1", "2"):
            put(f"B_{a}_{g}", f"residual max |meas - B|, A {a}", rho, [(f"nl0_A{a}", *ms[f"nl0_A{a}"]["B"])])
        for what, label in (("direct", "|direct check|"), ("skew", "|skewness|, t >= 1.5"), ("kurt", "|excess kurtosis|, t >= 1.5")):
            put(f"{what}_{g}", label, rho, [(k, *ms[k][what]) for k in ("nl0_A0.5", "nl0_A1", "nl0_A2")])
        put(f"p_{g}", "growth exponent p", rho, [(m[0], *ms[m[0]]["p"]) for m in MODES])
    o = A[60.0][MAIN]["modes"]
    put("lin_cont_60", "linear (A): max |amp-dev|, k 0-3", 60.0, [(k, *o[k]["cont"]) for k in ("lin0", "lin1", "lin2", "lin3")])
    for a in ("0.5", "1", "2"):
        put(f"nl_cont_{a}_60", f"nonlinear (A) against the pendulum, A {a}", 60.0, [(f"nl0_A{a}", *o[f"nl0_A{a}"]["cont"])])
    put("nl1_cont_60", "nonlinear (A) against sine-Gordon, A 1 cos x", 60.0, [("nl1_A1", *o["nl1_A1"]["cont"])])
    signs = []
    a30, a60 = A[30.0][MAIN]["modes"], A[60.0][MAIN]["modes"]
    for k in ("lin0", "lin1", "lin2", "lin3"):
        signs.append((f"convergence {k}", a30[k]["cont"][0] - a60[k]["cont"][0], math.hypot(a30[k]["cont"][1], a60[k]["cont"][1])))
    for m in MODES:
        r, se = sd_ratio(A[30.0], A[60.0], MAIN, m[0])
        signs.append((f"noise falls {m[0]}", 1 - r, se))
    for k in ("nl0_A1", "nl0_A2"):
        signs.append((f"renormalization falls {k}", a30[k]["renorm"] - a60[k]["renorm"], math.hypot(a30[k]["B"][1], a60[k]["B"][1])))
    return tol, table, signs, A


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=0)
    ap.add_argument("--derive3", type=str, default="", help="development JSON: print the tolerances of the final form")
    a = ap.parse_args()
    if a.derive3:
        res = json.load(open(a.derive3))
        tol, table, z, F0s = derive3(res)
        print(f"z (one-sided, family-wise 5 % over {len(table)} criteria) {z:.3f}")
        for k, m0, se, t in table:
            print(f"  {k}: M_dev {m0:.4f}, SE_M (bootstrap) {se:.4f} -> tolerance {t:g}")
        st, inst, pc, fc, info = verdict3(res, tol, F0s=F0s)
        print("the development itself:", st)
        print("TOL " + json.dumps(tol))
        return
    jobs = [(60.0, s, X0, T0, SEED0) for s in range(N_CARD[60.0])] + [(30.0, s, X0, T0, SEED0) for s in range(N_CARD[30.0])]
    t0 = time.time()
    res = parallel_map(cached, jobs, default_workers(a.workers), print,
                       lambda i, r, d, n: f"[C0 card b dev] {d}/{n}: rho {r['rho']:g} seed {r['seed']} N {r['N']} "
                                          f"{r['sec']:.0f}s ({time.time() - t0:.0f}s)")
    fn = os.path.join("diag", "out", f"c0_card_b_dev_{time.strftime('%Y%m%dT%H%M%S', time.gmtime())}.json")
    json.dump(res, open(fn, "w"))
    print("saved", fn)
    for rho in (30.0, 60.0):
        analyse([r for r in res if r["rho"] == rho])
    tol, table, signs, _ = derive(res)
    print("######## tolerances of revision b (rule of 10 Nov; development at the card's n)")
    print("criterion | rho | deciding row | development value | 2 SE | tolerance")
    for row in table:
        print(f"{row[0]} | {row[1]} | {row[2]} | {row[3]:.4f} | {row[4]:.4f} | {row[5]:g}")
    print("sign criteria (development): difference +- SE")
    for label, d, se in signs:
        print(f"{label}: {d:+.4f} +- {se:.4f}")
    print("TOL " + json.dumps(tol))



def verdict(res, tol, nse=2.0):
    """The criteria of revision b (MV8+ct) on runs at rho 30 and 60; returns (status, inst, passc, failc, A)."""
    rhos = (30.0, 60.0)
    A = {rho: card_stats([r for r in res if r["rho"] == rho]) for rho in rhos}
    inst, passc, failc = {}, {}, {}

    def t(store_p, store_f, key, val, se):
        store_p[key] = bool(val + nse * se < tol[key])
        if store_f is not None:
            store_f[key] = bool(val - nse * se > tol[key])

    def s(key, d, se):
        passc[key] = bool(d - nse * se > 0)
        failc[key] = bool(d + nse * se < 0)

    for rho in rhos:
        g = f"{rho:g}"
        o = A[rho][MAIN]
        ms = o["modes"]
        mean, se = shells([r for r in res if r["rho"] == rho], MAIN)
        for i in range(NRB):
            t(inst, None, f"eps_{g}_{i}", abs(mean[i]), se[i])
        t(inst, None, f"norm_{g}", abs(o["norm"][0]), o["norm"][1])
        q = max((ms[k]["A3"] for k in ("lin0", "lin1", "lin2", "lin3")), key=lambda z: z[0] + nse * z[1])
        t(inst, None, f"A3_{g}", *q)
        for a in ("0.5", "1", "2"):
            t(passc, failc, f"B_{a}_{g}", *ms[f"nl0_A{a}"]["B"])
        for what in ("direct", "skew", "kurt"):
            q = max((ms[k][what] for k in ("nl0_A0.5", "nl0_A1", "nl0_A2")), key=lambda z: z[0] + nse * z[1])
            t(passc, failc, f"{what}_{g}", *q)
        q = max((ms[m[0]]["p"] for m in MODES), key=lambda z: z[0] + nse * z[1])
        t(passc, failc, f"p_{g}", *q)
    o = A[60.0][MAIN]["modes"]
    q = max((o[k]["cont"] for k in ("lin0", "lin1", "lin2", "lin3")), key=lambda z: z[0] + nse * z[1])
    t(passc, failc, "lin_cont_60", *q)
    for a in ("0.5", "1", "2"):
        t(passc, failc, f"nl_cont_{a}_60", *o[f"nl0_A{a}"]["cont"])
    t(passc, failc, "nl1_cont_60", *o["nl1_A1"]["cont"])
    a30, a60 = A[30.0][MAIN]["modes"], A[60.0][MAIN]["modes"]
    for k in ("lin0", "lin1", "lin2", "lin3"):
        s(f"conv_{k}", a30[k]["cont"][0] - a60[k]["cont"][0], math.hypot(a30[k]["cont"][1], a60[k]["cont"][1]))
    for m in MODES:
        r, se = sd_ratio(A[30.0], A[60.0], MAIN, m[0])
        s(f"noise_{m[0]}", 1 - r, se)
    for k in ("nl0_A1", "nl0_A2"):
        s(f"renorm_{k}", a30[k]["renorm"] - a60[k]["renorm"], math.hypot(a30[k]["B"][1], a60[k]["B"][1]))
    if not all(inst.values()):
        st = "INSTRUMENT_FAILED"
    elif any(failc.values()):
        st = "FAILED"
    elif all(passc.values()):
        st = "PASSED"
    else:
        st = "UNDECIDED"
    return st, inst, passc, failc, A


# ---- revision b, statistical form (proposed 14 Nov, 15:00, before freezing; the author decides) ----------------------
# The instrument criteria compare the run with predictions that are exact in expectation ((A1) per shell, the
# normalisation, (A3) of the linear rows: docs/derivations/c0_mean_exactness.md). They become a significance test:
# INSTRUMENT_FAILED if any item deviates by more than T_CRIT x SE_eff, SE_eff^2 = SE^2 + E_FLOOR^2 (E_FLOOR covers the
# evaluation of the predictions and the first slabs, where the scatter between sprinklings is ~1e-4), family-wise
# alpha 1 % over all items of both densities (Bonferroni, Student t with n - 1 degrees of freedom).
# The physics tolerances keep the rule of 10 Nov with a statistical margin in place of the 25 %:
# tol = |dev| + (2 + Z_PHYS) SE (rounded up), Z_PHYS one-sided for family-wise 5 % over the physics criteria.
E_FLOOR = 0.002
ALPHA_INST = 0.01
ALPHA_PHYS = 0.05


def inst_items(res, rho):
    """(label, deviation, SE) of every instrument item at one density: 10 shells, the normalisation, 4 x slabs of (A3)."""
    from diag.c0_apex import dyn_pred
    R = [r for r in res if r["rho"] == rho]
    n = len(R)
    tb = np.array(R[0]["t"])
    name, w, ct = rows(rho)[0]
    items = []
    m, se = shells(R, name)
    items += [(f"eps shell {i}", float(m[i]), float(se[i])) for i in range(NRB)]
    late = tb >= 3.0
    nm = np.array([np.mean(np.array(r["k"][name]["norm"])[late]) for r in R])
    items.append(("norm", float(nm.mean()), float(nm.std(ddof=1) / math.sqrt(n))))
    for k in range(4):
        P = np.array(dyn_pred(w, rho, k, tb, X=R[0]["X"], ct=ct))
        D = np.array([r["k"][name]["modes"][f"lin{k}"] for r in R])[:, :, 0] - P
        for j in range(tb.size):
            items.append((f"lin{k} t {tb[j]:.2f}", float(D[:, j].mean()), float(D[:, j].std(ddof=1) / math.sqrt(n))))
    return items, n


def t_crit(nitems_total, n):
    from scipy import stats
    return float(stats.t.ppf(1 - ALPHA_INST / nitems_total / 2, n - 1))


def derive_stat(res):
    """Physics tolerances with the statistical margin; the instrument needs no tolerance."""
    from scipy import stats
    tol, table, signs, A = derive(res)
    keys = [k for k in tol if not (k.startswith("eps_") or k.startswith("norm_") or k.startswith("A3_"))]
    z = float(stats.norm.ppf(1 - ALPHA_PHYS / len(keys)))
    nd = {rho: len([r for r in res if r["rho"] == rho]) for rho in (30.0, 60.0)}
    out, tab = {}, []
    for label, rho, row, v, twose, t in table:
        key = [k for k in keys if (label, rho) == _label_of(k)]
        if not key:
            continue
        se = twose / 2
        out[key[0]] = ceil_to(v + (2 + z) * se)
        tab.append((label, rho, row, v, se, out[key[0]]))
    return out, tab, z, signs


def _label_of(key):
    parts = key.split("_")
    rho = parts[-1]
    base = "_".join(parts[:-1])
    names = {"B_0.5": "residual max |meas - B|, A 0.5", "B_1": "residual max |meas - B|, A 1",
             "B_2": "residual max |meas - B|, A 2", "direct": "|direct check|", "skew": "|skewness|, t >= 1.5",
             "kurt": "|excess kurtosis|, t >= 1.5", "p": "growth exponent p",
             "lin_cont": "linear (A): max |amp-dev|, k 0-3", "nl_cont_0.5": "nonlinear (A) against the pendulum, A 0.5",
             "nl_cont_1": "nonlinear (A) against the pendulum, A 1", "nl_cont_2": "nonlinear (A) against the pendulum, A 2",
             "nl1_cont": "nonlinear (A) against sine-Gordon, A 1 cos x"}
    return names.get(base), rho


def verdict_stat(res, tol, nse=2.0):
    """Revision b, statistical form: (status, instrument dict, passed, failed, worst instrument item)."""
    st0, inst0, passc, failc, A = verdict(res, {**tol, **{k: 1e9 for k in _inst_keys()}}, nse)
    I = {rho: inst_items(res, rho) for rho in (30.0, 60.0)}
    ntot = sum(len(v[0]) for v in I.values())
    inst, worst = {}, (None, 0.0)
    for rho, (items, n) in I.items():
        tc = t_crit(ntot, n)
        for label, d, se in items:
            z = abs(d) / math.sqrt(se * se + E_FLOOR * E_FLOOR)
            inst[f"{label} rho {rho:g}"] = bool(z < tc)
            if z / tc > worst[1]:
                worst = (f"{label} rho {rho:g}: |d| {abs(d):.4f}, SE {se:.4f}, z {z:.2f} (t_crit {tc:.2f})", z / tc)
    if not all(inst.values()):
        st = "INSTRUMENT_FAILED"
    elif any(failc.values()):
        st = "FAILED"
    elif all(passc.values()):
        st = "PASSED"
    else:
        st = "UNDECIDED"
    return st, inst, passc, failc, worst[0]


def _inst_keys():
    return [f"eps_{g}_{i}" for g in ("30", "60") for i in range(NRB)] + ["norm_30", "norm_60", "A3_30", "A3_60"]


# ---- physics criteria as max over slabs of (|d_j| + 2 SE_j) (proposed 14 Nov, 15:30, before freezing) --------------
def _slabavg(dev, tg, tb, X):
    out = []
    for b in tb:
        m = np.abs(tg - b) < 0.25
        wt = np.clip(X - tg[m] - 0.3, 0, None) ** 3
        out.append(float(np.sum(dev[m] * wt) / np.sum(wt)))
    return np.array(out)


def _B_curve(w, ct, rho, A, tb, sd, X, F0=None):
    from diag.c0_ct_renorm import mean_field
    from diag.c0_apex import pendulum, F_k
    tg = np.arange(0.0, tb[-1] + 0.3, 0.005)
    F0 = F_k(w, rho, 0, tg) if F0 is None else F0
    s = np.interp(tg, np.r_[0.0, tb], np.r_[0.0, sd * A])
    return _slabavg((mean_field(w, ct, rho, A, tg, s, 0.005, F0) - pendulum(A, tg)) / A, tg, tb, X)


def phys_stats(R, nse=2.0, F0cache=None):
    """Per density: for every physics criterion the per-slab (d, SE) and M = max_j (|d_j| + nse SE_j)."""
    from diag.c0_apex import dyn_pred, dyn_pred_nl0, F_k
    rho = R[0]["rho"]
    X = R[0]["X"]
    tb = np.array(R[0]["t"])
    n = len(R)
    name, w, ct = rows(rho)[0]
    tg = np.arange(0.0, tb[-1] + 0.3, 0.005)
    F0 = F_k(w, rho, 0, tg) if F0cache is None else F0cache
    out = {}

    def put(key, d, se, sl=slice(None)):
        d, se = np.asarray(d)[sl], np.asarray(se)[sl]
        out[key] = (float(np.max(np.abs(d) + nse * se)), d.tolist(), se.tolist())

    for a in ("0.5", "1", "2"):
        D = np.array([r["k"][name]["modes"][f"nl0_A{a}"] for r in R])
        sd = D[:, :, 1].mean(0)
        B = _B_curve(w, ct, rho, float(a), tb, sd, X, F0)
        res = D[:, :, 0] - B
        put(f"B_{a}", res.mean(0), res.std(0, ddof=1) / math.sqrt(n))
        put(f"nl_cont_{a}", D[:, :, 0].mean(0), D[:, :, 0].std(0, ddof=1) / math.sqrt(n))
        P = np.array(dyn_pred_nl0(w, rho, float(a), tb, X=X, ct=ct))
        out[f"renorm_{a}"] = float(np.max(np.abs(B - P)))
    late = tb >= 1.5
    for what, col, sl in (("skew", 2, late), ("kurt", 3, late), ("direct", 4, slice(None))):
        best = None
        for a in ("0.5", "1", "2"):
            D = np.array([r["k"][name]["modes"][f"nl0_A{a}"] for r in R])[:, :, col]
            put(f"{what}_{a}", D.mean(0), D.std(0, ddof=1) / math.sqrt(n), sl)
            if best is None or out[f"{what}_{a}"][0] > best[0]:
                best = out[f"{what}_{a}"]
        out[what] = best
    for k in range(4):
        D = np.array([r["k"][name]["modes"][f"lin{k}"] for r in R])[:, :, 0]
        put(f"lin_cont_{k}", D.mean(0), D.std(0, ddof=1) / math.sqrt(n))
    out["lin_cont"] = max((out[f"lin_cont_{k}"] for k in range(4)), key=lambda q: q[0])
    D = np.array([r["k"][name]["modes"]["nl1_A1"] for r in R])[:, :, 0]
    put("nl1_cont", D.mean(0), D.std(0, ddof=1) / math.sqrt(n))
    ps = []
    for m in MODES:
        Dm = np.array([r["k"][name]["modes"][m[0]] for r in R])
        pp = [growth_(tb, Dm[i, :, 1]) for i in range(n)]
        ps.append((float(np.mean(pp) + nse * np.std(pp, ddof=1) / math.sqrt(n)), m[0]))
    out["p"] = (max(ps)[0], max(ps)[1], None)
    out["sd"] = {m[0]: np.array([r["k"][name]["modes"][m[0]] for r in R])[:, :, 1].tolist() for m in MODES}
    return out


def growth_(tb, sd, lo=1.5, hi=5.5):
    m = (tb >= lo) & (tb <= hi) & (sd > 0)
    return float(np.polyfit(np.log(tb[m]), np.log(sd[m]), 1)[0])


PHYS_KEYS = [("B_0.5", (30, 60)), ("B_1", (30, 60)), ("B_2", (30, 60)), ("direct", (30, 60)), ("skew", (30, 60)),
             ("kurt", (30, 60)), ("p", (30, 60)), ("lin_cont", (60,)), ("nl_cont_0.5", (60,)), ("nl_cont_1", (60,)),
             ("nl_cont_2", (60,)), ("nl1_cont", (60,))]


def _jack_renorm(R, a, F0):
    """Jackknife SE of renorm max|B - (A3)| at amplitude a (B depends on the run only through the pooled sd)."""
    n = len(R)
    vals = []
    for i in range(n):
        S = R[:i] + R[i + 1:]
        vals.append(phys_stats_renorm(S, a, F0))
    vals = np.array(vals)
    return float(math.sqrt((n - 1) / n * np.sum((vals - vals.mean()) ** 2)))


def phys_stats_renorm(R, a, F0):
    from diag.c0_apex import dyn_pred_nl0
    rho = R[0]["rho"]
    tb = np.array(R[0]["t"])
    name, w, ct = rows(rho)[0]
    D = np.array([r["k"][name]["modes"][f"nl0_A{a}"] for r in R])
    B = _B_curve(w, ct, rho, float(a), tb, D[:, :, 1].mean(0), R[0]["X"], F0)
    P = np.array(dyn_pred_nl0(w, rho, float(a), tb, X=R[0]["X"], ct=ct))
    return float(np.max(np.abs(B - P)))


def verdict3(res, tol, nse=2.0, F0s=None, jack=True):
    """Revision b, final form proposed for freezing. Returns (status, inst, passc, failc, info)."""
    from diag.c0_apex import F_k
    rhos = (30.0, 60.0)
    Rs = {rho: [r for r in res if r["rho"] == rho] for rho in rhos}
    tb = np.array(res[0]["t"])
    tg = np.arange(0.0, tb[-1] + 0.3, 0.005)
    F0s = F0s or {rho: F_k(rows(rho)[0][1], rho, 0, tg) for rho in rhos}
    I = {rho: inst_items(res, rho) for rho in rhos}
    ntot = sum(len(v[0]) for v in I.values())
    inst, info = {}, {"inst": []}
    for rho, (items, n) in I.items():
        tc = t_crit(ntot, n)
        for label, d, se in items:
            z = abs(d) / math.sqrt(se * se + E_FLOOR * E_FLOOR)
            inst[f"{label} rho {rho:g}"] = bool(z < tc)
            info["inst"].append((f"{label} rho {rho:g}", d, se, z, tc))
    PS = {rho: phys_stats(Rs[rho], nse, F0s[rho]) for rho in rhos}
    passc, failc = {}, {}
    for key, rr in PHYS_KEYS:
        for rho in rr:
            o = PS[float(rho)][key]
            k = f"{key}_{rho}"
            if key == "p":
                passc[k] = bool(o[0] < tol[k])
                failc[k] = False
            else:
                lo = float(np.max(np.abs(np.array(o[1])) - nse * np.array(o[2])))
                passc[k] = bool(o[0] < tol[k])
                failc[k] = bool(lo > tol[k])
            info[k] = o[0]
    a30, a60 = PS[30.0], PS[60.0]
    for k in range(4):
        d30, s30 = np.array(a30[f"lin_cont_{k}"][1]), np.array(a30[f"lin_cont_{k}"][2])
        d60, s60 = np.array(a60[f"lin_cont_{k}"][1]), np.array(a60[f"lin_cont_{k}"][2])
        j30, j60 = int(np.argmax(np.abs(d30))), int(np.argmax(np.abs(d60)))
        d = abs(d30[j30]) - abs(d60[j60])
        se = math.hypot(s30[j30], s60[j60])
        passc[f"conv_lin{k}"], failc[f"conv_lin{k}"] = bool(d - nse * se > 0), bool(d + nse * se < 0)
        info[f"conv_lin{k}"] = (d, se)
    j = int(np.argmin(np.abs(tb - 4.75)))
    for m in MODES:
        a = np.array(a30["sd"][m[0]])[:, j]
        b = np.array(a60["sd"][m[0]])[:, j]
        r = b.mean() / a.mean()
        se = r * math.sqrt((a.std(ddof=1) / a.mean()) ** 2 / a.size + (b.std(ddof=1) / b.mean()) ** 2 / b.size)
        passc[f"noise_{m[0]}"], failc[f"noise_{m[0]}"] = bool(1 - r - nse * se > 0), bool(1 - r + nse * se < 0)
        info[f"noise_{m[0]}"] = (1 - r, se)
    for a in ("1", "2"):
        d = a30[f"renorm_{a}"] - a60[f"renorm_{a}"]
        se = math.hypot(_jack_renorm(Rs[30.0], a, F0s[30.0]), _jack_renorm(Rs[60.0], a, F0s[60.0])) if jack else 0.0
        passc[f"renorm_{a}"], failc[f"renorm_{a}"] = bool(d - nse * se > 0), bool(d + nse * se < 0)
        info[f"renorm_{a}"] = (d, se)
    if not all(inst.values()):
        st = "INSTRUMENT_FAILED"
    elif any(failc.values()):
        st = "FAILED"
    elif all(passc.values()):
        st = "PASSED"
    else:
        st = "UNDECIDED"
    return st, inst, passc, failc, info


def derive3(res, nboot=200, seed=1):
    """tol = ceil(M_dev + z SE_M), SE_M from a bootstrap of the development sprinklings (the statistic itself)."""
    from scipy import stats
    from diag.c0_apex import F_k
    rhos = (30.0, 60.0)
    tb = np.array(res[0]["t"])
    tg = np.arange(0.0, tb[-1] + 0.3, 0.005)
    F0s = {rho: F_k(rows(rho)[0][1], rho, 0, tg) for rho in rhos}
    Rs = {rho: [r for r in res if r["rho"] == rho] for rho in rhos}
    keys = [(key, rho) for key, rr in PHYS_KEYS for rho in rr]
    z = float(stats.norm.ppf(1 - ALPHA_PHYS / len(keys)))
    M0 = {rho: phys_stats(Rs[rho], 2.0, F0s[rho]) for rho in rhos}
    rng = np.random.default_rng(seed)
    boots = {f"{k}_{r}": [] for k, r in keys}
    for _ in range(nboot):
        for rho in rhos:
            R = Rs[rho]
            S = [R[i] for i in rng.integers(0, len(R), len(R))]
            P = phys_stats(S, 2.0, F0s[rho])
            for key, rr in PHYS_KEYS:
                if int(rho) in rr:
                    boots[f"{key}_{int(rho)}"].append(P[key][0])
    tol, table = {}, []
    for key, rho in keys:
        k = f"{key}_{rho}"
        m0 = M0[float(rho)][key][0]
        se = float(np.std(boots[k], ddof=1))
        tol[k] = ceil_to(m0 + z * se)
        table.append((k, m0, se, tol[k]))
    return tol, table, z, F0s


if __name__ == "__main__":
    main()
