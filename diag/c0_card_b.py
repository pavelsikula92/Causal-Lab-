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
    a = ap.parse_args()
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


if __name__ == "__main__":
    main()
