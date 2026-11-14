"""
C0: the tolerances of the card from the development rows of diag/c0_card.py, by the rule of 10 Nov (14 Nov)
==========================================================================================================
tolerance = ceil( 1.25 x (|development value| + 2 SE_card) ),  SE_card = SE_dev x sqrt(n_dev / n_card),
rounded up to 0.001 below 0.02, to 0.005 below 1, to 0.05 above. The development value is the largest over the modes a
criterion covers (and over the slabs, as the readout defines it); 1.25 is the margin. Sign criteria (convergence, the
noise falling with the density, the renormalization falling with the density) have no tolerance: difference - 2 SE > 0.

    python -m diag.c0_card_tol diag/out/c0_card_dev_<stamp>.json
"""
from __future__ import annotations

import json
import math
import sys

import numpy as np

from diag import c0_card as M

N_CARD = {30.0: 16, 60.0: 4}
MAIN = "MV8+ct"


def ceil_to(v):
    step = 0.001 if v < 0.02 else (0.005 if v < 1 else 0.05)
    return round(math.ceil(v / step - 1e-9) * step, 6)


def derive(res):
    A = {rho: M.card_stats([r for r in res if r["rho"] == rho]) for rho in (30.0, 60.0)}
    ndev = {rho: len([r for r in res if r["rho"] == rho]) for rho in (30.0, 60.0)}
    rows = []

    def tol(label, rho, items):
        f = math.sqrt(ndev[rho] / N_CARD[rho])
        best = max(items, key=lambda q: q[1] + 2 * q[2] * f)
        v = best[1] + 2 * best[2] * f
        t = ceil_to(1.25 * v)
        rows.append((label, rho, best[0], best[1], best[2] * f, t))
        return t

    out = {}
    for rho in (30.0, 60.0):
        o = A[rho][MAIN]
        ms = o["modes"]
        out[f"eps_{rho:g}"] = tol("eps: max_r |eps - (A1)|", rho, [("all r", *o["eps"])])
        out[f"norm_{rho:g}"] = tol("|normalisation - 1|", rho, [("t >= 3", abs(o["norm"][0]), o["norm"][1])])
        out[f"A3_{rho:g}"] = tol("max |amp-dev - (A3)|, linear k 0-3 and nl A 0.5", rho,
                                 [(k, *ms[k]["A3"]) for k in ("lin0", "lin1", "lin2", "lin3", "nl0_A0.5")])
        for A_ in ("0.5", "1", "2"):
            out[f"B_{A_}_{rho:g}"] = tol(f"residual max |meas - B|, A {A_}", rho, [(f"nl0_A{A_}", *ms[f"nl0_A{A_}"]["B"])])
        out[f"direct_{rho:g}"] = tol("|direct check|", rho, [(k, *ms[k]["direct"]) for k in ("nl0_A0.5", "nl0_A1", "nl0_A2")])
        out[f"skew_{rho:g}"] = tol("|skewness|, t >= 1.5", rho, [(k, *ms[k]["skew"]) for k in ("nl0_A0.5", "nl0_A1", "nl0_A2")])
        out[f"kurt_{rho:g}"] = tol("|excess kurtosis|, t >= 1.5", rho, [(k, *ms[k]["kurt"]) for k in ("nl0_A0.5", "nl0_A1", "nl0_A2")])
        out[f"p_{rho:g}"] = tol("growth exponent p", rho, [(m[0], *ms[m[0]]["p"]) for m in M.MODES])
    o = A[60.0][MAIN]["modes"]
    out["lin_cont_60"] = tol("linear (A): max |amp-dev|, k 0-3", 60.0, [(k, *o[k]["cont"]) for k in ("lin0", "lin1", "lin2", "lin3")])
    for A_ in ("0.5", "1", "2"):
        out[f"nl_cont_{A_}_60"] = tol(f"nonlinear (A): max |amp-dev| against the pendulum, A {A_}", 60.0,
                                      [(f"nl0_A{A_}", *o[f"nl0_A{A_}"]["cont"])])
    out["nl1_cont_60"] = tol("nonlinear (A): max |amp-dev| against sine-Gordon, A 1 cos x", 60.0, [("nl1_A1", *o["nl1_A1"]["cont"])])
    signs = []
    a30, a60 = A[30.0][MAIN]["modes"], A[60.0][MAIN]["modes"]
    for k in ("lin0", "lin1", "lin2", "lin3"):
        d = a30[k]["cont"][0] - a60[k]["cont"][0]
        se = math.hypot(a30[k]["cont"][1], a60[k]["cont"][1])
        signs.append((f"convergence {k}: max|A|(30) - max|A|(60)", d, se))
    for m in M.MODES:
        r, se = M.sd_ratio(A[30.0], A[60.0], MAIN, m[0])
        signs.append((f"noise falls {m[0]}: 1 - sd(60)/sd(30) at t 4.75", 1 - r, se))
    for k in ("nl0_A1", "nl0_A2"):
        d = a30[k]["renorm"] - a60[k]["renorm"]
        se = math.hypot(a30[k]["B"][1], a60[k]["B"][1])
        signs.append((f"renormalization falls {k}: renorm(30) - renorm(60)", d, se))
    return out, rows, signs


def main(fn):
    res = json.load(open(fn))
    out, rows, signs = derive(res)
    print("criterion | rho | deciding row | development value | 2 SE at the card's n | tolerance")
    for label, rho, k, v, se, t in rows:
        print(f"{label} | {rho:g} | {k} | {v:.4f} | {2 * se:.4f} | {t:g}")
    print("sign criteria (development): difference, SE at n_dev")
    for label, d, se in signs:
        print(f"{label}: {d:+.4f} +- {se:.4f}")
    print(json.dumps(out))


if __name__ == "__main__":
    main(sys.argv[1])
