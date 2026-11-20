"""
C2 diagnostic cl14 (decides nothing): is the residual of the kink at rest beyond the C1 mean field at t 1.0 real?
=================================================================================================================
Window 1 (seeds 0-63 per rho, base 21000, diag/out/c2_scan; xl 8, xr 10, L 5, T 2.45, the C0 field of diag/c2_dev on
the curved 3+1 order Pi^2 = (x + 24)/24 and on the flat order): the kink at rest moves, relative to the continuum's
gravity part xc_ex, by the C1 mean field (diag/c2_meanfield) plus a residual
    R(t) = < xc_sim - xc_mf > / < xc_ex >      (curved seeds; the same events for sim and mf)
R(1.0) = +0.068, +0.052, +0.041 at rho 15, 30, 60; pooled +0.048 +- 0.016 (z 3.1). It was found after looking at many
cells (t 0.5-1.5, rest and +-v, three rho), it does not fall with rho, and a second run with another domain (xl 13,
xr 18, L 7, rho 30, 16 seeds) gave about -0.01 +- 0.04 at t 1.0. Undecided.

This run repeats window 1 at rho 60 unchanged (the same geometry, readout and mean field) on fresh sprinklings
(seeds 1000.., base 21000, so no sprinkling is shared with window 1), 256 curved + 256 flat by default.

Pre-registered (19 Nov, before the run). Decisive: R(1.0) on the curved arm, with s its SE over seeds.
    replicated       R(1.0) - 2 s > 0
    not replicated   R(1.0) + 2 s < 0.048   and   R(1.0) - 2 s <= 0
    undecided        otherwise
Control: the flat kink at rest has expectation exactly 0 (mirror symmetry of the engine); its residual
<xc_sim - xc_mf> / <xc_ex(curved)> must lie within 2 SE of 0, otherwise the run is flagged and the verdict above is
not read. Reported, not judged: R(0.5), R(1.5); the gravity part of sim and mf separately. The fit
(phases/F3_structure/field1d.fit_kink) reads only events with |x - guess| < 7, inside the event sets of both readouts,
so sim and mf are read on the same events (checked in the quick run).

    python -m diag.cl14_c2_window --workers 8                 (full run, resumable: finished seeds are skipped)
    python -m diag.cl14_c2_window --n 2 --workers 2 --quick   (mechanical check; decides nothing)
Output: diag/out/c2_cl14/ (json and npy per seed; summary cl14_summary.{txt,json}).
"""
from __future__ import annotations

import argparse
import json
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

RHO, ELL, XL, XR, L, T = 60.0, 24.0, 8.0, 10.0, 5.0, 2.45
BS = (0.5, 1.0, 1.5)
R_WINDOW1 = 0.048
SEED0 = 1000


def _out(quick):
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "c2_cl14_quick" if quick else "c2_cl14")


def _job(args):
    ell, seed, out = args
    import diag.c2_dev as C
    C.OUT = out
    f = os.path.join(out, f"cl14_{RHO:g}_{ell:g}_{seed}.json")
    if os.path.exists(f) and os.path.exists(os.path.join(out, f"phi_rest_{RHO:g}_{ell:g}_{seed}.npy")):
        return f, None
    r = C.task(RHO, ell, seed, ["rest"], XL, XR, L, T)
    tmp = f + ".tmp"
    json.dump(r, open(tmp, "w"), indent=1)
    os.replace(tmp, f)
    return f, r["sec"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=256)
    ap.add_argument("--workers", type=int, default=os.cpu_count() or 1)
    ap.add_argument("--quick", action="store_true")
    a = ap.parse_args()
    out = _out(a.quick)
    os.makedirs(out, exist_ok=True)
    seeds = list(range(SEED0, SEED0 + a.n))
    jobs = [(ell, s, out) for s in seeds for ell in (ELL, float("inf"))]
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs = [ex.submit(_job, j) for j in jobs]
        for k, fu in enumerate(as_completed(futs), 1):
            f, sec = fu.result()
            print(f"{k}/{len(jobs)} {os.path.basename(f)} " + ("(exists)" if sec is None else f"{sec:.0f} s")
                  + f" ({time.time() - t0:.0f} s)", flush=True)

    import diag.c2_meanfield as M
    res = {}
    for ell, geo in ((ELL, "curved"), (float("inf"), "flat")):
        fs = [os.path.join(out, f"phi_rest_{RHO:g}_{ell:g}_{s}.npy") for s in seeds]
        mf = M.fit_readout(RHO, ell, 0.0, fs, bs=BS)
        for s, f, rows in zip(seeds, fs, mf):
            r = json.load(open(os.path.join(out, f"cl14_{RHO:g}_{ell:g}_{s}.json")))
            sim = {round(x["t"], 1): (x["fit"]["xc"], x["fit_ex"]["xc"]) for x in r["arms"]["rest"]}
            for b, xm, xe, wm, we in rows:
                b = round(b, 1)
                res[(geo, s, b)] = dict(sim=sim[b][0] - sim[b][1], mf=xm - xe, xex=sim[b][1])

    lines, summ = [], {"pre_registered": "R(1.0): replicated if R - 2s > 0; not replicated if R + 2s < 0.048 and "
                                         "R - 2s <= 0; undecided otherwise; flat control within 2 SE of 0",
                       "n": a.n, "quick": a.quick, "rows": {}}
    lines.append(f"cl14 (C2 diagnostic, decides nothing): rho {RHO:g}, ell {ELL:g} and flat, xl {XL:g} xr {XR:g} "
                 f"L {L:g} T {T:g}; seeds {seeds[0]}..{seeds[-1]} (n {a.n}); mean +- SE over seeds, "
                 f"in units of <xc_ex> (curved)")
    verdict, flag = None, None
    for b in BS:
        b = round(b, 1)
        C = [res[("curved", s, b)] for s in seeds]
        Fl = [res[("flat", s, b)] for s in seeds]
        den = np.mean([c["xex"] for c in C])
        def ms(v):
            v = np.asarray(v) / den
            return float(v.mean()), float(v.std(ddof=1) / np.sqrt(v.size)) if v.size > 1 else float("nan")
        R = ms([c["sim"] - c["mf"] for c in C])
        Gs, Gm = ms([c["sim"] for c in C]), ms([c["mf"] for c in C])
        Fr = ms([f["sim"] - f["mf"] for f in Fl])
        summ["rows"][str(b)] = dict(R=R, sim=Gs, mf=Gm, flat_res=Fr, xc_ex=den)
        lines.append(f"  t {b}: R {R[0]:+.4f} +- {R[1]:.4f}  "
                     f"sim {Gs[0]:+.4f} +- {Gs[1]:.4f}  mf {Gm[0]:+.4f} +- {Gm[1]:.4f}  "
                     f"flat res {Fr[0]:+.4f} +- {Fr[1]:.4f}  (<xc_ex> {den:+.5f})")
        if b == 1.0:
            flag = abs(Fr[0]) > 2 * Fr[1]
            if R[0] - 2 * R[1] > 0:
                verdict = "replicated"
            elif R[0] + 2 * R[1] < R_WINDOW1 and R[0] - 2 * R[1] <= 0:
                verdict = "not replicated"
            else:
                verdict = "undecided"
    lines.append(f"flat control at t 1.0 within 2 SE of 0: {'NO (flagged; the verdict is not read)' if flag else 'yes'}")
    lines.append(f"pre-registered reading of R(1.0) (window 1: +{R_WINDOW1}): {verdict}"
                 + ("  [quick run: mechanical only]" if a.quick else ""))
    lines.append(f"cost: {time.time() - t0:.0f} s")
    summ.update(verdict=verdict, flat_flag=bool(flag))
    txt = "\n".join(lines)
    print(txt)
    open(os.path.join(out, "cl14_summary.txt"), "w").write(txt + "\n")
    json.dump(summ, open(os.path.join(out, "cl14_summary.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
