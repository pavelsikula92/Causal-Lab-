"""
B2 development rows: walls in the prescribed curved order of the theory (decides nothing)
=========================================================================================
Order: diag/b2_curved.py (Pi^2 = (x + ell)/ell, exact relation, proper density rho, domain of dependence of the
rectangle -xl <= x <= xr, |y| <= Y). Field rule exactly as F3.15 (chain weight, variant B, frozen constants). The free
wave phi_0 and the exact reference come from the y-independent continuum solver on a grid (dx 0.005).
Arms (each on its own order: ell finite = curved, ell inf = flat):
  kink_rest   the kink wall at rest at x0 = 0 (continuum: exactly static; a test particle falls -t^2/(4 ell))
  kink_v0.3   launched at proper velocity 0.3 (continuum: constant proper speed, coordinate speed Pi^2 v)
  br0.7       breather wall omega 0.7 at maximal amplitude (falls; tension / energy 0.22)
Readouts on |y| < 2, slabs of width 1, centres 0.5 ... T - 1.5 (slab rule). Kink: fit centre, width, Q (as F3.15).
Breather: the centroid of the noise-debiased coordinate energy density over |x| < 7 (local regressions as F3.15).
Every readout is repeated on the exact continuum field at the same events (paired).

    python -m diag.b2_dev --rho 30 --seeds 4 --workers 2
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import time

import numpy as np

from core.order.sprinkle import derive_seed
from core.stats.parallel import parallel_map, default_workers
from phases.F3_structure import field1d as F
from phases.F3_structure.F3_15_wall_inertia_2p1 import _evolve, energy_db
from diag import b0_chain_field as B0
from diag.b2_curved import Continuum, sprinkle_curved, links_curved, tau_cut

CALIB = (0.2732723632276477, 0.21313211596642687, 1.2025922117299621)
GEOM = {"xl": 13.0, "xr": 22.0, "Y": 12.5, "T": 8.0}
XL_FLAT = 16.0                     # the flat domain of dependence closes faster on the left (development, 12 Nov)


def e_centroid(t, x, phi, tm, half=7.0, cell=0.5, step=0.25, x0=0.0):
    s = np.abs(t - tm) < 0.5 * cell
    ts, xs, ps = t[s], x[s], phi[s]
    o = np.argsort(xs)
    ts, xs, ps = ts[o], xs[o], ps[o]
    num = den = 0.0
    for xg in np.arange(x0 - half, x0 + half + 1e-9, step):
        i0, i1 = np.searchsorted(xs, [xg - 0.5 * cell, xg + 0.5 * cell])
        n = i1 - i0
        if n < 8:
            return float("nan"), float("nan")
        A = np.column_stack([np.ones(n), ts[i0:i1] - tm, xs[i0:i1] - xg])
        c, *_ = np.linalg.lstsq(A, ps[i0:i1], rcond=None)
        r = ps[i0:i1] - A @ c
        C = float(r @ r) / (n - 3) * np.linalg.inv(A.T @ A)
        e = 0.5 * (c[1] ** 2 - C[1, 1]) + 0.5 * (c[2] ** 2 - C[2, 2]) + (1.0 - math.cos(c[0]) * math.exp(0.5 * C[0, 0]))
        num += xg * e
        den += e
    return num / den, den * step


def task(rho, ell, k, arm, base):
    G = dict(GEOM, xl_flat=XL_FLAT)
    T = G["T"]
    sd = derive_seed(base, "B2dev", rho, ell if np.isfinite(ell) else -1.0, k)
    xl = G["xl"] if np.isfinite(ell) else G["xl_flat"]
    t, x, y = sprinkle_curved(rho, ell, xl, G["xr"], G["Y"], T, sd)
    om = float(arm[2:]) if arm.startswith("br") else None
    v = float(arm.split("_v")[1]) if "_v" in arm else 0.0
    a0 = time.time()
    Cx = Continuum(ell, 1.0, v, 0.0, T, -xl, G["xr"], omega=om)
    Cf = Continuum(ell, 1.0, v, 0.0, T, -xl, G["xr"], omega=om, mu_dyn=0.0)
    phi0, exact = Cf.at(t, x), Cx.at(t, x)
    sec_c = time.time() - a0
    ptr, idx = links_curved(t, x, y, ell if np.isfinite(ell) else -1.0, tau_cut(rho))
    sec_l = time.time() - a0 - sec_c
    phi = _evolve(ptr, idx, t, phi0, float(rho), 1.0, B0.weight_table("B", rho, CALIB), 0.0, 0.0, 0.0, 0.0, False, 0)
    sec = time.time() - a0
    band = np.abs(y) < 2.0
    tb, xb, pb, eb = t[band], x[band], phi[band], exact[band]
    rows = []
    for a in np.arange(0.0, T - 2.0 + 1e-9, 1.0):
        tm = a + 0.5
        row = {"t": tm}
        if om is None:
            xg = Cx.centre(tm)
            prev = rows[-1].get("xc") if rows else None
            guess = prev if (prev is not None and np.isfinite(prev)) else xg
            s = (tb >= a) & (tb < a + 1.0)
            r = F.fit_kink(tb[s], xb[s], pb[s], a, a + 1.0, guess, 1.0)
            rx = F.fit_kink(tb[s], xb[s], eb[s], a, a + 1.0, xg, 1.0)
            row.update({"xc": r["xc"] if r else float("nan"), "w": r["w"] if r else float("nan"),
                        "Q": r["Q"] if r else float("nan"), "xc_ex": rx["xc"] if rx else float("nan"),
                        "w_ex": rx["w"] if rx else float("nan"), "xc_grid": xg})
        else:
            c, e = e_centroid(tb, xb, pb, tm)
            cx, ex = e_centroid(tb, xb, eb, tm)
            c4, _ = e_centroid(tb, xb, pb, tm, half=4.0)
            cx4, _ = e_centroid(tb, xb, eb, tm, half=4.0)
            row.update({"c4": c4, "c4_ex": cx4})
            row.update({"c": c, "E": e, "c_ex": cx, "E_ex": ex, "amp": float(np.max(np.abs(pb[(np.abs(tb - tm) < 0.1)
                                                                                         & (np.abs(xb) < 4)])))})
        rows.append(row)
    os.makedirs(os.path.join("diag", "out", "b2_dev_fields"), exist_ok=True)
    np.savez_compressed(os.path.join("diag", "out", "b2_dev_fields", f"r{rho:g}_{ell}_{k}_{arm}.npz"), t=tb, x=xb,
                        y=y[band], phi=pb, exact=eb)
    return {"rho": rho, "ell": ell if np.isfinite(ell) else "inf", "k": k, "arm": arm, "N": int(t.size),
            "links": int(idx.size), "sec": sec, "sec_cont": sec_c, "sec_links": sec_l, "rows": rows}


def cached(rho, ell, k, arm, base):
    tag = [rho, str(ell), k, arm, base, GEOM, CALIB]
    if arm.startswith("br"):
        tag.append("v2")           # breather: half-window 4 added, fields saved, flat xl 16
    key = hashlib.sha256(json.dumps(tag).encode()).hexdigest()[:20]
    fn = os.path.join("diag", "out", "b2_dev_cache", key + ".json")
    if os.path.exists(fn):
        return json.load(open(fn))
    out = task(rho, ell, k, arm, base)
    os.makedirs(os.path.dirname(fn), exist_ok=True)
    json.dump(out, open(fn + ".tmp", "w"))
    os.replace(fn + ".tmp", fn)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rho", type=float, default=30.0)
    ap.add_argument("--seeds", type=int, default=4)
    ap.add_argument("--k0", type=int, default=0)
    ap.add_argument("--ell", type=float, default=24.0)
    ap.add_argument("--arms", type=str, default="kink_rest,br0.7,kink_v0.3,flat:br0.7")
    ap.add_argument("--workers", type=int, default=0)
    a = ap.parse_args()
    jobs = []
    for k in range(a.k0, a.k0 + a.seeds):
        for arm in a.arms.split(","):
            ell = float("inf") if arm.startswith("flat:") else a.ell
            jobs.append((a.rho, ell, k, arm.replace("flat:", ""), 20261112))
    t0 = time.time()
    res = parallel_map(cached, jobs, default_workers(a.workers), print,
                       lambda i, r, d, n: f"[B2 dev] {d}/{n}: {r['arm']} ell {r['ell']} k {r['k']} N {r['N']} "
                                          f"{r['sec']:.0f}s ({time.time() - t0:.0f}s)")
    os.makedirs(os.path.join("diag", "out"), exist_ok=True)
    fn = os.path.join("diag", "out", f"b2_dev_{time.strftime('%Y%m%dT%H%M%S', time.gmtime())}.json")
    json.dump({"args": vars(a), "geom": GEOM, "calib": CALIB, "results": res}, open(fn, "w"))
    print("saved", fn)


if __name__ == "__main__":
    main()
