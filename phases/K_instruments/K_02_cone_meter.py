"""
K.2 -- The causal cone as a meter of bending (from GTS Lab 3.x E47)
==================================================================

Calibrates the cone-opening instrument (core.geometry.cone_front). It does
NOT test whether hypothesis H1 holds; F7.1 uses the calibrated instrument.

Three sprinklings carry the same clock-rate field of a mass at the origin,
Pi = 1 - eps / sqrt(r^2 + soft^2) (weak: eps / b <= 0.04 on every ray):

    (a) weights    flat order, element weights Pi       relation: flat
    (b) conformal  density proportional to Pi^D          relation: flat
    (c) H1         ds^2 = Pi^2 dt^2 - Pi^-2 dx^2          local light speed Pi^2

A plane source slab at x = -L, a screen slab at x = +L, the mass between
them; impact parameter = screen coordinate y. Per screen cell the front
time of the causal future of the source is read with the representation's
relation and with the flat one on the same screen events; the difference
is the delay, its tilt the deflection.

GR prediction for a plane wave (weak field, gamma = 1):
    delay(y) = 2 eps int dx / sqrt(x^2 + y^2 + soft^2) over the path
    alpha(y) = -d delay / dy
"""

from __future__ import annotations

import time

import numpy as np

from core.order.sprinkle import Box, derive_seed
from core.order.relation import flat_future_of_set, h1_future_of_set
from core.geometry.cone_front import front_profile, local_tilt
from core.stats.stats import fmt
from core.testing import Result

TEST = {
    "id": "K.2",
    "phase": "K",
    "title": "The causal cone as an instrument of bending",
    "theory_part": "instrument for F7 (H1, Shapiro); no physical conclusion",
    "question": "Does the front of the causal future of a source measure the delay and the bending that GR predicts "
                "for the H1 metric?",
    "method": "a source slab and a screen slab, the field of a point mass between them; the front of the future in "
              "the cells of the screen; three representations of the field (weights, conformal density, H1 relation)",
    "prediction": "weights and conformal density: zero delay and zero bending (the relation is flat); H1: delay and "
                  "bending as in GR (gamma = 1), analytically from the integral along the path",
    "criterion": "(a) and (b): max |delay| and max |slope| <= zero_abs; (c): mean delay/GR and bending/GR within 1 "
                 "+- tol. Instrument: the H1 relation at eps = 0 reads an exactly flat front (otherwise "
                 "INSTRUMENT_FAILED).",
    "null_model": "representations (a) and (b), which must not give bending",
    "level": "calibration",
    "meters": [],
    "depends_on": [],
    "old_id": "E47",
    "params": {"tol": 0.2, "zero_abs": 1e-12, "eps": 0.003, "soft": 0.02, "L": 1.0, "w": 0.01, "t_src": 0.01,
               "y_lo": 0.08, "y_hi": 0.4, "y_pad": 0.3, "dy": 0.01, "t_extra": 0.08, "fit_half": 2,
               "n_src": 4000, "n_scr": 30000, "n_seeds": 8},
    "parts": {"d3": {"D": 3}, "d4": {"D": 4}},
    "quick": {"n_src": 600, "n_scr": 4000, "n_seeds": 1},
    "estimate": "10-20 min per part",
}


def gr_delay(y, eps, soft, x1, x2):
    """2 eps [asinh(x2/B) - asinh(x1/B)], B = sqrt(y^2 + soft^2)."""
    B = np.sqrt(np.asarray(y, float) ** 2 + soft * soft)
    return 2.0 * eps * (np.arcsinh(x2 / B) - np.arcsinh(x1 / B))


def gr_alpha(y, eps, soft, x1, x2):
    """-d gr_delay / dy (analytic)."""
    y = np.asarray(y, float)
    B2 = y * y + soft * soft
    return 2.0 * eps * y / B2 * (x2 / np.sqrt(x2 * x2 + B2) - x1 / np.sqrt(x1 * x1 + B2))


def slabs(D, p, rng):
    """Source and screen events (uniform), time-sorted."""
    L, w, dts = p["L"], p["w"], p["t_src"]
    y0, y1, pad = p["y_lo"], p["y_hi"], p["y_pad"]
    lo_s = [0.0, -L, y0 - pad] + [-w] * (D - 3)
    hi_s = [dts, -L + w, y1 + pad] + [w] * (D - 3)
    t_a, t_b = 2 * L - 2 * w - 0.01, 2 * L + p["t_extra"]
    lo_q = [t_a, L - w, y0] + [-w] * (D - 3)
    hi_q = [t_b, L, y1] + [w] * (D - 3)
    S = Box(D, lo_s, hi_s).propose(rng, int(p["n_src"]))
    Q = Box(D, lo_q, hi_q).propose(rng, int(p["n_scr"]))
    return S[np.argsort(S[:, 0])], Q[np.argsort(Q[:, 0])]


def run(p, ctx):
    D = int(p["D"])
    eps, soft, L = p["eps"], p["soft"], p["L"]
    edges = np.arange(p["y_lo"], p["y_hi"] + 1e-12, p["dy"])
    yc = 0.5 * (edges[:-1] + edges[1:])
    x1, x2 = -L + 0.5 * p["w"], L - 0.5 * p["w"]
    dpred, apred = gr_delay(yc, eps, soft, x1, x2), gr_alpha(yc, eps, soft, x1, x2)

    def Pi(P):
        return 1.0 - eps / np.sqrt((P ** 2).sum(axis=-1) + soft * soft)

    margin = 1.5 * float(gr_delay(np.array([0.0]), eps, soft, -2 * L, 2 * L)[0]) + 1e-6
    t0 = time.time()
    acc = {k: [] for k in ("a", "b", "c", "flat0")}
    for rep in range(int(p["n_seeds"])):
        rng = np.random.default_rng(derive_seed(ctx.seed, "K.2", D, rep))
        S, Q = slabs(D, p, rng)
        yq = Q[:, 2]
        base = front_profile(Q, flat_future_of_set(S, Q), yq, edges)
        acc["a"].append(front_profile(Q, flat_future_of_set(S, Q), yq, edges) - base)
        ks = rng.random(len(S)) < Pi(S[:, 1:]) ** D
        kq = rng.random(len(Q)) < Pi(Q[:, 1:]) ** D
        Sb, Qb = S[ks], Q[kq]
        acc["b"].append(front_profile(Qb, flat_future_of_set(Sb, Qb), Qb[:, 2], edges)
                        - front_profile(Qb, flat_future_of_set(Sb, Qb), Qb[:, 2], edges))
        acc["c"].append(front_profile(Q, h1_future_of_set(S, Q, Pi, margin), yq, edges) - base)
        f0 = h1_future_of_set(S, Q, lambda P: np.ones(P.shape[:-1]), 1e-9)
        acc["flat0"].append(front_profile(Q, f0, yq, edges) - base)
        ctx.progress(f"[K.2] {D-1}+1 seed {rep+1}/{p['n_seeds']} ({time.time()-t0:.0f}s)")
    prof = {k: np.nanmean(np.array(v), axis=0) for k, v in acc.items()}
    tilt = {k: local_tilt(yc, prof[k], int(p["fit_half"])) for k in prof}
    ok = np.isfinite(prof["c"]) & np.isfinite(tilt["c"])
    rd = float(np.mean(prof["c"][ok] / dpred[ok])) if ok.any() else float("nan")
    ra = float(np.mean(tilt["c"][ok] / apred[ok])) if ok.any() else float("nan")
    zero = {k: (float(np.nanmax(np.abs(prof[k]))), float(np.nanmax(np.abs(tilt[k])))) for k in ("a", "b", "flat0")}
    zero_ok = all(d <= p["zero_abs"] and a <= p["zero_abs"] for d, a in zero.values())
    c_ok = abs(rd - 1) <= p["tol"] and abs(ra - 1) <= p["tol"]
    inst_ok = zero["flat0"][0] == 0.0 and zero["flat0"][1] == 0.0
    status = "INSTRUMENT_FAILED" if not inst_ok else ("PASSED" if zero_ok and c_ok else "FAILED")
    k8 = max(1, len(yc) // 8)
    lines = [f"{D-1}+1: eps={eps} soft={soft} L={L}; {len(yc)} cells of {p['dy']} in y=[{p['y_lo']}, {p['y_hi']}]; seeds {p['n_seeds']}",
             "  y        " + " ".join(fmt(v, 2) for v in yc[::k8]),
             "  GR delay " + " ".join(fmt(v, 4) for v in dpred[::k8]),
             "  (c) H1   " + " ".join(fmt(v, 4) for v in prof["c"][::k8]),
             "  GR alpha " + " ".join(fmt(v, 4) for v in apred[::k8]),
             "  (c) tilt " + " ".join(fmt(v, 4) for v in tilt["c"][::k8]),
             f"  (c) delay / GR = {fmt(rd)}, deflection / GR = {fmt(ra)} (tol {p['tol']})",
             f"  (a) weights max |delay|, |tilt| = {fmt(zero['a'][0], 5)}, {fmt(zero['a'][1], 5)}; "
             f"(b) conformal {fmt(zero['b'][0], 5)}, {fmt(zero['b'][1], 5)} (zero within {p['zero_abs']})",
             f"  instrument: H1 relation at eps = 0 reads exactly flat: {'OK' if inst_ok else 'FAILS'}"]
    return Result(status, lines, {"y": yc, "gr_delay": dpred, "gr_alpha": apred, "profiles": prof,
                                  "tilts": tilt, "ratio_delay": rd, "ratio_alpha": ra, "zero": zero})
