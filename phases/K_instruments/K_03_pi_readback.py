"""
K.3 -- Reading Pi back from the density of chains (from GTS Lab 3.x E45c)
========================================================================

Calibrates the clock-rate reading core.geometry.chain_density in 2+1 and
3+1. In a sprinkling of uniform PROPER density, the longest chain of a
small diamond at rest of coordinate height H at a place with clock rate Pi
is proportional to its proper time Pi H (Brightwell-Gregory). So
Pi_order(r) = mean L(r) / mean over sites should follow Pi(r) / mean Pi.

A conformal sprinkling of the softened A5' profile (core.operators.
continuum.softened_point: 1/r in 3+1, log r in 2+1, Pi = 1 at |x| = R),
density proportional to Pi^D. At n_sites radii, n_diam small diamonds are
sprinkled (a Poisson process restricted to a subregion is the same process).
"""

from __future__ import annotations

import time

import numpy as np

from core.order.sprinkle import Interval, sprinkle_poisson, derive_seed
from core.geometry.chain_density import longest_chain_length, readback
from core.operators.continuum import softened_point
from core.stats.stats import fmt
from core.testing import Result

TEST = {
    "id": "K.3",
    "phase": "K",
    "title": "Reading Π back from the density of chains",
    "theory_part": "instrument of the clock rate Π for F0.5, F2 and F3; no physical conclusion",
    "question": "Does the length of the longest chain of small diamonds read back the profile of Π that was put into "
                "the density?",
    "method": "conformal sprinkling of a soft A5′ profile; n_diam diamonds of height H at n_sites radii; Pi_order = "
              "mean length / average over the sites",
    "prediction": "L is proportional to the proper time, so Pi_order(r) = Pi(r) / mean Pi",
    "criterion": "Pearson r >= r_min and the largest relative error <= err_max at every density of the plan. "
                 "Instrument: a flat sprinkling reads flat, |Pi_order - 1| <= flat_tol (otherwise "
                 "INSTRUMENT_FAILED). Null: the flat sprinkling compared with the field profile must not meet the "
                 "criterion (otherwise UNDECIDED).",
    "null_model": "the reading of a flat sprinkling compared with the field profile",
    "level": "calibration",
    "meters": [],
    "depends_on": [],
    "old_id": "E45c",
    "params": {"r_min": 0.95, "err_max": 0.10, "flat_tol": 0.05, "R": 0.5, "soft": 0.05, "H": 0.1,
               "n_sites": 9, "r_lo": 0.1, "r_hi": 0.42, "n_diam": 300},
    "parts": {"d3": {"D": 3, "eps": 0.2, "rhos": [1e5, 4e5, 1.6e6]},
              "d4": {"D": 4, "eps": 0.043, "rhos": [4e6, 1.6e7, 6.4e7]}},
    "quick": {"n_diam": 20},
    "quick_parts": {"d3": {"rhos": [1e5]}, "d4": {"rhos": [4e6]}},
    "estimate": "2+1: 10-30 min; 3+1: 1-3 h",
}


def site_lengths(D, rho, H, r, n_diam, Pi, rng):
    """Mean longest-chain length of n_diam diamonds of height H centred at radius r."""
    Ls = []
    for _ in range(n_diam):
        d = rng.normal(size=D - 1)
        d /= np.linalg.norm(d)
        reg = Interval(D, H, r * d)
        if Pi is None:
            X = sprinkle_poisson(reg, rho, rng)
        else:
            X = sprinkle_poisson(reg, rho, rng, density=lambda P: Pi(np.sqrt((P[:, 1:] ** 2).sum(axis=1))) ** D)
        Ls.append(longest_chain_length(X))
    return float(np.mean(Ls))


def run(p, ctx):
    D = int(p["D"])
    radii = [p["r_lo"] + (p["r_hi"] - p["r_lo"]) * k / (p["n_sites"] - 1) for k in range(int(p["n_sites"]))]
    Pi = softened_point(D - 1, p["eps"], p["soft"], p["R"])
    Pt = [float(Pi(r)) for r in radii]
    t0 = time.time()
    lines = [f"{D-1}+1: softened A5' profile eps={p['eps']} soft={p['soft']} R={p['R']}; diamonds H={p['H']}; "
             f"radii {' '.join(fmt(r, 2) for r in radii)}", f"  Pi(r): {' '.join(fmt(x) for x in Pt)}"]
    res, ok_all, null_any, flat_ok = {}, True, False, True
    for rho in p["rhos"]:
        rng = np.random.default_rng(derive_seed(ctx.seed, "K.3", D, rho))
        L = [site_lengths(D, rho, p["H"], r, int(p["n_diam"]), Pi, rng) for r in radii]
        Lf = [site_lengths(D, rho, p["H"], r, int(p["n_diam"]), None, rng) for r in radii]
        rb, nb = readback(L, Pt), readback(Lf, Pt)
        fdev = float(np.max(np.abs(np.asarray(Lf) / np.mean(Lf) - 1.0)))
        ok = rb["r"] >= p["r_min"] and rb["max_rel_err"] <= p["err_max"]
        nul = nb["r"] == nb["r"] and nb["r"] >= p["r_min"] and nb["max_rel_err"] <= p["err_max"]
        ok_all &= ok
        null_any |= nul
        flat_ok &= fdev <= p["flat_tol"]
        res[str(rho)] = {"L": L, "L_flat": Lf, "readback": rb, "null": nb, "flat_dev": fdev}
        lines.append(f"  rho={rho:.3g} (ell {rho ** (-1.0 / D):.4f}): r = {fmt(rb['r'])}, max err {fmt(rb['max_rel_err'])} "
                     f"-> {'meets' if ok else 'does NOT meet'}; flat dev {fmt(fdev)}; null r = {fmt(nb['r'])}")
        lines.append("     Pi_order " + " ".join(fmt(x) for x in rb["pi_order"]))
        ctx.progress(f"[K.3] {D-1}+1 rho={rho:.3g} r={fmt(rb['r'])} err={fmt(rb['max_rel_err'])} ({time.time()-t0:.0f}s)")
    if not flat_ok:
        status = "INSTRUMENT_FAILED"
    elif null_any:
        status = "UNDECIDED"
    else:
        status = "PASSED" if ok_all else "FAILED"
    lines.append(f"instrument (flat reads flat within {p['flat_tol']}): {'OK' if flat_ok else 'FAILS'}; "
                 f"null {'MEETS (not specific)' if null_any else 'fails (good)'}")
    return Result(status, lines, {"radii": radii, "Pi_true": Pt, "per_rho": res})
