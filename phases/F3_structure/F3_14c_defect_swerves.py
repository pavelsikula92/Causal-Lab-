"""
F3.14c -- Swerves of a field defect: the momentum diffusion predicted by the shot noise of the order (1+1)
=========================================================================================================

Why this test exists (10 Nov 2026; card docs/F3_14c_card_draft.md; derivation docs/derivations/swerves_of_a_defect.md)
---------------------------------------------------------------------------------------------------------------------
A kink of F3.14a at rest, with no force, wanders. The wandering study (findings, Discoveries) found a random walk in
momentum with D_p ~ 1/rho, rising steeply with the mass. The derivation: in 1+1 the retarded Green function is exactly
1/2 inside the cone, so the rule phi(x) = phi_0(x) - (1/2 rho) sum_{y < x} V'(phi(y)) is exactly the continuum field with
point sources at the events. Its difference from the continuous source is Poisson shot noise, white in spacetime, of
variance V'^2 / rho. The momentum balance dP/dt = int s phi_x dx then gives, at leading order (the kink keeps its
profile),
    D_p = (1 / 2 rho) int V'(phi)^2 phi_x^2 dx = (32/15) mu^(5/2) / rho   (sine-Gordon, V = mu (1 - cos phi)),
with no free constant. This test checks the prediction in magnitude and in its scaling with rho and mu.

Model (as F3.14a, unchanged rule)
---------------------------------
Flat 1+1 Poisson sprinkling of density rho in the domain of dependence |x| <= X - t, X = T + 10, T = 30; kink at rest at
x = 0 (exact initial data through the free wave phi_0); phi on all events by the Green function of the order (field1d,
Fenwick engine, numba). Cells: mu in {0.25, 1} (M 4, 8; w 2, 1) x rho in {400, 800}; reported: mu = 4 (M 16, w 0.5,
coarsely resolved, rho w^2 = 100-200) at rho 800.

Readout (frozen)
----------------
The kink fit of F3.14a per slab of width 1 (centres 0.5 ... 29.5), tracking the centre x_c. The velocity in window k:
v_k = slope of x_c over the 5 slabs of the window (6 windows, centres 2.5, 7.5, ..., 27.5). Per sprinkling the
increments dv_k = v_k - v_0 (this removes the initial velocity spread). Var(dv_k) over sprinklings against dt_k =
centre_k - centre_0 (5 points); its slope is 2 D_v (a constant readout floor drops out). D_p,meas = M^2 D_v, M = 8 sqrt(mu).
R = D_p,meas rho / ((32/15) mu^(5/2)). SE of R by bootstrap over sprinklings (1000 resamples, fixed seed).
Readout check (paired): on every sprinkling an exact kink is moved along a synthetic random walk in velocity with the
predicted D_v and sampled at the same events; the estimator applied to its readout must return the D_v of the true
synthetic path (ratio readout / truth, paired on the same paths).

Criteria (proposed in the card; frozen when the author approves)
-----------------------------------------------------------------
INSTRUMENT_FAILED if: in any decided cell more than 2 % of the sprinklings lose a fit; or the readout check
  |readout / truth - 1| + 2 SE >= 0.1 in any decided cell.
PASSED: |R - 1| + 2 SE < 0.3 in each of the 4 cells, and |R(800) - R(400)| - 2 SE < 0.3 for each mu (no extra density
  dependence beyond 1/rho).
FAILED: |R - 1| - 2 SE > 0.3 in any cell.
UNDECIDED otherwise. Every condition on its own; nothing averaged; nothing chosen after the run.
Reported: mu = 4 at rho 800; the exponent of Var x(t) between dt 9 and 29; the mean drift; D_p from the cubic growth of
Var x as a second estimator.

Reading declared now
--------------------
PASSED: the wandering of a bound defect on the order is the shot noise with which the order samples its own binding, with
the predicted constant and the predicted scaling (1/rho, mu^(5/2), i.e. M^5 at M w = 8). The swerves of causal set theory
are then derived for a bound object, set by its internal structure (int V'^2 (d phi)^2), and zero in the vacuum.
FAILED: the leading-order shot-noise picture is wrong or incomplete in magnitude (the feedback on the profile or the
radiation matter). Limits: 1+1 only (in 3+1 the Green function sits on links); the potential is chosen; no number for real
particles follows.
"""

from __future__ import annotations

import math
import time

import numpy as np

from core.order.sprinkle import derive_seed
from core.stats.parallel import parallel_map
from core.testing import Result
from phases.F3_structure import field1d as F

TEST = {
    "id": "F3.14c",
    "phase": "F3",
    "title": "Swerves of a field defect: momentum diffusion from the shot noise of the order (1+1)",
    "theory_part": "F3 emergence of structure; matter as a defect of a field of the order; discreteness effects "
                   "(swerves); docs/derivations/swerves_of_a_defect.md",
    "question": "Does a kink of F3.14a at rest make a random walk in momentum with D_p = (1/2 rho) int V'^2 phi_x^2 dx = "
                "(32/15) mu^(5/2) / rho, the shot noise of the order sampling its own binding, in magnitude and in its "
                "scaling with rho and mu?",
    "method": "flat 1+1, T 30; kink at rest; velocities over windows of 5 slabs; Var of the velocity increments against "
              "time; R = D_p rho / prediction for mu 0.25, 1 at rho 400, 800; readout check on the exact kink",
    "prediction": "R = 1 in every cell; R independent of rho",
    "criterion": "INSTRUMENT: fit losses <= 2 %; injected walk |readout / truth - 1| + 2 SE < 0.1. PASSED: |R - 1| + 2 SE < 0.3 "
                 "in each cell and |R(800) - R(400)| - 2 SE < 0.3 for each mu. FAILED: |R - 1| - 2 SE > 0.3 in any "
                 "cell. Otherwise UNDECIDED (card docs/F3_14c_card_draft.md)",
    "null_model": "an exact kink moved along a synthetic random walk of the predicted D_v at the same events (readout check)",
    "level": "mechanism",
    "meters": [],
    "depends_on": ["F3.14b"],
    "params": {"mus": [0.25, 1.0], "rhos": [400.0, 800.0], "mu_rep": 4.0, "rho_rep": 800.0, "n_seeds": 1024,
               "T": 30.0, "window": 5, "tol_R": 0.3, "tol_floor": 0.1, "max_lost": 0.02, "boot_B": 1000, "n_se": 2.0},
    "quick": {"rhos": [100.0, 200.0], "rho_rep": 200.0, "n_seeds": 8, "T": 20.0, "boot_B": 100},
    "estimate": "about 30-40 min with 8 processes (numba required): 4 cells x 1024 runs of 0.6-1.2 million events, plus "
                "1024 runs of mu 4 at rho 800",
}


def task(p, rho, mu, k, base):
    T = float(p["T"])
    X = T + 10.0
    t, x = F.sprinkle_trapezoid(rho, X, T, derive_seed(base, "F3.14c", rho, mu, k))
    phi = F.evolve(t, x, F.kink_free(t, x, 0.0, mu), rho, mu=mu)
    # readout check: an exact kink moved along a synthetic random walk in velocity with the predicted D_v
    rng = np.random.default_rng(derive_seed(base, "F3.14c-inject", rho, mu, k))
    d_v = (32.0 / 15.0) * mu ** 2.5 / rho / (64.0 * mu)
    h = 0.01
    tg = np.arange(0.0, T + h, h)
    vg = np.concatenate([[0.0], np.cumsum(rng.normal(0.0, math.sqrt(2.0 * d_v * h), tg.size - 1))])
    xg = np.concatenate([[0.0], np.cumsum(0.5 * (vg[1:] + vg[:-1]) * h)])
    ex = 4.0 * np.arctan(np.exp(np.clip(math.sqrt(mu) * (x - np.interp(t, tg, xg)), -50, 50)))
    true_c = np.interp(np.arange(0.5, T, 1.0), tg, xg).tolist()

    def track(field):
        xc = []
        for a in np.arange(0.0, T - 1.0 + 1e-9, 1.0):
            g = xc[-1] if (xc and np.isfinite(xc[-1])) else 0.0
            r = F.fit_kink(t, x, field, a, a + 1.0, g, mu)
            xc.append(r["xc"] if r else float("nan"))
        return xc

    return {"rho": rho, "mu": mu, "k": k, "N": int(t.size), "xc": track(phi), "xc_ex": track(ex), "true": true_c}


def run_task(*a):
    return task(*a)


def _windows(xc, W):
    tm = np.arange(0.5, len(xc), 1.0)
    v, c = [], []
    for s in range(0, len(xc) - W + 1, W):
        v.append(np.polyfit(tm[s:s + W], xc[s:s + W], 1)[0])
        c.append(tm[s:s + W].mean())
    return np.array(v), np.array(c)


def _slope_var(dV, dt):
    var = dV.var(axis=0, ddof=1)
    return float(np.polyfit(dt, var, 1)[0])


def analyse(R, W, mu, rho, B, seed):
    X = np.array([r["xc"] for r in R])
    Xe = np.array([r["xc_ex"] for r in R])
    Xt = np.array([r["true"] for r in R])
    good = np.all(np.isfinite(X), axis=1) & np.all(np.isfinite(Xe), axis=1)
    lost = 1.0 - good.mean()
    X, Xe, Xt = X[good], Xe[good], Xt[good]
    V = np.array([_windows(x, W)[0] for x in X])
    Ve = np.array([_windows(x, W)[0] for x in Xe])
    Vt = np.array([_windows(x, W)[0] for x in Xt])
    c = _windows(X[0], W)[1]
    dt = c[1:] - c[0]
    dV, dVe, dVt = (V - V[:, [0]])[:, 1:], (Ve - Ve[:, [0]])[:, 1:], (Vt - Vt[:, [0]])[:, 1:]
    pred = (32.0 / 15.0) * mu ** 2.5 / rho
    M2 = 64.0 * mu
    scale = 2.0 * pred / M2                      # predicted slope of Var(dv) = 2 D_v
    Rv = _slope_var(dV, dt) / scale
    fl = _slope_var(dVe, dt) / _slope_var(dVt, dt)          # readout of the injected walk / its truth (paired)
    rng = np.random.default_rng(seed)
    bsR, bsF = [], []
    n = len(dV)
    for _ in range(B):
        idx = rng.integers(0, n, n)
        bsR.append(_slope_var(dV[idx], dt) / scale)
        bsF.append(_slope_var(dVe[idx], dt) / _slope_var(dVt[idx], dt))
    D = X - X[:, [0]]
    var = D.var(axis=0, ddof=1)
    tt = np.arange(0.5, X.shape[1], 1.0) - 0.5
    i_a, i_b = int(np.argmin(np.abs(tt - 9))), len(tt) - 1
    ex = math.log(math.sqrt(var[i_b]) / math.sqrt(var[i_a])) / math.log(tt[i_b] / tt[i_a]) if i_a > 0 else float("nan")
    m = tt >= 5
    A = np.column_stack([np.ones(m.sum()), tt[m] ** 2, tt[m] ** 3])
    cc = np.linalg.lstsq(A, var[m], rcond=None)[0]
    R_cube = 1.5 * cc[2] * M2 / pred
    return {"R": Rv, "R_se": float(np.std(bsR, ddof=1)), "floor": fl, "floor_se": float(np.std(bsF, ddof=1)),
            "lost": lost, "n": int(n), "exponent": ex, "drift": float(D[:, -1].mean()),
            "drift_se": float(D[:, -1].std(ddof=1) / math.sqrt(n)), "R_cube": float(R_cube)}


def run(p, ctx):
    if F._nb is None:
        raise SystemExit("F3.14c needs numba (pip install numba)")
    t0 = time.time()
    ns, nse, W = int(p["n_seeds"]), float(p["n_se"]), int(p["window"])
    cells = [(rho, float(mu)) for mu in p["mus"] for rho in p["rhos"]]
    rep = (float(p["rho_rep"]), float(p["mu_rep"]))
    base = ctx.seed
    jobs = [(p, rho, mu, k, base) for (rho, mu) in cells + [rep] for k in range(ns)]
    jobs.sort(key=lambda j: -j[1])
    res = parallel_map(run_task, jobs, ctx.workers, ctx.progress,
                       lambda i, r, d, n: f"[F3.14c] {d}/{n}: rho {r['rho']:g} mu {r['mu']:g} N {r['N']} "
                                          f"({time.time() - t0:.0f}s)")
    lines = [f"flat 1+1, kink at rest, T {p['T']:g}; {ns} sprinklings per cell; velocity windows of {W} slabs; "
             f"prediction D_p rho = (32/15) mu^(5/2)"]
    inst, passc, failc = {}, {}, {}
    out = {}
    for (rho, mu) in cells + [rep]:
        R = [r for r in res if r["rho"] == rho and r["mu"] == mu]
        a = analyse(R, W, mu, rho, int(p["boot_B"]), derive_seed(base, "F3.14c-boot", rho, mu))
        out[(rho, mu)] = a
        tag = "" if (rho, mu) in cells else " (reported)"
        lines.append(f"  mu {mu:g} (M {8 * math.sqrt(mu):g}) rho {rho:g}{tag}: R = {a['R']:.3f} +- {a['R_se']:.3f}; "
                     f"readout check (injected walk, readout / truth) {a['floor']:.3f} +- {a['floor_se']:.3f}; lost {a['lost']:.3f} (n {a['n']}); "
                     f"R from the cubic growth {a['R_cube']:.3f}; exponent {a['exponent']:.2f}; drift {a['drift']:+.3f} "
                     f"+- {a['drift_se']:.3f}")
        if (rho, mu) in cells:
            key = f"mu {mu:g} rho {rho:g}"
            inst[f"fits {key}"] = a["lost"] <= float(p["max_lost"])
            inst[f"readout {key}"] = abs(a["floor"] - 1.0) + nse * a["floor_se"] < float(p["tol_floor"])
            passc[f"R {key}"] = abs(a["R"] - 1.0) + nse * a["R_se"] < float(p["tol_R"])
            failc[f"R {key}"] = abs(a["R"] - 1.0) - nse * a["R_se"] > float(p["tol_R"])
    rhos = [float(r) for r in p["rhos"]]
    for mu in p["mus"]:
        a1, a2 = out[(rhos[0], float(mu))], out[(rhos[-1], float(mu))]
        d = a2["R"] - a1["R"]
        se = math.sqrt(a1["R_se"] ** 2 + a2["R_se"] ** 2)
        lines.append(f"  density scaling mu {mu:g}: R({rhos[-1]:g}) - R({rhos[0]:g}) = {d:+.3f} +- {se:.3f}")
        passc[f"scaling mu {mu:g}"] = abs(d) - nse * se < float(p["tol_R"])
    lines.append("instrument: " + "; ".join(f"{k}: {'OK' if v else 'FAILS'}" for k, v in inst.items()))
    lines.append("passed conditions: " + "; ".join(f"{k}: {'yes' if v else 'no'}" for k, v in passc.items()))
    lines.append("failed conditions: " + "; ".join(f"{k}: {'YES' if v else 'no'}" for k, v in failc.items()))
    if not all(inst.values()):
        status = "INSTRUMENT_FAILED"
    elif any(failc.values()):
        status = "FAILED"
    elif all(passc.values()):
        status = "PASSED"
    else:
        status = "UNDECIDED"
    lines.append(f"({time.time() - t0:.0f}s)")
    return Result(status, lines, {"instrument": inst, "passed": passc, "failed": failc,
                                  "cells": {f"{k[0]}_{k[1]}": v for k, v in out.items()}, "results": res})
