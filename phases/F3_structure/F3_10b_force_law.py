"""
F3.10b -- Force law of two bundles: does the rule turn the partner's field into the force lambda F grad phi?
==========================================================================================================

Why this test exists (23 Oct 2026, after F3.10 law and the author's QA review)
------------------------------------------------------------------------------
F3.10 law (seed 0) ended UNDECIDED: every separation attracts, the centroid
exponent k = 0.77 +- 0.57 against 1.11 of a quasi-static point-line
prediction. The QA review found that this prediction is not valid where it
was used: with the calibrated lambda and mu the predicted point lines pass
through each other at d = 8 and close from 12 to 5.5 at d = 12; it puts the
current separation into a static-source formula and treats a bundle of
radius ~4-6 as a line. The momentum readout (a_P = dP/dt / M), about five
times more precise than the centroid, gave ratios to that approximation
of 0.61 / 0.76 / 0.91 / 0.98 for d = 8 / 12 / 16 / 20. F3.10 law stays
UNDECIDED (R132: its criteria are not changed); this is a new test with
its own criteria, declared before any of its runs.

What is tested, and what is not
-------------------------------
The 2+1 kernel 1/(2 pi tau) of the field and the sign of the coupling are
INPUT (F3.8). A law 1/r is therefore not derived here. The test asks
whether the discrete dynamics turns the field into force faithfully:
  force  for every d, the momentum the bundles receive equals the first-
         order prediction lambda sum F d_x phi_partner computed from the
         actual discrete field of the actual partner;
  law    the dependence of the received momentum on d has the exponent of
         that prediction (the shape of the force law follows from the kernel
         and the source, and the dynamics does not distort it).
The two hypotheses get separate verdicts: a constant wrong coefficient
(force FAILED, law PASSED) is a different finding from a distorted law.

Model and measurement (locked)
------------------------------
Exactly the F3.10 model and code (module F3_10_two_bundles): two species of
flux, each with its carried direction, its children choice and its retarded
field phi_i = sum F / (2 pi tau), tau >= 1; bending by the least-squares
gradient of ln Pi = -lambda sum phi over the parents (lag <= 2, >= 4
parents). epsilon = 0.1 locked; lambda = epsilon / phi_ref from the SAME
calibration as F3.10 (the same calibration streams and box, T = 48), so
lambda is the F3.10 value. Separations d = 8, 12, 16, 20, f = 1, 1; the box
ends at T = 2.5 d + 6, just after the window (the transport is causal, so
nothing later acts on the window), margin 12, radius 3, source time
[0, 4); independent sprinklings for each d (no pairing across d: the box
depends on d). 64 sprinklings per d: from the F3.10 law data (momentum
readout ~6-8 % per d with 24 sprinklings and a longer window), 32 would
leave 2 SE of the exponent near 0.3, above the tolerance 1/4.
Branches per sprinkling: A_L and A_R (each bundle alone), B (both, common
field). Instrument on every sprinkling (as F3.10 after the QA review; it
costs a fraction of a second): C = A bit for bit (choices, flux, directions), mirror x -> -x
and frame beta' = 0.6 (choices, flux bit for bit, directions mapped back
within 1e-9; a difference starting at a numerical tie, margin < 1e-9, is
reported as a tie).
Window (linear response, the same in units of d): [1.25 d, 2.5 d]. The
retarded build-up of the field of a line has the same shape in t/d, and
the displacement stays a small fraction of d. Validity is checked on the
data: the mean change of the separation over the window (B - A) must stay
within d/4, otherwise that d is outside linear response and no verdict
beyond UNDECIDED is given.
Measured: dP_i(t) = P_i(B) - P_i(A) per slab of 2, P = sum F u_x over the
carriers of species i; the slope of a linear fit over the window; the
received momentum rate of the pair Pdot = (Pdot_L - Pdot_R) / 2 (towards
each other positive).
Predicted (first order in the partner's field, from branch A only, no B
data): along the actual flow of a bundle alone (its choices, fluxes and
proper times of the chosen links), the kick of every link is lambda
d_x phi_partner(v) tau_link, with phi_partner the discrete retarded field
of the partner alone (sum over its flux events in the causal past, tau >= 1)
and d_x its symmetric finite difference with step h = 1 (the coarse-grained
gradient at the scale where the rule estimates it; the cut tau >= 1 moves
with the point, as in the continuum G_ret). The kicks are carried along the
flow (flux-weighted at mergers) to every event; P_pred per slab is the sum
of F times the carried kick; the same linear fit gives Pdot_pred. lambda is
the model's constant, the same in the simulation and in the prediction,
so the calibration noise of phi_ref does not enter the comparison.
Reported also: the continuum static-line exponent in the same windows
(k_line) and the kernel's 1 (the shape of the source and the retarded
build-up are what separates k_pred from 1).

Parts (verdicts; replication with --seed 1)
-------------------------------------------
force  ratio rho_d = mean Pdot / mean Pdot_pred for every d (bootstrap over
       sprinklings, measurement and prediction resampled together).
       PASSED if every d is resolved (Pdot - 2 SE > 0) and inside linear
       response and |rho_d - 1| + 2 SE <= 1/4 for every d; FAILED if
       |rho_d - 1| - 2 SE > 1/4 for some d inside linear response;
       UNDECIDED otherwise. Alternative: a coefficient off by a factor 2
       (rho = 1/2 or 2); a quarter of the distance to the nearer one.
law    exponent of Pdot against d: k_meas, and of Pdot_pred: k_pred (weighted
       least squares of ln on ln d over the four d); Delta k = k_meas -
       k_pred with the paired bootstrap. PASSED if every d is resolved and
       inside linear response and |Delta k| + 2 SE <= 1/4; FAILED if
       |Delta k| - 2 SE > 1/4; UNDECIDED otherwise. Alternative: a law one
       power off (Delta k = +-1); a quarter of that distance.
The output names the reason of an UNDECIDED: a d not resolved, a d outside
linear response, or a standard error too large (2 SE > 1/4 alone).
INSTRUMENT_FAILED if C, the mirror or the frame check fails.

Seen before this test was written: F3.10 attraction and law at seed 0
(above), the F3.10 quick runs, the QA review. Nothing of F3.10b has been run.
"""

from __future__ import annotations

import math
import time

import numpy as np

from core.order.sprinkle import derive_seed
from core.stats.parallel import parallel_map
from core.testing import Result
from phases.F3_structure.F3_10_two_bundles import (G_ret, Xm_of, boost_points, calibration_task, compare, lin_slope, make,
                                                   mink, species_series, transport_species, tube)

TEST = {
    "id": "F3.10b",
    "phase": "F3",
    "title": "Force law of two bundles: does the rule turn the partner's field into the force λ F ∇φ?",
    "theory_part": "fidelity of the conversion of a field into a force in the discrete dynamics; the shape of the "
                   "force law from the kernel of the field and the shape of the source (the kernel 1/(2πτ) is an "
                   "input, not an output)",
    "question": "Do the bundles at d = 8 / 12 / 16 / 20 receive a momentum equal to the first-order prediction from "
                "the field of the actual partner, and does the dependence on d have the exponent of that prediction?",
    "method": "model and code of F3.10; λ from the same calibration as F3.10; window of linear response [1.25 d, 2.5 "
              "d]; the measured rate of momentum transfer (B − A) against the prediction from branch A (the kick λ "
              "∂_x φ_partner τ_link carried along the actual flow of the lone bundle, gradient of the partner's "
              "discrete field by central difference h = 1)",
    "prediction": "ρ_d = 1 for every d (force) and k_meas = k_pred (law)",
    "criterion": "force: PASSED all d resolved, in linear response and |ρ_d − 1| + 2 SE <= 1/4; FAILED |ρ_d − 1| − 2 "
                 "SE > 1/4 for some d. law: PASSED |k_meas − k_pred| + 2 SE <= 1/4 (all d resolved and in linear "
                 "response); FAILED |Δk| − 2 SE > 1/4. Linear response: mean change of the separation in the window "
                 "<= d/4. INSTRUMENT: C, mirror, frame (every sprinkling). Otherwise UNDECIDED with the reason. "
                 "Replication --seed 1.",
    "null_model": "branch A (the bundles alone); the first-order prediction from branch A",
    "level": "mechanism",
    "meters": ["K.6"],
    "depends_on": ["F3.10"],
    "params": {"D": 3, "T": 56.0, "T_cal": 48.0, "margin": 12.0, "r_src": 3.0, "t_src": 4.0, "slab": 2.0, "tau_c": 1.0,
               "epsilon": 0.1, "L_g": 2.0, "min_parents": 4, "beta_check": 0.6, "cal_win": [12.0, 28.0], "n_cal": 8,
               "n_seeds": 64, "n_checked": 64, "n_se": 2.0, "boot_B": 4000, "seps": [8.0, 12.0, 16.0, 20.0],
               "w_lo": 1.25, "w_hi": 2.5, "h_fd": 1.0},
    "parts": {"force": {"mode": "force"}, "law": {"mode": "law"}},
    "quick": {"n_seeds": 3, "n_cal": 1, "n_checked": 1, "boot_B": 200, "T": 28.0, "T_cal": 24.0, "margin": 6.0,
              "r_src": 2.0, "t_src": 2.0, "cal_win": [6.0, 14.0], "seps": [6.0, 7.0, 8.0, 9.0]},
    "estimate": "force ~20–30 min, law ~20–30 min (10 processes; both parts compute the same, the verdicts are "
                "separate; the time goes into the sprinkling and the links; passes, checks and the prediction take a "
                "fraction of a second)",
}


# ============================================================================================== prediction

def field_at(points, P, F, tau_c):
    """Discrete retarded field sum F / (2 pi tau) of the events P (flux F) at the points (causal past, tau >= tau_c)."""
    out = np.zeros(len(points))
    for k0 in range(0, len(points), 256):
        Q = points[k0:k0 + 256]
        d = Q[:, None, :] - P[None, :, :]
        t2 = d[..., 0] ** 2 - (d[..., 1:] ** 2).sum(axis=-1)
        ok = (d[..., 0] > 0.0) & (t2 >= tau_c * tau_c)
        out[k0:k0 + 256] = (np.where(ok, F[None, :] / (2.0 * math.pi * np.sqrt(np.where(ok, t2, 1.0))), 0.0)).sum(axis=1)
    return out


def predicted_kicks(X, run_self, run_partner, lam, p):
    """
    First-order prediction (module notes): x-momentum per unit flux carried to every event of the bundle alone,
    from the kicks lambda d_x phi_partner(v) tau_link along its actual flow. Returns the carried kick K per event.
    """
    fl, dec = run_self["flux"][0], run_self["dec"][0]
    car = np.nonzero(fl > 0.0)[0]
    pf = run_partner["flux"][0]
    pe = np.nonzero(pf > 0.0)[0]
    h = float(p["h_fd"])
    ex = np.zeros(X.shape[1])
    ex[1] = h
    grad = (field_at(X[car] + ex, X[pe], pf[pe], float(p["tau_c"]))
            - field_at(X[car] - ex, X[pe], pf[pe], float(p["tau_c"]))) / (2.0 * h)
    g = dict(zip(car.tolist(), grad.tolist()))
    K = np.zeros(len(X))              # carried kick (x-momentum per unit flux) arriving at the event
    inK = np.zeros(len(X))            # flux-weighted sum of the incoming carried kicks
    for v in car:                      # increasing index = a linear extension of the causal order
        K[v] = inK[v] / fl[v]          # the source part carries no kick; fl[v] includes it
        k = int(dec[v])
        if k < 0:
            continue
        b = X[k] - X[v]
        tau = math.sqrt(max(float(mink(b, b)), 0.0))
        inK[k] += fl[v] * (K[v] + lam * g[int(v)] * tau)
    return K


def series_from(X, fl, val, p):
    """Per slab: sum over the carriers of F * val."""
    T, sl = float(p["T"]), float(p["slab"])
    edges = np.arange(0.0, T - sl + 1e-9, sl)
    t, s = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        m = (fl > 0.0) & (X[:, 0] >= a) & (X[:, 0] < b)
        t.append(0.5 * (a + b))
        s.append(float((fl[m] * val[m]).sum()))
    return np.array(t), np.array(s)


def line_slope(d, p):
    """Continuum static line (unit strength), the same window and fit: slope of the integrated gradient G_ret."""
    tc = float(p["tau_c"])
    ts = np.arange(0.0, float(p["T"]), 0.01)
    v = np.cumsum([G_ret(t, d, tc) for t in ts]) * 0.01
    lo, hi = window(d, p)
    return lin_slope(ts, v, lo, hi)


def window(d, p):
    return float(p["w_lo"]) * d, min(float(p["w_hi"]) * d, float(p["T"]) - 4.0)


# ============================================================================================== task

def force_task(p, seed, d, lam, checks):
    p = dict(p, T=float(p["w_hi"]) * d + 6.0)    # the box ends just after the window (the dynamics is causal)
    X, ptr, idx, cptr, cidx = make(p, seed, d)
    u0 = np.array([1.0, 0.0, 0.0])
    sL, sR = tube(X, -0.5 * d, p, 1.0), tube(X, 0.5 * d, p, 1.0)
    aL = transport_species(X, ptr, idx, cptr, cidx, [sL], [[0]], u0, lam, p)
    aR = transport_species(X, ptr, idx, cptr, cidx, [sR], [[0]], u0, lam, p)
    bb = transport_species(X, ptr, idx, cptr, cidx, [sL, sR], [[0, 1], [0, 1]], u0, lam, p)
    lo, hi = window(d, p)
    out = {}
    # measured: momentum received (B - A), towards the partner positive
    sAL, sAR = species_series(X, aL["flux"][0], aL["dir"][0], p), species_series(X, aR["flux"][0], aR["dir"][0], p)
    sBL, sBR = species_series(X, bb["flux"][0], bb["dir"][0], p), species_series(X, bb["flux"][1], bb["dir"][1], p)
    t = np.asarray(sAL["t"])
    dPL = np.asarray(sBL["P"]) - np.asarray(sAL["P"])
    dPR = np.asarray(sBR["P"]) - np.asarray(sAR["P"])
    out["meas"] = 0.5 * (lin_slope(t, dPL, lo, hi) - lin_slope(t, dPR, lo, hi))
    # linear-response validity: change of the separation over the window (B - A)
    dr = (np.asarray(sBR["x"]) - np.asarray(sBL["x"])) - (np.asarray(sAR["x"]) - np.asarray(sAL["x"]))
    w = (t >= lo) & (t <= hi)
    out["dr"] = float(dr[w][-1] - dr[w][0])
    # predicted from branch A only
    KL = predicted_kicks(X, aL, aR, lam, p)
    KR = predicted_kicks(X, aR, aL, lam, p)
    tL, PL = series_from(X, aL["flux"][0], KL, p)
    tR, PR = series_from(X, aR["flux"][0], KR, p)
    out["pred"] = 0.5 * (lin_slope(tL, PL, lo, hi) - lin_slope(tR, PR, lo, hi))
    if checks:
        both = [(0, 0), (1, 1)]
        cc = transport_species(X, ptr, idx, cptr, cidx, [sL, sR], [[0], [1]], u0, lam, p)
        cA = {k: [aL[k][0], aR[k][0]] for k in ("dec", "gap", "flux", "dir")}
        out["C"] = compare(cA, cc, both, lambda D_: D_, True)
        bm = transport_species(Xm_of(X), ptr, idx, cptr, cidx, [sL, sR], [[0, 1], [0, 1]], u0, lam, p)
        out["mirror"] = compare(bb, bm, both, lambda D_: D_ * np.array([1.0, -1.0, 1.0]), False)
        bc = float(p["beta_check"])
        ub = boost_points(u0[None, :], bc)[0]
        bf = transport_species(boost_points(X, bc), ptr, idx, cptr, cidx, [sL, sR], [[0, 1], [0, 1]], ub, lam, p)
        out["frame"] = compare(bb, bf, both, lambda D_: boost_points(D_, -bc), False)
    return out


# ============================================================================================== run

def kfit(ds, v):
    v = np.asarray(v, float)
    if np.any(v <= 0) or not np.all(np.isfinite(v)):
        return float("nan")
    return float(-np.polyfit(np.log(ds), np.log(v), 1)[0])


def run(p, ctx):
    mode = str(p["mode"])
    t0 = time.time()
    n, B, nse = int(p["n_seeds"]), int(p["boot_B"]), float(p["n_se"])
    pc = dict(p, T=float(p["T_cal"]))           # the F3.10 calibration: same streams, same box
    cal = parallel_map(calibration_task, [(pc, derive_seed(ctx.seed, "F3.10-cal", k)) for k in range(int(p["n_cal"]))],
                       ctx.workers, ctx.progress, lambda i, r_, dd, nn: f"[F3.10b] calibration {dd}/{nn}")
    phi_ref = float(np.median(np.concatenate([c[0] for c in cal])))
    lam = float(p["epsilon"]) / phi_ref
    ds = np.array([float(x) for x in p["seps"]])
    jobs, keys = [], []
    for d in ds:
        for k in range(n):
            jobs.append((p, derive_seed(ctx.seed, "F3.10b", float(d), k), float(d), lam, k < int(p["n_checked"])))
            keys.append((float(d), k))
    res = parallel_map(force_task, jobs, ctx.workers, ctx.progress,
                       lambda i, r_, dd, nn: f"[F3.10b] {mode} d = {keys[i][0]:g} sprinkling {keys[i][1]+1} ({dd}/{nn}, "
                                             f"{time.time()-t0:.0f}s)")
    chk = [r_ for r_ in res if "C" in r_]
    inst = {"C equals A bit for bit": all(r_["C"] == "identical" for r_ in chk),
            "mirror x -> -x": all(r_["mirror"] in ("identical", "tie") for r_ in chk),
            "frame beta' = 0.6": all(r_["frame"] in ("identical", "tie") for r_ in chk)}
    n_tie = sum(r_["mirror"] == "tie" for r_ in chk) + sum(r_["frame"] == "tie" for r_ in chk)
    meas = {d: np.array([r_["meas"] for r_, kk in zip(res, keys) if kk[0] == d]) for d in ds}
    pred = {d: np.array([r_["pred"] for r_, kk in zip(res, keys) if kk[0] == d]) for d in ds}
    drs = {d: np.array([r_["dr"] for r_, kk in zip(res, keys) if kk[0] == d]) for d in ds}
    rng = np.random.default_rng(derive_seed(ctx.seed, "F3.10b-boot", mode))
    boots = []                                    # paired bootstrap: the same sprinklings for measurement and prediction
    for _ in range(B):
        bm, bp = [], []
        for d in ds:
            j = rng.integers(0, n, n)
            bm.append(meas[d][j].mean())
            bp.append(pred[d][j].mean())
        boots.append((np.array(bm), np.array(bp)))
    mm = np.array([meas[d].mean() for d in ds])
    pm = np.array([pred[d].mean() for d in ds])
    se_m = np.array([np.std([b_[0][i] for b_ in boots], ddof=1) for i in range(ds.size)])
    rho = mm / pm
    se_rho = np.array([np.std([b_[0][i] / b_[1][i] for b_ in boots], ddof=1) for i in range(ds.size)])
    resolved = mm - nse * se_m > 0
    linear = np.array([abs(drs[d].mean()) <= d / 4.0 for d in ds])
    k_line = kfit(ds, [line_slope(d, p) for d in ds])
    lines = [f"2+1, rho = 1; F3.10 model; lambda = {lam:.5f} (F3.10 calibration, phi_ref = {phi_ref:.3f}); T = 2.5 d + 6; "
             f"window [{p['w_lo']} d, {p['w_hi']} d]; {n} sprinklings per d; received momentum rate of the pair "
             f"(B - A, towards each other positive) against the first-order prediction from branch A"]
    for i, d in enumerate(ds):
        lines.append(f"  d = {d:g}: measured {mm[i]:+.4f} +- {se_m[i]:.4f}, predicted {pm[i]:+.4f}, ratio rho = {rho[i]:.3f} +- "
                     f"{se_rho[i]:.3f}; change of the separation over the window {drs[d].mean():+.2f} (limit d/4 = {d/4:.1f}); "
                     f"resolved {bool(resolved[i])}, linear response {bool(linear[i])}")
    status = "UNDECIDED"
    why = []
    if not resolved.all():
        why.append("not resolved at d = " + ", ".join(f"{d:g}" for d, ok in zip(ds, resolved) if not ok))
    if not linear.all():
        why.append("outside linear response at d = " + ", ".join(f"{d:g}" for d, ok in zip(ds, linear) if not ok))
    if mode == "force":
        dev_hi = np.abs(rho - 1.0) + nse * se_rho
        dev_lo = np.abs(rho - 1.0) - nse * se_rho
        if np.any((dev_lo > 0.25) & linear):
            status = "FAILED"
        elif not why and np.all(dev_hi <= 0.25):
            status = "PASSED"
        elif not why:
            why.append("standard error too large or rho near the tolerance at d = "
                       + ", ".join(f"{d:g}" for d, x in zip(ds, dev_hi) if x > 0.25))
        lines.append("force: |rho - 1| + 2 SE = " + ", ".join(f"{x:.3f}" for x in dev_hi) + "; |rho - 1| - 2 SE = "
                     + ", ".join(f"{x:.3f}" for x in dev_lo) + " (tolerance 1/4)")
    else:
        k_m, k_p = kfit(ds, mm), kfit(ds, pm)
        dk = np.array([kfit(ds, b_[0]) - kfit(ds, b_[1]) for b_ in boots])
        ok = np.isfinite(dk)
        dk_se = float(np.std(dk[ok], ddof=1)) if ok.sum() > 2 else float("nan")
        dkm = k_m - k_p
        lines.append(f"law: k_meas = {k_m:.3f}, k_pred = {k_p:.3f} (discrete field of the actual partner), Delta k = {dkm:+.3f} "
                     f"+- {dk_se:.3f} (replicates without an exponent: {int((~ok).sum())} of {dk.size}); continuum static "
                     f"line in the same windows k_line = {k_line:.3f}; kernel 1")
        if np.isfinite(dkm) and np.isfinite(dk_se):
            if abs(dkm) - nse * dk_se > 0.25 and linear.all():
                status = "FAILED"
            elif not why and abs(dkm) + nse * dk_se <= 0.25:
                status = "PASSED"
            elif not why:
                why.append(f"2 SE of Delta k = {nse * dk_se:.3f}" + (" alone exceeds 1/4" if nse * dk_se > 0.25 else ""))
        else:
            why.append("the exponent is not defined (a non-positive mean)")
    if status == "UNDECIDED":
        lines.append("undecided because: " + "; ".join(why))
    if not all(inst.values()):
        status = "INSTRUMENT_FAILED"
    lines.append(f"instrument ({len(chk)} sprinklings): "
                 + "; ".join(f"{k}: {'OK' if v else 'FAILS'}" for k, v in inst.items()) + f"; numerical ties {n_tie}")
    data = {"lambda": lam, "phi_ref": phi_ref, "seps": ds.tolist(), "meas": {str(d): meas[d].tolist() for d in ds},
            "pred": {str(d): pred[d].tolist() for d in ds}, "dr": {str(d): drs[d].tolist() for d in ds},
            "rho": rho.tolist(), "se_rho": se_rho.tolist(), "k_line": k_line, "instrument": inst}
    return Result(status, lines, data)
