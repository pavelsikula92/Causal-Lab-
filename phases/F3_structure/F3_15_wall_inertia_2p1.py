"""
F3.15 -- Inertia of a domain wall of a field of the causal order in 2+1, intrinsic instrument (direction B, B1)
=============================================================================================================

Why this test exists (12 Nov 2026; card docs/B1_card_draft.md, revision 2, frozen by decision 21)
-------------------------------------------------------------------------------------------------
F3.14b (1+1, replicated) found that a sine-Gordon kink evolved by the retarded Green function of the order keeps its
velocity and contraction and has an inertial mass equal to its rest energy, while a relaxational dynamics and a
causality-blind sum do not. B1 asks whether this survives (i) the transverse dimension and (ii) the intrinsic weight
of 2+1: the Green function read from the order alone, w(L) of decision 17 (variant B, constants frozen below) with L
the longest chain between the two events (B0 closed, decision 20). The wall is a test body (decision 18); decision 16
(D4: the wall first). Development rows: docs/findings.md, Undecided, "B1 development rows" (12 Nov).

Model (frozen)
--------------
2+1 Poisson sprinkling of density rho in the pyramid of dependence {0 <= t <= T, |x| <= X - t, |y| <= Y - t}, X = T + 8,
Y = T + 2.5 (no wall enters the dynamics).
    phi(v) = phi_0(v) - (1/rho) sum over u < v (the order) of w(L_uv) [mu sin phi(u) + f chi(t_u)],
w(1) = a rho^(1/3), w(L) = c rho^(1/3) / (L + beta) for L >= 2 (a, c, beta frozen from diag/out/b0_weight_calib_
20261111T174914.json), L_uv the longest chain over the links (forward dynamic programme; links by diag.b0_chain_alg.links3);
mu = 1; phi_0 the free massless wave of the initial data of the sine-Gordon kink extended along y and launched at v
(d'Alembert in x; the launch is the only place where coordinates enter the dynamics).
Arms on the same sprinkling: rest, v = 0.3, v = 0.6, kick (from rest, f = 0.05 for 1 <= t < 3).
Controls (first n_ctrl sprinklings):
  telegraph null phi_tt + Gamma phi_t - box_x,y phi + sin phi = -f chi, Gamma = 1, 2, launched v = 0.3 and kick: evolved
     for psi = exp(Gamma t/2) phi with the same weights, the mass term -Gamma^2/4 moved into the source,
     J_psi = exp(Gamma t/2)(sin phi + f chi) - Gamma^2/4 psi; coordinate time of a medium frame (declared; a null model).
     Gamma = 1, 2 instead of the 2, 5 of F3.14b: in the smoke run (rho 10) the Gamma = 5 wall dissolved (Q 0.83 -> 0.3,
     width 1.4 -> 8) instead of stopping, most likely because the source growing as exp(Gamma t/2) weights the region
     near the cone, where the chain weight is pointwise least accurate (decision 17); at Gamma 2 the wall stayed intact;
  V = 0 (launched v = 0.3): the sum vanishes, phi = phi_0 (readout check);
  causality-blind (launched v = 0.3): each event v gives the same multiset of contributions w(L_vu) J(v) / rho as in the
     conservative rule, but to as many events u' drawn uniformly with replacement from all later events (index > v in the
     time order; the same time range); same normalisation, phi_0 and parameters. (The per-target form of F3.14b needs a
     unit weight; in 2+1 the weight depends on the pair, so the weights are kept and only the index set is randomised.)
Waves (n_wave own sprinklings, B0 readouts of diag.b0_chain_field): vacuum (8 random modes, total amplitude 0.01, box
T + 4) and plane waves k = (1, 0), (0, 1), (1, 1), (2, 1), amplitude 0.01.
Second density for I3: the rest arm at rho_i3 on n_i3 own sprinklings.

Readouts (frozen; coordinates, declared). Band |y| < 2. Slabs of width 1, centres 0.5 ... T - 1.5: the slab rule of the
card, the last slab (centre T - 0.5) is not read in any arm (at v = 0.6 the plateau window leaves the pyramid).
Per slab: kink fit (field1d.fit_kink, guess = the previous centre or v t_m; width, centre); Q from the plateaus
x_g +- (5 ... 7), x_g = v t_m; E = int [phi_t^2/2 + phi_x^2/2 + mu (1 - cos phi)] dx over |x - x_c| < 5 from local linear
regressions in cells 0.5 x 0.5, step 0.25, with the noise bias of the squared estimates removed (E[c^2] = c^2 + Var c,
E[cos c0] = cos c0 exp(-Var c0 / 2)); E_loc = the same over |x - x_g| < 3 with the model's own mu (independent of the
fit). The exact readouts Q_ex, w_ex, E_ex: the same readouts on the exact continuum solution at the same events.
Velocities: slopes of x_c over the read slabs (all, or t in [4, T - 1) for the kick and the telegraph); a slab whose fit
fails is skipped; fewer than n - 1 read slabs of a window gives nan. p = gamma v. A quantity that cannot be read is nan
and its condition is not met.

Criteria (card revision 2; SE = sd / sqrt(n) over sprinklings; "+ 2 SE" on the PASSED side, "- 2 SE" on the FAILED side)
INSTRUMENT_FAILED if any of:
  I1 vacuum: |mean lambda_eff| + 2 SE >= 0.02;  I2 each k: |omega/omega_c - 1| + 2 SE >= 0.015;
  I3 rest wall at rho_i3: |w/w_ex - 1| + 2 SE >= 0.05;  I4 forward/backward kernel, linear extension t + 0.3 x and the
  pyramid X against X - 2 on the common events: max relative difference >= 1e-12;
  rest (A): |Q - Q_ex| + 2 SE >= 0.01, |drift x_c(last read) - x_c(first)| + 2 SE >= 0.15, |w/w_ex - 1| + 2 SE >= 0.05,
  |E - E_ex| / 8 + 2 SE >= 0.05;
  telegraph, each Gamma and each of (v 0.3, kick): not [(p_cons - p_Gamma) - 2 SE > 0 and mean p_Gamma <= 0.5 mean p_cons],
  p over t in [4, T - 1), paired;
  V = 0: D_L = L_cons - L_V0 (L = E_loc(last read) / E_loc(first), paired with v 0.3) with D_L - 2 SE <= 0.2.
FAILED if any of (v = 0.3 and 0.6 separately): |v_meas - v| - 2 SE > 0.03; |C - 1| - 2 SE > 0.05 (C = w(v) gamma /
  w(rest), paired); |w/w_ex - 1| - 2 SE > 0.05; |e_B| - 2 SE > 0.05 (e_B = (E - E_ex) / (8 gamma)); |r_M - 1| - 2 SE > 0.10
  (r_M = 2 pi Q_kick f tau / mean Delta(gamma v) / M_rest, M_rest = mean E of the rest arm, SE by bootstrap over pairs).
PASSED if all: each v |v_meas - v| + 2 SE < 0.03, |C - 1| + 2 SE < 0.05, |w/w_ex - 1| + 2 SE < 0.05, |e_B| + 2 SE < 0.05,
  |Q - Q_ex| + 2 SE < 0.01; |r_M - 1| + 2 SE < 0.10; the blind control separated: D_L(blind) - 2 SE > 0.2.
Otherwise UNDECIDED. Every condition on its own; nothing averaged across conditions; nothing chosen after the run.

    python run.py F3.15 --workers N          (seed 0 decides; --seed 1 replicates)
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import time

import numpy as np

from core.order.sprinkle import derive_seed
from core.stats.parallel import parallel_map
from core.testing import Result
from phases.F3_structure import field1d as F
from diag.b0_green2p1 import sprinkle_pyramid
from diag.b0_chain_alg import links3
from diag import b0_chain_field as B0

TEST = {
    "id": "F3.15",
    "phase": "F3",
    "title": "Inertia of a domain wall of a field of the causal order in 2+1 (intrinsic Green weight, B1)",
    "theory_part": "F3 emergence of structure; direction B (2+1): the Green function from the order alone (longest chain, "
                   "decision 17, variant B); the wall as a test body (decisions 16, 18, 20)",
    "question": "Does a domain wall of a field on all events of a 2+1 causal order, evolved by the retarded Green function "
                "read from the order alone (the longest-chain weight), keep its identity at rest, keep its velocity and "
                "Lorentz contraction when launched, and answer a finite force with an inertial mass equal to its rest "
                "energy per unit length, while a relaxational dynamics and a causality-blind sum do not?",
    "method": "2+1 sprinkling in the pyramid of dependence; phi = phi_0 - (1/rho) sum over the causal past of w(L) [sin phi "
              "+ f chi]; sine-Gordon wall launched at v = 0, 0.3, 0.6 and kicked (f 0.05, tau 2); readouts on |y| < 2 by "
              "kink fits and noise-debiased local regressions per slab; I1-I4, rest, telegraph null (Gamma 1, 2), V = 0, "
              "causality-blind; deciding density 30, 32 sprinklings",
    "prediction": "continuum sine-Gordon: velocity kept, w gamma = 1, E = 8 gamma per unit length, M_in = M_rest = 8; the "
                  "telegraph wall stops; the blind sum loses localisation",
    "criterion": "INSTRUMENT: I1 |lambda| + 2 SE < 0.02; I2 |omega/omega_c - 1| + 2 SE < 0.015 each k; I3 |w/w_ex - 1| + 2 SE "
                 "< 0.05 at rho 60; I4 < 1e-12; rest Q 0.01, drift 0.15, w 0.05, E/8 0.05 (+ 2 SE); telegraph each Gamma "
                 "and arm separated and p_Gamma <= 0.5 p_cons; V0 D_L - 2 SE > 0.2. PASSED: each v |dv| + 2 SE < 0.03, "
                 "|C - 1| + 2 SE < 0.05, |w/w_ex - 1| + 2 SE < 0.05, |e_B| + 2 SE < 0.05, Q 0.01; |r_M - 1| + 2 SE < 0.10; "
                 "blind D_L - 2 SE > 0.2. FAILED: the same with - 2 SE above the tolerance (v, C, w, e_B, r_M). Otherwise "
                 "UNDECIDED (card docs/B1_card_draft.md rev. 2)",
    "null_model": "telegraph (overdamped) dynamics with the same potential, initial data, kick and weights; V = 0 (readout "
                  "check); causality-blind sum with the same weights and number of terms; exact continuum solutions",
    "level": "mechanism",
    "meters": [],
    "depends_on": [],
    "params": {"rho_dec": 30.0, "rho_i3": 60.0, "n_dec": 32, "n_ctrl": 8, "n_wave": 4, "n_i3": 8, "T": 8.0,
               "dX": 8.0, "dY": 2.5, "band": 2.0, "vs": [0.3, 0.6], "f": 0.05, "t1": 1.0, "t2": 3.0, "gammas": [1.0, 2.0],
               "calib": [0.2732723632276477, 0.21313211596642687, 1.2025922117299621],
               "waves": ["vacuum", "k1_0", "k0_1", "k1_1", "k2_1"],
               "i4_n": 2, "i4_rho": 30.0, "i4_T": 3.0, "i4_tol": 1e-12,
               "tol_I1": 0.02, "tol_I2": 0.015, "tol_I3": 0.05, "tol_Q": 0.01, "tol_drift": 0.15, "tol_w": 0.05,
               "tol_E": 0.05, "tol_v": 0.03, "tol_C": 0.05, "tol_M": 0.10, "drop": 0.5, "tol_DL": 0.2,
               "t_kick_read": 4.0, "n_se": 2.0, "boot_B": 2000},
    "quick": {"rho_dec": 10.0, "rho_i3": 15.0, "n_dec": 2, "n_ctrl": 2, "n_wave": 2, "n_i3": 2, "i4_n": 1,
              "i4_rho": 15.0, "boot_B": 200},
    "estimate": "long: ~80 000 events per wall arm at rho 30, ~7.5 min per arm on one core (O(N E_H)); 32 x 4 + 8 x 6 wall "
                "tasks, 20 wave tasks, 8 tasks at rho 60 (~30 min each): ~30 core-hours. numba required. Finished tasks "
                "are cached in diag/out/F3_15_cache (keyed by this file's hash), so an interrupted run resumes.",
}

_nb = F._nb

if _nb is not None:
    @_nb.njit(cache=True)
    def _evolve(ptr, idx, t, phi0, rho, mu, wL, f, t1, t2, gam, blind, seed):
        """The field rule with the chain weight. gam > 0: telegraph in psi = exp(gam t/2) phi (phi0 is psi_0); blind:
        each contribution of v goes to a uniformly drawn later event instead of u. Returns phi."""
        N = phi0.size
        if blind:
            np.random.seed(seed)
        acc = np.zeros(N)
        phi = np.empty(N)
        L = np.zeros(N, np.int64)
        stamp = np.full(N, -1, np.int64)
        for v in range(N):
            ps = phi0[v] - acc[v]
            if gam > 0.0:
                e = math.exp(0.5 * gam * t[v])
                ph = ps / e
            else:
                e = 1.0
                ph = ps
            phi[v] = ph
            J = mu * math.sin(ph)
            if t[v] >= t1 and t[v] < t2:
                J += f
            if gam > 0.0:
                J = e * J - 0.25 * gam * gam * ps
            coef = J / rho
            L[v] = 0
            stamp[v] = v
            for u in range(v + 1, N):
                b = -1
                for q in range(ptr[u], ptr[u + 1]):
                    p = idx[q]
                    if p >= v and stamp[p] == v and L[p] + 1 > b:
                        b = L[p] + 1
                if b > 0:
                    L[u] = b
                    stamp[u] = v
                    if blind:
                        acc[np.random.randint(v + 1, N)] += wL[b] * coef
                    else:
                        acc[u] += wL[b] * coef
        return phi


def energy_db(t, x, phi, tm, xc, mu=1.0, half=5.0, cell=0.5, step=0.25):
    """field1d.energy with the noise bias of the squared estimates removed (Var from the residuals of each regression)."""
    s = np.abs(t - tm) < 0.5 * cell
    ts, xs, ps = t[s], x[s], phi[s]
    o = np.argsort(xs)
    ts, xs, ps = ts[o], xs[o], ps[o]
    E, used = 0.0, 0
    for xg in np.arange(xc - half, xc + half + 1e-9, step):
        i0, i1 = np.searchsorted(xs, [xg - 0.5 * cell, xg + 0.5 * cell])
        n = i1 - i0
        if n < 8:
            continue
        A = np.column_stack([np.ones(n), ts[i0:i1] - tm, xs[i0:i1] - xg])
        c, *_ = np.linalg.lstsq(A, ps[i0:i1], rcond=None)
        r = ps[i0:i1] - A @ c
        s2 = float(r @ r) / (n - 3)
        C = s2 * np.linalg.inv(A.T @ A)
        E += (0.5 * (c[1] ** 2 - C[1, 1]) + 0.5 * (c[2] ** 2 - C[2, 2])
              + mu * (1.0 - math.cos(c[0]) * math.exp(0.5 * C[0, 0]))) * step
        used += 1
    n_all = int(round(2 * half / step)) + 1
    return E * n_all / used if used else float("nan")


def _weights(rho, calib):
    return B0.weight_table("B", rho, tuple(calib))


def _track(t, x, phi, v, T, mu, exact):
    rows = []
    for a in np.arange(0.0, T - 2.0 + 1e-9, 1.0):          # slab rule: the last slab [T - 1, T) is not read
        tm = a + 0.5
        xg = v * tm
        prev = rows[-1].get("xc") if rows else None
        guess = prev if (prev is not None and np.isfinite(prev)) else xg
        s = (t >= a) & (t < a + 1.0)
        hi = s & (x > xg + 5) & (x < xg + 7)
        lo = s & (x < xg - 5) & (x > xg - 7)
        row = {"t": float(tm), "E_loc": energy_db(t, x, phi, tm, xg, mu, half=3.0),
               "Q": float((np.mean(phi[hi]) - np.mean(phi[lo])) / (2 * math.pi)) if hi.any() and lo.any() else float("nan")}
        r = F.fit_kink(t[s], x[s], phi[s], a, a + 1.0, guess, 1.0)
        row["xc"] = r["xc"] if r else float("nan")
        row["w"] = r["w"] if r else float("nan")
        xe = row["xc"] if np.isfinite(row["xc"]) else xg
        row["E"] = energy_db(t, x, phi, tm, xe, mu)
        row["E_plain"] = F.energy(t, x, phi, tm, xe, mu)
        if exact is not None:
            rx = F.fit_kink(t[s], x[s], exact[s], a, a + 1.0, xg, 1.0)
            row["w_ex"] = rx["w"] if rx else float("nan")
            row["Q_ex"] = float((np.mean(exact[hi]) - np.mean(exact[lo])) / (2 * math.pi)) if hi.any() and lo.any() \
                else float("nan")
            row["E_ex"] = energy_db(t, x, exact, tm, xe, mu)
        rows.append(row)
    return rows


def task_wall(p, rho, k, arm, base):
    """arm: rest, v0.3, v0.6, kick; tele<G>_v0.3, tele<G>_kick; V0_v0.3; blind_v0.3."""
    T = float(p["T"])
    t, x, y = sprinkle_pyramid(rho, T + float(p["dX"]), T + float(p["dY"]), T, derive_seed(base, "F3.15", rho, k))
    variant, motion = arm.split("_") if "_" in arm else ("cons", arm)
    v = 0.0 if motion in ("rest", "kick") else float(motion[1:])
    f = float(p["f"]) if motion == "kick" else 0.0
    mu = 0.0 if variant == "V0" else 1.0
    t0 = time.time()
    if variant == "V0":
        phi = F.kink_free(t, x, v)
    else:
        ptr, idx = links3(t, x, y, (12.0 * 25.0 / (math.pi * rho)) ** (1.0 / 3.0))
        gam = float(variant[4:]) if variant.startswith("tele") else 0.0
        phi0 = F.kink_free_telegraph(t, x, v, gam) if gam > 0 else F.kink_free(t, x, v)
        sd = int(derive_seed(base, "F3.15-blind", rho, k) % (2 ** 31))
        phi = _evolve(ptr, idx, t, phi0, float(rho), mu, _weights(rho, p["calib"]), f, float(p["t1"]), float(p["t2"]),
                      gam, variant == "blind", sd)
    band = np.abs(y) < float(p["band"])
    tb, xb = t[band], x[band]
    exact = F.kink_exact(tb, xb, v) if (variant == "cons" and motion != "kick") else None
    rows = _track(tb, xb, phi[band], v, T, mu, exact)
    return {"kind": "wall", "rho": rho, "k": k, "arm": arm, "variant": variant, "motion": motion, "v": v,
            "N": int(t.size), "sec": time.time() - t0, "rows": rows}


def task_wave(p, rho, k, what, base):
    T = float(p["T"])
    X, Y = B0.box(what, T)
    sd = derive_seed(base, "F3.15-wave", rho, k)
    t, x, y = sprinkle_pyramid(rho, X, Y, T, sd)
    phi0, info = B0.initial(what, t, x, y, sd)
    t0 = time.time()
    ptr, idx = links3(t, x, y, (12.0 * 25.0 / (math.pi * rho)) ** (1.0 / 3.0))
    phi = _evolve(ptr, idx, t, phi0, float(rho), 1.0, _weights(rho, p["calib"]), 0.0, 0.0, 0.0, 0.0, False, 0)
    out = {"kind": "wave", "rho": rho, "k": k, "what": what, "N": int(t.size), "sec": time.time() - t0}
    out.update(B0.readout(what, t, x, y, phi, T, info))
    return out


def task_i4(p, k, base):
    """(1) forward kernel against the backward DP of diag.b0_weight_valid; (2) linear extension t + 0.3 x; (3) the pyramid
    X against X - 2 on the common events. Kicked wall at v 0.3 (f 0.05, t in [1, 3))."""
    from diag.b0_weight_valid import children, target_pass
    rho, T = float(p["i4_rho"]), float(p["i4_T"])
    calib = tuple(p["calib"])
    wL = _weights(rho, calib)
    t, x, y = sprinkle_pyramid(rho, T + 4.0, T + 3.0, T, derive_seed(base, "F3.15-i4", k))
    tc = (12.0 * 25.0 / (math.pi * rho)) ** (1.0 / 3.0)
    ptr, idx = links3(t, x, y, tc)
    # (1) linear Green application (J given) against the backward DP sum at 30 targets
    J = 1.0 / np.cosh((0.8 * x + 0.6 * y) / 0.7) ** 2
    lin, _, _ = B0.evolve_chain(ptr, idx, np.zeros(t.size), rho, 1.0, wL, 0.0, J, False)
    cptr, cidx = children(ptr, idx, t.size)
    rng = np.random.default_rng(derive_seed(base, "F3.15-i4-targets", k))
    tg = rng.choice(np.nonzero(t > T - 1)[0], 30, replace=False)
    sw = np.zeros((3, 120)); sg = np.zeros((3, 120)); cn = np.zeros((3, 120))
    d1 = 0.0
    for xi in tg:
        o = target_pass(int(xi), t, x, y, cptr, cidx, J[None, :], rho, *calib, 0.68 * rho ** (-1 / 3), sw, sg, cn)
        d1 = max(d1, abs(-lin[xi] - o[0, 2]) / max(abs(o[0, 2]), 1e-30))
    # (2) linear extension t + 0.3 x
    args = (float(p["f"]), float(p["t1"]), float(p["t2"]), 0.0, False, 0)
    phi0 = F.kink_free(t, x, 0.3)
    ref = _evolve(ptr, idx, t, phi0, rho, 1.0, wL, *args)
    sc = float(np.max(np.abs(ref)))
    o = np.argsort(t + 0.3 * x, kind="stable")
    inv = np.empty_like(o); inv[o] = np.arange(o.size)
    rows = [np.sort(inv[idx[ptr[j]:ptr[j + 1]]]) for j in o]
    ptr2 = np.zeros(t.size + 1, np.int64); ptr2[1:] = np.cumsum([r.size for r in rows])
    idx2 = np.concatenate(rows).astype(np.int64)
    alt = _evolve(ptr2, idx2, t[o], phi0[o], rho, 1.0, wL, *args)
    d2 = float(np.max(np.abs(alt - ref[o])) / sc)
    # (3) the smaller pyramid X - 2 (same events inside it; its own links)
    sub = np.nonzero(np.abs(x) <= (T + 2.0) - t)[0]
    ps, ix = links3(t[sub], x[sub], y[sub], tc)
    phs = _evolve(ps, ix, t[sub], phi0[sub], rho, 1.0, wL, *args)
    d3 = float(np.max(np.abs(phs - ref[sub])) / sc)
    return {"kind": "i4", "k": k, "N": int(t.size), "kernel": float(d1), "order": d2, "region": d3}


_SHA = None


def _sha():
    global _SHA
    if _SHA is None:
        with open(__file__, "rb") as fh:
            _SHA = hashlib.sha256(fh.read()).hexdigest()[:16]
    return _SHA


def run_task(kind, *a):
    """Runs a task; finished full-run tasks are cached on disk (keyed by this file's hash), so a run can resume."""
    p = a[0]
    key = json.dumps([kind] + [q for q in a[1:]], default=str)
    fn = os.path.join("diag", "out", "F3_15_cache", _sha(), hashlib.sha256(
        (key + json.dumps(p, sort_keys=True, default=str)).encode()).hexdigest()[:24] + ".json")
    if os.path.exists(fn):
        with open(fn) as fh:
            return json.load(fh)
    out = {"wall": task_wall, "wave": task_wave, "i4": task_i4}[kind](*a)
    os.makedirs(os.path.dirname(fn), exist_ok=True)
    with open(fn + ".tmp", "w") as fh:
        json.dump(out, fh)
    os.replace(fn + ".tmp", fn)
    return out


# ============================================================================================== analysis

def _col(rows, key):
    return np.array([r.get(key, float("nan")) for r in rows], float)


def _slope(rows, tmin=0.0):
    tm = _col(rows, "t")
    xc = _col(rows, "xc")
    m = tm >= tmin
    ok = m & np.isfinite(xc)
    if ok.sum() < max(2, m.sum() - 1):
        return float("nan")
    return float(np.polyfit(tm[ok], xc[ok], 1)[0])


def _gv(v):
    return v / math.sqrt(1 - v * v) if abs(v) < 1 else float("nan")


def _ms(a):
    a = np.asarray(a, float)
    if a.size < 2 or not np.all(np.isfinite(a)):
        return float("nan"), float("nan")
    return float(a.mean()), float(a.std(ddof=1) / math.sqrt(a.size))


def _fmt(m):
    return f"{m[0]:+.4f} +- {m[1]:.4f}"


def run(p, ctx):
    if _nb is None:
        raise SystemExit("F3.15 needs numba (pip install numba)")
    t0 = time.time()
    rd, ri3 = float(p["rho_dec"]), float(p["rho_i3"])
    n_dec, n_ctrl, n_wave, n_i3 = int(p["n_dec"]), int(p["n_ctrl"]), int(p["n_wave"]), int(p["n_i3"])
    nse, T = float(p["n_se"]), float(p["T"])
    v0 = float(p["vs"][0])
    base = ctx.seed
    arms = ["rest"] + [f"v{v:g}" for v in p["vs"]] + ["kick"]
    ctrl = [f"tele{g:g}_v{v0:g}" for g in p["gammas"]] + [f"tele{g:g}_kick" for g in p["gammas"]]
    ctrl += [f"V0_v{v0:g}", f"blind_v{v0:g}"]
    jobs = []
    for k in range(int(p["i4_n"])):
        jobs.append(("i4", p, k, base))
    for k in range(n_dec):
        for arm in arms + (ctrl if k < n_ctrl else []):
            jobs.append(("wall", p, rd, k, arm, base))
    for k in range(n_wave):
        for what in p["waves"]:
            jobs.append(("wave", p, rd, k, what, base))
    for k in range(n_i3):
        jobs.append(("wall", p, ri3, k, "rest", base))
    jobs.sort(key=lambda j: (0 if j[0] == "i4" else 1, -(j[2] if j[0] != "i4" else 0)))
    res = parallel_map(run_task, jobs, ctx.workers, ctx.progress,
                       lambda i, r, d, nn: f"[F3.15] {d}/{nn}: {r['kind']} {r.get('arm', r.get('what', ''))} "
                                           f"rho {r.get('rho', '')} ({time.time() - t0:.0f}s)")

    def W(rho, arm, m):
        out = {r["k"]: r for r in res if r["kind"] == "wall" and r["rho"] == rho and r["arm"] == arm}
        return [out[k] for k in range(m)]

    lines = [f"F3.15 (B1, 2+1, intrinsic weight variant B; a {p['calib'][0]:.4f}, c {p['calib'][1]:.4f}, beta "
             f"{p['calib'][2]:.4f}); deciding rho {rd:g}, T {T:g}, {n_dec} sprinklings; controls on {n_ctrl}; base seed "
             f"{base}; mean +- SE over sprinklings"]
    inst, passc, failc = {}, {}, {}

    # ---------------------------------------------------------------- instrument
    I4 = [r for r in res if r["kind"] == "i4"]
    worst = max(max(r["kernel"], r["order"], r["region"]) for r in I4)
    inst["I4 kernel, order, region"] = worst < float(p["i4_tol"])
    lines.append(f"  I4: forward/backward kernel {max(r['kernel'] for r in I4):.1e}, linear extension "
                 f"{max(r['order'] for r in I4):.1e}, region {max(r['region'] for r in I4):.1e}")
    WV = [r for r in res if r["kind"] == "wave"]
    vac = [r for r in WV if r["what"] == "vacuum"]
    tmv = np.arange(0.5, T, 1.0)
    lam = _ms([np.polyfit(tmv, np.log(np.array(r["A"]) / np.array(r["Ae"])), 1)[0] for r in vac])
    inst["I1 vacuum"] = abs(lam[0]) + nse * lam[1] < float(p["tol_I1"])
    lines.append(f"  I1 vacuum: lambda_eff {_fmt(lam)}")
    for what in p["waves"]:
        if what == "vacuum":
            continue
        R = [r for r in WV if r["what"] == what]
        d = _ms([r["omega"] / r["omega_c"] - 1 for r in R])
        inst[f"I2 {what}"] = abs(d[0]) + nse * d[1] < float(p["tol_I2"])
        lines.append(f"  I2 {what}: omega/omega_c - 1 {_fmt(d)}")

    def wrel(R):
        return [np.nanmean(_col(r["rows"], "w")) / np.nanmean(_col(r["rows"], "w_ex")) - 1 for r in R]

    def eps(R):
        return [np.nanmean(np.abs(_col(r["rows"], "w") - _col(r["rows"], "w_ex")) / _col(r["rows"], "w_ex")) for r in R]

    R3 = W(ri3, "rest", n_i3)
    d3 = _ms(wrel(R3))
    inst["I3 rest wall at rho_i3"] = abs(d3[0]) + nse * d3[1] < float(p["tol_I3"])
    rest = W(rd, "rest", n_dec)
    lines.append(f"  I3 rho {ri3:g}: w/w_ex - 1 {_fmt(d3)}; reported: eps (mean |w - w_ex|/w_ex) rho {rd:g} "
                 f"{_fmt(_ms(eps(rest)))}, rho {ri3:g} {_fmt(_ms(eps(R3)))}")

    Mrest = np.array([np.nanmean(_col(r["rows"], "E")) for r in rest])
    wrest = np.array([np.nanmean(_col(r["rows"], "w")) for r in rest])
    q = _ms([np.nanmean(_col(r["rows"], "Q")) - np.nanmean(_col(r["rows"], "Q_ex")) for r in rest])
    dr = _ms([_col(r["rows"], "xc")[-1] - _col(r["rows"], "xc")[0] for r in rest])
    wr = _ms(wrel(rest))
    er = _ms([(np.nanmean(_col(r["rows"], "E")) - np.nanmean(_col(r["rows"], "E_ex"))) / 8.0 for r in rest])
    inst["rest Q"] = abs(q[0]) + nse * q[1] < float(p["tol_Q"])
    inst["rest drift"] = abs(dr[0]) + nse * dr[1] < float(p["tol_drift"])
    inst["rest w"] = abs(wr[0]) + nse * wr[1] < float(p["tol_w"])
    inst["rest E"] = abs(er[0]) + nse * er[1] < float(p["tol_E"])
    ep = _ms([(np.nanmean(_col(r["rows"], "E_plain")) - 8.0) / 8.0 for r in rest])
    lines.append(f"  rest: Q - Q_ex {_fmt(q)}; drift {_fmt(dr)}; w/w_ex - 1 {_fmt(wr)}; (E - E_ex)/8 {_fmt(er)}; "
                 f"M_rest {Mrest.mean():.4f}; reported: plain E/8 - 1 {_fmt(ep)}")

    # ---------------------------------------------------------------- inertia (each v separately)
    for v in p["vs"]:
        g = 1 / math.sqrt(1 - v * v)
        R = W(rd, f"v{v:g}", n_dec)
        dv = _ms([_slope(r["rows"]) - v for r in R])
        C = _ms([np.nanmean(_col(r["rows"], "w")) * g / w0 - 1 for r, w0 in zip(R, wrest)])
        wv = _ms(wrel(R))
        eB = _ms([(np.nanmean(_col(r["rows"], "E")) - np.nanmean(_col(r["rows"], "E_ex"))) / (8 * g) for r in R])
        qv = _ms([np.nanmean(_col(r["rows"], "Q")) - np.nanmean(_col(r["rows"], "Q_ex")) for r in R])
        for lab, m, tol in (("velocity", dv, "tol_v"), ("contraction C", C, "tol_C"), ("width w/w_ex", wv, "tol_w"),
                            ("energy e_B", eB, "tol_E")):
            passc[f"v {v:g} {lab}"] = abs(m[0]) + nse * m[1] < float(p[tol])
            failc[f"v {v:g} {lab}"] = abs(m[0]) - nse * m[1] > float(p[tol])
        passc[f"v {v:g} Q"] = abs(qv[0]) + nse * qv[1] < float(p["tol_Q"])
        lines.append(f"  v {v:g}: v_meas - v {_fmt(dv)}; C - 1 {_fmt(C)}; w/w_ex - 1 {_fmt(wv)}; e_B {_fmt(eB)}; "
                     f"Q - Q_ex {_fmt(qv)}")

    # ---------------------------------------------------------------- masses
    K_ = W(rd, "kick", n_dec)
    tk = float(p["t_kick_read"])
    f_, tau = float(p["f"]), float(p["t2"]) - float(p["t1"])
    dp = np.array([_gv(_slope(k["rows"], tk)) - _gv(_slope(r["rows"], tk)) for k, r in zip(K_, rest)])
    Qk = float(np.nanmean([np.nanmean(_col(k["rows"], "Q")) for k in K_]))
    if np.all(np.isfinite(dp)) and np.all(np.isfinite(Mrest)):
        rM = 2 * math.pi * Qk * f_ * tau / dp.mean() / Mrest.mean()
        rng = np.random.default_rng(derive_seed(base, "F3.15-boot"))
        boot = []
        for _ in range(int(p["boot_B"])):
            i = rng.integers(0, dp.size, dp.size)
            boot.append(2 * math.pi * Qk * f_ * tau / dp[i].mean() / Mrest[i].mean())
        rMs = float(np.std(boot, ddof=1))
    else:
        rM, rMs = float("nan"), float("nan")
    passc["r_M"] = abs(rM - 1.0) + nse * rMs < float(p["tol_M"])
    failc["r_M"] = abs(rM - 1.0) - nse * rMs > float(p["tol_M"])
    after = _ms([_gv(_slope(k["rows"], tk)) for k in K_])
    lines.append(f"  kick: Delta(gamma v) {_fmt(_ms(dp))} (continuum {2 * math.pi * f_ * tau / 8:.4f}); Q_kick {Qk:.4f}; "
                 f"r_M = M_in / M_rest {rM:.4f} +- {rMs:.4f}; reported: gamma v after the kick {_fmt(after)}")

    # ---------------------------------------------------------------- controls
    cons_v = W(rd, f"v{v0:g}", n_ctrl)
    cons_k = K_[:n_ctrl]
    for g in p["gammas"]:
        for motion, cons in ((f"v{v0:g}", cons_v), ("kick", cons_k)):
            tel = W(rd, f"tele{g:g}_{motion}", n_ctrl)
            pc = np.array([_gv(_slope(r["rows"], tk)) for r in cons])
            pg = np.array([_gv(_slope(r["rows"], tk)) for r in tel])
            d = _ms(pc - pg)
            ok = (d[0] - nse * d[1] > 0.0) and (np.mean(pg) <= float(p["drop"]) * np.mean(pc))
            inst[f"telegraph Gamma {g:g} {motion}"] = bool(ok)
            lines.append(f"  telegraph Gamma {g:g} {motion}: p_Gamma {np.mean(pg):+.4f}, p_cons {np.mean(pc):+.4f}, "
                         f"p_cons - p_Gamma {_fmt(d)}")

    def Lr(r):
        e = _col(r["rows"], "E_loc")
        return e[-1] / e[0]

    Lc = np.array([Lr(r) for r in cons_v])
    for name, key in (("V = 0 readout check", f"V0_v{v0:g}"), ("blind control", f"blind_v{v0:g}")):
        ctl = W(rd, key, n_ctrl)
        Lx = np.array([Lr(r) for r in ctl])
        d = _ms(Lc - Lx)
        wl = float(np.nanmean([np.nanmean(_col(r["rows"], "w")[-3:]) for r in ctl]))
        vl = float(np.nanmean([_slope(r["rows"]) for r in ctl]))
        lines.append(f"  {name}: L {np.mean(Lx):.4f} (conservative {np.mean(Lc):.4f}); D_L {_fmt(d)}; reported: width "
                     f"(last 3 read slabs) {wl:.2f}, velocity {vl:+.3f}")
        sep = bool(d[0] - nse * d[1] > float(p["tol_DL"]))
        if key.startswith("V0"):
            inst["V = 0 readout check"] = sep
        else:
            passc["blind control separated"] = sep

    # ---------------------------------------------------------------- verdict
    inst = {k: bool(v) for k, v in inst.items()}
    passc = {k: bool(v) for k, v in passc.items()}
    failc = {k: bool(v) for k, v in failc.items()}
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
    lines.append(f"cost: {np.mean([r['sec'] for r in res if r['kind'] == 'wall' and r['rho'] == rd]):.0f} s per wall arm "
                 f"at rho {rd:g} ({time.time() - t0:.0f}s)")
    return Result(status, lines, {"instrument": inst, "passed": passc, "failed": failc, "results": res})
