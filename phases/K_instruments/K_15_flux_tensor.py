"""
K.15 -- Meter of the tensor of carried flux: energy density, momentum flux, pressure and shear from chains (3+1)
==============================================================================================================

Why this meter exists (23 Oct 2026, the author's design)
--------------------------------------------------------
The focusing principle (docs/derivations/spatial_part.md) needs the flux of
matter energy in every light direction, T_kk = T_ab k^a k^b. A count of links
crossing a cut is a scalar; pressure and shear are directional. The tensor
needs the carried 4-velocity u of the flux (F3.6, replicated). This meter
builds the candidate tensor from the chains, decomposes it for an observer
and checks it on populations whose answer is known. It also checks the
conservation of 4-momentum on cuts, which the focusing principle requires
(Bianchi identity), for the current merge rule and for the derived one.

Three layers (the author's design)
----------------------------------
  microscopic state  chains of the carried rule (F3.6): flux F per event,
                     carried 4-velocity u (upper index, signature + - - -);
  local measurement  projections on an observer n (unit timelike):
                     E = T_ab n^a n^b, q^a = h^a_c T^cb n_b, S^ab = h^a_c h^b_d T^cd,
                     h_ab = eta_ab - n_a n_b (spatial projector, + - - - so h is
                     negative definite on space; the spatial components below
                     are taken positive), p = (1/3) S^i_i, pi^ij = S^ij - p delta^ij;
  effective tensor   T^ab = (1 / dV4) sum over the links starting in a cell of
                     F u^a u^b tau_link   (the discrete form of the dust tensor
                     sum_a m_a integral u^a u^b d tau delta^4: every chain counted
                     by the proper time it spends in the cell; invariant: tau and
                     dV4 are Lorentz invariant, no cut, no double counting).
What F is: with the current merge rule (flux = the sum of the incoming
fluxes, direction P/|P|) F is a count of transfers and T is a FLUX tensor;
the 4-vector F u is not conserved at a merger of different directions (the
excess |P| - F is dropped, F3.10 unit). With the derived merge rule F = |P|
(rest mass grows at an inelastic merger, as in special relativity) F u is
the conserved 4-momentum and T is a candidate energy-momentum tensor. The
physical meaning of F (energy in which units) is the author's decision; the
meter reports both rules.

Populations (3+1, rho = 1, no field, no kicks: free transport)
--------------------------------------------------------------
A box t in [0, T], |x_i| <= half; source events in a central ball (radius
r_src = 4, about 80 source events) in a thin layer t < t_src = 0.3 emit F = 1 each with a direction
drawn per event (thin, so that a chain rarely passes another source event and
mixes with its emission: with t_src = 1 the first quick run lost 9 % of the
counter-stream pressure to such mixing at the sources):
  beam     all u = gamma (1, v, 0, 0), v = 0.6;
  counter  u = gamma (1, +-v, 0, 0) with probability 1/2 each;
  iso      u = gamma (1, v n), n uniform on the sphere, v = 0.6.
Analytic expectations for the directions as emitted (before any merger),
with rho = T^a_a (the trace, the invariant flux density):
  beam     in the rest frame of the beam (Landau frame: the timelike
           eigenvector of T^a_b): p = 0, pi = 0 (dust); in the lab
           E / rho = gamma^2, |q| / rho = gamma^2 v, S_xx / rho = gamma^2 v^2;
  counter  lab: q_x = 0, S_xx / E = v^2, S_yy = S_zz = 0, p / E = v^2 / 3,
           pi_xx / E = 2 v^2 / 3, pi_yy = pi_zz = -v^2 / 3 E;
  iso      lab: q_i = 0, pi_ij = 0, p / E = v^2 / 3.
Sampling noise: the directions are random per source event, so the measured
ratios fluctuate around the expectation; the meter averages over sprinklings
and tests components (q_x, pi_xy, ...), not magnitudes (|q|, max|pi| are
biased upwards by noise; changed after the first quick run; for the same
reason the emission ball is radius 4 and 64 sprinklings: a component of q/E of
an iso population fluctuates by ~0.07 per sprinkling in the quick runs). Mergers of
different directions (reported: the flux in merged events in the window and
the 'mixed' flux, whose speed is no longer v) lower the pressure of counter
and iso; in the quick runs p/E was 5-8 % below v^2/3 (counter 0.110, iso 0.111
against 0.12), inside the tolerance; it is physics of the population (inelastic
mergers cool it), not an error of the tensor, which the beam checks exactly.

Parts and criteria (meter: PASSED or INSTRUMENT_FAILED; a meter does not test
the theory)
-----------
beam, counter, iso  ratios over the measurement window t in [t_lo, t_hi]
         (all links starting there; the window ends before the flux reaches the
         box wall), for both merge rules. PASSED if every expected ratio holds
         within its tolerance: |measured - expected| + 2 SE <= 0.05 (in units of
         E; the ratios are of order 0.1-0.5; 0.05 is a tenth of the largest),
         for the derived merge rule F = |P|. The current rule is reported (it
         differs only where mergers of different directions happen).
         Reported also: the fraction of flux in merged events, the Landau-frame
         pressure (rest-frame pressure of the population), the ratio for cells
         of side 2 and 3 (spread over cells, stability with the cell size), and
         the frame check (T from boosted coordinates and directions equals the
         transformed tensor, relative 1e-9).
conserve the total 4-momentum on cuts t = const (sum of F u over the links
         crossing the cut) at several times in the window, for both rules, iso
         population (where mergers of different directions happen). PASSED if
         with the rule F = |P| the relative change is <= 1e-9 (exact
         conservation by construction, up to flux that ends at the box wall,
         which is reported and must be 0 in the window). The loss of the current
         rule is reported (no verdict: it is the known deficit).
Replication with --seed 1.

Limits (declared): the tensor is kinetic (flux only); a reservoir, a stress of
the order itself, or the energy of a field would need further terms, which
the model does not define. Nothing here tests the Einstein equations or the
focusing principle; the meter only provides T_kk for them.
"""

from __future__ import annotations

import math
import time

import numpy as np

from core.order.links import links_csr
from core.order.sprinkle import Box, derive_seed, sprinkle_poisson
from core.stats.parallel import parallel_map
from core.stats.stats import boot_se
from core.testing import Result

TEST = {
    "id": "K.15",
    "phase": "K",
    "title": "Instrument of the carried-flow tensor: energy density, momentum flux, pressure and shear from chains "
             "(3+1)",
    "theory_part": "the effective energy–momentum tensor of the carried flow (F3.6) for the focusing principle (T_kk)",
    "question": "Does T^ab = (1/ΔV₄) Σ F u^a u^b τ_link on populations with a known answer (a beam, "
                "counter-streaming beams, an isotropic gas) give the correct energy density, momentum flux, pressure "
                "and shear, and is the 4-momentum conserved on slices under the merge rule F = |P|?",
    "method": "3+1, free transport by the F3.6 rule without a field; sources in a central ball emit directions "
              "according to the population; tensor from links weighted by proper time, decomposition for an observer "
              "(E, q, p, π), Landau frame; two merge rules (current: flow = sum of flows; derived: flow = |P|); "
              "4-momentum on slices t = const",
    "prediction": "beam: dust in its rest frame (p = π = 0); counter-streaming: p/E = v²/3, π_xx/E = 2v²/3; "
                  "isotropic: π = 0, p/E = v²/3; with F = |P| the 4-momentum on slices is conserved exactly",
    "criterion": "beam/counter/iso: PASSED when every expected ratio holds |measured − expected| + 2 SE <= 0.05 (in "
                 "units of E) under the rule F = |P|, otherwise INSTRUMENT; conserve: PASSED when under F = |P| the "
                 "4-momentum on slices changes by at most 1e-9 relatively and no flow in the window ends at the "
                 "wall, otherwise INSTRUMENT. Replication --seed 1.",
    "null_model": "the current merge rule (flow = sum of flows), reported alongside",
    "level": "calibration",
    "meters": [],
    "depends_on": ["F3.6"],
    "params": {"D": 4, "T": 7.0, "half": 8.0, "r_src": 4.0, "t_src": 0.3, "v": 0.6, "win": [2.5, 4.5],
               "cuts": [2.5, 3.0, 3.5, 4.0, 4.5], "cells": [2.0, 3.0], "beta_check": 0.6, "n_seeds": 64, "n_se": 2.0,
               "boot_B": 2000, "tol": 0.05},
    "parts": {"beam": {"mode": "beam"}, "counter": {"mode": "counter"}, "iso": {"mode": "iso"},
              "conserve": {"mode": "conserve"}},
    "quick": {"n_seeds": 2, "boot_B": 200, "half": 6.0, "r_src": 3.0, "T": 5.5, "win": [2.0, 3.5], "cuts": [2.0, 2.75, 3.5]},
    "estimate": "each part ~5 min (10 processes; box ~29 000 events, ~40 s per sprinkling; 64 sprinklings)",
}

ETA = np.diag([1.0, -1.0, -1.0, -1.0])


def mink4(a, b):
    return a[..., 0] * b[..., 0] - (a[..., 1:] * b[..., 1:]).sum(axis=-1)


def boost_matrix(beta):
    g = 1.0 / math.sqrt(1.0 - beta * beta)
    B = np.eye(4)
    B[0, 0] = B[1, 1] = g
    B[0, 1] = B[1, 0] = -g * beta
    return B


# ============================================================================================== transport

def transport(X, ptr, idx, cptr, cidx, F_src, U_src, merge):
    """
    Free transport by the carried rule of F3.6 (no field). F_src, U_src: flux and 4-velocity emitted at every event
    (0 where none). merge = 'sum' (current: flux = sum of the incoming fluxes, direction P/|P|) or 'mass'
    (derived: flux = |P|, so that F u = P is conserved exactly). Returns flux, direction, chosen child, link tau,
    and the incoming |P| - sum F (the 4-momentum dropped by the 'sum' rule).
    """
    N = len(X)
    inF = np.zeros(N)                          # sum of the incoming fluxes
    inP = np.zeros((N, 4))                     # sum of the incoming 4-momenta F u
    flux = np.zeros(N)
    dirs = np.zeros((N, 4))
    dec = np.full(N, -1, dtype=np.int64)
    tau_l = np.zeros(N)
    excess = np.zeros(N)
    for v in range(N):                         # increasing index = a linear extension of the causal order
        Fin = inF[v] + F_src[v]
        if Fin <= 0.0:
            continue
        P = inP[v] + F_src[v] * U_src[v]
        m = math.sqrt(max(float(mink4(P, P)), 1e-300))   # invariant mass of the incoming 4-momentum
        u = P / m
        tot = m if merge == "mass" else Fin
        excess[v] = m - Fin                    # > 0 when different directions merge
        flux[v], dirs[v] = tot, u
        a, b = cptr[v], cptr[v + 1]
        if b == a:
            continue                           # no child: the flux ends here (box wall or end of time)
        ks = cidx[a:b]
        B = X[ks] - X[v]
        tau = np.sqrt(np.maximum(mink4(B, B), 1e-300))
        j = int(np.argmin(mink4(u[None, :], B) / tau))   # least rapidity relative to u (F3.6)
        k = int(ks[j])
        dec[v], tau_l[v] = k, float(tau[j])
        inF[k] += tot
        inP[k] += tot * u
    return {"flux": flux, "dir": dirs, "dec": dec, "tau": tau_l, "excess": excess}


# ============================================================================================== the tensor and its decomposition

def tensor(X, r, sel, dV4):
    """T^ab = (1 / dV4) sum over the links starting at the selected events of F u^a u^b tau_link."""
    m = sel & (r["flux"] > 0.0) & (r["dec"] >= 0)
    w = r["flux"][m] * r["tau"][m]
    U = r["dir"][m]
    return np.einsum("n,na,nb->ab", w, U, U) / dV4


def decompose(T, n):
    """E, q (3), S (3x3), p, pi (3x3) for the observer n (upper index), components in an orthonormal frame of n."""
    # orthonormal tetrad: e0 = n, e1..e3 spatial (Gram-Schmidt against the coordinate axes), all upper index
    E_ = [n / math.sqrt(float(n @ ETA @ n))]
    for k in range(1, 4):
        e = np.zeros(4)
        e[k] = 1.0
        for f in E_:
            e = e - (e @ ETA @ f) / (f @ ETA @ f) * f
        E_.append(e / math.sqrt(-float(e @ ETA @ e)))
    E_ = np.array(E_)
    Tl = ETA @ T @ ETA                                   # lower indices
    C = E_ @ Tl @ E_.T                                   # frame components T_(i)(j) = T_ab e_i^a e_j^b
    En = C[0, 0]
    q = C[0, 1:]                                         # momentum density (sign: flux along +e_i positive)
    S = C[1:, 1:]                                        # stress, positive for pressure
    p = float(np.trace(S)) / 3.0
    pi = S - p * np.eye(3)
    return En, q, S, p, pi


def landau_frame(T):
    """Timelike eigenvector of T^a_b (the rest frame of the population), normalised, future pointing."""
    w, V = np.linalg.eig(T @ ETA)                        # T^a_c eta_cb
    best = None
    for i in range(4):
        v = np.real(V[:, i])
        nn = float(v @ ETA @ v)
        if nn > 0 and (best is None or np.real(w[i]) > best[0]):
            best = (np.real(w[i]), v / math.sqrt(nn) * (1.0 if v[0] > 0 else -1.0))
    return best[1] if best else np.array([1.0, 0.0, 0.0, 0.0])


# ============================================================================================== task

def make(p, seed):
    rng = np.random.default_rng(seed)
    h, T = float(p["half"]), float(p["T"])
    X = sprinkle_poisson(Box(4, [0.0, -h, -h, -h], [T, h, h, h]), 1.0, rng)
    ptr, idx = links_csr(X, 1.0)
    child = np.repeat(np.arange(len(X)), np.diff(ptr))
    assert bool((idx < child).all())
    o = np.argsort(idx, kind="stable")
    cptr = np.concatenate([[0], np.cumsum(np.bincount(idx, minlength=len(X)))]).astype(np.int64)
    return X, ptr, idx, cptr, child[o], rng


def sources(X, p, mode, rng):
    """Emitted flux and 4-velocity per event for the population (module notes)."""
    v = float(p["v"])
    g = 1.0 / math.sqrt(1.0 - v * v)
    N = len(X)
    on = (X[:, 0] < float(p["t_src"])) & ((X[:, 1:] ** 2).sum(axis=1) <= float(p["r_src"]) ** 2)
    F = on.astype(float)
    U = np.zeros((N, 4))
    if mode == "beam":
        dirs = np.tile([1.0, 0.0, 0.0], (N, 1))
    elif mode == "counter":
        dirs = np.zeros((N, 3))
        dirs[:, 0] = np.where(rng.random(N) < 0.5, 1.0, -1.0)
    else:                                                 # iso (and conserve)
        dirs = rng.normal(size=(N, 3))
        dirs /= np.linalg.norm(dirs, axis=1)[:, None]
    U[:, 0] = g
    U[:, 1:] = g * v * dirs
    return F, U


def task(p, seed, mode):
    X, ptr, idx, cptr, cidx, rng = make(p, seed)
    pop = "iso" if mode == "conserve" else mode
    Fs, Us = sources(X, p, pop, rng)
    out = {"N": len(X)}
    lo, hi = p["win"]
    h = float(p["half"])
    sel = (X[:, 0] >= lo) & (X[:, 0] < hi)
    dV4 = (hi - lo) * (2.0 * h) ** 3
    for merge in ("mass", "sum"):
        r = transport(X, ptr, idx, cptr, cidx, Fs, Us, merge)
        T = tensor(X, r, sel, dV4)
        lab = decompose(T, np.array([1.0, 0.0, 0.0, 0.0]))
        nL = landau_frame(T)
        rest = decompose(T, nL)
        car = (r["flux"] > 0.0) & sel
        res = {"E": lab[0], "q": lab[1].tolist(), "S": lab[2].tolist(), "p": lab[3], "pi": lab[4].tolist(),
               "rho": float(np.trace(T @ ETA)), "p_rest": rest[3], "pi_rest": float(np.max(np.abs(rest[4]))),
               "E_rest": rest[0],
               "mixed": float(r["flux"][car & (np.abs(np.sqrt(np.maximum((r["dir"][:, 1:] ** 2).sum(axis=1), 0.0))
                                                         / np.maximum(r["dir"][:, 0], 1e-300) - float(p["v"])) > 1e-6)].sum()
                              / max(r["flux"][car].sum(), 1e-300)),
               "merged_frac": float(r["flux"][car & (r["excess"] > 1e-12)].sum() / max(r["flux"][car].sum(), 1e-300)),
               "wall": float(r["flux"][(np.diff(cptr) == 0) & (X[:, 0] < hi) & (r["flux"] > 0)].sum())}
        # p/E over cells of the given side (spread across cells and stability with the cell size)
        cells = {}
        for L in p["cells"]:
            ratios = []
            edges = np.arange(-h, h - 1e-9, float(L))
            for x0 in edges:
                for y0 in edges:
                    for z0 in edges:
                        m = sel & (X[:, 1] >= x0) & (X[:, 1] < x0 + L) & (X[:, 2] >= y0) & (X[:, 2] < y0 + L) \
                            & (X[:, 3] >= z0) & (X[:, 3] < z0 + L)
                        Tc = tensor(X, r, m, (hi - lo) * L ** 3)
                        if Tc[0, 0] > 0.2 * T[0, 0] * (2 * h / L) ** 3 / 8:      # cells with enough flux
                            dc = decompose(Tc, np.array([1.0, 0.0, 0.0, 0.0]))
                            ratios.append(dc[3] / dc[0])
            cells[str(L)] = [float(np.mean(ratios)) if ratios else float("nan"),
                             float(np.std(ratios)) if len(ratios) > 1 else float("nan"), len(ratios)]
        res["cells"] = cells
        # frame check: T from boosted coordinates and directions equals the transformed tensor (same links)
        Bm = boost_matrix(float(p["beta_check"]))
        rb = dict(r)
        rb["dir"] = r["dir"] @ Bm.T
        Tb = tensor(X @ Bm.T, rb, sel, dV4)
        res["frame_dev"] = float(np.max(np.abs(Tb - Bm @ T @ Bm.T)) / np.max(np.abs(T)))
        # 4-momentum on cuts t = c: sum of F u over the links crossing the cut
        Pc = []
        for c in p["cuts"]:
            dec = r["dec"]
            m = (r["flux"] > 0.0) & (dec >= 0) & (X[:, 0] < c)
            m[m] = X[dec[m], 0] >= c
            Pc.append((r["flux"][m][:, None] * r["dir"][m]).sum(axis=0).tolist())
        res["cuts"] = Pc
        out[merge] = res
    return out


# ============================================================================================== run

def expected(mode, v):
    """Expected ratios in units of E (module notes): name -> value."""
    g2 = 1.0 / (1.0 - v * v)
    if mode == "beam":
        return {"p_rest/E_rest": 0.0, "max|pi_rest|/E_rest": 0.0, "|q|/E": v, "S_xx/E": v * v, "E/rho": g2}
    if mode == "counter":
        return {"q_x/E": 0.0, "S_xx/E": v * v, "S_yy/E": 0.0, "p/E": v * v / 3.0, "pi_xx/E": 2.0 * v * v / 3.0}
    return {"q_x/E": 0.0, "q_y/E": 0.0, "q_z/E": 0.0, "pi_xx/E": 0.0, "pi_yy/E": 0.0, "pi_xy/E": 0.0,
            "p/E": v * v / 3.0}


def measured(mode, r):
    E = r["E"]
    S = np.array(r["S"])
    pi = np.array(r["pi"])
    q = np.array(r["q"])
    if mode == "beam":
        return {"p_rest/E_rest": r["p_rest"] / r["E_rest"], "max|pi_rest|/E_rest": r["pi_rest"] / r["E_rest"],
                "|q|/E": float(np.linalg.norm(q)) / E, "S_xx/E": S[0, 0] / E, "E/rho": E / r["rho"]}
    if mode == "counter":
        return {"q_x/E": q[0] / E, "S_xx/E": S[0, 0] / E, "S_yy/E": S[1, 1] / E, "p/E": r["p"] / E,
                "pi_xx/E": pi[0, 0] / E}
    return {"q_x/E": q[0] / E, "q_y/E": q[1] / E, "q_z/E": q[2] / E, "pi_xx/E": pi[0, 0] / E, "pi_yy/E": pi[1, 1] / E,
            "pi_xy/E": pi[0, 1] / E, "p/E": r["p"] / E}


def run(p, ctx):
    mode = str(p["mode"])
    t0 = time.time()
    n, B, nse, tol, v = int(p["n_seeds"]), int(p["boot_B"]), float(p["n_se"]), float(p["tol"]), float(p["v"])
    res = parallel_map(task, [(p, derive_seed(ctx.seed, "K.15", mode, k), mode) for k in range(n)], ctx.workers,
                       ctx.progress, lambda i, r_, dd, nn: f"[K.15] {mode} sprinkling {dd}/{nn} ({time.time()-t0:.0f}s)")
    lines = [f"3+1, rho = 1, box t [0, {p['T']:g}], |x_i| <= {p['half']:g} (~{np.mean([r_['N'] for r_ in res]):.0f} events); "
             f"sources in a ball of radius {p['r_src']:g}, t < {p['t_src']:g}, v = {v}; window t in {p['win']}; "
             f"T^ab = (1/dV4) sum F u^a u^b tau_link"]
    inst = {}
    if mode in ("beam", "counter", "iso"):
        exp_ = expected(mode, v)
        ok = True
        for merge in ("mass", "sum"):
            meas = [measured(mode, r_[merge]) for r_ in res]
            parts = []
            for k_, e_ in exp_.items():
                m_, s_ = boot_se([mm[k_] for mm in meas], B, derive_seed(ctx.seed, "K.15-b", mode, merge, k_))[:2]
                good = abs(m_ - e_) + nse * s_ <= tol
                if merge == "mass":
                    ok &= good
                parts.append(f"{k_} {m_:+.4f} +- {s_:.4f} (exp {e_:+.4f}{'' if merge != 'mass' else (', OK' if good else ', FAILS')})")
            mf = np.mean([r_[merge]["merged_frac"] for r_ in res])
            cl = "; ".join(f"cells L = {L}: p/E {np.nanmean([r_[merge]['cells'][L][0] for r_ in res]):.4f}, spread "
                           f"{np.nanmean([r_[merge]['cells'][L][1] for r_ in res]):.4f} over "
                           f"{np.mean([r_[merge]['cells'][L][2] for r_ in res]):.0f} cells" for L in map(str, p["cells"]))
            lines.append(f"  merge rule {'F = |P|' if merge == 'mass' else 'sum of fluxes'}: " + "; ".join(parts))
            lines.append(f"     flux in merged events {mf:.3f}, mixed flux (speed no longer v) "
                         f"{np.mean([r_[merge]['mixed'] for r_ in res]):.3f}; Landau-frame pressure p_rest/E_rest "
                         f"{np.mean([r_[merge]['p_rest'] / r_[merge]['E_rest'] for r_ in res]):+.4f}; {cl}")
        inst["expected ratios (merge rule F = |P|)"] = bool(ok)
        inst["frame check of the tensor (relative 1e-9)"] = max(r_[m]["frame_dev"] for r_ in res for m in ("mass", "sum")) <= 1e-9
    else:
        for merge in ("mass", "sum"):
            rel = []
            for r_ in res:
                Pc = np.array(r_[merge]["cuts"])
                rel.append(float(np.max(np.abs(Pc - Pc[0]))) / float(np.max(np.abs(Pc[0]))))
            dE = [(np.array(r_[merge]["cuts"])[-1, 0] - np.array(r_[merge]["cuts"])[0, 0]) / np.array(r_[merge]["cuts"])[0, 0]
                  for r_ in res]
            lines.append(f"  merge rule {'F = |P|' if merge == 'mass' else 'sum of fluxes'}: largest relative change of the "
                         f"4-momentum on the cuts {max(rel):.2e} (mean {np.mean(rel):.2e}); change of the energy component "
                         f"from the first to the last cut {np.mean(dE):+.2e}; flux ending at the wall in the window "
                         f"{max(r_[merge]['wall'] for r_ in res):.3g}; flux in merged events "
                         f"{np.mean([r_[merge]['merged_frac'] for r_ in res]):.3f}")
            if merge == "mass":
                inst["4-momentum conserved on the cuts (F = |P|, relative 1e-9)"] = max(rel) <= 1e-9
                inst["no flux ends at the wall in the window"] = max(r_[merge]["wall"] for r_ in res) == 0.0
    status = "PASSED" if all(inst.values()) else "INSTRUMENT_FAILED"
    lines.append("checks: " + "; ".join(f"{k_}: {'OK' if v_ else 'FAILS'}" for k_, v_ in inst.items()))
    return Result(status, lines, {"instrument": inst})
