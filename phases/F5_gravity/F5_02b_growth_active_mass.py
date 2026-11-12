"""
F5.2b -- Active gravitational mass from a causal growth rule of the order (1+1, scalar response)
===============================================================================================

Why this test exists (11 Nov 2026; card docs/F5_2b_card_draft.md; direction A approved by the author at 11:59)
-------------------------------------------------------------------------------------------------------------
F3.14b (inertia) and F5.2a3 (universality of free fall in an order of prescribed density) are PASSED and replicated. The
active half of H5 needs a dynamics of the order: matter must change the order so that other matter falls towards it.
Candidate (the author's corrections included): every candidate event x of a Poisson sprinkling of proper density rho
carries an independent mark u_x fixed at the sprinkling and is realised iff u_x < Pi(x)^2 = exp(2 L(x)), with
    L(x) = -(kappa / 2 rho) * sum over REALISED y < x of q(y),   q = 2 V_A(phi_A) + 2 V_B(phi_B)
(the load = trace of T per event; the normalisation of the retarded Green function of 1+1). The realised events have
coordinate density rho Pi^2, so the continuum law is (d_t^2 - d_x^2) L = -kappa Pi^2 q, i.e. box_g L = -kappa q in
ds^2 = Pi^2 (dt^2 - dx^2). In 1+1 a Poisson process of density rho Pi^2 is a sprinkling of that metric; the causal
relation stays Minkowski (conformal flatness), and the geometry is in the counts (proper times). Two separate fields
(A: mu_A, B: mu_B) evolve on the realised events with the unchanged rule of F3.14a and interact only through the order.
Development rows (findings, Discoveries, 11 Nov): B falls towards A as 1 : 1.97 : 3.86 for M_A 4 : 8 : 16, independently
of M_B, within 2-5 % of the continuum.

Model (frozen with the card)
----------------------------
rho 800, kappa 0.003, kinks at rest at x_A = -3, x_B = +3, T 20, trapezoid |x| <= T + 16 - t. Combinations (mu_A, mu_B):
(0.25, 1), (1, 1), (4, 1), (1, 0.25), i.e. (M_A, M_B) = (4, 8), (8, 8), (16, 8), (8, 4); and the boosted launch (A2):
(1, 1) with both kinks at v = 0.3. Every run is paired with kappa = 0 on the same candidates and marks (every candidate
realised). All combinations of one seed use the same candidates and marks. References: the coupled continuum equations
by finite differences, the same initial data (boosted for A2), centre by the pi crossing.

Readout (frozen)
----------------
Kink fit of F3.14a per slab of width 1 for each field, tracking the centre. dB = [x_B(19.5) - x_B(0.5)]_kappa -
[same]_kappa=0 (towards A is negative); dA likewise; for A2 the change of the separation x_B - x_A. r = mean(d) /
d_cont (ratio of means; SE by bootstrap over sprinklings). Universality: u = r_B(M_B 4) - r_B(M_B 8) at M_A 8, paired
by seed (bootstrap). Density: realised fraction / mean Pi^2 per x bin (t > 2T/3). Intervals: for 60 random realised
timelike pairs near A (|x - x_A| < 4, t > T/2, coordinate area > 1) on the first 8 sprinklings of (M_A 8, M_B 8): the
realised count inside / (rho x proper area int Pi^2) and / (rho x coordinate area).

Criteria (frozen with the card; 256 new sprinklings; 2 SE; every condition on its own)
--------------------------------------------------------------------------------------
INSTRUMENT_FAILED if any fails:
  A0  on the first 4 sprinklings, a change of the marks of the latest 10 % of candidates leaves every earlier event
      (realisation and fields) unchanged;
  A1  on the first 4 sprinklings, the linear extension t + 0.5 x gives the same realised set and fields within 1e-10;
  A3a density: |mean(realised fraction / Pi^2) - 1| + 2 SE < 0.01 over the bins of (M_A 8, M_B 8);
  A3b intervals: |mean(count / proper area) - 1| + 2 SE < 0.03;
  null: in the kappa = 0 runs, |mean B displacement| + 2 SE < 0.1 for every combination.
Level 1 (growth, PASSED requires): mean(count / coordinate area) + 2 SE < 0.9 (the counted proper times change).
Level 2 (active response): r_B for M_A = 4, 8, 16: PASSED |r - 1| + 2 SE < 0.15 each; FAILED |r - 1| - 2 SE > 0.15 any.
Level 3 (universality): |u| + 2 SE < 0.06 PASSED; |u| - 2 SE > 0.06 FAILED.
A2 (covariance): the boosted launch, r_sep = separation change / continuum: |r_sep - 1| + 2 SE < 0.15 PASSED;
  |r_sep - 1| - 2 SE > 0.15 FAILED.
PASSED if all PASSED conditions hold; FAILED if any FAILED condition; UNDECIDED otherwise.
Reported: r_A for every combination (A towards B), the third law M_A dA / (M_B dB), L along x at the end against the
continuum, the realised fraction.

Reading declared now
--------------------
PASSED: a causal growth rule of the order, using only the realised past, counts and independent marks, makes matter
thin its future order (slower counted clocks), and a test defect falls towards the source in proportion to the source's
energy and independently of its own mass, as the continuum response box_g L = -kappa q predicts, also for a moving
system. With F3.14b and F5.2a3, inertial, passive and active mass come from one field and one rule (1+1). A candidate for
the Axiom of the dynamics of the order. Limits: 1+1 (scalar response R2, no 1/r, the potential of a static source grows
in space and time); the causal relation does not change (conformal flatness; in 2+1 it must); kappa is chosen (G is not
derived); the potential is chosen.
"""

from __future__ import annotations

import math
import time

import numpy as np

from core.order.sprinkle import derive_seed
from core.stats.parallel import parallel_map
from core.testing import Result
from phases.F3_structure import field1d as F
from phases.F5_gravity import growth1d as G

TEST = {
    "id": "F5.2b",
    "phase": "F5",
    "title": "Active gravitational mass from a causal growth rule of the order (1+1, scalar response)",
    "theory_part": "F5 gravity; the Axiom of the dynamics of the order (candidate growth rule); H5 of "
                   "docs/cosmogenesis.md (active half); decision 15 (the chain field -> load -> response)",
    "question": "With candidates realised with probability Pi^2 computed from the realised causal past (ln Pi = -(kappa/2 rho) "
                "sum of the load 2V), does matter thin its future order, and does a test defect fall towards it in "
                "proportion to its energy and independently of its own mass, as the continuum response predicts?",
    "method": "1+1, rho 800, kappa 0.003, two separate fields (A, B) on the realised order, kinks at -3 and +3; M_A 4, 8, 16; "
              "M_B 4, 8; boosted launch v 0.3; paired kappa = 0; continuum references; 256 sprinklings",
    "prediction": "r = 1 for every M_A; u = 0; r_sep = 1; realised density = rho Pi^2; counts measure proper areas",
    "criterion": "INSTRUMENT: A0, A1 (1e-10), density |f/Pi^2 - 1| + 2 SE < 0.01, intervals |n/proper - 1| + 2 SE < 0.03, "
                 "null |dB| + 2 SE < 0.1. PASSED: n/coordinate + 2 SE < 0.9; |r_B - 1| + 2 SE < 0.15 (M_A 4, 8, 16); "
                 "|u| + 2 SE < 0.06; |r_sep - 1| + 2 SE < 0.15. FAILED: any |.| - 2 SE above its tolerance. Otherwise "
                 "UNDECIDED (card docs/F5_2b_card_draft.md)",
    "null_model": "kappa = 0 on the same candidates and marks (every candidate realised); continuum equations",
    "level": "mechanism",
    "meters": [],
    "depends_on": ["F3.14b", "F5.2a3"],
    "params": {"rho": 800.0, "kappa": 0.003, "d": 6.0, "T": 20.0, "n_seeds": 256, "v_boost": 0.3,
               "combos": [[0.25, 1.0], [1.0, 1.0], [4.0, 1.0], [1.0, 0.25]], "n_check": 4, "n_int": 8, "n_pairs": 60,
               "tol_dens": 0.01, "tol_int": 0.03, "tol_null": 0.1, "max_coord": 0.9, "tol_r": 0.15, "tol_u": 0.06,
               "tol_sep": 0.15, "boot_B": 1000, "n_se": 2.0},
    "quick": {"rho": 200.0, "n_seeds": 3, "T": 14.0, "n_check": 1, "n_int": 1, "n_pairs": 20, "boot_B": 100},
    "estimate": "about 2 h with 6 processes (numba required): 256 x 5 tasks of two evolutions of ~0.8 million candidates",
}


def task(p, k, muA, muB, v, base):
    rho, kappa, d, T = float(p["rho"]), float(p["kappa"]), float(p["d"]), float(p["T"])
    X = T + 16.0
    seed = derive_seed(base, "F5.2b", k)
    t, x = F.sprinkle_trapezoid(rho, X, T, seed)
    marks = np.random.default_rng(derive_seed(seed, "marks")).uniform(0.0, 1.0, t.size)
    xA, xB = -0.5 * d, 0.5 * d
    out = {"k": k, "muA": muA, "muB": muB, "v": v, "N": int(t.size)}
    keep = None
    for kap in (kappa, 0.0):
        phiA, phiB, L, real = G.evolve_growth(t, x, marks, muA, muB, xA, xB, rho, kap, v=v)
        tag = "k" if kap > 0 else "0"
        out[f"xA_{tag}"] = G.track(t, x, phiA, real, muA, xA, T)
        out[f"xB_{tag}"] = G.track(t, x, phiB, real, muB, xB, T)
        if kap > 0:
            keep = (phiA, phiB, L, real)
    phiA, phiB, L, real = keep
    if v == 0.0:
        bins = np.linspace(-12, 12, 13)
        s = t > (2.0 * T / 3.0)
        idx = np.digitize(x[s], bins) - 1
        out["dens"] = [float(real[s][idx == j].mean() / np.exp(2 * L[s][idx == j]).mean()) for j in range(12)
                       if np.any(idx == j)]
        s2 = t > T - 1.0
        out["L_end"] = [float(L[s2][np.abs(x[s2] - xg) < 0.5].mean()) for xg in np.linspace(-15, 15, 31)]
    if v == 0.0 and muA == 1.0 and muB == 1.0 and k < int(p["n_int"]):
        rng = np.random.default_rng(derive_seed(seed, "pairs"))
        U, V = (t + x) / math.sqrt(2.0), (t - x) / math.sqrt(2.0)
        ri = np.nonzero(real & (np.abs(x - xA) < 4.0) & (t > T / 2))[0]
        pairs = []
        tries = 0
        while len(pairs) < int(p["n_pairs"]) and ri.size > 10 and tries < 100000:
            tries += 1
            a, b = rng.choice(ri, 2, replace=False)
            if U[a] < U[b] and V[a] < V[b] and (U[b] - U[a]) * (V[b] - V[a]) > 1.0:
                pairs.append((a, b))
        pairs = np.array(pairs, np.int64)
        n = G._interval_counts(U, V, real, pairs)
        rp, rc = [], []
        for (a, b), c in zip(pairs, n):
            inside = (U > U[a]) & (U < U[b]) & (V > V[a]) & (V < V[b])
            vol = (U[b] - U[a]) * (V[b] - V[a])
            rp.append(c / (rho * vol * float(np.exp(2 * L[inside]).mean())))
            rc.append(c / (rho * vol))
        out["int_proper"], out["int_coord"] = rp, rc
    if v == 0.0 and muA == 1.0 and muB == 1.0 and k < int(p["n_check"]):
        a2, b2, L2, r2 = G.evolve_growth(t, x, marks, muA, muB, xA, xB, rho, kappa, order="boost")
        out["A1"] = bool(np.array_equal(r2, real)) and float(np.nanmax(np.abs(np.nan_to_num(a2 - phiA)) +
                                                                   np.abs(np.nan_to_num(b2 - phiB)))) < 1e-10
        m2 = marks.copy()
        late = t > np.quantile(t, 0.9)
        m2[late] = np.random.default_rng(derive_seed(seed, "late")).uniform(0, 1, late.sum())
        a3, b3, L3, r3 = G.evolve_growth(t, x, m2, muA, muB, xA, xB, rho, kappa)
        e = ~late
        out["A0"] = bool(np.array_equal(r3[e], real[e]) and np.array_equal(np.nan_to_num(a3[e]), np.nan_to_num(phiA[e]))
                         and np.array_equal(np.nan_to_num(b3[e]), np.nan_to_num(phiB[e])))
    return out


def run_task(*a):
    return task(*a)


def _disp(xc):
    return xc[-1] - xc[0]


def run(p, ctx):
    if G._nb is None:
        raise SystemExit("F5.2b needs numba (pip install numba)")
    t0 = time.time()
    n, nse, B = int(p["n_seeds"]), float(p["n_se"]), int(p["boot_B"])
    base = ctx.seed
    combos = [(float(a), float(b), 0.0) for a, b in p["combos"]] + [(1.0, 1.0, float(p["v_boost"]))]
    jobs = [(p, k, a, b, v, base) for k in range(n) for (a, b, v) in combos]
    res = parallel_map(run_task, jobs, ctx.workers, ctx.progress,
                       lambda i, r, d, nn: f"[F5.2b] {d}/{nn}: muA {r['muA']:g} muB {r['muB']:g} v {r['v']:g} "
                                           f"({time.time() - t0:.0f}s)")
    R = {(r["k"], r["muA"], r["muB"], r["v"]): r for r in res}
    T, d, kappa = float(p["T"]), float(p["d"]), float(p["kappa"])
    rng = np.random.default_rng(derive_seed(base, "F5.2b-boot"))
    boot_idx = [rng.integers(0, n, n) for _ in range(B)]

    def ratio(vals, ref):
        vals = np.asarray(vals, float)
        if not np.all(np.isfinite(vals)):
            return float("nan"), float("nan")
        bs = [vals[i].mean() / ref for i in boot_idx]
        return float(vals.mean() / ref), float(np.std(bs, ddof=1))

    def ms(vals):
        vals = np.asarray(vals, float)
        if not np.all(np.isfinite(vals)):
            return float("nan"), float("nan")
        return float(vals.mean()), float(vals.std(ddof=1) / math.sqrt(vals.size))

    lines = [f"1+1 growth rule, rho {p['rho']:g}, kappa {kappa:g}, kinks at {-d / 2:g} and {d / 2:g}, T {T:g}, {n} "
             f"sprinklings, paired with kappa = 0"]
    inst, passc, failc = {}, {}, {}
    chk = [r for r in res if "A1" in r]
    inst["A0 earlier events unchanged"] = all(r["A0"] for r in chk)
    inst["A1 order of computation"] = all(r["A1"] for r in chk)
    lines.append(f"  A0 {inst['A0 earlier events unchanged']}, A1 {inst['A1 order of computation']} ({len(chk)} sprinklings)")
    dens = [np.mean(R[(k, 1.0, 1.0, 0.0)]["dens"]) for k in range(n)]
    md = ms(dens)
    inst["A3a density"] = abs(md[0] - 1.0) + nse * md[1] < float(p["tol_dens"])
    ip = [np.mean(r["int_proper"]) for r in res if "int_proper" in r]
    ic = [np.mean(r["int_coord"]) for r in res if "int_coord" in r]
    mp, mc = ms(ip), ms(ic)
    inst["A3b intervals (proper)"] = abs(mp[0] - 1.0) + nse * mp[1] < float(p["tol_int"])
    passc["level 1: counted proper times change"] = mc[0] + nse * mc[1] < float(p["max_coord"])
    lines.append(f"  A3 density: realised / Pi^2 {md[0]:.4f} +- {md[1]:.4f}; intervals near A: count / proper area "
                 f"{mp[0]:.4f} +- {mp[1]:.4f}, count / coordinate area {mc[0]:.4f} +- {mc[1]:.4f}")
    out = {}
    for (a, b, v) in combos:
        cA, cB, Lc = G.continuum(a, b, -d / 2, d / 2, kappa, T, v=v)
        c0A, c0B, _ = G.continuum(a, b, -d / 2, d / 2, 0.0, T, v=v)
        rr = [R[(k, a, b, v)] for k in range(n)]
        dB = [_disp(r["xB_k"]) - _disp(r["xB_0"]) for r in rr]
        dA = [_disp(r["xA_k"]) - _disp(r["xA_0"]) for r in rr]
        nullB = ms([_disp(r["xB_0"]) - v * (T - 1.0) for r in rr])
        dBc = (cB[-1] - cB[0]) - (c0B[-1] - c0B[0])
        dAc = (cA[-1] - cA[0]) - (c0A[-1] - c0A[0])
        MA, MB = 8 * math.sqrt(a), 8 * math.sqrt(b)
        tag = f"M_A {MA:g}, M_B {MB:g}" + (f", launched v {v:g}" if v else "")
        inst[f"null {tag}"] = abs(nullB[0]) + nse * nullB[1] < float(p["tol_null"])
        rB, rA = ratio(dB, dBc), ratio(dA, dAc)
        out[(a, b, v)] = {"dB": dB, "dBc": dBc, "rB": rB, "rA": rA}
        third = (MA * np.mean(dA)) / (MB * -np.mean(dB)) if np.mean(dB) != 0 else float("nan")
        line = (f"  {tag}: B towards A {np.mean(dB):+.4f} (continuum {dBc:+.4f}), r_B {rB[0]:.4f} +- {rB[1]:.4f}; "
                f"A towards B {np.mean(dA):+.4f} (continuum {dAc:+.4f}), r_A {rA[0]:.4f} +- {rA[1]:.4f}; third law "
                f"M_A dA / (M_B |dB|) {third:.3f}; null B {nullB[0]:+.4f} +- {nullB[1]:.4f}")
        if v:
            sep = [(_disp(r["xB_k"]) - _disp(r["xA_k"])) - (_disp(r["xB_0"]) - _disp(r["xA_0"])) for r in rr]
            sepc = ((cB[-1] - cB[0]) - (cA[-1] - cA[0])) - ((c0B[-1] - c0B[0]) - (c0A[-1] - c0A[0]))
            rs = ratio(sep, sepc)
            line += f"; separation change {np.mean(sep):+.4f} (continuum {sepc:+.4f}), r_sep {rs[0]:.4f} +- {rs[1]:.4f}"
            passc["A2 covariance (launched)"] = abs(rs[0] - 1.0) + nse * rs[1] < float(p["tol_sep"])
            failc["A2 covariance (launched)"] = abs(rs[0] - 1.0) - nse * rs[1] > float(p["tol_sep"])
        elif b == 1.0:
            passc[f"level 2: r_B M_A {MA:g}"] = abs(rB[0] - 1.0) + nse * rB[1] < float(p["tol_r"])
            failc[f"level 2: r_B M_A {MA:g}"] = abs(rB[0] - 1.0) - nse * rB[1] > float(p["tol_r"])
        if not v:
            Le = np.mean([R[(k, a, b, v)]["L_end"] for k in range(n)], axis=0)
            line += "\n    L at the end x = -15, -9, -3, 3, 9, 15 (discrete / continuum): " + " ".join(
                f"{Le[j]:+.3f}/{Lc[j]:+.3f}" for j in (0, 6, 12, 18, 24, 30))
        lines.append(line)
    # universality, paired by seed
    b4 = np.array(out[(1.0, 0.25, 0.0)]["dB"]) / out[(1.0, 0.25, 0.0)]["dBc"]
    b8 = np.array(out[(1.0, 1.0, 0.0)]["dB"]) / out[(1.0, 1.0, 0.0)]["dBc"]
    u = float(b4.mean() - b8.mean())
    use = float(np.std([b4[i].mean() - b8[i].mean() for i in boot_idx], ddof=1))
    lines.append(f"  level 3 universality: u = r_B(M_B 4) - r_B(M_B 8) at M_A 8 = {u:+.4f} +- {use:.4f}")
    passc["level 3: universality"] = abs(u) + nse * use < float(p["tol_u"])
    failc["level 3: universality"] = abs(u) - nse * use > float(p["tol_u"])
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
    return Result(status, lines, {"instrument": inst, "passed": passc, "failed": failc, "results": res})
