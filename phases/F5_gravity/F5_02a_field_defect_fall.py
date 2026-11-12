"""
F5.2a -- Free fall of a field defect in a 1+1 order of nonuniform density (passive gravitational mass, WEP)
=========================================================================================================

Why this test exists (10 Nov 2026; card docs/F5_2a_card_draft.md, revision 2, frozen by the author on 10 Nov, 13:53)
------------------------------------------------------------------------------------------------------------------
F3.14a: a defect of a field on all events of a 1+1 order, evolved by the retarded Green function of the order, has
inertia with M_in = M_rest (seed 0 PASSED; replication INSTRUMENT_FAILED on the drift tolerance alone, all physical
conditions met twice). The gravitational half of H5 (docs/cosmogenesis.md) starts here with the passive mass; the
active mass is F5.2b (a back-reaction of the field on the order, 3+1). The planned F5.1i is carried over to this series.

Every 1+1 metric is conformally flat, ds^2 = Omega^2 (dt^2 - dx^2): the light cones and the causal order are those of
flat spacetime and only the number of events per coordinate area changes (rho Omega^2 for proper density rho). The rule
of F3.14a is NOT changed:
    phi(x) = phi_0(x) - (1 / 2 rho) * sum over events y < x (the order), t_y >= 0, of mu sin(phi(y));
on such an order it is the field equation in the curved spacetime (the massless operator is conformally invariant in
1+1). The rule knows nothing about the metric: a kink that falls falls because of the counting. In 1+1 the density and
the clock rate are the same function (Omega^2, Omega): the kink falls towards the sparser order; this sign is specific
to 1+1 and is not a test of the throughput principle of 3+1.

Model (frozen)
--------------
Phi(x) = g L tanh(x / L), g = 0.05, L = 10, Omega = exp(Phi) (curved: R ~ Phi'' != 0). Sprinkling of coordinate density
rho exp(2 Phi) in the trapezoid |x| <= 26 - t, T = 10 (domain of dependence of the initial slice). Kinks mu = 0.25, 1, 4
(M = 4, 8, 16, proper widths 2, 1, 0.5) at rest at x0 = 0, profile in proper distance s = int Omega dx, launched through
the free wave phi_0 (d'Alembert). Arms per seed: curved mu = 0.25, 1, 4 on the SAME curved sprinkling; flat (g = 0) the
same masses on the flat sprinkling of the same seed; counting null (mu = 1, curved sprinkling, every term divided by
Omega^2(y): the sum measures coordinate area, not the number of events; coordinates; a null model); up (mu = 1 launched
with coordinate velocity 0.3 against the field; reported). Readout: the kink fit of F3.14a per slab of width 1 (centres
0.5 ... 9.5), tracking the centre; D = x_c(9.5) - x_c(0.5). References: D_cont of the continuum field equation
phi_tt - phi_xx + Omega^2 mu sin phi = 0 (finite differences, the same initial data and centre readout) for each body;
the geodesic d2x/dt2 = -Phi'(1 - (dx/dt)^2) (reported). Helpers: phases/F5_gravity/field1d_curved.py.

Criteria (frozen with the card; deciding density 800, 128 new sprinklings, seed label "F5.2a"; SE = sd / sqrt(n), paired
where stated; every condition holds on its own; nothing averaged into a score, nothing chosen after the run)
-----------------------------------------------------------------------------------------------------------------------
INSTRUMENT_FAILED if any fails:
  I4 on the curved order: order (t vs t + 0.5 x), engine (Fenwick vs blocks), region (X vs X - 2): each < 1e-12;
  flat control, each mass: |mean D| + 2 SE < 0.2;
  counting null: |mean D| + 2 SE < 0.2;
  the fall is resolved, each mass: |mean D| - 2 SE > 1.0.
Decision: a_i = D_i / D_cont,i (each mass); eta_ij = a_i / a_j - 1, paired, pairs (4, 8), (8, 16), (4, 16).
  PASSED: |a_i - 1| + 2 SE < 0.05 for every mass and |eta_ij| + 2 SE < 0.06 for every pair.
  FAILED: |a_i - 1| - 2 SE > 0.05 or |eta_ij| - 2 SE > 0.06 for any.
  UNDECIDED otherwise.
Reported: raw eta (D_i / D_j - 1), D / D_geo, the launched arm against its continuum, rho = 400 (32 sprinklings), the
spread of D per mass (the wandering seen in F3.14a).

Reading declared now (card, section 6)
-------------------------------------
PASSED: in an order whose only change is the number of events, a bound defect falls as its continuum body in the
geometry the number defines, and defects of mass 4, 8 and 16 fall alike within 6 %: the weak equivalence principle
from counting, at the resolving power of this experiment (not a claim at the level of MICROSCOPE, eta ~ 1e-15); with
F3.14a, inertia and passive gravitational mass from the same rule. FAILED (fall): the counting does not give the
continuum fall at this scale. FAILED (WEP): defects of different mass fall differently by discreteness. Limits: 1+1 with
the opposite sign of density and clock rate to 3+1; the geometry is prescribed (no back-reaction: F5.2b); the potential
is chosen (H5-D); the launch, the readout and the counting null use coordinates (declared).
"""

from __future__ import annotations

import math
import time

import numpy as np

from core.order.sprinkle import derive_seed
from core.stats.parallel import parallel_map
from core.testing import Result
from phases.F3_structure import field1d as F
from phases.F5_gravity import field1d_curved as C

TEST = {
    "id": "F5.2a",
    "phase": "F5",
    "title": "Free fall of a field defect in a 1+1 order of nonuniform density (passive gravitational mass, WEP from "
             "counting)",
    "theory_part": "F5 gravity: equivalence principle; matter as a defect of a field of the order (H5 of "
                   "docs/cosmogenesis.md, gravitational half); order + number = geometry",
    "question": "In a 1+1 causal order whose only change is the number of events per coordinate area (a conformally flat "
                "metric), does the kink of F3.14a, under the unchanged rule, fall as its continuum body, and do kinks of "
                "mass 4, 8 and 16 fall alike, while it does not fall when the density is uniform or the counting is "
                "removed?",
    "method": "sprinkling of density rho exp(2 Phi), Phi = g L tanh(x/L), g 0.05, L 10; kinks mu 0.25, 1, 4 at rest on the "
              "same sprinkling; flat control; counting null (terms divided by Omega^2); continuum references by finite "
              "differences; 128 sprinklings at rho 800",
    "prediction": "a_i = D_i / D_cont,i = 1 for every mass; eta_ij = 0; flat and null do not fall",
    "criterion": "INSTRUMENT: I4 < 1e-12; flat and null |D| + 2 SE < 0.2; |D| - 2 SE > 1 each mass. PASSED: |a_i - 1| "
                 "+ 2 SE < 0.05 each mass and |eta_ij| + 2 SE < 0.06 each pair. FAILED: |a_i - 1| - 2 SE > 0.05 or "
                 "|eta_ij| - 2 SE > 0.06. Otherwise UNDECIDED (card docs/F5_2a_card_draft.md rev. 2, section 5b)",
    "null_model": "uniform density (flat); the same curved order with the counting removed (terms divided by Omega^2)",
    "level": "mechanism",
    "meters": [],
    "depends_on": ["F3.14a"],
    "params": {"rho_dec": 800.0, "rho_rep": 400.0, "n_seeds": 128, "n_rep": 32, "g": 0.05, "L": 10.0, "X": 26.0,
               "T": 10.0, "mus": [0.25, 1.0, 4.0], "v_up": 0.3, "tol_flat": 0.2, "tol_res": 1.0, "tol_a": 0.05,
               "tol_eta": 0.06, "i4_n": 2, "i4_rho": 100.0, "i4_tol": 1e-12, "n_se": 2.0},
    "quick": {"rho_dec": 100.0, "rho_rep": 50.0, "n_seeds": 2, "n_rep": 2, "i4_n": 1},
    "estimate": "20-30 min with 8 processes (numba required): ~400 000 events per arm at rho 800, 8 arms x 128 "
                "sprinklings, plus 32 at rho 400",
}


def task(p, rho, k, arm, base):
    X, T = float(p["X"]), float(p["T"])
    kind, mu = arm.split("_")
    mu = float(mu)
    g = 0.0 if kind == "flat" else float(p["g"])
    seed = derive_seed(base, "F5.2a", rho, k)
    t, x = C.sprinkle_curved(rho, X, T, g, seed)
    v0 = float(p["v_up"]) if kind == "up" else 0.0
    phi0 = C.free_wave(t, x, mu, g, v0)
    t0 = time.time()
    if kind == "null":
        U, V = (t + x) / math.sqrt(2.0), (t - x) / math.sqrt(2.0)
        phi = C._fenwick_weighted(U, V, np.arange(t.size, dtype=np.int64), phi0, np.exp(-2.0 * C.Phi(x, g)),
                                  float(rho), mu)
    else:
        phi = F.evolve(t, x, phi0, rho, mu=mu)
    xc = C.track(t, x, phi, mu, C.X0, T)
    return {"kind": "fall", "rho": rho, "k": k, "arm": arm, "N": int(t.size), "xc": xc, "sec": time.time() - t0}


def task_i4(p, k, base):
    rho, X, T, g = float(p["i4_rho"]), 12.0, 6.0, float(p["g"])
    t, x = C.sprinkle_curved(rho, X, T, g, derive_seed(base, "F5.2a-i4", k))
    phi0 = C.free_wave(t, x, 1.0, g)
    ref = F.evolve(t, x, phi0, rho)
    sc = float(np.max(np.abs(ref)))
    alt = F.evolve(t, x, phi0, rho, order="boost")
    blk = F.evolve(t, x, phi0, rho, engine="block", block=256)
    sub = np.abs(x) <= (X - 2.0) - t
    phs = F.evolve(t[sub], x[sub], phi0[sub], rho)
    return {"kind": "i4", "k": k, "N": int(t.size), "order": float(np.max(np.abs(alt - ref)) / sc),
            "engine": float(np.max(np.abs(blk - ref)) / sc), "region": float(np.max(np.abs(phs - ref[sub])) / sc)}


def run_task(kind, *a):
    return task(*a) if kind == "fall" else task_i4(*a)


def _ms(a):
    a = np.asarray(a, float)
    if a.size < 2 or not np.all(np.isfinite(a)):
        return float("nan"), float("nan"), float("nan")
    return float(a.mean()), float(a.std(ddof=1) / math.sqrt(a.size)), float(a.std(ddof=1))


def run(p, ctx):
    if F._nb is None:
        raise SystemExit("F5.2a needs numba (pip install numba)")
    if abs(float(p["g"]) - C.G) > 1e-15 or abs(float(p["L"]) - C.L) > 1e-15:
        raise SystemExit("g and L are frozen in phases/F5_gravity/field1d_curved.py")
    t0 = time.time()
    nse = float(p["n_se"])
    mus = [float(m) for m in p["mus"]]
    arms = [f"curved_{m:g}" for m in mus] + [f"flat_{m:g}" for m in mus] + ["null_1", f"up_1"]
    rd, rr = float(p["rho_dec"]), float(p["rho_rep"])
    base = ctx.seed
    jobs = [("fall", p, rd, k, arm, base) for k in range(int(p["n_seeds"])) for arm in arms]
    jobs += [("fall", p, rr, k, arm, base) for k in range(int(p["n_rep"])) for arm in arms]
    jobs += [("i4", p, k, base) for k in range(int(p["i4_n"]))]
    res = parallel_map(run_task, jobs, ctx.workers, ctx.progress,
                       lambda i, r, d, n: f"[F5.2a] {d}/{n}: {r['kind']} {r.get('arm', '')} rho {r.get('rho', '')} "
                                          f"N {r['N']} ({time.time() - t0:.0f}s)")
    tt, geo = C.geodesic()
    tm = np.arange(0.5, float(p["T"]), 1.0)
    gi = np.interp(tm, tt, geo)
    D_geo = float(gi[-1] - gi[0])
    cont = {m: C.continuum(m) for m in mus}
    D_cont = {m: float(cont[m][-1] - cont[m][0]) for m in mus}
    cu = C.continuum(1.0, v0=float(p["v_up"]))
    D_cont_up = float(cu[-1] - cu[0])
    lines = [f"1+1, Phi = g L tanh(x/L), g {p['g']}, L {p['L']}; trapezoid |x| <= {p['X']:g} - t, T {p['T']:g}; deciding "
             f"density {rd:g} with {p['n_seeds']} sprinklings; reported {rr:g} with {p['n_rep']}",
             f"references: D_geo {D_geo:+.4f}; D_cont " + ", ".join(f"mu {m:g}: {D_cont[m]:+.4f}" for m in mus) +
             f"; launched (mu 1, v 0.3) D_cont {D_cont_up:+.4f}"]
    inst, passc, failc = {}, {}, {}
    i4 = [r for r in res if r["kind"] == "i4"]
    for key in ("order", "engine", "region"):
        mx = max(r[key] for r in i4)
        inst[f"I4 {key}"] = mx < float(p["i4_tol"])
        lines.append(f"  I4 {key} (curved order): {mx:.2e}")
    for rho, n in ((rd, int(p["n_seeds"])), (rr, int(p["n_rep"]))):
        R = {(r["k"], r["arm"]): r for r in res if r["kind"] == "fall" and r["rho"] == rho}
        D = {arm: np.array([R[(k, arm)]["xc"][-1] - R[(k, arm)]["xc"][0] for k in range(n)], float) for arm in arms}
        tag = "DECIDING" if rho == rd else "reported"
        lines.append(f"--- rho {rho:g} ({tag})")
        a = {m: D[f"curved_{m:g}"] / D_cont[m] for m in mus}
        for m in mus:
            md = _ms(D[f"curved_{m:g}"])
            ma = _ms(a[m])
            mf = _ms(D[f"flat_{m:g}"])
            lines.append(f"  M {8 * math.sqrt(m):g} (mu {m:g}): D {md[0]:+.4f} +- {md[1]:.4f} (sd {md[2]:.3f}); a = D / D_cont "
                         f"{ma[0]:.4f} +- {ma[1]:.4f}; D / D_geo {md[0] / D_geo:.4f}; flat D {mf[0]:+.4f} +- {mf[1]:.4f} "
                         f"(sd {mf[2]:.3f})")
            if rho == rd:
                inst[f"flat M {8 * math.sqrt(m):g}"] = abs(mf[0]) + nse * mf[1] < float(p["tol_flat"])
                inst[f"fall resolved M {8 * math.sqrt(m):g}"] = abs(md[0]) - nse * md[1] > float(p["tol_res"])
                passc[f"fall M {8 * math.sqrt(m):g}"] = abs(ma[0] - 1.0) + nse * ma[1] < float(p["tol_a"])
                failc[f"fall M {8 * math.sqrt(m):g}"] = abs(ma[0] - 1.0) - nse * ma[1] > float(p["tol_a"])
        for i, j in ((0, 1), (1, 2), (0, 2)):
            mi, mj = mus[i], mus[j]
            eta = _ms(a[mi] / a[mj] - 1.0)
            raw = _ms(D[f"curved_{mi:g}"] / D[f"curved_{mj:g}"] - 1.0)
            name = f"WEP M {8 * math.sqrt(mi):g} vs {8 * math.sqrt(mj):g}"
            lines.append(f"  {name}: eta {eta[0]:+.4f} +- {eta[1]:.4f} (sd {eta[2]:.3f}); raw D_i / D_j - 1 "
                         f"{raw[0]:+.4f} +- {raw[1]:.4f}")
            if rho == rd:
                passc[name] = abs(eta[0]) + nse * eta[1] < float(p["tol_eta"])
                failc[name] = abs(eta[0]) - nse * eta[1] > float(p["tol_eta"])
        mn = _ms(D["null_1"])
        mu_ = _ms(D["up_1"])
        lines.append(f"  counting null (M 8): D {mn[0]:+.4f} +- {mn[1]:.4f} (curved {np.mean(D['curved_1']):+.4f})")
        lines.append(f"  launched v 0.3 against the field (reported): D {mu_[0]:+.4f} +- {mu_[1]:.4f}; continuum "
                     f"{D_cont_up:+.4f}; ratio {mu_[0] / D_cont_up:.4f}")
        if rho == rd:
            inst["counting null"] = abs(mn[0]) + nse * mn[1] < float(p["tol_flat"])
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
    return Result(status, lines, {"instrument": inst, "passed": passc, "failed": failc, "D_cont": D_cont,
                                  "D_geo": D_geo, "results": res})
