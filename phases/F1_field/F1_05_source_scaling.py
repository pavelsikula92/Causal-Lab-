"""
F1.5 -- Does the field scale with the source? (self-consistent field equation (a))
==================================================================================

Why this test is not trivial
----------------------------
With a fixed operator on a fixed (flat) causal set, A5' is linear: twice
the source gives exactly twice the field, and a test would only check
linear algebra. The author decided on 18 Oct 2026 that the field obeys

    (a)   B_g Pi = kappa rho_s,

B_g the operator of the causal set sprinkled into the geometry g = Pi^2 eta
that the field itself creates. Then the equation is nonlinear, and two
things become testable:

  P1 (linear limit): for a weak source the field is proportional to the
     source, with a coupling set by the curvature term. The local limit of
     B_g is box_g - R/2 (calibrated in K.11 local); for g = Pi^2 eta its
     linearization is D lap u (u = Pi - 1), D times the flat operator. So
     the weak field of (a) is the flat field divided by D: the -R/2 term
     changes the effective gravitational coupling to kappa / D.
     Alternatives: box_g alone or no back-reaction both give coupling 1.
  P2 (nonlinearity): for stronger sources the far-field strength grows
     more slowly than the source. In the local limit (a) reads
     lap Pi^p = (p / D) kappa Pi^(1+p) rho_s (p = 1 in 2+1), so the source
     is weighted by Pi^(1+p) < 1 where the field is deep. The predicted
     curve s(q) = m(q) / (q m_lin) (far-field strength per unit source,
     relative to the weak-field limit) comes from the continuum equation
     of (a) with the same source and boundary. Alternative: strict
     linearity, s = 1.

What is computed
----------------
Lengths in units of l_k (eps rho = 1). Gaussian source q exp(-r^2/2 sigma^2),
sigma = 8 l_k (wide compared with l_k, so that the nonlocal corrections of
the operator stay small: K.11 local gives 0.89 of the local limit at width
8 in 2+1, 0.95 at width 16). Cylinder radius R = 100, rows at t_x = 3R.

  * flat:   the flat expected operator (no back-reaction), linear: A_flat
            per unit source;
  * curved: the self-consistent solution of E[B_g(Pi)] Pi = q rho_hat
            (core.operators.bd_expect_curved.solve_selfconsistent, with the
            flat operator on the vacuum Pi = 1 subtracted, the difference
            form of F1.3) at a weak source q_lin = q_min / 100 and at q_list;
  * far-field strength: A from the least-squares fit u = A G(r) + B in the
            window [r_lo, r_hi] (G = ln r in 2+1, 1/r in 3+1), the SAME
            fit for the causal set and for the continuum candidates.
  * k = A_curved(q_lin) / (q_lin A_flat)        (coupling, P1)
  * s(q) = (A(q) / q) / (A(q_lin) / q_lin)        (nonlinearity, P2)

Quadrature: the coarsest level that K.11b far_d3 found within its
tolerance (volume-weighted far residue within 1 % of the source charge of
the finest level): [240, 8, 96, 8, 48, 1]; it also passes K.11 conv and
diamond. Set before the first run of F1.5 (18 Oct, after K.11b; F1.4 had
run at the coarser level [160, 8, 48, 6, 24, 1], which K.11b puts 1.3 %
outside that tolerance). To keep the cost of the finer quadrature, the
node step is 2 l_k (sigma / 4) with geometric ratio 1.12 beyond 4 sigma.
3+1: K.11b far_d4 FAILED (the far residue did not converge: Q_far/Q_src
3.5, 2.3, 1.9 over the three levels), so there is no 3+1 part.

Criteria (stated before the run; the quarter rule of F1.3)
----------------------------------------------------------
  P1: |k - 1/D| <= |1 - 1/D| / 4   (three times closer to (a) than to the
      alternatives with coupling 1);
  P2: RMS_q(s_cs - s_a) <= RMS_q(1 - s_a) / 4   (three times closer to the
      nonlinearity of (a) than to strict linearity).
Instrument: every self-consistent solve converged (largest change of Pi
in the last step < tol) and its residual RMS(E[B_g] Pi - rhs) <=
res_tol * max(rhs) (res_tol of F1.2b/F1.3). If not: INSTRUMENT_FAILED.
PASSED: P1 and P2. FAILED: P1 or P2 does not hold.

Strengths: chosen from the continuum equation of (a) before any causal-set
run, so that the central Pi spans ~0.9 to ~0.6 (2+1: q = 0.0025 ... 0.02
gives Pi(0) = 0.89 ... 0.58). Reported only: box_g alone where a static
solution exists, the central Pi of every solution, iterations.

Null model: no sampling enters. The competing hypotheses are the
alternative couplings and strict linearity; the criteria separate them.
"""

from __future__ import annotations

import time

import numpy as np

from core.operators.bd_expect import shape_nodes
from core.operators.bd_expect_curved import continuum_radial, solve_selfconsistent
from core.stats.parallel import parallel_map
from core.testing import Result

TEST = {
    "id": "F1.5",
    "phase": "F1",
    "title": "Scaling of the field with the strength of the source in the self-consistent equation (a)",
    "theory_part": "A5′ as the self-consistent equation B_g Π = κρ_s (author's decision of 18 Oct 2026)",
    "question": "Is the field of a weak source proportional to the source with the coupling κ/D (effect of the term "
                "−R/2), and does the field strength of a stronger source grow more slowly according to the "
                "non-linearity of equation (a)?",
    "method": "Gaussian source σ = 8 l_k; flat operator (no feedback) and the self-consistent solution E[B_g(Π)] Π = "
              "q ρ̂ at a weak q and four stronger q; strength of the far field by the fit u = A G(r) + B; comparison "
              "with the continuum equation (a)",
    "prediction": "k = A_curved/A_flat → 1/D in a weak field; s(q) = (A/q)/(A_lin/q_lin) according to the continuum "
                  "equation (a)",
    "criterion": "INSTRUMENT: all solutions converged and the residual <= res_tol·max(rhs). P1: |k − 1/D| <= |1 − "
                 "1/D|/4. P2: RMS(s_cs − s_a) <= RMS(1 − s_a)/4. PASSED: P1 and P2. FAILED: P1 or P2 does not hold.",
    "null_model": "no sampling; competing hypotheses: coupling 1 (□_g without −R/2, or no feedback) and strict "
                  "linearity s = 1",
    "level": "mechanism",
    "meters": ["K.4", "K.11", "K.11b"],
    "depends_on": ["K.11", "K.11b", "F1.3"],
    "params": {"sigma": 8.0, "R": 100.0, "dr": 2.0, "ratio": 1.12, "rows_per_element": 2,
               "r_lo": 24.0, "r_hi": 50.0, "q_lin_factor": 1e-2, "tol": 1e-5, "max_iter": 30, "res_tol": 0.05,
               "chunk": 4},
    "parts": {"d3": {"D": 3, "quad": [240, 8, 96, 8, 48, 1], "q_list": [0.0025, 0.005, 0.01, 0.02]}},
    "quick": {"sigma": 3.0, "R": 30.0, "dr": 1.0, "ratio": 1.2, "rows_per_element": 2, "r_lo": 9.0, "r_hi": 15.0,
              "tol": 1e-4, "max_iter": 12},
    "quick_parts": {"d3": {"quad": [40, 6, 12, 3, 8, 1], "q_list": [0.01, 0.04]}},
    "estimate": "2+1: 5 self-consistent solutions of ~5–12 iterations; finer quadrature (K.11b) and ~50 rows: of the "
                "order of 3–6 h with 16 processes",
}

_EPS, _RHO = 0.05, 20.0          # eps * rho = 1: lengths in units of l_k


def setup(p):
    """Nodes, row radii and the unit source profile at the rows."""
    R, sig = float(p["R"]), float(p["sigma"])
    nodes = shape_nodes(R, 2.0 * sig, float(p["dr"]), float(p["ratio"]))
    per = int(p["rows_per_element"])
    fr = (np.arange(per) + 0.5) / per
    rows = np.sort((nodes[:-1, None] + fr[None, :] * np.diff(nodes)[:, None]).ravel())
    src = lambda r: np.exp(-np.asarray(r, float) ** 2 / (2.0 * sig ** 2))
    return nodes, rows, src


def far_amplitude(r, u, window, D):
    """A of u = A G(r) + B in the window (G = ln r for D = 3, 1/r for D = 4)."""
    m = (r >= window[0]) & (r <= window[1])
    G = np.log(r[m]) if D == 3 else 1.0 / r[m]
    M = np.column_stack([G, np.ones(int(m.sum()))])
    coef, *_ = np.linalg.lstsq(M, u[m], rcond=None)
    return float(coef[0])


def run(p, ctx):
    D = int(p["D"])
    R = float(p["R"])
    nodes, rows, src = setup(p)
    window = [float(p["r_lo"]), float(p["r_hi"])]
    quad = [int(x) for x in p["quad"]]
    q_list = [float(x) for x in p["q_list"]]
    q_lin = float(p["q_lin_factor"]) * min(q_list)
    t0 = time.time()

    def solve(q, label, flat=False, Pi0=None):
        return solve_selfconsistent(D, R, nodes, rows, q * src(rows), quad, parallel_map, ctx.workers, ctx.progress,
                                    label=f"[F1.5] {D-1}+1 {label}", rho=_RHO, eps=_EPS, tol=float(p["tol"]),
                                    max_iter=int(p["max_iter"]), chunk=int(p["chunk"]), flat=flat, Pi0=Pi0)

    # flat operator, unit source (linear: the field per unit source)
    flat = solve(1.0, "flat operator", flat=True)
    A_flat = far_amplitude(nodes, flat["Pi"] - 1.0, window, D)
    sols = {"lin": solve(q_lin, f"curved q = {q_lin:g}")}
    prev_q, prev = q_lin, sols["lin"]["Pi"]
    for q in q_list:
        # continuation: start from the previous solution with u scaled to the new source (the fixed point is unchanged)
        sols[q] = solve(q, f"curved q = {q:g}", Pi0=1.0 + (prev - 1.0) * (q / prev_q))
        prev_q, prev = q, sols[q]["Pi"]
        ctx.progress(f"[F1.5] {D-1}+1 q = {q:g}: Pi(0) = {sols[q]['Pi'][0]:.4f}, {sols[q]['iterations']} iterations "
                     f"({time.time()-t0:.0f}s)")
    A = {k: far_amplitude(nodes, v["Pi"] - 1.0, window, D) for k, v in sols.items()}
    k_cs = (A["lin"] / q_lin) / A_flat
    s_cs = np.array([(A[q] / q) / (A["lin"] / q_lin) for q in q_list])

    # continuum candidates on a fine grid, read on the same nodes with the same fit
    def cont(model, q):
        c = continuum_radial(model, D, R, lambda r: q * src(r))
        return np.interp(nodes, c["r"], c["Pi"])
    Ac = {}
    for model in ("a", "flat", "box_g"):
        Ac[model] = {}
        for q in [q_lin] + q_list:
            try:
                Ac[model][q] = far_amplitude(nodes, cont(model, q) - 1.0, window, D)
            except ValueError:
                Ac[model][q] = float("nan")          # no static solution of that strength
    k_a = (Ac["a"][q_lin] / q_lin) / (Ac["flat"][q_lin] / q_lin)
    s_a = np.array([(Ac["a"][q] / q) / (Ac["a"][q_lin] / q_lin) for q in q_list])
    s_box = np.array([(Ac["box_g"][q] / q) / (Ac["box_g"][q_lin] / q_lin) for q in q_list])
    k_box = (Ac["box_g"][q_lin] / q_lin) / (Ac["flat"][q_lin] / q_lin)

    rhs_max = max(q_list)
    inst_conv = all(v["converged"] for v in sols.values())
    inst_res = all(v["residual"] <= float(p["res_tol"]) * max(q, q_lin) for q, v in
                   [(q_lin, sols["lin"])] + [(q, sols[q]) for q in q_list]) and \
        flat["residual"] <= float(p["res_tol"]) * 1.0
    inst = {"self-consistent solves converged": inst_conv, "residuals": inst_res}
    tolP1 = abs(1.0 - 1.0 / D) / 4.0
    P1 = abs(k_cs - 1.0 / D) <= tolP1
    dev_cs = float(np.sqrt(np.mean((s_cs - s_a) ** 2)))
    dev_lin = float(np.sqrt(np.mean((1.0 - s_a) ** 2)))
    P2 = dev_cs <= dev_lin / 4.0
    if not all(inst.values()):
        status = "INSTRUMENT_FAILED"
    elif P1 and P2:
        status = "PASSED"
    else:
        status = "FAILED"
    lines = [f"{D-1}+1 self-consistent field (a), lengths in l_k: Gaussian source sigma = {p['sigma']:g}, R = {R:g}, "
             f"{nodes.size} nodes, {rows.size} rows, quadrature {quad}; far-field window [{window[0]:g}, {window[1]:g}]",
             f"  flat operator: A_flat = {A_flat:.5f} per unit source (continuum flat {Ac['flat'][q_lin] / q_lin:.5f})",
             f"  P1 coupling k = A_curved(q_lin) / A_flat = {k_cs:.4f}; prediction 1/D = {1.0 / D:.4f} "
             f"(continuum (a) {k_a:.4f}, box_g {k_box:.4f}); |k - 1/D| = {abs(k_cs - 1.0 / D):.4f} (tol {tolP1:.4f}) -> "
             f"{'met' if P1 else 'NOT met'}",
             "  q            " + "  ".join(f"{q:9g}" for q in q_list),
             "  Pi(0) causal " + "  ".join(f"{sols[q]['Pi'][0]:9.4f}" for q in q_list),
             "  s causal set " + "  ".join(f"{x:9.4f}" for x in s_cs),
             "  s (a)        " + "  ".join(f"{x:9.4f}" for x in s_a),
             "  s box_g      " + "  ".join(f"{x:9.4f}" if np.isfinite(x) else "   no sol" for x in s_box),
             f"  P2 nonlinearity: RMS(s_cs - s_a) = {dev_cs:.4f} vs RMS(1 - s_a)/4 = {dev_lin / 4.0:.4f} -> "
             f"{'met' if P2 else 'NOT met'}",
             "  iterations   " + "  ".join(f"{sols[q]['iterations']:9d}" for q in q_list)
             + f"   (q_lin: {sols['lin']['iterations']})",
             "instrument: " + "; ".join(f"{k_}: {'OK' if v else 'FAILS'}" for k_, v in inst.items())
             + "; residuals " + ", ".join(f"{v['residual']:.1e}" for v in sols.values())]
    data = {"nodes": nodes, "rows": rows, "quad": quad, "q_lin": q_lin, "q_list": q_list, "A_flat": A_flat,
            "A_curved": {str(k_): v for k_, v in A.items()}, "k": k_cs, "s": s_cs, "s_a": s_a, "s_box": s_box,
            "k_a": k_a, "k_box": k_box, "continuum_A": {m: {str(q): v for q, v in d.items()} for m, d in Ac.items()},
            "Pi": {str(k_): v["Pi"] for k_, v in sols.items()}, "Pi_flat": flat["Pi"],
            "solves": {str(k_): {kk: v[kk] for kk in ("iterations", "converged", "residual", "steps")} for k_, v in sols.items()},
            "instrument": inst, "P1": P1, "P2": P2}
    return Result(status, lines, data)
