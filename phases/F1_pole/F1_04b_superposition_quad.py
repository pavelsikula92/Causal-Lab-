"""
F1.4b -- Superposition in (a), F1.4 with the quadrature calibrated by K.11b
===========================================================================

Revision of F1.4 (18 Oct 2026, the author's instruction): F1.4 PASSED
narrowly (RMS(Delta_cs - Delta_a) = 0.0101 against 0.0111) with the
quadrature [160, 8, 48, 6, 24, 1], which K.11b (run afterwards) puts 1.3 %
of the source charge outside its tolerance for the volume-weighted far
residue. F1.4b changes ONLY the quadrature, to the coarsest level within
the K.11b tolerance, [240, 8, 96, 8, 48, 1]. Sources (core sigma = 8,
shell r0 = 30, width 4, strengths), boundary, node grid, the definition of
the deficit, the evaluation range and the decision threshold are those of
F1.4. F1.4 stays recorded with its own result.

What it separates: if the result moves, the quadrature error mattered;
the second known systematic, the narrow shell (width 4 l_k, nonlocal
response ~0.77 of the local limit, K.11 local), is NOT removed here -- it
is diagnosed separately by F1.4c. If F1.4b stays just below the threshold,
the test is still not robust.

Original description of F1.4 (unchanged except as stated above):


Why this test is not trivial
----------------------------
For a linear equation on a fixed causal set the field of two sources is
exactly the sum of their fields; testing that would test linear algebra.
Under the author's choice (a), B_g Pi = kappa rho_s with B_g the operator of
the causal set in the geometry g = Pi^2 eta the field creates, the equation
is nonlinear and superposition fails at second order, by an amount the
theory predicts. In the local limit (calibrated in K.11 local) (a) reads

    lap Pi^p = (p / D) kappa Pi^(1+p) rho_s        (p = 1 + c_D / D; p = 1 in 2+1)

so outside the sources Pi^p is harmonic, but every source is weighted by
Pi^(1+p) at its own position, and that includes the field of the OTHER
source: each source sitting in the other's well is weakened. The
superposition deficit

    Delta(r) = u_12(r) - u_1(r) - u_2(r),       u = Pi - 1,

is therefore nonzero, of second order, and of a definite shape.
Alternatives: strict superposition (Delta = 0: any linear theory, e.g. no
back-reaction), and box_g alone (no -R/2) where a static solution exists
(reported).

Configuration (radial, so the calibrated radial machinery applies)
------------------------------------------------------------------
Two concentric sources: a Gaussian core q_c exp(-r^2 / 2 sigma^2),
sigma = 8 l_k, and a Gaussian shell q_s exp(-(r - r0)^2 / 2 w_s^2) at
r0 = 30 l_k with width w_s = 4 l_k. A shell is the classical second source
for this question: in a linear theory its field inside is constant, so
any change of the core's field by the shell is nonlinear. Strengths from
the continuum equation of (a) before any causal-set run: together the
central Pi is ~0.67 (core alone 0.71, shell alone 0.90) in 2+1. Cylinder
radius R = 100 l_k, rows at t_x = 3 R.

What is computed
----------------
Three self-consistent solutions of E[B_g(Pi)] Pi = rhs (core alone, shell
alone, both; core.operators.bd_expect_curved.solve_selfconsistent, vacuum
subtracted as in F1.5) at the quadrature level calibrated by K.11
(2+1: level 1), the deficit Delta_cs on the nodes with r <= r_eval, and the
same from the continuum equation of (a) on a fine grid, read on the same
nodes (Delta_a).

Criterion (stated before the run; the quarter rule of F1.3)
-----------------------------------------------------------
  S: RMS(Delta_cs - Delta_a) <= RMS(Delta_a) / 4 over the nodes with
     r <= r_eval (three times closer to the deficit of (a) than to strict
     superposition, Delta = 0).
Instrument: the three solves converged and their residuals are <= res_tol
times the largest source density. PASSED: S. FAILED: S does not hold.
Reported only: the deficit at the centre, the sign of the deficit
everywhere (u < 0 in a well, so Delta > 0 means the joint well is shallower
than the sum: weakening), box_g alone where it exists.

The 3+1 part follows after K.11b (see F1.5: far rows in 3+1 still carry a
quadrature residue whose volume-weighted effect has not been calibrated).

Null model: no sampling enters; the competing hypothesis is strict
superposition, which the criterion separates.
"""

from __future__ import annotations

import time

import numpy as np

from core.operators.bd_expect import shape_nodes
from core.operators.bd_expect_curved import continuum_radial, solve_selfconsistent
from core.stats.parallel import parallel_map
from core.testing import Result

TEST = {
    "id": "F1.4b",
    "phase": "F1",
    "title": "Superpozice v (a): F1.4 s kvadraturou podle K.11b",
    "theory_part": "A5′ jako samokonzistentní rovnice B_g Π = κρ_s (rozhodnutí autora 18. 10. 2026)",
    "question": "Sčítají se pole dvou zdrojů, nebo vzniká deficit superpozice druhého řádu přesně podle nelinearity rovnice (a)?",
    "method": "soustředné jádro (gauss σ = 8 l_k) a slupka (r0 = 30, šířka 4 l_k); tři samokonzistentní řešení (jádro, slupka, "
              "obojí); deficit Δ = u12 − u1 − u2; srovnání se spojitou rovnicí (a) se stejnými zdroji a okrajem",
    "prediction": "Δ_cs podle spojité rovnice (a) (zdroj ve studni druhého zdroje slábne), ne Δ = 0",
    "criterion": "MĚŘIDLO: tři řešení zkonvergovala a reziduum <= res_tol·max(rhs). S: RMS(Δ_cs − Δ_a) <= RMS(Δ_a)/4 "
                 "pro uzly r <= r_eval. PROŠEL: S. NEPROŠEL: S neplatí.",
    "null_model": "bez vzorkování; konkurenční hypotéza: přesná superpozice Δ = 0 (každá lineární teorie, i bez zpětné vazby)",
    "level": "mechanism",
    "meters": ["K.4", "K.11", "K.11b"],
    "depends_on": ["K.11", "K.11b", "F1.3"],
    "params": {"sigma": 8.0, "r0": 30.0, "w_s": 4.0, "R": 100.0, "dr": 1.0, "ratio": 1.08, "rows_per_element": 2,
               "r_eval": 50.0, "tol": 1e-5, "max_iter": 30, "res_tol": 0.05, "chunk": 4},
    "parts": {"d3": {"D": 3, "quad": [240, 8, 96, 8, 48, 1], "q_core": 0.01, "q_shell": 0.001}},
    "revision_of": "F1.4",
    "quick": {"sigma": 3.0, "r0": 10.0, "w_s": 2.0, "R": 30.0, "dr": 1.0, "ratio": 1.2, "r_eval": 15.0, "tol": 1e-4,
              "max_iter": 20},
    "quick_parts": {"d3": {"quad": [40, 6, 12, 3, 8, 1], "q_core": 0.04, "q_shell": 0.01}},
    "estimate": "2+1: jako F1.4, kvadratura ~8× dražší: ~13 h s 6 procesy, ~5 h s 16",
}

_EPS, _RHO = 0.05, 20.0          # eps * rho = 1: lengths in units of l_k


def run(p, ctx):
    D = int(p["D"])
    R, sig, r0, ws = float(p["R"]), float(p["sigma"]), float(p["r0"]), float(p["w_s"])
    qc, qs = float(p["q_core"]), float(p["q_shell"])
    # uniform nodes must cover the shell: shape_nodes is uniform up to 2a, here past r0 + 3 w_s
    nodes = shape_nodes(R, 0.5 * (r0 + 3.0 * ws), float(p["dr"]), float(p["ratio"]))
    per = int(p["rows_per_element"])
    fr = (np.arange(per) + 0.5) / per
    rows = np.sort((nodes[:-1, None] + fr[None, :] * np.diff(nodes)[:, None]).ravel())
    core = lambda r: qc * np.exp(-np.asarray(r, float) ** 2 / (2.0 * sig ** 2))
    shell = lambda r: qs * np.exp(-(np.asarray(r, float) - r0) ** 2 / (2.0 * ws ** 2))
    both = lambda r: core(r) + shell(r)
    quad = [int(x) for x in p["quad"]]
    t0 = time.time()
    sols = {}
    for name, f in (("core", core), ("shell", shell), ("both", both)):
        # "both" starts from core + shell (only the starting point; the fixed point does not depend on it)
        start = (sols["core"]["Pi"] + sols["shell"]["Pi"] - 1.0) if name == "both" else None
        sols[name] = solve_selfconsistent(D, R, nodes, rows, f(rows), quad, parallel_map, ctx.workers, ctx.progress,
                                          label=f"[F1.4b] {D-1}+1 {name}", rho=_RHO, eps=_EPS, tol=float(p["tol"]),
                                          max_iter=int(p["max_iter"]), chunk=int(p["chunk"]), Pi0=start)
        ctx.progress(f"[F1.4b] {D-1}+1 {name}: Pi(0) = {sols[name]['Pi'][0]:.4f}, {sols[name]['iterations']} iterations "
                     f"({time.time()-t0:.0f}s)")
    m = nodes <= float(p["r_eval"])
    u = {k: v["Pi"] - 1.0 for k, v in sols.items()}
    d_cs = u["both"] - u["core"] - u["shell"]

    def cont(model):
        out = {}
        for name, f in (("core", core), ("shell", shell), ("both", both)):
            c = continuum_radial(model, D, R, f)
            out[name] = np.interp(nodes, c["r"], c["Pi"]) - 1.0
        return out["both"] - out["core"] - out["shell"], out
    d_a, ua = cont("a")
    try:
        d_box, _ = cont("box_g")
    except ValueError:
        d_box = None                                  # no static solution of box_g at these strengths
    rms = lambda x: float(np.sqrt(np.mean(x[m] ** 2)))
    dev = rms(d_cs - d_a)
    ref = rms(d_a)
    S = dev <= ref / 4.0
    rhs_max = float(np.max(both(rows)))
    inst = {"self-consistent solves converged": all(v["converged"] for v in sols.values()),
            "residuals": all(v["residual"] <= float(p["res_tol"]) * rhs_max for v in sols.values())}
    if not all(inst.values()):
        status = "INSTRUMENT_FAILED"
    elif S:
        status = "PASSED"
    else:
        status = "FAILED"
    probe = [0.0, sig, r0, float(p["r_eval"])]
    at = lambda arr, x: float(np.interp(x, nodes, arr))
    lines = [f"{D-1}+1 self-consistent field (a), lengths in l_k: core sigma = {sig:g} (q = {qc:g}), shell r0 = {r0:g} "
             f"w = {ws:g} (q = {qs:g}); R = {R:g}, {nodes.size} nodes, {rows.size} rows, quadrature {quad}",
             "  Pi(0): core " + f"{sols['core']['Pi'][0]:.4f}, shell {sols['shell']['Pi'][0]:.4f}, both {sols['both']['Pi'][0]:.4f}"
             + f"   (continuum (a): {1 + ua['core'][0]:.4f}, {1 + ua['shell'][0]:.4f}, {1 + ua['both'][0]:.4f})",
             "  r                 " + "  ".join(f"{x:9g}" for x in probe),
             "  deficit causal    " + "  ".join(f"{at(d_cs, x):+9.5f}" for x in probe),
             "  deficit (a)       " + "  ".join(f"{at(d_a, x):+9.5f}" for x in probe),
             "  deficit box_g     " + ("  ".join(f"{at(d_box, x):+9.5f}" for x in probe) if d_box is not None
                                       else "no static solution at these strengths"),
             f"  S: RMS(Delta_cs - Delta_a) = {dev:.5f} vs RMS(Delta_a)/4 = {ref / 4.0:.5f} (r <= {p['r_eval']:g}) -> "
             f"{'met' if S else 'NOT met'}; sign of the deficit (causal set): "
             f"{'positive everywhere (the joint well is shallower than the sum: weakening)' if np.all(d_cs[m] >= 0) else 'negative everywhere (strengthening)' if np.all(d_cs[m] <= 0) else 'mixed'}",
             "  iterations: " + ", ".join(f"{k} {v['iterations']}" for k, v in sols.items()),
             "instrument: " + "; ".join(f"{k}: {'OK' if v else 'FAILS'}" for k, v in inst.items())
             + "; residuals " + ", ".join(f"{v['residual']:.1e}" for v in sols.values())]
    data = {"nodes": nodes, "rows": rows, "quad": quad, "deficit_cs": d_cs, "deficit_a": d_a,
            "deficit_box": d_box, "u": u, "u_continuum_a": ua, "rms_dev": dev, "rms_ref": ref,
            "solves": {k: {kk: v[kk] for kk in ("iterations", "converged", "residual", "steps")} for k, v in sols.items()},
            "instrument": inst, "S": S}
    return Result(status, lines, data)
