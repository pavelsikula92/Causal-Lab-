"""
K.11b -- Far rows of the curved operator: convergence of the volume-weighted residue
====================================================================================

Why this test exists
--------------------
K.11 conv passed under its criterion: between the two finest quadrature
levels every row of E[B_g] agrees within 1 % of the operator's response
inside the source. In 3+1 the far row r = 30 l_k nevertheless kept a
residue of a few 1e-4 that changed sign between levels (-8.3e-4, -4.7e-4,
+4.5e-4). For the self-consistent field of (a) (F1.4, F1.5) a residue in
the far rows acts like a distributed source, and its effect on the far
field grows with the volume it occupies (r^(n-1) dr out to tens of l_k):
a per-row criterion does not measure that. K.11 stays PASSED under its
own criterion; this test measures the quantity F1.4 and F1.5 depend on.

Development diagnosis (18 Oct, recorded, not used for the criterion): the
residue comes from rays whose intervals reach BEYOND the source
(rho_s > 40 l_k for a row at 30 l_k), mostly in directions close to the
inward one; it is the proper-volume (Phi) dependence of those long, nearly
null intervals. Concentrating the directions on the source did not remove
the level-to-level changes, so the instrument is left exactly as K.11
calibrated it.

What is computed
----------------
The field is the one F1.4 and F1.5 actually meet: the continuum solution
of (a) for their Gaussian source (sigma = 8 l_k, R = 100 l_k, central Pi
~0.7: q = 0.01 in 2+1, 0.04 in 3+1), on the node grid of F1.5. For this
field the continuum value of (a) is exactly the source inside and exactly
0 outside, so in the far rows the causal-set operator shows only its own
nonlocal response plus quadrature error.

Rows: r = 4 l_k (inside, the scale) and six far radii from 3 sigma = 24 to
R / 2 = 50 l_k (geometric). Quadrature levels: K.11 levels 1, 2, 3.

  Q_far(level) = integral over the far radii of E(r) omega_n r^(n-1) dr
                 (trapezoid on the six radii, the same rule at every level)
  Q_src        = integral of the source = q (2 pi sigma^2)^(n/2)

Criterion (one order below the 10 % field accuracy of K.3, as in K.11):
  |Q_far(L3) - Q_far(L2)| <= conv_tol * Q_src,  conv_tol = 0.01.
Reported: Q_far / Q_src at every level (if it converges to a nonzero
value, the causal set of (a) carries a genuine far response, a finding in
its own right), the coarsest level within tolerance of the finest (the
quadrature F1.4 and F1.5 will use in 3+1), and the inner row against the
source.

No null model: deterministic quadrature against refinement.
"""

from __future__ import annotations

import time

import numpy as np

from core.operators.bd_expect import shape_nodes
from core.operators.bd_expect_curved import continuum_radial, curved_radial_rows, omega_n
from core.stats.parallel import parallel_map
from core.testing import Result

TEST = {
    "id": "K.11b",
    "phase": "K",
    "title": "Far rows of the curved operator: the volume-weighted remainder",
    "theory_part": "instrument for the self-consistent field equation (a) in 3+1 (F1.4, F1.5); no physical conclusion",
    "question": "Does the volume-weighted remainder of the curved operator in the far rows (an effective spurious "
                "source) converge in the quadrature to the level that F1.4 and F1.5 need?",
    "method": "field = continuum solution of (a) for the source of F1.4/F1.5 (σ = 8 l_k, R = 100); a row inside and "
              "six far rows (24–50 l_k); quadrature levels 1–3 of K.11; Q_far = ∫ E ω r^(n−1) dr over the far rows",
    "prediction": "|Q_far(L3) − Q_far(L2)| <= 0.01 Q_src",
    "criterion": "PASSED: |Q_far(L3) − Q_far(L2)| <= conv_tol · Q_src. FAILED: otherwise.",
    "null_model": "—",
    "level": "calibration",
    "meters": ["K.11"],
    "depends_on": ["K.11"],
    "params": {"sigma": 8.0, "R": 100.0, "dr": 1.0, "ratio": 1.08, "r_in": 4.0, "r_far": [24.0, 50.0], "n_far": 6,
               "conv_tol": 0.01},
    "parts": {
        "far_d3": {"D": 3, "q": 0.01, "levels": [[160, 8, 48, 6, 24, 1], [240, 8, 96, 8, 48, 1],
                                                 [320, 10, 144, 10, 64, 1]]},
        "far_d4": {"D": 4, "q": 0.04, "levels": [[160, 8, 32, 6, 24, 6], [240, 8, 64, 8, 48, 8],
                                                 [320, 10, 96, 10, 64, 10]]},
    },
    "quick": {"sigma": 3.0, "R": 30.0, "dr": 1.0, "ratio": 1.2, "r_in": 2.0, "r_far": [9.0, 15.0], "n_far": 3},
    "quick_parts": {"far_d3": {"levels": [[40, 6, 12, 3, 8, 1], [60, 6, 16, 4, 12, 1]]},
                    "far_d4": {"levels": [[40, 6, 8, 3, 8, 3], [60, 6, 12, 4, 12, 4]]}},
    "estimate": "far_d3 ~10–20 min; far_d4 ~5 h on 1 core, with 16 processes ~40 min (the finest level dominates)",
}

_EPS, _RHO = 0.05, 20.0          # eps * rho = 1: lengths in units of l_k


def row_task(D, R, nodes, Pi_nodes, r_row, lev):
    """E[B_g Pi](r_row) at one quadrature level (module level: worker processes)."""
    row = curved_radial_rows(D, _RHO, _EPS, 3.0 * R, [r_row], R, nodes, Pi_nodes, *[int(x) for x in lev])[0]
    return float(row @ Pi_nodes)


def run(p, ctx):
    D, n = int(p["D"]), int(p["D"]) - 1
    R, sig, q = float(p["R"]), float(p["sigma"]), float(p["q"])
    nodes = shape_nodes(R, 2.0 * sig, float(p["dr"]), float(p["ratio"]))
    src = lambda r: q * np.exp(-np.asarray(r, float) ** 2 / (2.0 * sig ** 2))
    c = continuum_radial("a", D, R, src)
    Pin = np.interp(nodes, c["r"], c["Pi"])
    far = np.geomspace(float(p["r_far"][0]), float(p["r_far"][1]), int(p["n_far"]))
    radii = [float(p["r_in"])] + [float(x) for x in far]
    levels = p["levels"]
    args = [(D, R, nodes, Pin, r, lev) for lev in levels for r in radii]
    t0 = time.time()
    vals = parallel_map(row_task, args, ctx.workers, ctx.progress,
                        lambda i, r_, d, n_: f"[K.11b] {D-1}+1 level {i // len(radii)}, r = {radii[i % len(radii)]:.1f}: "
                                             f"{r_:+.4e} ({d}/{n_}, {time.time()-t0:.0f}s)")
    E = np.array(vals).reshape(len(levels), len(radii))
    w = omega_n(n) * far ** (n - 1)
    Qfar = np.array([float(np.trapezoid(E[li, 1:] * w, far)) for li in range(len(levels))])
    Qsrc = q * (2.0 * np.pi * sig ** 2) ** (n / 2.0)
    tol = float(p["conv_tol"]) * Qsrc
    ok = bool(abs(Qfar[-1] - Qfar[-2]) <= tol)
    within = [li for li in range(len(levels)) if abs(Qfar[li] - Qfar[-1]) <= tol]
    chosen = min(within) if within else None
    lines = [f"{D-1}+1: field = continuum solution of (a), Gaussian source sigma = {sig:g}, q = {q:g} (Pi(0) = {Pin[0]:.3f}), "
             f"R = {R:g}, {nodes.size} nodes; Q_src = {Qsrc:.4f}",
             "  r                  " + "  ".join(f"{r:10.1f}" for r in radii),
             f"  source (continuum) {src(radii[0]):+10.3e}" + "  ".join(f"{0.0:+10.3e}" for _ in far).rjust(12 * len(far)),
             ]
    for li, lev in enumerate(levels):
        lines.append(f"  level {li} " + "  ".join(f"{v:+10.3e}" for v in E[li]) + f"   Q_far/Q_src = {Qfar[li] / Qsrc:+.5f}  {list(lev)}")
    lines += [f"  |Q_far(L3) - Q_far(L2)| / Q_src = {abs(Qfar[-1] - Qfar[-2]) / Qsrc:.5f} (tol {p['conv_tol']}) -> "
              f"{'OK' if ok else 'FAILS'}",
              f"  coarsest level within tolerance of the finest: {chosen if chosen is not None else 'none'}"
              + (f" {list(levels[chosen])}" if chosen is not None else "")]
    data = {"radii": radii, "E": E, "Q_far": Qfar, "Q_src": Qsrc, "chosen_level": chosen, "levels": levels,
            "nodes": nodes, "Pi": Pin}
    return Result("PASSED" if ok else "FAILED", lines, data)
