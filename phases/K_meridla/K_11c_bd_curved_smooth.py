"""
K.11c -- Far rows of the curved operator with a C^2 (smooth) geometry
=====================================================================

Why this test exists
--------------------
K.11b far_d4 FAILED: in 3+1 the volume-weighted residue of the far rows of
E[B_g] did not converge in the quadrature (Q_far/Q_src = 3.5, 2.3, 1.9 at
the levels 1-3 of K.11, tolerance 0.01), so the self-consistent field
equation (a) could not be solved in 3+1 (F1.4, F1.5 have no 3+1 part).
K.11b stays FAILED under its own criterion. This test calibrates a changed
instrument with the SAME quantity and the SAME criterion.

Development diagnosis (19 Oct 2026; numbers on record, not used for the
criterion). The plan was to expand the contribution of every ray in the
thickness of its long, nearly null intervals and integrate it analytically.
Before that, the far row r = 40 l_k of the K.11b field (3+1) was computed
with each quadrature parameter doubled separately from level 1:

    geometry interpolated piecewise linearly between the nodes (K.11b):
        E(40) = +3.36e-3 (level 1); doubling n_rho, n_s, n_ang, n_mu, n_nu,
        n_phi: +4.46e-3, +3.36e-3, +1.32e-3, +3.83e-3, +4.78e-3, +3.17e-3
    geometry = the continuum solution itself (fine grid, smooth):
        E(40) = -4.15e-5 at level 1 and within 0.2 % of it for every doubling
    geometry = C^2 cubic spline through the same nodal values:
        E(40) = -4.152e-5 (smooth reference -4.150e-5)
    geometry = C^1 monotone spline (PCHIP) through the same values:
        E(40) = +4.4e-4 (wrong)

So the residue was not a property of long intervals as such: it came from
the piecewise-linear geometry. The curved operator responds to the
CURVATURE of the geometry (its local limit carries -R/2, i.e. lap Pi^D):
a long, nearly null interval is a thin needle whose proper volume changes
with its thickness in proportion to the transverse second derivative of
Pi^D, and a piecewise-linear Pi has that derivative as delta functions at
the nodes. In 3+1 the far rows weigh these needles ~r (2+1: ~1), which is
why 2+1 passed K.11b and 3+1 did not. The analytic expansion is therefore
not needed; the instrument change is the interpolation of the geometry:
C^2 cubic spline through the nodal values (core.operators.bd_expect_curved,
geom = 'spline'; the field ACTED ON stays in the hat basis, as in K.4).

What is computed (as K.11b, with geom = 'spline')
-------------------------------------------------
Field: the continuum solution of (a) for the Gaussian source of F1.4/F1.5
(sigma = 8 l_k, R = 100 l_k; q = 0.01 in 2+1, 0.04 in 3+1, central Pi
~0.7), read at the nodes of the part's grid. Rows: r = 4 l_k (inside) and
six far radii from 24 to 50 l_k (geometric). Quadrature levels 1-3 of K.11.
For this field the continuum value of (a) is the source itself (at 24 l_k
= 3 sigma still 1 % of its central value, printed with the rows); the
criterion compares levels, so that tail is common to all of them.

  Q_far(level) = trapezoid over the six far radii of E(r) omega_n r^(n-1)
  Q_src        = q (2 pi sigma^2)^(n/2)

Two node grids per dimension: g1 = step 1 l_k, ratio 1.08 (K.11b, F1.4,
F1.4b, F1.4c) and g2 = step 2 l_k, ratio 1.12 (current F1.5).

Criterion (K.11b, unchanged)
----------------------------
  |Q_far(L3) - Q_far(L2)| <= conv_tol * Q_src,   conv_tol = 0.01.
Instrument (else INSTRUMENT_FAILED): at level 2 the spline geometry agrees
with the smooth reference geometry (the continuum solution on its fine
grid, same field acted on): |Q_far(spline) - Q_far(reference)| <=
conv_tol * Q_src, i.e. the interpolation itself adds no residue above the
tolerance.
Reported: Q_far/Q_src at every level, the coarsest level within tolerance
of the finest (the quadrature 3+1 field tests will use), the inner row
against the source, and Q_far/Q_src with the piecewise-linear geometry at
level 2 (the K.11b instrument on the same grid): in 2+1 this says how much
the linear geometry of F1.4 and F1.5 differed from a smooth one.

No null model: deterministic quadrature against refinement and against a
smooth reference.
"""

from __future__ import annotations

import time

import numpy as np

from core.operators.bd_expect import shape_nodes
from core.operators.bd_expect_curved import continuum_radial, curved_radial_rows, omega_n
from core.stats.parallel import parallel_map
from core.testing import Result

TEST = {
    "id": "K.11c",
    "phase": "K",
    "title": "Vzdálené řádky zakřiveného operátoru s hladkou (C²) geometrií",
    "theory_part": "měřidlo pro samokonzistentní rovnici pole (a) ve 3+1 (F1.4, F1.5, F2.7); bez fyzikálního závěru",
    "question": "Konverguje objemově vážený zbytek vzdálených řádků zakřiveného operátoru ve 3+1, když je geometrie mezi "
                "uzly interpolována C² kubickým splajnem místo po částech lineárně?",
    "method": "jako K.11b (spojité řešení (a) pro zdroj σ = 8 l_k, R = 100; řádek uvnitř a šest vzdálených řádků 24–50 l_k; "
              "úrovně kvadratury 1–3 z K.11), geometrie C² splajnem; dvě sítě uzlů; reference s hladkou geometrií "
              "spojitého řešení; pro srovnání lineární geometrie na úrovni 2",
    "prediction": "|Q_far(L3) − Q_far(L2)| <= 0,01 Q_src (vývojový běh: E(40) ve 3+1 stabilní na 0,2 % při zdvojení každého "
                  "parametru kvadratury)",
    "criterion": "MĚŘIDLO: na úrovni 2 |Q_far(splajn) − Q_far(reference)| <= 0,01 Q_src. PROŠEL: |Q_far(L3) − Q_far(L2)| "
                 "<= 0,01 Q_src. NEPROŠEL: jinak.",
    "null_model": "—",
    "level": "calibration",
    "meters": ["K.11"],
    "depends_on": ["K.11", "K.11b"],
    "revision_of": "K.11b",
    "params": {"sigma": 8.0, "R": 100.0, "r_in": 4.0, "r_far": [24.0, 50.0], "n_far": 6, "conv_tol": 0.01,
               "ref_grid": 40001},
    "parts": {
        "far_d4_g1": {"D": 4, "q": 0.04, "dr": 1.0, "ratio": 1.08,
                      "levels": [[160, 8, 32, 6, 24, 6], [240, 8, 64, 8, 48, 8], [320, 10, 96, 10, 64, 10]]},
        "far_d4_g2": {"D": 4, "q": 0.04, "dr": 2.0, "ratio": 1.12,
                      "levels": [[160, 8, 32, 6, 24, 6], [240, 8, 64, 8, 48, 8], [320, 10, 96, 10, 64, 10]]},
        "far_d3_g1": {"D": 3, "q": 0.01, "dr": 1.0, "ratio": 1.08,
                      "levels": [[160, 8, 48, 6, 24, 1], [240, 8, 96, 8, 48, 1], [320, 10, 144, 10, 64, 1]]},
        "far_d3_g2": {"D": 3, "q": 0.01, "dr": 2.0, "ratio": 1.12,
                      "levels": [[160, 8, 48, 6, 24, 1], [240, 8, 96, 8, 48, 1], [320, 10, 144, 10, 64, 1]]},
    },
    "quick": {"sigma": 3.0, "R": 30.0, "r_in": 2.0, "r_far": [9.0, 15.0], "n_far": 3, "ref_grid": 6001},
    "quick_parts": {"far_d4_g1": {"levels": [[40, 6, 8, 3, 8, 3], [60, 6, 12, 4, 12, 4]]},
                    "far_d4_g2": {"levels": [[40, 6, 8, 3, 8, 3], [60, 6, 12, 4, 12, 4]]},
                    "far_d3_g1": {"levels": [[40, 6, 12, 3, 8, 1], [60, 6, 16, 4, 12, 1]]},
                    "far_d3_g2": {"levels": [[40, 6, 12, 3, 8, 1], [60, 6, 16, 4, 12, 1]]}},
    "estimate": "d3 části ~20–40 min; d4 části ~8 h na 1 jádře (nejjemnější úroveň dominuje), s 10–16 procesy ~1 h; "
                "spouštěj části zvlášť",
}

_EPS, _RHO = 0.05, 20.0          # eps * rho = 1: lengths in units of l_k


def row_task(D, R, nodes, Pi_nodes, r_row, lev, geom, g_nodes, g_values):
    """E[B_g Pi](r_row) at one quadrature level and geometry (module level: worker processes)."""
    row = curved_radial_rows(D, _RHO, _EPS, 3.0 * R, [r_row], R, nodes, Pi_nodes, *[int(x) for x in lev],
                             geom=geom, geom_nodes=g_nodes, geom_values=g_values)[0]
    return float(row @ Pi_nodes)


def run(p, ctx):
    D, n = int(p["D"]), int(p["D"]) - 1
    R, sig, q = float(p["R"]), float(p["sigma"]), float(p["q"])
    nodes = shape_nodes(R, 2.0 * sig, float(p["dr"]), float(p["ratio"]))
    src = lambda r: q * np.exp(-np.asarray(r, float) ** 2 / (2.0 * sig ** 2))
    c = continuum_radial("a", D, R, src, n_grid=int(p["ref_grid"]))
    Pin = np.interp(nodes, c["r"], c["Pi"])
    far = np.geomspace(float(p["r_far"][0]), float(p["r_far"][1]), int(p["n_far"]))
    radii = [float(p["r_in"])] + [float(x) for x in far]
    levels = p["levels"]
    # runs: the three levels with the spline geometry, then at level 2 (the second finest) the smooth reference
    # (continuum solution on its fine grid, interpolated linearly at spacing R / ref_grid) and the linear geometry
    runs = [("spline", lev, None, None) for lev in levels]
    runs += [("reference", levels[-2], c["r"], c["Pi"]), ("linear", levels[-2], None, None)]
    args = []
    for kind, lev, gn, gv in runs:
        geom = "spline" if kind == "spline" else "linear"
        for r in radii:
            args.append((D, R, nodes, Pin, r, lev, geom, gn, gv))
    t0 = time.time()
    labels = [f"{kind} {list(lev)}" for kind, lev, _, _ in runs]
    vals = parallel_map(row_task, args, ctx.workers, ctx.progress,
                        lambda i, r_, d, n_: f"[K.11c] {D-1}+1 {labels[i // len(radii)]}, r = {radii[i % len(radii)]:.1f}: "
                                             f"{r_:+.4e} ({d}/{n_}, {time.time()-t0:.0f}s)")
    E = np.array(vals).reshape(len(runs), len(radii))
    w = omega_n(n) * far ** (n - 1)
    Qfar = np.array([float(np.trapezoid(E[k, 1:] * w, far)) for k in range(len(runs))])
    Qsrc = q * (2.0 * np.pi * sig ** 2) ** (n / 2.0)
    tol = float(p["conv_tol"]) * Qsrc
    nl = len(levels)
    Qs, Qref, Qlin = Qfar[:nl], Qfar[nl], Qfar[nl + 1]
    conv_ok = bool(abs(Qs[-1] - Qs[-2]) <= tol)
    inst_ok = bool(abs(Qs[-2] - Qref) <= tol)
    within = [li for li in range(nl) if abs(Qs[li] - Qs[-1]) <= tol]
    chosen = min(within) if within else None
    status = "INSTRUMENT_FAILED" if not inst_ok else ("PASSED" if conv_ok else "FAILED")
    lines = [f"{D-1}+1: field = continuum solution of (a), Gaussian source sigma = {sig:g}, q = {q:g} (Pi(0) = {Pin[0]:.3f}), "
             f"R = {R:g}; grid step {p['dr']:g}, ratio {p['ratio']:g}, {nodes.size} nodes; Q_src = {Qsrc:.4f}",
             "  r                    " + "  ".join(f"{r:10.1f}" for r in radii),
             "  source (continuum)   " + "  ".join(f"{float(src(r)):+10.3e}" for r in radii)]
    for k, (kind, lev, _, _) in enumerate(runs):
        name = f"spline L{k + 1}" if kind == "spline" else f"{kind} L{nl - 1}"
        lines.append(f"  {name:20s} " + "  ".join(f"{v:+10.3e}" for v in E[k]) + f"   Q_far/Q_src = {Qfar[k] / Qsrc:+.5f}  "
                     f"{list(lev)}")
    lines += [f"  convergence: |Q_far(L3) - Q_far(L2)| / Q_src = {abs(Qs[-1] - Qs[-2]) / Qsrc:.5f} (tol {p['conv_tol']}) -> "
              f"{'OK' if conv_ok else 'FAILS'}",
              f"  instrument: |Q_far(spline L{nl - 1}) - Q_far(reference)| / Q_src = {abs(Qs[-2] - Qref) / Qsrc:.5f} -> "
              f"{'OK' if inst_ok else 'FAILS'}",
              f"  linear geometry (K.11b instrument) at L{nl - 1}: Q_far/Q_src = {Qlin / Qsrc:+.5f}; difference from the smooth "
              f"reference {(Qlin - Qref) / Qsrc:+.5f} (reported only)",
              f"  coarsest level within tolerance of the finest: {('L' + str(chosen + 1)) if chosen is not None else 'none'}"
              + (f" {list(levels[chosen])}" if chosen is not None else "")]
    data = {"radii": radii, "E": E, "runs": [(k, list(l)) for k, l, _, _ in runs], "Q_far": Qfar, "Q_src": Qsrc,
            "Q_far_spline": Qs, "Q_far_reference": Qref, "Q_far_linear": Qlin, "chosen_level": chosen,
            "levels": levels, "nodes": nodes, "Pi": Pin}
    return Result(status, lines, data)
