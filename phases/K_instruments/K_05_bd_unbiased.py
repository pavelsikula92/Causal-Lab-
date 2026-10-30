"""
K.5 -- The sampled BD operator is unbiased against the expected one
==================================================================

(BD-1 of GTS Lab 3.x E45b.) The sampled rows of the static BD operator
(core.operators.bd.bd_static_rows) at stratified elements of one
sprinkling, applied to a test vector (the continuum solution of a ball
source on the radial nodes), against the EXPECTED rows at the same
elements (core.operators.bd_expect, quadrature, no noise).

The rows of one sprinkling are strongly correlated (every row's past
contains the same source region), so the statistic is clustered: each of
n_sprinklings independent sprinklings contributes ONE mean of
(sampled - expected) . u, and the SE comes from those means. (A pooled
statistic read z = -4 in GTS Lab 3.1.0 smoke tests; the clustered one
+0.004 +- 0.008.)
"""

from __future__ import annotations

import time

import numpy as np

from core.order.sprinkle import sprinkle, Cylinder, derive_seed
from core.operators.bd import bd_static_rows, stratified_rows
from core.operators.bd_expect import expected_radial_rows
from core.operators.continuum import ball_source_solution
from core.stats.stats import mean_se, fmt
from core.stats.parallel import parallel_map
from core.testing import Result

TEST = {
    "id": "K.5",
    "phase": "K",
    "title": "Unbiasedness of the discrete BD operator",
    "theory_part": "instrument for A5′ (F1.2); no physical conclusion",
    "question": "Is the discrete BD operator of one sprinkling on average equal to the expected operator?",
    "method": "stratified rows in a time band; sampled minus expected rows applied to a test vector; one mean per "
              "sprinkling (clustering)",
    "prediction": "difference 0 within the statistical error",
    "criterion": "|z| <= z_max, z = mean of the means / their SE over independent sprinklings",
    "null_model": "—",
    "level": "calibration",
    "meters": ["K.4"],
    "depends_on": [],
    "old_id": "E45b",
    "params": {"z_max": 3.0, "R": 0.5, "T": 1.0, "a": 0.05, "q": 1.0, "n_nodes": 20, "bd_t0": 0.8, "bd_t1": 1.0,
               "bd_r_rows": 0.45, "rows_per_element": 2, "n_sprinklings": 16},
    "parts": {"d3_e0.2": {"D": 3, "N": 32000, "eps": 0.2}, "d3_e0.05": {"D": 3, "N": 32000, "eps": 0.05},
              "d4_e0.2": {"D": 4, "N": 32000, "eps": 0.2}},
    "quick": {"N": 2000, "n_sprinklings": 4, "rows_per_element": 1},
    "estimate": "a few minutes per part",
}


def one_sprinkling(p, seed):
    """One clustered unit: mean of (sampled - expected) . u_ref over the stratified rows."""
    D, R = int(p["D"]), float(p["R"])
    S = sprinkle(Cylinder(D, 0.0, float(p["T"]), R), int(p["N"]), seed)
    r = S.r()
    nodes = np.linspace(0.0, R, int(p["n_nodes"]) + 1)
    u_ref = ball_source_solution(D - 1, float(p["q"]), float(p["a"]), R)(nodes)
    rows = stratified_rows(S.X, r, float(p["bd_t0"]), float(p["bd_t1"]), nodes, float(p["bd_r_rows"]),
                           int(p["rows_per_element"]), np.random.default_rng(derive_seed(seed, "rows")))
    B = bd_static_rows(S.X, S.rho, float(p["eps"]), rows, nodes, R)
    Ae = np.vstack([expected_radial_rows(D, S.rho, float(p["eps"]), float(S.X[i, 0]), [float(r[i])], R, nodes) for i in rows])
    return {"mean": float(np.mean(B["A"] @ u_ref - Ae @ u_ref)), "scale": float(np.sqrt(np.mean((Ae @ u_ref) ** 2))),
            "rows": int(rows.size), "l_k": (float(p["eps"]) * S.rho) ** (-1.0 / D)}


def run(p, ctx):
    t0 = time.time()
    n = int(p["n_sprinklings"])
    args = [(p, derive_seed(ctx.seed, "K.5", int(p["D"]), int(p["N"]), float(p["eps"]), k)) for k in range(n)]
    res = parallel_map(one_sprinkling, args, ctx.workers, ctx.progress,
                       lambda i, r_, d, n_: f"[K.5] sprinkling {i+1} done ({d}/{n_}): mean diff {r_['mean']:+.4f} "
                                            f"({time.time()-t0:.0f}s)")
    m, se = mean_se([x["mean"] for x in res])
    z = m / se if se and se > 0 else float("nan")
    ok = np.isfinite(z) and abs(z) <= p["z_max"]
    lines = [f"{int(p['D'])-1}+1 N={p['N']} eps={p['eps']} (l_k = {res[0]['l_k']:.4f}); {n} sprinklings x {res[0]['rows']} rows",
             f"  sampled - expected = {m:+.4f} +- {se:.4f} (z {fmt(z, 2)}, max {p['z_max']}); scale of the rows {np.mean([x['scale'] for x in res]):.3f}"]
    return Result("PASSED" if ok else "FAILED", lines, {"mean": m, "se": se, "z": z, "units": res})
