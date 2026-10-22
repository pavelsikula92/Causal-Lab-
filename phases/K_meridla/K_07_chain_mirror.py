"""
K.7 -- Mirror identity of the chain functionals
==============================================

Gate of the motion tests (F2.x; in GTS Lab 3.x the E46 gate, R136). With no
field (eps = 0) every functional reading d of the mid-time displacement
must change sign exactly under the mirror M: x1 -> -x1 of the point set:
d(P) + d(MP) = 0 to round-off. It holds because links use squared
differences only and the dynamic programmes depend only on the relation
and the index order; a reading that breaks it carries a systematic
direction of its own and cannot measure bending. The readings must also be
non-trivial (spread over seeds > 0), and each functional's chain weights
must sum to 1.

Functionals (core.chains.readings): L, S, U (conformal representation),
H1, H2, G (weighted representation).
"""

from __future__ import annotations

import time

import numpy as np

from core.order.sprinkle import derive_seed
from core.chains.readings import mirror_residuals
from core.stats.parallel import parallel_map
from core.testing import Result

TEST = {
    "id": "K.7",
    "phase": "K",
    "title": "Zrcadlová identita funkcionálů řetězců",
    "theory_part": "měřidlo pro F2 (pohyb v poli); bez fyzikálního závěru",
    "question": "Mění čtení posunu každého funkcionálu znaménko přesně při zrcadlení bodové množiny (bez pole)?",
    "method": "eps = 0; čtení na P a na zrcadlové MP; funkcionály L, S, U, H1, H2, G",
    "prediction": "|d(P) + d(MP)| = 0 na zaokrouhlení; rozptyl čtení > 0; váhy řetězců dávají součet 1",
    "criterion": "nejhorší |d(P) + d(MP)| <= mirror_tol pro všechny funkcionály; sd čtení > 0; |součet vah - 1| <= weight_tol",
    "null_model": "—",
    "level": "calibration",
    "meters": [],
    "depends_on": [],
    "old_id": "E46",
    "params": {"mirror_tol": 1e-12, "weight_tol": 1e-6, "N0": 4000, "T": 1.0, "n_seeds": 8},
    "parts": {"d3": {"D": 3}, "d4": {"D": 4}},
    "quick": {"N0": 600, "n_seeds": 2},
    "estimate": "několik minut",
}


def run(p, ctx):
    D = int(p["D"])
    t0 = time.time()
    args = [(D, int(p["N0"]), derive_seed(ctx.seed, "K.7", D, k), float(p["T"])) for k in range(int(p["n_seeds"]))]
    res = parallel_map(mirror_residuals, args, ctx.workers, ctx.progress,
                       lambda i, r_, d, n_: f"[K.7] {D-1}+1 seed {i+1} done ({d}/{n_}) ({time.time()-t0:.0f}s)")
    keys = list(res[0]["residual"])
    worst = {k: max(r["residual"][k] for r in res) for k in keys}
    sd = {k: float(np.std([r["d"][k] for r in res], ddof=1)) if len(res) > 1 else float("nan") for k in keys}
    ww = max(r["weight_worst"] for r in res)
    ok = all(np.isfinite(worst[k]) and worst[k] <= p["mirror_tol"] for k in keys) \
        and all(np.isfinite(sd[k]) and sd[k] > 0 for k in keys) and ww <= p["weight_tol"]
    lines = [f"{D-1}+1 N0={p['N0']}, {len(res)} seeds",
             "  worst |d(P) + d(MP)|: " + ", ".join(f"{k} {worst[k]:.1e}" for k in keys) + f" (tol {p['mirror_tol']:g})",
             "  sd of the readings:  " + ", ".join(f"{k} {sd[k]:.3f}" for k in keys),
             f"  chain weights sum to 1 within {ww:.1e} (tol {p['weight_tol']:g})"]
    return Result("PASSED" if ok else "FAILED", lines, {"worst": worst, "sd": sd, "weight_worst": ww})
