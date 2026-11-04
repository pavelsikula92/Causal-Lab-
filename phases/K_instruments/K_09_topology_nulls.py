"""
K.9 -- Topology null models
==========================

The null models that rewire the link graph (core.flux.capacity) must do
what they claim.

  (a) the covering-relation check (redundant_edges) equals brute force
      (transitive closure) on a random DAG;
  (b) the original links of a sprinkling contain no redundant edge;
  (c) degree_time_DAG_null (rewire_links, reduce = False): every edge
      forward in time, no duplicates, in/out degrees exactly kept, and
      mixed (fraction of original edges kept <= kept_max); the fraction of
      non-covering edges is reported (it is NOT a causal set's link graph);
  (d) strict null (rewire_strict): no redundant edge after the swaps,
      degrees exactly kept; acceptance and mixing are reported (the strict
      null is interpretable only where it is mixed).
"""

from __future__ import annotations

import numpy as np

from core.order.sprinkle import sprinkle, Cylinder, derive_seed
from core.order.links import links_csr
from core.flux.capacity import (redundant_edges, redundant_sample, rewire_links, rewire_strict, is_time_dag)
from core.testing import Result

TEST = {
    "id": "K.9",
    "phase": "K",
    "title": "Topology null models",
    "theory_part": "instrument for F3 (null models of structure); no physical conclusion",
    "question": "Do the rewiring null models do what they declare (validity, degrees, mixing)?",
    "method": "a check of the covering relations against brute force; a DAG null with edge swaps; a strict null with "
              "a covering check",
    "prediction": "the covering check is exact; original links without redundant edges; the DAG null valid and "
                  "mixed; the strict null without redundant edges with preserved degrees",
    "criterion": "(a) exact agreement; (b) 0 redundant edges; (c) edges forward in time, no duplicates, degrees "
                 "exact, kept fraction <= kept_max; (d) 0 redundant edges and exact degrees",
    "null_model": "—",
    "level": "calibration",
    "meters": [],
    "depends_on": [],
    "old_id": "E55b",
    "params": {"kept_max": 0.5, "N_dag": 16000, "N_check": 4000, "N_strict": 1500, "swaps_per_edge": 10.0,
               "strict_trials": 300000, "local_ell": 1.0, "redundant_sample": 300, "T": 6.0, "R": 0.5},
    "parts": {"d3": {"D": 3}},
    "quick": {"N_dag": 3000, "N_check": 1000, "N_strict": 500, "strict_trials": 5000, "redundant_sample": 50},
    "estimate": "5-15 min",
}


def run(p, ctx):
    D = int(p["D"])
    lines, ok = [f"{D-1}+1"], True
    rng = np.random.default_rng(derive_seed(ctx.seed, "K.9", "dag"))
    n = 40
    A = np.triu(rng.random((n, n)) < 0.15, 1)
    vv, uu = np.nonzero(A.T)
    ptr = np.cumsum(np.r_[0, np.bincount(vv, minlength=n)]).astype(np.int64)
    red = redundant_edges(ptr, uu, n)
    M = A.astype(int)
    P = M.copy()
    for _ in range(n):
        P = ((P + P @ M) > 0).astype(int)
    a_ok = bool(np.array_equal(red, ((M @ P) > 0)[uu, vv]) and red.any())
    ok &= a_ok
    lines.append(f"  (a) covering check vs brute force on a random DAG: {'equal' if a_ok else 'DIFFERENT'} "
                 f"({int(red.sum())} redundant of {red.size})")
    reg = Cylinder(D, 0.0, float(p["T"]), float(p["R"]))
    S = sprinkle(reg, int(p["N_check"]), derive_seed(ctx.seed, "K.9", "check"))
    p2, i2 = links_csr(S.X, S.rho)
    r0 = int(redundant_edges(p2, i2, S.n).sum())
    ok &= r0 == 0
    lines.append(f"  (b) original links N={p['N_check']}: redundant edges {r0}")
    S = sprinkle(reg, int(p["N_dag"]), derive_seed(ctx.seed, "K.9", "null"))
    p3, i3 = links_csr(S.X, S.rho)
    tt = S.X[:, 0]
    nul = rewire_links(p3, i3, tt, float(p["swaps_per_edge"]), np.random.default_rng(derive_seed(ctx.seed, "K.9", "swap")),
                       reduce=False, full_check=False)
    dag = is_time_dag(nul["ptr"], nul["idx"], tt, S.n)
    rs = redundant_sample(nul["ptr"], nul["idx"], tt, int(p["redundant_sample"]), np.random.default_rng(1))
    c_ok = (dag["forward_in_time"] and dag["no_duplicates"] and nul["deg_in_L1"] == 0.0 and nul["deg_out_L1"] == 0.0
            and nul["original_edges_kept"] <= p["kept_max"])
    ok &= c_ok
    lines.append(f"  (c) DAG null N={p['N_dag']}: forward {dag['forward_in_time']}, no duplicates {dag['no_duplicates']}, "
                 f"degrees exact {nul['deg_in_L1'] == 0.0 and nul['deg_out_L1'] == 0.0}, original edges kept "
                 f"{nul['original_edges_kept']:.3f} (max {p['kept_max']}); non-covering edges {rs['fraction']:.3f} (reported)")
    S = sprinkle(reg, int(p["N_strict"]), derive_seed(ctx.seed, "K.9", "strict"))
    p4, i4 = links_csr(S.X, S.rho)
    st = rewire_strict(p4, i4, S.X[:, 0], int(p["strict_trials"]), np.random.default_rng(derive_seed(ctx.seed, "K.9", "st")),
                       X=S.X, local=float(p["local_ell"]) * S.ell)
    d_ok = st["redundant_after"] == 0 and st["deg_in_same"] and st["deg_out_same"]
    ok &= d_ok
    lines.append(f"  (d) strict null N={p['N_strict']}: redundant after {st['redundant_after']}, degrees exact "
                 f"{st['deg_in_same'] and st['deg_out_same']}; accepted {st['accepted_swaps']} of {st['tried']} trials, "
                 f"original edges kept {st['original_edges_kept']:.3f} (reported)")
    return Result("PASSED" if ok else "FAILED", lines,
                  {"dag_kept": nul["original_edges_kept"], "dag_noncovering": rs["fraction"], "strict_kept": st["original_edges_kept"],
                   "strict_accepted": st["accepted_swaps"], "strict_tried": st["tried"]})
