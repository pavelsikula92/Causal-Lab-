"""
K.8 -- The permeability probe
============================

The measured permeability Pi_v^meas = T_out / T_in of a marked probe with
proportional acceptance (core.flux.capacity.probe_permeability): the probe
is routed with the acceptance ratios of the joint routing (background +
probe). It must (a) obey its own balance (probe injected = stored +
escaped), (b) not depend on its size (probe of delta and of size_factor x
delta), (c) not shift the background (alpha without the probe vs Pi^meas).
Checked for an inserted source and for homogeneous driving, at a capacity
that saturates part of the events. It is NOT computed from a saturation
formula.
"""

from __future__ import annotations

import numpy as np

from core.order.sprinkle import sprinkle, Cylinder, derive_seed
from core.order.links import links_csr
from core.flux.capacity import route_capacity, probe_check
from core.testing import Result

TEST = {
    "id": "K.8",
    "phase": "K",
    "title": "Throughput probe Π_meas",
    "theory_part": "instrument for F3 (throughput Π, later Z = I/Π); no physical conclusion",
    "question": "Does a tagged probe measure the throughput T_out/T_in without dependence on its size and without "
                "shifting the background?",
    "method": "a tagged component routed with the acceptance ratios of the common routing; an embedded source and "
              "homogeneous driving; two probe sizes",
    "prediction": "balance of the tagged component to round-off; insensitivity to size; negligible shift of the "
                  "background",
    "criterion": "balance <= balance_tol; sensitivity to size <= size_tol; background shift <= shift_tol, for both "
                 "drivings",
    "null_model": "—",
    "level": "calibration",
    "meters": ["K.6"],
    "depends_on": [],
    "old_id": "E55",
    "params": {"balance_tol": 1e-9, "size_tol": 1e-3, "shift_tol": 1e-3, "probe_rel": 1e-6, "probe_factor": 10.0,
               "N": 16000, "T": 4.0, "R": 0.5, "R_abs": 0.4, "a": 0.05, "C_source": 0.5, "q_homog": 0.9},
    "parts": {"d3": {"D": 3}},
    "quick": {"N": 3000, "T": 2.0},
    "estimate": "1-2 min",
}


def run(p, ctx):
    D = int(p["D"])
    S = sprinkle(Cylinder(D, 0.0, float(p["T"]), float(p["R"])), int(p["N"]), derive_seed(ctx.seed, "K.8", D))
    ptr, idx = links_csr(S.X, S.rho)
    r = S.r()
    absorb = r > float(p["R_abs"])
    cases = {}
    src = (r < float(p["a"])).astype(float)
    cases["inserted source"] = (src, float(p["C_source"]))
    s1 = np.ones(S.n)
    F0 = route_capacity(ptr, idx, s1, np.inf, absorb)["offered"]
    cases["homogeneous driving"] = (s1, float(np.quantile(F0[~absorb], float(p["q_homog"]))))
    lines, data, ok = [f"{D-1}+1 N={p['N']}"], {}, True
    for name, (s, C) in cases.items():
        rc = route_capacity(ptr, idx, s, C, absorb)
        sat = float((rc["stored"] > 0).mean())
        pc = probe_check(ptr, idx, s, C, absorb, float(p["probe_rel"]), float(p["probe_factor"]))
        good = (pc["marked_balance"] <= p["balance_tol"] and pc["size_sensitivity"] <= p["size_tol"]
                and pc["background_shift"] <= p["shift_tol"])
        ok &= good
        data[name] = {k: v for k, v in pc.items() if k != "pi"} | {"C": C, "saturated": sat,
                                                                  "pi_min": float(np.min(pc["pi"]))}
        lines.append(f"  {name}: C = {C:.4g}, saturated {sat:.3f}, min Pi_meas {data[name]['pi_min']:.3f}; balance "
                     f"{pc['marked_balance']:.1e}, size sensitivity {pc['size_sensitivity']:.1e}, background shift "
                     f"{pc['background_shift']:.1e} -> {'OK' if good else 'FAILS'}")
    return Result("PASSED" if ok else "FAILED", lines, data)
