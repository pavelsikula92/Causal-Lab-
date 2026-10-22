"""
K.6 -- Identities of flux routing
================================

The routing kernels (core.flux.routing: uniform routing, the E41-A law;
core.flux.capacity: capacity-limited routing) must conserve flux exactly.

  (a) node-set balance and the closed-surface identity: for the region
      |x| < r, t in [t1, t2], sources inside = side net flux + top cap flux
      - bottom cap flux + sinks + discarded, at several radii, with and
      without the outflow boundary and the lossy null (keep < 1);
  (b) no source -> exactly no flux;
  (c) capacity routing: balance sum s = sum dQ + outflow exactly; C = inf
      reproduces the uniform routing exactly;
  (d) settling (reported): the impulse response of the routing decays
      exponentially; its rate sets how long a window must wait (F1.1).
"""

from __future__ import annotations

import numpy as np

from core.order.sprinkle import sprinkle, Cylinder, derive_seed
from core.order.links import links_csr
from core.flux.routing import route_uniform, link_endpoints, node_set_balance, net_side_flux, cap_flux, stock_profile
from core.flux.capacity import route_capacity, balance_residual
from core.stats.stats import fmt
from core.testing import Result

TEST = {
    "id": "K.6",
    "phase": "K",
    "title": "Identity toku (routing a kapacita)",
    "theory_part": "měřidlo pro axiom 3 a 4 (F1.1, F3); bez fyzikálního závěru",
    "question": "Zachovává routing tok přesně, s kapacitou i bez ní?",
    "method": "bilance množiny uzlů a identita uzavřené plochy; nulový zdroj; bilance kapacitního routingu; "
              "C = nekonečno proti uniformnímu routingu; impulzní odezva",
    "prediction": "všechny identity platí na strojovou přesnost",
    "criterion": "nejhorší reziduum identit <= tol; tok bez zdroje přesně 0; bilance kapacity <= tol; "
                 "rozdíl C = nekonečno proti uniformnímu routingu <= tol",
    "null_model": "—",
    "level": "calibration",
    "meters": [],
    "depends_on": [],
    "old_id": "E45a",
    "params": {"tol": 1e-9, "N": 16000, "T": 4.0, "R": 0.5, "R_abs": 0.4, "a": 0.05, "t1": 2.0, "t2": 3.5,
               "radii": [0.1, 0.2, 0.3], "keep_null": 0.9, "C_frac": 0.6, "n_seeds": 2},
    "parts": {"d3": {"D": 3}, "d4": {"D": 4, "N": 24000, "R": 0.4, "R_abs": 0.35}},
    "quick": {"N": 3000, "n_seeds": 1, "T": 2.0, "t1": 0.5, "t2": 1.5},
    "estimate": "několik minut",
}


def run(p, ctx):
    D = int(p["D"])
    worst, zero, cap_bal, cap_inf, rates = 0.0, 0.0, 0.0, 0.0, []
    for rep in range(int(p["n_seeds"])):
        S = sprinkle(Cylinder(D, 0.0, float(p["T"]), float(p["R"])), int(p["N"]), derive_seed(ctx.seed, "K.6", D, rep))
        ptr, idx = links_csr(S.X, S.rho)
        r, t = S.r(), S.X[:, 0]
        src = (r < float(p["a"])).astype(float)
        absorb = r > float(p["R_abs"])
        u, v = link_endpoints(ptr, idx)
        for keep, ab in ((1.0, None), (float(p["keep_null"]), absorb), (1.0, absorb)):
            F, Phi = route_uniform(ptr, idx, src, keep=keep, absorb=ab)
            for rr in p["radii"]:
                ins = (r < rr) & (t >= p["t1"]) & (t <= p["t2"])
                bal = node_set_balance(ptr, idx, src, F, Phi, ins, keep=keep, absorb=ab)
                side = net_side_flux(S.X, u, v, Phi, rr, p["t1"], p["t2"])["net"]
                caps = cap_flux(S.X, u, v, Phi, p["t2"], rr) - cap_flux(S.X, u, v, Phi, p["t1"], rr)
                res2 = bal["sources"] - (side + caps + bal["sinks"] + bal["discarded"])
                worst = max(worst, abs(bal["residual"]) / max(1.0, bal["sources"]), abs(res2) / max(1.0, bal["sources"]))
        F0, Phi0 = route_uniform(ptr, idx, np.zeros(S.n), absorb=absorb)
        zero = max(zero, float(np.abs(Phi0).max()) if Phi0.size else 0.0)
        s1 = np.ones(S.n)
        Fu, Phu = route_uniform(ptr, idx, s1, absorb=absorb)
        ri = route_capacity(ptr, idx, s1, np.inf, absorb)
        cap_inf = max(cap_inf, float(np.abs(ri["offered"] - Fu).max() + np.abs(ri["Phi"] - Phu).max()))
        C = float(np.quantile(ri["offered"], float(p["C_frac"])))
        rc = route_capacity(ptr, idx, s1, C, absorb)
        cap_bal = max(cap_bal, balance_residual(s1, rc))
        imp = src * (t < 0.1)
        if imp.sum() > 0:
            _, Phi_i = route_uniform(ptr, idx, imp, absorb=absorb)
            tg = np.arange(0.5, float(p["T"]) - 0.25 + 1e-9, 0.25)
            Pg = stock_profile(S.X, u, v, Phi_i, tg, [float(p["radii"][-1])])[:, 0] / imp.sum()
            m = (Pg > 1e-6) & (Pg < 0.5)
            if m.sum() >= 3:
                rates.append(float(-np.polyfit(tg[m], np.log(Pg[m]), 1)[0]))
        ctx.progress(f"[K.6] {D-1}+1 seed {rep+1}/{p['n_seeds']}: identity worst {worst:.1e}, capacity balance {cap_bal:.1e}")
    ok = worst <= p["tol"] and zero == 0.0 and cap_bal <= p["tol"] and cap_inf <= p["tol"]
    lines = [f"{D-1}+1 N={p['N']}: closed-surface and node-set identities worst {worst:.1e} (tol {p['tol']:g})",
             f"  no source -> flux max {zero:g}; capacity balance {cap_bal:.1e}; C = inf vs uniform routing {cap_inf:.1e}",
             f"  impulse response decay rate (reported) {fmt(np.mean(rates) if rates else float('nan'))} per unit t"]
    return Result("PASSED" if ok else "FAILED", lines, {"identity_worst": worst, "zero": zero, "capacity_balance": cap_bal,
                                                        "capacity_inf_diff": cap_inf, "decay_rates": rates})
