"""
Candidates for the local energy contribution E_v of the field on the order (development; decides nothing)
=====================================================================================================

Why (10 Nov 2026; docs/derivations/field_and_throughput.md, section 5; decision 15)
------------------------------------------------------------------------------------
The chain phi -> energy -> source of the response of the order -> clocks, cones, density -> resistance needs the load
of an event defined on the order: sigma_v = E_v[phi] / C_v. This diagnostic tries candidates on the kink of F3.14a
(and on the curved order of F5.2a2) and reads whether they are local, vanish in the vacuum, sum over a world tube to the
mass of the kink times its proper time, follow the kink, and count correctly in a nonuniform order.

The natural frame-free target (continuum, 1+1). The Lagrangian density L = (phi_t^2 - phi_x^2)/2 - V(phi) is a Lorentz
scalar; on a soliton, integrated over a region crossed by its world tube, int L d^2x = -M * (proper time of the kink in
the region). This is the action of a flow, S = -sum m dtau (Axiom of the action of a flow). The per-event contribution
l_v = L(v) / rho therefore sums to -M dtau, in any frame, and in a nonuniform order (density rho Omega^2) to
-M Omega dt = -M dtau automatically.

Candidates (per event; the sums are compared with -M dtau, M = 8 sqrt(mu)):
  V2   -2 V(phi_v) / rho. On a static soliton int (phi_x^2/2) dx = int V dx (virial / Derrick), so int L = -2 int V on the
       kink, in every frame (V is a scalar). Purely local, periodic in phi (relative phase), zero in the vacuum. It is NOT
       an energy for radiation (a wave has L ~ 0 on average but V > 0).
  grad [ (eps rho / 2) sum_{y < v} K_eps(n_vy) (phi_y - phi_v)^2 - V(phi_v) ] / rho: the kinetic term from the smeared BD
       operator in its difference form (shift invariant). Derivation: (1/2) phi B phi = (1/2) eps rho sum K (phi_y -
       phi_x)^2 + (1/4) B(phi^2) + O(1 + sum K); the middle term is a total derivative and is dropped, the last is noise.
  phiB [ (1/2) phi_v (B_eps phi)_v - V(phi_v) ] / rho: the kinetic term in its original form (not shift invariant).
The same estimators are applied to the exact continuum field at the same events (separates the estimator from the
dynamics). Vacuum check: the same sums over a region of the plateau away from the kink, per unit area.
Radiation check: a small plane wave (amplitude 0.3, k = 1): the time-averaged L is 0 in the continuum.

BD: the smeared 2D kernel of core.operators.bd, eps = 0.1 (nonlocality scale (eps rho)^(-1/2)), interval counts n_vy
exact by dominance counting among the candidates y with (u_v - u_y)(v_v - v_y) <= (1.6 n_cut + 60) / rho.

Run from the root of the laboratory (numba required):
    python -m diag.load_candidates --quick
    python -m diag.load_candidates --seeds 8 --workers 8
    python -m diag.load_candidates --seeds 16 --workers 8 --base 1     # new sprinklings (independent check)
"""

from __future__ import annotations

import argparse
import json
import math
import os
import time

import numpy as np

from core.operators.bd import smeared_kernel, kernel_cutoff_n
from core.order.sprinkle import derive_seed
from core.stats.parallel import parallel_map, default_workers
from phases.F3_structure import field1d as F
from phases.F5_gravity import field1d_curved as C

_nb = F._nb
EPS = 0.1
T_A, T_B = 4.0, 6.0          # the slab over which the sums are taken


if _nb is not None:
    @_nb.njit(cache=True)
    def _bd_sums(U, V, phi, ev_idx, rho, ncut, Ktab):
        """For each evaluated event: S0 = sum K, S1 = sum K phi_y, S2 = sum K (phi_y - phi_v)^2 over y < v."""
        N = U.size
        m_ev = ev_idx.size
        out = np.zeros((m_ev, 3))
        cfac = (ncut * 1.6 + 60.0) / rho
        cand = np.empty(N, np.int64)
        bit = np.zeros(N + 1, np.int64)
        for e in range(m_ev):
            i = ev_idx[e]
            m = 0
            for j in range(N):
                du = U[i] - U[j]
                dv = V[i] - V[j]
                if du > 0.0 and dv > 0.0 and du * dv <= cfac:
                    cand[m] = j
                    m += 1
            c = cand[:m]
            ou = np.argsort(-U[c])
            rk = np.argsort(np.argsort(V[c])) + 1
            for k in range(m + 1):
                bit[k] = 0
            s0 = 0.0
            s1 = 0.0
            s2 = 0.0
            for a in range(m):
                idx = ou[a]
                r = rk[idx]
                s = 0
                q = r
                while q > 0:
                    s += bit[q]
                    q -= q & (-q)
                n = a - s
                if n < ncut:
                    kk = Ktab[n]
                    py = phi[c[idx]]
                    s0 += kk
                    s1 += kk * py
                    s2 += kk * (py - phi[i]) ** 2
                q = r
                while q <= m:
                    bit[q] += 1
                    q += q & (-q)
            out[e, 0] = s0
            out[e, 1] = s1
            out[e, 2] = s2
        return out


def candidates(t, x, phi, mu, rho, sel):
    ncut = kernel_cutoff_n(EPS, 2, 1e-12)
    Ktab = smeared_kernel(np.arange(ncut), 2, EPS)
    U, V = (t + x) / math.sqrt(2.0), (t - x) / math.sqrt(2.0)
    ev = np.nonzero(sel)[0].astype(np.int64)
    S = _bd_sums(U, V, phi, ev, float(rho), int(ncut), Ktab)
    ph = phi[ev]
    Vv = mu * (1.0 - np.cos(ph))
    Bphi = -2.0 * EPS * rho * (ph + S[:, 1])
    l_V2 = -2.0 * Vv / rho
    l_grad = (0.5 * EPS * rho * S[:, 2] - Vv) / rho
    l_phiB = (0.5 * ph * Bphi - Vv) / rho
    return {"V2": l_V2, "grad": l_grad, "phiB": l_phiB, "sumK": S[:, 0]}, ev


def task(rho, seed, arm, quick):
    """arm: rest / v0.3 / v0.6 (mu 1), rest4 (mu 4), curved (F5.2a2 order, mu 1 at rest), wave (plane wave)."""
    X, T = (14.0, 7.0)
    t0 = time.time()
    mu = 4.0 if arm == "rest4" else 1.0
    if arm == "curved":
        t, x = C.sprinkle_curved2(rho, X, T, C.G, seed)
        phi = F.evolve(t, x, C.free_wave(t, x, mu, C.G), rho, mu=mu)
        ex = None
        v = 0.0
    elif arm == "wave":
        t, x = F.sprinkle_trapezoid(rho, X, T, seed)
        modes = [(0.3, 1.0, 0.0)]
        phi = F.evolve(t, x, F.modes_free(t, x, modes), rho, mu=mu)
        ex = F.modes_exact(t, x, modes)
        v = 0.0
    else:
        v = 0.0 if arm.startswith("rest") else float(arm[1:])
        t, x = F.sprinkle_trapezoid(rho, X, T, seed)
        phi = F.evolve(t, x, F.kink_free(t, x, v, mu), rho, mu=mu)
        ex = F.kink_exact(t, x, v, mu)
    xm = v * 0.5 * (T_A + T_B)
    slab = (t >= T_A) & (t < T_B) & (np.abs(x - xm) < 6.0)
    vac = (t >= T_A) & (t < T_B) & (x > xm + 6.0) & (x < X - T_B - 0.5)
    res = {"rho": rho, "seed": seed, "arm": arm, "mu": mu, "v": v, "N": int(t.size)}
    for name, field in (("disc", phi), ("exact", ex)):
        if field is None:
            continue
        lv, ev = candidates(t, x, field, mu, rho, slab | vac)
        tt, xx = t[ev], x[ev]
        in_slab = slab[ev]
        in_vac = vac[ev]
        area_vac = float((T_B - T_A) * max(1e-9, (X - T_B - 0.5) - (xm + 6.0)))
        for k in ("V2", "grad", "phiB"):
            res[f"{name}_{k}_slab"] = float(lv[k][in_slab].sum())
            res[f"{name}_{k}_vac_per_area"] = float(lv[k][in_vac].sum() / area_vac)
            w = -lv[k][in_slab]
            res[f"{name}_{k}_centroid"] = float(np.sum(w * xx[in_slab]) / np.sum(w)) if np.sum(w) != 0 else float("nan")
        res[f"{name}_sumK_mean"] = float(lv["sumK"].mean())
    if arm == "curved":
        r = F.fit_kink(t, x, phi, T_A, T_B, -0.6, mu)
        res["xc"] = r["xc"] if r else float("nan")
    res["sec"] = time.time() - t0
    return res


def _ms(a):
    a = np.asarray([v for v in a if np.isfinite(v)], float)
    return (float(a.mean()), float(a.std(ddof=1) / math.sqrt(a.size))) if a.size > 1 else (float("nan"), float("nan"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=8)
    ap.add_argument("--rho", type=float, default=400.0)
    ap.add_argument("--workers", type=int, default=0)
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--base", type=int, default=20261110, help="seed base; another value gives new sprinklings")
    a = ap.parse_args()
    if _nb is None:
        raise SystemExit("numba required")
    rho, ns = (100.0, 2) if a.quick else (a.rho, a.seeds)
    arms = ["rest", "v0.3", "v0.6", "rest4", "curved", "wave"]
    jobs = [(rho, derive_seed(a.base, "load", rho, k), arm, a.quick) for k in range(ns) for arm in arms]
    t0 = time.time()
    res = parallel_map(task, jobs, default_workers(a.workers), print,
                       lambda i, r, d, n: f"[load] {d}/{n}: {r['arm']} ({time.time() - t0:.0f}s)")
    L = [f"rho {rho:g}, eps {EPS} (l_k {(EPS * rho) ** -0.5:.3f}), slab t in [{T_A}, {T_B}), |x - x_kink| < 6; "
         f"target sum = -M dtau"]
    for arm in arms:
        R = [r for r in res if r["arm"] == arm]
        mu, v = R[0]["mu"], R[0]["v"]
        M = 8.0 * math.sqrt(mu)
        dtau = (T_B - T_A) * math.sqrt(1.0 - v * v)
        if arm == "curved":
            om = _ms([math.exp(C.Phi(np.array([r["xc"]]))[0]) for r in R])[0]
            dtau = (T_B - T_A) * om
            L.append(f"--- curved (F5.2a2 order, kink at rest falling; Omega at the centre {om:.3f})")
        elif arm == "wave":
            L.append("--- plane wave (amplitude 0.3, k 1): target time-averaged L = 0; per unit area")
        else:
            L.append(f"--- {arm} (mu {mu:g}, M {M:g}, v {v:g})")
        for name in ("disc", "exact"):
            if f"{name}_V2_slab" not in R[0]:
                continue
            for k in ("V2", "grad", "phiB"):
                s = _ms([r[f"{name}_{k}_slab"] for r in R])
                vac = _ms([r[f"{name}_{k}_vac_per_area"] for r in R])
                cen = _ms([r[f"{name}_{k}_centroid"] for r in R])
                if arm == "wave":
                    L.append(f"  {name:5s} {k:5s}: sum / area {s[0] / (12 * (T_B - T_A)):+.4f}")
                else:
                    L.append(f"  {name:5s} {k:5s}: sum / (-dtau) {-s[0] / dtau:+.4f} +- {s[1] / dtau:.4f} (target M {M:g}); "
                             f"vacuum per area {vac[0]:+.5f} +- {vac[1]:.5f}; centroid {cen[0]:+.3f} (kink "
                             f"{v * 0.5 * (T_A + T_B):+.3f})")
        L.append(f"  mean sum K per event {np.mean([r['disc_sumK_mean'] for r in R]):+.4f} (expectation -1)")
    print("\n".join(L))
    os.makedirs(os.path.join("diag", "out"), exist_ok=True)
    fn = os.path.join("diag", "out", f"load_candidates_{time.strftime('%Y%m%dT%H%M%S', time.gmtime())}.json")
    with open(fn, "w") as fh:
        json.dump({"args": vars(a), "summary": L, "results": res}, fh)
    print(f"raw records: {fn} ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
