"""
F3.12a -- Resistance of the substrate to rerouting: does it come from capacity or from topology?
=============================================================================================

Why this test exists (31 Oct 2026, the author's hypothesis for F3.12)
--------------------------------------------------------------------
A bound structure (a particle) needs something that resists compression without losing energy: F3.8-F3.11 showed
collapse, cooling and plateaus, never a stationary state. The author's proposal: where the cones narrow (throughput
principle, S = Pi^-2, c_eff = c Pi^2) the flow has to pass through a smaller set of causal paths. Three effects are
possible: capacity (events approach their maximal throughput), topology (fewer admissible paths to reroute the flow),
inertia (rearranging an existing configuration of relations costs). The working hypothesis is that the topological cost
of rerouting is the fundamental one:
    I_v = C(H_v -> H'_v | Phi, Q, K0),    Z_v = C / Pi_v,
a candidate microscopic definition, not a derived law. This test is the instrument that tells the mechanisms apart.

Definitions (the author's, 31 Oct)
---------------------------------
Original network H (links E) with the flow Phi; perturbed network H' (links E'). The rerouted flow is the smallest
change of the flow that the new network admits,
    dPhi_reroute = min over Phi' of 1/2 sum_{e in E u E'} |Phi'_e - Phi_e|,
with Phi_e = 0 on new links and Phi'_e = 0 on removed ones, subject to: the balance at every event that is not a sink
(in + source = out), the capacity of every link (0 <= Phi'_e <= K_e), the sources unchanged. Sinks (the top layer and
the side wall) absorb whatever arrives (an outflow boundary, as in the F1/F3 series). This is a linear program (an L1
distance to the original flow), not a shortest-path flow. Normalised: I_flow = dPhi_reroute / Phi_total.
The topological rebuilding is measured separately: created links |E' \\ E| and removed links |E \\ E'| (not combined:
lambda_+ and lambda_- are not known).
Split used for the reading (the laboratory's): I_flow = I_forced + I_cascade.
  I_forced: the flow that MUST move, because it ran on links that vanished or above the new capacity
            (sum over removed links of Phi_e + sum over capped links of max(0, Phi_e - K'_e), over Phi_total);
            a lower bound of I_flow, true by definition, no resistance in it;
  I_cascade = I_flow - I_forced: the change that the rest of the network has to make IN ADDITION, because the remaining
            paths do not suffice. A resistance of the substrate to rerouting would show here.

Model (locked before any run)
-----------------------------
2+1 box t in [0, 10], |x|, |y| <= 10, density 1, links of the flat sprinkling (core.order.links, rho_min = 1).
Source: unit source at every event with t < 1 and r < 2. Sinks: events with t > 9 or max(|x|, |y|) > 9. Original flow:
uniform routing (core.flux.routing.route_uniform, the F1 law) with the sinks absorbing; events without a future link
also absorb (as in route_uniform). Phi_total = the total flux carried by the links of H (sum of Phi_e), so that I_flow
is a fraction. A lens around the axis carries the perturbation, sigma = sigma0 exp(-r^2 / 2 w^2) for 3 <= t <= 7 (zero
elsewhere), w = 2.5: the source lies below it and the flow can pass through it or around it. Three configurations
per sigma0:
  T (topology): the cones narrow, a link is kept iff its speed <= exp(-sigma) at its parent (F7.1a); capacity unlimited.
               Fewer paths, the same capacity.
  K (capacity): the same links as H; the links whose parent lies in the lens (r < 2 w, 3 <= t <= 7) get one capacity K' = q, with q
               chosen per sprinkling so that I_forced(K) = I_forced(T) (I_forced = sum of max(0, Phi_e - q)). The same
               paths, less capacity, the same forced change.
  P (throughput, F7.1b): events added with density e^sigma - 1 (the original events kept, so H is a part of H'), links of
               the new sprinkling, cones narrowed as in T; capacity unlimited. New links appear (E' \\ E non-empty) and
               old ones vanish, also between old events, when a new event falls into their interval. The sinks of H stay
               sinks; added events at the boundary or without a future link are sinks too.
sigma0 = 0.1, 0.2, 0.3; 32 sprinklings.

Criteria (stated before any run)
--------------------------------
Instrument (INSTRUMENT_FAILED if any fails): with no perturbation (sigma0 = 0, capacity unlimited) the program returns
the original flow, I_flow <= 1e-9, in every sprinkling; the balance of every solution at non-sink events <= 1e-9 of
Phi_total; every program solved (no infeasible network).
Classification at sigma0 = 0.2 (the other strengths reported), with SE over sprinklings (bootstrap):
  topology dominates   I_cascade(T) - 2 SE > 2 (I_cascade(K) + 2 SE);
  capacity dominates   I_cascade(K) - 2 SE > 2 (I_cascade(T) + 2 SE);
  no resistance        I_cascade(T) + 2 SE < 0.01 and I_cascade(K) + 2 SE < 0.01 (rerouting costs nothing beyond the
                       forced change);
  composite            otherwise.
Status RECORDED with the classification (an instrument for F3.12b, not a test of a prediction). Reported: I_flow,
I_forced, I_cascade of T, K, P at all strengths; created and removed links of P and T; the classification of P
against K by the same rule.

Approved by the author on 31 Oct before any run (configurations T, K, P and the classification). Changed before any run
after the first quick run (mechanics, not criteria): the perturbation is a lens in time (3 <= t <= 7) instead of a
worldtube, because with the source inside the tube the capacity arm K was infeasible (the source could not emit through
the capped links); events without a future link absorb; Phi_total is the total link flux.

Reading declared now. Topology dominates: the cost of rerouting comes from the arrangement of the paths, as the author
proposes; F3.12b then tests whether a gradient of this cost deflects a carried flow (a force) without losing energy.
Capacity dominates: resistance is the capacity effect; the topological candidate is not supported as the primary one.
No resistance: losing paths alone costs nothing beyond moving the flow that used them; resistance has to come from
somewhere else (capacity in a stronger regime, or the inertial term). Composite: both, F3.12b needs both.
"""

from __future__ import annotations

import math
import time

import numpy as np
import scipy.sparse as sp
from scipy.optimize import linprog

from core.flux.routing import route_uniform
from core.order.links import links_csr
from core.order.sprinkle import derive_seed
from core.stats.parallel import parallel_map
from core.stats.stats import boot_se
from core.testing import Result

TEST = {
    "id": "F3.12a",
    "phase": "F3",
    "title": "Resistance of the substrate to rerouting: does it come from capacity or from topology?",
    "theory_part": "F3 emergence of structure; the impedance Z of the substrate as the cost of rerouting the flow "
                   "(the author's hypothesis of 31 Oct)",
    "question": "When the flow of a source has to pass a region where the cones narrow (fewer paths) or the capacity drops "
                "(same paths), how much of the flow has to change its paths beyond the flow that is forced to move "
                "(I_cascade), and which of the two perturbations makes it larger?",
    "method": "2+1, box t [0, 10], |x| <= 10; uniform routing from a source at the bottom; perturbations in a lens "
              "3 <= t <= 7, sigma = sigma0 exp(-r^2 / 2 w^2), w = 2.5: T narrowed cones, K capacity cut to the same forced change, P "
              "throughput (narrowed cones and added events); the smallest L1 change of the flow by linear programming; "
              "I_flow = I_forced + I_cascade; created and removed links",
    "prediction": "the author's working hypothesis: the topological cost of rerouting dominates (T above K)",
    "criterion": "at sigma0 = 0.2: topology dominates I_c(T) - 2 SE > 2 (I_c(K) + 2 SE); capacity dominates the reverse; "
                 "no resistance both I_c + 2 SE < 0.01; composite otherwise. INSTRUMENT: unperturbed I_flow <= 1e-9, "
                 "balance <= 1e-9, all programs solved. Status RECORDED with the classification.",
    "null_model": "unperturbed network (I_flow = 0 exactly)",
    "level": "calibration",
    "meters": ["K.6"],
    "depends_on": ["F7.1a", "F7.1b", "F1.1"],
    "params": {"T": 10.0, "half": 10.0, "width": 2.5, "sigmas": [0.1, 0.2, 0.3], "decide": 0.2, "window": [3.0, 7.0], "src_r": 2.0,
               "src_t": 1.0, "sink_t": 9.0, "sink_wall": 9.0, "n_seeds": 32, "n_se": 2.0, "boot_B": 4000,
               "zero_tol": 1e-9, "none_tol": 0.01},
    "parts": {"main": {}},
    "quick": {"n_seeds": 2, "T": 6.0, "half": 6.0, "sink_t": 5.0, "sink_wall": 5.0, "sigmas": [0.2], "decide": 0.2,
              "window": [2.0, 4.0], "boot_B": 200},
    "estimate": "minutes (32 sprinklings of ~4 000 events; ten linear programs of ~1e5 variables each)",
}


# ===================================================================================== networks

def sigma_at(Xe, s0, w, win):
    """The perturbation sigma at events Xe (t, x, y): a lens around the axis, active for win[0] <= t <= win[1]."""
    Xe = np.asarray(Xe)
    on = (Xe[..., 0] >= win[0]) & (Xe[..., 0] <= win[1])
    return np.where(on, s0 * np.exp(-(Xe[..., 1:] ** 2).sum(axis=-1) / (2.0 * w * w)), 0.0)


def edge_list(X, ids, keep_speed=None):
    """Links of the flat sprinkling X (time-sorted) as arrays of global ids (u earlier, v later); with keep_speed(s0, w)
    given, only the links whose speed is within exp(-sigma) at the parent are kept (F7.1a)."""
    ptr, idx = links_csr(X, 1.0)
    v = np.repeat(np.arange(len(X)), np.diff(ptr))
    u = idx
    if keep_speed is not None:
        d = X[v] - X[u]
        speed = np.hypot(d[:, 1], d[:, 2]) / d[:, 0]
        ok = speed <= np.exp(-keep_speed(X[u]))
        u, v = u[ok], v[ok]
    return ids[u], ids[v]


def maximal(X):
    """Events without a link to the future in the flat sprinkling (they absorb, as in route_uniform)."""
    ptr, idx = links_csr(X, 1.0)
    return np.bincount(idx, minlength=len(X)) == 0


def original_flow(X, src, sink):
    """Uniform routing (F1 law) on the flat links, sinks absorbing; returns (u, v, Phi) in local indices."""
    ptr, idx = links_csr(X, 1.0)
    _, Phi = route_uniform(ptr, idx, src, absorb=sink)
    v = np.repeat(np.arange(len(X)), np.diff(ptr))
    return idx, v, Phi


def reroute(n_nodes, sink, src, old, new, cap):
    """The smallest L1 change of the flow (linear program).
    old: dict (u, v) -> Phi_e of the original flow; new: list of (u, v) of the new network; cap: array of capacities
    aligned with new (np.inf allowed). Returns (dPhi_reroute, Phi', balance residual, status)."""
    m = len(new)
    nu = np.array([e[0] for e in new], dtype=np.int64)
    nv = np.array([e[1] for e in new], dtype=np.int64)
    phi_old = np.array([old.get(e, 0.0) for e in new])
    common = np.array([e in old for e in new])
    ci = np.flatnonzero(common)
    k = ci.size
    newset = set(new)
    removed = sum(val for e, val in old.items() if val != 0.0 and e not in newset)
    # variables: Phi'_e (m), d_e for common links (k)
    c = np.r_[np.where(common, 0.0, 0.5), np.full(k, 0.5)]
    # balance at non-sink nodes: in - out = -src
    rows = np.flatnonzero(~sink)
    rpos = -np.ones(n_nodes, dtype=np.int64)
    rpos[rows] = np.arange(rows.size)
    ri, cj, vals = [], [], []
    a = rpos[nv] >= 0                                   # inflow at v
    ri.append(rpos[nv][a]); cj.append(np.flatnonzero(a)); vals.append(np.ones(a.sum()))
    b = rpos[nu] >= 0                                   # outflow at u
    ri.append(rpos[nu][b]); cj.append(np.flatnonzero(b)); vals.append(-np.ones(b.sum()))
    A_eq = sp.csr_matrix((np.concatenate(vals), (np.concatenate(ri), np.concatenate(cj))), shape=(rows.size, m + k))
    b_eq = -src[rows]
    # |Phi' - Phi| <= d on common links
    I = np.arange(k)
    A_ub = sp.vstack([sp.csr_matrix((np.r_[np.ones(k), -np.ones(k)], (np.r_[I, I], np.r_[ci, m + I])), shape=(k, m + k)),
                      sp.csr_matrix((np.r_[-np.ones(k), -np.ones(k)], (np.r_[I, I], np.r_[ci, m + I])), shape=(k, m + k))])
    b_ub = np.r_[phi_old[ci], -phi_old[ci]]
    ub = np.where(sink[nu], 0.0, cap)                    # sinks forward nothing
    bounds = list(zip(np.zeros(m), np.where(np.isfinite(ub), ub, None))) + [(0.0, None)] * k
    res = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=bounds, method="highs")
    if res.status != 0:
        return float("nan"), None, float("nan"), res.status
    phi = res.x[:m]
    d = 0.5 * (np.abs(phi[common] - phi_old[common]).sum() + phi[~common].sum() + removed)
    bal = np.abs(A_eq @ res.x - b_eq).max() if rows.size else 0.0
    return float(d), phi, float(bal), 0


# ===================================================================================== one sprinkling

def seed_task(p, seed):
    rng = np.random.default_rng(seed)
    T, Hh, w = float(p["T"]), float(p["half"]), float(p["width"])
    n = rng.poisson(T * (2 * Hh) ** 2)
    X0 = np.column_stack([rng.uniform(0, T, n), rng.uniform(-Hh, Hh, n), rng.uniform(-Hh, Hh, n)])
    X0 = X0[np.argsort(X0[:, 0])]
    N0 = len(X0)
    r0 = np.hypot(X0[:, 1], X0[:, 2])
    src0 = ((X0[:, 0] < p["src_t"]) & (r0 < p["src_r"])).astype(float)
    sink0 = (X0[:, 0] > p["sink_t"]) | (np.abs(X0[:, 1:]).max(axis=1) > p["sink_wall"]) | maximal(X0)
    u, v, Phi = original_flow(X0, src0, sink0)
    total = float(Phi.sum())                             # the total flux carried by the links of H
    old = {(int(a), int(b)): float(f) for a, b, f in zip(u, v, Phi)}
    E = list(old)
    out = {"N": N0, "links": len(E), "total": total, "res": {}}
    # instrument: the unperturbed program returns the original flow
    d0, _, bal0, st0 = reroute(N0, sink0, src0, old, E, np.full(len(E), np.inf))
    out["zero"] = {"I_flow": d0 / total, "bal": bal0 / total, "status": st0}
    ok_all = st0 == 0
    for s0 in p["sigmas"]:
        win = p["window"]
        sig = lambda Xe: sigma_at(Xe, s0, w, win)
        r_ = {}
        # T: narrowed cones
        tu, tv = edge_list(X0, np.arange(N0), sig)
        newT = list(zip(tu.tolist(), tv.tolist()))
        setT = set(newT)
        forcedT = sum(f for e, f in old.items() if e not in setT) / total
        dT_, _, balT, stT = reroute(N0, sink0, src0, old, newT, np.full(len(newT), np.inf))
        r_["T"] = {"I_flow": dT_ / total, "I_forced": forcedT, "I_cascade": dT_ / total - forcedT, "bal": balT / total,
                   "status": stT, "removed": len(E) - len(setT & set(E)), "created": len(setT - set(E))}
        # K: the same links, the capacity in the disc cut to the same forced change
        disc = np.array([np.hypot(*X0[a, 1:]) < 2.0 * w and win[0] <= X0[a, 0] <= win[1] for a, _ in E])
        fl = np.array([old[e] for e in E])
        def forced_at(q):
            return float(np.maximum(fl[disc] - q, 0.0).sum()) / total
        lo, hi = 0.0, float(fl[disc].max()) if disc.any() else 0.0
        for _ in range(80):                              # bisection: forced_at is decreasing in q
            mid = 0.5 * (lo + hi)
            if forced_at(mid) > forcedT:
                lo = mid
            else:
                hi = mid
        q = 0.5 * (lo + hi)
        capK = np.where(disc, q, np.inf)
        forcedK = forced_at(q)
        dK, _, balK, stK = reroute(N0, sink0, src0, old, E, capK)
        r_["K"] = {"I_flow": dK / total, "I_forced": forcedK, "I_cascade": dK / total - forcedK, "bal": balK / total,
                   "status": stK, "q": q}
        # P: throughput, events added with density e^sigma - 1, cones narrowed
        smax = math.exp(s0)
        m_add = rng.poisson((smax - 1.0) * T * (2 * Hh) ** 2)
        Xa = np.column_stack([rng.uniform(0, T, m_add), rng.uniform(-Hh, Hh, m_add), rng.uniform(-Hh, Hh, m_add)])
        Xa = Xa[rng.random(m_add) < (np.exp(sig(Xa)) - 1.0) / (smax - 1.0)]
        XP = np.vstack([X0, Xa])
        idsP = np.arange(len(XP))                        # old events keep their ids 0..N0-1
        o = np.argsort(XP[:, 0], kind="stable")
        pu, pv = edge_list(XP[o], idsP[o], sig)
        newP = list(zip(pu.tolist(), pv.tolist()))
        setP = set(newP)
        NP = len(XP)
        srcP = np.r_[src0, np.zeros(len(Xa))]
        maxP = np.zeros(NP, dtype=bool)
        maxP[idsP[o]] = maximal(XP[o])
        sinkP = np.r_[sink0, (XP[N0:, 0] > p["sink_t"]) | (np.abs(XP[N0:, 1:]).max(axis=1) > p["sink_wall"]) | maxP[N0:]]
        forcedP = sum(f for e, f in old.items() if e not in setP) / total
        dP_, _, balP, stP = reroute(NP, sinkP, srcP, old, newP, np.full(len(newP), np.inf))
        r_["P"] = {"I_flow": dP_ / total, "I_forced": forcedP, "I_cascade": dP_ / total - forcedP, "bal": balP / total,
                   "status": stP, "removed": len(E) - len(setP & set(E)), "created": len(setP - set(E)), "added": len(Xa)}
        ok_all &= stT == 0 and stK == 0 and stP == 0
        out["res"][str(s0)] = r_
    out["ok"] = bool(ok_all)
    return out


# ===================================================================================== the test

def _classify(a, b, nse, tol):
    """a, b: (mean, se) of I_cascade of the 'topology' and the 'capacity' arm."""
    (ma, sa), (mb, sb) = a, b
    if ma - nse * sa > 2.0 * (mb + nse * sb):
        return "topology dominates"
    if mb - nse * sb > 2.0 * (ma + nse * sa):
        return "capacity dominates"
    if ma + nse * sa < tol and mb + nse * sb < tol:
        return "no resistance"
    return "composite"


def run(p, ctx):
    t0 = time.time()
    n, B, nse = int(p["n_seeds"]), int(p["boot_B"]), float(p["n_se"])
    jobs = [(p, derive_seed(ctx.seed, "F3.12a", ctx.part, k)) for k in range(n)]
    rows = parallel_map(seed_task, jobs, ctx.workers, ctx.progress,
                        lambda i, r_, d, nn: f"[F3.12a] {d}/{nn}: N {r_['N']} links {r_['links']} "
                                             f"({time.time() - t0:.0f}s)")
    zero = max(r_["zero"]["I_flow"] for r_ in rows)
    bal = max([r_["zero"]["bal"] for r_ in rows] +
              [r_["res"][k][c]["bal"] for r_ in rows for k in r_["res"] for c in ("T", "K", "P")])
    solved = all(r_["ok"] for r_ in rows)
    lines = [f"2+1 box t [0, {p['T']:g}], |x| <= {p['half']:g} (~{np.mean([r_['N'] for r_ in rows]):.0f} events, "
             f"~{np.mean([r_['links'] for r_ in rows]):.0f} links); source t < {p['src_t']:g}, r < {p['src_r']:g}; "
             f"lens {p['window'][0]:g} <= t <= {p['window'][1]:g}, w = {p['width']:g}; {n} sprinklings",
             f"  unperturbed: max I_flow {zero:.1e}; max balance residual {bal:.1e}; all programs solved: {solved}"]
    summ = {}
    for s0 in p["sigmas"]:
        k = str(s0)
        summ[k] = {}
        for cfg in ("T", "K", "P"):
            ok = [r_["res"][k][cfg] for r_ in rows if r_["res"][k][cfg]["status"] == 0]
            bs = lambda *tag: derive_seed(ctx.seed, "F3.12a-boot", ctx.part, s0, cfg, *tag)
            summ[k][cfg] = {q_: boot_se([o[q_] for o in ok], B, bs(q_))[:2] for q_ in ("I_flow", "I_forced", "I_cascade")}
            for q_ in ("removed", "created", "added"):
                if q_ in ok[0]:
                    summ[k][cfg][q_] = float(np.mean([o[q_] for o in ok]))
        lines.append(f"  sigma0 = {s0:g}:")
        for cfg, lab in (("T", "T topology "), ("K", "K capacity "), ("P", "P throughput")):
            s_ = summ[k][cfg]
            extra = "".join(f", {q_} {s_[q_]:.0f}" for q_ in ("removed", "created", "added") if q_ in s_)
            lines.append(f"    {lab}: I_flow {s_['I_flow'][0]:.4f} +- {s_['I_flow'][1]:.4f}; I_forced {s_['I_forced'][0]:.4f}; "
                         f"I_cascade {s_['I_cascade'][0]:+.4f} +- {s_['I_cascade'][1]:.4f}{extra}")
        lines.append(f"    reported: P against K: {_classify(summ[k]['P']['I_cascade'], summ[k]['K']['I_cascade'], nse, p['none_tol'])}")
    kd = str(p["decide"])
    cls = _classify(summ[kd]["T"]["I_cascade"], summ[kd]["K"]["I_cascade"], nse, float(p["none_tol"])) if kd in summ else "n/a"
    inst = {"unperturbed I_flow <= tol": zero <= p["zero_tol"], "balance <= tol": bal <= p["zero_tol"],
            "all programs solved": solved}
    lines.append("instrument: " + "; ".join(f"{k_}: {'OK' if v_ else 'FAILS'}" for k_, v_ in inst.items()))
    lines.append(f"classification at sigma0 = {p['decide']:g} (T against K): {cls}")
    status = "INSTRUMENT_FAILED" if not all(inst.values()) else "RECORDED"
    return Result(status, lines, {"summary": summ, "classification": cls, "rows": rows})
