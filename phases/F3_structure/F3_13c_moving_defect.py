"""
F3.13c -- A moving hedgehog: does the locked defect keep its charge and its cost when it moves, like a rest mass?
=================================================================================================================

Why this test exists (31 Oct 2026, after F3.13b2 PASSED; the author's choice of the next step)
--------------------------------------------------------------------------------------------
F3.13b2: a hedgehog of the carried orientation is locked under the alignment dynamics, and its core keeps a cost that a
twist without charge loses; the core re-forms after an erase with about the same cost. The author's candidate: rest mass
as the cost of a configuration. A rest mass has two properties a moving body must show:
  1. it belongs to the body, not to its state of motion: the cost of the defect read in its own rest frame is the same
     whether the defect rests or moves (the principle of relativity for the defect; the substrate is Lorentz invariant
     and the rules are written covariantly, but the transport, the merging, the alignment and the discreteness are a
     dynamics, and only a run shows that the locked structure survives motion with the same cost);
  2. inertia: the defect moves with the flow that carries it and keeps its charge while moving.
If the cost is carried with the flow like the energy density of dust, its energy in the laboratory is gamma times the
rest value (E = gamma m); this follows from 1 and is reported, not tested separately.

Model (locked before any run; everything not named is F3.13b2 with kappa = 0.3)
------------------------------------------------------------------------------
3+1 Poisson sprinkling, density 1, box t [0, 20], x [-7, 7 + 0.6 * 20], y, z [-7, 7] (one box for every velocity, so
that the arms of one sprinkling share the substrate). Source frame: boost by beta along x. Sources: events with
source-frame time 0 <= t_s < 2 and source-frame distance r_s < 5; unit flow with 4-velocity (1, 0, 0, 0) and
orientation (0, n_hat) in the source frame (n_hat of the configuration), both boosted to the laboratory. Transport,
merging and alignment exactly F3.13b (kappa = 0.3, rest-frame lag <= 2): written covariantly, they take the source
frame from the flow itself. One correction (found by the frame check of the first quick run, before any full run): F3.13b
pre-selected the neighbours by a window of coordinate time L_g, exact only for flows at rest in the coordinates; here
all processed carrying events are tested by the covariant conditions alone (identical to F3.13b at rest). Arms: defect at beta = 0, 0.3, 0.6; vacuum at beta = 0.6.
Readout in the source frame of each arm (coordinates, 4-velocities and orientations boosted by beta): slabs of
source-frame time of width 2 from 4 to 14 (the source-frame time that the box carries at beta = 0.6); per slab the
flow-weighted centroid, the degree on shells R = 2, 3, 4 (F3.13a readability rule), and the cost C = (sum of e over the
carrying events of the slab) / 2. Per sprinkling and arm: C_mean = the mean of C over slabs 6 <= t_s < 14; the ratio
r(beta) = C_mean(beta) / C_mean(0) on the same sprinkling; the centroid velocity v_c in the source frame (slope of the
centroid x over the slabs). Reported: the laboratory energy E(beta) = sum over carrying events with laboratory time in
[8, 12) of e (u^0)^2 / 4 against gamma E(0), and the core cost.
Frame check (first 4 sprinklings): the defect at beta = 0.6 computed again with all coordinates and vectors boosted by
beta' = 0.6 must hold flow at the same events and give the same costs (relative 1e-9).

Criteria (stated before any run)
--------------------------------
INSTRUMENT_FAILED if: vacuum (beta = 0.6) reads 0 on < 95 % of its readable shells; a defect arm has fewer than 50
readable shells; the frame check fails.
PASSED: for beta = 0.3 and 0.6 the defect reads +1 on >= 95 % of its readable shells, moves with its source
(|v_c| + 2 SE <= beta / 4) and keeps its cost: |r(beta) - 1| + 2 SE <= 0.15.
FAILED: for beta = 0.3 or 0.6 the defect reads +1 on < 80 % of its readable shells, or |v_c| - 2 SE > beta / 4, or
|r(beta) - 1| - 2 SE > 0.15.
UNDECIDED otherwise.
Development (31 Oct, before any full run; decides nothing). One full-size sprinkling (101 291 events, 1 135 s on one
core with the frame check): every readable shell read +1 (12 / 10 / 5 of 15 at beta = 0 / 0.3 / 0.6), vacuum 0 on 15 of
15; the centroid stays put in the source frame (x_c 0.68-1.11 at beta = 0.6); cost C per slab 7.6 7.0 6.2 5.6 4.5
(beta 0), 6.1 7.2 6.3 5.9 5.7 (0.3), 7.2 7.6 7.8 7.8 4.1 (0.6), i.e. r = 1.08 and 1.17 in this one row; frame check
identical. No parameter changed after it.

Reading declared now. PASSED: the locked defect is carried with its charge and its cost by a moving flow; its cost is a
property of the defect, not of its motion, as a rest mass must be (and its laboratory energy grows as gamma if the cost
is carried like dust). FAILED: motion changes the defect (the dynamics is not covariant in practice, or the moving flow
does not carry the configuration); the cost is not a rest mass. Limits: the cost has no unit yet (no conversion to
mass); cold flow, no field; one kappa; the hedgehog's cost still grows with time (Derrick), so rest and motion are
compared at the same source-frame time.
"""

from __future__ import annotations

import math
import time

import numpy as np

from core.order.links import links_csr
from core.order.sprinkle import Box, derive_seed, sprinkle_poisson
from core.stats.parallel import parallel_map
from core.stats.stats import boot_se
from core.testing import Result
from phases.F3_structure.F3_06c_hot_source import boost
from phases.F3_structure.F3_08_retarded_field import children_csr, mink
from phases.F3_structure.F3_13a_orientation_defect import EXPECT, orient, shell_value
from phases.F3_structure.F3_13b_orientation_locking import _proj_unit, transport_align  # noqa: F401

TEST = {
    "id": "F3.13c",
    "phase": "F3",
    "title": "A moving hedgehog: does the locked defect keep its charge and its cost when it moves, like a rest mass?",
    "theory_part": "F3 emergence of structure; the cost of a locked configuration as a rest mass (the author's candidate); "
                   "inertia of the defect",
    "question": "Carried by a flow moving with beta = 0.3 and 0.6, does the hedgehog of F3.13b2 keep its degree +1, move "
                "with its source and keep, in its own rest frame, the cost it has at rest?",
    "method": "3+1, one box for all velocities (t [0, 20], x [-7, 19], y, z [-7, 7]); F3.13b transport and alignment "
              "(kappa = 0.3); sources and orientations given in the source frame and boosted; readout in the source frame: "
              "degree on shells R = 2, 3, 4, centroid velocity, cost per source-time slab; ratio of the costs moving / at rest "
              "on the same sprinkling; frame check by a boost of 0.6",
    "prediction": "the defect keeps +1, moves with its source and keeps its rest-frame cost (r = 1)",
    "criterion": "INSTRUMENT: vacuum 0 >= 95 %, >= 50 readable shells per defect arm, frame check; PASSED (beta 0.3 and 0.6): "
                 "+1 >= 95 %, |v_c| + 2 SE <= beta/4, |r - 1| + 2 SE <= 0.15; FAILED: +1 < 80 %, |v_c| - 2 SE > beta/4 or "
                 "|r - 1| - 2 SE > 0.15; otherwise UNDECIDED",
    "null_model": "the defect at rest on the same sprinkling (r = 1); vacuum",
    "level": "mechanism",
    "meters": [],
    "depends_on": ["F3.13b2", "F3.13a", "F3.6"],
    "params": {"T": 20.0, "half": 7.0, "betas": [0.0, 0.3, 0.6], "r_src": 5.0, "t_src": 2.0, "slab": 2.0,
               "t_read": 4.0, "t_end": 14.0, "t_cost": 6.0, "shells": [2.0, 3.0, 4.0], "kappa": 0.3, "L_g": 2.0,
               "t_erase": 1e9, "a_erase": 0.0, "lab_window": [8.0, 12.0], "r_core": 1.5, "beta_check": 0.6,
               "n_frame": 4, "n_seeds": 16, "n_se": 2.0, "boot_B": 4000, "pass_frac": 0.95, "fail_frac": 0.8,
               "min_readable": 50, "tol_r": 0.15, "core_eps": 1e-6},
    "parts": {"main": {}},
    "quick": {"n_seeds": 2, "n_frame": 1, "T": 12.0, "half": 5.0, "r_src": 3.0, "t_end": 8.0, "t_cost": 4.0,
              "lab_window": [5.0, 8.0], "boot_B": 200, "min_readable": 1, "shells": [2.0, 3.0]},
    "estimate": "~20 min per order on one core (~100 000 events in 3+1, four arms, the frame check on 4); 16 orders; "
                "memory ~2 GB per process: --workers 6 on 16 GB, 10 on 32 GB",
}


# ============================================================================================== transport

def transport_align_cov(X, cptr, cidx, s, U0, N0, kappa, p, erase_seed=None):
    """transport_align of F3.13b with the neighbourhood searched among all processed carrying events by the covariant
    conditions only (timelike past, rest-frame lag <= L_g); F3.13b pre-selected them by a window of coordinate time L_g,
    exact only for flows at rest in the coordinates (all F3.13b had). Identical to F3.13b for flows at rest."""
    N, D = X.shape
    L_g, eps = float(p["L_g"]), float(p["core_eps"])
    te, ae = float(p["t_erase"]), float(p["a_erase"])
    Pin = np.zeros((N, D))
    nlist = [[] for _ in range(N)]
    nin = np.zeros(N, dtype=np.int64)
    Phold = np.zeros((N, D))
    Nhold = np.full((N, D), np.nan)
    cost = np.full(N, np.nan)
    core = np.zeros(N, dtype=bool)
    car_x, car_n, car_F, ncar = np.zeros((1024, D)), np.zeros((1024, D)), np.zeros(1024), 0   # processed carrying events
    for v in range(N):
        if nin[v] == 0 and s[v] <= 0.0:
            continue
        P = Pin[v] + s[v] * U0[v]
        F = math.sqrt(max(float(mink(P, P)), 1e-300))
        u = P / F
        items = nlist[v] + ([(s[v], N0[v])] if s[v] > 0.0 else [])
        acc = np.zeros(D)
        tot = 0.0
        for Fi, ni in items:
            acc += Fi * (ni - mink(ni, u) * u)
            tot += Fi
        Lc = math.sqrt(max(-float(mink(acc, acc)), 0.0))
        Phold[v] = P
        n_v = None
        if Lc >= eps * tot:
            n_v = acc / Lc
        else:
            core[v] = True
        t_v = float(X[v, 0])
        if n_v is not None and erase_seed is not None and te <= t_v < te + 1.0 and np.linalg.norm(X[v, 1:]) < ae:
            g = np.random.default_rng(derive_seed(erase_seed, v)).normal(size=D)
            n_v = _proj_unit(g, u)                      # a random orientation of the rest space (isotropic there)
        elif n_v is not None and ncar:
            W = car_x[:ncar]
            d = X[v][None, :] - W
            lag = mink(d, u[None, :])
            ok = (d[:, 0] > 0.0) & (mink(d, d) > 0.0) & (lag <= L_g)
            if ok.any():
                Nw = car_n[:ncar][ok]
                Fw = car_F[:ncar][ok]
                m = (Fw[:, None] * (Nw - mink(Nw, u[None, :])[:, None] * u[None, :])).sum(axis=0)
                mh = _proj_unit(m, u)
                if mh is not None:
                    cost[v] = 1.0 + float(mink(n_v, mh))
                    if kappa > 0.0:
                        nn = _proj_unit((1.0 - kappa) * n_v + kappa * mh, u)
                        if nn is not None:
                            n_v = nn
        if n_v is not None:
            Nhold[v] = n_v
            if ncar == car_x.shape[0]:                      # grow the carrier arrays
                car_x = np.vstack([car_x, np.zeros_like(car_x)])
                car_n = np.vstack([car_n, np.zeros_like(car_n)])
                car_F = np.concatenate([car_F, np.zeros_like(car_F)])
            car_x[ncar], car_n[ncar], car_F[ncar] = X[v], n_v, F
            ncar += 1
        a, b = cptr[v], cptr[v + 1]
        if b > a:
            ks = cidx[a:b]
            B = X[ks] - X[v]
            tau = np.sqrt(np.maximum(mink(B, B), 1e-300))
            k = int(ks[int(np.argmin(mink(u[None, :], B) / tau))])
            Pin[k] += P
            nin[k] += 1
            if n_v is not None:
                nlist[k].append((F, n_v))
    return Phold, Nhold, core, cost



# ============================================================================================== one sprinkling

def _source(X, beta, p):
    """Source events (source-frame time and radius) and their unit flows and orientations, boosted to the laboratory."""
    Xs = boost(X, [beta, 0.0, 0.0])
    src = np.flatnonzero((Xs[:, 0] >= 0.0) & (Xs[:, 0] < float(p["t_src"]))
                         & (np.linalg.norm(Xs[:, 1:], axis=1) < float(p["r_src"])))
    return Xs, src


def run_arm(X, cptr, cidx, beta, cfg, p, src_X=None, bc=0.0):
    """One arm: the transport of F3.13b with sources at rest in the frame moving with beta. For the frame check, X are
    the coordinates in the checked frame (boosted by bc), src_X the laboratory coordinates in which the sources are
    chosen, and the initial vectors are boosted by bc as well."""
    N = len(X)
    Xsel = X if src_X is None else src_X
    Xs, src = _source(Xsel, beta, p)
    s = np.zeros(N)
    s[src] = 1.0
    u_rest = np.zeros((src.size, 4))
    u_rest[:, 0] = 1.0
    n_rest = np.zeros((src.size, 4))
    n_rest[:, 1:] = orient(cfg, Xs[src, 1:], float(p["r_src"]))
    U0 = np.zeros((N, 4))
    U0[:, 0] = 1.0
    N0 = np.zeros((N, 4))
    U0[src] = boost(u_rest, [-beta, 0.0, 0.0])
    N0[src] = boost(n_rest, [-beta, 0.0, 0.0])
    if bc != 0.0:
        U0 = boost(U0, [bc, 0.0, 0.0])
        N0 = boost(N0, [bc, 0.0, 0.0])
    return transport_align_cov(X, cptr, cidx, s, U0, N0, float(p["kappa"]), p, None), src


def readout(X, Ph, Nh, cost, beta, p):
    """Source-frame readout (module notes)."""
    Xs = boost(X, [beta, 0.0, 0.0])
    held = Ph[:, 0] > 0.0
    Ps = boost(Ph, [beta, 0.0, 0.0])
    fin = np.isfinite(Nh[:, 0])
    Ns = np.full_like(Nh, np.nan)
    Ns[fin] = boost(Nh[fin], [beta, 0.0, 0.0])
    sl = float(p["slab"])
    edges = np.arange(float(p["t_read"]), float(p["t_end"]) + 1e-9, sl)
    out = {"t": [], "vals": [], "C": [], "core": [], "xc": []}
    for a, b in zip(edges[:-1], edges[1:]):
        m = held & (Xs[:, 0] >= a) & (Xs[:, 0] < b)
        tm = 0.5 * (a + b)
        out["t"].append(tm)
        if not m.any():
            out["C"].append(float("nan"))
            out["core"].append(float("nan"))
            out["xc"].append(float("nan"))
            continue
        w = Ps[m, 0]
        c = (w[:, None] * Xs[m, 1:]).sum(axis=0) / w.sum()
        out["xc"].append(float(c[0]))
        okn = m & fin
        rel = Xs[okn, 1:] - c
        rr = np.linalg.norm(rel, axis=1)
        nsp = Ns[okn, 1:] / np.maximum(np.linalg.norm(Ns[okn, 1:], axis=1), 1e-300)[:, None]
        for Rs in p["shells"]:
            sel = np.abs(rr - Rs) < 0.75
            v, ok = shell_value(rel[sel], nsp[sel], 3)
            out["vals"].append((tm, float(Rs), v, ok))
        cc = cost[okn]
        f = np.isfinite(cc)
        out["C"].append(float(cc[f].sum()) / sl)
        inc = (rr < float(p["r_core"])) & f
        out["core"].append(float(cc[inc].mean()) if inc.any() else float("nan"))
    # laboratory energy in a window of laboratory time
    lw = p["lab_window"]
    m = held & fin & np.isfinite(cost) & (X[:, 0] >= lw[0]) & (X[:, 0] < lw[1])
    u0 = Ph[m, 0] / np.sqrt(np.maximum(mink(Ph[m], Ph[m]), 1e-300))
    out["E_lab"] = float((cost[m] * u0 ** 2).sum()) / (lw[1] - lw[0])
    return out


def task(p, seed, frame):
    rng = np.random.default_rng(seed)
    T, H = float(p["T"]), float(p["half"])
    bmax = max(float(b) for b in p["betas"])
    X = sprinkle_poisson(Box(4, [0.0, -H, -H, -H], [T, H + bmax * T, H, H]), 1.0, rng)
    N = len(X)
    ptr, idx = links_csr(X, 1.0)
    cptr, cidx = children_csr(ptr, idx, N)
    out = {"N": N, "arms": {}}
    arms = [(f"defect_{float(b):g}", "defect", float(b)) for b in p["betas"]] + [(f"vacuum_{bmax:g}", "vacuum", bmax)]
    for key, cfg, beta in arms:
        (Ph, Nh, core, cost), src = run_arm(X, cptr, cidx, beta, cfg, p)
        rec = readout(X, Ph, Nh, cost, beta, p)
        rec["sources"] = int(src.size)
        rec["cores"] = int(core.sum())
        if frame and cfg == "defect" and beta == bmax:
            bc = float(p["beta_check"])
            Xb = boost(X, [bc, 0.0, 0.0])
            (Phb, Nhb, coreb, costb), srcb = run_arm(Xb, cptr, cidx, beta, cfg, p, src_X=X, bc=bc)
            same_held = bool(np.array_equal(Ph[:, 0] > 0.0, Phb[:, 0] > 0.0)) and np.array_equal(src, srcb)
            f1, f2 = np.isfinite(cost), np.isfinite(costb)
            same_cost = bool(np.array_equal(f1, f2) and np.allclose(cost[f1], costb[f2], rtol=1e-9, atol=1e-12))
            rec["frame"] = same_held and same_cost
        out["arms"][key] = rec
    return out


# ============================================================================================== run

def _frac(vals, expect):
    rd = [v for *_, v, ok in vals if ok]
    match = sum(1 for v in rd if abs(v - expect) < 1e-6)
    return (match / len(rd) if rd else float("nan")), len(rd), len(vals)


def run(p, ctx):
    t0 = time.time()
    n, B, nse = int(p["n_seeds"]), int(p["boot_B"]), float(p["n_se"])
    nfr = int(p["n_frame"])
    jobs = [(p, derive_seed(ctx.seed, "F3.13c", k), k < nfr) for k in range(n)]
    rows = parallel_map(task, jobs, ctx.workers, ctx.progress,
                        lambda i, r_, d, nn: f"[F3.13c] {d}/{nn}: N {r_['N']} ({time.time() - t0:.0f}s)")
    betas = [float(b) for b in p["betas"]]
    bmax = max(betas)
    tgrid = rows[0]["arms"][f"defect_{betas[0]:g}"]["t"]
    lines = [f"3+1 box t [0, {p['T']:g}], x [-{p['half']:g}, {p['half'] + bmax * p['T']:g}], y, z [-{p['half']:g}, "
             f"{p['half']:g}] (~{np.mean([r_['N'] for r_ in rows]):.0f} events); sources r_s < {p['r_src']:g}, t_s < "
             f"{p['t_src']:g} in the frame moving with beta; F3.13b transport and alignment kappa = {p['kappa']:g}; readout "
             f"in the source frame; {n} sprinklings, frame check on {nfr}",
             "        t_s          " + " ".join(f"{t:6.1f}" for t in tgrid)]
    tc = np.asarray(tgrid) >= float(p["t_cost"])
    Cm = {}
    res = {}
    for key in [f"defect_{b:g}" for b in betas] + [f"vacuum_{bmax:g}"]:
        vals = [x for r_ in rows for x in r_["arms"][key]["vals"]]
        cfg = key.split("_")[0]
        f, nr, nt = _frac(vals, EXPECT[cfg])
        C = np.array([r_["arms"][key]["C"] for r_ in rows], float)
        core = np.array([r_["arms"][key]["core"] for r_ in rows], float)
        with np.errstate(all="ignore"):
            Cm[key] = np.nanmean(C[:, tc], axis=1)
            Cser, cser = np.nanmean(C, axis=0), np.nanmean(core, axis=0)
        vcs = []
        for r_ in rows:
            xc = np.asarray(r_["arms"][key]["xc"], float)
            ok = np.isfinite(xc)
            vcs.append(float(np.polyfit(np.asarray(tgrid)[ok], xc[ok], 1)[0]) if ok.sum() >= 2 else float("nan"))
        vcs = np.array(vcs)
        vc = boot_se(vcs[np.isfinite(vcs)], B, derive_seed(ctx.seed, "F3.13c-boot", key, "vc"))[:2]
        E = np.array([r_["arms"][key]["E_lab"] for r_ in rows], float)
        fr = [r_["arms"][key]["frame"] for r_ in rows if "frame" in r_["arms"][key]]
        res[key] = {"f": f, "nr": nr, "vc": vc, "E": E, "frame": fr}
        lines.append(f"  {key:11s}: reads {EXPECT[cfg]:+d} on {f:.3f} of {nr} readable (of {nt}); centroid velocity in the "
                     f"source frame {vc[0]:+.4f} +- {vc[1]:.4f}; sources {np.mean([r_['arms'][key]['sources'] for r_ in rows]):.0f}"
                     f"; core events {np.mean([r_['arms'][key]['cores'] for r_ in rows]):.1f}"
                     + (f"; frame check {sum(fr)}/{len(fr)}" if fr else ""))
        lines.append("        cost C      " + " ".join(f"{x:6.1f}" if np.isfinite(x) else "   nan" for x in Cser))
        lines.append("        core cost   " + " ".join(f"{x:6.3f}" if np.isfinite(x) else "   nan" for x in cser))
    base = Cm[f"defect_{betas[0]:g}"]
    rr = {}
    for b in betas[1:]:
        r = Cm[f"defect_{b:g}"] / base
        r = r[np.isfinite(r)]
        rr[b] = boot_se(r, B, derive_seed(ctx.seed, "F3.13c-boot", b, "ratio"))[:2]
        g = 1.0 / math.sqrt(1.0 - b * b)
        e = res[f"defect_{b:g}"]["E"] / (g * res[f"defect_{betas[0]:g}"]["E"])
        e = e[np.isfinite(e)]
        eb = boot_se(e, B, derive_seed(ctx.seed, "F3.13c-boot", b, "E"))[:2]
        lines.append(f"  beta {b:g}: rest-frame cost ratio r = {rr[b][0]:.3f} +- {rr[b][1]:.3f}; reported: laboratory "
                     f"energy / (gamma E(0)) = {eb[0]:.3f} +- {eb[1]:.3f} (gamma {g:.3f})")
    pf, ff, mr, tol = float(p["pass_frac"]), float(p["fail_frac"]), int(p["min_readable"]), float(p["tol_r"])
    fr_all = res[f"defect_{bmax:g}"]["frame"]
    inst = {f"vacuum reads 0 on >= {pf:.0%}": res[f"vacuum_{bmax:g}"]["f"] >= pf,
            f">= {mr} readable shells per defect arm": all(res[f"defect_{b:g}"]["nr"] >= mr for b in betas),
            "frame check": bool(fr_all) and all(fr_all)}
    lines.append("instrument: " + "; ".join(f"{k}: {'OK' if v else 'FAILS'}" for k, v in inst.items()))
    mov = betas[1:]
    ok_pass = all(res[f"defect_{b:g}"]["f"] >= pf and abs(res[f"defect_{b:g}"]["vc"][0]) + nse * res[f"defect_{b:g}"]["vc"][1]
                  <= b / 4.0 and abs(rr[b][0] - 1.0) + nse * rr[b][1] <= tol for b in mov)
    bad = any(res[f"defect_{b:g}"]["f"] < ff or abs(res[f"defect_{b:g}"]["vc"][0]) - nse * res[f"defect_{b:g}"]["vc"][1]
              > b / 4.0 or abs(rr[b][0] - 1.0) - nse * rr[b][1] > tol for b in mov)
    if not all(inst.values()):
        status = "INSTRUMENT_FAILED"
    elif ok_pass:
        status = "PASSED"
    elif bad:
        status = "FAILED"
    else:
        status = "UNDECIDED"
    return Result(status, lines, {"res": {k: {kk: (vv.tolist() if isinstance(vv, np.ndarray) else vv) for kk, vv in v.items()}
                                          for k, v in res.items()}, "ratio": rr, "rows": rows})
