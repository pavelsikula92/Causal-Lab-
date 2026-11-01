"""
F3.12b -- Capacity as pressure: do streams that are reflected out of overloaded events hold a bundle in its own field?
=====================================================================================================================

Why this test exists (31 Oct 2026, after F3.12a2; the author's line towards a particle)
--------------------------------------------------------------------------------------
A bound structure needs something that resists compression without losing energy. The F3 series found the pull
(F3.8: the retarded field contracts a cold bundle) but no stationary state: capacity without an outward momentum gave
only a plateau (F3.9), internal motion held the bundle for a while but merging turned it into rest mass and it cooled
(F3.11), and flows that meet stick (F4.1). F3.12a2 then showed that the resistance of the substrate to rerouting is
the capacity effect, not the topology of the paths (I_cascade at sigma0 = 0.2: capacity 0.053, topology 0.013).
This test gives capacity the two properties F3.9 lacked:
  1. no merging: every stream keeps its own 4-momentum (and so its rest mass); streams may pass through the same
     event without fusing (the merging of F3.6c/F3.11 is the dissipative channel that cooled F3.11);
  2. an outward momentum: an event that carries more streams than its capacity C reflects every stream that moves
     towards higher load, elastically, in the rest frame of the load at that event (its energy there unchanged, its
     rest mass unchanged). This is the author's "the more something is compressed, the larger the resistance" in the
     form F3.12a2 supports: capacity, acting as a wall that returns momentum instead of absorbing it.
The question: is there a capacity C (and a temperature eta) at which a bundle in its own field keeps a constant width,
between the dispersal without a field and the collapse of a cold bundle?

Model (locked before any run; everything not named here is F3.11 unchanged)
---------------------------------------------------------------------------
2+1 Poisson sprinkling, rho = 1, box t [0, 48], half width 16; links (core.order.links). Source: world tube of radius 4,
time [0, 4), unit weight per source event along its 4-velocity (eta = 0 cold; eta > 0 the F3.6c 'sym' configuration).
Field and bending exactly F3.8 / F3.11: phi(x) = sum_{w < x, tau >= 1} m_w / (2 pi tau) with m_w the rest mass carried
through event w (the sum of the rest masses of the streams there: a scalar source, as F3.11), ln Pi = -lambda phi,
lambda = 0.1 / 44.15 locked; each stream is bent by the transverse gradient of ln Pi over the parents of its event
(rest-frame lag <= 2, at least 4 parents), its child chosen as the child of least rapidity relative to its unbent
direction (F3.6), the bent direction sent on (F3.11).
Streams (new): an event holding n streams p_i forwards each separately; nothing merges.
Capacity reflection (new): if n > C, let U be the unit vector of sum p_i (the rest frame of the load). The gradient of
the load (number of streams) over the parents of the event is fitted in the rest space of U exactly as the gradient of
ln Pi in F3.8 (lag <= 2 in that frame, at least 4 parents); nhat = the unit spacelike vector orthogonal to U pointing
down the gradient (towards lower load); a gradient below 1e-9 counts as none. Every stream with a = -p_i . nhat < 0 (moving towards higher load) is
reflected, p_i -> p_i - 2 a nhat: its rest mass and its energy in the frame U are unchanged, its momentum along nhat
reverses (the momentum goes to the substrate, as the bending momentum goes to the field). No reflection without a
fitted gradient.
Arms: free_0 (cold, no field, the F3.11 transport with merging: the reference of the classification, exactly as
F3.11), field_0 (cold, field, F3.11 transport: the positive control), free_cl_0.06 (streams, no field, no capacity:
the hot collisionless dispersal and the exact 4-momentum balance), cl_{eta} for eta in {0.03, 0.06, 0.1} (streams and
field, no capacity: does removing the merging alone bind?), and cap{C}_{eta} for C in {2, 4, 8}, eta in
{0, 0.03, 0.06, 0.1} (streams, field, capacity). A stream is reflected only if a < -1e-9 times its energy in the frame
U (a component at round-off is a numerical tie, not a direction), and no reflection when the fitted load gradient is
below 1e-9 (a flat load has no direction).

Readout (as F3.11, per slab of width 2 in the source frame)
-----------------------------------------------------------
Entries: events holding streams, weight = the summed energy of their streams in the source frame. R_rms about the
weighted centroid; centroid x; internal spread = energy-weighted mean rapidity of the streams relative to the source
frame; number of carrying events. Slopes of R_rms in W1 [12, 28] and W2 [28, 44]. Reported per arm: mean bending,
reflections per stream, the rest mass that leaves the box against the rest mass emitted (exact for every stream arm:
nothing merges, reflections and bending keep the rest mass).

Classification (exactly F3.8 / F3.11, reference free_0)
-------------------------------------------------------
Per window, tol = s(free_0) / 4: bound |s| + 2 SE <= tol; growing s - 2 SE > tol; collapsing s + 2 SE < -tol;
otherwise unresolved; a width above half / 2 in the window turns bound, unresolved or collapsing into wall-limited.

Parts and criteria (stated before any run)
------------------------------------------
search  RECORDED, at most one candidate: a stream arm with a field (cl_* or cap*) bound in W1 AND W2 with identical
        choices in the frame boosted by beta' = 0.6; the smallest |s_W1| + 2 SE, ties the smaller C, then the smaller
        eta. INSTRUMENT_FAILED if: an arm fails the frame check (done on the first n_frame sprinklings, all arms); the
        rest mass of any stream arm is not conserved to 1e-12 (relative); the 4-momentum balance of free_0 or
        free_cl_0.06 exceeds 1e-12; or the positive control fails (field_0 contracts against free_0 in W1 or W2:
        s(field_0) + 2 SE < s(free_0) in at least one window).
confirm --set arm=<candidate> (one of the search arms), new sprinklings, the source at rest and moving with beta = 0.5:
        PASSED  the arm bound at rest in W1 and W2, bound in motion in W1, inertia |v_c| + 2 SE <= beta/4 (centroid
                velocity in the source frame over [5, 28]), instruments OK;
        FAILED  the arm growing or collapsing at rest in W1 or W2, or |v_c| - 2 SE > beta/4;
        INSTRUMENT_FAILED frame check, rest mass or balances;
        UNDECIDED otherwise.
Development (31 Oct, before any full run; decides nothing). One full-size sprinkling (48 967 events, 188 source
events, all ten arms of the first grid eta in {0, 0.15}, frame check on; 172 s on one core): slopes s1 / s2
  free_0 +0.030 / +0.005; field_0 -0.052 / -0.054 (R 2.8 -> 1.5);
  free_cl_0.15 +0.111 / +0.153; cl_0.15 +0.086 / +0.128;
  cap2 / cap4 / cap8 at eta 0: -0.035 / -0.035 / -0.051 (W1), -0.028 / -0.029 / -0.031 (W2), reflections per stream
    0.41 / 0.38 / 0.31;
  cap2 / cap4 / cap8 at eta 0.15: +0.108 / +0.101 / +0.096 (W1), +0.146 / +0.130 / +0.123 (W2);
  rest mass conserved exactly in all stream arms; frame check identical in all arms except cap2_0.15.
Without merging nothing cools, so eta = 0.15 is far above what the locked field holds (cl_0.15 spreads almost as fast as
free_cl_0.15), while eta = 0 collapses, slower with capacity than without. The temperature grid was therefore moved
before any full run to eta in {0, 0.03, 0.06, 0.1} (free_cl at 0.06), with cl_eta at every hot eta as the control of
"merging was the only obstacle". The frame difference in cap2_0.15 was traced to the sign test of the reflection at
components of the order of round-off, and to an event whose parents all carried the same load: its fitted gradient is
1e-16, a direction made of round-off that differs between frames. Both are now ties: no reflection when the gradient of
the load is below 1e-9, and no reflection of a component below 1e-9 of the energy. Classification, criteria and
instrument unchanged.

Reading declared now. A bound arm: the first stationary structure of the laboratory, held by its own field against a
non-dissipative resistance; capacity then sets its size, and the next step is the relation between size, capacity
and mass (and the 3+1 version: in 2+1 the potential is logarithmic and confines). a cl arm bound alone: the merging was
what destroyed F3.11; capacity is not needed for binding. No bound arm: reflection off overloaded events does not
balance the field in this form; the plateaus and their lengths are reported. Declared limits: 2+1; the bending of the
F3 series is the conformal geometry (gamma = -1); the field is sourced by rest mass only; the reflection gives the
momentum to the substrate (the 4-momentum of matter is conserved only without field and capacity).
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
from phases.F3_structure.F3_06c_hot_source import boost, source_config
from phases.F3_structure.F3_08_retarded_field import (RetardedField, boost_towards, children_csr, classify, mink,
                                                       rest_basis, slope)
from phases.F3_structure.F3_11_hot_bundle_field import transport as transport_merge

TEST = {
    "id": "F3.12b",
    "phase": "F3",
    "title": "Capacity as pressure: do streams reflected out of overloaded events hold a bundle in its own field?",
    "theory_part": "binding of matter: the retarded field (F3.8) against the capacity of the substrate acting as an "
                   "elastic wall (F3.12a2: the resistance to rerouting is the capacity effect)",
    "question": "With streams that never merge and are reflected elastically out of events that carry more than C "
                "streams, is there a capacity C (and a temperature eta) at which a bundle in its own retarded field "
                "keeps a constant width?",
    "method": "2+1 as F3.11 (box t [0, 48], half 16, source radius 4, field and bending of F3.8 locked); streams carried "
              "separately; at an event with more than C streams every stream moving up the load gradient is reflected "
              "in the rest frame of the load; arms free_0, field_0 (F3.11 transport), free_cl_0.06, cl_{0.03,0.06,0.1}, "
              "cap{2,4,8}_{0,0.03,0.06,0.1}; width slopes in W1 [12, 28] and W2 [28, 44]; search -> at most one candidate -> confirm",
    "prediction": "not fixed; the hope: a capacity at which the width holds in both windows",
    "criterion": "as F3.11 (reference free_0): bound |s| + 2 SE <= s_free0/4 in W1 and W2; search: candidate = a stream arm "
                 "with field bound in both windows, frame check OK; INSTRUMENT: frame check, rest mass of stream arms "
                 "1e-12, balances of free_0 and free_cl_0.06 1e-12, field_0 contracts against free_0 in W1 or W2; "
                 "confirm: PASSED bound at rest in W1 and W2, bound in motion in W1, |v_c| + 2 SE <= beta/4; FAILED "
                 "growing or collapsing at rest, or |v_c| - 2 SE > beta/4; otherwise UNDECIDED",
    "null_model": "free_0 (cold, no field: the diffusion of the discreteness); free_cl_0.06 (hot streams, no field)",
    "level": "mechanism",
    "meters": ["K.6"],
    "depends_on": ["F3.11", "F3.12a2", "F3.8", "F3.6c"],
    "params": {"D": 3, "T": 48.0, "half": 16.0, "r_src": 4.0, "t_src": 4.0, "slab": 2.0, "W1": [12.0, 28.0],
               "W2": [28.0, 44.0], "win_move": [5.0, 28.0], "beta_check": 0.6, "beta_move": 0.5, "tau_c": 1.0,
               "epsilon": 0.1, "phi_ref": 44.15, "caps": [2, 4, 8], "etas": [0.0, 0.03, 0.06, 0.1], "eta_cl": 0.06,
               "L_g": 2.0, "min_parents": 4, "n_seeds": 32, "n_frame": 8, "n_se": 2.0, "boot_B": 4000,
               "tol_bal": 1e-12, "arm": None},
    "parts": {"search": {"mode": "search"}, "confirm": {"mode": "confirm"}},
    "settable": {"confirm": ["arm"]},
    "quick": {"n_seeds": 2, "n_frame": 1, "boot_B": 200, "T": 24.0, "half": 8.0, "r_src": 2.0, "t_src": 2.0,
              "W1": [6.0, 14.0], "W2": [14.0, 22.0], "win_move": [3.0, 14.0], "slab": 1.0, "caps": [2],
              "etas": [0.0, 0.06]},
    "estimate": "search ~20-40 min with 10 processes (32 sprinklings of ~49 000 events, 17 arms of ~12 s each, the frame "
                "check on 8); confirm ~10 min",
}


# ============================================================================================== streams

def _gradient(X, v, par, u, values, L_g, nmin):
    """The F3.8 fit: the gradient of `values` (given at the parents par of v) in the rest space of u, as a spacelike
    4-vector (G0 E0 + G1 E1); None if fewer than nmin parents lie within the rest-frame lag L_g."""
    d = X[par] - X[v]
    lag = -mink(d, u[None, :])
    sel = lag <= L_g
    if int(sel.sum()) < nmin:
        return None
    dd = d[sel]
    pp = dd - mink(dd, u[None, :])[:, None] * u[None, :]
    E = rest_basis(u)
    q = -np.stack([mink(pp, e[None, :]) for e in E], axis=1)
    A = np.column_stack([np.ones(q.shape[0]), q])
    coef, *_ = np.linalg.lstsq(A, np.asarray(values)[sel], rcond=None)
    G = coef[1:]
    gn = float(np.sqrt((G * G).sum()))
    if gn == 0.0:
        return None
    return G, E, gn


def transport_streams(X, ptr, idx, cptr, cidx, s, U0, lam, cap, p, Lsrc):
    """
    One causal pass with streams that never merge. lam > 0: the F3.8 field and bending per stream; cap: the capacity
    reflection (None: none). Lsrc: a function taking 4-vectors (rows) to the source frame (for the readout). Returns per event the summed
    4-momentum, the summed source-frame energy and energy x rapidity, the multiset of choices, and the bookkeeping.
    """
    N, D = X.shape
    L_g, nmin = float(p["L_g"]), int(p["min_parents"])
    field = RetardedField(X, float(p["tau_c"])) if lam > 0.0 else None
    inc = [[] for _ in range(N)]
    load = np.zeros(N)
    Psum = np.zeros((N, D))
    Esum = np.zeros(N)
    ESsum = np.zeros(N)
    choices = []
    bends, n_refl, n_stream = [], 0, 0
    endP, endm, srcm = np.zeros(D), 0.0, 0.0
    for v in range(N):
        streams = list(inc[v])
        if s[v] > 0.0:
            streams.append(s[v] * U0[v])
            srcm += float(s[v])
        if not streams:
            continue
        P = np.array(streams)
        m = np.sqrt(np.maximum(mink(P, P), 1e-300))
        load[v] = float(m.sum())
        par = idx[ptr[v]:ptr[v + 1]]
        if cap is not None and len(P) > cap and par.size:
            Pt = P.sum(axis=0)
            U = Pt / math.sqrt(max(float(mink(Pt, Pt)), 1e-300))
            g = _gradient(X, v, par, U, load[par], L_g, nmin)
            if g is not None and g[2] > 1e-9:                    # a flat load (gradient at round-off) has no direction
                G, E, gn = g
                nhat = -(G[0] * E[0] + G[1] * E[1]) / gn          # unit spacelike, orthogonal to U, towards lower load
                a = -mink(P, nhat[None, :])                       # component of each stream along nhat in the frame U
                r = a < -1e-9 * mink(P, U[None, :])               # a numerical tie (|a| ~ round-off) is not reflected
                if r.any():
                    P[r] = P[r] - 2.0 * a[r, None] * nhat[None, :]
                    n_refl += int(r.sum())
        Psum[v] = P.sum(axis=0)
        Ps = Lsrc(P)
        Es = Ps[:, 0]
        Esum[v] = float(Es.sum())
        ESsum[v] = float((Es * np.arccosh(np.maximum(Es / m, 1.0))).sum())
        a_, b_ = cptr[v], cptr[v + 1]
        if b_ == a_:
            endP += P.sum(axis=0)
            endm += float(m.sum())
        else:
            ks = cidx[a_:b_]
            B = X[ks] - X[v]
            tau = np.sqrt(np.maximum(mink(B, B), 1e-300))
            lnPi = None
            if field is not None and par.size:
                lnPi = np.array([-lam * field.value(int(w)) for w in par])
            for pi, mi in zip(P, m):
                n_stream += 1
                u = pi / mi
                j = int(np.argmin(mink(u[None, :], B) / tau))
                k = int(ks[j])
                choices.append((v, k))
                send = u
                if lnPi is not None:                             # F3.8 bending of this stream, unchanged
                    g = _gradient(X, v, par, u, lnPi, L_g, nmin)
                    if g is not None:
                        G, E, gn = g
                        e_dir = -(G[0] * E[0] + G[1] * E[1]) / gn
                        eta = gn * float(tau[j])
                        send = boost_towards(u, e_dir, eta)
                        bends.append(eta)
                inc[k].append(mi * send)
        if field is not None:
            field.add(v, load[v])                                # the rest mass carried through v sources the field
    choices.sort()
    return {"Psum": Psum, "Esum": Esum, "ESsum": ESsum, "choices": choices, "endP": endP, "endm": endm, "srcm": srcm,
            "bend": float(np.mean(bends)) if bends else float("nan"),
            "refl": n_refl / max(n_stream, 1)}


# ============================================================================================== readout

def readout(Xs, E, ES, p, beta):
    """Per slab of source time: R_rms about the energy-weighted centroid, centroid x, internal spread, carriers."""
    T, sl = float(p["T"]), float(p["slab"])
    t_end = T * math.sqrt(1.0 - beta * beta) - sl
    edges = np.arange(float(p["t_src"]), t_end + 1e-9, sl)
    car = E > 0.0
    out = {k: [] for k in ("t", "R", "xc", "spread", "n")}
    for a, b in zip(edges[:-1], edges[1:]):
        m = car & (Xs[:, 0] >= a) & (Xs[:, 0] < b)
        w = E[m]
        out["t"].append(0.5 * (a + b))
        out["n"].append(int(m.sum()))
        if w.sum() > 0:
            c = (w[:, None] * Xs[m, 1:]).sum(axis=0) / w.sum()
            out["R"].append(float(np.sqrt((w * ((Xs[m, 1:] - c) ** 2).sum(axis=1)).sum() / w.sum())))
            out["xc"].append(float(c[0]))
            out["spread"].append(float(ES[m].sum() / w.sum()))
        else:
            for k in ("R", "xc", "spread"):
                out[k].append(float("nan"))
    return out


# ============================================================================================== task

def arm_list(p):
    """(key, eta, kind, field on, cap): kind 'merge' = the F3.11 transport, 'streams' = this test's transport."""
    e = float(p["eta_cl"])
    arms = [("free_0", 0.0, "merge", False, None), ("field_0", 0.0, "merge", True, None),
            (f"free_cl_{e:g}", e, "streams", False, None)]
    for eta in p["etas"]:
        if float(eta) > 0.0:
            arms.append((f"cl_{float(eta):g}", float(eta), "streams", True, None))
    for c in p["caps"]:
        for eta in p["etas"]:
            arms.append((f"cap{int(c)}_{float(eta):g}", float(eta), "streams", True, int(c)))
    return arms


def run_task(p, seed, src_seed, beta, arms, frame):
    D, T, H = int(p["D"]), float(p["T"]), float(p["half"])
    rng = np.random.default_rng(seed)
    X = sprinkle_poisson(Box(D, [0.0, -H, -H], [T, H + beta * T, H]), 1.0, rng)
    N = len(X)
    ptr, idx = links_csr(X, 1.0)
    cptr, cidx = children_csr(ptr, idx, N)
    Xs = boost(X, [beta, 0.0])
    src = np.flatnonzero((Xs[:, 0] >= 0.0) & (Xs[:, 0] < float(p["t_src"]))
                         & ((Xs[:, 1:] ** 2).sum(axis=1) <= float(p["r_src"]) ** 2))
    bc = float(p["beta_check"])
    Xb = boost(X, [bc, 0.0])
    lam = float(p["epsilon"]) / float(p["phi_ref"])
    Ls = lambda P: boost(P, [beta, 0.0])                         # lab -> source frame
    Lsb = lambda P: boost(boost(P, [-bc, 0.0]), [beta, 0.0])      # checked frame -> lab -> source frame
    out = {"N": N, "sources": int(src.size), "arms": {}}
    for key, eta, kind, on, cap in arms:
        _, U_rest = source_config("cold" if eta == 0.0 else "sym", float(eta), src.size,
                                  np.random.default_rng(derive_seed(src_seed, "eta", float(eta))))
        s = np.zeros(N)
        U0 = np.zeros((N, D))
        s[src] = 1.0
        U0[src] = boost(U_rest, [-beta, 0.0])
        lam_a = lam if on else 0.0
        rec = {}
        if kind == "merge":
            r = transport_merge(X, ptr, idx, cptr, cidx, s, U0, lam_a, p)
            Ps = boost(r["P"], [beta, 0.0])
            E = np.where(Ps[:, 0] > 0.0, Ps[:, 0], 0.0)
            M = np.sqrt(np.maximum(mink(Ps, Ps), 1e-300))
            ES = np.where(E > 0.0, E * np.arccosh(np.maximum(E / M, 1.0)), 0.0)
            if frame:
                rb = transport_merge(Xb, ptr, idx, cptr, cidx, s, boost(U0, [bc, 0.0]), lam_a, p)
                rec["frame_identical"] = bool(np.array_equal(r["choice"], rb["choice"]))
            rec["bend"] = float(np.nanmean(r["bend"])) if np.isfinite(r["bend"]).any() else float("nan")
            if not on:
                held = r["P"][:, 0] > 0.0
                ends = held & (r["choice"] < 0)
                Psrc = U0[src].sum(axis=0)
                rec["balance"] = float(np.max(np.abs(Psrc - r["P"][ends].sum(axis=0))) / Psrc[0])
        else:
            r = transport_streams(X, ptr, idx, cptr, cidx, s, U0, lam_a, cap, p, Ls)
            E, ES = r["Esum"], r["ESsum"]
            if frame:
                rb = transport_streams(Xb, ptr, idx, cptr, cidx, s, boost(U0, [bc, 0.0]), lam_a, cap, p, Lsb)
                rec["frame_identical"] = r["choices"] == rb["choices"]
            rec["bend"], rec["refl"] = r["bend"], r["refl"]
            rec["mass"] = abs(r["endm"] - r["srcm"]) / max(r["srcm"], 1e-300)
            if not on and cap is None:
                Psrc = U0[src].sum(axis=0)
                rec["balance"] = float(np.max(np.abs(Psrc - r["endP"])) / Psrc[0])
        ser = readout(Xs, E, ES, p, beta)
        rec.update({"series": ser, "s1": slope(ser["t"], ser["R"], p["W1"]), "s2": slope(ser["t"], ser["R"], p["W2"]),
                    "vc": slope(ser["t"], ser["xc"], p["win_move"])})
        out["arms"][key] = rec
    return out


# ============================================================================================== statistics

def summarise(rows, keys, p, ctx, tag):
    B, nse = int(p["boot_B"]), float(p["n_se"])
    st = {}
    for key in keys:
        bs = lambda *t, key=key: derive_seed(ctx.seed, "F3.12b-boot", tag, key, *t)
        st[key] = {w: boot_se([r_["arms"][key][w] for r_ in rows], B, bs(w))[:2] for w in ("s1", "s2", "vc")}
        fr = [r_["arms"][key]["frame_identical"] for r_ in rows if "frame_identical" in r_["arms"][key]]
        st[key]["frame"] = all(fr) if fr else True
        st[key]["n_frame"] = len(fr)
        st[key]["balance"] = max((r_["arms"][key].get("balance", 0.0) for r_ in rows), default=0.0)
        st[key]["mass"] = max((r_["arms"][key].get("mass", 0.0) for r_ in rows), default=0.0)
        with np.errstate(all="ignore"):
            for k in ("R", "spread", "n"):
                st[key][k] = np.nanmean(np.array([r_["arms"][key]["series"][k] for r_ in rows], float), axis=0)
        st[key]["bend"] = float(np.nanmean([r_["arms"][key]["bend"] for r_ in rows]))
        st[key]["refl"] = float(np.nanmean([r_["arms"][key].get("refl", np.nan) for r_ in rows]))
    ref1, ref2 = st["free_0"]["s1"][0], st["free_0"]["s2"][0]
    t = np.asarray(rows[0]["arms"]["free_0"]["series"]["t"], float)
    wall = float(p["half"]) / 2.0
    for key in keys:
        st[key]["c1"] = classify(*st[key]["s1"], ref1, nse)
        st[key]["c2"] = classify(*st[key]["s2"], ref2, nse)
        for c, w in (("c1", p["W1"]), ("c2", p["W2"])):
            m = (t >= w[0]) & (t <= w[1])
            if m.any() and np.nanmax(st[key]["R"][m]) > wall and st[key][c] in ("bound", "unresolved", "collapsing"):
                st[key][c] = "wall-limited"
    return st


def lines_for(st, keys, t):
    out = ["        t       " + " ".join(f"{v:5.1f}" for v in t)]
    for key in keys:
        x = st[key]
        out.append(f"  {key:12s}: slope W1 {x['s1'][0]:+.4f} +- {x['s1'][1]:.4f} -> {x['c1']}; W2 {x['s2'][0]:+.4f} +- "
                   f"{x['s2'][1]:.4f} -> {x['c2']}; frame {'OK' if x['frame'] else 'FAILS'} ({x['n_frame']})"
                   + (f"; bending {x['bend']:.4f}" if np.isfinite(x["bend"]) else "")
                   + (f"; reflections per stream {x['refl']:.3f}" if np.isfinite(x["refl"]) else "")
                   + (f"; rest mass {x['mass']:.1e}" if x["mass"] else "")
                   + (f"; balance {x['balance']:.1e}" if key.startswith("free") else ""))
        out.append("        R_rms   " + " ".join(f"{v:5.2f}" if np.isfinite(v) else "  nan" for v in x["R"]))
        out.append("        spread  " + " ".join(f"{v:5.3f}" if np.isfinite(v) else "  nan" for v in x["spread"]))
        out.append("        carry   " + " ".join(f"{v:5.1f}" for v in x["n"]))
    return out


def _jsonable(st):
    return {k: {kk: (list(vv) if isinstance(vv, (tuple, np.ndarray)) else vv) for kk, vv in v.items()} for k, v in st.items()}


def _instrument(st, keys, nse, tol):
    stream_keys = [k for k in keys if k not in ("free_0", "field_0")]
    free_keys = [k for k in keys if k.startswith("free")]
    return {"frame check": all(st[k]["frame"] for k in keys),
            f"rest mass of stream arms <= {tol:g}": all(st[k]["mass"] <= tol for k in stream_keys),
            f"balances of free arms <= {tol:g}": all(st[k]["balance"] <= tol for k in free_keys)}


# ============================================================================================== run

def run(p, ctx):
    mode = str(p["mode"])
    t0 = time.time()
    n, nse, tol = int(p["n_seeds"]), float(p["n_se"]), float(p["tol_bal"])
    nfr = int(p["n_frame"])
    arms = arm_list(p)
    head = [f"2+1 sprinkling rho = 1, box t [0, {p['T']:g}], half width {p['half']:g}; source radius {p['r_src']:g}, time "
            f"[0, {p['t_src']:g}); field and bending of F3.8 (lambda = {p['epsilon']:g} / {p['phi_ref']:g}); streams never "
            f"merge; capacity reflection at events with more than C streams; windows W1 {p['W1']}, W2 {p['W2']}; "
            f"reference free_0; frame check beta' = {p['beta_check']} on the first {nfr} sprinklings"]
    if mode == "search":
        keys = [a[0] for a in arms]
        jobs = [(p, derive_seed(ctx.seed, "F3.12b-search", k), derive_seed(ctx.seed, "F3.12b-search-src", k), 0.0, arms,
                 k < nfr) for k in range(n)]
        rows = parallel_map(run_task, jobs, ctx.workers, ctx.progress,
                            lambda i, r_, d, nn: f"[F3.12b] search sprinkling {d}/{nn} ({time.time()-t0:.0f}s)")
        st = summarise(rows, keys, p, ctx, "search")
        lines = head + [f"search: {n} sprinklings, source at rest (~{np.mean([r_['N'] for r_ in rows]):.0f} events, "
                        f"{np.mean([r_['sources'] for r_ in rows]):.0f} source events)"]
        lines += lines_for(st, keys, rows[0]["arms"]["free_0"]["series"]["t"])
        inst = _instrument(st, keys, nse, tol)
        inst["positive control: field_0 contracts against free_0 in W1 or W2"] = bool(
            any(st["field_0"][w][0] + nse * st["field_0"][w][1] < st["free_0"][w][0] for w in ("s1", "s2")))
        lines.append("instrument: " + "; ".join(f"{k}: {'OK' if v else 'FAILS'}" for k, v in inst.items()))
        if not all(inst.values()):
            return Result("INSTRUMENT_FAILED", lines, {"stats": _jsonable(st), "instrument": inst})
        cands = []
        for key, eta, kind, on, cap in arms:
            if kind != "streams" or not on:
                continue
            x = st[key]
            if x["c1"] == "bound" and x["c2"] == "bound" and x["frame"]:
                cands.append((abs(x["s1"][0]) + nse * x["s1"][1], cap if cap is not None else 0, eta, key))
        cand = sorted(cands)[0][3] if cands else None
        lines.append(f"candidate (rule stated in the header): {cand if cand else 'none'}"
                     + (f"  -> confirm with: --set arm={cand}" if cand else ""))
        return Result("RECORDED", lines, {"candidate": cand, "stats": _jsonable(st), "instrument": inst})

    key = p.get("arm")
    allowed = {a[0]: a for a in arms if a[2] == "streams" and a[3]}
    if key not in allowed:
        raise SystemExit(f"F3.12b confirm: --set arm=<one of {sorted(allowed)}>")
    sel = [a for a in arms if a[0] in ("free_0", f"free_cl_{float(p['eta_cl']):g}", key)]
    keys = [a[0] for a in sel]
    bm = float(p["beta_move"])
    jr = [(p, derive_seed(ctx.seed, "F3.12b-confirm-rest", k), derive_seed(ctx.seed, "F3.12b-confirm-rest-src", k), 0.0,
           sel, k < nfr) for k in range(n)]
    jm = [(p, derive_seed(ctx.seed, "F3.12b-confirm-move", k), derive_seed(ctx.seed, "F3.12b-confirm-move-src", k), bm,
           sel, k < nfr) for k in range(n)]
    rows_r = parallel_map(run_task, jr, ctx.workers, ctx.progress,
                          lambda i, r_, d, nn: f"[F3.12b] confirm rest {d}/{nn} ({time.time()-t0:.0f}s)")
    rows_m = parallel_map(run_task, jm, ctx.workers, ctx.progress,
                          lambda i, r_, d, nn: f"[F3.12b] confirm move {d}/{nn} ({time.time()-t0:.0f}s)")
    sr = summarise(rows_r, keys, p, ctx, "confirm-rest")
    sm = summarise(rows_m, keys, p, ctx, "confirm-move")
    lines = head + [f"confirm: candidate {key} (CHOSEN FROM DATA); {n} sprinklings at rest, {n} moving (beta = {bm})",
                    "  at rest:"]
    lines += lines_for(sr, keys, rows_r[0]["arms"]["free_0"]["series"]["t"])
    lines.append(f"  moving source (beta = {bm}):")
    lines += lines_for(sm, keys, rows_m[0]["arms"]["free_0"]["series"]["t"])
    v, se = sm[key]["vc"]
    tolv = bm / 4.0
    lines.append(f"  inertia: centroid velocity of {key} in the source frame {v:+.4f} +- {se:.4f} vs beta/4 = {tolv:.4f}")
    inst = {**_instrument(sr, keys, nse, tol), **{k_ + " (moving)": v_ for k_, v_ in _instrument(sm, keys, nse, tol).items()}}
    lines.append("instrument: " + "; ".join(f"{k_}: {'OK' if v_ else 'FAILS'}" for k_, v_ in inst.items()))
    c = sr[key]
    if not all(inst.values()):
        status = "INSTRUMENT_FAILED"
    elif c["c1"] in ("growing", "collapsing") or c["c2"] in ("growing", "collapsing") or \
            (np.isfinite(v) and abs(v) - nse * se > tolv):
        status = "FAILED"
    elif c["c1"] == "bound" and c["c2"] == "bound" and sm[key]["c1"] == "bound" and np.isfinite(v) and abs(v) + nse * se <= tolv:
        status = "PASSED"
    else:
        status = "UNDECIDED"
    return Result(status, lines, {"candidate": key, "rest": _jsonable(sr), "move": _jsonable(sm)})
