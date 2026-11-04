"""
F3.11 -- A hot bundle in its own field: does internal motion balance the pull of the retarded field?
====================================================================================================

Why this test exists (25 Oct 2026, the author's question after F3.6c)
--------------------------------------------------------------------
Two measured halves of one equilibrium, never run together:
  F3.8  the retarded field of the bundle (the resistance of the substrate,
        slower clocks with range) pulls a COLD bundle in: epsilon = 0.1 stops
        and slowly collapses, epsilon = 0.5 collapses and bounces. The reading
        of F3.8 (written then): an equilibrium needs something that resists
        compression -- pressure, a velocity dispersion, or saturation.
  F3.6c a velocity dispersion with exact 4-momentum bookkeeping (merging sends
        on the whole P, flux F = |P|): alone it drives the bundle apart at about
        half the free speed; merging only cools.
The question is whether the two balance (virial equilibrium, as a star or a
cluster). No new mechanism: the field and the bending are F3.8, the hot
source and the merge rule are F3.6c. The strength epsilon = 0.1 is LOCKED
from F3.8 (the only strength that stopped the cold bundle without a bounce;
replicated in F3.9), with lambda = 0.1 / 44.15 (phi_ref of the F3.8 search,
taken as a number, not recalibrated). The scan is over the temperature eta.

Model (locked before any run)
-----------------------------
2+1 Poisson sprinkling, rho = 1 (lengths in l), links (core.order.links).
Source: world tube of radius 4, time [0, 4) in the source frame; source event i
emits weight 1 along its own unit 4-velocity u_i: eta = 0 -> u_i = u_src (cold);
eta > 0 -> the F3.6c 'sym' config (rapidity eta, angles phi_0 + 2 pi k / n over
a random permutation, spatial momenta cancel exactly in the source frame).
Transport (one causal pass, F3.6c rule B): an event holding streams has
P_v = sum of the incoming 4-momenta + its emission; flux F_v = |P_v| (invariant
mass), direction u_v = P_v / F_v; it sends everything to its child of smallest
rapidity relative to u_v. Field arms: the direction sent is bent exactly as
in F3.8 (transverse gradient of ln Pi = -lambda phi over the parents with
rest-frame lag <= 2, at least 4; boost by |G| tau_k towards -G); the child
receives F_v times the bent direction (bending keeps the rest mass; the
momentum change is the exchange with the field). The field is F3.8's
phi(x) = sum_{w < x, tau >= 1} F_w / (2 pi tau), with F_w = |P_w| (a Lorentz
scalar; for a cold bundle equal to the F3.8 flux). Declared limitation: the
kinetic energy of the hot streams sources the field only after merging has
turned it into rest mass (a scalar source; a tensor source would count it).
Arms (search): free_0 (cold, no field: the F3.6 transport; the reference of the
classification), field_0 (cold in its own field: the positive control), and for
eta in {0.15, 0.3, 0.45}: free_eta (no field) and field_eta.

Readout (per sprinkling, per arm; in the source frame)
------------------------------------------------------
Slabs of source time (width 2) from 4 up to where the source world line leaves
the box. Entries: events holding streams, weight = energy in the source frame
P . u_src. R_rms about the weighted centroid; centroid x; internal spread =
weighted mean rapidity of u_v relative to u_src; number of carriers. Slopes of
R_rms in W1 [12, 28] and W2 [28, 44]. Mean bending per step. Mergers. Balance
of the free arms: |sum_sources u_i - sum_ends P| / sum u_i^0 (implementation
check; field arms exchange momentum with the field and have no such balance).

Classification (exactly F3.8, reference free_0)
-----------------------------------------------
Per window, tol = s(free_0) / 4 in that window: bound |s| + 2 SE <= tol;
growing s - 2 SE > tol; collapsing s + 2 SE < -tol; otherwise unresolved; a
width above half / 2 in the window turns bound, unresolved or collapsing into
wall-limited.

Parts and criteria (stated before any run)
------------------------------------------
search  RECORDED, at most one candidate: a hot field arm bound in W1 AND W2 with
        identical choices in the frame boosted by beta' = 0.6; the smallest
        |s_W1| + 2 SE, ties the smaller eta. INSTRUMENT_FAILED if any arm fails
        the frame check, a free arm's balance exceeds 1e-12, or the positive
        control fails: field_0 must contract relative to free_0 in W1 or in W2
        (s(field_0) + 2 SE < s(free_0) in at least one window), i.e. the field of
        this run pulls.
confirm --set eta=<candidate> (one of the grid), new streams, the source at rest
        and moving with beta = 0.5 (directions boosted with the source):
        PASSED  field_eta bound at rest in W1 and W2, bound in motion in W1,
                inertia |v_c| + 2 SE <= beta/4 (centroid velocity in the source
                frame over [5, 28]), frame check and balances OK;
        FAILED  field_eta growing or collapsing at rest in W1 or W2, or
                |v_c| - 2 SE > beta/4;
        INSTRUMENT_FAILED frame check or balance;
        UNDECIDED otherwise.
Reading declared now. A bound hot bundle would be the first structure of the
laboratory held by its own field against its own internal motion. Expected
obstacle (F3.9 reading): merging is dissipative (it cools), so a plateau may
be followed by a slow contraction; W2 is there to see it. In 2+1 the potential
is logarithmic and confines; a PASSED asks for the 3+1 version before any
physical reading. The bending rule of the F3 series is the conformal geometry
(PPN gamma = -1, findings: Derived consequences): this tests structure, not light.

Development runs (after everything above was written; decide nothing). Quick run (3 small sprinklings): mechanics in
order (frame check identical, free balances 6e-16); the field is weak there (lambda is locked to the full-size
phi_ref, the quick source has 24 events). One full-size sprinkling (~49 000 events, 188 source events, ~75 s per
arm on one core, links dominate): free_0 R 2.9 -> 5.1 (s1 +0.037, s2 +0.022); field_0 contracts 3.4 -> 2.4 by
t ~ 29 and then rebounds to ~3.3 (s1 -0.050, s2 +0.052); field_0.3 R ~4.3 with s2 -0.003; mean bending 0.008 per
step. Changed after it, before any full run: the positive control originally asked for contraction of field_0 in
W2 only; a collapse in W1 followed by a rebound in W2 (seen here) is also a pulling field, so the control asks
for contraction in W1 OR W2. No criterion of the decision changed.
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

TEST = {
    "id": "F3.11",
    "phase": "F3",
    "title": "A hot bundle in its own field: does internal motion balance the pull of the retarded field?",
    "theory_part": "binding of matter: the resistance of a substrate with range (field F3.8) against internal motion "
                   "with conserved 4-momentum (F3.6c); a virial-type balance",
    "question": "Is there a source temperature η at which a bundle in its own retarded field (ε = 0.1 locked from "
                "F3.8) keeps a constant width, between dispersal without a field and the collapse of a cold bundle? "
                "Does such a bundle move with the source?",
    "method": "2+1 sprinkling (ρ = 1), source of radius 4 with streams of rapidity η (configuration sym of F3.6c), "
              "merging B (flow = |P|), bending and field exactly as F3.8 (λ = 0.1/44.15, field excited by |P|); "
              "branches free_0, field_0 (control) and free/field for η = 0.15, 0.3, 0.45; slope of the width in the "
              "windows [12, 28] and [28, 44]; protocol search → at most one candidate → confirm on new sprinklings "
              "including a moving source",
    "prediction": "not fixed; against binding stands the cooling by merging (F3.9: slow shrinking after a plateau)",
    "criterion": "as F3.8 (reference free_0): bound |s| + 2 SE <= s_free0/4 in W1 and W2; search: candidate = a hot "
                 "branch with a field bound in both windows with the same choice in the frame β′ = 0.6; INSTRUMENT: "
                 "frame check, balance of the free branches > 1e-12, or field_0 does not contract against free_0 in "
                 "either W1 or W2; confirm: PASSED bound at rest in W1 and W2, bound in motion in W1 and |v_c| + 2 "
                 "SE <= β/4; FAILED growing or collapse at rest, or |v_c| − 2 SE > β/4; otherwise UNDECIDED. Success "
                 "in 2+1 requires a 3+1 version before a physical reading.",
    "null_model": "free_η: the same hot bundle without a field; free_0: a cold bundle without a field (diffusion of "
                  "the discreteness)",
    "level": "mechanism",
    "meters": ["K.6"],
    "depends_on": ["F3.6c", "F3.8"],
    "params": {"D": 3, "T": 48.0, "half": 16.0, "r_src": 4.0, "t_src": 4.0, "slab": 2.0, "W1": [12.0, 28.0],
               "W2": [28.0, 44.0], "win_move": [5.0, 28.0], "beta_check": 0.6, "beta_move": 0.5, "tau_c": 1.0,
               "epsilon": 0.1, "phi_ref": 44.15, "etas": [0.15, 0.3, 0.45], "L_g": 2.0, "min_parents": 4,
               "n_seeds": 32, "n_se": 2.0, "boot_B": 4000, "tol_bal": 1e-12, "eta": None},
    "parts": {"search": {"mode": "search"}, "confirm": {"mode": "confirm"}},
    "settable": {"confirm": ["eta"]},
    "quick": {"n_seeds": 3, "boot_B": 200, "T": 24.0, "half": 8.0, "r_src": 2.0, "t_src": 2.0,
              "W1": [6.0, 14.0], "W2": [14.0, 22.0], "win_move": [3.0, 14.0], "slab": 1.0},
    "estimate": "search ~15–30 min with 10 processes (sprinkling ~49 000 events, links ~1 min, 8 branches also in "
                "another frame); confirm similar",
}


# ============================================================================================== transport

def transport(X, ptr, idx, cptr, cidx, s, U0, lam, p):
    """
    One causal pass of rule B with (lam > 0) or without the F3.8 bending. Returns per event: the 4-momentum held
    (P), the chosen child, the bending rapidity, the merger measure eta_eq, and the 4-momentum that leaves the box.
    """
    N, D = X.shape
    L_g, nmin = float(p["L_g"]), int(p["min_parents"])
    field = RetardedField(X, float(p["tau_c"])) if lam > 0.0 else None
    Pin = np.zeros((N, D))
    Fin = np.zeros(N)
    nin = np.zeros(N, dtype=np.int64)
    Phold = np.zeros((N, D))
    choice = np.full(N, -1, dtype=np.int64)
    bend = np.full(N, np.nan)
    eta_eq = np.full(N, np.nan)
    for v in range(N):
        if nin[v] == 0 and s[v] <= 0.0:
            continue
        P = Pin[v] + s[v] * U0[v]
        M = math.sqrt(max(float(mink(P, P)), 1e-300))
        u = P / M
        Phold[v] = P
        if nin[v] + (1 if s[v] > 0.0 else 0) >= 2:
            eta_eq[v] = 2.0 * math.acosh(max(M / (Fin[v] + s[v]), 1.0))
        a, b = cptr[v], cptr[v + 1]
        if b > a:
            ks = cidx[a:b]
            B = X[ks] - X[v]
            tau = np.sqrt(np.maximum(mink(B, B), 1e-300))
            j = int(np.argmin(mink(u[None, :], B) / tau))
            k = int(ks[j])
            choice[v] = k
            send = u
            if field is not None:                       # F3.8 bending, unchanged
                pa, pb = ptr[v], ptr[v + 1]
                par = idx[pa:pb]
                if par.size:
                    d = X[par] - X[v]
                    lag = -mink(d, u[None, :])
                    sel = lag <= L_g
                    if int(sel.sum()) >= nmin:
                        pars = par[sel]
                        dd = d[sel]
                        pp = dd - mink(dd, u[None, :])[:, None] * u[None, :]
                        E = rest_basis(u)
                        q = -np.stack([mink(pp, e[None, :]) for e in E], axis=1)
                        lnPi = np.array([-lam * field.value(int(w)) for w in pars])
                        A = np.column_stack([np.ones(q.shape[0]), q])
                        coef, *_ = np.linalg.lstsq(A, lnPi, rcond=None)
                        G = coef[1:]
                        gn = float(np.sqrt((G * G).sum()))
                        if gn > 0.0:
                            e_dir = -(G[0] * E[0] + G[1] * E[1]) / gn
                            eta = gn * float(tau[j])
                            send = boost_towards(u, e_dir, eta)
                            bend[v] = eta
            Pin[k] += M * send                          # the whole rest mass, along the (bent) direction
            Fin[k] += M
            nin[k] += 1
        if field is not None:
            field.add(v, M)                             # v is final: its invariant mass sources the field of its future
    return {"P": Phold, "choice": choice, "bend": bend, "eta_eq": eta_eq}


# ============================================================================================== readout

def readout(Xs, Ps, p, beta):
    """Per slab of source time: R_rms about the weighted centroid, centroid x, internal spread, carriers."""
    T, sl = float(p["T"]), float(p["slab"])
    t_end = T * math.sqrt(1.0 - beta * beta) - sl
    edges = np.arange(float(p["t_src"]), t_end + 1e-9, sl)
    E = Ps[:, 0]                                                   # energy in the source frame (u_src = (1, 0, 0) there)
    car = E > 0.0
    M = np.sqrt(np.maximum(mink(Ps, Ps), 1e-300))
    eta_rel = np.where(car, np.arccosh(np.maximum(E / M, 1.0)), 0.0)
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
            out["spread"].append(float((w * eta_rel[m]).sum() / w.sum()))
        else:
            for k in ("R", "xc", "spread"):
                out[k].append(float("nan"))
    return out


# ============================================================================================== task

def run_task(p, seed, src_seed, beta, arms):
    """One sprinkling: the arms (key, eta, field on/off), readout in the source frame, frame check, balances."""
    D, T, H = int(p["D"]), float(p["T"]), float(p["half"])
    rng = np.random.default_rng(seed)
    X = sprinkle_poisson(Box(D, [0.0, -H, -H], [T, H + beta * T, H]), 1.0, rng)
    N = len(X)
    ptr, idx = links_csr(X, 1.0)
    cptr, cidx = children_csr(ptr, idx, N)
    Xs = boost(X, [beta, 0.0])                                     # source frame
    src = np.flatnonzero((Xs[:, 0] >= 0.0) & (Xs[:, 0] < float(p["t_src"]))
                         & ((Xs[:, 1:] ** 2).sum(axis=1) <= float(p["r_src"]) ** 2))
    bc = float(p["beta_check"])
    Xb = boost(X, [bc, 0.0])
    lam = float(p["epsilon"]) / float(p["phi_ref"])
    out = {"N": N, "sources": int(src.size), "arms": {}}
    for key, eta, on in arms:
        _, U_rest = source_config("cold" if eta == 0.0 else "sym", float(eta), src.size,
                                  np.random.default_rng(derive_seed(src_seed, "eta", float(eta))))
        s = np.zeros(N)
        U0 = np.zeros((N, D))
        s[src] = 1.0
        U0[src] = boost(U_rest, [-beta, 0.0])                      # rest-frame directions carried with the source
        lam_a = lam if on else 0.0
        r = transport(X, ptr, idx, cptr, cidx, s, U0, lam_a, p)
        rb = transport(Xb, ptr, idx, cptr, cidx, s, boost(U0, [bc, 0.0]), lam_a, p)
        Ps = boost(r["P"], [beta, 0.0])
        ser = readout(Xs, Ps, p, beta)
        rec = {"series": ser, "s1": slope(ser["t"], ser["R"], p["W1"]), "s2": slope(ser["t"], ser["R"], p["W2"]),
               "vc": slope(ser["t"], ser["xc"], p["win_move"]),
               "frame_identical": bool(np.array_equal(r["choice"], rb["choice"])),
               "bend_mean": float(np.nanmean(r["bend"])) if np.isfinite(r["bend"]).any() else float("nan"),
               "n_merge": int(np.isfinite(r["eta_eq"]).sum()),
               "eta_merge": float(np.nanmean(r["eta_eq"])) if np.isfinite(r["eta_eq"]).any() else 0.0}
        if not on:
            held = r["P"][:, 0] > 0.0
            ends = held & (r["choice"] < 0)
            Psrc = U0[src].sum(axis=0)
            rec["balance"] = float(np.max(np.abs(Psrc - r["P"][ends].sum(axis=0))) / Psrc[0])
        out["arms"][key] = rec
    return out


# ============================================================================================== statistics

def summarise(rows, keys, p, ctx, tag):
    B, nse = int(p["boot_B"]), float(p["n_se"])
    st = {}
    for key in keys:
        bs = lambda *t, key=key: derive_seed(ctx.seed, "F3.11-boot", tag, key, *t)
        st[key] = {w: boot_se([r_["arms"][key][w] for r_ in rows], B, bs(w))[:2] for w in ("s1", "s2", "vc")}
        st[key]["frame"] = all(r_["arms"][key]["frame_identical"] for r_ in rows)
        st[key]["balance"] = max((r_["arms"][key].get("balance", 0.0) for r_ in rows), default=0.0)
        with np.errstate(all="ignore"):
            for k in ("R", "spread", "n"):
                st[key][k] = np.nanmean(np.array([r_["arms"][key]["series"][k] for r_ in rows], float), axis=0)
        st[key]["bend"] = float(np.nanmean([r_["arms"][key]["bend_mean"] for r_ in rows]))
        st[key]["n_merge"] = float(np.mean([r_["arms"][key]["n_merge"] for r_ in rows]))
        st[key]["eta_merge"] = float(np.mean([r_["arms"][key]["eta_merge"] for r_ in rows]))
    ref1, ref2 = st["free_0"]["s1"][0], st["free_0"]["s2"][0]
    t = np.asarray(rows[0]["arms"]["free_0"]["series"]["t"], float)
    wall = float(p["half"]) / 2.0
    for key in keys:
        st[key]["c1"] = classify(*st[key]["s1"], ref1, nse)
        st[key]["c2"] = classify(*st[key]["s2"], ref2, nse)
        for c, w in (("c1", p["W1"]), ("c2", p["W2"])):     # a width that reached the box cannot be read as bound
            m = (t >= w[0]) & (t <= w[1])
            if m.any() and np.nanmax(st[key]["R"][m]) > wall and st[key][c] in ("bound", "unresolved", "collapsing"):
                st[key][c] = "wall-limited"
    return st


def lines_for(st, keys, t):
    out = ["        t       " + " ".join(f"{v:5.1f}" for v in t)]
    for key in keys:
        x = st[key]
        out.append(f"  {key:11s}: slope W1 {x['s1'][0]:+.4f} +- {x['s1'][1]:.4f} -> {x['c1']}; W2 {x['s2'][0]:+.4f} +- "
                   f"{x['s2'][1]:.4f} -> {x['c2']}; frame {'OK' if x['frame'] else 'FAILS'}"
                   + (f"; mean bending {x['bend']:.4f}" if np.isfinite(x["bend"]) else "")
                   + f"; mergers {x['n_merge']:.0f} (eta_eq {x['eta_merge']:.3f})"
                   + (f"; balance {x['balance']:.1e}" if key.startswith("free") else ""))
        out.append("        R_rms   " + " ".join(f"{v:5.2f}" if np.isfinite(v) else "  nan" for v in x["R"]))
        out.append("        spread  " + " ".join(f"{v:5.3f}" if np.isfinite(v) else "  nan" for v in x["spread"]))
        out.append("        carry   " + " ".join(f"{v:5.1f}" for v in x["n"]))
    return out


def _jsonable(st):
    return {k: {kk: (list(vv) if isinstance(vv, (tuple, np.ndarray)) else vv) for kk, vv in v.items()} for k, v in st.items()}


# ============================================================================================== run

def run(p, ctx):
    mode = str(p["mode"])
    t0 = time.time()
    n, nse, tol_bal = int(p["n_seeds"]), float(p["n_se"]), float(p["tol_bal"])
    head = [f"2+1 sprinkling rho = 1, box t [0, {p['T']:g}], half width {p['half']:g}; source radius {p['r_src']:g}, time "
            f"[0, {p['t_src']:g}), proper streams of rapidity eta (sym); merge rule B (F = |P|); field phi = sum |P| / (2 pi tau), "
            f"tau >= {p['tau_c']:g}; ln Pi = -lambda phi, lambda = {p['epsilon']:g} / {p['phi_ref']:g} (locked from F3.8); windows "
            f"W1 {p['W1']}, W2 {p['W2']}; reference free_0; frame check beta' = {p['beta_check']}"]
    if mode == "search":
        etas = [float(e) for e in p["etas"]]
        arms = [("free_0", 0.0, False), ("field_0", 0.0, True)]
        for e in etas:
            arms += [(f"free_{e:g}", e, False), (f"field_{e:g}", e, True)]
        keys = [a[0] for a in arms]
        jobs = [(p, derive_seed(ctx.seed, "F3.11-search", k), derive_seed(ctx.seed, "F3.11-search-src", k), 0.0, arms)
                for k in range(n)]
        rows = parallel_map(run_task, jobs, ctx.workers, ctx.progress,
                            lambda i, r_, d, nn: f"[F3.11] search sprinkling {d}/{nn} ({time.time()-t0:.0f}s)")
        st = summarise(rows, keys, p, ctx, "search")
        lines = head + [f"search: {n} sprinklings, source at rest (~{np.mean([r_['N'] for r_ in rows]):.0f} events, "
                        f"{np.mean([r_['sources'] for r_ in rows]):.0f} source events)"]
        lines += lines_for(st, keys, rows[0]["arms"]["free_0"]["series"]["t"])
        inst = {"frame check (all arms)": all(st[k]["frame"] for k in keys),
                f"free balances <= {tol_bal:g}": all(st[k]["balance"] <= tol_bal for k in keys if k.startswith("free")),
                "positive control: field_0 contracts against free_0 in W1 or W2": bool(
                    any(st["field_0"][w][0] + nse * st["field_0"][w][1] < st["free_0"][w][0] for w in ("s1", "s2")))}
        lines.append("instrument: " + "; ".join(f"{k}: {'OK' if v else 'FAILS'}" for k, v in inst.items()))
        if not all(inst.values()):
            return Result("INSTRUMENT_FAILED", lines, {"stats": _jsonable(st), "instrument": inst})
        cands = []
        for e in etas:
            x = st[f"field_{e:g}"]
            if x["c1"] == "bound" and x["c2"] == "bound" and x["frame"]:
                cands.append((abs(x["s1"][0]) + nse * x["s1"][1], e))
        cand = sorted(cands)[0][1] if cands else None
        lines.append(f"candidate (rule stated in the header): {('eta = %g' % cand) if cand is not None else 'none'}"
                     + (f"  -> confirm with: --set eta={cand:g}" if cand is not None else ""))
        return Result("RECORDED", lines, {"candidate": cand, "stats": _jsonable(st), "instrument": inst})

    eta = p.get("eta")
    if eta is None:
        raise SystemExit("F3.11 confirm: give the candidate of the search, e.g. --set eta=0.3")
    eta = float(eta)
    if eta not in [float(e) for e in p["etas"]]:
        raise SystemExit(f"F3.11 confirm: eta must be one of the grid {p['etas']} (got {eta})")
    key = f"field_{eta:g}"
    arms = [("free_0", 0.0, False), (f"free_{eta:g}", eta, False), (key, eta, True)]
    keys = [a[0] for a in arms]
    bm = float(p["beta_move"])
    jr = [(p, derive_seed(ctx.seed, "F3.11-confirm-rest", k), derive_seed(ctx.seed, "F3.11-confirm-rest-src", k), 0.0, arms)
          for k in range(n)]
    jm = [(p, derive_seed(ctx.seed, "F3.11-confirm-move", k), derive_seed(ctx.seed, "F3.11-confirm-move-src", k), bm, arms)
          for k in range(n)]
    rows_r = parallel_map(run_task, jr, ctx.workers, ctx.progress,
                          lambda i, r_, d, nn: f"[F3.11] confirm rest {d}/{nn} ({time.time()-t0:.0f}s)")
    rows_m = parallel_map(run_task, jm, ctx.workers, ctx.progress,
                          lambda i, r_, d, nn: f"[F3.11] confirm move {d}/{nn} ({time.time()-t0:.0f}s)")
    sr = summarise(rows_r, keys, p, ctx, "confirm-rest")
    sm = summarise(rows_m, keys, p, ctx, "confirm-move")
    lines = head + [f"confirm: candidate eta = {eta:g} (chosen from the search, CHOSEN FROM DATA); {n} sprinklings at rest, "
                    f"{n} moving (beta = {bm})", "  at rest:"]
    lines += lines_for(sr, keys, rows_r[0]["arms"]["free_0"]["series"]["t"])
    lines.append(f"  moving source (beta = {bm}):")
    lines += lines_for(sm, keys, rows_m[0]["arms"]["free_0"]["series"]["t"])
    v, se = sm[key]["vc"]
    tol = bm / 4.0
    lines.append(f"  inertia: centroid velocity of {key} in the source frame {v:+.4f} +- {se:.4f} vs beta/4 = {tol:.4f}")
    inst_ok = all(s_[k]["frame"] for s_ in (sr, sm) for k in keys) and \
        all(s_[k]["balance"] <= tol_bal for s_ in (sr, sm) for k in keys if k.startswith("free"))
    c = sr[key]
    if not inst_ok:
        status = "INSTRUMENT_FAILED"
    elif c["c1"] in ("growing", "collapsing") or c["c2"] in ("growing", "collapsing") or \
            (np.isfinite(v) and abs(v) - nse * se > tol):
        status = "FAILED"
    elif c["c1"] == "bound" and c["c2"] == "bound" and sm[key]["c1"] == "bound" and np.isfinite(v) and abs(v) + nse * se <= tol:
        status = "PASSED"
    else:
        status = "UNDECIDED"
    lines.append(f"instrument: frame check and free balances {'OK' if inst_ok else 'FAILS'}")
    return Result(status, lines, {"candidate": eta, "rest": _jsonable(sr), "move": _jsonable(sm)})
