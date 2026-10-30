"""
F3.7 -- Self-lensing: does a flux bundle that slows its own clocks hold its width?
==================================================================================

Why this test exists (23 Oct 2026)
----------------------------------
F3.6 (rest and move PASSED): flux that carries its timelike direction as
a state stays on the world line of its source and moves with it, while a
direction taken from the last step decays (one link kicks the direction
by ~0.45 in rapidity). The bundle is not bound: its width grows slowly
(2.9 -> 4.0 in 24 l). F0.6 fixed (FAILED): cones narrowed in the
coordinates of the laboratory pick a preferred frame. The question now:
what holds the width? The author's core idea, made covariant with the
carried direction: matter slows clocks (axiom 5) and paths bend towards
slower clocks (a = -grad ln Pi, the acceleration of free fall in the
static metric -Pi^2 dt^2 + ...; c = 1). A narrowed cone alone is
symmetric about the direction and has no inward bias; only bending can
pull the edges in. The test compares both, with the protocol the author
set (23 Oct): reference load on calibration sprinklings, three strengths
fixed in advance, a search on its own sprinklings, at most one candidate
chosen by a rule stated now, its parameters locked and confirmed on new
sprinklings, every branch reported.

Model (definitions locked before any run)
-----------------------------------------
Substrate and source as F3.6: 2+1 Poisson sprinkling, rho = 1 (lengths in
l), links (core.order.links); source world tube radius 4, source time
[0, 4), every source event emits s = 1 with the source 4-velocity u_src.
Transport: one causal pass over the time-sorted events. Every event v
with flux F_v = inflow + s_v has a carried direction u_v = the normalised
flux-weighted sum of the directions sent to it (+ s_v u_src).
Clock of an event (axiom 5, the F3 response):
    Pi_v = 1 - (1 - Pi_sat) min(F_v / C, 1)^m,   Pi_sat = 0.5, m = 2;
an event without flux has Pi = 1 (vacuum).
Strength: C = kappa F_ref, F_ref = the median flux per carrying event of
the free arm in the window [12, 28] (source time) on the calibration
sprinklings of the part (own seed stream); kappa in {4, 1, 1/4} (weak:
typical load 1/4, Pi ~0.97; middle: Pi ~0.5 at a typical carrier; strong:
nearly every carrier saturated, Pi = 0.5 inside, the gradient only at the
edge). Fixed now, never changed after a run.
Arms:
  free    F3.6 "carried": v sends all its flux to its child of smallest
          rapidity relative to u_v; the child receives u_v. No clock acts.
  lens    as free, but the direction sent to the child is bent towards
          slower clocks. Transverse gradient of ln Pi at v, estimated from
          v's PARENTS only (causal: all of them are already evaluated;
          parents without flux count with ln Pi = 0): parents u whose
          time lag in v's rest frame, (x_v - x_u) . u_v, is <= L_g = 2
          (two discreteness lengths, the scale of one step, tau ~1.15);
          their offsets projected onto the space orthogonal to u_v,
          p = d - (d . u_v) u_v, d = x_u - x_v, expressed in an orthonormal
          basis of that space; least squares ln Pi_u = a + G . q_u (a free
          intercept; no update with fewer than 4 such parents). The direction
          sent along the step v -> k (proper time tau_k) is u_v boosted
          by the rapidity |G| tau_k towards -G (a = -grad_perp ln Pi; the
          projection keeps the bending transverse). No free coefficient.
  narrow  the author's narrowed cone, covariant: v splits its flux
          uniformly among its children whose speed relative to u_v,
          tanh(rapidity), is <= Pi_v^2 (c_eff = Pi^2 in the rest frame of
          the matter, H1 / F7); if none, all to the child of smallest
          rapidity. Children receive u_v. No bending.
Every quantity is a Lorentz scalar of the events, of their flux and of
the initial u_src; criterion (c) of the axiom audit is checked bit for bit
(below). As in F3.6 and the whole F3 series, choosing among children uses
the future of v (criterion (b)) and the rules use the metric at the scale
of a link, not only the order and the counts (criterion (a)).

Readout (source frame, as F3.6): flux-weighted R_rms of the carrying
events about the source world line per slab of source time 2; slope s of
R_rms over the window W1 = [12, 28] and over W2 = [28, 44] (the extension:
a bound bundle must stay bound when the run is prolonged; T = 48).
Classification of an arm (s_free = the slope of the free arm on the same
sprinklings and window, mean over sprinklings; bootstrap SE; quarter rule):
    bound       |s| + 2 SE <= s_free / 4
    growing     s - 2 SE > s_free / 4
    collapsing  s + 2 SE < -s_free / 4     (reported separately, not bound)
    unresolved  otherwise
    wall-limited  the mean R_rms in the window exceeds half the box half-width (the width is limited by the
                  box, not by the rule; never bound)

Parts
-----
search   (RECORDED, exploratory) calibration (8 sprinklings) -> F_ref;
         32 search sprinklings, source at rest: free, lens x 3 kappa,
         narrow x 3 kappa; every arm recomputed in a frame boosted by
         beta' = 0.6 (coordinates and u_src) must choose the same children
         at every event (a lens or narrow arm that does not is excluded).
         Candidate rule (stated now): among the (arm, kappa) that are
         bound in W1 and in W2 and passed the frame check, the one with the
         smallest |s_W1| + 2 SE; ties: lens before narrow, then the larger
         kappa (weaker). At most one candidate; "none" is a result.
confirm  (decides) the candidate given at the command line (settable arm,
         kappa: CHOSEN FROM DATA of search, recorded), new streams for
         calibration and sprinklings; 32 sprinklings at rest + 32 with the
         source moving (beta = 0.5, box x [-16, 40]).
         PASSED      the candidate is bound in W1 and W2 at rest, the frame
                     check holds, and the moving bundle stays with its
                     source: |v_c| + 2 SE <= beta / 4 (centroid velocity in
                     the source frame over [5, 28]) and bound in W1;
         FAILED      growing or collapsing in W1 or W2 at rest, or
                     |v_c| - 2 SE > beta / 4;
         INSTRUMENT_FAILED  the frame check fails;
         UNDECIDED   otherwise.

Seen before this test was written: F3.6 (rest, move; quick runs), F0.6
fixed, F6.3b A. None of the three rules of lens and narrow had been run.
Quick runs of this file (after the definitions and criteria above were
written; 3 sprinklings, box 24 x 16 x 16, source radius 2; decide nothing):
free grows 0.07 / 0.11 per unit time (W1 / W2); lens with kappa 1 and 1/4
grows FASTER (0.37, 0.33 in W1): the bending is strong and noisy (mean
0.21-0.30 rapidity per step, Pi down to 0.5) and heats the bundle instead
of binding it; lens kappa 4 is close to free; narrow spreads (0.14, 0.39)
until the box stops it. The box made two narrow arms look "bound" in W2;
the wall-limited class was added after these runs (it can only remove a
false "bound", never create one). The rules and the criteria were not
changed.
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

TEST = {
    "id": "F3.7",
    "phase": "F3",
    "title": "Self-lensing: does a bundle of flow that slows its own clocks keep its width?",
    "theory_part": "matter slows clocks and paths bend towards slower clocks (the core of the theory), covariantly "
                   "with respect to the carried direction of the flow (F3.6); the author's narrowing of cones, "
                   "covariantly; binding of matter",
    "question": "Does the spreading of a bundle with a carried direction stop when the direction bends along the "
                "transverse gradient of ln Π estimated from the parents (lens), or when the flow divides only inside "
                "a cone narrowed to Π² relative to the carried direction (narrow)? Does the rule stay independent of "
                "the frame and the bundle inertial?",
    "method": "F3.6 (carried direction) + clocks from the load Π = 1 − (1 − Π_sat) min(F/C, 1)^m, C = κ F_ref, κ ∈ "
              "{4, 1, 1/4}, F_ref from calibration sprinklings; branches free, lens, narrow; slope of R_rms in the "
              "windows [12, 28] and [28, 44]; search selects at most one candidate by a rule fixed in advance, "
              "confirm checks it on new sprinklings including a moving source",
    "prediction": "according to the hypothesis of the author and the laboratory: lens binds (the edge paths turn "
                  "towards the slower clocks inside); narrow does not bind (a symmetric cone has no inward pull)",
    "criterion": "bound |s| + 2 SE <= s_free/4, growing s − 2 SE > s_free/4, collapse s + 2 SE < −s_free/4. search: "
                 "candidate = bound in both windows and with the same choice in the frame β' = 0.6, the smallest |s| "
                 "+ 2 SE (lens before narrow, then the weaker κ); otherwise none. confirm: PASSED bound in both "
                 "windows, frame check, |v_c| + 2 SE <= β/4 and bound also in motion; FAILED growing or collapse, or "
                 "|v_c| − 2 SE > β/4; INSTRUMENT the frame check; otherwise UNDECIDED.",
    "null_model": "free (F3.6 without clocks): the natural diffusion of the bundle",
    "level": "mechanism",
    "meters": ["K.6"],
    "depends_on": ["F3.6", "F0.6"],
    "params": {"D": 3, "T": 48.0, "half": 16.0, "r_src": 4.0, "t_src": 4.0, "slab": 2.0, "W1": [12.0, 28.0],
               "W2": [28.0, 44.0], "win_move": [5.0, 28.0], "beta_check": 0.6, "beta_move": 0.5, "Pi_sat": 0.5, "m": 2.0,
               "kappas": [4.0, 1.0, 0.25], "L_g": 2.0, "min_parents": 4, "n_cal": 8, "n_seeds": 32, "n_se": 2.0,
               "boot_B": 4000, "arm": None, "kappa": None},
    "parts": {"search": {"mode": "search"}, "confirm": {"mode": "confirm"}},
    "settable": {"confirm": ["arm", "kappa"]},
    "quick": {"n_seeds": 3, "n_cal": 1, "boot_B": 200, "T": 24.0, "half": 8.0, "r_src": 2.0, "t_src": 2.0,
              "W1": [6.0, 14.0], "W2": [14.0, 22.0], "win_move": [3.0, 14.0], "slab": 1.0},
    "estimate": "search ~15–30 min with 10 processes (7 branches + recomputations in another frame; links ~49 000 "
                "events); confirm similar (rest + motion)",
}


# ============================================================================================== geometry

def boost_points(X: np.ndarray, beta: float) -> np.ndarray:
    """Coordinates (or 4-vectors) in the frame moving with velocity beta along x (c = 1)."""
    if beta == 0.0:
        return np.array(X, float, copy=True)
    g = 1.0 / math.sqrt(1.0 - beta * beta)
    Y = np.array(X, float, copy=True)
    Y[..., 0] = g * (X[..., 0] - beta * X[..., 1])
    Y[..., 1] = g * (X[..., 1] - beta * X[..., 0])
    return Y


def mink(a, b):
    """Minkowski product (+ - -) over the last axis."""
    return a[..., 0] * b[..., 0] - (a[..., 1:] * b[..., 1:]).sum(axis=-1)


def children_csr(ptr, idx, N):
    child = np.repeat(np.arange(N), np.diff(ptr))
    o = np.argsort(idx, kind="stable")
    cptr = np.concatenate([[0], np.cumsum(np.bincount(idx, minlength=N))]).astype(np.int64)
    return cptr, child[o]


def rest_basis(u: np.ndarray) -> np.ndarray:
    """Orthonormal spacelike basis (rows e_i, e_i . e_i = -1, e_i . u = 0) of the rest space of the unit vector u (Gram-Schmidt)."""
    D = u.size
    out = []
    for k in range(1, D):
        e = np.zeros(D)
        e[k] = 1.0
        e = e - mink(e, u) * u                      # remove the component along u (u . u = 1)
        for f in out:
            e = e + mink(e, f) * f                  # f . f = -1: e - (e.f)/(f.f) f = e + (e.f) f
        n = math.sqrt(max(-float(mink(e, e)), 0.0))
        out.append(e / n)
    return np.array(out)


def boost_towards(u: np.ndarray, e: np.ndarray, eta: float) -> np.ndarray:
    """u boosted by rapidity eta along the unit spacelike vector e orthogonal to u: cosh(eta) u + sinh(eta) e."""
    return math.cosh(eta) * u + math.sinh(eta) * e


# ============================================================================================== transport

def clock(F, C, Pi_sat, m):
    return 1.0 - (1.0 - Pi_sat) * min(F / C, 1.0) ** m if F > 0.0 else 1.0


def transport(X, ptr, idx, cptr, cidx, s, u_src, arm, C, p):
    """
    One causal pass (module notes). arm: 'free', 'lens', 'narrow'. Returns the flux of every event, the decisions
    (for the frame check: the chosen child, or the number of allowed children for narrow), the clocks and the
    bending rapidities of the lens arm.
    """
    N, D = X.shape
    Pi_sat, m = float(p["Pi_sat"]), float(p["m"])
    L_g, nmin = float(p["L_g"]), int(p["min_parents"])
    inflow = np.zeros(N)
    U = np.zeros((N, D))
    flux = np.zeros(N)
    lnPi = np.zeros(N)
    decision = np.full(N, -1, dtype=np.int64)
    bend = np.full(N, np.nan)
    for v in range(N):
        tot = inflow[v] + s[v]
        if tot <= 0.0:
            continue
        P = U[v] + s[v] * u_src
        u = P / math.sqrt(max(float(mink(P, P)), 1e-300))
        flux[v] = tot
        Pi_v = clock(tot, C, Pi_sat, m) if arm != "free" else 1.0
        lnPi[v] = math.log(Pi_v)
        a, b = cptr[v], cptr[v + 1]
        if b == a:
            continue
        ks = cidx[a:b]
        B = X[ks] - X[v]
        tau = np.sqrt(np.maximum(mink(B, B), 1e-300))
        ch = mink(u[None, :], B) / tau                       # cosh rapidity of each step relative to u
        if arm == "narrow":
            speed = np.sqrt(np.maximum(1.0 - 1.0 / (ch * ch), 0.0))
            ok = speed <= Pi_v ** 2
            if ok.any():
                sel = ks[ok]
                share = tot / sel.size
                np.add.at(inflow, sel, share)
                np.add.at(U, sel, share * u)
                decision[v] = int(sel.size)
            else:
                j = int(np.argmin(ch))
                inflow[ks[j]] += tot
                U[ks[j]] += tot * u
                decision[v] = -2 - int(ks[j])                # fallback, encoded
            continue
        j = int(np.argmin(ch))
        k = int(ks[j])
        decision[v] = k
        send = u
        if arm == "lens":
            pa, pb = ptr[v], ptr[v + 1]
            par = idx[pa:pb]
            if par.size:
                d = X[par] - X[v]                             # past-directed offsets
                lag = -mink(d, u[None, :])                    # time lag of each parent in v's rest frame (> 0)
                sel = lag <= L_g
                if int(sel.sum()) >= nmin:
                    dd = d[sel]
                    pp = dd - mink(dd, u[None, :])[:, None] * u[None, :]     # projection orthogonal to u
                    E = rest_basis(u)
                    q = -np.stack([mink(pp, e[None, :]) for e in E], axis=1)  # coordinates in the rest space
                    A = np.column_stack([np.ones(q.shape[0]), q])
                    coef, *_ = np.linalg.lstsq(A, lnPi[par[sel]], rcond=None)
                    G = coef[1:]                              # gradient of ln Pi in the rest space of u
                    gn = float(np.sqrt((G * G).sum()))
                    if gn > 0.0:
                        e_dir = -(G[0] * E[0] + (G[1] * E[1] if E.shape[0] > 1 else 0.0)) / gn
                        eta = gn * float(tau[j])
                        send = boost_towards(u, e_dir, eta)
                        bend[v] = eta
        inflow[k] += tot
        U[k] += tot * send
    return {"flux": flux, "decision": decision, "lnPi": lnPi, "bend": bend}


# ============================================================================================== readout

def readout(Xs, flux, p, beta):
    T, sl = float(p["T"]), float(p["slab"])
    t_end = T * math.sqrt(1.0 - beta * beta) - sl
    edges = np.arange(float(p["t_src"]), t_end + 1e-9, sl)
    out = {k: [] for k in ("t", "R", "xc", "n", "F")}
    car = flux > 0.0
    for a, b in zip(edges[:-1], edges[1:]):
        m = car & (Xs[:, 0] >= a) & (Xs[:, 0] < b)
        w = flux[m]
        out["t"].append(0.5 * (a + b))
        out["n"].append(int(m.sum()))
        out["F"].append(float(w.sum()))
        if w.sum() > 0:
            out["R"].append(float(np.sqrt((w * (Xs[m, 1:] ** 2).sum(axis=1)).sum() / w.sum())))
            out["xc"].append(float((w * Xs[m, 1]).sum() / w.sum()))
        else:
            out["R"].append(float("nan"))
            out["xc"].append(float("nan"))
    return out


def slope(t, y, win):
    t, y = np.asarray(t, float), np.asarray(y, float)
    m = (t >= win[0]) & (t <= win[1]) & np.isfinite(y)
    return float(np.polyfit(t[m], y[m], 1)[0]) if m.sum() >= 2 else float("nan")


# ============================================================================================== tasks

def setup(p, seed, beta):
    D, T, H = int(p["D"]), float(p["T"]), float(p["half"])
    rng = np.random.default_rng(seed)
    X = sprinkle_poisson(Box(D, [0.0, -H, -H], [T, H + beta * T, H]), 1.0, rng)
    ptr, idx = links_csr(X, 1.0)
    cptr, cidx = children_csr(ptr, idx, len(X))
    g = 1.0 / math.sqrt(1.0 - beta * beta)
    u_src = np.array([g, g * beta, 0.0])
    Xs = boost_points(X, beta)
    tube = (Xs[:, 0] >= 0.0) & (Xs[:, 0] < float(p["t_src"])) & ((Xs[:, 1:] ** 2).sum(axis=1) <= float(p["r_src"]) ** 2)
    return X, ptr, idx, cptr, cidx, Xs, tube.astype(float), u_src


def calibration_task(p, seed):
    """Flux per carrying event of the free arm in W1 (source at rest): the reference load."""
    X, ptr, idx, cptr, cidx, Xs, s, u_src = setup(p, seed, 0.0)
    r = transport(X, ptr, idx, cptr, cidx, s, u_src, "free", 1.0, p)
    w1 = p["W1"]
    m = (r["flux"] > 0) & (Xs[:, 0] >= w1[0]) & (Xs[:, 0] <= w1[1])
    return r["flux"][m]


def run_task(p, seed, beta, arms, F_ref):
    """One sprinkling: every (arm, kappa) in arms, readout in the source frame, frame check (beta')."""
    X, ptr, idx, cptr, cidx, Xs, s, u_src = setup(p, seed, beta)
    bc = float(p["beta_check"])
    Xb = boost_points(X, bc)
    ub = boost_points(u_src[None, :], bc)[0]
    out = {"N": len(X), "sources": int(s.sum()), "arms": {}}
    for arm, kappa in arms:
        C = float(kappa) * F_ref if kappa is not None else 1.0
        r = transport(X, ptr, idx, cptr, cidx, s, u_src, arm, C, p)
        rb = transport(Xb, ptr, idx, cptr, cidx, s, ub, arm, C, p)
        ser = readout(Xs, r["flux"], p, beta)
        key = arm if kappa is None else f"{arm}_{kappa:g}"
        out["arms"][key] = {"series": ser, "s1": slope(ser["t"], ser["R"], p["W1"]), "s2": slope(ser["t"], ser["R"], p["W2"]),
                            "vc": slope(ser["t"], ser["xc"], p["win_move"]),
                            "frame_identical": bool(np.array_equal(r["decision"], rb["decision"])),
                            "bend_mean": float(np.nanmean(r["bend"])) if np.isfinite(r["bend"]).any() else float("nan"),
                            "Pi_min": float(np.exp(r["lnPi"].min()))}
    return out


# ============================================================================================== statistics

def classify(sv, se, sfree, nse):
    if not all(np.isfinite([sv, se, sfree])):
        return "unresolved"
    tol = sfree / 4.0
    if abs(sv) + nse * se <= tol:
        return "bound"
    if sv - nse * se > tol:
        return "growing"
    if sv + nse * se < -tol:
        return "collapsing"
    return "unresolved"


def summarise(rows, keys, p, ctx, tag):
    B, nse = int(p["boot_B"]), float(p["n_se"])
    st = {}
    for key in keys:
        bs = lambda *t, key=key: derive_seed(ctx.seed, "F3.7-boot", tag, key, *t)
        st[key] = {w: boot_se([r_["arms"][key][w] for r_ in rows], B, bs(w))[:2] for w in ("s1", "s2", "vc")}
        st[key]["frame"] = all(r_["arms"][key]["frame_identical"] for r_ in rows)
        with np.errstate(all="ignore"):
            st[key]["R"] = np.nanmean(np.array([r_["arms"][key]["series"]["R"] for r_ in rows], float), axis=0)
            st[key]["n"] = np.nanmean(np.array([r_["arms"][key]["series"]["n"] for r_ in rows], float), axis=0)
            st[key]["F"] = np.nanmean(np.array([r_["arms"][key]["series"]["F"] for r_ in rows], float), axis=0)
            st[key]["bend"] = float(np.nanmean([r_["arms"][key]["bend_mean"] for r_ in rows]))
            st[key]["Pi_min"] = float(np.nanmin([r_["arms"][key]["Pi_min"] for r_ in rows]))
    sfree1, sfree2 = st["free"]["s1"][0], st["free"]["s2"][0]
    t = np.asarray(rows[0]["arms"]["free"]["series"]["t"], float)
    wall = float(p["half"]) / 2.0
    for key in keys:
        st[key]["c1"] = classify(*st[key]["s1"], sfree1, nse)
        st[key]["c2"] = classify(*st[key]["s2"], sfree2, nse)
        for c, w in (("c1", p["W1"]), ("c2", p["W2"])):     # a width that reached the box cannot be read as bound
            m = (t >= w[0]) & (t <= w[1])
            if m.any() and np.nanmax(st[key]["R"][m]) > wall and st[key][c] in ("bound", "unresolved", "collapsing"):
                st[key][c] = "wall-limited"
    return st


def lines_for(st, keys, t):
    out = []
    for key in keys:
        x = st[key]
        out.append(f"  {key:12s}: slope W1 {x['s1'][0]:+.4f} +- {x['s1'][1]:.4f} -> {x['c1']}; W2 {x['s2'][0]:+.4f} +- "
                   f"{x['s2'][1]:.4f} -> {x['c2']}; frame check {'OK' if x['frame'] else 'FAILS'}"
                   + (f"; mean bending {x['bend']:.4f}" if np.isfinite(x["bend"]) else "")
                   + (f"; min Pi {x['Pi_min']:.3f}" if key != "free" else ""))
        out.append("        R_rms  " + " ".join(f"{v:5.2f}" if np.isfinite(v) else "  nan" for v in x["R"]))
        out.append("        carry  " + " ".join(f"{v:5.1f}" for v in x["n"]))
    out.insert(0, "        t      " + " ".join(f"{v:5.1f}" for v in t))
    return out


# ============================================================================================== run

def run(p, ctx):
    mode = str(p["mode"])
    t0 = time.time()
    n = int(p["n_seeds"])
    if mode == "confirm":
        arm, kappa = p.get("arm"), p.get("kappa")
        if arm not in ("lens", "narrow") or kappa is None:
            raise SystemExit("F3.7 confirm: give the candidate of the search, e.g. --set arm=lens --set kappa=1")
        kappa = float(kappa)
        if kappa not in [float(k) for k in p["kappas"]]:
            raise SystemExit(f"F3.7 confirm: kappa must be one of the grid {p['kappas']} (got {kappa})")
    stream = "F3.7-" + mode
    cal = parallel_map(calibration_task, [(p, derive_seed(ctx.seed, stream + "-cal", k)) for k in range(int(p["n_cal"]))],
                       ctx.workers, ctx.progress, lambda i, r_, d, nn: f"[F3.7] calibration {d}/{nn} ({time.time()-t0:.0f}s)")
    F_ref = float(np.median(np.concatenate(cal)))
    head = [f"2+1 sprinkling rho = 1, box t [0, {p['T']:g}], half width {p['half']:g}; source radius {p['r_src']:g}, time "
            f"[0, {p['t_src']:g}); clocks Pi = 1 - (1 - {p['Pi_sat']}) min(F/C, 1)^{p['m']:g}, C = kappa F_ref, F_ref = "
            f"{F_ref:.3f} (median flux per carrier of the free arm in W1, {p['n_cal']} calibration sprinklings); lens: "
            f"parents with rest-frame lag <= {p['L_g']:g}, >= {p['min_parents']} parents; windows W1 {p['W1']}, W2 {p['W2']}; "
            f"frame check beta' = {p['beta_check']}"]
    if mode == "search":
        arms = [("free", None)] + [(a, float(k)) for a in ("lens", "narrow") for k in p["kappas"]]
        keys = ["free"] + [f"{a}_{float(k):g}" for a in ("lens", "narrow") for k in p["kappas"]]
        seeds = [derive_seed(ctx.seed, stream, k) for k in range(n)]
        rows = parallel_map(run_task, [(p, sd, 0.0, arms, F_ref) for sd in seeds], ctx.workers, ctx.progress,
                            lambda i, r_, d, nn: f"[F3.7] search sprinkling {d}/{nn} ({time.time()-t0:.0f}s)")
        st = summarise(rows, keys, p, ctx, "search")
        lines = head + [f"search: {n} sprinklings, source at rest (~{np.mean([r_['N'] for r_ in rows]):.0f} events, "
                        f"{np.mean([r_['sources'] for r_ in rows]):.0f} source events)"]
        lines += lines_for(st, keys, rows[0]["arms"]["free"]["series"]["t"])
        nse = float(p["n_se"])
        cands = []
        for key in keys[1:]:
            x = st[key]
            if x["c1"] == "bound" and x["c2"] == "bound" and x["frame"]:
                arm, kap = key.split("_")
                cands.append((abs(x["s1"][0]) + nse * x["s1"][1], 0 if arm == "lens" else 1, -float(kap), key))
        cand = sorted(cands)[0][3] if cands else None
        lines.append(f"candidate (rule stated in the header): {cand if cand else 'none'}"
                     + (f"  -> confirm with: --set arm={cand.split('_')[0]} --set kappa={cand.split('_')[1]}" if cand else ""))
        return Result("RECORDED", lines, {"F_ref": F_ref, "candidate": cand, "stats": _jsonable(st)})

    # confirm
    arms = [("free", None), (arm, kappa)]
    key = f"{arm}_{kappa:g}"
    keys = ["free", key]
    seeds_r = [derive_seed(ctx.seed, stream + "-rest", k) for k in range(n)]
    seeds_m = [derive_seed(ctx.seed, stream + "-move", k) for k in range(n)]
    bm = float(p["beta_move"])
    rows_r = parallel_map(run_task, [(p, sd, 0.0, arms, F_ref) for sd in seeds_r], ctx.workers, ctx.progress,
                          lambda i, r_, d, nn: f"[F3.7] confirm rest {d}/{nn} ({time.time()-t0:.0f}s)")
    rows_m = parallel_map(run_task, [(p, sd, bm, arms, F_ref) for sd in seeds_m], ctx.workers, ctx.progress,
                          lambda i, r_, d, nn: f"[F3.7] confirm move {d}/{nn} ({time.time()-t0:.0f}s)")
    sr = summarise(rows_r, keys, p, ctx, "confirm-rest")
    sm = summarise(rows_m, keys, p, ctx, "confirm-move")
    nse = float(p["n_se"])
    lines = head + [f"confirm: candidate {key} (chosen from the search, CHOSEN FROM DATA); {n} sprinklings at rest, {n} moving "
                    f"(beta = {bm})", "  at rest:"]
    lines += lines_for(sr, keys, rows_r[0]["arms"]["free"]["series"]["t"])
    lines.append(f"  moving source (beta = {bm}):")
    lines += lines_for(sm, keys, rows_m[0]["arms"]["free"]["series"]["t"])
    v, se = sm[key]["vc"]
    tol = bm / 4.0
    lines.append(f"  inertia: centroid velocity of {key} in the source frame {v:+.4f} +- {se:.4f} vs beta/4 = {tol:.4f}")
    frame_ok = sr[key]["frame"] and sm[key]["frame"]
    c = sr[key]
    if not frame_ok:
        status = "INSTRUMENT_FAILED"
    elif c["c1"] in ("growing", "collapsing") or c["c2"] in ("growing", "collapsing") or \
            (np.isfinite(v) and abs(v) - nse * se > tol):
        status = "FAILED"
    elif c["c1"] == "bound" and c["c2"] == "bound" and sm[key]["c1"] == "bound" and np.isfinite(v) and abs(v) + nse * se <= tol:
        status = "PASSED"
    else:
        status = "UNDECIDED"
    lines.append(f"instrument: frame check (beta' = {p['beta_check']}) {'OK' if frame_ok else 'FAILS'}")
    return Result(status, lines, {"F_ref": F_ref, "candidate": key, "rest": _jsonable(sr), "move": _jsonable(sm)})


def _jsonable(st):
    return {k: {kk: (list(vv) if isinstance(vv, (tuple, np.ndarray)) else vv) for kk, vv in v.items()} for k, v in st.items()}
