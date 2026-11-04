"""
F3.6 -- Direction memory: does flux that carries its own timelike direction (momentum) stay localized?
=====================================================================================================

Why this test exists (22-23 Oct 2026)
-------------------------------------
Every rule of the F3 series that uses only the order and no memory of a
direction (uniform split of the flux among the children, the copied
pulse) spreads the influence of a saturated region over the whole future
cone, at the speed of light (F3.3g A). The only compact domain the lab has
seen (F6.3b A, 23 Oct) came from cones narrowed in the coordinates of the
laboratory, i.e. from a timelike direction borrowed from the lab frame,
and the quick runs of F0.6 show it is frame-dependent. In a Lorentz-
invariant substrate no event has a preferred direction of its own (the
children of an event are distributed alike in every frame; Bombelli,
Henson, Sorkin 2006: no covariant local choice of a direction). A
localized object is a world tube, so it must know which timelike direction
it goes, and the only covariant source of that is its own past. The
hypothesis of the laboratory (23 Oct): a localized structure in a Lorentz-
invariant substrate needs flux that carries the memory of its timelike
direction; that memory IS momentum, so matter and momentum are one
problem, not two.

Model (pure transport: no capacity, no pulse, no threshold, no drive
after the ignition; the question is the spreading itself)
--------------------------------------------------------------------------
Fixed causal set: a Poisson sprinkling of a 2+1 box (units: rho = 1, so
lengths are in the discreteness length l), links by core.order.links. A
source (the "ignition"): every event of the source world tube (radius 4,
source proper time 0-4) emits s = 1 with the source 4-velocity u_src.
Three arms, all parameter-free:
  uniform  every event splits its flux equally among its children (the F3
           routing, core.flux.capacity.route_capacity with C = infinity);
  step     every event v sends all its flux to the ONE child k whose step
           v -> k has the smallest rapidity relative to v's direction
           u_v; the direction of k is the flux-weighted sum of the unit
           step vectors arriving at it (momentum = the displacement of the
           last step);
  carried  the same choice of the child, but the direction is a carried
           state: u_k = the flux-weighted sum of the directions u_v of the
           parents that send to k (momentum carried with the flux,
           combined when streams merge; a source event adds s u_src).
Rapidity between a unit timelike vector u and a step b (both future-
timelike): cosh eta = (u . b) / tau_b with the Minkowski product
u . b = u_0 b_0 - u_vec . b_vec and tau_b = sqrt(b . b). Every quantity of
the step and carried rules is a Lorentz scalar of the events and of the
initial u_src (a covariant initial condition: it transforms with the
source). They use the embedding (the metric at the scale of one link), not
only the order and the counts: criterion (a) of the axiom audit is not
met, criterion (c) is (checked bit for bit below). Choosing among the
children uses the future of v (the children exist later), as the uniform
split of the whole F3 series does; criterion (b) is not met by any arm.

Development, stated before the criteria (22-23 Oct 2026, one sprinkling,
decide nothing): an event has ~48 children (links) in 2+1 at rho = 1; the
child of smallest rapidity relative to a fixed frame has rapidity 0.45
(median 0.44) and proper time ~1.15. So a single link step kicks the
direction by ~0.45 in rapidity. Prediction from that number: in the step
arm the direction random-walks (~0.45 per step), becomes relativistic
within a handful of steps and the flux spreads at nearly the speed of
light; in the carried arm the direction does not change and the position
jitters around the world line of u_src (diffusive growth of the width).

Quick runs of this file (after the criteria were written; 3 sprinklings, box
16 x 16 x 16, source radius 2, decide nothing): uniform spreads at 0.84 and
0.76; step 0.37 +- 0.13 and 0.61 +- 0.03; carried 0.20 +- 0.06 and
0.16 +- 0.07 (R_rms stays ~1.5-2.5 while uniform grows to the walls); the
carrying events coalesce (from ~12 to ~5 per slab); the moving carried
bundle stays with the source (v_c -0.02 +- 0.04); the boosted
recomputation chose identically. Changed after them: the readout slabs end
where the source world line leaves the box (T / g - 2), because late slabs
of the moving part mixed in the box edge; no criterion changed.

Parts
-----
rest  the source at rest (u_src = (1, 0, 0)); box t [0, 32], x, y [-16, 16]
      (~33 000 events).
move  the source moving with beta = 0.5 along x (u_src = g (1, 0.5, 0),
      the source world tube boosted); box x [-16, 32]. Readout in the
      source frame (events Lorentz-transformed by the source velocity).
Instrument in both parts: the step and carried arms computed again with
all coordinates and u_src boosted by beta' = 0.6 (a different frame of the
computation) must choose the same child at every event (rapidities are
invariant; a mismatch would mean a frame-dependent rule or a numerical
near-tie); the uniform arm must spread (below), otherwise the readout
cannot tell spreading from localization.

Observables (per sprinkling, per arm, slabs of source proper time 2 from
4 up to T / g - 2, where the source world line leaves the sprinkled box): flux-weighted R_rms of the carrying events about the source
world line (in the source frame), flux-weighted centroid, number of
carrying events, flux left in the box. Spreading rate v_s = least-squares
slope of R_rms against source time over the window [5, 11] (before the
uniform arm reaches the box walls). Mean step rapidity relative to the
carried direction (step and carried arms).

Criterion (stated before any full run)
--------------------------------------
Alternatives for an arm: spreading at the speed of light (v_s ~ 1, the
uniform behaviour) against localization (v_s ~ 0). Quarter rule:
  localized   v_s + 2 SE <= 1/4;    spreading   v_s - 2 SE > 1/4.
INSTRUMENT_FAILED if the uniform arm is not spreading or the boosted
recomputation chooses differently at any event.
rest:  PASSED if the carried arm is localized; FAILED if it is spreading;
       UNDECIDED otherwise. The step arm is reported (prediction above).
move:  inertia of the carried arm: the centroid velocity v_c in the source
       frame (slope of the centroid x' over [5, 24]); alternatives: moves
       with the source (v_c = 0) or stays at rest in the lab (v_c = -beta).
       PASSED if the carried arm is localized and |v_c| + 2 SE <= beta / 4;
       FAILED if it is spreading or |v_c| - 2 SE > beta / 4; UNDECIDED otherwise.

Reading declared now. A localized carried arm shows that a covariant rule
CAN hold flux on a world line, provided momentum is a carried state, and
a spreading step arm that momentum cannot be the displacement between
neighbouring events (the discreteness kicks it by ~0.45 per link). It does
NOT show a self-bound object: the carried bundle persists because the flux
is conserved and keeps its direction; nothing yet holds its width (the
next question: what binds it). If the carried arm spreads too, direction
memory is not enough for localization in this substrate.
"""

from __future__ import annotations

import math
import time

import numpy as np

from core.flux.capacity import route_capacity
from core.order.links import links_csr
from core.order.sprinkle import Box, derive_seed, sprinkle_poisson
from core.stats.parallel import parallel_map
from core.stats.stats import boot_se
from core.testing import Result

TEST = {
    "id": "F3.6",
    "phase": "F3",
    "title": "Direction memory: does a flow that carries its own timelike direction (momentum) stay localised?",
    "theory_part": "matter and momentum (laboratory hypothesis of 23 Oct: localisation in a Lorentz-invariant "
                   "substrate requires a flow with a memory of direction); existence gate, inertia",
    "question": "Does the flow from a source spread at the speed of light (uniform division), or does it stay at the "
                "world line of the source when every event sends the flow to the descendant nearest to its direction "
                "and the direction is (a) taken from the last step, (b) carried as a state of the flow? Does the "
                "bundle move with the source (inertia)?",
    "method": "2+1 sprinkling (ρ = 1), links; source (radius 4, time 0–4) with 4-velocity u_src; three branches "
              "without parameters: uniform division, step (direction = last step), carried direction (direction = "
              "flow-weighted sum of the parents' directions); spreading speed of R_rms in the window [5, 11]; moving "
              "source β = 0.5: velocity of the centroid in the frame of the source; check: recomputation in a frame "
              "β' = 0.6 must choose the same descendants",
    "prediction": "uniform and step: spreading almost at the speed of light (a kick of ~0.45 in rapidity per link); "
                  "carried direction: localisation with diffusive spreading and motion with the source",
    "criterion": "localised branch v_s + 2 SE <= 1/4, spreading v_s − 2 SE > 1/4. INSTRUMENT: the uniform branch "
                 "spreads and the recomputation in another frame chooses the same. rest: PASSED carried direction "
                 "localised, FAILED spreads, otherwise UNDECIDED. move: PASSED carried direction localised and |v_c| "
                 "+ 2 SE <= β/4, FAILED spreads or |v_c| − 2 SE > β/4, otherwise UNDECIDED.",
    "null_model": "uniform division of the flow (F3 series) and a memory of the last step only",
    "level": "mechanism",
    "meters": ["K.6"],
    "depends_on": ["F3.3g", "F0.6"],
    "params": {"D": 3, "T": 32.0, "half": 16.0, "r_src": 4.0, "t_src": 4.0, "slab": 2.0, "win": [5.0, 11.0],
               "win_move": [5.0, 24.0], "beta_check": 0.6, "n_se": 2.0, "boot_B": 4000, "n_seeds": 32},
    "parts": {"rest": {"mode": "rest", "beta": 0.0}, "move": {"mode": "move", "beta": 0.5}},
    "quick": {"n_seeds": 3, "boot_B": 200, "T": 16.0, "half": 8.0, "r_src": 2.0, "t_src": 2.0, "win": [3.0, 6.0],
              "win_move": [3.0, 12.0], "slab": 1.0},
    "estimate": "rest ~3–5 min, move ~5–8 min with 10 processes (one sprinkling ~1–1.5 min, the links dominate; the "
                "branches are cheap)",
}

ARMS = ("uniform", "step", "carried")


# ============================================================================================== geometry

def boost_points(X: np.ndarray, beta: float) -> np.ndarray:
    """Coordinates in the frame moving with velocity beta along x (c = 1); also used for 4-vectors."""
    if beta == 0.0:
        return np.array(X, float, copy=True)
    g = 1.0 / math.sqrt(1.0 - beta * beta)
    Y = np.array(X, float, copy=True)
    Y[..., 0] = g * (X[..., 0] - beta * X[..., 1])
    Y[..., 1] = g * (X[..., 1] - beta * X[..., 0])
    return Y


def mink(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Minkowski product a . b = a_0 b_0 - a_vec . b_vec (signature + - -), broadcast over leading axes."""
    return a[..., 0] * b[..., 0] - (a[..., 1:] * b[..., 1:]).sum(axis=-1)


def children_csr(ptr: np.ndarray, idx: np.ndarray, N: int):
    """Children lists from the parents CSR of core.order.links: (cptr, cidx)."""
    child = np.repeat(np.arange(N), np.diff(ptr))
    o = np.argsort(idx, kind="stable")
    cptr = np.concatenate([[0], np.cumsum(np.bincount(idx, minlength=N))]).astype(np.int64)
    return cptr, child[o]


# ============================================================================================== rules

def route_direction(X: np.ndarray, cptr: np.ndarray, cidx: np.ndarray, s: np.ndarray, u_src: np.ndarray, mode: str):
    """
    The step / carried rule (module notes) in one causal pass over the time-sorted events (a linear extension).
    Returns the flux of every event, the chosen child (-1 if none or no flux), the direction of every carrying
    event and the rapidity of every chosen step relative to the direction.
    """
    N, D = X.shape
    inflow = np.zeros(N)
    U = np.zeros((N, D))                        # flux-weighted momentum sum arriving at each event
    flux = np.zeros(N)
    choice = np.full(N, -1, dtype=np.int64)
    eta = np.full(N, np.nan)
    udir = np.zeros((N, D))
    for v in range(N):
        tot = inflow[v] + s[v]
        if tot <= 0.0:
            continue
        P = U[v] + s[v] * u_src
        norm = math.sqrt(max(float(mink(P, P)), 0.0))
        u = P / norm
        udir[v] = u
        flux[v] = tot
        a, b = cptr[v], cptr[v + 1]
        if b == a:
            continue                              # no child: the flux leaves the sprinkled region
        ks = cidx[a:b]
        B = X[ks] - X[v]
        tau = np.sqrt(np.maximum(mink(B, B), 1e-300))
        ch = mink(u[None, :], B) / tau            # cosh of the rapidity of each step relative to u
        j = int(np.argmin(ch))
        k = int(ks[j])
        choice[v] = k
        eta[v] = math.acosh(max(float(ch[j]), 1.0))
        inflow[k] += tot
        U[k] += tot * (u if mode == "carried" else B[j] / tau[j])
    return {"flux": flux, "choice": choice, "dir": udir, "eta": eta}


# ============================================================================================== readout

def slab_readout(Xs: np.ndarray, flux: np.ndarray, p):
    """Per slab of source time (columns of Xs in the source frame): flux-weighted R_rms, centroid x, carriers, flux."""
    T, sl, beta = float(p["T"]), float(p["slab"]), float(p["beta"])
    # the source world line stays inside the sprinkled box up to source time T / g; the last slab before that
    # is dropped too (events at the top of the box have no children and keep their flux)
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
            r2 = (Xs[m, 1:] ** 2).sum(axis=1)
            out["R"].append(float(np.sqrt((w * r2).sum() / w.sum())))
            out["xc"].append(float((w * Xs[m, 1]).sum() / w.sum()))
        else:
            out["R"].append(float("nan"))
            out["xc"].append(float("nan"))
    return out


def slope(t, y, win):
    """Least-squares slope of y(t) over the window (finite points only)."""
    t, y = np.asarray(t, float), np.asarray(y, float)
    m = (t >= win[0]) & (t <= win[1]) & np.isfinite(y)
    if m.sum() < 2:
        return float("nan")
    return float(np.polyfit(t[m], y[m], 1)[0])


# ============================================================================================== task

def seed_task(p, seed):
    """One sprinkling: the three arms, the readout in the source frame, the boosted recomputation check."""
    D, T, H = int(p["D"]), float(p["T"]), float(p["half"])
    beta = float(p["beta"])
    g = 1.0 / math.sqrt(1.0 - beta * beta)
    u_src = np.array([g, g * beta, 0.0])
    rng = np.random.default_rng(seed)
    X = sprinkle_poisson(Box(D, [0.0, -H, -H], [T, H + beta * T, H]), 1.0, rng)
    N = len(X)
    ptr, idx = links_csr(X, 1.0)
    cptr, cidx = children_csr(ptr, idx, N)
    Xs = boost_points(X, beta)                                    # source frame (readout)
    in_tube = (Xs[:, 0] >= 0.0) & (Xs[:, 0] < float(p["t_src"])) & ((Xs[:, 1:] ** 2).sum(axis=1) <= float(p["r_src"]) ** 2)
    s = in_tube.astype(float)
    out = {"N": N, "links": int(idx.size), "sources": int(in_tube.sum()), "arms": {}, "check_identical": True}
    uni = route_capacity(ptr, idx, s, np.inf, np.zeros(N, dtype=bool))
    res = {"uniform": {"flux": uni["accepted"]}}
    bc = float(p["beta_check"])
    Xb = boost_points(X, bc)
    ub = boost_points(u_src[None, :], bc)[0]
    for mode in ("step", "carried"):
        r = route_direction(X, cptr, cidx, s, u_src, mode)
        rb = route_direction(Xb, cptr, cidx, s, ub, mode)            # same events, another frame of the computation
        out["check_identical"] &= bool(np.array_equal(r["choice"], rb["choice"]))
        res[mode] = r
    for arm in ARMS:
        ser = slab_readout(Xs, res[arm]["flux"], p)
        rec = {"series": ser, "v_s": slope(ser["t"], ser["R"], p["win"]),
               "v_c": slope(ser["t"], ser["xc"], p["win_move"])}
        if arm != "uniform":
            e = res[arm]["eta"]
            rec["eta_mean"] = float(np.nanmean(e)) if np.isfinite(e).any() else float("nan")
        out["arms"][arm] = rec
    return out


# ============================================================================================== run

def classify(v, se, nse):
    if not (np.isfinite(v) and np.isfinite(se)):
        return "undefined"
    if v + nse * se <= 0.25:
        return "localized"
    if v - nse * se > 0.25:
        return "spreading"
    return "inconclusive"


def run(p, ctx):
    t0 = time.time()
    n, B, nse = int(p["n_seeds"]), int(p["boot_B"]), float(p["n_se"])
    seeds = [derive_seed(ctx.seed, "F3.6", ctx.part, k) for k in range(n)]
    rows = parallel_map(seed_task, [(p, sd) for sd in seeds], ctx.workers, ctx.progress,
                        lambda i, r_, d, nn: f"[F3.6] {ctx.part} sprinkling {d}/{nn}: v_s uniform "
                                             f"{r_['arms']['uniform']['v_s']:.3f} step {r_['arms']['step']['v_s']:.3f} "
                                             f"carried {r_['arms']['carried']['v_s']:.3f} ({time.time()-t0:.0f}s)")
    beta = float(p["beta"])
    lines = [f"2+1 sprinkling rho = 1 (lengths in l), box t [0, {p['T']:g}], x [-{p['half']:g}, {p['half'] + beta * p['T']:g}], "
             f"y [-{p['half']:g}, {p['half']:g}]; source radius {p['r_src']:g}, source time [0, {p['t_src']:g}), 4-velocity "
             f"beta = {beta:g}; {n} sprinklings (~{np.mean([r_['N'] for r_ in rows]):.0f} events, "
             f"{np.mean([r_['sources'] for r_ in rows]):.0f} source events); readout in the source frame"]
    stats = {}
    for arm in ARMS:
        bs = lambda *tag, arm=arm: derive_seed(ctx.seed, "F3.6-boot", ctx.part, arm, *tag)
        vs = boot_se([r_["arms"][arm]["v_s"] for r_ in rows], B, bs("vs"))
        vc = boot_se([r_["arms"][arm]["v_c"] for r_ in rows], B, bs("vc"))
        cls = classify(vs[0], vs[1], nse)
        with np.errstate(all="ignore"):
            R = np.nanmean(np.array([r_["arms"][arm]["series"]["R"] for r_ in rows], float), axis=0)
            nc = np.nanmean(np.array([r_["arms"][arm]["series"]["n"] for r_ in rows], float), axis=0)
            F = np.nanmean(np.array([r_["arms"][arm]["series"]["F"] for r_ in rows], float), axis=0)
            xc = np.nanmean(np.array([r_["arms"][arm]["series"]["xc"] for r_ in rows], float), axis=0)
        eta = float(np.nanmean([r_["arms"][arm].get("eta_mean", np.nan) for r_ in rows])) if arm != "uniform" else float("nan")
        stats[arm] = {"v_s": vs[:2], "v_c": vc[:2], "class": cls, "R": R, "n": nc, "F": F, "xc": xc, "eta": eta}
        tser = rows[0]["arms"][arm]["series"]["t"]
        lines.append(f"  {arm:8s}: spreading rate v_s = {vs[0]:+.4f} +- {vs[1]:.4f} -> {cls}; centroid velocity (source frame) "
                     f"v_c = {vc[0]:+.4f} +- {vc[1]:.4f}" + (f"; mean step rapidity {eta:.3f}" if np.isfinite(eta) else ""))
        lines.append("        t      " + " ".join(f"{x:6.1f}" for x in tser))
        lines.append("        R_rms  " + " ".join(f"{x:6.2f}" if np.isfinite(x) else "   nan" for x in R))
        lines.append("        x_c    " + " ".join(f"{x:+6.2f}" if np.isfinite(x) else "   nan" for x in xc))
        lines.append("        carry  " + " ".join(f"{x:6.1f}" for x in nc))
        lines.append("        flux   " + " ".join(f"{x:6.1f}" for x in F))
    inst = {"uniform arm spreads": stats["uniform"]["class"] == "spreading",
            "boosted recomputation chooses the same children": all(r_["check_identical"] for r_ in rows)}
    c = stats["carried"]
    if not all(inst.values()):
        status = "INSTRUMENT_FAILED"
    elif p["mode"] == "rest":
        status = {"localized": "PASSED", "spreading": "FAILED"}.get(c["class"], "UNDECIDED")
    else:
        v, se = c["v_c"]
        tol = beta / 4.0
        if c["class"] == "spreading" or (np.isfinite(v) and abs(v) - nse * se > tol):
            status = "FAILED"
        elif c["class"] == "localized" and np.isfinite(v) and abs(v) + nse * se <= tol:
            status = "PASSED"
        else:
            status = "UNDECIDED"
        lines.append(f"  inertia: carried centroid velocity in the source frame {v:+.4f} +- {se:.4f} vs tolerance beta/4 = {tol:.4f} "
                     f"(alternative: at rest in the lab, v_c = {-beta:+.2f})")
    lines.append("instrument: " + "; ".join(f"{k}: {'OK' if v else 'FAILS'}" for k, v in inst.items()))
    data = {"stats": {a: {k: (list(v) if isinstance(v, (tuple, np.ndarray)) else v) for k, v in s_.items()}
                      for a, s_ in stats.items()},
            "instrument": inst, "rows": [{"N": r_["N"], "links": r_["links"], "sources": r_["sources"],
                                          "arms": {a: {k: v for k, v in r_["arms"][a].items()} for a in ARMS}} for r_ in rows]}
    return Result(status, lines, data)
