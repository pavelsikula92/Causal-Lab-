"""
F3.6c -- Hot source: does inelastic merging with exact 4-momentum (F = |P|) hold together a bundle with internal motion?
=======================================================================================================================

Why this test exists (25 Oct 2026, author's design)
---------------------------------------------------
F3.6 showed that flux carrying its timelike direction as a state stays
localized and moves with its source. It left open what holds the width
of the bundle: it persists only because the flux is conserved and keeps
its direction. The derivation docs/derivations/spatial_part.md (23 Oct)
requires matter that conserves 4-momentum (integrability of the focusing
principle). The F3.6 carried rule does not: when streams merge it sends on
F u, with F the scalar sum of the incoming fluxes and u = P / |P|, and
|P| >= F (reverse triangle inequality for future timelike vectors), so the
excess |P| - F is thrown away. The corrected rule sends on the whole
4-momentum P; the scalar flux of the merged stream is its invariant mass
F = |P| (rest mass grows in an inelastic merger, as in special relativity).

In F3.6 itself the correction is empty: every source event emits along the
same u_src, all streams carry exactly u_src and |P| = F at every merger
(checked 25 Oct on one sprinkling: largest deviation of a carried direction
from u_src 0 at rest and 2e-16 moving). The rules differ only when streams
of different directions meet, so the source must be hot.

Model (no new mechanism: the F3.6 routing with three merge rules)
-----------------------------------------------------------------
Fixed causal set: Poisson sprinkling of a 2+1 box, rho = 1 (lengths in l),
links by core.order.links; signature (+ - -), so for a future timelike
vector the invariant mass is |P| = sqrt(P . P) = sqrt(P_0^2 - |P_vec|^2)
(the same as sqrt(-P^mu P_mu) in (- + + +)). Source: the world tube of
radius r_src at rest in the laboratory, time [0, t_src); every source event
i emits weight s_i along its own unit 4-velocity u_i (the source config,
below). Every event that holds a stream sends ALL of it to the one child
whose step has the smallest rapidity relative to the stream direction u
(cosh eta = u . b / tau_b, exactly F3.6). Three arms, identical initial
events, weights, directions and sprinkling:
  A  old       the F3.6 carried rule: F_v = sum of incoming F + s_v,
               u_v = P_v / |P_v| with P_v = sum of incoming F u + s_v u_v;
               the child receives F_v u_v (excess |P_v| - F_v discarded).
  B  |P|       the child receives P_v itself; F_v = |P_v|. 4-momentum is
               conserved exactly by construction; the scalar flux is not
               (the rest mass grows when directions differ).
  C  free      no merging: every source event's stream is followed alone
               along its own chain with its fixed direction, even through
               events that other streams also visit. The free-dispersal
               reference, not a physical alternative to merging.
In A and B a stream changes direction only at a merger; an unmerged stream
moves exactly as in C. A and B differ in a first-generation merger only in
the weight sent on, so their directions differ only from the second
generation on (B weights heavier, already merged streams more).

Source configs (directions in the source rest frame = laboratory)
------------------------------------------------------------------
  cold   s_i = 1, u_i = (1, 0, 0): the F3.6 case, control (A and B must coincide).
  sym    s_i = 1, rapidity eta, angles phi_0 + 2 pi k / n over a random
         permutation of the source events: the spatial momenta cancel exactly.
  asym   s_i uniform in [0.5, 1.5], rapidity uniform in [0.5 eta, 1.5 eta],
         angle uniform in (-pi/2, pi/2) with probability 1/2, otherwise uniform
         in (-pi, pi]: net momentum, unequal weights, no special symmetry.
Parts: eta = 0.3 and 0.6 (each with all three configs, on the same sprinklings).

Readout (per sprinkling, arm, config)
-------------------------------------
The CM frame of the total source 4-momentum P_src (velocity P_vec / P_0), its
4-velocity U = P_src / |P_src|. An entry is a stream at an event (A, B: every
event holding flux; C: every stream at every event of its chain), with the
4-momentum it sends on. Weight = its CM energy P . U (a Lorentz scalar).
Per slab of CM time (width 2, from t_src up to where the CM world line leaves the box):
  - R_rms of the entries about their weighted centroid (not the source line);
  - internal directional spread: weighted mean of acosh(u . U), the rapidity
    of the stream relative to the common motion;
  - number of entries (active streams), their total weight.
Spreading rate v_s = least-squares slope of R_rms over CM time [6, 18].
Mergers (A, B): events receiving two or more streams (a source emission
counts as one); at each, the equivalent relative rapidity
eta_eq = 2 acosh(|P| / sum F_in) (exact for two equal streams), weighted by |P|.
Balances: Delta P = sum_sources s u - sum_ends P (ends = entries with no
child, the flux leaving the box), every component, relative to the CM
energy of the source. Rest-mass gain sum_ends F / sum s.

Criterion (stated before any full or quick run)
-----------------------------------------------
INSTRUMENT_FAILED if any of:
  (i)   B or C: |Delta P^mu| / E_src > 1e-12 for some component, sprinkling or config;
  (ii)  the recomputation in the frame boosted by beta' = 0.6 chooses a
        different child anywhere, in any arm and config;
  (iii) cold: A and B differ in any choice, or their fluxes by more than
        1e-12 relative (the correction must be empty there);
  (iv)  a hot config is not discriminating: the weighted mean eta_eq at the
        mergers of B (mean over sprinklings) is below eta / 4 (the streams
        would meet nearly parallel, the case already known to be empty);
  (v)   C does not spread in a hot config: v_s(C) - 2 SE < v_free / 2, with
        v_free the CM-energy-weighted RMS velocity of the source streams in the
        CM frame (the window could not tell free dispersal from binding).
Per hot config, against free dispersal (quarter rule of F3.6, scaled to C):
  bound     v_s(B) + 2 SE <= v_s(C) / 4;
  spreading v_s(B) - 2 SE >  v_s(C) / 4;
  inconclusive otherwise.
Part (one eta): PASSED if B is bound in sym AND asym; FAILED if B spreads in
both; UNDECIDED otherwise.
Reported, not decided: the same classification for A; the paired difference
v_s(B) - v_s(A) (bootstrap over sprinklings; "no measurable difference" if
|Delta| <= 2 SE); the time series of R_rms and of the internal spread for all
arms; the internal spread late (last three slabs) for A, B, C; the number of
mergers and of active streams; the rest-mass gain of B; the energy deficit of
A (what the old rule loses: a stability of A could be partly that loss).

Reading declared now (the author's decision rules, 25 Oct):
  B breaks the balance            -> fix the implementation, read no physics (INSTRUMENT_FAILED (i));
  B conserves, the bundle widens  -> the merge rule keeps the books, it does not hold the width (FAILED);
  B narrows the internal spread and holds the width -> interesting, needs checks that the stability
                                     is not an artificial loss of internal energy (excluded by (i): in B
                                     the relative kinetic energy becomes rest mass, measured as the gain)
                                     or another numerical artefact (PASSED, followed by replication);
  B does not differ from A        -> the new rule has no measurable effect in this regime (reported).
A localized B is not yet a self-bound object of the theory: the routing
still uses the embedding at the scale of one link and the future of each
event (audit criteria (a), (b), as in the whole F3 series).

Quick runs (after the criteria were written; decide nothing): mechanics in order. B balances to 1e-16, the cold
control is identical bit for bit, the boosted recomputation chooses identically. The old rule A loses 2-3 % of the
energy at eta = 0.3 and 7 % at 0.6; B gains the same as rest mass (gain 1.03 and 1.08). One full-size sprinkling
(eta = 0.6, ~61 000 events, 222 source events, ~140 mergers, 106 s on one core): v_s of A / B / C 0.29 / 0.27 / 0.50
(sym) and 0.23 / 0.22 / 0.42 (asym); internal spread late 0.25 / 0.23 / 0.60 (sym). No criterion changed.
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
    "id": "F3.6c",
    "phase": "F3",
    "title": "Hot source: does inelastic merging with exact 4-momentum (F = |P|) hold together a bundle with "
             "internal motion?",
    "theory_part": "matter and momentum: conservation of the 4-momentum of matter (a condition of the focusing "
                   "principle, docs/derivations/spatial_part.md); what holds the width of a bundle (open after F3.6)",
    "question": "A source emits streams of different directions. Does a merging that passes on the whole 4-momentum "
                "(flow = invariant mass |P|) hold the bundle together against free dispersal, and does it differ "
                "from the old rule, which discards the excess |P| − F?",
    "method": "2+1 sprinkling (ρ = 1), links, F3.6 routing (the descendant of the smallest rapidity relative to the "
              "direction of the stream); branches A old rule, B F = |P|, C free motion without merging, with the "
              "same events, weights and directions; source configurations cold (control), sym (spatial momenta "
              "cancelled exactly), asym (net momentum, unequal weights); parts η = 0.3 and 0.6; readout in the "
              "centre-of-momentum frame: R_rms about the centroid, internal spread of directions, merging, "
              "4-momentum balance",
    "prediction": "not fixed in advance; B differs from A only from the second generation of merging on (weights |P| "
                  "against F)",
    "criterion": "INSTRUMENT: balance of B or C > 1e-12, the recomputation in the frame β' = 0.6 chooses "
                 "differently, cold A ≠ B, merging in B with a mean relative rapidity < η/4, or C does not disperse "
                 "(v_s(C) − 2 SE < v_free/2). Configuration: B bound v_s(B) + 2 SE <= v_s(C)/4, dispersing v_s(B) − "
                 "2 SE > v_s(C)/4. Part: PASSED bound in sym and asym, FAILED dispersing in both, otherwise "
                 "UNDECIDED.",
    "null_model": "free motion of the streams without merging (C); the old merging rule without conservation of "
                  "4-momentum (A)",
    "level": "mechanism",
    "meters": ["K.6"],
    "depends_on": ["F3.6"],
    "params": {"D": 3, "T": 32.0, "half": 20.0, "x_extra": 8.0, "r_src": 4.0, "t_src": 4.0, "slab": 2.0,
               "win": [6.0, 18.0], "beta_check": 0.6, "n_se": 2.0, "boot_B": 4000, "n_seeds": 32,
               "tol_bal": 1e-12, "configs": ["cold", "sym", "asym"]},
    "parts": {"eta03": {"eta": 0.3}, "eta06": {"eta": 0.6}},
    "quick": {"n_seeds": 3, "boot_B": 200, "T": 16.0, "half": 10.0, "x_extra": 4.0, "r_src": 2.0, "t_src": 2.0,
              "slab": 1.0, "win": [3.0, 9.0]},
    "estimate": "per part ~6–10 min with 10 processes (one sprinkling of ~61 000 events ~2 min per core, the links "
                "dominate); both parts use the same sprinklings",
}

ARMS = ("old", "absP", "free")
ARM_CS = {"old": "A old", "absP": "B |P|", "free": "C free"}


# ============================================================================================== geometry

def mink(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Minkowski product a . b = a_0 b_0 - a_vec . b_vec (signature + - -), broadcast over leading axes."""
    return a[..., 0] * b[..., 0] - (a[..., 1:] * b[..., 1:]).sum(axis=-1)


def boost(X: np.ndarray, v) -> np.ndarray:
    """Coordinates (or 4-vectors) in the frame moving with spatial velocity v (2 components, |v| < 1, c = 1)."""
    v = np.asarray(v, float)
    v2 = float(v @ v)
    if v2 == 0.0:
        return np.array(X, float, copy=True)
    g = 1.0 / math.sqrt(1.0 - v2)
    X = np.asarray(X, float)
    vx = X[..., 1:] @ v                                            # v . x
    Y = np.empty_like(X)
    Y[..., 0] = g * (X[..., 0] - vx)
    Y[..., 1:] = X[..., 1:] + (((g - 1.0) * vx / v2) - g * X[..., 0])[..., None] * v
    return Y


def children_csr(ptr: np.ndarray, idx: np.ndarray, N: int):
    """Children lists from the parents CSR of core.order.links: (cptr, cidx)."""
    child = np.repeat(np.arange(N), np.diff(ptr))
    o = np.argsort(idx, kind="stable")
    cptr = np.concatenate([[0], np.cumsum(np.bincount(idx, minlength=N))]).astype(np.int64)
    return cptr, child[o]


def best_child(X, cptr, cidx, v, u):
    """The child of v whose step has the smallest rapidity relative to the direction u (-1 if v has no child)."""
    a, b = cptr[v], cptr[v + 1]
    if b == a:
        return -1
    ks = cidx[a:b]
    B = X[ks] - X[v]
    tau = np.sqrt(np.maximum(mink(B, B), 1e-300))
    return int(ks[int(np.argmin(mink(u[None, :], B) / tau))])


# ============================================================================================== source configs

def source_config(cfg: str, eta: float, n: int, rng: np.random.Generator):
    """Weights s and unit 4-velocities U (n x 3) of the n source events (module notes)."""
    if cfg == "cold":
        return np.ones(n), np.tile([1.0, 0.0, 0.0], (n, 1))
    if cfg == "sym":
        s = np.ones(n)
        phi = rng.uniform(0.0, 2.0 * math.pi) + 2.0 * math.pi * rng.permutation(n) / n
        r = np.full(n, eta)
    elif cfg == "asym":
        s = rng.uniform(0.5, 1.5, n)
        r = rng.uniform(0.5 * eta, 1.5 * eta, n)
        fwd = rng.random(n) < 0.5
        phi = np.where(fwd, rng.uniform(-0.5 * math.pi, 0.5 * math.pi, n), rng.uniform(-math.pi, math.pi, n))
    else:
        raise ValueError(cfg)
    U = np.stack([np.cosh(r), np.sinh(r) * np.cos(phi), np.sinh(r) * np.sin(phi)], axis=1)
    return s, U


# ============================================================================================== rules

def route_merge(X, cptr, cidx, s, U0, mode: str):
    """
    Arms A (mode 'old') and B (mode 'absP') in one causal pass over the time-sorted events.
    Returns per event: the 4-momentum it sends on (Pout), its scalar flux F, the chosen child, the number of incoming
    streams (a source emission counts as one), and for mergers the equivalent relative rapidity and |P|.
    """
    N, D = X.shape
    Pin = np.zeros((N, D))                       # 4-momentum arriving (A: sum F u; B: sum P)
    Fin = np.zeros(N)                            # scalar flux arriving (A: sum F; B: sum |P|)
    nin = np.zeros(N, dtype=np.int64)
    Pout = np.zeros((N, D))
    F = np.zeros(N)
    choice = np.full(N, -1, dtype=np.int64)
    eta_eq = np.full(N, np.nan)
    for v in range(N):
        if nin[v] == 0 and s[v] <= 0.0:
            continue
        P = Pin[v] + s[v] * U0[v]
        Fsum = Fin[v] + s[v]
        M = math.sqrt(max(float(mink(P, P)), 0.0))
        u = P / M
        n_streams = nin[v] + (1 if s[v] > 0.0 else 0)
        if n_streams >= 2:
            eta_eq[v] = 2.0 * math.acosh(max(M / Fsum, 1.0))
        if mode == "old":
            F[v] = Fsum
            Pout[v] = Fsum * u                    # the excess |P| - F is discarded
        else:
            F[v] = M
            Pout[v] = P                           # the whole 4-momentum is sent on
        k = best_child(X, cptr, cidx, v, u)
        choice[v] = k
        if k >= 0:
            Pin[k] += Pout[v]
            Fin[k] += F[v]
            nin[k] += 1
    return {"Pout": Pout, "F": F, "choice": choice, "nin": nin, "eta_eq": eta_eq}


def route_free(X, cptr, cidx, s, U0, src):
    """Arm C: every source event's stream alone along its chain. Returns entries (event, 4-momentum, end flag) and chains."""
    ev, P, end, chains = [], [], [], []
    for i in src:
        p = s[i] * U0[i]
        u = U0[i]
        v = int(i)
        chain = []
        while True:
            chain.append(v)
            k = best_child(X, cptr, cidx, v, u)
            ev.append(v)
            P.append(p)
            end.append(k < 0)
            if k < 0:
                break
            v = k
        chains.append(chain)
    return {"ev": np.array(ev, dtype=np.int64), "P": np.array(P, float).reshape(-1, X.shape[1]),
            "end": np.array(end, bool), "chains": chains}


# ============================================================================================== readout

def slab_readout(Xc, P, Ucm, p, beta_cm):
    """Per slab of CM time: R_rms about the weighted centroid, internal spread, entries, total weight."""
    T, sl = float(p["T"]), float(p["slab"])
    t_end = T * math.sqrt(1.0 - beta_cm * beta_cm) - sl
    edges = np.arange(float(p["t_src"]), t_end + 1e-9, sl)
    w_all = mink(P, Ucm[None, :])                                  # CM energy of each entry
    m_all = np.sqrt(np.maximum(mink(P, P), 1e-300))
    eta_all = np.arccosh(np.maximum(w_all / m_all, 1.0))           # rapidity relative to the common motion
    out = {k: [] for k in ("t", "R", "spread", "n", "W")}
    for a, b in zip(edges[:-1], edges[1:]):
        m = (Xc[:, 0] >= a) & (Xc[:, 0] < b)
        w = w_all[m]
        out["t"].append(0.5 * (a + b))
        out["n"].append(int(m.sum()))
        out["W"].append(float(w.sum()))
        if w.sum() > 0:
            c = (w[:, None] * Xc[m, 1:]).sum(axis=0) / w.sum()
            r2 = ((Xc[m, 1:] - c) ** 2).sum(axis=1)
            out["R"].append(float(np.sqrt((w * r2).sum() / w.sum())))
            out["spread"].append(float((w * eta_all[m]).sum() / w.sum()))
        else:
            out["R"].append(float("nan"))
            out["spread"].append(float("nan"))
    return out


def slope(t, y, win):
    """Least-squares slope of y(t) over the window (finite points only)."""
    t, y = np.asarray(t, float), np.asarray(y, float)
    m = (t >= win[0]) & (t <= win[1]) & np.isfinite(y)
    if m.sum() < 2:
        return float("nan")
    return float(np.polyfit(t[m], y[m], 1)[0])


# ============================================================================================== task

def seed_task(p, seed, cfg_seeds):
    """One sprinkling: for every source config the three arms, the readout, balances and the boosted recomputation."""
    D, T, H = int(p["D"]), float(p["T"]), float(p["half"])
    eta = float(p["eta"])
    rng = np.random.default_rng(seed)
    X = sprinkle_poisson(Box(D, [0.0, -H, -H], [T, H + float(p["x_extra"]), H]), 1.0, rng)
    N = len(X)
    ptr, idx = links_csr(X, 1.0)
    cptr, cidx = children_csr(ptr, idx, N)
    src = np.flatnonzero((X[:, 0] < float(p["t_src"])) & ((X[:, 1:] ** 2).sum(axis=1) <= float(p["r_src"]) ** 2))
    bc = float(p["beta_check"])
    Xb = boost(X, [bc, 0.0])
    out = {"N": N, "links": int(idx.size), "sources": int(src.size), "configs": {}}
    for cfg in p["configs"]:
        s_src, U_src = source_config(cfg, eta, src.size, np.random.default_rng(cfg_seeds[cfg]))
        s = np.zeros(N)
        U0 = np.zeros((N, D))
        s[src], U0[src] = s_src, U_src
        P_src = (s_src[:, None] * U_src).sum(axis=0)
        M_src = math.sqrt(float(mink(P_src, P_src)))
        Ucm = P_src / M_src
        v_cm = P_src[1:] / P_src[0]
        beta_cm = float(np.hypot(*v_cm))
        Xc = boost(X, v_cm)
        U_src_cm = boost(U_src, v_cm)                              # free-stream velocities in the CM frame
        vel2 = (U_src_cm[:, 1:] ** 2).sum(axis=1) / U_src_cm[:, 0] ** 2
        wE = s_src * U_src_cm[:, 0]
        v_free = float(math.sqrt((wE * vel2).sum() / wE.sum()))
        U0b = boost(U0, [bc, 0.0])
        rec = {"beta_cm": beta_cm, "v_free": v_free, "E_src": M_src, "sum_s": float(s_src.sum()), "arms": {},
               "check_identical": True}
        res = {}
        for mode in ("old", "absP"):
            r = route_merge(X, cptr, cidx, s, U0, mode)
            rb = route_merge(Xb, cptr, cidx, s, U0b, mode)
            rec["check_identical"] &= bool(np.array_equal(r["choice"], rb["choice"]))
            res[mode] = r
        fr = route_free(X, cptr, cidx, s, U0, src)
        frb = route_free(Xb, cptr, cidx, s, U0b, src)
        rec["check_identical"] &= (fr["chains"] == frb["chains"])
        if cfg == "cold":
            a, b = res["old"], res["absP"]
            rec["cold_same_choice"] = bool(np.array_equal(a["choice"], b["choice"]))
            mF = a["F"] > 0
            rec["cold_flux_rel"] = float(np.max(np.abs(a["F"][mF] - b["F"][mF]) / a["F"][mF])) if mF.any() else 0.0
        for arm in ARMS:
            if arm == "free":
                ev, P, end = fr["ev"], fr["P"], fr["end"]
            else:
                r = res[arm]
                ev = np.flatnonzero(r["F"] > 0)
                P = r["Pout"][ev]
                end = r["choice"][ev] < 0
            ser = slab_readout(Xc[ev], P, Ucm, p, beta_cm)
            dP = P_src - P[end].sum(axis=0)
            a = {"series": ser, "v_s": slope(ser["t"], ser["R"], p["win"]),
                 "spread_late": float(np.nanmean(ser["spread"][-3:])),
                 "bal_rel": float(np.max(np.abs(dP)) / M_src),
                 "E_deficit": float(mink(dP, Ucm) / M_src),
                 "mass_gain": float(np.sqrt(np.maximum(mink(P[end], P[end]), 0.0)).sum() / s_src.sum())}
            if arm != "free":
                r = res[arm]
                mm = np.isfinite(r["eta_eq"])
                a["n_merge"] = int(mm.sum())
                a["eta_merge"] = float((r["eta_eq"][mm] * r["F"][mm]).sum() / r["F"][mm].sum()) if mm.any() else 0.0
            rec["arms"][arm] = a
        out["configs"][cfg] = rec
    return out


# ============================================================================================== run

def classify(v, se, ref, nse):
    """bound / spreading / inconclusive against a quarter of the free spreading rate."""
    if not all(np.isfinite(x) for x in (v, se, ref)):
        return "undefined"
    if v + nse * se <= ref / 4.0:
        return "bound"
    if v - nse * se > ref / 4.0:
        return "spreading"
    return "inconclusive"


def run(p, ctx):
    t0 = time.time()
    n, B, nse, eta = int(p["n_seeds"]), int(p["boot_B"]), float(p["n_se"]), float(p["eta"])
    tol_bal = float(p["tol_bal"])
    seeds = [derive_seed(ctx.seed, "F3.6c", k) for k in range(n)]           # the same sprinklings in every part
    jobs = [(p, sd, {c: derive_seed(ctx.seed, "F3.6c-src", ctx.part, c, k) for c in p["configs"]})
            for k, sd in enumerate(seeds)]
    rows = parallel_map(seed_task, jobs, ctx.workers, ctx.progress,
                        lambda i, r_, d, nn: f"[F3.6c] {ctx.part} sprinkling {d}/{nn}: v_s sym A/B/C "
                                             + "/".join(f"{r_['configs']['sym']['arms'][a]['v_s']:.3f}" for a in ARMS)
                                             + f" ({time.time()-t0:.0f}s)")
    lines = [f"2+1 sprinkling rho = 1 (lengths in l), box t [0, {p['T']:g}], x [-{p['half']:g}, {p['half'] + p['x_extra']:g}], "
             f"y [-{p['half']:g}, {p['half']:g}]; source tube radius {p['r_src']:g}, time [0, {p['t_src']:g}) at rest in the lab; "
             f"eta = {eta:g}; {n} sprinklings (~{np.mean([r_['N'] for r_ in rows]):.0f} events, "
             f"{np.mean([r_['sources'] for r_ in rows]):.0f} source events); readout in the CM frame, window {list(p['win'])}"]
    inst = {"balance B and C <= %g" % tol_bal: True, "boosted recomputation identical": True,
            "cold: A = B": True}
    verdicts, summary = {}, {}
    for cfg in p["configs"]:
        R = [r_["configs"][cfg] for r_ in rows]
        bs = lambda *tag, cfg=cfg: derive_seed(ctx.seed, "F3.6c-boot", ctx.part, cfg, *tag)
        inst["boosted recomputation identical"] &= all(x["check_identical"] for x in R)
        for arm in ("absP", "free"):
            inst["balance B and C <= %g" % tol_bal] &= all(x["arms"][arm]["bal_rel"] <= tol_bal for x in R)
        st = {}
        for arm in ARMS:
            vs = boot_se([x["arms"][arm]["v_s"] for x in R], B, bs(arm, "vs"))
            sp = boot_se([x["arms"][arm]["spread_late"] for x in R], B, bs(arm, "sp"))
            with np.errstate(all="ignore"):
                Rm = np.nanmean(np.array([x["arms"][arm]["series"]["R"] for x in R], float), axis=0)
                Sm = np.nanmean(np.array([x["arms"][arm]["series"]["spread"] for x in R], float), axis=0)
                Nm = np.nanmean(np.array([x["arms"][arm]["series"]["n"] for x in R], float), axis=0)
            st[arm] = {"v_s": vs[:2], "spread_late": sp[:2], "R": Rm, "spread": Sm, "n": Nm,
                       "bal": max(x["arms"][arm]["bal_rel"] for x in R),
                       "E_def": float(np.mean([x["arms"][arm]["E_deficit"] for x in R])),
                       "gain": float(np.mean([x["arms"][arm]["mass_gain"] for x in R]))}
            if arm != "free":
                st[arm]["n_merge"] = float(np.mean([x["arms"][arm]["n_merge"] for x in R]))
                st[arm]["eta_merge"] = float(np.mean([x["arms"][arm]["eta_merge"] for x in R]))
        dBA = boot_se([x["arms"]["absP"]["v_s"] - x["arms"]["old"]["v_s"] for x in R], B, bs("dBA"))
        v_free = float(np.mean([x["v_free"] for x in R]))
        vC = st["free"]["v_s"]
        lines.append(f"  config {cfg}: CM velocity {np.mean([x['beta_cm'] for x in R]):.3f}; free-stream RMS velocity "
                     f"v_free {v_free:.3f}")
        tser = R[0]["arms"]["free"]["series"]["t"]
        for arm in ARMS:
            a = st[arm]
            cls = classify(a["v_s"][0], a["v_s"][1], vC[0], nse) if arm != "free" else "reference"
            a["class"] = cls
            extra = (f"; mergers {a['n_merge']:.0f}, mean eta_eq {a['eta_merge']:.3f}" if arm != "free" else "")
            lines.append(f"    {ARM_CS[arm]:9s}: v_s = {a['v_s'][0]:+.4f} +- {a['v_s'][1]:.4f} ({cls}); internal spread late "
                         f"{a['spread_late'][0]:.4f} +- {a['spread_late'][1]:.4f}; max |dP|/E {a['bal']:.1e}; energy deficit "
                         f"{a['E_def']:+.2e}; rest-mass gain {a['gain']:.4f}{extra}")
            lines.append("        t       " + " ".join(f"{x:6.1f}" for x in tser))
            lines.append("        R_rms   " + " ".join(f"{x:6.2f}" if np.isfinite(x) else "   nan" for x in a["R"]))
            lines.append("        spread  " + " ".join(f"{x:6.3f}" if np.isfinite(x) else "   nan" for x in a["spread"]))
            lines.append("        streams " + " ".join(f"{x:6.1f}" for x in a["n"]))
        diff = "no measurable difference" if abs(dBA[0]) <= nse * dBA[1] else "differ"
        lines.append(f"    v_s(B) - v_s(A) = {dBA[0]:+.5f} +- {dBA[1]:.5f} ({diff})")
        if cfg == "cold":
            same = all(x["cold_same_choice"] for x in R)
            frel = max(x["cold_flux_rel"] for x in R)
            inst["cold: A = B"] = same and frel <= tol_bal
            lines.append(f"    cold control: A and B choose identically: {'yes' if same else 'NO'}; max flux difference {frel:.1e}")
        else:
            disc = st["absP"]["eta_merge"] >= eta / 4.0
            spreads = np.isfinite(vC[0]) and vC[0] - nse * vC[1] >= v_free / 2.0
            inst[f"{cfg}: mergers discriminate (mean eta_eq >= {eta / 4:g})"] = disc
            inst[f"{cfg}: C spreads (v_s - 2 SE >= v_free / 2 = {v_free / 2:.3f})"] = bool(spreads)
            verdicts[cfg] = st["absP"]["class"]
        summary[cfg] = {"v_free": v_free, "dBA": dBA[:2],
                        "arms": {a: {k: (list(v) if isinstance(v, (tuple, np.ndarray)) else v) for k, v in st[a].items()}
                                 for a in ARMS}}
    hot = [c for c in p["configs"] if c != "cold"]
    if not all(inst.values()):
        status = "INSTRUMENT_FAILED"
    elif all(verdicts[c] == "bound" for c in hot):
        status = "PASSED"
    elif all(verdicts[c] == "spreading" for c in hot):
        status = "FAILED"
    else:
        status = "UNDECIDED"
    lines.append("verdicts of B against free dispersal: " + "; ".join(f"{c}: {verdicts[c]}" for c in hot))
    lines.append("instrument: " + "; ".join(f"{k}: {'OK' if v else 'FAILS'}" for k, v in inst.items()))
    data = {"summary": summary, "instrument": inst, "verdicts": verdicts,
            "rows": [{"N": r_["N"], "links": r_["links"], "sources": r_["sources"],
                      "configs": {c: {k: v for k, v in x.items() if k != "arms"} | {"arms": x["arms"]}
                                  for c, x in r_["configs"].items()}} for r_ in rows]}
    return Result(status, lines, data)
