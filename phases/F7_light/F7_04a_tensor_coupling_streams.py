"""
F7.4a -- Tensor coupling from the response of the order: do two counter-streaming flows gravitate as T_ab?
==========================================================================================================

Why this test exists (28 Oct 2026, author's decision 13: route A)
-----------------------------------------------------------------
In 3+1 the throughput principle (F7.1c: gamma ~ 1) is, in the linear order, the same thing as coupling to the full
energy-momentum tensor: the weak field of a dust source is hbar_ab = 4U u_a u_b (docs/derivations/tensor_coupling.md).
Route A takes the principle as the law of the plateau and general relativity as its nonlinear completion. What is
still open is how the SUBSTRATE combines the saturation of flows that move differently. Two equally strong streams
with velocities +v and -v (Gamma_v = 1/sqrt(1 - v^2), not to be confused with the PPN gamma) have zero momentum
density but nonzero energy and spatial stress, so the candidate rules give different clock fields.

F7.4 (23 Oct) was set aside because its tensor part only verified inserted algebra (the trace reversal of the
linearised Einstein equation and a Green function). Nothing of that kind enters here: no field equation, no Green
function, no trace reversal. The tensor structure, if it appears, comes from the response rule of the order that
F7.1c established for one flow, applied per flow in its own rest frame.

The substrate rule under test (input, from the ontology)
--------------------------------------------------------
  per-event, own frame  every saturated flow narrows the capacity cone in ITS OWN rest frame (F7.1c rule, Lorentz
                        covariant; a flow carries its 4-velocity as a state, F3.6);
  additive capacity     flows that load the same link share its finite capacity: their cone deficits add (series
                        impedance, book 5.3);
  throughput            the event density rises by exp(sigma_tot) (events are Lorentz scalars).
For a link d from a parent event with stream saturations sigma_i (rest frame) and 4-velocities u_i (lower index,
u_a d^a = Gamma_i(-d_t + v_i d_x)):
    kept  iff  d_t > 0  and  -d_t^2 + |d_x|^2 + sum_i a_i (u_i . d)^2 <= 0,
    a_i = (sigma_i / sigma_tot) (1 - exp(-2 sigma_tot)),
so that for streams at rest this is exactly the F7.1c rule with sigma_tot (cones exp(-sigma_tot)). The relation is
reachability through kept links; virtual diamond endpoints relate directly to the window within the cone of the
endpoint (as F7.1a-c). No metric is sprinkled: the order is built from a flat sprinkling and per-link judgements.

Normalisation (frozen before any run; the author's point of 28 Oct)
--------------------------------------------------------------------
Fixed when v changes: the REST-FRAME saturation of each stream, sigma_i = sigma0 exp(-rho^2 / 2 w^2) with
rho = sqrt(y^2 + z^2) (the streams fill a cylinder along x and are uniform along x), i.e. the proper density of
each stream; the transverse profile in the lab; the box. A boost along x keeps rho and a uniform-in-x profile, so
this normalisation is exact. With it the linear predictions for the clock field, relative to v = 0, are
    tensor (this rule)                          Gamma_v^2 (1 + v^2)     2.125 at v = 0.6
    energy density in the mean-flow frame       Gamma_v^2               1.5625
    lab-frame density                           Gamma_v                 1.25
    rest density (event count)                  1                       1.000
(The energy-density reading is what the F3 merge rule F = |P| produces when counter-streams merge into rest dust;
it is listed because it differs from the tensor only by the stress term.) The three alternatives are "mean-flow
scalar" rules: isotropic cones in the local mean-flow frame narrowed by f(v) sigma_tot, density exp(f sigma_tot);
their clock field is exactly f(v) times the v = 0 field (the F7.1c rule is exact for the clocks: phi = -sigma/2).

Why the amplitude of the near field and not the 1/r tail: for a compact, closed, stationary system the virial
(von Laue) theorem gives int T_ij d^3x = 0, so the 1/r amplitude is the total energy for the tensor AND for the
energy-density rule; the stress term shows only along a direction where the system is not closed. Hence uniform
streams along x and the transverse profile of the clock field.

Observables and predictions (continuum, exact for this rule; computed by predict() below before any run)
--------------------------------------------------------------------------------------------------------
Configurations (16 sprinklings each, 3+1 box t in [0, 9], |x|, |y|, |z| <= 9, w = 4):
  flat        no streams (the null);
  two_v0      two streams at rest, sigma0 = 0.1 each (= F7.1c with sigma_tot = 0.2 in a cylinder);
  two_v06     two streams, v = +0.6 and -0.6, sigma0 = 0.1 each;
  single_v06  one stream, v = 0.6, sigma0 = 0.2 (reported: the boosted F7.1c order and the dragging of the cones).
(1) PRIMARY, amplitude: phi_clock per transverse shell (diamonds at rest in the lab, height 3). A = least-squares
    amplitude against the predicted v = 0 profile; R = A(two_v06) / A(two_v0). The exact prediction of the rule
    R_T is computed by projecting the predicted profiles in the same way (nonlinear: above 2.125 at sigma0 = 0.1).
(2) SECOND, anisotropy (normalisation free): the cone deficit along the streams against across them,
    K = sum_shells s_par / sum_shells s_perp (inner three shells); tensor: K_T > 1 (1 + v^2 = 1.36 in the linear
    order), every mean-flow scalar rule: K = 1 exactly (isotropic cones).
(3) reported: single_v06 cones forward (+x) and backward (-x): the tensor rule tilts the cone (g_0x, dragging of the
    cones by the flow; linear ratio of the deficits (1 + v)^2 / (1 - v)^2 = 16 at v = 0.6); its amplitude R.

Readers (blind: events and relation)
------------------------------------
Cylindrical shells rho in [0, 1.5, 3, 4.5, 6]; anchors and diamonds with |x| <= half - 3 and 1 <= t <= T - 1.
  density  ln(events per coordinate volume);
  cones    300 anchors per shell, pairs with dt in [1, 2], |dx| <= 1.5 dt, sorted by direction n = dx/|dx|:
           forward n_x >= cos 35 deg, backward n_x <= -cos 35 deg, across |n_x| <= sin 25 deg, along = forward and
           backward; pooled boundary (K.16b) per class; s = ln c;
  clocks   200 diamonds at rest per shell (random x, angle and start), height 3; phi = ln(tau / 3),
           tau = (24 <N> / pi)^{1/4} (K.16c).
Predictions per class: the solid-angle average of ln c(n) over the class (c(n) the null speed of the rule's cone).

Criteria (stated before any run)
--------------------------------
INSTRUMENT_FAILED if (i) flat: some reading of some shell |mean| > 0.01 + 3 SE; or (ii) static amplitude:
|A(two_v0) / A_pred(two_v0) - 1| - 2 SE > 0.25; or (iii) the directional reader sees anisotropy in isotropic
cones: |K(two_v0) - 1| - 2 SE > 0.1.
Amplitude:   PASSED |R - R_T| + 2 SE <= 0.25;  FAILED |R - R_T| - 2 SE > 0.25;  UNDECIDED otherwise.
Anisotropy:  tol_K = (K_T - 1) / 2;  PASSED |K - K_T| + 2 SE <= tol_K;  FAILED |K - K_T| - 2 SE > tol_K;  else UNDECIDED.
Overall: PASSED if both PASSED; FAILED if either FAILED; UNDECIDED otherwise. Reported: R against every candidate
(the nearest named), K, the profiles against the predictions, single_v06.

What is derived, what is tested, what is input
----------------------------------------------
  input:   the per-event own-frame rule and additive capacity sharing (ontology, not derived from the axioms);
  derived: its continuum clock field Gamma_v^2 (1 + v^2) in the linear order, i.e. coupling to T_00 + T_ii, and the
           anisotropy K_T (analytic, docs/derivations/tensor_coupling.md);
  tested:  that the discrete order built from a flat sprinkling by these per-link judgements, read blind, realises
           the tensor field and is told apart from the mean-flow scalar rules.
The test cannot tell which rule NATURE uses (that pressure and motion gravitate is known observationally); it tells
whether the substrate rule that route A needs works on the causal set.

Quick run (2 sprinklings, box |x| <= 6; decides nothing): the first attempt stopped on an operator-precedence error in
the diamond preselection (fixed, no result had been produced); the second is mechanically in order: links kept 28 %
(two_v0) and 16 % (two_v06), the static clock profile follows the prediction within ~10 %, the directional reader
reads isotropic cones as isotropic; the flat check fails with 2 sprinklings (SE from two values, as in K.16b). No
criterion, tolerance or parameter changed after it.
"""

from __future__ import annotations

import math
import time

import numpy as np

from core.order.links import links_csr
from core.order.sprinkle import derive_seed
from core.stats.parallel import parallel_map
from core.testing import Result
from phases.F7_light.F7_01a_capacity_cones import CapacityOrder
from phases.K_instruments import K_16c_gamma_reader_3d as k16c

TEST = {
    "id": "F7.4a",
    "phase": "F7",
    "title": "Tensor coupling from the response of the order: do counter-streaming flows gravitate according to T_ab?",
    "theory_part": "the source of gravity (decision 13, path A): how the substrate composes the saturation of flows "
                   "with different 4-velocities",
    "question": "When each flow narrows the cone in its own rest frame and the capacity deficits add, does the blind "
                "instrument read for two counter-streaming flows (±0.6c) a clock field Γ_v²(1 + v²) relative to rest "
                "(tensor) and anisotropic cones, and not Γ_v², Γ_v or 1 (scalar rules of the mean frame)?",
    "method": "flat 3+1 sprinkling, density exp(σ_tot), link kept when −dt² + |dx|² + Σ a_i (u_i·d)² <= 0 at the "
              "parent, a_i = σ_i/σ_tot (1 − e^{−2σ_tot}); streams along x, σ_i = σ0·exp(−ρ²/32) in the rest frame of "
              "the stream (fixed proper density); relation = reachability; readings: clocks (diamonds at rest) and "
              "cones by direction in radial shells around the x axis; configurations flat, two_v0, two_v06, "
              "single_v06; 16 sprinklings",
    "prediction": "amplitude R = A(v = 0.6)/A(0) = Γ_v²(1 + v²) (linearly 2.125, exactly for σ0 = 0.1 higher, see "
                  "predict); anisotropy K = deficit along / across > 1 (linearly 1 + v²); alternatives: Γ_v² 1.5625, "
                  "Γ_v 1.25, 1",
    "criterion": "INSTRUMENT: the flat configuration reads non-zero, static amplitude |A/A_pred − 1| − 2 SE > 0.25, "
                 "or anisotropy in isotropic cones |K − 1| − 2 SE > 0.1. Amplitude PASSED |R − R_T| + 2 SE <= 0.25, "
                 "FAILED |R − R_T| − 2 SE > 0.25; anisotropy with the tolerance (K_T − 1)/2; overall PASSED both, "
                 "FAILED either, otherwise UNDECIDED.",
    "null_model": "no streams (uniform events, all links)",
    "level": "mechanism",
    "meters": ["K.16c"],
    "depends_on": ["F7.1c", "K.16c"],
    "params": {"T": 9.0, "half": 9.0, "width": 4.0, "bins": [0, 1.5, 3.0, 4.5, 6.0], "x_margin": 3.0,
               "configs": [["flat", []], ["two_v0", [[0.1, 0.0], [0.1, 0.0]]],
                           ["two_v06", [[0.1, 0.6], [0.1, -0.6]]], ["single_v06", [[0.2, 0.6]]]],
               "n_seeds": 16, "anchors": 300, "diamonds": 200, "T_d": 3.0, "cone_dt": [1.0, 2.0], "cone_q": 1.5,
               "cos_par": math.cos(math.radians(35.0)), "sin_perp": math.sin(math.radians(25.0)), "k_shells": 3,
               "tol_R": 0.25, "tol_static": 0.25, "tol_iso": 0.1, "null_abs": 0.01, "null_se": 3.0, "n_se": 2.0},
    "parts": {"main": {}},
    "quick": {"n_seeds": 2, "anchors": 60, "diamonds": 40, "half": 6.0, "T": 7.0, "bins": [0, 1.5, 3.0, 4.5]},
    "estimate": "~25–45 min with 8–10 processes (64 3+1 orders of ~55 000 events, links ~2 min per order, cone "
                "readings 300 anchors per shell); memory ~0.5 GB per process",
}

CLASSES = ("fwd", "bwd", "perp", "par")


# ============================================================================================== the substrate rule

class Streams:
    """Rest-frame saturation sigma_i(rho) and 4-velocity (along x) of each stream; the per-link cone rule."""

    def __init__(self, streams, width):
        self.s0 = np.array([float(s) for s, _ in streams])
        self.v = np.array([float(v) for _, v in streams])
        self.G = 1.0 / np.sqrt(1.0 - self.v ** 2) if self.v.size else self.v
        self.w = float(width)

    def sig(self, P):
        """(n_streams, n) rest-frame saturations at spatial points P (n, 3)."""
        P = np.atleast_2d(P)
        prof = np.exp(-(P[:, 1] ** 2 + P[:, 2] ** 2) / (2.0 * self.w ** 2))
        return self.s0[:, None] * prof[None, :]

    def total(self, P):
        return self.sig(P).sum(axis=0) if self.s0.size else np.zeros(len(np.atleast_2d(P)))

    def in_cone(self, P, D):
        """Is the displacement D (n, 4) from the parent at spatial P (n, 3) within the rule's cone (future only)?"""
        D = np.atleast_2d(D)
        q = -D[:, 0] ** 2 + (D[:, 1:] ** 2).sum(axis=1)
        if self.s0.size:
            s = self.sig(P)
            st = s.sum(axis=0)
            safe = np.where(st > 0, st, 1.0)
            a = np.where(st > 0, s / safe * (1.0 - np.exp(-2.0 * st)), 0.0)          # (n_streams, n)
            ud = self.G[:, None] * (-D[None, :, 0] + self.v[:, None] * D[None, :, 1])  # u_a d^a per stream
            q = q + (a * ud ** 2).sum(axis=0)
        return (D[:, 0] > 0) & (q <= 0)


class StreamOrder(CapacityOrder):
    """Flat sprinkling; the link v -> k is kept iff it lies within the stream cone at the parent v."""

    def __init__(self, X, streams):
        self.X, self.N, self.st = X, len(X), streams
        ptr, idx = links_csr(X, 1.0)
        child = np.repeat(np.arange(self.N), np.diff(ptr))
        parent = idx
        keep = np.empty(parent.size, dtype=bool)
        for a in range(0, parent.size, 1_000_000):              # in chunks: ~7e6 links in 3+1
            pa, ch = parent[a:a + 1_000_000], child[a:a + 1_000_000]
            keep[a:a + 1_000_000] = streams.in_cone(X[pa, 1:], X[ch] - X[pa])
        self.n_links, self.n_kept = int(keep.size), int(keep.sum())
        parent, child = parent[keep], child[keep]
        o = np.argsort(parent, kind="stable")
        self.cptr = np.concatenate([[0], np.cumsum(np.bincount(parent, minlength=self.N))]).astype(np.int64)
        self.cidx = child[o]
        o = np.argsort(child, kind="stable")
        self.pptr = np.concatenate([[0], np.cumsum(np.bincount(child, minlength=self.N))]).astype(np.int64)
        self.pidx = parent[o]


def sprinkle_streams(streams, T, H, rng):
    """Poisson events of coordinate density exp(sigma_tot(rho)) in 3+1 (thinning), sorted by time."""
    smax = math.exp(float(streams.s0.sum()))
    n = rng.poisson(smax * T * (2 * H) ** 3)
    X = np.column_stack([rng.uniform(0, T, n)] + [rng.uniform(-H, H, n) for _ in range(3)])
    X = X[rng.random(n) < np.exp(streams.total(X[:, 1:])) / smax]
    return X[np.argsort(X[:, 0])]


# ============================================================================================== continuum predictions

def _metric(streams, P):
    """The cone metric G = eta + sum a_i u_i u_i at one spatial point, and sigma_tot."""
    G = np.diag([-1.0, 1.0, 1.0, 1.0])
    st = float(streams.total(P)[0])
    if st > 0:
        s = streams.sig(P)[:, 0]
        for si, gi, vi in zip(s, streams.G, streams.v):
            u = np.array([-gi, gi * vi, 0.0, 0.0])
            G = G + si / st * (1.0 - math.exp(-2.0 * st)) * np.outer(u, u)
    return G, st


def _sphere(n=20000):
    """Deterministic, nearly uniform directions (Fibonacci sphere)."""
    k = np.arange(n) + 0.5
    z = 1.0 - 2.0 * k / n
    ph = math.pi * (1.0 + math.sqrt(5.0)) * k
    r = np.sqrt(1.0 - z ** 2)
    return np.column_stack([z, r * np.cos(ph), r * np.sin(ph)])  # first component = n_x


def predict(streams, centres, p):
    """Per shell: d = sigma_tot; phi_clock = (1/4)[sigma_tot + 2 ln(-G_00) - (1/2) ln(-det G)] (the coordinate volume
    of a diamond of the cone metric G, times the density); s per direction class = solid-angle mean of ln c(n)."""
    n = _sphere()
    nx = n[:, 0]
    masks = {"fwd": nx >= p["cos_par"], "bwd": nx <= -p["cos_par"], "perp": np.abs(nx) <= p["sin_perp"]}
    masks["par"] = masks["fwd"] | masks["bwd"]
    out = {"d": [], "phic": []} | {c: [] for c in CLASSES}
    for rc in centres:
        P = np.array([[0.0, rc, 0.0]])
        G, st = _metric(streams, P)
        out["d"].append(st)
        out["phic"].append(0.25 * (st + 2.0 * math.log(-G[0, 0]) - 0.5 * math.log(-np.linalg.det(G))))
        A = np.einsum("ki,ij,kj->k", n, G[1:, 1:], n)
        B = 2.0 * n @ G[0, 1:]
        c = (-B + np.sqrt(B * B - 4.0 * A * G[0, 0])) / (2.0 * A)
        for k_, m in masks.items():
            out[k_].append(float(np.mean(np.log(c[m]))))
    return {k_: np.array(v_) for k_, v_ in out.items()}


# ============================================================================================== readers

def _rho(X):
    return np.hypot(X[:, 2], X[:, 3])


def read_density_cyl(X, edges, T, xcut):
    rho = _rho(X)
    m = (X[:, 0] >= 1.0) & (X[:, 0] <= T - 1.0) & (np.abs(X[:, 1]) <= xcut)
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        k = int(((rho >= lo) & (rho < hi) & m).sum())
        vol = math.pi * (hi ** 2 - lo ** 2) * 2.0 * xcut * (T - 2.0)
        out.append(math.log(k / vol) if k > 0 else float("nan"))
    return np.array(out)


def read_cones_dir(order, edges, T, xcut, p, rng):
    """Pooled boundary (K.16b) per shell and direction class; the relation is reachability through kept links."""
    X = order.X
    rho = _rho(X)
    dt_lo, dt_hi = p["cone_dt"]
    out = {c: [] for c in CLASSES}
    for lo, hi in zip(edges[:-1], edges[1:]):
        idx = np.flatnonzero((rho >= lo) & (rho < hi) & (np.abs(X[:, 1]) <= xcut)
                             & (X[:, 0] >= 1.0) & (X[:, 0] <= T - 1.0 - dt_hi))
        pool = {c: ([], []) for c in CLASSES}
        if idx.size:
            for i in rng.choice(idx, min(int(p["anchors"]), idx.size), replace=False):
                d = X - X[i]
                dt = d[:, 0]
                dx = np.sqrt((d[:, 1:] ** 2).sum(axis=1))
                cand = np.flatnonzero((dt >= dt_lo) & (dt <= dt_hi) & (dx <= float(p["cone_q"]) * dt))
                if cand.size == 0:
                    continue
                fut = order.future(order.cidx[order.cptr[i]:order.cptr[i + 1]], X[i, 0] + dt_hi)[cand]
                q = dx[cand] / dt[cand]
                nx = d[cand, 1] / np.maximum(dx[cand], 1e-12)
                cls = {"fwd": nx >= p["cos_par"], "bwd": nx <= -p["cos_par"], "perp": np.abs(nx) <= p["sin_perp"]}
                cls["par"] = cls["fwd"] | cls["bwd"]
                for c, m in cls.items():
                    pool[c][0].append(q[m])
                    pool[c][1].append(fut[m])
        for c in CLASSES:
            qs, oks = pool[c]
            qs = np.concatenate(qs) if qs else np.zeros(0)
            oks = np.concatenate(oks) if oks else np.zeros(0, bool)
            out[c].append(math.log(k16c.pooled_boundary(qs, oks)) if qs.size > 1 else float("nan"))
    return {c: np.array(v) for c, v in out.items()}


def read_clocks_cyl(order, edges, T, xcut, p, rng):
    """Diamonds at rest in the lab; virtual endpoints relate directly within the rule's cone (as F7.1a-c)."""
    X, st = order.X, order.st
    T_d = float(p["T_d"])
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        rc = 0.5 * (lo + hi)
        counts = []
        for _ in range(int(p["diamonds"])):
            ang = rng.uniform(0.0, 2.0 * math.pi)
            x0 = np.array([rng.uniform(-xcut, xcut), rc * math.cos(ang), rc * math.sin(ang)])
            t0 = rng.uniform(1.0, T - 1.0 - T_d)
            pp, qq = np.concatenate([[t0], x0]), np.concatenate([[t0 + T_d], x0])
            near = np.flatnonzero((X[:, 0] > t0) & (X[:, 0] < t0 + T_d)
                                  & (((X[:, 1:] - x0) ** 2).sum(axis=1) <= (0.5 * T_d) ** 2))
            if near.size == 0:
                counts.append(0)
                continue
            E = X[near]
            first = near[st.in_cone(np.repeat(x0[None, :], near.size, axis=0), E - pp)]
            last = near[st.in_cone(E[:, 1:], qq - E)]
            if first.size == 0 or last.size == 0:
                counts.append(0)
                continue
            counts.append(int((order.future(first, qq[0]) & order.past(last, pp[0])).sum()))
        tau = (24.0 * float(np.mean(counts)) / math.pi) ** 0.25
        out.append(math.log(tau / T_d) if tau > 0 else float("nan"))
    return np.array(out)


# ============================================================================================== task and evaluation

def order_task(p, streams_list, seed):
    rng = np.random.default_rng(seed)
    T, H = float(p["T"]), float(p["half"])
    xcut = H - float(p["x_margin"])
    st = Streams(streams_list, p["width"])
    X = sprinkle_streams(st, T, H, rng)
    order = StreamOrder(X, st)
    edges = np.array(p["bins"], float)
    out = {"N": len(X), "kept": order.n_kept / max(order.n_links, 1),
           "d": read_density_cyl(X, edges, T, xcut).tolist(),
           "phic": read_clocks_cyl(order, edges, T, xcut, p, rng).tolist()}
    for c, v in read_cones_dir(order, edges, T, xcut, p, rng).items():
        out[c] = v.tolist()
    return out


def _amp(phi, shape):
    phi, shape = np.asarray(phi, float), np.asarray(shape, float)
    m = np.isfinite(phi)
    return float((phi[m] * shape[m]).sum() / (shape[m] ** 2).sum())


def _aniso(row, k):
    par, perp = np.asarray(row["par"][:k], float), np.asarray(row["perp"][:k], float)
    return float(par.sum() / perp.sum())


def _ms(vals):
    v = np.asarray([x for x in vals if np.isfinite(x)], float)
    return float(v.mean()), float(v.std(ddof=1) / math.sqrt(v.size)) if v.size > 1 else float("nan")


def run(p, ctx):
    t0 = time.time()
    n = int(p["n_seeds"])
    jobs, keys = [], []
    for name, streams_list in p["configs"]:
        for k in range(n):
            jobs.append((p, streams_list, derive_seed(ctx.seed, "F7.4a", name, k)))
            keys.append(name)
    res = parallel_map(order_task, jobs, ctx.workers, ctx.progress,
                       lambda i, r_, dd, nn: f"[F7.4a] {keys[i]} sprinkling {dd}/{nn} ({time.time() - t0:.0f}s)")
    edges = np.array(p["bins"], float)
    centres = 0.5 * (edges[:-1] + edges[1:])
    rows = {name: [r_ for r_, k in zip(res, keys) if k == name] for name, _ in p["configs"]}
    pred = {name: predict(Streams(sl, p["width"]), centres, p) for name, sl in p["configs"]}
    nse, kk = float(p["n_se"]), int(p["k_shells"])
    shape = pred["two_v0"]["phic"]

    lines = [f"3+1, t [0, {p['T']:g}], |x|, |y|, |z| <= {p['half']:g}; streams along x, rest-frame saturation "
             f"sigma0 exp(-rho^2/(2*{p['width']:g}^2)); density exp(sigma_tot); link kept iff -dt^2 + |dx|^2 + "
             f"sum a_i (u_i.d)^2 <= 0 at the parent; reachability; {n} sprinklings per configuration"]
    st = {}
    for name, sl in p["configs"]:
        R_ = rows[name]
        x = {k_: (np.nanmean([r_[k_] for r_ in R_], axis=0), np.nanstd([r_[k_] for r_ in R_], axis=0, ddof=1)
                  / math.sqrt(len(R_))) for k_ in ("d", "phic") + CLASSES}
        x["A"] = _ms([_amp(r_["phic"], shape) for r_ in R_])
        x["K"] = _ms([_aniso(r_, kk) for r_ in R_])
        st[name] = x
        lines.append(f"  {name} {sl}: ~{np.mean([r_['N'] for r_ in R_]):.0f} events, links kept "
                     f"{np.mean([r_['kept'] for r_ in R_]):.3f}; amplitude A {x['A'][0]:+.4f} ± {x['A'][1]:.4f}; "
                     f"K (along/across) {x['K'][0]:.3f} ± {x['K'][1]:.3f}")
        lines.append("        rho            " + " ".join(f"{c:7.2f}" for c in centres))
        for k_ in ("d", "phic", "perp", "par", "fwd", "bwd"):
            lines.append(f"        {k_:5s} predicted " + " ".join(f"{v:+7.4f}" for v in pred[name][k_]))
            lines.append(f"        {k_:5s} read      " + " ".join(f"{v:+7.4f}" for v in x[k_][0]))

    # instrument
    worst, inst_null = -np.inf, True
    for k_ in ("d", "phic") + CLASSES:
        m, s = st["flat"][k_]
        excess = np.abs(m) - float(p["null_se"]) * s
        worst = max(worst, float(np.nanmax(excess)))
        inst_null &= bool(np.all(np.abs(m) <= float(p["null_abs"]) + float(p["null_se"]) * s))
    A0p = _amp(pred["two_v0"]["phic"], shape)
    (A0, sA0) = st["two_v0"]["A"]
    static_dev = A0 / A0p - 1.0
    inst_static = abs(static_dev) - nse * sA0 / abs(A0p) <= float(p["tol_static"])
    (K0, sK0) = st["two_v0"]["K"]
    inst_iso = abs(K0 - 1.0) - nse * sK0 <= float(p["tol_iso"])
    lines.append(f"  instrument (i) flat: every reading |mean| <= {p['null_abs']} + {p['null_se']:g} SE -> "
                 f"{'OK' if inst_null else 'FAILS'} (largest |mean| - {p['null_se']:g} SE = {worst:+.4f})")
    lines.append(f"  instrument (ii) static amplitude: A(two_v0) / predicted - 1 = {static_dev:+.3f} ± {sA0 / abs(A0p):.3f} -> "
                 f"{'OK' if inst_static else 'FAILS'} (tolerance {p['tol_static']})")
    lines.append(f"  instrument (iii) isotropic cones: K(two_v0) = {K0:.3f} ± {sK0:.3f} -> {'OK' if inst_iso else 'FAILS'} "
                 f"(needs |K - 1| - 2 SE <= {p['tol_iso']})")

    # amplitude
    (A6, sA6) = st["two_v06"]["A"]
    R = A6 / A0
    sR = abs(R) * math.hypot(sA6 / A6, sA0 / A0)
    R_T = _amp(pred["two_v06"]["phic"], shape) / A0p
    G6 = 1.0 / math.sqrt(1.0 - 0.6 ** 2)
    cands = {"tensor (exact, this rule)": R_T, "tensor (linear)": G6 ** 2 * (1 + 0.36), "energy density": G6 ** 2,
             "lab density": G6, "rest density": 1.0}
    nearest = min(cands, key=lambda c: abs(R - cands[c]))
    tol = float(p["tol_R"])
    amp = ("PASSED" if abs(R - R_T) + nse * sR <= tol else "FAILED" if abs(R - R_T) - nse * sR > tol else "UNDECIDED")
    lines.append(f"  AMPLITUDE R = A(v = 0.6) / A(0) = {R:.3f} ± {sR:.3f}; candidates: "
                 + ", ".join(f"{c} {v:.3f}" for c, v in cands.items()) + f"; nearest: {nearest}; against R_T: "
                 f"{R - R_T:+.3f} ± {sR:.3f} -> {amp}")

    # anisotropy
    K_T = float(pred["two_v06"]["par"][:kk].sum() / pred["two_v06"]["perp"][:kk].sum())
    (K6, sK6) = st["two_v06"]["K"]
    tol_K = (K_T - 1.0) / 2.0
    ani = ("PASSED" if abs(K6 - K_T) + nse * sK6 <= tol_K else "FAILED" if abs(K6 - K_T) - nse * sK6 > tol_K
           else "UNDECIDED")
    lines.append(f"  ANISOTROPY K = along / across = {K6:.3f} ± {sK6:.3f}; tensor K_T {K_T:.3f}, every mean-flow scalar "
                 f"rule 1; against K_T: {K6 - K_T:+.3f} ± {sK6:.3f} (tolerance {tol_K:.3f}) -> {ani}")

    # reported: the single boosted stream
    y = st["single_v06"]
    R1 = y["A"][0] / A0
    R1p = _amp(pred["single_v06"]["phic"], shape) / A0p
    fwd_r, bwd_r = float(np.nansum(y["fwd"][0][:kk])), float(np.nansum(y["bwd"][0][:kk]))
    fwd_p, bwd_p = float(pred["single_v06"]["fwd"][:kk].sum()), float(pred["single_v06"]["bwd"][:kk].sum())
    lines.append(f"  single stream v = 0.6 (reported): R = {R1:.3f} (predicted {R1p:.3f}); cone deficit backward / "
                 f"forward (inner shells) read {bwd_r / fwd_r if fwd_r else float('nan'):.2f}, predicted "
                 f"{bwd_p / fwd_p:.2f} (dragging of the cones; mean-flow scalar rules: 1)")

    if not (inst_null and inst_static and inst_iso):
        status = "INSTRUMENT_FAILED"
    elif amp == "PASSED" and ani == "PASSED":
        status = "PASSED"
    elif "FAILED" in (amp, ani):
        status = "FAILED"
    else:
        status = "UNDECIDED"
    data = {"R": [R, sR], "R_T": R_T, "candidates": cands, "K": [K6, sK6], "K_T": K_T, "amplitude": amp,
            "anisotropy": ani, "single_R": [R1, R1p],
            "configs": {name: {k_: [list(map(float, v[0])), list(map(float, v[1]))] if k_ not in ("A", "K") else list(v)
                               for k_, v in x_.items()} for name, x_ in st.items()},
            "predicted": {name: {k_: list(map(float, v)) for k_, v in pr.items()} for name, pr in pred.items()}}
    return Result(status, lines, data)
