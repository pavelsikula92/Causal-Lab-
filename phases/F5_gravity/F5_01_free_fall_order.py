"""
F5.1 -- Free fall in an order shaped by matter: does a flow fall without any bending rule?
==========================================================================================

Why this test exists (29 Oct 2026, the author's go-ahead for the F5 phase)
--------------------------------------------------------------------------
Gravity in Causal Physics changes the causal order itself: near matter the cones narrow and events are added
(throughput principle, F7.1b/c; tensor rule F7.4a). Matter is a flow that carries its 4-velocity (F3.6) and
merges on contact (F4.1). Falling has so far been tested only with an INSERTED bending rule (F3.8-F3.10: the
carried direction turned by the gradient of a scalar clock field, the geometry of gamma = -1). Nobody has asked
whether a flow falls in an order that matter has shaped, with no bending rule at all. The carried rule of F3.6
chooses the child of least rapidity relative to the carried direction, computed with the FLAT Minkowski product
of the sprinkling coordinates: it reads nothing from the order except which children exist. If such a flow does
not fall, the matter rule and the gravity of the order do not fit together, and the rule must be reformulated
in terms of the order (principle M3, which the F3.6 rule does not meet, as its own header says).

The order (locked; F7.1b in 2+1)
--------------------------------
2+1 box t [0, T], |x|, |y| <= H. Events: Poisson with coordinate density S(r) = exp(sigma(r)),
sigma = sigma0 exp(-r^2 / 2 w^2) (a static source at the centre); links of this sprinkling; the link v -> k is
kept iff its speed is <= c_eff(r_v) = exp(-sigma(r_v)) (phases/F7_light/F7_01a, CapacityOrder). The relation is
reachability through kept links. Continuum geometry of this order: clocks phi = -sigma/3, space psi = -2 sigma/3
(F7.1b; its blind reading gave gamma ~ 1.85 against 2). A slow body falls with the coordinate acceleration
    d^2 r / dt^2 = -d phi / dr = -(sigma0 r / 3 w^2) exp(-r^2 / 2 w^2)     (towards the centre).

Three arms on the same order (test chains of zero mass: one chain each, nothing merges)
-------------------------------------------------------------------------------------
  carried    the F3.6 rule: from v go to the kept child k of least rapidity relative to the carried direction u
             (flat Minkowski product); u is carried unchanged (one chain, no merging): u = (1, 0, 0) at rest.
  intrinsic  a rule that uses only the order and counts: the chain carries its ORIGIN event p as its state (in
             place of a direction). From v it goes to a kept child k that lies on a longest chain from p: the least
             longest-chain length L(p, k) among the children (normally L(p, v) + 1); among those, the one with most
             events in I(p, k); ties: the earliest child. Every prefix of the walk is then a longest chain from its
             origin, i.e. a discrete geodesic from p (F2: maximal chains follow geodesics): the local, order-only form
             of "matter follows the longest histories".
  maximal    control: the longest chain (most links) of the kept-link order between a start event p on the ring
             and the event q closest to the same spatial point a time dT later; a geodesic between two events at the
             same place rises away from the mass and comes back: mid-time sag OUTWARDS, +a dT^2/8 for constant a.
Start events: n_ring events nearest to points on the ring r = r0 at time t0, at equally spaced angles. The
origin of an intrinsic chain is the event of the start's past nearest to the same place a time anchor_dt earlier
(the initial condition "at rest"; the geodesic from the origin through the start already falls, which changes the
initial velocity, not the fitted acceleration).

Development, before the criteria (29 Oct, container, 2 sprinklings per arm, 12 chains, decides nothing): carried:
loaded R = 0.36 (r(t) rises by ~0.5 and stays), only 6.5 of 12 chains reach the end (some events have no kept child);
flat drifts by -0.5 over 22 (12 chains). A first order-only rule (maximise I(p, k) with p six steps back) did not keep
a velocity: in the flat order its chains left their place at ~0.4 c (the counts of short intervals are too noisy to
select a direction); a second one (go to the child that maximises |I(p, k)| from the origin p) jumped along the
long, nearly null links of the sprinkling (one link of time 11): maximising the interval rewards long steps, not
straight ones. Both were replaced, before any full run, by the longest-chain shell rule above. With it (1 sprinkling,
8 chains) the chains keep their place to +-1-2 over 11 time units; the readout by |r| still drifted outwards in both
orders (random walk, |r| bias), so the readout became the radial projection above. Maximal chains: loaded R_max
1.05 and 1.5, flat -0.26 +- 0.20 and -0.59 +- 0.34 (6 chains).
With the projection (2 sprinklings x 24 chains per arm and order): carried loaded R = 0.23 / 0.30, flat 0.07 / -0.33;
intrinsic loaded -0.03 / 0.80, flat -1.06 / 0.36 (the order-only walk is noisy: ~+-1 in R per sprinkling); maximal
chains loaded R_max 1.35 / 1.45, flat 1.95 / 0.12 (8 chains). Loaded chains often end early (an event with no kept
child): 14-18 of 24 reach the end. Changed after these runs, before any full run: 32 sprinklings, 40 chains, 40
maximal chains; chains that end count up to their end (no selection by survival); instrument tolerances 1/2 for the
flat nulls and [1/2, 2] for R_max (with 1/4 and [1/2, 3/2] the nulls fail on noise alone at this resolution).

Readout and predictions
-----------------------
For carried and intrinsic: the displacement of every chain from its start, projected on the outward radial unit
vector at the start (a projection: a transverse random walk would bias |r| outwards), against coordinate time; the
mean over chains is fitted by v0 (t - t0) + a/2 (t - t0)^2 over [t0, t0 + span]; the same fit applied to the numerical
solution of the fall equation above gives a_pred; the ratio R = a / a_pred (inwards positive for both). For
maximal: the mean radial deviation at mid-time of the longest chains, against the geodesic sag computed from the
same equation (numerically: the path at the same place in t0 and t0 + dT that maximises proper time), ratio R_max.
A flat order (sigma0 = 0) with the same machinery gives the null.

Criteria (stated before any full run)
-------------------------------------
Classification of an arm by its ratio R over sprinklings (bootstrap SE): falls  R - 2 SE > 1/2;  does not fall
R + 2 SE < 1/2; inconclusive otherwise.
Instrument (INSTRUMENT_FAILED if any fails): (i) the flat order: |R| + 2 SE <= 1/2 for the arm and for the maximal
chains (no drift or sag without matter at the resolution of the classification); (ii) the maximal chains read the
order's geometry: R_max - 2 SE > 1/2 and R_max + 2 SE < 2 (the discrete order realises the continuum geometry only
approximately: F7.1b read gamma ~ 1.85 against 2; the arms are classified against the continuum prediction and the
ratio to the maximal chains R / R_max is reported).
Status, per part:
  carried    PASSED if the carried arm falls; FAILED if it does not fall; UNDECIDED otherwise.
  intrinsic  PASSED if the intrinsic arm falls; FAILED if it does not fall; UNDECIDED otherwise.
Reading declared now. carried FAILED and intrinsic PASSED would mean: matter as defined so far (direction in
flat coordinates) does not feel the gravity of the order, and the rule that does is "follow the longest history
from your own past", which reads the geometry from counts (M3). Both PASSED: the flat-coordinate rule falls
because the order restricts its children; both FAILED: neither rule falls and the gravity of the order does not
reach matter of this kind.
"""

from __future__ import annotations

import math
import time

import numpy as np

from core.order.sprinkle import derive_seed
from core.stats.parallel import parallel_map
from core.stats.stats import boot_se
from core.testing import Result
from phases.F7_light.F7_01a_capacity_cones import CapacityOrder
from phases.F7_light.F7_01b_throughput_cones import sprinkle_saturated

TEST = {
    "id": "F5.1",
    "phase": "F5",
    "title": "Free fall in an order shaped by matter: does a flow fall without any bending rule?",
    "theory_part": "F5 gravity and equivalence; Axiom 6 (flows carry momentum) against Axiom 9 (the order responds to matter); M3",
    "question": "In an order shaped by a static source (throughput rule: narrower cones and more events), does a test flow fall "
                "towards the source (a) with the carried rule of F3.6 (direction in flat coordinates), (b) with a rule that "
                "reads only the order (go to the child that maximises the interval from the flow's own past)?",
    "method": "2+1, the F7.1b order (density exp(sigma), cones exp(-sigma), sigma = sigma0 exp(-r^2/2w^2)); test chains from a "
              "ring at rest; radial fall fitted by a quadratic and compared with the fall equation of the order's geometry "
              "(phi = -sigma/3); control: longest chains between two events at the same place (outward sag); flat null",
    "prediction": "geodesic fall with a = -dphi/dr for a rule that reads the order; the flat-coordinate rule may not fall",
    "criterion": "falls R - 2 SE > 1/2, does not fall R + 2 SE < 1/2. INSTRUMENT: flat |R| + 2 SE <= 1/2 (arm and maximal "
                 "chains); loaded maximal chains 1/2 < R_max -+ 2 SE < 2. carried / intrinsic part: PASSED its arm falls, FAILED it does not fall, "
                 "UNDECIDED otherwise.",
    "null_model": "flat order (sigma0 = 0): no fall",
    "level": "mechanism",
    "meters": ["K.16b"],
    "depends_on": ["F3.6", "F7.1b", "F2.5"],
    "params": {"T": 32.0, "half": 18.0, "width": 8.0, "sigma0": 0.45, "r0": 8.0, "t0": 7.0, "span": 22.0,
               "n_ring": 40, "anchor_dt": 6.0, "dT": 22.0, "n_max": 40, "n_seeds": 32, "n_se": 2.0,
               "boot_B": 4000},
    "parts": {"carried": {"arm": "carried"}, "intrinsic": {"arm": "intrinsic"}},
    "quick": {"n_seeds": 2, "n_ring": 12, "n_max": 6, "boot_B": 200},
    "estimate": "carried ~10-20 min, intrinsic ~30-60 min with 10 processes (64 orders of ~48 000 events; bitset closures of "
                "the intrinsic rule; 40 longest chains per order)",
}


# ===================================================================================== geometry and prediction

def accel(r, p):
    """Coordinate acceleration (radial, outward positive) of a slow body: -dphi/dr with phi = -sigma/3."""
    s0, w = float(p["sigma0"]), float(p["width"])
    return -(s0 * r / (3.0 * w * w)) * math.exp(-r * r / (2.0 * w * w))


def fall_curve(p, ts):
    """Numerical fall from rest at r0 (RK4), r(t) - r0 at the times ts - t0."""
    r, v, t, dt = float(p["r0"]), 0.0, 0.0, 0.01
    out, k = [], 0
    ts = np.asarray(ts, float) - float(p["t0"])
    for tt in ts:
        while t < tt - 1e-12:
            h = min(dt, tt - t)
            k1v, k1x = accel(r, p), v
            k2v, k2x = accel(r + 0.5 * h * k1x, p), v + 0.5 * h * k1v
            k3v, k3x = accel(r + 0.5 * h * k2x, p), v + 0.5 * h * k2v
            k4v, k4x = accel(r + h * k3x, p), v + h * k3v
            r += h * (k1x + 2 * k2x + 2 * k3x + k4x) / 6.0
            v += h * (k1v + 2 * k2v + 2 * k3v + k4v) / 6.0
            t += h
        out.append(r - float(p["r0"]))
    return np.array(out)


def sag_prediction(p):
    """Mid-time radial deviation of the proper-time-maximising path between (t0, r0) and (t0 + dT, r0): the trajectory
    launched upward with the v0 that returns to r0 at dT (shooting), its deviation at dT/2 (outward positive)."""
    dT = float(p["dT"])
    q = dict(p)

    def end(v0):
        r, v, t, h = float(p["r0"]), v0, 0.0, 0.01
        mid = None
        while t < dT - 1e-12:
            a = accel(r, q)
            v += a * h
            r += v * h
            t += h
            if mid is None and t >= dT / 2:
                mid = r
        return r - float(p["r0"]), mid - float(p["r0"])
    lo, hi = 0.0, 1.0
    for _ in range(60):
        m = 0.5 * (lo + hi)
        if end(m)[0] > 0:
            hi = m
        else:
            lo = m
    return end(0.5 * (lo + hi))[1]


# ===================================================================================== rules

def _mink(a, b):
    return a[..., 0] * b[..., 0] - (a[..., 1:] * b[..., 1:]).sum(axis=-1)


def walk_carried(order, v0, t_end):
    """F3.6 carried rule for one chain at rest: the kept child of least rapidity relative to u = (1, 0, 0)."""
    X, cptr, cidx = order.X, order.cptr, order.cidx
    u = np.array([1.0, 0.0, 0.0])
    path, v = [v0], v0
    while X[v, 0] < t_end:
        a, b = cptr[v], cptr[v + 1]
        if b == a:
            break
        ks = cidx[a:b]
        B = X[ks] - X[v]
        tau = np.sqrt(np.maximum(_mink(B, B), 1e-300))
        v = int(ks[int(np.argmin(_mink(u[None, :], B) / tau))])
        path.append(v)
    return path


def past_closures(order, origin, t_end):
    """For every event of the kept-link future of the origin (up to t_end), the bitset of the events of that future
    that lie in its causal past, itself included (dynamic programming in time order over the parents). The number of
    events strictly between the origin and k is then popcount(C[k]) - 1."""
    X = order.X
    fut = order.future([origin], t_end)
    fut[origin] = False
    ids = np.where(fut)[0]
    ids = ids[np.argsort(X[ids, 0], kind="stable")]
    pos = np.full(order.N, -1, dtype=np.int64)
    pos[ids] = np.arange(ids.size)
    nb = (ids.size + 7) // 8
    C = np.zeros((ids.size, nb), dtype=np.uint8)
    L = np.zeros(ids.size, dtype=np.int64)               # longest chain (links) from the origin
    for j, w in enumerate(ids):
        par = order.pidx[order.pptr[w]:order.pptr[w + 1]]
        pj = pos[par]
        pj = pj[pj >= 0]
        if pj.size:
            C[j] = np.bitwise_or.reduce(C[pj], axis=0)
            L[j] = 1 + int(L[pj].max())
        else:
            L[j] = 1                                     # a child of the origin itself
        C[j, j >> 3] |= np.uint8(1 << (7 - (j & 7)))
    return pos, C, L


def walk_intrinsic(order, v0, origin, t_end):
    """Order-only rule: the chain carries its ORIGIN event; at every step it goes to a kept child that lies on a longest
    chain from the origin (the least longest-chain length among the children, normally L(v) + 1), among those to the
    one with most events between the origin and itself; ties: the earliest child (module notes)."""
    X, cptr, cidx = order.X, order.cptr, order.cidx
    pos, C, L = past_closures(order, origin, t_end + 2.0)
    path, v = [v0], v0
    while X[v, 0] < t_end:
        a, b = cptr[v], cptr[v + 1]
        if b == a:
            break
        ks = cidx[a:b]
        pk = pos[ks]
        ok = pk >= 0
        if not ok.any():
            break
        ks, pk = ks[ok], pk[ok]
        Lk = L[pk]
        sel = Lk == Lk.min()                             # the next shell of the longest chains from the origin
        ks, pk = ks[sel], pk[sel]
        n = np.unpackbits(C[pk], axis=1).sum(axis=1)
        cand = ks[n == n.max()]
        v = int(cand[np.argmin(X[cand, 0])])
        path.append(v)
    return path


def longest_chain(order, p, q):
    """The chain p -> q with most kept links (DP over the interval in time order)."""
    X, N = order.X, order.N
    fut = order.future([p], float(X[q, 0]))
    past = order._bfs([q], order.pptr, order.pidx, X[:, 0] >= float(X[p, 0]))
    inside = np.where(fut & past)[0]
    inside = inside[np.argsort(X[inside, 0])]
    dist = {p: 0}
    prev = {}
    for w in inside:
        if w == p:
            continue
        best, arg = -1, -1
        for u in order.pidx[order.pptr[w]:order.pptr[w + 1]]:
            du = dist.get(int(u))
            if du is not None and du > best:
                best, arg = du, int(u)
        if arg >= 0:
            dist[int(w)] = best + 1
            prev[int(w)] = arg
    if q not in dist:
        return None
    path, w = [q], q
    while w != p:
        w = prev[w]
        path.append(w)
    return path[::-1]


# ===================================================================================== one sprinkling

def ring_starts(X, p, n):
    """The events nearest to n equally spaced points of the ring r = r0 at time t0."""
    r0, t0 = float(p["r0"]), float(p["t0"])
    out = []
    for j in range(n):
        ang = 2.0 * math.pi * j / n
        target = np.array([t0, r0 * math.cos(ang), r0 * math.sin(ang)])
        out.append(int(np.argmin(((X - target) ** 2).sum(axis=1))))
    return out


def radial_series(X, paths, p, grid):
    """Mean over chains of the displacement from the start projected on the outward radial unit vector at the start
    (a projection, not |r|: a transverse random walk then averages to zero instead of biasing r outwards)."""
    cols, full = [], 0
    for path in paths:
        P = X[path]
        rhat = P[0, 1:] / np.linalg.norm(P[0, 1:])
        proj = (P[:, 1:] - P[0, 1:]) @ rhat
        y = np.interp(grid, P[:, 0], proj)
        y[grid > P[-1, 0]] = np.nan                      # a chain that ended (no kept child) counts up to its end
        cols.append(y)
        full += int(P[-1, 0] >= grid[-1])
    with np.errstate(all="ignore"):
        return (np.nanmean(cols, axis=0) if cols else np.full(len(grid), np.nan)), full


def fit_a(t, y, t0):
    """a from y = v0 (t - t0) + a/2 (t - t0)^2 (least squares, no constant: y(t0) = 0)."""
    s = np.asarray(t) - t0
    A = np.c_[s, 0.5 * s * s]
    c, *_ = np.linalg.lstsq(A, np.asarray(y), rcond=None)
    return float(c[1]), float(c[0])


def seed_task(p, sigma0, arm, seed):
    rng = np.random.default_rng(seed)
    T, H, w = float(p["T"]), float(p["half"]), float(p["width"])
    X = sprinkle_saturated(sigma0, w, T, H, rng)
    order = CapacityOrder(X, sigma0, w)
    t0, span = float(p["t0"]), float(p["span"])
    grid = np.linspace(t0, t0 + span, 23)
    starts = ring_starts(X, p, int(p["n_ring"]))
    out = {"N": len(X), "sigma0": sigma0}
    if arm == "carried":
        paths = [walk_carried(order, s_, t0 + span + 1.0) for s_ in starts]
    else:
        paths = []
        for s_ in starts:
            ta = X[s_, 0] - float(p["anchor_dt"])
            past = order._bfs([s_], order.pptr, order.pidx, X[:, 0] >= ta - 1.0)
            cand = np.where(past & (np.abs(X[:, 0] - ta) < 1.0))[0]
            anchor = int(cand[np.argmin(((X[cand, 1:] - X[s_, 1:]) ** 2).sum(axis=1))]) if cand.size else s_
            paths.append(walk_intrinsic(order, s_, anchor, t0 + span + 1.0))
    y, n_ok = radial_series(X, paths, p, grid)
    ok = np.isfinite(y)
    a, v0 = fit_a(grid[ok], y[ok], t0)
    out.update({"a": a, "v0": v0, "n_chains": n_ok, "series": y.tolist(), "grid": grid.tolist()})
    # maximal chains (control)
    dT = float(p["dT"])
    devs = []
    for s_ in starts[: int(p["n_max"])]:
        tq = X[s_, 0] + dT
        target = np.array([tq, X[s_, 1], X[s_, 2]])
        fut = order.future([s_], tq + 1.0)
        cand = np.where(fut & (np.abs(X[:, 0] - tq) < 1.0))[0]
        if cand.size == 0:
            continue
        q = int(cand[np.argmin(((X[cand] - target) ** 2).sum(axis=1))])
        path = longest_chain(order, s_, q)
        if path is None:
            continue
        P = X[path]
        rhat = P[0, 1:] / np.linalg.norm(P[0, 1:])
        proj = (P[:, 1:] - P[0, 1:]) @ rhat                 # outward radial projection, as for the walks
        tm = 0.5 * (P[0, 0] + P[-1, 0])
        rline = np.interp(tm, [P[0, 0], P[-1, 0]], [proj[0], proj[-1]])
        devs.append(float(np.interp(tm, P[:, 0], proj) - rline))
    out["sag"] = float(np.mean(devs)) if devs else float("nan")
    out["n_max"] = len(devs)
    return out


# ===================================================================================== run

def classify(R, se, nse):
    if not (np.isfinite(R) and np.isfinite(se)):
        return "undefined"
    if R - nse * se > 0.5:
        return "falls"
    if R + nse * se < 0.5:
        return "does not fall"
    return "inconclusive"


def run(p, ctx):
    t0 = time.time()
    arm = p["arm"]
    n, B, nse = int(p["n_seeds"]), int(p["boot_B"]), float(p["n_se"])
    grid = np.linspace(float(p["t0"]), float(p["t0"]) + float(p["span"]), 23)
    a_pred, _ = fit_a(grid, fall_curve(p, grid), float(p["t0"]))
    sag_pred = sag_prediction(p)
    jobs = []
    for s0 in (float(p["sigma0"]), 0.0):
        jobs += [(p, s0, arm, derive_seed(ctx.seed, "F5.1", ctx.part, s0, k)) for k in range(n)]
    rows = parallel_map(seed_task, jobs, ctx.workers, ctx.progress,
                        lambda i, r_, d, nn: f"[F5.1] {arm} {d}/{nn}: sigma0 {r_['sigma0']:g} a {r_['a']:+.4f} sag "
                                             f"{r_['sag']:+.3f} ({time.time() - t0:.0f}s)")
    res = {}
    for s0 in (float(p["sigma0"]), 0.0):
        rr = [r_ for r_ in rows if r_["sigma0"] == s0]
        bs = lambda *tag: derive_seed(ctx.seed, "F5.1-boot", ctx.part, s0, *tag)
        R = boot_se([r_["a"] / a_pred for r_ in rr], B, bs("R"))[:2]
        Rm = boot_se([r_["sag"] / sag_pred for r_ in rr if np.isfinite(r_["sag"])], B, bs("Rm"))[:2]
        v0 = boot_se([r_["v0"] for r_ in rr], B, bs("v0"))[:2]
        res[s0] = {"R": R, "R_max": Rm, "v0": v0, "series": np.nanmean([r_["series"] for r_ in rr], axis=0).tolist(),
                   "N": float(np.mean([r_["N"] for r_ in rr])), "chains": float(np.mean([r_["n_chains"] for r_ in rr]))}
    s0 = float(p["sigma0"])
    lines = [f"2+1 order of F7.1b: density exp(sigma), cones exp(-sigma), sigma0 = {s0:g}, w = {p['width']:g}; box t [0, {p['T']:g}], "
             f"|x|, |y| <= {p['half']:g} (~{res[s0]['N']:.0f} events); ring r0 = {p['r0']:g} at t0 = {p['t0']:g}, {p['n_ring']} "
             f"chains, fit over {p['span']:g}; arm {arm}",
             f"  prediction (phi = -sigma/3): a_pred = {a_pred:+.5f} (fit of the fall curve; inwards negative); maximal-chain "
             f"sag at mid-time {sag_pred:+.4f} (outwards positive)"]
    for key, lab in ((s0, "loaded"), (0.0, "flat  ")):
        r_ = res[key]
        lines.append(f"  {lab}: R = a / a_pred = {r_['R'][0]:+.3f} +- {r_['R'][1]:.3f} ({classify(*r_['R'], nse)}); v0 "
                     f"{r_['v0'][0]:+.4f} +- {r_['v0'][1]:.4f}; maximal chains R_max = {r_['R_max'][0]:+.3f} +- "
                     f"{r_['R_max'][1]:.3f}; chains used {r_['chains']:.1f}")
        lines.append("          r(t) - r0: " + " ".join(f"{x:+.2f}" for x in r_["series"][::2]))
    f_ = res[0.0]
    inst = {"flat: no drift (|R| + 2 SE <= 1/2)": abs(f_["R"][0]) + nse * f_["R"][1] <= 0.5,
            "flat maximal chains: no sag (|R_max| + 2 SE <= 1/2)": abs(f_["R_max"][0]) + nse * f_["R_max"][1] <= 0.5,
            "loaded maximal chains read the geometry (1/2 < R_max -+ 2 SE < 2)":
                res[s0]["R_max"][0] - nse * res[s0]["R_max"][1] > 0.5 and res[s0]["R_max"][0] + nse * res[s0]["R_max"][1] < 2.0}
    lines.append("instrument: " + "; ".join(f"{k}: {'OK' if ok else 'FAILS'}" for k, ok in inst.items()))
    cls = classify(*res[s0]["R"], nse)
    if not all(inst.values()):
        status = "INSTRUMENT_FAILED"
    else:
        status = {"falls": "PASSED", "does not fall": "FAILED"}.get(cls, "UNDECIDED")
    data = {"a_pred": a_pred, "sag_pred": sag_pred, "res": {str(k): v for k, v in res.items()}, "rows": rows}
    return Result(status, lines, data)
