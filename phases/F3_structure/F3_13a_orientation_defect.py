"""
F3.13a -- A topological defect of a carried orientation: is its winding number an invariant of the flow?
========================================================================================================

Why this test exists (31 Oct 2026, the author's decision after F3.12b)
---------------------------------------------------------------------
F3.12b closed the search for a particle as a bundle held by its own field: neither capacity nor internal motion gives
a stationary width. The author's new principle: a particle as a topologically locked defect of the causal structure,
stable not because forces balance but because there is no cheap local way back to the vacuum. His condition, adopted
here: first an invariant that tells the defect from the vacuum, only then the cost of erasing it (F3.13b, F3.13c later).
Of the two candidates discussed (a conical defect of the geometry; a defect of a direction field carried by the flow)
the author chose the second: it is a point defect in 3+1 (a "hedgehog"), and its charge is an integer.
The carried direction of F3.6 is the direction of motion (momentum); a hedgehog in it would be a radial outflow, an
explosion, not a particle at rest. The author therefore introduced (31 Oct) a NEW carried state: an internal
orientation n, a unit spacelike vector orthogonal to the 4-velocity u of the flow, carried unchanged from event to
event (like the direction in F3.6). A hedgehog of n is static: the flow rests, the orientation points outwards.

Model (locked before any run)
-----------------------------
Orientation: every flow carries (u, n) with u . u = 1, n . n = -1, n . u = 0. Transport (cold, no field): the F3.6
rule, an event sends everything to its child of least rapidity relative to u; flows that meet merge as in F3.6c rule B
(P = sum of the 4-momenta, F = |P|, u = P / F). The merge rule for n (the author's choice, 31 Oct): the flow-weighted
sum of the incoming orientations, each projected onto the rest space of the merged u, normalised:
    n = N[ sum_i F_i (n_i - (n_i . u) u) ],
and if that sum is shorter than 1e-6 of sum F_i the event has no orientation (a core event: the singular point of the
defect). An event with one incoming flow keeps its n.
Invariant (the degree of the map from a shell around the centre to the sphere of orientations):
  2+1: the events of a ring R - 3/4 <= r < R + 3/4 around the centre, sorted by polar angle; the winding number is the
       sum of the wrapped differences of the polar angle of n between neighbours, over 2 pi;
  3+1: the events of a shell R - 3/4 <= r < R + 3/4, their directions from the centre triangulated by the convex hull
       (outward orientation); the degree is the sum over triangles of the signed solid angle of the image triangle
       (n_a, n_b, n_c), Omega = 2 atan2(n_a . (n_b x n_c), 1 + n_a . n_b + n_b . n_c + n_c . n_a), over 4 pi.
  A ring or shell is readable if it has at least 6 (2+1) / 12 (3+1) events with an orientation and no angular gap of
  the positions larger than 2 pi / 3 (2+1), or no hull facet closer to the centre than 0.3 (3+1), and if the field is
  resolved on it: neighbouring orientations (angular neighbours on a ring, hull edges on a shell) differ by less than
  pi / 2, so that no wrap is ambiguous; core events are left out. The value is an integer by construction; the test
  reads which integer.
Configurations of n at the sources (spatial components in the source frame; the flows start at rest):
  vacuum   n = e_x everywhere (degree 0);
  defect   n = r_hat (2+1: a vortex, winding +1; 3+1: a hedgehog, degree +1);
  anti     2+1 (x, -y)/r, 3+1 (-x, y, z)/r (degree -1);
  double   2+1 angle 2 theta, 3+1 (sin T cos 2 P, sin T sin 2 P, cos T) in spherical angles (degree 2; reported);
  twist    n = e_x rotated by alpha(r) = pi (1 - r / r_src) in a plane containing e_x: strongly non-uniform but of
           degree 0 on every shell (the control: a distorted region without an invariant).
Parts:
  static_d3, static_d4: one slab of a Poisson sprinkling (density 1) in a ball of radius 6, n assigned directly;
           shells R = 3, 4, 5 (R = 2 holds too few events of a ring in 2+1); noise: each n rotated by a random angle uniform in [0, theta] about a random axis,
           theta = 0, pi/8, pi/4, 3 pi/8, pi/2 (theta = 0 decides, the rest reported); 32 sprinklings.
  flow_d3: 2+1 box t [0, 24], half 10, links (core.order.links); sources: events with t < 2, r < 5, unit flow at rest
           with n of the configuration; readout per slab of width 2 from t = 4 to 22, rings R = 2, 3, 4 around the
           flow-weighted centroid of the slab; 32 sprinklings.
  flow_d4: 3+1 box t [0, 12], half 7, sources t < 2, r < 5; slabs from t = 4 to 10; shells R = 2, 3, 4; 16 sprinklings.

Criteria (stated before any run)
--------------------------------
static (instrument): at theta = 0, vacuum 0, defect +1, anti -1, double +2, twist 0 on every readable shell of every
  sprinkling, and at least 90 % of the shells readable for vacuum, defect, anti and twist (double, whose field turns
  twice as fast, has its readability reported); PASSED if so, INSTRUMENT_FAILED otherwise. The value of theta
  up to which the defect still reads +1 on >= 95 % of shells is reported.
flow (the question): over all readable (sprinkling, slab, shell) triples of an arm, the fraction that reads the degree
  of the configuration. INSTRUMENT_FAILED if fewer than 50 % of the triples of the defect arm are readable, or vacuum
  or twist read their 0 on fewer than 95 %. PASSED: defect and anti read +1 and -1 on >= 95 % of their readable
  triples; FAILED: defect or anti below 80 %; UNDECIDED otherwise. The double arm is reported (a degree-2 defect may
  split into two unit defects: shells around both still read 2).
Development (31 Oct, before any full run; decides nothing). Quick runs showed two faults of the first readout, fixed
before any full run: rings of width 1 in 2+1 were too sparse at R = 2 (now width 1.5 as the shells), and with sparse
carriers a winding could be counted wrongly when neighbouring orientations differ by more than pi/2 (now such a ring or
shell is unreadable instead: the field-resolution rule above). One full-size row per flow part:
  flow_d3 (9 631 events, 168 sources, 7 s): readable 23 / 16 / 16 / 2 / 16 of 27 (vacuum, defect, anti, double,
    twist); every readable value is the degree of its configuration; no core event.
  flow_d4 (32 799 events, 1 104 sources, 104 s): readable 9 / 5 / 5 / 0 / 2 of 9; every readable value is the degree.
The cold flow moves almost straight up and merges only neighbouring streams, so in this set-up the charge is mainly
carried, not challenged: a PASSED here is the necessary condition (the transport of the laboratory, merging included,
does not destroy the charge), not yet the locking against a perturbation that tries to unwind it (F3.13b).

Reading declared now. PASSED in flow: the integer charge of a defect of the carried orientation survives the
transport of the lab, merging and the rebuilding of the core included; the defect is a locked configuration of the
flow, and F3.13b may measure the cost of erasing it and its motion. FAILED: merging (or the sparse carriers of a cold
bundle) destroys the charge; the orientation would need another merge rule or a carrier that does not merge.
Declared limits: the orientation is a new assumption of the theory (the author's decision of 31 Oct); the defect is
put in by hand (whether the dynamics creates one is a later question); no field, no interaction between orientations.
"""

from __future__ import annotations

import math
import time

import numpy as np
from scipy.spatial import ConvexHull

from core.order.links import links_csr
from core.order.sprinkle import Box, derive_seed, sprinkle_poisson
from core.stats.parallel import parallel_map
from core.testing import Result
from phases.F3_structure.F3_08_retarded_field import children_csr, mink

TEST = {
    "id": "F3.13a",
    "phase": "F3",
    "title": "A topological defect of a carried orientation: is its winding number an invariant of the flow?",
    "theory_part": "F3 emergence of structure; a particle as a topologically locked defect (the author's principle of "
                   "31 Oct); a new carried state: the internal orientation n of a flow",
    "question": "Does a vortex (2+1) or hedgehog (3+1) of the carried orientation keep its integer degree on shells "
                "around the centre while the cold flow of the laboratory carries it, merges and rebuilds its core?",
    "method": "static parts: degree of the map shell -> sphere of orientations on sprinkled shells, configurations "
              "vacuum / defect / anti / double / twist, noise sweep; flow parts: F3.6 cold transport with merging "
              "(rule B) and the flow-weighted merge of n, readout per slab on rings (2+1) and shells (3+1)",
    "prediction": "the defect keeps its degree +1 on readable shells (topological charge of the flow)",
    "criterion": "static: theta = 0 reads 0 / +1 / -1 / +2 / 0 on every readable shell, >= 90 % readable (else "
                 "INSTRUMENT_FAILED); flow: PASSED defect and anti >= 95 % of readable triples, FAILED < 80 %, "
                 "INSTRUMENT_FAILED readable < 50 % or vacuum / twist < 95 %, otherwise UNDECIDED",
    "null_model": "vacuum (uniform n) and twist (distorted, degree 0)",
    "level": "mechanism",
    "meters": [],
    "depends_on": ["F3.6", "F3.6c", "F3.12b"],
    "params": {"r_ball": 6.0, "shells": [3.0, 4.0, 5.0], "thetas": [0.0, 0.3927, 0.7854, 1.1781, 1.5708],
               "n_static": 32, "T3": 24.0, "half3": 10.0, "T4": 12.0, "half4": 7.0, "r_src": 5.0, "t_src": 2.0,
               "slab": 2.0, "t_read": 4.0, "rings": [2.0, 3.0, 4.0], "n_flow3": 32, "n_flow4": 16,
               "read_static": 0.9, "read_flow": 0.5, "pass_frac": 0.95, "fail_frac": 0.8, "core_eps": 1e-6},
    "parts": {"static_d3": {"mode": "static", "D": 3}, "static_d4": {"mode": "static", "D": 4},
              "flow_d3": {"mode": "flow", "D": 3}, "flow_d4": {"mode": "flow", "D": 4}},
    "quick": {"n_static": 2, "n_flow3": 1, "n_flow4": 1, "T3": 12.0, "half3": 7.0, "T4": 8.0, "half4": 6.0,
              "r_src": 4.0, "rings": [2.0, 3.0]},
    "estimate": "static: seconds; flow_d3 ~1 min; flow_d4 ~5-10 min with 8 processes (16 orders of ~33 000 events in 3+1)",
}

CONFIGS = ["vacuum", "defect", "anti", "double", "twist"]
EXPECT = {"vacuum": 0, "defect": 1, "anti": -1, "double": 2, "twist": 0}


# ============================================================================================== orientations

def orient(cfg, P, r_src):
    """Spatial orientation (unit vectors, rows) of configuration cfg at spatial points P (relative to the centre)."""
    P = np.asarray(P, float)
    d = P.shape[1]
    r = np.linalg.norm(P, axis=1)
    rs = np.maximum(r, 1e-12)
    ex = np.zeros_like(P)
    ex[:, 0] = 1.0
    if cfg == "vacuum":
        return ex
    if cfg == "defect":
        return P / rs[:, None]
    if cfg == "anti":
        Q = P.copy()
        if d == 2:
            Q[:, 1] = -Q[:, 1]
        else:
            Q[:, 0] = -Q[:, 0]
        return Q / rs[:, None]
    if cfg == "double":
        if d == 2:
            a = 2.0 * np.arctan2(P[:, 1], P[:, 0])
            return np.column_stack([np.cos(a), np.sin(a)])
        th = np.arccos(np.clip(P[:, 2] / rs, -1.0, 1.0))
        ph = 2.0 * np.arctan2(P[:, 1], P[:, 0])
        return np.column_stack([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)])
    if cfg == "twist":
        a = math.pi * np.clip(1.0 - r / r_src, 0.0, 1.0)
        out = np.zeros_like(P)
        out[:, 0] = np.cos(a)
        out[:, 1] = np.sin(a)                         # rotation in the (e_x, e_y) plane: constant on every shell
        return out
    raise ValueError(cfg)


def noise(n, theta, rng):
    """Each unit vector rotated by a random angle uniform in [0, theta] about a random axis orthogonal to it."""
    if theta <= 0.0:
        return n
    d = n.shape[1]
    ang = rng.uniform(0.0, theta, len(n))
    w = rng.normal(size=n.shape)
    w -= (w * n).sum(axis=1)[:, None] * n
    w /= np.maximum(np.linalg.norm(w, axis=1), 1e-12)[:, None]
    return np.cos(ang)[:, None] * n + np.sin(ang)[:, None] * w


# ============================================================================================== the invariant

def winding2(pos, n):
    """Winding number on a ring (2+1): (value, readable)."""
    if len(pos) < 6:
        return float("nan"), False
    phi = np.arctan2(pos[:, 1], pos[:, 0])
    o = np.argsort(phi)
    gaps = np.diff(np.r_[phi[o], phi[o][0] + 2.0 * math.pi])
    if gaps.max() > 2.0 * math.pi / 3.0:
        return float("nan"), False
    a = np.arctan2(n[o, 1], n[o, 0])
    da = np.diff(np.r_[a, a[0]])
    da = (da + math.pi) % (2.0 * math.pi) - math.pi
    if np.abs(da).max() >= 0.5 * math.pi:              # the field is not resolved: a wrap could be ambiguous
        return float("nan"), False
    return float(da.sum() / (2.0 * math.pi)), True


def degree3(pos, n):
    """Degree of the map shell -> S^2 (3+1): (value, readable)."""
    if len(pos) < 12:
        return float("nan"), False
    p = pos / np.linalg.norm(pos, axis=1)[:, None]
    try:
        hull = ConvexHull(p)
    except Exception:
        return float("nan"), False
    if hull.equations[:, 3].max() > -0.3:             # a facet closer to the centre than 0.3: a hole in the shell
        return float("nan"), False
    S = hull.simplices
    for e0, e1 in ((0, 1), (1, 2), (2, 0)):            # the field is not resolved if neighbours differ by >= pi/2
        if (np.einsum("ij,ij->i", n[S[:, e0]], n[S[:, e1]]) <= 0.0).any():
            return float("nan"), False
    tot = 0.0
    for i, j, k in hull.simplices:
        if np.dot(np.cross(p[j] - p[i], p[k] - p[i]), p[i] + p[j] + p[k]) < 0.0:
            j, k = k, j
        a, b, c = n[i], n[j], n[k]
        num = float(np.dot(a, np.cross(b, c)))
        den = 1.0 + float(np.dot(a, b) + np.dot(b, c) + np.dot(c, a))
        tot += 2.0 * math.atan2(num, den)
    return float(tot / (4.0 * math.pi)), True


def shell_value(pos, n, d):
    return winding2(pos, n) if d == 2 else degree3(pos, n)


# ============================================================================================== static part

def static_task(p, D, seed):
    rng = np.random.default_rng(seed)
    d = D - 1
    R = float(p["r_ball"])
    P = rng.uniform(-R, R, size=(rng.poisson((2 * R) ** d), d))   # Poisson, density 1, in the cube; the ball kept
    P = P[np.linalg.norm(P, axis=1) < R]
    r = np.linalg.norm(P, axis=1)
    half = 0.75
    out = {}
    for cfg in CONFIGS:
        base = orient(cfg, P, R)
        for th in p["thetas"]:
            nn = noise(base, float(th), np.random.default_rng(derive_seed(seed, cfg, float(th))))
            vals = []
            for Rs in p["shells"]:
                sel = np.abs(r - Rs) < half
                v, ok = shell_value(P[sel], nn[sel], d)
                vals.append((float(Rs), v, ok))
            out[f"{cfg}|{float(th):g}"] = vals
    return out


# ============================================================================================== flow part

def transport_orient(X, cptr, cidx, s, U0, N0, eps):
    """F3.6 cold transport with merging (rule B) and the carried orientation (module notes)."""
    N, D = X.shape
    Pin = np.zeros((N, D))
    Nin = np.zeros((N, D))                            # sum of F_i n_i (projected at the event)
    Fin = np.zeros(N)
    nin = np.zeros(N, dtype=np.int64)
    Phold = np.zeros((N, D))
    Nhold = np.full((N, D), np.nan)
    core = np.zeros(N, dtype=bool)
    nlist = [[] for _ in range(N)]
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
            proj = ni - mink(ni, u) * u
            acc += Fi * proj
            tot += Fi
        L = math.sqrt(max(-float(mink(acc, acc)), 0.0))
        Phold[v] = P
        if L < eps * tot:
            core[v] = True
        else:
            Nhold[v] = acc / L
        a, b = cptr[v], cptr[v + 1]
        if b > a:
            ks = cidx[a:b]
            B = X[ks] - X[v]
            tau = np.sqrt(np.maximum(mink(B, B), 1e-300))
            k = int(ks[int(np.argmin(mink(u[None, :], B) / tau))])
            Pin[k] += P
            Fin[k] += F
            nin[k] += 1
            if not core[v]:
                nlist[k].append((F, Nhold[v]))
    return Phold, Nhold, core


def flow_task(p, D, seed):
    rng = np.random.default_rng(seed)
    T = float(p["T3"] if D == 3 else p["T4"])
    H = float(p["half3"] if D == 3 else p["half4"])
    X = sprinkle_poisson(Box(D, [0.0] + [-H] * (D - 1), [T] + [H] * (D - 1)), 1.0, rng)
    N = len(X)
    ptr, idx = links_csr(X, 1.0)
    cptr, cidx = children_csr(ptr, idx, N)
    rs, ts = float(p["r_src"]), float(p["t_src"])
    src = np.flatnonzero((X[:, 0] < ts) & (np.linalg.norm(X[:, 1:], axis=1) < rs))
    s = np.zeros(N)
    s[src] = 1.0
    U0 = np.zeros((N, D))
    U0[:, 0] = 1.0
    d = D - 1
    half = 0.75
    sl, t0 = float(p["slab"]), float(p["t_read"])
    edges = np.arange(t0, T - sl + 1e-9, sl)
    out = {"N": N, "sources": int(src.size), "arms": {}}
    for cfg in CONFIGS:
        N0 = np.zeros((N, D))
        N0[src, 1:] = orient(cfg, X[src, 1:], rs)
        Ph, Nh, core = transport_orient(X, cptr, cidx, s, U0, N0, float(p["core_eps"]))
        held = Ph[:, 0] > 0.0
        rec = {"vals": [], "core": int(core.sum()), "held": int(held.sum())}
        for a, b in zip(edges[:-1], edges[1:]):
            m = held & (X[:, 0] >= a) & (X[:, 0] < b)
            if not m.any():
                continue
            w = Ph[m, 0]
            c = (w[:, None] * X[m, 1:]).sum(axis=0) / w.sum()
            ok_n = m & ~core
            rel = X[ok_n, 1:] - c
            rr = np.linalg.norm(rel, axis=1)
            nsp = Nh[ok_n, 1:]
            nsp = nsp / np.maximum(np.linalg.norm(nsp, axis=1), 1e-300)[:, None]
            for Rs in p["rings"]:
                sel = np.abs(rr - Rs) < half
                v, ok = shell_value(rel[sel], nsp[sel], d)
                rec["vals"].append((0.5 * (a + b), float(Rs), v, ok))
        out["arms"][cfg] = rec
    return out


# ============================================================================================== run

def _frac(vals, expect):
    rd = [v for *_, v, ok in vals if ok]
    allv = len(vals)
    match = sum(1 for v in rd if abs(v - expect) < 1e-6)
    return (match / len(rd) if rd else float("nan")), (len(rd) / allv if allv else 0.0), len(rd)


def run(p, ctx):
    t0 = time.time()
    mode, D = str(p["mode"]), int(p["D"])
    if mode == "static":
        n = int(p["n_static"])
        jobs = [(p, D, derive_seed(ctx.seed, "F3.13a-static", D, k)) for k in range(n)]
        rows = parallel_map(static_task, jobs, ctx.workers, ctx.progress,
                            lambda i, r_, dd, nn: f"[F3.13a] static {dd}/{nn} ({time.time() - t0:.0f}s)")
        lines = [f"static {D - 1}+1: ball radius {p['r_ball']:g}, shells {p['shells']}, {n} sprinklings; "
                 f"value read on each readable shell (the degree of shell -> sphere of orientations)"]
        ok_all, read_all = True, []
        robust = {}
        for cfg in CONFIGS:
            for th in p["thetas"]:
                key = f"{cfg}|{float(th):g}"
                vals = [(R, v, ok) for r_ in rows for (R, v, ok) in r_[key]]
                f, rd, nr = _frac(vals, EXPECT[cfg])
                if float(th) == 0.0:
                    ok_all &= (f == 1.0)
                    if cfg != "double":                    # double: a finer field, readability reported only
                        read_all.append(rd)
                lines.append(f"  {cfg:7s} theta {float(th):.3f}: reads {EXPECT[cfg]:+d} on {f:.3f} of {nr} readable shells "
                             f"(readable {rd:.2f})")
                if cfg == "defect":
                    robust[float(th)] = f
        thr = max([th for th, f in robust.items() if f >= float(p["pass_frac"])], default=0.0)
        rdmin = min(read_all) if read_all else 0.0
        lines.append(f"  reported: the defect reads +1 on >= {p['pass_frac']:.0%} of shells up to theta = {thr:.3f} rad")
        good = ok_all and rdmin >= float(p["read_static"])
        lines.append(f"instrument: every configuration reads its degree at theta = 0: {'OK' if ok_all else 'FAILS'}; "
                     f"readable >= {p['read_static']:.0%}: {'OK' if rdmin >= float(p['read_static']) else 'FAILS'} "
                     f"(min {rdmin:.2f})")
        return Result("PASSED" if good else "INSTRUMENT_FAILED", lines, {"rows": rows, "robust": robust})

    n = int(p["n_flow3"] if D == 3 else p["n_flow4"])
    jobs = [(p, D, derive_seed(ctx.seed, "F3.13a-flow", D, k)) for k in range(n)]
    rows = parallel_map(flow_task, jobs, ctx.workers, ctx.progress,
                        lambda i, r_, dd, nn: f"[F3.13a] flow {dd}/{nn}: N {r_['N']} ({time.time() - t0:.0f}s)")
    T = p["T3"] if D == 3 else p["T4"]
    lines = [f"flow {D - 1}+1: box t [0, {T:g}], ~{np.mean([r_['N'] for r_ in rows]):.0f} events, "
             f"{np.mean([r_['sources'] for r_ in rows]):.0f} sources (r < {p['r_src']:g}, t < {p['t_src']:g}), cold F3.6 "
             f"transport with merging and the flow-weighted merge of n; rings/shells {p['rings']} around the centroid of "
             f"each slab; {n} sprinklings"]
    res = {}
    for cfg in CONFIGS:
        vals = [(t, R, v, ok) for r_ in rows for (t, R, v, ok) in r_["arms"][cfg]["vals"]]
        f, rd, nr = _frac([(t, v, ok) for t, R, v, ok in vals], EXPECT[cfg])
        res[cfg] = (f, rd, nr)
        rdv = [v for t, R, v, ok in vals if ok]
        hist = {}
        for v in rdv:
            k = int(round(v))
            hist[k] = hist.get(k, 0) + 1
        core = np.mean([r_["arms"][cfg]["core"] for r_ in rows])
        lines.append(f"  {cfg:7s}: reads {EXPECT[cfg]:+d} on {f:.3f} of {nr} readable triples (readable {rd:.2f}); "
                     f"values {dict(sorted(hist.items()))}; core events {core:.1f}")
        by_t = {}
        for t, R, v, ok in vals:
            if ok:
                by_t.setdefault(t, []).append(abs(v - EXPECT[cfg]) < 1e-6)
        lines.append("           match by slab: " + " ".join(f"t{t:g}:{np.mean(m):.2f}" for t, m in sorted(by_t.items())))
    pf, ff = float(p["pass_frac"]), float(p["fail_frac"])
    inst = {f"defect readable >= {p['read_flow']:.0%}": res["defect"][1] >= float(p["read_flow"]),
            f"vacuum reads 0 on >= {pf:.0%}": res["vacuum"][0] >= pf,
            f"twist reads 0 on >= {pf:.0%}": res["twist"][0] >= pf}
    lines.append("instrument: " + "; ".join(f"{k}: {'OK' if v else 'FAILS'}" for k, v in inst.items()))
    if not all(inst.values()):
        status = "INSTRUMENT_FAILED"
    elif res["defect"][0] >= pf and res["anti"][0] >= pf:
        status = "PASSED"
    elif res["defect"][0] < ff or res["anti"][0] < ff:
        status = "FAILED"
    else:
        status = "UNDECIDED"
    return Result(status, lines, {"res": res, "rows": rows})
