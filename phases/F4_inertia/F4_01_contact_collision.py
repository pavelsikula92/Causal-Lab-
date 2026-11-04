"""
F4.1 -- Contact collision of two flows: do flows that meet exchange momentum, and how (stick, pass, bounce)?
===========================================================================================================

Why this test exists (29 Oct 2026, the author's choice to start the F4/F5 phases)
---------------------------------------------------------------------------------
Matter in Causal Physics is a flow that carries its 4-velocity as a state (F3.6, Axiom 6); when chains of flow
meet at an event they merge into one with F = |P| (K.15: 4-momentum conserved exactly, rest mass grows at an
inelastic merger). F3.6 showed that a flow with carried direction stays localised and moves with its source; it
did not ask what happens when two such flows MEET. Each event sends its whole flow to ONE child (the child of
least rapidity relative to the carried direction), so two flows interact by contact only if their chains choose
the same events. Whether that happens often (contact interaction, "matter can push matter") or almost never
(flows pass through each other and interact only through the order) was not known.

Development, before the criteria (29 Oct, container, 2 sprinklings each, decides nothing): two balls of flow
(radius 3, emitted in a thin layer) at 0 and moving at 0.6 towards it, offset 8 (2+1) / 6 (3+1). Fraction of
the projectile flow found on events that also carry target flow at the end: 0.94 and 0.42 (2+1), 0.46 and 0.69
(3+1); of the target flow 0.75 / 0.64 and 0.27 / 0.41. Flows meet and merge a lot, in 3+1 too.

What the rule implies (derived before the run)
----------------------------------------------
A merger is irreversible and conserves 4-momentum: two chains with 4-momenta m u_A and m u_B become ONE chain with
P = m (u_A + u_B). There is no elastic channel: a chain cannot come out of a contact with its velocity reversed.
For a symmetric head-on contact in the centre-of-momentum frame (u_A = gamma (1, -v), u_B = gamma (1, +v)) the
merged chain is at rest and its rest mass is 2 m gamma: the kinetic energy becomes rest mass (E = m c^2). So the
alternatives for the bundles are
  stick   merged flow at rest in the CM frame, rest mass up by gamma per merged pair;
  pass    the flows cross without merging and keep +-v;
  bounce  flow reverses its velocity (impossible by the rule; its absence is reported, not judged).

Model (locked)
--------------
Poisson sprinkling rho = 1 (lengths in l), links (core.order.links). Two flows, A and B, each a source tube of
radius r_src moving with its velocity during lab time [0, t_src); every source event emits F = 1 with the
4-velocity of its tube. Configurations on the SAME sprinkling:
  head  B starts at x = -d/2 with +v, A at x = +d/2 with -v (head-on, the lab is the CM frame);
  miss  the same, with A and B displaced to y = +-b/2, b = 2 r_src + 2 (no geometric contact; instrument).
Transport: the carried rule of F3.6 with the mass merge of K.15 (F = |P|); the species fractions of every chain
(the parts of its rest mass that came from A and from B) are carried along for the readout only.

Readout (lab time window after the contact, flux-weighted over the carrying events)
-----------------------------------------------------------------------------------
  f_mix   the fraction of the flow on events that carry both species (contact happened);
  v_mix   the velocity P_x / P_0 of the mixed flow (stick: 0);
  v_A, v_B the velocities of the pure flows (pass: -v, +v);
  bounce  the fraction of pure flow moving against its initial direction;
  mass    on mixed events, the rest mass m against the sum of the species parts f_A + f_B: the prediction of the
          rule for parts that kept their directions is gamma sqrt((f_A + f_B)^2 - v^2 (f_A - f_B)^2) (instrument:
          4-momentum conservation, built in, checked to 1e-9); the mean ratio m / (f_A + f_B) is the energy that
          became rest mass (reported).

Criteria (stated before any full run)
-------------------------------------
Instrument (INSTRUMENT_FAILED if any fails): (i) miss: f_mix <= 0.05 on average (flows merge only on contact);
(ii) the head configuration computed again with all coordinates and velocities boosted by beta' = 0.6 chooses the
same child at every event; (iii) 4-momentum conservation on mixed events, |m - m_pred| / m <= 1e-9.
head, per part (quarter rule on f_mix over sprinklings, bootstrap SE):
  PASSED     f_mix - 2 SE > 1/4 and |v_mix| + 2 SE <= v/4: flows that meet stick: contact transfers momentum and the
             merged flow comes to rest in the CM frame;
  FAILED     f_mix + 2 SE < 1/4: flows pass through each other; contact does not transfer momentum;
  UNDECIDED  otherwise.
Reading declared now. PASSED shows that contact between flows is a perfectly inelastic collision at the level of
chains: momentum is shared, kinetic energy becomes rest mass, and the merged flow's inertia is its energy. It
also means that matter in this model cannot scatter elastically by contact: an elastic channel would have to come
from fields or from a bound structure, which the model does not have (open). FAILED would mean that flows are
transparent to each other and every interaction of matter goes through the order (gravity only).

Parts
-----
d3  2+1, box t [0, 24], x [-20, 20], y [-10, 10]; r_src = 3, t_src = 2, d = 8, v = 0.6; window t [18, 22];
    32 sprinklings.
d4  3+1, box t [0, 16], x [-13, 13], y, z [-7, 7]; r_src = 3, t_src = 2, d = 6, v = 0.6; window t [13, 15];
    16 sprinklings.
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
    "id": "F4.1",
    "phase": "F4",
    "title": "Contact collision of two flows: do flows that meet exchange momentum (stick, pass through, bounce)?",
    "theory_part": "Axiom 6 (flows carry momentum; merging F = |P|); inertia and E = mc^2 of merged flow; matter interacting by contact",
    "question": "When two flows with carried 4-velocity meet head-on, do their chains merge (contact transfers momentum, the merged "
                "flow comes to rest in the CM frame and gains rest mass), or do the flows pass through each other?",
    "method": "Poisson sprinkling (2+1 and 3+1, rho = 1), two source tubes moving at +-0.6 towards each other, carried rule of F3.6 "
              "with the mass merge of K.15, species fractions tracked; fraction of mixed flow, velocities and rest mass after the "
              "contact; null: the same flows displaced so that they do not touch",
    "prediction": "the rule has no elastic channel; development runs show substantial merging (2+1 and 3+1)",
    "criterion": "INSTRUMENT: miss f_mix <= 0.05; boosted recomputation chooses the same children; 4-momentum conservation on mixed "
                 "events <= 1e-9. head: PASSED f_mix - 2 SE > 1/4 and |v_mix| + 2 SE <= v/4; FAILED f_mix + 2 SE < 1/4; "
                 "UNDECIDED otherwise.",
    "null_model": "flows displaced by b = 2 r_src + 2 (no contact): no merging",
    "level": "mechanism",
    "meters": ["K.15"],
    "depends_on": ["F3.6", "K.15"],
    "params": {"v": 0.6, "r_src": 3.0, "t_src": 2.0, "beta_check": 0.6, "n_se": 2.0, "boot_B": 4000,
               "tol_miss": 0.05, "tol_cons": 1e-9},
    "parts": {"d3": {"D": 3, "T": 24.0, "hx": 20.0, "hy": 10.0, "d": 8.0, "win": [18.0, 22.0], "n_seeds": 32},
              "d4": {"D": 4, "T": 16.0, "hx": 13.0, "hy": 7.0, "d": 6.0, "win": [13.0, 15.0], "n_seeds": 16}},
    "quick": {"n_seeds": 3, "boot_B": 200},
    "estimate": "d3 ~2-4 min, d4 ~30-60 min with 10 processes (links of ~80 000 events dominate)",
}


def mink(a, b):
    """Minkowski product (+ - - ...), broadcast over leading axes."""
    return a[..., 0] * b[..., 0] - (a[..., 1:] * b[..., 1:]).sum(axis=-1)


def boost(X, beta):
    """Coordinates (or 4-vectors) in the frame moving with beta along x."""
    g = 1.0 / math.sqrt(1.0 - beta * beta)
    Y = np.array(X, float, copy=True)
    Y[..., 0] = g * (X[..., 0] - beta * X[..., 1])
    Y[..., 1] = g * (X[..., 1] - beta * X[..., 0])
    return Y


def transport(X, cptr, cidx, FA, FB, uA, uB):
    """
    Carried rule (F3.6) with the mass merge (K.15): every event sends its whole flow to the child of least rapidity
    relative to its direction; at an event the incoming 4-momenta add and the flow is |P|. fA, fB: the rest-mass
    parts that came from the sources of A and B (bookkeeping for the readout, not used by the rule).
    """
    N, D = X.shape
    inP = np.zeros((N, D))
    inA = np.zeros(N)
    inB = np.zeros(N)
    flux = np.zeros(N)
    P_ev = np.zeros((N, D))
    fa = np.zeros(N)
    fb = np.zeros(N)
    choice = np.full(N, -1, dtype=np.int64)
    for k in range(N):                                   # increasing time = a linear extension of the order
        a_ = inA[k] + FA[k]
        b_ = inB[k] + FB[k]
        if a_ + b_ <= 0.0:
            continue
        P = inP[k] + FA[k] * uA + FB[k] * uB
        m = math.sqrt(max(float(mink(P, P)), 1e-300))
        u = P / m
        flux[k], P_ev[k], fa[k], fb[k] = m, P, a_, b_
        s, e = cptr[k], cptr[k + 1]
        if e == s:
            continue                                     # no child: the flow leaves the sprinkled region
        ks = cidx[s:e]
        B = X[ks] - X[k]
        tau = np.sqrt(np.maximum(mink(B, B), 1e-300))
        j = int(np.argmin(mink(u[None, :], B) / tau))
        c = int(ks[j])
        choice[k] = c
        inP[c] += P
        inA[c] += a_
        inB[c] += b_
    return {"flux": flux, "P": P_ev, "fa": fa, "fb": fb, "choice": choice}


def sources(X, p, conf):
    """Source events of A and B and their 4-velocities (module notes)."""
    D = X.shape[1]
    v, r, ts, d = float(p["v"]), float(p["r_src"]), float(p["t_src"]), float(p["d"])
    g = 1.0 / math.sqrt(1.0 - v * v)
    b = (2.0 * r + 2.0) if conf == "miss" else 0.0
    t = X[:, 0]

    def tube(x0, y0, vx):
        dx = X[:, 1] - (x0 + vx * t)
        dy = X[:, 2] - y0
        rr = dx * dx + dy * dy + ((X[:, 3:] ** 2).sum(axis=1) if D > 3 else 0.0)
        return ((t >= 0.0) & (t < ts) & (rr <= r * r)).astype(float)
    FA = tube(+d / 2.0, +b / 2.0, -v)
    FB = tube(-d / 2.0, -b / 2.0, +v)
    uA = np.zeros(D); uA[0], uA[1] = g, -g * v
    uB = np.zeros(D); uB[0], uB[1] = g, +g * v
    return FA, FB, uA, uB


def readout(X, r, p, uA, uB):
    """Window observables (module notes)."""
    v = float(p["v"])
    g = 1.0 / math.sqrt(1.0 - v * v)
    lo, hi = p["win"]
    car = (r["flux"] > 0.0) & (X[:, 0] >= lo) & (X[:, 0] < hi)
    fa, fb, F, P = r["fa"], r["fb"], r["flux"], r["P"]
    mix = car & (fa > 0.0) & (fb > 0.0)
    pa = car & (fa > 0.0) & (fb == 0.0)
    pb = car & (fb > 0.0) & (fa == 0.0)
    tot = float(F[car].sum())
    vel = lambda m: float(P[m, 1].sum() / P[m, 0].sum()) if m.any() else float("nan")
    out = {"f_mix": float(F[mix].sum()) / tot if tot > 0 else float("nan"),
           "v_mix": vel(mix), "v_A": vel(pa), "v_B": vel(pb), "n_car": int(car.sum()), "flux": tot}
    vx = np.where(car, P[:, 1] / np.maximum(P[:, 0], 1e-300), 0.0)
    pure = float(F[pa].sum() + F[pb].sum())
    out["bounce"] = float(F[pa & (vx > 0.0)].sum() + F[pb & (vx < 0.0)].sum()) / pure if pure > 0 else float("nan")
    if mix.any():
        s_, dlt = fa[mix] + fb[mix], fa[mix] - fb[mix]
        m_pred = g * np.sqrt(np.maximum(s_ * s_ - v * v * dlt * dlt, 0.0))
        out["cons"] = float(np.max(np.abs(F[mix] - m_pred) / F[mix]))
        out["mass_ratio"] = float(F[mix].sum() / s_.sum())
    else:
        out["cons"], out["mass_ratio"] = 0.0, float("nan")
    return out


def seed_task(p, seed):
    """One sprinkling: head and miss configurations, the boosted recomputation of head."""
    D, T, hx, hy = int(p["D"]), float(p["T"]), float(p["hx"]), float(p["hy"])
    rng = np.random.default_rng(seed)
    X = sprinkle_poisson(Box(D, [0.0, -hx] + [-hy] * (D - 2), [T, hx] + [hy] * (D - 2)), 1.0, rng)
    N = len(X)
    ptr, idx = links_csr(X, 1.0)
    child = np.repeat(np.arange(N), np.diff(ptr))
    o = np.argsort(idx, kind="stable")
    cptr = np.concatenate([[0], np.cumsum(np.bincount(idx, minlength=N))]).astype(np.int64)
    cidx = child[o]
    out = {"N": N}
    for conf in ("head", "miss"):
        FA, FB, uA, uB = sources(X, p, conf)
        r = transport(X, cptr, cidx, FA, FB, uA, uB)
        out[conf] = readout(X, r, p, uA, uB)
        out[conf]["sources"] = [int(FA.sum()), int(FB.sum())]
        if conf == "head":
            bc = float(p["beta_check"])
            rb = transport(boost(X, bc), cptr, cidx, FA, FB, boost(uA[None, :], bc)[0], boost(uB[None, :], bc)[0])
            out["same_choice"] = bool(np.array_equal(r["choice"], rb["choice"]))
    return out


def run(p, ctx):
    t0 = time.time()
    n, B, nse = int(p["n_seeds"]), int(p["boot_B"]), float(p["n_se"])
    seeds = [derive_seed(ctx.seed, "F4.1", ctx.part, k) for k in range(n)]
    rows = parallel_map(seed_task, [(p, sd) for sd in seeds], ctx.workers, ctx.progress,
                        lambda i, r_, d, nn: f"[F4.1] {ctx.part} sprinkling {d}/{nn}: f_mix head {r_['head']['f_mix']:.3f} "
                                             f"miss {r_['miss']['f_mix']:.3f} ({time.time() - t0:.0f}s)")
    v = float(p["v"])
    g = 1.0 / math.sqrt(1.0 - v * v)
    bs = lambda *tag: derive_seed(ctx.seed, "F4.1-boot", ctx.part, *tag)
    fin = lambda xs: [x for x in xs if np.isfinite(x)]
    stats = {}
    for conf in ("head", "miss"):
        stats[conf] = {k: boot_se(fin([r_[conf][k] for r_ in rows]), B, bs(conf, k))[:2]
                       for k in ("f_mix", "v_mix", "v_A", "v_B", "bounce", "mass_ratio")}
    lines = [f"{int(p['D']) - 1}+1 sprinkling rho = 1, box t [0, {p['T']:g}], x [-{p['hx']:g}, {p['hx']:g}], other axes "
             f"[-{p['hy']:g}, {p['hy']:g}]; tubes radius {p['r_src']:g}, emission t [0, {p['t_src']:g}), d = {p['d']:g}, "
             f"v = +-{v:g}; {n} sprinklings (~{np.mean([r_['N'] for r_ in rows]):.0f} events, sources "
             f"{np.mean([r_['head']['sources'][0] for r_ in rows]):.0f} + {np.mean([r_['head']['sources'][1] for r_ in rows]):.0f}); "
             f"window t [{p['win'][0]:g}, {p['win'][1]:g}]"]
    for conf in ("head", "miss"):
        s_ = stats[conf]
        lines.append(f"  {conf:4s}: f_mix {s_['f_mix'][0]:.4f} +- {s_['f_mix'][1]:.4f}; v_mix {s_['v_mix'][0]:+.4f} +- "
                     f"{s_['v_mix'][1]:.4f}; v_A {s_['v_A'][0]:+.4f} +- {s_['v_A'][1]:.4f}; v_B {s_['v_B'][0]:+.4f} +- "
                     f"{s_['v_B'][1]:.4f}; bounce {s_['bounce'][0]:.4f}; mixed rest mass / parts {s_['mass_ratio'][0]:.4f} "
                     f"(gamma = {g:.4f} for equal parts)")
    cons = max(max(r_["head"]["cons"], r_["miss"]["cons"]) for r_ in rows)
    inst = {f"miss f_mix <= {p['tol_miss']:g}": stats["miss"]["f_mix"][0] <= float(p["tol_miss"]),
            "boosted recomputation chooses the same children": all(r_["same_choice"] for r_ in rows),
            f"4-momentum conservation on mixed events (max {cons:.1e})": cons <= float(p["tol_cons"])}
    lines.append("instrument: " + "; ".join(f"{k}: {'OK' if ok else 'FAILS'}" for k, ok in inst.items()))
    f, fse = stats["head"]["f_mix"]
    vm, vse = stats["head"]["v_mix"]
    if not all(inst.values()):
        status = "INSTRUMENT_FAILED"
    elif f - nse * fse > 0.25 and np.isfinite(vm) and abs(vm) + nse * vse <= v / 4.0:
        status = "PASSED"
    elif f + nse * fse < 0.25:
        status = "FAILED"
    else:
        status = "UNDECIDED"
    lines.append(f"  head: f_mix - 2 SE = {f - nse * fse:.4f} (> 1/4 stick), f_mix + 2 SE = {f + nse * fse:.4f} (< 1/4 pass); "
                 f"|v_mix| + 2 SE = {abs(vm) + nse * vse:.4f} vs v/4 = {v / 4:.3f}")
    data = {"stats": {c: {k: list(x) for k, x in s_.items()} for c, s_ in stats.items()}, "instrument": inst,
            "rows": rows}
    return Result(status, lines, data)
