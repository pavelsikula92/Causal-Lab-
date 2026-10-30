"""
F7.1a -- The capacity candidate: loaded events keep only their slower links. Does that give gamma = 1?
====================================================================================================

Why this test exists (27 Oct 2026, author's decision)
-----------------------------------------------------
gamma = 1 needs the substrate to change BOTH the clocks and space, from one dynamics
(the real problem, docs/findings.md 27 Oct). The first candidate coupling "Pi -> causal
relations" comes from axiom 4 (capacity): an event under load cannot sustain all its
links; it keeps only those whose step is slow enough, so the cones narrow where the load
is. The blind reader K.16b (PASSED) decides what gamma such an order has.

Derivation stated before any run (2+1; the dictionary of K.16)
--------------------------------------------------------------
For the static isotropic weak field ds^2 = -e^{2 phi} dt^2 + e^{-2 psi} dx^2 a substrate of
proper density 1 has coordinate density e^{phi - 2 psi} and light speed e^{phi + psi}. A
rule that only REMOVES links and keeps every event keeps the coordinate density uniform,
so phi - 2 psi = 0: psi = phi / 2, gamma = 1/2 (3+1: phi - 3 psi = 0, gamma = 1/3), for ANY
strength function c(sigma). The narrowing of the cones then fixes phi = (2/3) ln c (the
clocks follow from the count of a diamond, N ~ c^2 T^3 = tau^3). General relativity
(gamma = 1) would need the substrate ALSO to add events near the matter: coordinate density
e^{-phi} in 2+1 (e^{-2 phi} in 3+1) together with the light speed e^{2 phi}.
Prediction of this test therefore: gamma = 1/2, not 1. Estimate A (cones + density) equals
1/2 by construction once the density is uniform; estimate B (clocks + cones) is the real
check that an order built from capacity-limited links is metric-like (its diamonds count
like the cones predict).

Model (locked before any run)
-----------------------------
Flat 2+1 Poisson sprinkling of density 1, box t in [0, 14], |x|, |y| <= 18 (the generator of
K.16 with a = 0), links by core.order.links. Load of a static source, prescribed (the source
of the load is not modelled): sigma(r) = sigma0 exp(-r^2 / 128) (the width of K.16b, inside
the reader's calibrated range). Capacity rule: the link v -> k is kept iff its coordinate
speed |dx|/dt <= c(sigma(r_v)) = exp(-sigma(r_v)), judged at the loaded event v, in the rest
frame of the source (the frame of the matter, F3.6). The relation of the substrate is the
reachability through kept links (a chain of steps no faster than c cannot be faster than c).
Virtual diamond endpoints p, q: p relates directly to every event of the diamond window within
its capacity cone (speed <= c(r_p)), q is reached directly from every event e of the window with
speed <= c(r_e); the rest of the diamond follows by reachability through kept links.
Configurations: flat (sigma0 = 0: all links kept, the instrument control), capacity
(sigma0 = 0.2: |phi| ~ 0.13 at the centre if gamma = 1/2), and sigma0 = 0.3 (strength
independence, reported only). 16 sprinklings each.

Readers: those of K.16b (pooled cone boundary, 100 anchors; diamonds of height 3, 200 per
bin; density), the cones and diamonds reading the reachability above. gamma per sprinkling
by least squares through the origin over the bins (as K.16b).

Criteria (stated before any run)
--------------------------------
INSTRUMENT_FAILED if in the flat configuration some reading of some bin has
|mean| > 0.01 + 3 SE (the link reachability must reproduce the flat order).
UNDECIDED (not metric-like, the reader does not apply) if for sigma0 = 0.2
|gamma_A - gamma_B| - 2 SE > 0.25.
PASSED (the candidate gives general relativity) if |gamma - 1| + 2 SE <= 0.25 for both
estimates; FAILED if |gamma - 1| - 2 SE > 0.25 for both; UNDECIDED otherwise.
Reported: the distance from the predicted 1/2, the read light speed against the nominal
c(r) (discreteness: chains through kept links may reach less than c), sigma0 = 0.3.

Development (quick run, decide nothing; changed before any full run): the first version let the virtual endpoints
step only to events within dt <= 1 and reach the rest through links. Without any load that read phi_clock = -0.53:
an event near the light cone of p is generally not above any event of the first unit of time, so the diamond lost
its rim. The endpoints now relate directly to the whole window (first_dt = T_d), as a real event would. The cone
reader (real anchors) read s = 0.000 without load and the nominal ln c with it. No criterion changed.
"""

from __future__ import annotations

import math
import time

import numpy as np

from core.order.links import links_csr
from core.order.sprinkle import derive_seed
from core.stats.parallel import parallel_map
from core.testing import Result
from phases.K_instruments import K_16_gamma_reader as k16

TEST = {
    "id": "F7.1a",
    "phase": "F7",
    "title": "A candidate from capacity: loaded events keep only slower links. Does it give γ = 1?",
    "theory_part": "why γ = 1: the coupling Π → causal relation from axiom 4 (capacity); measured by the blind "
                   "instrument K.16b",
    "question": "When a loaded event keeps only the links with speed <= c(σ) (the cones narrow, the events stay), "
                "what γ does the blind instrument read? Prediction from the derivation: 1/2 for any function c(σ), "
                "because the density does not change.",
    "method": "flat 2+1 sprinkling (density 1), links; load of a static source σ = σ0·exp(−r²/128), link v → k kept "
              "when |Δx|/Δt <= exp(−σ(r_v)) in the rest frame of the source; relation = reachability through the "
              "kept links; reading K.16b (cones, diamonds, density); σ0 = 0 (control), 0.2, 0.3; 16 sprinklings",
    "prediction": "γ = 1/2 (2+1; in 3+1 1/3) for any c(σ): removing links without adding events does not change the "
                  "density, so φ − 2ψ = 0",
    "criterion": "INSTRUMENT: σ0 = 0 some bin reading |mean| > 0.01 + 3 SE. UNDECIDED (not metric): |γ_A − γ_B| − 2 "
                 "SE > 0.25. PASSED (the candidate gives GR): |γ − 1| + 2 SE <= 0.25 both paths; FAILED: |γ − 1| − 2 "
                 "SE > 0.25 both paths; otherwise UNDECIDED.",
    "null_model": "σ0 = 0 (all links kept: flat order)",
    "level": "mechanism",
    "meters": ["K.16b"],
    "depends_on": ["K.16b"],
    "params": {"T": 14.0, "half": 18.0, "width": 8.0, "bins": [0, 2, 4, 6, 8, 10, 12, 14],
               "configs": [["flat", 0.0], ["capacity", 0.2], ["capacity_strong", 0.3]],
               "n_seeds": 16, "anchors": 100, "diamonds": 200, "T_d": 3.0, "cone_dt": [1.0, 2.0], "cone_q": 1.5,
               "first_dt": 3.0, "tol_gamma": 0.25, "null_abs": 0.01, "null_se": 3.0, "n_se": 2.0},
    "parts": {"main": {}},
    "quick": {"n_seeds": 2, "anchors": 30, "diamonds": 40, "half": 12.0, "bins": [0, 2, 4, 6, 8, 10]},
    "estimate": "~5–15 min with 8–10 processes (the links of 48 sprinklings of ~18 000 events dominate)",
}


# ============================================================================================== the substrate

class CapacityOrder:
    """A flat sprinkling whose links are kept only if their speed is within the capacity cone of the parent event."""

    def __init__(self, X, sigma0, width):
        self.X = X
        self.N = len(X)
        self.sigma0, self.width = float(sigma0), float(width)
        ptr, idx = links_csr(X, 1.0)                              # parents CSR
        child = np.repeat(np.arange(self.N), np.diff(ptr))        # the later event of each link
        parent = idx
        d = X[child] - X[parent]
        speed = np.hypot(d[:, 1], d[:, 2]) / d[:, 0]
        keep = speed <= self.c_at(X[parent, 1:])
        self.n_links, self.n_kept = int(keep.size), int(keep.sum())
        parent, child = parent[keep], child[keep]
        o = np.argsort(parent, kind="stable")                     # children CSR of the kept links
        self.cptr = np.concatenate([[0], np.cumsum(np.bincount(parent, minlength=self.N))]).astype(np.int64)
        self.cidx = child[o]
        o = np.argsort(child, kind="stable")                      # parents CSR of the kept links
        self.pptr = np.concatenate([[0], np.cumsum(np.bincount(child, minlength=self.N))]).astype(np.int64)
        self.pidx = parent[o]

    def c_at(self, P):
        """Light speed allowed by the capacity at spatial points P: exp(-sigma(r))."""
        r2 = (np.asarray(P) ** 2).sum(axis=-1)
        return np.exp(-self.sigma0 * np.exp(-r2 / (2.0 * self.width ** 2)))

    def _bfs(self, start, ptr, idx, inside):
        """Events reachable from the start set along the CSR, restricted to the events where inside is True."""
        seen = np.zeros(self.N, dtype=bool)
        front = np.unique(np.asarray(start, dtype=np.int64))
        front = front[inside[front]]
        seen[front] = True
        while front.size:
            nxt = np.concatenate([idx[ptr[v]:ptr[v + 1]] for v in front]) if front.size else np.zeros(0, np.int64)
            if nxt.size == 0:
                break
            nxt = np.unique(nxt)
            nxt = nxt[inside[nxt] & ~seen[nxt]]
            seen[nxt] = True
            front = nxt
        return seen

    def future(self, start, t_max):
        return self._bfs(start, self.cptr, self.cidx, self.X[:, 0] <= t_max)

    def past(self, start, t_min):
        return self._bfs(start, self.pptr, self.pidx, self.X[:, 0] >= t_min)


# ============================================================================================== readers on the substrate

def read_cones(order, edges, T, n_anchor, dt_rng, qmax, rng):
    """Pooled cone boundary per bin (K.16b), the relation being reachability through kept links."""
    X = order.X
    r = np.hypot(X[:, 1], X[:, 2])
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        idx = np.flatnonzero((r >= lo) & (r < hi) & (X[:, 0] >= 1.0) & (X[:, 0] <= T - 1.0 - dt_rng[1]))
        if idx.size == 0:
            out.append(float("nan"))
            continue
        qs, oks = [], []
        for a in rng.choice(idx, min(n_anchor, idx.size), replace=False):
            d = X - X[a]
            dt = d[:, 0]
            dx = np.hypot(d[:, 1], d[:, 2])
            cand = np.flatnonzero((dt >= dt_rng[0]) & (dt <= dt_rng[1]) & (dx <= qmax * dt))
            if cand.size == 0:
                continue
            fut = order.future(order.cidx[order.cptr[a]:order.cptr[a + 1]], X[a, 0] + dt_rng[1])
            qs.append(dx[cand] / dt[cand])
            oks.append(fut[cand])
        if not qs:
            out.append(float("nan"))
            continue
        q, ok = np.concatenate(qs), np.concatenate(oks)
        o = np.argsort(q)
        q, ok = q[o], ok[o]
        rel_after = ok[::-1].cumsum()[::-1]
        unrel_upto = (~ok).cumsum()
        err = np.append(rel_after[1:], 0) + unrel_upto
        k = int(np.argmin(err))
        c = 0.5 * (q[k] + q[k + 1]) if k + 1 < q.size else q[k]
        out.append(math.log(c))
    return np.array(out)


def read_clocks(order, edges, T, n_diam, T_d, first_dt, rng):
    """phi_b = ln(tau / T_d) from the mean count of diamonds at rest between virtual endpoints (K.16b volume law)."""
    X = order.X
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        rc = 0.5 * (lo + hi)
        counts = []
        for _ in range(n_diam):
            ang = rng.uniform(0, 2 * math.pi)
            x0, y0 = rc * math.cos(ang), rc * math.sin(ang)
            t0 = rng.uniform(1.0, T - 1.0 - T_d)
            t1 = t0 + T_d
            dxp = np.hypot(X[:, 1] - x0, X[:, 2] - y0)
            c_p = float(order.c_at(np.array([[x0, y0]]))[0])
            first = np.flatnonzero((X[:, 0] > t0) & (X[:, 0] <= t0 + first_dt) & (dxp <= c_p * (X[:, 0] - t0)))
            last = np.flatnonzero((X[:, 0] < t1) & (X[:, 0] >= t1 - first_dt)
                                  & (dxp <= order.c_at(X[:, 1:]) * (t1 - X[:, 0])))
            if first.size == 0 or last.size == 0:
                counts.append(0)
                continue
            fut = order.future(first, t1)
            pas = order.past(last, t0)
            counts.append(int((fut & pas).sum()))
        tau = (12.0 * float(np.mean(counts)) / math.pi) ** (1.0 / 3.0)
        out.append(math.log(tau / T_d) if tau > 0 else float("nan"))
    return np.array(out)


def seed_task(p, sigma0, seed):
    rng = np.random.default_rng(seed)
    T, H = float(p["T"]), float(p["half"])
    X = k16.sprinkle(k16.Field(0.0, 0.0, p["width"]), T, H, rng)
    order = CapacityOrder(X, sigma0, p["width"])
    edges = np.array(p["bins"], float)
    d = k16.read_density(X, edges, T)
    s = read_cones(order, edges, T, int(p["anchors"]), p["cone_dt"], float(p["cone_q"]), rng)
    pc = read_clocks(order, edges, T, int(p["diamonds"]), float(p["T_d"]), float(p["first_dt"]), rng)
    phiA, psiA = (2 * s + d) / 3.0, (s - d) / 3.0
    phiB, psiB = pc, s - pc
    return {"N": len(X), "kept": order.n_kept / max(order.n_links, 1), "d": d.tolist(), "s": s.tolist(),
            "phic": pc.tolist(), "gA": k16.gamma_fit(phiA, psiA), "gB": k16.gamma_fit(phiB, psiB),
            "phiA": phiA.tolist(), "psiA": psiA.tolist(), "phiB": phiB.tolist(), "psiB": psiB.tolist()}


def _mean_se(rows, key):
    a = np.array([r_[key] for r_ in rows], float)
    if a.ndim == 1:
        a = a[np.isfinite(a)]
        return float(a.mean()), float(a.std(ddof=1) / math.sqrt(a.size)) if a.size > 1 else float("nan")
    return np.nanmean(a, axis=0), np.nanstd(a, axis=0, ddof=1) / math.sqrt(a.shape[0])


def run(p, ctx):
    t0 = time.time()
    n = int(p["n_seeds"])
    jobs, keys = [], []
    for name, s0 in p["configs"]:
        for k in range(n):
            jobs.append((p, float(s0), derive_seed(ctx.seed, "F7.1a", name, k)))
            keys.append(name)
    res = parallel_map(seed_task, jobs, ctx.workers, ctx.progress,
                       lambda i, r_, dd, nn: f"[F7.1a] {keys[i]} sprinkling {dd}/{nn} ({time.time() - t0:.0f}s)")
    edges = np.array(p["bins"], float)
    centres = 0.5 * (edges[:-1] + edges[1:])
    lines = [f"2+1 flat sprinkling, density 1, t [0, {p['T']:g}], |x|, |y| <= {p['half']:g}; load sigma0 exp(-r^2/(2*{p['width']:g}^2)); "
             f"link kept iff speed <= exp(-sigma) at the parent; relation = reachability through kept links; readers of K.16b; "
             f"{n} sprinklings per configuration (~{np.mean([r_['N'] for r_ in res]):.0f} events)"]
    st = {}
    for name, s0 in p["configs"]:
        rows = [r_ for r_, k in zip(res, keys) if k == name]
        st[name] = {k: _mean_se(rows, k) for k in ("d", "s", "phic", "phiA", "psiA", "phiB", "psiB", "gA", "gB", "kept")}
        x = st[name]
        nominal = -s0 * np.exp(-centres ** 2 / (2 * p["width"] ** 2))       # ln c of the rule (an upper bound of s)
        lines.append(f"  {name} (sigma0 = {s0:g}): links kept {x['kept'][0]:.3f}; gamma A (cones + density) "
                     f"{x['gA'][0]:+.3f} ± {x['gA'][1]:.3f}; gamma B (clocks + cones) {x['gB'][0]:+.3f} ± {x['gB'][1]:.3f}")
        lines.append("        r          " + " ".join(f"{c:7.1f}" for c in centres))
        lines.append("        ln c rule  " + " ".join(f"{v:+7.4f}" for v in nominal))
        lines.append("        s read     " + " ".join(f"{v:+7.4f}" for v in x["s"][0]))
        lines.append("        d read     " + " ".join(f"{v:+7.4f}" for v in x["d"][0]))
        lines.append("        phi clock  " + " ".join(f"{v:+7.4f}" for v in x["phic"][0]))
        lines.append("        phi A      " + " ".join(f"{v:+7.4f}" for v in x["phiA"][0]))
        lines.append("        psi A      " + " ".join(f"{v:+7.4f}" for v in x["psiA"][0]))
        lines.append("        phi B      " + " ".join(f"{v:+7.4f}" for v in x["phiB"][0]))
        lines.append("        psi B      " + " ".join(f"{v:+7.4f}" for v in x["psiB"][0]))
    # instrument: the flat link reachability must read nothing
    inst_ok, worst = True, -math.inf
    for f in ("d", "s", "phic"):
        m, se = st["flat"][f]
        for mb, sb in zip(m, se):
            bad = not (np.isfinite(mb) and np.isfinite(sb)) or abs(mb) > float(p["null_abs"]) + float(p["null_se"]) * sb
            inst_ok &= not bad
            if np.isfinite(mb) and np.isfinite(sb):
                worst = max(worst, abs(mb) - float(p["null_se"]) * sb)
    lines.append(f"  instrument (flat, all links): every reading |mean| <= {p['null_abs']} + {p['null_se']:g} SE -> "
                 f"{'OK' if inst_ok else 'FAILS'} (largest |mean| - {p['null_se']:g} SE = {worst:+.4f})")
    x = st["capacity"]
    (gA, sA), (gB, sB) = x["gA"], x["gB"]
    nse, tol = float(p["n_se"]), float(p["tol_gamma"])
    metric_like = abs(gA - gB) - nse * math.hypot(sA, sB) <= tol
    gr = all(abs(g - 1.0) + nse * s_ <= tol for g, s_ in ((gA, sA), (gB, sB)))
    not_gr = all(abs(g - 1.0) - nse * s_ > tol for g, s_ in ((gA, sA), (gB, sB)))
    lines.append(f"  capacity (sigma0 = 0.2): |gamma_A - gamma_B| = {abs(gA - gB):.3f} ± {math.hypot(sA, sB):.3f} -> "
                 f"{'metric-like' if metric_like else 'NOT metric-like'}; distance from the predicted 1/2: A {gA - 0.5:+.3f}, "
                 f"B {gB - 0.5:+.3f}; from GR (1): A {gA - 1:+.3f}, B {gB - 1:+.3f}")
    xs = st["capacity_strong"]
    lines.append(f"  strength independence (sigma0 = 0.3, reported): gamma A {xs['gA'][0]:+.3f} ± {xs['gA'][1]:.3f}, "
                 f"B {xs['gB'][0]:+.3f} ± {xs['gB'][1]:.3f}")
    if not inst_ok:
        status = "INSTRUMENT_FAILED"
    elif not metric_like:
        status = "UNDECIDED"
    elif gr:
        status = "PASSED"
    elif not_gr:
        status = "FAILED"
    else:
        status = "UNDECIDED"
    data = {name: {k: (list(map(float, v[0])), list(map(float, v[1]))) if isinstance(v[0], np.ndarray) else list(v)
                   for k, v in x_.items()} for name, x_ in st.items()}
    return Result(status, lines, data)
