"""
K.16 -- A blind reader of gamma: clocks, cones and density read from an order, without knowing the metric
=========================================================================================================

Why this instrument exists (27 Oct 2026, author's decision)
-----------------------------------------------------------
The real problem of the theory is WHY gamma = 1: the weak-field response of the
substrate must slow the clocks AND change space, and the equality psi = phi must
come out of one dynamics, not be put in (H1 puts it in). Before any candidate
coupling "Pi -> causal relations" can be tested, the laboratory needs an
instrument that reads phi and psi from an order without being told the metric.

Dictionary (an order plus a count of elements is a metric): for the static weak
field ds^2 = -e^{2 phi} dt^2 + e^{-2 psi} dx^2 in 2+1 (lengths in the discreteness
length, proper density rho = 1),
    clocks   proper time of a diamond at rest       tau = e^{phi} T           -> phi
    cones    coordinate light speed (the relation)  c   = e^{phi + psi}       -> s = phi + psi
    density  elements per coordinate volume        sqrt(-g) = e^{phi - 2 psi} -> d = phi - 2 psi
(3+1: d = phi - 3 psi). Any two determine (phi, psi); the third is a check. gamma = psi/phi:
the conformal geometry of a fixed order (A6 on a sprinkling, the F3 rule) has gamma = -1
(s = 0), H1 and general relativity gamma = +1 (s = 2 phi, d = -phi in 2+1).

Readers (blind: they get the events, the relation of the substrate as a function of
pairs, and the declared proper density rho = 1; never phi, psi or gamma)
---------------------------------------------------------------------------------
Radial bins of width 2 around the source (r < 14), events with 1 <= t <= T - 1.
  density  d_b = ln( events in the bin / coordinate volume ).
  cones    anchors = events of the bin; for each, the events with dt in [1, 2] and
           |dx| <= 1.5 dt; q = |dx| / dt; the relation splits them at the local light
           speed: c = (largest q of a related event + smallest larger q of an unrelated
           one) / 2; s_b = ln(mean c over the anchors of the bin).
  clocks   diamonds at rest at the bin centre (random angle and start), coordinate
           height T_d = 5: the events p < e < q counted with the relation; proper time
           from the count, tau = (12 <N> / (pi rho))^{1/3} (2+1 diamond volume pi tau^3 / 12);
           phi_b = ln(tau / T_d).
Two estimates of gamma per sprinkling, least squares through the origin over the bins:
  A (cones + density):  phi = (2 s + d)/3, psi = (s - d)/3;
  B (clocks + cones):   phi = phi_clock,   psi = s - phi_clock.
Mean and SE over the sprinklings.

Calibration sprinklings (the generator; the readers do not see it)
-------------------------------------------------------------------
2+1 box t in [0, 14], x, y in [-18, 18]; field phi(r) = -a exp(-r^2 / (2 * 8^2)),
psi = gamma phi. Elements: Poisson of proper density 1 (coordinate rate
e^{phi - 2 psi}, by thinning). Relation p < q iff t_q - t_p >= int e^{-(phi + psi)} dl
along the straight segment (core.order.relation.optical_length; the bending of the
path is second order). Configurations: null a = 0; a = 0.1 with gamma = -1, 0, +1;
a = 0.03 with gamma = +1 (resolution, reported only). 16 sprinklings each.

Criteria (stated before any run)
--------------------------------
INSTRUMENT_FAILED (the readers are biased without a field) if in the null configuration
some reading of some bin has |mean| > 0.01 + 3 SE.
PASSED if for every gamma_true in {-1, 0, +1} (a = 0.1) both estimates satisfy
|gamma_est - gamma_true| + 2 SE <= 0.25; FAILED otherwise.
Reported, not decided: phi, psi per bin against the true profile (smoothing of the
field over the diamond and the cone segment), the a = 0.03 resolution, and the
agreement of the two estimates.
Declared limits: the readers use coordinates (a declared readout, audit of axioms) and
the proper density; they serve substrates built on an embedding (every current
candidate of the laboratory). A coordinate-free reader (radar echo, Shapiro delay, as
K.2 does for one geometry) is a later step.

After the full run (27 Oct): the evaluation was moved unchanged into evaluate() so that K.16b can share it; no
logic changed. Result: FAILED (estimate B for gamma = +1: |dev| + 2 SE = 0.253); see K.16b.
"""

from __future__ import annotations

import math
import time

import numpy as np

from core.order.relation import optical_length
from core.order.sprinkle import derive_seed
from core.stats.parallel import parallel_map
from core.testing import Result

TEST = {
    "id": "K.16",
    "phase": "K",
    "title": "Slepé čtení γ: hodiny, kužele a hustota z uspořádání bez znalosti metriky",
    "theory_part": "proč γ = 1: časová a prostorová odezva substrátu musí vyjít z jedné dynamiky; měřidlo pro kandidátní "
                   "vazby Π → kauzální relace",
    "question": "Určí měřidlo z událostí, relace a hustoty (bez φ, ψ a γ) γ = −1, 0 a +1 slabého pole ve 2+1 dvěma "
                "nezávislými cestami (kužele + hustota, hodiny + kužele)?",
    "method": "sprinklingy 2+1 s vlastní hustotou 1 do −e^{2φ}dt² + e^{−2ψ}dx², φ = −a·exp(−r²/128), ψ = γφ; čtení v "
              "radiálních binech: hustota → φ − 2ψ, hranice relace → φ + ψ, počet prvků diamantu → φ; γ jako fit přes biny; "
              "16 sprinklingů na konfiguraci, nulová kontrola a = 0",
    "prediction": "měřidlo vrátí vložené γ (kalibrace)",
    "criterion": "MĚŘIDLO: v nulové konfiguraci některé čtení některého binu |průměr| > 0,01 + 3 SE. PROŠEL: pro γ = −1, 0, +1 "
                 "(a = 0,1) obě cesty |γ_odhad − γ| + 2 SE <= 0,25; NEPROŠEL jinak.",
    "null_model": "a = 0 (bez pole)",
    "level": "calibration",
    "meters": ["K.2", "K.3"],
    "depends_on": [],
    "params": {"T": 14.0, "half": 18.0, "width": 8.0, "bins": [0, 2, 4, 6, 8, 10, 12, 14],
               "configs": [["null", 0.0, 0.0], ["g-1", 0.1, -1.0], ["g0", 0.1, 0.0], ["g+1", 0.1, 1.0],
                           ["g+1_weak", 0.03, 1.0]],
               "n_seeds": 16, "anchors": 30, "diamonds": 30, "T_d": 5.0, "cone_dt": [1.0, 2.0], "cone_q": 1.5,
               "tol_gamma": 0.25, "null_abs": 0.01, "null_se": 3.0, "n_se": 2.0},
    "parts": {"main": {}},
    "quick": {"n_seeds": 3, "anchors": 12, "diamonds": 12},
    "estimate": "~5–10 min s 10 procesy (80 sprinklingů po ~18 000 událostech; čtení je levné)",
}


# ============================================================================================== generator (hidden)

class Field:
    """phi(r) = -a exp(-r^2 / 2 w^2), psi = gamma phi; the generator's metric (never given to the readers)."""

    def __init__(self, a, gamma, w):
        self.a, self.g, self.w = float(a), float(gamma), float(w)

    def phi(self, P):
        r2 = (np.asarray(P) ** 2).sum(axis=-1)
        return -self.a * np.exp(-r2 / (2.0 * self.w ** 2))

    def inv_speed(self, P):                                       # 1 / c = e^{-(phi + psi)}
        return np.exp(-(1.0 + self.g) * self.phi(P))

    def volume_factor(self, P):                                   # sqrt(-g) in 2+1 = e^{phi - 2 psi}
        return np.exp((1.0 - 2.0 * self.g) * self.phi(P))


def sprinkle(field, T, H, rng):
    """Poisson process of proper density 1: coordinate rate sqrt(-g), by thinning; sorted by time."""
    r = np.linspace(0.0, H * math.sqrt(2.0), 400)
    vmax = float(np.max(field.volume_factor(np.stack([r, np.zeros_like(r)], axis=1))))
    n = rng.poisson(vmax * T * (2 * H) ** 2)
    X = np.column_stack([rng.uniform(0, T, n), rng.uniform(-H, H, n), rng.uniform(-H, H, n)])
    keep = rng.random(n) < field.volume_factor(X[:, 1:]) / vmax
    X = X[keep]
    return X[np.argsort(X[:, 0])]


def make_relation(field):
    """The substrate's relation as a function of paired arrays (A_i < B_i), straight-line optical length."""
    def rel(A, B):
        dt = B[:, 0] - A[:, 0]
        return (dt > 0) & (dt >= optical_length(A[:, 1:], B[:, 1:], field.inv_speed))
    return rel


# ============================================================================================== readers (blind)

def read_density(X, edges, T):
    """d_b = ln(events per coordinate volume) in radial bins, events with 1 <= t <= T - 1."""
    r = np.hypot(X[:, 1], X[:, 2])
    m = (X[:, 0] >= 1.0) & (X[:, 0] <= T - 1.0)
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        n = int(((r >= lo) & (r < hi) & m).sum())
        vol = math.pi * (hi * hi - lo * lo) * (T - 2.0)
        out.append(math.log(n / vol) if n > 0 else float("nan"))
    return np.array(out)


def read_cones(X, rel, edges, T, n_anchor, dt_rng, qmax, rng):
    """s_b = ln(local light speed) from the split of related / unrelated events around anchors of each bin."""
    r = np.hypot(X[:, 1], X[:, 2])
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        idx = np.flatnonzero((r >= lo) & (r < hi) & (X[:, 0] >= 1.0) & (X[:, 0] <= T - 1.0 - dt_rng[1]))
        if idx.size == 0:
            out.append(float("nan"))
            continue
        cs = []
        for i in rng.choice(idx, min(n_anchor, idx.size), replace=False):
            d = X - X[i]
            dt = d[:, 0]
            dx = np.hypot(d[:, 1], d[:, 2])
            cand = np.flatnonzero((dt >= dt_rng[0]) & (dt <= dt_rng[1]) & (dx <= qmax * dt))
            if cand.size < 4:
                continue
            q = dx[cand] / dt[cand]
            ok = rel(np.repeat(X[i][None, :], cand.size, axis=0), X[cand])
            if not ok.any() or ok.all():
                continue
            c_lo = q[ok].max()
            above = q[(~ok) & (q > c_lo)]
            if above.size == 0:
                continue
            cs.append(0.5 * (c_lo + above.min()))
        out.append(math.log(float(np.mean(cs))) if cs else float("nan"))
    return np.array(out)


def read_clocks(X, rel, edges, T, n_diam, T_d, rng):
    """phi_b = ln(tau / T_d), tau from the mean count of diamonds at rest at the bin centre (2+1 volume pi tau^3/12)."""
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        rc = 0.5 * (lo + hi)
        counts = []
        for _ in range(n_diam):
            ang = rng.uniform(0, 2 * math.pi)
            x0, y0 = rc * math.cos(ang), rc * math.sin(ang)
            t0 = rng.uniform(1.0, T - 1.0 - T_d)
            p, q = np.array([t0, x0, y0]), np.array([t0 + T_d, x0, y0])
            cand = np.flatnonzero((X[:, 0] > t0) & (X[:, 0] < t0 + T_d)
                                  & (np.hypot(X[:, 1] - x0, X[:, 2] - y0) <= 0.75 * T_d))
            if cand.size == 0:
                counts.append(0)
                continue
            E = X[cand]
            ok = rel(np.repeat(p[None, :], cand.size, axis=0), E) & rel(E, np.repeat(q[None, :], cand.size, axis=0))
            counts.append(int(ok.sum()))
        tau = (12.0 * float(np.mean(counts)) / math.pi) ** (1.0 / 3.0)
        out.append(math.log(tau / T_d) if tau > 0 else float("nan"))
    return np.array(out)


def gamma_fit(phi, psi):
    """Least squares psi = gamma phi through the origin over the finite bins."""
    m = np.isfinite(phi) & np.isfinite(psi)
    den = float((phi[m] ** 2).sum())
    return float((phi[m] * psi[m]).sum() / den) if den > 0 else float("nan")


# ============================================================================================== task

def seed_task(p, a, gamma, seed):
    rng = np.random.default_rng(seed)
    field = Field(a, gamma, p["width"])
    T, H = float(p["T"]), float(p["half"])
    X = sprinkle(field, T, H, rng)
    rel = make_relation(field)
    edges = np.array(p["bins"], float)
    d = read_density(X, edges, T)
    s = read_cones(X, rel, edges, T, int(p["anchors"]), p["cone_dt"], float(p["cone_q"]), rng)
    pc = read_clocks(X, rel, edges, T, int(p["diamonds"]), float(p["T_d"]), rng)
    phiA, psiA = (2 * s + d) / 3.0, (s - d) / 3.0
    phiB, psiB = pc, s - pc
    return {"N": len(X), "d": d.tolist(), "s": s.tolist(), "phic": pc.tolist(),
            "gA": gamma_fit(phiA, psiA), "gB": gamma_fit(phiB, psiB),
            "phiA": phiA.tolist(), "psiA": psiA.tolist(), "phiB": phiB.tolist(), "psiB": psiB.tolist()}


def run(p, ctx):
    t0 = time.time()
    n = int(p["n_seeds"])
    jobs, keys = [], []
    for name, a, g in p["configs"]:
        for k in range(n):
            jobs.append((p, float(a), float(g), derive_seed(ctx.seed, "K.16", name, k)))
            keys.append(name)
    res = parallel_map(seed_task, jobs, ctx.workers, ctx.progress,
                       lambda i, r_, dd, nn: f"[K.16] {keys[i]} sprinkling {dd}/{nn} ({time.time() - t0:.0f}s)")
    return evaluate(p, res, keys, n)


def evaluate(p, res, keys, n):
    """Readings, criteria and verdict from the per-sprinkling results (shared with K.16b; unchanged since K.16 ran)."""
    edges = np.array(p["bins"], float)
    centres = 0.5 * (edges[:-1] + edges[1:])
    lines = [f"2+1, proper density 1, box t [0, {p['T']:g}], |x|, |y| <= {p['half']:g}; phi = -a exp(-r^2/(2*{p['width']:g}^2)), "
             f"psi = gamma phi; {n} sprinklings per configuration (~{np.mean([r_['N'] for r_ in res]):.0f} events); "
             f"bins r = {list(p['bins'])}; readers blind (events, relation, proper density)"]
    by = {}
    for name, a, g in p["configs"]:
        rows = [r_ for r_, k in zip(res, keys) if k == name]
        st = {}
        for f in ("d", "s", "phic", "phiA", "psiA", "phiB", "psiB"):
            arr = np.array([r_[f] for r_ in rows], float)
            st[f] = (np.nanmean(arr, axis=0), np.nanstd(arr, axis=0, ddof=1) / math.sqrt(len(rows)))
        for f in ("gA", "gB"):
            v = np.array([r_[f] for r_ in rows], float)
            v = v[np.isfinite(v)]
            st[f] = (float(v.mean()) if v.size else float("nan"),
                     float(v.std(ddof=1) / math.sqrt(v.size)) if v.size > 1 else float("nan"))
        by[name] = (a, g, st)
        true_phi = -a * np.exp(-centres ** 2 / (2 * p["width"] ** 2))
        lines.append(f"  {name} (a = {a:g}, gamma = {g:+g}): gamma A (cones + density) {st['gA'][0]:+.3f} ± {st['gA'][1]:.3f}; "
                     f"gamma B (clocks + cones) {st['gB'][0]:+.3f} ± {st['gB'][1]:.3f}")
        lines.append("        r       " + " ".join(f"{c:7.1f}" for c in centres))
        lines.append("        phi true" + " ".join(f"{v:+7.4f}" for v in true_phi))
        lines.append("        phi A   " + " ".join(f"{v:+7.4f}" for v in st["phiA"][0]))
        lines.append("        phi B   " + " ".join(f"{v:+7.4f}" for v in st["phiB"][0]))
        lines.append("        psi true" + " ".join(f"{v:+7.4f}" for v in g * true_phi))
        lines.append("        psi A   " + " ".join(f"{v:+7.4f}" for v in st["psiA"][0]))
        lines.append("        psi B   " + " ".join(f"{v:+7.4f}" for v in st["psiB"][0]))
        lines.append("        d, s, phi_clock ± SE (first bin): "
                     f"{st['d'][0][0]:+.4f} ± {st['d'][1][0]:.4f}, {st['s'][0][0]:+.4f} ± {st['s'][1][0]:.4f}, "
                     f"{st['phic'][0][0]:+.4f} ± {st['phic'][1][0]:.4f}")
    # instrument: no field, no reading
    _, _, sn = by["null"]
    worst = 0.0
    inst_ok = True
    for f in ("d", "s", "phic"):
        m, se = sn[f]
        for mb, sb in zip(m, se):
            if not (np.isfinite(mb) and np.isfinite(sb)) or abs(mb) > float(p["null_abs"]) + float(p["null_se"]) * sb:
                inst_ok = False
            worst = max(worst, abs(mb) - float(p["null_se"]) * sb) if np.isfinite(mb) else worst
    lines.append(f"  instrument (null, a = 0): every reading |mean| <= {p['null_abs']} + {p['null_se']:g} SE -> "
                 f"{'OK' if inst_ok else 'FAILS'} (largest |mean| - {p['null_se']:g} SE = {worst:+.4f})")
    tol, nse = float(p["tol_gamma"]), float(p["n_se"])
    ok_all = True
    for name in ("g-1", "g0", "g+1"):
        a, g, st = by[name]
        for f in ("gA", "gB"):
            m, se = st[f]
            ok = np.isfinite(m) and np.isfinite(se) and abs(m - g) + nse * se <= tol
            ok_all &= bool(ok)
        lines.append(f"  {name}: |gamma - {g:+g}| + 2 SE: A {abs(st['gA'][0] - g) + nse * st['gA'][1]:.3f}, "
                     f"B {abs(st['gB'][0] - g) + nse * st['gB'][1]:.3f} (tolerance {tol})")
    status = "INSTRUMENT_FAILED" if not inst_ok else ("PASSED" if ok_all else "FAILED")
    data = {name: {"a": a, "gamma": g, **{f: ([list(map(float, v[0])), list(map(float, v[1]))]
                                             if isinstance(v[0], np.ndarray) else list(v)) for f, v in st.items()}}
            for name, (a, g, st) in by.items()}
    return Result(status, lines, data)
