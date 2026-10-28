"""
K.16c -- The blind gamma reader in 3+1: clocks, cones and density read from an order, in a box small enough to run
==================================================================================================================

Why this instrument exists (27 Oct 2026, author's decision)
-----------------------------------------------------------
The principle of conserved causal throughput predicts gamma = 2/(d_s - 1): 2 in 2+1 (F7.1b: ~1.85 in two runs, the
alternatives 1/2 and 1 far away) and 1 in 3+1, the dimension of our world. F7.1c tests it in 3+1; it needs the
reader of K.16b in 3+1, calibrated on exact metrics in a box that 3+1 can afford.

Dictionary in 3+1 (ds^2 = -e^{2 phi} dt^2 + e^{-2 psi} dx^2, proper density 1)
------------------------------------------------------------------------------
    clocks   diamond at rest, coordinate height T_d:  tau = e^{phi} T_d, volume pi tau^4 / 24 -> phi
    cones    light speed (the relation)              c = e^{phi + psi}                     -> s = phi + psi
    density  events per coordinate volume            sqrt(-g) = e^{phi - 3 psi}             -> d = phi - 3 psi
Estimates of (phi, psi) per radial bin:
    A (cones + density):  phi = (3 s + d)/4,  psi = (s - d)/4;
    B (clocks + cones):   phi = phi_clock,    psi = s - phi_clock.
gamma per sprinkling: least squares psi = gamma phi through the origin over the bins (as K.16b); mean and SE.
Alternatives F7.1c must tell apart in 3+1: fixed order -1, narrowing cones only 1/3, conserved throughput 1.

Readers (as K.16b, in three spatial dimensions; blind: events, relation, proper density 1)
-----------------------------------------------------------------------------------------
Radial shells of width 1.5 (r < 7.5), events with 1 <= t <= T - 1.
  density  ln(events per coordinate volume of the shell);
  cones    100 anchors per shell, events with dt in [1, 2] and |dx| <= 1.5 dt, the boundary c that misclassifies
           the fewest pooled (q = |dx|/dt, related) pairs; s = ln c;
  clocks   200 diamonds at rest per shell (random direction and start), height T_d = 3, events p < e < q counted
           with the relation; tau = (24 <N> / pi)^{1/4}; phi = ln(tau / T_d).

Calibration sprinklings (hidden generator)
------------------------------------------
3+1 box t in [0, 9], |x|, |y|, |z| <= 9 (~52 000 events); phi(r) = -a exp(-r^2 / (2 * 4^2)), psi = gamma phi;
elements of proper density 1 (coordinate rate e^{phi - 3 psi}, thinning); relation t_q - t_p >= int e^{-(phi + psi)} dl
along the straight segment. Configurations: null (a = 0), a = 0.1 with gamma = -1, 1/3, +1; 16 sprinklings each.
(gamma = 1/3 replaces 0 of K.16b: it is the 3+1 value of the cones-only alternative that F7.1c must exclude.)

Criteria (stated before any run; those of K.16b)
------------------------------------------------
INSTRUMENT_FAILED if in the null configuration some reading of some shell has |mean| > 0.01 + 3 SE.
PASSED if for gamma_true = -1, 1/3, +1 both estimates satisfy |gamma_est - gamma_true| + 2 SE <= 0.25; FAILED
otherwise. Reported: the per-shell profiles against the truth.
"""

from __future__ import annotations

import math
import time

import numpy as np

from core.order.relation import optical_length
from core.order.sprinkle import derive_seed
from core.stats.parallel import parallel_map
from core.testing import Result
from phases.K_meridla import K_16_gamma_reader as k16

TEST = {
    "id": "K.16c",
    "phase": "K",
    "title": "Slepé čtení γ ve 3+1: hodiny, kužele a hustota v boxu, který 3+1 unese",
    "theory_part": "proč γ = 1: měřidlo pro test principu propustnosti ve 3+1 (F7.1c)",
    "question": "Určí měřidlo K.16b převedené do 3+1 (hustota → φ − 3ψ, diamant πτ⁴/24) γ = −1, 1/3 a +1 slabého pole "
                "naslepo oběma cestami v boxu t ∈ [0, 9], |x| <= 9?",
    "method": "sprinklingy 3+1 s vlastní hustotou 1 do −e^{2φ}dt² + e^{−2ψ}dx², φ = −a·exp(−r²/32), ψ = γφ; čtení "
              "v radiálních slupkách šířky 1,5: hustota, sdružená hranice relace, diamanty výšky 3; 16 sprinklingů na "
              "konfiguraci, nulová kontrola a = 0",
    "prediction": "měřidlo vrátí vložené γ (kalibrace)",
    "criterion": "MĚŘIDLO: v nulové konfiguraci některé čtení slupky |průměr| > 0,01 + 3 SE. PROŠEL: pro γ = −1, 1/3, +1 "
                 "obě cesty |γ_odhad − γ| + 2 SE <= 0,25; NEPROŠEL jinak.",
    "null_model": "a = 0 (bez pole)",
    "level": "calibration",
    "meters": ["K.16b"],
    "depends_on": ["K.16b"],
    "params": {"T": 9.0, "half": 9.0, "width": 4.0, "bins": [0, 1.5, 3.0, 4.5, 6.0, 7.5],
               "configs": [["null", 0.0, 0.0], ["g-1", 0.1, -1.0], ["g1/3", 0.1, 1.0 / 3.0], ["g+1", 0.1, 1.0]],
               "n_seeds": 16, "anchors": 100, "diamonds": 200, "T_d": 3.0, "cone_dt": [1.0, 2.0], "cone_q": 1.5,
               "tol_gamma": 0.25, "null_abs": 0.01, "null_se": 3.0, "n_se": 2.0},
    "parts": {"main": {}},
    "quick": {"n_seeds": 2, "anchors": 30, "diamonds": 40, "half": 7.0, "T": 8.0, "bins": [0, 1.5, 3.0, 4.5, 6.0]},
    "estimate": "~10–30 min s 8–10 procesy (64 sprinklingů po ~52 000 událostech; čtení bez vazeb)",
}


# ============================================================================================== generator (hidden)

class Field3:
    """phi(r) = -a exp(-r^2 / 2 w^2), psi = gamma phi, 3+1."""

    def __init__(self, a, gamma, w):
        self.a, self.g, self.w = float(a), float(gamma), float(w)

    def phi(self, P):
        return -self.a * np.exp(-(np.asarray(P) ** 2).sum(axis=-1) / (2.0 * self.w ** 2))

    def inv_speed(self, P):
        return np.exp(-(1.0 + self.g) * self.phi(P))

    def volume_factor(self, P):                                   # sqrt(-g) = e^{phi - 3 psi}
        return np.exp((1.0 - 3.0 * self.g) * self.phi(P))


def sprinkle3(field, T, H, rng):
    r = np.linspace(0.0, H * math.sqrt(3.0), 400)
    vmax = float(np.max(field.volume_factor(np.stack([r, np.zeros_like(r), np.zeros_like(r)], axis=1))))
    n = rng.poisson(vmax * T * (2 * H) ** 3)
    X = np.column_stack([rng.uniform(0, T, n)] + [rng.uniform(-H, H, n) for _ in range(3)])
    X = X[rng.random(n) < field.volume_factor(X[:, 1:]) / vmax]
    return X[np.argsort(X[:, 0])]


def make_relation(field):
    def rel(A, B):
        dt = B[:, 0] - A[:, 0]
        return (dt > 0) & (dt >= optical_length(A[:, 1:], B[:, 1:], field.inv_speed))
    return rel


# ============================================================================================== readers (blind, 3+1)

def read_density3(X, edges, T):
    r = np.sqrt((X[:, 1:] ** 2).sum(axis=1))
    m = (X[:, 0] >= 1.0) & (X[:, 0] <= T - 1.0)
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        n = int(((r >= lo) & (r < hi) & m).sum())
        vol = 4.0 / 3.0 * math.pi * (hi ** 3 - lo ** 3) * (T - 2.0)
        out.append(math.log(n / vol) if n > 0 else float("nan"))
    return np.array(out)


def pooled_boundary(q, ok):
    """The c that misclassifies the fewest (q, related) pairs (K.16b)."""
    o = np.argsort(q)
    q, ok = q[o], ok[o]
    rel_after = ok[::-1].cumsum()[::-1]
    unrel_upto = (~ok).cumsum()
    err = np.append(rel_after[1:], 0) + unrel_upto
    k = int(np.argmin(err))
    return 0.5 * (q[k] + q[k + 1]) if k + 1 < q.size else q[k]


def read_cones3(X, related_from, edges, T, n_anchor, dt_rng, qmax, rng):
    """related_from(i, cand) -> bool array: is event cand related to (in the future of) anchor i."""
    r = np.sqrt((X[:, 1:] ** 2).sum(axis=1))
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        idx = np.flatnonzero((r >= lo) & (r < hi) & (X[:, 0] >= 1.0) & (X[:, 0] <= T - 1.0 - dt_rng[1]))
        if idx.size == 0:
            out.append(float("nan"))
            continue
        qs, oks = [], []
        for i in rng.choice(idx, min(n_anchor, idx.size), replace=False):
            d = X - X[i]
            dt = d[:, 0]
            dx = np.sqrt((d[:, 1:] ** 2).sum(axis=1))
            cand = np.flatnonzero((dt >= dt_rng[0]) & (dt <= dt_rng[1]) & (dx <= qmax * dt))
            if cand.size == 0:
                continue
            qs.append(dx[cand] / dt[cand])
            oks.append(related_from(i, cand))
        out.append(math.log(pooled_boundary(np.concatenate(qs), np.concatenate(oks))) if qs else float("nan"))
    return np.array(out)


def read_clocks3(X, diamond_count, edges, T, n_diam, T_d, rng):
    """diamond_count(p, q) -> number of events p < e < q; phi = ln(tau / T_d), tau = (24 <N> / pi)^{1/4}."""
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        rc = 0.5 * (lo + hi)
        counts = []
        for _ in range(n_diam):
            u = rng.normal(size=3)
            u /= np.linalg.norm(u)
            x0 = rc * u
            t0 = rng.uniform(1.0, T - 1.0 - T_d)
            counts.append(diamond_count(np.concatenate([[t0], x0]), np.concatenate([[t0 + T_d], x0])))
        tau = (24.0 * float(np.mean(counts)) / math.pi) ** 0.25
        out.append(math.log(tau / T_d) if tau > 0 else float("nan"))
    return np.array(out)


def estimates3(d, s, pc):
    phiA, psiA = (3 * s + d) / 4.0, (s - d) / 4.0
    phiB, psiB = pc, s - pc
    return {"d": d.tolist(), "s": s.tolist(), "phic": pc.tolist(),
            "gA": k16.gamma_fit(phiA, psiA), "gB": k16.gamma_fit(phiB, psiB),
            "phiA": phiA.tolist(), "psiA": psiA.tolist(), "phiB": phiB.tolist(), "psiB": psiB.tolist()}


# ============================================================================================== task

def seed_task(p, a, gamma, seed):
    rng = np.random.default_rng(seed)
    field = Field3(a, gamma, p["width"])
    T, H = float(p["T"]), float(p["half"])
    X = sprinkle3(field, T, H, rng)
    rel = make_relation(field)
    edges = np.array(p["bins"], float)
    T_d = float(p["T_d"])

    def related_from(i, cand):
        return rel(np.repeat(X[i][None, :], cand.size, axis=0), X[cand])

    def diamond_count(pp, qq):
        dxp = np.sqrt(((X[:, 1:] - pp[1:]) ** 2).sum(axis=1))
        cand = np.flatnonzero((X[:, 0] > pp[0]) & (X[:, 0] < qq[0]) & (dxp <= 0.75 * T_d))
        if cand.size == 0:
            return 0
        E = X[cand]
        ok = rel(np.repeat(pp[None, :], cand.size, axis=0), E) & rel(E, np.repeat(qq[None, :], cand.size, axis=0))
        return int(ok.sum())

    d = read_density3(X, edges, T)
    s = read_cones3(X, related_from, edges, T, int(p["anchors"]), p["cone_dt"], float(p["cone_q"]), rng)
    pc = read_clocks3(X, diamond_count, edges, T, int(p["diamonds"]), T_d, rng)
    out = estimates3(d, s, pc)
    out["N"] = len(X)
    return out


def mean_se(rows, key):
    a = np.array([r_[key] for r_ in rows], float)
    if a.ndim == 1:
        a = a[np.isfinite(a)]
        return float(a.mean()), float(a.std(ddof=1) / math.sqrt(a.size)) if a.size > 1 else float("nan")
    return np.nanmean(a, axis=0), np.nanstd(a, axis=0, ddof=1) / math.sqrt(a.shape[0])


def null_check(st, p):
    """True if every reading of every shell of the null configuration is compatible with zero (K.16b rule)."""
    ok, worst = True, -math.inf
    for f in ("d", "s", "phic"):
        m, se = st[f]
        for mb, sb in zip(m, se):
            bad = not (np.isfinite(mb) and np.isfinite(sb)) or abs(mb) > float(p["null_abs"]) + float(p["null_se"]) * sb
            ok &= not bad
            if np.isfinite(mb) and np.isfinite(sb):
                worst = max(worst, abs(mb) - float(p["null_se"]) * sb)
    return ok, worst


def run(p, ctx):
    t0 = time.time()
    n = int(p["n_seeds"])
    jobs, keys = [], []
    for name, a, g in p["configs"]:
        for k in range(n):
            jobs.append((p, float(a), float(g), derive_seed(ctx.seed, "K.16c", name, k)))
            keys.append(name)
    res = parallel_map(seed_task, jobs, ctx.workers, ctx.progress,
                       lambda i, r_, dd, nn: f"[K.16c] {keys[i]} sprinkling {dd}/{nn} ({time.time() - t0:.0f}s)")
    edges = np.array(p["bins"], float)
    centres = 0.5 * (edges[:-1] + edges[1:])
    lines = [f"3+1, proper density 1, box t [0, {p['T']:g}], |x|, |y|, |z| <= {p['half']:g}; phi = -a exp(-r^2/(2*{p['width']:g}^2)), "
             f"psi = gamma phi; {n} sprinklings per configuration (~{np.mean([r_['N'] for r_ in res]):.0f} events); "
             f"shells r = {list(p['bins'])}; readers blind"]
    st = {}
    for name, a, g in p["configs"]:
        rows = [r_ for r_, k in zip(res, keys) if k == name]
        st[name] = {k: mean_se(rows, k) for k in ("d", "s", "phic", "phiA", "psiA", "phiB", "psiB", "gA", "gB")}
        x = st[name]
        tp = -a * np.exp(-centres ** 2 / (2 * p["width"] ** 2))
        lines.append(f"  {name} (a = {a:g}, gamma = {g:+.3f}): gamma A (cones + density) {x['gA'][0]:+.3f} ± {x['gA'][1]:.3f}; "
                     f"gamma B (clocks + cones) {x['gB'][0]:+.3f} ± {x['gB'][1]:.3f}")
        lines.append("        r        " + " ".join(f"{c:7.2f}" for c in centres))
        lines.append("        phi true " + " ".join(f"{v:+7.4f}" for v in tp))
        lines.append("        phi A    " + " ".join(f"{v:+7.4f}" for v in x["phiA"][0]))
        lines.append("        phi B    " + " ".join(f"{v:+7.4f}" for v in x["phiB"][0]))
        lines.append("        psi true " + " ".join(f"{v:+7.4f}" for v in g * tp))
        lines.append("        psi A    " + " ".join(f"{v:+7.4f}" for v in x["psiA"][0]))
        lines.append("        psi B    " + " ".join(f"{v:+7.4f}" for v in x["psiB"][0]))
    inst_ok, worst = null_check(st["null"], p)
    lines.append(f"  instrument (null): every reading |mean| <= {p['null_abs']} + {p['null_se']:g} SE -> "
                 f"{'OK' if inst_ok else 'FAILS'} (largest |mean| - {p['null_se']:g} SE = {worst:+.4f})")
    tol, nse = float(p["tol_gamma"]), float(p["n_se"])
    ok_all = True
    for name, a, g in p["configs"][1:]:
        x = st[name]
        vals = []
        for e in ("gA", "gB"):
            m, se = x[e]
            v = abs(m - g) + nse * se if np.isfinite(m) and np.isfinite(se) else float("inf")
            ok_all &= v <= tol
            vals.append(v)
        lines.append(f"  {name}: |gamma - {g:+.3f}| + 2 SE: A {vals[0]:.3f}, B {vals[1]:.3f} (tolerance {tol})")
    status = "INSTRUMENT_FAILED" if not inst_ok else ("PASSED" if ok_all else "FAILED")
    data = {name: {k: (list(map(float, v[0])), list(map(float, v[1]))) if isinstance(v[0], np.ndarray) else list(v)
                   for k, v in x_.items()} for name, x_ in st.items()}
    return Result(status, lines, data)
