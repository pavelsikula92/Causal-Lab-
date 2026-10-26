"""
F1.9a -- beta of the discrete field equation B_g ln Pi = kappa rho_s in 3+1 (H1): the second-order response outside the source
=============================================================================================================================

Why this test exists (author's decision, 25 Oct 2026)
-----------------------------------------------------
The continuum equation B_g ln Pi = kappa rho_s in the H1 geometry
(ds^2 = -Pi^2 dt^2 + Pi^-2 dx^2, gamma = 1) predicts the PPN parameter
beta = 1 (F1.7, sympy). Is that a property of the DISCRETE model (the
expected Benincasa-Dowker operator of a causal set sprinkled into H1), or
only of its continuum limit? K.13a verified the first-order response of the
discrete operator for Gaussian test fields near the source (rows 0, sigma,
2 sigma). beta, however, is read OUTSIDE the source, in the 1/r field. This
test measures the discrete response exactly there.

Derivation (docs/derivations/k12_dno.md does not cover it; stated here)
------------------------------------------------------------------------
Expand L = ln Pi = e L1 + e^2 L2. Order e: B_0[L1] = kappa rho_s (the flat
operator, K.4, F1.3). Order e^2:

    B_0[L2] = - dB_{L1}[L1] + (source terms),

dB_u[w] the first-order change of the operator acting on w when the geometry
changes by u (core.operators.bd_expect_h1, the instrument of K.13a). Every
source term is proportional to rho_s (a volume factor sqrt(-g) = e^((1-3 gamma) L),
proper or coordinate density, the normalisation of rho_s): it lives inside the
source, its effect outside is harmonic (1/r) and only renormalises the active
mass. The coefficient beta is the 1/r^2 part ~ L1^2 of g_tt and does NOT
depend on the source convention. Outside the source:

    lap L2 = -S,   S = dB_{L1}[L1],   S = c |grad L1|^2  =>  L2 = -(c/2) L1^2,  beta = 1 - c/2.

Continuum: dB_u[w] = 2 gamma u lap w + (1 - gamma (n-2)) grad u.grad w + (1 - gamma (n-1)) w lap u
(n = 3); outside the source lap L1 = 0, so c = 1 - gamma:
    H1 (gamma = 1):          c = 0  (beta = 1)        -- the prediction tested here;
    conformal (gamma = -1):  c = 2                     -- positive control from the SAME quadrature.
The alternatives that define the decision band are the neighbouring
candidates f = (Pi^k - 1)/k, k = +-1, of the field equation (beta = 1 - k/2,
F1.7): they differ from ln Pi by Delta c = +-1.

The field
---------
u = w = L1-shaped field of a Gaussian source of width sigma:
    (sigma / r) erf(r / sqrt2 sigma) = sqrt(2/pi) int_0^1 exp(-r^2 t^2 / 2 sigma^2) dt,
represented EXACTLY as the finite mixture
    L(r) = sum_i a_i exp(-r^2 t_i^2 / 2 sigma^2),  t_i, a_i = sqrt(2/pi) x Gauss-Legendre on [t_min, 1],
with t_min = sigma / R (the widest component has the width of the box). The mixture is the field of
the test (its continuum response is computed exactly from it), so the representation error does not
enter the comparison; it only makes the field slightly less harmonic than 1/r, which is reported
(continuum c of the mixture, local exponent of |grad L|). The instrument is linear in u (the geometry)
and accepts any w (it needs only w(y)), so the discrete response of the mixture is the exact sum of the
responses of its Gaussian components, each computed by the calibrated K.13a code.

Protocol (stated before any full run)
-------------------------------------
Rows outside the source only: r = 3, 4, 6 sigma (the Gaussian source density is below 1.1 % of its centre
at r = 3 sigma); no row is excluded (|grad L|^2 is monotone and far from 0 at these rows). sigma = 4, 8, 16
l_k (the node step of the radial quadrature is fixed in l_k, as K.13a), R = 12 sigma, t_x = 3 R. Every row
at the main level [160, 16, 24, 6, 16, 8] and the check level [80, 12, 8, 3, 8, 4] (K.13a). Wall control:
sigma = 4 with R = 24 sigma (the same mixture, a box twice as large), main level.
For each (sigma, row) and gamma = +1, -1:  c_disc = dE/|grad L|^2, c_cont = dB_cont/|grad L|^2, dev = |c_disc - c_cont|.
    INSTRUMENT_FAILED  |c_main - c_check| > 1/16 for some (sigma, row, gamma)   (quarter of the tolerance),
                       or the wall control moves c by more than 1/16,
                       or the positive control (gamma = -1) has dev > 1/4 at sigma = 16 in some row
                       (the instrument would not reproduce a known nonzero response);
    PASSED             gamma = +1: dev <= 1/4 at sigma = 16 in every row (a quarter of the band Delta c = 1),
                       and dev(16) <= dev(4) + 1/16 in every row (no growth with sigma / l_k);
    FAILED             otherwise (the discrete operator gives beta != 1 outside the source).
Reported, not decided: the profile c(r) of both gammas at every sigma (not averaged), beta_disc = 1 - c/2,
the density and cone parts against their limits, the linear extrapolation c_inf = 2 c(16) - c(8) (first order
in l_k / sigma; K.13a found local exponents 1.56-1.65), the continuum c of the mixture, and d ln|grad L| / d ln r.

Review of 25 Oct 2026 (before any full result was seen; two changes, the bands are unchanged)
-------------------------------------------------------------------------------------------
(1) A Gaussian source has no compact support, so a row is "outside the source" only if the source-convention
    terms (q kappa rho_s L1, |q| <= |1 - 3 gamma| = 2) cannot move c noticeably. Their bound
    c_src = 2 |L1 lap L1| / |grad L1|^2 for the exact L1 is 0.51 at r = 3 sigma, 0.034 at 4 sigma and 5e-6 at
    6 sigma. Only rows with c_src <= 1/16 (the instrument tolerance) enter the beta decision (4 and 6 sigma); the
    row 3 sigma is computed and reported as a check of the operator against the continuum of the same field.
    The quadrature, wall and control checks of the instrument still apply to every row.
(2) The instrument gate accepts only the LATEST full K.13a lin record, which must be PASSED and carry the same
    core hash as the current run (a later failure or any change of core closes the gate).

Crash resilience: every task (one Gaussian component at one row and level) is cached in
results/F1.9a/cache, keyed by its inputs and the hash of the instrument file, so a run interrupted by a
crash continues where it stopped (quick runs are not cached).

Seen before this test was written: K.13a (first order near the source, Gaussian fields), the continuum
derivation above. No discrete exterior value of c is known.
"""

from __future__ import annotations

import hashlib
import json
import math
import time
from pathlib import Path

import numpy as np

from core import provenance
from core.operators.bd_expect_h1 import GaussianField, first_order_row, response
from core.stats.parallel import parallel_map
from core.testing import Result

TEST = {
    "id": "F1.9a",
    "phase": "F1",
    "title": "β diskrétní rovnice B_g ln Π = κρ_s ve 3+1 (H1): odezva druhého řádu mimo zdroj",
    "theory_part": "rovnice pole B_g ln Π = κρ_s v geometrii H1 (rozhodnutí 1, Path A); PPN β diskrétního modelu",
    "question": "Zachová diskrétní BD operátor v poli 1/r mimo zdroj nulovou odezvu H1 (c = 0, β = 1) a nenulovou "
                "odezvu konformní kontroly (c = 2), a konverguje k nim s σ/ℓ_k?",
    "method": "druhý řád rozvoje L = ln Π: S = dB_{L1}[L1] mimo zdroj přístrojem K.13a (první řád operátoru v H1), pole "
              "L1 tvaru erf/r jako přesná směs gaussovských složek; c = S/|∇L1|² v řádcích r = 3, 4, 6 σ, σ = 4, 8, 16 ℓ_k, "
              "R = 12σ, γ = +1 (H1) a γ = −1 (kontrola) ze stejné kvadratury, dvě úrovně kvadratury, kontrola stěny",
    "prediction": "c = 0 pro H1 (β = 1), c = 2 pro γ = −1 (spojitá limita)",
    "criterion": "MĚŘIDLO: |c_hlavní − c_kontrolní| > 1/16, stěna mění c o více než 1/16, nebo kontrola γ = −1 odchylka "
                 "> 1/4 při σ = 16. PROŠEL: H1 |c_disc − c_spoj| <= 1/4 při σ = 16 v každém řádku a odchylka neroste s σ "
                 "(dev(16) <= dev(4) + 1/16). NEPROŠEL: jinak. Pásmo Δc = 1 (sousední kandidáti f, β = 1 ∓ 1/2).",
    "null_model": "alternativy f = (Π^k − 1)/k, k = ±1: c = ±1 (β = 1/2, 3/2); konformní geometrie c = 2",
    "level": "mechanism",
    "meters": ["K.4", "K.13a"],
    "depends_on": ["K.13a", "F1.7"],
    "params": {"D": 4, "sigmas": [4.0, 8.0, 16.0], "row_fracs": [3.0, 4.0, 6.0], "R_factor": 12.0,
               "wall_R_factor": 24.0, "wall_sigma": 4.0, "n_comp": 12,
               "quad": [160, 16, 24, 6, 16, 8], "quad_check": [80, 12, 8, 3, 8, 4],
               "tol_c": 0.25, "tol_inst": 0.0625, "cache": True},
    # added 25 Oct before any full result (review): rows enter the decision only if the source-term bound c_src <= tol_inst
    "parts": {"main": {}},
    "quick": {"sigmas": [2.0, 4.0], "row_fracs": [3.0, 6.0], "R_factor": 10.0, "wall_R_factor": 16.0,
              "wall_sigma": 2.0, "n_comp": 6, "quad": [40, 10, 8, 3, 6, 4], "quad_check": [30, 10, 6, 2, 4, 4],
              "cache": False},
    "estimate": "~250 úloh (složka × řádek × úroveň); hlavní úroveň ~10–20 min na úlohu (1 jádro, jako řádek K.13a), "
                "kontrolní ~1–2 min; s 10 procesy odhadem 4–6 h; přerušený běh pokračuje z mezipaměti",
}

_EPS, _RHO = 0.05, 20.0          # eps * rho = 1: lengths in units of l_k (as K.13a)
_N = 3                           # spatial dimensions (D = 4)
_INSTRUMENT = Path(__file__).resolve().parents[2] / "core" / "operators" / "bd_expect_h1.py"
_CACHE = Path(__file__).resolve().parents[2] / "results" / "F1.9a" / "cache"


# =============================================================================
# the field: an exact Gaussian mixture of the shape (sigma/r) erf(r/sqrt2 sigma)
# =============================================================================

def mixture(sigma: float, R: float, n_comp: int):
    """Component widths s_i = sigma / t_i and weights a_i of L(r) = sum a_i exp(-r^2 / 2 s_i^2) (module notes)."""
    t_min = sigma / R
    x, w = np.polynomial.legendre.leggauss(int(n_comp))
    t = t_min + 0.5 * (1.0 - t_min) * (x + 1.0)
    a = math.sqrt(2.0 / math.pi) * 0.5 * (1.0 - t_min) * w
    return [(float(sigma / ti), float(ai)) for ti, ai in zip(t, a)]


class MixtureField:
    """Sum of weighted Gaussians: value (what the instrument needs from w) and exact radial derivatives."""

    def __init__(self, comps):
        self.comps = [(float(s), float(a)) for s, a in comps]

    def value(self, r2):
        r2 = np.asarray(r2, float)
        return sum(a * np.exp(-r2 / (2.0 * s * s)) for s, a in self.comps)

    def radial(self, r, n=_N):
        """(L, L'(r), lap L) at radius r in n spatial dimensions."""
        r = np.asarray(r, float)
        L = np.zeros_like(r)
        dL = np.zeros_like(r)
        lap = np.zeros_like(r)
        for s, a in self.comps:
            g = a * np.exp(-r * r / (2.0 * s * s))
            L += g
            dL += -g * r / (s * s)
            lap += g * (r * r / s ** 4 - n / (s * s))
        return L, dL, lap


def source_bound(r: float, sigma: float) -> float:
    """
    Upper bound of what the source-convention terms can contribute to c at radius r (in units of |grad L1|^2):
    they are q kappa rho_s L1 with |q| <= |1 - 3 gamma| = 2 for H1 (sqrt(-g) of the metric), and kappa rho_s = -lap L1.
    For the exact L1 = (sigma/r) erf(r / sqrt2 sigma): lap L1 = -sqrt(2/pi) exp(-r^2 / 2 sigma^2) / sigma^2.
    """
    x = r / sigma
    L = math.erf(x / math.sqrt(2.0)) / x
    dL = (-math.erf(x / math.sqrt(2.0)) / (x * x) + math.sqrt(2.0 / math.pi) * math.exp(-0.5 * x * x) / x) / sigma
    lap = -math.sqrt(2.0 / math.pi) * math.exp(-0.5 * x * x) / sigma ** 2
    return 2.0 * abs(L * lap) / (dL * dL)


def continuum_S(field: MixtureField, r: float, gamma: float) -> float:
    """Continuum dB_u[w] for u = w = field: 2 g u lap u + (1 - g (n-2)) |grad u|^2 + (1 - g (n-1)) u lap u."""
    L, dL, lap = (float(v[0]) for v in field.radial(np.array([r])))
    return 2.0 * gamma * L * lap + (1.0 - gamma * (_N - 2)) * dL * dL + (1.0 - gamma * (_N - 1)) * L * lap


def continuum_PQ(field: MixtureField, r: float):
    """Local limits of the density part P and the cone part Q (core.operators.bd_expect_h1.continuum_parts)."""
    L, dL, lap = (float(v[0]) for v in field.radial(np.array([r])))
    D = _N + 1
    P = (-2.0 * L * lap + (_N - 1) * dL * dL + _N * L * lap) / D
    Q = (2.0 * L * lap + 2.0 * dL * dL + L * lap) / D
    return P, Q


# =============================================================================
# one task: one Gaussian component of the geometry u at one row and level (cached)
# =============================================================================

def _key(args) -> str:
    D, t_x, r, R, s_i, comps, quad = args
    blob = json.dumps({"D": D, "t_x": t_x, "r": r, "R": R, "s": s_i, "w": comps, "quad": list(quad),
                       "eps": _EPS, "rho": _RHO, "instr": provenance.sha256_file(_INSTRUMENT)}, sort_keys=True)
    return hashlib.sha256(blob.encode()).hexdigest()[:32]


def component_task(D, t_x, r, R, s_i, comps, quad, use_cache):
    """first_order_row with u = the Gaussian component of width s_i and w = the whole mixture."""
    args = (D, t_x, r, R, s_i, comps, list(quad))
    path = _CACHE / f"{_key(args)}.json"
    if use_cache and path.exists():
        try:
            return json.loads(path.read_text())
        except (OSError, ValueError):
            pass                                          # a damaged cache file is recomputed
    out = first_order_row(D, _RHO, _EPS, t_x, r, R, GaussianField(s_i), MixtureField(comps), quad)
    out = {k: float(v) for k, v in out.items()}
    if use_cache:
        try:
            _CACHE.mkdir(parents=True, exist_ok=True)
            tmp = path.with_suffix(".tmp")
            tmp.write_text(json.dumps(out))
            tmp.replace(path)                             # atomic: a crash never leaves a half-written entry
        except OSError:
            pass
    return out


def combine(parts_list, comps):
    """Parts of the mixture geometry: the instrument is linear in u, so P, Q are the weighted sums."""
    out = {k: 0.0 for k in ("P1", "Q1", "Pw", "Qw")}
    for pt, (_, a) in zip(parts_list, comps):
        for k in out:
            out[k] += a * pt[k]
    out["w_x"] = parts_list[0]["w_x"]
    return out


# =============================================================================
# the test
# =============================================================================

def _gate(ctx):
    """
    The instrument must be calibrated in its CURRENT version: the latest full K.13a lin record (by time) must be
    PASSED and carry the same core hash as now (a later failure, or any change of core since, closes the gate).
    Returns (ok, message).
    """
    recs = [r for r in provenance.load_all("K.13a") if r.get("part") == "lin"]
    if not recs:
        return False, "no full K.13a lin record"
    last = max(recs, key=lambda r: r.get("time_utc", ""))
    if last.get("status") != "PASSED":
        return False, f"the latest K.13a lin record ({last.get('time_utc')}) is {last.get('status')}"
    if last.get("core_sha256") != provenance.core_sha256():
        return False, (f"core changed since the latest K.13a lin PASSED ({last.get('time_utc')}): re-run "
                       f"'python run.py K.13a --part lin' first")
    return True, f"K.13a lin PASSED {last.get('time_utc')} with the current core"


def run(p, ctx):
    D = int(p["D"])
    gate_ok, gate_msg = (True, "quick run: gate not checked") if ctx.quick else _gate(ctx)
    if not gate_ok:
        return Result("INSTRUMENT_FAILED", [f"instrument gate: {gate_msg}; nothing computed"], {"gate": gate_msg})
    sigmas = [float(s) for s in p["sigmas"]]
    fracs = [float(f) for f in p["row_fracs"]]
    nc = int(p["n_comp"])
    use_cache = bool(p["cache"]) and not ctx.quick
    # configurations: (label, sigma, R, level)
    confs = []
    for level in ("check", "main"):
        for sig in sigmas:
            confs.append(("ladder", sig, float(p["R_factor"]) * sig, level))
    ws = float(p["wall_sigma"])
    confs.append(("wall", ws, float(p["wall_R_factor"]) * ws, "main"))
    jobs, keys = [], []
    for lab, sig, R, level in confs:
        # the SAME mixture for the ladder box and the wall box at the wall sigma (t_min from the ladder R)
        comps = mixture(sig, float(p["R_factor"]) * sig, nc)
        quad = [int(q) for q in p["quad" if level == "main" else "quad_check"]]
        for f in fracs:
            for i, (s_i, _) in enumerate(comps):
                jobs.append((D, 3.0 * R, f * sig, R, s_i, comps, quad, use_cache))
                keys.append((lab, sig, R, level, f, i))
    t0 = time.time()
    res = parallel_map(component_task, jobs, ctx.workers, ctx.progress,
                       lambda i, r_, d, nn: f"[F1.9a] tasks {d}/{nn}: {keys[i][0]} {keys[i][3]} sigma {keys[i][1]:g} "
                                            f"R {keys[i][2]:g} r/sigma {keys[i][4]:g} component {keys[i][5]} "
                                            f"({time.time()-t0:.0f}s)")
    by = {}
    for k, v in zip(keys, res):
        by.setdefault(k[:5], [None] * nc)[k[5]] = v

    def measure(lab, sig, R, level, f):
        comps = mixture(sig, float(p["R_factor"]) * sig, nc)
        fld = MixtureField(comps)
        r = f * sig
        parts = combine(by[(lab, sig, R, level, f)], comps)
        _, dL, lap = (float(v[0]) for v in fld.radial(np.array([r])))
        g2 = dL * dL
        Pc, Qc = continuum_PQ(fld, r)
        out = {"grad2": g2, "lapL": lap}
        for gam, tag in ((1.0, "h1"), (-1.0, "conf")):
            resp = response(parts, gam, D, _RHO, _EPS, "w")
            cd = resp["total"] / g2
            cc = continuum_S(fld, r, gam) / g2
            out[tag] = {"c_disc": cd, "c_cont": cc, "dev": abs(cd - cc), "beta_disc": 1.0 - cd / 2.0}
        resp = response(parts, 1.0, D, _RHO, _EPS, "w")
        out["P_ratio"] = resp["P"] / Pc
        out["Q_ratio"] = resp["Q"] / Qc
        # local exponent of |grad L| (ideal 1/r^2: -2)
        h = 0.05 * r
        g_p = abs(float(fld.radial(np.array([r + h]))[1][0]))
        g_m = abs(float(fld.radial(np.array([r - h]))[1][0]))
        out["grad_exponent"] = math.log(g_p / g_m) / math.log((r + h) / (r - h))
        return out

    M = {}
    for lab, sig, R, level in confs:
        for f in fracs:
            M[(lab, sig, level, f)] = measure(lab, sig, R, level, f)
    lines = [f"3+1, lengths in l_k; u = w = Gaussian mixture of (sigma/r) erf(r/sqrt2 sigma), {nc} components, t_min = "
             f"1/{p['R_factor']:g}; rows r/sigma = {fracs}; R = {p['R_factor']:g} sigma, t_x = 3R; main {list(p['quad'])}, "
             f"check {list(p['quad_check'])}; c = dB_L[L] / |grad L|^2, beta = 1 - c/2"]
    tol, tol_i = float(p["tol_c"]), float(p["tol_inst"])
    inst_ok, control_ok, h1_ok, growth_ok = True, True, True, True
    s_max, s_min = sigmas[-1], sigmas[0]
    table = []
    cb = {f: source_bound(f, 1.0) for f in fracs}                 # scale-free: depends on r / sigma only
    clean = [f for f in fracs if cb[f] <= tol_i]
    lines.append(f"  instrument gate: {gate_msg}")
    lines.append("  source-term bound c_src (|1 - 3 gamma| kappa rho_s L1 / |grad L1|^2): "
                 + "; ".join(f"r = {f:g} sigma: {cb[f]:.2e}" for f in fracs)
                 + f" -> rows in the decision (c_src <= {tol_i:g}): {clean}; the others are reported only")
    for f in fracs:
        for tag, name in (("h1", "H1 gamma = +1 (prediction c = 0)"), ("conf", "conformal gamma = -1 (control, c = 2)")):
            seq = [M[("ladder", s, "main", f)][tag] for s in sigmas]
            chk = [M[("ladder", s, "check", f)][tag] for s in sigmas]
            dq = [abs(a["c_disc"] - b["c_disc"]) for a, b in zip(seq, chk)]
            inst_ok &= all(d <= tol_i for d in dq)                # the instrument checks hold in EVERY row
            if f in clean:                                        # the beta decision only outside the source
                if tag == "conf":
                    control_ok &= seq[-1]["dev"] <= tol
                else:
                    h1_ok &= seq[-1]["dev"] <= tol
                growth_ok &= seq[-1]["dev"] <= seq[0]["dev"] + tol_i
            c_inf = 2.0 * seq[-1]["c_disc"] - seq[-2]["c_disc"] if len(seq) > 1 else float("nan")
            table.append({"r_frac": f, "gamma": tag, "c_disc": [x["c_disc"] for x in seq],
                          "c_cont": [x["c_cont"] for x in seq], "dev": [x["dev"] for x in seq], "quad_diff": dq,
                          "c_inf_linear": c_inf, "beta_disc": [x["beta_disc"] for x in seq]})
            lines.append(f"  r = {f:g} sigma{'' if f in clean else ' (reported only: source tail)'}, {name}: c_disc(sigma = {', '.join(f'{s:g}' for s in sigmas)}) = "
                         + " ".join(f"{x['c_disc']:+.4f}" for x in seq) + "; continuum of the mixture "
                         + " ".join(f"{x['c_cont']:+.4f}" for x in seq) + "; |dev| "
                         + " ".join(f"{x['dev']:.4f}" for x in seq) + f"; |main - check| "
                         + " ".join(f"{d:.4f}" for d in dq) + f"; linear extrapolation c_inf {c_inf:+.4f}")
            if tag == "h1":
                lines.append("        beta_disc = 1 - c/2: " + " ".join(f"{x['beta_disc']:.4f}" for x in seq)
                             + "; density part P / limit: "
                             + " ".join(f"{M[('ladder', s, 'main', f)]['P_ratio']:.4f}" for s in sigmas)
                             + "; cone part Q / limit: "
                             + " ".join(f"{M[('ladder', s, 'main', f)]['Q_ratio']:.4f}" for s in sigmas)
                             + "; d ln|grad L| / d ln r: "
                             + " ".join(f"{M[('ladder', s, 'main', f)]['grad_exponent']:.3f}" for s in sigmas))
    # wall control
    wall_ok = True
    for f in fracs:
        for tag in ("h1", "conf"):
            a = M[("ladder", ws, "main", f)][tag]["c_disc"]
            b = M[("wall", ws, "main", f)][tag]["c_disc"]
            wall_ok &= abs(a - b) <= tol_i
            lines.append(f"  wall control sigma {ws:g}, r = {f:g} sigma, {tag}: R = {p['R_factor']:g} sigma {a:+.4f}, "
                         f"R = {p['wall_R_factor']:g} sigma {b:+.4f} (difference {b - a:+.4f}, tolerance {tol_i:g})")
    if not (inst_ok and wall_ok and control_ok):
        status = "INSTRUMENT_FAILED"
    elif h1_ok and growth_ok:
        status = "PASSED"
    else:
        status = "FAILED"
    lines.append(f"decision: quadrature within {tol_i:g}: {'OK' if inst_ok else 'FAILS'}; wall within {tol_i:g}: "
                 f"{'OK' if wall_ok else 'FAILS'}; control gamma = -1 within {tol:g} at sigma {s_max:g}: "
                 f"{'OK' if control_ok else 'FAILS'}; H1 within {tol:g} at sigma {s_max:g} in every decision row: "
                 f"{'yes' if h1_ok else 'NO'}; no growth from sigma {s_min:g}: {'yes' if growth_ok else 'NO'}")
    if not clean:
        status = "INSTRUMENT_FAILED"
        lines.append("no row is outside the source by the bound: nothing to decide")
    data = {"source_bound": cb, "decision_rows": clean, "gate": gate_msg, "table": table, "measures": [{"conf": k[0], "sigma": k[1], "level": k[2], "r_frac": k[3], **v}
                                         for k, v in M.items()],
            "mixture_ladder": {str(s): mixture(s, float(p["R_factor"]) * s, nc) for s in sigmas}}
    return Result(status, lines, data)
