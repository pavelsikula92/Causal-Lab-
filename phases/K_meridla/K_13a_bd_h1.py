"""
K.13a -- The expected BD operator of a causal set in the H1 geometry (first order, fixed test fields)
====================================================================================================

Why this test exists (decisions of the author, 22 Oct 2026)
-----------------------------------------------------------
The field equation keeps the Benincasa-Dowker operator and changes its
argument, B_g f(Pi) = rho_s (decision 1). Before any nonlinear solution,
the operator itself has to be checked in the geometry it is meant for:
H1, ds^2 = -Pi^2 dt^2 + Pi^-2 dx^2 (gamma = 1), where the cones are no
longer flat (speed Pi^2). Its continuum limit should be box_g - R/2 in
every dimension (Benincasa-Dowker 2010, Belenchia et al. 2016), so for a
fixed static test field w and Pi = 1 + e u

    E_g[B_g w] = B_0[w] + e dB_u[w] + O(e^2),
    dB_u[w] -> 2 gamma u lap w + (1 - gamma(n-1)) w lap u + (1 - gamma(n-2)) grad u.grad w

(sympy, C6 of docs/derivations/local_form_symbolic.py; n = d_s = 3). The
test measures dB_u[w] of the DISCRETE construction (Poisson sprinkling of
the H1 geometry, exact layers and smearing of K.4/K.5) and compares it with
that limit. With w = 1 it is the curvature term alone, -R/2 = (1 - gamma(n-1))
lap u = -lap u for H1 in 3+1: the term that cancels the linear response of
equation (a) (decision 1, lab note).

What this test does NOT do (stated with the test, 22 Oct 2026): it does not
compare candidate functions f. The BD operator is linear in the field, so
E_g[f(Pi)] = e f'(1) B_0[u] + O(e^2) and the difference between ln Pi and
Pi - 1 at second order is -(1/2) e^2 B_0[u^2], the FLAT operator on u^2; the
H1 geometry enters only at e^3. The planned K.13b (amplitude dependence of
B_g ln Pi against B_g (Pi - 1)) was dropped for this reason; f must come from
a principle (docs/derivations/ln_pi_potencial.md) and be tested in the full
nonlinear equation (F1.7 N).

Instrument (core.operators.bd_expect_h1, definitions in its docstring)
---------------------------------------------------------------------
Causal order of H1: y < x iff t_x - t_y >= T(y, x), T the optical time
(Fermat, int Pi^-(1+gamma) dl), first order T = L - (1+gamma) e U (straight
segment). Sprinkling: Poisson, rho per proper volume (sqrt(-g) = Pi^(1-gamma n)).
Interval volume: V = int sqrt(-g) [dt - T(y,z) - T(z,x)]_+ d^n z. Expected
operator with the fundamental rho in the prefactor. First order in e by
quadrature, split into the density part (1 - gamma n) P (flat cones, the K.11
instrument) and the cone part (1 + gamma) Q; their limits (C7) are
P -> (-2 u lap w + (n-1) grad u.grad w + n w lap u)/D, Q -> (2 u lap w + 2 grad u.grad w + w lap u)/D.
Lengths in l_k (eps rho = 1), u = exp(-r^2/2 sigma^2), w = 1 and
w = exp(-r^2/2 (2 sigma)^2), rows at r = 0, sigma, 2 sigma, t_x = 3R.

Parts
-----
geom -- the H1 geometry the instrument assumes (independent code paths).
  (a) Fermat: the optical time of a numerically optimised polyline path
      (12 interior vertices, the exact Pi^-(1+gamma), gamma = 1 and 0) against
      the first-order straight-line T, at e = 0.04, 0.02, 0.01, for random
      pairs. If the first-order T is right, the difference is O(e^2):
      log2 of its ratio between successive e must be >= 1.75 for every pair
      (an error of first order would give 1; 1.75 is three quarters of the
      way to 2).
  (b) Interval volume: A1 = int u l0 and A2 = int [U(y,z) + U(z,x) - (L/dt) U(y,x)]
      over flat diamonds, spheroid quadrature at the main level against a
      Monte Carlo integral over the box around the diamond with the plain
      erf segment formula (random diamonds, L up to 3 sigma, tau up to sigma).
      Criterion: max |z| <= z_max (z_max = 3.5, as K.11).
  (c) Long thin intervals: the cancellation-free S = U(y,z) + U(z,x) - U(y,x)
      (GaussianField.broken_excess) against mpmath at 50 digits, random
      needles of length up to 160 sigma. Criterion: relative error <= 1e-6
      (or absolute <= 1e-15 sigma, one rounding of a single U).
  PASSED if (a), (b), (c) hold.

valid -- the instrument against the two operators already calibrated.
  (a) Pi = 1: the flat value from the same nodes against
      core.operators.bd_expect.expected_apply_general (K.4); relative <= 1e-9.
  (b) gamma = -1 (conformal geometry, no cone change): the first-order
      response against a central difference (h = 1e-3) of the K.11 instrument
      core.operators.bd_expect_curved.curved_radial_rows (C^2 spline geometry
      on a fine grid, the same quadrature), sigma = 4, all rows, both fields.
      Criterion: |difference| <= |continuum(gamma = -1)| / 64 (the quadrature
      tolerance of K.12).
  PASSED if (a) and (b) hold.

lin -- the measurement (gamma = 1, H1). Runs only after full geom and valid
  runs of THIS test file PASSED (re-read from the records, same file hash).
  sigma = 4, 8, 16, 32 (l_k / sigma from 1/4 to 1/32), R = 10 sigma, every
  row and field at the main quadrature level and at a coarser check level.
  rho(sigma) = dE (discrete) / dB (continuum, gamma = 1) for each (field, row).
  Alternatives (continuum, gamma = 1): "flat cones" (only the density
  change, Q dropped: the conformal K.11 operator with the H1 density) and
  "box_g without -R/2". d_alt = the distance of the NEARER alternative from 1
  in units of the target (w = 1: flat cones 1.5, box_g 0 -> d_alt = 0.5).
  Criterion (stated before the full runs):
    INSTRUMENT_FAILED  at some sigma, |rho_main - rho_check| > d_alt / 16
                       (a quarter of the decision band, the quarter rule);
    PASSED   (i) |1 - rho| decreases strictly along the sigma ladder for every
             (field, row) and (ii) at the largest sigma |1 - rho| <= d_alt / 4
             for every (field, row);
    FAILED   otherwise.
  Reported separately, not averaged into the verdict: the finite-scale
  systematic |1 - rho(sigma_max)| and the local exponent
  log2(|1 - rho(sigma_max/2)| / |1 - rho(sigma_max)|) for each (field, row);
  the coefficient of lap u for w = 1 (expected 1 - gamma(n-1) = -1); the
  density and cone parts separately against their limits; gamma = 0 and
  gamma = -1 from the same parts.

Seen before this test was written (development runs, 22 Oct 2026, coarse
quadrature [80, 12, 8, 3, 8, 4], decide nothing): rho for w = 1 at r = 0:
0.743, 0.894, 0.961, 0.987 for sigma = 4, 8, 16, 32 (deficit falling roughly
as sigma^-1.5 towards sigma^-2); P and Q parts with nearly the same ratio at
each sigma; at sigma = 64 the coarse level is not converged. Levels compared
at sigma = 4 (r = 2 sigma), 32 and 64 (r = sigma): [120, 12, 16, 4, 12, 6] and the
main level [160, 16, 24, 6, 16, 8] agree to <= 5e-4 in every ratio at sigma 4
and 32; at sigma = 64 the cone part for w = 1 still moves by ~1 %, so the
ladder stops at 32. The development
also found and fixed two precision losses of the instrument (needle volumes:
sigma_max - 1 computed as dt/L - 1; S as a difference of three integrals of
size sigma), both now computed without cancellation and checked by geom (c).
The criteria use only the alternatives and the quarter rule.
"""

from __future__ import annotations

import math
import time
from pathlib import Path

import numpy as np

from core import provenance
from core.operators.bd_expect import expected_apply_general
from core.operators.bd_expect_curved import curved_radial_rows
from core.operators.bd_expect_h1 import (GaussianField, _spheroid_A, continuum_parts, first_order_row, flat_value,
                                         response)
from core.order.sprinkle import derive_seed, interval_volume
from core.stats.parallel import parallel_map
from core.testing import Result

TEST = {
    "id": "K.13a",
    "phase": "K",
    "title": "Očekávaný BD operátor v geometrii H1 (první řád, pevná testovací pole)",
    "theory_part": "měřidlo operátoru B_g pro rovnici pole B_g f(Π) = ρ_s v geometrii H1 (rozhodnutí 1); bez volby f "
                   "a bez fyzikálního závěru",
    "question": "Reprodukuje diskrétní konstrukce (Poissonův sprinkling geometrie −Π²dt² + Π⁻²dx², vrstvy a rozmazání "
                "K.4/K.5) v prvním řádu spojitou limitu □_g − R/2 na pevných testovacích polích, včetně členu −R/2 "
                "= −∇²u, a konverguje k ní s ℓ_k/σ → 0?",
    "method": "první řád střední hodnoty operátoru kvadraturou: optický čas (Fermat), objem intervalu se zakřivenými "
              "kužely, hustota √−g; rozklad na hustotní a kuželovou část; pole w = 1 a gaussovské w, řádky r = 0, σ, 2σ, "
              "σ = 4, 8, 16, 32 ℓ_k; kontroly geometrie (optimalizovaná dráha, Monte Carlo objemů, mpmath), Π = 1 proti "
              "K.4 a γ = −1 proti přístroji K.11",
    "prediction": "poměr k □_g − R/2 se blíží 1 s rostoucím σ; vývoj (hrubá kvadratura): w = 1, r = 0: 0,743 / 0,894 / "
                  "0,961 / 0,987",
    "criterion": "geom: Fermat log2-poměr odchylek >= 1,75; Monte Carlo max |z| <= 3,5; tenké intervaly rel. chyba <= 1e-6. "
                 "valid: Π = 1 rel. <= 1e-9; γ = −1 proti K.11 |rozdíl| <= |spojitá hodnota|/64. "
                 "lin (jen po PROŠLÝCH geom a valid téhož souboru): MĚŘIDLO, pokud |ρ_hlavní − ρ_kontrolní| > d_alt/16; "
                 "PROŠEL, pokud |1 − ρ| ostře klesá s σ pro každé pole a řádek a při největším σ |1 − ρ| <= d_alt/4; "
                 "jinak NEPROŠEL. d_alt = vzdálenost bližší alternativy (ploché kužely, □_g bez −R/2) od 1.",
    "null_model": "alternativy: ploché kužely s hustotou H1 (bez kuželové části) a □_g bez −R/2",
    "level": "calibration",
    "meters": ["K.4", "K.11"],
    "depends_on": ["K.11"],
    "params": {"D": 4, "w_factor": 2.0, "row_fracs": [0.0, 1.0, 2.0], "R_factor": 10.0,
               "quad": [160, 16, 24, 6, 16, 8], "quad_check": [80, 12, 8, 3, 8, 4],
               "sigmas": [4.0, 8.0, 16.0, 32.0],
               # geom
               "geom_sigma": 4.0, "fermat_pairs": 8, "fermat_eps": [0.04, 0.02, 0.01], "fermat_vertices": 12,
               "fermat_slope": 1.75, "mc_diamonds": 20, "mc_samples": 400000, "z_max": 3.5,
               "needles": 40, "needle_tol": 1e-6,
               # valid
               "valid_sigma": 4.0, "fd_h": 1e-3, "flat_tol": 1e-9, "valid_frac": 64.0},
    "parts": {"geom": {"mode": "geom"}, "valid": {"mode": "valid"}, "lin": {"mode": "lin"}},
    "quick": {"quad": [40, 10, 8, 3, 6, 4], "quad_check": [30, 10, 6, 2, 4, 4], "sigmas": [2.0, 4.0],
              "fermat_pairs": 2, "mc_diamonds": 3, "mc_samples": 40000, "needles": 6, "valid_sigma": 3.0,
              "R_factor": 8.0},
    "estimate": "geom: minuty; valid: minuty; lin: 24 řádkových úloh, hlavní úroveň ~10–20 min na řádek mimo osu "
                "(1 jádro; vývoj 22. 10.: ~20 min při sdíleném jádře), řádky na ose ~1 min; s 10 procesy ~30–40 min",
}

_EPS, _RHO = 0.05, 20.0          # eps * rho = 1: lengths in units of l_k
_FIELDS = ("one", "gauss")


# ============================================================================================== geom

def _fermat_time(Y, X, u, e, gamma, K, gl):
    """Optical time int Pi^-(1+gamma) dl of the best polyline from Y to X with K interior vertices (exact Pi)."""
    from scipy.optimize import minimize
    xi, wi = gl
    tq = 0.5 * (xi + 1.0)
    wq = 0.5 * wi

    def T_of(flat):
        P = np.vstack([Y, flat.reshape(K, 3), X])
        A, B = P[:-1], P[1:]
        seg = np.linalg.norm(B - A, axis=1)
        pts = A[:, None, :] + tq[None, :, None] * (B - A)[:, None, :]
        Pi = 1.0 + e * u.value((pts ** 2).sum(axis=2))
        return float((seg[:, None] * wq[None, :] * Pi ** (-(1.0 + gamma))).sum())

    straight = (Y[None, :] + (np.arange(1, K + 1) / (K + 1.0))[:, None] * (X - Y)[None, :]).ravel()
    res = minimize(T_of, straight, method="BFGS", options={"gtol": 1e-11, "maxiter": 2000})
    return min(float(res.fun), T_of(straight))


def fermat_check(p, rng):
    """Part geom (a): |T_opt - T_lin| = O(e^2) for random pairs (module notes)."""
    sig = float(p["geom_sigma"])
    u = GaussianField(sig)
    gl = np.polynomial.legendre.leggauss(16)
    out = []
    for k in range(int(p["fermat_pairs"])):
        Y = rng.normal(size=3) * sig
        X = rng.normal(size=3) * sig
        for gamma in (1.0, 0.0):
            U = float(u.segment([np.array(c) for c in Y], [np.array(c) for c in X]))
            L = float(np.linalg.norm(X - Y))
            devs = []
            for e in p["fermat_eps"]:
                T_opt = _fermat_time(Y, X, u, float(e), gamma, int(p["fermat_vertices"]), gl)
                devs.append(abs(T_opt - (L - (1.0 + gamma) * float(e) * U)))
            slopes = [math.log2(devs[i] / devs[i + 1]) for i in range(len(devs) - 1)]
            out.append({"pair": k, "gamma": gamma, "L": L, "U": U, "dev": devs, "slopes": slopes})
    return out


def mc_check(p, rng, quad):
    """Part geom (b): spheroid quadrature of A1, A2 against Monte Carlo in the box around the diamond."""
    sig = float(p["geom_sigma"])
    u = GaussianField(sig)
    n = int(p["D"]) - 1
    out = []
    for k in range(int(p["mc_diamonds"])):
        r_x = rng.uniform(0.0, 2.0 * sig)
        CO = rng.uniform(-1.0, 1.0)
        L = rng.uniform(0.1, 3.0) * sig
        tau = rng.uniform(0.1, 1.0) * sig
        dt = math.sqrt(L * L + tau * tau)
        ex = tau * tau / (dt + L)
        A1q, A2q, _ = _spheroid_A(u, n, r_x, CO, np.array([[L]]), np.array([[dt]]), np.array([[ex]]),
                                  int(quad[3]), int(quad[4]), int(quad[5]), 0.0)
        SN = math.sqrt(1.0 - CO * CO)
        x = np.array([r_x, 0.0, 0.0])
        nh = np.array([CO, SN, 0.0])
        y = x + L * nh
        e1 = np.array([SN, -CO, 0.0])
        e2 = np.array([0.0, 0.0, 1.0])
        a, b = 0.5 * dt, 0.5 * tau                           # semi-axes of the spheroid
        M = int(p["mc_samples"])
        s_ax = rng.uniform(-a, a, M)
        s1 = rng.uniform(-b, b, M)
        s2 = rng.uniform(-b, b, M)
        Z = 0.5 * (x + y)[None, :] + s_ax[:, None] * nh[None, :] + s1[:, None] * e1[None, :] + s2[:, None] * e2[None, :]
        ell = dt - np.linalg.norm(Z - y, axis=1) - np.linalg.norm(Z - x, axis=1)
        inside = ell > 0
        box = (2 * a) * (2 * b) ** 2
        zc = [Z[:, i] for i in range(3)]
        yc = [np.full(M, y[i]) for i in range(3)]
        xc = [np.full(M, x[i]) for i in range(3)]
        f1 = np.where(inside, u.value((Z ** 2).sum(axis=1)) * ell, 0.0)
        Uyx = float(u.segment([np.array(c) for c in y], [np.array(c) for c in x]))
        f2 = np.where(inside, u.segment(yc, zc) + u.segment(zc, xc) - (L / dt) * Uyx, 0.0)
        rec = {"diamond": k, "r_x": r_x, "cos": CO, "L": L, "tau": tau}
        for name, f, q in (("A1", f1, float(A1q[0, 0])), ("A2", f2, float(A2q[0, 0]))):
            m, se = box * float(f.mean()), box * float(f.std(ddof=1)) / math.sqrt(M)
            rec[name] = {"quad": q, "mc": m, "se": se, "z": (q - m) / se if se > 0 else 0.0}
        out.append(rec)
    return out


def needle_check(p, rng):
    """Part geom (c): broken_excess against mpmath (50 digits) on long thin intervals."""
    import mpmath as mp
    mp.mp.dps = 50
    sig = float(p["geom_sigma"])
    u = GaussianField(sig)

    def U_mp(a, b):
        a = [mp.mpf(float(v)) for v in a]
        b = [mp.mpf(float(v)) for v in b]
        d = [bb - aa for aa, bb in zip(a, b)]
        Ln = mp.sqrt(sum(v * v for v in d))
        dh = [v / Ln for v in d]
        ad = sum(v * w_ for v, w_ in zip(a, dh))
        bd = sum(v * w_ for v, w_ in zip(b, dh))
        m2 = sum(v * v for v in a) - ad ** 2
        s = mp.mpf(sig)
        k = 1 / (mp.sqrt(2) * s)
        return s * mp.sqrt(mp.pi / 2) * mp.e ** (-m2 / (2 * s * s)) * (mp.erf(bd * k) - mp.erf(ad * k))

    out = []
    for k in range(int(p["needles"])):
        L = sig * 10.0 ** rng.uniform(-2.0, math.log10(160.0))
        x = rng.normal(size=3) * sig * rng.uniform(0.0, 2.0)
        d = rng.normal(size=3)
        d /= np.linalg.norm(d)
        y = x + L * d
        perp = rng.normal(size=3)
        perp -= perp.dot(d) * d
        perp *= rng.uniform(0.0, 1.0) * min(L, 2.0) / max(np.linalg.norm(perp), 1e-30) * 10.0 ** rng.uniform(-3.0, 0.0)
        z = y + rng.uniform() * (x - y) + perp
        exact = float(U_mp(y, z) + U_mp(z, x) - U_mp(y, x))
        got = float(u.broken_excess([np.array(c) for c in y], [np.array(c) for c in z], [np.array(c) for c in x]))
        err = abs(got - exact)
        ok = err <= float(p["needle_tol"]) * abs(exact) or err <= 1e-15 * sig
        out.append({"L": L, "exact": exact, "got": got, "abs_err": err, "ok": bool(ok)})
    return out


def run_geom(p, ctx):
    rng = np.random.default_rng(derive_seed(ctx.seed, "K.13a-geom"))
    t0 = time.time()
    fer = fermat_check(p, rng)
    ctx.progress(f"[K.13a] geom (a) Fermat done ({time.time()-t0:.0f}s)")
    mc = mc_check(p, rng, p["quad"])
    ctx.progress(f"[K.13a] geom (b) Monte Carlo done ({time.time()-t0:.0f}s)")
    nd = needle_check(p, rng)
    ctx.progress(f"[K.13a] geom (c) needles done ({time.time()-t0:.0f}s)")
    min_slope = min(min(f["slopes"]) for f in fer)
    zmax = max(max(abs(r["A1"]["z"]), abs(r["A2"]["z"])) for r in mc)
    ok = {"(a) Fermat": min_slope >= float(p["fermat_slope"]),
          "(b) interval volumes": zmax <= float(p["z_max"]),
          "(c) thin intervals": all(r["ok"] for r in nd)}
    lines = [f"H1 geometry, sigma = {p['geom_sigma']:g} l_k, 3+1",
             f"(a) Fermat: {len(fer)} pairs x gamma (1, 0), e = {p['fermat_eps']}: smallest log2 ratio of |T_opt - T_lin| "
             f"{min_slope:.3f} (needs >= {p['fermat_slope']})"]
    for f in fer:
        lines.append(f"     pair {f['pair']} gamma {f['gamma']:+.0f} L {f['L']:.2f} U {f['U']:.3f}: |T_opt - T_lin| "
                     + " ".join(f"{d:.2e}" for d in f["dev"]) + "  slopes " + " ".join(f"{s:.2f}" for s in f["slopes"]))
    lines.append(f"(b) interval volumes: {len(mc)} diamonds, quadrature {list(p['quad'][3:])} vs Monte Carlo "
                 f"({p['mc_samples']} samples): max |z| {zmax:.2f} (needs <= {p['z_max']})")
    for r in mc:
        lines.append(f"     L {r['L']:.2f} tau {r['tau']:.2f} r_x {r['r_x']:.2f}: A1 {r['A1']['quad']:.5g} vs {r['A1']['mc']:.5g} "
                     f"(z {r['A1']['z']:+.2f}); A2 {r['A2']['quad']:.5g} vs {r['A2']['mc']:.5g} (z {r['A2']['z']:+.2f})")
    worst = max(nd, key=lambda r: r["abs_err"] / max(abs(r["exact"]), 1e-300))
    lines.append(f"(c) thin intervals: {sum(r['ok'] for r in nd)}/{len(nd)} within tolerance; worst relative error "
                 f"{worst['abs_err'] / max(abs(worst['exact']), 1e-300):.2e} (L = {worst['L']:.1f})")
    lines.append("checks: " + "; ".join(f"{k}: {'OK' if v else 'FAILS'}" for k, v in ok.items()))
    return Result("PASSED" if all(ok.values()) else "FAILED", lines,
                  {"fermat": fer, "mc": mc, "needles": nd, "checks": ok})


# ============================================================================================== valid

def run_valid(p, ctx):
    D, n = int(p["D"]), int(p["D"]) - 1
    sig = float(p["valid_sigma"])
    R = float(p["R_factor"]) * sig
    t_x = 3.0 * R
    u, w = GaussianField(sig), GaussianField(float(p["w_factor"]) * sig)
    quad = [int(q) for q in p["quad_check"]]
    rows = [f * sig for f in p["row_fracs"]]
    t0 = time.time()
    parts = parallel_map(first_order_row, [(D, _RHO, _EPS, t_x, r, R, u, w, quad) for r in rows], ctx.workers,
                         ctx.progress, lambda i, r_, d, nn: f"[K.13a] valid: H1 rows {d}/{nn} ({time.time()-t0:.0f}s)")
    # (a) Pi = 1 against the flat instrument of K.4
    flat_rows = []
    for r, pt in zip(rows, parts):
        mine = flat_value(pt, D, _RHO, _EPS)
        ref = expected_apply_general(D, _RHO, _EPS, t_x, r, R, lambda d_, rr, tt: w.value(rr * rr), pt["w_x"],
                                     n_rho=quad[0], n_s=quad[1], n_ang=quad[2])
        flat_rows.append({"r": r, "mine": mine, "ref": ref, "rel": abs(mine - ref) / abs(ref)})
    # (b) gamma = -1 against a central difference of the K.11 instrument
    h = float(p["fd_h"])
    dr = sig / 16.0
    nodes = np.arange(0.0, R + 0.5 * dr, dr)
    nodes[-1] = max(nodes[-1], R)
    gnodes = np.arange(0.0, 1.25 * R + 0.5 * dr / 2, dr / 2)
    wn = w.value(nodes ** 2)
    k11 = parallel_map(curved_radial_rows,
                       [(D, _RHO, _EPS, t_x, rows, R, nodes, np.ones_like(nodes), quad[0], quad[1], quad[2], quad[3],
                         quad[4], quad[5], "spline", gnodes, 1.0 + sgn * h * u.value(gnodes ** 2)) for sgn in (1.0, -1.0)],
                       ctx.workers, ctx.progress, lambda i, r_, d, nn: f"[K.13a] valid: K.11 rows {d}/{nn} ({time.time()-t0:.0f}s)")
    Ap, Am = k11
    fd = {"one": (Ap.sum(axis=1) - Am.sum(axis=1)) / (2.0 * h), "gauss": (Ap @ wn - Am @ wn) / (2.0 * h)}
    comp = []
    for i, (r, pt) in enumerate(zip(rows, parts)):
        for fld, ww, key in (("one", None, "one"), ("gauss", w, "w")):
            mine = response(pt, -1.0, D, _RHO, _EPS, key)["total"]
            cont = float(continuum_parts(n, np.array([r]), u, ww)["total"](-1.0)[0])
            tol = abs(cont) / float(p["valid_frac"])
            comp.append({"r": r, "field": fld, "mine": mine, "k11": float(fd[fld][i]), "cont": cont,
                         "diff": abs(mine - float(fd[fld][i])), "tol": tol, "ok": abs(mine - float(fd[fld][i])) <= tol})
    ok = {"(a) Pi = 1 equals K.4": all(f["rel"] <= float(p["flat_tol"]) for f in flat_rows),
          "(b) gamma = -1 equals K.11": all(c["ok"] for c in comp)}
    lines = [f"3+1, sigma = {sig:g} l_k, R = {R:g}, rows r = {', '.join(f'{r:g}' for r in rows)}, quadrature {quad}"]
    for f in flat_rows:
        lines.append(f"(a) r {f['r']:g}: flat E[B w] {f['mine']:+.8e} vs K.4 {f['ref']:+.8e} (rel {f['rel']:.1e})")
    for c in comp:
        lines.append(f"(b) r {c['r']:g} {c['field']:5s}: first order gamma = -1 {c['mine']:+.6e}, K.11 central difference "
                     f"{c['k11']:+.6e}, |diff| {c['diff']:.2e} vs tol {c['tol']:.2e}; continuum {c['cont']:+.6e}")
    lines.append("checks: " + "; ".join(f"{k}: {'OK' if v else 'FAILS'}" for k, v in ok.items()))
    return Result("PASSED" if all(ok.values()) else "FAILED", lines,
                  {"flat": flat_rows, "gamma_m1": comp, "checks": ok, "parts": parts})


# ============================================================================================== lin

def _gates(ctx):
    """lin runs only after full geom and valid runs of THIS test file PASSED."""
    here = provenance.sha256_file(Path(__file__))
    recs = [r for r in provenance.load_all("K.13a") if r.get("test_sha256") == here and r.get("status") == "PASSED"]
    return {v: any(r.get("part") == v for r in recs) for v in ("geom", "valid")}


def _row_task(D, t_x, r, R, sig, wf, quad):
    """One row at one level (module-level for the worker processes)."""
    return first_order_row(D, _RHO, _EPS, t_x, r, R, GaussianField(sig), GaussianField(wf * sig), quad)


def _ratios(pt, D, r, sig, wf, gamma):
    """rho for both fields at gamma, with the density and cone parts against their limits."""
    n = D - 1
    u = GaussianField(sig)
    out = {}
    for fld, ww, key in (("one", None, "one"), ("gauss", GaussianField(wf * sig), "w")):
        c = continuum_parts(n, np.array([r]), u, ww)
        resp = response(pt, gamma, D, _RHO, _EPS, key)
        tot = float(c["total"](gamma)[0])
        out[fld] = {"dE": resp["total"], "cont": tot, "rho": resp["total"] / tot,
                    "rho_P": resp["P"] / float(c["P"][0]), "rho_Q": resp["Q"] / float(c["Q"][0]),
                    "alt_flat_cones": (1.0 - gamma * n) * float(c["P"][0]) / tot,
                    "alt_box_only": float(c["box_only"](gamma)[0]) / tot}
    return out


def run_lin(p, ctx):
    D, n = int(p["D"]), int(p["D"]) - 1
    gates = {"geom": True, "valid": True} if ctx.quick else _gates(ctx)
    if not all(gates.values()):
        missing = [k for k, v in gates.items() if not v]
        return Result("INSTRUMENT_FAILED", [f"lin needs PASSED full runs of {missing} of this file first; nothing computed"],
                      {"gates": gates})
    sigmas = [float(s) for s in p["sigmas"]]
    wf = float(p["w_factor"])
    t0 = time.time()
    jobs, keys = [], []
    for level in ("check", "main"):
        quad = [int(q) for q in p["quad" if level == "main" else "quad_check"]]
        for sig in sigmas:
            R = float(p["R_factor"]) * sig
            for f in p["row_fracs"]:
                jobs.append((D, 3.0 * R, f * sig, R, sig, wf, quad))
                keys.append((level, sig, float(f)))
    res = parallel_map(_row_task, jobs, ctx.workers, ctx.progress,
                       lambda i, r_, d, nn: f"[K.13a] lin rows {d}/{nn}: {keys[i][0]} sigma {keys[i][1]:g} r/sigma "
                                            f"{keys[i][2]:g} ({time.time()-t0:.0f}s)")
    parts = {k: v for k, v in zip(keys, res)}
    R1 = {k: _ratios(v, D, k[2] * k[1], k[1], wf, 1.0) for k, v in parts.items()}
    lines = [f"H1 (gamma = 1), 3+1, lengths in l_k; u = exp(-r^2/2 sigma^2), w = 1 and exp(-r^2/2 ({wf:g} sigma)^2); "
             f"R = {p['R_factor']:g} sigma, rows r/sigma = {p['row_fracs']}; main {list(p['quad'])}, check {list(p['quad_check'])}"]
    inst_ok, monotone, final_ok = True, True, True
    table = []
    for fld in _FIELDS:
        for f in p["row_fracs"]:
            f = float(f)
            seq = [R1[("main", s, f)][fld] for s in sigmas]
            alt = (seq[0]["alt_flat_cones"], seq[0]["alt_box_only"])
            d_alt = min(abs(a - 1.0) for a in alt)
            devs = [abs(1.0 - x["rho"]) for x in seq]
            dq = [abs(R1[("main", s, f)][fld]["rho"] - R1[("check", s, f)][fld]["rho"]) for s in sigmas]
            ok_q = all(d <= d_alt / 16.0 for d in dq)
            ok_m = all(devs[i + 1] < devs[i] for i in range(len(devs) - 1))
            ok_f = devs[-1] <= d_alt / 4.0
            inst_ok &= ok_q
            monotone &= ok_m
            final_ok &= ok_f
            expo = math.log2(devs[-2] / devs[-1]) if len(devs) > 1 and devs[-1] > 0 else float("nan")
            table.append({"field": fld, "r_frac": f, "rho": [x["rho"] for x in seq], "dev": devs, "d_alt": d_alt,
                          "alt": alt, "quad_diff": dq, "ok_quad": ok_q, "ok_monotone": ok_m, "ok_final": ok_f,
                          "systematic": devs[-1], "local_exponent": expo,
                          "rho_P": [x["rho_P"] for x in seq], "rho_Q": [x["rho_Q"] for x in seq],
                          "dE": [x["dE"] for x in seq], "cont": [x["cont"] for x in seq]})
            lines.append(f"  {fld:5s} r = {f:g} sigma: rho(sigma = {', '.join(f'{s:g}' for s in sigmas)}) = "
                         + " ".join(f"{x['rho']:.4f}" for x in seq)
                         + f"; alternatives flat cones {alt[0]:.3f}, box_g {alt[1]:.3f} -> d_alt {d_alt:.3f}; "
                           f"|1 - rho| at sigma {sigmas[-1]:g}: {devs[-1]:.4f} vs d_alt/4 = {d_alt/4:.4f}")
            lines.append("        density part P / limit: " + " ".join(f"{x['rho_P']:.4f}" for x in seq)
                         + "; cone part Q / limit: " + " ".join(f"{x['rho_Q']:.4f}" for x in seq)
                         + f"; quadrature |main - check|: " + " ".join(f"{d:.4f}" for d in dq) + f" (tol {d_alt/16:.4f})")
            lines.append(f"        finite-scale systematic at sigma {sigmas[-1]:g}: {devs[-1]:.4f}, local exponent "
                         f"log2(dev(sigma/2)/dev(sigma)) = {expo:.2f}; strictly decreasing: {'yes' if ok_m else 'NO'}")
    # the coefficient of lap u for w = 1 (expected 1 - gamma (n - 1) = -1)
    coef = []
    for s in sigmas:
        for f in p["row_fracs"]:
            r = float(f) * s
            lu = float(GaussianField(s).radial(np.array([r]), n)[2][0])
            coef.append({"sigma": s, "r_frac": float(f), "coef": R1[("main", s, float(f))]["one"]["dE"] / lu})
    lines.append("  coefficient of lap u for w = 1 (expected -1): " + "; ".join(
        f"sigma {c['sigma']:g} r/sigma {c['r_frac']:g}: {c['coef']:+.4f}" for c in coef))
    for g in (0.0, -1.0):
        Rg = {k: _ratios(v, D, k[2] * k[1], k[1], wf, g) for k, v in parts.items() if k[0] == "main"}
        lines.append(f"  gamma = {g:+g} (same parts): " + "; ".join(
            f"{fld} r/sigma {float(f):g}: " + " ".join(f"{Rg[('main', s, float(f))][fld]['rho']:.4f}" for s in sigmas)
            for fld in _FIELDS for f in p["row_fracs"]))
    if not inst_ok:
        status = "INSTRUMENT_FAILED"
    elif monotone and final_ok:
        status = "PASSED"
    else:
        status = "FAILED"
    lines.append(f"decision: quadrature within d_alt/16: {'OK' if inst_ok else 'FAILS'}; |1 - rho| strictly decreasing: "
                 f"{'yes' if monotone else 'NO'}; at sigma {sigmas[-1]:g} within d_alt/4: {'yes' if final_ok else 'NO'}")
    data = {"table": table, "coef_lap_u": coef, "gates": gates,
            "parts": [{"level": k[0], "sigma": k[1], "r_frac": k[2], **v} for k, v in parts.items()]}
    return Result(status, lines, data)


def run(p, ctx):
    mode = str(p["mode"])
    if mode == "geom":
        return run_geom(p, ctx)
    if mode == "valid":
        return run_valid(p, ctx)
    if mode == "lin":
        return run_lin(p, ctx)
    raise ValueError(mode)
