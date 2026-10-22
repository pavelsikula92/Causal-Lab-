"""
F1.3 -- Asymptotic shape and strength of the field of a localized source (BD-3)
==============================================================================

Question. A5' in kernel form says that a static source of total strength Q
produces the field u = Pi - Pi_vac that solves lap u = Q delta (smeared).
Far from the source the field must therefore be the Green function of the
Laplacian with the Gauss-law amplitude:

    3+1:  u(r) = -Q / (4 pi) * 1/r + const      (Newtonian potential)
    2+1:  u(r) =  Q / (2 pi) * ln r + const     (the planar analogue)

In the causal set the operator is not the Laplacian but the EXPECTED
smeared Benincasa-Dowker operator E[B], nonlocal on the scale
l_k = (eps rho)^(-1/d) (core.operators.bd_expect, calibrated by K.4).
BD-3 claims that the nonlocality only reshapes the field near the source:
for r >> l_k the field is again Green-function shaped, with the amplitude
fixed by Q. That is the claim tested here. The answer is not guaranteed:
the BD kernel reaches far along the light cone, and a nonlocal operator
can leave power-law corrections that change the far field.

What is computed (deterministic quadrature, no sprinkling)
-----------------------------------------------------------
All lengths in units of l_k (eps rho = 1): the expected operator depends
on eps and rho only through l_k, so one computation per dimension covers
every N and eps.

  * Source: a Gaussian q exp(-r^2 / 2 sigma^2), sigma = 1 l_k (the source
    width of F1.2b in these units), total strength Q = q (2 pi sigma^2)^(n/2),
    n = D - 1 spatial dimensions.
  * Field: rows of E[B] on a radial hat basis (core.operators.bd_expect.
    expected_radial_rows), placed rows_per_element times per element at
    t_x = 3 R (the bottom of the region never enters; the past is cut only
    by the wall |x| = R, as in a tall cylinder); least squares for the
    nodal values with u(R) = 0.
  * Shape: the best Box-Cox exponent g of u(r) in the window [r_lo, r_hi]
    (core.operators.solve.fit_profile: u = A (r^g - 1)/g + B, g = 0 is the
    logarithm). The Green function has g = -1 in 3+1 and g = 0 in 2+1.
  * Strength: the amplitude A of the predicted shape with g FIXED at its
    predicted value (linear least squares of u on [boxcox(r, g_pred), 1]
    in the window), compared with A_pred = Q / omega_n (omega_2 = 2 pi,
    omega_3 = 4 pi). The constant B stays free: it carries the boundary
    condition and the wall, not the source.

Why these choices (stated before the run; nothing is fitted on the result)
--------------------------------------------------------------------------
  * Window [10, 40] l_k: far from the source (r >= 10 l_k, ten nonlocality
    lengths) and at most R / 2 from the axis (the wall at R = 80 l_k cuts
    the past of every row near it).
  * g_tol = 0.25: the purpose of F1.3 in the plan is to tell the candidate
    laws apart -- 1/r, 1/r^2, 1/r^3 and ln r, i.e. g = -1, -2, -3, 0, spaced
    by 1. A law counts as identified when the fitted g lies within a
    quarter of that spacing, i.e. at least three times closer to the
    predicted law than to any alternative.
  * amp_tol = 0.10: the accuracy of the Pi readback from chain density
    (K.3, err <= 0.10). Later tests read the field through K.3; a field
    strength known better than that is not needed, a worse one would
    dominate their error budget.
  * Convergence and residual tolerances are those of K.4 (conv_g 0.05,
    conv relative 0.02) and F1.2b (res_tol 0.05 q).

Seen before this test was written (disclosed, not used): K.4 "shape"
reported g = -1.05 / -1.04 (3+1) and g = -0.12 (2+1) for a BALL source
(radius 1 l_k) with R = 40 l_k in the window [5, 20] l_k. F1.3 uses a
different source, a doubled domain and a window that starts twice as far
out; the prediction comes from the theory (Green function, Gauss law),
not from those numbers.

Instrument (no conclusion if it fails)
---------------------------------------
  * two resolutions (radial basis AND quadrature nodes) agree:
    |g1 - g2| <= conv_g and |A1 / A2 - 1| <= conv_amp;
  * the field solves the equation: RMS(E[B] u - rhs) / q <= res_tol at
    both resolutions.

Decision (finer resolution res2)
--------------------------------
  S (shape):    |g - g_pred| <= g_tol;
  M (strength): |A_fix / A_pred - 1| <= amp_tol.
PASSED: S and M. FAILED: S or M does not hold (the fitted exponent and
the nearest candidate law are reported either way).

Null model: no sampling enters, so there is no random null. The competing
hypotheses are the alternative laws (g = -3, -2, -1, 0 other than the
predicted one); the shape criterion is set so that it separates them.

Reported only (not part of the decision): g in the two halves of the
window (a trend toward or away from the Green function), the relative
deviation from the continuum solution with the same source and u(R) = 0,
and the local relative deviation at a few radii.
"""

from __future__ import annotations

import math
import time

import numpy as np

from core.operators.bd_expect import expected_radial_rows, shape_nodes
from core.operators.solve import fit_profile, boxcox
from core.stats.stats import fmt
from core.stats.parallel import parallel_map
from core.testing import Result

TEST = {
    "id": "F1.3",
    "phase": "F1",
    "title": "Asymptotika pole lokalizovaného zdroje (BD-3)",
    "theory_part": "A5′ v jádrovém tvaru a jeho nelokální realizace očekávaným BD operátorem (BD-3)",
    "question": "Má pole zdroje řešené očekávaným nelokálním BD operátorem pro r ≫ l_k tvar Greenovy funkce Laplaceova "
                "operátoru (1/r ve 3+1, ln r ve 2+1) s amplitudou danou celkovým zdrojem (Gaussův zákon)?",
    "method": "gaussovský zdroj šířky 1 l_k; řádky očekávaného operátoru na radiální bázi (kvadratura, jednotky l_k, R = 80); "
              "pole nejmenšími čtverci s u(R) = 0; Box-Coxův exponent g v okně [10, 40] l_k; amplituda při g pevném na předpovědi; "
              "dvě rozlišení báze i kvadratury",
    "prediction": "g = −1 (3+1), g = 0 (2+1); amplituda A = Q/ω_n (Q/4π ve 3+1, Q/2π ve 2+1), bez fitovaného měřítka",
    "criterion": "MĚŘIDLO: dvě rozlišení shodná (|Δg| <= conv_g, |A1/A2 − 1| <= conv_amp) a reziduum rovnice <= res_tol q. "
                 "S (tvar): |g − g_pred| <= g_tol (čtvrtina rozestupu kandidátních zákonů). "
                 "M (síla): |A/A_pred − 1| <= amp_tol (přesnost čtení Π, K.3). PROŠEL: S a M. NEPROŠEL: S nebo M neplatí.",
    "null_model": "deterministická kvadratura bez vzorkování; konkurenční hypotézy jsou zákony 1/r, 1/r², 1/r³ a ln r "
                  "(g = −1, −2, −3, 0), kritérium S je odděluje",
    "level": "mechanism",
    "meters": ["K.4"],
    "depends_on": ["K.4"],
    "params": {"sigma": 1.0, "q": 1.0, "R": 80.0, "r_lo": 10.0, "r_hi": 40.0,
               "g_tol": 0.25, "amp_tol": 0.10, "conv_g": 0.05, "conv_amp": 0.02, "res_tol": 0.05,
               # resolution = [inner step dr, geometric ratio, rows per element, n_rho, n_ang]
               "res1": [0.25, 1.06, 3, 300, 32], "res2": [0.2, 1.04, 4, 400, 40],
               "probe_r": [10.0, 20.0, 40.0], "chunk": 16},
    "parts": {"d3": {"D": 3}, "d4": {"D": 4}},
    "quick": {"R": 24.0, "r_lo": 5.0, "r_hi": 12.0, "res1": [0.5, 1.15, 2, 150, 16], "res2": [0.4, 1.12, 2, 200, 20],
              "probe_r": [5.0, 8.0, 12.0]},
    "estimate": "na dimenzi ~7 min na 1 jádře (res1 ~1,5 min, res2 ~5 min); s více procesy úměrně méně",
}

# fixed units: eps * rho = 1, so l_k = 1 (any eps works; this one matches K.4)
_EPS = 0.05
_RHO = 1.0 / _EPS


def omega(n: int) -> float:
    """Surface area of the unit sphere in n spatial dimensions (2 pi for n = 2, 4 pi for n = 3)."""
    return 2.0 * math.pi ** (n / 2.0) / math.gamma(n / 2.0)


def total_source(n: int, q: float, sigma: float) -> float:
    """Q = integral of q exp(-r^2 / 2 sigma^2) over n-dimensional space = q (2 pi sigma^2)^(n/2)."""
    return q * (2.0 * math.pi * sigma * sigma) ** (n / 2.0)


def continuum_gaussian(n: int, q: float, sigma: float, R: float, r: np.ndarray) -> np.ndarray:
    """
    The continuum reference: lap u = q exp(-r^2/2 sigma^2), u(R) = 0, regular at 0.
    Gauss law per shell: u'(r) = Q(r) / (omega_n r^(n-1)), Q(r) the source inside r;
    u(r) = -integral_r^R u'(s) ds, by the trapezoid rule on a fine grid.
    """
    s = np.linspace(0.0, R, 200001)
    dens = q * np.exp(-s ** 2 / (2.0 * sigma ** 2)) * omega(n) * s ** (n - 1)
    Qs = np.concatenate([[0.0], np.cumsum(0.5 * (dens[1:] + dens[:-1]) * np.diff(s))])
    du = np.zeros_like(s)
    du[1:] = Qs[1:] / (omega(n) * s[1:] ** (n - 1))
    U = np.concatenate([[0.0], np.cumsum(0.5 * (du[1:] + du[:-1]) * np.diff(s))])
    u_grid = U - U[-1]                                  # u(R) = 0
    return np.interp(np.asarray(r, float), s, u_grid)


def rows_chunk(D, t_x, rows, R, nodes, n_rho, n_ang):
    """Rows of the expected operator for a chunk of row radii (module level: runs in worker processes)."""
    return expected_radial_rows(D, _RHO, _EPS, t_x, rows, R, nodes, n_rho=n_rho, n_ang=n_ang)


def solve_field(p, res, ctx, label):
    """Nodal field u for one resolution; returns nodes, u and the RMS residual in units of q."""
    D, R, sig, q = int(p["D"]), float(p["R"]), float(p["sigma"]), float(p["q"])
    dr, ratio, per, n_rho, n_ang = float(res[0]), float(res[1]), int(res[2]), int(res[3]), int(res[4])
    nodes = shape_nodes(R, 2.0 * sig, dr, ratio)        # uniform up to max(4, 4 sigma), then geometric
    fr = (np.arange(per) + 0.5) / per
    rows = np.sort((nodes[:-1, None] + fr[None, :] * np.diff(nodes)[:, None]).ravel())
    c = int(p["chunk"])
    chunks = [rows[i:i + c] for i in range(0, rows.size, c)]
    t0 = time.time()
    parts = parallel_map(rows_chunk, [(D, 3.0 * R, ch, R, nodes, n_rho, n_ang) for ch in chunks], ctx.workers,
                         ctx.progress, lambda i, r_, d, n_: f"[F1.3] {D-1}+1 {label}: rows {d}/{n_} chunks ({time.time()-t0:.0f}s)")
    A = np.vstack(parts)
    rhs = q * np.exp(-rows ** 2 / (2.0 * sig ** 2))
    sol, *_ = np.linalg.lstsq(A[:, :-1], rhs, rcond=None)
    u = np.r_[sol, 0.0]
    resid = float(np.sqrt(np.mean((A @ u - rhs) ** 2))) / q
    return nodes, u, resid, int(rows.size)


def amplitude_fixed(nodes, u, window, g):
    """Amplitude A and offset B of u = A boxcox(r, g) + B in the window (g fixed)."""
    m = (nodes >= window[0]) & (nodes <= window[1])
    M = np.column_stack([boxcox(nodes[m], g), np.ones(int(m.sum()))])
    coef, *_ = np.linalg.lstsq(M, u[m], rcond=None)
    return float(coef[0]), float(coef[1])


def run(p, ctx):
    D = int(p["D"])
    n = D - 1
    R, sig, q = float(p["R"]), float(p["sigma"]), float(p["q"])
    window = [float(p["r_lo"]), float(p["r_hi"])]
    g_pred = -1.0 if n == 3 else 0.0                    # Green function: 1/r (n = 3), ln r (n = 2)
    Q = total_source(n, q, sig)
    A_pred = Q / omega(n)
    out = {}
    for label in ("res1", "res2"):
        t0 = time.time()
        nodes, u, resid, nrows = solve_field(p, p[label], ctx, label)
        fit = fit_profile(nodes, u, window)
        A_fix, B_fix = amplitude_fixed(nodes, u, window, g_pred)
        half = math.sqrt(window[0] * window[1])         # geometric midpoint of the window
        g_in = fit_profile(nodes, u, [window[0], half])["g"]
        g_out = fit_profile(nodes, u, [half, window[1]])["g"]
        uc = continuum_gaussian(n, q, sig, R, nodes)
        wn = (nodes >= window[0]) & (nodes <= window[1])
        dev = float(np.sqrt(((u[wn] - uc[wn]) ** 2).sum() / (uc[wn] ** 2).sum()))
        out[label] = {"nodes": nodes, "u": u, "residual": resid, "rows": nrows, "g": fit["g"], "A_fix": A_fix,
                      "B_fix": B_fix, "g_inner": g_in, "g_outer": g_out, "dev_continuum": dev, "u_continuum": uc}
        ctx.progress(f"[F1.3] {n}+1 {label}: g = {fit['g']:+.3f}, A/A_pred = {A_fix / A_pred:.4f}, "
                     f"residual {resid:.1e} q ({time.time()-t0:.0f}s)")
    r1, r2 = out["res1"], out["res2"]
    dg = abs(r1["g"] - r2["g"])
    damp = abs(r1["A_fix"] / r2["A_fix"] - 1.0)
    inst = {"convergence of g": dg <= p["conv_g"],
            "convergence of amplitude": damp <= p["conv_amp"],
            "field residual": max(r1["residual"], r2["residual"]) <= p["res_tol"]}
    g, ratio = r2["g"], r2["A_fix"] / A_pred
    S_ok = abs(g - g_pred) <= p["g_tol"]
    M_ok = abs(ratio - 1.0) <= p["amp_tol"]
    candidates = {0.0: "ln r", -1.0: "1/r", -2.0: "1/r^2", -3.0: "1/r^3"}
    nearest = candidates[min(candidates, key=lambda c: abs(g - c))]
    if not all(inst.values()):
        status = "INSTRUMENT_FAILED"
    elif S_ok and M_ok:
        status = "PASSED"
    else:
        status = "FAILED"
    # local relative deviation from the continuum at a few radii (reported only)
    probes = [float(x) for x in p["probe_r"]]
    loc = [float((np.interp(x, r2["nodes"], r2["u"]) - np.interp(x, r2["nodes"], r2["u_continuum"]))
                 / abs(np.interp(x, r2["nodes"], r2["u_continuum"]))) for x in probes]
    lines = [f"{n}+1 expected nonlocal operator, lengths in l_k: Gaussian source sigma = {sig:g}, q = {q:g}, Q = {Q:.4f}; "
             f"R = {R:g}; window [{window[0]:g}, {window[1]:g}]",
             f"  prediction: g = {g_pred:+.0f} ({candidates[g_pred]}), amplitude A_pred = Q / omega_{n} = {A_pred:.5f}",
             f"  res1 ({r1['rows']} rows): g = {r1['g']:+.3f}, A/A_pred = {r1['A_fix'] / A_pred:.4f}, residual {r1['residual']:.1e} q",
             f"  res2 ({r2['rows']} rows): g = {r2['g']:+.3f}, A/A_pred = {ratio:.4f}, residual {r2['residual']:.1e} q",
             f"  convergence |dg| = {dg:.3f} (tol {p['conv_g']}), |A1/A2 - 1| = {damp:.4f} (tol {p['conv_amp']})",
             f"  S shape: |g - g_pred| = {abs(g - g_pred):.3f} (tol {p['g_tol']}) -> {'met' if S_ok else 'NOT met'}; "
             f"nearest candidate law: {nearest}",
             f"  M strength: |A/A_pred - 1| = {abs(ratio - 1.0):.4f} (tol {p['amp_tol']}) -> {'met' if M_ok else 'NOT met'}",
             f"  reported only: g in [{window[0]:g}, {half:.1f}] = {r2['g_inner']:+.3f}, in [{half:.1f}, {window[1]:g}] = "
             f"{r2['g_outer']:+.3f}; deviation from the continuum solution (u(R) = 0) {r2['dev_continuum']:.4f}",
             "  local (u - u_continuum) / |u_continuum| at r = " + ", ".join(f"{x:g}: {v:+.4f}" for x, v in zip(probes, loc)),
             "instrument: " + "; ".join(f"{k}: {'OK' if v else 'FAILS'}" for k, v in inst.items())]
    data = {"Q": Q, "A_pred": A_pred, "g_pred": g_pred, "window": window, "instrument": inst,
            "shape_met": S_ok, "strength_met": M_ok, "nearest_law": nearest, "local_deviation": dict(zip(map(str, probes), loc)),
            "res1": r1, "res2": r2}
    return Result(status, lines, data)
