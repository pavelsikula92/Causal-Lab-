"""
K.18 -- The needle volume of the optical sector at finite amplitude (scaling limit): calibration (DRAFT, NOT FROZEN)
====================================================================================================================

Why this test exists (16 Nov 2026; the author's order of 22:49, point 1)
-----------------------------------------------------------------------
F6.8f (the optical response of the expected BD action at finite amplitude) stands on a new instrument: the volume of
the causal interval of a long, thin pair (a needle) in the optical sector n = exp(b cos(k.w + theta)), at finite
eps = b (kT)^2, in the scaling limit kappa = kT -> infinity, b -> 0 at fixed eps and p = k_par T, at q = k sigma. By the
roadmap every instrument used in a phase needs a calibration record in the range in which it is used. The working
criteria the author set for step 1b of F6.8f-dev (16 Nov, 17:42) are calibration criteria; they are stated here.

The instrument (code in diag, frozen with this card by md5: diag/f68f_1b.py, diag/f68f_1c.py, diag/f68f_1b_exact.py,
diag/f68f_1d.py, diag/f68f_layer.py, diag/f68f_common.py; derivations in docs/derivations/f68f_dev.md, steps 1b-1d)
  * the first-arrival ray (driven pendulum psi'' = -eps sin(p s + theta + psi), all rays by the shooting map with every
    fold resolved) and the plane-wave volume A_0 = 6 int sqrt(s(1-s) E1 E2 / E1(1)) ds (Jacobi fields, RK4);
  * the q^2 series A = A_0 + A_1 q^2 + A_2 q^4 + A_3 q^6 from the Taylor series of the first-arrival actions to order 8
    (truncated power series on RK4, reversion, Lagrange inversion);
  * the exact needle (lower envelopes of the ray families, no expansion) as the independent check of the series;
  * the brute force at finite kappa and b in the real 3+1 medium (diag/f68f_1d.py: exact rays of n(z) as graphs, lower
    envelopes, the interval volume integrated directly, no expansion in b, q or 1 / kappa) as the check of the
    scaling limit itself.

Range: eps = 0.1, 0.25, 0.5, 1, 2, 4 (the author's grid; below the first caustic, one ray at every phase), p on the
panels of F6.8e3 ([0, 24]), 128 phases theta_j = 2 pi j / 128, q = 1/2, 1/4, 1/8.

Parts and criteria (the author's working criteria of 16 Nov 17:42 where given; [proposal] where added)
--------------------------------------------------------------------------------------------------
  single   one ray x -> y at every phase of the range (the series needs it; the conjugate-point criterion
           |Delta eps_c| <= 0.01 of the author has no object below the first caustic): INSTRUMENT_FAILED otherwise.
  jacobi   delta_J = max over the range |A_0(N) - A_0(2N)| (N = 4000 RK steps) <= 1e-4, and the series A_0 against the
           plane-wave A_0 of step 1b <= 1e-9 [proposal]: INSTRUMENT_FAILED otherwise.
  series   the series against the exact needle at q = 1/2 and 1/4 on 12 phases (eps 1, 2, 4; p 0, 2, 4, 8; theta 0.7,
           2.0, 4.0 cycling): |A_exact - A_series| <= 1e-6 [proposal; the exact needle is accurate to ~3e-7]; a phase
           where the exact needle's window search fails is reported and not counted: INSTRUMENT_FAILED if more than
           2 of 12 fail or the tolerance is exceeded.
  qconv    delta_A = |mean A(1/8) - mean A(1/4)| / max(1, ...) <= 1e-3 at every (eps, p) of the range (series; the
           pairwise distribution over the phases is reported): FAILED otherwise.
  e3       b -> 0 against the F6.8e3 instrument (an independent code): the F6.8e3 bulk density per direction at
           k = 0.12, n = 1536 against 4 mu_3 k^5 n^5 C2(p) / l2^4 from the series (C2 from eps 0.1, 0.2 by
           Richardson): relative difference of the p-integrated density <= 5e-3 [proposal]: FAILED otherwise. Needs an
           80-bit long double (Linux / WSL). Development (docs/derivations/f68f_dev.md, step 1c): +2.9e-3 at n = 384,
           -5.5e-4 at n = 1536; a part c / n^2 (+3.7e-3 at n = 384) and a remainder of about -8e-4 that does not
           vanish with n (k^2 or quadrature, not separated). n = 1536 keeps the 1/n^2 part at 2e-4.
  brute    the author's criterion "agreement with brute force at moderate eps (2, 4) <= 1e-3" with the levels
           kappa = 10 sqrt(eps), 20 sqrt(eps) (b = 0.01, 0.0025): the finite-kappa volume A_kappa (diag/f68f_1d.py) at
           the 8 points of series_points with eps 2 and 4, q = 1/4, 1/8. The finite-kappa corrections are physical,
           O(b) = O(1 / kappa^2), so the comparable configuration is the extrapolation A_ext = (4 A_20 - A_10) / 3:
           |A_ext - A_series| <= 1e-3 [the author's number; the extrapolation is a proposal], and the O(1 / kappa^2)
           form is checked: (A_10 - A_series) / (A_20 - A_series) in [3.5, 4.5] where |A_10 - A_series| > 1e-4
           [proposal]: FAILED otherwise. The raw differences at both levels are reported. Development (step 1d): 10
           points at eps 2, 4, q 1/8 - 1/2: |A_ext - A_series| <= 3.7e-6, ratios 3.97-4.07 above 1e-4, raw differences up to
           3.7e-3 (kappa = 20, eps 4, p 4) and 9.4e-4 (kappa = 40).

Verdict: INSTRUMENT_FAILED if single, jacobi or series fails; PASSED if qconv, e3 and brute pass; FAILED otherwise.
"""

from __future__ import annotations

import math
import time

import numpy as np

from core.stats.parallel import parallel_map
from core.testing import Result

TEST = {
    "id": "K.18",
    "phase": "K",
    "title": "Needle volume of the optical sector at finite amplitude (scaling limit): calibration (draft, not frozen)",
    "theory_part": "instrument of F6.8f: the interval volume of near-cone pairs at finite eps = b (kT)^2",
    "question": "Does the needle-volume instrument (plane wave, q^2 series to q^6, exact needle) give the interval volume "
                "of the optical sector to the working accuracy in the range eps <= 4 that F6.8f uses?",
    "method": "first-arrival rays of the scaling limit (all rays, folds resolved); plane-wave volume from Jacobi fields; "
              "q^2 series from Taylor series of the first-arrival actions (order 8); exact needle by lower envelopes; "
              "eps 0.1-4, p on the F6.8e3 panels, 128 phases, q = 1/2, 1/4, 1/8; b -> 0 against F6.8e3 at n = 1536; "
              "brute force at finite kappa = 10, 20 sqrt(eps) in the real medium at eps 2, 4",
    "prediction": "the instrument returns the volume to the stated accuracy (calibration)",
    "criterion": "INSTRUMENT: one ray per phase; delta_J <= 1e-4 and series A_0 vs plane wave <= 1e-9; series vs exact "
                 "needle <= 1e-6 (at most 2 of 12 window failures). PASSED: delta_A(1/4 -> 1/8) <= 1e-3 at every point "
                 "the b -> 0 density against F6.8e3 within 5e-3, and the finite-kappa brute force extrapolated in 1/kappa^2 "
                 "within 1e-3 of the series (ratio of the level differences 3.5-4.5); FAILED otherwise.",
    "null_model": "none (calibration)",
    "level": "calibration",
    "meters": [],
    "depends_on": ["F6.8e3"],
    "params": {"eps": [0.1, 0.25, 0.5, 1.0, 2.0, 4.0], "p_panels": [[0.0, 12.0, 10], [12.0, 24.0, 6]], "n_theta": 128,
               "qs": [0.5, 0.25, 0.125], "tol_J": 1e-4, "tol_A0": 1e-9, "tol_series": 1e-6, "max_window_fail": 2,
               "tol_qconv": 1e-3, "tol_e3": 5e-3, "e3_k": 0.12, "e3_n": 1536.0, "e3_eps": [0.1, 0.2],
               "series_points": [[1.0, 0.0, 0.7], [1.0, 2.0, 2.0], [1.0, 4.0, 4.0], [1.0, 8.0, 0.7], [2.0, 0.0, 2.0],
                                 [2.0, 2.0, 4.0], [2.0, 4.0, 0.7], [2.0, 8.0, 2.0], [4.0, 0.0, 4.0], [4.0, 2.0, 0.7],
                                 [4.0, 4.0, 2.0], [4.0, 8.0, 4.0]],
               "brute_levels": [10.0, 20.0], "brute_qs": [0.25, 0.125], "tol_brute": 1e-3, "brute_ratio": [3.5, 4.5],
               "brute_min_diff": 1e-4},
    "parts": {"main": {}},
    "quick": {"eps": [1.0, 4.0], "p_panels": [[0.0, 12.0, 2], [12.0, 24.0, 1]], "n_theta": 8,
              "series_points": [[4.0, 2.0, 2.0]], "e3_eps": [0.1, 0.2], "e3_n": 384.0, "brute_qs": [0.25]},
    "estimate": "~5-7 h with 2 processes (the series at 96 points of 128 phases, the plane-wave check, F6.8e3 at 16 "
                "directions in 80-bit at n = 1536 ~1 h, the brute force 32 volumes ~15 min); needs an 80-bit long double "
                "for part e3 (Linux / WSL)",
}


def _jacobi(eps, p, n_th):
    from diag import f68f_1b as B
    r = B.point(eps, p, 4000, n_th)
    return {"dA": r["dA"], "nbr_max": r["nbr_max"], "A": r["A"]}


def _series(eps, p, n_th):
    from diag import f68f_1c as S
    ths = 2 * math.pi * np.arange(n_th) / n_th
    r = S.series_coefficients(eps, p, ths)
    from diag import f68f_1b as B
    v, _, _ = S.rays(eps, p, ths)
    A0pw, _, _, _ = B.volume_A(v, eps, p, ths, 4000)
    return {"A": r["A"].tolist(), "nray": r["nray"], "lin": r["lin"],
            "dA0": float(np.max(np.abs(r["A"][:, 0] - A0pw)))}


def _exact(eps, p, th):
    from diag import f68f_1b_exact as X
    try:
        return {"A": X.phase_value(eps, p, th)["A"], "ok": True}
    except RuntimeError as err:
        return {"A": None, "ok": False, "why": str(err)}


def _brute(eps, p, th, kap, q):
    from diag import f68f_1d as D
    return D.volume(eps, p, th, kap, q)


def _e3(k, n, p):
    from diag import f68f_layer as L
    return L.task(k, n, p, 16, True)


def run(p, ctx):
    t0 = time.time()
    from phases.F6_order_dynamics.F6_08e3_bd_action_optical_asymptotic_rev2 import _p_nodes
    from diag import f68f_1c as S
    pn, wp = _p_nodes(p["p_panels"])
    pts = [(float(e), float(x)) for e in p["eps"] for x in pn]
    nt = int(p["n_theta"])
    lines = [f"K.18 (draft): needle volume; eps {p['eps']}, {len(pn)} directions, {nt} phases"]
    ser = parallel_map(_series, [(e, x, nt) for e, x in pts], ctx.workers, ctx.progress,
                       lambda i, r, d, n: f"[K.18] series {d}/{n}")
    jac = parallel_map(_jacobi, [(e, x, nt) for e, x in pts], ctx.workers, ctx.progress,
                       lambda i, r, d, n: f"[K.18] plane wave {d}/{n}")
    ex = parallel_map(_exact, [tuple(t) for t in p["series_points"]], ctx.workers, ctx.progress,
                      lambda i, r, d, n: f"[K.18] exact needle {d}/{n}")
    sp = [tuple(t) for t in p["series_points"]]
    ser_at = [S.series_coefficients(e, x, [th]) for (e, x, th) in sp]          # the series at the exact phases
    qs = p["qs"]
    # single, jacobi
    nray = max(r["nray"] for r in ser)
    dJ = max(r["dA"] for r in jac)
    dA0 = max(r["dA0"] for r in ser)
    ok_single = nray == 1
    ok_jac = dJ <= p["tol_J"] and dA0 <= p["tol_A0"]
    lines.append(f"  single: largest number of rays {nray} ({'ok' if ok_single else 'FAIL'})")
    lines.append(f"  jacobi: delta_J {dJ:.2e} (<= {p['tol_J']:g}); series A0 vs plane wave {dA0:.2e} (<= {p['tol_A0']:g})")
    # series vs exact
    worst, fails = 0.0, 0
    for (e, x, th), r, s_ in zip(sp, ex, ser_at):
        A = s_["A"][0]
        if not r["ok"]:
            fails += 1
            lines.append(f"    exact needle failed at eps {e} p {x} theta {th}: {r['why']}")
            continue
        for q, a in zip(qs[:2], r["A"][:2]):
            d = abs(a - (A[0] + A[1] * q ** 2 + A[2] * q ** 4 + A[3] * q ** 6))
            worst = max(worst, d)
    ok_series = worst <= p["tol_series"] and fails <= p["max_window_fail"]
    lines.append(f"  series vs exact needle: largest difference {worst:.2e} (<= {p['tol_series']:g}), window failures "
                 f"{fails} (<= {p['max_window_fail']})")
    # qconv
    worst_q, worst_pair = 0.0, 0.0
    for (e, x), r in zip(pts, ser):
        A = np.array(r["A"])
        a4 = A[:, 0] + A[:, 1] * qs[1] ** 2 + A[:, 2] * qs[1] ** 4 + A[:, 3] * qs[1] ** 6
        a8 = A[:, 0] + A[:, 1] * qs[2] ** 2 + A[:, 2] * qs[2] ** 4 + A[:, 3] * qs[2] ** 6
        dq = abs(a8.mean() - a4.mean()) / max(1.0, abs(a8.mean()), abs(a4.mean()))
        worst_q = max(worst_q, dq)
        worst_pair = max(worst_pair, float(np.max(np.abs(a8 - a4))))
    ok_q = worst_q <= p["tol_qconv"]
    lines.append(f"  qconv: largest delta_A {worst_q:.2e} (<= {p['tol_qconv']:g}); largest pairwise {worst_pair:.2e}")
    # e3
    k, n = float(p["e3_k"]), float(p["e3_n"])
    ok_e3, rel = True, float("nan")
    if np.finfo(np.longdouble).nmant < 63:
        lines.append("  e3: no 80-bit long double on this platform: not run (Linux / WSL)")
        ok_e3 = False
    else:
        bulk = parallel_map(_e3, [(k, n, float(x)) for x in pn], ctx.workers, ctx.progress,
                            lambda i, r, d, nn: f"[K.18] F6.8e3 bulk {d}/{nn}")
        e1, e2 = (float(v) for v in p["e3_eps"])
        small = parallel_map(_series, [(e_, float(x), nt) for x in pn for e_ in (e1, e2)], ctx.workers, ctx.progress,
                             lambda i, r, d, nn: f"[K.18] small eps {d}/{nn}")
        prof = []
        for i in range(len(pn)):
            c = []
            for r_, e_ in ((small[2 * i], e1), (small[2 * i + 1], e2)):
                A = np.array(r_["A"])
                c.append(float(np.mean(np.array([S.s4hat(a) for a in A]) * A[:, 0] ** -2.0)) / e_ ** 2)
            prof.append((4 * c[0] - c[1]) / 3 if abs(e2 - 2 * e1) < 1e-12 else c[0])
        prof = np.array(prof)
        L2, MU3 = S.L2, S.MU3
        pred = 4 * MU3 * k ** 5 * n ** 5 * prof / L2 ** 4
        meas = np.array([b["bulk"] for b in bulk])
        rel = abs(float((wp * meas).sum()) / float((wp * pred).sum()) - 1.0)
        ok_e3 = rel <= p["tol_e3"]
        lines.append(f"  e3: p-integrated F6.8e3 bulk / series prediction - 1 = {rel:+.2e} (<= {p['tol_e3']:g})")
        for x, a, b in zip(pn, meas, pred):
            lines.append(f"     p {x:6.3f}: F6.8e3 {a:+.6e}  series {b:+.6e}  ratio {a / b:.6f}")
    # brute force at finite kappa (the series points with eps 2 and 4)
    bp = [(e, x, th) for (e, x, th) in sp if e in (2.0, 4.0)]
    levs = [float(v) for v in p["brute_levels"]]
    bjobs = [(e, x, th, L * math.sqrt(e), float(q)) for (e, x, th) in bp for q in p["brute_qs"] for L in levs]
    bres = parallel_map(_brute, bjobs, ctx.workers, ctx.progress, lambda i, r, d, nn: f"[K.18] brute force {d}/{nn}")
    ok_b, worst_b, rat_bad = True, 0.0, []
    brute_rows = []
    for i in range(0, len(bjobs), 2):
        e, x, th, _, q = bjobs[i]
        A10, A20 = bres[i], bres[i + 1]
        s_ = S.series_coefficients(e, x, [th])["A"][0]
        ser_v = s_[0] + s_[1] * q ** 2 + s_[2] * q ** 4 + s_[3] * q ** 6
        ext = (4 * A20 - A10) / 3
        d10, d20 = A10 - ser_v, A20 - ser_v
        ratio = d10 / d20 if d20 != 0 else float("inf")
        worst_b = max(worst_b, abs(ext - ser_v))
        if abs(d10) > p["brute_min_diff"] and not (p["brute_ratio"][0] <= ratio <= p["brute_ratio"][1]):
            rat_bad.append((e, x, th, q, ratio))
        brute_rows.append({"eps": e, "p": x, "theta": th, "q": q, "A10": A10, "A20": A20, "series": ser_v, "ext": ext})
        lines.append(f"     brute eps {e:g} p {x:g} theta {th:g} q {q:g}: kappa-levels - series {d10:+.2e} {d20:+.2e} "
                     f"(ratio {ratio:.3f}); extrapolated - series {ext - ser_v:+.2e}")
    ok_b = worst_b <= p["tol_brute"] and not rat_bad
    lines.append(f"  brute: largest |A_ext - A_series| {worst_b:.2e} (<= {p['tol_brute']:g}); ratios outside "
                 f"{p['brute_ratio']}: {len(rat_bad)}")
    if not (ok_single and ok_jac and ok_series):
        status = "INSTRUMENT_FAILED"
    elif ok_q and ok_e3 and ok_b:
        status = "PASSED"
    else:
        status = "FAILED"
    lines.append(f"  cost {time.time() - t0:.0f} s")
    return Result(status, lines, {"points": pts, "series": ser, "jacobi": jac, "exact": ex, "nray": nray, "dJ": dJ,
                                  "dA0": dA0, "worst_series": worst, "fails": fails, "worst_q": worst_q, "e3_rel": rel,
                                  "brute": brute_rows, "worst_brute": worst_b})
