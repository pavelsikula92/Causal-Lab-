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
diag/f68f_common.py; derivations in docs/derivations/f68f_dev.md, steps 1b and 1c)
  * the first-arrival ray (driven pendulum psi'' = -eps sin(p s + theta + psi), all rays by the shooting map with every
    fold resolved) and the plane-wave volume A_0 = 6 int sqrt(s(1-s) E1 E2 / E1(1)) ds (Jacobi fields, RK4);
  * the q^2 series A = A_0 + A_1 q^2 + A_2 q^4 + A_3 q^6 from the Taylor series of the first-arrival actions to order 8
    (truncated power series on RK4, reversion, Lagrange inversion);
  * the exact needle (lower envelopes of the ray families, no expansion) as the independent check of the series.

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
           k = 0.12, n = 384 against 2 mu_3 k^5 n^5 C2(p) / l2^4 from the series (C2 from eps 0.1, 0.2 by
           Richardson): relative difference of the p-integrated density <= 5e-3 [proposal; expected O(1/n^2) = 7e-6 and
           O(k^2) ~ 1e-3]: FAILED otherwise. Needs an 80-bit long double (Linux / WSL).
NOT INCLUDED, blocking the freeze: the author's criterion "agreement with brute force at moderate eps (2, 4) <= 1e-3"
and the two levels kappa = 10 sqrt(eps), 20 sqrt(eps). Both need the exact volume at finite kappa and b in the real
medium; that instrument is not built (a draft was removed untested on 16 Nov). The author decides whether the e3 part
(an independent code, but only at b -> 0) replaces it or the card waits for it.

Verdict: INSTRUMENT_FAILED if single, jacobi or series fails; PASSED if qconv and e3 pass; FAILED otherwise.
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
              "eps 0.1-4, p on the F6.8e3 panels, 128 phases, q = 1/2, 1/4, 1/8; b -> 0 against F6.8e3 at n = 384",
    "prediction": "the instrument returns the volume to the stated accuracy (calibration)",
    "criterion": "INSTRUMENT: one ray per phase; delta_J <= 1e-4 and series A_0 vs plane wave <= 1e-9; series vs exact "
                 "needle <= 1e-6 (at most 2 of 12 window failures). PASSED: delta_A(1/4 -> 1/8) <= 1e-3 at every point "
                 "and the b -> 0 density against F6.8e3 within 5e-3; FAILED otherwise.",
    "null_model": "none (calibration)",
    "level": "calibration",
    "meters": [],
    "depends_on": ["F6.8e3"],
    "params": {"eps": [0.1, 0.25, 0.5, 1.0, 2.0, 4.0], "p_panels": [[0.0, 12.0, 10], [12.0, 24.0, 6]], "n_theta": 128,
               "qs": [0.5, 0.25, 0.125], "tol_J": 1e-4, "tol_A0": 1e-9, "tol_series": 1e-6, "max_window_fail": 2,
               "tol_qconv": 1e-3, "tol_e3": 5e-3, "e3_k": 0.12, "e3_n": 384.0, "e3_eps": [0.1, 0.2],
               "series_points": [[1.0, 0.0, 0.7], [1.0, 2.0, 2.0], [1.0, 4.0, 4.0], [1.0, 8.0, 0.7], [2.0, 0.0, 2.0],
                                 [2.0, 2.0, 4.0], [2.0, 4.0, 0.7], [2.0, 8.0, 2.0], [4.0, 0.0, 4.0], [4.0, 2.0, 0.7],
                                 [4.0, 4.0, 2.0], [4.0, 8.0, 4.0]]},
    "parts": {"main": {}},
    "quick": {"eps": [1.0, 4.0], "p_panels": [[0.0, 12.0, 2], [12.0, 24.0, 1]], "n_theta": 8,
              "series_points": [[4.0, 2.0, 2.0]], "e3_eps": [0.1, 0.2], "e3_n": 96.0},
    "estimate": "~4-6 h with 2 processes (the series at 96 points of 128 phases, the plane-wave check, F6.8e3 at 16 "
                "directions in 80-bit); needs an 80-bit long double for part e3 (Linux / WSL)",
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
        pred = 2 * MU3 * k ** 5 * n ** 5 * prof / L2 ** 4
        meas = np.array([b["bulk"] for b in bulk])
        rel = abs(float((wp * meas).sum()) / float((wp * pred).sum()) - 1.0)
        ok_e3 = rel <= p["tol_e3"]
        lines.append(f"  e3: p-integrated F6.8e3 bulk / series prediction - 1 = {rel:+.2e} (<= {p['tol_e3']:g})")
        for x, a, b in zip(pn, meas, pred):
            lines.append(f"     p {x:6.3f}: F6.8e3 {a:+.6e}  series {b:+.6e}  ratio {a / b:.6f}")
    if not (ok_single and ok_jac and ok_series):
        status = "INSTRUMENT_FAILED"
    elif ok_q and ok_e3:
        status = "PASSED"
    else:
        status = "FAILED"
    lines.append(f"  cost {time.time() - t0:.0f} s")
    return Result(status, lines, {"points": pts, "series": ser, "jacobi": jac, "exact": ex, "nray": nray, "dJ": dJ,
                                  "dA0": dA0, "worst_series": worst, "fails": fails, "worst_q": worst_q, "e3_rel": rel})
