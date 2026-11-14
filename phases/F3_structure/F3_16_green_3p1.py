"""
F3.16 -- The retarded Green function of a 3+1 order and the field it carries (C0)
================================================================================

Why this test exists (14 Nov 2026; card docs/C0_card_draft.md, revision 1; decisions 25 and 26)
---------------------------------------------------------------------------------------------
Every 3+1 step of the programme (C1 curved orders, B4 localized defects) needs an instrument that computes a field on
a 3+1 causal order from the order alone. In 3+1 the retarded Green function of the wave operator lies on the light
cone, G = delta(t - r)/(4 pi r). C0 asks whether a fixed layer kernel w(n) (n = the number of events in the interval of
a pair) with a local counterterm, both derived from the continuum and not from any measurement, carries it and the
linear and nonlinear field of the continuum, with deviations that are the predicted finite-density effects and fall
with the density, and with noise that has no growing mode.

Model (frozen; diag/c0_card.py, diag/c0_apex.py, diag/c0_ct.py, diag/c0_ct_renorm.py, diag/c0_minvar.py,
diag/c0_links4.py; md5 in the decision that freezes the card)
-----------------------------------------------------------------------------------------------------------
3+1 pyramids |x|, |y|, |z| <= X - t, X 7, T 6, Poisson density rho (30 and 60). Pairs with rho V <= 22, classes
n <= 16. The rule, event by event in time order:
    phi(v) - delta J(phi(v)) = phi0(v) - c sum_{u < v} w(n_uv) J(phi(u)),   c = 1/(2 pi sqrt(6 rho)),
    delta = (a / (8 sqrt pi)) sum_n w_n H_n,  a = sqrt(24/(pi rho))     (solved per event; Newton for J = mu sin phi),
w = MV8 (the minimum-variance kernel with sum g_n w_n = Gamma(1/2), sum w_n = 0, n <= 8); MV2+ct on the same pairs as a
reference (reported, not deciding). mu 1. Data phi0 = A cos(k x) cos(k t).

Readouts (coordinates only here)
-------------------------------
  eps      max over shells r < 2.5 of |eps(r) - (A1)|;
  norm     the constant-source sum against its expectation at each event's own t, t >= 3;
  modes    amplitude deviation from the exact solution (projection on cos(kx), slabs |t - b| < 0.25 of the inner events):
           linear k = 0..3; nonlinear homogeneous A = 0.5, 1, 2 (pendulum); nonlinear A = 1 cos x (1+1 sine-Gordon);
           against A (the continuum: amp-dev), (A3) (the kernel's finite-density expectation, no noise) and B
           (J_eff = mu exp(-s^2/2) sin f with the measured s); residual = meas - B; renormalization = B - (A3);
  gauss    skewness, excess kurtosis (t >= 1.5) of the event noise and the direct check
           <sin phi> - exp(-s^2/2) <sin(phi - eta)>, nonlinear homogeneous rows;
  noise    the growth exponent p of the sd (t 1.5 .. 5.5); sd(rho 60)/sd(rho 30) at t 4.75.

Criteria (card revision 1, section 4; MV8+ct; every condition on its own; SE over sprinklings)
INSTRUMENT_FAILED if any of, at each rho: eps + 2 SE >= tol_eps; |norm| + 2 SE >= tol_norm; max |amp-dev - (A3)|
  (linear k 0..3 and nonlinear A 0.5) + 2 SE >= tol_A3.
FAILED if any of: a tolerance criterion with value - 2 SE > tol; a sign criterion with difference + 2 SE < 0.
PASSED if all tolerance criteria have value + 2 SE < tol and all sign criteria difference - 2 SE > 0.
  tolerance criteria: linear (A) at rho 60; nonlinear (A) at rho 60 for A 0.5, 1, 2 and A 1 cos x; residual against B
  for A 0.5, 1, 2 at each rho; |direct check|, |skewness|, |excess kurtosis| at each rho; p of every mode at each rho;
  sign criteria: convergence max|A|(30) - max|A|(60), k 0..3; 1 - sd(60)/sd(30) at t 4.75, every mode;
  renormalization(30) - renormalization(60), A 1 and 2.
Otherwise UNDECIDED.

    python run.py F3.16 --workers N          (seed 0 decides)
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import time

import numpy as np

from core.order.sprinkle import derive_seed
from core.stats.parallel import parallel_map
from core.testing import Result
from phases.F3_structure import field1d as F
from diag import c0_card as C

TEST = {
    "id": "F3.16",
    "phase": "F3",
    "title": "The retarded Green function of a 3+1 order (layer kernel MV8 + local counterterm) and the field it carries",
    "theory_part": "F3 structure, C0 (decisions 25, 26): the 3+1 Green function from the counts in intervals of the "
                   "order, with a local term fixed by the continuum; instrument for C1 and B4",
    "question": "Does a layer kernel with a local counterterm, both derived from the continuum, carry the retarded Green "
                "function of 3+1 and the linear and nonlinear field on it, with deviations that are the predicted "
                "finite-density effects and fall with the density, the nonlinear softening explained by the fluctuations, "
                "and noise without a growing mode?",
    "method": "3+1 sprinklings X 7, T 6, rho 30 (16) and 60 (4); MV8+ct (MV2+ct reported); linear modes k 0-3, nonlinear "
              "homogeneous A 0.5/1/2 and A 1 cos x; comparisons A (continuum), (A3) (kernel expectation), B (J_eff)",
    "prediction": "eps and the linear modes as (A1)/(A3); the nonlinear modes as B; errors falling with rho; p ~ 1",
    "criterion": "card docs/C0_card_draft.md revision 1, section 4 (tolerances in params)",
    "null_model": "the continuum solutions (Klein-Gordon, pendulum, 1+1 sine-Gordon) and the kernel's expectation",
    "level": "mechanism",
    "meters": [],
    "depends_on": [],
    "params": {"rhos": [30.0, 60.0], "n": {"30": 16, "60": 4}, "X": 7.0, "T": 6.0, "n_se": 2.0,
               "tol": {"eps_30": 0.013, "norm_30": 0.011, "A3_30": 0.065, "B_0.5_30": 0.02, "B_1_30": 0.025,
                       "B_2_30": 0.065, "direct_30": 0.017, "skew_30": 0.34, "kurt_30": 0.69, "p_30": 1.45,
                       "eps_60": 0.008, "norm_60": 0.005, "A3_60": 0.045, "B_0.5_60": 0.015, "B_1_60": 0.025,
                       "B_2_60": 0.04, "direct_60": 0.012, "skew_60": 0.385, "kurt_60": 0.715, "p_60": 1.45,
                       "lin_cont_60": 0.09, "nl_cont_0.5_60": 0.085, "nl_cont_1_60": 0.105, "nl_cont_2_60": 0.285,
                       "nl1_cont_60": 0.09}},
    "quick": {"rhos": [8.0, 12.0], "n": {"8": 2, "12": 2}, "X": 4.0, "T": 3.0},
    "estimate": "~320 s per sprinkling at rho 30, ~2600 s at rho 60 (4 processes on the author's machine); about 1.1 h "
                "with 4 workers. numba required. Finished tasks are cached in diag/out/F3_16_cache (keyed by this "
                "file's hash).",
}

_SHA = None


def _sha():
    global _SHA
    if _SHA is None:
        with open(__file__, "rb") as fh:
            _SHA = hashlib.sha256(fh.read()).hexdigest()[:16]
    return _SHA


def run_task(p, rho, k, base):
    key = json.dumps([str(rho), k, base, p["X"], p["T"]])
    fn = os.path.join("diag", "out", "F3_16_cache", _sha(), hashlib.sha256(key.encode()).hexdigest()[:24] + ".json")
    if os.path.exists(fn):
        with open(fn) as fh:
            return json.load(fh)
    out = C.task(float(rho), int(derive_seed(base, "F3.16", float(rho), int(k))), float(p["X"]), float(p["T"]), 0)
    out["kk"] = int(k)
    os.makedirs(os.path.dirname(fn), exist_ok=True)
    with open(fn + ".tmp", "w") as fh:
        json.dump(out, fh)
    os.replace(fn + ".tmp", fn)
    return out


def run(p, ctx):
    if F._nb is None:
        raise SystemExit("F3.16 needs numba (pip install numba)")
    t0 = time.time()
    base, nse = ctx.seed, float(p.get("n_se", 2.0))
    rhos = sorted(float(r) for r in p["rhos"])
    lo, hi = rhos
    tol = p["tol"]
    jobs = [(p, rho, k, base) for rho in sorted(rhos, reverse=True) for k in range(int(p["n"][f"{rho:g}"]))]
    res = parallel_map(run_task, jobs, ctx.workers, ctx.progress,
                       lambda i, r, d, nn: f"[F3.16] {d}/{nn}: rho {r['rho']:g} k {r['kk']} N {r['N']} "
                                           f"{r['sec']:.0f}s ({time.time() - t0:.0f}s)")
    A = {rho: C.card_stats([r for r in res if r["rho"] == rho]) for rho in rhos}
    main = "MV8+ct"
    lines = [f"F3.16 (C0: 3+1 Green function, MV8 + counterterm; MV2+ct reported); X {p['X']:g}, T {p['T']:g}; rho "
             + ", ".join(f"{r:g} ({p['n'][f'{r:g}']} sprinklings)" for r in rhos) + f"; base seed {base}; mean +- SE"]
    inst, passc, failc = {}, {}, {}

    def tcrit(store_pass, store_fail, label, val, se, t):
        store_pass[label] = val + nse * se < t
        if store_fail is not None:
            store_fail[label] = val - nse * se > t
        lines.append(f"  {label}: {val:.4f} +- {se:.4f} (tolerance {t:g})")

    def scrit(label, d, se):
        passc[label] = d - nse * se > 0
        failc[label] = d + nse * se < 0
        lines.append(f"  {label}: {d:+.4f} +- {se:.4f} (sign)")

    for rho in rhos:
        g = f"{rho:g}" if rho in (30.0, 60.0) else ("30" if rho == lo else "60")
        o = A[rho][main]
        ms = o["modes"]
        lines.append(f"rho {rho:g}, {main}:")
        tcrit(inst, None, f"eps rho {rho:g}", *o["eps"], tol[f"eps_{g}"])
        tcrit(inst, None, f"|norm - 1| rho {rho:g}", abs(o["norm"][0]), o["norm"][1], tol[f"norm_{g}"])
        a3 = max(((k, *ms[k]["A3"]) for k in ("lin0", "lin1", "lin2", "lin3", "nl0_A0.5")),
                 key=lambda q: q[1] + nse * q[2])
        tcrit(inst, None, f"max |amp-dev - (A3)| rho {rho:g} ({a3[0]})", a3[1], a3[2], tol[f"A3_{g}"])
        for A_ in ("0.5", "1", "2"):
            tcrit(passc, failc, f"residual |meas - B| A {A_} rho {rho:g}", *ms[f"nl0_A{A_}"]["B"], tol[f"B_{A_}_{g}"])
        for what in ("direct", "skew", "kurt"):
            q = max(((k, *ms[k][what]) for k in ("nl0_A0.5", "nl0_A1", "nl0_A2")), key=lambda z: z[1] + nse * z[2])
            tcrit(passc, failc, f"|{what}| rho {rho:g} ({q[0]})", q[1], q[2], tol[f"{what}_{g}"])
        q = max(((m[0], *ms[m[0]]["p"]) for m in C.MODES), key=lambda z: z[1] + nse * z[2])
        tcrit(passc, failc, f"growth exponent p rho {rho:g} ({q[0]})", q[1], q[2], tol[f"p_{g}"])
        for m in C.MODES:
            mo = ms[m[0]]
            lines.append(f"     {m[0]:9s} max|A| {mo['cont'][0]:.3f}" + (f"  max|meas-(A3)| {mo['A3'][0]:.3f}" if "A3" in mo else "")
                         + (f"  max|meas-B| {mo['B'][0]:.3f}  renorm {mo['renorm']:.3f}" if "B" in mo else "")
                         + f"  p {mo['p'][0]:.2f}  sd " + " ".join(f"{v:.2f}" for v in mo["sd"][1::2]))
        r2 = A[rho]["MV2+ct"]
        lines.append("     reported, MV2+ct: max|A| " + ", ".join(f"{m[0]} {r2['modes'][m[0]]['cont'][0]:.3f}" for m in C.MODES)
                     + f"; eps {r2['eps'][0]:.4f}; norm {r2['norm'][0]:+.4f}")
    o = A[hi][main]["modes"]
    g = "60"
    q = max(((k, *o[k]["cont"]) for k in ("lin0", "lin1", "lin2", "lin3")), key=lambda z: z[1] + nse * z[2])
    tcrit(passc, failc, f"linear (A) max |amp-dev| rho {hi:g} ({q[0]})", q[1], q[2], tol[f"lin_cont_{g}"])
    for A_ in ("0.5", "1", "2"):
        tcrit(passc, failc, f"nonlinear (A) A {A_} rho {hi:g}", *o[f"nl0_A{A_}"]["cont"], tol[f"nl_cont_{A_}_{g}"])
    tcrit(passc, failc, f"nonlinear (A) A 1 cos x rho {hi:g}", *o["nl1_A1"]["cont"], tol[f"nl1_cont_{g}"])
    a30, a60 = A[lo][main]["modes"], A[hi][main]["modes"]
    for k in ("lin0", "lin1", "lin2", "lin3"):
        scrit(f"convergence {k}", a30[k]["cont"][0] - a60[k]["cont"][0], math.hypot(a30[k]["cont"][1], a60[k]["cont"][1]))
    for m in C.MODES:
        r, se = C.sd_ratio(A[lo], A[hi], main, m[0], tb=np.array(res[0]["t"]))
        scrit(f"noise falls {m[0]}", 1 - r, se)
    for k in ("nl0_A1", "nl0_A2"):
        scrit(f"renormalization falls {k}", a30[k]["renorm"] - a60[k]["renorm"], math.hypot(a30[k]["B"][1], a60[k]["B"][1]))
    inst = {k: bool(v) for k, v in inst.items()}
    passc = {k: bool(v) for k, v in passc.items()}
    failc = {k: bool(v) for k, v in failc.items()}
    lines.append("instrument: " + "; ".join(f"{k}: {'OK' if v else 'FAILS'}" for k, v in inst.items()))
    lines.append("passed conditions: " + "; ".join(f"{k}: {'yes' if v else 'no'}" for k, v in passc.items()))
    lines.append("failed conditions: " + "; ".join(f"{k}: {'YES' if v else 'no'}" for k, v in failc.items()))
    if not all(inst.values()):
        status = "INSTRUMENT_FAILED"
    elif any(failc.values()):
        status = "FAILED"
    elif all(passc.values()):
        status = "PASSED"
    else:
        status = "UNDECIDED"
    lines.append(f"cost: {np.mean([r['sec'] for r in res]):.0f} s per sprinkling ({time.time() - t0:.0f}s)")
    return Result(status, lines, {"instrument": inst, "passed": passc, "failed": failc, "results": res})
