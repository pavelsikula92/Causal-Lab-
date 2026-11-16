"""
F6.8f -- The optical age term of the expected BD action at finite amplitude (DRAFT, NOT FROZEN)
===============================================================================================

Why this test exists (16 Nov 2026; the author's development branch of 09:59 and the order of 22:49, point 3)
---------------------------------------------------------------------------------------------------------
F6.8e3 (PASSED) found the optical age term of the second variation growing as (kT)^5. The second variation holds only
for eps^2 = b^2 (kT)^4 << 1; static fields of physical interest have eps^2 ~ 1e40-1e55. Whether the expected BD action
gives a usable static response for an old field is decided by the full dependence on eps, not by the b^2 term.

What the development found (docs/derivations/f68f_dev.md, steps 0-1c; decides nothing)
-------------------------------------------------------------------------------------
  * dE[S]/dT at finite b is the integral over the events of the past of x at time distance T (first arrival <= T).
    In the scaling limit (kappa = kT -> infinity, b -> 0 at fixed eps, p) only near-cone pairs (needles) count, and per
    direction the density is (T / 2 l2) int ds <g0(s^2 A(eps, p, theta; q))>_theta, A the needle volume (K.18).
  * The BD moments mu_0 = mu_1 = mu_2 = 0 remove every order below q^6 (q = k sigma); the survivor is
        dE/dn ~ k^5 Phi(eps),   Phi(eps) = int dp < s4hat A_0^-2 >_theta,
    s4hat from the q^2 series of A to q^6 (K.18). At small eps Phi = C2 eps^2: the (kT)^5 law of F6.8e3.
  * E(n) ~ (1 / 2 sqrt b) int_0^{b n^2} eps^-1/2 Phi(eps) d eps: an old field has a finite (static) age term if and only
    if eps^-1/2 Phi(eps) is integrable at infinity, i.e. Phi decays faster than eps^-1/2.

What this card can and cannot decide
------------------------------------
The instrument (K.18) is calibrated for eps <= 4 (one ray per phase; the first caustic is at eps = pi^2 for p = 0).
On that range the card decides whether the second-order extrapolation still describes the age term (the transition
region). It cannot decide the static limit (readings A, B, C of the author's message of 09:59), which is a property
of Phi at eps -> infinity; that needs the multi-ray range (an extension K.18b of the instrument with the cut loci,
and a card F6.8g). The readings are recorded as indications only.

Question, quantities, criteria (PROPOSALS for the author's decision; nothing is frozen)
-------------------------------------------------------------------------------------
  N(eps) = Phi(eps) / (C2 eps^2), C2 from eps = 0.1, 0.2 by Richardson; judged eps = 0.25, 1, 2, 4 (development rows
  were computed at eps = 0.1, 0.2, 0.5 and 3 only, so the judged rows are not seen).
  INSTRUMENT_FAILED: K.18 not PASSED; or the b -> 0 check against F6.8e3 (k = 0.12, n = 384, p-integrated density)
  beyond 5e-3; or more than one ray at some phase; or the linear term of Phi above 1e-8.
  Hypothesis H2 (the second order holds on the transition range): PASSED if |N(eps) - 1| <= 0.05 at every judged eps;
  FAILED if |N(4) - 1| > 0.2 or Phi changes sign; UNDECIDED otherwise.
  Recorded (no verdict): the local exponent a(eps) = d ln|Phi| / d ln eps between judged neighbours, as an indication
  for the static limit (a < 2: weaker than the second order; a < -1/2 would be needed asymptotically for reading A).
  Checks of the author's list (09:59): 1 b -> 0 against F6.8e3 (above); 2 zero age (Phi(eps -> 0) -> 0); 3 k = 0
  (A = 1, Phi = 0); 4 both signs of b (Phi even in eps; exact: A(-eps, theta) = A(eps, theta + pi)); 5 conformal and
  mixed sectors (not computed here; unchanged); 6 the stable instrument (K.18).
Level: mechanism (the age term of the expected BD action in a declared model: a static plane-wave optical index).
"""

from __future__ import annotations

import math
import time

import numpy as np

from core.stats.parallel import parallel_map
from core.testing import Result

TEST = {
    "id": "F6.8f",
    "phase": "F6",
    "title": "Optical age term of the expected BD action at finite amplitude eps = b (kT)^2 (draft, not frozen)",
    "theory_part": "the static response of the order for old fields (B6-full); input of the static limit and of C.1",
    "question": "Does the second-order extrapolation describe the optical age term of the expected BD action up to "
                "eps = b (kT)^2 = 4, and how does the full term depart from it?",
    "method": "scaling limit of the age derivative: needle volumes A(eps, p, theta; q) as a q^2 series to q^6 (K.18), the "
              "BD moments select the q^6 survivor; Phi(eps) = int dp <s4hat A0^-2>_theta on the p panels of F6.8e3, 128 "
              "phases; N(eps) = Phi / (C2 eps^2) at eps = 0.25, 1, 2, 4",
    "prediction": "none fixed: the development computed eps = 0.1, 0.2, 0.5, 3 only; the null model is the second order",
    "criterion": "INSTRUMENT: K.18 PASSED; b -> 0 against F6.8e3 within 5e-3; one ray per phase; linear term <= 1e-8. "
                 "PASSED |N(eps) - 1| <= 0.05 at every judged eps; FAILED |N(4) - 1| > 0.2 or a sign change of Phi; "
                 "UNDECIDED otherwise (proposal, not frozen)",
    "null_model": "the second-order extrapolation N(eps) = 1 (F6.8e3 extended by eps^2)",
    "level": "mechanism",
    "meters": ["K.18"],
    "depends_on": ["F6.8e3", "K.18"],
    "params": {"eps": [0.25, 1.0, 2.0, 4.0], "c2_eps": [0.1, 0.2], "p_panels": [[0.0, 12.0, 10], [12.0, 24.0, 6]],
               "n_theta": 128, "tol_pass": 0.05, "tol_fail": 0.2, "tol_lin": 1e-8},
    "parts": {"main": {}},
    "quick": {"eps": [0.5, 3.0], "p_panels": [[0.0, 12.0, 2], [12.0, 24.0, 1]], "n_theta": 8},   # development eps only
    "estimate": "~2 h with 2 processes (the series at 6 eps x 16 directions x 128 phases)",
}


def _point(eps, p, n_th):
    from diag import f68f_1c as S
    r = S.point(eps, p, n_th)
    return {"phi": r["phi"], "nray": r["nray"], "lin": r["lin"], "A0": r["A0"]}


def run(p, ctx):
    t0 = time.time()
    from phases.F6_order_dynamics.F6_08e3_bd_action_optical_asymptotic_rev2 import _p_nodes
    pn, wp = _p_nodes(p["p_panels"])
    nt = int(p["n_theta"])
    epss = [float(e) for e in p["c2_eps"]] + [float(e) for e in p["eps"]]
    jobs = [(e, float(x), nt) for e in epss for x in pn]
    res = parallel_map(_point, jobs, ctx.workers, ctx.progress, lambda i, r, d, n: f"[F6.8f] {d}/{n}")
    Phi = {}
    nray, lin = 1, 0.0
    for e in epss:
        R = [r for j, r in zip(jobs, res) if j[0] == e]
        Phi[e] = float((wp * np.array([r["phi"] for r in R])).sum())
        nray = max(nray, max(r["nray"] for r in R))
        lin = max(lin, max(r["lin"] for r in R))
    e1, e2 = (float(v) for v in p["c2_eps"])
    c1, c2_ = Phi[e1] / e1 ** 2, Phi[e2] / e2 ** 2
    C2 = (4 * c1 - c2_) / 3 if abs(e2 - 2 * e1) < 1e-12 else c1
    lines = [f"F6.8f (draft): optical age term at finite eps; C2 = {C2:+.10e} (from eps {e1:g}, {e2:g})",
             "     eps      Phi(eps)          N(eps)"]
    N = {}
    for e in p["eps"]:
        e = float(e)
        N[e] = Phi[e] / (C2 * e * e)
        lines.append(f"  {e:6.3f}  {Phi[e]:+.10e}  {N[e]:+.6f}")
    js = sorted(N)
    for a, b in zip(js[:-1], js[1:]):
        if Phi[a] * Phi[b] > 0:
            lines.append(f"  local exponent {a:g} -> {b:g}: {math.log(abs(Phi[b] / Phi[a])) / math.log(b / a):+.4f}")
    inst_ok = nray == 1 and lin <= p["tol_lin"]
    lines.append(f"  instrument: largest number of rays {nray}; largest linear term {lin:.1e}; (K.18 and the F6.8e3 "
                 f"check are separate records)")
    signs = {np.sign(Phi[float(e)]) for e in p["eps"]} | {np.sign(C2)}
    if not inst_ok:
        status = "INSTRUMENT_FAILED"
    elif len(signs) > 1 or abs(N[max(N)] - 1) > p["tol_fail"]:
        status = "FAILED"
    elif all(abs(v - 1) <= p["tol_pass"] for v in N.values()):
        status = "PASSED"
    else:
        status = "UNDECIDED"
    lines.append(f"  cost {time.time() - t0:.0f} s")
    return Result(status, lines, {"Phi": {str(k): v for k, v in Phi.items()}, "C2": C2, "N": {str(k): v for k, v in N.items()},
                                  "nray": nray, "lin": lin})
