"""
F3.16b -- The retarded Green function of a 3+1 order and the field it carries, revision b (C0)
==============================================================================================

Why this revision exists (14 Nov 2026; card docs/C0_card_b.md)
-------------------------------------------------------------
F3.16 ended INSTRUMENT_FAILED by design errors of its criteria, not by the instrument (docs/findings.md):
- the largest of 10 shells of eps, taken with the SE of the chosen shell, was compared with a tolerance from the SE of
  another shell;
- a nonlinear row sat in the instrument criterion against the noise-free (A3);
- 4 sprinklings at rho 60;
- tolerances of the form 1.25 (|dev| + 2 SE) on about 45 quantities whose true value is near zero.

A bootstrap of the development rows showed that the last form passes data identical to the development in about 1 %
of cases. The late deviation of F3.16 was shown to be a fluctuation (diag/c0_late.py; the mean of the linear rule
obeys (A3) exactly, docs/derivations/c0_mean_exactness.md).

Model, kernel, readouts, predictions unchanged (F3.16; diag/c0_card_b.task = diag/c0_card.task through the bit-identical
diag/c0_fast.evolve_fast). rho 30 (16 sprinklings) and 60 (12).

Criteria (card docs/C0_card_b.md, section 4; MV8+ct)
---------------------------------------------------
INSTRUMENT (the predictions exact in expectation): every item, ie the eps of 10 shells, the normalisation and (A3) of
  the linear rows k 0..3 in every slab, at both densities, has |d| / sqrt(SE^2 + 0.002^2) < t_crit. t_crit is Student
  t with n - 1 degrees of freedom at family-wise alpha 1 % (Bonferroni over all items). Any item above:
  INSTRUMENT_FAILED.
PHYSICS, tolerance criteria: the statistic is M = max over slabs (and amplitudes for the conditions of B) of
  |d_j| + 2 SE_j. PASSED side M < tol; FAILED side max_j (|d_j| - 2 SE_j) > tol. The growth exponent: mean + 2 SE < tol.
  The tolerances (params) are tol = M_dev + z SE_M, with SE_M from a bootstrap of the development sprinklings
  (diag/c0_card_b.derive3) and z one-sided at family-wise 5 % over the 19 criteria.
PHYSICS, sign criteria: convergence k 0..3, the noise falling in every mode, the renormalization falling at A 1 and 2
  (SE by jackknife over sprinklings): difference - 2 SE > 0 PASSED, + 2 SE < 0 FAILED.
Verdict order: INSTRUMENT_FAILED -> FAILED -> PASSED -> UNDECIDED.

    python run.py F3.16b --workers N          (seed 0 decides)
"""
from __future__ import annotations

import hashlib
import json
import os
import time

import numpy as np

from core.order.sprinkle import derive_seed
from core.stats.parallel import parallel_map
from core.testing import Result
from phases.F3_structure import field1d as F
from diag import c0_card_b as B

TEST = {
    "id": "F3.16b",
    "phase": "F3",
    "title": "The retarded Green function of a 3+1 order (MV8 + counterterm), revision b: statistical criteria",
    "theory_part": "F3 structure, C0 (decisions 25, 26; revision after F3.16 INSTRUMENT_FAILED by design errors)",
    "question": "Does a layer kernel with a local counterterm, both derived from the continuum, carry the retarded Green "
                "function of 3+1 and the linear and nonlinear field on it: the instrument as its exact predictions, the "
                "field within fixed tolerances of the continuum, the nonlinear softening explained by the fluctuations, "
                "errors and noise falling with the density?",
    "method": "as F3.16; rho 30 (16) and 60 (12); instrument as a significance test, physics with bootstrap margins",
    "prediction": "no instrument item significant; the physics statistics below their tolerances; all signs positive",
    "criterion": "card docs/C0_card_b.md section 4 (tolerances in params)",
    "null_model": "(A1), (A3) exact in expectation; the continuum solutions; B",
    "level": "mechanism",
    "meters": [],
    "depends_on": [],
    "params": {"rhos": [30.0, 60.0], "n": {"30": 16, "60": 12}, "X": 7.0, "T": 6.0, "n_se": 2.0,
               "tol": {"B_0.5_30": 0.055, "B_0.5_60": 0.045, "B_1_30": 0.04, "B_1_60": 0.04, "B_2_30": 0.065,
                       "B_2_60": 0.04, "direct_30": 0.014, "direct_60": 0.011, "skew_30": 0.41, "skew_60": 0.31,
                       "kurt_30": 0.6, "kurt_60": 0.615, "p_30": 1.2, "p_60": 1.15, "lin_cont_60": 0.08,
                       "nl_cont_0.5_60": 0.075, "nl_cont_1_60": 0.1, "nl_cont_2_60": 0.25, "nl1_cont_60": 0.09}},
    "quick": {"rhos": [30.0, 60.0], "n": {"30": 2, "60": 2}, "X": 4.0, "T": 6.0},
    "estimate": "~145 s per sprinkling at rho 30 and ~650 s at rho 60 with 8 processes on the author's machine; about "
                "45 min with 8 workers. numba required. Cached in diag/out/F3_16b_cache (keyed by this file's hash).",
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
    fn = os.path.join("diag", "out", "F3_16b_cache", _sha(), hashlib.sha256(key.encode()).hexdigest()[:24] + ".json")
    if os.path.exists(fn):
        with open(fn) as fh:
            return json.load(fh)
    out = B.task(float(rho), int(derive_seed(base, "F3.16b", float(rho), int(k))), float(p["X"]), float(p["T"]), 0)
    out["kk"] = int(k)
    os.makedirs(os.path.dirname(fn), exist_ok=True)
    with open(fn + ".tmp", "w") as fh:
        json.dump(out, fh)
    os.replace(fn + ".tmp", fn)
    return out


def run(p, ctx):
    if F._nb is None:
        raise SystemExit("F3.16b needs numba (pip install numba)")
    t0 = time.time()
    base, nse = ctx.seed, float(p.get("n_se", 2.0))
    rhos = sorted(float(r) for r in p["rhos"])
    jobs = [(p, rho, k, base) for rho in sorted(rhos, reverse=True) for k in range(int(p["n"][f"{rho:g}"]))]
    res = parallel_map(run_task, jobs, ctx.workers, ctx.progress,
                       lambda i, r, d, nn: f"[F3.16b] {d}/{nn}: rho {r['rho']:g} k {r['kk']} N {r['N']} "
                                           f"{r['sec']:.0f}s ({time.time() - t0:.0f}s)")
    st, inst, passc, failc, info = B.verdict3(res, p["tol"], nse)
    lines = [f"F3.16b (C0 revision b: 3+1 Green function, MV8 + counterterm); X {p['X']:g}, T {p['T']:g}; rho "
             + ", ".join(f"{r:g} ({p['n'][f'{r:g}']} sprinklings)" for r in rhos) + f"; base seed {base}"]
    worst = sorted(info["inst"], key=lambda q: -q[3] / q[4])[:6]
    lines.append(f"instrument: {len(inst)} items, t_crit {worst[0][4]:.2f} / {worst[-1][4]:.2f} (rho 30 / 60 df); "
                 f"floor 0.002; the six largest |d|/SE_eff:")
    for label, d, se, z, tc in worst:
        lines.append(f"  {label}: d {d:+.4f}, SE {se:.4f}, z {z:.2f} (t_crit {tc:.2f}) {'OK' if z < tc else 'FAILS'}")
    lines.append("physics (M = max_j |d_j| + 2 SE_j; p = mean + 2 SE):")
    for key, rr in B.PHYS_KEYS:
        for rho in rr:
            k = f"{key}_{rho}"
            lines.append(f"  {k}: M {info[k]:.4f} (tolerance {p['tol'][k]:g}) passed {'yes' if passc[k] else 'no'}"
                         f"{', FAILED' if failc[k] else ''}")
    lines.append("sign criteria (difference +- SE):")
    for k, v in info.items():
        if isinstance(v, tuple):
            lines.append(f"  {k}: {v[0]:+.4f} +- {v[1]:.4f} passed {'yes' if passc[k] else 'no'}{', FAILED' if failc[k] else ''}")
    lines.append("instrument failures: " + (", ".join(k for k, v in inst.items() if not v) or "none"))
    lines.append("not passed: " + (", ".join(k for k, v in passc.items() if not v) or "none"))
    lines.append("failed: " + (", ".join(k for k, v in failc.items() if v) or "none"))
    lines.append(f"cost: {np.mean([r['sec'] for r in res]):.0f} s per sprinkling ({time.time() - t0:.0f}s)")
    return Result(st, lines, {"instrument": inst, "passed": passc, "failed": failc,
                              "info": {k: v for k, v in info.items() if k != "inst"}, "inst_items": info["inst"],
                              "results": res})
