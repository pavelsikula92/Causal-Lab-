"""
F5.3 -- The field of the C0 instrument on a curved 3+1 order: the exact mean, not a local field equation (C1)
==========================================================================================================

Why this test exists (14 Nov 2026; card docs/C1_card.md)
-------------------------------------------------------
The C0 instrument (MV8 + counterterm; F3.16b) reads only the order and the counts in intervals. On a curved order the
counts measure the curved volume. Development (docs/findings.md, "C1 development rows") found:
- on a conformally flat planar well the field follows the exact mean of the rule (Mecke's formula, the diamond-averaged
  curved volume; diag/c1_apex.py), with no free parameter;
- that mean converges with the density to a Green function in which a near-null pair weighs rho Omega_u^4
  <Omega^4>_I^(-1/2). This is not the conformally coupled (or any local xi R) Green function: the curvature along the
  whole null segment enters.

This test fixes the claim on fresh sprinklings, against both predictions computed before the run by frozen code.

Model (frozen; diag/c1_conf.py, diag/c1_apex.py, diag/c0_fast.py, diag/c0_apex.py, diag/c0_minvar.py, diag/c0_links4.py)
-------------------------------------------------------------------------------------------------------------------
g = Omega(x)^2 eta, Omega = 1 - 0.15 exp(-x^2/(2 1.5^2)). The flat cones make the pyramid |x|, |y|, |z| <= 7 - t past-complete,
T 6. A flat sprinkling of proper density rho is paired with its curved child, thinned with probability Omega^4. The C0
rule runs unchanged on both: MV8 kernel, c and the counterterm delta at the proper rho, linear J = mu phi, mu 1,
data phi0 = A. rho 30 and 60, 8 pairs each.

Readouts: D = <phi>_curved - <phi>_flat in bins of t and x (width 0.5) over the inner events, mean and SE over the
pairs.
- P_exact: the exact mean (diag/c1_apex.mean_field, mode 0).
- P_local: the conformally coupled Green function Omega_u^-1 Omega_v^-1 G_eta, smeared by the same flat kernel
  (mode 1). This is the local continuum theory with the same finite-density smearing.

Criteria (card docs/C1_card.md; decision 28 form; Bonferroni family-wise alpha 1 %, Student t with n - 1 dof)
INSTRUMENT: the flat parents per time bin (t >= 0.5, |x| <= 2.5) against the flat (A3): |d| / sqrt(SE^2 + 0.002^2) < t_crit.
P1 (the exact mean): every bin t >= 1, |x| <= 2.5, both densities: |D - P_exact| / sqrt(SE^2 + 0.005^2) < t_crit. Any bin
  above: FAILED.
P2 (not the local theory): the centre window |x| < 0.75, t >= 4.5 at each density. Delta = mean over the pairs of
  (D - P_local). PASSED side Delta - 2 SE > 0, the direction and size predicted by P_exact - P_local; FAILED side
  Delta + 2 SE < 0.
Verdict order: INSTRUMENT_FAILED -> FAILED -> PASSED -> UNDECIDED.

    python run.py F5.3 --workers N          (seed 0 decides)
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import time

import numpy as np
from scipy import stats

from core.order.sprinkle import derive_seed
from core.stats.parallel import parallel_map
from core.testing import Result
from phases.F3_structure import field1d as F
from diag import c1_apex as AX
from diag import c1_conf as C1

TEST = {
    "id": "F5.3",
    "phase": "F5",
    "title": "The field of the C0 instrument on a curved 3+1 order: the exact mean, not a local field equation",
    "theory_part": "F5 gravity, C1: the 3+1 field rule (MV8 + counterterm) on a conformally flat curved order",
    "question": "On a curved 3+1 order (conformally flat planar well), does the field of the order-only rule follow its "
                "exact mean (the curved volume along each pair), and does it differ, in the predicted direction, from "
                "the local conformally coupled field with the same finite-density smearing?",
    "method": "paired flat/curved sprinklings, rho 30 and 60 (8 pairs each), X 7, T 6; D = <phi>_curved - <phi>_flat by "
              "(t, x); predictions by diag/c1_apex.mean_field (mode 0 exact, mode 1 local conformal)",
    "prediction": "every bin as the exact mean; the centre at t >= 4.5 above the local model by about +0.14",
    "criterion": "card docs/C1_card.md section 4",
    "null_model": "the local conformally coupled field with the same smearing (mode 1)",
    "level": "mechanism",
    "meters": [],
    "depends_on": ["F3.16b"],
    "params": {"rhos": [30.0, 60.0], "n": {"30": 8, "60": 8}, "X": 7.0, "T": 6.0, "alpha": 0.01, "floor": 0.005,
               "flat_floor": 0.002, "n_se": 2.0},
    "quick": {"rhos": [8.0, 12.0], "n": {"8": 2, "12": 2}, "X": 4.0, "T": 3.0},
    "estimate": "~2000 s per pair at rho 60 and ~460 s at rho 30 with 8 processes on the author's machine; about 45 min "
                "with 8 workers. numba required. Cached in diag/out/F5_3_cache (keyed by this file's hash).",
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
    fn = os.path.join("diag", "out", "F5_3_cache", _sha(), hashlib.sha256(key.encode()).hexdigest()[:24] + ".json")
    if os.path.exists(fn):
        with open(fn) as fh:
            return json.load(fh)
    s1 = int(derive_seed(base, "F5.3", float(rho), int(k)))
    s2 = int(derive_seed(base, "F5.3-thin", float(rho), int(k)))
    out = C1.task_seeds(float(rho), s1, s2, float(p["X"]), float(p["T"]), int(k))
    out["kk"] = int(k)
    os.makedirs(os.path.dirname(fn), exist_ok=True)
    with open(fn + ".tmp", "w") as fh:
        json.dump(out, fh)
    os.replace(fn + ".tmp", fn)
    return out


def run(p, ctx):
    if F._nb is None:
        raise SystemExit("F5.3 needs numba (pip install numba)")
    t0 = time.time()
    base, nse = ctx.seed, float(p.get("n_se", 2.0))
    rhos = sorted(float(r) for r in p["rhos"])
    jobs = [(p, rho, k, base) for rho in sorted(rhos, reverse=True) for k in range(int(p["n"][f"{rho:g}"]))]
    res = parallel_map(run_task, jobs, ctx.workers, ctx.progress,
                       lambda i, r, d, nn: f"[F5.3] {d}/{nn}: rho {r['rho']:g} k {r['kk']} "
                                           f"{r['flat']['sec'] + r['curved']['sec']:.0f}s ({time.time() - t0:.0f}s)")
    S = {rho: AX.card_stats([r for r in res if r["rho"] == rho]) for rho in rhos}
    lines = [f"F5.3 (C1: the C0 instrument on a conformally flat well, a {C1.A_WELL}, L {C1.L_WELL}); X {p['X']:g}, "
             f"T {p['T']:g}; rho " + ", ".join(f"{r:g} ({S[r]['n']} pairs)" for r in rhos) + f"; base seed {base}"]
    nflat = sum(len(S[r]["flat"]) for r in rhos)
    nbins = sum(len(S[r]["bins"]) for r in rhos)
    inst, passc, failc = {}, {}, {}
    worst_flat, worst_p1 = [], []
    for rho in rhos:
        n = S[rho]["n"]
        tcf = float(stats.t.ppf(1 - p["alpha"] / (2 * nflat), n - 1))
        tcb = float(stats.t.ppf(1 - p["alpha"] / (2 * nbins), n - 1))
        for t, d, se in S[rho]["flat"]:
            z = abs(d) / math.hypot(se, p["flat_floor"])
            inst[f"flat rho {rho:g} t {t:.2f}"] = bool(z < tcf)
            worst_flat.append((z / tcf, f"flat rho {rho:g} t {t:.2f}: d {d:+.4f}, SE {se:.4f}, z {z:.2f} (t_crit {tcf:.2f})"))
        ok = True
        for t, x, D, se, Pe, Pl in S[rho]["bins"]:
            z = abs(D - Pe) / math.hypot(se, p["floor"])
            if z >= tcb:
                ok = False
            worst_p1.append((z / tcb, f"rho {rho:g} t {t:.2f} x {x:+.2f}: D {D:+.3f} +- {se:.3f}, exact {Pe:+.3f}, "
                                      f"local {Pl:+.3f}; z {z:.2f} (t_crit {tcb:.2f})"))
        passc[f"P1 exact mean rho {rho:g}"] = ok
        failc[f"P1 exact mean rho {rho:g}"] = not ok
        dl, sl = S[rho]["centre_local"]
        de, sd = S[rho]["centre_exact"]
        passc[f"P2 not local rho {rho:g}"] = bool(dl - nse * sl > 0)
        failc[f"P2 not local rho {rho:g}"] = bool(dl + nse * sl < 0)
        lines.append(f"rho {rho:g}: centre t >= 4.5: D - local {dl:+.4f} +- {sl:.4f} (predicted {S[rho]['centre_pred']:+.4f}); "
                     f"D - exact {de:+.4f} +- {sd:.4f}")
        b = S[rho]["bins"]
        c2e = np.mean([((D - Pe) / math.hypot(se, p["floor"])) ** 2 for t, x, D, se, Pe, Pl in b])
        c2l = np.mean([((D - Pl) / math.hypot(se, p["floor"])) ** 2 for t, x, D, se, Pe, Pl in b])
        lines.append(f"   {len(b)} bins: chi2/bin against the exact mean {c2e:.2f}, against the local model {c2l:.2f}")
    lines.append("instrument (flat parents), the four largest:")
    lines += ["  " + s for _, s in sorted(worst_flat, reverse=True)[:4]]
    lines.append("P1, the four largest bins against the exact mean:")
    lines += ["  " + s for _, s in sorted(worst_p1, reverse=True)[:4]]
    lines.append("instrument failures: " + (", ".join(k for k, v in inst.items() if not v) or "none"))
    lines.append("passed: " + "; ".join(f"{k}: {'yes' if v else 'no'}" for k, v in passc.items()))
    lines.append("failed: " + (", ".join(k for k, v in failc.items() if v) or "none"))
    if not all(inst.values()):
        st = "INSTRUMENT_FAILED"
    elif any(failc.values()):
        st = "FAILED"
    elif all(passc.values()):
        st = "PASSED"
    else:
        st = "UNDECIDED"
    lines.append(f"cost: {time.time() - t0:.0f}s")
    return Result(st, lines, {"instrument": inst, "passed": passc, "failed": failc,
                              "stats": {str(k): v for k, v in S.items()}})
