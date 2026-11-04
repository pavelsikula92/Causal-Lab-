"""
F3.12a2 -- Resistance of the substrate to rerouting, revision of F3.12a (instrument tolerance fitted to the solver)
=================================================================================================================

Why this revision exists (31 Oct 2026; F3.12a INSTRUMENT_FAILED, its verdict and criteria stay, R132)
---------------------------------------------------------------------------------------------------
F3.12a failed only in its instrument: the unperturbed program returned the original flow to 1.3e-7 of the total link
flux (limit 1e-9) and the balance held to 4.6e-9 (limit 1e-9). Those limits were a design error: the linear program is
solved by HiGHS, whose primal feasibility tolerance is 1e-7; a solution exact to 1e-9 cannot be demanded of it. The
effects the test reads are 1e-3 to 1e-1 of the total flux. The reported numbers of F3.12a (seen, they decide nothing
here): I_cascade at sigma0 = 0.2: T 0.0129 +- 0.0004, K 0.0539 +- 0.0011, P 0.0115 +- 0.0004.
This revision changes only the two tolerances of the instrument, to 1e-6 of the total link flux (ten times the solver
tolerance, four orders below the smallest effect read), and uses new sprinklings (a new seed label). Model,
configurations T, K, P, strengths, the classification and its thresholds are those of F3.12a unchanged (see that file).

Criteria (stated before any run)
--------------------------------
Instrument (INSTRUMENT_FAILED if any fails): unperturbed I_flow <= 1e-6; balance of every solution <= 1e-6 of the total
link flux; every program solved.
Classification at sigma0 = 0.2 as F3.12a: topology dominates I_c(T) - 2 SE > 2 (I_c(K) + 2 SE); capacity dominates the
reverse; no resistance both I_c + 2 SE < 0.01; composite otherwise. Status RECORDED with the classification.
Reading as declared in F3.12a.
"""

from __future__ import annotations

import time

import numpy as np

from core.order.sprinkle import derive_seed
from core.stats.parallel import parallel_map
from core.stats.stats import boot_se
from core.testing import Result
from phases.F3_structure.F3_12a_rerouting_resistance import TEST as _T0, _classify, seed_task

TEST = dict(_T0)
TEST.update({
    "id": "F3.12a2",
    "title": "Resistance of the substrate to rerouting (revision of F3.12a: instrument tolerance fitted to the solver)",
    "criterion": "as F3.12a; INSTRUMENT: unperturbed I_flow <= 1e-6, balance <= 1e-6 of the total link flux, all programs "
                 "solved. Status RECORDED with the classification.",
    "depends_on": ["F3.12a", "F7.1a", "F7.1b", "F1.1"],
    "params": dict(_T0["params"], zero_tol=1e-6, bal_tol=1e-6),
    "quick": dict(_T0["quick"]),
})


def run(p, ctx):
    t0 = time.time()
    n, B, nse = int(p["n_seeds"]), int(p["boot_B"]), float(p["n_se"])
    jobs = [(p, derive_seed(ctx.seed, "F3.12a2", ctx.part, k)) for k in range(n)]
    rows = parallel_map(seed_task, jobs, ctx.workers, ctx.progress,
                        lambda i, r_, d, nn: f"[F3.12a2] {d}/{nn}: N {r_['N']} links {r_['links']} "
                                             f"({time.time() - t0:.0f}s)")
    zero = max(r_["zero"]["I_flow"] for r_ in rows)
    bal = max([r_["zero"]["bal"] for r_ in rows] +
              [r_["res"][k][c]["bal"] for r_ in rows for k in r_["res"] for c in ("T", "K", "P")])
    solved = all(r_["ok"] for r_ in rows)
    lines = [f"2+1 box t [0, {p['T']:g}], |x| <= {p['half']:g} (~{np.mean([r_['N'] for r_ in rows]):.0f} events, "
             f"~{np.mean([r_['links'] for r_ in rows]):.0f} links); source t < {p['src_t']:g}, r < {p['src_r']:g}; "
             f"lens {p['window'][0]:g} <= t <= {p['window'][1]:g}, w = {p['width']:g}; {n} sprinklings",
             f"  unperturbed: max I_flow {zero:.1e}; max balance residual {bal:.1e}; all programs solved: {solved}"]
    summ = {}
    for s0 in p["sigmas"]:
        k = str(s0)
        summ[k] = {}
        for cfg in ("T", "K", "P"):
            ok = [r_["res"][k][cfg] for r_ in rows if r_["res"][k][cfg]["status"] == 0]
            bs = lambda *tag: derive_seed(ctx.seed, "F3.12a2-boot", ctx.part, s0, cfg, *tag)
            summ[k][cfg] = {q_: boot_se([o[q_] for o in ok], B, bs(q_))[:2] for q_ in ("I_flow", "I_forced", "I_cascade")}
            for q_ in ("removed", "created", "added"):
                if q_ in ok[0]:
                    summ[k][cfg][q_] = float(np.mean([o[q_] for o in ok]))
        lines.append(f"  sigma0 = {s0:g}:")
        for cfg, lab in (("T", "T topology "), ("K", "K capacity "), ("P", "P throughput")):
            s_ = summ[k][cfg]
            extra = "".join(f", {q_} {s_[q_]:.0f}" for q_ in ("removed", "created", "added") if q_ in s_)
            lines.append(f"    {lab}: I_flow {s_['I_flow'][0]:.4f} +- {s_['I_flow'][1]:.4f}; I_forced {s_['I_forced'][0]:.4f}; "
                         f"I_cascade {s_['I_cascade'][0]:+.4f} +- {s_['I_cascade'][1]:.4f}{extra}")
        lines.append(f"    reported: P against K: {_classify(summ[k]['P']['I_cascade'], summ[k]['K']['I_cascade'], nse, p['none_tol'])}")
    kd = str(p["decide"])
    cls = _classify(summ[kd]["T"]["I_cascade"], summ[kd]["K"]["I_cascade"], nse, float(p["none_tol"])) if kd in summ else "n/a"
    inst = {"unperturbed I_flow <= tol": zero <= p["zero_tol"], "balance <= tol": bal <= p["bal_tol"],
            "all programs solved": solved}
    lines.append("instrument: " + "; ".join(f"{k_}: {'OK' if v_ else 'FAILS'}" for k_, v_ in inst.items()))
    lines.append(f"classification at sigma0 = {p['decide']:g} (T against K): {cls}")
    status = "INSTRUMENT_FAILED" if not all(inst.values()) else "RECORDED"
    return Result(status, lines, {"summary": summ, "classification": cls, "rows": rows})
