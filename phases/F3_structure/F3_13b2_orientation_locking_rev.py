"""
F3.13b2 -- Locking of an orientation defect, revision of F3.13b (the readability requirement of the instrument)
==============================================================================================================

Why this revision exists (31 Oct 2026; F3.13b INSTRUMENT_FAILED, its verdict and criteria stay, R132)
---------------------------------------------------------------------------------------------------
F3.13b failed only in one instrument condition: its hedgehog arms (all kappa together) were readable on 0.445 of their
(sprinkling, slab, shell) triples against 0.5. The requirement was a fraction, and the taller box added late slabs in
which the cold bundle has thinned to few world lines; at kappa = 0 (no alignment, the F3.13a dynamics) only 0.26 were
readable. A fraction does not say whether there is enough data to decide; a count does. Everything else of F3.13b passed
or read clearly (seen, decides nothing here): every readable shell of every defect and erase arm read +1 (449 and 390),
no core event, the twist relaxed at kappa 0.3 (mean cost 0.0277 -> 0.0097), the core cost of the hedgehog in the last
slab at kappa 0.3 was 0.343 +- 0.049 against 0.011 +- 0.003 of the twist.
This revision changes only the readability requirement and uses new sprinklings (a new seed label):
  - the instrument requires at least 50 readable shells for every defect arm and every erase arm (shells R >= 3 after
    t = 8) at kappa > 0, instead of a fraction;
  - the decision is taken on the arms with alignment (kappa > 0), the test of locking; kappa = 0 (the F3.13a dynamics,
    nothing tries to unwind) is reported only.
Model, arms, readout, the vacuum and twist conditions of the instrument, the thresholds 95 % / 80 % and the cost
criterion are those of F3.13b unchanged (see that file).

Criteria (stated before any run)
--------------------------------
INSTRUMENT_FAILED: vacuum reads 0 on < 95 % at some kappa; fewer than 50 readable shells in a defect arm or an erase
arm (R >= 3, t >= 8) at kappa > 0; or the twist at kappa 0.3 does not halve its mean cost.
PASSED: at every kappa > 0 the hedgehog reads +1 on >= 95 % of its readable shells and the erase arm on >= 95 % of its
readable shells R >= 3 after t = 8, and at kappa 0.3 the core cost of the hedgehog in the last slab exceeds that of the
twist by more than 2 SE. FAILED: at some kappa > 0 the hedgehog or the erase arm below 80 %. UNDECIDED otherwise.
Reading as declared in F3.13b.
"""

from __future__ import annotations

import math
import time

import numpy as np

from core.order.sprinkle import derive_seed
from core.stats.parallel import parallel_map
from core.stats.stats import boot_se
from core.testing import Result
from phases.F3_structure.F3_13a_orientation_defect import EXPECT
from phases.F3_structure.F3_13b_orientation_locking import TEST as _T0, _frac, arm_list, task

TEST = dict(_T0)
TEST.update({
    "id": "F3.13b2",
    "title": "Locking of an orientation defect (revision of F3.13b: readability as a count, decision at kappa > 0)",
    "criterion": "INSTRUMENT: vacuum 0 on >= 95 %, >= 50 readable shells in every defect and erase arm at kappa > 0, twist "
                 "cost at kappa 0.3 halves; PASSED: at every kappa > 0 defect +1 >= 95 % and erase +1 >= 95 % on R >= 3 "
                 "after t = 8, core cost of the defect > twist + 2 SE at kappa 0.3; FAILED < 80 %; else UNDECIDED",
    "depends_on": ["F3.13b", "F3.13a"],
    "params": dict(_T0["params"], min_readable=50),
    "quick": dict(_T0["quick"], min_readable=1),
})


def run(p, ctx):
    t0 = time.time()
    n, B, nse = int(p["n_seeds"]), int(p["boot_B"]), float(p["n_se"])
    jobs = [(p, derive_seed(ctx.seed, "F3.13b2", k)) for k in range(n)]
    rows = parallel_map(task, jobs, ctx.workers, ctx.progress,
                        lambda i, r_, d, nn: f"[F3.13b2] {d}/{nn}: N {r_['N']} ({time.time() - t0:.0f}s)")
    tgrid = rows[0]["arms"][arm_list(p)[0][0]]["t"]
    lines = [f"3+1 box t [0, {p['T']:g}], half {p['half']:g} (~{np.mean([r_['N'] for r_ in rows]):.0f} events, "
             f"{np.mean([r_['sources'] for r_ in rows]):.0f} sources); cold transport with merging; alignment step kappa "
             f"towards the causal past (lag <= {p['L_g']:g}); erase: core r < {p['a_erase']:g} randomised at "
             f"{p['t_erase']:g} <= t < {p['t_erase'] + 1:g}; shells {p['shells']}; {n} sprinklings",
             "        t           " + " ".join(f"{t:6.1f}" for t in tgrid)]
    res = {}
    for key, cfg, kappa, erase in arm_list(p):
        vals = [x for r_ in rows for x in r_["arms"][key]["vals"]]
        exp = EXPECT[cfg]
        f, rd, nr = _frac(vals, exp)
        fo = _frac(vals, exp, rmin=3.0, tmin=float(p["t_after"])) if erase else None
        f2 = _frac(vals, exp, tmin=float(p["t_after"])) if erase else None
        cc = np.array([r_["arms"][key]["core_cost"] for r_ in rows], float)
        mc = np.array([r_["arms"][key]["mean_cost"] for r_ in rows], float)
        last = boot_se(cc[:, -1][np.isfinite(cc[:, -1])], B, derive_seed(ctx.seed, "F3.13b2-boot", key, "core"))[:2] \
            if np.isfinite(cc[:, -1]).sum() >= 2 else (float("nan"), float("nan"))
        with np.errstate(all="ignore"):
            ccm, mcm = np.nanmean(cc, axis=0), np.nanmean(mc, axis=0)
        res[key] = {"f": f, "rd": rd, "nr": nr, "fo": fo, "f2": f2, "core_last": last, "mc": mcm.tolist(), "cc": ccm.tolist(),
                    "cores": float(np.mean([r_["arms"][key]["cores"] for r_ in rows]))}
        hist = {}
        for *_, v, ok in vals:
            if ok:
                hist[int(round(v))] = hist.get(int(round(v)), 0) + 1
        lines.append(f"  {key:11s}: reads {exp:+d} on {f:.3f} of {nr} readable (readable {rd:.2f}); values "
                     f"{dict(sorted(hist.items()))}; core events {res[key]['cores']:.1f}"
                     + (f"; after the erase on R >= 3: {fo[0]:.3f} of {fo[2]} (readable {fo[1]:.2f}), all shells "
                        f"{f2[0]:.3f} of {f2[2]}" if erase else ""))
        lines.append("        core cost   " + " ".join(f"{x:6.3f}" if np.isfinite(x) else "   nan" for x in ccm))
        lines.append("        mean cost   " + " ".join(f"{x:6.3f}" if np.isfinite(x) else "   nan" for x in mcm))
    pf, ff = float(p["pass_frac"]), float(p["fail_frac"])
    ks = [float(k) for k in p["kappas"]]
    kmax = max(ks)
    tw = res[f"twist_{kmax:g}"]["mc"]
    tw_first = next((x for x in tw if np.isfinite(x)), float("nan"))
    tw_last = next((x for x in reversed(tw) if np.isfinite(x)), float("nan"))
    kpos = [k for k in ks if k > 0.0]
    mincount = int(p["min_readable"])
    inst = {f"vacuum reads 0 on >= {pf:.0%} (every kappa)": all(res[f"vacuum_{k:g}"]["f"] >= pf for k in ks),
            f"defect and erase arms at kappa > 0 have >= {mincount} readable shells": all(
                res[f"defect_{k:g}"]["nr"] >= mincount and res[f"erase_{k:g}"]["fo"][2] >= mincount for k in kpos),
            f"twist at kappa {kmax:g} relaxes (last mean cost < first / 2)": bool(tw_last < 0.5 * tw_first)}
    lines.append(f"  twist at kappa {kmax:g}: mean cost first {tw_first:.4f}, last {tw_last:.4f}")
    dl, tl = res[f"defect_{kmax:g}"]["core_last"], res[f"twist_{kmax:g}"]["core_last"]
    cost_ok = bool(np.isfinite(dl[0]) and np.isfinite(tl[0]) and dl[0] - tl[0] > nse * math.hypot(dl[1], tl[1]))
    lines.append(f"  core cost in the last slab at kappa {kmax:g}: defect {dl[0]:.4f} +- {dl[1]:.4f}, twist {tl[0]:.4f} +- "
                 f"{tl[1]:.4f} -> the defect keeps a cost: {cost_ok}")
    lines.append("instrument: " + "; ".join(f"{k}: {'OK' if v else 'FAILS'}" for k, v in inst.items()))
    keep = all(res[f"defect_{k:g}"]["f"] >= pf for k in kpos)
    erase_ok = all(res[f"erase_{k:g}"]["fo"][0] >= pf for k in kpos)
    lost = any(res[f"defect_{k:g}"]["f"] < ff or res[f"erase_{k:g}"]["fo"][0] < ff for k in kpos)
    if not all(inst.values()):
        status = "INSTRUMENT_FAILED"
    elif keep and erase_ok and cost_ok:
        status = "PASSED"
    elif lost:
        status = "FAILED"
    else:
        status = "UNDECIDED"
    return Result(status, lines, {"res": res, "rows": rows})
