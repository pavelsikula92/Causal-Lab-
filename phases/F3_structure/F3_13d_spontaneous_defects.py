"""
F3.13d -- Spontaneous defects: does the alignment of disordered orientations create locked defects by itself?
=============================================================================================================

Why this test exists (31 Oct 2026, the author's objection after F3.13c: "the hedgehog is a predetermined concept")
-----------------------------------------------------------------------------------------------------------------
F3.13a-b2 showed a conditional fact: if the flow carries an orientation that aligns with its causal neighbourhood, then a
hedgehog put in by hand cannot be removed. The shape of the defect was an input. Here no shape is given: every source
event receives an independent random orientation, and the question is whether the alignment of this disorder leaves
behind defects with a non-zero charge that persist (the Kibble-Zurek mechanism: domains that align independently meet in
knots that cannot be untied locally). The disorder is an input; the defects, their places and their charges are not.
The orientation and its alignment remain the two assumptions added in F3.13a / F3.13b.

Model (approved by the author on 31 Oct before any run; kappa and T changed on 1 Nov, see below; everything not named
is F3.13b2)
-------------------------------------------------------------------------------------------------------------------
3+1 Poisson sprinkling, density 1, box t [0, 30], |x|, |y|, |z| <= 7 (~82 000 events). Sources: events with t < 2 and
r < 5, unit flow at rest. Orientation of each source event: an independent isotropic unit vector of its rest space (a
seed of its own per sprinkling and arm). Transport, merging and alignment of F3.13b (the coordinate window of F3.13b is
exact here: every flow is at rest, since every source is at rest and P is carried and summed). Arms:
  random_0.9 (decides), random_0.3 (reported), random_0 (no alignment: transport and merging only; reported),
  vacuum_0.9 (a uniform orientation: nothing may arise from nothing; instrument).
Readout. Slabs of width 2 with centres t = 5, 7, ..., 27. Candidate centres: the points of the integer grid with |c| <= 3
(123 points). For each slab and centre the degree of the shell R = 2 around c (carrying events of the slab with
||x - c| - 2| < 0.75), read by the F3.13a rule (an unresolved field is unreadable, so noise cannot give a wrong integer).
A defect site is a cluster of centres of one slab that read the same non-zero integer and are connected by steps of
distance <= 1.5; its place is the mean of its centres, its charge that integer.
  S        per sprinkling: the mean number of sites per slab over the late slabs (centre t >= 11); S = mean over
           sprinklings, SE by bootstrap over sprinklings;
  p        persistence: of all sites in slabs with centre 7 <= t <= 25 (pooled over sprinklings), the fraction that have a
           site of the same charge within distance 1.5 in the next slab (t + 2);
reported: S per slab, the readable fraction, the charges of the sites, the net charge on the shell R = 4 around the
origin, the mean cost per carrying event per slab.

Criteria (stated before any run, approved by the author)
--------------------------------------------------------
INSTRUMENT_FAILED: vacuum_0.9 reads 0 on < 95 % of its readable shells; or random_0.9 has fewer than 200 readable
shells (all slabs and centres, pooled).
PASSED (at kappa = 0.9): S - 2 SE > 0.25 and p >= 0.8: defects arise by themselves and are locked.
FAILED: S + 2 SE < 0.1 (practically no defect arises), or p < 0.5 (what arises flickers, it is not locked).
UNDECIDED otherwise.

Reading declared now. PASSED: locked defects are a natural product of the dynamics of a carried orientation; the
hedgehog is not a predetermined shape but what a disordered flow leaves behind; the open question moves to the
orientation itself (can it be read from the order instead of being added?). FAILED: the branch of a carried orientation
does not create particles by itself; a serious negative result for the line F3.13. Limits: the orientation and its
alignment are assumptions; cold flow at rest, no field; one size of the region; the grid sees defects at the scale R = 2
only.
A development sprinkling in the container before the full run checks the readability only; if random_0.3 had too few
readable shells there, the author is asked before anything is changed.

Development and the change of 1 Nov (before any full run; decides nothing)
-------------------------------------------------------------------------
The approved set-up (kappa 0.3 deciding, T = 20) could not be read: one full-size sprinkling in the container gave
random_0.3 no readable shell of 861 (mean cost 0.49 -> 0.17 by t = 17; vacuum 859 of 861). Two diagnostic quick runs
(4 sprinklings, T = 30, on the author's computer): kappa 0.3 / 0.1 / 0 gave 255 / 0 / 0 readable shells of 5904
(readability rising to 0.12 from t = 23; mean cost at kappa 0.3 falling about as 1/t, 0.485 -> 0.099); kappa 0.6 / 0.9
gave 588 / 1308 (readability up to 0.20 / 0.39; mean cost levelling at ~0.05 for both, i.e. limited by the causal
reach of the alignment, not by kappa). In every arm every readable shell read 0: no site at any kappa. Approved by the
author on 1 Nov: the deciding arm is kappa = 0.9 and T = 30 (late slabs t = 11 to 27), for readability only; since no
kappa showed a site, the choice does not favour an outcome. Arms kappa 0.3 and 0 are reported. All thresholds are
unchanged; new sprinklings (seed label F3.13d-main, not the label of the diagnostic runs). Expected from the development: FAILED (S = 0).
"""

from __future__ import annotations

import math
import time

import numpy as np

from core.order.links import links_csr
from core.order.sprinkle import Box, derive_seed, sprinkle_poisson
from core.stats.parallel import parallel_map
from core.stats.stats import boot_se
from core.testing import Result
from phases.F3_structure.F3_08_retarded_field import children_csr
from phases.F3_structure.F3_13a_orientation_defect import shell_value
from phases.F3_structure.F3_13b_orientation_locking import transport_align

TEST = {
    "id": "F3.13d",
    "phase": "F3",
    "title": "Spontaneous defects: does the alignment of disordered orientations create locked defects by itself?",
    "theory_part": "F3 emergence of structure; a particle as a topologically locked defect that the dynamics creates "
                   "(the author's objection of 31 Oct: the hedgehog must not be put in by hand)",
    "question": "With independent random orientations at the sources and the alignment of F3.13b (kappa = 0.9), do defects "
                "with a non-zero degree arise by themselves and persist?",
    "method": "3+1 as F3.13b2 with a taller box (t [0, 30], half 7); random isotropic orientation per source event; arms "
              "random at kappa 0.9 / 0.3 / 0, vacuum at 0.9; degree on shells R = 2 around the points of the integer grid |c| <= 3 in "
              "slabs of width 2; sites = clusters of centres with the same non-zero degree; S = sites per late slab, "
              "p = persistence of a site to the next slab",
    "prediction": "the Kibble-Zurek mechanism: defects arise from the disorder and persist (some annihilate in pairs)",
    "criterion": "INSTRUMENT: vacuum 0 >= 95 %, random_0.9 >= 200 readable shells; PASSED (kappa 0.9): S - 2 SE > 0.25 "
                 "and p >= 0.8; FAILED: S + 2 SE < 0.1 or p < 0.5; otherwise UNDECIDED",
    "null_model": "vacuum (uniform orientation); random without alignment (kappa = 0, reported)",
    "level": "mechanism",
    "meters": [],
    "depends_on": ["F3.13b2", "F3.13a"],
    "params": {"T": 30.0, "half": 7.0, "r_src": 5.0, "t_src": 2.0, "slab": 2.0, "t_read": 4.0, "R": 2.0, "R_net": 4.0,
               "grid_r": 3.0, "grid_step": 1.0, "link": 1.5, "t_late": 11.0, "t_pers": 7.0, "kappas": [0.9, 0.3, 0.0], "kappa_dec": 0.9,
               "L_g": 2.0, "t_erase": 1e9, "a_erase": 0.0, "core_eps": 1e-6, "n_seeds": 16, "n_se": 2.0,
               "boot_B": 4000, "pass_frac": 0.95, "min_readable": 200, "S_pass": 0.25, "S_fail": 0.1,
               "p_pass": 0.8, "p_fail": 0.5},
    "parts": {"main": {}},
    "quick": {"n_seeds": 1, "T": 10.0, "half": 6.0, "r_src": 4.0, "grid_r": 2.0, "t_late": 7.0, "t_pers": 5.0,
              "boot_B": 200, "min_readable": 1},
    "estimate": "~7.5 min per order on one core (~82 000 events in 3+1, four arms, ~1 500 shells per arm); 16 orders; "
                "memory ~1.5 GB per process: --workers 8 on 16 GB",
}


# ============================================================================================== orientations and grid

def random_orient(n, seed):
    """n independent isotropic unit 3-vectors (rows)."""
    g = np.random.default_rng(seed).normal(size=(n, 3))
    return g / np.maximum(np.linalg.norm(g, axis=1), 1e-300)[:, None]


def grid(p):
    """Candidate centres: points of the grid with step grid_step inside the ball |c| <= grid_r."""
    r, h = float(p["grid_r"]), float(p["grid_step"])
    k = np.arange(-math.floor(r / h), math.floor(r / h) + 1) * h
    G = np.array(np.meshgrid(k, k, k, indexing="ij")).reshape(3, -1).T
    return G[np.linalg.norm(G, axis=1) <= r + 1e-9]


def sites(G, vals, link):
    """Clusters of centres with the same non-zero integer value, connected by steps <= link: [(place, charge, size)]."""
    out = []
    nz = [i for i, v in enumerate(vals) if v is not None and v != 0]
    seen = set()
    for i in nz:
        if i in seen:
            continue
        comp, stack = [], [i]
        seen.add(i)
        while stack:                                  # flood fill over centres with the same value
            j = stack.pop()
            comp.append(j)
            for k in nz:
                if k not in seen and vals[k] == vals[i] and np.linalg.norm(G[j] - G[k]) <= link + 1e-9:
                    seen.add(k)
                    stack.append(k)
        out.append((G[comp].mean(axis=0).tolist(), int(vals[i]), len(comp)))
    return out


# ============================================================================================== one sprinkling

def task(p, seed):
    rng = np.random.default_rng(seed)
    T, H = float(p["T"]), float(p["half"])
    X = sprinkle_poisson(Box(4, [0.0, -H, -H, -H], [T, H, H, H]), 1.0, rng)
    N = len(X)
    ptr, idx = links_csr(X, 1.0)
    cptr, cidx = children_csr(ptr, idx, N)
    rs, ts = float(p["r_src"]), float(p["t_src"])
    src = np.flatnonzero((X[:, 0] < ts) & (np.linalg.norm(X[:, 1:], axis=1) < rs))
    s = np.zeros(N)
    s[src] = 1.0
    U0 = np.zeros((N, 4))
    U0[:, 0] = 1.0
    G = grid(p)
    sl, R, Rn = float(p["slab"]), float(p["R"]), float(p["R_net"])
    edges = np.arange(float(p["t_read"]), T - sl + 1e-9, sl)
    arms = [(f"random_{float(k):g}", "random", float(k)) for k in p["kappas"]] + [(f"vacuum_{float(p['kappa_dec']):g}", "vacuum", float(p["kappa_dec"]))]
    out = {"N": N, "sources": int(src.size), "centres": int(len(G)), "arms": {}}
    for key, cfg, kappa in arms:
        N0 = np.zeros((N, 4))
        if cfg == "random":
            N0[src, 1:] = random_orient(src.size, derive_seed(seed, "F3.13d-orient", key))
        else:
            N0[src, 1] = 1.0
        Ph, Nh, core, cost = transport_align(X, cptr, cidx, s, U0, N0, kappa, p, None)
        held = Ph[:, 0] > 0.0
        rec = {"t": [], "n_read": [], "n_zero": [], "sites": [], "net": [], "mean_cost": [], "cores": int(core.sum())}
        for a, b in zip(edges[:-1], edges[1:]):
            okn = held & np.isfinite(Nh[:, 0]) & (X[:, 0] >= a) & (X[:, 0] < b)
            P = X[okn, 1:]
            nsp = Nh[okn, 1:] / np.maximum(np.linalg.norm(Nh[okn, 1:], axis=1), 1e-300)[:, None]
            vals = []
            for c in G:                               # degree of the shell R around every candidate centre
                rel = P - c
                sel = np.abs(np.linalg.norm(rel, axis=1) - R) < 0.75
                v, ok = shell_value(rel[sel], nsp[sel], 3)
                vals.append(int(round(v)) if ok and abs(v - round(v)) < 1e-6 else None)
            rec["t"].append(0.5 * (a + b))
            rec["n_read"].append(sum(v is not None for v in vals))
            rec["n_zero"].append(sum(v == 0 for v in vals))
            rec["sites"].append(sites(G, vals, float(p["link"])))
            sel = np.abs(np.linalg.norm(P, axis=1) - Rn) < 0.75
            v, ok = shell_value(P[sel], nsp[sel], 3)
            rec["net"].append(int(round(v)) if ok and abs(v - round(v)) < 1e-6 else None)
            cc = cost[okn]
            f = np.isfinite(cc)
            rec["mean_cost"].append(float(cc[f].mean()) if f.any() else float("nan"))
        out["arms"][key] = rec
    return out


# ============================================================================================== run

def persistence(rec, t_min, link):
    """(persisting sites, sites) over slabs with centre t_min <= t <= the last but one."""
    hit = tot = 0
    for i in range(len(rec["t"]) - 1):
        if rec["t"][i] < t_min - 1e-9:
            continue
        nxt = rec["sites"][i + 1]
        for place, q, _ in rec["sites"][i]:
            tot += 1
            if any(q2 == q and np.linalg.norm(np.subtract(place, pl2)) <= link + 1e-9 for pl2, q2, _ in nxt):
                hit += 1
    return hit, tot


def run(p, ctx):
    t0 = time.time()
    n, B, nse = int(p["n_seeds"]), int(p["boot_B"]), float(p["n_se"])
    # seed label "F3.13d-main" (1 Nov): the diagnostic quick runs used "F3.13d" with the same box, so the full run must
    # not reuse their sprinklings
    jobs = [(p, derive_seed(ctx.seed, "F3.13d-main", k)) for k in range(n)]
    rows = parallel_map(task, jobs, ctx.workers, ctx.progress,
                        lambda i, r_, d, nn: f"[F3.13d] {d}/{nn}: N {r_['N']} ({time.time() - t0:.0f}s)")
    vac = f"vacuum_{float(p['kappa_dec']):g}"
    tg = rows[0]["arms"][vac]["t"]
    lines = [f"3+1 box t [0, {p['T']:g}], half {p['half']:g} (~{np.mean([r_['N'] for r_ in rows]):.0f} events, "
             f"{np.mean([r_['sources'] for r_ in rows]):.0f} sources); random isotropic orientation per source; F3.13b "
             f"transport and alignment; degree on shells R = {p['R']:g} around {rows[0]['centres']} grid centres "
             f"(|c| <= {p['grid_r']:g}); {n} sprinklings",
             "        t                " + " ".join(f"{t:6.1f}" for t in tg)]
    late = np.asarray(tg) >= float(p["t_late"]) - 1e-9
    res = {}
    for key in rows[0]["arms"]:
        recs = [r_["arms"][key] for r_ in rows]
        nread = int(sum(sum(r["n_read"]) for r in recs))
        nzero = int(sum(sum(r["n_zero"]) for r in recs))
        ntot = int(len(recs) * len(tg) * rows[0]["centres"])
        cnt = np.array([[len(s_) for s_ in r["sites"]] for r in recs], float)
        S_k = cnt[:, late].mean(axis=1)
        S = boot_se(S_k, B, derive_seed(ctx.seed, "F3.13d-boot", key))[:2]
        hit, tot = 0, 0
        for r in recs:
            h_, t_ = persistence(r, float(p["t_pers"]), float(p["link"]))
            hit, tot = hit + h_, tot + t_
        pers = hit / tot if tot else float("nan")
        charges = {}
        for r in recs:
            for s_ in r["sites"]:
                for _, q, _ in s_:
                    charges[q] = charges.get(q, 0) + 1
        nets = {}
        for r in recs:
            for t_, v in zip(r["t"], r["net"]):
                if v is not None and t_ >= float(p["t_late"]) - 1e-9:
                    nets[v] = nets.get(v, 0) + 1
        mc = np.nanmean(np.array([r["mean_cost"] for r in recs], float), axis=0)
        res[key] = {"nread": nread, "zero_frac": (nzero / nread if nread else float("nan")), "S": S, "p": pers,
                    "p_n": tot, "charges": charges, "net": nets}
        lines.append(f"  {key:11s}: readable {nread} of {ntot} shells ({nread / ntot:.2f}); zero on {res[key]['zero_frac']:.3f}"
                     f" of readable; S (late) = {S[0]:.3f} +- {S[1]:.3f}; persistence {hit}/{tot} = {pers:.3f}; site "
                     f"charges {dict(sorted(charges.items()))}; net charge R = {p['R_net']:g} (late) {dict(sorted(nets.items()))}"
                     f"; core events {np.mean([r['cores'] for r in recs]):.1f}")
        lines.append("        sites per slab   " + " ".join(f"{x:6.2f}" for x in cnt.mean(axis=0)))
        lines.append("        readable         " + " ".join(f"{x:6.2f}" for x in
                                                         np.mean([r["n_read"] for r in recs], axis=0) / rows[0]["centres"]))
        lines.append("        mean cost        " + " ".join(f"{x:6.3f}" if np.isfinite(x) else "   nan" for x in mc))
    pf, mr = float(p["pass_frac"]), int(p["min_readable"])
    dec = f"random_{float(p['kappa_dec']):g}"       # the deciding arm (kappa 0.9)
    if dec not in res:                                # only in a quick diagnostic with other kappas: decides nothing
        dec = [k for k in res if k.startswith("random_")][0]
        lines.append(f"note: kappa_dec = {p['kappa_dec']:g} was not run; the status below is computed on {dec}")
    d = res[dec]
    inst = {f"vacuum reads 0 on >= {pf:.0%}": res[vac]["zero_frac"] >= pf,
            f"{dec} has >= {mr} readable shells": d["nread"] >= mr}
    lines.append("instrument: " + "; ".join(f"{k}: {'OK' if v else 'FAILS'}" for k, v in inst.items()))
    S, se = d["S"]
    pv = d["p"]
    lines.append(f"decision ({dec}): S - {nse:g} SE = {S - nse * se:.3f} (pass > {p['S_pass']:g}), S + {nse:g} SE = "
                 f"{S + nse * se:.3f} (fail < {p['S_fail']:g}); p = {pv:.3f} (pass >= {p['p_pass']:g}, fail < {p['p_fail']:g})")
    if not all(inst.values()):
        status = "INSTRUMENT_FAILED"
    elif S - nse * se > float(p["S_pass"]) and np.isfinite(pv) and pv >= float(p["p_pass"]):
        status = "PASSED"
    elif S + nse * se < float(p["S_fail"]) or (np.isfinite(pv) and pv < float(p["p_fail"])):
        status = "FAILED"
    else:
        status = "UNDECIDED"
    return Result(status, lines, {"res": res, "rows": rows})
