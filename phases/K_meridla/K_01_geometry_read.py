"""
K.1 -- Calibration of the geometry reading of a causal order (3+1), ported from E54
==================================================================================

Why: F6 (dynamics of order) and later F0 must decide whether an order
produced by a dynamics is a manifold-like geometry and of which dimension.
The instrument is the E38-B reading of the 2.x archive (gate G2 of the
E43 certificate): interval dimensions (Myrheim-Meyer and chain scaling),
homogeneity over the height, the smeared Benincasa-Dowker action, and
three matched null models (core.geometry.order_read, module notes). In
3+1 the reading sits at the edge of the instrument at small sizes
(E42a: D_MM 3.56 vs D_chain 4.34 at 8 000), so the size from which it
reads a 3+1 sprinkling as GEOMETRIC(d* = 4) must be calibrated. GTS Lab 3.x
did that in E54 (history: N* = 32 000, replicated); this is the port into
Causality Lab, with new seeds and no archive code.

The port (core.geometry.order_read) reproduces the archive statement by
statement with the same random streams. Checked against the archive on
19 Oct 2026 (2+1 and 3+1 sprinklings of 1 500 elements): identical link
structure, identical readings (interval dimensions, homogeneity, dS for
d = 2, 3, 4, links, relation density), identical matched nulls
(tp, seq, layered) and their readings, identical decisions. The two code
paths of the port (bitsets for general orders, coordinates for large
sprinklings) agree element by element (core check `order_read_paths`).

Subject: the 3+1 sprinkling of the unit diamond with the genesis element
(core.order.sprinkle.diamond_2x, the archive's point sets), read through
the coordinate path at sizes up to 128 000 (n_read sprinklings per size).
Matched nulls (tp, seq, layered) are built from the order of the
subject's sprinklings and read through the bitset path; they need memory
~N^2/4 bytes, so they are built at sizes <= null_cap = 16 000. Above it
the null reading of the largest size <= null_cap is used: a null that is
not manifold-like at 16 000 is assumed not to become 4-dimensional at
128 000 (assumption of E54, stated). G3 (bulk link growth) of E54 was a
reported number only and is not ported.

Instrument check: the 2+1 sprinkling at inst_sizes reads GEOMETRIC with
d* = 3 (2+1 is calibrated in 2.x), with its own matched nulls.
Null model: originary percolation p = 0.01 (the E43 negative control),
read through the bitset path at min(null_cap, largest size): it must not
be GEOMETRIC.
Criterion (E54): the last `consecutive` = 3 sizes all GEOMETRIC with
d* = 4 -> PASSED (calibrated) with N* = the first size of that trailing
run. Fewer -> FAILED (N* not reached in the plan). Percolation GEOMETRIC
-> UNDECIDED (not specific). Instrument check failing -> INSTRUMENT_FAILED.
"""

from __future__ import annotations

import time

import numpy as np

from core.order.sprinkle import diamond_2x, derive_seed
from core.geometry import order_read as OR
from core.stats.parallel import parallel_map
from core.testing import Result

TEST = {
    "id": "K.1",
    "phase": "K",
    "title": "Čtení geometrie z uspořádání (3+1)",
    "theory_part": "měřidlo pro F6 a F0: je uspořádání varieta a jaké dimenze; bez fyzikálního závěru",
    "question": "Od jaké velikosti čte čtení E38-B (intervalové dimenze, homogenita, BD akce, párované nulové modely) "
                "sprinkling 3+1 jako GEOMETRIC s d* = 4?",
    "method": "sprinklingy diamantu (body archivu) při 8 000 až 128 000 přes souřadnicovou cestu; párované nulové modely "
              "(tp, seq, layered) přes bitovou cestu do 16 000; kontrola 2+1; perkolace p = 0,01 jako negativní kontrola",
    "prediction": "GEOMETRIC(d* = 4) od N* = 32 000 nebo 64 000 (E54 v historii: od 32 000)",
    "criterion": "MĚŘIDLO: 2+1 při inst_sizes GEOMETRIC s d* = 3. PROŠEL: poslední 3 velikosti GEOMETRIC s d* = 4, N* = první "
                 "z nich. NEPROŠEL: méně. NEROZHODNUTO: perkolace GEOMETRIC (nespecifické).",
    "null_model": "párované nulové modely tp, seq, layered pro každý sprinkling (do 16 000); perkolace p = 0,01",
    "level": "calibration",
    "meters": [],
    "depends_on": [],
    "old_id": "E54",
    "params": {"D": 4, "sizes": [8000, 16000, 32000, 64000, 128000], "n_read": 2, "n_intervals": 60, "n_action": 60,
               "null_cap": 16000, "inst_sizes": [4000, 8000], "consecutive": 3, "perc_p": 0.01},
    "parts": {"all": {}},
    "quick": {"sizes": [600, 1000, 1500], "n_read": 1, "n_intervals": 15, "n_action": 5, "null_cap": 1000,
              "inst_sizes": [800]},
    "estimate": "mnoho hodin (E54): souřadnicová čtení až 128 000 a nulové modely 16 000 v čistém Pythonu (bitové množiny); "
                "paměť ~2–4 GB na proces u 128 000; úlohy běží paralelně",
}


def task(kind, D, n, seed, n_intervals, n_action, perc_p):
    """One unit of work (module level: worker processes)."""
    if kind == "read":                                   # subject / instrument: coordinate path
        S = diamond_2x(D, n, seed)
        return OR.read_coords(S.X, S.rho, seed, n_intervals, n_action)
    if kind == "nulls":                                  # matched nulls of the sprinkling's order: bitset path
        S = diamond_2x(D, n, seed)
        P = OR.poset_from_points(S.X)
        return {k: OR.read_poset(NP, seed, n_intervals, n_action) for k, NP in OR.matched_nulls(P, seed).items()}
    if kind == "perc":
        return OR.read_poset(OR.transitive_percolation_np(n, perc_p, seed=seed), seed, n_intervals, n_action)
    raise ValueError(kind)


def run(p, ctx):
    D = int(p["D"])
    sizes = [int(x) for x in p["sizes"]]
    inst_sizes = [int(x) for x in p["inst_sizes"]]
    nread, nint, nact, cap = int(p["n_read"]), int(p["n_intervals"]), int(p["n_action"]), int(p["null_cap"])
    pp = float(p["perc_p"])
    sd = lambda d, n, rep: derive_seed(ctx.seed, "K.1", d, n, rep) % (2 ** 31)
    null_size = lambda n: n if n <= cap else max([m for m in sizes if m <= cap] or [min(sizes)])
    jobs = []
    for n in inst_sizes:
        for rep in range(nread):
            jobs.append(("read", 3, n, sd(3, n, rep)))
        for rep in range(nread):
            jobs.append(("nulls", 3, min(n, cap), sd(3, min(n, cap), rep)))
    for n in sizes:
        for rep in range(nread):
            jobs.append(("read", D, n, sd(D, n, rep)))
    for nn in sorted({null_size(n) for n in sizes}):
        for rep in range(nread):
            jobs.append(("nulls", D, nn, sd(D, nn, rep)))
    n_perc = min(cap, sizes[-1])
    for rep in range(nread):
        jobs.append(("perc", 0, n_perc, derive_seed(ctx.seed, "K.1-perc", n_perc, rep) % (2 ** 31)))
    t0 = time.time()
    res = parallel_map(task, [(k, d, n, s, nint, nact, pp) for k, d, n, s in jobs], ctx.workers, ctx.progress,
                       lambda i, r_, dn, nn: f"[K.1] {jobs[i][0]} {jobs[i][1]-1 if jobs[i][1] else '-'}+1 N={jobs[i][2]} "
                                             f"done ({dn}/{nn}, {time.time()-t0:.0f}s)")
    reads, nulls, perc = {}, {}, []
    for (k, d, n, s), r_ in zip(jobs, res):
        if k == "read":
            reads.setdefault((d, n), []).append(r_)
        elif k == "nulls":
            nulls.setdefault((d, n), []).append(r_)
        else:
            perc.append(r_)

    def null_subjects(d, n):
        out = {}
        for key in ("tp", "seq", "layered"):
            sb = OR.subject({n: [r_[key] for r_ in nulls[(d, n)]]})
            out[key] = {"label": sb["label"], "numbers": sb["numbers"]}
        return out

    def decide_at(d, n, n_null):
        sub = OR.subject({1: reads[(d, n)]})
        return OR.decide(sub, null_subjects(d, n_null)), sub

    inst = []
    for n in inst_sizes:
        (dec, dstar), sub = decide_at(3, n, min(n, cap))
        inst.append({"n": n, "decision": dec, "d_star": dstar})
    inst_ok = all(x["decision"].startswith("GEOMETRIC") and x["d_star"] == 3 for x in inst)
    per = []
    for n in sizes:
        (dec, dstar), sub = decide_at(D, n, null_size(n))
        nl = null_subjects(D, null_size(n))
        per.append({"n": n, "decision": dec, "d_star": dstar, "null_size": null_size(n),
                    "D_MM": sub["numbers"].get("D_MM"), "D_chain": sub["numbers"].get("D_chain"),
                    "D_sd": sub["numbers"].get("D_sd"), "homog_spread": sub["homog_spread"], "dS4": sub["dS"].get(4),
                    "nulls": {k: v["label"] for k, v in nl.items()}})
    geo = [x["decision"].startswith("GEOMETRIC") and x["d_star"] == D for x in per]
    run_len = 0
    for g in reversed(geo):
        if g:
            run_len += 1
        else:
            break
    n_star = per[len(per) - run_len]["n"] if run_len else None
    (pdec, _), _ = (OR.decide(OR.subject({n_perc: perc}), {}), None)
    null_geo = pdec.startswith("GEOMETRIC")
    if not inst_ok:
        status = "INSTRUMENT_FAILED"
    elif null_geo:
        status = "UNDECIDED"
    elif run_len >= int(p["consecutive"]):
        status = "PASSED"
    else:
        status = "FAILED"
    fmt = lambda x, nd=2: "-" if x is None or (isinstance(x, float) and np.isnan(x)) else f"{x:.{nd}f}"
    lines = ["instrument (2+1, coordinate path, own matched nulls): " + "; ".join(
        f"N={x['n']} {x['decision'][:40]}" for x in inst) + f" -> {'OK' if inst_ok else 'FAILS'}"]
    for x in per:
        lines.append(f"  {D-1}+1 N={x['n']:<7d} {x['decision'][:44]:44s} D_MM {fmt(x['D_MM'])} D_chain {fmt(x['D_chain'])} "
                     f"D_sd {fmt(x['D_sd'])} homog {fmt(x['homog_spread'])} dS_4 {fmt(x['dS4'])} (nulls at {x['null_size']}: "
                     + ", ".join(f"{k} {v[:24]}" for k, v in x["nulls"].items()) + ")")
    lines += [f"  trailing GEOMETRIC(4) sizes: {run_len} (need {p['consecutive']}); N* = {n_star}",
              f"  null: percolation p = {pp:g} at N = {n_perc}: {pdec[:60]}"]
    data = {"instrument": inst, "per_size": per, "n_star": n_star, "trailing": run_len, "percolation": pdec,
            "reads": {f"{d}_{n}": v for (d, n), v in reads.items()}}
    return Result(status, lines, data)
