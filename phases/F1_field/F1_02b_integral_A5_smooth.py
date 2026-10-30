"""
F1.2b -- The integral field equation A5' in one ordering, smooth source (revision of F1.2)
=========================================================================================

Why a revision (18 Oct 2026). F1.2 ended INSTRUMENT_FAILED: its instrument
check "residual of the field solution <= res_tol" read 0.0875 q. Cause: the
source of F1.2 is a ball with a sharp edge (radius a = 0.05), while the
expected operator is nonlocal on l_k = 0.05 (2+1) and 0.095 (3+1). Its image
is smooth on that scale, so no field solves E[B] u = q 1(r < a) at the edge;
a finer radial basis does not help (40 nodes: 0.062 q). The tolerance was
fixed without measuring the residual on the full grid (the quick runs had a
square system, residual 0) -- a design error, not a property of the theory.
F1.2 stays recorded as INSTRUMENT_FAILED under its own criterion.

Change: the source is a Gaussian q exp(-r^2 / 2 sigma^2), sigma = 0.05 (the
size of the F1.2 source, >= l_k / 2): the measured residual on the full grid
is 0.003 q (2+1) and 0.008 q (3+1), below the UNCHANGED tolerance 0.05.
Everything else -- question, aggregates, calibration, null model and all
thresholds -- is unchanged. New seed streams ("F1.2b", "F1.2b-cal"), so no
sprinkling of F1.2 is reused. Seen before this revision: the F1.2 run (A met,
max |z| 1.83; B met, p_test 1.010 +- 0.161 vs p_cal 0.737 +- 0.246; null
rejected, max |z| 21); it plays no role in F1.2b's decision.

Original description of F1.2 (with the ball source replaced by the Gaussian):

Register R133-R135 (H-N1, H-N2, A5' as a mesoscopic law). The discrete BD
operator of ONE sprinkling is noise-dominated element by element (F1.6,
E45b: INSTRUMENT_FAILED), and its expectation is a nonlocal operator on the
scale l_k = (eps rho)^(-1/d) (R133). A5' is therefore read as a law for
aggregated quantities: averaged over a region of scale L, the operator of
one ordering applied to the field should approach the EXPECTED operator,
with an error that falls as (L / ell)^(-p) (H-N2, R134).

What is computed
----------------
Field. u = Pi - Pi_vac for a Gaussian source of width sigma and peak density q (A5' in
kernel units, kappa = 1): u_e is the solution of the EXPECTED nonlocal
equation E[B] u = q exp(-r^2 / 2 sigma^2) with u(R) = 0, from expected rows on a fixed
grid of (t, r) points of the time band (deterministic quadrature, the same
geometry as the rows; no sprinkling enters).

Rows. In each sprinkling, stratified elements of the band [t0, t1]
(rows_per_element per radial element below r_rows). For each row x: the
sampled operator applied to the field, (B u)(x), and the expected operator at
the same element, (E[B] u)(x).

Aggregate. For a region |x| < L of the band, the VOLUME average of the
row differences d(x) = (B u)(x) - (E[B] u)(x), each row weighted by the
volume of its radial shell divided by the number of rows in it:
delta(L) = <d>_L, in units of q. One value per sprinkling and L.

Error. epsilon(L) = root mean square of delta(L) over sprinklings; its
slope against L / ell on a log-log scale gives the aggregation exponent
p = -d log epsilon / d log(L / ell).

Calibration (non-circular, independent of the decision data):
  * the same procedure on n_cal INDEPENDENT sprinklings (seed stream
    "F1.2b-cal") with a SOURCE-FREE smooth field u_cal = A ((r/R)^2 - 1),
    A = max |u_e| (same amplitude), gives p_cal and checks that the
    calibration aggregates are unbiased (instrument);
  * the expected operator (u_e, expected rows) is computed by quadrature
    (K.4), independently of every sprinkling;
  * nothing is fitted on the decision data: no normalisation, no kappa,
    no exponent.
If the aggregation exponent depended on the shape of the field, the
decision would show it as a mismatch of p: that is a failure of the
universality claimed by H-N2, not an instrument fault.

Null model. The same rows with the kernel weights K_eps(n_xy) of each row
randomly permuted among its candidate ancestors (core.operators.bd,
null_rng): same weight distribution, but the causal structure no longer
says which y gets which weight. It must NOT reproduce the expected
operator in aggregate.

Not part of this test (separate levels, user's decision of 18 Oct 2026):
the asymptotic shape Pi - Pi_vac ~ -1/r for r >> l_k (F1.3, derived in
BD-3), superposition (F1.4) and scaling with the source (F1.5).
"""

from __future__ import annotations

import time

import numpy as np

from core.order.sprinkle import sprinkle, Cylinder, derive_seed
from core.operators.bd import bd_static_rows, stratified_rows
from core.operators.bd_expect import expected_radial_rows
from core.operators.solve import solve_radial_ls
from core.stats.stats import mean_se, slope_loglog, fmt
from core.stats.parallel import parallel_map
from core.testing import Result

TEST = {
    "id": "F1.2b",
    "phase": "F1",
    "title": "Integral A5′ in one order, smooth source",
    "theory_part": "axiom 5 and A5′ as a mesoscopic law (R133–R135, H-N2)",
    "question": "Does the aggregated BD operator of one order applied to the field of a source approach the expected "
                "operator, with an error falling as (L/ℓ)^(−p) and the exponent p from an independent calibration?",
    "method": "Gaussian source of width sigma; field u_e from the expected non-local equation (quadrature); "
              "stratified rows in a time band; volume average of the difference sampled − expected operator over "
              "regions |x| < L; RMS over sprinklings; exponent from the slope; calibration on independent "
              "sprinklings with a field without a source",
    "prediction": "H-N2: the aggregates are unbiased (mean difference 0) at every L; the error falls with the "
                  "exponent p_cal of the calibration; the permuted null does not reproduce the expected operator",
    "criterion": "INSTRUMENT: residual of the solution u_e <= res_tol (in units of q); calibration aggregates "
                 "unbiased (max |z| <= z_max); p_cal >= p_min. A (agreement): max_L |z_test(L)| <= z_max. B "
                 "(scaling): p_test >= p_min and |p_test - p_cal| <= max(p_tol, n_se * sqrt(SE_test^2 + SE_cal^2)), "
                 "SE by bootstrap over sprinklings. PASSED: A and B and the null does not meet A. FAILED: A or B "
                 "does not hold. UNDECIDED: the null meets A (unspecific).",
    "null_model": "the same rows with the kernel weights randomly permuted among the candidate ancestors of each row",
    "level": "mechanism",
    "meters": ["K.4", "K.5"],
    "depends_on": ["F1.1"],
    "revision_of": "F1.2",
    "old_id": "F1.2 rev.",
    "params": {"res_tol": 0.05, "z_max": 3.0, "p_min": 0.1, "p_tol": 0.15, "n_se": 2.0, "boot_B": 2000,
               "T": 1.0, "R": 0.5, "sigma": 0.05, "q": 1.0, "eps": 0.2, "n_nodes": 20, "t0": 0.8, "t1": 1.0,
               "r_rows": 0.45, "L_min": 0.1, "n_L": 6, "grid_t": [0.85, 0.9, 0.95], "grid_per_element": 3},
    "parts": {"d3": {"D": 3, "N": 32000, "rows_per_element": 20, "n_cal": 16, "n_test": 16},
              "d4": {"D": 4, "N": 32000, "n_nodes": 12, "rows_per_element": 6, "n_cal": 12, "n_test": 12}},
    "quick": {"N": 2500, "rows_per_element": 2, "n_cal": 3, "n_test": 3, "boot_B": 200, "grid_per_element": 1,
              "grid_t": [0.9]},
    "estimate": "2+1: ~1.7 h on 1 core (32 sprinklings of ~3 min); 3+1: ~35 min on 1 core; with more processes "
                "proportionally less; run K.4 and K.5 first",
}


def expected_field(p, rho):
    """u_e on the radial nodes from expected rows on the fixed (t, r) grid; residual in units of q."""
    D, R = int(p["D"]), float(p["R"])
    nodes = np.linspace(0.0, R, int(p["n_nodes"]) + 1)
    k = int(p["grid_per_element"])
    fr = (np.arange(k) + 0.5) / k
    rr = (nodes[:-1, None] + fr[None, :] * np.diff(nodes)[:, None]).ravel()
    rr = rr[rr < float(p["r_rows"])]
    A = np.vstack([expected_radial_rows(D, rho, float(p["eps"]), float(tg), rr, R, nodes) for tg in p["grid_t"]])
    rg = np.tile(rr, len(p["grid_t"]))
    rhs = float(p["q"]) * np.exp(-rg ** 2 / (2.0 * float(p["sigma"]) ** 2))
    u = solve_radial_ls(A, rhs)
    res = float(np.sqrt(np.mean((A @ u - rhs) ** 2))) / float(p["q"])
    return nodes, u, res


def one_sprinkling(p, seed, field, with_null):
    """Volume-averaged row differences delta(L) (and the null's) for one sprinkling and one field vector."""
    D, R = int(p["D"]), float(p["R"])
    S = sprinkle(Cylinder(D, 0.0, float(p["T"]), R), int(p["N"]), seed)
    r = S.r()
    nodes = np.linspace(0.0, R, int(p["n_nodes"]) + 1)
    rows = stratified_rows(S.X, r, float(p["t0"]), float(p["t1"]), nodes, float(p["r_rows"]),
                           int(p["rows_per_element"]), np.random.default_rng(derive_seed(seed, "rows")))
    B = bd_static_rows(S.X, S.rho, float(p["eps"]), rows, nodes, R,
                       null_rng=np.random.default_rng(derive_seed(seed, "null")) if with_null else None)
    Ae = np.vstack([expected_radial_rows(D, S.rho, float(p["eps"]), float(S.X[i, 0]), [float(r[i])], R, nodes)
                    for i in rows])
    f = np.asarray(field, float)
    e_val = Ae @ f
    d = B["A"] @ f - e_val
    dn = (B["A_null"] @ f - e_val) if with_null else None
    # volume weights: shell volume of the row's radial element / number of rows in it
    rr = r[rows]
    k = np.clip(np.searchsorted(nodes, rr, side="right") - 1, 0, len(nodes) - 2)
    shell = nodes[k + 1] ** (D - 1) - nodes[k] ** (D - 1)
    cnt = np.bincount(k, minlength=len(nodes) - 1)[k]
    w = shell / np.maximum(cnt, 1)
    Ls = np.geomspace(float(p["L_min"]), float(p["r_rows"]), int(p["n_L"]))
    out = {"delta": [], "delta_null": [], "expected": [], "rows": int(rows.size), "ell": S.ell}
    for L in Ls:
        m = rr < L
        sw = w[m].sum()
        out["delta"].append(float((w[m] * d[m]).sum() / sw) if sw > 0 else float("nan"))
        out["expected"].append(float((w[m] * e_val[m]).sum() / sw) if sw > 0 else float("nan"))
        if with_null:
            out["delta_null"].append(float((w[m] * dn[m]).sum() / sw) if sw > 0 else float("nan"))
    return out


def stats(deltas: np.ndarray, Ls: np.ndarray, ell: float, B: int, seed: int):
    """Bias, SE and z per L; RMS error per L; exponent p and its bootstrap SE (resampling sprinklings)."""
    n = deltas.shape[0]
    bias = deltas.mean(axis=0)
    se = deltas.std(axis=0, ddof=1) / np.sqrt(n) if n > 1 else np.full(deltas.shape[1], np.nan)
    z = bias / se
    eps_L = np.sqrt((deltas ** 2).mean(axis=0))
    p = -slope_loglog(Ls / ell, eps_L)
    rng = np.random.default_rng(seed)
    ps = []
    for _ in range(B):
        idx = rng.integers(0, n, size=n)
        ps.append(-slope_loglog(Ls / ell, np.sqrt((deltas[idx] ** 2).mean(axis=0))))
    ps = np.asarray([x for x in ps if np.isfinite(x)])
    return {"bias": bias, "se": se, "z": z, "eps": eps_L, "p": p, "p_se": float(ps.std(ddof=1)) if ps.size > 1 else float("nan")}


def run(p, ctx):
    D = int(p["D"])
    t0 = time.time()
    reg = Cylinder(D, 0.0, float(p["T"]), float(p["R"]))
    rho = int(p["N"]) / reg.volume()
    ell = rho ** (-1.0 / D)
    l_k = (float(p["eps"]) * rho) ** (-1.0 / D)
    nodes, u_e, res = expected_field(p, rho)
    amp = float(np.max(np.abs(u_e)))
    u_cal = amp * ((nodes / float(p["R"])) ** 2 - 1.0)
    ctx.progress(f"[F1.2b] {D-1}+1 N={p['N']}: expected field solved (residual {res:.2e} q) ({time.time()-t0:.0f}s)")
    n_cal, n_test = int(p["n_cal"]), int(p["n_test"])
    args = [(p, derive_seed(ctx.seed, "F1.2b-cal", D, int(p["N"]), k), u_cal, False) for k in range(n_cal)]
    args += [(p, derive_seed(ctx.seed, "F1.2b", D, int(p["N"]), k), u_e, True) for k in range(n_test)]
    res_all = parallel_map(one_sprinkling, args, ctx.workers, ctx.progress,
                           lambda i, r_, d, n_: f"[F1.2b] {'calibration' if i < n_cal else 'test'} sprinkling "
                                                f"{(i if i < n_cal else i - n_cal) + 1} done ({d}/{n_}) ({time.time()-t0:.0f}s)")
    Ls = np.geomspace(float(p["L_min"]), float(p["r_rows"]), int(p["n_L"]))
    cal = np.array([r_["delta"] for r_ in res_all[:n_cal]])
    tst = np.array([r_["delta"] for r_ in res_all[n_cal:]])
    nul = np.array([r_["delta_null"] for r_ in res_all[n_cal:]])
    exp_test = np.mean([r_["expected"] for r_ in res_all[n_cal:]], axis=0)
    B = int(p["boot_B"])
    sc = stats(cal, Ls, ell, B, derive_seed(ctx.seed, "F1.2b-boot", "cal"))
    st = stats(tst, Ls, ell, B, derive_seed(ctx.seed, "F1.2b-boot", "test"))
    sn = stats(nul, Ls, ell, B, derive_seed(ctx.seed, "F1.2b-boot", "null"))
    zmax = float(p["z_max"])
    inst = {"field residual": res <= p["res_tol"],
            "calibration unbiased": bool(np.all(np.isfinite(sc["z"])) and np.max(np.abs(sc["z"])) <= zmax),
            "calibration exponent defined": bool(np.isfinite(sc["p"]) and sc["p"] >= p["p_min"])}
    A_ok = bool(np.all(np.isfinite(st["z"])) and np.max(np.abs(st["z"])) <= zmax)
    dp = abs(st["p"] - sc["p"])
    tol_p = max(float(p["p_tol"]), float(p["n_se"]) * float(np.sqrt(st["p_se"] ** 2 + sc["p_se"] ** 2)))
    B_ok = bool(np.isfinite(st["p"]) and st["p"] >= p["p_min"] and dp <= tol_p)
    null_A = bool(np.all(np.isfinite(sn["z"])) and np.max(np.abs(sn["z"])) <= zmax)
    if not all(inst.values()):
        status = "INSTRUMENT_FAILED"
    elif not (A_ok and B_ok):
        status = "FAILED"
    elif null_A:
        status = "UNDECIDED"
    else:
        status = "PASSED"
    row = lambda v, nd=4: " ".join(fmt(x, nd) for x in v)
    lines = [f"{D-1}+1 N={p['N']} eps={p['eps']}: ell = {ell:.4f}, l_k = {l_k:.4f}; Gaussian source sigma = {p['sigma']}, q = {p['q']}; "
             f"band t in [{p['t0']}, {p['t1']}]; {res_all[0]['rows']} rows per sprinkling; {n_cal} calibration + {n_test} test sprinklings",
             f"  L            {row(Ls, 3)}",
             f"  L / ell      {row(Ls / ell, 2)}",
             f"  expected <E[B]u_e>_L (units of q) {row(exp_test)}",
             f"  calibration (field without source): bias {row(sc['bias'])}",
             f"               z {row(sc['z'], 2)};  RMS error {row(sc['eps'])};  p_cal = {fmt(sc['p'])} +- {fmt(sc['p_se'])}",
             f"  test (source field u_e):            bias {row(st['bias'])}",
             f"               z {row(st['z'], 2)};  RMS error {row(st['eps'])};  p_test = {fmt(st['p'])} +- {fmt(st['p_se'])}",
             f"  null (permuted kernel weights):     bias {row(sn['bias'])}",
             f"               z {row(sn['z'], 2)}",
             f"  A agreement (max |z_test| <= {zmax}): {'met' if A_ok else 'NOT met'} (max |z| {fmt(np.nanmax(np.abs(st['z'])), 2)})",
             f"  B scaling: p_test {fmt(st['p'])} vs p_cal {fmt(sc['p'])}, |dp| = {fmt(dp)} (allowed {fmt(tol_p)}), "
             f"p_min {p['p_min']} -> {'met' if B_ok else 'NOT met'}",
             f"  null meets A: {'YES (not specific)' if null_A else 'no (good)'} (max |z_null| {fmt(np.nanmax(np.abs(sn['z'])), 2)})",
             "instrument: " + "; ".join(f"{k}: {'OK' if v else 'FAILS'}" for k, v in inst.items())
             + f" (field residual {res:.2e} q, tol {p['res_tol']})"]
    data = {"ell": ell, "l_k": l_k, "rho": rho, "L": Ls, "nodes": nodes, "u_e": u_e, "u_cal": u_cal, "field_residual": res,
            "expected": exp_test, "calibration": sc, "test": st, "null": sn, "instrument": inst,
            "agreement": A_ok, "scaling": B_ok, "null_meets": null_A,
            "deltas": {"calibration": cal, "test": tst, "null": nul}}
    return Result(status, lines, data)
