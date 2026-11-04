"""
K.11 -- The expected BD operator of a causal set in the geometry Pi^2 eta
=========================================================================

Instrument for the self-consistent field equation (a) chosen by the author
on 18 Oct 2026: B_g Pi = kappa rho_s, the operator of the causal set
sprinkled into the geometry g = Pi^2 eta that the field itself creates.
F1.4 (superposition) and F1.5 (scaling with the source) are built on it and
may only run with a quadrature this test has calibrated.

The instrument (core.operators.bd_expect_curved) computes E[B_g] by
quadrature: the causal order of g is flat, only the density changes to
rho Pi^D, so the expected kernel is evaluated at the PROPER volume of each
interval, V_g = V_flat * Phi, Phi = <Pi^D> over the flat diamond (computed
on the exact spatial projection of the diamond, a prolate spheroid).
With Pi = 1 it is exactly the flat operator of K.4 (core check
`curved_rows_flat_limit`).

Parts
-----
diamond_d3, diamond_d4 -- the diamond average Phi.
    Random diamonds (end point radius, spatial length up to 40 l_k,
    direction, proper time) in the test field Pi = 1 - 0.3 exp(-r^2/18)
    (width 3 l_k, the source width of F1.4/F1.5). Phi from the spheroid
    rule at the finest quadrature level against a Monte Carlo integral over
    the diamond sampled uniformly in its REST frame (an independent code
    path). Criterion: max |z| <= z_max over all diamonds (z = difference /
    Monte Carlo SE; z_max = 3.5 for 40 comparisons).

conv_d3, conv_d4 -- convergence of the operator in the quadrature.
    The curved operator applied to the same test field at r = 2 (inside
    the source), 10 and 30 l_k (outside), at four levels of quadrature
    refined together (directions, radial nodes, s nodes and the three
    spheroid directions). Background: the far rows are a small residue of
    large contributions that cancel along every ray; the residue depends on
    how well the proper volume of the rays crossing the source is resolved,
    and at the coarsest levels it moved by ~1e-3 when any of these
    parameters was doubled (development runs, 18 Oct; numbers reported here
    only so that they are on record).
    Criterion: between the two finest levels every row agrees within
    conv_tol times the operator's response inside the source,
    |E_L3(r) - E_L2(r)| <= conv_tol |E_L3(2)|, conv_tol = 0.01 (one order
    below the 10 % field accuracy of K.3, the accuracy the field tests need).
    Reported: the coarsest level within that tolerance of the finest
    (the quadrature F1.4 and F1.5 will use), and the continuum values of (a)
    (D Pi^-2 lap Pi + c_D Pi^-3 |grad Pi|^2), of box_g alone and of the flat
    Laplacian at the same radii.

local_d3, local_d4 -- the local limit box_g - R/2.
    Averaged causal-set d'Alembertians tend to box_g - R/2 in every
    dimension (Belenchia et al. 2016). For g = Pi^2 eta and static Pi that
    is (a)'s local form D Pi^-2 lap Pi + c_D Pi^-3 |grad Pi|^2 (c_3 = 0,
    c_4 = 2); box_g alone would give Pi^-2 lap Pi + (D - 2) Pi^-3 |grad Pi|^2.
    Field Pi = 1 - 0.3 exp(-r^2 / 2 w^2) for widths w = 4, 8, 16 l_k, row at
    r = w, quadrature level local_level. Ratio rho(w) = E / (a)-continuum.
    Criterion: (i) |1 - rho| decreases with w (approach to the limit) and
    (ii) at the largest w the ratio is at least three times closer to (a)
    than to box_g: |rho - 1| <= |rho_box - 1| / 4 (the quarter rule used for
    candidate laws in F1.3). Seen in development (coarse and medium
    quadrature): 2+1 rho = 0.77, 0.89, 0.95; 3+1 rho = 0.87, 0.955, 0.968
    (box_g: 0.41 and 0.29).

No null model: deterministic quadrature against independent references
(Monte Carlo, refinement, analytic local limit).
"""

from __future__ import annotations

import time

import numpy as np

from core.operators.bd_expect import shape_nodes
from core.operators.bd_expect_curved import curved_radial_rows, diamond_phi_spheroid, model_constants
from core.order.sprinkle import derive_seed
from core.stats.stats import fmt
from core.stats.parallel import parallel_map
from core.testing import Result

TEST = {
    "id": "K.11",
    "phase": "K",
    "title": "The expected BD operator in the geometry Π²η",
    "theory_part": "instrument for the self-consistent field equation (a) B_g Π = κρ_s (F1.4, F1.5); no physical "
                   "conclusion",
    "question": "Does the quadrature compute the expected BD operator of a causal set in the geometry Π²η correctly: "
                "the proper volume of intervals, convergence of the quadrature and the local limit □_g − R/2?",
    "method": "Φ = ⟨Π^D⟩ over a diamond against Monte Carlo in the rest frame of the diamond; the operator on a test "
              "field at four quadrature levels at r = 2, 10, 30 l_k; ratio to the continuum limit (a) for widths 4, "
              "8, 16 l_k",
    "prediction": "Φ agrees with Monte Carlo; the two finest levels agree within 1 % of the response inside the "
                  "source; the ratio to (a) approaches 1 and is closer to (a) than to □_g",
    "criterion": "diamond: max |z| <= z_max. conv: |E_L3(r) − E_L2(r)| <= conv_tol |E_L3(2)| for all r. local: |1 − "
                 "ρ| decreases with w and at the largest w |ρ − 1| <= |ρ_box − 1| / 4.",
    "null_model": "—",
    "level": "calibration",
    "meters": ["K.4"],
    "depends_on": ["K.4"],
    "params": {"amp": -0.3, "sigma": 3.0, "R": 80.0, "dr": 0.5, "ratio": 1.06, "rows": [2.0, 10.0, 30.0],
               "conv_tol": 0.01, "z_max": 3.5, "n_pairs": 40, "mc_samples": 200000,
               "widths": [4.0, 8.0, 16.0], "local_level": 2, "local_nodes": 241},
    "parts": {
        "diamond_d3": {"mode": "diamond", "D": 3, "levels": [[160, 8, 24, 4, 12, 1], [160, 8, 48, 6, 24, 1],
                                                              [240, 8, 96, 8, 48, 1], [320, 10, 144, 10, 64, 1]]},
        "diamond_d4": {"mode": "diamond", "D": 4, "levels": [[160, 8, 16, 4, 12, 4], [160, 8, 32, 6, 24, 6],
                                                              [240, 8, 64, 8, 48, 8], [320, 10, 96, 10, 64, 10]]},
        "conv_d3": {"mode": "conv", "D": 3, "levels": [[160, 8, 24, 4, 12, 1], [160, 8, 48, 6, 24, 1],
                                                        [240, 8, 96, 8, 48, 1], [320, 10, 144, 10, 64, 1]]},
        "conv_d4": {"mode": "conv", "D": 4, "levels": [[160, 8, 16, 4, 12, 4], [160, 8, 32, 6, 24, 6],
                                                        [240, 8, 64, 8, 48, 8], [320, 10, 96, 10, 64, 10]]},
        "local_d3": {"mode": "local", "D": 3, "levels": [[160, 8, 24, 4, 12, 1], [160, 8, 48, 6, 24, 1],
                                                          [240, 8, 96, 8, 48, 1], [320, 10, 144, 10, 64, 1]]},
        "local_d4": {"mode": "local", "D": 4, "levels": [[160, 8, 16, 4, 12, 4], [160, 8, 32, 6, 24, 6],
                                                          [240, 8, 64, 8, 48, 8], [320, 10, 96, 10, 64, 10]]},
    },
    "quick": {"n_pairs": 6, "mc_samples": 20000, "rows": [2.0, 10.0], "R": 30.0, "widths": [3.0, 5.0],
              "local_level": 1, "local_nodes": 61},
    "quick_parts": {
        "diamond_d3": {"levels": [[40, 6, 12, 3, 8, 1], [60, 6, 16, 4, 12, 1]]},
        "diamond_d4": {"levels": [[40, 6, 8, 3, 8, 3], [60, 6, 12, 4, 12, 4]]},
        "conv_d3": {"levels": [[40, 6, 12, 3, 8, 1], [60, 6, 16, 4, 12, 1]]},
        "conv_d4": {"levels": [[40, 6, 8, 3, 8, 3], [60, 6, 12, 4, 12, 4]]},
        "local_d3": {"levels": [[40, 6, 12, 3, 8, 1], [60, 6, 16, 4, 12, 1]]},
        "local_d4": {"levels": [[40, 6, 8, 3, 8, 3], [60, 6, 12, 4, 12, 4]]},
    },
    "estimate": "diamond: minutes; conv_d3 ~1 h, conv_d4 several hours on 1 core (the finest level dominates; 12 "
                "jobs in parallel); local: tens of minutes; run the parts separately with --workers",
}

# lengths in units of l_k: eps * rho = 1
_EPS, _RHO = 0.05, 20.0


def test_field(amp, sigma):
    """Pi = 1 + amp exp(-r^2 / 2 sigma^2) and its radial Laplacian and squared gradient (n spatial dimensions)."""
    def Pi(r):
        return 1.0 + amp * np.exp(-np.asarray(r, float) ** 2 / (2.0 * sigma ** 2))

    def lap(r, n):
        g = amp * np.exp(-r ** 2 / (2.0 * sigma ** 2))
        return g * (r ** 2 / sigma ** 4 - n / sigma ** 2)

    def grad2(r):
        g = amp * np.exp(-r ** 2 / (2.0 * sigma ** 2))
        return (g * r / sigma ** 2) ** 2
    return Pi, lap, grad2


def continuum_values(D, P, lap, g2):
    """Local operator values at one point: (a), box_g alone, flat Laplacian."""
    al, c, _ = model_constants("a", D)
    return {"a": al * lap / P ** 2 + c * g2 / P ** 3,
            "box_g": lap / P ** 2 + (D - 2) * g2 / P ** 3,
            "flat": lap}


# ---------------------------------------------------------------- diamond part

def mc_phi(D, r_x, L, CO, tau, amp, sigma, M, seed):
    """Monte Carlo <Pi^D> over the diamond, sampled uniformly in its rest frame (independent of the spheroid rule)."""
    rng = np.random.default_rng(seed)
    Pi, _, _ = test_field(amp, sigma)
    n = D - 1
    dt = np.sqrt(L * L + tau * tau)
    nvec = np.zeros(n)
    nvec[0], nvec[1] = CO, np.sqrt(max(1.0 - CO * CO, 0.0))
    xs = np.zeros(n)
    xs[0] = r_x
    ys = xs + L * nvec
    u = -nvec                                       # spatial direction from y to x
    perp = np.linalg.svd(u[None, :])[2][1:]         # rest-frame transverse axes (orthogonal to u)
    vals = []
    left = M
    while left > 0:
        m = min(left, 200000)
        T = tau * rng.random(m)
        X = (rng.random((m, n)) * 2.0 - 1.0) * tau / 2.0
        ok = np.linalg.norm(X, axis=1) <= np.minimum(T, tau - T)
        T, X = T[ok], X[ok]
        # rest-frame point y + (T/tau) Delta + X_L e_L + X_T e_T; spatial part of e_L is (dt/tau) u
        pos = ys + (T / tau)[:, None] * (L * u)[None, :] + X[:, :1] * (dt / tau) * u[None, :] + X[:, 1:] @ perp
        vals.append(Pi(np.linalg.norm(pos, axis=1)) ** D)
        left -= m
    v = np.concatenate(vals)
    return float(v.mean()), float(v.std(ddof=1) / np.sqrt(v.size))


def _diamond(p, ctx):
    D = int(p["D"])
    rng = np.random.default_rng(derive_seed(ctx.seed, "K.11-diamonds", D))
    k = int(p["n_pairs"])
    pairs = np.column_stack([rng.uniform(0.5, 12.0, k), np.exp(rng.uniform(np.log(0.1), np.log(40.0), k)),
                             rng.uniform(-1.0, 1.0, k), np.exp(rng.uniform(np.log(0.05), np.log(6.0), k))])
    args = [(D, float(a), float(b), float(c), float(d), float(p["amp"]), float(p["sigma"]), int(p["mc_samples"]),
             derive_seed(ctx.seed, "K.11-mc", D, i)) for i, (a, b, c, d) in enumerate(pairs)]
    t0 = time.time()
    mc = parallel_map(mc_phi, args, ctx.workers, ctx.progress,
                      lambda i, r_, d, n_: f"[K.11] {D-1}+1 Monte Carlo diamond {d}/{n_} ({time.time()-t0:.0f}s)")
    Pi, _, _ = test_field(float(p["amp"]), float(p["sigma"]))
    lines = [f"{D-1}+1: {k} random diamonds in Pi = 1 {p['amp']:+g} exp(-r^2/2*{p['sigma']:g}^2); Monte Carlo "
             f"{p['mc_samples']} rest-frame samples each"]
    zs_by_level = []
    for li, lev in enumerate(p["levels"]):
        n_mu, n_nu, n_phi = int(lev[3]), int(lev[4]), int(lev[5])
        z = []
        for (r_x, L, CO, tau), (m, se) in zip(pairs, mc):
            dt = np.sqrt(L * L + tau * tau)
            ph = float(diamond_phi_spheroid(Pi, D, r_x, np.array([[L]]), CO, np.array([[dt]]), np.array([[tau]]),
                                            n_mu, n_nu, n_phi)[0, 0])
            z.append((ph - m) / se)
        z = np.array(z)
        zs_by_level.append(z)
        lines.append(f"  level {li} (n_mu {n_mu}, n_nu {n_nu}, n_phi {n_phi}): max |z| = {np.max(np.abs(z)):.2f}, "
                     f"RMS z = {np.sqrt(np.mean(z ** 2)):.2f}")
    zf = zs_by_level[-1]
    ok = bool(np.max(np.abs(zf)) <= float(p["z_max"]))
    lines.append(f"  finest level: max |z| {np.max(np.abs(zf)):.2f} (z_max {p['z_max']}) -> {'OK' if ok else 'FAILS'}")
    return ok, lines, {"pairs": pairs, "mc": mc, "z": zs_by_level}


# ------------------------------------------------------------ convergence part

def conv_task(D, R, nodes, Pi_nodes, r_row, lev):
    """One operator value E[B_g Pi](r_row) at one quadrature level (module level: worker processes)."""
    row = curved_radial_rows(D, _RHO, _EPS, 3.0 * R, [r_row], R, nodes, Pi_nodes, *[int(x) for x in lev])[0]
    return float(row @ Pi_nodes)


def _conv(p, ctx):
    D, n = int(p["D"]), int(p["D"]) - 1
    R, sig, amp = float(p["R"]), float(p["sigma"]), float(p["amp"])
    nodes = shape_nodes(R, 2.0 * sig, float(p["dr"]), float(p["ratio"]))
    Pi, lap, g2 = test_field(amp, sig)
    Pin = Pi(nodes)
    rows = [float(x) for x in p["rows"]]
    levels = p["levels"]
    args = [(D, R, nodes, Pin, r, lev) for lev in levels for r in rows]
    t0 = time.time()
    vals = parallel_map(conv_task, args, ctx.workers, ctx.progress,
                        lambda i, r_, d, n_: f"[K.11] {D-1}+1 level {i // len(rows)}, r = {rows[i % len(rows)]:g}: "
                                             f"{r_:+.5e} ({d}/{n_}, {time.time()-t0:.0f}s)")
    E = np.array(vals).reshape(len(levels), len(rows))
    scale = abs(E[-1, 0])
    tol = float(p["conv_tol"]) * scale
    ok = bool(np.all(np.abs(E[-1] - E[-2]) <= tol))
    within = [li for li in range(len(levels)) if np.all(np.abs(E[li] - E[-1]) <= tol)]
    chosen = min(within) if within else None
    cont = [continuum_values(D, float(Pi(r)), float(lap(r, n)), float(g2(r))) for r in rows]
    lines = [f"{D-1}+1: E[B_g Pi] for Pi = 1 {amp:+g} exp(-r^2/2*{sig:g}^2), R = {R:g}, {nodes.size} nodes; rows r = "
             + ", ".join(f"{r:g}" for r in rows)]
    for li, lev in enumerate(levels):
        lines.append(f"  level {li} {list(lev)}: " + "  ".join(f"{v:+.5e}" for v in E[li]))
    lines += ["  continuum (a):       " + "  ".join(f"{c['a']:+.5e}" for c in cont),
              "  continuum box_g:     " + "  ".join(f"{c['box_g']:+.5e}" for c in cont),
              "  continuum flat lap:  " + "  ".join(f"{c['flat']:+.5e}" for c in cont),
              f"  |E_L3 - E_L2| = " + ", ".join(f"{abs(a - b):.2e}" for a, b in zip(E[-1], E[-2]))
              + f" (tol {p['conv_tol']} x |E(r={rows[0]:g})| = {tol:.2e}) -> {'OK' if ok else 'FAILS'}",
              f"  coarsest level within tolerance of the finest: {chosen if chosen is not None else 'none'}"
              + (f" {list(levels[chosen])}" if chosen is not None else "")]
    return ok, lines, {"E": E, "rows": rows, "levels": levels, "chosen_level": chosen, "continuum": cont,
                       "tolerance": tol}


# ------------------------------------------------------------ local limit part

def local_task(D, w, amp, n_nodes, lev):
    """Ratio of the curved operator to the continuum values at r = w for the test field of width w."""
    R = 6.0 * w
    nodes = np.linspace(0.0, R, int(n_nodes))
    Pi, lap, g2 = test_field(amp, w)
    row = curved_radial_rows(D, _RHO, _EPS, 3.0 * R, [w], R, nodes, Pi(nodes), *[int(x) for x in lev])[0]
    E = float(row @ Pi(nodes))
    c = continuum_values(D, float(Pi(w)), float(lap(w, D - 1)), float(g2(w)))
    return {"E": E, "a": c["a"], "box_g": c["box_g"], "flat": c["flat"]}


def _local(p, ctx):
    D = int(p["D"])
    lev = p["levels"][min(int(p["local_level"]), len(p["levels"]) - 1)]
    widths = [float(w) for w in p["widths"]]
    t0 = time.time()
    res = parallel_map(local_task, [(D, w, float(p["amp"]), int(p["local_nodes"]), lev) for w in widths], ctx.workers,
                       ctx.progress, lambda i, r_, d, n_: f"[K.11] {D-1}+1 local limit w = {widths[i]:g} done "
                                                          f"({d}/{n_}, {time.time()-t0:.0f}s)")
    rho = np.array([r_["E"] / r_["a"] for r_ in res])
    rho_box = np.array([r_["box_g"] / r_["a"] for r_ in res])
    dev = np.abs(1.0 - rho)
    mono = bool(np.all(np.diff(dev) < 0))
    closer = bool(dev[-1] <= abs(rho_box[-1] - 1.0) / 4.0)
    lines = [f"{D-1}+1: Pi = 1 {p['amp']:+g} exp(-r^2/2w^2), row at r = w, quadrature level {list(lev)}",
             "  w            " + "  ".join(f"{w:8g}" for w in widths),
             "  E / (a)      " + "  ".join(f"{x:8.4f}" for x in rho),
             "  box_g / (a)  " + "  ".join(f"{x:8.4f}" for x in rho_box),
             f"  |1 - rho| decreases with w: {'yes' if mono else 'NO'}; at w = {widths[-1]:g}: |rho - 1| = {dev[-1]:.4f} "
             f"<= |rho_box - 1|/4 = {abs(rho_box[-1] - 1.0) / 4.0:.4f}: {'yes' if closer else 'NO'}"]
    return mono and closer, lines, {"widths": widths, "results": res, "rho": rho, "rho_box": rho_box}


def run(p, ctx):
    mode = p["mode"]
    if mode == "diamond":
        ok, lines, data = _diamond(p, ctx)
    elif mode == "conv":
        ok, lines, data = _conv(p, ctx)
    else:
        ok, lines, data = _local(p, ctx)
    return Result("PASSED" if ok else "FAILED", lines, data)
