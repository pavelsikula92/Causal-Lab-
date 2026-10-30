"""
F1.7 -- Benchmark of the field equation B_g f(Pi) = kappa rho_s in 3+1 (continuum, static, spherical)
======================================================================================================

Why this test exists
--------------------
Equation (a), B_g Pi = kappa rho_s with B_g = box_g - R/2 (the continuum limit
of the Benincasa-Dowker operator, K.11), has no Newtonian response in 3+1 in
the geometry H1 (gamma = 1): in the weak field (Pi = 1 + u)
    box_g Pi -> + lap u,     -(R/2) Pi -> (1 - gamma (d_s - 1)) lap u,
    lambda_1 = 2 - gamma (d_s - 1) = 0   for d_s = 3, gamma = 1
(d_s spatial dimensions, D = d_s + 1). The author's decision (22 Oct 2026,
docs/author_decisions.md, point 1): keep the BD operator and gamma = 1,
change what the operator acts on: B_g f(Pi) = kappa rho_s, where
    B_g f(Pi) = f'(Pi) box_g Pi + f''(Pi) (grad Pi)^2 - (R/2) f(Pi).
Every f with f(1) = 0 and f'(1) = 1 has lambda_1 = 1, so the weak field does not
tell them apart:
    primary candidate   f = ln Pi        (a candidate, not an axiom)
    linear control      f = Pi - 1
    nonlinear controls  f = (Pi^k - 1)/k,  k = -1, 2, 3
    negative control    f = Pi           (equation (a): lambda_1 = 0)
This test is the author's benchmark for them: the same source, the same
boundary conditions, the same weak-field normalisation (kappa = 1).

Geometry and the local form
---------------------------
Metric ds^2 = -Pi^2 dt^2 + Pi^(-2 gamma) dx^2, gamma = 1 (H1), static, d_s = 3.
The local forms of box_g Pi and R are derived symbolically at the start of
every run (sympy, docs/derivations/local_form_symbolic.py) and printed:
    box_g Pi = Pi^(2g) [ lap Pi - (g - 1) Pi^-1 |grad Pi|^2 ],
    R        = Pi^(2g) [ 2 (2g - 1) Pi^-1 lap Pi - 2 g (g + 1) Pi^-2 |grad Pi|^2 ],
    (grad Pi)^2 = g^ij d_i Pi d_j Pi = Pi^(2g) |grad Pi|^2,
so the equation is Pi^(2g) [ P(Pi) lap Pi + Q(Pi) |grad Pi|^2 ] = kappa rho_s with
    P = f' - (2g - 1) f / Pi,    Q = -(g - 1) f' / Pi + f'' + g (g + 1) f / Pi^2.
For f = Pi, g = 1: P = 0 and the equation is 2 Pi |grad Pi|^2 = kappa rho_s (C5).
Radial: lap Pi = Pi'' + (2/r) Pi', |grad Pi|^2 = Pi'^2.

Source and boundary conditions
------------------------------
Gaussian rho_s = q exp(-r^2 / 2 sigma^2). Strength is labelled by the central
value of the LINEAR solution u_lin(0) = -kappa M sqrt(2/pi) / (4 pi sigma),
M = q (2 pi)^(3/2) sigma^3 (the same q for every candidate at a label).
Pi'(r_min) = 0 at r_min = 1e-4 sigma; at r = R = 50 sigma the exterior
condition Pi(R) - 1 + R Pi'(R) = 0 (u proportional to 1/r). scipy solve_bvp,
tolerance 1e-8, continuation in strength from the weak field.

Parts and criteria (stated before the runs)
-------------------------------------------
V  validation of the instrument:
   (i)  the separated symbolic forms reproduce C5 ((box_g - R/2) Pi = 2 Pi
        |grad Pi|^2 at g = 1) and lambda_1(g) = 2 - g (d_s - 1);
   (ii) the solver reproduces the analytic linear solution
        u = -kappa M erf(r / sqrt(2) sigma) / (4 pi r) of lap u = kappa rho_s:
        max |u - u_an| <= 1e-6 max |u_an|;
   (iii) the continuum has no length scale: ln Pi at u_lin(0) = -0.3 for
        sigma = 4 and sigma = 16 (q scaled) gives the same Pi(r / sigma)
        within 1e-6.
   PASSED if (i)-(iii) hold.
W  weak field (u_lin(0) = -1e-4): every normalised candidate agrees with the
   linear solution, max |u - u_lin| <= 0.01 max |u_lin|; the negative control
   f = Pi has a zero linear coefficient Pi^(2g) P(Pi) at Pi = 1 (no weak-field
   response). PASSED if both.
N  nonlinear (comparison, RECORDED): strengths u_lin(0) = -0.03, -0.1, -0.2,
   -0.3, -0.5, -0.7, -0.9 for every normalised candidate. Existence is
   decided by SHOOTING, not by a boundary-value solver (a solver that stops
   converging cannot tell "no solution" from "numerical trouble", and a
   branch of solutions can turn back at a largest source strength): for
   each strength, Pi(0) is scanned over (0, 1) (Pi'(0) = 0, the regular
   start Pi''(0) = kappa rho(0) / (3 Pi0^2g P(Pi0))), the radial equation is
   integrated to R (stopping if Pi <= 0 or P <= 0: no regular solution from
   that start), and the roots of the exterior condition Pi(R) - 1 + R Pi'(R)
   are counted. A strength has a regular solution if a root exists whose
   profile has Pi > 0, P > 0 and Pi non-decreasing; the weak-field branch is
   the root with the largest Pi(0). Reported per strength: the number of
   regular solutions, Pi(0) of the weak-field branch, the exterior charge
   R^2 Pi'(R) against the linear kappa M / 4 pi; per candidate the largest
   declared strength with a regular solution; the differences between the
   candidates.
   No single verdict: the decision between the candidates is the author's
   (and needs check 3, below).

Benchmark points of the author and what this test covers
--------------------------------------------------------
(1) weak field: part W; (2) a regular nonlinear profile returning to Pi -> 1:
part N; (4) acceptability when the size and the strength change: the
continuum equation has no length scale (part V iii shows it), so the size
enters only through the strength u_lin(0), which part N scans.
(3) the second potential psi and gamma = 1 is NOT tested here: with the
metric given by the ansatz H1, psi = phi and gamma = 1 are input, so the
check would be a tautology. A real test needs both potentials read
independently (clocks for phi; spatial distances or cone widths of the
dynamic order, F6.3b, for psi).

Seen before this test was written: the analytic split of lambda_1 (22 Oct
2026) and the symbolic forms above (derived for this file). No run of the
nonlinear equations.
"""

from __future__ import annotations

import math
import time

import numpy as np
from scipy.integrate import solve_bvp
from scipy.special import erf

from core.testing import Result

TEST = {
    "id": "F1.7",
    "phase": "F1",
    "title": "Benchmark of the field equation B_g f(Π) = κρ_s in 3+1 (continuum): ln Π against controls",
    "theory_part": "field equation: the BD operator B_g = □_g − R/2 acting on f(Π) (author's decision of 22 Oct), H1 "
                   "geometry",
    "question": "Do the candidates f = ln Π, Π − 1 and (Π^k − 1)/k in 3+1 in the H1 geometry give the same "
                "weak-field response and a regular non-linear solution for a static spherically symmetric source, "
                "and how do they differ non-linearly? (f = Π is a negative control with zero linear response.)",
    "method": "local form of □_g Π and R from the metric −Π²dt² + Π^(−2γ)dx² (γ = 1, d_s = 3) symbolically; radial "
              "problem Π²[P lap Π + Q |∇Π|²] = κρ_s, Gaussian source, Π'(0) = 0, outer condition Π − 1 + rΠ' = 0 at "
              "R = 50σ, solve_bvp with continuation in the strength; strength by the centre of the linear solution "
              "u_lin(0)",
    "prediction": "weak field: all normalised candidates agree with the linear solution, f = Π without a response; "
                  "non-linearly no prediction",
    "criterion": "V PASSED: the symbolic forms give C5 and λ₁(γ) = 2 − γ(d_s − 1); solver = analytic linear solution "
                 "to 1e-6; scaling σ = 4 and 16 agrees to 1e-6. W PASSED: normalised candidates at u_lin(0) = −1e-4 "
                 "deviate from the linear solution by <= 1 %, f = Π linear coefficient 0. N RECORDED: strengths "
                 "−0.03 to −0.9, for every candidate regularity (a solution exists, Π > 0, monotonic, P > 0) and a "
                 "comparison of Π(0) and of the outer charge; the author decides.",
    "null_model": "f = Π (equation (a)): no weak-field response in 3+1 at γ = 1",
    "level": "mechanism",
    "meters": [],
    "depends_on": ["F1.5"],
    "params": {"gamma": 1.0, "d_s": 3, "kappa": 1.0, "sigma": 8.0, "R_over_sigma": 50.0, "rmin_over_sigma": 1e-4,
               "n_shoot": 160, "shoot_rtol": 1e-10,
               "tol": 1e-8, "max_nodes": 200000, "lin_tol": 1e-6, "scale_tol": 1e-6, "weak_u0": -1e-4, "weak_tol": 0.01,
               "strengths": [-0.03, -0.1, -0.2, -0.3, -0.5, -0.7, -0.9], "ks": [-1.0, 2.0, 3.0], "scale_sigmas": [4.0, 16.0],
               "scale_u0": -0.3},
    "parts": {"V": {"mode": "V"}, "W": {"mode": "W"}, "N": {"mode": "N"}},
    "quick": {"tol": 1e-6, "strengths": [-0.03, -0.3], "lin_tol": 1e-4, "scale_tol": 1e-4, "n_shoot": 40},
    "estimate": "seconds to minutes (no causal set)",
}


# ------------------------------------------------------------------------------------------------- local form

def local_forms(d_s: int):
    """Symbolic (box_g Pi, R) coefficients: box = Pi^2g [a lap + b |grad|^2], R = Pi^2g [c lap + e |grad|^2] (functions of p, g)."""
    import sys
    from pathlib import Path
    import sympy as sp
    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "docs" / "derivations"))
    import local_form_symbolic as L
    gam = sp.Symbol("gamma")
    coords, xs, Pi, g, ginv, sqrtg = L.metric_objects(d_s, gam)
    B = L.box(coords, ginv, sqrtg, Pi)
    R = L.ricci_scalar(coords, g, ginv)
    p = sp.Symbol("p", positive=True)
    d1, d2 = sp.symbols(f"d0:{d_s}"), sp.symbols(f"e0:{d_s}")
    sub = {sp.Derivative(Pi, x, 2): d2[i] for i, x in enumerate(xs)}
    sub.update({sp.Derivative(Pi, x): d1[i] for i, x in enumerate(xs)})

    def split(E):
        E = sp.simplify(E.subs(sub).subs(Pi, p))
        A = sp.simplify(sp.diff(E, d2[0]) / p ** (2 * gam))
        Bq = sp.simplify(sp.diff(E, d1[0], 2) / 2 / p ** (2 * gam))
        rest = sp.simplify(E - p ** (2 * gam) * (A * sum(d2) + Bq * sum(d ** 2 for d in d1)))
        return A, Bq, rest

    a, b, rb = split(B)
    c, e, rr = split(R)
    # checks: C5 at gamma = 1 and lambda_1(gamma) in the weak field
    u = sp.Symbol("u")
    lap, grad2 = sp.symbols("lap grad2")
    full = p ** (2 * gam) * ((a - p * c / 2) * lap + (b - p * e / 2) * grad2)   # (box - R/2) Pi
    c5 = sp.simplify(full.subs(gam, 1) - 2 * p * grad2) == 0 if d_s == 3 else True
    lam1 = sp.simplify((a - p * c / 2).subs(p, 1))
    lam_ok = sp.simplify(lam1 - (2 - gam * (d_s - 1))) == 0
    return {"a": a, "b": b, "c": c, "e": e, "rest_zero": rb == 0 and rr == 0, "c5": bool(c5), "lambda1": lam1,
            "lambda1_ok": bool(lam_ok), "gamma": gam, "p": p, "text": f"box_g Pi = Pi^(2g)[({a}) lap + ({b}) grad^2]; "
                                                                       f"R = Pi^(2g)[({c}) lap + ({e}) grad^2]"}


def candidates(ks):
    """(name, f, f', f'') for the normalised candidates and the negative control."""
    out = [("ln Pi", np.log, lambda P: 1.0 / P, lambda P: -1.0 / P ** 2),
           ("Pi - 1", lambda P: P - 1.0, lambda P: np.ones_like(P), lambda P: np.zeros_like(P))]
    for k in ks:
        out.append((f"(Pi^{k:g} - 1)/{k:g}", (lambda P, k=k: (P ** k - 1.0) / k), (lambda P, k=k: P ** (k - 1.0)),
                    (lambda P, k=k: (k - 1.0) * P ** (k - 2.0))))
    return out


def PQ_functions(forms, gamma: float):
    """Numeric a(p), b(p), c(p), e(p) at the given gamma."""
    import sympy as sp
    fn = {k: sp.lambdify(forms["p"], forms[k].subs(forms["gamma"], gamma), "numpy") for k in ("a", "b", "c", "e")}
    return {k: (lambda P, f=f: np.broadcast_to(np.asarray(f(P), float), np.shape(P)).astype(float)) for k, f in fn.items()}


# ------------------------------------------------------------------------------------------------- the radial problem

def source(p, q, sigma):
    return lambda r: q * np.exp(-np.asarray(r, float) ** 2 / (2.0 * sigma ** 2))


def q_of(u0, sigma, kappa):
    """Source amplitude whose LINEAR solution has central value u0."""
    return -u0 * 4.0 * math.pi * sigma / (kappa * math.sqrt(2.0 / math.pi) * (2.0 * math.pi) ** 1.5 * sigma ** 3)


def linear_solution(r, q, sigma, kappa):
    M = q * (2.0 * math.pi) ** 1.5 * sigma ** 3
    r = np.maximum(np.asarray(r, float), 1e-300)
    return -kappa * M * erf(r / (math.sqrt(2.0) * sigma)) / (4.0 * math.pi * r)


def solve_radial(fs, fun, g, rho, R, rmin, tol, max_nodes, guess):
    """Pi^(2g)[P lap Pi + Q Pi'^2] = kappa rho on [rmin, R]; returns the scipy result."""
    name, f, f1, f2 = fun

    def P_Q(Pi):
        a, b, c, e = fs["a"](Pi), fs["b"](Pi), fs["c"](Pi), fs["e"](Pi)
        P = f1(Pi) * a - 0.5 * f(Pi) * c
        Q = f1(Pi) * b + f2(Pi) * Pi ** (2.0 * g) / Pi ** (2.0 * g) - 0.5 * f(Pi) * e   # (grad Pi)^2 = Pi^2g |grad|^2
        return P, Q

    def ode(r, y):
        Pi, dPi = np.maximum(y[0], 1e-12), y[1]
        P, Q = P_Q(Pi)
        d2 = (rho(r) / Pi ** (2.0 * g) - Q * dPi ** 2) / P - 2.0 * dPi / r
        return np.vstack([dPi, d2])

    def bc(ya, yb):
        return np.array([ya[1], yb[0] - 1.0 + R * yb[1]])

    res = solve_bvp(ode, bc, guess[0], guess[1], tol=tol, max_nodes=max_nodes)
    return res, P_Q


def linear_ode_solve(rho, R, rmin, tol, max_nodes, r0):
    """The linear problem lap u = rho (validation of the solver)."""
    def ode(r, y):
        return np.vstack([y[1], rho(r) - 2.0 * y[1] / r])

    def bc(ya, yb):
        return np.array([ya[1], yb[0] + R * yb[1]])
    return solve_bvp(ode, bc, r0, np.zeros((2, r0.size)), tol=tol, max_nodes=max_nodes)


def grid(rmin, R, sigma, n=4000):
    return np.unique(np.r_[np.linspace(rmin, 6.0 * sigma, n // 2), np.geomspace(6.0 * sigma, R, n // 2)])


def continuation(fs, fun, p, sigma, strengths):
    """Solve from the weak field to each strength in turn; per strength the regularity data."""
    g, kappa = float(p["gamma"]), float(p["kappa"])
    R, rmin = float(p["R_over_sigma"]) * sigma, float(p["rmin_over_sigma"]) * sigma
    r = grid(rmin, R, sigma)
    q0 = q_of(float(p["weak_u0"]), sigma, kappa)
    y = np.vstack([1.0 + linear_solution(r, q0, sigma, kappa), np.gradient(linear_solution(r, q0, sigma, kappa), r)])
    out, x = [], r
    path = [float(p["weak_u0"])] + list(strengths)
    for u0 in path:
        # finer steps in between keep the continuation on the branch connected to the weak field
        last = out[-1]["u0"] if out else float(p["weak_u0"])
        steps = np.linspace(last, u0, 6)[1:] if out else [u0]
        ok = True
        for us in steps:
            q = q_of(us, sigma, kappa)
            rho = source(p, kappa * q, sigma)
            res, P_Q = solve_radial(fs, fun, g, rho, R, rmin, float(p["tol"]), int(p["max_nodes"]), (x, y))
            if not res.success:
                ok = False
                break
            x, y = res.x, res.y
        Pi = y[0]
        P, _ = P_Q(np.maximum(Pi, 1e-12))
        lin = linear_solution(x, q_of(u0, sigma, kappa), sigma, kappa)
        rec = {"u0": u0, "success": bool(ok), "Pi0": float(Pi[0]), "Pi_min": float(Pi.min()), "Pi_R": float(Pi[-1]),
               "monotone": bool(np.all(np.diff(Pi) >= -1e-10)), "P_min": float(np.min(P)), "x": x, "Pi": Pi,
               "charge": float(x[-1] ** 2 * y[1][-1]), "charge_lin": float(kappa * q_of(u0, sigma, kappa) * (2 * math.pi) ** 1.5
                                                                              * sigma ** 3 / (4 * math.pi)),
               "dev_lin": float(np.max(np.abs((Pi - 1.0) - np.interp(x, x, lin))) / max(np.max(np.abs(lin)), 1e-300))}
        rec["regular"] = rec["success"] and rec["Pi_min"] > 0 and rec["monotone"] and rec["P_min"] > 0
        out.append(rec)
        if not rec["regular"]:
            break
    return out


def shoot(fs, fun, g, rho, rho0, R, rmin, Pi0, rtol):
    """Integrate from the centre with Pi(0) = Pi0; returns (mismatch of the exterior condition, profile) or (None, reason)."""
    from scipy.integrate import solve_ivp
    name, f, f1, f2 = fun

    def PQ(P_):
        a, b, c, e = fs["a"](P_), fs["b"](P_), fs["c"](P_), fs["e"](P_)
        return f1(P_) * a - 0.5 * f(P_) * c, f1(P_) * b + f2(P_) - 0.5 * f(P_) * e

    P0, _ = PQ(np.array(Pi0))
    if not P0 > 0:
        return None, "P <= 0 at the centre"
    d2c = rho0 / (3.0 * Pi0 ** (2.0 * g) * float(P0))

    def ode(r, y):
        Pi = y[0]
        P, Q = PQ(np.array(Pi))
        return [y[1], (rho(r) / Pi ** (2.0 * g) - float(Q) * y[1] ** 2) / float(P) - 2.0 * y[1] / r]

    def ev_zero(r, y):
        return y[0] - 1e-6
    ev_zero.terminal = True

    def ev_sing(r, y):
        return float(PQ(np.array(max(y[0], 1e-9)))[0]) - 1e-9
    ev_sing.terminal = True
    sol = solve_ivp(ode, (rmin, R), [Pi0 + 0.5 * d2c * rmin ** 2, d2c * rmin], method="LSODA", rtol=rtol, atol=1e-13,
                    events=(ev_zero, ev_sing), dense_output=False)
    if sol.status != 0 or sol.t[-1] < R * (1 - 1e-12):
        return None, "Pi reached 0 or P reached 0 before R" if sol.status == 1 else "integration failed"
    Pi, dPi = sol.y[0], sol.y[1]
    return float(Pi[-1] - 1.0 + R * dPi[-1]), {"r": sol.t, "Pi": Pi, "dPi": dPi,
                                                "P_min": float(np.min(PQ(Pi)[0])), "monotone": bool(np.all(dPi >= -1e-12))}


def shoot_strength(fs, fun, p, sigma, u0):
    """All regular solutions at one strength (roots of the exterior condition over Pi(0) in (0, 1))."""
    from scipy.optimize import brentq
    g, kappa = float(p["gamma"]), float(p["kappa"])
    R, rmin = float(p["R_over_sigma"]) * sigma, float(p["rmin_over_sigma"]) * sigma
    q = q_of(u0, sigma, kappa)
    rho = source(p, kappa * q, sigma)
    rho0 = kappa * q
    rt = float(p["shoot_rtol"])
    grid0 = 1.0 - np.geomspace(1e-6, 0.999, int(p["n_shoot"]))            # Pi(0) from just below 1 down to 0.001
    vals = []
    for P0 in grid0:
        F, prof = shoot(fs, fun, g, rho, rho0, R, rmin, float(P0), rt)
        vals.append((float(P0), F))
    sols = []
    for (a, Fa), (b, Fb) in zip(vals[:-1], vals[1:]):
        if Fa is None or Fb is None or Fa * Fb > 0:
            continue
        root = brentq(lambda x: shoot(fs, fun, g, rho, rho0, R, rmin, x, rt)[0], b, a, xtol=1e-12, rtol=1e-12)
        F, prof = shoot(fs, fun, g, rho, rho0, R, rmin, root, rt)
        if prof is None or isinstance(prof, str):
            continue
        reg = prof["P_min"] > 0 and prof["monotone"] and float(np.min(prof["Pi"])) > 0
        sols.append({"Pi0": root, "regular": bool(reg), "Pi_R": float(prof["Pi"][-1]),
                     "charge": float(R ** 2 * prof["dPi"][-1]),
                     "charge_lin": float(kappa * q * (2 * math.pi) ** 1.5 * sigma ** 3 / (4 * math.pi)),
                     "P_min": prof["P_min"], "monotone": prof["monotone"]})
    sols.sort(key=lambda s: -s["Pi0"])
    return sols


# ------------------------------------------------------------------------------------------------- run

def run(p, ctx):
    t0 = time.time()
    mode = str(p["mode"])
    g, d_s, kappa, sigma = float(p["gamma"]), int(p["d_s"]), float(p["kappa"]), float(p["sigma"])
    forms = local_forms(d_s)
    fs = PQ_functions(forms, g)
    cands = candidates([float(k) for k in p["ks"]])
    lines = [f"3+1 (d_s = {d_s}), metric -Pi^2 dt^2 + Pi^(-2 gamma) dx^2, gamma = {g:g}; {forms['text']}",
             f"  weak-field coefficient of (box_g - R/2) Pi: lambda_1(gamma) = {forms['lambda1']}"]
    if mode == "V":
        ok_sym = forms["rest_zero"] and forms["c5"] and forms["lambda1_ok"]
        R, rmin = float(p["R_over_sigma"]) * sigma, float(p["rmin_over_sigma"]) * sigma
        q = q_of(-0.3, sigma, kappa)
        rho = source(p, kappa * q, sigma)
        res = linear_ode_solve(rho, R, rmin, float(p["tol"]) * 1e-2, int(p["max_nodes"]), grid(rmin, R, sigma))
        an = linear_solution(res.x, q, sigma, kappa)
        # the analytic solution obeys u(R) = -kappa M erf/(4 pi R), the exterior condition is exact up to erf - 1
        dev = float(np.max(np.abs(res.y[0] - an)) / np.max(np.abs(an)))
        ok_lin = bool(res.success) and dev <= float(p["lin_tol"])
        ln = cands[0]
        sc = []
        for s in p["scale_sigmas"]:
            sp_ = dict(p, sigma=s)
            rec = continuation(fs, ln, sp_, float(s), [float(p["scale_u0"])])[-1]
            sc.append((float(s), rec))
        xs = np.linspace(0.0, 10.0, 200)
        prof = [np.interp(xs * s, rec["x"], rec["Pi"]) for s, rec in sc]
        dsc = float(np.max(np.abs(prof[0] - prof[1])))
        ok_sc = all(rec["regular"] for _, rec in sc) and dsc <= float(p["scale_tol"])
        lines += [f"V (i) symbolic: forms complete {forms['rest_zero']}; C5 reproduced {forms['c5']}; lambda_1(gamma) = 2 - "
                  f"gamma (d_s - 1) {forms['lambda1_ok']}",
                  f"V (ii) linear solver against the analytic solution: max relative deviation {dev:.2e} (allowed {p['lin_tol']})",
                  f"V (iii) scaling: ln Pi at u_lin(0) = {p['scale_u0']} for sigma = {p['scale_sigmas']}: max |Pi_a(r/sigma) - "
                  f"Pi_b(r/sigma)| = {dsc:.2e} (allowed {p['scale_tol']}); Pi(0) = {[round(rec['Pi0'], 6) for _, rec in sc]}"]
        status = "PASSED" if (ok_sym and ok_lin and ok_sc) else "FAILED"
        return Result(status, lines, {"symbolic_ok": ok_sym, "lin_dev": dev, "scale_dev": dsc})
    if mode == "W":
        rows, ok = [], True
        for fun in cands:
            rec = continuation(fs, fun, p, sigma, [])[0]
            okc = rec["regular"] and rec["dev_lin"] <= float(p["weak_tol"])
            ok &= okc
            rows.append((fun[0], rec))
            lines.append(f"W {fun[0]:18s}: u_lin(0) = {p['weak_u0']:g}: Pi(0) - 1 = {rec['Pi0'] - 1:+.4e}; max |u - u_lin| / "
                         f"max |u_lin| = {rec['dev_lin']:.2e} (allowed {p['weak_tol']}) -> {'OK' if okc else 'FAILS'}")
        # negative control f = Pi: the linear coefficient Pi^2g P at Pi = 1
        a1, c1 = float(fs["a"](np.array(1.0))), float(fs["c"](np.array(1.0)))
        coef = 1.0 * a1 - 0.5 * 1.0 * c1
        okn = abs(coef) < 1e-12
        lines.append(f"W negative control f = Pi (equation (a)): linear coefficient {coef:+.3e} -> "
                     f"{'no weak-field response, as derived' if okn else 'UNEXPECTED response'}")
        status = "PASSED" if (ok and okn) else "FAILED"
        return Result(status, lines, {"rows": [(n, {k: v for k, v in r_.items() if k not in ('x', 'Pi')}) for n, r_ in rows],
                                      "control_coefficient": coef})
    # part N: existence by shooting
    strengths = [float(s) for s in p["strengths"]]
    table = {}
    for fun in cands:
        rows_ = []
        for u0 in strengths:
            sols = shoot_strength(fs, fun, p, sigma, u0)
            reg = [s_ for s_ in sols if s_["regular"]]
            rows_.append({"u0": u0, "n_regular": len(reg), "n_roots": len(sols), "branch": reg[0] if reg else None,
                          "all": [round(s_["Pi0"], 6) for s_ in sols]})
        table[fun[0]] = rows_
        last_ok = max((r_["u0"] for r_ in rows_ if r_["n_regular"]), key=abs, default=None)
        first_bad = next((r_["u0"] for r_ in rows_ if not r_["n_regular"]), None)
        lines.append(f"N {fun[0]:18s}: " + ("regular solution at every declared strength" if first_bad is None else
                     f"regular up to u_lin(0) = {last_ok if last_ok is not None else 'none'}; none from {first_bad:g}"))
        for r_ in rows_:
            b = r_["branch"]
            lines.append(f"     u_lin(0) {r_['u0']:+.2f}: regular solutions {r_['n_regular']} (roots {r_['all']})"
                         + (f"; weak-field branch Pi(0) {b['Pi0']:.5f}, Pi(R) {b['Pi_R']:.6f}, min P {b['P_min']:+.3f}, "
                            f"exterior charge / linear {b['charge'] / b['charge_lin']:.4f}" if b else ""))
    names = list(table)
    for i, s in enumerate(strengths):
        vals = {n: (table[n][i]["branch"]["Pi0"] if table[n][i]["branch"] else float("nan")) for n in names}
        lines.append(f"  Pi(0) of the weak-field branch at u_lin(0) {s:+.2f}: " + "; ".join(f"{n} {v:.5f}" for n, v in vals.items()))
    lines.append(f"  (check 3, psi and gamma, not tested here: with the H1 ansatz it is input; {time.time() - t0:.0f} s)")
    return Result("RECORDED", lines, {n: [{k: v for k, v in r_.items()} for r_ in rows_] for n, rows_ in table.items()})
