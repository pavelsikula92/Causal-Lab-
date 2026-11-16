"""
F6.8f-dev, step 2: the assembly of K_eff from the age term Phi(eps) (development, decides nothing)
===================================================================================================

docs/derivations/f68f_dev.md, step 2 (the author's order of 22:49, point 2; the goal of 09:59). In the scaling limit
(step 1c) the optical age term per unit age n = kT is dE/dn = c Phi(b n^2) with c fixed by F6.8e3 (c C2 = the b^2
coefficient). At fixed k and age n:

    E(b, n)  = c int_0^n Phi(b m^2) dm,
    E_opt    = dE/db = c int_0^n m^2 Phi'(b m^2) dm = (c / 2 b^{3/2}) M(eps),   M(eps) = int_0^eps sqrt(u) Phi'(u) du,
    K_eff    = (E_opt(b) - E_opt(0)) / b = (c / 2 b^{5/2}) M(eps)              (E_opt(0) = 0: Phi is even in eps),
    K_2      = lim_{b -> 0} K_eff = (2/5) c C2 n^5                             (the second order, F6.8e3: (kT)^5),
    R(eps)   = K_eff / K_2 = M(eps) / ((4/5) C2 eps^{5/2}).

R = 1 is the null model (the second-order extrapolation). At fixed b the age dependence is K_eff(T) = K_2(T) R(b k^2 T^2),
with the local exponent gamma_T = d ln|K_eff| / d ln T = 2 d ln|M| / d ln eps (5 for the second order; a static limit
needs gamma_T -> 0). M is computed by parts without Phi':

    M(eps) = sqrt(eps) Phi(eps) - (1/2) int_0^eps u^{-1/2} Phi(u) du = sqrt(eps) Phi(eps) - int_0^{sqrt eps} Phi(t^2) dt,

the last integral by Gauss-Legendre on panels in t = sqrt(u) (Phi(t^2) is smooth and even in t). Phi(eps) =
int dp <s4hat A0^-2>_theta from diag.f68f_1c (its per-job cache is reused), p by Gauss-Legendre on panels as F6.8e3.

Checks of the author's list (09:59): 1 b -> 0 (R and N at the smallest eps; the F6.8e3 comparison is diag.f68f_1c_e3);
2 zero age (M ~ eps^{5/2} -> 0; N at the smallest quadrature node); 3 k = 0 (eps = 0: Phi = 0 at every p, computed);
4 both signs of b (exact: n(-b, theta) = n(b, theta + pi), and an even number of phases maps the theta grid onto itself,
so Phi(-eps) = Phi(eps) and K_eff(-b) = K_eff(b) identically); 5 conformal and mixed sectors (no code path here);
6 stable instrument (one ray per phase, stationarity of the ray, G against G - 1 nodes per t panel, n_theta against
2 n_theta at the ends). --selftest checks the assembly on Phi = u^2 (R = 1) and Phi = u^3 (M = 6/7 eps^{7/2}).

Run from the repository root (one job = one (eps, p, n_theta), cached in diag/out/f68f_1c_cache/):
    python -m diag.f68f_keff --ends 0.5,3 --g 4 --n-theta 16 --p-panels 0:12:3,12:24:1 --workers 2
    python -m diag.f68f_keff --selftest
Output diag/out/f68f_keff_<UTC time>.json and .txt.
"""
from __future__ import annotations

import argparse
import math
import os
import time

import numpy as np

from core.stats.parallel import default_workers, parallel_map
from diag import f68f_1b as B
from diag import f68f_1c as S
from diag import f68f_common as C

NAME = "f68f_keff"
CODE = C.code_md5([os.path.abspath(__file__), S.__file__, B.__file__])


def gl(a, b, n):
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5 * (b - a) * x + 0.5 * (b + a), 0.5 * (b - a) * w


def p_nodes(spec):
    """'0:12:3,12:24:1' -> Gauss-Legendre nodes and weights on the panels (the rule of F6.8e3)."""
    xs, ws = [], []
    for item in spec.split(","):
        a, b, n = item.split(":")
        x, w = gl(float(a), float(b), int(n))
        xs.append(x)
        ws.append(w)
    return np.concatenate(xs), np.concatenate(ws)


def t_nodes(ends, g):
    """Gauss-Legendre nodes in t on [0, sqrt e1], [sqrt e1, sqrt e2], ...; per panel (nodes, weights)."""
    edges = [0.0] + [math.sqrt(e) for e in ends]
    return [gl(a, b, g) for a, b in zip(edges[:-1], edges[1:])]


def assemble(ends, panels, phi):
    """M(eps) at the ends from Phi (a callable of eps) by parts; cumulative over the t panels."""
    M, I = [], 0.0
    for e, (t, w) in zip(ends, panels):
        I += float(np.sum(w * np.array([phi(x * x) for x in t])))
        M.append(math.sqrt(e) * phi(e) - I)
    return M


def selftest():
    ends = [0.25, 1.0, 2.0, 4.0]
    lines = ["selftest of the assembly (G = 4 nodes per t panel)"]
    for name, f, exact in [("Phi = u^2", lambda u: u * u, lambda e: 0.8 * e ** 2.5),
                           ("Phi = u^3", lambda u: u ** 3, lambda e: 6.0 / 7.0 * e ** 3.5)]:
        M = assemble(ends, t_nodes(ends, 4), f)
        err = max(abs(m / exact(e) - 1) for m, e in zip(M, ends))
        lines.append(f"  {name}: largest relative error of M {err:.1e}")
    return lines


def job(eps, p, n_th):
    r = S.job(eps, p, n_th)
    return {"phi": r["phi"], "nray": r["nray"], "lin": r["lin"]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ends", type=str, default="0.5,3", help="eps at which R is assembled (development points)")
    ap.add_argument("--g", type=int, default=4, help="Gauss nodes per t panel")
    ap.add_argument("--n-theta", type=int, default=16)
    ap.add_argument("--p-panels", type=str, default="0:12:3,12:24:1")
    ap.add_argument("--c2", type=str, default="0.1,0.2")
    ap.add_argument("--no-checks", action="store_true", help="skip G - 1, 2 n_theta and k = 0")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--workers", type=int, default=0)
    a = ap.parse_args()
    t0 = time.time()
    if a.selftest:
        lines = selftest()
        stem = C.save(NAME + "_selftest", vars(a), CODE, {}, lines, t0)
        print("\n".join(lines))
        print("saved", stem + ".json")
        return
    ends, c2e = C.floats(a.ends), C.floats(a.c2)
    pn, wp = p_nodes(a.p_panels)
    nt = a.n_theta
    pan = t_nodes(ends, a.g)
    pan_c = t_nodes(ends, a.g - 1)
    eps_main = sorted(set(c2e + ends + [float(x * x) for t, _ in pan for x in t]))
    jobs = [(e, float(x), nt) for e in eps_main for x in pn]
    if not a.no_checks:
        eps_c = sorted(set(float(x * x) for t, _ in pan_c for x in t) - set(eps_main))
        jobs += [(e, float(x), nt) for e in eps_c for x in pn]
        jobs += [(e, float(x), 2 * nt) for e in ends for x in pn]
        jobs += [(0.0, float(x), nt) for x in pn]
    res = parallel_map(job, jobs, default_workers(a.workers), print,
                       lambda i, r, d, n: f"[f68f_keff] {d}/{n}: eps {jobs[i][0]:.4g} p {jobs[i][1]:.3f} "
                                          f"n_theta {jobs[i][2]} ({time.time() - t0:.0f}s)")
    tab = {(j[0], j[1], j[2]): r for j, r in zip(jobs, res)}

    def Phi(e, n=nt):
        return float(np.sum(wp * np.array([tab[(e, float(x), n)]["phi"] for x in pn])))

    def key(e):                        # node values were stored under the float of x * x
        for k_ in {j[0] for j in jobs}:
            if abs(k_ - e) <= 1e-15 * max(1.0, e):
                return k_
        raise KeyError(e)

    phi = lambda e: Phi(key(e))
    e1, e2 = c2e
    C2 = (4 * Phi(e1) / e1 ** 2 - Phi(e2) / e2 ** 2) / 3 if abs(e2 - 2 * e1) < 1e-12 else Phi(e1) / e1 ** 2
    M = assemble(ends, pan, phi)
    R = [m / (0.8 * C2 * e ** 2.5) for m, e in zip(M, ends)]
    N = [Phi(e) / (C2 * e * e) for e in ends]
    nray = max(r["nray"] for r in res)
    lin = max(r["lin"] for r in res)
    lines = [f"F6.8f-dev step 2: K_eff / K_2 = R(eps) from the age term (scaling limit); development, decides nothing",
             f"  p panels {a.p_panels} ({pn.size} directions), {nt} phases, G = {a.g} nodes per t panel; "
             f"C2 = {C2:+.10e} (from eps {e1:g}, {e2:g})", "",
             "     eps      Phi(eps)           N = Phi/(C2 eps^2)   M(eps)             R = K_eff/K_2"]
    for e, ph, n_, m, r in zip(ends, [Phi(e) for e in ends], N, M, R):
        lines.append(f"  {e:6.3f}  {ph:+.10e}  {n_:+.6f}            {m:+.10e}  {r:+.6f}")
    gam = []
    prev = (None, None)
    for e, m in zip(ends, M):
        if prev[0] is not None and m * prev[1] > 0:
            g_ = 2 * math.log(m / prev[1]) / math.log(e / prev[0])
            gam.append((prev[0], e, g_))
            lines.append(f"  gamma_T = d ln K_eff / d ln T between eps {prev[0]:g} and {e:g}: {g_:+.4f} "
                         f"(second order: 5)")
        prev = (e, m)
    lines += ["", "checks (the author's list of 09:59)",
              f"  1 b -> 0: N(eps {e2:g}) = {Phi(e2) / (C2 * e2 * e2):+.6f}, R at eps {ends[0]:g} = {R[0]:+.6f} "
              "(F6.8e3 itself: diag.f68f_1c_e3)"]
    t_small = float(pan[0][0][0] ** 2)
    lines.append(f"  2 zero age: N at the smallest node eps {t_small:.2e} = {phi(t_small) / (C2 * t_small ** 2):+.6f}; "
                 "M ~ eps^5/2 -> 0 by construction")
    data = {"ends": ends, "Phi": [Phi(e) for e in ends], "N": N, "M": M, "R": R, "C2": C2, "gamma_T": gam,
            "nray": nray, "lin": lin, "p": pn, "wp": wp}
    if not a.no_checks:
        k0 = max(abs(tab[(0.0, float(x), nt)]["phi"]) for x in pn)
        Mc = assemble(ends, pan_c, lambda e: Phi(key(e)))
        Rc = [m / (0.8 * C2 * e ** 2.5) for m, e in zip(Mc, ends)]
        dth = [abs(Phi(e, 2 * nt) / Phi(e) - 1) for e in ends]
        lines += [f"  3 k = 0 (eps = 0): largest |<s4hat A0^-2>_theta| over the directions {k0:.1e}",
                  "  4 both signs of b: exact (n(-b, theta) = n(b, theta + pi); the theta grid is mapped onto itself)",
                  "  5 conformal and mixed sectors: not computed here (no code path)",
                  f"  6 stable instrument: rays per phase {nray} (1 required); stationarity of the ray {lin:.1e}; "
                  f"R with G - 1 = {a.g - 1}: " + ", ".join(f"{abs(x - y):.1e}" for x, y in zip(Rc, R))
                  + f" (|difference|); Phi with {2 * nt} phases against {nt}: "
                  + ", ".join(f"{d:.1e}" for d in dth) + " (relative)"]
        data.update({"k0": k0, "R_g_minus_1": Rc, "dphi_theta": dth})
    lines.append(f"cost {time.time() - t0:.0f} s")
    stem = C.save(NAME, vars(a), CODE, data, lines, t0)
    print("\n".join(lines))
    print("saved", stem + ".json")


if __name__ == "__main__":
    main()
