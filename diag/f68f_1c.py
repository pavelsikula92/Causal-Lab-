"""
F6.8f-dev, step 1c (development, decides nothing): the needle volume as a series in q^2 = (k sigma)^2, and the age
term of the expected BD action at finite eps
==================================================================================================================

docs/derivations/f68f_dev.md, step 1c. Scaling limit of step 1b (kappa = kT -> infinity, b -> 0 at fixed eps = b kappa^2
and p = k_par T). In the range eps < pi^2 there is one ray x -> y for every phase (step 1b), so the needle volume is
analytic in q^2:

    A(eps, p, theta; q) = A_0 + A_1 q^2 + A_2 q^4 + A_3 q^6 + O(q^8),       A_0 = A_PW of step 1b.

The coefficients come from the Taylor series of the first-arrival actions G1(s, psi) and G2(s, psi) around the ray up to
(psi - psi_r)^8: the families of rays from x and from y are integrated as truncated power series in the initial slope
(degree 8, RK4 on the coefficients), reverted and composed at every node s. With Phi = G1 + G2 - S_min = K x^2 / 2 (1 +
d_1 x + ... + d_6 x^6) and the variable y, K y^2 / 2 = Phi, x = y + x_2 y^2 + ... + x_7 y^7 (Lagrange inversion):

    int dx [c - Phi]_+^{3/2} = sum_j (2j + 1) x_{2j+1} (2 / K)^{j + 1/2} B_j c^{j + 2},   B_j = Gamma(j + 1/2) Gamma(5/2) / Gamma(j + 3),
    A_j = (32 sqrt 2 / pi) (1/2)^{j + 2} B_j (2j + 1) int_0^1 ds sqrt(s (1 - s)) x_{2j+1}(s) (2 / K(s))^{j + 1/2}.

The age term. With the age cut at coordinate time T, dE[S]/dT is the integral over the events y of the past of x at
time distance exactly T (first arrival <= T). In the scaling limit only pairs near the cone contribute (the kernel),
r = T, the factors sqrt(-g) and the measure tend to constants, and per direction p

    rho(p) = (T / (2 l2)) int_0^inf ds < g0(s^2 A(s)) >_theta,     A(s) = A(q^2 = k^2 s / l2),

s = l2 sigma^2 the variable of F6.8c. With u = s (A(s) / A_0)^{1/2} the BD moments mu_0 = mu_1 = mu_2 = 0 remove every
term below s^3; the first survivor (mu_3 = 1/2) gives

    rho(p) - rho_flat(p) = (T / (2 l2)) 4 mu_3 (k^2 / l2)^3 < s4hat A_0^{-2} >_theta + O(k^8),
    s4hat = -e4 + 5 e2 e3 - 5 e2^3,  e2 = a1/2, e3 = a2/2 - a1^2/8, e4 = a3/2 - a1 a2/4 + a1^3/16,  a_j = A_j / A_0.

At small eps, < s4hat A_0^{-2} > = C2(p) eps^2, and the bulk density of F6.8e3 at large n must be
    bulk(p) = 4 mu_3 k^5 n^5 C2(p) / l2^4      (no free constant; corrections O(1/n^2) and O(k^2)); the bulk of F6.8e3 is
    the full second derivative in b (d^2/db^2, as K_hat_bb = (d^2 E/db^2) / (2 k^2) in F6.8c), hence 4, not 2.
That is the check of this file against F6.8e3. The finite-eps age term is Phi(eps) = int dp < s4hat A_0^{-2} >_theta
(per unit of the prefactor above); Phi(eps) / (Phi_2 eps^2) = 1 is the second-order extrapolation (the null model).

Run from the repository root (one job = one point (eps, p), cached in diag/out/f68f_1c_cache/):
    python -m diag.f68f_1c --eps 0.01,0.1,0.25,0.5,1,2,4 --workers 2
    python -m diag.f68f_1c --check                  series against the exact needle and against diag.f68f_1b_q
Output diag/out/f68f_1c_<UTC time>.json and .txt.
"""
from __future__ import annotations

import argparse
import math
import os
import time

import numpy as np

from core.stats.parallel import default_workers, parallel_map
from diag import f68f_1b as B
from diag import f68f_common as C

D = 8                      # degree of the truncated power series
N_RK = 2000
L2 = 0.3618006272791338    # l2 of F6.8a
MU3 = 0.5                  # int s^3 g0(s^2) ds

# ----------------------------------------------------------------------------- truncated power series (axis 0)
_M = np.zeros((D + 1, D + 1, D + 1))
for _i in range(D + 1):
    for _j in range(D + 1 - _i):
        _M[_i, _j, _i + _j] = 1.0


def tmul(a, b):
    return np.einsum("i...,j...,ijk->k...", a, b, _M)


def tpow_list(a, n):
    """[a^0, a^1, ..., a^n] for a series a."""
    one = np.zeros_like(a)
    one[0] = 1.0
    out = [one, a]
    for _ in range(n - 1):
        out.append(tmul(out[-1], a))
    return out


def tsincos(x):
    """sin and cos of a series x WITH constant term: sin(x0 + xt), cos(x0 + xt)."""
    x0 = x[0]
    xt = x.copy()
    xt[0] = 0.0
    pw = tpow_list(xt, D)
    S = np.zeros_like(x)
    Cc = np.zeros_like(x)
    f = 1.0
    for m in range(D + 1):
        if m > 0:
            f *= m
        term = pw[m] / f
        if m % 4 == 0:
            Cc += term
        elif m % 4 == 1:
            S += term
        elif m % 4 == 2:
            Cc -= term
        else:
            S -= term
    s0, c0 = np.sin(x0), np.cos(x0)
    return s0 * Cc + c0 * S, c0 * Cc - s0 * S


def compose(f, g):
    """f(g) for series f (in t) and g (in x) with g[0] = 0."""
    pw = tpow_list(g, D)
    out = np.zeros_like(g)
    for m in range(D + 1):
        out += f[m] * pw[m]
    return out


def revert(g):
    """x(y) with g(x(y)) = y, for g[0] = 0, g[1] != 0."""
    a1 = g[1]
    h = g / a1                                  # h(x) = x + ...
    x = np.zeros_like(g)
    x[1] = 1.0
    for _ in range(D):
        # x = y - sum_{n >= 2} h_n x^n
        pw = tpow_list(x, D)
        nx = np.zeros_like(g)
        nx[1] = 1.0
        for n in range(2, D + 1):
            nx -= h[n] * pw[n]
        x = nx
    # x(y) for h; for g: x(y) = x_h(y / a1)
    out = x.copy()
    for n in range(D + 1):
        out[n] = x[n] / a1 ** n
    return out


def family_series(eps, p, th, v0, backward=False):
    """Series in delta (the change of the slope at the start) of psi(s) and of the action S(s) along the rays of the
    family, at all N_RK + 1 nodes of [0, 1] (increasing s). From x: psi(0) = 0, psi'(0) = v0 + delta, S from 0. From y
    (backward): psi(1) = 0, psi'(1) = v0 + delta, S = action from s to 1."""
    n = np.size(v0)
    th = np.broadcast_to(np.asarray(th, float), (n,))
    psi = np.zeros((D + 1, n))
    dpsi = np.zeros((D + 1, n))
    dpsi[0] = v0
    dpsi[1] = 1.0
    S = np.zeros((D + 1, n))
    h = (-1.0 if backward else 1.0) / N_RK

    def f(s, psi, dpsi):
        sn, cs = tsincos(psi + (p * s + th)[None, :] * np.eye(D + 1)[:, :1])
        return dpsi, -eps * sn, 0.5 * tmul(dpsi, dpsi) + eps * cs

    out_psi = np.empty((N_RK + 1, D + 1, n))
    out_S = np.empty((N_RK + 1, D + 1, n))
    out_dpsi = np.empty((N_RK + 1, D + 1, n))
    s = 1.0 if backward else 0.0
    out_psi[0], out_S[0], out_dpsi[0] = psi, S, dpsi
    for i in range(N_RK):
        k1 = f(s, psi, dpsi)
        k2 = f(s + h / 2, psi + h / 2 * k1[0], dpsi + h / 2 * k1[1])
        k3 = f(s + h / 2, psi + h / 2 * k2[0], dpsi + h / 2 * k2[1])
        k4 = f(s + h, psi + h * k3[0], dpsi + h * k3[1])
        psi = psi + h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        dpsi = dpsi + h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        S = S + h / 6 * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2])
        s += h
        out_psi[i + 1], out_S[i + 1], out_dpsi[i + 1] = psi, S, dpsi
    if backward:
        out_psi, out_S, out_dpsi = out_psi[::-1], -out_S[::-1], out_dpsi[::-1]
    return out_psi, out_S, out_dpsi


def rays(eps, p, ths):
    """The first-arrival ray x -> y at every phase: slope at x, slope at y, number of rays."""
    ths = np.asarray(ths, float)
    roots = B.branches(eps, p, ths, 1000)
    v = np.empty(ths.size)
    nr = np.empty(ths.size, int)
    for i, th in enumerate(ths):
        vv = B.newton(roots[i], eps, p, np.full(roots[i].size, th), 4000)
        vv = np.unique(np.round(vv, 10))
        _, _, S = B.end_state(vv, eps, p, np.full(vv.size, th), 4000)
        v[i] = vv[np.argmin(S)]
        nr[i] = vv.size
    out = B.integrate(v, eps, p, ths, 4000)
    return v, out[-1, 1], nr


def series_coefficients(eps, p, ths):
    """A_0 ... A_3 for every phase in ths, and diagnostics (max |linear term of Phi|, min K, number of rays)."""
    ths = np.atleast_1d(np.asarray(ths, float))
    v, w, nray = rays(eps, p, ths)
    fp, fS, _ = family_series(eps, p, ths, v)
    bp, bS, _ = family_series(eps, p, ths, w, backward=True)
    idx = np.arange(1, N_RK)                         # interior nodes
    nt = ths.size
    s = np.repeat(idx / N_RK, nt)
    flat = lambda X: X[idx].transpose(1, 0, 2).reshape(D + 1, -1).copy()     # (D + 1, nodes * phases)
    Pf, Sf, Pb, Sb = flat(fp), flat(fS), flat(bp), flat(bS)
    gf = Pf.copy()
    gf[0] = 0.0
    gb = Pb.copy()
    gb[0] = 0.0
    G1 = compose(Sf, revert(gf))
    G2 = compose(Sb, revert(gb))
    Phi = G1 + G2
    lin = float(np.max(np.abs(Phi[1])))
    Phi[0] = 0.0
    Phi[1] = 0.0
    K = 2.0 * Phi[2]
    dser = np.zeros_like(Phi)
    dser[0] = 1.0
    dser[1:D - 1] = Phi[3:D + 1] / Phi[2]
    t = dser.copy()
    t[0] = 0.0
    pw = tpow_list(t, D)
    sq = np.zeros_like(t)
    coef = 1.0
    for m in range(D + 1):
        sq += coef * pw[m]
        coef *= (0.5 - m) / (m + 1)
    yx = np.zeros_like(t)
    yx[1:] = sq[:-1]                                 # y = x * sqrt(1 + sum d_m x^m)
    xy = revert(yx)
    w_ = np.ones(N_RK - 1)
    w_[0::2] = 4.0
    w_[1::2] = 2.0                                   # Simpson over the interior nodes (end values are 0)
    w_ = np.repeat(w_ / (3.0 * N_RK), nt)
    A = np.empty((nt, 4))
    for j in range(4):
        Bj = math.gamma(j + 0.5) * math.gamma(2.5) / math.gamma(j + 3.0)
        integrand = np.sqrt(s * (1 - s)) * (2 * j + 1) * xy[2 * j + 1] * (2.0 / K) ** (j + 0.5)
        A[:, j] = (32 * math.sqrt(2) / math.pi) * 0.5 ** (j + 2) * Bj * (w_ * integrand).reshape(-1, nt).sum(0)
    return {"A": A, "lin": lin, "minK": float(K.min()), "nray": int(nray.max())}


def s4hat(A):
    a1, a2, a3 = A[1] / A[0], A[2] / A[0], A[3] / A[0]
    e2 = a1 / 2
    e3 = a2 / 2 - a1 ** 2 / 8
    e4 = a3 / 2 - a1 * a2 / 4 + a1 ** 3 / 16
    return -e4 + 5 * e2 * e3 - 5 * e2 ** 3


def point(eps, p, n_th):
    ths = 2 * math.pi * np.arange(n_th) / n_th
    r = series_coefficients(eps, p, ths)
    A = r["A"]
    sh = np.array([s4hat(a) for a in A])
    return {"A": A, "s4hat": sh, "phi": float(np.mean(sh * A[:, 0] ** -2.0)), "A0": float(A[:, 0].mean()),
            "lin": r["lin"], "minK": r["minK"], "nray": r["nray"]}


NAME = "f68f_1c"
CODE = C.code_md5([os.path.abspath(__file__), B.__file__])


def job(eps, p, n_th):
    return C.cached(NAME, CODE, point, (eps, p, n_th))


def check():
    """The series against the exact needle (diag.f68f_1b_exact) and against the O(q^2) expansion (diag.f68f_1b_q)."""
    from diag import f68f_1b_exact as X
    from diag import f68f_1b_q as Q
    lines = ["check of the q^2 series against the exact needle and the O(q^2) expansion (single phases)", ""]
    for (eps, p, th) in [(4.0, 2.0, 2.0), (4.0, 0.0, 0.7), (2.0, 1.0, 4.0), (1.0, 4.0, 1.0)]:
        sc = series_coefficients(eps, p, [th])
        A = sc["A"][0]
        try:
            ex = X.phase_value(eps, p, th)["A"]
        except RuntimeError as err:                    # the exact needle's window search failed (instrument limit)
            ex = [float("nan")] * len(X.QS)
            lines.append(f"  exact needle failed at eps {eps} p {p} theta {th}: {err}")
        v = float(rays(eps, p, [th])[0][0])
        Apw, _, _, _ = B.volume_A(np.array([v]), eps, p, np.array([th]), 4000)
        S_, K_, c3_, c4_, _, _ = Q.needle_terms(np.array([v]), eps, p, np.array([th]), 4000)
        a_q = [float(Q.A_expansion(K_, c3_, c4_, np.ones(1), q, 4000)[0][0]) for q in X.QS]
        A1q = (a_q[2] - float(Apw[0])) / X.QS[2] ** 2
        lines.append(f"eps {eps} p {p} theta {th}: A0 {A[0]:.12f} (step 1b {float(Apw[0]):.12f}); A1 {A[1]:+.10e} "
                     f"(O(q^2) expansion {A1q:+.10e}); A2 {A[2]:+.6e}; A3 {A[3]:+.6e}; max linear term {sc['lin']:.1e}")
        for q, e in zip(X.QS, ex):
            ser = A[0] + A[1] * q ** 2 + A[2] * q ** 4 + A[3] * q ** 6
            lines.append(f"    q {q:.3f}: exact {e:.10f}  series {ser:.10f}  difference {e - ser:+.2e}  (q^8 = {q ** 8:.1e})")
    return lines


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--eps", type=str, default="0.01,0.1,0.25,0.5,1,2,4")
    ap.add_argument("--p", type=str, default="")
    ap.add_argument("--n-theta", type=int, default=128)
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--workers", type=int, default=0)
    a = ap.parse_args()
    t0 = time.time()
    if a.check:
        lines = check()
        lines.append(f"cost {time.time() - t0:.0f} s")
        stem = C.save(NAME + "_check", vars(a), CODE, {}, lines, t0)
        print("\n".join(lines))
        print("saved", stem + ".txt")
        return
    from phases.F6_order_dynamics.F6_08e3_bd_action_optical_asymptotic_rev2 import _p_nodes
    pn, wp = _p_nodes([[0.0, 12.0, 10], [12.0, 24.0, 6]])
    ps = C.floats(a.p) if a.p else [float(x) for x in pn]
    epss = C.floats(a.eps)
    jobs = [(e, p, a.n_theta) for e in epss for p in ps]
    res = parallel_map(job, jobs, default_workers(a.workers), print,
                       lambda i, r, d, n: f"[f68f_1c] {d}/{n}: eps {jobs[i][0]:g} p {jobs[i][1]:.3f} ({time.time() - t0:.0f}s)")
    lines = ["F6.8f-dev step 1c: needle volume as a series in q^2 and the age term at finite eps (scaling limit); "
             "development, decides nothing", "",
             "per eps: Phi(eps) = int dp <s4hat A0^-2>_theta (p panels of F6.8e3, [0, 24]); N(eps) = Phi / (eps^2 C2), "
             "C2 = Phi / eps^2 at the smallest eps (N = 1 is the second-order extrapolation)", ""]
    data = []
    phis = {}
    for e in epss:
        R = [r for j, r in zip(jobs, res) if j[0] == e]
        prof = np.array([r["phi"] for r in R])
        Phi = float((wp * prof).sum()) if not a.p else float("nan")
        phis[e] = Phi
        data.append({"eps": e, "p": ps, "profile": prof, "Phi": Phi, "A0": [r["A0"] for r in R],
                     "lin": max(r["lin"] for r in R), "minK": min(r["minK"] for r in R),
                     "nray": max(r["nray"] for r in R)})
    e0 = min(epss)
    C2 = phis[e0] / e0 ** 2
    lines.append("     eps      Phi(eps)          N(eps)      rays  min K      max |linear|")
    for d in data:
        lines.append(f"  {d['eps']:6.3f}  {d['Phi']:+.10e}  {d['Phi'] / (C2 * d['eps'] ** 2):+.6f}   {d['nray']}    "
                     f"{d['minK']:.3e}  {d['lin']:.1e}")
    lines.append("")
    lines.append("profile <s4hat A0^-2>_theta / eps^2 per p:")
    for d in data:
        lines.append(f"  eps {d['eps']:g}: " + " ".join(f"{v / d['eps'] ** 2:+.4e}" for v in d["profile"]))
    lines.append(f"cost {time.time() - t0:.0f} s")
    stem = C.save(NAME, vars(a), CODE, {"points": data, "C2": C2}, lines, t0)
    print("\n".join(lines))
    print("saved", stem + ".json")


if __name__ == "__main__":
    main()
