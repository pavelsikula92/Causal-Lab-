"""
F6.8f-dev, step 1b at finite q = k sigma, b -> 0 (development, decides nothing)
==============================================================================

docs/derivations/f68f_dev.md, step 1b. Scaling limit kappa -> infinity at fixed eps, p (b -> 0) but finite q = k sigma.
With psi = k delta (in-plane transverse phase offset) and eta = k Y (out of the plane of k), the interval is

    { (s, psi, eta) :  Phi(s, psi) + eta^2 / (2 s (1 - s)) <= q^2 / 2 },
    Phi = G1(s, psi) + G2(s, psi) - S_min,

G1 (G2) the first-arrival action from x to (s, psi) (from (s, psi) to y), S the action of step 1b. The eta direction is
exactly quadratic (k has no component there) and is integrated in closed form:

    A(q) = 32 sqrt(2) / (pi q^4)  int_0^1 ds sqrt(s (1 - s))  int dpsi [q^2/2 - Phi]_+^{3/2}      (A(0) = A_PW).

The psi integral by expansion: around each ray b, Phi = S_b - S_min + K x^2/2 + c3 x^3 + c4 x^4 + ..., x = psi - psi_b,
    with K, c3, c4 from the variational equations of order 1-3 along the ray (Hamilton-Jacobi derivatives of G1, G2).
    Then  int dx [c - Phi_b]^{3/2} = (3 pi / 8) c^2 sqrt(2/K) [1 + c (5 c3^2 / (2 K^3) - c4 / K^2) + O(c^2)],
    c = q^2/2 - (S_b - S_min). Every ray with S_b - S_min < q^2/2 contributes its own needle (it is in the past of y
    at that sigma). A ray with K <= 0 somewhere inside (a saddle of the action) is flagged, not integrated.
The expansion is not valid near caustics (column corr); there the exact computation is diag.f68f_1b_exact.
(The 16 Nov version also had an exact check by boundary-value problems, output diag/out/f68f_1b_q_check.txt; it
stopped on a fixed search window and was removed; its two completed rows stay as recorded.)

Run from the repository root (one job = one point (eps, p), cached in diag/out/f68f_1b_q_cache/):
    python -m diag.f68f_1b_q --workers 8
    python -m diag.f68f_1b_q --eps 4 --p 2
Output diag/out/f68f_1b_q_<UTC time>.json and .txt. The run of 16 Nov before this form is diag/out/f68f_1b_q.txt.
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

QS = [0.5, 0.25, 0.125]


def variational(v0, eps, p, th, N, backward=False):
    """Ray from x (slope v0 at s = 0) or, backward, from y (slope v0 at s = 1); with the derivatives of (psi, psi')
    with respect to the slope up to third order. Returns arrays over the N + 1 nodes s = 0 ... 1 (always in
    increasing s): psi, P = psi', S (action from the start), and E, F, H (d psi/dv, d2, d3), E', F', H'."""
    v0 = np.atleast_1d(np.asarray(v0, float))
    th = np.broadcast_to(np.asarray(th, float), v0.shape)
    y = np.zeros((9,) + v0.shape)        # psi, psi', E, E', F, F', H, H', S
    y[1] = v0
    y[3] = 1.0
    h = (-1.0 if backward else 1.0) / N

    def f(s, y):
        ph = p * s + th + y[0]
        sn, cs = np.sin(ph), np.cos(ph)
        f0, f1, f2, f3 = -eps * sn, -eps * cs, eps * sn, eps * cs
        E, F, H = y[2], y[4], y[6]
        return np.stack([y[1], f0, y[3], f1 * E, y[5], f1 * F + f2 * E * E, y[7],
                         f1 * H + 3 * f2 * E * F + f3 * E ** 3, 0.5 * y[1] ** 2 + eps * cs])

    out = np.empty((N + 1,) + y.shape)
    s = 1.0 if backward else 0.0
    out[0] = y
    for i in range(N):
        k1 = f(s, y)
        k2 = f(s + h / 2, y + h / 2 * k1)
        k3 = f(s + h / 2, y + h / 2 * k2)
        k4 = f(s + h, y + h * k3)
        y = y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        s += h
        out[i + 1] = y
    if backward:
        out = out[::-1]
        out[:, 8] = -out[:, 8]           # action from s to 1 (integrated with ds < 0)
    return out


def hj_derivs(out, sign):
    """d^n G / dpsi^n, n = 2, 3, 4, of the action G of the family, from the slope derivatives; sign = +1 from x,
    -1 from y (dG2/dpsi = -psi')."""
    E, dE, F, dF, H, dH = out[:, 2], out[:, 3], out[:, 4], out[:, 5], out[:, 6], out[:, 7]
    psv, psvv, psvvv = E, F, H
    Pv, Pvv, Pvvv = sign * dE, sign * dF, sign * dH
    with np.errstate(divide="ignore", invalid="ignore"):
        G2 = Pv / psv
        Nn = Pvv * psv - Pv * psvv
        G3 = Nn / psv ** 3
        Nv = Pvvv * psv - Pv * psvvv
        G4 = (Nv * psv - 3 * Nn * psvv) / psv ** 5
    return G2, G3, G4


def needle_terms(v, eps, p, th, N):
    """For rays v (one per entry of th): S, and per node s: K, c3, c4 (arrays (N + 1, n))."""
    fo = variational(v, eps, p, th, N)
    w = fo[-1, 1]                                        # slope at y
    bo = variational(w, eps, p, th, N, backward=True)
    a2, a3, a4 = hj_derivs(fo, 1.0)
    b2, b3, b4 = hj_derivs(bo, -1.0)
    K = a2 + b2
    c3 = (a3 + b3) / 6.0
    c4 = (a4 + b4) / 24.0
    return fo[-1, 8], K, c3, c4, fo, bo


def A_expansion(K, c3, c4, cfac, q, N):
    """A contribution of one needle: c = cfac * q^2/2 (cfac = 1 for the first arrival)."""
    s = np.linspace(0.0, 1.0, N + 1)[:, None]
    c = cfac * q * q / 2.0
    inner = np.zeros_like(K)
    ok = (s > 0) & (s < 1)
    ok = np.broadcast_to(ok, K.shape)
    with np.errstate(divide="ignore", invalid="ignore"):
        corr = 1.0 + c * (2.5 * c3 ** 2 / K ** 3 - c4 / K ** 2)
        val = np.sqrt(s * (1 - s) / K) * corr * cfac ** 2
    inner[ok] = val[ok]
    wS = np.ones(N + 1)
    wS[1:-1:2] = 4.0
    wS[2:-1:2] = 2.0
    I = (wS[:, None] * inner).sum(0) / (3.0 * N)
    return 6.0 * I, np.nanmax(np.abs(np.where(ok, c * (2.5 * c3 ** 2 / K ** 3 - c4 / K ** 2), 0.0)), 0)


def point_arrays(eps, p, N=4000, n_th=B.N_TH):
    ths = 2 * math.pi * np.arange(n_th) / n_th
    res = {q: np.zeros(len(ths)) for q in QS}
    flags = {"saddle": 0, "extra": {q: 0 for q in QS}, "maxcorr": {q: 0.0 for q in QS}}
    if eps == 0.0:
        for q in QS:
            res[q][:] = 1.0
        return res, flags, np.ones(len(ths))
    roots = B.branches(eps, p, ths, N // 4)
    ti = np.concatenate([np.full(len(r), i) for i, r in enumerate(roots)])
    vv = B.newton(np.concatenate(roots), eps, p, ths[ti], N)
    keep = np.ones(vv.size, bool)
    for i in range(len(ths)):                             # remove duplicates per phase
        idx = np.nonzero(ti == i)[0]
        _, first = np.unique(np.round(vv[idx], 9), return_index=True)
        m = np.zeros(idx.size, bool)
        m[first] = True
        keep[idx] = m
    ti, vv = ti[keep], vv[keep]
    S, K, c3, c4, fo, bo = needle_terms(vv, eps, p, ths[ti], N)
    Smin = np.array([S[ti == i].min() for i in range(len(ths))])
    dS = S - Smin[ti]
    inside = (np.arange(N + 1) > 0) & (np.arange(N + 1) < N)
    Kmin = K[inside].min(0)
    apw = np.zeros(len(ths))
    for q in QS:
        c = q * q / 2.0
        sel = dS < c
        good = sel & (Kmin > 0)
        flags["saddle"] = max(flags["saddle"], int(np.sum(sel & ~(Kmin > 0))))
        cf = (c - dS[good]) / c
        Ab, mc = A_expansion(K[:, good], c3[:, good], c4[:, good], cf, q, N)
        np.add.at(res[q], ti[good], Ab)
        flags["extra"][q] = int(np.sum(good & (dS > 0)))
        flags["maxcorr"][q] = float(mc.max()) if mc.size else 0.0
    m0 = dS == 0
    A0, _ = A_expansion(K[:, m0], c3[:, m0], c4[:, m0], np.ones(int(m0.sum())), 0.0, N)
    apw[ti[m0]] = A0
    return res, flags, apw


def point(eps, p, N=4000, n_th=B.N_TH):
    res, fl, apw = point_arrays(eps, p, N, n_th)
    return {"A": [res[q] for q in QS], "apw": apw, "saddle": fl["saddle"], "extra": [fl["extra"][q] for q in QS],
            "maxcorr": [fl["maxcorr"][q] for q in QS]}


NAME = "f68f_1b_q"
CODE = C.code_md5([os.path.abspath(__file__), B.__file__])


def job(eps, p, n_th):
    return C.cached(NAME, CODE, point, (eps, p, 4000, n_th))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--eps", type=str, default=",".join(f"{e:g}" for e in B.EPS_GRID))
    ap.add_argument("--p", type=str, default=",".join(f"{x:g}" for x in B.P_GRID))
    ap.add_argument("--n-theta", type=int, default=B.N_TH)
    ap.add_argument("--workers", type=int, default=0)
    a = ap.parse_args()
    t0 = time.time()
    epss, ps = C.floats(a.eps), C.floats(a.p)
    jobs = [(e, p, a.n_theta) for p in ps for e in epss]
    res = parallel_map(job, jobs, default_workers(a.workers), print,
                       lambda i, r, d, n: f"[f68f_1b_q] {d}/{n}: eps {jobs[i][0]:g} p {jobs[i][1]:g} ({time.time() - t0:.0f}s)")
    lines = ["F6.8f-dev step 1b at finite q = k sigma, b -> 0: mean over the phases of A(q) (expansion to O(q^2) "
             "per needle, every ray within q^2/2 of the minimum action included); development, decides nothing", "",
             "delta_A = |A(1/8) - A(1/4)| / max(1, |A(1/8)|, |A(1/4)|) (criterion 1e-3); pair = largest |A_j(1/8) - "
             "A_j(1/4)| over phases; extra = needles of other rays included at q = 1/2, 1/4, 1/8; saddle = rays within "
             "q^2/2 with K <= 0 somewhere (not integrated); corr = largest relative q^2 correction at q = 1/4 (the "
             "expansion is not valid where it exceeds ~0.1)", ""]
    data = []
    for p in ps:
        lines.append(f"p = {p:g}")
        lines.append("     eps     A_PW        A(1/2)      A(1/4)      A(1/8)      delta_A    pair       extra     saddle  corr")
        for e in epss:
            r = res[jobs.index((e, p, a.n_theta))]
            A = np.array(r["A"])
            m = A.mean(1)
            dA = abs(m[2] - m[1]) / max(1.0, abs(m[2]), abs(m[1]))
            pair = float(np.max(np.abs(A[2] - A[1])))
            ex = "/".join(str(x) for x in r["extra"])
            lines.append(f"  {e:6.2f}  {np.mean(r['apw']):.7f}  {m[0]:.7f}  {m[1]:.7f}  {m[2]:.7f}  {dA:.1e}    "
                         f"{pair:.1e}    {ex:9s} {r['saddle']:3d}     {r['maxcorr'][1]:.1e}")
            data.append({"eps": e, "p": p, "q": QS, **r, "mean": m, "delta_A": dA, "pair": pair})
        lines.append("")
    lines.append(f"cost {time.time() - t0:.0f} s")
    stem = C.save(NAME, vars(a), CODE, data, lines, t0)
    print("\n".join(lines))
    print("saved", stem + ".json")


if __name__ == "__main__":
    main()
