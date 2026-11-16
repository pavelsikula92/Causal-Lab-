"""
F6.8f-dev, step 1b at finite q, b -> 0: the exact needle with all branches (development, decides nothing)
=======================================================================================================

docs/derivations/f68f_dev.md, step 1b. For the points where the expansion of diag/f68f_1b_q.py is not valid (near the
caustics), A(q) without any expansion:

    A(q) = 32 sqrt(2) / (pi q^4)  int_0^1 ds sqrt(s (1 - s))  int dpsi [q^2/2 - Phi(s, psi)]_+^{3/2},
    Phi = G1 + G2 - S_min,

with G1(s, .) the TRUE first-arrival action from x (lower envelope over the whole family of rays from x) and G2 the
same from y. The families are integrated once per phase on a dense grid of slopes (RK4, with the slope derivatives
d psi/dv, d2 psi/dv2 and the momentum), and at every Gauss node s the envelope is taken over all monotone pieces of
the curve v -> (psi(s; v), S(s; v)), with quintic Hermite interpolation in v (dS/dv = P dpsi/dv). No ray is selected
in advance: every ray, every fold and every cut locus is in the envelope. The psi integral is done on a fine grid
over a window that is widened until Phi > q^2/2 at both ends.

Run from the repository root (one job = one phase; jobs cached in diag/out/f68f_1b_exact_cache/, so an interrupted
run resumes):
    python -m diag.f68f_1b_exact --points 10:1,9.87:1,16:4 --workers 8
    python -m diag.f68f_1b_exact --points 4:2 --n-theta 4
Output diag/out/f68f_1b_exact_<UTC time>.json (per-phase A for every q) and .txt. The runs of 16 Nov before this form
are diag/out/f68f_1b_exact.txt and diag/out/f68f_1b_exact_e*_p*.npy.
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
N_RK = 2000
N_V = 3001
N_S = 24
N_PSI = 4001


def family(eps, p, th, vmax, sn, backward=False):
    """Rays from x (or from y, backward) with slopes v in [-vmax, vmax]; at the nodes sn returns psi, E = dpsi/dv,
    F = d2psi/dv2, P (momentum: +psi' from x, -psi' from y so that dG/dpsi = P) and S (action from the start)."""
    v = np.linspace(-vmax, vmax, N_V)
    y = np.zeros((7, N_V))          # psi, psi', E, E', F, F', S
    y[1] = v
    y[3] = 1.0
    h = (-1.0 if backward else 1.0) / N_RK

    def f(s, y):
        ph = p * s + th + y[0]
        sn_, cs = np.sin(ph), np.cos(ph)
        f1 = -eps * cs
        return np.stack([y[1], -eps * sn_, y[3], f1 * y[2], y[5], f1 * y[4] + eps * sn_ * y[2] ** 2,
                         0.5 * y[1] ** 2 + eps * cs])

    grid = np.linspace(0, 1, N_RK + 1)
    if backward:
        grid = grid[::-1]
    out = np.empty((N_RK + 1, 7, N_V))
    out[0] = y
    s = grid[0]
    for i in range(N_RK):
        k1 = f(s, y)
        k2 = f(s + h / 2, y + h / 2 * k1)
        k3 = f(s + h / 2, y + h / 2 * k2)
        k4 = f(s + h, y + h * k3)
        y = y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        s += h
        out[i + 1] = y
    if backward:
        out = out[::-1]
        out[:, 6] = -out[:, 6]
        out[:, 1] = -out[:, 1]        # momentum of G2: dG2/dpsi = -psi'
    sg = np.linspace(0, 1, N_RK + 1)
    # values at the nodes sn by local quintic interpolation in s would cost; use the RK grid nodes nearest and RK
    # sub-steps instead: integrate exactly to sn by choosing N_RK so that sn are grid points is not possible for
    # Gauss nodes; we interpolate with 7-point Lagrange in s (RK grid spacing 1e-3, error ~ h^7)
    res = {}
    idx = np.clip(np.searchsorted(sg, sn) - 3, 0, N_RK - 6)
    for name, j in (("psi", 0), ("P", 1), ("E", 2), ("F", 4), ("S", 6)):
        vals = np.empty((len(sn), N_V))
        for m, (s0, i0) in enumerate(zip(sn, idx)):
            xs = sg[i0:i0 + 7]
            ys = out[i0:i0 + 7, j]
            w = np.ones(7)
            for a in range(7):
                for b in range(7):
                    if a != b:
                        w[a] *= (s0 - xs[b]) / (xs[a] - xs[b])
            vals[m] = (w[:, None] * ys).sum(0)
        res[name] = vals
    res["v"] = v
    return res


def envelope(psi, S, P, E, F, x):
    """Lower envelope G(x) = min over the curve v -> (psi(v), S(v)) at the points x (1-d arrays over v)."""
    G = np.full(x.shape, np.inf)
    sgn = np.sign(E)
    cuts = np.nonzero(sgn[1:] != sgn[:-1])[0] + 1
    bounds = np.r_[0, cuts, len(psi)]
    for a, b in zip(bounds[:-1], bounds[1:]):
        if b - a < 2:
            continue
        ps = psi[a:b]
        if ps[0] > ps[-1]:
            ps = ps[::-1]
            Ss, Ps, Es, Fs = S[a:b][::-1], P[a:b][::-1], E[a:b][::-1], F[a:b][::-1]
        else:
            Ss, Ps, Es, Fs = S[a:b], P[a:b], E[a:b], F[a:b]
        m = (x >= ps[0]) & (x <= ps[-1])
        if not np.any(m):
            continue
        xm = x[m]
        j = np.clip(np.searchsorted(ps, xm) - 1, 0, len(ps) - 2)
        # Hermite in psi: S as a function of psi on this piece has dS/dpsi = P and d2S/dpsi2 = dP/dpsi = (dP/dv)/E;
        # dP/dv is not stored, so use cubic Hermite with (S, P) at both ends (error ~ dpsi^4 S'''')
        x0, x1 = ps[j], ps[j + 1]
        hh = x1 - x0
        t = (xm - x0) / hh
        h00 = 2 * t ** 3 - 3 * t ** 2 + 1
        h10 = t ** 3 - 2 * t ** 2 + t
        h01 = -2 * t ** 3 + 3 * t ** 2
        h11 = t ** 3 - t ** 2
        val = h00 * Ss[j] + h10 * hh * Ps[j] + h01 * Ss[j + 1] + h11 * hh * Ps[j + 1]
        G[m] = np.minimum(G[m], val)
    return G


def phase_A(eps, p, th, Smin, psi_rays, qs=QS):
    """A(q) for one phase by the exact envelope. psi_rays: the rays x -> y within q_max^2/2 (to place the window)."""
    xs, wsg = np.polynomial.legendre.leggauss(N_S)
    sn = 0.5 * (xs + 1.0)
    ws = 0.5 * wsg
    vmax = eps + 2.0 * math.sqrt(eps) + 4.0 + 2.0 / math.sqrt(sn.min())
    fx = family(eps, p, th, vmax, sn)
    fy = family(eps, p, th, vmax, sn, backward=True)
    cmax = max(qs) ** 2 / 2
    tot = {q: 0.0 for q in qs}
    worst = 0.0
    for m in range(N_S):
        s = sn[m]
        lo = min(r[m] for r in psi_rays)
        hi = max(r[m] for r in psi_rays)
        W = 3.0 * math.sqrt(2 * cmax * s * (1 - s)) + 0.05
        for _ in range(8):
            x = np.linspace(lo - W, hi + W, N_PSI)
            G1 = envelope(fx["psi"][m], fx["S"][m], fx["P"][m], fx["E"][m], fx["F"][m], x)
            G2 = envelope(fy["psi"][m], fy["S"][m], fy["P"][m], fy["E"][m], fy["F"][m], x)
            Phi = G1 + G2 - Smin
            if np.isfinite(Phi[0]) and np.isfinite(Phi[-1]) and Phi[0] > cmax and Phi[-1] > cmax:
                break
            W *= 2.0
        else:
            raise RuntimeError("window")
        worst = min(worst, float(np.nanmin(Phi)))
        dx = x[1] - x[0]
        for q in qs:
            c = q * q / 2
            g = np.maximum(c - Phi, 0.0) ** 1.5
            g[~np.isfinite(g)] = 0.0
            tot[q] += ws[m] * math.sqrt(s * (1 - s)) * dx * (g.sum() - 0.5 * (g[0] + g[-1]))
    return {q: 32 * math.sqrt(2) / (math.pi * q ** 4) * tot[q] for q in qs}, worst


def phase_value(eps, p, th):
    """A(q) for every q of QS at one phase, and the smallest Phi found (must not be below -1e-9)."""
    roots = B.branches(eps, p, np.array([th]), 1000)[0]
    vv = B.newton(roots, eps, p, np.full(roots.size, th), 4000)
    _, _, S = B.end_state(vv, eps, p, np.full(vv.size, th), 4000)
    Smin = float(S.min())
    near = vv[S - Smin < max(QS) ** 2 / 2]
    paths = B.integrate(near, eps, p, np.full(near.size, th), N_RK)[:, 0, :]   # (N_RK + 1, n)
    xs, _ = np.polynomial.legendre.leggauss(N_S)
    sn = 0.5 * (xs + 1.0)
    sg = np.linspace(0, 1, N_RK + 1)
    psi_rays = [np.interp(sn, sg, paths[:, j]) for j in range(paths.shape[1])]
    A, worst = phase_A(eps, p, th, Smin, psi_rays)
    return {"A": [A[q] for q in QS], "minPhi": worst}


NAME = "f68f_1b_exact"
CODE = C.code_md5([os.path.abspath(__file__), B.__file__])


def job(eps, p, th):
    return C.cached(NAME, CODE, phase_value, (eps, p, th))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--points", type=str, default="4:2", help="eps:p,eps:p,...")
    ap.add_argument("--n-theta", type=int, default=B.N_TH)
    ap.add_argument("--workers", type=int, default=0)
    a = ap.parse_args()
    t0 = time.time()
    pts = C.pairs(a.points)
    ths = [2 * math.pi * j / a.n_theta for j in range(a.n_theta)]
    jobs = [(e, p, th) for (e, p) in pts for th in ths]
    res = parallel_map(job, jobs, default_workers(a.workers), print,
                       lambda i, r, d, n: f"[f68f_1b_exact] {d}/{n}: eps {jobs[i][0]:g} p {jobs[i][1]:g} "
                                          f"theta {jobs[i][2]:.4f} ({time.time() - t0:.0f}s)")
    lines = [f"F6.8f-dev step 1b at finite q, b -> 0: exact needle with all branches (lower envelopes), mean over "
             f"{a.n_theta} phases; development, decides nothing", "",
             "delta_A = |A(1/8) - A(1/4)| / max(1, ...); pair = largest |A_j(1/8) - A_j(1/4)|; minPhi = smallest Phi "
             "found (must be >= -1e-9: the envelope is consistent with S_min)", "",
             "     eps    p     A(1/2)       A(1/4)       A(1/8)       delta_A    pair       minPhi"]
    data = []
    for k, (e, p) in enumerate(pts):
        R = res[k * a.n_theta:(k + 1) * a.n_theta]
        Aq = np.array([r["A"] for r in R])            # (n_theta, 3)
        m = Aq.mean(0)
        dA = abs(m[2] - m[1]) / max(1.0, abs(m[2]), abs(m[1]))
        pair = float(np.max(np.abs(Aq[:, 2] - Aq[:, 1])))
        worst = min(r["minPhi"] for r in R)
        lines.append(f"  {e:6.2f}  {p:3.1f}  {m[0]:.8f}  {m[1]:.8f}  {m[2]:.8f}  {dA:.1e}    {pair:.1e}    {worst:+.1e}")
        data.append({"eps": e, "p": p, "q": QS, "theta": ths, "A": Aq, "mean": m, "delta_A": dA, "pair": pair,
                     "minPhi": worst})
    lines.append(f"cost {time.time() - t0:.0f} s")
    stem = C.save(NAME, vars(a), CODE, data, lines, t0)
    print("\n".join(lines))
    print("saved", stem + ".json")


if __name__ == "__main__":
    main()
