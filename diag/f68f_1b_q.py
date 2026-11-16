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

Two ways of computing the psi integral:
  * expansion (all 128 phases): around each ray b, Phi = S_b - S_min + K x^2/2 + c3 x^3 + c4 x^4 + ..., x = psi - psi_b,
    with K, c3, c4 from the variational equations of order 1-3 along the ray (Hamilton-Jacobi derivatives of G1, G2).
    Then  int dx [c - Phi_b]^{3/2} = (3 pi / 8) c^2 sqrt(2/K) [1 + c (5 c3^2 / (2 K^3) - c4 / K^2) + O(c^2)],
    c = q^2/2 - (S_b - S_min). Every ray with S_b - S_min < q^2/2 contributes its own needle (it is in the past of y
    at that sigma). A ray with K <= 0 somewhere inside (a saddle of the action) is flagged, not integrated.
  * exact (check on a few phases): Phi from the boundary-value problems x -> (s, psi) and (s, psi) -> y solved by
    Newton for every node (s, psi), first-arrival branch by continuation from the ray, psi integral on the exact
    support. No expansion.

Run: python diag/f68f_1b_q.py [check]   (output diag/out/f68f_1b_q.txt, or f68f_1b_q_check.txt)
"""
import math
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import f68f_1b as B  # noqa: E402

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


def point(eps, p, N=4000):
    ths = 2 * math.pi * np.arange(B.N_TH) / B.N_TH
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


# ----------------------------------------------------------------------------- exact check (no expansion)

def bvp_action(sa, psa, sb, psb, v0, eps, p, th, N=1024, it=25):
    """Action of the ray from (sa, psa) to (sb, psb), Newton on the slope from v0 (continuation); vectorised."""
    v = np.array(v0, float)
    L = sb - sa
    for _ in range(it):
        y = np.zeros((5,) + v.shape)
        y[0] = psa
        y[1] = v
        y[3] = 1.0
        h = 1.0 / N

        def f(u, y):
            s = sa + L * u
            ph = p * s + th + y[0]
            c = np.cos(ph)
            # derivatives in u: psi_u = L psi_s
            return np.stack([L * y[1], -L * eps * np.sin(ph), L * y[3], -L * eps * c * y[2], L * (0.5 * y[1] ** 2 + eps * c)])
        u = 0.0
        for i in range(N):
            k1 = f(u, y)
            k2 = f(u + h / 2, y + h / 2 * k1)
            k3 = f(u + h / 2, y + h / 2 * k2)
            k4 = f(u + h, y + h * k3)
            y = y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            u += h
        F = y[0] - psb
        step = F / y[2]
        v = v - np.clip(step, -0.5, 0.5)
        if np.max(np.abs(step)) < 1e-13:
            break
    return y[4], v


def exact_check(eps, p, th, qs=QS, n_s=24, n_psi=24):
    """A(q) by the exact needle for one phase (first-arrival ray only, no other branch)."""
    roots = B.branches(eps, p, np.array([th]), 1000)[0]
    vv = B.newton(roots, eps, p, np.full(roots.size, th), 4000)
    _, _, Sr = B.end_state(vv, eps, p, np.full(vv.size, th), 4000)
    v = vv[np.argmin(Sr)]
    Smin = float(np.min(Sr))
    xs, wsg = np.polynomial.legendre.leggauss(n_s)
    sn = 0.5 * (xs + 1.0)
    ws = 0.5 * wsg
    # the ray and its slope at the nodes
    out = B.integrate(np.array([v]), eps, p, th, 4000)[:, :, 0]
    sg = np.linspace(0, 1, 4001)
    psi_r = np.interp(sn, sg, out[:, 0])
    dpsi_r = np.interp(sn, sg, out[:, 1])
    _, K, c3, c4, fo, bo = needle_terms(np.array([v]), eps, p, np.array([th]), 4000)
    Kn = np.interp(sn, sg, K[:, 0])
    res = {}
    for q in qs:
        c = q * q / 2.0
        hw = 1.6 * math.sqrt(2 * c / 1.0) * np.sqrt(1.0 / Kn)          # search half-width (paraxial x 1.6)
        t = np.cos(math.pi * (np.arange(n_psi) + 0.5) / n_psi)          # Chebyshev nodes on [-1, 1]
        SN = np.repeat(sn, n_psi)
        X = (hw[:, None] * t[None, :]).ravel()
        PS = np.repeat(psi_r, n_psi) + X
        # G1: from (0, 0) to (s, psi); warm start v + x / E1(s) with E1 from the ray
        E1n = np.interp(sn, sg, fo[:, 2, 0])
        G1, _ = bvp_action(0.0, 0.0, SN, PS, v + X / np.repeat(E1n, n_psi), eps, p, th)
        # G2: from (s, psi) to (1, 0); slope at s warm start dpsi_r + x * (d slope / d psi) ~ -K2
        G2, _ = bvp_action(SN, PS, 1.0, 0.0, np.repeat(dpsi_r, n_psi), eps, p, th)
        Phi = (G1 + G2 - Smin).reshape(n_s, n_psi)
        tot = 0.0
        for j in range(n_s):
            # interpolate Phi in x by Chebyshev, integrate [c - Phi]^{3/2} on a fine grid with the exact support
            co = np.polynomial.chebyshev.chebfit(t, Phi[j], n_psi - 1)
            tt = np.linspace(-1, 1, 20001)
            ph = np.polynomial.chebyshev.chebval(tt, co)
            g = np.maximum(c - ph, 0.0) ** 1.5
            if g[0] > 0 or g[-1] > 0:
                raise RuntimeError("support wider than the search window")
            tot += ws[j] * math.sqrt(sn[j] * (1 - sn[j])) * hw[j] * np.trapezoid(g, tt)
        res[q] = 32 * math.sqrt(2) / (math.pi * q ** 4) * tot
    return res


def main():
    t0 = time.time()
    if len(sys.argv) > 1 and sys.argv[1] == "check":
        lines = ["F6.8f-dev step 1b, finite q: exact needle against the expansion, single phases (development)", ""]
        for (eps, p, th) in [(4.0, 0.0, 0.7), (4.0, 2.0, 2.0), (12.0, 1.0, 0.3), (16.0, 4.0, 4.0), (40.0, 8.0, 1.0)]:
            ex = exact_check(eps, p, th)
            ths = np.array([th])
            roots = B.branches(eps, p, ths, 1000)[0]
            vv = B.newton(roots, eps, p, np.full(roots.size, th), 4000)
            S, K, c3, c4, _, _ = needle_terms(vv, eps, p, np.full(vv.size, th), 4000)
            j = int(np.argmin(S))
            row = [f"eps {eps:5.1f} p {p:3.1f} theta {th:.2f}:"]
            for q in QS:
                Ae, _ = A_expansion(K[:, j:j + 1], c3[:, j:j + 1], c4[:, j:j + 1], np.ones(1), q, 4000)
                row.append(f"q {q:.3f} exact {ex[q]:.7f} expansion {Ae[0]:.7f} diff {ex[q] - Ae[0]:+.1e}")
            lines.append("  ".join(row))
            print(lines[-1], flush=True)
        lines.append(f"cost {time.time() - t0:.0f} s")
        name = "f68f_1b_q_check.txt"
    else:
        lines = ["F6.8f-dev step 1b at finite q = k sigma, b -> 0: mean over 128 phases of A(q) (expansion to O(q^2) "
                 "per needle, every ray within q^2/2 of the minimum action included); development, decides nothing", "",
                 "delta_A = |A(1/8) - A(1/4)| / max(1, |A(1/8)|, |A(1/4)|) (criterion 1e-3); pair = largest |A_j(1/8) - "
                 "A_j(1/4)| over phases; extra = needles of other rays included at q = 1/2, 1/4, 1/8; saddle = rays within "
                 "q^2/2 with K <= 0 somewhere (not integrated); corr = largest relative q^2 correction at q = 1/4", ""]
        for p in B.P_GRID:
            lines.append(f"p = {p}")
            lines.append("     eps     A_PW        A(1/2)      A(1/4)      A(1/8)      delta_A    pair       extra     saddle  corr")
            for eps in B.EPS_GRID:
                res, fl, apw = point(eps, p)
                a = {q: float(res[q].mean()) for q in QS}
                dA = abs(a[0.125] - a[0.25]) / max(1.0, abs(a[0.125]), abs(a[0.25]))
                pair = float(np.max(np.abs(res[0.125] - res[0.25])))
                ex = "/".join(str(fl["extra"][q]) for q in QS)
                lines.append(f"  {eps:6.2f}  {apw.mean():.7f}  {a[0.5]:.7f}  {a[0.25]:.7f}  {a[0.125]:.7f}  {dA:.1e}    "
                             f"{pair:.1e}    {ex:9s} {fl['saddle']:3d}     {fl['maxcorr'][0.25]:.1e}")
                print(lines[-1], flush=True)
            lines.append("")
        lines.append(f"cost {time.time() - t0:.0f} s")
        name = "f68f_1b_q.txt"
    txt = "\n".join(lines)
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", name), "w") as fh:
        fh.write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
