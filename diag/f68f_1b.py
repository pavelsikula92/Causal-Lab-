"""
F6.8f-dev, step 1b (development, decides nothing): the plane-wave (Penrose) limit of the interval volume
=====================================================================================================

docs/derivations/f68f_dev.md, step 1b. Optical sector n = exp(b cos(k.w + theta)). A pair x -> y at distance L, with
kappa = kL -> infinity at fixed eps = b kappa^2 and p = k_par L, and sigma / L -> 0 with k sigma -> 0 (the needle of
the interval is thin compared with 1/k). In this limit (s = l / L in [0, 1], psi = k delta the transverse
displacement of the ray along k_perp, in units of the wavelength):

  * the ray (first arrival) minimises  S[psi] = int_0^1 [ psi'^2 / 2 + eps cos(p s + theta + psi) ] ds,
    (T - L = L S / kappa^2), with psi(0) = psi(1) = 0; stationary rays solve psi'' = -eps sin(p s + theta + psi)
    (a driven pendulum); there can be several (branches);
  * the transverse Jacobi field in the plane of k:  J'' = -eps cos(p s + theta + psi(s)) J;  out of that plane J'' = 0;
  * the interval is the needle around the first-arrival ray: a plane wave with that profile. Its volume relative to
    flat space at the same sigma is (derivation in the dev log)

        A = 6 int_0^1 sqrt( s (1 - s) E1(s) E2(s) / E1(1) ) ds,

    E1(0) = 0, E1'(0) = 1;  E2(1) = 0, E2'(1) = -1  (flat: E1 = s, E2 = 1 - s, A = 1).
    E1(1) = d psi(1) / d psi'(0) is the slope of the shooting map: E1(1) = 0 is a conjugate point of y, i.e. a fold of
    the family of rays from x (two branches are born or die there). The first-arrival ray has no interior conjugate
    point (it minimises S), so E1 > 0 on (0, 1] and E2 > 0 on [0, 1) are checked, not assumed.

The branch structure is found independently of the Jacobi fields: all roots of the shooting map psi(1; v) = 0,
v = psi'(0), in |v| <= eps + 2 sqrt(eps) + 1 (outside it S > eps >= S[0] >= min S, so no first arrival), on a grid in v
with every interval where E1(1) changes sign examined for a pair of close roots.

Numerics: classical RK4 vectorised over (theta, v), N steps on [0, 1]; the numerical error of A is the difference
between N and 2N (reported, criterion 1e-4 of the author). Newton on v for the roots.

Not here: the error of the Penrose approximation itself (finite k sigma and finite b), which needs the exact medium
(the brute force of 1a on long thin intervals): see the dev log.

Run: python diag/f68f_1b.py   (output diag/out/f68f_1b.txt)
"""
import math
import os
import sys
import time

import numpy as np

EPS_GRID = [0.0, 1.0, 2.0, 4.0, 8.0, 9.0, 9.5, 9.8, 9.87, 9.95, 10.0, 10.5, 12.0, 16.0, 40.0]
P_GRID = [0.0, 1.0, 2.0, 4.0, 8.0]
N_TH = 128
N_V = 400


def integrate(v, eps, p, th, N):
    """RK4 for (psi, psi', E, E', C, C', S) from s = 0, vectorised; returns the states at all N + 1 nodes."""
    v = np.asarray(v, float)
    th = np.broadcast_to(np.asarray(th, float), v.shape)
    y = np.zeros((7,) + v.shape)
    y[1] = v
    y[3] = 1.0
    y[4] = 1.0
    h = 1.0 / N
    out = np.empty((N + 1,) + y.shape)
    out[0] = y

    def f(s, y):
        ph = p * s + th + y[0]
        c = np.cos(ph)
        g = -eps * c
        return np.stack([y[1], -eps * np.sin(ph), y[3], g * y[2], y[5], g * y[4], 0.5 * y[1] ** 2 + eps * c])

    s = 0.0
    for i in range(N):
        k1 = f(s, y)
        k2 = f(s + h / 2, y + h / 2 * k1)
        k3 = f(s + h / 2, y + h / 2 * k2)
        k4 = f(s + h, y + h * k3)
        y = y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        s += h
        out[i + 1] = y
    return out


def end_state(v, eps, p, th, N):
    """psi(1), E1(1), S at s = 1 without storing the path."""
    v = np.asarray(v, float)
    th = np.broadcast_to(np.asarray(th, float), v.shape)
    y = np.zeros((5,) + v.shape)        # psi, psi', E, E', S
    y[1] = v
    y[3] = 1.0
    h = 1.0 / N

    def f(s, y):
        ph = p * s + th + y[0]
        c = np.cos(ph)
        return np.stack([y[1], -eps * np.sin(ph), y[3], -eps * c * y[2], 0.5 * y[1] ** 2 + eps * c])

    s = 0.0
    for i in range(N):
        k1 = f(s, y)
        k2 = f(s + h / 2, y + h / 2 * k1)
        k3 = f(s + h / 2, y + h / 2 * k2)
        k4 = f(s + h, y + h * k3)
        y = y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        s += h
    return y[0], y[2], y[4]


def newton(v, eps, p, th, N, it=30):
    v = np.array(v, float)
    for _ in range(it):
        F, dF, _ = end_state(v, eps, p, th, N)
        step = F / np.where(dF == 0, 1e-300, dF)
        v = v - np.clip(step, -0.5, 0.5)
        if np.max(np.abs(step)) < 1e-14:
            break
    return v


def _bisect(fun, th, a, b, it=55):
    """Vectorised bisection on intervals [a, b] (one per entry of th) for the sign change of fun(v, th)."""
    fa = fun(a, th)
    for _ in range(it):
        m = 0.5 * (a + b)
        fm = fun(m, th)
        left = fm * fa > 0
        a = np.where(left, m, a)
        fa = np.where(left, fm, fa)
        b = np.where(left, b, m)
    return 0.5 * (a + b)


def branches(eps, p, ths, N):
    """All roots v of psi(1; v) = 0 in the bounded window, per theta (list of arrays). Sign changes of psi(1) on a v
    grid; every grid interval where E1(1) = d psi(1)/dv changes sign without a sign change of psi(1) is resolved: the
    fold v* (E1(1) = 0) is found and, if psi(1; v*) has the opposite sign, the two roots on either side are added."""
    vm = eps + 2.0 * math.sqrt(eps) + 1.0
    vg = np.linspace(-vm, vm, N_V)
    TH, VG = np.meshgrid(ths, vg, indexing="ij")
    F, D, _ = end_state(VG, eps, p, TH, N)
    fpsi = lambda v, th: end_state(v, eps, p, th, N)[0]
    fE = lambda v, th: end_state(v, eps, p, th, N)[1]
    sc = (F[:, :-1] * F[:, 1:] < 0) | (F[:, :-1] == 0)
    fo = (~sc) & (D[:, :-1] * D[:, 1:] < 0)
    I, J = np.nonzero(sc)
    ti, aa, bb = [I], [vg[J]], [vg[J + 1]]
    I2, J2 = np.nonzero(fo)
    if I2.size:
        vs = _bisect(fE, ths[I2], vg[J2], vg[J2 + 1])
        Fs = fpsi(vs, ths[I2])
        two = Fs * F[I2, J2] < 0
        I2, J2, vs = I2[two], J2[two], vs[two]
        ti += [I2, I2]
        aa += [vg[J2], vs]
        bb += [vs, vg[J2 + 1]]
    ti = np.concatenate(ti)
    aa = np.concatenate(aa)
    bb = np.concatenate(bb)
    roots = [np.array([]) for _ in ths]
    if ti.size:
        r = _bisect(fpsi, ths[ti], aa, bb)
        for i in range(len(ths)):
            x = np.sort(r[ti == i])
            roots[i] = x
    return roots


def volume_A(v, eps, p, th, N):
    """A of the branches v (one per theta in th), and the minima of E1 on (0, 1] and E2 on [0, 1) (the conjugate-point
    check) and E1(1)."""
    v = np.atleast_1d(np.asarray(v, float))
    th = np.broadcast_to(np.asarray(th, float), v.shape)
    out = integrate(v, eps, p, th, N)                       # (N + 1, 7, n)
    s = np.linspace(0.0, 1.0, N + 1)[:, None]
    E1, dE1, C, dC = out[:, 2], out[:, 3], out[:, 4], out[:, 5]
    det = E1[-1] * dC[-1] - C[-1] * dE1[-1]
    a = C[-1] / det                                        # solves [E1 C; E1' C'] (a, c) = (0, -1)
    c = -E1[-1] / det
    E2 = a * E1 + c * C
    W = E1[-1]
    integ = np.sqrt(np.maximum(s * (1 - s) * E1 * E2 / W, 0.0))
    wS = np.ones(N + 1)
    wS[1:-1:2] = 4.0
    wS[2:-1:2] = 2.0
    A = 6.0 * (wS[:, None] * integ).sum(0) / (3.0 * N)
    return A, np.min(E1[1:], 0), np.min(E2[:-1], 0), W


def point(eps, p, N=4000):
    ths = 2 * math.pi * np.arange(N_TH) / N_TH
    if eps == 0.0:
        return {"A": 1.0, "dA": 0.0, "nbr_max": 1, "frac_multi": 0.0, "jumps": 0, "minW": 1.0, "Amin": 1.0,
                "Amax": 1.0, "conj_viol": 0}
    roots = branches(eps, p, ths, N // 4)
    ti = np.concatenate([np.full(len(r), i) for i, r in enumerate(roots)])
    vv = newton(np.concatenate(roots), eps, p, ths[ti], N) if ti.size else np.array([])
    _, _, S = end_state(vv, eps, p, ths[ti], N)
    vf = np.empty(len(ths))
    nbr = np.zeros(len(ths), int)
    for i in range(len(ths)):
        sel = ti == i
        nbr[i] = len(np.unique(np.round(vv[sel], 9)))
        vf[i] = vv[sel][np.argmin(S[sel])]
    A, m1, m2, W = volume_A(vf, eps, p, ths, N)
    vf2 = newton(vf, eps, p, ths, 2 * N)
    A2, _, _, _ = volume_A(vf2, eps, p, ths, 2 * N)
    jumps = int(np.sum(np.abs(np.diff(np.r_[vf, vf[0]])) > 0.5))
    return {"A": float(A.mean()), "dA": float(np.max(np.abs(A - A2))), "nbr_max": int(nbr.max()),
            "frac_multi": float(np.mean(nbr > 1)), "jumps": jumps, "minW": float(W.min()),
            "Amin": float(A.min()), "Amax": float(A.max()), "conj_viol": int(np.sum((m1 <= 0) | (m2 <= 0))),
            "A2": float(A2.mean())}


def eps_branch_onset(p, N=1000, lo=0.5, hi=40.0):
    """Smallest eps at which some theta has more than one ray (a fold of the family from x reaches y)."""
    ths = 2 * math.pi * np.arange(N_TH) / N_TH

    def multi(e):
        return max(len(r) for r in branches(e, p, ths, N)) > 1
    if not multi(hi):
        return None
    for _ in range(14):
        m = 0.5 * (lo + hi)
        if multi(m):
            hi = m
        else:
            lo = m
    return 0.5 * (lo + hi), hi - lo


def main():
    t0 = time.time()
    lines = ["F6.8f-dev step 1b: plane-wave (Penrose) limit of the interval volume; first-arrival ray; "
             "development, decides nothing", "",
             "A = theta average (128 phases) of V / V0 at fixed sigma; dA = max over theta |A(N) - A(2N)|, N = 4000 "
             "(numerical error of the Jacobi solution); nbr = largest number of rays x -> y in the window; multi = "
             "fraction of theta with more than one ray; jumps = switches of the first-arrival branch around theta "
             "(cut locus); minE1(1) = smallest E1(1) of the first-arrival ray; conj = phases with E1 or E2 <= 0 "
             "inside (must be 0)", ""]
    for p in P_GRID:
        lines.append(f"p = {p}")
        lines.append("     eps      A            Amin        Amax        dA         nbr  multi   jumps  minE1(1)    conj")
        for eps in EPS_GRID:
            r = point(eps, p)
            lines.append(f"  {eps:6.2f}  {r['A']:.8f}  {r['Amin']:.6f}  {r['Amax']:.6f}  {r['dA']:.1e}   {r['nbr_max']:3d}  "
                         f"{r['frac_multi']:.3f}   {r['jumps']:3d}   {r['minW']:+.3e}  {r['conj_viol']}")
            print(lines[-1], flush=True)
        ob = eps_branch_onset(p)
        lines.append(f"  onset of a second ray (fold reaches y), smallest over theta: eps = "
                     f"{'none below 40' if ob is None else f'{ob[0]:.4f} +- {ob[1]:.1e}'}")
        lines.append("")
        print(lines[-2], flush=True)
    lines.append(f"cost {time.time() - t0:.0f} s")
    txt = "\n".join(lines)
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "f68f_1b.txt"), "w") as fh:
        fh.write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
