"""
F6.8i-dev, D2b-nf (development, decides nothing): the age term of the cusp normal form. docs/derivations/f68i_dev.md
section 6.2.

D2a (section 6.1) found that near p_c the q^6 coefficient <s4hat A0^-2>_theta grows as |p - p_c|^-3.5 on both sides,
with min K ~ Delta, the 3-ray half-width w ~ Delta^1.5 and the window slope alpha ~ Delta^0.5. These are the exponents
of the cusp catastrophe: near p_c the needle well at every s is (with a the amplitude of the soft Jacobi mode)

    V(a) = a^4/4 + sig a^2/2 + tau a,     sig ~ (p - p_c),  tau ~ (theta - theta_c),

and the needle volume is A = N' w(c; sig, tau) / c^2 with w(c) = int da [c - (V(a) - Vmin)]_+^{3/2}. Then
A(c; sig, tau) = c^(-1/4) F(sig c^(-1/2), tau c^(-3/4)) exactly, and the age term integrated over the unfolding plane,

    I(gamma) = int dsig dtau J(sig, tau; gamma),     J = int_0^inf ds g0(s^2 A(gamma s; sig, tau)),

is homogeneous: I(gamma) = C gamma^(11/7) (s ~ gamma^(1/7), c ~ gamma^(8/7), sig ~ c^(1/2), tau ~ c^(3/4)). The
smooth part of the age term is O(gamma^3); a cusp term of order gamma^(11/7) would dominate it at small gamma. The BD
moments mu_0..mu_2 vanish (they remove the gamma^0..gamma^2 terms of an analytic needle) but they do not remove a
non-integer power: in the exchanged order C = (4/7) G(9/7) FP int F^(-9/7), G(z) = int v^(z-1) g0(v) dv.

This file computes C (canonical form: N' = 1, lambda = 1), in the correct order (s innermost), with checks:
  1. w(c) against the closed form at the cusp, the paraxial limit and brute-force quadrature of a tilted double well;
  2. J against doubled s nodes;
  3. I at two node sets (sig and tau nodes doubled), and the sig tails (power fit of int dtau J over |sig| in
     [Sig/2, Sig], added analytically);
  4. homogeneity: I(gamma) at gamma = 1 and 1/4 must have the ratio 4^(-11/7) (the grids are fixed, the structures
     move by factors 2^(8/7) in sig and 2^(12/7) in tau, so this tests the quadrature, not the scaling argument).

    python -m diag.f68i_d2b_nf --workers 8
Output diag/out/f68i_d2b_nf_<UTC time>.json and .txt.
"""
from __future__ import annotations

import argparse
import math
import os
import time
from multiprocessing import Pool

import numpy as np
from scipy.integrate import quad
from scipy.special import roots_jacobi

from diag import f68f_common as Cm
from phases.F6_order_dynamics.F6_08a_bd_action_conformal import BD_COEFFS, D, kernel_derivs

NJ = 64
XJ, WJ = roots_jacobi(NJ, 1.5, 1.5)
SIG_MAX = 16.0
SIG_BRK = [0.0, 1.0, 2.5, 5.0, 10.0, 16.0]
# t = tau / max(|sig|, 1)^1.5, up to t = 10: beyond, J (~ tau^(-10/3)) is below its rounding floor (~1e-11, the O(1)
# terms of the s integral cancel); the part beyond t = 10 is reported as a bound, not added
T_BRK_POS = [0.0, 0.1, 0.3, 0.6, 1.0, 2.0, 4.0, 10.0]
T_BRK_NEG = [0.1, 0.3, 2 / 3 ** 1.5, 0.6, 1.0, 2.0, 4.0, 10.0]  # beyond the window; 2/3^1.5: the fold line


def g0(w):
    return kernel_derivs(w)[0]


def G_mellin(z):
    return sum(c / math.factorial(i) * math.gamma(i + z) for i, c in enumerate(BD_COEFFS[D]))


# ------------------------------------------------------------------------------------------------ the normal form

def extrema(sig, tau):
    r = np.roots([1.0, 0.0, sig, tau])
    r = np.sort(r[np.abs(r.imag) < 1e-9 * (1 + np.abs(r.real))].real)
    for _ in range(3):
        d = 3 * r * r + sig
        ok = np.abs(d) > 1e-12
        r = np.where(ok, r - (r ** 3 + sig * r + tau) / np.where(ok, d, 1.0), r)
    return r, r ** 4 / 4 + sig * r * r / 2 + tau * r


def wells(sig, tau):
    """Vmin, the depth of the second minimum above Vmin and the barrier height above Vmin (None for one well)."""
    r, V = extrema(sig, tau)
    vmin = float(V.min())
    if r.size == 3:
        return vmin, float(max(V[0], V[2]) - vmin), float(V[1] - vmin)
    return vmin, None, None


def w_of_c(sig, tau, c, vmin=None):
    """w(c) = int da [c - (V - Vmin)]_+^{3/2}, vectorised over c: the real roots of V - Vmin = c bound one or two
    intervals; on each, c - U = (1/4)(a - a1)(a2 - a) Q(a) with Q > 0, Gauss-Jacobi (3/2, 3/2) in a."""
    c = np.atleast_1d(np.asarray(c, float))
    if vmin is None:
        vmin = wells(sig, tau)[0]
    n = c.size
    M = np.zeros((n, 4, 4))
    M[:, 0, 1] = -2 * sig
    M[:, 0, 2] = -4 * tau
    M[:, 0, 3] = 4 * (vmin + c)
    M[:, 1, 0] = M[:, 2, 1] = M[:, 3, 2] = 1.0
    z = np.linalg.eigvals(M)
    isr = np.abs(z.imag) < 1e-7 * (1 + np.abs(z.real))
    re = np.sort(np.where(isr, z.real, np.inf), axis=1)
    k = isr.sum(1)
    fin = np.isfinite(re)
    x = np.where(fin, re, 0.0)
    cc = c[:, None]
    for _ in range(2):                       # Newton polish of the real roots
        f = x ** 4 / 4 + sig * x * x / 2 + tau * x - vmin - cc
        df = x ** 3 + sig * x + tau
        ok = fin & (np.abs(df) > 1e-12)
        x = np.where(ok, x - f / np.where(ok, df, 1.0), x)
    out = np.zeros(n)
    for lo, hi, need in ((0, 1, k >= 2), (2, 3, k == 4)):
        if not np.any(need):
            continue
        a1, a2 = x[need, lo], x[need, hi]
        m, h = 0.5 * (a1 + a2), 0.5 * (a2 - a1)
        a = m[:, None] + h[:, None] * XJ[None, :]
        Q = a * a + (a1 + a2)[:, None] * a + (a1 * a1 + a1 * a2 + a2 * a2 + 2 * sig)[:, None]
        out[need] += h ** 4 / 8.0 * (np.maximum(Q, 0.0) ** 1.5 * WJ[None, :]).sum(1)
    return out


def J(sig, tau, gam=1.0, ngl=48):
    """int_0^inf ds g0(s^2 A(gam s)), A = w / c^2; s split at the well onsets (second well, barrier top),
    graded towards 0, and ended where s^2 A > 70."""
    vmin, d2, bar = wells(sig, tau)

    def A(s):
        cc = gam * np.asarray(s, float)
        return w_of_c(sig, tau, cc, vmin) / cc ** 2

    smax = 0.5
    while smax * smax * A([smax])[0] <= 70.0:
        smax *= 1.5
    brk = [b / gam for b in (d2, bar) if b is not None and 0 < b / gam < smax]
    first = min(brk + [smax])
    pts = sorted(set([0.0, smax] + brk + [first * 4.0 ** (-j) for j in range(1, 8)]))
    x, wt = np.polynomial.legendre.leggauss(ngl)
    tot = 0.0
    for a, b in zip(pts[:-1], pts[1:]):
        s = 0.5 * (b - a) * (x + 1) + a
        tot += float(np.sum(0.5 * (b - a) * wt * g0(s * s * A(s))))
    return tot


# ------------------------------------------------------------------------------------------------ the plane

def tau_nodes(sig, nt, gam):
    """Nodes and weights in tau >= 0 (J is even in tau), and the tail exponent region end."""
    s3 = max(abs(sig), 1.0) ** 1.5
    x, w = np.polynomial.legendre.leggauss(nt)
    segs = []
    if sig < 0:
        tw = 2.5 * abs(sig) ** -0.25 * gam ** 0.5 * max(1.0, gam ** 0.5)   # the window (c ~ gamma s, s ~ |sig|^1/4)
        tw = min(tw, 0.1 * s3)
        segs.append((0.0, tw, 2 * nt))
        tb = [tw] + [t * s3 for t in T_BRK_NEG if t * s3 > tw]
    else:
        tb = [t * s3 for t in T_BRK_POS]
    segs += [(a, b, nt) for a, b in zip(tb[:-1], tb[1:])]
    T, W = [], []
    for a, b, n in segs:
        xx, ww = np.polynomial.legendre.leggauss(n)
        T.append(0.5 * (b - a) * (xx + 1) + a)
        W.append(0.5 * (b - a) * ww)
    return np.concatenate(T), np.concatenate(W), tb[-1]


def row_job(args):
    """2 int_0^T dtau J(sig, tau) at one sig, T = 10 max(|sig|, 1)^1.5; the tail beyond as a bound |J(T)| T / (7/3)."""
    sig, nt, gam = args
    T, W, tend = tau_nodes(sig, nt, gam)
    vals = np.array([J(sig, t, gam) for t in T])
    jend = J(sig, tend, gam)
    tail = abs(jend) * tend / (7.0 / 3.0)
    return {"sig": sig, "int": 2.0 * float(np.sum(W * vals)), "tail_bound": 2.0 * tail,
            "maxabs": float(np.max(np.abs(vals)))}


def plane(pool, ns, nt, gam):
    brk = sorted(set([-b for b in SIG_BRK] + SIG_BRK))
    x, w = np.polynomial.legendre.leggauss(ns)
    S, Wt = [], []
    for a, b in zip(brk[:-1], brk[1:]):
        S += list(0.5 * (b - a) * (x + 1) + a)
        Wt += list(0.5 * (b - a) * w)
    rows = pool.map(row_job, [(s, nt, gam) for s in S])
    core = float(sum(wt * r["int"] for wt, r in zip(Wt, rows)))
    # sig tails: fit |int| ~ B |sig|^-p on [Sig/2, Sig] on each side, add B Sig^(1-p)/(p-1)
    tails = {}
    for side in (-1, 1):
        pts = [(abs(r["sig"]), r["int"]) for r in rows if side * r["sig"] >= SIG_MAX / 2]
        xs = np.log([p[0] for p in pts])
        ys = [p[1] for p in pts]
        sgn = np.sign(ys[-1])
        same = all(np.sign(y) == sgn for y in ys)
        if same:
            pfit = np.polyfit(xs, np.log(np.abs(ys)), 1)
            p = -pfit[0]
            B = sgn * math.exp(pfit[1])
            tails[side] = {"p": float(p), "B": float(B),
                           "tail": float(B * SIG_MAX ** (1 - p) / (p - 1)) if p > 1 else float("nan")}
        else:
            tails[side] = {"p": float("nan"), "B": float("nan"), "tail": float("nan"), "note": "sign changes"}
    total = core + sum(t["tail"] for t in tails.values())
    tb = float(sum(wt * r["tail_bound"] for wt, r in zip(Wt, rows)))
    return {"core": core, "tails": tails, "total": total, "tau_tail_bound": tb, "rows": rows}


# ------------------------------------------------------------------------------------------------ checks

def checks(lines):
    out = {}
    B = quad(lambda x: (1 - x ** 4) ** 1.5, -1, 1, epsabs=1e-14, epsrel=1e-14)[0]
    r1 = [abs(w_of_c(0, 0, c)[0] / ((4 * c) ** 0.25 * c ** 1.5 * B) - 1) for c in (1e-4, 1e-2, 1.0, 10.0)]
    r2 = abs(w_of_c(5.0, 0, 1e-7)[0] / 1e-14 / (3 * math.pi / 8 * math.sqrt(2 / 5)) - 1)
    r3 = []
    for sig, tau in ((-2.0, 0.3), (-5.0, 0.05), (1.0, -0.7)):
        vmin = wells(sig, tau)[0]
        r = extrema(sig, tau)[0]
        for c in (0.2, 1.5, 6.0):
            bf = quad(lambda a: max(c - (a ** 4 / 4 + sig * a * a / 2 + tau * a - vmin), 0) ** 1.5, -6, 6,
                      limit=800, points=sorted(set(np.round(r, 6))), epsabs=1e-13, epsrel=1e-12)[0]
            r3.append(abs(w_of_c(sig, tau, c)[0] / bf - 1))
    js = [(sig, tau, J(sig, tau), J(sig, tau, ngl=96)) for sig, tau in ((0, 0), (3, 1), (-10, 0.1), (30, 0))]
    r4 = max(abs(a - b) / max(abs(b), 1e-300) for _, _, a, b in js)
    out.update(cusp_closed=max(r1), paraxial=r2, brute_double_well=max(r3), J_48_96=r4, B=B)
    lines += ["Checks of the normal form:",
              f"  w at the cusp against (4c)^(1/4) c^(3/2) int(1 - x^4)^(3/2): max rel {max(r1):.1e}",
              f"  paraxial w -> c^2 (3 pi / 8) sqrt(2 / sig) at sig = 5, c = 1e-7: rel {r2:.1e}",
              f"  w against adaptive quadrature (3 tilted wells x 3 c): max rel {max(r3):.1e}",
              "  J at 48 against 96 s nodes per segment: " + ", ".join(
                  f"({s:g}, {t:g}) {a:+.6e}/{b:+.6e}" for s, t, a, b in js) + f"; max rel {r4:.1e}", ""]
    return out


NAME = "f68i_d2b_nf"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ns", type=str, default="12,20", help="Gauss-Legendre nodes per sig segment (two sets)")
    ap.add_argument("--nt", type=str, default="20,32", help="nodes per tau segment (window: twice)")
    ap.add_argument("--gammas", type=str, default="1,0.25")
    ap.add_argument("--workers", type=int, default=0)
    a = ap.parse_args()
    t0 = time.time()
    import phases.F6_order_dynamics.F6_08a_bd_action_conformal as KER
    code = Cm.code_md5([os.path.abspath(__file__), KER.__file__])
    lines = ["F6.8i-dev D2b-nf: the age term of the cusp normal form V = a^4/4 + sig a^2/2 + tau a, A = w(c)/c^2 "
             "(development, decides nothing; docs/derivations/f68i_dev.md section 6.2)", ""]
    data = {"G_9_7": G_mellin(9 / 7), "moments_G": {str(z): G_mellin(z) for z in (0.5, 1.0, 1.5, 2.0)}}
    lines.append(f"G(z) = int v^(z-1) g0(v) dv: G(1/2) {G_mellin(.5):+.1e}, G(1) {G_mellin(1):+.1e}, "
                 f"G(3/2) {G_mellin(1.5):+.1e}, G(2) {G_mellin(2):+.6f}, G(9/7) {G_mellin(9 / 7):+.6f}")
    lines.append("")
    data["checks"] = checks(lines)
    print("\n".join(lines))
    nss = [int(x) for x in a.ns.split(",")]
    nts = [int(x) for x in a.nt.split(",")]
    gams = [float(x) for x in a.gammas.split(",")]
    res = []
    nw = a.workers or os.cpu_count()
    with Pool(nw) as pool:
        for gam in gams:
            for ns, nt in zip(nss, nts):
                t1 = time.time()
                r = plane(pool, ns, nt, gam)
                r.update(gamma=gam, ns=ns, nt=nt, C=r["total"] / gam ** (11 / 7))
                res.append(r)
                tl = "; ".join(f"{'sig<0' if s < 0 else 'sig>0'} tail p {t['p']:.2f} {t['tail']:+.3e}"
                               for s, t in r["tails"].items())
                ln = (f"gamma {gam:g}, nodes sig {ns} tau {nt}: core {r['core']:+.8e}  ({tl}; tau beyond t = 10: "
                      f"|.| < {r['tau_tail_bound']:.1e})  I {r['total']:+.8e}"
                      f"  C = I / gamma^(11/7) = {r['C']:+.6e}  ({time.time() - t1:.0f} s)")
                lines.append(ln)
                print(ln, flush=True)
    lines.append("")
    byg = {}
    for r in res:
        byg.setdefault(r["gamma"], []).append(r)
    for g, rr in byg.items():
        if len(rr) > 1:
            lines.append(f"gamma {g:g}: node sets differ by {abs(rr[-1]['total'] - rr[-2]['total']) / abs(rr[-1]['total']):.1e} (rel)")
    if len(byg) > 1:
        gs = sorted(byg)
        I0, I1 = byg[gs[0]][-1]["total"], byg[gs[-1]][-1]["total"]
        expo = math.log(I0 / I1) / math.log(gs[0] / gs[-1]) if I0 * I1 > 0 else float("nan")
        lines.append(f"homogeneity: I({gs[0]:g}) / I({gs[-1]:g}) gives the exponent {expo:.5f} (11/7 = {11 / 7:.5f})")
        data["exponent"] = expo
    data["runs"] = [{k: v for k, v in r.items()} for r in res]
    lines.append(f"total {time.time() - t0:.0f} s")
    stem = Cm.save(NAME, vars(a), code, data, lines, t0)
    print("\n".join(lines[-4:]))
    print("saved", stem + ".json")


if __name__ == "__main__":
    main()
