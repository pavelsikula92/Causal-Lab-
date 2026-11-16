"""
F6.8f-dev, step 0 (development, decides nothing): the exact first-arrival time in the optical sector at finite b
==============================================================================================================

The author's development branch of 16 Nov 2026 (docs/derivations/substrate_and_action_audit.md, section 8): the
optical response of the expected BD action at finite amplitude b. The second-order instrument of F6.8e uses the
first-order deviation of the rays (z = w + b delta) and cannot be continued to finite b. This file builds the piece
that replaces it: the exact optical length (first-arrival time) between two points.

Medium. The optical sector of F6.8c: g = -dt^2 + n^2 dx^2, nu = ln n = b cos(k z + th). The medium is stratified
(n depends on z only), so along a ray the transverse slowness q = n sin(theta) (theta from the z-axis) is conserved.
In units k = 1 (u = k z), between heights u1 and u2 with transverse distance X:
    X = int q / sqrt(n^2 - q^2) du,     t = q X + W,     W = int sqrt(n^2 - q^2) du,
the integrals taken along the ray (turning points where n(u) = q). The first arrival is the smallest t over all ray
families that reach the target, so caustics and multiple arrivals are included exactly, with no expansion in b.
d t / d b at fixed end points (Fermat) = int cos(u) n dl = int cos(u) n^2 / sqrt(n^2 - q^2) du.

Families. q < n_min = e^-|b|: monotone in u (no turning). q in [n_min, min(n(u1), n(u2))]: the ray is confined to
the allowed interval [m - a, m + a] around the maximum m of n containing u1 and u2 (cos a = ln q / b) and bounces;
unfolding its phase gives the families (start up/down, end up/down, j half-period pairs).
Turning-point integrals use u = m - a cos(psi), which removes the square-root end singularities, and
cos u - cos a = 2 sin(a cos^2(psi/2)) sin(a sin^2(psi/2)) (no cancellation).

Run: python diag/f68f_dev.py   (checks T1-T5; output diag/out/f68f_step0.txt)
"""
import math
import sys
import os
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

QOPT = dict(epsabs=1e-15, epsrel=1e-13, limit=400)


# ----------------------------------------------------------------------------- integrands (k = 1, b > 0)

def _D(psi, a):
    """cos u - cos a for u = m - a cos psi."""
    return 2.0 * math.sin(a * math.cos(0.5 * psi) ** 2) * math.sin(a * math.sin(0.5 * psi) ** 2)


def _turn_integrals(a, q, b, psi_hi):
    """F, G, H from the lower turning point (psi = 0) to psi_hi on the interval of half-width a."""
    def fF(psi):
        e = math.expm1(2.0 * b * _D(psi, a))
        return a * math.sin(psi) / math.sqrt(e) if e > 0 else 0.0

    def fG(psi):
        e = math.expm1(2.0 * b * _D(psi, a))
        return q * a * math.sin(psi) * math.sqrt(max(e, 0.0))

    def fH(psi):
        D = _D(psi, a)
        e = math.expm1(2.0 * b * D)
        if e <= 0:
            return 0.0
        cu = math.cos(a * math.cos(psi))           # cos(u - m) = cos u (m a multiple of 2 pi)
        return cu * q * math.exp(2.0 * b * D) / math.sqrt(e) * a * math.sin(psi)
    F = quad(fF, 0.0, psi_hi, **QOPT)[0]
    G = quad(fG, 0.0, psi_hi, **QOPT)[0]
    H = quad(fH, 0.0, psi_hi, **QOPT)[0]
    return F, G, H


def _mono_integrals(u1, u2, lam, b):
    """Monotone ray from u1 to u2 (u2 >= u1) with q = e^(-b - lam), lam > 0 (no turning point).
    n^2 - q^2 = q^2 expm1(2 b (1 + cos u) + 2 lam), written without cancellation."""
    q = math.exp(-b - lam)
    ee = lambda u: math.expm1(4.0 * b * math.cos(0.5 * u) ** 2 + 2.0 * lam)
    fF = lambda u: 1.0 / math.sqrt(ee(u))
    fG = lambda u: q * math.sqrt(ee(u))
    fH = lambda u: math.cos(u) * q * (1.0 + ee(u)) / math.sqrt(ee(u))
    pts = [p for p in np.arange(math.ceil((u1 - math.pi) / (2 * math.pi)), math.floor((u2 - math.pi) / (2 * math.pi)) + 1)
           * 2 * math.pi + math.pi if u1 < p < u2]
    kw = dict(QOPT, points=pts) if pts else QOPT
    return (quad(fF, u1, u2, **kw)[0], quad(fG, u1, u2, **kw)[0], quad(fH, u1, u2, **kw)[0]), q


def _vert(u1, u2, b):
    return quad(lambda u: math.exp(b * math.cos(u)), u1, u2, **QOPT)[0]


def _vertH(u1, u2, b):
    return quad(lambda u: math.cos(u) * math.exp(b * math.cos(u)), u1, u2, **QOPT)[0]


# ----------------------------------------------------------------------------- arrivals

def arrivals(u1, u2, X, b, n_grid=240, jmax=None):
    """All ray arrivals (t, dt/db, q, family) from height u1 to u2 at transverse distance X >= 0; k = 1."""
    if b == 0.0:
        return [(math.hypot(X, u2 - u1), 0.0, None, "flat")]
    sgn = 1.0
    if b < 0:                                   # n(u; b) = n(u + pi; -b): work with b > 0
        u1, u2, b, sgn = u1 + math.pi, u2 + math.pi, -b, -1.0
    if u2 < u1:
        u1, u2 = u2, u1
    out = []
    nmin = math.exp(-b)
    # monotone family
    if u2 > u1:
        def Xm(lam):
            return _mono_integrals(u1, u2, lam, b)[0][0] - X
        lam_lo = 1e-13
        if X == 0.0:
            out.append((_vert(u1, u2, b), sgn * _vertH(u1, u2, b), 0.0, "mono"))
        elif 0.0 <= Xm(lam_lo):
            lam = brentq(Xm, lam_lo, 300.0, xtol=1e-16, rtol=1e-15)
            (F, G, H), q = _mono_integrals(u1, u2, lam, b)
            out.append((q * X + G, sgn * H, q, "mono"))
    elif X == 0.0:
        out.append((0.0, 0.0, None, "point"))
    # turning families: both ends in the same allowed interval
    m1 = 2 * math.pi * round(u1 / (2 * math.pi))
    m2 = 2 * math.pi * round(u2 / (2 * math.pi))
    if m1 != m2:
        return out
    m = m1
    n1, n2 = math.exp(b * math.cos(u1)), math.exp(b * math.cos(u2))
    qhi = min(n1, n2)
    if u1 == u2 == m:                           # the ray along the layer of largest index
        out.append((qhi * X, sgn * qhi * X, qhi, "axis"))

    def fam_vals(a):
        q = math.exp(b * math.cos(a))
        ps = lambda u: math.acos(max(-1.0, min(1.0, (m - u) / a)))
        F1, G1, H1 = _turn_integrals(a, q, b, ps(u1))
        F2, G2, H2 = _turn_integrals(a, q, b, ps(u2))
        P, PG, PH = _turn_integrals(a, q, b, math.pi)
        return (F1, G1, H1), (F2, G2, H2), (P, PG, PH)

    def families(v, J):
        (F1, G1, H1), (F2, G2, H2), (P, PG, PH) = v
        starts = {"up": (F1, G1, H1), "down": (2 * P - F1, 2 * PG - G1, 2 * PH - H1)}
        res = {}
        for sk, s in starts.items():
            for j in range(-1, J + 1):
                for ek, e in (("up", (F2 + 2 * P * j, G2 + 2 * PG * j, H2 + 2 * PH * j)),
                              ("down", (2 * P - F2 + 2 * P * j, 2 * PG - G2 + 2 * PG * j, 2 * PH - H2 + 2 * PH * j))):
                    res[(sk, ek, j)] = tuple(e[i] - s[i] for i in range(3))
        return res

    # a = half-width of the allowed interval: from max(|u1 - m|, |u2 - m|) (q = qhi) to pi (q = n_min)
    alo = max(abs(u1 - m), abs(u2 - m), 1e-9)
    ahi = math.pi * (1.0 - 1e-12)
    if ahi <= alo:
        return out
    xs = np.linspace(0.0, math.pi, n_grid + 1)
    qs = alo + (ahi - alo) * 0.5 * (1.0 - np.cos(xs))     # grid in a (named qs for the loop below)
    vals = [fam_vals(a) for a in qs]
    Pmin = min(v[2][0] for v in vals)
    J = jmax if jmax is not None else int(X / (2 * Pmin)) + 1
    fams = [families(v, J) for v in vals]
    for key in fams[0]:
        g = np.array([f[key][0] - X for f in fams])
        for i in range(len(qs) - 1):
            if g[i] == 0.0 or g[i] * g[i + 1] < 0:
                fx = lambda q: families(fam_vals(q), J)[key][0] - X
                aa = brentq(fx, qs[i], qs[i + 1], xtol=1e-16, rtol=1e-15)
                q = math.exp(b * math.cos(aa))
                Xf, W, H = families(fam_vals(aa), J)[key]
                if Xf >= -1e-12 and W >= -1e-12:
                    out.append((q * X + W, sgn * H, q, key))
    return out


def first_arrival(z1, z2, X, b, k=1.0, th=0.0):
    """Exact optical length (first arrival) between (z1, X = 0) and (z2, X) in n = exp(b cos(k z + th)); returns
    (t, dt/db)."""
    arr = arrivals(k * z1 + th, k * z2 + th, k * X, b)
    t, h = min((a[0], a[1]) for a in arr)
    return t / k, h / k


# ----------------------------------------------------------------------------- checks

def _f68c():
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
    from phases.F6_order_dynamics import F6_08c_bd_action_optical as c
    return c


def check_T1():
    """b -> 0: flat."""
    rows = []
    for (z1, z2, X) in ((0.3, 2.1, 1.7), (0.0, 0.0, 4.0), (-1.0, 3.0, 0.0)):
        t, _ = first_arrival(z1, z2, X, 1e-9, k=0.8)
        rows.append((z1, z2, X, t, math.hypot(X, z2 - z1), abs(t - math.hypot(X, z2 - z1))))
    return rows


def check_T2():
    """The phase-averaged second-order coefficient against the closed form t2_mean of F6.8c."""
    c = _f68c()
    rows = []
    k = 0.9
    for d, ang in ((3.0, 0.4), (5.0, 1.1), (2.0, 1.5), (40.0, 1.2)):
        dz, X = d * math.cos(ang), d * math.sin(ang)
        nph = 16
        acc = 0.0
        for i in range(nph):
            th = 2 * math.pi * i / nph
            eps = 1e-4
            tp = first_arrival(0.0, dz, X, eps, k, th)[0]
            tm = first_arrival(0.0, dz, X, -eps, k, th)[0]
            acc += (tp + tm - 2 * d) / (2 * eps * eps)
        num = acc / nph
        ana = float(c.t2_mean(k, np.array(d), np.array(k * dz)))
        rows.append((d, ang, num, ana, abs(num - ana) / abs(ana)))
    return rows


def check_T3():
    """Moderate b (finite, no caustic) against the polyline Fermat minimisation of F6.8c (Richardson in n)."""
    c = _f68c()
    rows = []
    for (b, k, d, ang, th) in ((0.3, 1.0, 5.0, 1.0, 0.4), (0.6, 0.7, 6.0, 0.6, 2.0)):
        A = np.zeros(3)
        kvec = k * np.array([1.0, 0.0, 0.0])                 # z-axis of this file = x-axis there
        B = d * np.array([math.cos(ang), math.sin(ang), 0.0])
        tf = [c._fermat_T(A, B, kvec, th, b, n=nn) for nn in (200, 400)]
        tr = (4 * tf[1] - tf[0]) / 3
        t, _ = first_arrival(0.0, d * math.cos(ang), d * math.sin(ang), b, k, th)
        rows.append((b, k, d, ang, th, t, tf[1], tr, abs(t - tr) / t))
    return rows


def check_T4(b=1e-2):
    """The first caustic: from the layer of largest index the bounced family exists from X = pi/(k sqrt b) on and
    arrives before the axis ray."""
    rows = []
    Xc = math.pi / math.sqrt(b)
    for f in (0.9, 0.99, 1.01, 1.1, 1.5):
        arr = arrivals(0.0, 0.0, f * Xc, b)
        axis = [a[0] for a in arr if a[3] == "axis"][0]
        other = sorted(a[0] for a in arr if a[3] != "axis")
        rows.append((f, len(arr), axis, other[0] if other else float("nan"),
                     (other[0] - axis) if other else float("nan")))
    return Xc, rows


def check_T5():
    """dt/db along the ray (Fermat) against a central difference in b."""
    rows = []
    for (z1, z2, X, b, k, th) in ((0.0, 3.0, 2.0, 0.5, 1.0, 0.3), (0.2, 0.9, 12.0, 0.2, 0.8, 1.0)):
        t, h = first_arrival(z1, z2, X, b, k, th)
        e = 1e-5
        fd = (first_arrival(z1, z2, X, b + e, k, th)[0] - first_arrival(z1, z2, X, b - e, k, th)[0]) / (2 * e)
        rows.append((z1, z2, X, b, k, th, h, fd, abs(h - fd) / abs(fd)))
    return rows


def main():
    lines = ["F6.8f-dev step 0: exact first arrival in n = exp(b cos(k z + th)) (development, decides nothing)", ""]
    lines.append("T1 b -> 0 (b = 1e-9, k = 0.8; first order ~1e-9): z1, z2, X, t, flat, |diff|")
    for r in check_T1():
        lines.append("   %6.2f %6.2f %6.2f  %.15f  %.15f  %.1e" % r)
    lines.append("T2 phase-averaged second-order coefficient (16 phases, b = +-1e-4) against t2_mean of F6.8c:"
                 " d, angle, num, closed form, rel")
    for r in check_T2():
        lines.append("   %6.1f %5.2f  %+.9e  %+.9e  %.1e" % r)
    lines.append("T3 finite b against the polyline minimisation of F6.8c (n 200, 400, Richardson):"
                 " b, k, d, angle, th, exact, poly400, Richardson, rel")
    for r in check_T3():
        lines.append("   %.1f %.1f %.1f %.1f %.1f  %.12f  %.12f  %.12f  %.1e" % r)
    Xc, rows = check_T4()
    lines.append("T4 first caustic from the layer of largest index, b = 1e-2, k = 1: pi/sqrt(b) = %.6f;"
                 " X/Xc, number of arrivals, axis time, earliest other, other - axis" % Xc)
    for r in rows:
        lines.append("   %.2f  %d  %.12f  %.12f  %+.3e" % r)
    lines.append("T5 dt/db along the ray against a central difference: z1, z2, X, b, k, th, ray, difference, rel")
    for r in check_T5():
        lines.append("   %.1f %.1f %.1f %.1f %.1f %.1f  %+.10e  %+.10e  %.1e" % r)
    txt = "\n".join(lines)
    print(txt)
    os.makedirs(os.path.join(os.path.dirname(os.path.abspath(__file__)), "out"), exist_ok=True)
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "f68f_step0.txt"), "w") as fh:
        fh.write(txt + "\n")


if __name__ == "__main__":
    main()
