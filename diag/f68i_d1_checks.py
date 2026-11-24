"""
F6.8i-dev, D1 checks 2 and 3 (development, decides nothing): the cut-locus window of the needle volume, integrated over
the phase theta, at finite q and in the age term. docs/derivations/f68i_dev.md section 5 (5.2 for this file).

Near a cut locus theta_c the first-arrival rays a and b tie. On each side of theta_c (delta = |theta - theta_c|) let
"first" be the first arrival and "second" the other ray, d(delta) = S_second - S_first. The two-well needle volume
(check 1, diag/f68i_d1.py) is

    A(delta; c) = sum_j A_j^first(delta) (2c)^j + [c > d] sum_j A_j^second(delta) 2^j (c - d)^(j + 2) / c^2,   c = q^2/2.

Model. On each side the rays, d and the series coefficients A_0..A_3 of both rays are computed at --nodes values of
delta in (0, delta_max] and fitted by polynomials in delta (d without a constant term). Model checks: the fit residual
at the nodes, a held-out node, and continuity across theta_c (first on one side = second on the other at delta = 0).

Check 2 (finite q). The window integral (1/2pi) int dtheta [A_exact - A_first-series] over both sides, with the exact
needle (diag/f68f_1b_exact: lower envelopes over all rays, no ray chosen) at Gauss-Legendre nodes inside the window
d < c, against
  - the leading term dA = (A_a + A_b) c / (6 pi alpha) of f68f_dev section 6 (the relative difference must be O(c)),
  - the same integral of the two-well model (the difference is the truncation of the series, check 1).

Check 3 (the age term). In the age term c = gamma s (gamma = k^2 / (2 l2)), rho = (T/(2 l2)) int ds <g0(s^2 A)>_theta.
  3a: for one ray (a phase near theta_c) F(gamma) = int ds g0(s^2 A(gamma s)) / (4 mu_3 (2 gamma)^3) must tend to
      s4hat A_0^-2 of diag/f68f_1c (the moment selection of the smooth part, with this file's quadrature).
  3b: the window part W(gamma) = int ds (1/2pi) sum_sides int ddelta [g0(s^2 A_two-well) - g0(s^2 A_first)] must have
      no gamma^1 and gamma^2 terms (BD moments mu_1 = mu_2 = 0, f68i_dev 5): W / gamma^3 tends to a constant W3. Its
      contribution to <s4hat A_0^-2>_theta is W3 / (4 mu_3 8) per cut locus (reported with the smooth part of the same
      (eps, p) on equal phases for scale only: equal phases across the jump at theta_c converge as 1/N).

Cut loci: --points eps:p:theta_c (theta_c optional; default the symmetric tie theta_c = -p/2 mod 2 pi, refined).

    python -m diag.f68i_d1_checks --points 12:1,16:1,12:0 --workers 8
Output diag/out/f68i_d1_checks_<UTC time>.json and .txt.
"""
from __future__ import annotations

import argparse
import math
import os
import time
from multiprocessing import Pool

import numpy as np

from diag import f68f_1b as B
from diag import f68f_1b_cut as CUT
from diag import f68f_1b_exact as EX
from diag import f68f_1c as SER
from diag import f68f_common as C
from diag.f68i_d1 import series_for_ray
from phases.F6_order_dynamics.F6_08a_bd_action_conformal import kernel_derivs

MU3 = 0.5                     # int s^3 g0(s^2) ds (f68f_dev section 3)
QS = EX.QS                    # 1/2, 1/4, 1/8
C_MAX = max(QS) ** 2 / 2      # the largest c the model must cover
GAMMAS = [0.0008, 0.0011, 0.0016, 0.0023, 0.0032, 0.0045, 0.0064, 0.009]   # check 3b (the window part)
GAMMAS_3A = [0.005, 0.007, 0.01, 0.014, 0.019, 0.027, 0.038]   # check 3a: F is O(gamma^3) of O(1) pieces; small gamma loses digits
DEG_D, DEG_A = 5, 5


# ------------------------------------------------------------------------------------------------ jobs (in workers)

def node_job(args):
    """Rays at one phase: d to the second and third ray, series A_0..A_3 of the first and second ray."""
    eps, p, th = args
    vv, S = CUT.rays(eps, p, th)
    o = np.argsort(S)
    gaps = S[o] - S[o[0]]
    return {"th": th, "d": float(gaps[1]), "d3": float(gaps[2]) if gaps.size > 2 else float("inf"),
            "A1": list(map(float, series_for_ray(eps, p, th, vv[o[0]]))),
            "A2": list(map(float, series_for_ray(eps, p, th, vv[o[1]]))), "nray": int(vv.size)}


def exact_job(args):
    eps, p, th = args
    r = EX.phase_value(eps, p, th)
    return {"th": th, "A": [float(a) for a in r["A"]], "minPhi": float(r["minPhi"])}


def smooth_job(args):
    eps, p, n_th = args
    r = SER.point(eps, p, n_th)
    return {"phi": r["phi"], "nray": r["nray"], "minK": r["minK"]}


# ------------------------------------------------------------------------------------------------ the model

def refine_tie(eps, p, th0, h=1e-4):
    """theta where the two lowest rays tie: secant on the signed gap, continuing the two rays."""
    vv, S = CUT.rays(eps, p, th0)
    o = np.argsort(S)
    a, b = vv[o[0]], vv[o[1]]

    def f(th):
        return CUT.S_of(CUT.follow(b, eps, p, th), eps, p, th) - CUT.S_of(CUT.follow(a, eps, p, th), eps, p, th)

    x0, x1 = th0 - h, th0 + h
    f0, f1 = f(x0), f(x1)
    for _ in range(30):
        if f1 == f0:
            break
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        x0, f0, x1, f1 = x1, f1, x2, f(x2)
        if abs(x1 - x0) < 1e-13:
            break
    return float(x1)


def fit_side(rows, deg_d=DEG_D, deg_a=DEG_A):
    dl = np.array([r["delta"] for r in rows])
    d = np.array([r["d"] for r in rows])
    V = np.vander(dl, deg_d + 1, increasing=True)[:, 1:]
    cd = np.linalg.lstsq(V, d, rcond=None)[0]
    pd = np.polynomial.Polynomial(np.concatenate([[0.0], cd]))
    P1, P2 = [], []
    for key, out in (("A1", P1), ("A2", P2)):
        Y = np.array([r[key] for r in rows])
        for j in range(4):
            out.append(np.polynomial.Polynomial.fit(dl, Y[:, j], deg_a, domain=[0.0, dl.max()]).convert())
    res_d = float(np.max(np.abs(pd(dl) - d)))
    res_a = float(max(np.max(np.abs(P(dl) - np.array([r[k] for r in rows])[:, j]))
                      for k, PP in (("A1", P1), ("A2", P2)) for j, P in enumerate(PP)))
    return {"d": pd, "dd": pd.deriv(), "P1": P1, "P2": P2, "res_d": res_d, "res_a": res_a,
            "alpha": float(cd[0]), "dmax": float(dl.max())}


def delta_of(M, target, alpha):
    """Solve d(delta) = target (array) on one side by Newton from target / alpha."""
    x = np.asarray(target, float) / alpha
    for _ in range(60):
        f = M["d"](x) - target
        x = x - f / M["dd"](x)
        if np.max(np.abs(f)) < 1e-15:
            break
    return x


def g0(w):
    return kernel_derivs(w)[0]


# ------------------------------------------------------------------------------------------------ check 3

def window_W(sides, gam, ns=240, nz=32, absval=False):
    """(1/2pi) sum_sides int ds int ddelta [g0(s^2 A2w) - g0(s^2 A1)] at c = gam s."""
    xs, ws = np.polynomial.legendre.leggauss(ns)
    zx, zw = np.polynomial.legendre.leggauss(nz)
    z = 0.5 * (zx + 1.0)
    wz = 0.5 * zw
    total = 0.0
    for M in sides:
        amin = min(float(M["P1"][0](0.0)), float(M["P1"][0](M["dmax"])))
        smax = math.sqrt(80.0 / max(amin, 0.2))
        smax = min(smax, C_MAX / gam)
        s = 0.5 * smax * (xs + 1.0)
        wss = 0.5 * smax * ws
        c = gam * s[:, None]
        dl = delta_of(M, c * z[None, :], M["alpha"])
        if np.any(dl > M["dmax"] * 1.0001):
            raise ValueError("window beyond the model range")
        A1 = sum(M["P1"][j](dl) * (2 * c) ** j for j in range(4))
        X = sum(M["P2"][j](dl) * 2 ** j * c ** j * (1 - z[None, :]) ** (j + 2) for j in range(4))
        jac = c / M["dd"](dl)
        S2 = (s * s)[:, None]
        f = (g0(S2 * (A1 + X)) - g0(S2 * A1)) * jac
        if absval:
            f = np.abs(f)
        total += float((wss[:, None] * wz[None, :] * f).sum())
    return total / (2 * math.pi)


def smooth_F(A, gam, ns=240):
    """int ds g0(s^2 A(gam s)) / (4 mu_3 (2 gam)^3) for the series A_0..A_3 of one ray."""
    xs, ws = np.polynomial.legendre.leggauss(ns)
    smax = math.sqrt(80.0 / A[0])
    s = 0.5 * smax * (xs + 1.0)
    w = 0.5 * smax * ws
    Aq = sum(A[j] * (2 * gam * s) ** j for j in range(4))
    return float((w * g0(s * s * Aq)).sum()) / (4 * MU3 * (2 * gam) ** 3)


# ------------------------------------------------------------------------------------------------ one cut locus

def run_point(pool, eps, p, thc0, n_nodes, n_gl, n_smooth, lines):
    t0 = time.time()
    thc = refine_tie(eps, p, thc0)
    probe = pool.map(node_job, [(eps, p, thc - 1e-3), (eps, p, thc + 1e-3)])
    alpha0 = 0.5 * (probe[0]["d"] + probe[1]["d"]) / 1e-3
    dmax = 1.25 * C_MAX / alpha0
    u = 0.5 * (1 - np.cos(math.pi * (np.arange(n_nodes) + 0.5) / n_nodes))     # Chebyshev points in (0, 1)
    deltas = list(dmax * u / u.max())
    held = float(dmax * 0.5 * (u[n_nodes // 2] + u[n_nodes // 2 + 1]) / u.max())  # a held-out node between two
    jobs = [(eps, p, thc + sg * dl) for sg in (-1, 1) for dl in deltas + [held]]
    res = pool.map(node_job, jobs, chunksize=1)
    sides, raw = [], {}
    k = 0
    for sg in (-1, 1):
        rows = []
        for dl in deltas + [held]:
            r = dict(res[k])
            r["delta"] = dl
            rows.append(r)
            k += 1
        fit_rows, held_row = rows[:-1], rows[-1]
        M = fit_side(fit_rows)
        M["held_d"] = abs(float(M["d"](held_row["delta"])) - held_row["d"])
        M["held_a"] = max(abs(float(M["P" + str(i)][j](held_row["delta"])) - held_row["A" + str(i)][j])
                          for i in (1, 2) for j in range(4))
        M["min_d3"] = min(r["d3"] for r in rows)
        M["nray"] = max(r["nray"] for r in rows)
        sides.append(M)
        raw[sg] = rows
    cont = max(abs(float(sides[0]["P1"][j](0.0) - sides[1]["P2"][j](0.0))) for j in range(4))
    cont = max(cont, max(abs(float(sides[1]["P1"][j](0.0) - sides[0]["P2"][j](0.0))) for j in range(4)))
    Aa, Ab = float(sides[1]["P2"][0](0.0)), float(sides[1]["P1"][0](0.0))
    alpha = 0.5 * (sides[0]["alpha"] + sides[1]["alpha"])
    lines.append(f"eps {eps:g} p {p:g}: theta_c {thc:.9f}  alpha {sides[0]['alpha']:.5f} / {sides[1]['alpha']:.5f} "
                 f"(left / right)  A_a {Aa:.6f}  A_b {Ab:.6f}  rays {max(M['nray'] for M in sides)}  "
                 f"min gap to the third ray {min(M['min_d3'] for M in sides):.3f} (c_max {C_MAX})")
    lines.append(f"  model: {n_nodes} nodes per side in (0, {dmax:.4f}]; fit residual d {max(M['res_d'] for M in sides):.1e},"
                 f" A {max(M['res_a'] for M in sides):.1e}; held-out node d {max(M['held_d'] for M in sides):.1e}, "
                 f"A {max(M['held_a'] for M in sides):.1e}; continuity across theta_c {cont:.1e}")

    # ---- check 2: exact needle inside the window
    xg, wg = np.polynomial.legendre.leggauss(n_gl)
    ex_jobs, meta = [], []
    for iq, q in enumerate(QS):
        c = q * q / 2
        for si, sg in enumerate((-1, 1)):
            M = sides[si]
            de = float(delta_of(M, np.array([c]), M["alpha"])[0])
            for x, w in zip(xg, wg):
                dl = 0.5 * de * (x + 1.0)
                ex_jobs.append((eps, p, thc + sg * dl))
                meta.append((iq, si, dl, 0.5 * de * w))
    ex = pool.map(exact_job, ex_jobs, chunksize=1)
    I_ex = np.zeros(len(QS))
    I_mo = np.zeros(len(QS))
    minphi = min(e["minPhi"] for e in ex)
    for (iq, si, dl, w), e in zip(meta, ex):
        q = QS[iq]
        c = q * q / 2
        M = sides[si]
        first = sum(float(M["P1"][j](dl)) * (2 * c) ** j for j in range(4))
        d = float(M["d"](dl))
        second = sum(float(M["P2"][j](dl)) * 2 ** j * max(c - d, 0.0) ** (j + 2) / c ** 2 for j in range(4))
        I_ex[iq] += w * (e["A"][iq] - first)
        I_mo[iq] += w * second
    I_ex /= 2 * math.pi
    I_mo /= 2 * math.pi
    dA = np.array([(Aa + Ab) * (q * q / 2) / (6 * math.pi * alpha) for q in QS])
    lines.append(f"  check 2 (window integral over theta, exact needle at {n_gl} GL nodes per side and q; min Phi {minphi:.1e}):")
    for iq, q in enumerate(QS):
        lines.append(f"    q {q:<6g} c {q * q / 2:.5f}: exact {I_ex[iq]:.6e}  two-well {I_mo[iq]:.6e} "
                     f"(diff {I_ex[iq] - I_mo[iq]:+.1e})  dA {dA[iq]:.6e}  exact/dA - 1 {I_ex[iq] / dA[iq] - 1:+.4e}")
    r = (I_ex / dA - 1)
    lines.append(f"    (exact/dA - 1) falls by {r[0] / r[1]:.2f} and {r[1] / r[2]:.2f} per halving of q "
                 f"(O(c): 4)")

    # ---- check 3a: the moment selection of one ray
    ref = raw[1][n_nodes // 2]
    A = np.array(ref["A1"])
    target = SER.s4hat(A) * A[0] ** -2.0
    Fs = [smooth_F(A, g) for g in GAMMAS_3A]
    lines.append(f"  check 3a (one ray at delta {ref['delta']:.4f}): s4hat A0^-2 = {target:.6e};  F(gamma) = " +
                 ", ".join(f"{f:.6e}" for f in Fs[:3]) + f" at gamma {GAMMAS_3A[0]}, {GAMMAS_3A[1]}, {GAMMAS_3A[2]}")
    lin = np.polyfit(GAMMAS_3A, Fs, 3)
    lines.append(f"    F(gamma -> 0) by a cubic fit over gamma {GAMMAS_3A[0]}..{GAMMAS_3A[-1]}: {lin[-1]:.6e} "
                 f"(relative to s4hat A0^-2 {lin[-1] / target - 1:+.1e}; limited by rounding where s4hat A0^-2 is small)")

    # ---- check 3b: the window term in the age term
    W = np.array([window_W(sides, g) for g in GAMMAS])
    N = np.array([window_W(sides, g, absval=True) for g in GAMMAS])
    W2x = window_W(sides, GAMMAS[2], ns=480, nz=64)
    G = np.array(GAMMAS)

    def fit(powers):
        V = np.stack([G ** m for m in powers], 1)
        c = np.linalg.lstsq(V, W, rcond=None)[0]
        return c, float(np.max(np.abs(V @ c - W) / np.abs(W)))

    coef, res_all = fit((1, 2, 3, 4, 5, 6))
    coef3, res3 = fit((3, 4, 5, 6))
    N1 = float(np.median(N / G))
    lines.append(f"  check 3b (window part of the age term, c = gamma s; N = the same with |.|, the size of the gamma^1 "
                 f"piece that the moments remove):")
    for g, w, n in zip(GAMMAS, W, N):
        lines.append(f"    gamma {g:<7g} W {w:+.6e}  N {n:.4e}  W/N {w / n:+.2e}  W/gamma^3 {w / g ** 3:+.6e}")
    lines.append(f"    fit W = sum_m W_m gamma^m, m 1..6: W1/N1 {coef[0] / N1:+.1e}  W2/N1 {coef[1] / N1:+.1e}  "
                 f"W3 {coef[2]:+.6e}  W4 {coef[3]:+.3e}  (N1 = N/gamma {N1:.4e}; max rel. residual {res_all:.1e})")
    lines.append(f"    fit with m 3..6 only: W3 {coef3[0]:+.6e} (max rel. residual {res3:.1e}); the two W3 differ by "
                 f"{abs(coef3[0] / coef[2] - 1):.1e}")
    lines.append(f"    quadrature: doubling s and z nodes at gamma {GAMMAS[2]} changes W by {abs(W2x - W[2]) / abs(W[2]):.1e} "
                 f"(relative)")
    w_eq = coef[2] / (4 * MU3 * 8)
    out = {"eps": eps, "p": p, "theta_c": thc, "alpha": [sides[0]["alpha"], sides[1]["alpha"]], "A_a": Aa, "A_b": Ab,
           "model": {"nodes": n_nodes, "delta_max": dmax, "res_d": max(M["res_d"] for M in sides),
                     "res_a": max(M["res_a"] for M in sides), "held_d": max(M["held_d"] for M in sides),
                     "held_a": max(M["held_a"] for M in sides), "continuity": cont,
                     "min_gap_third": min(M["min_d3"] for M in sides)},
           "check2": {"q": QS, "exact": I_ex.tolist(), "two_well": I_mo.tolist(), "dA": dA.tolist(), "minPhi": minphi},
           "check3a": {"target": target, "gamma": GAMMAS_3A, "F": Fs, "F0_fit": float(lin[-1])},
           "check3b": {"gamma": GAMMAS, "W": W.tolist(), "N": N.tolist(), "coef": coef.tolist(), "N1": N1,
                       "coef_3to6": coef3.tolist(), "res_all": res_all, "res_3to6": res3,
                       "W_doubled": W2x, "window_s4hat_equiv": w_eq},
           "nodes": {str(k): v for k, v in raw.items()}}
    if n_smooth:
        sm = pool.map(smooth_job, [(eps, p, n_smooth)])[0]
        out["smooth"] = sm
        lines.append(f"  the window term per cut locus in <s4hat A0^-2>_theta: W3 / (4 mu_3 8) = {w_eq:+.4e};  smooth part "
                     f"on {n_smooth} equal phases (scale only, across the jump) {sm['phi']:+.4e}")
    else:
        lines.append(f"  the window term per cut locus in <s4hat A0^-2>_theta: W3 / (4 mu_3 8) = {w_eq:+.4e}")
    lines.append(f"  ({time.time() - t0:.0f} s)")
    return out


NAME = "f68i_d1_checks"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--points", default="12:1,16:1,12:0", help="eps:p[:theta_c],...")
    ap.add_argument("--nodes", type=int, default=9, help="model nodes per side")
    ap.add_argument("--gl", type=int, default=6, help="exact-needle Gauss-Legendre nodes per side and q (check 2)")
    ap.add_argument("--smooth", type=int, default=64, help="equal phases for the smooth part (0: skip)")
    ap.add_argument("--workers", type=int, default=2)
    a = ap.parse_args()
    t0 = time.time()
    code = C.code_md5([os.path.abspath(__file__), B.__file__, CUT.__file__, EX.__file__, SER.__file__,
                       os.path.join(os.path.dirname(os.path.abspath(__file__)), "f68i_d1.py")])
    lines = ["F6.8i-dev D1 checks 2 and 3: the cut-locus window integrated over theta (development, decides nothing; "
             "docs/derivations/f68i_dev.md section 5.2)", ""]
    data = []
    with Pool(a.workers) as pool:
        for item in a.points.split(","):
            f = [float(x) for x in item.split(":")]
            eps, p = f[0], f[1]
            thc0 = f[2] if len(f) > 2 else (-p / 2) % (2 * math.pi)
            data.append(run_point(pool, eps, p, thc0, a.nodes, a.gl, a.smooth, lines))
            print("\n".join(lines[-14:]), flush=True)
    lines.append(f"total {time.time() - t0:.0f} s")
    stem = C.save(NAME, vars(a), code, data, lines, t0)
    print("\n".join(lines))
    print("saved", stem + ".json")


if __name__ == "__main__":
    main()
