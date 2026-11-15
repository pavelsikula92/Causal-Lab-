"""
F6.8d -- Where the mixed sector departs from covariance: the age term of the mixed (omega x nu) sector in 3+1
============================================================================================================

Why this test exists (15 Nov 2026; the author's order "B, then C, then A, then D", 08:30; reformulated by the author's
choice of 13:0x, "Přeformulovat")
-----------------------------------------------------------------------------------------------------------
F6.8b2 found two departures of the mixed sector from the conformal one that a covariant quadratic action forbids:
Delta = K_loc,ab - K_loc ~ 0.46 k^2 and c_T,ab ~ 0.195 against c_T ~ 0.110 (ratio 1.78), both read from linear fits over
three ages at kT = 8-16. The planned question was "from the age cut or from the action?". The preparation
(docs/derivations/f68d_sectors_beyond_leading_order.md) found:
  * any covariant quadratic form R + R F0 R + C F2 C, with any age-dependent F0, gives B/8 = A/6 (the Weyl tensor of the
    test waves contains nu only), so the ratio 1.78 means the age part is not of covariant form;
  * the explicit boundary term of the b-dependent age cut is <= 2.4e-3 in K_hat units and carries neither departure
    (the record of F6.8b2 refitted without it);
  * the age term is the accumulation of long, nearly null pairs; no covariant weight suppresses them, so any finite age
    picks a frame (the rest frame of the field's history);
  * development rows: the angle-integrated mixed bulk density grows as r^2 (local exponent 1.50, 1.87, 1.96 at
    r = 200, 400, 800, k = 0.2; double precision), and at fixed phase advance p = 0 the long-double density grows as r^3
    (exponents 2.87, 2.97, 2.99, 3.00 at r = 400 ... 3200, stable between levels to 1e-5): the mixed age term grows as
    (kT)^3, not linearly. The linear reading of F6.8b2 is pre-asymptotic (its verdict stands, R132). With F6.8e (optical
    density ~ r^4, age term ~ (kT)^5) the pattern is: each factor of the optical index brings (kL)^2.

The question decided here: does the age term of the mixed sector grow as (kT)^3, i.e. is the departure of the mixed
sector from the conformal one at large ages the accumulation of long pairs focused by the optical index (age-made),
rather than a linear term with a non-covariant coefficient?

Prediction (stated before this computation)
-------------------------------------------
The first-order volume change of a long interval by the optical index has a tau-independent part growing as (kL)^2
(first-order focusing). The BD moments remove its tau^0, tau^2, tau^4 parts; the surviving tau^6 part carries the same
factor. Hence the bulk pair density at fixed p, I_ab(r, p) ~ k^5 (kr)^3 psi(p), the angle-integrated density
rho_ab ~ r^2, and the mixed age part of K_hat_ab grows as (kT)^3: the local exponent alpha_ab(n) = 1 + d ln rho_ab / d ln n
tends to 3. Not decided here, reported: Delta at fixed n (the T-independent departure cannot be separated from the
age part at small n, where the cut is within a wavelength of x; F6.8b2 development: Delta = 0.006 at k = 0.1, kT = 2).

What is computed
----------------
As F6.8e for the mixed sector: the s-integrated bulk pair density of F6.8b2 (_pair_terms, full s range) with the
first-order spheroid pieces Dm, Qm computed in long double with sigma - 1 formed from tau^2 (F6.8b2 forms it as
dt / L - 1, which loses digits for long thin intervals) and tip-clustered axial nodes; k = 0.15 (profile) at
n = kT = 96, 192, 384, 768, 1536, 3072 (r = n / k); p on Gauss nodes in [0, 12] (10) and [12, 24] (6), doubled by the
symmetry p -> -p; rho_ab = (2 pi / k r) int dp I_ab; F'_ab(n) = -a_4 rho_ab / (2 k^5) (the mixed age part of K_hat_ab
per unit n in units of k^2). The development reached n = 640 (k = 0.2, p = 0 only); the judged pair n 1536 -> 3072 at
k = 0.15 was not computed before this test.

Instrument (INSTRUMENT_FAILED if any fails)
-------------------------------------------
  I0 long double eps <= 1e-18 (Linux; numpy on Windows has no 80-bit long double);
  I1 the long-double pieces Dm, Qm against F6.8b2.spheroid_pieces (converged) at three diamonds (k 0.2: L 2, CO 0.4,
     tau 1.2; L 20, CO 0 and 0.3, tau 1): relative difference <= 1e-9;
  I2 two quadrature levels at p = 0, smallest and largest n: relative difference <= 1e-3;
  I3 the conformal control on the same nodes: c_T within 3 % of 0.1035;
  I4 tail |I_ab(r, 24)| / |I_ab(r, 0)| <= 1e-3 at the largest n;
  I5 the form k^2 F(kT): I_ab(0.15, n/0.15) / I_ab(0.25, n/0.25) = (0.15/0.25)^5 at p = 0 for n = 192 and 1536, <= 1e-3.
Status (alpha_ab between n = 1536 and 3072 at k = 0.15)
-------------------------------------------------------
  PASSED     |alpha_ab - 3| <= 0.1: the mixed age term grows as (kT)^3; the departure of the mixed sector from the
             conformal one at large ages is age-made (the focusing of long pairs), not a linear non-covariant term;
  FAILED     |alpha_ab - 3| > 0.5 (in particular alpha_ab near 1: a linear age term as read in F6.8b2);
  UNDECIDED  otherwise.
Reported, not judged: alpha_ab at every pair of n; c3 = F'_ab / (3 n^2) at the largest n; the ratio of the mixed to the
conformal age slope at each n (F6.8b2 read 1.78 at n = 8-16).

What a pass means and what it does not
--------------------------------------
It means the mixed sector's age dependence is not linear: the "c_T,ab = 0.195, ratio 1.78" of F6.8b2 is a pre-asymptotic
reading of a cubic growth, and the planned question "cut or action" has the answer "the age": the long pairs of a field
with a finite history in its rest frame. It does not settle whether the T-independent part (Delta) is covariant, and
it is the bulk of the past only (the explicit boundary term was small at kT <= 16; not computed at large n).
"""

from __future__ import annotations

import math
import time

import numpy as np

from core.stats.parallel import parallel_map
from core.testing import Result
from phases.F6_order_dynamics.F6_08a_bd_action_conformal import _A4, _L2, _gl, kernel_derivs, mean_exp_diamond
from phases.F6_order_dynamics.F6_08b2_bd_action_mixed_rev import _s_nodes
from phases.F6_order_dynamics.F6_08e_bd_action_optical_asymptotic import LD, PI, _Ef, _te_nodes

TEST = {
    "id": "F6.8d",
    "phase": "F6",
    "title": "Age term of the mixed (omega x nu) sector of the expected BD action: linear, or (kT)^3 from the focusing of "
             "long pairs?",
    "theory_part": "where the static sectors depart from covariance beyond leading order; the reading of F6.8b2",
    "question": "Does the age term of the mixed sector of the expected discrete BD action grow as (kT)^3 for kT -> infinity, "
                "so that its departure from the conformal sector is made by the age (long pairs focused by the optical "
                "index) and is not a linear term with a non-covariant coefficient?",
    "method": "the bulk mixed pair density of F6.8b2 with long-double first-order spheroid pieces (sigma - 1 from tau^2, "
              "tip-clustered nodes); k = 0.15, n = kT = 96 ... 3072; profile in the phase advance p; local exponent",
    "prediction": "alpha_ab -> 3 (rho_ab ~ r^2), from the first-order (kL)^2 focusing that survives the BD moments",
    "criterion": "INSTRUMENT: I0 long double; I1 pieces vs F6.8b2 <= 1e-9 rel.; I2 levels <= 1e-3; I3 conformal c_T within "
                 "3 % of 0.1035; I4 tail <= 1e-3; I5 k^5 scaling <= 1e-3. alpha_ab (n 1536 -> 3072, k 0.15): PASSED "
                 "|alpha_ab - 3| <= 0.1; FAILED > 0.5; UNDECIDED otherwise.",
    "null_model": "a linear mixed age term (alpha_ab = 1), as read by F6.8b2",
    "level": "chain",
    "meters": [],
    "depends_on": ["F6.8a", "F6.8b2", "F6.8e"],
    "params": {"ks": [0.15, 0.25], "profile_ks": [0.15], "ns": [96.0, 192.0, 384.0, 768.0, 1536.0, 3072.0],
               "p_panels": [[0.0, 12.0, 10], [12.0, 24.0, 6]], "lev1": [16, 16, 16, 16], "lev2": [24, 24, 24, 16],
               "scale_ns": [192.0, 1536.0], "c_T_ref": 0.1035,
               "tol_pair": 1e-9, "tol_level": 1e-3, "tol_cT": 0.03, "tol_tail": 1e-3, "tol_scale": 1e-3,
               "tol_pass": 0.1, "tol_fail": 0.5, "alpha_pred": 3.0},
    "parts": {"main": {}},
    "quick": {"ks": [0.15], "profile_ks": [0.15], "ns": [48.0, 96.0], "p_panels": [[0.0, 12.0, 4], [12.0, 24.0, 2]],
              "scale_ns": [], "lev1": [16, 8, 8, 8], "lev2": [16, 16, 8, 8]},
    "estimate": "~125 pair-density evaluations of ~30-60 s (long double): ~45 min with 2 processes; Linux only",
}


def pieces_ld(k, CO, L, tau2s, n_sig=16, n_te=16, n_phi=16, npan=None):
    """Dm, Qm of F6.8b2.spheroid_pieces for one L and the proper times tau2s, in long double; and D_xy."""
    k, CO, L = LD(k), LD(CO), LD(L)
    tau2 = np.asarray(tau2s, LD)[:, None, None]
    dt = np.sqrt(L * L + tau2)
    kpar = k * CO
    kperp = k * np.sqrt(max(LD(1) - CO * CO, LD(0)))
    npan = npan or int(8 + math.ceil(float(abs(kpar) * L)))
    c = L / 2
    em1 = tau2 / (L * (dt + L))                         # sigma_max - 1, from tau^2
    E_xy, _ = _Ef(kpar * L)
    Dxy = L * E_xy
    LDx = L / dt * Dxy
    xs, ws = _gl(0.0, 1.0, n_sig)
    te, wt = _te_nodes(n_te, npan)
    te, wt, xs, ws = (np.asarray(v, LD) for v in (te, wt, xs, ws))
    ph = 2 * PI * (np.arange(n_phi, dtype=LD) + LD(0.5)) / n_phi
    TE, PH = np.meshgrid(te, ph, indexing="ij")
    TE, PH = TE[None], PH[None]
    WT = ((wt[:, None] * np.ones_like(PH[0])) * (2 * PI / n_phi))[None]
    cph = np.cos(PH)
    eky = np.exp(1j * kpar * L)
    numD = np.zeros(tau2.shape[0], np.clongdouble)
    numQ = np.zeros(tau2.shape[0], np.clongdouble)
    den = np.zeros(tau2.shape[0], LD)
    for x, wx in zip(xs, ws):
        sm1 = em1 * x
        sig = 1 + sm1
        lw = c * (1 + sig * TE)                         # axial distance from x
        rp = c * np.sqrt(np.maximum(sm1 * (2 + sm1) * (1 - TE * TE), LD(0)))
        d1 = c * (sig - TE)                             # to y
        d2 = c * (sig + TE)                             # to x
        kw = kpar * lw + kperp * rp * cph               # k . (z - x)
        E1, _ = _Ef(kw - kpar * L)
        E2, _ = _Ef(kw)
        dlt = LDx - (d1 * eky * E1 + d2 * E2)
        wgt = (em1 * wx) * WT * c ** 3 * (sig * sig - TE * TE)
        numD += (wgt * dlt).sum(axis=(1, 2))
        numQ += (wgt * np.exp(1j * kw) * np.conj(dlt)).sum(axis=(1, 2))
        den += wgt.sum(axis=(1, 2))
    return (numD / den).astype(complex), (numQ / den).astype(complex), complex(Dxy)


def pair_density(k, CO, r, n_s=16, lev=(16, 16, 16)):
    """s-integrated bulk pair density at distance r, direction CO (full s range): mixed and conformal."""
    s, ws = _s_nodes(np.array([1e12]), n_s, True)
    s, ws = s[0], ws[0]
    keep = ws > 0
    s, ws = s[keep], ws[keep]
    tau2 = s / _L2
    dt = np.sqrt(r * r + tau2)
    g0, g1, g2 = kernel_derivs(s * s)
    ph = np.exp(1j * k * CO * r)
    one = 1.0 + ph
    M = mean_exp_diamond(k, np.full((1, s.size), r), CO, dt[None, :])[0]
    Dm, Qm, Dxy = pieces_ld(k, CO, r, tau2, *lev)
    AV = 4.0 * dt / tau2
    alpha = 4.0 * M
    beta = 3.0 * M + AV * Dm
    gam = 6.0 + 2.0 * AV * np.real(Qm)
    mu = -r * Dxy / (dt * dt)
    mixed = (12.0 * (1.0 + ph.real) * g0 + 2.0 * g1 * np.real(one * np.conj(beta))
             + 1.5 * g1 * np.real(alpha * np.conj(one)) + 0.5 * g2 * np.real(alpha * np.conj(beta)) + g1 * gam
             + 0.5 * np.real((4.0 * one * g0 + g1 * alpha) * np.conj(mu)))
    conf = (1.0 + ph.real) * g0 + np.real(one * np.conj(M)) * g1 + 0.5 * np.abs(M) ** 2 * g2 + 0.5 * g1
    base = ws * r * r / (_L2 * 2.0 * dt)
    return float((base * mixed).sum()), float((base * conf).sum())


def _p_nodes(panels):
    ps, wp = [], []
    for a, b, n in panels:
        x, w = _gl(float(a), float(b), int(n))
        ps.append(x)
        wp.append(w)
    return np.concatenate(ps), np.concatenate(wp)


def _dispatch(kind, args):
    if kind == "dens":
        k, r, p, lev = args
        b, c = pair_density(k, p / (k * r), r, lev[0], (lev[1], lev[2], lev[3]))
        return {"k": k, "r": r, "p": p, "ab": b, "conf": c}
    if kind == "pair":
        from phases.F6_order_dynamics.F6_08b2_bd_action_mixed_rev import spheroid_pieces
        k, CO, L, tau = args
        dt = math.sqrt(L * L + tau * tau)
        Dm, Qm, _ = spheroid_pieces(k, CO, np.array([[L]]), np.array([[dt]]), 64, 256, 64)
        D2, Q2, _ = pieces_ld(k, CO, L, [tau * tau], 32, 32, 32)
        rel = max(abs(complex(Dm[0, 0]) - D2[0]) / abs(D2[0]), abs(complex(Qm[0, 0]) - Q2[0]) / abs(Q2[0]))
        return {"k": k, "CO": CO, "L": L, "tau": tau, "rel": rel}
    raise ValueError(kind)


def run(p, ctx):
    t0 = time.time()
    lines = ["3+1, eps rho = 1 (l_k = 1); mixed (omega x nu) bulk pair density of F6.8b2 with long-double first-order "
             "pieces; p = k r cos(n, k)"]
    eps = float(np.finfo(np.longdouble).eps)
    lines.append(f"  I0 long double eps {eps:.2e} (need <= 1e-18)")
    if eps > 1e-18:
        lines.append("  no 80-bit long double on this platform (numpy on Windows): run under Linux")
        return Result("INSTRUMENT_FAILED", lines, {"eps": eps})
    ks = [float(k) for k in p["ks"]]
    pks = [float(k) for k in p.get("profile_ks", ks)]
    ns = [float(n) for n in p["ns"]]
    pn, wp = _p_nodes(p["p_panels"])
    lev1, lev2 = tuple(p["lev1"]), tuple(p["lev2"])
    jobs = [("dens", (k, n / k, float(pp), lev1)) for k in pks for n in ns for pp in pn]
    nm = len(jobs)
    zn = sorted({ns[0], ns[-1]} | {float(n) for n in p.get("scale_ns", [])})
    zero = [("dens", (k, n / k, 0.0, lev1)) for k in ks for n in zn]
    tail = [("dens", (k, ns[-1] / k, 24.0, lev1)) for k in pks]
    lv = [("dens", (k, n / k, 0.0, lev2)) for k in pks for n in sorted({ns[0], ns[-1]})]
    pairs = [("pair", (0.2, 0.4, 2.0, 1.2)), ("pair", (0.2, 0.0, 20.0, 1.0)), ("pair", (0.2, 0.3, 20.0, 1.0))]
    res = parallel_map(_dispatch, jobs + zero + tail + lv + pairs, ctx.workers, ctx.progress,
                       lambda i, r_, d, n: f"[F6.8d] task {d}/{n} ({time.time() - t0:.0f}s)")
    i = 0
    main = res[i:i + nm]; i += nm
    rz = res[i:i + len(zero)]; i += len(zero)
    rt = res[i:i + len(tail)]; i += len(tail)
    rl = res[i:i + len(lv)]; i += len(lv)
    rp = res[i:]

    def find(rr, k, r, pp):
        for x in rr:
            if abs(x["k"] - k) < 1e-12 and abs(x["r"] - r) < 1e-9 * r and abs(x["p"] - pp) < 1e-12:
                return x
        return None

    ok1 = all(x["rel"] <= float(p["tol_pair"]) for x in rp)
    lines.append("  I1 long-double pieces vs F6.8b2: " + "; ".join(f"L {x['L']:g} CO {x['CO']:g}: rel {x['rel']:.1e}"
                                                                for x in rp) + f" (tol {p['tol_pair']:g})")
    ok2 = True
    for x in rl:
        y = find(rz, x["k"], x["r"], 0.0)
        d = abs(x["ab"] - y["ab"]) / abs(x["ab"])
        ok2 &= d <= float(p["tol_level"])
        lines.append(f"  I2 k {x['k']:g} n {x['r'] * x['k']:g}: I_ab(r, 0) level 1 {y['ab']:+.6e}, level 2 {x['ab']:+.6e} "
                     f"(rel {d:.1e}, tol {p['tol_level']:g})")
    prof = {}
    for k in pks:
        for n in ns:
            r = n / k
            vb = np.array([find(main, k, r, float(pp))["ab"] for pp in pn])
            vc = np.array([find(main, k, r, float(pp))["conf"] for pp in pn])
            pref = 2.0 * math.pi / (k * r)
            prof[(k, n)] = {"rho_ab": pref * 2.0 * float((wp * vb).sum()), "rho_c": pref * 2.0 * float((wp * vc).sum())}
    ok3 = True
    for k in pks:
        cT = -(16.0 / 3.0) * _A4 * prof[(k, ns[-1])]["rho_c"] / k ** 5
        d = abs(cT - float(p["c_T_ref"])) / float(p["c_T_ref"])
        ok3 &= d <= float(p["tol_cT"])
        lines.append(f"  I3 conformal c_T at k {k:g}, n {ns[-1]:g}: {cT:.4f} (F6.8a 0.1035; rel {d:.1e}, tol {p['tol_cT']:g})")
    ok4 = True
    for x in rt:
        y = find(rz, x["k"], x["r"], 0.0)
        q = abs(x["ab"]) / abs(y["ab"])
        ok4 &= q <= float(p["tol_tail"])
        lines.append(f"  I4 tail k {x['k']:g} n {x['r'] * x['k']:g}: |I(24)| / |I(0)| = {q:.1e} (tol {p['tol_tail']:g})")
    ok5 = True
    if len(ks) >= 2:
        for n in p.get("scale_ns", []):
            a_, b_ = find(rz, ks[0], n / ks[0], 0.0), find(rz, ks[1], n / ks[1], 0.0)
            want = (ks[0] / ks[1]) ** 5
            d = abs(a_["ab"] / b_["ab"] / want - 1.0)
            ok5 &= d <= float(p["tol_scale"])
            lines.append(f"  I5 n {n:g}: I_ab(k {ks[0]:g}) / I_ab(k {ks[1]:g}) = {a_['ab'] / b_['ab']:.6e} vs {want:.6e} "
                         f"(rel {d:.1e}, tol {p['tol_scale']:g})")
    inst_ok = ok1 and ok2 and ok3 and ok4 and ok5
    lines.append("  k      n      r          rho_ab          F'_ab(n) = -a4 rho_ab / (2 k^5)   slope ratio mixed/conf   "
                 "alpha_ab")
    alphas = {}
    for k in pks:
        prev = None
        for n in ns:
            pr = prof[(k, n)]
            Fab = -_A4 * pr["rho_ab"] / (2.0 * k ** 5)
            Fc = -(16.0 / 3.0) * _A4 * pr["rho_c"] / k ** 5
            al = ""
            if prev is not None:
                a_ = 1.0 + math.log(pr["rho_ab"] / prev[1]) / math.log(n / prev[0])
                alphas[(k, n)] = a_
                al = f"{a_:.4f}"
            lines.append(f"  {k:<6g} {n:<6g} {n / k:<10.5g} {pr['rho_ab']:+.6e}   {Fab:+.6e}                      "
                         f"{Fab / Fc:+.4f}                   {al}")
            prev = (n, pr["rho_ab"])
        Fab = -_A4 * prof[(k, ns[-1])]["rho_ab"] / (2.0 * k ** 5)
        lines.append(f"  reported k {k:g}: c3 = F'_ab(n) / (3 n^2) at n {ns[-1]:g} = {Fab / (3.0 * ns[-1] ** 2):.4e}")
    if not inst_ok:
        status = "INSTRUMENT_FAILED"
    elif len(ns) < 2 or (not ctx.quick and len(ks) < 2):
        status = "UNDECIDED"
    else:
        dev_ = [abs(alphas[(k, ns[-1])] - float(p["alpha_pred"])) for k in pks]
        if all(d <= float(p["tol_pass"]) for d in dev_):
            status = "PASSED"
        elif any(d > float(p["tol_fail"]) for d in dev_):
            status = "FAILED"
        else:
            status = "UNDECIDED"
    data = {"profiles": {f"{k}_{n}": v for (k, n), v in prof.items()},
            "alphas": {f"{k}_{n}": v for (k, n), v in alphas.items()}, "levels": rl, "zero": rz, "tail": rt, "pairs": rp,
            "p_nodes": pn.tolist(), "p_weights": wp.tolist()}
    return Result(status, lines, data)
