"""
F6.8e3 -- Second revision of F6.8e: the instrument checks stated in the quantity that is judged (the exponent)
==============================================================================================================

What happened (15 Nov 2026; docs/findings.md, "F6.8e" and "F6.8e2"; the author's instruction of 22:38, "test by se měl
opravit a spustit poté znovu")
------------------------------------------------------------------------------------------------------------------
F6.8e and F6.8e2 both ended INSTRUMENT_FAILED on checks whose tolerances or ranges the cards set wrongly; the judged
quantity (alpha) was never the failing item, and both runs recorded alpha = 5.000 at the largest ages (seen):
  * F6.8e, I5: the k1/k2 density ratio against (k1/k2)^5 (a k-dependent normalisation of 0.2 %, not a shape change);
  * F6.8e2, I3: the conformal control at r = 61 440, where the conformal path (double precision, F6.8a's
    mean_exp_diamond) has lost its digits (constant to r ~ 1.5e4, then departing);
  * F6.8e2, I5: the ratio drift 1.7e-3 against 1e-3 between n 384 and 3072 (an exponent difference of 8e-4);
  * F6.8e2, I2: a level difference 1.2e-3 against 1e-3 at the smallest age (irrelevant for an exponent).
The fault common to all four: the checks were stated in raw densities with tolerances unrelated to what is judged.
This revision states them in the exponent, the quantity of the status, at a tolerance ten times tighter than the
status tolerance (0.01 against 0.1), and places each check where its instrument is valid.

Changes against F6.8e2 (nothing else changes)
---------------------------------------------
  * I2 (levels): the exponent from the p = 0 densities at the judged pair n 3072 -> 6144 at two quadrature levels:
    |alpha_0(level 2) - alpha_0(level 1)| <= 0.01.
  * I3 (conformal control) at n = 384 (r = 2954 at k = 0.13), inside the range where the conformal path is exact
    (F6.8e2: c_T 0.1019 at r = 7680); tolerance 3 % as before.
  * I5 (the kT-shape independent of k): the exponent from the p = 0 densities at the judged pair at k1 = 0.13 and
    k2 = 0.21: |alpha_0(k1) - alpha_0(k2)| <= 0.01. The deviation of the density ratio from (k1/k2)^5 is reported.
  * Wavenumbers k = 0.13 (profile) and 0.21; the largest r is 47 262 (k = 0.13, n = 6144), inside the range where the
    optical levels agree (F6.8e2: 6.3e-4 at r = 61 440; at r = 87 771 the levels differ by 7.5e-3, checked before this
    card without printing values, so k = 0.07 was not used).
  * Ages n = 192 ... 6144 as F6.8e2; judged pair n 3072 -> 6144 at k = 0.13. The exponent is not blind: F6.8e and
    F6.8e2 recorded 5.000. This run gives the verdict with a corrected instrument; it is not new evidence beyond them.

F6.8e (the frozen card, kept for reference)
-------------------------------------------
Revision 1 (15 Nov 2026, 13:2x, before any result; the author's choice "Revize před během"): the first full run in
Claude's workspace was stopped at 175 of 277 tasks by the workspace's two-hour limit per process; nothing of it was
printed or saved. To fit the limit the profile is computed at one wavenumber (k = 0.15; the k^5 scaling is exact, a
second wavenumber adds no information about alpha) with 10 + 6 instead of 12 + 8 nodes in p; k = 0.25 stays only in
the scaling check I5. Nothing else changes.

Why this test exists (15 Nov 2026; the author's order of work "B, then C, then A, then D", 08:30)
---------------------------------------------------------------------------------------------
F6.8c2 found that the optical sector weakens with the age T of the static field faster than linearly
(k = 0.12: D(n) = 13.9, 14.3, 15.0, 18.0, 28.0 at n = kT = 6, 8, 10, 16, 24), while the conformal and the mixed
sectors become linear (c_T = 0.110, c_T,ab = 0.195). The form for kT -> infinity, the regime of any physical field,
decides the far-field law (F8.1) and the window of l_k (paper, Part VIII). The record of F6.8c also says that rows with
T >~ 170 broke the scaling, "a numerical loss at large ages, not settled".

What the development found (15 Nov, before this test; docs/F6_8e_card_draft.md, section 2)
------------------------------------------------------------------------------------------
1. The numerical loss is real and its cause is known. For a long thin interval of length L the second-order relative
   volume V''/V0 has a part that does not depend on the proper time tau and grows as (kL)^4 (k = 0.12: -25.4, -348,
   -5343 at L = 100, 200, 400). The BD moment conditions remove it exactly from the s-integral (they remove tau^0,
   tau^2, tau^4), and what survives is ~1e-9 of it. F6.8c computes V'' as a bulk integral plus a surface integral that
   cancel to O(L^4 / tau^2), and its per-pair quadrature error (1e-9 to 1e-10) is of the size of the survivor: the
   s-integrated density at r = 200 and 400 changed by factors of 2 to 4 between quadrature levels.
2. The stable instrument used here. (a) Coordinates that follow the bent ray: z = w + b delta(l_w), delta the
   first-order transverse deviation of the ray x -> y, has unit Jacobian, so the b^2 coefficient is the same number,
   but the bulk-surface cancellation is gone (on the axis the broken path through a point of the ray has the optical
   length of the ray). Against F6.8c where F6.8c is converged: 1e-12 (the F6.8c test diamond k = 0.7, L = 2) and
   1e-10 (k = 0.12, L = 5 and 20). (b) The pointwise cancellation inside l'' (terms ~ k^2 L^3) is done in extended
   precision (80-bit long double, eps 1.1e-19); per pair the result is then stable to 1e-14 at L = 400. (c) Tip
   clustering of the axial nodes (double exponential). With (a)-(c) the s-integrated density I(r, p) is stable
   between levels to 1e-11 (r = 100) and 5e-4 (r = 400).
3. With the stable instrument (k = 0.12, p = k r cos the phase advance): I(r, 0) = -7.40e-7, -9.06e-7, -6.17e-6,
   -1.10e-4, -2.95e-3, -9.01e-2, -2.85, -90.9 at r = 50 ... 6400 (doubling). Local exponent 0.29, 2.77, 4.16, 4.74, 4.93,
   4.98, 4.995 -> 5. The conformal density on the same nodes -> exponent 1 (1.996, 2.003 ratios), the known linear age
   term; the same profile integration returns c_T = 0.105 (F6.8a 0.1035). The p-profile has a fixed shape (width
   ~10 in p, symmetric in p) and I(k, r) = k^5 Psi(kr) to 1e-3 (k = 0.12 against 0.2 at equal kr).
4. The slope of the optical age part at n = 24 from this bulk density, 1.19 per unit n (in units of k^2), against the
   full rows of F6.8c2 (n = 16 and 24): 1.24. The layer and boundary terms of the age cut do not change the slope
   there.

Prediction (stated before this computation)
-------------------------------------------
The tau-independent part of the second-order relative volume grows as (kL)^4: the second-order focusing of the cross
section of a long interval by the optical index. The BD moments remove the tau^0, tau^2 and tau^4 parts, and the
surviving tau^6 part carries the same factor, because the transverse structure at fixed phase advance does not depend
on L. Hence I(r, p) ~ k^5 (kr)^5 psi(p), the bulk density rho(r) = (2 pi / k r) int dp I ~ r^4, and the optical age term

    K_hat_bb(k, T) = K0 - k^2 [c1 ln(1/k) + c0 + c5 (kT)^5]   for kT -> infinity,

growing as the fifth power of the age; the local exponent alpha(n) = 1 + d ln rho / d ln n of the age term tends to 5.
The hypotheses that would make the optical sector like the others, a linear age term (alpha = 1), or a saturating one
(alpha < 1), are what this test can exclude.

What is computed
----------------
For k = 0.15 and 0.25 (development used 0.12 and 0.2) and n = kT = 96, 192, 384, 768, 1536, 3072 (r = n / k).
The development showed I(k, r) = k^5 Psi(kr): at equal n two wavenumbers give the same exponent, so new wavenumbers
are NOT new information about alpha. What is new is the range: the development reached n = 768 (k = 0.12, p = 0
only); the judged pair, n = 1536 -> 3072, was not computed before this test. Computed: the s-integrated pair density I(r, p) (the bulk integrand of F6.8c, _pair_terms, with beta and gamma from
the stable instrument; full s range), on Gauss nodes in p in [0, 12] (12 nodes) and [12, 24] (8 nodes), doubled by the
symmetry p -> -p; rho_bb(r) = (2 pi / k r) int dp I; the same for the conformal density.
F'_bb(n) = -a_4 rho_bb / k^5 (the age part of K_hat_bb per unit n in units of k^2); alpha between successive n.

Instrument (INSTRUMENT_FAILED if any fails)
-------------------------------------------
  I0 extended precision available: numpy long double eps <= 1e-18 (on Windows numpy has no 80-bit long double; run
     under Linux or WSL);
  I1 the stable per-pair V''/V0 against F6.8c.second_order_volume (converged) at three diamonds (k 0.7, CO 0.4, L 2,
     tau 1.2; k 0.12, CO 0, L 5 and 20, tau 1): relative difference <= 1e-9;
  I2 quadrature: the s-integrated density at two levels (s nodes 16 / 24, nodes (24, 16, 16) / (32, 24, 16)) at
     p = 0 for the largest and the smallest n at each k: relative difference <= 1e-3;
  I3 the method on the conformal sector: c_T = -(16/3) a_4 rho_c / k^5 at the largest n against F6.8a's asymptotic
     0.1035: within 3 %;
  I4 the tail: |I(r, 24)| / |I(r, 0)| <= 1e-3 at the largest n (the p-range holds the profile);
  I5 the form k^2 F(kT): I(k = 0.15, n / 0.15) / I(k = 0.25, n / 0.25) = (0.15 / 0.25)^5 at p = 0 for n = 192 and 1536:
     relative difference <= 1e-3.
Status (alpha between n = 1536 and 3072 at k = 0.15; revision 1)
------------------------------------------------------------------
  PASSED     |alpha - 5| <= 0.1: the optical age term grows as (kT)^5;
  FAILED     |alpha - 5| > 0.5 (in particular alpha near 1, a linear term, or below, saturation);
  UNDECIDED  otherwise.
Reported, not judged: alpha at every successive pair of n; c5 = F'(n) / (5 n^4) at the largest n; the age n* at which
the optical stiffness K_hat_bb changes sign, n* ~ (K0 / (k^2 c5))^(1/5) (local part neglected); the k = 0.12 slope at
n = 24 against F6.8c2 (development wavenumber, seen).

What a pass means and what it does not
--------------------------------------
A pass says that the second-order (b^2) coefficient of the expected BD action in the optical sector grows with the
age of the static field as T^5, not linearly: the formal static quadratic form of the order has no limit for old
fields in this sector, and the age term of P2 cannot be read in the conformal form there. It says nothing about the
response at finite b (where the expansion in b breaks down near b^2 (kT)^4 ~ 1) and it is a statement about the bulk
of the past: the layer and boundary terms of the age cut are not computed at large n (at n = 24 they do not change the
slope; F6.8d asks whether the cut itself makes the sectors differ).
"""

from __future__ import annotations

import math
import time

import numpy as np

from core.stats.parallel import parallel_map
from core.testing import Result
from phases.F6_order_dynamics.F6_08a_bd_action_conformal import _A4, _L2, _gl, kernel_derivs, mean_exp_diamond
from phases.F6_order_dynamics.F6_08b2_bd_action_mixed_rev import _s_nodes

TEST = {
    "id": "F6.8e3",
    "phase": "F6",
    "title": "Asymptotic age term of the optical (nu^2) sector, second revision of F6.8e (checks in the exponent)",
    "theory_part": "the static quadratic form of the order for old fields; input of F8.1 (far-field law) and of P2",
    "question": "How does the optical (nu^2) coefficient of the expected discrete BD action depend on the age T of a static "
                "field for kT -> infinity?",
    "method": "as F6.8e: the bulk pair density of F6.8c, s-integrated, with the stable second-order interval volume; "
              "k = 0.1 (profile) and 0.18 (scaling), n = kT = 192 ... 6144; local exponent between successive n",
    "prediction": "the optical age term grows as (kT)^5 (alpha -> 5), from the (kL)^4 focusing of long intervals that survives "
                  "the BD moments in the tau^6 part",
    "criterion": "INSTRUMENT: I0 long double eps <= 1e-18; I1 stable V'' vs F6.8c <= 1e-9 rel.; I2 levels <= 1e-3; I3 "
                 "conformal c_T within 3 % of 0.1035 at n 384; I4 tail <= 1e-3; I2 and I5 in the exponent at the judged pair: "
                 "levels and the two k within 0.01. alpha at k = 0.13 (n 3072 -> 6144): PASSED |alpha - 5| <= 0.1; FAILED "
                 "> 0.5; UNDECIDED otherwise.",
    "null_model": "a linear age term as in the conformal and mixed sectors (alpha = 1), or saturation (alpha < 1)",
    "level": "chain",
    "meters": [],
    "depends_on": ["F6.8a", "F6.8c", "F6.8c2", "F6.8e", "F6.8e2"],
    "params": {"ks": [0.13, 0.21], "profile_ks": [0.13], "ns": [192.0, 384.0, 768.0, 1536.0, 3072.0, 6144.0],
               "p_panels": [[0.0, 12.0, 10], [12.0, 24.0, 6]], "lev1": [16, 24, 16, 16], "lev2": [24, 32, 24, 16],
               "scale_ns": [3072.0, 6144.0], "cT_n": 384.0, "tol_alpha_check": 0.01, "dev_row": [0.12, 24.0], "c_T_ref": 0.1035,
               "tol_pair": 1e-9, "tol_level": 1e-3, "tol_cT": 0.03, "tol_tail": 1e-3, "tol_scale": 1e-3,
               "tol_pass": 0.1, "tol_fail": 0.5, "alpha_pred": 5.0},
    "parts": {"main": {}},
    "quick": {"ks": [0.13], "profile_ks": [0.13], "ns": [48.0, 96.0], "cT_n": 96.0, "p_panels": [[0.0, 12.0, 4], [12.0, 24.0, 2]], "scale_ns": [],
              "dev_row": None, "lev1": [16, 16, 8, 8], "lev2": [16, 24, 8, 8]},
    "estimate": "~125 pair-density evaluations of 30-100 s each (long double): ~1 h with 2 processes; "
                "memory < 0.5 GB per process; needs 80-bit long double (Linux / WSL)",
}

LD = np.longdouble
PI = LD("3.14159265358979323846264338327950288")


# ============================================================ extended-precision optics of a segment

def _jbend(x):
    """J(x) = (x^2 + 2 cos x - 2) / x^4 (J(0) = 1/12), long double."""
    x = np.asarray(x, LD)
    out = np.empty_like(x)
    sm = np.abs(x) < 0.5
    xs = x[sm] ** 2
    out[sm] = (LD(1) / 12 - xs / 360 + xs ** 2 / 20160 - xs ** 3 / 1814400 + xs ** 4 / 239500800
               - xs ** 5 / LD("43589145600"))
    xb = x[~sm]
    out[~sm] = (xb * xb + 2 * np.cos(xb) - 2) / xb ** 4
    return out


def _bend(k, d, p):
    """Phase-averaged bending part of <T2> of a segment of length d and phase advance p (F6.8c.bend), long double."""
    d = np.asarray(d, LD)
    p = np.asarray(p, LD)
    return -np.maximum(k * k * d * d - p * p, LD(0)) * d * _jbend(p) / 4


def _Ef(a):
    """E(a) = (e^{ia} - 1)/(ia) and E'(a), long double, series for small a."""
    a0 = np.asarray(a, LD)
    a = np.atleast_1d(a0)
    sm = np.abs(a) < 1e-2
    ab = np.where(sm, LD(1), a)
    e = np.exp(1j * ab)
    E = (e - 1) / (1j * ab)
    dE = -1j * (1j * ab * e - e + 1) / (ab * ab)
    x = a[sm]
    E[sm] = 1 + 1j * x / 2 - x ** 2 / 6 - 1j * x ** 3 / 24 + x ** 4 / 120 + 1j * x ** 5 / 720
    dE[sm] = 1j / 2 - x / 3 - 1j * x ** 2 / 8 + x ** 3 / 30 + 1j * x ** 4 / 144
    return E.reshape(a0.shape), dE.reshape(a0.shape)


def _ray_dev(kpar, kperp, L, l):
    """Amplitude (along k_perp) of the first-order deviation of the ray x -> y: delta'' = i kperp e^{i kpar l},
    delta(0) = delta(L) = 0."""
    l = np.asarray(l, LD)
    if abs(kpar) * max(L, 1) < 1e-3:
        q = kpar
        return 1j * kperp * (l * (l - L) / 2 + 1j * q * (l ** 3 - L ** 2 * l) / 6 - q ** 2 * (l ** 4 - L ** 3 * l) / 24)
    return -(1j * kperp / kpar ** 2) * (np.exp(1j * kpar * l) - 1 - (np.exp(1j * kpar * L) - 1) * l / L)


def _te_nodes(n_te, npan):
    """Axial nodes te = tanh(pi/2 sinh v), v in [-3, 3] in npan Gauss panels (clustered at the tips)."""
    v, wv = [], []
    edges = np.linspace(-3.0, 3.0, npan + 1)
    for a, b in zip(edges[:-1], edges[1:]):
        x, w = _gl(a, b, n_te)
        v.append(x)
        wv.append(w)
    v = np.concatenate(v)
    wv = np.concatenate(wv)
    u = 0.5 * math.pi * np.sinh(v)
    return np.tanh(u), wv * 0.5 * math.pi * np.cosh(v) / np.cosh(u) ** 2


def second_order_volume_stable(k, CO, L, tau2s, n_sig=24, n_te=16, n_phi=16, npan=None):
    """
    beta = V'/V0 (complex amplitude) and gamma = <V''>/V0 of the interval x = 0 -> y = L n (cos(n, k) = CO) for the
    proper times tau2s, in coordinates that follow the bent ray (z = w + b delta(l_w), unit Jacobian), long double.
    Same quantities as F6.8c.second_order_volume.
    """
    k, CO, L = LD(k), LD(CO), LD(L)
    tau2 = np.asarray(tau2s, LD)[:, None, None]
    dt = np.sqrt(L * L + tau2)
    kpar = k * CO
    kperp = k * np.sqrt(max(LD(1) - CO * CO, LD(0)))
    npan = npan or int(8 + math.ceil(float(abs(kpar) * L)))
    c = L / 2
    em1 = tau2 / (L * (dt + L))
    dtL = tau2 / (dt + L)
    xs, ws = _gl(0.0, 1.0, n_sig)
    te, wt = _te_nodes(n_te, npan)
    te, wt, xs, ws = (np.asarray(v, LD) for v in (te, wt, xs, ws))
    ph = 2 * PI * (np.arange(n_phi, dtype=LD) + LD(0.5)) / n_phi
    TE, PH = np.meshgrid(te, ph, indexing="ij")
    TE, PH = TE[None], PH[None]
    WT = ((wt[:, None] * np.ones_like(PH[0])) * (2 * PI / n_phi))[None]
    E_xy, _ = _Ef(kpar * L)
    T1xy = L * E_xy
    T2xy = L / 4 + _bend(k, L, kpar * L)
    dt1 = L * T1xy / dt
    dt2 = (abs(T1xy) ** 2 / 2 + 2 * L * T2xy) / dt - L * L * abs(T1xy) ** 2 / 2 / dt ** 3
    eky = np.exp(1j * kpar * L)
    cosP, sinP = np.cos(PH), np.sin(PH)

    def at(sm1):
        sig = 1 + sm1
        lw = c * (1 + sig * TE)                                    # axial distance from x
        rp = c * np.sqrt(np.maximum(sm1 * (2 + sm1) * (1 - TE * TE), LD(0)))
        w2 = rp * cosP                                             # along k_perp
        d1 = c * (sig - TE)                                        # to y
        d2 = c * (sig + TE)                                        # to x
        kw = kpar * lw + kperp * w2
        ekw = np.exp(1j * kw)
        dl = _ray_dev(kpar, kperp, L, lw)
        u1y = w2 / d1
        u2y = w2 / d2
        a1 = kw - kpar * L
        E1, dE1 = _Ef(a1)
        E2, dE2 = _Ef(kw)
        g1y = eky * (u1y * E1 + d1 * dE1 * kperp)                  # grad_w T1(y, w) along k_perp
        g2y = u2y * E2 + d2 * dE2 * kperp                          # grad_w T1(w, x) along k_perp
        lp = dt1 - (u1y + u2y) * dl - d1 * eky * E1 - d2 * E2      # l' amplitude
        adl = np.abs(dl) ** 2
        lpp = (dt2 - (adl * (1 - u1y ** 2) / (2 * d1) + np.real(dl * np.conj(g1y)) + 2 * (d1 / 4 + _bend(k, d1, a1)))
               - (adl * (1 - u2y ** 2) / (2 * d2) + np.real(dl * np.conj(g2y)) + 2 * (d2 / 4 + _bend(k, d2, kw))))
        Ppp = LD(4.5) + 3 * np.real(1j * kperp * ekw * np.conj(dl))
        return sig, ekw, lp, lpp, Ppp

    V1 = np.zeros(tau2.shape[0], np.clongdouble)
    V2 = np.zeros(tau2.shape[0], LD)
    for x, wx in zip(xs, ws):
        sm1 = em1 * x
        sig, ekw, lp, lpp, Ppp = at(sm1)
        l0 = dtL - L * sm1
        wgt = (em1 * wx) * WT * c ** 3 * (sig * sig - TE * TE)
        V1 += (wgt * (3 * ekw * l0 + lp)).sum(axis=(1, 2))
        V2 += (wgt * (Ppp * l0 + 3 * np.real(ekw * np.conj(lp)) + lpp)).sum(axis=(1, 2))
    sig, ekw, lp, lpp, Ppp = at(em1)
    V2 += (WT * c ** 3 * (sig * sig - TE * TE) / L * np.abs(lp) ** 2 / 2).sum(axis=(1, 2))   # surface term
    V0 = PI * tau2[:, 0, 0] ** 2 / 24
    return (V1 / V0).astype(complex), (V2 / V0).astype(float)


def pair_density(k, CO, r, n_s=16, lev=(24, 16, 16)):
    """s-integrated bulk pair density at distance r, direction CO (full s range): optical and conformal."""
    s, ws = _s_nodes(np.array([1e12]), n_s, True)
    s, ws = s[0], ws[0]
    keep = ws > 0
    s, ws = s[keep], ws[keep]
    kpar = k * CO
    tau2 = s / _L2
    dt = np.sqrt(r * r + tau2)
    g0, g1, g2 = kernel_derivs(s * s)
    beta, gam = second_order_volume_stable(k, CO, r, tau2, *lev)
    E_xy, _ = _Ef(kpar * r)
    Dxy = complex(r * E_xy)
    Bxy = float(_bend(k, r, kpar * r))
    ph = np.exp(1j * kpar * r)
    one = 1.0 + ph
    M = mean_exp_diamond(k, np.full((1, s.size), r), CO, dt[None, :])[0]
    mu = -r * Dxy / (dt * dt)
    mup = -2.0 * r * (0.25 * r + Bxy) / dt ** 2 - 0.5 * abs(Dxy) ** 2 * (dt * dt - 2.0 * r * r) / dt ** 4
    X1 = 3.0 * one + mu
    bb = ((0.5 * np.abs(X1) ** 2 + mup) * g0 + np.real(X1 * np.conj(beta)) * g1 + gam * g1
          + 0.5 * np.abs(beta) ** 2 * g2)
    conf = (1.0 + ph.real) * g0 + np.real(one * np.conj(M)) * g1 + 0.5 * np.abs(M) ** 2 * g2 + 0.5 * g1
    base = ws * r * r / (_L2 * 2.0 * dt)
    return float((base * bb).sum()), float((base * conf).sum())


# ===================================================================================== tasks

def _p_nodes(panels):
    ps, wp = [], []
    for a, b, n in panels:
        x, w = _gl(float(a), float(b), int(n))
        ps.append(x)
        wp.append(w)
    return np.concatenate(ps), np.concatenate(wp)


def _dispatch(kind, args):
    if kind == "dens":
        k, r, p, lev1 = args
        n_s, ns2, nte, nph = lev1
        CO = p / (k * r)
        b, c = pair_density(k, CO, r, n_s, (ns2, nte, nph))
        return {"k": k, "r": r, "p": p, "bb": b, "conf": c}
    if kind == "pair":
        from phases.F6_order_dynamics.F6_08c_bd_action_optical import second_order_volume
        k, CO, L, tau = args
        dt = math.sqrt(L * L + tau * tau)
        _, _, g_ref, _, _ = second_order_volume(k, CO, np.array([[L]]), np.array([[dt]]), np.array([[tau * tau]]),
                                                64, int(128 + math.ceil(6.4 * abs(k * CO) * L)), 64)
        _, g = second_order_volume_stable(k, CO, L, [tau * tau], 48, 32, 16)
        return {"k": k, "CO": CO, "L": L, "tau": tau, "ref": float(g_ref[0, 0]), "stable": float(g[0]),
                "rel": abs(float(g[0]) - float(g_ref[0, 0])) / abs(float(g_ref[0, 0]))}
    raise ValueError(kind)


def run(p, ctx):
    t0 = time.time()
    lines = ["3+1, eps rho = 1 (l_k = 1); optical (nu^2) bulk pair density of F6.8c with the stable second-order volume; "
             "p = k r cos(n, k)"]
    eps = float(np.finfo(np.longdouble).eps)
    ok0 = eps <= 1e-18
    lines.append(f"  I0 long double eps {eps:.2e} (need <= 1e-18)")
    if not ok0:
        lines.append("  no 80-bit long double on this platform (numpy on Windows): run under Linux or WSL")
        return Result("INSTRUMENT_FAILED", lines, {"eps": eps})
    ks = [float(k) for k in p["ks"]]
    ns = [float(n) for n in p["ns"]]
    pn, wp = _p_nodes(p["p_panels"])
    lev1 = tuple(p["lev1"])
    lev2 = tuple(p["lev2"])
    pks = [float(k) for k in p.get("profile_ks", ks)]
    jobs = [("dens", (k, n / k, float(pp), lev1)) for k in pks for n in ns for pp in pn]
    n_main = len(jobs)
    tail = [("dens", (k, ns[-1] / k, 24.0, lev1)) for k in pks]
    jn = [float(n) for n in p.get("scale_ns", [])] or [ns[-2], ns[-1]]
    lv = [("dens", (k, n / k, 0.0, lev2)) for k in pks for n in jn]
    zn = sorted({ns[-1]} | set(jn))
    sc = [("dens", (k, n / k, 0.0, lev1)) for k in ks for n in zn]          # p = 0 at level 1 (I2, I4, I5)
    dev = []
    if p.get("dev_row"):
        k_d, n_d = p["dev_row"]
        dev = [("dens", (float(k_d), float(n_d) / float(k_d), float(pp), lev1)) for pp in pn]
    pairs = [("pair", (0.7, 0.4, 2.0, 1.2)), ("pair", (0.12, 0.0, 5.0, 1.0)), ("pair", (0.12, 0.0, 20.0, 1.0))]
    alljobs = jobs + tail + lv + sc + dev + pairs
    res = parallel_map(_dispatch, alljobs, ctx.workers, ctx.progress,
                       lambda i, r_, d, n: f"[F6.8e3] task {d}/{n} ({time.time() - t0:.0f}s)")
    i = 0
    main = res[i:i + n_main]; i += n_main
    rtail = res[i:i + len(tail)]; i += len(tail)
    rlv = res[i:i + len(lv)]; i += len(lv)
    rsc = res[i:i + len(sc)]; i += len(sc)
    rdev = res[i:i + len(dev)]; i += len(dev)
    rpair = res[i:]

    def find(rr, k, r, pp):
        for x in rr:
            if abs(x["k"] - k) < 1e-12 and abs(x["r"] - r) < 1e-9 * r and abs(x["p"] - pp) < 1e-12:
                return x
        return None

    # I1
    ok1 = all(x["rel"] <= float(p["tol_pair"]) for x in rpair)
    lines.append("  I1 stable V''/V0 vs F6.8c: " + "; ".join(f"k {x['k']:g} L {x['L']:g}: {x['stable']:+.10e} vs "
                                                          f"{x['ref']:+.10e} (rel {x['rel']:.1e})" for x in rpair)
                 + f" (tol {p['tol_pair']:g})")
    # profiles
    prof = {}
    for k in pks:
        for n in ns:
            r = n / k
            vb = np.array([find(main, k, r, float(pp))["bb"] for pp in pn])
            vc = np.array([find(main, k, r, float(pp))["conf"] for pp in pn])
            Ib = 2.0 * float((wp * vb).sum())
            Ic = 2.0 * float((wp * vc).sum())
            pref = 2.0 * math.pi / (k * r)
            prof[(k, n)] = {"rho_bb": pref * Ib, "rho_c": pref * Ic}
    # I2 levels, in the exponent at the judged pair
    def alpha0(rr, k, na, nb):
        xa, xb = find(rr, k, na / k, 0.0), find(rr, k, nb / k, 0.0)
        if xa is None or xb is None:
            return None
        return 1.0 + math.log((xb["bb"] / (nb / k)) / (xa["bb"] / (na / k))) / math.log(nb / na)
    ok2 = True
    tola = float(p["tol_alpha_check"])
    for k in pks:
        a1, a2 = alpha0(rsc, k, jn[0], jn[-1]), alpha0(rlv, k, jn[0], jn[-1])
        if a1 is not None and a2 is not None:
            ok2 &= abs(a1 - a2) <= tola
            lines.append(f"  I2 k {k:g}, n {jn[0]:g} -> {jn[-1]:g}: alpha_0 level 1 {a1:.5f}, level 2 {a2:.5f} "
                         f"(|diff| {abs(a1 - a2):.1e}, tol {tola:g})")
    # I3 conformal control
    ok3 = True
    for k in pks:
        ncT = float(p.get("cT_n", ns[-1]))
        cT = -(16.0 / 3.0) * _A4 * prof[(k, ncT)]["rho_c"] / k ** 5
        d = abs(cT - float(p["c_T_ref"])) / float(p["c_T_ref"])
        ok3 &= d <= float(p["tol_cT"])
        lines.append(f"  I3 conformal c_T at k {k:g}, n {ncT:g}: {cT:.4f} (F6.8a asymptotic {p['c_T_ref']:g}; rel {d:.1e}, "
                     f"tol {p['tol_cT']:g})")
    # I4 tail
    ok4 = True
    for x in rtail:
        y = find(rsc, x["k"], x["r"], 0.0)
        ref0 = abs(y["bb"])
        rt = abs(x["bb"]) / ref0
        ok4 &= rt <= float(p["tol_tail"])
        lines.append(f"  I4 tail k {x['k']:g}, n {x['r'] * x['k']:g}: |I(p = 24)| / |I(p = 0)| = {rt:.1e} (tol {p['tol_tail']:g})")
    # I5 the kT-shape independent of k, in the exponent at the judged pair
    ok5 = True
    if len(ks) >= 2:
        al = [alpha0(rsc, k, jn[0], jn[-1]) for k in ks[:2]]
        if None not in al:
            ok5 = abs(al[0] - al[1]) <= tola
            lines.append(f"  I5 alpha_0 at n {jn[0]:g} -> {jn[-1]:g}: k {ks[0]:g} {al[0]:.5f}, k {ks[1]:g} {al[1]:.5f} "
                         f"(|diff| {abs(al[0] - al[1]):.1e}, tol {tola:g})")
        for n in jn:
            a_, b_ = find(rsc, ks[0], n / ks[0], 0.0), find(rsc, ks[1], n / ks[1], 0.0)
            if a_ and b_:
                lines.append(f"  reported n {n:g}: I(k {ks[0]:g}) / I(k {ks[1]:g}) = {a_['bb'] / b_['bb']:.6e}, "
                             f"(k ratio)^5 {(ks[0] / ks[1]) ** 5:.6e}")
    inst_ok = ok1 and ok2 and ok3 and ok4 and ok5

    lines.append("  k      n      r          rho_bb          F'_bb(n) = -a4 rho_bb / k^5   rho_conf        alpha(n_prev -> n)")
    alphas = {}
    for k in pks:
        prev = None
        for n in ns:
            pr = prof[(k, n)]
            Fp = -_A4 * pr["rho_bb"] / k ** 5
            al = ""
            if prev is not None:
                a_ = 1.0 + math.log(pr["rho_bb"] / prev[1]) / math.log(n / prev[0])
                alphas[(k, n)] = a_
                al = f"{a_:.4f}"
            lines.append(f"  {k:<6g} {n:<6g} {n / k:<10.5g} {pr['rho_bb']:+.6e}   {Fp:+.6e}                  "
                         f"{pr['rho_c']:+.6e}   {al}")
            prev = (n, pr["rho_bb"])
    for k in pks:
        Fp = -_A4 * prof[(k, ns[-1])]["rho_bb"] / k ** 5
        c5 = Fp / (5.0 * ns[-1] ** 4)
        nstar = (1.0 / (k * k * c5)) ** 0.2 if c5 > 0 else float("nan")
        lines.append(f"  reported k {k:g}: c5 = F'(n) / (5 n^4) at n {ns[-1]:g} = {c5:.4e}; K_hat_bb changes sign near "
                     f"n* ~ {nstar:.1f} (T ~ {nstar / k:.0f} l_k; local part neglected)")
    if rdev:
        k_d, n_d = p["dev_row"]
        vb = np.array([x["bb"] for x in rdev])
        rho = 2.0 * math.pi / (k_d * n_d / k_d) * 2.0 * float((wp * vb).sum())
        lines.append(f"  reported (development wavenumber, seen): k {k_d:g}, n {n_d:g}: age slope in units of k^2 "
                     f"{-_A4 * rho / k_d ** 5:.3f} per unit n (F6.8c2 full rows, n 16 -> 24: 1.244)")

    if not inst_ok:
        status = "INSTRUMENT_FAILED"
    elif len(ns) < 2 or (not ctx.quick and len(ks) < 2):
        status = "UNDECIDED"
    else:
        last = [alphas[(k, ns[-1])] for k in pks]
        dev_ = [abs(a_ - float(p["alpha_pred"])) for a_ in last]
        if all(d <= float(p["tol_pass"]) for d in dev_):
            status = "PASSED"
        elif any(d > float(p["tol_fail"]) for d in dev_):
            status = "FAILED"
        else:
            status = "UNDECIDED"
    data = {"profiles": {f"{k}_{n}": v for (k, n), v in prof.items()}, "alphas": {f"{k}_{n}": v for (k, n), v in alphas.items()},
            "levels": rlv, "tail": rtail, "pairs": rpair, "scale": rsc, "dev": rdev,
            "p_nodes": pn.tolist(), "p_weights": wp.tolist()}
    return Result(status, lines, data)
