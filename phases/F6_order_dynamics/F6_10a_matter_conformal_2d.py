"""
F6.10a -- The response of the pairwise matter action to a conformal change of a 2D order (OV2 T1)
================================================================================================

FROZEN by the author, 21 Nov 2026 09:17 (decision 40; card docs/OV2_coupling_cards_draft.md rev. 3, T1). Prepared by
Claude on 21 Nov at the author's request.

Why this test exists
--------------------
OV2 (docs/derivations/order_variation_OV2.md, decision 39) gives matter an action that depends on the order only:
S_m = S_kin + S_pot, S_kin = -eps^2 sum_pairs f(n, eps) (phi_x - phi_y)^2 (the pairwise, Dirichlet form: exactly
invariant under phi -> phi + 2 pi k), S_pot = -(1/rho) sum V(phi_v). Matter acts on the order only if its action
responds to the metric as the continuum does: delta S_m = (1/2) int T^{ab} delta g_ab sqrt|g| at fixed matter (T_ab
off shell, decision 39). In 2D every metric is locally conformal to flat; for g -> e^{2 sigma} g the response is
int sigma T^a_a sqrt|g| with T^a_a = -2V for sine-Gordon: the kinetic part is conformally invariant, the potential
part is the measure.

What is computed (ensemble level, deterministic; no sprinkling)
---------------------------------------------------------------
For a Poisson sprinkling, E[f(n, eps)] with n ~ Poisson(Lam) is exactly g(eps Lam) = e^{-z}(1 - 2z + z^2/2), so by
the Mecke formula
    <S_kin>[sigma] = -eps^2 int int_{x<y} lam_x lam_y g(eps Lam_xy) (phi_x - phi_y)^2,   lam = N e^{2 sigma},
on the light-cone square [0, 1]^2 (ds^2 = -du dv, proper density 2N). sigma = a sh, sh = exp(-((u-u0)^2 + (v-v0)^2)
/ (2 s^2)); the first order in a is the analytic derivative under the integral (exact rectangle integrals of the
Gaussian). The kink phi = 4 arctan(e^{m x}), x = (v - u)/2. The target int sigma T^a_a sqrt|g| = -a int sh V du dv
equals the first order of <S_pot> exactly, so the judged quantity is the kinetic share
    R = d<S_kin>/da / (-int sh V du dv).
Variables per x: t = ln du (lower cut 1e-10 / (eps N): the near-null pairs), w = eps N du dv; outer nodes clustered at
the future edges (development: OV2 sec. 10, where a fixed cut 1e-7 was found to be an error growing as N^2).

Parameters (fresh; development used m = 12, s = 0.12, centre (0.5, 0.5)): eps = 0.12, N = 66 560 and 266 240 (T/w = 32
and 64 of OV1), m = 9, s = 0.09, centre (0.45, 0.55).

Criterion
---------
Per density, with Delta the change of R between the two quadrature levels:
  PASSED if |R| + Delta <= 0.03 at both densities; FAILED if |R| - Delta > 0.06 at either; UNDECIDED otherwise.
Instrument (INSTRUMENT_FAILED if any fails):
  I1 the flat <S_kin> at N = 266 240 within 1 % of the continuum int phi_u phi_v du dv;
  I2 quadrature: Delta <= 0.01 at both densities (absolute, in units of the target);
  I3 the derivative of <S_pot> computed on the quadrature grid against the closed form -int sh V within 1e-6.
Reported: R per density and level, the flat <S_kin> against the continuum, the development rows (OV2 sec. 10).
Not part of this test: the sprinkled cross-check (card I4); decision 40: not required.

Grid choice (development, 21 Nov, Claude): R at levels (64,48,48)/(96,64,64)/(128,96,64)/(192,128,96)/
(256,160,128) = +0.551/-0.006/-0.0139/-0.0139/-0.0139 at N 66 560 and -8.49/-0.341/-0.0048/-0.0047/-0.0047 at
N 266 240; the judged levels are the last two. Since the potential part is exact, |delta_1 <S_m>/target - 1| = |R|.
"""

from __future__ import annotations

import math
import time

import numba as nb
import numpy as np

from core.testing import Result

TEST = {
    "id": "F6.10a",
    "phase": "F6",
    "title": "Response of the pairwise matter action to a conformal change of a 2D order (OV2 T1)",
    "theory_part": "OV2 (docs/derivations/order_variation_OV2.md): the coupling of matter to the order; T_ab from "
                   "the first variation of the matter action at fixed matter",
    "question": "Is the first-order response of the expected pairwise matter action of a kink to a conformal bump "
                "int sigma T^a_a sqrt|g| with T^a_a = -2V, the kinetic part within 3 % of the target?",
    "method": "Mecke formula with the exact kernel average g(eps Lam); analytic derivative under the integral; "
              "light-cone square, eps 0.12, N 66560 and 266240; two quadrature levels",
    "prediction": "kinetic share R -> 0 (conformal invariance of the kinetic term in 2D); development: R = 0.015 at "
                  "T/w 32, 0.005 at 64 (m = 12)",
    "criterion": "PASSED if |R| + Delta <= 0.03 at both densities, FAILED if |R| - Delta > 0.06 at either, "
                 "else UNDECIDED; INSTRUMENT_FAILED if I1 (flat S_kin within 1 %), I2 (Delta <= 0.01) or I3 "
                 "(potential derivative to 1e-6) fails",
    "null_model": "the potential part (the measure), whose first order equals the target exactly",
    "level": "mechanism",
    "meters": [],
    "depends_on": ["F6.9a"],
    "params": {"eps": 0.12, "n_exp": [66560, 266240], "m": 9.0, "s": 0.09, "u0": 0.45, "v0": 0.55,
               "levels": [[192, 128, 96], [256, 160, 128]], "tol": 0.03, "fail": 0.06, "tol_flat": 0.01,
               "tol_delta": 0.01},
    "parts": {"main": {}},
    "quick": {"n_exp": [266240], "levels": [[128, 96, 64], [192, 128, 96]]},
    "estimate": "~4 min on one core (deterministic quadrature)",
}


@nb.njit(cache=True, inline="always")
def _phi(u, v, m):
    return 4.0 * np.arctan(np.exp(m * 0.5 * (v - u)))


@nb.njit(cache=True, inline="always")
def _gint(a1, a2, c, s):
    """int_{a1}^{a2} exp(-(x - c)^2 / (2 s^2)) dx."""
    d = a2 - a1
    if d < 1e-7 * s:
        xm = 0.5 * (a1 + a2) - c
        return d * np.exp(-xm * xm / (2.0 * s * s))
    r = 1.0 / (s * np.sqrt(2.0))
    return s * np.sqrt(np.pi / 2.0) * (math.erf((a2 - c) * r) - math.erf((a1 - c) * r))


@nb.njit(cache=True)
def _quad(n_exp, eps, m, s, u0, v0, xo, wo, xt, wt, xw, ww, pmin, wcap):
    """(flat sum, d/da at a = 0) of <sum f(n) (phi_x - phi_y)^2> (without -eps^2)."""
    flat = 0.0
    d1 = 0.0
    for i in range(xo.size):
        u1 = xo[i]
        for j in range(xo.size):
            v1 = xo[j]
            shx = np.exp(-((u1 - u0) ** 2 + (v1 - v0) ** 2) / (2.0 * s * s))
            px = _phi(u1, v1, m)
            tlo = np.log(pmin)
            thi = np.log(1.0 - u1)
            if thi <= tlo:
                continue
            af = 0.0
            ad = 0.0
            for k in range(xt.size):
                t = tlo + (thi - tlo) * xt[k]
                p = np.exp(t)
                dp = p * (thi - tlo) * wt[k]
                wmax = min(wcap, eps * n_exp * p * (1.0 - v1))
                if wmax <= 0.0:
                    continue
                Iu = _gint(u1, u1 + p, u0, s)
                for l in range(xw.size):
                    w = wmax * xw[l]
                    q = w / (eps * n_exp * p)
                    dq = wmax * ww[l] / (eps * n_exp * p)
                    u2 = u1 + p
                    v2 = v1 + q
                    ez = np.exp(-w)
                    g = ez * (1.0 - 2.0 * w + 0.5 * w * w)
                    gp = ez * (-3.0 + 3.0 * w - 0.5 * w * w)
                    shy = np.exp(-((u2 - u0) ** 2 + (v2 - v0) ** 2) / (2.0 * s * s))
                    dL = 2.0 * n_exp * Iu * _gint(v1, v2, v0, s)
                    d = _phi(u2, v2, m) - px
                    af += g * d * d * dp * dq
                    ad += (2.0 * (shx + shy) * g + gp * eps * dL) * d * d * dp * dq
            flat += n_exp * n_exp * af * wo[i] * wo[j]
            d1 += n_exp * n_exp * ad * wo[i] * wo[j]
    return flat, d1


def _gl01(n):
    x, w = np.polynomial.legendre.leggauss(int(n))
    return 0.5 * (x + 1), 0.5 * w


def _outer(no, n_exp, eps):
    tmax = math.log(1e3 * eps * n_exp)
    k = max(1, int(no) // 8)
    x8, w8 = np.polynomial.legendre.leggauss(8)
    xs, ws = [], []
    for i in range(k):
        a_, b_ = tmax * i / k, tmax * (i + 1) / k
        t = a_ + (b_ - a_) * 0.5 * (x8 + 1)
        xs.append(1 - np.exp(-t))
        ws.append(np.exp(-t) * (b_ - a_) * 0.5 * w8)
    return np.concatenate(xs), np.concatenate(ws)


def _continuum(m, s, u0, v0, n=1600):
    g = (np.arange(n) + 0.5) / n
    U, V = np.meshgrid(g, g, indexing="ij")
    x = 0.5 * (V - U)
    ph_x = 4 * m * np.exp(m * x) / (1 + np.exp(2 * m * x))
    kin = float(np.sum(-0.25 * ph_x * ph_x)) / n ** 2                 # int phi_u phi_v du dv
    Vp = m * m * (1 - np.cos(4 * np.arctan(np.exp(m * x))))
    sh = np.exp(-((U - u0) ** 2 + (V - v0) ** 2) / (2 * s * s))
    return kin, float(-np.sum(sh * Vp)) / n ** 2


def _pot_d1_grid(m, s, u0, v0, no=192):
    """I3: -int sh V du dv on the outer quadrature nodes (Gauss-Legendre on [0, 1]^2)."""
    x, w = _gl01(no)
    U, V = np.meshgrid(x, x, indexing="ij")
    W = np.outer(w, w)
    xx = 0.5 * (V - U)
    Vp = m * m * (1 - np.cos(4 * np.arctan(np.exp(m * xx))))
    sh = np.exp(-((U - u0) ** 2 + (V - v0) ** 2) / (2 * s * s))
    return float(-np.sum(W * sh * Vp))


def run(p, ctx):
    t0 = time.time()
    eps, m, s, u0, v0 = float(p["eps"]), float(p["m"]), float(p["s"]), float(p["u0"]), float(p["v0"])
    kin_c, target = _continuum(m, s, u0, v0)
    lines = [f"2D light-cone square, eps {eps:g}, kink m {m:g}, bump s {s:g} at ({u0:g}, {v0:g}); target "
             f"int sigma T^a_a sqrt|g| / a = {target:+.5f}; continuum flat S_kin {kin_c:+.5f}"]
    rows = {}
    ok = True
    i2 = True
    flat_last = None
    for ne in p["n_exp"]:
        Rs = []
        for lev in p["levels"]:
            no, nt, nw = lev
            xo, wo = _outer(no, ne, eps)
            xt, wt = _gl01(nt)
            xw, ww = _gl01(nw)
            fl, d1 = _quad(float(ne), eps, m, s, u0, v0, xo, wo, xt, wt, xw, ww, 1e-10 / (eps * ne), 60.0)
            Sk, dSk = -eps ** 2 * fl, -eps ** 2 * d1
            R = dSk / target
            Rs.append(R)
            lines.append(f"  N {ne} level {tuple(lev)}: flat <S_kin> {Sk:+.4f} ({Sk / kin_c:.4f} of the continuum); "
                         f"d<S_kin>/da {dSk:+.5f}; R = {R:+.4f}")
            flat_last = Sk
            ctx.progress(f"[F6.10a] N {ne} level {lev} ({time.time() - t0:.0f}s)")
        delta = abs(Rs[-1] - Rs[-2])
        i2 &= delta <= float(p["tol_delta"])
        R = Rs[-1]
        v = ("PASSED" if abs(R) + delta <= float(p["tol"]) else
             "FAILED" if abs(R) - delta > float(p["fail"]) else "UNDECIDED")
        rows[ne] = {"R": R, "delta": delta, "verdict": v}
        lines.append(f"  N {ne}: R = {R:+.4f}, Delta {delta:.4f}; |R| + Delta = {abs(R) + delta:.4f} -> {v}")
    i1 = abs(flat_last / kin_c - 1) <= float(p["tol_flat"])
    pg = _pot_d1_grid(m, s, u0, v0)
    i3 = abs(pg - target) <= 1e-6 * max(1.0, abs(target))
    inst = {"I1": i1, "I2": i2, "I3": i3}
    lines.append(f"  I1 flat S_kin at N {p['n_exp'][-1]}: {flat_last / kin_c:.4f} of the continuum (tol "
                 f"{p['tol_flat']}) -> {'OK' if i1 else 'FAILS'}; I2 quadrature -> {'OK' if i2 else 'FAILS'}; "
                 f"I3 potential derivative {pg:+.7f} vs {target:+.7f} -> {'OK' if i3 else 'FAILS'}")
    if not all(inst.values()):
        status = "INSTRUMENT_FAILED"
    elif all(r["verdict"] == "PASSED" for r in rows.values()):
        status = "PASSED"
    elif any(r["verdict"] == "FAILED" for r in rows.values()):
        status = "FAILED"
    else:
        status = "UNDECIDED"
    lines.append(f"cost: {time.time() - t0:.0f} s")
    return Result(status, lines, {"rows": {str(k): v for k, v in rows.items()}, "instrument": inst,
                                  "target": target, "kin_continuum": kin_c})
