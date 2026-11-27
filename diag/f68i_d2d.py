"""
F6.8i-dev, step D2d (development, decides nothing): the direct check of the cusp term with the exact needle.
docs/derivations/f68i_dev.md sections 6.2, 6.5 and 7.

Per phase (eps, p, theta) the age term of the needle is (f68f_dev, step 1c; f68i_dev section 6.2)

    J(gamma) = int_0^inf ds g0(s^2 A(c = gamma s)),     gamma = k^2 / (2 l2),

with A(c) the exact needle volume (f68f_1b_exact: the envelopes of the ray families from x and from y). For an analytic
needle the BD moments leave J = 16 gamma^3 s4hat A0^-2 + O(gamma^4) (f68f_1c: 4 mu_3 (k^2/l2)^3, k^2/l2 = 2 gamma).
At a cusp of the first arrival the normal form (section 6.2) predicts, over a patch of the (p, theta) plane around it,

    Q(gamma) = int dp int dtheta/(2 pi) J = C Psi (2 pi)^-1 gamma^(11/7) + O(gamma^3),

with C = -8.8324e-3 (diag/out/f68i_d2b_nf_20261125T101212Z) and Psi = (sig' tau')^-1 N'^(-9/7) from D2c (section 6.5).
This file computes J from the exact needle only; the normal form is used nowhere except to place the quadrature nodes
and to state the prediction.

Modes:
  --check     the needle of this file against f68f_1b_exact.phase_A (one window, the q of EX.QS): must agree to rounding.
  --control   phases without a cusp: J(gamma) / gamma^3 against 16 s4hat A0^-2 from the series code f68f_1c
              (independent code); exponent 3 and the coefficient.

The needle A(c) at many c: the families of rays are integrated once (for the largest c); the window in psi and the
integral are then done per decade of c, each decade with its own window (phase_A uses one window for all q). Guide rays
psi_ray +- GUIDE W(s) widen the window where the valley of Phi is flat (f68i_d2c revision 2.2).

Revision 3 (27 Nov): MAX_WIDEN 5 -> 12 (rows 14, 15 of the eps 12 patch at gamma 6.25e-5 refused the window near the
cusp); a failing patch row is reported in the record instead of stopping the batch. Phases that passed under
revision 2 take the same path (the widening loop stops at the same step), so their values are unchanged.

Revision 2 (26 Nov): the s nodes of J get a panel edge at every s = gap / gamma, gap the action of a later ray above
S_min (rays_row now returns the gaps). Without it the three-ray rows of the patch carried a ripple in theta of ~1.5 %
of J (Gauss nodes across the jump of A''), which the tau quadrature turned into +-4 % of Q (record 20261126T171747Z,
rows at sig_hat <= -5). Single-ray phases (the control, the check) are unchanged.
"""
from __future__ import annotations

import argparse
import math
import os
import time

import numpy as np

from diag import f68f_common as Cm
from diag import f68f_1b as B
from diag import f68f_1b_exact as EX
from diag import f68f_1c as SER
from diag import f68i_d2b_nf as NF
from diag.f68i_d2b_nf import g0

NAME = "f68i_d2d"
PATCH_CODE = None             # set below (code md5 of this file and its imports; the cache key of the patch rows)
C_CUSP = -8.8324e-3          # f68i_d2b_nf_20261125T101212Z (section 6.2)
S_CAP = 16.0                 # s range of the J integral; s^2 A must exceed V_END there
V_END = 45.0                 # g0(45) ~ 1e-15
N_GL = 10
GUIDE = 3.0                  # guide rays psi +- GUIDE W for a refused window (f68i_d2c revision 2.2)
NV_RAYS = 6400               # B.N_V for the ray search (as f68i_d2c revision 2: the three rays near a cusp)
NPSI = 16001                 # points of the psi window per decade of c (EX uses 4001 for all q at once)
MAX_WIDEN = 12               # doublings of the window (EX: 2); revision 3: 5 -> 12 (the quartic valley at the cusp
                             # is ~c^(1/4) wide against the c^(1/2) of W, ~90 W at gamma 6.25e-5)
NV_FAM = 8001                # rays per family (EX.N_V = 4001); the quintic Hermite envelope between rays is the
                             # noise floor of J: ~4e-9 with 4001, ~2e-10 with 8001 and 16001 (eps 8, p 2, 26 Nov)


def s_panels(s_end, breaks=()):
    """Gauss-Legendre nodes on [0, s_end]: graded geometrically towards 0 (A ~ c^(-1/4) at a cusp), then panels of 1.
    breaks (revision 2): the s where c = gamma s reaches the action gap of a later ray; there A gains a (c - gap)_+^2
    term (a jump of A''), and a node set that ignores it gives J a ripple in theta (~1.5 % of J at gamma 2.5e-4 in the
    three-ray region, 20261126 dense profile). Each break is a panel edge, graded geometrically on both sides."""
    pts = [0.0] + [4.0 ** (-j) for j in range(6, 0, -1)] + list(np.arange(1.0, math.ceil(s_end) + 1e-9, 1.0))
    end = pts[-1]
    for b in breaks:
        if 0.0 < b < end:
            pts += [b] + [b + sg * 0.5 * 4.0 ** (-j) for j in range(0, 6) for sg in (-1.0, 1.0)]
    pts = sorted(x for x in pts if 0.0 <= x <= end)
    pts = [x for i, x in enumerate(pts) if i == 0 or x - pts[i - 1] > 1e-12]
    x, w = np.polynomial.legendre.leggauss(N_GL)
    S, W = [], []
    for a, b in zip(pts[:-1], pts[1:]):
        S.append(0.5 * (b - a) * (x + 1) + a)
        W.append(0.5 * (b - a) * w)
    return np.concatenate(S), np.concatenate(W)


def rays_row(eps, p, ths, cmax):
    """For every theta of a row (fixed p): S_min, the rays within cmax of it at the nodes s_m, and the ray count.
    One vectorised call per row (the per-phase calls cost ~17 s in Python overhead)."""
    B.N_V = NV_RAYS
    ths = np.atleast_1d(np.asarray(ths, float))
    roots = B.branches(eps, p, ths, 1000)
    v = np.concatenate(roots)
    tt = np.concatenate([np.full(r.size, t) for r, t in zip(roots, ths)])
    vv = B.newton(v, eps, p, tt, 4000)
    _, _, S = B.end_state(vv, eps, p, tt, 4000)
    xs, _ = np.polynomial.legendre.leggauss(EX.N_S)
    sn = 0.5 * (xs + 1.0)
    sg = np.linspace(0, 1, EX.N_RK + 1)
    out = []
    k = 0
    for r, t in zip(roots, ths):
        Sv, vr = S[k:k + r.size], vv[k:k + r.size]
        k += r.size
        Smin = float(Sv.min())
        near = vr[Sv - Smin < cmax]
        gaps = sorted(float(x) for x in Sv[Sv - Smin < cmax] - Smin)
        paths = B.integrate(near, eps, p, np.full(near.size, t), EX.N_RK)[:, 0, :]
        out.append((Smin, [np.interp(sn, sg, paths[:, j]) for j in range(paths.shape[1])], int(r.size), gaps))
    return out


class Needle:
    """The exact needle at one phase: rays and the two ray families are computed once (for c up to cmax); A(c) is
    then evaluated for any set of c, per decade of c with its own psi window."""

    def __init__(self, eps, p, th, cmax, guide=0.0, pre=None, nv=NV_FAM):
        EX.N_V = nv
        self.eps, self.p, self.th, self.cmax, self.guide = eps, p, th, cmax, guide
        self.Smin, self.rays, self.nray, self.gaps = pre if pre is not None else rays_row(eps, p, [th], cmax)[0]
        xs, wsg = np.polynomial.legendre.leggauss(EX.N_S)
        self.sn = 0.5 * (xs + 1.0)
        self.ws = 0.5 * wsg
        sn = self.sn
        Wtop = (3.0 * math.sqrt(2 * cmax) + 0.1) * np.sqrt(sn * (1 - sn))
        allr = self.rays + ([r + k * guide * Wtop for r in self.rays for k in (-1.0, 1.0)] if guide > 0 else [])
        vmax = eps + 2.0 * math.sqrt(eps) + 2.0 + max((max(abs(r[m]) for r in allr) + 4.0 * Wtop[m])
                                                       / min(sn[m], 1 - sn[m]) for m in range(EX.N_S))
        self.fx = EX.family(eps, p, th, vmax, sn)
        self.fy = EX.family(eps, p, th, vmax, sn, backward=True)
        self.worst = 0.0
        self.widen = 0
        self.shift = 0.0

    def A(self, cs, decades=True, npsi=None, renorm=True, wfloor=0.0):
        """A(c) for every c in cs (c <= cmax). decades=False, npsi=EX.N_PSI, renorm=False, wfloor=0.1 is the window
        and the trapezoid of EX.phase_A (one window for all c)."""
        npsi = npsi or NPSI
        cs = np.asarray(cs, float)
        cref = self.cmax
        if decades:
            k = np.floor(np.log10(cref / cs) + 1e-12).astype(int)     # c in (cmax 10^-(k+1), cmax 10^-k]
        else:
            k = np.zeros(cs.size, int)
        sn, ws, fx, fy = self.sn, self.ws, self.fx, self.fy
        A = np.empty(cs.size)
        for kk in np.unique(k):
            idx = np.nonzero(k == kk)[0]
            cg = cref * 10.0 ** (-kk) if decades else float(cs.max())
            # no 0.1 floor (EX has one): at small c it leaves ~1e3 points across the needle and a trapezoid edge
            # error ~1e-8 (the floor of J grew as gamma fell, control run 20261126T130444Z); the widening loop
            # guards the edges
            Wg = (3.0 * math.sqrt(2 * cg) + wfloor) * np.sqrt(sn * (1 - sn))
            tot = np.zeros(idx.size)
            cc = cs[idx]
            for m in range(EX.N_S):
                s = sn[m]
                lo = min(r[m] for r in self.rays) - self.guide * Wg[m]
                hi = max(r[m] for r in self.rays) + self.guide * Wg[m]
                W = Wg[m]
                for t in range(MAX_WIDEN + 1):
                    x = np.linspace(lo - W, hi + W, npsi)
                    G1 = EX.envelope(fx["psi"][m], fx["S"][m], fx["P"][m], fx["E"][m], fx["dP"][m], x)
                    G2 = EX.envelope(fy["psi"][m], fy["S"][m], fy["P"][m], fy["E"][m], fy["dP"][m], x)
                    Phi = G1 + G2 - self.Smin
                    if not (np.isfinite(Phi[0]) and np.isfinite(Phi[-1])):
                        raise RuntimeError("window outside the families (coverage)")
                    if Phi[0] > cg and Phi[-1] > cg:
                        break
                    W *= 2.0
                    self.widen += 1
                else:
                    raise RuntimeError(f"window: Phi <= c at the edge after widening {MAX_WIDEN} times")
                self.worst = min(self.worst, float(np.nanmin(Phi)))
                if renorm:
                    # per-node zero of Phi: the minimum of G1 + G2 over psi is S_min exactly; the families (RK4,
                    # N_RK steps from x and from y) and S_min (RK4, 4000 steps along the ray) differ by ~1e-12, which
                    # gives A a relative error ~ 2.5e-12 / c, i.e. a floor of J that grows as gamma falls (control
                    # 20261126T130444Z). The parabola through the grid minimum gives the zero to ~h^3.
                    i = int(np.nanargmin(Phi))
                    if 0 < i < Phi.size - 1:
                        y0, y1, y2 = Phi[i - 1], Phi[i], Phi[i + 1]
                        den = y0 - 2 * y1 + y2
                        pmin = y1 - 0.125 * (y2 - y0) ** 2 / den if den > 0 else y1
                    else:
                        pmin = Phi[i]
                    Phi = Phi - pmin
                    self.shift = max(self.shift, abs(float(pmin)))
                dx = x[1] - x[0]
                g = np.maximum(cc[:, None] - Phi[None, :], 0.0) ** 1.5
                g[~np.isfinite(g)] = 0.0
                tot += ws[m] * math.sqrt(s * (1 - s)) * dx * (g.sum(1) - 0.5 * (g[:, 0] + g[:, -1]))
            A[idx] = 32 * math.sqrt(2) / (math.pi * 4.0 * cc ** 2) * tot
        return A


def J_phase(eps, p, th, gammas, guide=0.0, pre=None):
    """J(gamma) for every gamma. s_end per gamma where s^2 A(gamma s) = V_END (two probes), capped at S_CAP."""
    gmax = max(gammas)
    N = Needle(eps, p, th, gmax * S_CAP, guide, pre)
    ends = []
    for g in gammas:
        se = 6.0
        for _ in range(2):
            a = N.A([g * se])[0]
            se = min(S_CAP, 1.1 * math.sqrt(V_END / a))
        ends.append(se)
    nodes = [s_panels(se, [gp / g for gp in N.gaps if gp > 0.0]) for g, se in zip(gammas, ends)]
    cs = np.concatenate([g * sw[0] for g, sw in zip(gammas, nodes)])
    A = N.A(cs)
    J, vend, k = [], [], 0
    for g, (s, w) in zip(gammas, nodes):
        a = A[k:k + s.size]
        k += s.size
        J.append(float(np.sum(w * g0(s * s * a))))
        vend.append(float(s[-1] ** 2 * a[-1]))
    return {"J": J, "s_end": ends, "v_end": vend, "minPhi": N.worst, "nray": N.nray, "widen": N.widen,
            "zero_shift": N.shift}


# ------------------------------------------------------------------------------------------------ modes

def check_job(arg):
    eps, p, th = arg
    qs = list(EX.QS)
    cmax = max(qs) ** 2 / 2
    pre = rays_row(eps, p, [th], cmax)[0]
    EX.N_V = 4001
    ref, _ = EX.phase_A(eps, p, th, pre[0], pre[1], qs=qs)
    cs = np.array([q * q / 2 for q in qs])
    mine = Needle(eps, p, th, cmax, pre=pre, nv=4001).A(cs, decades=False, npsi=EX.N_PSI, renorm=False, wfloor=0.1)
    dec = Needle(eps, p, th, cmax, pre=pre).A(cs)
    r = np.array([ref[q] for q in qs])
    return {"eps": eps, "p": p, "th": th, "q": qs, "ref": r.tolist(), "one_window": mine.tolist(),
            "decades": dec.tolist(), "dev_one": float(np.max(np.abs(mine / r - 1))),
            "dev_dec": float(np.max(np.abs(dec / r - 1)))}


def control_job(arg):
    eps, p, th, gammas = arg
    t = time.time()
    sr = SER.series_coefficients(eps, p, np.array([th]))
    A4 = sr["A"][0]
    pred = 16.0 * SER.s4hat(A4) * A4[0] ** -2
    out = J_phase(eps, p, th, gammas)
    out |= {"eps": eps, "p": p, "th": th, "gammas": gammas, "series_coef": float(pred), "series_A": A4.tolist(),
            "minK": sr["minK"], "seconds": time.time() - t}
    return out


# ------------------------------------------------------------------------------------------------ the patch

# cusp constants from f68i_d2c revision 2 (f68i_dev section 6.5): p_c, sig', tau', N', Psi
CUSPS = {12.0: (2.449046, 0.9703, 6.8218, 1.1963, 0.11998),
         20.0: (4.3813227, 1.7901, 6.0318, 1.3264, 0.06441)}
CUSPS_ALL = {9.9: 0.3135853, 12.0: 2.449046, 20.0: 4.3813227, 35.0: 5.427158}
COUNT_P = [0.25, 0.5, 0.9, 0.98, 1.02, 1.1, 1.5, 2.5]
NS, NT = 3, 6                # Gauss nodes per sig panel and per tau segment (f68i_d2b_nf layout; on the normal form
                             # this quadrature is 0.56 % off its converged value, the same factor at every gamma)


def patch_nodes(eps, gam, ns=NS, nt=NT):
    """The nodes of f68i_d2b_nf.plane (canonical sig, tau at gamma = 1; tau >= 0), scaled to gamma_nf = gamma / sqrt(N')
    and mapped to the physical plane: p = p_c + sig / sig', theta = -p/2 + tau / tau'. Returns the rows."""
    pc, sp, tp, Np, _ = CUSPS[eps]
    gn = gam / math.sqrt(Np)
    a, b = gn ** (4 / 7), gn ** (6 / 7)
    brk = sorted(set([-x for x in NF.SIG_BRK] + NF.SIG_BRK))
    x, w = np.polynomial.legendre.leggauss(ns)
    rows = []
    for lo, hi in zip(brk[:-1], brk[1:]):
        for xi, wi in zip(0.5 * (hi - lo) * (x + 1) + lo, 0.5 * (hi - lo) * w):
            T, W, _ = NF.tau_nodes(float(xi), nt, 1.0)
            rows.append({"sig_hat": float(xi), "w_sig": float(wi), "tau_hat": T.tolist(), "w_tau": W.tolist(),
                         "p": pc + xi * a / sp, "sig": xi * a, "tau": (T * b).tolist(),
                         "jac": a * b / (sp * tp)})
    return rows


def row_job(eps, gam, row_index):
    """J at the nodes of one row (fixed p): exact needle (guide rays if the plain window is refused) and the normal form
    at the same (sig, tau). Cached per row."""
    pc, sp, tp, Np, _ = CUSPS[eps]
    row = patch_nodes(eps, gam)[row_index]
    p = row["p"]
    ths = [(-p / 2 + t / tp) % (2 * math.pi) for t in row["tau"]]
    t0 = time.time()
    pre = rays_row(eps, p, ths, gam * S_CAP)
    Jx, guides, nray = [], [], []
    for th, pr in zip(ths, pre):
        try:
            r = J_phase(eps, p, th, [gam], pre=pr)
            guides.append(0)
        except RuntimeError:
            r = J_phase(eps, p, th, [gam], guide=GUIDE, pre=pr)
            guides.append(1)
        Jx.append(r["J"][0])
        nray.append(r["nray"])
    gn = gam / math.sqrt(Np)
    Jn = [Np ** -0.5 * NF.J(row["sig"], t, gn) for t in row["tau"]]
    return {"J": Jx, "J_nf": Jn, "guide": guides, "nray": nray, "seconds": time.time() - t0}


def symmetry_job(arg):
    """J at theta_c + t and theta_c - t (the patch uses tau >= 0 only)."""
    eps, gam, row_index, k = arg
    pc, sp, tp, Np, _ = CUSPS[eps]
    row = patch_nodes(eps, gam)[row_index]
    p = row["p"]
    t = row["tau"][k] / tp
    out = []
    try:
        for th in ((-p / 2 + t) % (2 * math.pi), (-p / 2 - t) % (2 * math.pi)):
            out.append(J_phase(eps, p, th, [gam])["J"][0])
    except RuntimeError as err:
        return {"eps": eps, "gam": gam, "p": p, "dtheta": t, "error": str(err)}
    return {"eps": eps, "gam": gam, "p": p, "dtheta": t, "J_plus": out[0], "J_minus": out[1]}


def patch_sum(rows, key, sig_max=16.0):
    """(1/2pi) int dp int dtheta J over the nodes (tau >= 0 doubled), plus the sig tails beyond |sig_hat| = sig_max from a
    power fit of the row integrals on |sig_hat| >= sig_max / 2 (as f68i_d2b_nf.plane)."""
    core = 0.0
    ints = []
    for r in rows:
        ri = 2.0 * float(np.dot(r["w_tau"], r[key]))
        ints.append((r["sig_hat"], ri))
        core += r["w_sig"] * ri * r["jac"]
    jac = rows[0]["jac"]
    tails = {}
    for side in (-1, 1):
        pts = [(abs(sh), v) for sh, v in ints if side * sh >= sig_max / 2]
        ys = np.array([v for _, v in pts])
        if len(pts) >= 2 and (np.all(ys > 0) or np.all(ys < 0)):
            fit = np.polyfit(np.log([q for q, _ in pts]), np.log(np.abs(ys)), 1)
            pw = -fit[0]
            Bc = np.sign(ys[0]) * math.exp(fit[1])
            tails[side] = {"p": float(pw), "tail": float(Bc * sig_max ** (1 - pw) / (pw - 1) * jac)}
        else:
            tails[side] = {"p": float("nan"), "tail": 0.0, "note": "no fit (sign change or too few rows)"}
    total = core + sum(t["tail"] for t in tails.values())
    return {"core": core / (2 * math.pi), "tails": tails, "total": total / (2 * math.pi),
            "abs": sum(r["w_sig"] * 2.0 * float(np.dot(r["w_tau"], np.abs(r[key]))) * r["jac"] for r in rows) / (2 * math.pi)}


PATCH_CODE = Cm.code_md5([os.path.abspath(__file__), B.__file__, EX.__file__, NF.__file__])


def patch_job(arg):
    """One cached row; a failing row is reported and not cached, so that the other rows of the batch finish (rev 3)."""
    eps, gam, i = arg
    try:
        return Cm.cached(NAME + "_patch", PATCH_CODE, row_job, (eps, gam, i))
    except RuntimeError as err:
        return {"error": str(err), "args": [eps, gam, i]}


# ------------------------------------------------------------------------------------------------ the cusp count

N_TH_COUNT = 720             # phases of the full theta scan per p
NV_COUNT = 800               # B.N_V for the scan


def first_arrival_scan(eps, p, n_th=N_TH_COUNT):
    """Over theta in [0, 2 pi): the slope v and action S of the first-arrival ray, the gap to the next ray, the number
    of rays. Cut loci (ties) are where the first-arrival ray jumps; each jump is refined by bisection in theta, and the
    gap of the two rays at the refined point must vanish."""
    B.N_V = NV_COUNT
    ths = 2 * math.pi * np.arange(n_th) / n_th

    def lowest(thv):
        roots = B.branches(eps, p, thv, 1000)
        v = np.concatenate(roots)
        tt = np.concatenate([np.full(r.size, t) for r, t in zip(roots, thv)])
        vv = B.newton(v, eps, p, tt, 4000)
        _, _, S = B.end_state(vv, eps, p, tt, 4000)
        out, k = [], 0
        for r in roots:
            Sv, vr = S[k:k + r.size], vv[k:k + r.size]
            k += r.size
            o = np.argsort(Sv)
            gap = float(Sv[o[1]] - Sv[o[0]]) if r.size > 1 else float("inf")
            out.append((float(vr[o[0]]), float(Sv[o[0]]), gap, int(r.size), vr[o[:2]].tolist()))
        return out

    res = []
    for k in range(0, n_th, 90):
        res += lowest(ths[k:k + 90])
    v = np.array([r[0] for r in res])
    dv = np.abs(np.diff(np.r_[v, v[0]]))
    med = float(np.median(dv))
    jumps = np.nonzero(dv > max(20 * med, 1e-3))[0]
    ties = []
    for j in jumps:
        a, b = ths[j], ths[j] + 2 * math.pi / n_th
        va = res[j][0]
        for _ in range(30):
            m = 0.5 * (a + b)
            vm = lowest(np.array([m % (2 * math.pi)]))[0][0]
            if abs(vm - va) < 0.5 * dv[j]:
                a = m
            else:
                b = m
        r = lowest(np.array([0.5 * (a + b) % (2 * math.pi)]))[0]
        ties.append({"theta": 0.5 * (a + b) % (2 * math.pi), "dv": float(dv[j]), "gap_at_tie": r[2],
                     "theta_minus_sym": float(((0.5 * (a + b) + p / 2 + math.pi) % (2 * math.pi)) - math.pi)})
    return {"eps": eps, "p": p, "n_th": n_th, "max_rays": max(r[3] for r in res), "median_dv": med,
            "ties": ties, "Smin": [r[1] for r in res]}


def count_job(arg):
    eps, p = arg
    return first_arrival_scan(eps, p)


def run(jobs, fn, workers):
    if workers > 1:
        from multiprocessing import Pool
        with Pool(min(workers, len(jobs))) as pool:
            return list(pool.map(fn, jobs))
    return [fn(j) for j in jobs]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", type=str, default="", help="eps:p:theta,... (needle against EX.phase_A)")
    ap.add_argument("--control", type=str, default="", help="eps:p:theta,... (J against the series)")
    ap.add_argument("--gammas", type=str, default="0.008,0.004,0.002,0.001")
    ap.add_argument("--patch", type=float, default=0.0, help="eps of the cusp patch (a key of CUSPS)")
    ap.add_argument("--pgammas", type=str, default="0.004,0.001,0.00025")
    ap.add_argument("--rows", type=str, default="", help="only these row indices (a partial run fills the cache)")
    ap.add_argument("--count", type=str, default="", help="eps,... : full theta scan of the first arrival at p/p_c in "
                                                           "COUNT_P and at -p_c/2 (mirror)")
    ap.add_argument("--workers", type=int, default=0)
    a = ap.parse_args()
    t0 = time.time()
    code = Cm.code_md5([os.path.abspath(__file__), B.__file__, EX.__file__, SER.__file__])
    gammas = [float(x) for x in a.gammas.split(",")]
    L = ["F6.8i-dev D2d: the age term J from the exact needle (development, decides nothing; "
         "docs/derivations/f68i_dev.md section 7)", ""]
    data = {}
    if a.check:
        pts = [tuple(float(v) for v in x.split(":")) for x in a.check.split(",")]
        res = run(pts, check_job, a.workers)
        data["check"] = res
        L.append("needle of this file against f68f_1b_exact.phase_A (q of EX.QS):")
        for r in res:
            L.append(f"  eps {r['eps']:g} p {r['p']:g} theta {r['th']:.6f}: one window max rel dev {r['dev_one']:.1e}; "
                     f"per-decade windows {r['dev_dec']:.1e}")
    if a.control:
        pts = [tuple(float(v) for v in x.split(":")) + (gammas,) for x in a.control.split(",")]
        res = run(pts, control_job, a.workers)
        data["control"] = res
        L.append("control without a cusp: J(gamma) against 16 gamma^3 s4hat A0^-2 (series, f68f_1c):")
        for r in res:
            L.append(f"  eps {r['eps']:g} p {r['p']:g} theta {r['th']:.6f} (min K {r['minK']:.3f}, {r['nray']} rays, "
                     f"{r['seconds']:.0f} s): series coefficient {r['series_coef']:+.6e}")
            g = np.array(r["gammas"])
            J = np.array(r["J"])
            for i in range(len(g)):
                L.append(f"    gamma {g[i]:.5f}: J {J[i]:+.6e}  J/gamma^3 {J[i] / g[i] ** 3:+.6e}  "
                         f"ratio to series {J[i] / g[i] ** 3 / r['series_coef']:.6f}  (s^2 A at s_end {r['v_end'][i]:.0f})")
            sl = np.diff(np.log(np.abs(J))) / np.diff(np.log(g))
            L.append("    local exponents d ln|J| / d ln gamma: " + ", ".join(f"{x:.4f}" for x in sl))
    if a.count:
        jobs = []
        for e in [float(x) for x in a.count.split(",")]:
            pc = CUSPS_ALL[e]
            jobs += [(e, f * pc) for f in COUNT_P] + [(e, -0.5 * pc)]
        res = run(jobs, count_job, a.workers)
        data["count"] = [{k: v for k, v in r.items() if k != "Smin"} for r in res]
        L.append(f"cusp count: first-arrival ties over theta (scan of {N_TH_COUNT} phases, B.N_V {NV_COUNT}):")
        for r in res:
            pc = CUSPS_ALL[r["eps"]]
            tl = "; ".join(f"theta {t['theta']:.5f} (theta + p/2 {t['theta_minus_sym']:+.2e}, gap {t['gap_at_tie']:.1e})"
                           for t in r["ties"])
            L.append(f"  eps {r['eps']:g} p {r['p']:+.5f} (p/p_c {r['p'] / pc:+.2f}), up to {r['max_rays']} rays: "
                     f"{len(r['ties'])} tie(s)" + (f": {tl}" if tl else ""))
        for e in sorted(set(r["eps"] for r in res)):
            a1 = [r for r in res if r["eps"] == e and abs(r["p"] / CUSPS_ALL[e] - 0.5) < 1e-9][0]
            a2 = [r for r in res if r["eps"] == e and abs(r["p"] / CUSPS_ALL[e] + 0.5) < 1e-9][0]
            S1 = np.array(a1["Smin"])
            S2 = np.array(a2["Smin"])
            n = S1.size
            mirror = S2[(-np.arange(n)) % n]          # theta -> -theta
            L.append(f"  eps {e:g}: mirror S_min(-p, -theta) against S_min(p, theta) at p = p_c/2: "
                     f"max |diff| {float(np.max(np.abs(mirror - S1))):.1e}")
    if a.patch:
        eps = a.patch
        pc, sp, tp, Np, Psi = CUSPS[eps]
        pg = [float(x) for x in a.pgammas.split(",")]
        nrow = len(patch_nodes(eps, pg[0]))
        sel = [int(x) for x in a.rows.split(",")] if a.rows else list(range(nrow))
        jobs = [(eps, g, i) for g in pg for i in sel]
        res = run(jobs, patch_job, a.workers)
        sym = run([(eps, pg[0], nrow // 2 - 2, 3), (eps, pg[0], nrow // 2 - 2, 8), (eps, pg[-1], nrow // 2 + 1, 4)],
                  symmetry_job, a.workers)
        data["patch"] = {"eps": eps, "cusp": CUSPS[eps], "gammas": pg, "rows": sel, "jobs": res, "symmetry": sym}
        L.append(f"patch around the cusp at eps {eps:g} (p_c {pc}, sig' {sp}, tau' {tp}, N' {Np}, Psi {Psi}); "
                 f"{nrow} rows of the f68i_d2b_nf layout (ns {NS}, nt {NT}), tau >= 0 doubled:")
        for sy in sym:
            if "error" in sy:
                L.append(f"  symmetry: gamma {sy['gam']:g} p {sy['p']:.5f} dtheta {sy['dtheta']:.3e}: FAILED {sy['error']}")
                continue
            L.append(f"  symmetry: gamma {sy['gam']:g} p {sy['p']:.5f} dtheta {sy['dtheta']:.3e}: J(+) {sy['J_plus']:+.9e} "
                     f"J(-) {sy['J_minus']:+.9e}")
        failed = [jr for jr in res if "error" in jr]
        for jr in failed:
            L.append(f"  row failed: eps {jr['args'][0]:g} gamma {jr['args'][1]:g} row {jr['args'][2]}: {jr['error']}")
        if failed:
            L.append(f"  {len(failed)} rows failed (not cached); no sums")
        elif len(sel) == nrow:
            out = []
            for g in pg:
                rows = patch_nodes(eps, g)
                for i, r in enumerate(rows):
                    jr = res[pg.index(g) * nrow + i]
                    r["J"], r["J_nf"] = jr["J"], jr["J_nf"]
                ex = patch_sum(rows, "J")
                nf = patch_sum(rows, "J_nf")
                pred = C_CUSP * Psi * g ** (11 / 7) / (2 * math.pi)
                ng = sum(sum(jr["guide"]) for jr in res[pg.index(g) * nrow:(pg.index(g) + 1) * nrow])
                out.append({"gamma": g, "exact": ex, "nf_same_nodes": nf, "prediction": pred, "guided": ng})
                L.append(f"  gamma {g:g}: Q exact {ex['total']:+.6e} (core {ex['core']:+.6e}, sig tails "
                         f"{ex['tails'][-1]['tail'] / (2 * math.pi):+.2e} / {ex['tails'][1]['tail'] / (2 * math.pi):+.2e}; "
                         f"int|J| {ex['abs']:.3e}; {ng} phases with guide rays)")
                L.append(f"      prediction C Psi gamma^(11/7) / 2pi {pred:+.6e}: ratio {ex['total'] / pred:.5f}; "
                         f"normal form on the same nodes {nf['total']:+.6e}: ratio exact / nf {ex['total'] / nf['total']:.5f}")
            for i in range(len(out) - 1):
                e1 = math.log(out[i + 1]["exact"]["total"] / out[i]["exact"]["total"]) / math.log(pg[i + 1] / pg[i])
                L.append(f"  exponent between gamma {pg[i]:g} and {pg[i + 1]:g}: {e1:.4f} (11/7 = {11 / 7:.4f}; smooth: 3)")
            data["patch"]["summary"] = out
        else:
            L.append(f"  partial run: {len(sel)} of {nrow} rows cached; no sums")
    L.append(f"total {time.time() - t0:.0f} s")
    print("\n".join(L))
    stem = Cm.save(NAME, vars(a), code, data, L, t0)
    print("saved", stem + ".json")


if __name__ == "__main__":
    main()
