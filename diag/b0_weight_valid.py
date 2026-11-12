"""
B0, development stage 5: validation of the calibrated Green weight on pyramid sprinklings (decision 17; decides nothing)
=======================================================================================================================
The weight w(L; rho) = rho^(1/3) w~(L), w~(1) = a, w~(L >= 2) = c / (L + beta), with (a, c, beta) from
diag/b0_weight_calib.py (read from its newest output), is checked on configurations it was not derived from:
  1. pairs of real pyramid sprinklings, L by the DP over the links (diag/b0_chain_alg.py), at three densities: the ratio
     E[w(L)] / E[G(tau)] per bin of s = rho^(1/3) tau, compared with the calibration curve (density scaling), and split
     by the rapidity of the pair (boost independence; the interval is never cut in a causally convex region);
  2. the Green integral I(x) = int G J d^3y of prescribed sources against its exact value, discrete
     I_hat(x) = sum over y < x of w(L_xy) J(y) / rho. Sources:
       homogeneous J = 1 (I = t^2 / 2);
       walls J = sech^2((n.r - u0) / l), l = 1 and 0.5, normals at 0 and 35 degrees (exact in closed form);
       a localised blob J = exp(-|r - r0|^2 / l^2), l = 1 (exact by quadrature);
       outside the declared class (diagnostic only): a wall of width 0.1 and a step J = theta(n.r).
     Three variants on the same pairs: coordinate tau with the cap 0.68 rho^(-1/3) (reference, development only),
     intrinsic A (single alpha_3 = 1.2: G = alpha_3 rho^(1/3) / (2 pi L)), intrinsic B (the table).
     Reported: the ratio of means sum I_hat / sum I over the targets, with a SE over sprinklings, and the rms of the
     per-target relative error.
Coordinates enter to build the relation of the flat sprinkling, to evaluate the prescribed sources and the exact
integrals, and (for the coordinate variant and the binning) the true tau.

    python -m diag.b0_weight_valid --seeds 8 --workers 1
"""

from __future__ import annotations

import argparse
import glob
import json
import math
import os
import time

import numpy as np

from core.order.sprinkle import derive_seed
from core.stats.parallel import parallel_map, default_workers
from phases.F3_structure import field1d as F
from diag.b0_green2p1 import sprinkle_pyramid
from diag.b0_chain_alg import links3

_nb = F._nb
ALPHA3 = 1.2
SBIN = 0.25
NSB = 120
ETA_EDGES = (0.5, 1.5)

if _nb is not None:
    @_nb.njit(cache=True)
    def children(ptr, idx, N):
        cnt = np.zeros(N + 1, np.int64)
        for a in range(idx.size):
            cnt[idx[a] + 1] += 1
        cptr = np.cumsum(cnt)
        cidx = np.empty(idx.size, np.int64)
        fill = cptr[:-1].copy()
        for v in range(N):
            for a in range(ptr[v], ptr[v + 1]):
                u = idx[a]
                cidx[fill[u]] = v
                fill[u] += 1
        return cptr, cidx

    @_nb.njit(cache=True)
    def target_pass(xi, t, x, y, cptr, cidx, J, rho, a, c, beta, tau_c, sw, sg, cn):
        """Backward DP from target xi: L(u) = longest chain u -> xi. Returns the sums of w J / rho for every source row
        of J and the three variants (coordinate, A, B), shape (nsrc, 3); accumulates the pair bins sw, sg, cn."""
        N = t.size
        L = np.full(xi + 1, -1, np.int64)
        L[xi] = 0
        r13 = rho ** (1.0 / 3.0)
        ns = J.shape[0]
        out = np.zeros((ns, 3))
        inv2pi = 1.0 / (2.0 * math.pi)
        for u in range(xi - 1, -1, -1):
            b = -1
            for q in range(cptr[u], cptr[u + 1]):
                v = cidx[q]
                if v <= xi and L[v] >= 0 and L[v] + 1 > b:
                    b = L[v] + 1
            L[u] = b
            if b < 1:
                continue
            dt = t[xi] - t[u]
            dx = x[xi] - x[u]
            dy = y[xi] - y[u]
            tau = math.sqrt(max(dt * dt - dx * dx - dy * dy, 1e-300))
            wB = a if b == 1 else c / (b + beta)
            wB *= r13
            wA = ALPHA3 * r13 * inv2pi / b
            wC = inv2pi / max(tau, tau_c)
            for k in range(ns):
                jv = J[k, u] / rho
                out[k, 0] += wC * jv
                out[k, 1] += wA * jv
                out[k, 2] += wB * jv
            s = r13 * tau
            sb = int(s / SBIN)
            if sb < NSB:
                eta = math.acosh(max(dt / tau, 1.0))
                e = 0 if eta < ETA_EDGES[0] else (1 if eta < ETA_EDGES[1] else 2)
                sw[e, sb] += wB / r13
                sg[e, sb] += inv2pi / s
                cn[e, sb] += 1
        return out


def logcosh(z):
    z = np.abs(z)
    return z + np.log1p(np.exp(-2 * z)) - math.log(2.0)


def exact_wall(t, d, l):
    return 0.5 * l * l * (logcosh((d + t) / l) + logcosh((d - t) / l) - 2 * logcosh(d / l))


def exact_step(t, d):
    # 1/2 int_0^t |[d - t', d + t'] intersect [0, inf)| dt'
    tp = np.linspace(0, 1, 2001)[None, :] * t[:, None]
    f = np.clip(d[:, None] + tp, 0, None) - np.clip(d[:, None] - tp, 0, None)
    trap = getattr(np, "trapezoid", None) or getattr(np, "trapz")
    return 0.5 * trap(f, tp, axis=1)


_U, _WU = np.polynomial.legendre.leggauss(96)
_U = 0.5 * (_U + 1.0)
_WU = 0.5 * _WU
_KU = _U * np.arccosh(1.0 / _U)
_PH = np.linspace(0, 2 * math.pi, 129)[:-1]


def exact_blob(t, px, py, cx, cy, l):
    """(1/2pi) int_{|r'-r|<t} J(r') arccosh(t/|r'-r|) d^2r' = (t^2/2pi) int_0^1 u arccosh(1/u) du int dphi J."""
    out = np.empty(t.size)
    for i in range(t.size):
        r = t[i] * _U[:, None]
        X = px[i] + r * np.cos(_PH)[None, :] - cx
        Y = py[i] + r * np.sin(_PH)[None, :] - cy
        Jv = np.exp(-(X * X + Y * Y) / (l * l))
        out[i] = t[i] ** 2 / (2 * math.pi) * np.sum(_WU * _KU * Jv.mean(axis=1)) * 2 * math.pi
    return out


def sources(t, x, y):
    def n(deg):
        th = math.radians(deg)
        return math.cos(th), math.sin(th)
    S = []
    S.append(("homogeneous", np.ones_like(t), lambda tt, px, py: tt ** 2 / 2, True))
    for l, deg in ((1.0, 0.0), (1.0, 35.0), (0.5, 0.0), (0.5, 35.0)):
        nx, ny = n(deg)
        S.append((f"wall l {l:g} at {deg:g} deg", 1.0 / np.cosh((nx * x + ny * y) / l) ** 2,
                  (lambda l_, nx_, ny_: (lambda tt, px, py: exact_wall(tt, nx_ * px + ny_ * py, l_)))(l, nx, ny), True))
    S.append(("blob l 1", np.exp(-(x ** 2 + y ** 2)), lambda tt, px, py: exact_blob(tt, px, py, 0.0, 0.0, 1.0), True))
    nx, ny = n(35.0)
    S.append(("wall l 0.1 at 35 deg (outside the class)", 1.0 / np.cosh((nx * x + ny * y) / 0.1) ** 2,
              lambda tt, px, py: exact_wall(tt, nx * px + ny * py, 0.1), False))
    S.append(("step at 35 deg (outside the class)", (nx * x + ny * y > 0).astype(float),
              lambda tt, px, py: exact_step(tt, nx * px + ny * py), False))
    return S


def task(rho, seed, T, ntarg, a, c, beta):
    X = T + 2.0
    t, x, y = sprinkle_pyramid(rho, X, X, T, seed)
    lam = 25.0
    tau_cut = (12.0 * lam / (math.pi * rho)) ** (1.0 / 3.0)
    t0 = time.time()
    ptr, idx = links3(t, x, y, tau_cut)
    cptr, cidx = children(ptr, idx, t.size)
    tl = time.time() - t0
    S = sources(t, x, y)
    J = np.vstack([s[1] for s in S])
    rng = np.random.default_rng(derive_seed(seed, "targets"))
    cand = np.nonzero((t > T - 1.5) & (t < T - 0.5) & (np.abs(x) < 1.0) & (np.abs(y) < 1.0))[0]
    tg = np.sort(rng.choice(cand, min(ntarg, cand.size), replace=False))
    sw = np.zeros((3, NSB))
    sg = np.zeros((3, NSB))
    cn = np.zeros((3, NSB))
    Ih = np.zeros((tg.size, len(S), 3))
    t1 = time.time()
    for i, xi in enumerate(tg):
        Ih[i] = target_pass(int(xi), t, x, y, cptr, cidx, J, float(rho), a, c, beta, 0.68 * rho ** (-1 / 3), sw, sg, cn)
    td = time.time() - t1
    Ie = np.column_stack([s[2](t[tg], x[tg], y[tg]) for s in S])
    return {"rho": rho, "seed": seed, "N": int(t.size), "links": int(idx.size), "sec_links": tl,
            "sec_dp_per_target": td / max(tg.size, 1), "names": [s[0] for s in S], "in_class": [s[3] for s in S],
            "Ih": Ih.tolist(), "Ie": Ie.tolist(), "sw": sw.tolist(), "sg": sg.tolist(), "cn": cn.tolist()}


def newest_calib():
    fs = sorted(glob.glob(os.path.join("diag", "out", "b0_weight_calib_*.json")))
    if not fs:
        raise SystemExit("run diag.b0_weight_calib first")
    d = json.load(open(fs[-1]))
    return d["a"], d["c"], d["beta"], fs[-1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=8)
    ap.add_argument("--rhos", type=str, default="50,100,200")
    ap.add_argument("--T", type=float, default=4.0)
    ap.add_argument("--targets", type=int, default=150)
    ap.add_argument("--workers", type=int, default=0)
    a_ = ap.parse_args()
    if _nb is None:
        raise SystemExit("numba required")
    a, c, beta, fn_cal = newest_calib()
    rhos = [float(r) for r in a_.rhos.split(",")]
    jobs = [(rho, derive_seed(20261111, "B0w-valid", rho, k), a_.T, a_.targets, a, c, beta)
            for rho in rhos for k in range(a_.seeds)]
    jobs.sort(key=lambda j: -j[0])
    t0 = time.time()
    res = parallel_map(task, jobs, default_workers(a_.workers), print,
                       lambda i, r, d, n: f"[valid] {d}/{n}: rho {r['rho']:g} N {r['N']} ({time.time() - t0:.0f}s)")
    L_ = [f"calibration {fn_cal}: a {a:.4f}, c {c:.4f}, beta {beta:.3f}; pyramid |x|, |y| <= T + 2 - t, T {a_.T:g}; "
          f"targets t in [T-1.5, T-0.5], |x|, |y| < 1"]
    for rho in rhos:
        R = [r for r in res if r["rho"] == rho]
        L_.append(f"=== rho {rho:g} ({len(R)} sprinklings, N ~{int(np.mean([r['N'] for r in R]))}, links per event "
                  f"{np.mean([r['links'] / r['N'] for r in R]):.1f}, links {np.mean([r['sec_links'] for r in R]):.1f} s, "
                  f"DP per target {1e3 * np.mean([r['sec_dp_per_target'] for r in R]):.1f} ms)")
        sw = np.sum([r["sw"] for r in R], axis=0)
        sg = np.sum([r["sg"] for r in R], axis=0)
        cn = np.sum([r["cn"] for r in R], axis=0)
        L_.append("  pairs: E[w~] / E[G] per bin of s (all rapidities | eta < 0.5 | 0.5-1.5 | > 1.5)")
        for lo, hi in ((0, 0.5), (0.5, 1), (1, 1.5), (1.5, 2), (2, 2.5), (2.5, 3), (3, 4), (4, 6), (6, 10), (10, 30)):
            i0, i1 = int(lo / SBIN), min(int(hi / SBIN), NSB)
            tot = sw[:, i0:i1].sum() / sg[:, i0:i1].sum() if cn[:, i0:i1].sum() > 0 else float("nan")
            parts = [sw[e, i0:i1].sum() / sg[e, i0:i1].sum() if cn[e, i0:i1].sum() > 100 else float("nan")
                     for e in range(3)]
            L_.append(f"    s [{lo:g}, {hi:g}): {tot:.3f} | " + " | ".join(f"{p:.3f}" for p in parts)
                      + f"  (pairs {int(cn[:, i0:i1].sum())})")
        names = R[0]["names"]
        L_.append("  Green integral, ratio of means sum I_hat / sum I (SE over sprinklings); rms per-target rel. error"
                  " -- coordinate | intrinsic A | intrinsic B")
        for k, nm in enumerate(names):
            parts = []
            for v in range(3):
                rat = [np.sum(np.array(r["Ih"])[:, k, v]) / np.sum(np.array(r["Ie"])[:, k]) for r in R]
                rel = np.concatenate([np.array(r["Ih"])[:, k, v] / np.array(r["Ie"])[:, k] - 1 for r in R])
                parts.append(f"{np.mean(rat):.4f} +- {np.std(rat, ddof=1) / math.sqrt(len(rat)):.4f} "
                             f"(rms {np.sqrt(np.mean(rel ** 2)):.3f})")
            L_.append(f"    {nm}: " + " | ".join(parts))
    print("\n".join(L_))
    os.makedirs(os.path.join("diag", "out"), exist_ok=True)
    fn = os.path.join("diag", "out", f"b0_weight_valid_{time.strftime('%Y%m%dT%H%M%S', time.gmtime())}.json")
    with open(fn, "w") as fh:
        json.dump({"args": vars(a_), "calib": fn_cal, "summary": L_, "results": res}, fh)
    print(f"raw records: {fn} ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
