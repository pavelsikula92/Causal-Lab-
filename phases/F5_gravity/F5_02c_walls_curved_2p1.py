"""
F5.2c -- Walls of a field in a prescribed curved 2+1 order of the theory, intrinsic instrument (direction B, B2)
===============================================================================================================

Why this test exists (12 Nov 2026; card docs/B2_card_draft.md, revision 2, frozen by decision 22)
-------------------------------------------------------------------------------------------------
B1 (F3.15, PASSED) showed that a domain wall of a field on a flat 2+1 order, evolved by the Green weight read from the
order alone, is an inertial body. B2 puts walls into a curved order of the theory, the exact sprinkling of the metric
of the theory ds^2 = -Pi^2 dt^2 + Pi^-2 (dx^2 + dy^2) (narrower cones and more events; decisions 12, 16), and asks
whether the same rule, reading only the order, moves them as the continuum field theory in that metric does.
The continuum makes a sharp prediction (docs/derivations/walls_in_curved_order.md): the passive weight of a wall
uniform along y is E - sigma (energy minus tension). The kink (sigma = E) does not fall at all, it is an exact static
solution at every position; the breather wall (omega 0.7, sigma/E = 0.22) falls; a launched kink keeps its proper
speed, so its coordinate speed follows the cones (Pi^2 v). The order must give all three: no spurious force from the
gradient of the density of events (1.7x across the domain), and the right fall of an object with little tension.

Model (frozen)
--------------
Order: Poisson sprinkling of density rho per proper volume in the domain of dependence of the rectangle
-xl <= x <= xr, |y| <= Y of the metric with Pi^2 = (x + ell)/ell; the relation in closed form (optical metric = the
hyperbolic plane of radius ell): u < v iff t_v - t_u > d_h, cosh(d_h/ell) = 1 + ((x_u - x_v)^2 + (y_u - y_v)^2) /
(2 (x_u + ell)(x_v + ell)); links by the shell method with the exact relation (diag/b2_curved.py). The order is
prescribed (B2), not grown (B3). Flat control: ell = inf, the same rectangle.
Field rule exactly as F3.15: phi(v) = phi_0(v) - (1/rho) sum over u < v of w(L_uv) mu sin phi(u), mu = 1, the chain
weight of decision 17 (variant B, frozen constants), L the longest chain over the links.
phi_0: the free wave (mu = 0) of the same initial data in the same metric, from the y-independent continuum solver
(diag/b2_curved.Continuum, leapfrog dx 0.005); the initial data are the only place where coordinates enter the dynamics.
Initial data at x0 = 0 (Pi = 1): kink 4 atan(exp(gamma s)) launched at proper velocity v (phi_t = -Pi gamma v 2/cosh),
or the breather 4 atan((k/omega)/cosh(k s)) at maximal amplitude (phi_t = 0), k = sqrt(1 - omega^2); s the proper
distance from x0 along x.
Arms: on each curved sprinkling: kink_rest, kink_v0.3, br0.7 (the same events); on the first n_ctrl flat sprinklings:
br0.7 (the breather on the flat order; continuum: no motion).

Readouts (frozen; coordinates, declared). Band |y| < 2; slabs of width 1, centres 0.5 ... T - 1.5 (slab rule of B1).
Kink: the fit of field1d.fit_kink (centre, width, Q), guess = the previous centre or the continuum centre. Breather:
the centre xb of a fit of the profile 4 atan(A(t)/cosh(k (x - xb))), A quadratic in t over the slab, over |x| < 5
(fit_breather); reported: the centroid c of the noise-debiased coordinate energy density over |x| < 7 (as F3.15).
The flat control uses the rectangle -16 <= x <= 22 (its domain of dependence closes faster). Every readout is repeated on the exact continuum field at the same events (paired, suffix _ex).
  D_K   = [xc(last) - xc(first)] - [xc_ex(last) - xc_ex(first)]   (kink at rest; continuum 0, a particle -0.44)
  dv_L  = slope(xc) - slope(xc_ex)                                (launched kink; slopes over the read slabs)
  cones = slope(xc) - v                                           (the launched kink against flat cones; continuum +0.013)
  r_B   = sum_t d(t) d_ex(t) / sum_t d_ex(t)^2, d = xb - xb(first) (breather fall projected on the continuum fall)
  D_0   = xb(last) - xb(first) on the flat order, divided by the curved continuum fall d_ex(last) (no fall without
          curvature)

Criteria (card revision 2, section 4; SE = sd / sqrt(n) over sprinklings)
INSTRUMENT_FAILED if any of: I1 the links differ from the brute-force transitive reduction (curved or flat, small
  domain); I4 linear extension t + 0.2 x or the smaller rectangle differ by >= 1e-12; I0 |N / (rho V_proper) - 1| + 2 SE
  >= 0.01; the kink at rest |Q - Q_ex| + 2 SE >= 0.01 or |w/w_ex - 1| + 2 SE >= 0.05; the flat breather |D_0| + 2 SE
  >= 0.25.
FAILED if any of: |D_K| - 2 SE > 0.10; |dv_L| - 2 SE > 0.02; |r_B - 1| - 2 SE > 0.15.
PASSED if all: |D_K| + 2 SE < 0.10; |dv_L| + 2 SE < 0.02; cones - 2 SE > 0; |r_B - 1| + 2 SE < 0.15.
Otherwise UNDECIDED. Each condition on its own; nothing averaged; nothing chosen after the run.

    python run.py F5.2c --workers N          (seed 0 decides; --seed 1 replicates)
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import time

import numpy as np

from core.order.sprinkle import derive_seed
from core.stats.parallel import parallel_map
from core.testing import Result
from phases.F3_structure import field1d as F
from diag import b0_chain_field as B0
from diag.b2_curved import Continuum, sprinkle_curved, links_curved, links_brute, tau_cut, proper_volume

TEST = {
    "id": "F5.2c",
    "phase": "F5",
    "title": "Walls of a field in a prescribed curved 2+1 order of the theory (intrinsic Green weight, B2)",
    "theory_part": "F5 gravity; direction B (2+1): the order of the theory (narrower cones, more events) prescribed; "
                   "the Green function from the order alone (decisions 16, 17, 20); the passive response of walls",
    "question": "In an exact sprinkling of the 2+1 metric of the theory, does the field rule that reads only the order "
                "move walls as the continuum does: the kink (tension = energy) stays at rest, a launched kink follows the "
                "cones of the order, and a breather wall (little tension) falls as the continuum predicts?",
    "method": "2+1 sprinkling of Pi^2 = (x + ell)/ell, ell 24, in the domain of dependence; exact relation and links; "
              "phi = phi_0 - (1/rho) sum w(L) sin phi with the chain weight of F3.15; kink at rest, kink launched at 0.3, "
              "breather omega 0.7; breather on a flat order as a control; readouts paired with the continuum solution",
    "prediction": "continuum: kink exactly static (D_K = 0), launched kink with constant proper speed, breather falls "
                  "(about 0.78 of a test particle); on the flat order the breather stays",
    "criterion": "INSTRUMENT: links exact, I4 < 1e-12, density 0.01, kink Q 0.01, w 0.05, flat breather |D_0| 0.25 (+ 2 SE). "
                 "PASSED: |D_K| + 2 SE < 0.10, |dv_L| + 2 SE < 0.02, cones - 2 SE > 0, |r_B - 1| + 2 SE < 0.15. FAILED: "
                 "|D_K| - 2 SE > 0.10 or |dv_L| - 2 SE > 0.02 or |r_B - 1| - 2 SE > 0.15. Otherwise UNDECIDED (card "
                 "docs/B2_card_draft.md rev. 2)",
    "null_model": "the flat order (ell = inf) for the breather; the test-particle fall for the kink; exact continuum "
                  "solutions at the same events",
    "level": "mechanism",
    "meters": [],
    "depends_on": [],
    "params": {"rho": 30.0, "ell": 24.0, "n_dec": 32, "n_ctrl": 8, "xl": 13.0, "xl_flat": 16.0, "xr": 22.0, "Y": 12.5, "T": 8.0,
               "band": 2.0, "v": 0.3, "omega": 0.7,
               "calib": [0.2732723632276477, 0.21313211596642687, 1.2025922117299621],
               "i_rho": 30.0, "i_box": 3.0, "i_T": 2.5, "i4_T": 3.0, "i4_tol": 1e-12,
               "tol_N": 0.01, "tol_Q": 0.01, "tol_w": 0.05, "tol_K": 0.10, "tol_L": 0.02, "tol_B": 0.15,
               "tol_0": 0.25, "n_se": 2.0},
    "quick": {"rho": 12.0, "n_dec": 2, "n_ctrl": 1, "i_rho": 15.0},
    "estimate": "long: ~110 000 events per curved arm and ~127 000 per flat arm at rho 30, 15-20 min per arm on one core "
                "(the cloud container; ~0.8 of that on the author's machine); 3 x 32 + 8 = 104 arms, ~26 core-hours. numba "
                "required. Finished tasks are cached in diag/out/F5_2c_cache (keyed by this file's hash), so an "
                "interrupted run resumes.",
}

_nb = F._nb

if _nb is not None:
    @_nb.njit(cache=True)
    def _evolve(ptr, idx, phi0, rho, mu, wL):
        """The field rule with the chain weight (as F3.15, conservative arm). Returns phi."""
        N = phi0.size
        acc = np.zeros(N)
        phi = np.empty(N)
        L = np.zeros(N, np.int64)
        stamp = np.full(N, -1, np.int64)
        for v in range(N):
            ph = phi0[v] - acc[v]
            phi[v] = ph
            coef = mu * math.sin(ph) / rho
            L[v] = 0
            stamp[v] = v
            for u in range(v + 1, N):
                b = -1
                for q in range(ptr[u], ptr[u + 1]):
                    p = idx[q]
                    if p >= v and stamp[p] == v and L[p] + 1 > b:
                        b = L[p] + 1
                if b > 0:
                    L[u] = b
                    stamp[u] = v
                    acc[u] += wL[b] * coef
        return phi


def e_centroid(t, x, phi, tm, half=7.0, cell=0.5, step=0.25, x0=0.0):
    """Centroid of the noise-debiased coordinate energy density over |x - x0| < half at time tm (local regressions)."""
    s = np.abs(t - tm) < 0.5 * cell
    ts, xs, ps = t[s], x[s], phi[s]
    o = np.argsort(xs)
    ts, xs, ps = ts[o], xs[o], ps[o]
    num = den = 0.0
    for xg in np.arange(x0 - half, x0 + half + 1e-9, step):
        i0, i1 = np.searchsorted(xs, [xg - 0.5 * cell, xg + 0.5 * cell])
        n = i1 - i0
        if n < 8:
            return float("nan"), float("nan")
        A = np.column_stack([np.ones(n), ts[i0:i1] - tm, xs[i0:i1] - xg])
        c, *_ = np.linalg.lstsq(A, ps[i0:i1], rcond=None)
        r = ps[i0:i1] - A @ c
        C = float(r @ r) / (n - 3) * np.linalg.inv(A.T @ A)
        e = 0.5 * (c[1] ** 2 - C[1, 1]) + 0.5 * (c[2] ** 2 - C[2, 2]) + (1.0 - math.cos(c[0]) * math.exp(0.5 * C[0, 0]))
        num += xg * e
        den += e
    return num / den, den * step


def fit_breather(t, x, phi, a, b, guess=0.0, k0=0.714):
    """Centre of the breather wall in the slab a <= t < b: least squares of 4 atan(A(t) / cosh(k (x - xc))) with
    A(t) = A0 + A1 (t - tm) + A2 (t - tm)^2 over |x - guess| < 5. Returns xc, k, A0."""
    from scipy.optimize import least_squares
    s = (t >= a) & (t < b) & (np.abs(x - guess) < 5.0)
    tt, xx, pp = t[s] - 0.5 * (a + b), x[s], phi[s]
    A0 = float(np.tan(np.clip(np.mean(pp[np.abs(xx - guess) < 0.5]), -6, 6) / 4.0))

    def f(q):
        xc, lk, a0, a1, a2 = q
        return 4.0 * np.arctan((a0 + a1 * tt + a2 * tt * tt) / np.cosh(math.exp(lk) * (xx - xc))) - pp

    r = least_squares(f, [guess, math.log(k0), A0, 0.0, 0.0])
    return float(r.x[0]), float(math.exp(r.x[1])), float(r.x[2])


def _ell(e):
    return e if np.isfinite(e) else -1.0


def task_arm(p, ell, k, arm, base):
    T = float(p["T"])
    rho = float(p["rho"])
    sd = derive_seed(base, "F5.2c", rho, _ell(ell), k)
    xl = float(p["xl"]) if np.isfinite(ell) else float(p["xl_flat"])
    t, x, y = sprinkle_curved(rho, ell, xl, float(p["xr"]), float(p["Y"]), T, sd)
    om = float(p["omega"]) if arm.startswith("br") else None
    v = float(p["v"]) if arm == "kink_v" else 0.0
    t0 = time.time()
    Cx = Continuum(ell, 1.0, v, 0.0, T, -xl, float(p["xr"]), omega=om)
    Cf = Continuum(ell, 1.0, v, 0.0, T, -xl, float(p["xr"]), omega=om, mu_dyn=0.0)
    phi0, exact = Cf.at(t, x), Cx.at(t, x)
    ptr, idx = links_curved(t, x, y, _ell(ell), tau_cut(rho))
    phi = _evolve(ptr, idx, phi0, rho, 1.0, B0.weight_table("B", rho, tuple(p["calib"])))
    band = np.abs(y) < float(p["band"])
    tb, xb, pb, eb = t[band], x[band], phi[band], exact[band]
    rows = []
    for a in np.arange(0.0, T - 2.0 + 1e-9, 1.0):            # slab rule: the last slab is not read
        tm = a + 0.5
        row = {"t": tm}
        if om is None:
            xg = Cx.centre(tm)
            prev = rows[-1].get("xc") if rows else None
            guess = prev if (prev is not None and np.isfinite(prev)) else xg
            s = (tb >= a) & (tb < a + 1.0)
            r = F.fit_kink(tb[s], xb[s], pb[s], a, a + 1.0, guess, 1.0)
            rx = F.fit_kink(tb[s], xb[s], eb[s], a, a + 1.0, xg, 1.0)
            row.update({"xc": r["xc"] if r else float("nan"), "w": r["w"] if r else float("nan"),
                        "Q": r["Q"] if r else float("nan"), "xc_ex": rx["xc"] if rx else float("nan"),
                        "w_ex": rx["w"] if rx else float("nan"), "Q_ex": rx["Q"] if rx else float("nan")})
        else:
            c, e = e_centroid(tb, xb, pb, tm)
            cx, ex = e_centroid(tb, xb, eb, tm)
            fb, fbx = fit_breather(tb, xb, pb, a, a + 1.0), fit_breather(tb, xb, eb, a, a + 1.0)
            row.update({"c": c, "E": e, "c_ex": cx, "E_ex": ex, "xb": fb[0], "kb": fb[1], "xb_ex": fbx[0],
                        "kb_ex": fbx[1]})
        rows.append(row)
    return {"kind": "arm", "ell": ell if np.isfinite(ell) else "inf", "k": k, "arm": arm, "N": int(t.size),
            "V": proper_volume(ell, xl, float(p["xr"]), float(p["Y"]), T) if arm != "kink_v" else None,
            "sec": time.time() - t0, "rows": rows}


def task_check(p, base):
    """I1 links against the brute-force transitive reduction (curved and flat); I4 linear extension t + 0.2 x and the
    smaller rectangle (xr - 2) on a small curved domain (kink at rest)."""
    out = {"kind": "check"}
    rho = float(p["i_rho"])
    b = float(p["i_box"])
    ok = True
    for j, ell in enumerate((float(p["ell"]), float("inf"))):
        t, x, y = sprinkle_curved(rho, ell, b, b, b, float(p["i_T"]), derive_seed(base, "F5.2c-links", j))
        ptr, idx = links_curved(t, x, y, _ell(ell), tau_cut(rho))
        B, cnt = links_brute(t, x, y, _ell(ell))
        same = cnt == idx.size and all(B[kk, jj] for jj in range(t.size) for kk in idx[ptr[jj]:ptr[jj + 1]])
        ok = ok and bool(same)
        out[f"links_{j}"] = [int(t.size), int(idx.size), int(cnt)]
    out["links_ok"] = ok
    ell, T = float(p["ell"]), float(p["i4_T"])
    t, x, y = sprinkle_curved(rho, ell, 6.0, 8.0, 6.0, T, derive_seed(base, "F5.2c-i4"))
    C0 = Continuum(ell, 1.0, 0.0, 0.0, T, -6.0, 8.0, mu_dyn=0.0)
    phi0 = C0.at(t, x)
    wL = B0.weight_table("B", rho, tuple(p["calib"]))
    ptr, idx = links_curved(t, x, y, ell, tau_cut(rho))
    ref = _evolve(ptr, idx, phi0, rho, 1.0, wL)
    sc = float(np.max(np.abs(ref)))
    o = np.argsort(t + 0.2 * x, kind="stable")
    inv = np.empty_like(o)
    inv[o] = np.arange(o.size)
    rows = [np.sort(inv[idx[ptr[j]:ptr[j + 1]]]) for j in o]
    ptr2 = np.zeros(t.size + 1, np.int64)
    ptr2[1:] = np.cumsum([r.size for r in rows])
    idx2 = np.concatenate(rows).astype(np.int64)
    alt = _evolve(ptr2, idx2, phi0[o], rho, 1.0, wL)
    out["order"] = float(np.max(np.abs(alt - ref[o])) / sc)
    from diag.b2_curved import dist_boundary
    sub = np.nonzero(t <= dist_boundary(x, y, ell, 6.0, 6.0, 6.0))[0]
    ps, ix = links_curved(t[sub], x[sub], y[sub], ell, tau_cut(rho))
    phs = _evolve(ps, ix, phi0[sub], rho, 1.0, wL)
    out["region"] = float(np.max(np.abs(phs - ref[sub])) / sc)
    out["N_i4"] = int(t.size)
    return out


_SHA = None


def _sha():
    global _SHA
    if _SHA is None:
        with open(__file__, "rb") as fh:
            _SHA = hashlib.sha256(fh.read()).hexdigest()[:16]
    return _SHA


def run_task(kind, *a):
    """Runs a task; finished tasks are cached on disk (keyed by this file's hash), so a run can resume."""
    p = a[0]
    key = json.dumps([kind] + [str(q) for q in a[1:]])
    fn = os.path.join("diag", "out", "F5_2c_cache", _sha(), hashlib.sha256(
        (key + json.dumps(p, sort_keys=True, default=str)).encode()).hexdigest()[:24] + ".json")
    if os.path.exists(fn):
        with open(fn) as fh:
            return json.load(fh)
    out = {"arm": task_arm, "check": task_check}[kind](*a)
    os.makedirs(os.path.dirname(fn), exist_ok=True)
    with open(fn + ".tmp", "w") as fh:
        json.dump(out, fh)
    os.replace(fn + ".tmp", fn)
    return out


# ============================================================================================== analysis

def _col(rows, key):
    return np.array([r.get(key, float("nan")) for r in rows], float)


def _ms(a):
    a = np.asarray(a, float)
    if a.size < 2 or not np.all(np.isfinite(a)):
        return float("nan"), float("nan")
    return float(a.mean()), float(a.std(ddof=1) / math.sqrt(a.size))


def _fmt(m):
    return f"{m[0]:+.4f} +- {m[1]:.4f}"


def _slope(rows, key):
    tm, xc = _col(rows, "t"), _col(rows, key)
    ok = np.isfinite(xc)
    if ok.sum() < tm.size - 1:
        return float("nan")
    return float(np.polyfit(tm[ok], xc[ok], 1)[0])


def _rB(r, key="xb"):
    c, cx = _col(r["rows"], key), _col(r["rows"], key + "_ex")
    d, dx = c[1:] - c[0], cx[1:] - cx[0]
    return float(np.sum(d * dx) / np.sum(dx * dx))


def run(p, ctx):
    if _nb is None:
        raise SystemExit("F5.2c needs numba (pip install numba)")
    t0 = time.time()
    ell = float(p["ell"])
    n_dec, n_ctrl, nse = int(p["n_dec"]), int(p["n_ctrl"]), float(p["n_se"])
    base = ctx.seed
    jobs = [("check", p, base)]
    for k in range(n_dec):
        for arm in ("kink_rest", "kink_v", "br"):
            jobs.append(("arm", p, ell, k, arm, base))
    for k in range(n_ctrl):
        jobs.append(("arm", p, float("inf"), k, "br", base))
    res = parallel_map(run_task, jobs, ctx.workers, ctx.progress,
                       lambda i, r, d, nn: f"[F5.2c] {d}/{nn}: {r['kind']} {r.get('arm', '')} ell {r.get('ell', '')} "
                                           f"k {r.get('k', '')} ({time.time() - t0:.0f}s)")

    def A(arm, ell_, m):
        out = {r["k"]: r for r in res if r["kind"] == "arm" and r["arm"] == arm and str(r["ell"]) == str(ell_)}
        return [out[k] for k in range(m)]

    lines = [f"F5.2c (B2, 2+1, prescribed order of the theory Pi^2 = (x + {ell:g})/{ell:g}; intrinsic weight variant B); "
             f"rho {p['rho']:g}, T {p['T']:g}, {n_dec} sprinklings; flat control on {n_ctrl}; base seed {base}; "
             f"mean +- SE over sprinklings"]
    inst, passc, failc = {}, {}, {}
    ck = [r for r in res if r["kind"] == "check"][0]
    inst["I1 links exact"] = bool(ck["links_ok"])
    inst["I4 order, region"] = max(ck["order"], ck["region"]) < float(p["i4_tol"])
    lines.append(f"  I1 links against the brute-force transitive reduction: {ck['links_0']} (curved), {ck['links_1']} "
                 f"(flat) [N, links, brute]: {'identical' if ck['links_ok'] else 'DIFFERENT'}")
    lines.append(f"  I4 linear extension t + 0.2 x {ck['order']:.1e}, region {ck['region']:.1e} (N {ck['N_i4']})")
    KR, KV, BR = A("kink_rest", ell, n_dec), A("kink_v", ell, n_dec), A("br", ell, n_dec)
    B0f = A("br", "inf", n_ctrl)
    dens = _ms([r["N"] / (float(p["rho"]) * r["V"]) - 1.0 for r in KR])
    inst["I0 density"] = abs(dens[0]) + nse * dens[1] < float(p["tol_N"])
    lines.append(f"  I0 density: N / (rho V_proper) - 1 {_fmt(dens)}")
    q = _ms([np.nanmean(_col(r["rows"], "Q")) - np.nanmean(_col(r["rows"], "Q_ex")) for r in KR])
    w = _ms([np.nanmean(_col(r["rows"], "w")) / np.nanmean(_col(r["rows"], "w_ex")) - 1 for r in KR])
    inst["kink Q"] = abs(q[0]) + nse * q[1] < float(p["tol_Q"])
    inst["kink w"] = abs(w[0]) + nse * w[1] < float(p["tol_w"])
    lines.append(f"  kink at rest: Q - Q_ex {_fmt(q)}; w/w_ex - 1 {_fmt(w)}")
    # flat breather: no fall without curvature (relative to the curved continuum fall)
    fall_ex = float(np.mean([_col(r["rows"], "xb_ex")[-1] - _col(r["rows"], "xb_ex")[0] for r in BR]))
    d0 = _ms([(_col(r["rows"], "xb")[-1] - _col(r["rows"], "xb")[0]) / fall_ex for r in B0f])
    inst["flat breather stays"] = abs(d0[0]) + nse * d0[1] < float(p["tol_0"])
    lines.append(f"  flat order, breather: D_0 = (xb(last) - xb(first)) / curved continuum fall ({fall_ex:+.4f}) "
                 f"{_fmt(d0)}")
    # decision
    part = -(_col(KR[0]["rows"], "t")[-1] ** 2 - _col(KR[0]["rows"], "t")[0] ** 2) / (4 * ell)
    DK = _ms([(_col(r["rows"], "xc")[-1] - _col(r["rows"], "xc")[0])
              - (_col(r["rows"], "xc_ex")[-1] - _col(r["rows"], "xc_ex")[0]) for r in KR])
    dv = _ms([_slope(r["rows"], "xc") - _slope(r["rows"], "xc_ex") for r in KV])
    rB = _ms([_rB(r) for r in BR])
    for lab, m, tol, ref in (("kink stays: D_K", DK, "tol_K", 0.0), ("launched kink: dv_L", dv, "tol_L", 0.0),
                             ("breather falls: r_B", rB, "tol_B", 1.0)):
        passc[lab] = abs(m[0] - ref) + nse * m[1] < float(p[tol])
        failc[lab] = abs(m[0] - ref) - nse * m[1] > float(p[tol])
    cones = _ms([_slope(r["rows"], "xc") - float(p["v"]) for r in KV])
    passc["launched kink: faster than flat cones"] = cones[0] - nse * cones[1] > 0.0
    sl_ex = float(np.mean([_slope(r["rows"], "xc_ex") for r in KV]))
    lines.append(f"  kink at rest: D_K {_fmt(DK)} (a test particle: {part:+.4f})")
    lines.append(f"  launched kink v {p['v']:g}: slope - continuum slope {_fmt(dv)} (continuum slope {sl_ex:+.4f}; "
                 f"flat cones would give {p['v']:g}); slope - v {_fmt(cones)}")
    cfall = float(np.mean([_col(r["rows"], "c_ex")[-1] - _col(r["rows"], "c_ex")[0] for r in BR]))
    lines.append(f"  breather omega {p['omega']:g}: r_B {_fmt(rB)} (continuum fall of the profile centre {fall_ex:+.4f}); "
                 f"reported: the energy centroid r {_fmt(_ms([_rB(r, 'c') for r in BR]))} (continuum {cfall:+.4f})")
    inst = {k: bool(v) for k, v in inst.items()}
    passc = {k: bool(v) for k, v in passc.items()}
    failc = {k: bool(v) for k, v in failc.items()}
    lines.append("instrument: " + "; ".join(f"{k}: {'OK' if v else 'FAILS'}" for k, v in inst.items()))
    lines.append("passed conditions: " + "; ".join(f"{k}: {'yes' if v else 'no'}" for k, v in passc.items()))
    lines.append("failed conditions: " + "; ".join(f"{k}: {'YES' if v else 'no'}" for k, v in failc.items()))
    if not all(inst.values()):
        status = "INSTRUMENT_FAILED"
    elif any(failc.values()):
        status = "FAILED"
    elif all(passc.values()):
        status = "PASSED"
    else:
        status = "UNDECIDED"
    lines.append(f"cost: {np.mean([r['sec'] for r in res if r['kind'] == 'arm']):.0f} s per arm ({time.time() - t0:.0f}s)")
    return Result(status, lines, {"instrument": inst, "passed": passc, "failed": failc, "results": res})
