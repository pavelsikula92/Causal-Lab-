"""
B0, development stage 6: the field rule in 2+1 with the intrinsic Green weight of the longest chain (decides nothing)
====================================================================================================================
Decisions 16 (D1) and 17. The field rule
    phi(v) = phi_0(v) - (1 / rho) * sum over u < v of w(L_uv) J(u),   J = V'(phi) = mu sin(phi)
is evaluated with the longest chain L_uv from the order alone: the links of the sprinkling (shell method,
diag/b0_chain_alg.links3) and, for every source u in time order, a forward dynamic programme over the links of its
future (L_u(u) = 0, L_u(v) = 1 + max over links p -> v with p in the future of u of L_u(p)); the contribution
w(L) J(u) / rho is added to every v in the future of u. Events are processed in time order, so phi(v) is final before v
becomes a source. Cost O(sum over sources of the links of their future) ~ S * E_H / 2.

General source (the author's requirement): the kernel takes either the self-consistent source mu sin(phi) (the field
evolution) or a prescribed array J (linear application of the Green function, e.g. J_g = T or 2V for the clock response).
Threshold: sources with |J| <= eps_J are skipped. eps_J is a declared numerical approximation; the default is 0 (none).

Variants on the same sprinkling and the same initial data (the free wave phi_0):
  coord  coordinate tau with the cap 0.68 rho^(-1/3) (diag/b0_green2p1._evolve3; reference, development only);
  A      intrinsic, single alpha_3 = 1.2: w(L) = alpha_3 rho^(1/3) / (2 pi L);
  B      intrinsic, the table of decision 17: w(1) = a rho^(1/3), w(L >= 2) = c rho^(1/3) / (L + beta), from the newest
         diag/out/b0_weight_calib_*.json.
Arms (readouts as diag/b0_green2p1.py and diag/b0_vacuum_long.py):
  rest, v0.3   the wall (1+1 kink along y): velocity, rms against the exact kink per slab;
  k1_0, k0_1, k1_1, k2_1   plane waves of amplitude 0.01: omega against sqrt(k^2 + 1), and the amplitude A(t) / A(0)
               (the energy of a single linear wave is proportional to A^2: the energy proxy);
  vacuum       8 random waves, total amplitude 0.01, pyramid |x|, |y| <= T + 4 - t: A / A_exact per slab, lambda_eff.

    python -m diag.b0_chain_field --check
    python -m diag.b0_chain_field --rhos 50 --seeds 4 --T 6 --arms rest,v0.3,k1_0,k1_1 --variants coord,A,B --workers 6
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
from diag.b0_green2p1 import _evolve3, sprinkle_pyramid, modes_free_2d
from diag.b0_chain_alg import links3

_nb = F._nb
ALPHA3 = 1.2
LMAX = 4096

if _nb is not None:
    @_nb.njit(cache=True)
    def evolve_chain(ptr, idx, phi0, rho, mu, wL, eps_J, Jext, selfcons):
        """selfcons: J = mu sin(phi) (field evolution); else J = Jext (linear Green application; phi = phi0 - acc).
        Returns phi, the number of sources used, and the largest L met."""
        N = phi0.size
        acc = np.zeros(N)
        phi = np.empty(N)
        L = np.zeros(N, np.int64)
        stamp = np.full(N, -1, np.int64)
        nsrc = 0
        lmax = 0
        for v in range(N):
            ph = phi0[v] - acc[v]
            phi[v] = ph
            J = mu * math.sin(ph) if selfcons else Jext[v]
            if abs(J) <= eps_J:
                continue
            nsrc += 1
            coef = J / rho
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
                    if b > lmax:
                        lmax = b
                    acc[u] += wL[b] * coef
        return phi, nsrc, lmax


def newest_calib():
    fs = sorted(glob.glob(os.path.join("diag", "out", "b0_weight_calib_*.json")))
    if not fs:
        raise SystemExit("run diag.b0_weight_calib first")
    d = json.load(open(fs[-1]))
    return d["a"], d["c"], d["beta"], fs[-1]


def weight_table(variant, rho, calib):
    r13 = rho ** (1.0 / 3.0)
    Ls = np.arange(LMAX, dtype=float)
    w = np.zeros(LMAX)
    if variant == "A":
        w[1:] = ALPHA3 * r13 / (2 * math.pi * Ls[1:])
    elif variant == "B":
        a, c, beta = calib
        w[1] = a * r13
        w[2:] = c * r13 / (Ls[2:] + beta)
    else:
        raise ValueError(variant)
    return w


def box(arm, T):
    return (T + 4.0, T + 4.0) if arm == "vacuum" else (T + 8.0, T + 3.0)


def initial(arm, t, x, y, seed):
    if arm in ("rest", "v0.3"):
        v = 0.0 if arm == "rest" else 0.3
        return F.kink_free(t, x, v), {"v": v}
    if arm.startswith("k"):
        kx, ky = (float(u) for u in arm[1:].split("_"))
        return modes_free_2d(t, x, y, [(0.01, kx, ky, 0.0)]), {"kx": kx, "ky": ky}
    if arm == "vacuum":
        rng = np.random.default_rng(derive_seed(seed, "modes"))
        modes = []
        for _ in range(8):
            k = rng.uniform(0.3, 3.0)
            th = rng.uniform(0, 2 * math.pi)
            modes.append((0.01 / 2.0, k * math.cos(th), k * math.sin(th), rng.uniform(0, 2 * math.pi)))
        return modes_free_2d(t, x, y, modes), {"modes": modes}
    raise ValueError(arm)


def readout(arm, t, x, y, phi, T, info):
    out = {}
    if arm in ("rest", "v0.3"):
        v = info["v"]
        ex = F.kink_exact(t, x, v)
        c = np.abs(y) < 2.0
        xc, rms = [], []
        for a in np.arange(0.0, T - 1.0 + 1e-9, 1.0):
            g = xc[-1] if (xc and np.isfinite(xc[-1])) else v * (a + 0.5)
            s = c & (t >= a) & (t < a + 1.0)
            r = F.fit_kink(t[s], x[s], phi[s], a, a + 1.0, g, 1.0)
            xc.append(r["xc"] if r else float("nan"))
            sn = s & (np.abs(x - v * (a + 0.5)) < 4.0)
            rms.append(float(np.sqrt(np.mean((phi[sn] - ex[sn]) ** 2))))
        out.update({"xc": xc, "rms": rms, "v": v})
    elif arm.startswith("k"):
        kx, ky = info["kx"], info["ky"]
        ph, amp, tm = [], [], []
        c = (np.abs(x) < 6.0) & (np.abs(y) < 4.0)
        for a in np.arange(0.0, T - 0.25 + 1e-9, 0.25):
            s = c & (t >= a) & (t < a + 0.25)
            arg = kx * x[s] + ky * y[s]
            M = np.column_stack([np.cos(arg), np.sin(arg), np.ones(s.sum())])
            cf, *_ = np.linalg.lstsq(M, phi[s], rcond=None)
            ph.append(math.atan2(cf[1], cf[0]))
            amp.append(math.hypot(cf[0], cf[1]))
            tm.append(a + 0.125)
        om = float(np.polyfit(np.array(tm), np.unwrap(np.array(ph)), 1)[0])
        out.update({"omega": om, "omega_c": math.sqrt(kx * kx + ky * ky + 1.0), "amp": amp, "tm": tm})
    elif arm == "vacuum":
        ex = np.zeros_like(x)
        for a_, kx, ky, th in info["modes"]:
            ex += a_ * np.cos(kx * x + ky * y - math.sqrt(kx * kx + ky * ky + 1.0) * t + th)
        c = (np.abs(x) < 4.0) & (np.abs(y) < 4.0)
        A, Ae, E = [], [], []
        for s0 in np.arange(0.0, T - 1.0 + 1e-9, 1.0):
            m = c & (t >= s0) & (t < s0 + 1.0)
            A.append(float(np.sqrt(np.mean(phi[m] ** 2))))
            Ae.append(float(np.sqrt(np.mean(ex[m] ** 2))))
            E.append(float(np.sqrt(np.mean((phi[m] - ex[m]) ** 2))))
        out.update({"A": A, "Ae": Ae, "E": E})
    return out


def task(rho, seed, arm, variant, T, calib, eps_J):
    X, Y = box(arm, T)
    t, x, y = sprinkle_pyramid(rho, X, Y, T, seed)
    phi0, info = initial(arm, t, x, y, seed)
    t0 = time.time()
    extra = {}
    if variant == "coord":
        phi = _evolve3(t, x, y, phi0, float(rho), 1.0, 0.68 * rho ** (-1.0 / 3.0))
    else:
        tau_cut = (12.0 * 25.0 / (math.pi * rho)) ** (1.0 / 3.0)
        ptr, idx = links3(t, x, y, tau_cut)
        extra["sec_links"] = time.time() - t0
        extra["links"] = int(idx.size)
        phi, nsrc, lmax = evolve_chain(ptr, idx, phi0, float(rho), 1.0, weight_table(variant, rho, calib), eps_J,
                                       np.zeros(1), True)
        if lmax >= LMAX - 1:
            raise RuntimeError("L beyond the weight table")
        extra.update({"nsrc": int(nsrc), "lmax": int(lmax)})
    out = {"rho": rho, "seed": seed, "arm": arm, "variant": variant, "N": int(t.size), "sec": time.time() - t0}
    out.update(extra)
    out.update(readout(arm, t, x, y, phi, T, info))
    return out


def check():
    """Instrument checks: (1) the linear Green application equals the backward-DP sum of diag/b0_weight_valid on the same
    sprinkling and source; (2) a linear extension other than time gives the same field (order independence is not
    available here: the kernel needs a time order of indices; checked instead: a relabelling of events with equal
    order gives the same field)."""
    from diag.b0_weight_valid import children, target_pass
    calib = newest_calib()[:3]
    rho, T = 50.0, 3.0
    t, x, y = sprinkle_pyramid(rho, T + 2, T + 2, T, 3)
    tau_cut = (12.0 * 25.0 / (math.pi * rho)) ** (1.0 / 3.0)
    ptr, idx = links3(t, x, y, tau_cut)
    J = 1.0 / np.cosh((0.8 * x + 0.6 * y) / 0.7) ** 2
    phi, nsrc, lmax = evolve_chain(ptr, idx, np.zeros(t.size), rho, 1.0, weight_table("B", rho, calib), 0.0, J, False)
    cptr, cidx = children(ptr, idx, t.size)
    rng = np.random.default_rng(1)
    tg = rng.choice(np.nonzero(t > T - 1)[0], 30, replace=False)
    sw = np.zeros((3, 120)); sg = np.zeros((3, 120)); cn = np.zeros((3, 120))
    dev = 0.0
    for xi in tg:
        o = target_pass(int(xi), t, x, y, cptr, cidx, J[None, :], rho, *calib, 0.68 * rho ** (-1 / 3), sw, sg, cn)
        dev = max(dev, abs(-phi[xi] - o[0, 2]) / max(abs(o[0, 2]), 1e-30))
    print(f"(1) forward kernel against the backward DP of the validation, 30 targets: max relative difference {dev:.2e}")
    # (2) the same order under another time-compatible labelling: t' = t + 0.3 x (spacelike-tilted slicing is a valid
    # linear extension only if it respects the order; t + 0.3 x does, since |dx| < dt for every relation)
    key = t + 0.3 * x
    o = np.argsort(key, kind="stable")
    inv = np.empty_like(o); inv[o] = np.arange(o.size)
    t2, x2, y2 = t[o], x[o], y[o]
    # links in the new labelling: map the parents and resort
    ptr2 = np.zeros(t.size + 1, np.int64)
    rows = [np.sort(inv[idx[ptr[j]:ptr[j + 1]]]) for j in o]
    ptr2[1:] = np.cumsum([r.size for r in rows])
    idx2 = np.concatenate(rows).astype(np.int64)
    phi0 = F.kink_free(t, x, 0.3)
    p1, *_ = evolve_chain(ptr, idx, phi0, rho, 1.0, weight_table("B", rho, calib), 0.0, np.zeros(1), True)
    p2, *_ = evolve_chain(ptr2, idx2, phi0[o], rho, 1.0, weight_table("B", rho, calib), 0.0, np.zeros(1), True)
    print(f"(2) field under the linear extension t + 0.3 x against time order: max |difference| "
          f"{np.max(np.abs(p2 - p1[o])):.2e}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=4)
    ap.add_argument("--rhos", type=str, default="50")
    ap.add_argument("--T", type=float, default=6.0)
    ap.add_argument("--arms", type=str, default="rest,v0.3,k1_0,k0_1,k1_1,k2_1,vacuum")
    ap.add_argument("--variants", type=str, default="coord,A,B")
    ap.add_argument("--eps", type=float, default=0.0)
    ap.add_argument("--workers", type=int, default=0)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    if _nb is None:
        raise SystemExit("numba required")
    if a.check:
        check()
        return
    ca, cc, cb, fn_cal = newest_calib()
    calib = (ca, cc, cb)
    rhos = [float(r) for r in a.rhos.split(",")]
    arms = a.arms.split(",")
    variants = a.variants.split(",")
    jobs = [(rho, derive_seed(20261111, "B0", rho, k), arm, var, a.T, calib, a.eps)
            for rho in rhos for k in range(a.seeds) for arm in arms for var in variants]
    jobs.sort(key=lambda j: (-j[0], j[3] == "coord"))
    t0 = time.time()
    res = parallel_map(task, jobs, default_workers(a.workers), print,
                       lambda i, r, d, n: f"[B0 chain] {d}/{n}: {r['arm']} {r['variant']} rho {r['rho']:g} N {r['N']} "
                                          f"{r['sec']:.0f}s ({time.time() - t0:.0f}s)")
    L_ = [f"2+1 field rule, variants {variants} on the same sprinklings and initial data; calibration {fn_cal} "
          f"(a {ca:.4f}, c {cc:.4f}, beta {cb:.3f}); alpha_3 {ALPHA3}; eps_J {a.eps:g}; T {a.T:g}"]
    for rho in rhos:
        L_.append(f"=== rho {rho:g}")
        for arm in arms:
            for var in variants:
                R = [r for r in res if r["rho"] == rho and r["arm"] == arm and r["variant"] == var]
                if not R:
                    continue
                tag = f"  {arm:7s} {var:5s}:"
                if arm in ("rest", "v0.3"):
                    X = np.array([r["xc"] for r in R])
                    tm = np.arange(0.5, a.T, 1.0)
                    sl = [np.polyfit(tm, xc, 1)[0] for xc in X if np.all(np.isfinite(xc))]
                    rms = np.mean([r["rms"] for r in R], axis=0)
                    se = np.std(sl, ddof=1) / math.sqrt(len(sl)) if len(sl) > 1 else float("nan")
                    L_.append(f"{tag} velocity {np.mean(sl):+.4f} +- {se:.4f} (launch {R[0]['v']}); rms against the kink "
                              + " ".join(f"{q:.3f}" for q in rms))
                elif arm.startswith("k"):
                    om = [r["omega"] for r in R]
                    oc = R[0]["omega_c"]
                    amp = np.mean([np.array(r["amp"]) for r in R], axis=0)
                    se = np.std(om, ddof=1) / math.sqrt(len(om)) if len(om) > 1 else float("nan")
                    n4 = len(amp) // 4
                    L_.append(f"{tag} omega/omega_c {np.mean(om) / oc:.4f} +- {se / oc:.4f}; amplitude first quarter "
                              f"{np.mean(amp[:n4]):.5f}, last quarter {np.mean(amp[-n4:]):.5f} (ratio "
                              f"{np.mean(amp[-n4:]) / np.mean(amp[:n4]):.3f})")
                elif arm == "vacuum":
                    tm = np.arange(0.5, a.T, 1.0)
                    ratio = np.mean([np.array(r["A"]) / np.array(r["Ae"]) for r in R], axis=0)
                    lam = [np.polyfit(tm, np.log(np.array(r["A"]) / np.array(r["Ae"])), 1)[0] for r in R]
                    se = np.std(lam, ddof=1) / math.sqrt(len(lam)) if len(lam) > 1 else float("nan")
                    L_.append(f"{tag} A/A_exact " + " ".join(f"{v:.3f}" for v in ratio)
                              + f"; lambda_eff {np.mean(lam):+.4f} +- {se:.4f}")
                if var != "coord":
                    L_[-1] += (f"  [links/event {np.mean([r['links'] / r['N'] for r in R]):.1f}, max L "
                               f"{max(r['lmax'] for r in R)}, s/run {np.mean([r['sec'] for r in R]):.0f}]")
                else:
                    L_[-1] += f"  [s/run {np.mean([r['sec'] for r in R]):.0f}]"
    print("\n".join(L_))
    os.makedirs(os.path.join("diag", "out"), exist_ok=True)
    fn = os.path.join("diag", "out", f"b0_chain_field_{time.strftime('%Y%m%dT%H%M%S', time.gmtime())}.json")
    with open(fn, "w") as fh:
        json.dump({"args": vars(a), "calib": fn_cal, "summary": L_, "results": res}, fh)
    print(f"raw records: {fn} ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
