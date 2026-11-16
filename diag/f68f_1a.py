"""
F6.8f-dev, step 1a (development, decides nothing): the exact interval volume of short intervals at finite b
===========================================================================================================

docs/derivations/f68f_dev.md, step 1a (the author's order of 16 Nov, 10:34). Optical sector g = -dt^2 + n^2 dx^2,
n = exp(b cos(k.w + theta)), 3+1. For the pair x = 0, y = (dt, L e) with cos(e, k) = CO:

    V(b, theta) = int d^3w n(w)^3 [dt - T(x, w) - T(w, y)]_+ ,

with T the optical distance (first arrival) computed exactly by Fermat, not by an expansion in b:
the ray between two points lies in the plane of the segment and k; its deviation from the chord is
delta(u) = u (1 - u) sum_j c_j T_j(2u - 1) (Chebyshev, spectral for the analytic ray), and the optical length
int n dl is minimised by Newton (analytic gradient and Hessian, Gauss nodes). Checked against the exact
Snell first arrival of step 0 (diag/f68f_dev.py) on random segments.

The volume: prolate coordinates around the chord, w = (c (1 + sigma te), rp cos phi, rp sin phi), c = L/2,
d^3w = c^3 (sigma^2 - te^2) dsigma dte dphi; on every ray (te, phi) the root sigma_r of dt - T(x,w) - T(w,y)
is found numerically (secant), and sigma in [1, sigma_r] is integrated by Gauss. te nodes of F6.8e3
(clustered at the tips), phi trapezoid, theta average on n_th equally spaced phases.

Two parametrisations of the pair:
    fixed sigma (as F6.8c / F6.8e3): dt = sqrt(T(x, y)^2 + tau^2), T(x, y) exact;
    fixed tau:                       dt = sqrt(L^2 + tau^2).
Reported per row and b:
    R(b) = 2 (<V(b)>_theta - V0) / (b^2 V0)   -> gamma as b -> 0,
    B(b) = 2 <V(b) e^{-i theta}>_theta / (b V0) -> beta,
    delta(tau^2)/tau^2 = max_theta |L^2 - T(x,y)^2| / tau^2 (the difference of the two parametrisations).
Compared with gamma, beta of F6.8e3 (second_order_volume_stable, fixed sigma) and with the same second-order
instrument at fixed tau (dt' = dt'' = 0, the small-b prediction for fixed tau).

Run: python diag/f68f_1a.py [workers]   (output diag/out/f68f_1a.txt)
"""
import inspect
import math
import os
import sys
import time
from multiprocessing import Pool

import numpy as np
from numpy.polynomial import chebyshev as C

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from phases.F6_order_dynamics import F6_08e3_bd_action_optical_asymptotic_rev2 as E3  # noqa: E402

M_RAY = 14
Q_RAY = 56


def _basis():
    x, w = np.polynomial.legendre.leggauss(Q_RAY)
    u = 0.5 * (x + 1.0)
    w = 0.5 * w
    Tm = np.stack([C.chebval(x, np.eye(M_RAY)[j]) for j in range(M_RAY)], 1)            # T_j(x)
    dTm = np.stack([C.chebval(x, C.chebder(np.eye(M_RAY)[j])) for j in range(M_RAY)], 1)  # T_j'(x)
    phi = (u * (1 - u))[:, None] * Tm
    dphi = (1 - 2 * u)[:, None] * Tm + (u * (1 - u))[:, None] * 2.0 * dTm
    return u, w, phi, dphi


_U, _W, _PHI, _DPHI = _basis()


def optical_length(zA, zB, d, b, k, th, iters=30):
    """Exact optical length of the rays between points with coordinates zA, zB along k (unit vector) and distance d
    (arrays of equal shape); n = exp(b cos(k z + th)). Fermat minimum by Newton in the Chebyshev basis."""
    shp = np.shape(d)
    zA, zB, d = (np.ravel(np.asarray(v, float)) for v in (zA, zB, d))
    d = np.maximum(d, 1e-300)
    pz = np.sqrt(np.maximum(1.0 - ((zB - zA) / d) ** 2, 0.0))
    z0 = zA[:, None] + _U[None, :] * (zB - zA)[:, None]
    c = np.zeros((d.size, M_RAY))
    for it in range(iters):
        dl = c @ _PHI.T
        dv = c @ _DPHI.T
        z = z0 + dl * pz[:, None]
        a = k * z + th
        n = np.exp(b * np.cos(a))
        n1 = -b * k * np.sin(a) * n
        n2 = (b * b * k * k * np.sin(a) ** 2 - b * k * k * np.cos(a)) * n
        Q = np.sqrt(d[:, None] ** 2 + dv * dv)
        fz = n1 * Q
        fv = n * dv / Q
        fzz = n2 * Q
        fzv = n1 * dv / Q
        fvv = n * d[:, None] ** 2 / Q ** 3
        P = pz[:, None]
        g = ((_W * fz * P)[:, :, None] * _PHI[None] + (_W * fv)[:, :, None] * _DPHI[None]).sum(1)
        A = _W * fzz * P * P
        Bm = _W * fzv * P
        Cc = _W * fvv
        Y = A[:, :, None] * _PHI[None] + Bm[:, :, None] * _DPHI[None]        # n x q x j
        Z = Bm[:, :, None] * _PHI[None] + Cc[:, :, None] * _DPHI[None]
        H = np.matmul(Y.transpose(0, 2, 1), _PHI) + np.matmul(Z.transpose(0, 2, 1), _DPHI)
        step = np.linalg.solve(H, g[:, :, None])[:, :, 0]
        c = c - step
        if np.max(np.abs(step)) < 1e-15 * max(1.0, float(np.max(d))):
            break
    dl = c @ _PHI.T
    dv = c @ _DPHI.T
    z = z0 + dl * pz[:, None]
    n = np.exp(b * np.cos(k * z + th))
    S = (_W * n * np.sqrt(d[:, None] ** 2 + dv * dv)).sum(1)
    return S.reshape(shp)


def volume(k, CO, L, dt, b, th, n_sig=16, n_te=12, n_phi=16, npan=None):
    """V(b, theta) of the pair x = 0, y = (dt, L e1), k = k (CO, sqrt(1 - CO^2), 0)."""
    kx, ky = CO, math.sqrt(max(1.0 - CO * CO, 0.0))
    c = L / 2.0
    npan = npan or int(8 + math.ceil(abs(k * CO) * L))
    te, wt = E3._te_nodes(n_te, npan)
    te, wt = np.asarray(te, float), np.asarray(wt, float)
    ph = 2 * math.pi * (np.arange(n_phi) + 0.5) / n_phi
    TE, PH = np.meshgrid(te, ph, indexing="ij")
    WT = (wt[:, None] * np.ones_like(PH)) * (2 * math.pi / n_phi)
    TE, PH, WT = TE.ravel(), PH.ravel(), WT.ravel()
    zy = kx * L

    def geom(sig):
        lw = c * (1 + sig * TE)
        rp = c * np.sqrt(np.maximum((sig * sig - 1) * (1 - TE * TE), 0.0))
        zw = kx * lw + ky * rp * np.cos(PH)
        d2 = c * (sig + TE)            # |w - x|
        d1 = c * (sig - TE)            # |y - w|
        return zw, d1, d2

    def ell(sig):
        zw, d1, d2 = geom(sig)
        return dt - optical_length(np.zeros_like(zw), zw, d2, b, k, th) - optical_length(zw, np.full_like(zw, zy), d1, b, k, th)

    # safeguarded secant on a bracket [lo, hi]: l(1) > 0 on the chord, l(hi) < 0 (hi expanded from the flat root)
    lo = np.ones(TE.shape)
    flo = ell(lo + 1e-12)
    if np.any(flo <= 0):
        return float("nan")                      # the chord leaves the diamond: prolate coordinates do not apply
    hi = np.full(TE.shape, dt / L)
    fhi = ell(hi)
    for _ in range(60):
        bad = fhi >= 0
        if not np.any(bad):
            break
        hi = np.where(bad, 1.0 + 2.0 * (hi - 1.0), hi)
        fhi = np.where(bad, ell(hi), fhi)
    s0, f0 = hi, fhi
    s1 = hi - 0.01 * (hi - 1.0)
    f1 = ell(s1)
    for it in range(80):
        lo = np.where(f1 > 0, np.maximum(lo, s1), lo)
        hi = np.where(f1 <= 0, np.minimum(hi, s1), hi)
        den = f1 - f0
        ok = den != 0
        s2 = np.where(ok, s1 - f1 * (s1 - s0) / np.where(ok, den, 1.0), s1)
        out = ~((s2 > lo) & (s2 < hi)) & ok
        s2 = np.where(out, 0.5 * (lo + hi), s2)
        s0, f0 = s1, f1
        s1 = s2
        f1 = ell(s1)
        if np.max(np.abs(s1 - s0)) < 1e-15 * float(np.max(s1)):
            break
    sr = s1
    xg, wg = np.polynomial.legendre.leggauss(n_sig)
    V = 0.0
    for x, wx in zip(xg, wg):
        sig = 1.0 + (sr - 1.0) * 0.5 * (x + 1.0)
        zw, d1, d2 = geom(sig)
        lval = dt - optical_length(np.zeros_like(zw), zw, d2, b, k, th) - optical_length(zw, np.full_like(zw, zy), d1, b, k, th)
        P = np.exp(3 * b * np.cos(k * zw + th))
        V += float((WT * 0.5 * wx * (sr - 1.0) * c ** 3 * (sig * sig - TE * TE) * P * lval).sum())
    return V


def pair_T(k, CO, L, b, th):
    return float(optical_length(np.array([0.0]), np.array([CO * L]), np.array([L]), b, k, th)[0])


def _fixed_tau_instrument():
    """second_order_volume_stable of F6.8e3 with dt' = dt'' = 0 (fixed tau). The source is copied and the two lines
    that set dt' and dt'' are replaced; each replacement must occur exactly once."""
    src = inspect.getsource(E3.second_order_volume_stable)
    for old, new in (("    dt1 = L * T1xy / dt\n", "    dt1 = 0 * T1xy\n"),
                     ("    dt2 = (abs(T1xy) ** 2 / 2 + 2 * L * T2xy) / dt - L * L * abs(T1xy) ** 2 / 2 / dt ** 3\n",
                      "    dt2 = 0 * T2xy\n")):
        assert src.count(old) == 1, old
        src = src.replace(old, new)
    src = src.replace("def second_order_volume_stable(", "def second_order_volume_fixed_tau(")
    ns = dict(vars(E3))
    exec(src, ns)
    return ns["second_order_volume_fixed_tau"]


def task(args):
    k, CO, L, tau, b, th, mode, lev = args
    if mode == "sigma":
        T = pair_T(k, CO, L, b, th) if b != 0 else L
        dt = math.sqrt(T * T + tau * tau)
    else:
        dt = math.sqrt(L * L + tau * tau)
        T = pair_T(k, CO, L, b, th) if b != 0 else L
        if T >= dt:
            return args, float("nan"), T          # the pair is not timelike at fixed tau: no interval
    V = volume(k, CO, L, dt, b, th, *lev)
    return args, V, T


def check_fermat():
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import f68f_dev as S0
    rows = []
    rng = np.random.default_rng(1)
    for _ in range(6):
        z1, z2 = rng.uniform(-2, 2, 2)
        X = rng.uniform(0.1, 3.0)
        b = rng.choice([0.01, 0.05, 0.16])
        th = rng.uniform(0, 2 * math.pi)
        k = 1.0
        t_ex, _ = S0.first_arrival(z1, z2, X, b, k=k, th=th)
        d = math.hypot(X, z2 - z1)
        t_f = float(optical_length(np.array([z1]), np.array([z2]), np.array([d]), b, k, th)[0])
        rows.append((z1, z2, X, b, t_ex, t_f, abs(t_f - t_ex) / t_ex))
    return rows


ROWS = [  # k, CO, L, tau
    (0.7, 0.4, 2.0, 1.2),     # the test diamond of F6.8c I1(b)
    (1.0, 0.0, 2.0, 1.0),
    (1.0, 1.0, 2.0, 1.0),
    (1.0, 0.7, 4.0, 1.5),
]
BS = [0.01, 0.02, 0.04, 0.08, 0.16]
N_TH = 8


def main():
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    t0 = time.time()
    lines = ["F6.8f-dev step 1a: exact interval volume of short intervals at finite b (development, decides nothing)", ""]
    lines.append("Fermat instrument against the exact Snell first arrival of step 0 (k = 1):")
    for r in check_fermat():
        lines.append("  z1 %+.3f z2 %+.3f X %.3f b %.2f: Snell %.15f  Fermat %.15f  rel %.1e" % r)
    lines.append("")
    lev = (16, 12, 16, None)
    ths = [2 * math.pi * j / N_TH for j in range(N_TH)]
    jobs = []
    for (k, CO, L, tau) in ROWS:
        jobs.append((k, CO, L, tau, 0.0, 0.0, "sigma", lev))
        for b in BS:
            for th in ths:
                for mode in ("sigma", "tau"):
                    jobs.append((k, CO, L, tau, b, th, mode, lev))
    lev2 = (24, 16, 24, None)
    k0, CO0, L0, tau0 = ROWS[0]
    jobs.append((k0, CO0, L0, tau0, 0.0, 0.0, "sigma", lev2))
    for b in (0.01, 0.04):
        for th in ths:
            jobs.append((k0, CO0, L0, tau0, b, th, "sigma", lev2))
    cdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "f68f_1a_cache")
    os.makedirs(cdir, exist_ok=True)

    def key(a):
        return "_".join(repr(x) for x in a[:7]) + "_" + "-".join(repr(x) for x in a[7][:3])

    done = {}
    for fn in os.listdir(cdir):
        with open(os.path.join(cdir, fn)) as fh:
            kk, V, T = fh.read().split()
            done[kk] = (float(V), float(T))
    todo = [a for a in jobs if key(a) not in done]
    print(f"jobs {len(jobs)}, cached {len(jobs) - len(todo)}, to compute {len(todo)}", flush=True)
    with Pool(workers) as pool:
        for a, V, T in pool.imap_unordered(task, todo, chunksize=1):
            kk = key(a)
            with open(os.path.join(cdir, f"{abs(hash(kk)) % 10 ** 12:012d}_{len(done)}.txt"), "w") as fh:
                fh.write(f"{kk} {V!r} {T!r}\n")
            done[kk] = (V, T)
    res_all = [(a,) + done[key(a)] for a in jobs]
    res = [x for x in res_all if x[0][7] == lev]
    res2 = [x for x in res_all if x[0][7] == lev2]
    fixed_tau = _fixed_tau_instrument()
    for (k, CO, L, tau) in ROWS:
        V0n = [V for (a, V, T) in res if a[:4] == (k, CO, L, tau) and a[4] == 0.0][0]
        V0 = math.pi * tau ** 4 / 24
        beta_s, gam_s = E3.second_order_volume_stable(k, CO, L, [tau * tau], 24, 16, 16)
        beta_t, gam_t = fixed_tau(k, CO, L, [tau * tau], 24, 16, 16)
        lines.append(f"row k = {k}, CO = {CO}, L = {L}, tau = {tau} (kL = {k * L:.2f}, k tau = {k * tau:.2f}); "
                     f"V0 numeric / exact - 1 = {V0n / V0 - 1:+.2e}")
        lines.append(f"  second order (F6.8e3 instrument): fixed sigma gamma {gam_s[0]:+.6e}, beta {beta_s[0]:.6e}; "
                     f"fixed tau gamma {gam_t[0]:+.6e}, beta {beta_t[0]:.6e}")
        lines.append("     b      eps=b(kL)^2  max|d(tau2)|/tau2   R_sigma(b)      R_tau(b)        B_sigma(b)                 B_tau(b)")
        for b in BS:
            out = {}
            for mode in ("sigma", "tau"):
                rr = [(a[5], V, T) for (a, V, T) in res if a[:4] == (k, CO, L, tau) and a[4] == b and a[6] == mode]
                Vs = np.array([V for (_, V, _) in rr])
                ths_ = np.array([t for (t, _, _) in rr])
                Ts = np.array([T for (_, _, T) in rr])
                R = 2 * (Vs.mean() / V0n - 1) / b ** 2
                B = 2 * (Vs * np.exp(-1j * ths_)).mean() / (b * V0n)
                out[mode] = (R, B, np.max(np.abs(L * L - Ts * Ts)) / tau ** 2)
            lines.append(f"  {b:6.3f}  {b * (k * L) ** 2:9.4f}   {out['sigma'][2]:.3e}        {out['sigma'][0]:+.6e}  "
                         f"{out['tau'][0]:+.6e}  {out['sigma'][1]:.6e}  {out['tau'][1]:.6e}")
        lines.append("")
    V02 = [V for (a, V, T) in res2 if a[4] == 0.0][0]
    for b in (0.01, 0.04):
        V1 = np.array([V for (a, V, T) in res if a[:4] == ROWS[0] and a[4] == b and a[6] == "sigma"])
        V2 = np.array([V for (a, V, T) in res2 if a[4] == b])
        V01 = [V for (a, V, T) in res if a[:4] == ROWS[0] and a[4] == 0.0][0]
        lines.append(f"level check, first row, fixed sigma, b = {b}: R_sigma {2 * (V1.mean() / V01 - 1) / b ** 2:+.6e} "
                     f"(levels {lev[:3]}) against {2 * (V2.mean() / V02 - 1) / b ** 2:+.6e} (levels {lev2[:3]})")
    lines.append(f"levels: n_sig, n_te, n_phi = {lev[:3]}, ray basis M = {M_RAY}, Q = {Q_RAY}, n_theta = {N_TH}; "
                 f"cost {time.time() - t0:.0f} s, {workers} processes")
    txt = "\n".join(lines)
    print(txt)
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "f68f_1a.txt"), "w") as fh:
        fh.write(txt + "\n")


if __name__ == "__main__":
    main()
