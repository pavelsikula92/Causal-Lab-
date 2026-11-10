"""
F3.14a -- Inertia of a localised defect of a field of the causal order (1+1, hypothesis H5)
==========================================================================================

Why this test exists (10 Nov 2026; card docs/H5_card_draft.md, revision 4, frozen by the author on 10 Nov, 12:05)
---------------------------------------------------------------------------------------------------------------
The author's choice (10 Nov): a locked defect must move as a whole through its own internal binding; no carrier flow
enters the base hypothesis (H5 of docs/cosmogenesis.md). F3.13c did not test this: its hedgehog moved as dust, every
world line kept its own velocity by the coordinate rule of F3.6, and the alignment never touched the motion. A
relaxational binding gives no inertia (defects of an overdamped dynamics move with v ~ force and stop with it), so the
author decided (D1-D4, D2'): a field phi on ALL events of the order, valued in an internal angle; a conservative,
second-order dynamics; in 1+1 implemented through the retarded Green function of the order (1/2 x the causal relation);
the sine-Gordon kink as the analytic reference. The BD operator in retarded form amplified the Poisson noise (probe
P-a: a constant field diverged within ~2 units) and stays a separate instrument. The potential is chosen here; its
derivation from the finite capacity of the flow is H5-D, later and separate (K0 plays no part).
A new branch: a defect of a field of the order (F3.13 tested a defect of an orientation carried by a flow).

Model (frozen)
--------------
1+1 Poisson sprinkling of density rho in the domain of dependence {0 <= t <= T, |x| <= X - t} of the initial slice, so
every computed value is exact for its order and no wall enters the dynamics.
    phi(x) = phi_0(x) - (1 / 2 rho) * sum over events y < x (the order), t_y >= 0, of [mu sin(phi(y)) + f chi(t_y)],
chi = 1 for t1 <= t < t2 (the kick), else 0; mu = 1; phi_0 the free wave of the initial data on t = 0 (d'Alembert's
formula with the exact kink of velocity v: the launch, the only place where coordinates enter the dynamics). Events are
computed in a linear extension of the order (helpers: phases/F3_structure/field1d.py; Fenwick engine, numba required).
Arms (all kink arms of one density share the sprinkling of that density and seed; the waves have their own):
  conservative: rest, v = 0.3, v = 0.6, kick (f = 0.05, 1 <= t < 3, from rest) -- at rho = 100, 200, 400;
  telegraph null (rho = 400): phi_tt + Gamma phi_t - phi_xx + sin phi = -f chi, Gamma = 2, 5, launched v = 0.3 and kick;
     evolved for psi = exp(Gamma t / 2) phi with the same sum over the causal past; uses the coordinate time of a medium
     frame (declared; a null model, not a candidate);
  V = 0 (rho = 400, launched v = 0.3): the sum vanishes, phi is the free wave of the initial data: a check of the readout;
  causality-blind (rho = 400, launched v = 0.3): every event sums over m_x events drawn uniformly with replacement from all
     earlier events (the same time range 0 <= t_y < t_x), m_x = the size of its causal past counted by the same engine
     with unit weights; the same normalisation, phi_0 and parameters;
  waves (all densities): vacuum (8 modes, k uniform in (0.25, 4], total amplitude 0.01) and plane waves k = 0.5, 1, 2, 4
     (amplitude 0.01), exact linear solutions with omega = sqrt(k^2 + 1).

Readouts (frozen; coordinates, declared): slabs of width 1 with centres 0.5 ... 9.5; kink fit
phi = 2 pi n0 + 4 arctan exp((x - x_c - v_s (t - t_m)) / w) on the events with |x - x_g| < 4, x_g = v_launch t_m;
Q = (mean phi on [x_g + 5, x_g + 7] - mean phi on [x_g - 7, x_g - 5]) / 2 pi; E = int [phi_t^2/2 + phi_x^2/2 + mu (1 -
cos phi)] dx over |x - x_c| < 5 from local regressions in cells 0.5 x 0.5, step 0.25; E_loc = the same over
|x - x_g| < 3 with the model's own mu (independent of the fit); velocities = slopes of x_c over the named slabs, p = gamma v;
the exact readouts (Q_ex, w_ex, E_ex) are the same readouts applied to the exact continuum solution at the same events.
A quantity that cannot be read (a failed fit in an arm where it is required) is nan and its condition is not met.

Criteria (frozen with the card, revision 4; deciding density rho = 400; 32 new sprinklings per density; SE = sd / sqrt(n)
over sprinklings, paired where stated; the SE of r_M by bootstrap over sprinklings, 2000 resamples; every condition
holds on its own, nothing is averaged into a score, nothing is chosen after the run)
-----------------------------------------------------------------------------------------------------------------------
INSTRUMENT_FAILED if any fails:
  I1  delta1 = R / R_c - 1 (R = rms phi last slab / first slab, R_c the same of the exact solution): |mean| + 2 SE < 0.1;
      mean rms(phi - exact) / rms(exact) in the last slab < 0.2;
  I2  for each k: delta2 = omega_meas / sqrt(k^2 + 1) - 1 (phase of the projection on cos kx, sin kx, slabs of 0.25):
      |mean| + 2 SE < 0.005;
  I3  eps = mean over the 10 slabs of |w - w_ex| / w_ex (rest arm): eps(100) - eps(400) > 2 sqrt(SE_100^2 + SE_400^2);
  I4  same sprinkling, initial data and parameters (mu 1, f 0.05 on [1, 3), v 0.3; 4 sprinklings, rho 100, X 12, T 6):
      linear extension t vs t + 0.5 x, Fenwick vs block engine, trapezoid X vs X - 2 on common events: each < 1e-12;
  H5-A (rest): |mean(Q - Q_ex)| + 2 SE < 0.02; drift x_c(slab 10) - x_c(slab 1): |mean| + 2 SE < 0.1;
      |mean(w/w_ex - 1)| + 2 SE < 0.1; |mean(E - E_ex)| / 8 + 2 SE < 0.05 (per-sprinkling slab means); M_rest = mean E;
  H5-E for each Gamma in {2, 5} and each arm (launched v 0.3, kick): p over slabs t in [5, 10) paired with the
      conservative arm: (p_cons - p_Gamma) - 2 SE > 0 and mean p_Gamma <= 0.5 mean p_cons;
  V = 0: L = E_loc(slab 10) / E_loc(slab 1); D_L = L_cons(v 0.3) - L_V0: mean - 2 SE > 0.2.
FAILED if any: H5-B (v 0.3 or 0.6): |mean v_meas - v| - 2 SE > 0.03 (v_meas over all 10 slabs), |mean C - 1| - 2 SE >
  0.1 (C = w(v) gamma / w(rest), paired), |mean e_B| - 2 SE > 0.05 (e_B = (E(v) - E_ex(v)) / (8 gamma)); H5-C: |r_M - 1|
  - 2 SE > 0.1, r_M = M_in / M_rest, M_in = 2 pi Q_kick f tau / mean(Delta p), Delta p = gamma v(kick) - gamma v(rest) over
  slabs t in [4, 10), paired.
PASSED if all: H5-B at both v: |mean v_meas - v| + 2 SE < 0.03, |mean C - 1| + 2 SE < 0.1, |mean e_B| + 2 SE < 0.05,
  |mean(Q - Q_ex)| + 2 SE < 0.02; H5-C: |r_M - 1| + 2 SE < 0.1; blind control separated: D_L = L_cons - L_blind,
  mean - 2 SE > 0.2.
UNDECIDED otherwise (if the blind control does not separate while H5-B and H5-C pass: the sum, not the causal structure,
  carries the kink).
Reported only: rho = 100 and 200 (r_M trend, eps(200)), the free motion after the kick (gamma v over [4, 7) vs [7, 10)),
the field momentum P, widths and velocities of the controls.

Reading declared now (card, section 9)
-------------------------------------
PASSED (1+1): the order's own retarded Green function gives a defect with inertia whose inertial mass equals the energy
of its configuration, and a relaxational dynamics does not; matter's inertia is then a property of a conservative
dynamics of the order, not of a carrier flow. A limit stated before the run: in a Lorentz-invariant sprinkling with a
covariant rule a moving kink is in distribution a boosted resting kink, so H5-B is largely fixed by symmetry once H5-A
holds; not fixed by symmetry: that the discrete sum is stable (BD is not), the size of the discreteness corrections, the
separation from relaxation. FAILED: the discreteness breaks inertia or the equality of the masses at this scale. Limits:
1+1 (no transverse directions), the potential chosen (H5-D derives it), the internal target space is an input (D1),
the launch uses coordinates, nothing about gravity (F5.1i).
"""

from __future__ import annotations

import math
import time

import numpy as np

from core.order.sprinkle import derive_seed
from core.stats.parallel import parallel_map
from core.testing import Result
from phases.F3_structure import field1d as F

TEST = {
    "id": "F3.14a",
    "phase": "F3",
    "title": "Inertia of a localised defect of a field of the causal order (1+1, sine-Gordon kink, Green function of "
             "the order)",
    "theory_part": "F3 emergence of structure; matter as a locked defect with inertia from its own binding (H5 of "
                   "docs/cosmogenesis.md); a field of the order on all events; conservative dynamics (decisions D1-D4, D2')",
    "question": "Does a topological defect of a field on all events of a 1+1 causal order, evolved by the retarded Green "
                "function of the order (1/2 x the causal relation) without any source after its launch, keep its "
                "identity at rest, keep its velocity and Lorentz contraction when launched, and answer a finite force with "
                "an inertial mass equal to its rest energy, while a relaxational dynamics and a causality-blind sum do not?",
    "method": "1+1 sprinkling in the domain of dependence of the initial slice; phi(x) = phi_0(x) - (1/2 rho) sum over the "
              "causal past of [sin phi + f chi]; sine-Gordon kink launched at v = 0, 0.3, 0.6 and kicked (f 0.05, tau 2); "
              "readouts by kink fits and local regressions per slab; instrument I1-I4, H5-A, telegraph null (Gamma 2, 5), "
              "V = 0 readout check; causality-blind control; 32 sprinklings per density, deciding density 400",
    "prediction": "continuum sine-Gordon: velocity kept, w gamma = 1, E = 8 gamma, M_in = M_rest = 8; the telegraph kink "
                  "stops (dP/dt = -Gamma P); the blind sum loses localisation",
    "criterion": "INSTRUMENT: I1 |delta1| + 2 SE < 0.1 and err < 0.2; I2 |delta2| + 2 SE < 0.005 each k; I3 eps(100) - "
                 "eps(400) > 2 SE; I4 < 1e-12; H5-A Q 0.02, drift 0.1, w 0.1, E/8 0.05 (+ 2 SE); H5-E each Gamma and arm "
                 "separated and p_Gamma <= 0.5 p_cons; V0 D_L - 2 SE > 0.2. PASSED: H5-B both v |dv| + 2 SE < 0.03, "
                 "|C - 1| + 2 SE < 0.1, |e_B| + 2 SE < 0.05, Q; H5-C |r_M - 1| + 2 SE < 0.1; blind D_L - 2 SE > 0.2. "
                 "FAILED: |dv| - 2 SE > 0.03, |C - 1| - 2 SE > 0.1, |e_B| - 2 SE > 0.05 or |r_M - 1| - 2 SE > 0.1. "
                 "Otherwise UNDECIDED (card docs/H5_card_draft.md rev. 4, section 8)",
    "null_model": "telegraph (overdamped) dynamics with the same potential, initial data and kick; V = 0 (readout check); "
                  "causality-blind sum with the same number of terms; exact continuum solutions",
    "level": "mechanism",
    "meters": [],
    "depends_on": [],
    "params": {"rhos": [100.0, 200.0, 400.0], "rho_dec": 400.0, "rho_lo": 100.0, "n_seeds": 32, "X_kink": 26.0,
               "X_wave": 20.0, "T": 10.0, "vs": [0.3, 0.6], "f": 0.05, "t1": 1.0, "t2": 3.0, "gammas": [2.0, 5.0],
               "ks": [0.5, 1.0, 2.0, 4.0], "n_modes": 8, "amp": 0.01, "k_lo": 0.25, "k_hi": 4.0,
               "i4_n": 4, "i4_rho": 100.0, "i4_X": 12.0, "i4_T": 6.0, "i4_v": 0.3, "i4_dX": 2.0, "i4_tol": 1e-12,
               "tol_I1": 0.1, "tol_I1_err": 0.2, "tol_I2": 0.005, "tol_Q": 0.02, "tol_drift": 0.1, "tol_w": 0.1,
               "tol_E": 0.05, "tol_v": 0.03, "tol_C": 0.1, "tol_M": 0.1, "drop": 0.5, "tol_DL": 0.2,
               "t_kick_read": 4.0, "t_null_read": 5.0, "n_se": 2.0, "boot_B": 2000},
    "quick": {"rhos": [50.0, 100.0], "rho_dec": 100.0, "rho_lo": 50.0, "n_seeds": 2, "i4_n": 1, "boot_B": 200},
    "estimate": "minutes with numba (pip install numba; required): ~170 000 events per kink arm at rho 400; the "
                "causality-blind arm is the slowest (~1 min each); 32 x (3 x 4 + 6) kink tasks and 32 x 3 x 5 wave tasks",
}


# ============================================================================================== tasks

def _momentum(t, x, phi, tm, xc, half=5.0, cell=0.5, step=0.25):
    """P = -int phi_t phi_x dx with the regressions of F.energy (reported only)."""
    s = np.abs(t - tm) < 0.5 * cell
    ts, xs, ps = t[s], x[s], phi[s]
    o = np.argsort(xs)
    ts, xs, ps = ts[o], xs[o], ps[o]
    P, used = 0.0, 0
    for xg in np.arange(xc - half, xc + half + 1e-9, step):
        i0, i1 = np.searchsorted(xs, [xg - 0.5 * cell, xg + 0.5 * cell])
        if i1 - i0 < 6:
            continue
        A = np.column_stack([np.ones(i1 - i0), ts[i0:i1] - tm, xs[i0:i1] - xg])
        c, *_ = np.linalg.lstsq(A, ps[i0:i1], rcond=None)
        P += -c[1] * c[2] * step
        used += 1
    n_all = int(round(2 * half / step)) + 1
    return P * n_all / used if used else float("nan")


def _track(t, x, phi, v, T, mu, exact):
    rows = []
    for a in np.arange(0.0, T - 1.0 + 1e-9, 1.0):
        tm = a + 0.5
        guess = v * tm
        row = {"t": float(tm), "E_loc": F.energy(t, x, phi, tm, guess, mu, half=3.0)}
        r = F.fit_kink(t, x, phi, a, a + 1.0, guess, mu)
        if r is None:
            row["fail"] = True
            rows.append(row)
            continue
        row.update({k: r[k] for k in ("xc", "vs", "w", "res", "Q")})
        row["E"] = F.energy(t, x, phi, tm, r["xc"], mu)
        row["P"] = _momentum(t, x, phi, tm, r["xc"])
        if exact is not None:
            rx = F.fit_kink(t, x, exact, a, a + 1.0, guess, mu)
            row["Q_ex"] = rx["Q"] if rx else float("nan")
            row["w_ex"] = rx["w"] if rx else float("nan")
            row["E_ex"] = F.energy(t, x, exact, tm, r["xc"], mu)
        rows.append(row)
    return rows


def task_kink(p, rho, k, arm, base):
    """arm: rest, v0.3, v0.6, kick; tele<G>_v0.3, tele<G>_kick; V0_v0.3; blind_v0.3."""
    X, T = float(p["X_kink"]), float(p["T"])
    t, x = F.sprinkle_trapezoid(rho, X, T, derive_seed(base, "F3.14a", rho, k))
    variant, motion = arm.split("_") if "_" in arm else ("cons", arm)
    v = 0.0 if motion in ("rest", "kick") else float(motion[1:])
    kw = dict(f=float(p["f"]), t1=float(p["t1"]), t2=float(p["t2"])) if motion == "kick" else {}
    mu = 0.0 if variant == "V0" else F.MU
    t0 = time.time()
    if variant.startswith("tele"):
        gam = float(variant[4:])
        phi = F.evolve(t, x, F.kink_free_telegraph(t, x, v, gam), rho, mu=mu, gam=gam, **kw)
    elif variant == "blind":
        m = F.evolve(t, x, np.zeros_like(t), rho, unit=True)
        phi = F._blind_evolve(t, F.kink_free(t, x, v), np.rint(m).astype(np.int64), float(rho), mu,
                              float(kw.get("f", 0.0)), float(kw.get("t1", 0.0)), float(kw.get("t2", 0.0)),
                              int(derive_seed(base, "F3.14a-blind", rho, k) % (2 ** 31)))
    else:
        phi = F.evolve(t, x, F.kink_free(t, x, v), rho, mu=mu, **kw)
    exact = F.kink_exact(t, x, v) if (variant == "cons" and motion != "kick") else None
    rows = _track(t, x, phi, v, T, mu, exact)
    return {"kind": "kink", "rho": rho, "k": k, "arm": arm, "variant": variant, "motion": motion, "v": v,
            "N": int(t.size), "sec": time.time() - t0, "rows": rows}


def task_wave(p, rho, k, what, base):
    X, T = float(p["X_wave"]), float(p["T"])
    sd = derive_seed(base, "F3.14a-wave", rho, k)
    t, x = F.sprinkle_trapezoid(rho, X, T, sd)
    out = {"kind": "wave", "rho": rho, "k": k, "what": what, "N": int(t.size)}
    if what == "vacuum":
        rng = np.random.default_rng(derive_seed(sd, "modes"))
        n = int(p["n_modes"])
        ks = rng.uniform(float(p["k_lo"]), float(p["k_hi"]), n)
        modes = [(float(p["amp"]) / math.sqrt(n / 2.0), float(kk), float(rng.uniform(0, 2 * math.pi))) for kk in ks]
        ex = F.modes_exact(t, x, modes)
        phi = F.evolve(t, x, F.modes_free(t, x, modes), rho)
        def slab(a):
            s = (t >= a) & (t < a + 1.0)
            return (float(np.sqrt(np.mean(phi[s] ** 2))), float(np.sqrt(np.mean(ex[s] ** 2))),
                    float(np.sqrt(np.mean((phi[s] - ex[s]) ** 2))))
        f0, f1 = slab(0.0), slab(T - 1.0)
        out.update({"R": f1[0] / f0[0], "R_c": f1[1] / f0[1], "err_last": f1[2] / f1[1]})
    else:
        kk = float(what[1:])
        phi = F.evolve(t, x, F.modes_free(t, x, [(float(p["amp"]), kk, 0.0)]), rho)
        ph, tm = [], []
        for a in np.arange(0.0, T - 0.25 + 1e-9, 0.25):
            ph.append(F.project(t, x, phi, kk, a, a + 0.25)[0])
            tm.append(a + 0.125)
        om = float(np.polyfit(np.array(tm), np.unwrap(np.array(ph)), 1)[0])
        out.update({"kw": kk, "omega": om, "omega_c": math.sqrt(kk * kk + F.MU)})
    return out


def task_i4(p, k, base):
    rho, X, T, v = float(p["i4_rho"]), float(p["i4_X"]), float(p["i4_T"]), float(p["i4_v"])
    t, x = F.sprinkle_trapezoid(rho, X, T, derive_seed(base, "F3.14a-i4", k))
    kw = dict(f=float(p["f"]), t1=float(p["t1"]), t2=float(p["t2"]))
    phi0 = F.kink_free(t, x, v)
    ref = F.evolve(t, x, phi0, rho, **kw)
    sc = float(np.max(np.abs(ref)))
    alt = F.evolve(t, x, phi0, rho, order="boost", **kw)
    blk = F.evolve(t, x, phi0, rho, engine="block", block=256, **kw)
    sub = np.abs(x) <= (X - float(p["i4_dX"])) - t
    phs = F.evolve(t[sub], x[sub], phi0[sub], rho, **kw)
    return {"kind": "i4", "k": k, "N": int(t.size), "order": float(np.max(np.abs(alt - ref)) / sc),
            "engine": float(np.max(np.abs(blk - ref)) / sc), "region": float(np.max(np.abs(phs - ref[sub])) / sc)}


def run_task(kind, *a):
    return {"kink": task_kink, "wave": task_wave, "i4": task_i4}[kind](*a)


# ============================================================================================== analysis

def _good(r):
    return not r.get("fail", False)


def _col(rows, key):
    return np.array([r[key] if _good(r) and key in r else np.nan for r in rows], float)


def _slope(rows, tmin=0.0, tmax=1e9):
    tt = np.array([r["t"] for r in rows])
    xc = _col(rows, "xc")
    m = (tt >= tmin) & (tt < tmax)
    if m.sum() < 2 or not np.all(np.isfinite(xc[m])):
        return float("nan")
    return float(np.polyfit(tt[m], xc[m], 1)[0])


def _gv(v):
    return v / math.sqrt(1.0 - v * v) if np.isfinite(v) and abs(v) < 1.0 else float("nan")


def _ms(a):
    a = np.asarray(a, float)
    if a.size < 2 or not np.all(np.isfinite(a)):
        return float("nan"), float("nan")
    return float(a.mean()), float(a.std(ddof=1) / math.sqrt(a.size))


def _fmt(m):
    return f"{m[0]:+.4f} +- {m[1]:.4f}"


def run(p, ctx):
    if F._nb is None:
        raise SystemExit("F3.14a needs numba (pip install numba): the telegraph null, the past counts and the speed")
    t0 = time.time()
    n, nse = int(p["n_seeds"]), float(p["n_se"])
    rhos, rd, rlo = [float(r) for r in p["rhos"]], float(p["rho_dec"]), float(p["rho_lo"])
    base = ctx.seed
    jobs = []
    for rho in rhos:
        arms = ["rest"] + [f"v{v:g}" for v in p["vs"]] + ["kick"]
        if rho == rd:
            arms += [f"tele{g:g}_v{p['vs'][0]:g}" for g in p["gammas"]] + [f"tele{g:g}_kick" for g in p["gammas"]]
            arms += [f"V0_v{p['vs'][0]:g}", f"blind_v{p['vs'][0]:g}"]
        for k in range(n):
            for arm in arms:
                jobs.append(("kink", p, rho, k, arm, base))
            for what in ["vacuum"] + [f"k{kk:g}" for kk in p["ks"]]:
                jobs.append(("wave", p, rho, k, what, base))
    for k in range(int(p["i4_n"])):
        jobs.append(("i4", p, k, base))
    jobs.sort(key=lambda j: (j[0] != "kink" or "blind" not in j[4], -(j[2] if j[0] != "i4" else 0)))
    res = parallel_map(run_task, jobs, ctx.workers, ctx.progress,
                       lambda i, r, d, nn: f"[F3.14a] {d}/{nn}: {r['kind']} {r.get('arm', r.get('what', ''))} rho "
                                           f"{r.get('rho', '')} ({time.time() - t0:.0f}s)")

    def K(rho, arm):
        out = {r["k"]: r for r in res if r["kind"] == "kink" and r["rho"] == rho and r["arm"] == arm}
        return [out[k] for k in range(n)]

    def W(rho, what):
        out = {r["k"]: r for r in res if r["kind"] == "wave" and r["rho"] == rho and r["what"] == what}
        return [out[k] for k in range(n)]

    lines = [f"1+1, trapezoid |x| <= X - t, T {p['T']:g}; {n} sprinklings per density {rhos}; deciding density {rd:g}; "
             f"Green function of the order; engine Fenwick (numba)"]
    inst, passc, failc = {}, {}, {}

    # ---------------------------------------------------------------- I4
    i4 = [r for r in res if r["kind"] == "i4"]
    for key in ("order", "engine", "region"):
        mx = max(r[key] for r in i4)
        inst[f"I4 {key} < {p['i4_tol']:g}"] = mx < float(p["i4_tol"])
        lines.append(f"  I4 {key}: largest relative difference {mx:.2e} ({len(i4)} sprinklings, N {i4[0]['N']})")

    # ---------------------------------------------------------------- waves (I1, I2)
    for rho in rhos:
        vac = W(rho, "vacuum")
        d1 = _ms([r["R"] / r["R_c"] - 1.0 for r in vac])
        err = float(np.mean([r["err_last"] for r in vac]))
        lines.append(f"  rho {rho:g} I1: delta1 {_fmt(d1)}; error last slab {err:.3f}")
        if rho == rd:
            inst["I1 amplitude ratio"] = abs(d1[0]) + nse * d1[1] < float(p["tol_I1"])
            inst["I1 error last slab"] = err < float(p["tol_I1_err"])
        for kk in p["ks"]:
            ws = W(rho, f"k{kk:g}")
            d2 = _ms([r["omega"] / r["omega_c"] - 1.0 for r in ws])
            lines.append(f"  rho {rho:g} I2 k {kk:g}: delta2 {d2[0]:+.5f} +- {d2[1]:.5f}")
            if rho == rd:
                inst[f"I2 k {kk:g}"] = abs(d2[0]) + nse * d2[1] < float(p["tol_I2"])

    # ---------------------------------------------------------------- per density: conservative arms
    def eps_of(rho):
        return [float(np.nanmean(np.abs(_col(r["rows"], "w") - _col(r["rows"], "w_ex")) / _col(r["rows"], "w_ex")))
                if np.all(np.isfinite(_col(r["rows"], "w"))) else float("nan") for r in K(rho, "rest")]

    eps = {rho: _ms(eps_of(rho)) for rho in rhos}
    lines.append("  I3 eps(rho): " + ", ".join(f"{rho:g}: {_fmt(eps[rho])}" for rho in rhos))
    inst["I3 convergence"] = (eps[rlo][0] - eps[rd][0]) > nse * math.sqrt(eps[rlo][1] ** 2 + eps[rd][1] ** 2)

    f_, tau = float(p["f"]), float(p["t2"]) - float(p["t1"])
    tk, tn = float(p["t_kick_read"]), float(p["t_null_read"])
    rM = {}
    out = {}
    for rho in rhos:
        rest = K(rho, "rest")
        Qd = _ms([np.mean(_col(r["rows"], "Q") - _col(r["rows"], "Q_ex")) for r in rest])
        drift = _ms([_col(r["rows"], "xc")[-1] - _col(r["rows"], "xc")[0] for r in rest])
        wr = _ms([np.mean(_col(r["rows"], "w")) / np.mean(_col(r["rows"], "w_ex")) - 1.0 for r in rest])
        Ed = _ms([(np.mean(_col(r["rows"], "E")) - np.mean(_col(r["rows"], "E_ex"))) / 8.0 for r in rest])
        Erest = np.array([np.mean(_col(r["rows"], "E")) for r in rest])
        lines.append(f"  rho {rho:g} H5-A rest: Q - Q_ex {_fmt(Qd)}; drift {_fmt(drift)}; w/w_ex - 1 {_fmt(wr)}; "
                     f"(E - E_ex)/8 {_fmt(Ed)}; M_rest = E {np.mean(Erest):.4f}")
        if rho == rd:
            inst["H5-A Q"] = abs(Qd[0]) + nse * Qd[1] < float(p["tol_Q"])
            inst["H5-A drift"] = abs(drift[0]) + nse * drift[1] < float(p["tol_drift"])
            inst["H5-A width"] = abs(wr[0]) + nse * wr[1] < float(p["tol_w"])
            inst["H5-A energy"] = abs(Ed[0]) + nse * Ed[1] < float(p["tol_E"])
        wrest = np.array([np.mean(_col(r["rows"], "w")) for r in rest])
        for v in p["vs"]:
            g = 1.0 / math.sqrt(1.0 - v * v)
            arm = K(rho, f"v{v:g}")
            dv = _ms([_slope(r["rows"]) - v for r in arm])
            C = _ms([np.mean(_col(r["rows"], "w")) * g / wrest[i] - 1.0 for i, r in enumerate(arm)])
            eB = _ms([(np.mean(_col(r["rows"], "E")) - np.mean(_col(r["rows"], "E_ex"))) / (8.0 * g) for r in arm])
            Qv = _ms([np.mean(_col(r["rows"], "Q") - _col(r["rows"], "Q_ex")) for r in arm])
            lines.append(f"  rho {rho:g} H5-B v {v:g}: v_meas - v {_fmt(dv)}; C - 1 {_fmt(C)}; e_B {_fmt(eB)}; "
                         f"Q - Q_ex {_fmt(Qv)}")
            if rho == rd:
                passc[f"H5-B v {v:g} velocity"] = abs(dv[0]) + nse * dv[1] < float(p["tol_v"])
                passc[f"H5-B v {v:g} contraction"] = abs(C[0]) + nse * C[1] < float(p["tol_C"])
                passc[f"H5-B v {v:g} energy"] = abs(eB[0]) + nse * eB[1] < float(p["tol_E"])
                passc[f"H5-B v {v:g} identity Q"] = abs(Qv[0]) + nse * Qv[1] < float(p["tol_Q"])
                failc[f"H5-B v {v:g} velocity"] = abs(dv[0]) - nse * dv[1] > float(p["tol_v"])
                failc[f"H5-B v {v:g} contraction"] = abs(C[0]) - nse * C[1] > float(p["tol_C"])
                failc[f"H5-B v {v:g} energy"] = abs(eB[0]) - nse * eB[1] > float(p["tol_E"])
        kick = K(rho, "kick")
        dp = np.array([_gv(_slope(a["rows"], tk)) - _gv(_slope(b["rows"], tk)) for a, b in zip(kick, rest)])
        Qk = np.array([np.mean(_col(r["rows"], "Q")) for r in kick])

        def rm(idx):
            return 2.0 * math.pi * np.mean(Qk[idx]) * f_ * tau / np.mean(dp[idx]) / np.mean(Erest[idx])

        r0 = rm(np.arange(n)) if np.all(np.isfinite(dp)) and np.all(np.isfinite(Erest)) else float("nan")
        if np.isfinite(r0):
            rng = np.random.default_rng(derive_seed(base, "F3.14a-boot", rho))
            bs = [rm(rng.integers(0, n, n)) for _ in range(int(p["boot_B"]))]
            se = float(np.std(bs, ddof=1))
        else:
            se = float("nan")
        rM[rho] = (r0, se)
        free = _ms([_gv(_slope(r["rows"], tk, 7.0)) - _gv(_slope(r["rows"], 7.0)) for r in kick])
        P_after = float(np.nanmean([np.nanmean(_col(r["rows"], "P")[int(tk):]) for r in kick]))
        lines.append(f"  rho {rho:g} H5-C kick: Delta p {_fmt(_ms(dp))} (continuum M = 8: {2 * math.pi * f_ * tau / 8:.4f}); "
                     f"Q_kick {np.mean(Qk):.4f}; M_in {2 * math.pi * np.mean(Qk) * f_ * tau / np.mean(dp):.4f}; "
                     f"r_M = M_in / M_rest {r0:.4f} +- {se:.4f}; reported: gamma v [4,7) - [7,10) {_fmt(free)}, "
                     f"field momentum P after the kick {P_after:+.4f} (continuum {2 * math.pi * f_ * tau:.4f})")
        out[rho] = {"eps": eps[rho], "rM": rM[rho]}
    passc["H5-C r_M"] = abs(rM[rd][0] - 1.0) + nse * rM[rd][1] < float(p["tol_M"])
    failc["H5-C r_M"] = abs(rM[rd][0] - 1.0) - nse * rM[rd][1] > float(p["tol_M"])

    # ---------------------------------------------------------------- controls at the deciding density
    v0 = float(p["vs"][0])
    cons_v = K(rd, f"v{v0:g}")
    cons_k = K(rd, "kick")
    for g in p["gammas"]:
        for motion, cons in ((f"v{v0:g}", cons_v), ("kick", cons_k)):
            tel = K(rd, f"tele{g:g}_{motion}")
            pc = np.array([_gv(_slope(r["rows"], tn)) for r in cons])
            pg = np.array([_gv(_slope(r["rows"], tn)) for r in tel])
            d = _ms(pc - pg)
            ok = (d[0] - nse * d[1] > 0.0) and (np.mean(pg) <= float(p["drop"]) * np.mean(pc))
            inst[f"H5-E Gamma {g:g} {motion}"] = bool(ok)
            lines.append(f"  H5-E Gamma {g:g} {motion}: p_Gamma {np.mean(pg):+.4f}, p_cons {np.mean(pc):+.4f}, "
                         f"p_cons - p_Gamma {_fmt(d)}")

    def L(r):
        e = np.array([rr["E_loc"] for rr in r["rows"]], float)
        return e[-1] / e[0]

    Lc = np.array([L(r) for r in cons_v])
    for name, key in (("V = 0 readout check", f"V0_v{v0:g}"), ("blind control", f"blind_v{v0:g}")):
        ctl = K(rd, key)
        Lx = np.array([L(r) for r in ctl])
        d = _ms(Lc - Lx)
        wl = float(np.nanmean([np.nanmean(_col(r["rows"], "w")[-3:]) for r in ctl]))
        vl = float(np.nanmean([_slope(r["rows"]) for r in ctl]))
        lines.append(f"  {name}: L {np.mean(Lx):.4f} (conservative {np.mean(Lc):.4f}); D_L {_fmt(d)}; reported: width "
                     f"(last 3 slabs) {wl:.2f}, velocity {vl:+.3f}")
        sep = d[0] - nse * d[1] > float(p["tol_DL"])
        if key.startswith("V0"):
            inst["V = 0 readout check"] = bool(sep)
        else:
            passc["blind control separated"] = bool(sep)

    # ---------------------------------------------------------------- verdict
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
    lines.append(f"({time.time() - t0:.0f}s)")
    return Result(status, lines, {"instrument": inst, "passed": passc, "failed": failc, "per_density": out,
                                  "results": res})
