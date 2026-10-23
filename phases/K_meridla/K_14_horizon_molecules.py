"""
K.14 -- Horizon molecules: does counting links across a horizon measure its area?
==================================================================================

Why this test exists (author's decision of 22 Oct 2026)
-------------------------------------------------------
The Jacobson route (docs/derivations/jacobson_cesta.md): if every local
causal horizon carries an entropy proportional to its area, with one
universal constant, and heat crossing it obeys dQ = T dS (T the Unruh
temperature), then Einstein's equations follow as an equation of state
(Jacobson 1995). Gravity would then be derived, not put in: no field
equation, no f(Pi), no H1 ansatz, and gamma = 1 comes with GR. The lab
can supply the first premise from the causal set itself: the entropy of a
horizon as the number of HORIZON MOLECULES (Dou-Sorkin; the 1-molecule of
Barton et al. 2019, core.order.horizon). This test calibrates that
instrument:
  (1) in flat space the count per unit area equals the constant a^(d)
      computed independently by quadrature, at every density (the scaling
      rho^((d-2)/d) is exact there);
  (2) in a curved (conformally flat) geometry it counts the PROPER area of
      the horizon section, not the coordinate area and not the volume: the
      universality the Jacobson argument needs.

Instrument (core.order.horizon)
-------------------------------
Sigma = {t = 0}, horizon H = the null plane t = x (outside: x - t > 0).
A molecule: p- outside, p+ inside, both below Sigma, p+ the ONLY element
below Sigma in the future of p-. Sprinkling: a true Poisson process
(core.order.sprinkle.sprinkle_poisson) in a box below Sigma, extended
beyond the counting window by the height h_cut below which the chance of a
lower element is < e^-30 (so the future cone of every candidate is fully
sprinkled). a^(d) by adaptive quadrature of the module formula: a^(2) = 1/3,
a^(3) = 0.218853, a^(4) = 0.173205 (= sqrt 3 / 10).

Parts
-----
flat_d2, flat_d3, flat_d4 -- flat space at two densities rho = 1 and 4 (per
    coordinate volume; the window is the same in coordinates, so the
    expected count changes by 4^((d-2)/d)). d = 2: J is a point, many small
    sprinklings; d = 3, 4: a transverse window of side L.
    Criterion per density: z = (mean count - a^(d) rho^((d-2)/d) A) / SE
    (SE over the sprinklings); PASSED if |z| <= z_max = 3.5 for both
    densities, otherwise FAILED (an instrument fault, no physics).
conf_d3, conf_d4 -- conformally flat geometry Pi(y)^2 eta with
    Pi(y) = 0.75 + 0.25 cos(2 pi y / lambda), lambda = 40 (in units of the
    discreteness length at Pi = 1; Pi from 0.5 to 1, so the local
    discreteness length l / Pi stays <= 2, twenty times below lambda). The
    order is flat (null plane and Sigma unchanged, the horizon section has
    no expansion); the density is rho Pi^d per coordinate volume. The proper
    area element of J is Pi^(d-2) dA. Counts in 16 bins of y over two
    periods; model count_bin = C int_bin Pi^k dA, k and C by least squares
    on the bin means (weights from the spread over sprinklings), SE of k by
    bootstrap over sprinklings.
    Alternatives: k = 0 (coordinate area) and k = d (the volume density);
    both at distance 2 from k = d - 2 in d = 4, and at 1 and 2 in d = 3.
    Criterion (quarter rule): PASSED if |k - (d-2)| + 2 SE <= d_alt / 4,
    FAILED if |k - (d-2)| - 2 SE > d_alt / 4, UNDECIDED otherwise, d_alt the
    distance to the nearer alternative (d = 4: 2; d = 3: 1). Reported, not
    decided: C / a^(d) (the parameter-free leading-order prediction is 1;
    finite-rho corrections are a derivative expansion in l / lambda) and the
    chi^2 of the parameter-free model.

Seen before this test was written (development runs, 22 Oct 2026, decide
nothing): flat d = 2, 3, 4 at rho = 1: 0.321 +- 0.009, 0.213 +- 0.014,
0.175 +- 0.003 against a^(d) = 0.333, 0.219, 0.173; conformal d = 4 (8
sprinklings): the model k = 2 with C = a fits (chi^2 19 for 16 bins), k = 0
and k = 4 do not (chi^2 1465 and 172).

Level: calibration (an instrument for the entropy of a horizon; no claim
about gravity until the Jacobson step is tested, docs/derivations/jacobson_cesta.md).
"""

from __future__ import annotations

import math
import time

import numpy as np

from core.order.horizon import count_molecules, h_cut, molecule_constant
from core.order.sprinkle import Box, derive_seed, sprinkle_poisson
from core.stats.parallel import parallel_map
from core.testing import Result

TEST = {
    "id": "K.14",
    "phase": "K",
    "title": "Molekuly horizontu: měří počet vazeb přes horizont jeho plochu?",
    "theory_part": "měřidlo entropie horizontu pro Jacobsonovu cestu (gravitace jako stavová rovnice); bez závěru o gravitaci",
    "question": "Je počet molekul horizontu (Barton et al. 2019) v plochém prostoročasu roven konstantě a^(d) krát plocha "
                "v jednotkách diskrétnosti a v konformně plochém prostoročasu úměrný VLASTNÍ ploše řezu horizontu?",
    "method": "Poissonův sprinkling pod plochou Σ = {t = 0}, horizont nadrovina t = x; molekula = p− vně, p+ uvnitř, p+ jediný "
              "prvek pod Σ v budoucnosti p−; plochý prostor ve 2, 3, 4 rozměrech při dvou hustotách proti a^(d) z kvadratury; "
              "konformní Π(y) = 0,75 + 0,25 cos(2πy/40): sklon k počtu v pásech proti Π",
    "prediction": "plochý: a^(2) = 1/3, a^(3) = 0,2189, a^(4) = 0,1732 na jednotku plochy v jednotkách ℓ^(d−2); konformní: k = d − 2",
    "criterion": "flat: |z| <= 3,5 při obou hustotách (jinak NEPROŠEL = vada měřidla). conf: PROŠEL |k − (d−2)| + 2 SE <= d_alt/4, "
                 "NEPROŠEL |k − (d−2)| − 2 SE > d_alt/4, jinak NEROZHODNUTO; alternativy k = 0 (souřadnicová plocha) a k = d "
                 "(objemová hustota), d_alt = vzdálenost bližší (d = 4: 2, d = 3: 1). C/a^(d) a χ² bez volných parametrů hlášeny.",
    "null_model": "alternativy k = 0 a k = d",
    "level": "calibration",
    "meters": [],
    "depends_on": [],
    "params": {"densities": [1.0, 4.0], "z_max": 3.5, "lam": 40.0, "Pi_mid": 0.75, "Pi_amp": 0.25, "periods": 2,
               "n_bins": 16, "L_z": 20.0, "boot_B": 2000, "n_se": 2.0, "tail": 30.0, "chunk": 1},
    "parts": {"flat_d2": {"mode": "flat", "d": 2, "n_spr": 20000, "L": 0.0, "chunk": 500},
              "flat_d3": {"mode": "flat", "d": 3, "n_spr": 64, "L": 60.0},
              "flat_d4": {"mode": "flat", "d": 4, "n_spr": 48, "L": 40.0},
              "conf_d3": {"mode": "conf", "d": 3, "n_spr": 400},
              "conf_d4": {"mode": "conf", "d": 4, "n_spr": 64}},
    "quick": {"densities": [1.0, 4.0], "boot_B": 200},
    "quick_parts": {"flat_d2": {"n_spr": 400, "chunk": 100}, "flat_d3": {"n_spr": 4, "L": 20.0},
                    "flat_d4": {"n_spr": 3, "L": 12.0}, "conf_d3": {"n_spr": 20, "lam": 20.0},
                    "conf_d4": {"n_spr": 3, "lam": 16.0, "L_z": 8.0}},
    "estimate": "flat části po minutách; conf_d4 ~6 s na sprinkling (64 sprinklingů), conf_d3 rychlejší; vše s --workers",
}


# ============================================================================================== tasks

def flat_task(p, seeds, rho):
    """Molecule counts of flat sprinklings at density rho (one per seed); window side L (d >= 3)."""
    d, L = int(p["d"]), float(p["L"])
    hm = h_cut(d, rho, float(p["tail"]))
    out = []
    for sd in seeds:
        rng = np.random.default_rng(sd)
        lo = [-hm, -2.0 * hm] + [0.0] * (d - 2)
        hi = [0.0, 2.0 * hm] + [L + 2.0 * hm] * (d - 2)
        X = sprinkle_poisson(Box(d, lo, hi), rho, rng)
        X = X[X[:, 0] < 0.0]
        r = count_molecules(X, hm, [hm] * (d - 2), [L + hm] * (d - 2))
        out.append(int(r["lower"].size))
    return out


def Pi_of(p):
    lam, mid, amp = float(p["lam"]), float(p["Pi_mid"]), float(p["Pi_amp"])
    return lambda y: mid + amp * np.cos(2.0 * np.pi * np.asarray(y, float) / lam)


def conf_task(p, seed):
    """Molecule counts per y bin in the conformally flat geometry Pi(y)^2 eta (density Pi^d, rho = 1 at Pi = 1)."""
    d = int(p["d"])
    Pi = Pi_of(p)
    Pmin = float(p["Pi_mid"]) - float(p["Pi_amp"])
    hm = h_cut(d, Pmin ** d, float(p["tail"]))                   # the lowest density sets the cut
    Ly = float(p["periods"]) * float(p["lam"])
    Lz = float(p["L_z"])
    lo = [-hm, -2.0 * hm, -hm] + ([-hm] if d == 4 else [])
    hi = [0.0, 2.0 * hm, Ly + hm] + ([Lz + hm] if d == 4 else [])
    rng = np.random.default_rng(seed)
    X = sprinkle_poisson(Box(d, lo, hi), 1.0, rng, density=lambda P: Pi(P[:, 2]) ** d)
    X = X[X[:, 0] < 0.0]
    r = count_molecules(X, hm, [0.0] + ([0.0] if d == 4 else []), [Ly] + ([Lz] if d == 4 else []))
    edges = np.linspace(0.0, Ly, int(p["n_bins"]) + 1)
    return np.histogram(r["transverse"][:, 0], edges)[0].astype(float)


# ============================================================================================== conformal fit

def bin_integrals(p, k):
    """int_bin Pi(y)^k dA for every bin (dA = dy in d = 3, dy * L_z in d = 4), Gauss-Legendre per bin."""
    d = int(p["d"])
    Pi = Pi_of(p)
    Ly = float(p["periods"]) * float(p["lam"])
    edges = np.linspace(0.0, Ly, int(p["n_bins"]) + 1)
    xg, wg = np.polynomial.legendre.leggauss(32)
    out = []
    for a, b in zip(edges[:-1], edges[1:]):
        y = 0.5 * (a + b) + 0.5 * (b - a) * xg
        out.append(0.5 * (b - a) * float((wg * Pi(y) ** k).sum()))
    return np.asarray(out) * (float(p["L_z"]) if d == 4 else 1.0)


def fit_k(mean, se, p):
    """Weighted least squares of mean = C int Pi^k over the bins: (k, C, chi2)."""
    from scipy.optimize import minimize_scalar
    w = 1.0 / np.maximum(se, 1e-12) ** 2

    def best_C(k):
        f = bin_integrals(p, k)
        C = float((w * f * mean).sum() / (w * f * f).sum())
        return C, float((w * (mean - C * f) ** 2).sum())

    res = minimize_scalar(lambda k: best_C(k)[1], bounds=(-2.0, int(p["d"]) + 2.0), method="bounded",
                          options={"xatol": 1e-6})
    C, chi2 = best_C(float(res.x))
    return float(res.x), C, chi2


# ============================================================================================== run

def run(p, ctx):
    d, mode = int(p["d"]), str(p["mode"])
    t0 = time.time()
    a = molecule_constant(d)
    lines = []
    if mode == "flat":
        n, chunk = int(p["n_spr"]), int(p["chunk"])
        A = float(p["L"]) ** (d - 2) if d > 2 else 1.0
        rows, ok_all = [], True
        for rho in [float(r) for r in p["densities"]]:
            seeds = [derive_seed(ctx.seed, "K.14", ctx.part, rho, k) for k in range(n)]
            groups = [seeds[i:i + chunk] for i in range(0, n, chunk)]
            res = parallel_map(flat_task, [(p, g, rho) for g in groups], ctx.workers, ctx.progress,
                               lambda i, r_, dd, nn: f"[K.14] {ctx.part} rho {rho:g}: {dd}/{nn} ({time.time()-t0:.0f}s)")
            c = np.asarray([v for g in res for v in g], float)
            pred = a * rho ** ((d - 2.0) / d) * A
            m, se = float(c.mean()), float(c.std(ddof=1) / math.sqrt(c.size))
            z = (m - pred) / se if se > 0 else float("inf")
            ok = abs(z) <= float(p["z_max"])
            ok_all &= ok
            rows.append({"rho": rho, "mean": m, "se": se, "pred": pred, "z": z, "n": int(c.size),
                         "per_area": m / A, "per_area_se": se / A})
            lines.append(f"  rho {rho:g}: molecules per sprinkling {m:.4f} +- {se:.4f} (n {c.size}), prediction "
                         f"a^({d}) rho^({d - 2}/{d}) A = {pred:.4f}; per unit area {m / A:.5f} +- {se / A:.5f} "
                         f"(a^({d}) rho^((d-2)/d) = {a * rho ** ((d - 2.0) / d):.5f}); z = {z:+.2f}")
        head = (f"flat {d - 1}+1: Sigma t = 0, horizon t = x; a^({d}) = {a:.6f} (quadrature); window "
                + (f"side L = {p['L']:g}" if d > 2 else "a point (d = 2)") + f", tail e^-{p['tail']:g}")
        lines.insert(0, head)
        return Result("PASSED" if ok_all else "FAILED", lines, {"a": a, "rows": rows})

    # conformal
    n = int(p["n_spr"])
    seeds = [derive_seed(ctx.seed, "K.14", ctx.part, k) for k in range(n)]
    cnt = np.asarray(parallel_map(conf_task, [(p, sd) for sd in seeds], ctx.workers, ctx.progress,
                                  lambda i, r_, dd, nn: f"[K.14] {ctx.part}: {dd}/{nn} ({time.time()-t0:.0f}s)"), float)
    mean = cnt.mean(axis=0)
    se = cnt.std(axis=0, ddof=1) / math.sqrt(n)
    k, C, chi2 = fit_k(mean, se, p)
    rng = np.random.default_rng(derive_seed(ctx.seed, "K.14-boot", ctx.part))
    ks = []
    for _ in range(int(p["boot_B"])):
        b = cnt[rng.integers(0, n, n)]
        ks.append(fit_k(b.mean(axis=0), b.std(axis=0, ddof=1) / math.sqrt(n), p)[0])
    k_se = float(np.std(ks, ddof=1))
    target = d - 2.0
    d_alt = min(abs(target - 0.0), abs(target - d))
    tol = d_alt / 4.0
    nse = float(p["n_se"])
    if abs(k - target) + nse * k_se <= tol:
        status = "PASSED"
    elif abs(k - target) - nse * k_se > tol:
        status = "FAILED"
    else:
        status = "UNDECIDED"
    f_true = bin_integrals(p, target)
    chi2_free = float((((mean - a * f_true) / np.maximum(se, 1e-12)) ** 2).sum())
    chi2_alt = {kk: float((((mean - a * bin_integrals(p, kk)) / np.maximum(se, 1e-12)) ** 2).sum()) for kk in (0.0, float(d))}
    lines += [f"conformal {d - 1}+1: Pi(y) = {p['Pi_mid']} + {p['Pi_amp']} cos(2 pi y / {p['lam']}), {p['periods']} periods, "
              f"{p['n_bins']} bins" + (f", L_z = {p['L_z']}" if d == 4 else "") + f"; {n} sprinklings; a^({d}) = {a:.6f}",
              f"  fitted exponent k = {k:.4f} +- {k_se:.4f} (bootstrap), target d - 2 = {target:g}, alternatives 0 and {d}; "
              f"tolerance d_alt/4 = {tol:.3f}",
              f"  amplitude C / a^({d}) = {C / a:.4f} (leading-order prediction 1); chi^2 of the fit {chi2:.1f} for "
              f"{len(mean) - 2} dof",
              f"  parameter-free model (k = {target:g}, C = a): chi^2 {chi2_free:.1f} for {len(mean)} bins; with k = 0: "
              f"{chi2_alt[0.0]:.1f}; with k = {d}: {chi2_alt[float(d)]:.1f}",
              "  bins (mean count, prediction a int Pi^(d-2)): " + "; ".join(f"{m_:.2f}/{f_:.2f}" for m_, f_ in
                                                                             zip(mean, a * f_true))]
    return Result(status, lines, {"a": a, "k": k, "k_se": k_se, "C": C, "chi2": chi2, "chi2_free": chi2_free,
                                  "chi2_alt": {str(k_): v for k_, v in chi2_alt.items()}, "mean": mean, "se": se,
                                  "counts": cnt})
