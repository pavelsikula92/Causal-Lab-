"""
F8.1 -- Prediction P1: the static field of a point mass from the three sectors of the BD action
==============================================================================================

Why this test exists (29 Oct 2026; paper Part VIII, P1 and the window of l_k)
----------------------------------------------------------------------------
The static field of a point mass follows from the three coefficients A, B, C of the static quadratic
form of the action of the order (omega = conformal factor, nu = optical index; EH: 6, 8, 2), WITHOUT
assuming covariance. Matter couples to omega only; eliminating nu gives

    Phi(k)/Phi_GR(k) = -2 / (A - B^2/4C),        gamma(k) = B/(2C) - 1.

F6.8a / F6.8b2 / F6.8c computed, on the same grid of k and ages T,
    A = 6 [K_loc(k) - c_T k^3 T],  B = 8 [K_loc,ab(k) - c_T,ab k^3 T],  C = 2 [K_loc,bb(k) - c_T,bb k^3 T]
(l_k = 1). At the distances where gravity is measured k l_k is tiny, so K_loc -> 1 (corrections
(k l_k)^2 ln(1/k l_k), unobservable for l_k <~ 1e-19 m), but the age term is not suppressed: with
x = k^3 l_k^2 T and a, b, c the three age coefficients, to first order in x

    Phi/Phi_GR = 1 + rho x,  rho = -3a + 8b - 4c;      gamma = 1 + 2 (c - b) x;
    light (Phi + Psi) relative to GR: 1 + (rho + c - b) x.

Fourier inversion (int d^3k/(2pi)^3 e^{ik.r} |k| = -1/(pi^2 r^4) for r != 0) turns the k^3 term into a
power law in the far field r >> L_*:

    V(r) = -G M m / r [1 + beta_4 (1 mm / r)^3],   beta_4 = -(4/pi) rho l_k^2 T / (1 mm)^3,

exactly the form that torsion-balance experiments constrain (Adelberger et al. 2007, Table I:
|beta_4| <= 4.9e-5 at 68% CL). The shape is fixed by the counting of intervals; only the product l_k^2 T is
free. The linear expansion breaks down at L_* = (rho l_k^2 T)^(1/3) (A - B^2/4C = 0, the static pole), and
what happens there is open; the far-field law does not depend on it.

What is computed (nothing is fitted to gravity data)
---------------------------------------------------
From the latest full runs of F6.8b2 and F6.8c (recorded per-k linear fits in T on the same grid):
the age coefficients a, b, c on the plateau k = 0.07-0.2 (declared before the run; their spread over the
plateau is the stated uncertainty; k = 0.05, the longest ages, and k = 0.3, where the linear model fails,
are reported separately), rho, the coefficients of the far-field laws for Phi and for light, the
long-wavelength gamma from the local parts, and, for each reading of T, the bound on l_k from |beta_4| and
the length L_*.

Instruments (INSTRUMENT_FAILED if any fails)
--------------------------------------------
  I1 inputs: full runs of F6.8b2 and F6.8c exist, their instruments passed (status not INSTRUMENT_FAILED),
     and they share the grid (same k and same ages);
  I2 the conformal K_hat computed independently inside F6.8b2 and F6.8c on the same nodes agrees row by row
     (|diff| <= 1e-8): the two runs used the same quadrature;
  I3 the Fourier identity behind the power law: (1/(2 pi^2 r)) int k^2 sin(kr) e^{-eps k} dk by quadrature
     against the closed form Im[2/(eps - i r)^3]/(2 pi^2 r) (relative <= 1e-6), whose eps -> 0 limit is
     -1/(pi^2 r^4).
Status: RECORDED (a prediction: numbers fixed before any comparison; the comparison with data gives a
bound on l_k^2 T, conditional on what T is). If F6.8c did not pass (C(0) != 2) the lines say so: then gamma != 1
at large scales and the far-field numbers are not a prediction of a viable theory.
"""

from __future__ import annotations

import math

import numpy as np

from core import provenance
from core.testing import Result

TEST = {
    "id": "F8.1",
    "phase": "F8",
    "title": "Prediction P1: far-field law of a point mass from the three sectors of the BD action (age term)",
    "theory_part": "paper Part VIII (P1, P2, window of l_k); docs/derivations/path_to_predictions.md",
    "question": "What deviation from Newton's law and from gamma = 1 does the action of the order predict for the "
                "static field of a point mass, with the coefficients of all three sectors computed and nothing "
                "fitted to gravity data, and what does the torsion-balance bound on power-law potentials imply "
                "for l_k for each reading of the age T?",
    "method": "per-k linear fits in T recorded by F6.8b2 (conformal and mixed) and F6.8c (optical) on the same grid; "
              "age coefficients a, b, c on the plateau k = 0.07-0.2; Phi/Phi_GR = -2/(A - B^2/4C) and gamma = B/2C - 1 "
              "to first order in x = k^3 l_k^2 T; Fourier inversion of the k^3 term to a 1/r^4 potential; comparison "
              "with |beta_4| <= 4.9e-5 (Adelberger et al. 2007, 68% CL)",
    "prediction": "V(r) = -GMm/r [1 + beta_4 (1 mm/r)^3] with beta_4 = -(4/pi) rho l_k^2 T/(1 mm)^3, "
                  "rho = -3a + 8b - 4c from the counting of intervals",
    "criterion": "INSTRUMENT: I1 inputs are full runs with passed instruments on the same grid; I2 conformal K_hat of "
                 "F6.8b2 and F6.8c agree row by row (<= 1e-8); I3 Fourier identity by quadrature (<= 1e-6 rel.). "
                 "Status RECORDED (prediction numbers and the conditional bounds on l_k).",
    "null_model": "Einstein-Hilbert (a = b = c = 0): no deviation, beta_4 = 0",
    "level": "prediction",
    "meters": [],
    "depends_on": ["F6.8a", "F6.8b2", "F6.8c"],
    "params": {"plateau_ks": [0.07, 0.1, 0.14, 0.2], "beta4_bound": 4.9e-5, "r0_m": 1e-3,
               "T_readings_m": {"age of the universe (13.8 Gyr)": 1.3057e26,
                                "age of the field of the Earth (4.5 Gyr)": 4.258e25,
                                "one year": 9.461e15, "one day": 2.590e13},
               "l_k_examples_m": [1.616e-35, 1e-20, 1e-19],
               "tol_conf": 1e-8, "tol_ft": 1e-6},
    "parts": {"main": {}},
    "quick": {"tol_ft": 1e-6},
    "estimate": "seconds (reads recorded results)",
}


def _latest_full(tid):
    """Newest full (non-quick) record of part main, or None."""
    recs = [r for r in provenance.load_all(tid) if r.get("part") == "main" and not r.get("quick")]
    return recs[-1] if recs else None


def check_fourier(r=1.7, eps=0.05):
    """I3: (1/(2 pi^2 r)) int_0^inf k^2 sin(kr) e^{-eps k} dk by quadrature vs Im[2/(eps - i r)^3]/(2 pi^2 r)."""
    kmax = 60.0 / eps
    n = 400000
    k = np.linspace(0.0, kmax, n + 1)
    f = k * k * np.sin(k * r) * np.exp(-eps * k)
    w = np.full(n + 1, 2.0)
    w[1:-1:2] = 4.0
    w[0] = w[-1] = 1.0
    num = float((w * f).sum() * (kmax / n) / 3.0) / (2.0 * math.pi ** 2 * r)
    ana = float(np.imag(2.0 / (eps - 1j * r) ** 3)) / (2.0 * math.pi ** 2 * r)
    lim = -1.0 / (math.pi ** 2 * r ** 4)
    return {"num": num, "ana": ana, "lim": lim, "rel": abs(num - ana) / abs(ana)}


def run(p, ctx):
    lines = ["Static field of a point mass from the three sectors (A, B, C) of the expected BD action; far field r >> L_*; "
             "units c = 1, l_k = 1 in the sectors"]
    rb, rc = _latest_full("F6.8b2"), _latest_full("F6.8c")
    ok1 = rb is not None and rc is not None
    if ok1:
        ok1 = rb["status"] != "INSTRUMENT_FAILED" and rc["status"] != "INSTRUMENT_FAILED"
    grid_ok = False
    if rb is not None and rc is not None:
        gb = sorted((round(r_["k"], 6), round(r_["T"], 4)) for r_ in rb["data"]["rows"])
        gc = sorted((round(r_["k"], 6), round(r_["T"], 4)) for r_ in rc["data"]["rows"])
        grid_ok = gb == gc
    ok1 = ok1 and grid_ok
    lines.append(f"  I1 inputs: F6.8b2 {rb['status'] if rb else 'missing'} ({rb['_file'] if rb else '-'}), "
                 f"F6.8c {rc['status'] if rc else 'missing'} ({rc['_file'] if rc else '-'}); same grid: {grid_ok}")
    if rb is None or rc is None or not grid_ok:
        return Result("INSTRUMENT_FAILED", lines, {"inputs": [bool(rb), bool(rc), grid_ok]})

    # I2: conformal K_hat of both runs, row by row
    cb = {(round(r_["k"], 6), round(r_["T"], 4)): r_["Kc"] for r_ in rb["data"]["rows"]}
    dmax = max(abs(r_["Kc"] - cb[(round(r_["k"], 6), round(r_["T"], 4))]) for r_ in rc["data"]["rows"])
    ok2 = dmax <= float(p["tol_conf"])
    lines.append(f"  I2 conformal K_hat of F6.8b2 vs F6.8c on the same rows: max |diff| {dmax:.1e} (tol {p['tol_conf']:g})")
    ft = check_fourier()
    ok3 = ft["rel"] <= float(p["tol_ft"])
    lines.append(f"  I3 Fourier identity: quadrature {ft['num']:+.9e} vs closed form {ft['ana']:+.9e} (rel {ft['rel']:.1e}; "
                 f"eps -> 0 limit {ft['lim']:+.6e}) (tol {p['tol_ft']:g})")
    if not (ok1 and ok2 and ok3):
        return Result("INSTRUMENT_FAILED", lines, {"conf_diff": dmax, "ft": ft})

    pb, pc = rb["data"]["per_k"], rc["data"]["per_k"]
    ks = sorted(float(k) for k in pc)
    lines.append("  per k (linear fits in T, same grid): a = c_T (conformal), b = c_T,ab (mixed), c = c_T,bb (optical); "
                 "rho = -3a + 8b - 4c; gamma_loc = 2 K_loc,ab / K_loc,bb - 1")
    per = {}
    for k in ks:
        kb, kc = pb[str(k)] if str(k) in pb else pb[repr(k)], pc[str(k)]
        a, b, c = kc["cT"], kb["cT_ab"], kc["cT_bb"]
        rho = -3.0 * a + 8.0 * b - 4.0 * c
        gl = 2.0 * kb["Kloc_ab"] / kc["Kloc_bb"] - 1.0
        per[k] = {"a": a, "b": b, "c": c, "rho": rho, "gamma_loc": gl, "Kloc": kc["Kloc"], "Kloc_ab": kb["Kloc_ab"],
                  "Kloc_bb": kc["Kloc_bb"]}
        tag = "" if k in p["plateau_ks"] else "  (outside the plateau, reported only)"
        lines.append(f"    k = {k:<5g}: a {a:+.4f}  b {b:+.4f}  c {c:+.4f}  rho {rho:+.4f}  gamma_loc {gl:+.5f}{tag}")
    pl = [per[k] for k in p["plateau_ks"] if k in per]
    arr = lambda key: np.array([q[key] for q in pl])
    a_m, b_m, c_m = float(arr("a").mean()), float(arr("b").mean()), float(arr("c").mean())
    rho_v = arr("rho")
    rho_m, rho_s = float(rho_v.mean()), float(0.5 * (rho_v.max() - rho_v.min()))
    lt_v = rho_v + arr("c") - arr("b")
    lt_m, lt_s = float(lt_v.mean()), float(0.5 * (lt_v.max() - lt_v.min()))
    g_m = 2.0 * (c_m - b_m)
    lines.append(f"  plateau k = {', '.join(f'{k:g}' for k in p['plateau_ks'])}: a = {a_m:.4f}, b = {b_m:.4f}, c = {c_m:.4f}; "
                 f"rho = {rho_m:+.4f} +- {rho_s:.4f} (half spread)")
    lines.append(f"  far field (r >> L_*): Phi/Phi_N = 1 - (4/pi) rho l_k^2 T / r^3 = 1 {-4 / math.pi * rho_m:+.4f} l_k^2 T / r^3; "
                 f"light (Phi + Psi): 1 {-4 / math.pi * lt_m:+.4f} (+- {4 / math.pi * lt_s:.4f}) l_k^2 T / r^3; "
                 f"gamma(k) = 1 + {g_m:+.4f} k^3 l_k^2 T")
    r0 = float(p["r0_m"])
    bnd = float(p["beta4_bound"])
    lines.append(f"  torsion balance, |beta_4| <= {bnd:g} (r0 = 1 mm): l_k^2 T <= {bnd * r0 ** 3 * math.pi / (4 * abs(rho_m)):.3e} m^3")
    bounds = {}
    for name, T in p["T_readings_m"].items():
        lmax = math.sqrt(bnd * r0 ** 3 * math.pi / (4.0 * abs(rho_m) * T)) if rho_m != 0 else float("inf")
        Ls = {f"{l:g}": (abs(rho_m) * l * l * T) ** (1.0 / 3.0) for l in p["l_k_examples_m"]}
        b4 = {f"{l:g}": -(4.0 / math.pi) * rho_m * l * l * T / r0 ** 3 for l in p["l_k_examples_m"]}
        bounds[name] = {"T_m": T, "l_k_max_m": lmax, "L_star_m": Ls, "beta4": b4}
        lines.append(f"    T = {name} ({T:.3g} m): l_k <= {lmax:.2e} m; L_* = "
                     + ", ".join(f"{v:.2e} m (l_k = {l})" for l, v in Ls.items())
                     + "; beta_4 = " + ", ".join(f"{v:+.2e}" for v in b4.values()))
    fb, fc = rb["data"].get("fit") or {}, rc["data"].get("fit") or {}
    if "K0_ab" in fb and "K0_bb" in fc:
        g0 = 2.0 * fb["K0_ab"] / fc["K0_bb"] - 1.0
        lines.append(f"  gamma at long wavelength from the local fits: 2 K0_ab / K0_bb - 1 = 2 x {fb['K0_ab']:.5f} / "
                     f"{fc['K0_bb']:.5f} - 1 = {g0:+.5f} (EH 1; the fits are good to ~1e-3, far coarser than Cassini's 2.3e-5: a consistency check, not that bound)")
    if rc["status"] != "PASSED":
        lines.append(f"  NOTE: F6.8c status {rc['status']}: C(0) = 2 is not established; the far-field numbers are "
                     f"conditional on it")
    data = {"inputs": {"F6.8b2": rb["_file"], "F6.8c": rc["_file"], "F6.8c_status": rc["status"]},
            "per_k": {str(k): v for k, v in per.items()}, "plateau": {"a": a_m, "b": b_m, "c": c_m, "rho": rho_m,
            "rho_halfspread": rho_s, "light": lt_m, "light_halfspread": lt_s, "gamma_coef": g_m},
            "bounds": bounds, "ft": ft, "conf_diff": dmax}
    return Result("RECORDED", lines, data)
