"""
F5.2d -- Active gravitational mass from the causal growth of a 2+1 order whose density and cones follow the potential
         of matter (direction B, B3a: the active side)
=====================================================================================================================

Why this test exists (13 Nov 2026; card docs/B3_card_draft.md, revision 2, frozen by decision 23)
-----------------------------------------------------------------------------------------------
B1 (F3.15) and B2 (F5.2c) are PASSED: on a 2+1 order the field rule that reads only the order gives inertia, and a
prescribed curved order of the theory moves walls as the continuum metric does. B3 asks the active question: when the
order GROWS, event by event, with its density and its cones set by a potential that matter produces from its realised
causal past, does the order become the order of the theory around the source, as much as the source's load says?
The author chose (13 Nov) B3a: the cones are cut in the rest frame of the static source (the coordinate frame, a
declared proxy as F7.1b), and only the active side is tested; the fall of a test body moves to B4 (walls do not fall in
the order of the theory, B2; a breather's fall is swamped by its redshifted phase).

Model (frozen; diag/b3_growth.py, md5 in decision 23)
----------------------------------------------------
Candidates: Poisson in the pyramid |x| <= T + 6 - t, |y| <= T + 2.5 - t, T 9, density rho e^{-Lmin}, independent marks.
In time order, every candidate v
  1. finds its parents among the REALISED earlier events by the cone rule of B3a: u < v iff |dx|^2 + |dy|^2 <
     (Pi_m^2 dt)^2, Pi_m^2 the square of Pi at the coordinate midpoint of the pair from the binned potential of the
     realised events (bins dt 0.25, dx 0.5, dy 2); links = the maximal ones in the shell tau^2 <= tau_cut^2;
  2. gets its longest chains from every earlier realised event (backward dynamic programme over the links);
  3. phi(v) = phi0(v) - (1/rho) sum w(L_uv) mu sin phi(u)  (the chain weight of B0/B1, variant B, frozen constants);
     L(v) = -(kappa/rho) sum w(L_uv) q(u) - (kappa/2) s(t_v) f(v), f = phi (2 pi - phi)/2,
     q = s [3V + (pi - phi) mu sin(phi)/2] - s'' f/2   (the load = trace of T by the identity of decision 19 with a
     quadratic vanishing on both vacua; gravity switched on by the smoothstep s from t 3 to 4);
  4. is realised iff its mark < exp(Lmin - L(v)): coordinate density rho Pi^-1 (the metric of the theory
     -Pi^2 dt^2 + Pi^-2 (dx^2 + dy^2)), with Pi = e^L.
phi0 is the flat free wave of the kink data (exact: gravity starts after the fronts of the free wave have left the
wall). Source: one kink wall at rest at x = 0, mu 1 (E 8, load 2E = 16) and mu 0.25 (E 4, load 8); kappa 0.01. Every
sprinkling is run at kappa and at kappa = 0 on the same candidates and marks.
Continuum reference: box_g phi = mu sin phi, box_g L = kappa s q with the full trace (diag/b3_continuum, exact).

Readouts (frozen; coordinates declared)
--------------------------------------
  r_L    = sum L L_c / sum L_c^2 over the realised events of the band |y| < 2, |x| < 6, t in [5, 8) (the potential the
           order computed, projected on the continuum potential);
  u      = r_L(mu 1) - r_L(mu 0.25), paired by sprinkling (the potential in proportion to the load: active mass);
  clocks = [(near/far) of the longest-chain length]_kappa / [same]_0, near = |x| < 0.5, far = |x - 4.5| < 0.5, from 12
           realised events in |y| < 1, t in [4.5, 5) each, to the events within 0.4 of straight above at t + [2.5, 3);
           against the continuum ratio of the proper times int Pi dt;
  cones  = the same ratio for the largest |dy| reached in the causal future (reachability through the links) at
           t + [2.5, 3), |dx| < 0.4; against the continuum ratio of Pi^2 at the midpoint;
  density = realised / sum of the acceptance probabilities (t > 4).
Reported: the kink drift and width, L by x in the last slab, r_L per slab.

Criteria (card revision 2, section 4; n = n_dec sprinklings; SE = sd/sqrt(n); every condition on its own)
INSTRUMENT_FAILED if any of: A0 a change of the marks of the latest 10 % of candidates changes any earlier event
  (2 sprinklings); overflow (L < Lmin) on any event; |density - 1| + 2 SE >= 0.01 (each mu); |cones/cont - 1| + 2 SE
  >= 0.03 (each mu; the rule realises the cones it was given).
FAILED if any of: |r_L - 1| - 2 SE > tol_L (mu 1: 0.06; mu 0.25: 0.12); |u| - 2 SE > 0.12; |clocks/cont - 1| - 2 SE
  > 0.06 (each mu).
PASSED if all: |r_L - 1| + 2 SE < tol_L (each mu); |u| + 2 SE < 0.12; |clocks/cont - 1| + 2 SE < 0.06 (each mu);
  (1 - clocks) - 2 SE > 0 (each mu: the clocks near the source are slowed).
Otherwise UNDECIDED.

    python run.py F5.2d --workers N          (seed 0 decides)
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
from diag import b3_growth as G

TEST = {
    "id": "F5.2d",
    "phase": "F5",
    "title": "Active gravitational mass from the causal growth of a 2+1 order (density and cones from the potential, B3a)",
    "theory_part": "F5 gravity; direction B (2+1), B3a (the author's choice of 13 Nov): the order grows with density "
                   "Pi^-1 and cones Pi^2 from a potential that matter produces from its realised past (decisions 15, 16, "
                   "19, 23); the Green weight from the order (decision 17, variant B)",
    "question": "When a 2+1 order grows event by event with density and cones set by the potential L = ln Pi that a static "
                "source wall produces from its realised past through the chain weight, does the order become the order of "
                "the theory around the source, as much as the source's load (2E) says: the potential, the slowed clocks "
                "and the narrowed cones of the continuum response box_g L = kappa q?",
    "method": "growth of candidates in time order (rho 30, kappa 0.01, T 9); cones from the binned potential at the "
              "midpoint (B3a proxy); field and potential by the chain weight; load = trace of T by the identity of "
              "decision 19; kink source mu 1 and 0.25; each paired with kappa = 0; continuum references",
    "prediction": "r_L = 1 for both sources; u = 0; clocks near/far as int Pi dt (0.878 and 0.954); cones as Pi^2",
    "criterion": "INSTRUMENT: A0 exact, no overflow, density 0.01, cones 0.03 (+ 2 SE). PASSED: |r_L - 1| + 2 SE < 0.06 "
                 "(mu 1) and 0.12 (mu 0.25); |u| + 2 SE < 0.12; |clocks/cont - 1| + 2 SE < 0.06 and (1 - clocks) - 2 SE > 0 "
                 "(each mu). FAILED: any |.| - 2 SE above its tolerance. Otherwise UNDECIDED (card docs/B3_card_draft.md)",
    "null_model": "kappa = 0 on the same candidates and marks (a flat order of density rho); continuum equations",
    "level": "mechanism",
    "meters": [],
    "depends_on": ["F3.15", "F5.2c"],
    "params": {"rho": 30.0, "kappa": 0.01, "mus": [1.0, 0.25], "n_dec": 32, "n_a0": 2, "n_probe": 12,
               "growth": {"T": 9.0, "dX": 6.0, "dY": 2.5, "t1": 3.0, "t2": 4.0, "Lmin": -0.7, "band": 2.0,
                          "bt": 0.25, "bx": 0.5, "by": 2.0, "mid": True},
               "calib": [0.2732723632276477, 0.21313211596642687, 1.2025922117299621],
               "tol_dens": 0.01, "tol_cones": 0.03, "tol_L": [0.06, 0.12], "tol_u": 0.12, "tol_clock": 0.06,
               "n_se": 2.0},
    "quick": {"rho": 12.0, "n_dec": 2, "n_a0": 1, "n_probe": 6},
    "estimate": "~90 000 realised events (175 000 candidates) per run at rho 30, ~200 s per run on one core of the "
                "container; 4 runs per sprinkling (2 sources x kappa, 0) + A0; ~7 core-hours for n_dec 32. numba "
                "required. Finished tasks are cached in diag/out/F5_2d_cache (keyed by this file's hash).",
}


def task_run(p, mu, kappa, k, base):
    P = dict(p["growth"])
    if tuple(p["calib"]) != G.CALIB:
        raise SystemExit("F5.2d: calibration differs from diag/b3_growth.CALIB")
    out = G.task(float(p["rho"]), float(mu), float(kappa), int(k), base, P, label="F5.2d", n_probe=int(p["n_probe"]))
    out["kind"] = "run"
    return out


def task_a0(p, mu, k, base):
    """A0: the marks of the latest 10 % of candidates are redrawn; every earlier event must be unchanged."""
    P = dict(p["growth"])
    rho = float(p["rho"])
    T = P["T"]
    X, Y = T + P["dX"], T + P["dY"]
    sd = derive_seed(base, "F5.2d", rho, k)
    t, x, y = G.sprinkle_pyramid(rho * math.exp(-P["Lmin"]), X, Y, T, sd)
    mark = np.random.default_rng(derive_seed(sd, "marks")).uniform(0, 1, t.size)
    f = lambda z: 4.0 * np.arctan(np.exp(np.clip(math.sqrt(mu) * z, -50, 50)))
    phi0 = 0.5 * (f(x - t) + f(x + t))
    wL = B0.weight_table("B", rho, tuple(p["calib"]))
    tc = 1.15 * (12.0 * 25.0 / (math.pi * rho)) ** (1.0 / 3.0)
    args = (phi0, rho, float(mu), float(p["kappa"]), wL, tc * tc, P["Lmin"], P["t1"], P["t2"], X, Y, P["bt"], P["bx"],
            P["by"], P["mid"])
    r1 = G.grow(t, x, y, mark, *args)
    m2 = mark.copy()
    cut = int(0.9 * t.size)
    m2[cut:] = np.random.default_rng(derive_seed(sd, "marks2")).uniform(0, 1, t.size - cut)
    r2 = G.grow(t, x, y, m2, *args)
    same_real = bool(np.array_equal(r1[0][:cut], r2[0][:cut]))
    dphi = float(np.max(np.abs(r1[1][:cut] - r2[1][:cut])))
    dL = float(np.max(np.abs(r1[2][:cut] - r2[2][:cut])))
    changed = int(np.sum(r1[0][cut:] != r2[0][cut:]))
    return {"kind": "a0", "mu": mu, "k": k, "same_real": same_real, "dphi": dphi, "dL": dL, "changed_late": changed}


_SHA = None


def _sha():
    global _SHA
    if _SHA is None:
        with open(__file__, "rb") as fh:
            _SHA = hashlib.sha256(fh.read()).hexdigest()[:16]
    return _SHA


def run_task(kind, *a):
    p = a[0]
    key = json.dumps([kind] + [str(q) for q in a[1:]])
    fn = os.path.join("diag", "out", "F5_2d_cache", _sha(), hashlib.sha256(
        (key + json.dumps(p, sort_keys=True, default=str)).encode()).hexdigest()[:24] + ".json")
    if os.path.exists(fn):
        with open(fn) as fh:
            return json.load(fh)
    out = {"run": task_run, "a0": task_a0}[kind](*a)
    os.makedirs(os.path.dirname(fn), exist_ok=True)
    with open(fn + ".tmp", "w") as fh:
        json.dump(out, fh)
    os.replace(fn + ".tmp", fn)
    return out


def _ms(a):
    a = np.asarray(a, float)
    if a.size < 2 or not np.all(np.isfinite(a)):
        return float("nan"), float("nan")
    return float(a.mean()), float(a.std(ddof=1) / math.sqrt(a.size))


def _fmt(m):
    return f"{m[0]:+.4f} +- {m[1]:.4f}"


def run(p, ctx):
    if F._nb is None:
        raise SystemExit("F5.2d needs numba (pip install numba)")
    t0 = time.time()
    n, base, nse = int(p["n_dec"]), ctx.seed, float(p["n_se"])
    mus = [float(m) for m in p["mus"]]
    kap = float(p["kappa"])
    jobs = [("a0", p, mus[0], k, base) for k in range(int(p["n_a0"]))]
    for k in range(n):
        for mu in mus:
            for kp in (kap, 0.0):
                jobs.append(("run", p, mu, kp, k, base))
    res = parallel_map(run_task, jobs, ctx.workers, ctx.progress,
                       lambda i, r, d, nn: f"[F5.2d] {d}/{nn}: {r['kind']} mu {r.get('muA', r.get('mu'))} "
                                           f"kappa {r.get('kappa', '')} k {r['k']} ({time.time() - t0:.0f}s)")
    runs = [r for r in res if r["kind"] == "run"]
    a0 = [r for r in res if r["kind"] == "a0"]

    def get(mu, kp):
        d = {r["k"]: r for r in runs if r["muA"] == mu and r["kappa"] == kp}
        return [d[k] for k in range(n)]

    lines = [f"F5.2d (B3a, 2+1 growth: density Pi^-1, cones Pi^2 at the midpoint, potential by the chain weight); rho "
             f"{p['rho']:g}, kappa {kap:g}, T {p['growth']['T']:g}, {n} sprinklings x (mu {mus[0]:g}, {mus[1]:g}) x "
             f"(kappa, 0); base seed {base}; mean +- SE over sprinklings"]
    inst, passc, failc = {}, {}, {}
    a0ok = all(r["same_real"] and r["dphi"] == 0.0 and r["dL"] == 0.0 for r in a0)
    inst["A0 earlier events unchanged"] = a0ok
    lines.append("  A0: " + "; ".join(f"k {r['k']}: realisation same {r['same_real']}, max |dphi| {r['dphi']:.1e}, "
                                     f"max |dL| {r['dL']:.1e} (late candidates changed {r['changed_late']})" for r in a0))
    over = sum(r["over"] for r in runs)
    inst["no overflow"] = over == 0
    lines.append(f"  overflow (L < Lmin {p['growth']['Lmin']}): {over}; lowest L {min(r['Lmin_seen'] for r in runs):+.3f}")
    rL = {}
    late = None
    for j, mu in enumerate(mus):
        K, Z = get(mu, kap), get(mu, 0.0)
        late = [i for i, w in enumerate(K[0]["rows"]) if 5.0 <= w["t"] < 8.0]
        dens = _ms([r["dens"][0] / r["dens"][1] - 1 for r in K])
        inst[f"density mu {mu:g}"] = abs(dens[0]) + nse * dens[1] < float(p["tol_dens"])
        rL[mu] = [sum(r["rows"][i]["LL"] for i in late) / sum(r["rows"][i]["LcLc"] for i in late) for r in K]
        m = _ms(rL[mu])
        tl = float(p["tol_L"][j])
        passc[f"potential r_L mu {mu:g}"] = abs(m[0] - 1) + nse * m[1] < tl
        failc[f"potential r_L mu {mu:g}"] = abs(m[0] - 1) - nse * m[1] > tl
        probes = {}
        for what, key, ck in (("clocks", "clen", "tau_cont"), ("cones", "reach", "Pi2_mid")):
            rat, cont = [], []
            for a, z in zip(K, Z):
                pa, pz = a["probe"], z["probe"]
                rat.append((np.mean(pa["near"][key]) / np.mean(pa["far"][key]))
                           / (np.mean(pz["near"][key]) / np.mean(pz["far"][key])))
                cont.append(np.mean(pa["near"][ck]) / np.mean(pa["far"][ck]))
            c = float(np.mean(cont))
            probes[what] = (_ms(np.array(rat) / c - 1), _ms(rat), c)
        cm, craw, cc = probes["cones"]
        inst[f"cones mu {mu:g}"] = abs(cm[0]) + nse * cm[1] < float(p["tol_cones"])
        km, kraw, kc = probes["clocks"]
        passc[f"clocks mu {mu:g}"] = abs(km[0]) + nse * km[1] < float(p["tol_clock"])
        failc[f"clocks mu {mu:g}"] = abs(km[0]) - nse * km[1] > float(p["tol_clock"])
        passc[f"clocks slowed mu {mu:g}"] = (1 - kraw[0]) - nse * kraw[1] > 0
        b = np.array([[bb[:2] for bb in r["rows"][late[-1]]["bins"]] for r in K])
        drift = _ms([r["rows"][late[-1]]["xc"] - r["rows"][0]["xc"] for r in K])
        drift0 = _ms([r["rows"][late[-1]]["xc"] - r["rows"][0]["xc"] for r in Z])
        wr = _ms([r["rows"][late[-1]]["w"] / r["rows"][late[-1]]["w_ex"] - 1 for r in K])
        lines.append(f"  source mu {mu:g} (load {16 * math.sqrt(mu):g}): r_L {_fmt(m)} (per slab "
                     + ", ".join(_fmt(_ms([r['rows'][i]['LL'] / r['rows'][i]['LcLc'] for r in K])) for i in late)
                     + f"); density - 1 {_fmt(dens)}")
        lines.append(f"     clocks near/far {_fmt(kraw)} against the continuum {kc:.4f} (ratio - 1 {_fmt(km)}); cones "
                     f"{_fmt(craw)} against {cc:.4f} (ratio - 1 {_fmt(cm)})")
        lines.append("     L by x in the last read slab, order: " + " ".join(f"{v:+.3f}" for v in np.nanmean(b[:, :, 0], 0)))
        lines.append("                               continuum: " + " ".join(f"{v:+.3f}" for v in np.nanmean(b[:, :, 1], 0)))
        lines.append(f"     reported: kink drift {_fmt(drift)} (kappa 0: {_fmt(drift0)}); width w/w_ex - 1 in the last "
                     f"slab {_fmt(wr)}")
    u = _ms(np.array(rL[mus[0]]) - np.array(rL[mus[1]]))
    passc["active mass: u"] = abs(u[0]) + nse * u[1] < float(p["tol_u"])
    failc["active mass: u"] = abs(u[0]) - nse * u[1] > float(p["tol_u"])
    lines.append(f"  active mass: u = r_L(mu {mus[0]:g}) - r_L(mu {mus[1]:g}) {_fmt(u)}")
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
    lines.append(f"cost: {np.mean([r['sec'] for r in runs]):.0f} s per growth ({time.time() - t0:.0f}s)")
    return Result(status, lines, {"instrument": inst, "passed": passc, "failed": failc, "results": res})
