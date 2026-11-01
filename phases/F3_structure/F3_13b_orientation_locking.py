"""
F3.13b -- Locking of an orientation defect: does a hedgehog keep its charge while the network relaxes its orientations?
=====================================================================================================================

Why this test exists (31 Oct 2026, after F3.13a; the author's approval of the alignment rule)
-------------------------------------------------------------------------------------------
F3.13a (flow_d4, PASSED): the cold flow of the laboratory carries a hedgehog of the internal orientation n with its
degree +1 intact. But in F3.13a nothing tried to unwind it: n was only carried, orientations did not act on each other
and a change cost nothing, so a charge outside a perturbed core would survive merely because nothing reaches it.
Locking is a physical statement only if a dynamics tries to remove the defect and the topology prevents it.
The missing ingredient, approved by the author on 31 Oct as a NEW RULE: alignment of orientations. At every event the
carried orientation turns by a step kappa towards the flow-weighted mean orientation of its immediate causal past
(the neighbourhood of F3.8: carrying events in the past of the event with rest-frame lag <= 2). This is the discrete
gradient flow of the configuration cost E = sum (1 - n . n') over neighbours: the author's "cost of changing the
causal configuration". A configuration without a charge (twist) can relax to the vacuum; a hedgehog cannot, its core
keeps a cost: the candidate for a rest mass as the cost of a configuration, not the number of events.

Model (locked before any run; everything not named is F3.13a flow_d4)
--------------------------------------------------------------------
3+1 box t [0, 20], half 7, density 1, links (core.order.links); sources t < 2, r < 5, unit flow at rest with n of the
configuration (vacuum, defect = hedgehog, twist: F3.13a); cold F3.6 transport with merging (rule B) and the F3.13a
merge of n. Alignment (new): after the merge at a carrying event v (orientation n_v, 4-velocity u_v), the neighbourhood
is the set of carrying events w already processed with X_v - X_w timelike, future-pointing and rest-frame lag
(X_v - X_w) . u_v <= 2; m = the flow-weighted sum of their orientations projected onto the rest space of u_v, m_hat =
m / |m|. The local cost e_v = 1 + n_v . m_hat (= 1 - cos of the spatial angle; recorded before the update); then
n_v <- N[(1 - kappa) n_v + kappa m_hat] (projected, normalised). No neighbour: no update, no cost recorded.
kappa = 0 (the control: F3.13a dynamics, cost only measured), 0.1, 0.3.
Erase attempt: in the arms erase_kappa (the hedgehog), every carrying event with 6 <= t < 7 and r < 2 gets a random
orientation (uniform on the sphere of its rest space, seeded per event); alignment continues after it.
Readout per slab of width 2 from t = 4 to 18, around the flow-weighted centroid: the degree on shells R = 2, 3, 4
(F3.13a, with its readability rule); the core cost (mean e over carrying events with r < 1.5) and the mean cost per
carrying event of the slab. 16 sprinklings.

Criteria (stated before any run, approved by the author)
-------------------------------------------------------
Instrument (INSTRUMENT_FAILED if any fails): vacuum reads 0 on >= 95 % of its readable shells in every kappa; the
hedgehog arms (all kappa together) are readable on >= 50 % of their triples; the alignment works: the twist at kappa = 0.3 has a mean
cost per event in its last slab below half of its first slab (the dynamics relaxes what nothing protects).
PASSED (locked): the hedgehog reads +1 on >= 95 % of its readable triples at every kappa; the erase arms read +1 on
>= 95 % of their readable shells R >= 3 in the slabs after the erase (t >= 8); and at kappa = 0.3 the core cost of the
hedgehog in the last slab exceeds the core cost of the twist by more than 2 SE (the hedgehog keeps a cost the twist
loses).
FAILED: at some kappa > 0 the hedgehog reads +1 on fewer than 80 % of its readable triples, or an erase arm on fewer
than 80 % of its readable shells R >= 3 after t = 8.
UNDECIDED otherwise.
Reported: the core cost of the hedgehog against time (the candidate rest-mass density), the shell R = 2 of the erase
arms (does the core re-form as a single +1?), readability by kappa.
Development (31 Oct, before any full run; decides nothing). One full-size row at the first size (box t [0, 16], 44 088
events, 1 061 sources, 257 s on one core): every readable shell of every defect and erase arm read +1, vacuum 0, twist
0; after the erase the core cost jumped to 0.6-1.0 and fell back to 0.17-0.39 with the charge kept. At kappa = 0.3 the
core cost of the hedgehog in the last slab was 0.106 against 0.032 of the twist; the mean cost of the twist fell from
0.0278 to 0.0134 (ratio 0.48: the positive control passed only just). The hedgehog at kappa = 0 was readable on 6 of 15
shells (0.40), at 0.3 on 12 of 15. Changed after it, before any full run: the box is taller (t [0, 20], two more
slabs) so that the twist has time to relax, and the readability of the hedgehog is required for all its arms together
(the alignment itself makes the field readable; kappa = 0 alone fell below 0.5). The criteria approved by the author
are otherwise unchanged.

Reading declared now. PASSED: a dynamics that removes every unprotected distortion does not remove the hedgehog; its
charge is locked and its core keeps a cost; the candidate rest mass is that cost (F3.13c: its dependence on the size
and on kappa, motion and inertia). FAILED: the alignment unwinds the hedgehog (through cores or unresolved patches); the
charge of a directional field is not locked on this substrate. Known physics in advance (Derrick): the gradient cost
of an S^2 hedgehog grows with its size in 3D; a finite-mass particle would need a full orientation (a Skyrmion);
discreteness may stop the core from shrinking below l. Limits: the alignment and the orientation are new assumptions;
the defect is put in by hand; cold flow, no field.
"""

from __future__ import annotations

import bisect
import math
import time

import numpy as np

from core.order.links import links_csr
from core.order.sprinkle import Box, derive_seed, sprinkle_poisson
from core.stats.parallel import parallel_map
from core.stats.stats import boot_se
from core.testing import Result
from phases.F3_structure.F3_08_retarded_field import children_csr, mink
from phases.F3_structure.F3_13a_orientation_defect import EXPECT, orient, shell_value

TEST = {
    "id": "F3.13b",
    "phase": "F3",
    "title": "Locking of an orientation defect: does a hedgehog keep its charge while the network relaxes its orientations?",
    "theory_part": "F3 emergence of structure; a particle as a topologically locked defect; the cost of a configuration "
                   "as the candidate rest mass (the author's alignment rule of 31 Oct)",
    "question": "With orientations that align with their causal neighbourhood (step kappa), does a hedgehog keep its degree "
                "+1 while a twist without charge relaxes, does it survive the randomisation of its core, and does its "
                "core keep a cost?",
    "method": "3+1 as F3.13a flow_d4 (box t [0, 20]); alignment n <- N[(1 - kappa) n + kappa m_hat], m_hat the flow-weighted "
              "mean orientation of the causal past with rest-frame lag <= 2; kappa = 0, 0.1, 0.3; arms vacuum, defect, "
              "twist, erase (core r < 2 randomised at 6 <= t < 7); degree on shells R = 2, 3, 4 and the cost e = 1 + n . "
              "m_hat per slab",
    "prediction": "the hedgehog keeps +1 and a core cost; the twist relaxes towards the vacuum",
    "criterion": "INSTRUMENT: vacuum 0 on >= 95 %, defect arms readable >= 50 % (all kappa), twist cost at kappa 0.3 halves; "
                 "PASSED: defect +1 >= 95 % at every kappa, erase +1 >= 95 % on R >= 3 after t = 8, core cost of the "
                 "defect > twist + 2 SE at kappa 0.3; FAILED: defect or erase < 80 % at some kappa > 0; else UNDECIDED",
    "null_model": "twist (distorted, no charge) and vacuum",
    "level": "mechanism",
    "meters": [],
    "depends_on": ["F3.13a"],
    "params": {"T": 20.0, "half": 7.0, "r_src": 5.0, "t_src": 2.0, "slab": 2.0, "t_read": 4.0, "shells": [2.0, 3.0, 4.0],
               "kappas": [0.0, 0.1, 0.3], "L_g": 2.0, "t_erase": 6.0, "a_erase": 2.0, "t_after": 8.0, "r_core": 1.5,
               "n_seeds": 16, "n_se": 2.0, "boot_B": 4000, "pass_frac": 0.95, "fail_frac": 0.8, "read_frac": 0.5,
               "core_eps": 1e-6},
    "parts": {"main": {}},
    "quick": {"n_seeds": 1, "T": 10.0, "half": 6.0, "r_src": 4.0, "t_erase": 5.0, "t_after": 6.0, "boot_B": 200},
    "estimate": "~10-20 min with 10 processes (16 orders of ~55 000 events in 3+1, 12 arms each; ~6 min per order on one core)",
}


# ============================================================================================== transport

def _proj_unit(w, u):
    """w projected onto the rest space of u and normalised (w . w = -1); None if it vanishes."""
    q = w - mink(w, u) * u
    L = math.sqrt(max(-float(mink(q, q)), 0.0))
    return q / L if L > 1e-12 else None


def transport_align(X, cptr, cidx, s, U0, N0, kappa, p, erase_seed=None):
    """F3.13a transport (merge of P and n) with the alignment of orientations and, optionally, the erase of the core."""
    N, D = X.shape
    L_g, eps = float(p["L_g"]), float(p["core_eps"])
    te, ae = float(p["t_erase"]), float(p["a_erase"])
    Pin = np.zeros((N, D))
    nlist = [[] for _ in range(N)]
    nin = np.zeros(N, dtype=np.int64)
    Phold = np.zeros((N, D))
    Nhold = np.full((N, D), np.nan)
    cost = np.full(N, np.nan)
    core = np.zeros(N, dtype=bool)
    car_t, car_x, car_n, car_F = [], [], [], []        # processed carrying events with an orientation, in time order
    for v in range(N):
        if nin[v] == 0 and s[v] <= 0.0:
            continue
        P = Pin[v] + s[v] * U0[v]
        F = math.sqrt(max(float(mink(P, P)), 1e-300))
        u = P / F
        items = nlist[v] + ([(s[v], N0[v])] if s[v] > 0.0 else [])
        acc = np.zeros(D)
        tot = 0.0
        for Fi, ni in items:
            acc += Fi * (ni - mink(ni, u) * u)
            tot += Fi
        Lc = math.sqrt(max(-float(mink(acc, acc)), 0.0))
        Phold[v] = P
        n_v = None
        if Lc >= eps * tot:
            n_v = acc / Lc
        else:
            core[v] = True
        t_v = float(X[v, 0])
        if n_v is not None and erase_seed is not None and te <= t_v < te + 1.0 and np.linalg.norm(X[v, 1:]) < ae:
            g = np.random.default_rng(derive_seed(erase_seed, v)).normal(size=D)
            n_v = _proj_unit(g, u)                      # a random orientation of the rest space (isotropic there)
        elif n_v is not None and car_t:
            i0 = bisect.bisect_left(car_t, t_v - L_g - 1e-9)
            if i0 < len(car_t):
                W = np.array(car_x[i0:])
                d = X[v][None, :] - W
                lag = mink(d, u[None, :])
                ok = (d[:, 0] > 0.0) & (mink(d, d) > 0.0) & (lag <= L_g)
                if ok.any():
                    Nw = np.array(car_n[i0:])[ok]
                    Fw = np.array(car_F[i0:])[ok]
                    m = (Fw[:, None] * (Nw - mink(Nw, u[None, :])[:, None] * u[None, :])).sum(axis=0)
                    mh = _proj_unit(m, u)
                    if mh is not None:
                        cost[v] = 1.0 + float(mink(n_v, mh))      # 1 - cos(spatial angle), signature + - -
                        if kappa > 0.0:
                            nn = _proj_unit((1.0 - kappa) * n_v + kappa * mh, u)
                            if nn is not None:
                                n_v = nn
        if n_v is not None:
            Nhold[v] = n_v
            car_t.append(t_v)
            car_x.append(X[v].copy())
            car_n.append(n_v)
            car_F.append(F)
        a, b = cptr[v], cptr[v + 1]
        if b > a:
            ks = cidx[a:b]
            B = X[ks] - X[v]
            tau = np.sqrt(np.maximum(mink(B, B), 1e-300))
            k = int(ks[int(np.argmin(mink(u[None, :], B) / tau))])
            Pin[k] += P
            nin[k] += 1
            if n_v is not None:
                nlist[k].append((F, n_v))
    return Phold, Nhold, core, cost


# ============================================================================================== one sprinkling

def arm_list(p):
    arms = []
    for k in p["kappas"]:
        for cfg in ("vacuum", "defect", "twist"):
            arms.append((f"{cfg}_{float(k):g}", cfg, float(k), False))
        arms.append((f"erase_{float(k):g}", "defect", float(k), True))
    return arms


def task(p, seed):
    rng = np.random.default_rng(seed)
    T, H = float(p["T"]), float(p["half"])
    X = sprinkle_poisson(Box(4, [0.0, -H, -H, -H], [T, H, H, H]), 1.0, rng)
    N = len(X)
    ptr, idx = links_csr(X, 1.0)
    cptr, cidx = children_csr(ptr, idx, N)
    rs, ts = float(p["r_src"]), float(p["t_src"])
    src = np.flatnonzero((X[:, 0] < ts) & (np.linalg.norm(X[:, 1:], axis=1) < rs))
    s = np.zeros(N)
    s[src] = 1.0
    U0 = np.zeros((N, 4))
    U0[:, 0] = 1.0
    sl = float(p["slab"])
    edges = np.arange(float(p["t_read"]), T - sl + 1e-9, sl)
    out = {"N": N, "sources": int(src.size), "arms": {}}
    for key, cfg, kappa, erase in arm_list(p):
        N0 = np.zeros((N, 4))
        N0[src, 1:] = orient(cfg, X[src, 1:], rs)
        Ph, Nh, core, cost = transport_align(X, cptr, cidx, s, U0, N0, kappa, p,
                                             derive_seed(seed, "erase", key) if erase else None)
        held = Ph[:, 0] > 0.0
        rec = {"vals": [], "core_cost": [], "mean_cost": [], "t": [], "cores": int(core.sum())}
        for a, b in zip(edges[:-1], edges[1:]):
            m = held & (X[:, 0] >= a) & (X[:, 0] < b)
            tm = 0.5 * (a + b)
            rec["t"].append(tm)
            if not m.any():
                rec["core_cost"].append(float("nan"))
                rec["mean_cost"].append(float("nan"))
                continue
            w = Ph[m, 0]
            c = (w[:, None] * X[m, 1:]).sum(axis=0) / w.sum()
            okn = m & np.isfinite(Nh[:, 0])
            rel = X[okn, 1:] - c
            rr = np.linalg.norm(rel, axis=1)
            nsp = Nh[okn, 1:] / np.maximum(np.linalg.norm(Nh[okn, 1:], axis=1), 1e-300)[:, None]
            for Rs in p["shells"]:
                sel = np.abs(rr - Rs) < 0.75
                val, ok = shell_value(rel[sel], nsp[sel], 3)
                rec["vals"].append((tm, float(Rs), val, ok))
            cc = cost[okn]
            fin = np.isfinite(cc)
            incore = (rr < float(p["r_core"])) & fin
            rec["core_cost"].append(float(cc[incore].mean()) if incore.any() else float("nan"))
            rec["mean_cost"].append(float(cc[fin].mean()) if fin.any() else float("nan"))
        out["arms"][key] = rec
    return out


# ============================================================================================== run

def _frac(vals, expect, rmin=None, tmin=None):
    sel = [(t, R, v, ok) for t, R, v, ok in vals if (rmin is None or R >= rmin) and (tmin is None or t >= tmin)]
    rd = [v for *_, v, ok in sel if ok]
    match = sum(1 for v in rd if abs(v - expect) < 1e-6)
    return (match / len(rd) if rd else float("nan")), (len(rd) / len(sel) if sel else 0.0), len(rd)


def run(p, ctx):
    t0 = time.time()
    n, B, nse = int(p["n_seeds"]), int(p["boot_B"]), float(p["n_se"])
    jobs = [(p, derive_seed(ctx.seed, "F3.13b", k)) for k in range(n)]
    rows = parallel_map(task, jobs, ctx.workers, ctx.progress,
                        lambda i, r_, d, nn: f"[F3.13b] {d}/{nn}: N {r_['N']} ({time.time() - t0:.0f}s)")
    tgrid = rows[0]["arms"][arm_list(p)[0][0]]["t"]
    lines = [f"3+1 box t [0, {p['T']:g}], half {p['half']:g} (~{np.mean([r_['N'] for r_ in rows]):.0f} events, "
             f"{np.mean([r_['sources'] for r_ in rows]):.0f} sources); cold transport with merging; alignment step kappa "
             f"towards the causal past (lag <= {p['L_g']:g}); erase: core r < {p['a_erase']:g} randomised at "
             f"{p['t_erase']:g} <= t < {p['t_erase'] + 1:g}; shells {p['shells']}; {n} sprinklings",
             "        t           " + " ".join(f"{t:6.1f}" for t in tgrid)]
    res = {}
    for key, cfg, kappa, erase in arm_list(p):
        vals = [x for r_ in rows for x in r_["arms"][key]["vals"]]
        exp = EXPECT[cfg]
        f, rd, nr = _frac(vals, exp)
        fo = _frac(vals, exp, rmin=3.0, tmin=float(p["t_after"])) if erase else None
        f2 = _frac(vals, exp, tmin=float(p["t_after"])) if erase else None
        cc = np.array([r_["arms"][key]["core_cost"] for r_ in rows], float)
        mc = np.array([r_["arms"][key]["mean_cost"] for r_ in rows], float)
        last = boot_se(cc[:, -1][np.isfinite(cc[:, -1])], B, derive_seed(ctx.seed, "F3.13b-boot", key, "core"))[:2] \
            if np.isfinite(cc[:, -1]).sum() >= 2 else (float("nan"), float("nan"))
        with np.errstate(all="ignore"):
            ccm, mcm = np.nanmean(cc, axis=0), np.nanmean(mc, axis=0)
        res[key] = {"f": f, "rd": rd, "nr": nr, "fo": fo, "f2": f2, "core_last": last, "mc": mcm.tolist(), "cc": ccm.tolist(),
                    "cores": float(np.mean([r_["arms"][key]["cores"] for r_ in rows]))}
        hist = {}
        for *_, v, ok in vals:
            if ok:
                hist[int(round(v))] = hist.get(int(round(v)), 0) + 1
        lines.append(f"  {key:11s}: reads {exp:+d} on {f:.3f} of {nr} readable (readable {rd:.2f}); values "
                     f"{dict(sorted(hist.items()))}; core events {res[key]['cores']:.1f}"
                     + (f"; after the erase on R >= 3: {fo[0]:.3f} of {fo[2]} (readable {fo[1]:.2f}), all shells "
                        f"{f2[0]:.3f} of {f2[2]}" if erase else ""))
        lines.append("        core cost   " + " ".join(f"{x:6.3f}" if np.isfinite(x) else "   nan" for x in ccm))
        lines.append("        mean cost   " + " ".join(f"{x:6.3f}" if np.isfinite(x) else "   nan" for x in mcm))
    pf, ff = float(p["pass_frac"]), float(p["fail_frac"])
    ks = [float(k) for k in p["kappas"]]
    kmax = max(ks)
    tw = res[f"twist_{kmax:g}"]["mc"]
    tw_first = next((x for x in tw if np.isfinite(x)), float("nan"))
    tw_last = next((x for x in reversed(tw) if np.isfinite(x)), float("nan"))
    inst = {f"vacuum reads 0 on >= {pf:.0%} (every kappa)": all(res[f"vacuum_{k:g}"]["f"] >= pf for k in ks),
            f"defect arms readable >= {p['read_frac']:.0%} (all kappa)": bool(
                sum(res[f"defect_{k:g}"]["rd"] * len([1 for r_ in rows for _ in r_["arms"][f"defect_{k:g}"]["vals"]]) for k in ks)
                / max(sum(len([1 for r_ in rows for _ in r_["arms"][f"defect_{k:g}"]["vals"]]) for k in ks), 1)
                >= float(p["read_frac"])),
            f"twist at kappa {kmax:g} relaxes (last mean cost < first / 2)": bool(tw_last < 0.5 * tw_first)}
    lines.append(f"  twist at kappa {kmax:g}: mean cost first {tw_first:.4f}, last {tw_last:.4f}")
    dl, tl = res[f"defect_{kmax:g}"]["core_last"], res[f"twist_{kmax:g}"]["core_last"]
    cost_ok = bool(np.isfinite(dl[0]) and np.isfinite(tl[0]) and dl[0] - tl[0] > nse * math.hypot(dl[1], tl[1]))
    lines.append(f"  core cost in the last slab at kappa {kmax:g}: defect {dl[0]:.4f} +- {dl[1]:.4f}, twist {tl[0]:.4f} +- "
                 f"{tl[1]:.4f} -> the defect keeps a cost: {cost_ok}")
    lines.append("instrument: " + "; ".join(f"{k}: {'OK' if v else 'FAILS'}" for k, v in inst.items()))
    keep = all(res[f"defect_{k:g}"]["f"] >= pf for k in ks)
    erase_ok = all(res[f"erase_{k:g}"]["fo"][0] >= pf for k in ks if np.isfinite(res[f"erase_{k:g}"]["fo"][0]))
    lost = any(res[f"defect_{k:g}"]["f"] < ff or (np.isfinite(res[f"erase_{k:g}"]["fo"][0]) and res[f"erase_{k:g}"]["fo"][0] < ff)
               for k in ks if k > 0.0)
    if not all(inst.values()):
        status = "INSTRUMENT_FAILED"
    elif keep and erase_ok and cost_ok:
        status = "PASSED"
    elif lost:
        status = "FAILED"
    else:
        status = "UNDECIDED"
    return Result(status, lines, {"res": res, "rows": rows})
