"""
K.17 -- Feasibility: can the quadratic form of the Benincasa-Dowker action be measured on a 3+1 sprinkling?
=========================================================================================================

Why this test exists (28 Oct 2026, docs/derivations/pravidlo_z_akce.md)
-----------------------------------------------------------------------
The substrate rule of F7.4a ("every flow narrows the cone in its own rest frame, deficits add") follows in the
linear order from the action of the flows (A6, proper time counted by the order) and from the stiffness of the order.
Lorentz invariance leaves one number, gamma = -c_PhiPsi / (2 c_PsiPsi), fixed by the static quadratic form of the
action of the order,
    L = c_PhiPhi |grad Phi|^2 + c_PhiPsi grad Phi . grad Psi + c_PsiPsi |grad Psi|^2
(ds^2 = -(1 + 2 Phi) dt^2 + (1 - 2 Psi) dx^2). The Einstein-Hilbert action, the continuum limit of the BD action,
has c_PhiPhi = 0, c_PhiPsi = -4, c_PsiPsi = 2 in 3+1, hence gamma = 1. Measuring these coefficients from the
discrete BD action of the causal set (planned F6.8) would derive the substrate rule on the substrate instead of
inserting it. The signal is second order in the field and the BD action fluctuates strongly; this test only decides
whether such a measurement is possible at the sizes the laboratory can afford.

What is measured
----------------
Static Psi-sector metric of exponential form (the K.16c family with phi = 0):
    ds^2 = -dt^2 + e^{2 omega} dx^2,   omega = a f(r),   f = exp(-r^2 / (2 w^2)),
exactly sqrt(-g) R = e^omega (-4 Lap omega - 2 |grad omega|^2) (spatially conformally flat, static).
The even part in a removes the first-order (total-derivative) term:
    Q = [S(+a) + S(-a)] / 2 - S(0),   S = -(2/rho) sum_{x in rows} (B_eps 1)(x)  ->  int sqrt(-g) R over the rows,
    Q_pred = a^2 int_band int_{r <= R_rows} f^2 (12/w^2 - 6 r^2/w^4) d^3x dt   (continuum, Einstein-Hilbert).
B_eps is the smeared 4D BD operator of the laboratory (core.operators.bd: kernel coefficients -1, 9, -16, 8,
a_4 = 4/sqrt(6)); (B_eps 1)(x) = -a_4 (eps rho)^(1/2) [1 + sum_{y < x} K_eps(n_xy)], n_xy the number of elements
strictly between y and x; E[B_eps 1] -> -R/2.

Correlated sampling: one Poisson superset of density rho_max = rho e^{3|a|} per sprinkling; every configuration
(0, +a, -a) keeps a point with the same uniform number u < sqrt(-g)(x) / rho_max (sqrt(-g) = e^{3 omega}), so the
three causal sets share most elements. Relation: y < x iff t_x - t_y >= int e^{omega} dl along the straight segment
(closed form with erf to third order in a; pairs decided by the bounds e^{+-|a| f_max} without it). n_xy by bitsets:
popcount(future(y) AND past(x)).

Declared limits: the straight segment is exact only to first order in the field (Fermat), so it biases Q at order
a^2, the order of the signal itself; nonlocality (l_k = eps^(-1/4)) and the finite box bias it too. The MEAN of Q is
therefore reported, not judged. The test judges only the NOISE, which these biases do not change materially.

Model (locked before any run)
-----------------------------
3+1 box t in [0, 8], |x|, |y|, |z| <= 8, rho = 1 (~33 000 elements flat); w = 2.5, a = 0.2; rows: every element
with t in [6.5, 8] and r <= 6 in each configuration; eps = 0.3 and 0.1 (both from the same bitsets); 48 sprinklings.

Criteria (stated before any full run; changed after the first quick run, see below)
-----------------------------------------------------------------------------------
INSTRUMENT_FAILED if the perturbed code path does not reduce to the flat one: on the first sprinkling, the
configuration with a = 1e-9 must reproduce S(0) to |S(1e-9) - S(0)| <= 1e-6 (|S(0)| + 1) for every eps.
For each eps: sigma_Q = sd of Q over sprinklings; M_needed = (sigma_Q / (0.1 |Q_pred|))^2, the number of sprinklings
that would measure the Psi-sector coefficient to 10 %.
PASSED (feasible) if M_needed <= 1000 for at least one eps; FAILED (not feasible at this size) otherwise.
Reported: sigma of S(0) alone, corr(S(0), S(+a)) (how much the correlated sampling helps), mean Q against Q_pred,
the elements and rows per configuration, time per sprinkling, and the scaling hint sigma_Q / Q_pred.

Quick run (3 sprinklings, box 5, decides nothing): mechanics in order, 5 s per sprinkling of ~5000 elements. It
showed corr(S(0), S(+a)) = 0.88 and sd(Q) ~ 60 x Q_pred. The first instrument criterion (corr >= 0.9) mixed the
premise of the method with its outcome: a low correlation is a property of the BD action (the field reshuffles the
interval counts of near-null pairs), not a fault of the code. Before the full run it was replaced by the exact
null check above, which tests the code; the correlation is now reported and enters sigma_Q. eps = 0.03 was added
(the same bitsets, only another kernel). Nothing else changed.
"""

from __future__ import annotations

import math
import time

import numpy as np
from scipy.special import erf

from core.operators.bd import BD_NORM, kernel_cutoff_n, smeared_kernel
from core.order.sprinkle import derive_seed, interval_volume
from core.stats.parallel import parallel_map
from core.testing import Result

TEST = {
    "id": "K.17",
    "phase": "K",
    "title": "Proveditelnost: lze na sprinklingu 3+1 změřit kvadratickou formu BD akce?",
    "theory_part": "odvození pravidla substrátu z akce řádu (pravidlo_z_akce.md); předpoklad pro F6.8",
    "question": "Je šum diskrétní BD akce při korelovaném vzorkování dost malý, aby šlo v boxu laboratoře změřit "
                "kvadratickou (Einsteinovu–Hilbertovu) odezvu na statické pole Ψ s přesností 10 % do 1000 sprinklingů?",
    "method": "3+1 box t [0, 8], |x| <= 8, ρ = 1; metrika −dt² + e^{2af}dx², f = exp(−r²/2w²), w = 2,5, a = 0,2; "
              "jedna nadmnožina bodů pro konfigurace 0, +a, −a (ztenčení e^{3af}); relace přes optickou délku úsečky; "
              "n_xy z bitsetů; S = −(2/ρ) Σ (B_ε 1)(x) přes prvky s t ∈ [6,5; 8], r <= 6; Q = (S₊ + S₋)/2 − S₀; "
              "ε = 0,3, 0,1 a 0,03; 48 sprinklingů",
    "prediction": "Q_pred = a² ∫ f² (12/w² − 6r²/w⁴) d³x dt (Einstein–Hilbert); šum rozhodne o proveditelnosti",
    "criterion": "MĚŘIDLO: konfigurace a = 1e−9 nereprodukuje S(0) (|ΔS| > 1e−6 (|S0| + 1)). PROŠEL (proveditelné): M_needed = (σ_Q / 0,1 |Q_pred|)² <= 1000 "
                 "pro aspoň jedno ε; NEPROŠEL jinak. Střední Q proti Q_pred jen hlášeno (úsečka a nelokalita).",
    "null_model": "a = 0 (konfigurace S(0))",
    "level": "calibration",
    "meters": [],
    "depends_on": ["K.13a"],
    "params": {"T": 8.0, "half": 8.0, "rho": 1.0, "w": 2.5, "a": 0.2, "band": [6.5, 8.0], "R_rows": 6.0,
               "eps": [0.3, 0.1, 0.03], "n_seeds": 48, "null_a": 1e-9, "null_tol": 1e-6, "rel_prec": 0.1, "M_max": 1000, "block": 256},
    "parts": {"main": {}},
    "quick": {"n_seeds": 3, "T": 5.0, "half": 5.0, "band": [4.0, 5.0], "R_rows": 3.0, "w": 1.5, "eps": [0.3]},
    "estimate": "~30–60 min s 8 procesy (na sprinkling 3 konfigurace po ~33 000 prvcích: relace N² ~ 1 min, bitsety "
                "~140 MB, popcount řádků); paměť ~0,6 GB na proces",
}


# ============================================================================================== metric and relation

def line_integrals(P, Q, w, orders=(1, 2, 3)):
    """int_0^L f^k ds along the straight segments P -> Q (arrays (n, 3)), f = exp(-r^2 / 2 w^2), closed form."""
    d = Q - P
    L = np.sqrt((d * d).sum(axis=1))
    Ls = np.maximum(L, 1e-300)
    n = d / Ls[:, None]
    s0 = (P * n).sum(axis=1)
    b2 = np.maximum((P * P).sum(axis=1) - s0 * s0, 0.0)
    out = {}
    for k in orders:
        c = math.sqrt(k / 2.0) / w
        out[k] = np.exp(-k * b2 / (2 * w * w)) * (w * math.sqrt(math.pi / (2.0 * k))) * (erf((s0 + L) * c)
                                                                                         - erf(s0 * c))
    return L, out


def seg_fmax(P, Q, w):
    """Largest f on the segment: exp(-dmin^2 / 2 w^2), dmin the distance of the segment from the origin."""
    d = Q - P
    L2 = np.maximum((d * d).sum(axis=1), 1e-300)
    s = np.clip(-(P * d).sum(axis=1) / L2, 0.0, 1.0)
    c = P + s[:, None] * d
    return np.exp(-(c * c).sum(axis=1) / (2 * w * w))


def relation_rows(X, i0, i1, a, w):
    """Boolean (i1 - i0, N): X[i] < X[j] for i in [i0, i1) (time-sorted X; only j > i can be related)."""
    N = len(X)
    A = X[i0:i1]
    dt = X[None, :, 0] - A[:, None, 0]
    dx = X[None, :, 1:] - A[:, None, 1:]
    L = np.sqrt((dx * dx).sum(axis=2))
    rel = (dt > 0) & (dt >= L)
    if a == 0.0:
        return rel
    # bounds of the inverse light speed e^{a f} on the segment: e^{-|a| fmax} .. e^{|a| fmax}
    ii, jj = np.nonzero((dt > 0) & (dt >= L * math.exp(-abs(a))) & (dt <= L * math.exp(abs(a))))
    if ii.size:
        P, Qp = A[ii, 1:], X[jj, 1:]
        fm = seg_fmax(P, Qp, w)
        Lij, dtij = L[ii, jj], dt[ii, jj]
        hi = Lij * np.exp(abs(a) * fm)
        lo = Lij * np.exp(-abs(a) * fm)
        sure_yes = dtij >= hi
        sure_no = dtij < lo
        amb = ~(sure_yes | sure_no)
        res = sure_yes.copy()
        if amb.any():
            _, I = line_integrals(P[amb], Qp[amb], w)
            opt = Lij[amb] + a * I[1] + a * a / 2.0 * I[2] + a ** 3 / 6.0 * I[3]
            res[amb] = dtij[amb] >= opt
        rel[ii, jj] = res
    return rel


def future_bitsets(X, a, w, block):
    """future(y) as packed bit rows (uint64), one per element."""
    N = len(X)
    nw = (N + 63) // 64
    bits = np.zeros((N, nw), dtype=np.uint64)
    for i0 in range(0, N, block):
        i1 = min(N, i0 + block)
        rel = relation_rows(X, i0, i1, a, w)
        pk = np.packbits(rel, axis=1, bitorder="little")
        pad = np.zeros((i1 - i0, nw * 8), dtype=np.uint8)
        pad[:, :pk.shape[1]] = pk
        bits[i0:i1] = pad.view(np.uint64)
    return bits


def column(bits, x):
    """Boolean (N,): y < x, read from the future bitsets (bit x of every row)."""
    word, bit = divmod(x, 64)
    return ((bits[:, word] >> np.uint64(bit)) & np.uint64(1)).astype(bool)


# ============================================================================================== action

def action_sums(X, bits, rows, eps_list, rho):
    """For each eps: sum over the rows of (1 + sum_{y < x} K_eps(n_xy)); n by popcount."""
    D = 4
    nmax = max(kernel_cutoff_n(e, D) for e in eps_list)
    tc2 = ((1.5 * nmax + 50.0) / (rho * interval_volume(D, 1.0))) ** 0.5 * 1.3  # proper-time^2 cutoff, with slack
    N = len(X)
    nw = bits.shape[1]
    sums = np.zeros(len(eps_list))
    for x in rows:
        past = column(bits, x)
        past[x] = False
        cand = np.flatnonzero(past)
        if cand.size == 0:
            sums += 1.0
            continue
        d = X[x] - X[cand]
        tau2 = d[:, 0] ** 2 - (d[:, 1:] ** 2).sum(axis=1)
        cand = cand[tau2 <= tc2]
        pk = np.packbits(past, bitorder="little")
        pad = np.zeros(nw * 8, dtype=np.uint8)
        pad[:pk.size] = pk
        pw = pad.view(np.uint64)
        n = np.zeros(cand.size)
        for c0 in range(0, cand.size, 1024):
            n[c0:c0 + 1024] = np.bitwise_count(bits[cand[c0:c0 + 1024]] & pw[None, :]).sum(axis=1)
        for k, e in enumerate(eps_list):
            sums[k] += 1.0 + float(smeared_kernel(n, D, e).sum())
    return sums


def sprinkle_task(p, seed, null_check=False):
    t0 = time.time()
    rng = np.random.default_rng(seed)
    T, H, rho, w, a = float(p["T"]), float(p["half"]), float(p["rho"]), float(p["w"]), float(p["a"])
    rho_max = rho * math.exp(3 * abs(a))
    m = rng.poisson(rho_max * T * (2 * H) ** 3)
    P = np.column_stack([rng.uniform(0, T, m)] + [rng.uniform(-H, H, m) for _ in range(3)])
    u = rng.random(m)
    r = np.sqrt((P[:, 1:] ** 2).sum(axis=1))
    f = np.exp(-r * r / (2 * w * w))
    out = {"N": {}, "rows": {}, "S": {}}
    configs = [("0", 0.0), ("+", a), ("-", -a)] + ([("null", float(p["null_a"]))] if null_check else [])
    for name, aa in configs:
        keep = u < np.exp(3 * aa * f) * rho / rho_max
        X = P[keep]
        o = np.argsort(X[:, 0], kind="stable")
        X = X[o]
        rr = np.sqrt((X[:, 1:] ** 2).sum(axis=1))
        rows = np.flatnonzero((X[:, 0] >= p["band"][0]) & (X[:, 0] <= p["band"][1]) & (rr <= p["R_rows"]))
        bits = future_bitsets(X, aa, w, int(p["block"]))
        sums = action_sums(X, bits, rows, p["eps"], rho)
        # S = -(2/rho) sum (B 1) = (2/rho) a_4 sqrt(eps rho) sum (1 + sum K)
        out["S"][name] = [2.0 / rho * BD_NORM[4] * math.sqrt(e * rho) * s for e, s in zip(p["eps"], sums)]
        out["N"][name], out["rows"][name] = len(X), len(rows)
        del bits
    out["time"] = time.time() - t0
    return out


def q_pred(p):
    """Continuum Einstein-Hilbert even part over the rows: a^2 int f^2 (12/w^2 - 6 r^2/w^4) 4 pi r^2 dr x band."""
    w, a, R = float(p["w"]), float(p["a"]), float(p["R_rows"])
    r = np.linspace(0.0, R, 20001)
    g = np.exp(-r * r / (w * w)) * (12 / w ** 2 - 6 * r * r / w ** 4) * 4 * math.pi * r * r
    return a * a * float(getattr(np, "trapezoid", getattr(np, "trapz", None))(g, r)) * (p["band"][1] - p["band"][0])


def run(p, ctx):
    t0 = time.time()
    n = int(p["n_seeds"])
    jobs = [(p, derive_seed(ctx.seed, "K.17", k), k == 0) for k in range(n)]
    res = parallel_map(sprinkle_task, jobs, ctx.workers, ctx.progress,
                       lambda i, r_, dd, nn: f"[K.17] sprinkling {dd}/{nn} ({r_['time']:.0f}s, total "
                                             f"{time.time() - t0:.0f}s)")
    Qp = q_pred(p)
    lines = [f"3+1 box t [0, {p['T']:g}], |x| <= {p['half']:g}, rho = {p['rho']:g}; metric -dt^2 + e^(2 a f) dx^2, "
             f"w = {p['w']:g}, a = {p['a']:g}; rows t in {p['band']}, r <= {p['R_rows']:g}; {n} sprinklings; "
             f"elements (0/+/-) ~{np.mean([r_['N']['0'] for r_ in res]):.0f} / {np.mean([r_['N']['+'] for r_ in res]):.0f}"
             f" / {np.mean([r_['N']['-'] for r_ in res]):.0f}, rows ~{np.mean([r_['rows']['0'] for r_ in res]):.0f}; "
             f"{np.mean([r_['time'] for r_ in res]):.0f} s per sprinkling",
             f"continuum prediction Q_pred = {Qp:.4f} (Einstein-Hilbert even part over the rows)"]
    feas, inst_ok = [], True
    first = res[0]
    for k, e in enumerate(p["eps"]):
        s0, sn = first["S"]["0"][k], first["S"]["null"][k]
        ok = abs(sn - s0) <= float(p["null_tol"]) * (abs(s0) + 1.0)
        inst_ok &= ok
        lines.append(f"  instrument: eps = {e:g}, S(a = {p['null_a']:g}) - S(0) = {sn - s0:+.3e} -> {'OK' if ok else 'FAILS'}")
    data = {"Q_pred": Qp, "eps": {}}
    for k, e in enumerate(p["eps"]):
        S0 = np.array([r_["S"]["0"][k] for r_ in res])
        Sp = np.array([r_["S"]["+"][k] for r_ in res])
        Sm = np.array([r_["S"]["-"][k] for r_ in res])
        Q = 0.5 * (Sp + Sm) - S0
        corr = float(np.corrcoef(S0, Sp)[0, 1]) if n > 2 else float("nan")
        sQ = float(Q.std(ddof=1)) if n > 1 else float("nan")
        M_need = (sQ / (float(p["rel_prec"]) * abs(Qp))) ** 2 if Qp != 0 else float("inf")
        feas.append(M_need <= float(p["M_max"]))
        lines.append(f"  eps = {e:g} (l_k = {e ** -0.25:.2f}): S(0) = {S0.mean():+.2f} ± {S0.std(ddof=1):.2f} (sd); "
                     f"corr(S0, S+) = {corr:.4f}; Q = {Q.mean():+.4f} ± {sQ / math.sqrt(n):.4f} (sd {sQ:.4f}) against "
                     f"Q_pred {Qp:+.4f}; M_needed for 10 % = {M_need:.3g}")
        data["eps"][str(e)] = {"S0": S0.tolist(), "Sp": Sp.tolist(), "Sm": Sm.tolist(), "corr": corr, "sd_Q": sQ,
                               "M_needed": M_need}
    lines.append("  the mean of Q is biased at the order of the signal by the straight-segment relation, the "
                 "nonlocality and the box; only the noise is judged")
    if not inst_ok:
        status = "INSTRUMENT_FAILED"
    elif any(feas):
        status = "PASSED"
    else:
        status = "FAILED"
    return Result(status, lines, data)
