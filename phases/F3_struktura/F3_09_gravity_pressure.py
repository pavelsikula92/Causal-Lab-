"""
F3.9 -- Gravity and pressure: does the field with range plus a finite channel capacity give a bound bundle?
===========================================================================================================

Why this test exists (23 Oct 2026)
----------------------------------
F3.8 search (no candidate, recorded under "Objevy"): the retarded field
of the bundle attracts it for the first time. epsilon = 0.1 stops the
growth (slope -0.004 +- 0.005 in W1) and then collapses (-0.036 in W2);
epsilon = 0.5 collapses (R 3.0 -> 2.1) and rebounds to the walls. The
bundle behaves like cold dust: nothing resists compression. New hypothesis
after that failure (own criteria; nothing of F3.8 is changed): an
equilibrium needs pressure, and the theory has it: the finite capacity of
a channel (axiom 4; F1.8: capacity limits the flux, not the number of
links, so this is NOT the refuted degree bound K0). When the bundle
contracts, the flux concentrates into fewer carriers (73 -> 15 in F3.8);
with a finite capacity the excess must take less aligned channels and the
bundle resists compression.

Model (locked before any run; protocol of the author as in F3.7/F3.8)
---------------------------------------------------------------------
Everything as F3.8 (2+1 sprinkling rho = 1, source tube radius 4 / time
[0, 4) with u_src; carried direction; retarded field
phi = sum F / (2 pi tau) over the past flux with tau >= 1; ln Pi = -lambda phi,
lambda = epsilon / phi_ref; bending of the carried direction by the
transverse gradient of ln Pi from the parents) plus:
Capacity of a channel: an event v with flux F_v <= C sends everything to
its child of smallest rapidity (as F3.8); with F_v > C it fills its
children in order of increasing rapidity relative to its direction, at
most C each (the last one takes the remainder; if C times the number of
children is not enough, the rest goes to the best child). All of them
receive v's (bent) direction: momentum is conserved, the excess only takes
less aligned paths. C = kappa F_ref, F_ref = the median flux per carrying
event of the free arm in W1 (calibration sprinklings, own stream).
Grid fixed now: epsilon in {0.01, 0.1, 0.5} (the grid of F3.8, kept, not
chosen from its results) x kappa in {4, 1, 1/4} (4: the capacity acts
only where the flux concentrates; 1: at a typical carrier; 1/4: almost
everywhere).
Arms: free (no field, no capacity); field_eps (F3.8, no capacity; a
replication on new sprinklings, reported, not a candidate); cap_eps_kappa
(field + capacity; the candidates).
Covariance checked bit for bit in a frame boosted by beta' = 0.6 (the
candidates; the decisions include the number of channels used). As in
F3.6-F3.8: the rules use the metric at the scale of a link, and choosing
among the children uses the future of v (criteria (a), (b) of the audit).
Remark declared now: in 2+1 the potential is logarithmic and confines
every orbit; a PASSED here asks for the 3+1 version before any physical
reading.

Classification and protocol as F3.7/F3.8: bound |s| + 2 SE <= s_free/4 in
W1 [12, 28] and W2 [28, 44]; growing, collapsing, unresolved, wall-limited;
if the free arm does not grow in a window (s_free <= 0) nothing is
classified there (unresolved).
search RECORDED; candidate among the cap arms: bound in W1 and W2 with the
frame check, the smallest |s_W1| + 2 SE, ties: weaker epsilon, then larger
kappa; at most one, "none" is a result. confirm (settable epsilon, kappa:
CHOSEN FROM DATA) on new streams, rest and moving source (beta = 0.5):
PASSED bound at rest in W1 and W2, frame check, |v_c| + 2 SE <= beta/4 and
bound in motion; FAILED growing or collapsing at rest, or
|v_c| - 2 SE > beta/4; INSTRUMENT_FAILED frame check; else UNDECIDED.

Seen before this test was written: F3.6, F3.7 search, F3.8 search. The cap
arms had not been run. After everything above was written: quick runs
(3 sprinklings, box 24 x 16 x 16; the free arm itself does not grow
there in W2, so the quick classification is meaningless: the guard
"s_free <= 0 -> unresolved" was added after them, it can only remove a
verdict) and one full-size timing sprinkling (slopes W1 / W2: free
0.087 / 0.079, cap_0.1_1 0.013 / 0.002, cap_0.5_0.25 -0.076 / 0.021; one
sprinkling, decides nothing). Nothing else was changed.
"""

from __future__ import annotations

import math
import time

import numpy as np

from core.order.links import links_csr
from core.order.sprinkle import Box, derive_seed, sprinkle_poisson
from core.stats.parallel import parallel_map
from core.stats.stats import boot_se
from core.testing import Result

TEST = {
    "id": "F3.9",
    "phase": "F3",
    "title": "Gravitace a tlak: dá pole s dosahem spolu s konečnou kapacitou kanálu vázaný svazek?",
    "theory_part": "hmota jako rovnováha gravitace pole (A5′, F3.8) a tlaku konečné kapacity toku (axiom 4; ne K₀); vázanost "
                   "hmoty, setrvačnost",
    "question": "Ustálí se šířka svazku s neseným směrem, když ho k sobě táhne jeho retardované pole (F3.8) a konečná kapacita "
                "kanálu nutí přebytek toku do méně souběžných potomků?",
    "method": "F3.8 + kapacita kanálu C = κ F_ref (událost s tokem nad C plní potomky podle rapidity, nejvýš C každému); "
              "mřížka ε ∈ {0,01; 0,1; 0,5} × κ ∈ {4, 1, 1/4}; protokol search → nejvýš jeden kandidát → confirm "
              "na nových sprinklinzích včetně pohyblivého zdroje",
    "prediction": "podle hypotézy po F3.8: gravitace pole a tlak kapacity dají ustálenou šířku (vázaný svazek)",
    "criterion": "jako F3.8: vázaný |s| + 2 SE <= s_free/4 v oknech [12, 28] a [28, 44]; search: kandidát mezi větvemi "
                 "s kapacitou, vázaný v obou oknech se shodnou volbou v soustavě β′ = 0,6, nejmenší |s| + 2 SE (slabší ε, pak "
                 "větší κ); confirm: PROŠEL vázaný v obou oknech, kontrola soustavy, |v_c| + 2 SE <= β/4 a vázaný v pohybu; "
                 "NEPROŠEL rostoucí nebo kolaps, nebo |v_c| − 2 SE > β/4; MĚŘIDLO kontrola soustavy; jinak NEROZHODNUTO. "
                 "Úspěch ve 2+1 vyžaduje verzi 3+1 před fyzikálním čtením.",
    "null_model": "free (bez pole a kapacity) a field (F3.8, bez kapacity)",
    "level": "mechanism",
    "meters": ["K.6"],
    "depends_on": ["F3.6", "F3.8"],
    "params": {"D": 3, "T": 48.0, "half": 16.0, "r_src": 4.0, "t_src": 4.0, "slab": 2.0, "W1": [12.0, 28.0],
               "W2": [28.0, 44.0], "win_move": [5.0, 28.0], "beta_check": 0.6, "beta_move": 0.5, "tau_c": 1.0,
               "epsilons": [0.01, 0.1, 0.5], "kappas": [4.0, 1.0, 0.25], "L_g": 2.0, "min_parents": 4, "n_cal": 8,
               "n_seeds": 32, "n_se": 2.0, "boot_B": 4000, "epsilon": None, "kappa": None},
    "parts": {"search": {"mode": "search"}, "confirm": {"mode": "confirm"}},
    "settable": {"confirm": ["epsilon", "kappa"]},
    "quick": {"n_seeds": 3, "n_cal": 1, "boot_B": 200, "T": 24.0, "half": 8.0, "r_src": 2.0, "t_src": 2.0,
              "W1": [6.0, 14.0], "W2": [14.0, 22.0], "win_move": [3.0, 14.0], "slab": 1.0},
    "estimate": "search ~25–40 min s 10 procesy (13 větví, přepočet v jiné soustavě pro 9 větví s kapacitou; jeden "
                "sprinkling ~6–8 min); confirm ~10–15 min",
}

# ============================================================================================== geometry

def boost_points(X: np.ndarray, beta: float) -> np.ndarray:
    """Coordinates (or 4-vectors) in the frame moving with velocity beta along x (c = 1)."""
    if beta == 0.0:
        return np.array(X, float, copy=True)
    g = 1.0 / math.sqrt(1.0 - beta * beta)
    Y = np.array(X, float, copy=True)
    Y[..., 0] = g * (X[..., 0] - beta * X[..., 1])
    Y[..., 1] = g * (X[..., 1] - beta * X[..., 0])
    return Y


def mink(a, b):
    """Minkowski product (+ - -) over the last axis."""
    return a[..., 0] * b[..., 0] - (a[..., 1:] * b[..., 1:]).sum(axis=-1)


def children_csr(ptr, idx, N):
    child = np.repeat(np.arange(N), np.diff(ptr))
    o = np.argsort(idx, kind="stable")
    cptr = np.concatenate([[0], np.cumsum(np.bincount(idx, minlength=N))]).astype(np.int64)
    return cptr, child[o]


def rest_basis(u: np.ndarray) -> np.ndarray:
    """Orthonormal spacelike basis (rows e_i, e_i . e_i = -1, e_i . u = 0) of the rest space of the unit vector u (Gram-Schmidt)."""
    D = u.size
    out = []
    for k in range(1, D):
        e = np.zeros(D)
        e[k] = 1.0
        e = e - mink(e, u) * u                      # remove the component along u (u . u = 1)
        for f in out:
            e = e + mink(e, f) * f                  # f . f = -1: e - (e.f)/(f.f) f = e + (e.f) f
        n = math.sqrt(max(-float(mink(e, e)), 0.0))
        out.append(e / n)
    return np.array(out)


def boost_towards(u: np.ndarray, e: np.ndarray, eta: float) -> np.ndarray:
    """u boosted by rapidity eta along the unit spacelike vector e orthogonal to u: cosh(eta) u + sinh(eta) e."""
    return math.cosh(eta) * u + math.sinh(eta) * e



# ============================================================================================== field and transport

class RetardedField:
    """phi(x) = sum over the flux events w in the causal past of x with tau >= tau_c of F_w / (2 pi tau) (rho = 1)."""

    def __init__(self, X, tau_c):
        self.X = X
        self.tau_c2 = float(tau_c) ** 2
        self.pos = []                               # finished flux events (positions) and their flux
        self.flx = []
        self._P = np.zeros((0, X.shape[1]))
        self._F = np.zeros(0)
        self._dirty = False
        self.cache = {}

    def add(self, w, F):
        self.pos.append(self.X[w])
        self.flx.append(F)
        self._dirty = True

    def value(self, x_idx):
        if x_idx in self.cache:
            return self.cache[x_idx]
        if self._dirty:
            self._P = np.array(self.pos)
            self._F = np.array(self.flx)
            self._dirty = False
        if self._F.size == 0:
            return 0.0
        d = self.X[x_idx][None, :] - self._P
        t2 = mink(d, d)
        ok = (d[:, 0] > 0.0) & (t2 >= self.tau_c2)
        val = float((self._F[ok] / (2.0 * math.pi * np.sqrt(t2[ok]))).sum())
        # safe to cache: the field is asked for only at events whose whole past is already processed
        self.cache[x_idx] = val
        return val


def transport(X, ptr, idx, cptr, cidx, s, u_src, arm, lam, p, want_phi=False, C=np.inf):
    """
    One causal pass. arm: 'free' or 'field'. Events are processed in time order; an event's flux is final when it
    is processed, and it then joins the field sources (so every field value uses only the finished flux of its
    past: all of x's past is processed before x). The field at the parents of v is evaluated when v is processed.
    Returns flux, chosen children, bending rapidities, and (want_phi) the field at the carrying events.
    """
    N, D = X.shape
    L_g, nmin = float(p["L_g"]), int(p["min_parents"])
    field = RetardedField(X, float(p["tau_c"]))
    inflow = np.zeros(N)
    U = np.zeros((N, D))
    flux = np.zeros(N)
    decision = np.full(N, -1, dtype=np.int64)
    bend = np.full(N, np.nan)
    phi_car = np.full(N, np.nan)
    order_t = X[:, 0]
    for v in range(N):
        tot = inflow[v] + s[v]
        if tot <= 0.0:
            continue
        # the field value of v and of its parents use only events before v; flux events with t < t_v are final
        P = U[v] + s[v] * u_src
        u = P / math.sqrt(max(float(mink(P, P)), 1e-300))
        flux[v] = tot
        if want_phi:
            phi_car[v] = field.value(v)
        a, b = cptr[v], cptr[v + 1]
        send = u
        if b > a:
            ks = cidx[a:b]
            B = X[ks] - X[v]
            tau = np.sqrt(np.maximum(mink(B, B), 1e-300))
            ch = mink(u[None, :], B) / tau
            j = int(np.argmin(ch))
            k = int(ks[j])
            decision[v] = k
            if arm == "field" and lam > 0.0:
                pa, pb = ptr[v], ptr[v + 1]
                par = idx[pa:pb]
                if par.size:
                    d = X[par] - X[v]
                    lag = -mink(d, u[None, :])
                    sel = lag <= L_g
                    if int(sel.sum()) >= nmin:
                        pars = par[sel]
                        dd = d[sel]
                        pp = dd - mink(dd, u[None, :])[:, None] * u[None, :]
                        E = rest_basis(u)
                        q = -np.stack([mink(pp, e[None, :]) for e in E], axis=1)
                        lnPi = np.array([-lam * field.value(int(w)) for w in pars])
                        A = np.column_stack([np.ones(q.shape[0]), q])
                        coef, *_ = np.linalg.lstsq(A, lnPi, rcond=None)
                        G = coef[1:]
                        gn = float(np.sqrt((G * G).sum()))
                        if gn > 0.0:
                            e_dir = -(G[0] * E[0] + G[1] * E[1]) / gn
                            eta = gn * float(tau[j])
                            send = boost_towards(u, e_dir, eta)
                            bend[v] = eta
            if tot <= C:
                inflow[k] += tot
                U[k] += tot * send
                n_used = 1
            else:
                # finite capacity of a channel (axiom 4): the most aligned children first, at most C each
                # (water-filling; the last one takes the remainder); all carry the (bent) direction of v
                order = np.argsort(ch, kind="stable")
                rest, n_used = tot, 0
                for jj in order:
                    give = min(C, rest)
                    kk = int(ks[jj])
                    inflow[kk] += give
                    U[kk] += give * send
                    rest -= give
                    n_used += 1
                    if rest <= 0.0:
                        break
                if rest > 0.0:                           # more flux than C times the children: the rest to the best child
                    inflow[k] += rest
                    U[k] += rest * send
            decision[v] = k * 100000 + n_used
        field.add(v, tot)                              # v's flux is final: it sources the field of its future
    return {"flux": flux, "decision": decision, "bend": bend, "phi": phi_car}

# ============================================================================================== readout

def readout(Xs, flux, p, beta):
    T, sl = float(p["T"]), float(p["slab"])
    t_end = T * math.sqrt(1.0 - beta * beta) - sl
    edges = np.arange(float(p["t_src"]), t_end + 1e-9, sl)
    out = {k: [] for k in ("t", "R", "xc", "n", "F")}
    car = flux > 0.0
    for a, b in zip(edges[:-1], edges[1:]):
        m = car & (Xs[:, 0] >= a) & (Xs[:, 0] < b)
        w = flux[m]
        out["t"].append(0.5 * (a + b))
        out["n"].append(int(m.sum()))
        out["F"].append(float(w.sum()))
        if w.sum() > 0:
            out["R"].append(float(np.sqrt((w * (Xs[m, 1:] ** 2).sum(axis=1)).sum() / w.sum())))
            out["xc"].append(float((w * Xs[m, 1]).sum() / w.sum()))
        else:
            out["R"].append(float("nan"))
            out["xc"].append(float("nan"))
    return out


def slope(t, y, win):
    t, y = np.asarray(t, float), np.asarray(y, float)
    m = (t >= win[0]) & (t <= win[1]) & np.isfinite(y)
    return float(np.polyfit(t[m], y[m], 1)[0]) if m.sum() >= 2 else float("nan")



# ============================================================================================== tasks

def setup(p, seed, beta):
    D, T, H = int(p["D"]), float(p["T"]), float(p["half"])
    rng = np.random.default_rng(seed)
    X = sprinkle_poisson(Box(D, [0.0, -H, -H], [T, H + beta * T, H]), 1.0, rng)
    ptr, idx = links_csr(X, 1.0)
    cptr, cidx = children_csr(ptr, idx, len(X))
    g = 1.0 / math.sqrt(1.0 - beta * beta)
    u_src = np.array([g, g * beta, 0.0])
    Xs = boost_points(X, beta)
    tube = (Xs[:, 0] >= 0.0) & (Xs[:, 0] < float(p["t_src"])) & ((Xs[:, 1:] ** 2).sum(axis=1) <= float(p["r_src"]) ** 2)
    return X, ptr, idx, cptr, cidx, Xs, tube.astype(float), u_src


def calibration_task(p, seed):
    """Flux per carrier and field at the carriers of the free arm in W1 (no back-reaction): F_ref and phi_ref."""
    X, ptr, idx, cptr, cidx, Xs, s, u_src = setup(p, seed, 0.0)
    r = transport(X, ptr, idx, cptr, cidx, s, u_src, "free", 0.0, p, want_phi=True)
    w1 = p["W1"]
    m = (r["flux"] > 0) & (Xs[:, 0] >= w1[0]) & (Xs[:, 0] <= w1[1]) & np.isfinite(r["phi"])
    return r["flux"][m], r["phi"][m]


def run_task(p, seed, beta, arms, phi_ref, F_ref, check):
    """One sprinkling: the arms (key, field strength epsilon or None, capacity kappa or None), readout, frame check."""
    X, ptr, idx, cptr, cidx, Xs, s, u_src = setup(p, seed, beta)
    bc = float(p["beta_check"])
    Xb = boost_points(X, bc)
    ub = boost_points(u_src[None, :], bc)[0]
    out = {"N": len(X), "sources": int(s.sum()), "arms": {}}
    for key, eps, kap in arms:
        lam = float(eps) / phi_ref if eps is not None else 0.0
        C = float(kap) * F_ref if kap is not None else np.inf
        arm = "field" if eps is not None else "free"
        r = transport(X, ptr, idx, cptr, cidx, s, u_src, arm, lam, p, C=C)
        if key in check:
            rb = transport(Xb, ptr, idx, cptr, cidx, s, ub, arm, lam, p, C=C)
            same = bool(np.array_equal(r["decision"], rb["decision"]))
        else:
            same = True
        ser = readout(Xs, r["flux"], p, beta)
        out["arms"][key] = {"series": ser, "s1": slope(ser["t"], ser["R"], p["W1"]), "s2": slope(ser["t"], ser["R"], p["W2"]),
                            "vc": slope(ser["t"], ser["xc"], p["win_move"]), "frame_identical": same,
                            "bend_mean": float(np.nanmean(r["bend"])) if np.isfinite(r["bend"]).any() else float("nan"),
                            "Pi_min": float("nan")}
    return out

# ============================================================================================== statistics

def classify(sv, se, sfree, nse):
    if not all(np.isfinite([sv, se, sfree])) or sfree <= 0.0:
        return "unresolved"                     # no reference growth to compare with (the free arm must grow)
    tol = sfree / 4.0
    if abs(sv) + nse * se <= tol:
        return "bound"
    if sv - nse * se > tol:
        return "growing"
    if sv + nse * se < -tol:
        return "collapsing"
    return "unresolved"


def summarise(rows, keys, p, ctx, tag):
    B, nse = int(p["boot_B"]), float(p["n_se"])
    st = {}
    for key in keys:
        bs = lambda *t, key=key: derive_seed(ctx.seed, "F3.7-boot", tag, key, *t)
        st[key] = {w: boot_se([r_["arms"][key][w] for r_ in rows], B, bs(w))[:2] for w in ("s1", "s2", "vc")}
        st[key]["frame"] = all(r_["arms"][key]["frame_identical"] for r_ in rows)
        with np.errstate(all="ignore"):
            st[key]["R"] = np.nanmean(np.array([r_["arms"][key]["series"]["R"] for r_ in rows], float), axis=0)
            st[key]["n"] = np.nanmean(np.array([r_["arms"][key]["series"]["n"] for r_ in rows], float), axis=0)
            st[key]["F"] = np.nanmean(np.array([r_["arms"][key]["series"]["F"] for r_ in rows], float), axis=0)
            st[key]["bend"] = float(np.nanmean([r_["arms"][key]["bend_mean"] for r_ in rows]))
            st[key]["Pi_min"] = float(np.nanmin([r_["arms"][key]["Pi_min"] for r_ in rows]))
    sfree1, sfree2 = st["free"]["s1"][0], st["free"]["s2"][0]
    t = np.asarray(rows[0]["arms"]["free"]["series"]["t"], float)
    wall = float(p["half"]) / 2.0
    for key in keys:
        st[key]["c1"] = classify(*st[key]["s1"], sfree1, nse)
        st[key]["c2"] = classify(*st[key]["s2"], sfree2, nse)
        for c, w in (("c1", p["W1"]), ("c2", p["W2"])):     # a width that reached the box cannot be read as bound
            m = (t >= w[0]) & (t <= w[1])
            if m.any() and np.nanmax(st[key]["R"][m]) > wall and st[key][c] in ("bound", "unresolved", "collapsing"):
                st[key][c] = "wall-limited"
    return st


def lines_for(st, keys, t):
    out = []
    for key in keys:
        x = st[key]
        out.append(f"  {key:12s}: slope W1 {x['s1'][0]:+.4f} +- {x['s1'][1]:.4f} -> {x['c1']}; W2 {x['s2'][0]:+.4f} +- "
                   f"{x['s2'][1]:.4f} -> {x['c2']}; frame check {'OK' if x['frame'] else 'FAILS'}"
                   + (f"; mean bending {x['bend']:.4f}" if np.isfinite(x["bend"]) else "")
                   )
        out.append("        R_rms  " + " ".join(f"{v:5.2f}" if np.isfinite(v) else "  nan" for v in x["R"]))
        out.append("        carry  " + " ".join(f"{v:5.1f}" for v in x["n"]))
    out.insert(0, "        t      " + " ".join(f"{v:5.1f}" for v in t))
    return out



# ============================================================================================== run

def run(p, ctx):
    mode = str(p["mode"])
    t0 = time.time()
    n = int(p["n_seeds"])
    if mode == "confirm":
        eps, kap = p.get("epsilon"), p.get("kappa")
        if eps is None or kap is None:
            raise SystemExit("F3.9 confirm: give the candidate of the search, e.g. --set epsilon=0.1 --set kappa=1")
        eps, kap = float(eps), float(kap)
        if eps not in [float(e) for e in p["epsilons"]] or kap not in [float(k) for k in p["kappas"]]:
            raise SystemExit(f"F3.9 confirm: epsilon and kappa must be on the grid {p['epsilons']} x {p['kappas']}")
    stream = "F3.9-" + mode
    cal = parallel_map(calibration_task, [(p, derive_seed(ctx.seed, stream + "-cal", k)) for k in range(int(p["n_cal"]))],
                       ctx.workers, ctx.progress, lambda i, r_, d, nn: f"[F3.9] calibration {d}/{nn} ({time.time()-t0:.0f}s)")
    F_ref = float(np.median(np.concatenate([c[0] for c in cal])))
    phi_ref = float(np.median(np.concatenate([c[1] for c in cal])))
    head = [f"2+1 sprinkling rho = 1, box t [0, {p['T']:g}], half width {p['half']:g}; source radius {p['r_src']:g}, time "
            f"[0, {p['t_src']:g}); field phi = sum F / (2 pi tau), tau >= {p['tau_c']:g}; ln Pi = -lambda phi, lambda = "
            f"epsilon / phi_ref, phi_ref = {phi_ref:.3f}; channel capacity C = kappa F_ref, F_ref = {F_ref:.3f} (free arm, W1, "
            f"{p['n_cal']} calibration sprinklings); windows W1 {p['W1']}, W2 {p['W2']}; frame check beta' = {p['beta_check']}"]
    nse = float(p["n_se"])
    if mode == "search":
        arms = [("free", None, None)] + [(f"field_{float(e):g}", float(e), None) for e in p["epsilons"]] + \
               [(f"cap_{float(e):g}_{float(k):g}", float(e), float(k)) for e in p["epsilons"] for k in p["kappas"]]
        keys = [a[0] for a in arms]
        check = [k for k in keys if k.startswith("cap_")]
        seeds = [derive_seed(ctx.seed, stream, k) for k in range(n)]
        rows = parallel_map(run_task, [(p, sd, 0.0, arms, phi_ref, F_ref, check) for sd in seeds], ctx.workers, ctx.progress,
                            lambda i, r_, d, nn: f"[F3.9] search sprinkling {d}/{nn} ({time.time()-t0:.0f}s)")
        st = summarise(rows, keys, p, ctx, "search")
        lines = head + [f"search: {n} sprinklings, source at rest (~{np.mean([r_['N'] for r_ in rows]):.0f} events, "
                        f"{np.mean([r_['sources'] for r_ in rows]):.0f} source events); field_* = F3.8 replication (no capacity)"]
        lines += lines_for(st, keys, rows[0]["arms"]["free"]["series"]["t"])
        cands = []
        for key in check:
            x = st[key]
            if x["c1"] == "bound" and x["c2"] == "bound" and x["frame"]:
                _, e, k = key.split("_")
                cands.append((abs(x["s1"][0]) + nse * x["s1"][1], float(e), -float(k), key))
        cand = sorted(cands)[0][3] if cands else None
        lines.append(f"candidate (rule stated in the header): {cand if cand else 'none'}"
                     + (f"  -> confirm with: --set epsilon={cand.split('_')[1]} --set kappa={cand.split('_')[2]}" if cand else ""))
        return Result("RECORDED", lines, {"phi_ref": phi_ref, "F_ref": F_ref, "candidate": cand, "stats": _jsonable(st)})

    key = f"cap_{eps:g}_{kap:g}"
    arms = [("free", None, None), (key, eps, kap)]
    keys = ["free", key]
    bm = float(p["beta_move"])
    rows_r = parallel_map(run_task, [(p, derive_seed(ctx.seed, stream + "-rest", k), 0.0, arms, phi_ref, F_ref, [key])
                                     for k in range(n)],
                          ctx.workers, ctx.progress, lambda i, r_, d, nn: f"[F3.9] confirm rest {d}/{nn} ({time.time()-t0:.0f}s)")
    rows_m = parallel_map(run_task, [(p, derive_seed(ctx.seed, stream + "-move", k), bm, arms, phi_ref, F_ref, [key])
                                     for k in range(n)],
                          ctx.workers, ctx.progress, lambda i, r_, d, nn: f"[F3.9] confirm move {d}/{nn} ({time.time()-t0:.0f}s)")
    sr = summarise(rows_r, keys, p, ctx, "confirm-rest")
    sm = summarise(rows_m, keys, p, ctx, "confirm-move")
    lines = head + [f"confirm: candidate {key} (chosen from the search, CHOSEN FROM DATA); {n} sprinklings at rest, {n} moving "
                    f"(beta = {bm})", "  at rest:"]
    lines += lines_for(sr, keys, rows_r[0]["arms"]["free"]["series"]["t"])
    lines.append(f"  moving source (beta = {bm}):")
    lines += lines_for(sm, keys, rows_m[0]["arms"]["free"]["series"]["t"])
    v, se = sm[key]["vc"]
    tol = bm / 4.0
    lines.append(f"  inertia: centroid velocity of {key} in the source frame {v:+.4f} +- {se:.4f} vs beta/4 = {tol:.4f}")
    frame_ok = sr[key]["frame"] and sm[key]["frame"]
    c = sr[key]
    if not frame_ok:
        status = "INSTRUMENT_FAILED"
    elif c["c1"] in ("growing", "collapsing") or c["c2"] in ("growing", "collapsing") or \
            (np.isfinite(v) and abs(v) - nse * se > tol):
        status = "FAILED"
    elif c["c1"] == "bound" and c["c2"] == "bound" and sm[key]["c1"] == "bound" and np.isfinite(v) and abs(v) + nse * se <= tol:
        status = "PASSED"
    else:
        status = "UNDECIDED"
    lines.append(f"instrument: frame check (beta' = {p['beta_check']}) {'OK' if frame_ok else 'FAILS'}")
    return Result(status, lines, {"phi_ref": phi_ref, "F_ref": F_ref, "candidate": key, "rest": _jsonable(sr),
                                  "move": _jsonable(sm)})


def _jsonable(st):
    return {k: {kk: (list(vv) if isinstance(vv, (tuple, np.ndarray)) else vv) for kk, vv in v.items()} for k, v in st.items()}
