"""
F3.8 -- A field with range: does a flux bundle bind itself through its retarded field?
======================================================================================

Why this test exists (23 Oct 2026)
----------------------------------
F3.7 search (no candidate): bending towards slower clocks, with the clocks
taken from the LOCAL load of each event, pulls the bundle in at first
(R 3.0 -> 2.5) and then heats it to the walls: that field is rough and
short-ranged (it exists only inside the bundle) and its gradient is mostly
the noise of the discrete load. Gravity is not like that: the field is
smooth and long-ranged because the source acts through a field equation
and every event feels the sum of many sources. New hypothesis after that
failure (with its own criteria; nothing of F3.7 is changed): the bundle
can bind itself only through a field with range, produced by the flux via
the retarded Green function of the wave operator, i.e. a causal sum over
the flux in the past of every event. This joins the two branches of the
laboratory: matter with momentum (F3.6) and the field (A5', BD operator).

Model (locked before any run; protocol of the author as in F3.7)
----------------------------------------------------------------
Substrate, source, transport and readout as F3.7 (2+1 Poisson
sprinkling, rho = 1, lengths in l; source tube radius 4, time [0, 4),
u_src; carried direction; each event sends all its flux to its child of
smallest rapidity relative to its direction).
Field: the retarded Green function of the 2+1 wave operator,
G(x, w) = 1 / (2 pi tau(x, w)) for w in the causal past of x (tau the
proper time; the Green function of box phi = -source in 2+1), summed over
the flux of the past events:
    phi(x) = sum_{w < x, tau >= tau_c} F_w / (2 pi tau(x, w))     (rho = 1: each event carries volume 1),
tau_c = 1 (the discreteness length; below it the continuum Green function
has no meaning for a causal set; it also removes the integrable 1/tau peak
on the light cone, whose discrete sum would be noise). Every event has a
field (events without flux too); it uses only the flux of its past, which
the causal pass has finished.
Clock: ln Pi = -lambda phi (weak-field form of the slow clock; Pi < 1 where
the field is strong). Strength: lambda = epsilon / phi_ref, phi_ref = the
median field at the carrying events of the free arm in W1 (calibration
sprinklings, own stream, field computed without back-reaction); epsilon in
{0.01, 0.1, 0.5}: the depth of ln Pi at a typical carrier, weak / medium /
strong. Fixed now.
Bending (as F3.7 lens, with this clock): transverse gradient of ln Pi at v
by least squares over v's parents with rest-frame lag <= L_g = 2 (at least
4), projected orthogonal to the carried direction; the direction sent
along the step (proper time tau_k) is boosted by |G| tau_k towards -G.
Arms: free (no field acts; reference) and field_epsilon x 3.
Covariance (criterion (c)) is checked bit for bit in a frame boosted by
beta' = 0.6; as in F3.6/F3.7 the rules use the metric at the scale of a
link and the choice of the child uses the future (criteria (a), (b)).
Remark declared now: in 2+1 the Newtonian potential of a source is
logarithmic and confines every bound orbit; binding in 2+1 is therefore
weaker evidence than in 3+1 (1/r, with an escape velocity). A PASSED here
asks for the 3+1 version before any physical reading.

Classification, parts, candidate rule and confirm criteria: exactly those
of F3.7 (bound |s| + 2 SE <= s_free/4 in W1 [12, 28] and W2 [28, 44],
growing, collapsing, unresolved, wall-limited; search RECORDED with at
most one candidate: bound in W1 and W2 and frame check, smallest
|s_W1| + 2 SE, ties the weaker epsilon; confirm on new streams with the
moving source: PASSED bound at rest in W1 and W2, frame check, inertia
|v_c| + 2 SE <= beta/4 and bound in motion; FAILED growing or collapsing,
or |v_c| - 2 SE > beta/4; INSTRUMENT_FAILED frame check; else UNDECIDED).

Seen before this test was written: F3.6, F3.7 search (above). The field
arms had not been run. Quick runs of this file (after everything above
was written; 3 sprinklings, box 24 x 16 x 16, source radius 2; decide
nothing): free grows 0.07 per unit time; field epsilon 0.01 like free;
epsilon 0.1 keeps R_rms ~2.0-2.5 (slope +0.009 +- 0.031 in W1, -0.025 +-
0.038 in W2, unresolved with 3 sprinklings); epsilon 0.5 noisy (bending
0.075 rapidity per step). Nothing was changed after them.
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
    "id": "F3.8",
    "phase": "F3",
    "title": "Pole s dosahem: sváže svazek toku sám sebe přes své retardované pole?",
    "theory_part": "hmota budí pole s dosahem (A5′, retardovaná Greenova funkce), pole zpomaluje hodiny, dráhy se ohýbají "
                   "k pomalejším hodinám; vázanost hmoty; spojení větve hmoty (F3.6) a pole",
    "question": "Zastaví se rozšiřování svazku s neseným směrem, když hodiny určuje hladké pole φ = Σ F/(2πτ) přes tok "
                "v minulosti každé události a směr se ohýbá podle příčného gradientu ln Π = −λφ? Zůstane pravidlo "
                "nezávislé na soustavě a svazek setrvačný?",
    "method": "F3.7 lens s hodinami z retardovaného pole místo místní zátěže; λ = ε/φ_ref, ε ∈ {0,01; 0,1; 0,5}, φ_ref "
              "z kalibračních sprinklingů; protokol search → nejvýš jeden kandidát → confirm na nových sprinklinzích "
              "včetně pohyblivého zdroje",
    "prediction": "podle hypotézy po F3.7: hladké pole s dosahem svazek sváže (ve 2+1 logaritmický potenciál); proti "
                  "tomu stojí zahřívání šumem odhadu gradientu",
    "criterion": "jako F3.7: vázaný |s| + 2 SE <= s_free/4 v oknech [12, 28] a [28, 44]; search: kandidát vázaný v obou "
                 "oknech se shodnou volbou v soustavě β′ = 0,6, nejmenší |s| + 2 SE (slabší ε při shodě); confirm: PROŠEL "
                 "vázaný v obou oknech, kontrola soustavy, |v_c| + 2 SE <= β/4 a vázaný v pohybu; NEPROŠEL rostoucí nebo "
                 "kolaps, nebo |v_c| − 2 SE > β/4; MĚŘIDLO kontrola soustavy; jinak NEROZHODNUTO. Úspěch ve 2+1 "
                 "vyžaduje verzi 3+1 před fyzikálním čtením.",
    "null_model": "free (F3.6 bez pole): přirozená difuze svazku",
    "level": "mechanism",
    "meters": ["K.6"],
    "depends_on": ["F3.6", "F3.7"],
    "params": {"D": 3, "T": 48.0, "half": 16.0, "r_src": 4.0, "t_src": 4.0, "slab": 2.0, "W1": [12.0, 28.0],
               "W2": [28.0, 44.0], "win_move": [5.0, 28.0], "beta_check": 0.6, "beta_move": 0.5, "tau_c": 1.0,
               "epsilons": [0.01, 0.1, 0.5], "L_g": 2.0, "min_parents": 4, "n_cal": 8, "n_seeds": 32, "n_se": 2.0,
               "boot_B": 4000, "epsilon": None},
    "parts": {"search": {"mode": "search"}, "confirm": {"mode": "confirm"}},
    "settable": {"confirm": ["epsilon"]},
    "quick": {"n_seeds": 3, "n_cal": 1, "boot_B": 200, "T": 24.0, "half": 8.0, "r_src": 2.0, "t_src": 2.0,
              "W1": [6.0, 14.0], "W2": [14.0, 22.0], "win_move": [3.0, 14.0], "slab": 1.0},
    "estimate": "search ~15–30 min s 10 procesy (jeden sprinkling ~2,5 min) (4 větve + přepočty v jiné soustavě; pole = součet přes tok v minulosti); "
                "confirm podobně",
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


def transport(X, ptr, idx, cptr, cidx, s, u_src, arm, lam, p, want_phi=False):
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
            inflow[k] += tot
            U[k] += tot * send
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
    """Field at the carrying events of the free arm in W1 (no back-reaction): the reference field."""
    X, ptr, idx, cptr, cidx, Xs, s, u_src = setup(p, seed, 0.0)
    r = transport(X, ptr, idx, cptr, cidx, s, u_src, "free", 0.0, p, want_phi=True)
    w1 = p["W1"]
    m = (r["flux"] > 0) & (Xs[:, 0] >= w1[0]) & (Xs[:, 0] <= w1[1]) & np.isfinite(r["phi"])
    return r["phi"][m]


def run_task(p, seed, beta, arms, phi_ref):
    """One sprinkling: free and field arms, readout in the source frame, frame check (beta')."""
    X, ptr, idx, cptr, cidx, Xs, s, u_src = setup(p, seed, beta)
    bc = float(p["beta_check"])
    Xb = boost_points(X, bc)
    ub = boost_points(u_src[None, :], bc)[0]
    out = {"N": len(X), "sources": int(s.sum()), "arms": {}}
    for arm, eps in arms:
        lam = float(eps) / phi_ref if eps is not None else 0.0
        r = transport(X, ptr, idx, cptr, cidx, s, u_src, arm, lam, p)
        rb = transport(Xb, ptr, idx, cptr, cidx, s, ub, arm, lam, p)
        ser = readout(Xs, r["flux"], p, beta)
        key = arm if eps is None else f"{arm}_{eps:g}"
        out["arms"][key] = {"series": ser, "s1": slope(ser["t"], ser["R"], p["W1"]), "s2": slope(ser["t"], ser["R"], p["W2"]),
                            "vc": slope(ser["t"], ser["xc"], p["win_move"]),
                            "frame_identical": bool(np.array_equal(r["decision"], rb["decision"])),
                            "bend_mean": float(np.nanmean(r["bend"])) if np.isfinite(r["bend"]).any() else float("nan"),
                            "Pi_min": float("nan")}
    return out

# ============================================================================================== statistics

def classify(sv, se, sfree, nse):
    if not all(np.isfinite([sv, se, sfree])):
        return "unresolved"
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
        eps = p.get("epsilon")
        if eps is None:
            raise SystemExit("F3.8 confirm: give the candidate of the search, e.g. --set epsilon=0.1")
        eps = float(eps)
        if eps not in [float(e) for e in p["epsilons"]]:
            raise SystemExit(f"F3.8 confirm: epsilon must be one of the grid {p['epsilons']} (got {eps})")
    stream = "F3.8-" + mode
    cal = parallel_map(calibration_task, [(p, derive_seed(ctx.seed, stream + "-cal", k)) for k in range(int(p["n_cal"]))],
                       ctx.workers, ctx.progress, lambda i, r_, d, nn: f"[F3.8] calibration {d}/{nn} ({time.time()-t0:.0f}s)")
    phi_ref = float(np.median(np.concatenate(cal)))
    head = [f"2+1 sprinkling rho = 1, box t [0, {p['T']:g}], half width {p['half']:g}; source radius {p['r_src']:g}, time "
            f"[0, {p['t_src']:g}); field phi = sum F / (2 pi tau) over the past flux, tau >= {p['tau_c']:g}; ln Pi = -lambda phi, "
            f"lambda = epsilon / phi_ref, phi_ref = {phi_ref:.3f} (median field at the carriers of the free arm in W1, "
            f"{p['n_cal']} calibration sprinklings); gradient from parents with rest-frame lag <= {p['L_g']:g}; windows "
            f"W1 {p['W1']}, W2 {p['W2']}; frame check beta' = {p['beta_check']}"]
    if mode == "search":
        arms = [("free", None)] + [("field", float(e)) for e in p["epsilons"]]
        keys = ["free"] + [f"field_{float(e):g}" for e in p["epsilons"]]
        seeds = [derive_seed(ctx.seed, stream, k) for k in range(n)]
        rows = parallel_map(run_task, [(p, sd, 0.0, arms, phi_ref) for sd in seeds], ctx.workers, ctx.progress,
                            lambda i, r_, d, nn: f"[F3.8] search sprinkling {d}/{nn} ({time.time()-t0:.0f}s)")
        st = summarise(rows, keys, p, ctx, "search")
        lines = head + [f"search: {n} sprinklings, source at rest (~{np.mean([r_['N'] for r_ in rows]):.0f} events, "
                        f"{np.mean([r_['sources'] for r_ in rows]):.0f} source events)"]
        lines += lines_for(st, keys, rows[0]["arms"]["free"]["series"]["t"])
        nse = float(p["n_se"])
        cands = []
        for key in keys[1:]:
            x = st[key]
            if x["c1"] == "bound" and x["c2"] == "bound" and x["frame"]:
                cands.append((abs(x["s1"][0]) + nse * x["s1"][1], float(key.split("_")[1]), key))
        cand = sorted(cands)[0][2] if cands else None
        lines.append(f"candidate (rule stated in the header): {cand if cand else 'none'}"
                     + (f"  -> confirm with: --set epsilon={cand.split('_')[1]}" if cand else ""))
        return Result("RECORDED", lines, {"phi_ref": phi_ref, "candidate": cand, "stats": _jsonable(st)})

    arms = [("free", None), ("field", eps)]
    key = f"field_{eps:g}"
    keys = ["free", key]
    bm = float(p["beta_move"])
    rows_r = parallel_map(run_task, [(p, derive_seed(ctx.seed, stream + "-rest", k), 0.0, arms, phi_ref) for k in range(n)],
                          ctx.workers, ctx.progress, lambda i, r_, d, nn: f"[F3.8] confirm rest {d}/{nn} ({time.time()-t0:.0f}s)")
    rows_m = parallel_map(run_task, [(p, derive_seed(ctx.seed, stream + "-move", k), bm, arms, phi_ref) for k in range(n)],
                          ctx.workers, ctx.progress, lambda i, r_, d, nn: f"[F3.8] confirm move {d}/{nn} ({time.time()-t0:.0f}s)")
    sr = summarise(rows_r, keys, p, ctx, "confirm-rest")
    sm = summarise(rows_m, keys, p, ctx, "confirm-move")
    nse = float(p["n_se"])
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
    return Result(status, lines, {"phi_ref": phi_ref, "candidate": key, "rest": _jsonable(sr), "move": _jsonable(sm)})


def _jsonable(st):
    return {k: {kk: (list(vv) if isinstance(vv, (tuple, np.ndarray)) else vv) for kk, vv in v.items()} for k, v in st.items()}
