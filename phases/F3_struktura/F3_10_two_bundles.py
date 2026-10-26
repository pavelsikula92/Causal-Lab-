"""
F3.10 -- Two bundles: do separate flux bundles attract each other through their retarded fields?
================================================================================================

Why this test exists (23 Oct 2026, the author's design)
-------------------------------------------------------
F3.7-F3.9 searched for a self-bound bundle and found none; each step added
a mechanism after a failure (F3.9 would have been followed by a fourth,
"pressure with momentum"). The author stopped that: no further mechanism;
the question changes. F3.6 gave localization and inertia of a bundle that
carries its direction (momentum) as a state; F3.8 and its replication in
F3.9 gave a reproducible response of a bundle to its OWN retarded field.
Mutual gravity, the force law and the universality of free fall are still
hypotheses. Two separate bundles test attraction without asking gravity to
stabilise a particle. No new mechanism enters: the dynamics is exactly the
F3.8 field arm.

Model (locked before any run)
-----------------------------
2+1 Poisson sprinkling, rho = 1 (lengths in l), links (core.order.links).
Two source tubes ("bundles") at rest, centred at x = -d/2 (left, L) and
x = +d/2 (right, R), radius 3, source time [0, 4); every source event of a
bundle emits s = f (its flux factor) with the 4-velocity (1, 0, 0). Each
bundle is a separate SPECIES of flux: it has its own carried direction and
its own children choice (F3.6 carried rule) and its own retarded field
    phi_i(x) = sum over the species-i flux events w in the causal past of x,
               tau >= 1, of F_w / (2 pi tau)          (F3.8).
The species interact ONLY through the field: a species responds to the
sum of the fields of the species in its response set (ln Pi = -lambda
sum phi_j; bending of its carried direction by the transverse gradient of
ln Pi estimated from the parents, rest-frame lag <= 2, >= 4 parents, as
F3.8). lambda = epsilon / phi_ref, phi_ref = the median own field at the
carriers of a single bundle (f = 1) in the window [12, 28] on calibration
sprinklings (own stream). epsilon = 0.1, CHOSEN FROM DATA of F3.8 / F3.9:
the strength at which the self-field visibly acts on a bundle without the
collapse-and-rebound of epsilon = 0.5; locked here, not tuned in F3.10.
Branches on the SAME sprinkling and sources:
  A  each bundle alone (only its own field): the reference drift;
  B  both bundles, each responds to the sum of both fields (decides);
  C  both bundles, each responds only to its own field: must equal A
     choice for choice (instrument: the only coupling is the field);
  D  the mirror: B on the sprinkling reflected x -> -x with the reflected
     sources (positions and identities swap; the causal order and links are
     unchanged). The rule uses only Minkowski products, so the reflection is
     an exact symmetry: every choice must be identical; any difference is a
     preferred coordinate direction or tie-breaking of the instrument
     (instrument, every part; noticed in the first quick run, which gave a
     difference of exactly 0 for a statistical version of this check);
  frame check: B recomputed with coordinates and source velocity boosted
     by beta' = 0.6 must choose identically (instrument).
Readout (per slab of time 2, flux-weighted over the carriers of a
species): centroid x_i(t); x-momentum P_i(t) = sum F u_x. Relative
separation correction delta r(t) = [x_R - x_L](B) - [x_R - x_L](A);
acceleration a = the coefficient of (t - t_on)^2 / 2 in a least-squares
quadratic over the window [d + 2, T - 4] (the field of the partner needs
the time d to arrive), t_on = d + 2; the same fit for delta x_L and for
delta P_i (slope, linear fit).
Diagnostic without a verdict: the acceleration from the momentum,
a_P,i = (dP_i/dt) / M_i with M_i the mean flux of species i per slab in
the window (alone branch); if the bundle moved as a rigid body carrying
its momentum, a_P = the centroid acceleration. A difference means that
the centroid moves also by deformation (one side of the bundle bent more
than the other), not only by a change of the carried momentum.
Prediction (independent of the data that test it; a QUASI-STATIC POINT-LINE
APPROXIMATION, see the QA review below): two point world lines
under the continuum retarded field of a static line source of strength mu
(flux per unit time) in 2+1 with the same cut tau >= 1:
    a_i = lambda mu_j G(t, r),   G = (1/2 pi) [ t / (r sqrt(t^2 - r^2)) - 1 / (r^2 sqrt(1 + 1/r^2)) ]
for t >= sqrt(r^2 + 1) (0 before), towards the partner; mu from the
calibration (flux per unit time of a single bundle with f = 1, times f);
integrated from rest and fitted with the same quadratic over the same
window: a_pred. For t >> r this is the static 1/r force of 2+1; at the
separations used the retarded build-up still matters, so the law is
tested against the full retarded G, not against 1/r alone.

Parts (each hypothesis its own verdict; replication with --seed 1 required)
--------------------------------------------------------------------------
attraction  d = 12, f = 1, 1. INSTRUMENT_FAILED if C != A, the frame
            check or the mirror check fails. PASSED if a + 2 SE < 0
            (the bundles accelerate towards each other); FAILED if
            a - 2 SE >= 0; UNDECIDED otherwise. Reported: a / a_pred.
law         d = 8, 12, 16, 20, f = 1, 1. Exponent k of |a| against d
            (weighted least squares of ln|a| on ln d; bootstrap SE over
            sprinklings) against k_pred of the predicted a_pred(d) through
            the same fit; alternatives k_pred +- 1 (a force law one power
            off). PASSED if every d attracts (a + 2 SE < 0) and
            |k - k_pred| + 2 SE <= 1/4; FAILED if |k - k_pred| - 2 SE > 1/4;
            UNDECIDED otherwise. Reported: a / a_pred per d.
momentum    d = 12, f = 1 (L), 2 (R). The model has no field momentum, so
            the balance must hold between the bundles: dP_L/dt + dP_R/dt
            (each relative to its alone reference). Both changes must be
            resolved (|dP_i/dt| - 2 SE > 0), else UNDECIDED. PASSED if
            |dP_L + dP_R| + 2 SE <= (|dP_L| + |dP_R|) / 8 (a quarter of the
            mean single change); FAILED if |dP_L + dP_R| - 2 SE > that;
            UNDECIDED otherwise. If it fails, the theory does not define a
            conserved total momentum (a defect of the theory, not of the test).
universality d = 12, two configurations on the same sprinklings: L with
            f = 1 and with f = 2, R with f = 2 in both. (Changed by the QA
            review before this part was run: R is an EXTERNAL source, it
            responds only to its own field, so the field acting on L is the
            same field in both configurations, checked bit for bit; the
            acceleration is the one from the carried momentum a_P,L =
            (dP_L/dt) / M_L, in the linear-response window.) Acceleration of L
            towards R (+x, positive), a_P,L(1) and a_P,L(2), each resolved
            (a_P,L - 2 SE > 0) else UNDECIDED. Alternative: a_L proportional to
            1 / f_L (twice as large for f = 1): the difference is then
            (2/3) |a_mean|; a quarter of that distance: PASSED if
            |a_L(1) - a_L(2)| + 2 SE <= |a_mean| / 6; FAILED if
            |a_L(1) - a_L(2)| - 2 SE > |a_mean| / 6. Declared now: in this
            model the bending does not depend on the carrier's own flux, so
            universality of the response is built into the rule; the test
            checks that the self-field (stronger for f = 2), merging and the
            discreteness do not spoil it. It is a consistency check, not an
            independent confirmation of the equivalence principle.
Predefined failure (the author): attraction that comes only with a large
loss of total momentum (momentum FAILED) or with frame dependence (frame
check) is not a confirmation of the gravitational mechanism.

QA review (the author, 23 Oct, after attraction and law had run at seed 0)
------------------------------------------------------------------------
Checked in the code and recorded; nothing of it changes a verdict that has
been reached.
- The quasi-static point-line prediction is not valid at small d: with the
  calibrated lambda and mu its separation collapses inside the window (at
  d = 8 the point lines pass through each other, r_min -4.7; at d = 12 it
  falls to 5.5). It uses the current separation in the static-source
  formula, not the real past of a moving source. Its ratios are reported
  as an approximation only; the law part is superseded by F3.10b (a
  first-order prediction from the actual discrete field of the actual
  partner, in a linear-response window) and is not replicated.
- momentum and universality (not yet run) use the linear-response window
  [1.25 d, 2.5 d]: the displacement stays a small fraction of d and the
  retarded build-up has the same shape in units of d. universality has an
  external source (above): with a responding partner, f_L = 2 would pull
  R closer and enlarge the force on L, a trajectory effect, not a response.
- The processing order is asserted to be a linear extension of the causal
  order (every parent before its child); a boost keeps the causal order, so
  the same order serves every frame, and the retarded field at an event
  uses exactly the finished events of its causal past (tau >= 1, a Lorentz-
  invariant set).
- The checks compare more than the choices: C against A bit for bit
  (choices, flux, directions); the mirror and the frame check also the
  flux (bit for bit) and the directions mapped back (1e-9). A difference
  whose first event is a numerical tie (margin of the choice < 1e-9) is
  reported as a tie, not as a frame effect.
- part unit: deterministic checks on small sprinklings: lambda = 0 gives the
  separate bundles bit for bit; the scalar flux is conserved (sources = flux
  at the events without children); a bundle responding only to the
  partner's field is kicked towards it; C, mirror and frame as above.
- Theory, reported without a verdict: the model conserves the scalar flux,
  not the 4-momentum. At a merger of different directions the incoming
  4-vector has |P| > F and the outgoing one has norm F: the excess is
  dropped (measured: ~2-3e-4 of the flux, ~0.03-0.05 % of the x-momentum per
  event, small against the effects). The field carries no momentum. The
  momentum part tests whether the mutual kicks balance anyway.
- lambda is the model's constant, used identically in the simulation and in
  any prediction, so its calibration noise does not enter a comparison;
  mu enters only the approximation.

Remarks declared now: 2+1 only (logarithmic potential); the 3+1 version
(a ~ 1/r^2 from the 3+1 Green function) follows only after this one. As in
F3.6-F3.9, the rules use the metric at the scale of a link and the choice
of the child uses the future of v (criteria (a), (b) of the axiom audit).

Seen before this test was written: F3.6-F3.9. Nothing of F3.10 has been run.
"""

from __future__ import annotations

import hashlib
import math
import time

import numpy as np

from core.order.links import links_csr
from core.order.sprinkle import Box, derive_seed, sprinkle_poisson
from core.stats.parallel import parallel_map
from core.stats.stats import boot_se
from core.testing import Result

TEST = {
    "id": "F3.10",
    "phase": "F3",
    "title": "Dva svazky: přitahují se oddělené svazky toku přes svá retardovaná pole?",
    "theory_part": "vzájemná gravitace hmoty s hybností (F3.6) přes pole s dosahem (F3.8); zákon síly, zachování hybnosti, "
                   "univerzalita volného pádu; bez nového mechanismu",
    "question": "Zrychlují dva oddělené svazky k sobě (proti referenci, kde každý cítí jen své pole), podle zákona "
                "retardovaného pole 2+1 přes několik vzdáleností, se zachováním celkové hybnosti a se zrychlením nezávislým "
                "na vlastním toku?",
    "method": "dva druhy toku (každý svazek svůj nesený směr, volbu potomka a pole φ_i = Σ F/(2πτ)); interakce jen polem; "
              "větve A (sám), B (společné pole), C (jen vlastní pole, musí = A), D (zrcadlo x → −x), kontrola soustavy; "
              "ε = 0,1 z F3.8/F3.9 zamčené; předpověď z retardovaného pole statické úsečky ve 2+1 se stejným řezem",
    "prediction": "přitahování se zrychlením podle retardovaného zákona 2+1 (k_pred z téhož fitu), bilance hybnosti, "
                  "univerzalita (zabudovaná do pravidla, kontrola konzistence)",
    "criterion": "attraction: PROŠEL a + 2 SE < 0, NEPROŠEL a − 2 SE >= 0. law: PROŠEL všechna d přitahují a "
                 "|k − k_pred| + 2 SE <= 1/4, NEPROŠEL |k − k_pred| − 2 SE > 1/4. momentum: PROŠEL |Ṗ_L + Ṗ_R| + 2 SE <= "
                 "(|Ṗ_L| + |Ṗ_R|)/8, NEPROŠEL − 2 SE > tolerance. universality: PROŠEL |a_L(1) − a_L(2)| + 2 SE <= |ā|/6, "
                 "NEPROŠEL − 2 SE > |ā|/6. MĚŘIDLO: C ≠ A, kontrola soustavy, zrcadlo (volby B na odraženém sprinklingu ≠ volby B). "
                 "Jinak NEROZHODNUTO. unit: deterministické kontroly implementace (λ = 0, bilance toku, znaménko kopnutí, "
                 "C/zrcadlo/soustava i pro tok a směry), PROŠEL / MĚŘIDLO. Replikace --seed 1. law nahrazena F3.10b "
                 "(QA 23. 10.: kvazistatická předpověď neplatí při d = 8 a 12).",
    "null_model": "větev A (každý svazek jen se svým polem) a C (obě, bez vzájemné vazby)",
    "level": "mechanism",
    "meters": ["K.6"],
    "depends_on": ["F3.6", "F3.8", "F3.9"],
    "params": {"D": 3, "T": 48.0, "margin": 12.0, "r_src": 3.0, "t_src": 4.0, "slab": 2.0, "tau_c": 1.0,
               "epsilon": 0.1, "L_g": 2.0, "min_parents": 4, "beta_check": 0.6, "cal_win": [12.0, 28.0],
               "n_cal": 8, "n_seeds": 32, "n_se": 2.0, "boot_B": 4000, "d_mass": 12.0},
    "parts": {"unit": {"mode": "unit", "n_seeds": 4, "n_cal": 1, "T": 24.0, "margin": 6.0, "r_src": 2.0, "t_src": 2.0,
                       "cal_win": [6.0, 14.0], "d_unit": 8.0},
              "attraction": {"mode": "attraction", "seps": [12.0]},
              "law": {"mode": "law", "seps": [8.0, 12.0, 16.0, 20.0], "n_seeds": 24},
              "momentum": {"mode": "momentum", "window": "lr"},
              "universality": {"mode": "universality", "window": "lr", "external_R": True}},
    "quick": {"n_seeds": 3, "n_cal": 1, "boot_B": 200, "T": 24.0, "margin": 6.0, "r_src": 2.0, "t_src": 2.0,
              "cal_win": [6.0, 14.0], "d_mass": 8.0},
    "quick_parts": {"unit": {"n_seeds": 2}, "attraction": {"seps": [8.0]}, "law": {"seps": [6.0, 7.0, 8.0, 9.0]}},
    "estimate": "unit < 1 min, attraction ~5 min, law ~15 min (nahrazena F3.10b), momentum ~5 min, universality ~10 min "
                "(10 procesů; jeden sprinkling 70–100 s na jednom jádře: 6 průchodů s polem, kontroly C, soustavy a zrcadla "
                "na každém)",
}


# ============================================================================================== geometry

def boost_points(X, beta):
    if beta == 0.0:
        return np.array(X, float, copy=True)
    g = 1.0 / math.sqrt(1.0 - beta * beta)
    Y = np.array(X, float, copy=True)
    Y[..., 0] = g * (X[..., 0] - beta * X[..., 1])
    Y[..., 1] = g * (X[..., 1] - beta * X[..., 0])
    return Y


def mink(a, b):
    return a[..., 0] * b[..., 0] - (a[..., 1:] * b[..., 1:]).sum(axis=-1)


def children_csr(ptr, idx, N):
    child = np.repeat(np.arange(N), np.diff(ptr))
    o = np.argsort(idx, kind="stable")
    cptr = np.concatenate([[0], np.cumsum(np.bincount(idx, minlength=N))]).astype(np.int64)
    return cptr, child[o]


def rest_basis(u):
    D = u.size
    out = []
    for k in range(1, D):
        e = np.zeros(D)
        e[k] = 1.0
        e = e - mink(e, u) * u
        for f in out:
            e = e + mink(e, f) * f
        out.append(e / math.sqrt(max(-float(mink(e, e)), 0.0)))
    return np.array(out)


# ============================================================================================== field and transport

class RetardedField:
    """phi(x) = sum over the finished flux events w of one species in the causal past of x, tau >= tau_c, of F_w / (2 pi tau)."""

    def __init__(self, X, tau_c):
        self.X = X
        self.tau_c2 = float(tau_c) ** 2
        self.pos, self.flx = [], []
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
        self.cache[x_idx] = val                         # asked only at events whose whole past is processed
        return val


def transport_species(X, ptr, idx, cptr, cidx, sources, resp, u_src, lam, p):
    """
    One causal pass with several species of flux (module notes). sources: list of source arrays (one per species);
    resp: list of lists, the species whose fields species i responds to. Returns per species: flux, direction,
    decision (chosen child).
    """
    N, D = X.shape
    S = len(sources)
    L_g, nmin = float(p["L_g"]), int(p["min_parents"])
    fields = [RetardedField(X, float(p["tau_c"])) for _ in range(S)]
    inflow = [np.zeros(N) for _ in range(S)]
    U = [np.zeros((N, D)) for _ in range(S)]
    flux = [np.zeros(N) for _ in range(S)]
    dirs = [np.zeros((N, D)) for _ in range(S)]
    dec = [np.full(N, -1, dtype=np.int64) for _ in range(S)]
    gap = [np.full(N, np.inf) for _ in range(S)]        # margin between the best and the second-best child (ties)
    kick = [np.zeros(N) for _ in range(S)]              # x-momentum given by the bending at v: tot (send - u)_x
    normP = [np.zeros(N) for _ in range(S)]             # Minkowski norm of the incoming 4-vector (>= tot; > at mergers)
    for v in range(N):
        active = [i for i in range(S) if inflow[i][v] + sources[i][v] > 0.0]
        if not active:
            continue
        a, b = cptr[v], cptr[v + 1]
        ks = cidx[a:b]
        B = X[ks] - X[v] if b > a else None
        tau = np.sqrt(np.maximum(mink(B, B), 1e-300)) if b > a else None
        par = idx[ptr[v]:ptr[v + 1]]
        for i in active:
            tot = inflow[i][v] + sources[i][v]
            P = U[i][v] + sources[i][v] * u_src
            u = P / math.sqrt(max(float(mink(P, P)), 1e-300))
            flux[i][v] = tot
            dirs[i][v] = u
            normP[i][v] = math.sqrt(max(float(mink(P, P)), 0.0))
            if b == a:
                continue
            ch = mink(u[None, :], B) / tau
            j = int(np.argmin(ch))
            k = int(ks[j])
            dec[i][v] = k
            if ch.size > 1:
                gap[i][v] = float(np.partition(ch, 1)[1] - ch[j])
            send = u
            if lam > 0.0 and par.size:
                d = X[par] - X[v]
                lag = -mink(d, u[None, :])
                sel = lag <= L_g
                if int(sel.sum()) >= nmin:
                    pars = par[sel]
                    dd = d[sel]
                    pp = dd - mink(dd, u[None, :])[:, None] * u[None, :]
                    E = rest_basis(u)
                    q = -np.stack([mink(pp, e[None, :]) for e in E], axis=1)
                    lnPi = np.array([-lam * sum(fields[jj].value(int(w)) for jj in resp[i]) for w in pars])
                    A = np.column_stack([np.ones(q.shape[0]), q])
                    coef, *_ = np.linalg.lstsq(A, lnPi, rcond=None)
                    G = coef[1:]
                    gn = float(np.sqrt((G * G).sum()))
                    if gn > 0.0:
                        e_dir = -(G[0] * E[0] + G[1] * E[1]) / gn
                        eta = gn * float(tau[j])
                        send = math.cosh(eta) * u + math.sinh(eta) * e_dir
            kick[i][v] = tot * float(send[1] - u[1])
            inflow[i][k] += tot
            U[i][k] += tot * send
        for i in active:
            fields[i].add(v, flux[i][v])
    return {"flux": flux, "dir": dirs, "dec": dec, "gap": gap, "kick": kick, "normP": normP}


# ============================================================================================== readout and fits

def species_series(X, flux, dirs, p):
    """Per slab: flux-weighted centroid x, x-momentum sum F u_x, total flux."""
    T, sl = float(p["T"]), float(p["slab"])
    edges = np.arange(0.0, T - sl + 1e-9, sl)
    out = {"t": [], "x": [], "P": [], "F": []}
    car = flux > 0.0
    for a, b in zip(edges[:-1], edges[1:]):
        m = car & (X[:, 0] >= a) & (X[:, 0] < b)
        w = flux[m]
        out["t"].append(0.5 * (a + b))
        out["F"].append(float(w.sum()))
        out["P"].append(float((w * dirs[m, 1]).sum()))
        out["x"].append(float((w * X[m, 1]).sum() / w.sum()) if w.sum() > 0 else float("nan"))
    return out


def quad_acc(t, y, t_on, t_end):
    """Coefficient a of y = c0 + c1 (t - t_on) + a (t - t_on)^2 / 2 over [t_on, t_end] (finite points)."""
    t, y = np.asarray(t, float), np.asarray(y, float)
    m = (t >= t_on) & (t <= t_end) & np.isfinite(y)
    if m.sum() < 4:
        return float("nan")
    s = t[m] - t_on
    A = np.column_stack([np.ones(s.size), s, 0.5 * s * s])
    return float(np.linalg.lstsq(A, y[m], rcond=None)[0][2])


def lin_slope(t, y, t_on, t_end):
    t, y = np.asarray(t, float), np.asarray(y, float)
    m = (t >= t_on) & (t <= t_end) & np.isfinite(y)
    return float(np.polyfit(t[m], y[m], 1)[0]) if m.sum() >= 2 else float("nan")


def G_ret(t, r, tau_c):
    """Magnitude of the gradient of the retarded field of a unit static line source started at t = 0 (2+1, cut tau >= tau_c)."""
    if r <= 0.0 or t * t < r * r + tau_c * tau_c:
        return 0.0
    c = math.sqrt(1.0 + (tau_c / r) ** 2)
    return (t / (r * math.sqrt(t * t - r * r)) - tau_c / (r * r * c)) / (2.0 * math.pi)


def predict(d, mu_L, mu_R, lam, p, dt=0.01):
    """Two point world lines from rest at -d/2, +d/2 under the mutual retarded fields: x_L(t), x_R(t) at the slab centres."""
    T, sl, tc = float(p["T"]), float(p["slab"]), float(p["tau_c"])
    xL, xR, vL, vR = -0.5 * d, 0.5 * d, 0.0, 0.0
    ts = np.arange(0.0, T - sl + 1e-9, sl)
    centres = 0.5 * (ts[:-1] + ts[1:])
    out_L, out_R = [], []
    t, ci = 0.0, 0
    while ci < centres.size:
        if t >= centres[ci]:
            out_L.append(xL)
            out_R.append(xR)
            ci += 1
            continue
        r = xR - xL
        g = G_ret(t, r, tc)
        aL, aR = lam * mu_R * g, -lam * mu_L * g
        vL += aL * dt
        vR += aR * dt
        xL += vL * dt
        xR += vR * dt
        t += dt
    return centres, np.array(out_L), np.array(out_R)


# ============================================================================================== tasks

def make(p, seed, d):
    """Sprinkling of the box for separation d, links, children."""
    T, M = float(p["T"]), float(p["margin"])
    rng = np.random.default_rng(seed)
    half = 0.5 * d + M
    X = sprinkle_poisson(Box(int(p["D"]), [0.0, -half, -M], [T, half, M]), 1.0, rng)
    ptr, idx = links_csr(X, 1.0)
    # the pass `for v in range(N)` needs every parent before its child; the points are time-sorted, and a
    # boost keeps the causal order, so the same index order stays a linear extension in every frame
    child = np.repeat(np.arange(len(X)), np.diff(ptr))
    assert bool((idx < child).all()), "processing order is not a linear extension of the causal order"
    cptr, cidx = children_csr(ptr, idx, len(X))
    return X, ptr, idx, cptr, cidx


def tube(X, x0, p, f):
    """Source of a bundle centred at x = x0: f on every event of the tube (radius r_src, time [0, t_src))."""
    r2 = (X[:, 1] - x0) ** 2 + X[:, 2] ** 2
    return f * ((X[:, 0] < float(p["t_src"])) & (r2 <= float(p["r_src"]) ** 2)).astype(float)


def calibration_task(p, seed):
    """A single bundle (f = 1): its own field at its carriers in cal_win and its flux per unit time."""
    X, ptr, idx, cptr, cidx = make(p, seed, 0.0)
    u0 = np.array([1.0, 0.0, 0.0])
    r = transport_species(X, ptr, idx, cptr, cidx, [tube(X, 0.0, p, 1.0)], [[0]], u0, 0.0, p)
    fl = r["flux"][0]
    w = p["cal_win"]
    m = (fl > 0) & (X[:, 0] >= w[0]) & (X[:, 0] <= w[1])
    fld = RetardedField(X, float(p["tau_c"]))
    order = np.nonzero(fl > 0)[0]
    phis = []
    for v in order:                                   # own field at the carriers in the window (no back-reaction)
        if m[v]:
            phis.append(fld.value(int(v)))
        fld.add(int(v), fl[v])
    mu = float(fl[m].sum() / (w[1] - w[0]))
    return np.array(phis), mu


def Xm_of(X):
    """Reflection x -> -x: same order, same links; the source arrays (indexed by event) are reflected with it."""
    Xm = np.array(X, copy=True)
    Xm[:, 1] = -Xm[:, 1]
    return Xm


TIE = 1e-9          # a decision margin below this is a numerical tie: a different choice there is not a frame effect


def compare(ref, oth, pairs, map_dir, exact):
    """
    Compare two runs of the transport, species pairs [(i_ref, i_oth), ...]. Returns 'identical', 'tie' or 'differs'.
    identical: the same child at every event, the same flux bit for bit, and the same directions (bit for bit if
    exact, else within 1e-9 after mapping the other run's directions back with map_dir). tie: the first event (in
    the processing order) where a choice differs has a margin < TIE in either run, so the difference starts at a
    numerical tie and says nothing about the rule. differs: anything else.
    """
    first = None
    for ir, io in pairs:
        dif = np.nonzero(ref["dec"][ir] != oth["dec"][io])[0]
        if dif.size:
            v = int(dif[0])
            if first is None or v < first[0]:
                first = (v, ir, io)
    if first is not None:
        v, ir, io = first
        return "tie" if min(ref["gap"][ir][v], oth["gap"][io][v]) < TIE else "differs"
    for ir, io in pairs:
        if not np.array_equal(ref["flux"][ir], oth["flux"][io]):
            return "differs"
        a, b = ref["dir"][ir], map_dir(oth["dir"][io])
        if exact:
            if not np.array_equal(a, b):
                return "differs"
        elif float(np.max(np.abs(a - b))) > 1e-9:
            return "differs"
    return "identical"


def window(d, p):
    """Fit window: 'full' [d + 2, T - 4] (attraction and the superseded law part, as run on 23 Oct); 'lr' the linear-
    response window [1.25 d, 2.5 d] (capped at T - 4), declared 23 Oct after the QA review, before momentum and
    universality were run: the displacement stays a small fraction of d and the window is the same in units of d."""
    if str(p.get("window", "full")) == "lr":
        return 1.25 * d, min(2.5 * d, float(p["T"]) - 4.0)
    return d + 2.0, float(p["T"]) - 4.0


def pair_task(p, seed, d, fL, fR, lam, checks=True):
    """All branches for one sprinkling and one configuration (module notes)."""
    X, ptr, idx, cptr, cidx = make(p, seed, d)
    u0 = np.array([1.0, 0.0, 0.0])
    out = {"N": len(X)}
    resp_B = [[0, 1], [1]] if p.get("external_R") else [[0, 1], [0, 1]]   # universality: R is an external source
    sL, sR = tube(X, -0.5 * d, p, fL), tube(X, 0.5 * d, p, fR)
    aL = transport_species(X, ptr, idx, cptr, cidx, [sL], [[0]], u0, lam, p)
    aR = transport_species(X, ptr, idx, cptr, cidx, [sR], [[0]], u0, lam, p)
    bb = transport_species(X, ptr, idx, cptr, cidx, [sL, sR], resp_B, u0, lam, p)
    out["main"] = {"A_L": species_series(X, aL["flux"][0], aL["dir"][0], p),
                   "A_R": species_series(X, aR["flux"][0], aR["dir"][0], p),
                   "B_L": species_series(X, bb["flux"][0], bb["dir"][0], p),
                   "B_R": species_series(X, bb["flux"][1], bb["dir"][1], p)}
    out["R_hash"] = hashlib.sha256(bb["dec"][1].tobytes()).hexdigest()   # the external field is the same field
    fl = bb["flux"][0] + bb["flux"][1]
    out["merge_excess"] = float((bb["normP"][0] + bb["normP"][1] - fl).sum() / fl.sum())
    if checks:
        both = [(0, 0), (1, 1)]
        cc = transport_species(X, ptr, idx, cptr, cidx, [sL, sR], [[0], [1]], u0, lam, p)
        # C against A: the species of C against the separate runs; exact (the same arithmetic)
        cA = {k: [aL[k][0], aR[k][0]] for k in ("dec", "gap", "flux", "dir")}
        out["C_equals_A"] = compare(cA, cc, both, lambda D_: D_, True)
        bm = transport_species(Xm_of(X), ptr, idx, cptr, cidx, [sL, sR], resp_B, u0, lam, p)
        out["mirror"] = compare(bb, bm, both, lambda D_: D_ * np.array([1.0, -1.0, 1.0]), False)
        bc = float(p["beta_check"])
        ub = boost_points(u0[None, :], bc)[0]
        bf = transport_species(boost_points(X, bc), ptr, idx, cptr, cidx, [sL, sR], resp_B, ub, lam, p)
        out["frame"] = compare(bb, bf, both, lambda D_: boost_points(D_, -bc), False)
    return out


def measures(res, d, p):
    """Relative acceleration, left acceleration towards R, momentum rates (each relative to the alone reference)."""
    t = np.asarray(res["B_L"]["t"])
    t_on, t_end = window(d, p)
    rB = np.asarray(res["B_R"]["x"]) - np.asarray(res["B_L"]["x"])
    rA = np.asarray(res["A_R"]["x"]) - np.asarray(res["A_L"]["x"])
    dxL = np.asarray(res["B_L"]["x"]) - np.asarray(res["A_L"]["x"])
    dPL = np.asarray(res["B_L"]["P"]) - np.asarray(res["A_L"]["P"])
    dPR = np.asarray(res["B_R"]["P"]) - np.asarray(res["A_R"]["P"])
    win = (t >= t_on) & (t <= t_end)
    ML = float(np.mean(np.asarray(res["A_L"]["F"])[win]))
    MR = float(np.mean(np.asarray(res["A_R"]["F"])[win]))
    sL, sR = lin_slope(t, dPL, t_on, t_end), lin_slope(t, dPR, t_on, t_end)
    dr = rB - rA
    return {"a_rel": quad_acc(t, dr, t_on, t_end), "a_L": quad_acc(t, dxL, t_on, t_end),
            "dPL": sL, "dPR": sR, "aP_L": sL / ML if ML > 0 else float("nan"),
            "aP_R": sR / MR if MR > 0 else float("nan"),
            "dr_end": float(dr[win][-1] - dr[win][0]) if win.any() else float("nan")}


def unit_task(p, seed, lam):
    """Small deterministic checks of the implementation on one sprinkling (module notes, part unit)."""
    d = float(p["d_unit"])
    X, ptr, idx, cptr, cidx = make(p, seed, d)
    u0 = np.array([1.0, 0.0, 0.0])
    sL, sR = tube(X, -0.5 * d, p, 1.0), tube(X, 0.5 * d, p, 1.0)
    out = {}
    # (1) no coupling strength: the pair is the two separate bundles, bit for bit
    z0 = transport_species(X, ptr, idx, cptr, cidx, [sL, sR], [[0, 1], [0, 1]], u0, 0.0, p)
    zL = transport_species(X, ptr, idx, cptr, cidx, [sL], [[0]], u0, 0.0, p)
    zR = transport_species(X, ptr, idx, cptr, cidx, [sR], [[0]], u0, 0.0, p)
    zA = {k: [zL[k][0], zR[k][0]] for k in ("dec", "gap", "flux", "dir")}
    out["lambda0_equals_alone"] = compare(zA, z0, [(0, 0), (1, 1)], lambda D_: D_, True) == "identical"
    # (2) scalar flux conservation: what the sources emit ends at the events without children
    bb = transport_species(X, ptr, idx, cptr, cidx, [sL, sR], [[0, 1], [0, 1]], u0, lam, p)
    sink = np.diff(cptr) == 0
    cons = []
    for i, s_ in enumerate((sL, sR)):
        cons.append(abs(float(bb["flux"][i][sink].sum()) - float(s_.sum())) <= 1e-9 * float(s_.sum()))
    out["flux_conserved"] = all(cons)
    # (3) sign: each bundle responding only to the partner's field is kicked towards the partner
    pp_ = transport_species(X, ptr, idx, cptr, cidx, [sL, sR], [[1], [0]], u0, lam, p)
    late = X[:, 0] > d + 1.0
    kL, kR = float(pp_["kick"][0][late].sum()), float(pp_["kick"][1][late].sum())
    nzL = late & (pp_["kick"][0] != 0.0)
    out["kick_L"], out["kick_R"] = kL, kR
    out["kick_L_frac_towards"] = float((pp_["kick"][0][nzL] > 0).mean()) if nzL.any() else float("nan")
    out["kick_sign"] = kL > 0.0 and kR < 0.0
    # (4) the 4-momentum deficit of the model (reported, no verdict): |P| - tot at the mergers
    fl = bb["flux"][0] + bb["flux"][1]
    out["merge_excess"] = float((bb["normP"][0] + bb["normP"][1] - fl).sum() / fl.sum())
    out["n_merge"] = int(((bb["normP"][0] - bb["flux"][0]) > 1e-12).sum() + ((bb["normP"][1] - bb["flux"][1]) > 1e-12).sum())
    return out


def predicted(d, fL, fR, mu, lam, p):
    """Quasi-static point-line prediction (module notes): accelerations through the same fits, and the smallest
    predicted separation in the window (the approximation needs it to stay close to d)."""
    t, xL, xR = predict(d, fL * mu, fR * mu, lam, p)
    t_on, t_end = window(d, p)
    m = (t >= t_on) & (t <= t_end)
    vL = np.gradient(xL, t)
    return {"a_rel": quad_acc(t, (xR - xL) - d, t_on, t_end), "a_L": quad_acc(t, xL + 0.5 * d, t_on, t_end),
            "aP_L": lin_slope(t, vL, t_on, t_end), "r_min": float((xR - xL)[m].min()) if m.any() else float("nan")}


# ============================================================================================== run

def run(p, ctx):
    mode = str(p["mode"])
    t0 = time.time()
    n, B, nse = int(p["n_seeds"]), int(p["boot_B"]), float(p["n_se"])
    cal = parallel_map(calibration_task, [(p, derive_seed(ctx.seed, "F3.10-cal", k)) for k in range(int(p["n_cal"]))],
                       ctx.workers, ctx.progress, lambda i, r_, dd, nn: f"[F3.10] calibration {dd}/{nn} ({time.time()-t0:.0f}s)")
    phi_ref = float(np.median(np.concatenate([c[0] for c in cal])))
    mu = float(np.mean([c[1] for c in cal]))
    lam = float(p["epsilon"]) / phi_ref
    lines = [f"2+1, rho = 1; bundles radius {p['r_src']:g}, source time [0, {p['t_src']:g}), box margin {p['margin']:g}, "
             f"T = {p['T']:g}; field phi_i = sum F / (2 pi tau), tau >= {p['tau_c']:g}; epsilon = {p['epsilon']} (from F3.8/F3.9), "
             f"phi_ref = {phi_ref:.3f}, lambda = {lam:.5f}; flux per unit time of a single bundle mu = {mu:.2f}; window "
             + ("[1.25 d, 2.5 d] (linear response)" if str(p.get("window", "full")) == "lr" else "[d + 2, T - 4]")]

    if mode == "unit":
        seeds = [derive_seed(ctx.seed, "F3.10-unit", k) for k in range(n)]
        U_ = parallel_map(unit_task, [(p, sd, lam) for sd in seeds], ctx.workers, ctx.progress,
                          lambda i, r_, dd, nn: f"[F3.10] unit {dd}/{nn}")
        P_ = parallel_map(pair_task, [(p, sd, float(p["d_unit"]), 1.0, 1.0, lam, True) for sd in seeds], ctx.workers,
                          ctx.progress, lambda i, r_, dd, nn: f"[F3.10] unit checks {dd}/{nn}")
        checks = {"lambda = 0: the pair equals the two separate bundles bit for bit": all(u["lambda0_equals_alone"] for u in U_),
                  "scalar flux conserved (sources = flux at the events without children)": all(u["flux_conserved"] for u in U_),
                  "sign: a bundle responding only to the partner's field is kicked towards it": all(u["kick_sign"] for u in U_),
                  "C equals A bit for bit (choices, flux, directions)": all(q["C_equals_A"] == "identical" for q in P_),
                  "mirror x -> -x (choices, flux, directions)": all(q["mirror"] in ("identical", "tie") for q in P_),
                  "frame beta' = 0.6 (choices, flux, directions boosted back)": all(q["frame"] in ("identical", "tie") for q in P_)}
        for k_, v_ in checks.items():
            lines.append(f"  {k_}: {'OK' if v_ else 'FAILS'}")
        lines.append("  kick towards the partner: fraction of the kicked events of L "
                     + ", ".join(f"{u['kick_L_frac_towards']:.3f}" for u in U_)
                     + "; total kick L " + ", ".join(f"{u['kick_L']:+.3f}" for u in U_)
                     + ", R " + ", ".join(f"{u['kick_R']:+.3f}" for u in U_))
        lines.append("  numerical ties (first differing choice at a margin < 1e-9): mirror "
                     f"{sum(q['mirror'] == 'tie' for q in P_)}, frame {sum(q['frame'] == 'tie' for q in P_)} of {len(P_)}")
        lines.append("  theory (no verdict): the model conserves the scalar flux, not the 4-momentum: at a merger of different "
                     "directions |P| > F and the excess is dropped; (|P| - F) summed / flux summed = "
                     + ", ".join(f"{u['merge_excess']:.2e}" for u in U_) + "; merger events "
                     + ", ".join(str(u["n_merge"]) for u in U_))
        status = "PASSED" if all(checks.values()) else "INSTRUMENT_FAILED"
        return Result(status, lines, {"checks": checks, "unit": U_, "lambda": lam})

    if mode in ("attraction", "law"):
        configs = [(float(d), 1.0, 1.0) for d in p["seps"]]
    elif mode == "momentum":
        configs = [(float(p["d_mass"]), 1.0, 2.0)]
    else:
        configs = [(float(p["d_mass"]), 1.0, 2.0), (float(p["d_mass"]), 2.0, 2.0)]
    jobs, keys = [], []
    for ci, (d, fL, fR) in enumerate(configs):
        for k in range(n):
            # the universality configurations share the sprinklings (paired); the others have their own
            sd = derive_seed(ctx.seed, "F3.10", mode, d if mode != "universality" else "u", k)
            jobs.append((p, sd, d, fL, fR, lam, True))
            keys.append((ci, k))
    res = parallel_map(pair_task, jobs, ctx.workers, ctx.progress,
                       lambda i, r_, dd, nn: f"[F3.10] {mode} config {keys[i][0]} sprinkling {keys[i][1]+1} ({dd}/{nn}, "
                                             f"{time.time()-t0:.0f}s)")
    inst = {"C equals A bit for bit (coupling only through the field)": all(r_["C_equals_A"] == "identical" for r_ in res),
            "frame check beta' = 0.6 (choices, flux, directions)": all(r_["frame"] in ("identical", "tie") for r_ in res),
            "mirror x -> -x (choices, flux, directions)": all(r_["mirror"] in ("identical", "tie") for r_ in res)}
    if mode == "universality":
        by = {}
        for r_, kk in zip(res, keys):
            by.setdefault(kk[1], []).append(r_["R_hash"])
        inst["the external source R is the same in both configurations"] = all(len(set(h)) == 1 for h in by.values())
    n_tie = sum(r_["frame"] == "tie" for r_ in res) + sum(r_["mirror"] == "tie" for r_ in res)
    per = {}
    for ci, (d, fL, fR) in enumerate(configs):
        rows = [r_ for r_, kk in zip(res, keys) if kk[0] == ci]
        M = [measures(r_["main"], d, p) for r_ in rows]
        bs = lambda *tag, ci=ci: derive_seed(ctx.seed, "F3.10-boot", mode, ci, *tag)
        st = {q: boot_se([m_[q] for m_ in M], B, bs(q))[:2] for q in ("a_rel", "a_L", "dPL", "dPR", "aP_L", "aP_R", "dr_end")}
        st["sum_dP"] = boot_se([m_["dPL"] + m_["dPR"] for m_ in M], B, bs("sum"))[:2]
        st["pred"] = predicted(d, fL, fR, mu, lam, p)
        st["raw"] = M
        per[ci] = st
        pr = st["pred"]
        lines.append(f"  d = {d:g}, f_L = {fL:g}, f_R = {fR:g}: relative acceleration a = {st['a_rel'][0]:+.5f} +- {st['a_rel'][1]:.5f} "
                     f"(quasi-static point-line approximation {pr['a_rel']:+.5f}, ratio "
                     f"{st['a_rel'][0] / pr['a_rel'] if pr['a_rel'] else float('nan'):.3f}; its smallest separation in the window "
                     f"{pr['r_min']:.2f}); change of the separation over the window {st['dr_end'][0]:+.2f} +- {st['dr_end'][1]:.2f}; "
                     f"a_L towards R {st['a_L'][0]:+.5f} +- {st['a_L'][1]:.5f}; dP_L/dt {st['dPL'][0]:+.4f} +- {st['dPL'][1]:.4f}, "
                     f"dP_R/dt {st['dPR'][0]:+.4f} +- {st['dPR'][1]:.4f}, sum {st['sum_dP'][0]:+.4f} +- {st['sum_dP'][1]:.4f}; "
                     f"acceleration from momentum a_P,L {st['aP_L'][0]:+.5f} +- {st['aP_L'][1]:.5f} (approximation "
                     f"{pr['aP_L']:+.5f}), a_P,R {st['aP_R'][0]:+.5f} +- {st['aP_R'][1]:.5f}")
    lines.append(f"  4-momentum dropped at mergers (theory, no verdict): (|P| - F) / F = "
                 f"{np.mean([r_['merge_excess'] for r_ in res]):.2e}; numerical ties in the checks: {n_tie}")
    status = "UNDECIDED"
    if mode == "attraction":
        a, se = per[0]["a_rel"]
        status = "PASSED" if a + nse * se < 0 else ("FAILED" if a - nse * se >= 0 else "UNDECIDED")
        lines.append(f"attraction: a + 2 SE = {a + nse * se:+.5f} (needs < 0); a - 2 SE = {a - nse * se:+.5f}")
    elif mode == "law":
        ds = np.array([c[0] for c in configs])
        a_m = np.array([per[ci]["a_rel"][0] for ci in range(len(configs))])
        a_p = np.array([per[ci]["pred"]["a_rel"] for ci in range(len(configs))])
        attract = all(per[ci]["a_rel"][0] + nse * per[ci]["a_rel"][1] < 0 for ci in per)

        def kfit(av):
            av = np.asarray(av, float)
            if np.any(av >= 0) or not np.all(np.isfinite(av)):
                return float("nan")
            return float(-np.polyfit(np.log(ds), np.log(-av), 1)[0])
        k_meas, k_pred = kfit(a_m), kfit(a_p)
        rng = np.random.default_rng(derive_seed(ctx.seed, "F3.10-kboot"))
        raws = [per[ci]["raw"] for ci in range(len(configs))]
        ks = []
        for _ in range(B // 4):
            means = [np.mean([raws[ci][j]["a_rel"] for j in rng.integers(0, len(raws[ci]), len(raws[ci]))])
                     for ci in range(len(configs))]
            ks.append(kfit(means))
        ks = np.asarray(ks)
        k_se = float(np.nanstd(ks, ddof=1)) if np.isfinite(ks).sum() > 2 else float("nan")
        lines.append(f"law (SUPERSEDED by F3.10b; the quasi-static prediction is not valid where its separation collapses): "
                     f"exponent of |a| against d: measured k = {k_meas:.3f} +- {k_se:.3f}, approximation k_pred = {k_pred:.3f}; "
                     f"every d attracts: {attract}; bootstrap replicates without a defined exponent "
                     f"{int((~np.isfinite(ks)).sum())} of {ks.size}; ratios a / a_pred: "
                     + ", ".join(f"d {d_:g}: {m_ / p_:.3f}" for d_, m_, p_ in zip(ds, a_m, a_p)))
        if np.isfinite(k_meas) and np.isfinite(k_se) and np.isfinite(k_pred):
            if attract and abs(k_meas - k_pred) + nse * k_se <= 0.25:
                status = "PASSED"
            elif abs(k_meas - k_pred) - nse * k_se > 0.25:
                status = "FAILED"
    elif mode == "momentum":
        st = per[0]
        (pl, sl_), (pr_, sr_), (sm, ss) = st["dPL"], st["dPR"], st["sum_dP"]
        tol = (abs(pl) + abs(pr_)) / 8.0
        resolved = abs(pl) - nse * sl_ > 0 and abs(pr_) - nse * sr_ > 0
        lines.append(f"momentum: dP_L + dP_R = {sm:+.4f} +- {ss:.4f} vs tolerance (|dP_L| + |dP_R|)/8 = {tol:.4f}; both "
                     f"changes resolved: {resolved}")
        if resolved and abs(sm) + nse * ss <= tol:
            status = "PASSED"
        elif abs(sm) - nse * ss > tol:
            status = "FAILED"
    else:
        (a1, s1), (a2, s2) = per[0]["aP_L"], per[1]["aP_L"]
        diffs = [m1["aP_L"] - m2["aP_L"] for m1, m2 in zip(per[0]["raw"], per[1]["raw"])]
        dm, dse = boot_se(diffs, B, derive_seed(ctx.seed, "F3.10-uni"))[:2]
        amean = 0.5 * (a1 + a2)
        tol = abs(amean) / 6.0
        resolved = a1 - nse * s1 > 0 and a2 - nse * s2 > 0          # towards R is positive for the left bundle
        lines.append(f"universality (acceleration from the carried momentum, external source R): a_P,L(f = 1) = {a1:+.5f} +- "
                     f"{s1:.5f}, a_P,L(f = 2) = {a2:+.5f} +- {s2:.5f}, paired difference {dm:+.5f} +- {dse:.5f} vs tolerance "
                     f"|a|/6 = {tol:.5f}; both resolved towards the partner: {resolved}; centroid a_L "
                     f"{per[0]['a_L'][0]:+.5f} +- {per[0]['a_L'][1]:.5f} and {per[1]['a_L'][0]:+.5f} +- {per[1]['a_L'][1]:.5f}")
        if resolved and abs(dm) + nse * dse <= tol:
            status = "PASSED"
        elif abs(dm) - nse * dse > tol:
            status = "FAILED"
    if not all(inst.values()):
        status = "INSTRUMENT_FAILED"
    lines.append("instrument: " + "; ".join(f"{k}: {'OK' if v else 'FAILS'}" for k, v in inst.items()))
    data = {"phi_ref": phi_ref, "mu": mu, "lambda": lam, "instrument": inst,
            "per": {str(ci): {k: (list(v) if isinstance(v, tuple) else v) for k, v in st.items()} for ci, st in per.items()}}
    return Result(status, lines, data)
