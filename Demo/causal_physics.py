"""
Causal Physics: an interactive demonstration
=============================================

A guided tour of Causal Physics (P. Sikula, "Physics of Causality", doi:10.5281/zenodo.24427611) in seven chapters.
Every picture in chapters 1-5 is computed live from a causal order of events: a Poisson sprinkling, its causal
relation, and nothing else. Chapter 6 shows the field equation from the laboratory's frozen records, and chapter 7
the status of every claim.

    1  The world grows         events appear one by one; time is the order growing; space and dimension are read
                               from the order; a sprinkling has no rest frame, a lattice has one
    2  The order of events     past and future of an event, the links, the longest chain as proper time
    3  Gravity as a response   clocks, cones and density near matter; the three responses and gamma; a thrown
                               ball whose path is the longest chain of the order
    4  Matter as a defect      a kink of a field carried by the retarded response of the order: identity, motion,
                               Lorentz contraction, inertia equal to rest energy, swerves
    5  The field equation      the action of the order plus the action of matter: Einstein's equation at first
                               order, with a form factor counted from intervals
    6  What is known           tested, derived, adopted and open, item by item

Each chapter has two layers of text: "In plain words" and "For physicists" (numbers, test IDs, formulas).
Demonstrations are not laboratory tests: where a chapter shows the same physics as a test, the test's ID and its
frozen result are given next to it.

Requirements: Python 3.9+, numpy, tkinter.
    python Demo/causal_physics.py                 the tour
    python Demo/causal_physics.py --chapter 4     start at a chapter
    python Demo/causal_physics.py --selftest      numerical checks without a window
"""
from __future__ import annotations

import argparse
import math
import sys
import time

import numpy as np

# ===================================================================================================== numerics


def sprinkle_diamond(rho, T, rng, dim=2):
    """Poisson sprinkling of the causal diamond between (t, x) = (0, 0) and (T, 0) in 1+1 (dim 2) or 2+1 (dim 3).
    Returns t, x (and y for dim 3), sorted by t."""
    if dim == 2:
        n = rng.poisson(rho * T * T / 2.0)
        u, w = rng.uniform(0, T, n), rng.uniform(0, T, n)
        t, x = 0.5 * (u + w), 0.5 * (w - u)
        o = np.argsort(t)
        return t[o], x[o], None
    vol = 2.0 * math.pi * (T / 2.0) ** 3 / 3.0
    n = rng.poisson(rho * vol)
    t = np.empty(0)
    x = np.empty(0)
    y = np.empty(0)
    while t.size < n:
        m = 4 * n + 100
        tt = rng.uniform(0, T, m)
        xx = rng.uniform(-T / 2, T / 2, m)
        yy = rng.uniform(-T / 2, T / 2, m)
        ok = np.hypot(xx, yy) <= np.minimum(tt, T - tt)
        t, x, y = np.r_[t, tt[ok]], np.r_[x, xx[ok]], np.r_[y, yy[ok]]
    t, x, y = t[:n], x[:n], y[:n]
    o = np.argsort(t)
    return t[o], x[o], y[o]


def relation(t, x, y=None):
    """R[i, j] = event i precedes event j (strictly inside the future light cone)."""
    dt = t[None, :] - t[:, None]
    if y is None:
        d = np.abs(x[None, :] - x[:, None])
    else:
        d = np.hypot(x[None, :] - x[:, None], y[None, :] - y[:, None])
    return dt > d


def links(R):
    """Links: i precedes j with nothing in between."""
    Rf = R.astype(np.float32)
    return R & ~((Rf @ Rf) > 0.5)


def mm_fraction(d):
    """Myrheim-Meyer: the expected ordering fraction of an interval in d spacetime dimensions."""
    return math.gamma(d + 1.0) * math.gamma(d / 2.0) / (2.0 * math.gamma(1.5 * d))


def mm_dimension(R):
    """The dimension read from the order: solve mm_fraction(d) = 2 R / (N (N - 1))."""
    n = R.shape[0]
    if n < 3:
        return float("nan"), float("nan")
    f = 2.0 * float(R.sum()) / (n * (n - 1.0))
    lo, hi = 1.0, 10.0
    if not (mm_fraction(hi) < f < mm_fraction(lo)):
        return float("nan"), f
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if mm_fraction(mid) > f:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi), f


def maximal(idx, u, w):
    """The maximal elements (nothing above them) of a 2D dominance set: u, w are the light-cone coordinates."""
    if len(idx) == 0:
        return np.array([], int)
    idx = np.asarray(idx)
    o = idx[np.argsort(-u[idx], kind="stable")]
    keep, best = [], -np.inf
    for i in o:
        if w[i] > best:
            keep.append(i)
            best = w[i]
    return np.array(keep, int)


def longest_chain(R, a, b, t):
    """The longest chain from event a to event b in the order R (a precedes b). Returns the list of events."""
    inside = np.flatnonzero(R[a] & R[:, b])
    nodes = np.r_[a, inside[np.argsort(t[inside])], b]
    sub = R[np.ix_(nodes, nodes)]
    L = np.zeros(nodes.size, int)
    prev = -np.ones(nodes.size, int)
    for j in range(1, nodes.size):
        cand = np.flatnonzero(sub[:j, j])
        if cand.size:
            k = cand[np.argmax(L[cand])]
            L[j], prev[j] = L[k] + 1, k
    path, j = [], nodes.size - 1
    while j >= 0:
        path.append(int(nodes[j]))
        j = prev[j]
    return path[::-1]


# ---------------------------------------------------------------------------------------- gravity (chapter 3)

def responses(phi, mode, ds=3):
    """Spatial potential psi, coordinate cone speed c and coordinate density rho_c for the three responses of the
    order to a clock potential phi (metric -e^{2 phi} dt^2 + e^{-2 psi} dx^2), with gamma = psi / phi."""
    gam = {"fixed": -1.0, "cones": 1.0 / ds, "throughput": 2.0 / (ds - 1.0)}[mode]
    psi = gam * phi
    return psi, np.exp(phi + psi), np.exp(phi - ds * psi), gam


def ball_order(g, T, rho, rng, xmin=-2.0, xmax=None):
    """A 1+1 order with the clock field of uniform gravity, ds^2 = -e^{2 g x} dt^2 + dx^2, sprinkled with proper
    density rho over the causal diamond from (0, 0) to (T, 0). Causal relation: |X_j - X_i| < t_j - t_i with the
    optical coordinate X = (1 - e^{-g x}) / g. Returns t, x, R and the indices of the two ends."""
    if xmax is None:
        xmax = 0.75 * T
    X = lambda xx: (1.0 - np.exp(-g * xx)) / g
    area_box = T * (xmax - xmin)
    dmax = math.exp(g * xmax)
    n = rng.poisson(rho * area_box * dmax)
    t = rng.uniform(0, T, n)
    x = rng.uniform(xmin, xmax, n)
    keep = rng.uniform(0, dmax, n) < np.exp(g * x)                    # proper density: sqrt(-g) = e^{g x}
    t, x = t[keep], x[keep]
    Xo = X(x)
    inside = np.abs(Xo) < np.minimum(t, T - t)                         # the interval between the two ends
    t, x, Xo = t[inside], x[inside], Xo[inside]
    t = np.r_[0.0, t, T]
    x = np.r_[0.0, x, 0.0]
    Xo = np.r_[0.0, Xo, 0.0]
    o = np.argsort(t, kind="stable")
    t, x, Xo = t[o], x[o], Xo[o]
    R = (t[None, :] - t[:, None]) > np.abs(Xo[None, :] - Xo[:, None])
    a, b = int(np.flatnonzero(t == 0.0)[0]), int(np.flatnonzero(t == T)[0])
    return t, x, R, a, b


def ball_geodesic(g, T, n=2000):
    """The timelike geodesic of -e^{2 g x} dt^2 + dx^2 from (0, 0) to (T, 0) (it maximises the proper time):
    shooting on the initial slope; returns t, x."""
    def shoot(v0):
        dt = T / n
        x, v = 0.0, v0
        xs = [0.0]
        for _ in range(n):
            # x'' = -g (e^{2 g x} - 2 x'^2): Euler-Lagrange of the proper time sqrt(e^{2 g x} - x'^2) dt
            k1x, k1v = v, -g * (math.exp(2 * g * x) - 2 * v * v)
            k2x, k2v = v + 0.5 * dt * k1v, -g * (math.exp(2 * g * (x + 0.5 * dt * k1x)) - 2 * (v + 0.5 * dt * k1v) ** 2)
            k3x, k3v = v + 0.5 * dt * k2v, -g * (math.exp(2 * g * (x + 0.5 * dt * k2x)) - 2 * (v + 0.5 * dt * k2v) ** 2)
            k4x, k4v = v + dt * k3v, -g * (math.exp(2 * g * (x + dt * k3x)) - 2 * (v + dt * k3v) ** 2)
            x += dt * (k1x + 2 * k2x + 2 * k3x + k4x) / 6.0
            v += dt * (k1v + 2 * k2v + 2 * k3v + k4v) / 6.0
            xs.append(x)
        return np.array(xs)
    lo, hi = 0.0, 0.99
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if shoot(mid)[-1] > 0.0:
            hi = mid
        else:
            lo = mid
    return np.linspace(0.0, T, n + 1), shoot(0.5 * (lo + hi))


# ---------------------------------------------------------------------------------------- the defect (chapter 4)

def _gd(y):
    """The Gudermannian, the integral of sech."""
    return 2.0 * np.arctan(np.tanh(0.5 * y))


def kink_region(rho, X, T, rng):
    """Events of the domain of dependence |x| <= X - t, 0 < t < T (1+1), Poisson with density rho."""
    n = rng.poisson(rho * (2.0 * X * T - T * T))
    t = np.empty(0)
    x = np.empty(0)
    while t.size < n:
        m = 2 * n + 100
        tt, xx = rng.uniform(0, T, m), rng.uniform(-X, X, m)
        ok = np.abs(xx) <= X - tt
        t, x = np.r_[t, tt[ok]], np.r_[x, xx[ok]]
    return t[:n], x[:n]


def kink_free(t, x, v, x0=0.0):
    """The free wave (d'Alembert) of the initial data of a sine-Gordon kink (mu = 1) at x0 moving with velocity v."""
    g = 1.0 / math.sqrt(1.0 - v * v)
    f = lambda s: 4.0 * np.arctan(np.exp(g * (s - x0)))
    return 0.5 * (f(x - t) + f(x + t)) - v * (_gd(g * (x + t - x0)) - _gd(g * (x - t - x0)))


class Fenwick:
    def __init__(self, n):
        self.n = n
        self.a = [0.0] * (n + 1)

    def add(self, i, val):
        i += 1
        a, n = self.a, self.n
        while i <= n:
            a[i] += val
            i += i & -i

    def prefix(self, i):
        """Sum of entries 0 .. i - 1."""
        s, a = 0.0, self.a
        while i > 0:
            s += a[i]
            i -= i & -i
        return s


class KinkRun:
    """The field rule of the order in 1+1 (the rule C0 of the laboratory):
        phi(v) = phi_0(v) - (1/rho) sum_{u < v} G(u, v) [sin phi(u) + f chi(t_u)],   G = 1/2 on the causal relation,
    phi_0 the free wave of the initial data on t = 0. The sum over the causal past is a 2D dominance sum in the
    light-cone coordinates; events are processed in the order of w = t + x (every event of the past comes earlier)
    and the sum over u = t - x is a Fenwick tree. Exact, O(N log N)."""

    def __init__(self, rho, X, T, v, rng, f=0.0, t1=1.0, t2=3.0, events=None):
        self.rho, self.X, self.T, self.v, self.f, self.t1, self.t2 = rho, X, T, v, f, t1, t2
        self.t, self.x = events if events is not None else kink_region(rho, X, T, rng)
        self.phi = np.zeros(self.t.size)
        self.phi0 = kink_free(self.t, self.x, v)
        u, w = self.t - self.x, self.t + self.x
        self.order = np.argsort(w, kind="stable")
        self.urank = np.empty(self.t.size, int)
        self.urank[np.argsort(u, kind="stable")] = np.arange(self.t.size)
        self.tree = Fenwick(self.t.size)
        self.done = 0

    def step(self, nmax=4000):
        """Process up to nmax events; returns the fraction done."""
        o, ur, tr = self.order, self.urank, self.tree
        phi, phi0, t = self.phi, self.phi0, self.t
        c = 0.5 / self.rho
        end = min(self.done + nmax, o.size)
        for k in range(self.done, end):
            i = o[k]
            p = phi0[i] - c * tr.prefix(int(ur[i]))
            phi[i] = p
            src = math.sin(p)
            if self.f and self.t1 <= t[i] < self.t2:
                src += self.f
            tr.add(int(ur[i]), src)
        self.done = end
        return end / o.size

    def track(self, dt=0.5):
        """The kink by time slabs: centre and width from a fit of the kink profile (fit_kink_slab)."""
        out = []
        for a in np.arange(0.0, self.T - 1e-9, dt):
            m = (self.t >= a) & (self.t < a + dt)
            xs, ph = self.x[m], self.phi[m]
            if xs.size < 10:
                continue
            core = np.abs(np.mod(ph, 2 * np.pi) - np.pi) < 0.5 * np.pi       # events near the middle of the turn
            guess = float(np.median(xs[core])) if core.any() else 0.0
            fit = fit_kink_slab(xs, ph, guess)
            if fit is not None:
                out.append((a + 0.5 * dt, fit[0], fit[1]))
        return np.array(out)


def fit_kink_slab(x, phi, guess, half=4.0, iters=12):
    """Least-squares fit of a kink profile 4 arctan(exp((x - xc) / w)) + 2 pi n0 to the events of a slab within
    |x - guess| < half (Gauss-Newton on xc and log w; the laboratory's readout fit_kink also fits a profile).
    Returns (xc, w) or None. The width of a continuum kink is w = 1 / gamma (mu = 1)."""
    m = np.abs(x - guess) < half
    if m.sum() < 15:
        return None
    xx, pp = x[m], phi[m]
    n0 = round(float(np.mean(pp[xx < guess - 2.5])) / (2 * math.pi)) if (xx < guess - 2.5).any() else 0
    xc, lw = guess, 0.0
    for _ in range(iters):
        w = math.exp(lw)
        z = np.clip((xx - xc) / w, -40, 40)
        model = 2 * math.pi * n0 + 4.0 * np.arctan(np.exp(z))
        dm = 2.0 / np.cosh(z)                                       # d model / dz
        J = np.stack([-dm / w, -dm * z], axis=1)                    # d/dxc, d/dlw
        r = pp - model
        try:
            step = np.linalg.lstsq(J, r, rcond=None)[0]
        except np.linalg.LinAlgError:
            return None
        xc += float(np.clip(step[0], -1, 1))
        lw += float(np.clip(step[1], -0.5, 0.5))
        if abs(step[0]) < 1e-7 and abs(step[1]) < 1e-7:
            break
    return xc, math.exp(lw)


def fit_velocity(tr, ta, tb):
    m = (tr[:, 0] >= ta) & (tr[:, 0] < tb)
    if m.sum() < 3:
        return float("nan")
    return float(np.polyfit(tr[m, 0], tr[m, 1], 1)[0])


def width_pred(v):
    """The width w of the profile 4 arctan(exp(x / w)) of a continuum kink moving with velocity v: w = 1/gamma."""
    return math.sqrt(1.0 - v * v)


# ------------------------------------------------------------------------------- the field equation (chapter 5)

# From the frozen records of the laboratory (numbers as recorded; see docs/findings.md):
F68G = {   # F6.8g, PASSED: K_hat at n = kT = 4; results/F6_order_dynamics/F6.8g/20261121T183904Z_main_seed0_86ec1e05.json
    "k": [0.045, 0.06, 0.08, 0.10, 0.12],
    "Kc": [0.97095, 0.95245, 0.92275, 0.88810, 0.84918],
    "Kab": [0.97093, 0.95241, 0.92268, 0.88799, 0.84902],
    "Kbb": [0.96127, 0.93697, 0.89766, 0.85175, 0.80017],
    "fit_c": (0.99993, 3.92202, 2.15337), "fit_ab": (0.99995, 3.92630, 2.15610), "fit_bb": (0.99943, 5.08718, 3.05198),
}
F610B = {  # F6.10b, PASSED: source form factors; results/F6_order_dynamics/F6.10b/20261121T091000Z_main_seed0_14801b85.json
    "k": [0.0, 0.03, 0.06, 0.12],
    "Fw": [0.99654, 0.98465, 0.95447, 0.85466],
    "Fn": [1.00601, 0.98980, 0.94949, 0.81791],
}


def ov3_solve(K, Kab, Kbb, Fw, Fn, tau, sig):
    """Static linear response per k (OV3): alpha k^2 [6K w + 4Kab n] = Fw tau, alpha k^2 [4Kab w + 2Kbb n] = Fn sig,
    with alpha k^2 = 1 (only ratios are used). Phi = w, Psi = -(w + n)."""
    A = np.array([[6.0 * K, 4.0 * Kab], [4.0 * Kab, 2.0 * Kbb]])
    w, n = np.linalg.solve(A, [Fw * tau, Fn * sig])
    return w, -(w + n)


def ov3_response(K, Kab, Kbb, Fw, Fn):
    """Dust (tau = -rho, sigma = 0) and a static BPS wall (tau = -3 rho_E, sigma = -2 rho_E), relative to GR."""
    out = {}
    for name, tau, sig in (("dust", -1.0, 0.0), ("wall", -3.0, -2.0)):
        P, S = ov3_solve(K, Kab, Kbb, Fw, Fn, tau, sig)
        P0, S0 = ov3_solve(1, 1, 1, 1, 1, tau, sig)
        out[name] = {"Phi": P, "Psi": S, "Phi/GR": P / P0, "gamma": S / P}
    out["wall/dust"] = out["wall"]["Phi"] / out["dust"]["Phi"]
    return out


def kfit(k, K0, c1, c0):
    k = np.asarray(k, float)
    return K0 - k * k * (c1 * np.log(1.0 / k) + c0)


# ======================================================================================================= texts

TEXTS = {
    "intro": (
        "Causal Physics in one minute",
        """The world is a growing order of events. An event is something that happens; nothing more is assumed about it. Between two events there is either influence (one can affect the other) or none. That relation, "can influence", is the only primitive.

Everything else is read from it:
- time is the order itself, growing;
- space is the set of events that cannot influence one another;
- distances and durations are counted: how many events lie between two events;
- matter is a pattern that persists in the order;
- gravity is the way the order responds to matter: clocks slow down, cones of influence narrow, and events become denser.

The chapters on the left show this step by step. In chapters 1-5 nothing is drawn from a formula of space or time: every point is an event of a random order, and every line is read from the relation "can influence".

Each chapter has two layers. "In plain words" (this tab) tells the story. "For physicists" gives the numbers, the laboratory tests that check each claim, and what is still open.""",
        """Causal Physics (P. Šikula, "Physics of Causality", v2, doi:10.5281/zenodo.24427611) takes the causal order of events as the primitive: order + number = geometry (Sorkin). Every quantitative step is checked by a pre-registered test of the Causality Lab (verdicts PASSED / FAILED / UNDECIDED / INSTRUMENT_FAILED; criteria never changed after a run).

The chain shown in this demonstration:
1. Substrate: a Lorentz-invariant causal set (Poisson sprinkling). Finite-degree substrates with a rest frame are not supported (F0.3, F1.8, F6.1, F6.2; Bombelli-Henson-Sorkin theorem).
2. Geometry is read blindly from the order (K.1, K.16b, K.16c).
3. Gravity is the response of the order: γ = −1 (fixed order), 1/d_s (cones only), 2/(d_s − 1) (throughput conserved), which is 1 in 3+1 (F7.1a/b/c).
4. Matter: flows that carry their four-velocity (F3.6, K.15), and defects of a field of the order with inertial mass = rest energy (F3.14b, F3.15), swerves (F3.14c) and fall (F5.2a3, F5.2b, F5.2c2).
5. The expected Benincasa–Dowker action has the Einstein–Hilbert stiffness in all three static sectors (F6.8a, F6.8b2, F6.8c2, F6.8g); a matter action of the order sources it through T_ab (F6.10a/b); together: the linearised Einstein equation with a non-local spin-0 form factor.

What a demonstration is: chapters 1–5 compute small orders live (thousands of events, seconds). They show the mechanism. The laboratory tests use larger orders, many sprinklings, and criteria fixed in advance; their numbers are quoted next to each demonstration."""),
    "grow": (
        "1  The world grows",
        """Press "Grow". Events appear one by one, each in the future of some earlier ones. Nothing else exists: no space in which they sit, no clock that ticks. The picture places them only so that you can see them.

The thin lines are links: the direct influences, with nothing in between. The bright events at the top are "now": the events that nothing has yet followed. The present is the growing edge of the order.

Time is this growth. Space is what lies side by side: events that cannot influence each other.

Dimension is not put in. It is counted: in a world with more directions of space, fewer pairs of events can influence each other. The readout "dimension from the order" counts the related pairs and returns 2 for this picture (one time and one space direction), and 3 if you switch to 2+1.

The last switch shows why the events must be random. A regular lattice looks different when you move past it: it has a rest frame. A random sprinkling looks the same to every observer. That is how the order can carry special relativity without a preferred frame.""",
        """Sprinkling: Poisson process of density ρ in the causal diamond between (0,0) and (T,0) (1+1) or the double cone (2+1). Relation: x ≺ y iff y is in the future light cone of x. Links: x ≺ y with no z, x ≺ z ≺ y (Hasse diagram).

"Now" = the maximal elements of the grown part (nothing above them yet). The growth shown is in time order of the sprinkling, one realisation of a linear extension; the laboratory requires results independent of the order of realisation (principle M4, the counterpart of causal invariance).

Dimension (Myrheim–Meyer): for an interval in d dimensions the expected ordering fraction is r(d) = Γ(d+1)Γ(d/2) / (2Γ(3d/2)): r(2) = 1/2, r(3) = 0.229, r(4) = 1/10. The readout inverts r = 2R/(N(N−1)).

Lorentz invariance: a Poisson sprinkling is statistically invariant under boosts (the boosted set is again a Poisson process of the same density); a lattice is not, and a substrate with finite degree and a rest frame is excluded by the tests together with the theorem of Bombelli, Henson and Sorkin (Paper Sec. 2; F0.3 degree grows with N, F6.1: 0/7 local growth rules manifold-like, F6.2 no 3D).

Laboratory: K.1 (geometry readable from the order in 3+1 from N ≈ 3·10⁴), K.16b / K.16c (blind instrument reads γ from events and relation only)."""),
    "order": (
        "2  The order of events",
        """Click an event. Its past turns blue and its future orange. These are not light cones drawn by hand: they are just the events this one can be influenced by, and the events it can influence. The cone shape appears by itself.

Click a second event in the orange future. The program finds the longest chain of events between the two: the longest sequence in which each event influences the next. Its length is the proper time between them: the time a clock would show if it travelled from the first to the second.

The readout compares the count with the true proper time of the picture. They agree better the denser the events: geometry is order plus number.""",
        """Interval I(a, b) = {z : a ≺ z ≺ b}. In 1+1 its volume is V = τ²/2 and the expected number of events is ρV.

Longest chain: for a Poisson sprinkling of an interval with n = ρV events, the longest chain is the longest increasing subsequence of n random points (Ulam's problem): L ≈ 2√n − 1.77 n^{1/6} (Baik–Deift–Johansson). Hence τ_chain = L/√(2ρ) → τ as ρ → ∞; the relative fluctuations fall as n^{−1/3}.

In the laboratory, the length of the longest chain is the proper time of a flow (Axiom "Action of a flow", S_m = −Σ m Δτ); K.3 calibrates the clock rate read from chains."""),
    "gravity": (
        "3  Gravity as the response of the order",
        """Near matter the order changes in three ways at once: clocks run slower, the cones of influence narrow (influence travels more slowly), and events become denser. The upper panel shows the three, for three possible ways the order could respond.

- A fixed order: matter changes only the clocks. Light would not bend at all.
- Cones that only narrow: light bends, but too little.
- An order that keeps its throughput (the number of events a light front crosses per unit time stays the same, like cars that bunch up where they slow down): light bends exactly as Einstein predicted, twice Newton's value.

The lower panel is an experiment. A ball is thrown up and caught. In Causal Physics nothing pushes it: its path is the longest chain of events between the throw and the catch. Press "Throw": the program sprinkles an order in which clocks run faster higher up, finds the longest chain, and draws it. It rises and falls like a ball, along the path a physicist calculates from gravity.""",
        """Metric of a static field, ds² = −e^{2φ}dt² + e^{−2ψ}dx², γ = ψ/φ. Dictionary of the order: clocks ↔ φ, cone speed c = e^{φ+ψ}, coordinate density ρ_c = ρ e^{φ−d_sψ}.
- fixed order (c = 1): ψ = −φ, γ = −1 (Nordström class);
- cones only (ρ_c fixed): γ = 1/d_s (1/3 in 3+1);
- throughput ρ_c c = const: γ = 2/(d_s − 1) (= 1 in 3+1, the linearised coupling of GR).
Laboratory: F7.1a (cones only, 2+1: γ ≈ 0.47 against 1/2; 3+1 0.35 against 1/3), F7.1b (throughput 2+1: ≈ 1.85 against 2, UNDECIDED), F7.1c (3+1: γ = 1.000 ± 0.045, UNDECIDED because the control path read low). F7.4a: the tensor rule (each flow narrows the cone in its own rest frame), counter-streaming beams R = 2.401 ± 0.053 against 2.361. The upper panel draws the continuum profiles for d_s = 3 on a sprinkled slab (an illustration).

The ball: 1+1 order with ds² = −e^{2gx}dt² + dx² (uniform field g), proper density ρ, relation from the optical coordinate X = (1 − e^{−gx})/g. The longest chain between (0,0) and (T,0) follows the timelike geodesic, which maximises proper time; g = 0.05, T = 20: the geodesic rises to 2.88 (Newton's slow-motion value gT²/8 = 2.5; the throw is not slow, v₀ ≈ 0.5). A single longest chain fluctuates transversally (KPZ class), so the demonstration averages chains. At finite density the chains lie on average somewhat above the geodesic: in our development runs (16 orders each) the height at the top was 3.52 ± 0.36 at density 30 and 2.95 ± 0.34 at density 60, against 2.78; the bias shrinks as the density grows. This 1+1 picture is a demonstration. The laboratory's test is F5.1d (2+1, order of F7.1b): a history along longest chains from the apex of a thrown ball falls with R = a/a_pred = 1.069 ± 0.051 (replicated 1.101); in 3+1 the same rule does not fall (F5.1e, FAILED), so the motion of matter in 3+1 is open."""),
    "defect": (
        "4  Matter as a defect of a field",
        """Give every event a value, an angle. Each event takes its value from its causal past: the free wave that arrives there, minus the pull of every earlier event it can be influenced by. Nothing else; there is no grid and no time step. The colours show the angle at each event.

Start with a kink: on the left the angle is 0, on the right it has turned once around, 2π. No local change can undo a full turn, so the kink cannot disappear. It is a particle.

Try it:
- At rest it stays (the white track is its centre).
- Launched, it keeps its speed and becomes narrower: the Lorentz contraction of special relativity, without any formula of relativity put in.
- Kicked by a uniform force, it gains exactly as much speed as a body whose mass is its energy: inertia equals rest energy, E = mc².
- Watched long enough at rest, it trembles: the order is random, and its randomness shakes the kink. This is the "swerve" of causal set theory, here derived for a body with structure.""",
        """The rule (C0 in 1+1): φ(v) = φ₀(v) − (1/ρ) Σ_{u≺v} G(u,v) [V′(φ(u)) + f χ(t_u)], G = 1/2 on the causal relation (the retarded Green function of □ in 1+1), V = μ²(1 − cos φ), μ = 1, φ₀ the free wave of the initial data on t = 0 (d'Alembert). Region: the domain of dependence |x| ≤ X − t. Computed exactly with a Fenwick tree in light-cone coordinates.

Readouts: centre x_c and width w from a least-squares fit of the profile 4 arctan(e^{(x−x_c)/w}) per slab of Δt = 0.5 (as the laboratory's fit_kink); velocity from a linear fit of x_c; a continuum kink has w = 1/γ; kick f during 1 ≤ t < 3, read as the difference kicked − at rest on the same events (paired, as F3.14b), gives Δ(γv) = 2πfτ/M in the continuum, M = 8μ, so r_M = M_in/M_rest = (2πfτ/8)/Δ(γv). A single small order is noisy and, at finite density, the order drags the kink slightly (velocity and r_M approach the continuum as the density grows: try 100, 400, 800).

Laboratory (larger orders, frozen criteria):
- F3.14b (1+1, ρ 800, 128 sprinklings): identity, velocity to 0.002, contraction to 0.001, r_M = 0.9988 ± 0.0077 (replication 0.9895 ± 0.0079): PASSED, replicated.
- F3.15 (2+1 wall, Green weight from the longest chain): r_M = 0.978 ± 0.012 and 0.952 ± 0.011: PASSED, replicated.
- F3.14c (swerves): D_p = (1/2ρ)∫V′²φ_x² dx = (32/15)μ^{5/2}/ρ, no free constant; measured ratio 1.016 / 0.991 / 0.882 / 1.008 (± 0.05): PASSED (seed 0).
- F5.2a3, F5.2b, F5.2c2: the defect falls; under a growth rule it falls in proportion to the source energy and independently of its own mass.
Limits: the potential is chosen; 1+1 and 2+1 only; no 3+1 defect dynamics yet."""),
    "fieldeq": (
        "5  The field equation of the order",
        """The order has an action: a count over its intervals (the Benincasa–Dowker action of causal set theory). Matter that is a field of the order has an action too, counted over pairs of related events. The order takes the shape in which the two together are stationary. That condition is the field equation.

The laboratory computed both, with nothing from Einstein's theory put in:
- the stiffness of the order (left plot) comes out as Einstein's at long wavelength in all three directions in which a static field can bend the order;
- the way matter loads the order (middle) comes out as Einstein's source, the energy and the stresses of matter;
- together (right) they give Einstein's equation at first order. At wavelengths of tens of the discreteness scale the order departs from it: the field of dust is weaker, light bends relatively more, and a domain wall pushes instead of pulls.

The first-order equation is the first step. The second order (which decides the perihelion of Mercury), fields that have existed for a long time, and the value of Newton's constant are still open.""",
        """Γ = ⟨S_BD⟩ + ⟨S_m⟩, δΓ/δg_ab = 0 (ensemble mean in the manifold-like phase at linear order; beyond it the balance of action and entropy, Γ_λ = −(1/λ) ln Σ e^{−λS}).

Gravity side, F6.8g (PASSED, n = kT = 4): K̂ = stiffness / EH for the conformal (A/6), mixed (B/8) and optical (C/2) sectors; K₀ = 0.99993, 0.99995, 0.99943. Quadratic action ∫√−g[R + R F₀(□) R + C F₂(□) C], A = 6 + 9φ₀, B = 8 + 12φ₀, C = 2 + 4φ₀ + φ₂, φ₀ = 4k²F₀ = −k²ℓ_k²(2.63 ln(1/kℓ_k) + 1.41).

Matter side, F6.10b (PASSED): S_m = −ε²Σ f(n)(Δχ)² − (1/ρ)ΣV, response (1/2)∫T^ab δg_ab × F_ω (trace) and F_ν (Σ T_ii); F6.10a (2D): T^a_a = −2V, kinetic share −0.014 / −0.005.

Response (OV3-lin-a, derived from the two frozen records, no test of its own): α k² [A ω + B ν/2] = F_ω T^a_a, α k² [B ω/2 + C ν] = F_ν Σ T_ii; Φ = ω, Ψ = −(ω + ν). At kℓ_k = 0.06 / 0.12: dust Φ/Φ_GR = 0.940 / 0.809, γ(k) = 1.033 / 1.122; a static BPS wall (ρ + Σp = −ρ_E) repels, wall/dust lapse −1.045 / −1.062.

The field equation to first order (young static field): G_ab + 2(η_ab□ − ∂_a∂_b)F₀(□)R = 8πG T_ab^eff.
Open: the second order (β; OV3 §10: β − 1 = −κ₃/8, κ₃ the cubic action along the exponential metric; card F6.12 in development, stage β1 conformal cubic → 0.9995); old fields (age growth (kT)⁵, F6.8e3; no saturation at finite amplitude, F6.8f); time-dependent and tensor sectors; the value of G."""),
    "status": (
        "6  What is known, and what is not",
        """Causal Physics is a research programme with a laboratory. Some of its claims are tested, some are derived, some are adopted from known physics, and some are open. The table on the canvas says which is which.

The honest summary: from events and their order, the laboratory reaches Einstein's gravity at first order, matter that has inertia equal to its energy, and a field equation counted from the order itself. It has not yet reached particles in our three dimensions of space, gravity of fields that have existed for a long time, or the quantum.""",
        """Status labels as in the paper: tested (a frozen laboratory test passed), derived (follows from tested or adopted pieces), adopted (taken from known physics, derivation open), open.

The full list of tests with their latest verdicts: docs/test_list.md (python run.py report); the verdicts and their reading: docs/findings.md; the decisions of the author: docs/author_decisions.md. Negative results are kept: e.g. F6.8c (optical sector read linearly in age, FAILED), F6.8f (no saturation at finite amplitude, FAILED), F5.1e (longest chains do not fall in 3+1, FAILED), F3.6c (merging does not bind a bundle, FAILED)."""),
}

STATUS_ROWS = [
    ("Time is the order; no global now", "derived", "Ax. 1"),
    ("Geometry read from order and number", "tested", "K.1, K.16b, K.16c"),
    ("Lorentz-invariant substrate (no rest frame)", "tested", "F0.3, F6.1, F6.2 + theorem"),
    ("Dimension counted from the order", "tested", "K.1"),
    ("Gauss's law for a source", "tested", "F1.1"),
    ("γ = 1 from conserved throughput (3+1)", "derived; raw 1.000 ± 0.045", "F7.1c (undecided)"),
    ("Pressure and motion gravitate; cones tilt", "tested", "F7.4a"),
    ("Einstein–Hilbert stiffness of the order", "tested, 3 sectors", "F6.8a, b2, c2, g"),
    ("Source T_ab from an order-only matter action", "tested", "F6.10a, F6.10b"),
    ("Linear field equation (young static field)", "derived", "OV3-lin-a"),
    ("Defect: inertial mass = rest energy", "tested, 1+1 and 2+1", "F3.14b, F3.15"),
    ("Defect: swerves with no free constant", "tested, 1+1", "F3.14c"),
    ("Defect falls; equivalence", "tested, 1+1 and 2+1", "F5.2a3, F5.2b, F5.2c2"),
    ("Free fall along longest chains", "2+1 passed, 3+1 failed", "F5.1d, F5.1e"),
    ("Non-linear completion = general relativity", "adopted", "Jacobson route open"),
    ("Second order of the field equation (β)", "open (in development)", "F6.12 draft"),
    ("Static response of old fields", "open", "F6.8e3, F6.8f"),
    ("The value of G", "open", ""),
    ("Self-bound matter; particles in 3+1", "open", "F3.8, F3.11, F3.12b"),
    ("Quantum superposition and measurement", "interpretation only", ""),
]


# ========================================================================================================= GUI

BG, PANEL, INK, INK2, GRID = "#fbfaf7", "#f1efe9", "#1c1b19", "#6b6964", "#e3e1db"
C_PAST, C_FUT, C_SEL, C_LINK, C_NOW, C_CHAIN = "#2a78d6", "#eb6834", "#111111", "#c9c6bd", "#d4a017", "#1baf7a"


class View:
    """Maps world coordinates (h horizontal, v vertical, v upwards) onto a canvas box."""

    def __init__(self, x0, y0, w, h, hmin, hmax, vmin, vmax):
        self.x0, self.y0, self.w, self.h = x0, y0, w, h
        self.hmin, self.hmax, self.vmin, self.vmax = hmin, hmax, vmin, vmax

    def p(self, hh, vv):
        return (self.x0 + (hh - self.hmin) / (self.hmax - self.hmin) * self.w,
                self.y0 + self.h - (vv - self.vmin) / (self.vmax - self.vmin) * self.h)

    def inv(self, X, Y):
        return (self.hmin + (X - self.x0) / self.w * (self.hmax - self.hmin),
                self.vmin + (self.y0 + self.h - Y) / self.h * (self.vmax - self.vmin))


def hue_color(a):
    """A cyclic colour for an angle a (radians)."""
    a = (a / (2 * math.pi)) % 1.0
    r = 0.5 + 0.5 * math.cos(2 * math.pi * a)
    g = 0.5 + 0.5 * math.cos(2 * math.pi * (a - 1 / 3))
    b = 0.5 + 0.5 * math.cos(2 * math.pi * (a - 2 / 3))
    return "#%02x%02x%02x" % (int(40 + 200 * r), int(40 + 200 * g), int(40 + 200 * b))


class Chapter:
    key = ""

    def __init__(self, app):
        self.app = app
        self.cv = app.canvas
        self.rng = np.random.default_rng(1)

    def controls(self, frame):
        pass

    def draw(self):
        pass

    def click(self, X, Y):
        pass

    def stop(self):
        pass

    def readout(self, s):
        self.app.readout.set(s)


class Intro(Chapter):
    key = "intro"

    def draw(self):
        cv = self.cv
        cv.delete("all")
        W, H = self.app.CW, self.app.CH
        boxes = [
            ("events", "something happens", 0.5, 0.08),
            ("the order", "\"can influence\"", 0.5, 0.24),
            ("time", "the order growing", 0.17, 0.42),
            ("space, dimension", "events side by side; counted", 0.5, 0.42),
            ("proper time", "the longest chain", 0.83, 0.42),
            ("matter", "flows; defects of a field", 0.25, 0.62),
            ("gravity", "the response of the order", 0.75, 0.62),
            ("the field equation", "action of the order + action of matter", 0.5, 0.82),
        ]
        pos = {}
        for name, sub, fx, fy in boxes:
            x, y = fx * W, fy * H
            cv.create_rectangle(x - 120, y - 30, x + 120, y + 30, fill=PANEL, outline=INK2, width=1.2)
            cv.create_text(x, y - 9, text=name, font=("Helvetica", 13, "bold"), fill=INK)
            cv.create_text(x, y + 12, text=sub, font=("Helvetica", 10), fill=INK2)
            pos[name] = (x, y)
        for a, b in [("events", "the order"), ("the order", "time"), ("the order", "space, dimension"),
                     ("the order", "proper time"), ("space, dimension", "matter"), ("proper time", "gravity"),
                     ("time", "matter"), ("matter", "gravity"), ("matter", "the field equation"),
                     ("gravity", "the field equation")]:
            (x0, y0), (x1, y1) = pos[a], pos[b]
            cv.create_line(x0, y0 + 30, x1, y1 - 30, arrow="last", fill=INK2, width=1.1)
        cv.create_text(W / 2, H - 22, text="Choose a chapter on the left. Chapters 1–4 compute random orders live.",
                       font=("Helvetica", 11, "italic"), fill=INK2)
        self.readout("Causal Physics v2 · doi:10.5281/zenodo.24427611 · Causality Lab")


class Grow(Chapter):
    key = "grow"

    def __init__(self, app):
        super().__init__(app)
        self.dim = 2
        self.boost = 0.0
        self.running = False
        self.new_sprinkling()

    def new_sprinkling(self):
        T = 10.0
        rho = 9.0 if self.dim == 2 else 1.6
        self.t, self.x, self.y = sprinkle_diamond(rho, T, self.rng, self.dim)
        self.T = T
        self.n_shown = 0
        self.R = relation(self.t, self.x, self.y)
        self.L = links(self.R)
        self.u, self.w = self.t - self.x, self.t + self.x

    def controls(self, fr):
        import tkinter as tk
        from tkinter import ttk
        ttk.Button(fr, text="Grow", command=self.start).pack(fill="x", pady=2)
        ttk.Button(fr, text="Pause", command=self.stop).pack(fill="x", pady=2)
        ttk.Button(fr, text="All at once", command=self.all).pack(fill="x", pady=2)
        ttk.Button(fr, text="New order", command=self.reset).pack(fill="x", pady=2)
        self.dimvar = tk.StringVar(value="1+1")
        ttk.Label(fr, text="Dimensions").pack(anchor="w", pady=(10, 0))
        for d in ("1+1", "2+1"):
            ttk.Radiobutton(fr, text=d, value=d, variable=self.dimvar, command=self.set_dim).pack(anchor="w")
        ttk.Label(fr, text="Move past it (boost, rapidity)").pack(anchor="w", pady=(10, 0))
        self.bvar = tk.DoubleVar(value=0.0)
        ttk.Scale(fr, from_=0.0, to=1.2, variable=self.bvar, command=lambda e: self.set_boost()).pack(fill="x")

    def set_dim(self):
        self.stop()
        self.dim = 2 if self.dimvar.get() == "1+1" else 3
        self.new_sprinkling()
        self.draw()

    def set_boost(self):
        self.boost = float(self.bvar.get())
        self.draw()

    def reset(self):
        self.stop()
        self.new_sprinkling()
        self.draw()

    def all(self):
        self.stop()
        self.n_shown = self.t.size
        self.draw()

    def start(self):
        if self.n_shown >= self.t.size:
            self.n_shown = 0
        self.running = True
        self.tick()

    def stop(self):
        self.running = False

    def tick(self):
        if not self.running or self.app.current is not self:
            return
        self.n_shown = min(self.t.size, self.n_shown + max(1, self.t.size // 150))
        self.draw()
        if self.n_shown < self.t.size:
            self.app.root.after(40, self.tick)
        else:
            self.running = False

    def draw(self):
        cv = self.cv
        cv.delete("all")
        W, H = self.app.CW, self.app.CH
        T = self.T
        if self.boost > 0:
            self.draw_boost()
            return
        v = View(40, 20, W - 80, H - 60, -T / 2 - 0.3, T / 2 + 0.3, -0.3, T + 0.3)
        n = self.n_shown
        sh = np.arange(n)
        r = 2.5 if self.dim == 2 else 2.0
        if self.dim == 2:
            Ls = self.L[:n, :n]
            ii, jj = np.nonzero(Ls)
            for i, j in zip(ii, jj):
                cv.create_line(*v.p(self.x[i], self.t[i]), *v.p(self.x[j], self.t[j]), fill=C_LINK)
            now = maximal(sh, self.u, self.w) if n else []
        else:
            now = [i for i in sh if not self.R[i, :n].any()]
        nowset = set(int(i) for i in now)
        for i in sh:
            X, Y = v.p(self.x[i], self.t[i])
            col = C_NOW if i in nowset else INK
            rr = r + 1.5 if i in nowset else r
            cv.create_oval(X - rr, Y - rr, X + rr, Y + rr, fill=col, outline="")
        cv.create_text(48, 30, anchor="nw", text="time ↑ (the order grows)", fill=INK2, font=("Helvetica", 10))
        if n >= 10:
            d, f = mm_dimension(self.R[:n, :n])
            self.readout(f"events {n}   related pairs {int(self.R[:n, :n].sum())}   ordering fraction {f:.3f}   "
                         f"dimension from the order {d:.2f} (true {self.dim})   now: {len(nowset)} events")
        else:
            self.readout(f"events {n}")

    def draw_boost(self):
        cv = self.cv
        W, H = self.app.CW, self.app.CH
        eta = self.boost
        ch, sh = math.cosh(eta), math.sinh(eta)
        half = (W - 60) / 2
        for k, title in enumerate(("random sprinkling", "regular lattice")):
            x0 = 20 + k * (half + 20)
            v = View(x0, 40, half, H - 80, -4, 4, 2, 10)
            cv.create_rectangle(x0, 40, x0 + half, H - 40, outline=GRID)
            cv.create_text(x0 + half / 2, 25, text=f"{title}, seen moving at v = {math.tanh(eta):.2f}",
                           font=("Helvetica", 11, "bold"), fill=INK)
            if k == 0:
                rng = np.random.default_rng(7)
                n = rng.poisson(4.0 * 40 * 40)
                t, x = rng.uniform(-14, 26, n), rng.uniform(-20, 20, n)
            else:
                g = np.arange(-40, 41) * 0.5
                t, x = [a.ravel() for a in np.meshgrid(g + 6, g)]
            tb, xb = ch * t - sh * x, ch * x - sh * t
            m = (xb > -4) & (xb < 4) & (tb > 2) & (tb < 10)
            for a, b in zip(xb[m], tb[m]):
                X, Y = v.p(a, b)
                cv.create_oval(X - 1.6, Y - 1.6, X + 1.6, Y + 1.6, fill=INK, outline="")
        self.readout("A lattice shows its rest frame when you move past it; a Poisson sprinkling looks the same in "
                     "every frame (statistically). Set the slider to 0 to return.")


class Order(Chapter):
    key = "order"

    def __init__(self, app):
        super().__init__(app)
        self.rho = 12.0
        self.make()

    def make(self):
        self.T = 10.0
        self.t, self.x, _ = sprinkle_diamond(self.rho, self.T, self.rng)
        self.R = relation(self.t, self.x)
        self.L = None
        self.a = self.b = None
        self.chain = None

    def controls(self, fr):
        import tkinter as tk
        from tkinter import ttk
        ttk.Label(fr, text="Density of events").pack(anchor="w")
        self.rvar = tk.StringVar(value="12")
        for r in ("4", "12", "30"):
            ttk.Radiobutton(fr, text=r, value=r, variable=self.rvar, command=self.set_rho).pack(anchor="w")
        self.showlinks = tk.BooleanVar(value=False)
        ttk.Checkbutton(fr, text="Show links", variable=self.showlinks, command=self.draw).pack(anchor="w", pady=6)
        ttk.Button(fr, text="New order", command=lambda: (self.make(), self.draw())).pack(fill="x", pady=2)
        ttk.Button(fr, text="Clear selection", command=self.clear).pack(fill="x", pady=2)
        ttk.Label(fr, text="Click an event, then an event\nin its future.", foreground=INK2).pack(anchor="w", pady=8)

    def set_rho(self):
        self.rho = float(self.rvar.get())
        self.make()
        self.draw()

    def clear(self):
        self.a = self.b = None
        self.chain = None
        self.draw()

    def view(self):
        W, H = self.app.CW, self.app.CH
        return View(40, 20, W - 80, H - 60, -self.T / 2 - 0.3, self.T / 2 + 0.3, -0.3, self.T + 0.3)

    def click(self, X, Y):
        v = self.view()
        xs, ts = v.inv(X, Y)
        i = int(np.argmin((self.x - xs) ** 2 + (self.t - ts) ** 2))
        if self.a is None or self.b is not None:
            self.a, self.b, self.chain = i, None, None
        elif self.R[self.a, i]:
            self.b = i
            self.chain = longest_chain(self.R, self.a, self.b, self.t)
        else:
            self.a, self.chain = i, None
        self.draw()

    def draw(self):
        cv = self.cv
        cv.delete("all")
        v = self.view()
        if self.showlinks.get() if hasattr(self, "showlinks") else False:
            if self.L is None:
                self.L = links(self.R)
            ii, jj = np.nonzero(self.L)
            for i, j in zip(ii, jj):
                cv.create_line(*v.p(self.x[i], self.t[i]), *v.p(self.x[j], self.t[j]), fill=C_LINK)
        col = np.array([INK] * self.t.size, dtype=object)
        if self.a is not None:
            col[self.R[:, self.a]] = C_PAST
            col[self.R[self.a]] = C_FUT
        r = 2.6 if self.rho < 20 else 1.8
        for i in range(self.t.size):
            X, Y = v.p(self.x[i], self.t[i])
            cv.create_oval(X - r, Y - r, X + r, Y + r, fill=col[i], outline="")
        if self.chain:
            pts = [c for i in self.chain for c in v.p(self.x[i], self.t[i])]
            cv.create_line(*pts, fill=C_CHAIN, width=2.5)
        for i in (self.a, self.b):
            if i is not None:
                X, Y = v.p(self.x[i], self.t[i])
                cv.create_oval(X - 6, Y - 6, X + 6, Y + 6, outline=C_SEL, width=2)
        if self.a is not None and self.b is not None:
            dt, dx = self.t[self.b] - self.t[self.a], self.x[self.b] - self.x[self.a]
            tau = math.sqrt(dt * dt - dx * dx)
            n_in = int((self.R[self.a] & self.R[:, self.b]).sum())
            Lc = len(self.chain) - 1
            self.readout(f"proper time τ = {tau:.2f}   events in the interval {n_in} (expected ρτ²/2 = "
                         f"{self.rho * tau * tau / 2:.1f})   longest chain {Lc} links → τ ≈ L/√(2ρ) = "
                         f"{Lc / math.sqrt(2 * self.rho):.2f}  (Ulam: {(2 * math.sqrt(max(n_in, 1)) - 1.77 * max(n_in, 1) ** (1 / 6)) / math.sqrt(2 * self.rho):.2f})")
        elif self.a is not None:
            self.readout(f"past: {int(self.R[:, self.a].sum())} events   future: {int(self.R[self.a].sum())} events   "
                         f"(the cones are not drawn: they are the order)")
        else:
            self.readout(f"{self.t.size} events, density ρ = {self.rho:g}. Click an event.")


class Gravity(Chapter):
    key = "gravity"

    def __init__(self, app):
        super().__init__(app)
        self.mode = "throughput"
        self.A = 0.25
        self.ball = None
        self.slab = None

    def controls(self, fr):
        import tkinter as tk
        from tkinter import ttk
        ttk.Label(fr, text="How the order responds").pack(anchor="w")
        self.mvar = tk.StringVar(value=self.mode)
        for m, lab in (("fixed", "fixed order"), ("cones", "cones only"), ("throughput", "throughput kept")):
            ttk.Radiobutton(fr, text=lab, value=m, variable=self.mvar, command=self.set_mode).pack(anchor="w")
        ttk.Label(fr, text="Strength of the field").pack(anchor="w", pady=(8, 0))
        self.avar = tk.DoubleVar(value=self.A)
        ttk.Scale(fr, from_=0.05, to=0.45, variable=self.avar, command=lambda e: self.set_A()).pack(fill="x")
        ttk.Separator(fr).pack(fill="x", pady=10)
        ttk.Label(fr, text="Density of events (ball)").pack(anchor="w")
        self.bvar = tk.StringVar(value="30")
        for r in ("30", "60"):
            ttk.Radiobutton(fr, text=r + (" (slower)" if r == "60" else ""), value=r, variable=self.bvar).pack(anchor="w")
        ttk.Button(fr, text="Throw a ball", command=self.throw).pack(fill="x", pady=2)
        ttk.Button(fr, text="Throw 6 balls (average)", command=lambda: self.throw(6)).pack(fill="x", pady=2)

    def set_mode(self):
        self.mode = self.mvar.get()
        self.slab = None
        self.draw()

    def set_A(self):
        self.A = float(self.avar.get())
        self.slab = None
        self.draw()

    def throw(self, n=1):
        g, T, rho = 0.05, 20.0, float(self.bvar.get())
        paths = []
        self.readout("sprinkling the order and finding the longest chain ...")
        self.app.root.update()
        t0 = time.time()
        for _ in range(n):
            t, x, R, a, b = ball_order(g, T, rho, self.rng)
            ch = longest_chain(R, a, b, t)
            paths.append((t[ch], x[ch], t, x))
        tg, xg = ball_geodesic(g, T)
        self.ball = (paths, tg, xg, g, T, rho, time.time() - t0)
        self.draw()

    def draw(self):
        cv = self.cv
        cv.delete("all")
        W, H = self.app.CW, self.app.CH
        top = View(60, 30, W - 90, H * 0.42, -12, 12, 0, 1)
        # upper panel: profiles across a mass at x = 0
        xs = np.linspace(-12, 12, 241)
        phi = -self.A * np.exp(-xs ** 2 / (2 * 3.0 ** 2))
        psi, c, rc, gam = responses(phi, self.mode)
        clock = np.exp(phi)
        cv.create_text(60, 14, anchor="w", text=f"Near a mass at the centre — response: {self.mode};  γ = {gam:.3g} "
                       f"(3+1)", font=("Helvetica", 11, "bold"), fill=INK)
        series = [("clock rate e^φ", clock, C_PAST), ("cone speed c", c, C_FUT), ("event density ρ_c", rc, C_CHAIN)]
        lo = min(s.min() for _, s, _ in series)
        hi = max(s.max() for _, s, _ in series)
        lo, hi = min(lo, 0.6), max(hi, 1.6)
        top.vmin, top.vmax = lo, hi
        X1, Y1 = top.p(-12, 1.0)
        X2, _ = top.p(12, 1.0)
        cv.create_line(X1, Y1, X2, Y1, fill=GRID, dash=(3, 3))
        for k, (lab, s, col) in enumerate(series):
            pts = [q for a, b in zip(xs, s) for q in top.p(a, b)]
            cv.create_line(*pts, fill=col, width=2)
            cv.create_text(W - 40, 40 + 16 * k, anchor="e", text=lab, fill=col, font=("Helvetica", 10, "bold"))
        # little cones along x
        for xc in (-9, -4.5, 0, 4.5, 9):
            cc = float(np.interp(xc, xs, c))
            X, Y = top.p(xc, lo + 0.08 * (hi - lo))
            L = 26
            cv.create_line(X, Y, X - L * cc / math.sqrt(1 + cc * cc), Y - L / math.sqrt(1 + cc * cc), fill=INK2)
            cv.create_line(X, Y, X + L * cc / math.sqrt(1 + cc * cc), Y - L / math.sqrt(1 + cc * cc), fill=INK2)
        # lower panel: the ball
        bot = View(60, H * 0.52, W - 90, H * 0.42, -0.5, 20.5, -1.5, 5.5)
        cv.create_text(60, H * 0.50, anchor="w", text="A thrown ball = the longest chain of the order (time →, "
                       "height ↑; clocks run faster higher up)", font=("Helvetica", 11, "bold"), fill=INK)
        Xa, Ya = bot.p(0, -1)
        Xb, Yb = bot.p(20, -1)
        cv.create_line(Xa, Ya, Xb, Yb, fill=INK2)
        if self.ball:
            paths, tg, xg, g, T, rho, sec = self.ball
            t_all, x_all = paths[-1][2], paths[-1][3]
            for a, b in zip(t_all[::3], x_all[::3]):
                if not (bot.vmin <= b <= bot.vmax):
                    continue
                X, Y = bot.p(a, b)
                cv.create_oval(X - 1, Y - 1, X + 1, Y + 1, fill="#bdbab2", outline="")
            apex = []
            for tc, xc, _, _ in paths:
                pts = [q for a, b in zip(tc, xc) for q in bot.p(a, b)]
                cv.create_line(*pts, fill=C_CHAIN, width=1.6)
                apex.append(float(np.mean(xc[np.abs(tc - T / 2) < 3.0])))
            geo_mid = float(np.mean(xg[np.abs(tg - T / 2) < 3.0]))
            pts = [q for a, b in zip(tg[::20], xg[::20]) for q in bot.p(a, b)]
            cv.create_line(*pts, fill=INK, width=2, dash=(6, 4))
            Xs, Ys = bot.p(0, 0)
            Xe, Ye = bot.p(T, 0)
            cv.create_line(Xs, Ys, Xe, Ye, fill=INK2, dash=(2, 3))
            ap = float(np.mean(apex))
            se = float(np.std(apex, ddof=1) / math.sqrt(len(apex))) if len(apex) > 1 else float("nan")
            self.readout(f"{len(paths)} order(s) of ~{t_all.size} events (density {rho:g}): height of the longest "
                         f"chain at the top (|t − T/2| < 3) {ap:.2f}" + (f" ± {se:.2f}" if len(apex) > 1 else "") +
                         f";  the geodesic of general relativity (dashed) {geo_mid:.2f};  without gravity 0.  "
                         f"A single chain fluctuates by about ±0.8; at density 30 the chains also lie a little high (finite density). "
                         f"({sec:.1f} s)")
        else:
            self.readout("Choose a response above; press \"Throw a ball\" below.")


class Defect(Chapter):
    key = "defect"

    def __init__(self, app):
        super().__init__(app)
        self.run = None
        self.pair = None
        self.mode = "launched v = 0.6"
        self.rho = 400.0
        self.swerve = None
        self.busy = False

    def controls(self, fr):
        import tkinter as tk
        from tkinter import ttk
        ttk.Label(fr, text="Experiment").pack(anchor="w")
        self.mvar = tk.StringVar(value=self.mode)
        for m in ("at rest", "launched v = 0.3", "launched v = 0.6", "kicked (force f = 0.05, 1 ≤ t < 3)"):
            ttk.Radiobutton(fr, text=m, value=m, variable=self.mvar).pack(anchor="w")
        ttk.Label(fr, text="Density of events").pack(anchor="w", pady=(8, 0))
        self.rvar = tk.StringVar(value="400")
        for r in ("100", "400", "800"):
            ttk.Radiobutton(fr, text=r, value=r, variable=self.rvar).pack(anchor="w")
        ttk.Button(fr, text="Run", command=self.start).pack(fill="x", pady=(10, 2))
        ttk.Separator(fr).pack(fill="x", pady=8)
        ttk.Button(fr, text="Swerves: 8 orders at rest", command=self.swerves).pack(fill="x", pady=2)

    def start(self):
        if self.busy:
            return
        self.mode = self.mvar.get()
        self.rho = float(self.rvar.get())
        v = 0.3 if "0.3" in self.mode else 0.6 if "0.6" in self.mode else 0.0
        f = 0.05 if self.mode.startswith("kick") else 0.0
        X = 20.0 if v > 0.5 else 16.0
        self.run = KinkRun(self.rho, X, 10.0, v, self.rng, f=f)
        # the kick is read against the same kink at rest on the same events (paired, as in F3.14b): the swerves of
        # the order then cancel
        self.pair = KinkRun(self.rho, X, 10.0, 0.0, self.rng, events=(self.run.t, self.run.x)) if f else None
        self.swerve = None
        self.busy = True
        self.t0 = time.time()
        self.tick()

    def tick(self):
        if self.app.current is not self:
            self.busy = False
            return
        frac = self.run.step(3000)
        if self.pair is not None:
            frac = 0.5 * (frac + self.pair.step(3000))
        self.readout(f"computing the field on {self.run.t.size} events: {100 * frac:.0f} %")
        if frac < 1.0:
            self.app.root.after(1, self.tick)
        else:
            self.busy = False
            self.draw()

    def swerves(self):
        if self.busy:
            return
        self.busy = True
        rho, n = 100.0, 8
        vs = []
        tracks = []
        t0 = time.time()
        for k in range(n):
            self.readout(f"swerves: order {k + 1} of {n} ...")
            self.app.root.update()
            r = KinkRun(rho, 15.0, 10.0, 0.0, self.rng)
            r.step(10 ** 9)
            tr = r.track()
            tracks.append(tr)
            vs.append(fit_velocity(tr, 6.0, 10.0))
        self.busy = False
        self.swerve = (rho, tracks, np.array(vs), time.time() - t0)
        self.run = None
        self.draw()

    def draw(self):
        cv = self.cv
        cv.delete("all")
        W, H = self.app.CW, self.app.CH
        if self.swerve:
            rho, tracks, vs, sec = self.swerve
            span = max(0.5, max(float(np.max(np.abs(tr[:, 1]))) for tr in tracks) * 1.3)
            v = View(60, 30, W - 100, H - 80, -span, span, 0, 10)
            cv.create_text(60, 14, anchor="w", text=f"The centre of a kink at rest in 8 random orders (density {rho:g}; "
                           f"horizontal scale ±{span:.2f})",
                           font=("Helvetica", 11, "bold"), fill=INK)
            X0, Y0 = v.p(0, 0)
            X1, Y1 = v.p(0, 10)
            cv.create_line(X0, Y0, X1, Y1, fill=GRID, dash=(3, 3))
            for k, tr in enumerate(tracks):
                pts = [q for a, b in zip(tr[:, 1], tr[:, 0]) for q in v.p(a, b)]
                cv.create_line(*pts, fill=hue_color(k * 0.8), width=2)
            Dp = (32.0 / 15.0) / rho
            # velocity variance at t ~ 8 (the middle of the fit window) for a random walk in momentum, M = 8
            pred = math.sqrt(2.0 * Dp * 8.0) / 8.0
            self.readout(f"spread of the final velocities {np.std(vs, ddof=1):.3f} (8 orders);  shot noise of the order: "
                         f"D_p = (32/15)/ρ = {Dp:.4f} → σ_v(t≈8) ≈ √(2 D_p t)/M = {pred:.3f}   ({sec:.1f} s; a rough "
                         f"demonstration: F3.14c measures it with 1024 orders per cell)")
            return
        if not self.run or self.run.done < self.run.t.size:
            cv.create_text(W / 2, H / 2, text="Choose an experiment and press Run.", font=("Helvetica", 13), fill=INK2)
            return
        r = self.run
        v = View(50, 30, W - 80, H - 70, -r.X, r.X, 0, r.T)
        cv.create_text(50, 14, anchor="w", text=f"The field on {r.t.size} events (colour = angle; time ↑)   "
                       f"— {self.mode}", font=("Helvetica", 11, "bold"), fill=INK)
        step = max(1, r.t.size // 9000)
        for i in range(0, r.t.size, step):
            X, Y = v.p(r.x[i], r.t[i])
            cv.create_rectangle(X - 1.4, Y - 1.4, X + 1.4, Y + 1.4, fill=hue_color(r.phi[i]), outline="")
        tr = r.track()
        pts = [q for a, b in zip(tr[:, 1], tr[:, 0]) for q in v.p(a, b)]
        cv.create_line(*pts, fill="white", width=4)
        cv.create_line(*pts, fill=INK, width=1.5)
        vel = r.v
        g = 1 / math.sqrt(1 - vel * vel)
        if r.f:
            trp = self.pair.track()
            pts = [q for a, b in zip(trp[:, 1], trp[:, 0]) for q in v.p(a, b)]
            cv.create_line(*pts, fill=INK2, width=1.2, dash=(4, 3))
            va, vr = fit_velocity(tr, 4.0, 10.0), fit_velocity(trp, 4.0, 10.0)
            gv = va / math.sqrt(1 - va ** 2) - vr / math.sqrt(1 - vr ** 2)
            pred = 2 * math.pi * r.f * (r.t2 - r.t1) / 8.0
            self.readout(f"Δ(γv) kicked − at rest (dashed), same events: {gv:+.4f}; a body whose mass is its rest energy 8 "
                         f"gains {pred:.4f}  →  M_in/M_rest = {pred / gv:.3f}   (one order of {r.t.size} events; "
                         f"F3.14b with 128 orders at density 800: 0.9988 ± 0.0077)")
        else:
            vm = fit_velocity(tr, 1.0, 10.0)
            wid = float(np.median(tr[2:, 2]))
            self.readout(f"velocity {vm:+.3f} (launched {vel:g});  width {wid:.3f}: at rest a continuum kink has "
                         f"{width_pred(0.0):.3f}, Lorentz-contracted at this speed {width_pred(vel):.3f}   (one order of "
                         f"{r.t.size} events; at low density the order drags the kink a little: compare 100, 400, 800)")


class FieldEq(Chapter):
    key = "fieldeq"

    def draw(self):
        cv = self.cv
        cv.delete("all")
        W, H = self.app.CW, self.app.CH
        pw = (W - 100) / 3
        plots = []
        for i in range(3):
            plots.append(View(50 + i * (pw + 25), 60, pw - 10, H * 0.5, 0, 0.13, 0.76, 1.16 if i == 2 else 1.02))
        titles = ("the order: stiffness / EH (F6.8g, kT = 4)", "matter: source / T_ab (F6.10b)",
                  "together: response / GR (derived)")
        for v, tt in zip(plots, titles):
            cv.create_rectangle(v.x0, v.y0, v.x0 + v.w, v.y0 + v.h, outline=GRID)
            cv.create_text(v.x0, v.y0 - 12, anchor="w", text=tt, font=("Helvetica", 9, "bold"), fill=INK)
            X1, Y1 = v.p(0, 1.0)
            X2, _ = v.p(0.13, 1.0)
            cv.create_line(X1, Y1, X2, Y1, fill=INK2, dash=(2, 3))
            for kk in (0.0, 0.05, 0.10):
                X, Y = v.p(kk, v.vmin)
                cv.create_text(X, Y + 10, text=f"{kk:g}", font=("Helvetica", 8), fill=INK2)
            cv.create_text(v.x0 + v.w / 2, v.y0 + v.h + 24, text="k ℓ_k", font=("Helvetica", 9), fill=INK2)
            for yy in (0.8, 0.9, 1.0, 1.1):
                if v.vmin <= yy <= v.vmax:
                    X, Y = v.p(0, yy)
                    cv.create_text(X - 4, Y, anchor="e", text=f"{yy:.1f}", font=("Helvetica", 8), fill=INK2)
        kk = np.linspace(1e-4, 0.13, 100)
        v = plots[0]
        for key, fk, col, lab in (("Kc", "fit_c", C_PAST, "conformal"), ("Kab", "fit_ab", C_FUT, "mixed"),
                                  ("Kbb", "fit_bb", C_CHAIN, "optical")):
            pts = [q for a, b in zip(kk, kfit(kk, *F68G[fk])) for q in v.p(a, b)]
            cv.create_line(*pts, fill=col, width=1.5)
            for a, b in zip(F68G["k"], F68G[key]):
                X, Y = v.p(a, b)
                cv.create_oval(X - 3, Y - 3, X + 3, Y + 3, fill=col, outline="white")
        for k, (lab, col) in enumerate((("conformal A/6", C_PAST), ("mixed B/8", C_FUT), ("optical C/2", C_CHAIN))):
            cv.create_text(v.x0 + 6, v.y0 + v.h - 44 + 14 * k, anchor="w", text=lab, fill=col, font=("Helvetica", 8, "bold"))
        v = plots[1]
        for key, col, lab in (("Fw", C_PAST, "F_ω (trace)"), ("Fn", C_CHAIN, "F_ν (Σ T_ii)")):
            pts = [q for a, b in zip(F610B["k"], F610B[key]) for q in v.p(a, b)]
            cv.create_line(*pts, fill=col, width=1.5)
            for a, b in zip(F610B["k"], F610B[key]):
                X, Y = v.p(a, b)
                cv.create_oval(X - 3, Y - 3, X + 3, Y + 3, fill=col, outline="white")
        cv.create_text(v.x0 + 6, v.y0 + v.h - 30, anchor="w", text="F_ω (trace)", fill=C_PAST, font=("Helvetica", 8, "bold"))
        cv.create_text(v.x0 + 6, v.y0 + v.h - 16, anchor="w", text="F_ν (Σ T_ii, cones)", fill=C_CHAIN,
                       font=("Helvetica", 8, "bold"))
        v = plots[2]
        ks = [0.0, 0.06, 0.12]
        res = [ov3_response(1, 1, 1, 1, 1)]
        for k in (0.06, 0.12):
            i = F68G["k"].index(k)
            j = F610B["k"].index(k)
            res.append(ov3_response(F68G["Kc"][i], F68G["Kab"][i], F68G["Kbb"][i], F610B["Fw"][j], F610B["Fn"][j]))
        series = (("dust Φ/Φ_GR", [r["dust"]["Phi/GR"] for r in res], C_PAST),
                  ("dust γ(k)", [r["dust"]["gamma"] for r in res], C_CHAIN),
                  ("−wall/dust (repels)", [-r["wall/dust"] for r in res], C_FUT))
        for k, (lab, ys, col) in enumerate(series):
            pts = [q for a, b in zip(ks, ys) for q in v.p(a, b)]
            cv.create_line(*pts, fill=col, width=1.5)
            for a, b in zip(ks, ys):
                X, Y = v.p(a, b)
                cv.create_oval(X - 3, Y - 3, X + 3, Y + 3, fill=col, outline="white")
            cv.create_text(v.x0 + 6, v.y0 + v.h - 44 + 14 * k, anchor="w", text=lab, fill=col, font=("Helvetica", 8, "bold"))
        y = 60 + H * 0.5 + 60
        cv.create_text(W / 2, y, text="G_ab + 2(η_ab□ − ∂_a∂_b) F₀(□) R  =  8πG T_ab^eff", font=("Helvetica", 15, "bold"),
                       fill=INK)
        cv.create_text(W / 2, y + 26, text="the linearised field equation of the order (young static field); "
                       "F₀ and the source factors are counted from intervals", font=("Helvetica", 10), fill=INK2)
        cv.create_text(W / 2, y + 56, text="open: second order (β) · old fields · time-dependent and tensor sectors · "
                       "the value of G", font=("Helvetica", 10, "italic"), fill=C_FUT)
        r12 = res[2]
        self.readout(f"kℓ_k = 0.12: dust Φ/Φ_GR {r12['dust']['Phi/GR']:.3f}, γ {r12['dust']['gamma']:.3f}, wall/dust "
                     f"{r12['wall/dust']:+.3f}  (OV3-lin-a, from the frozen records of F6.8g and F6.10b)")


class Status(Chapter):
    key = "status"

    def draw(self):
        cv = self.cv
        cv.delete("all")
        W, H = self.app.CW, self.app.CH
        y = 24
        cv.create_text(20, y, anchor="w", text="claim", font=("Helvetica", 11, "bold"), fill=INK)
        cv.create_text(W * 0.55, y, anchor="w", text="status", font=("Helvetica", 11, "bold"), fill=INK)
        cv.create_text(W * 0.78, y, anchor="w", text="laboratory", font=("Helvetica", 11, "bold"), fill=INK)
        dy = (H - 60) / len(STATUS_ROWS)
        for k, (claim, st, lab) in enumerate(STATUS_ROWS):
            y = 50 + dy / 2 + k * dy
            if k % 2 == 0:
                cv.create_rectangle(10, y - dy / 2, W - 10, y + dy / 2, fill=PANEL, outline="")
            col = (C_CHAIN if st.startswith("tested") else C_PAST if st.startswith("derived") else
                   INK2 if st.startswith(("adopted", "interpretation")) else C_FUT)
            cv.create_text(20, y, anchor="w", text=claim, font=("Helvetica", 10), fill=INK)
            cv.create_text(W * 0.55, y, anchor="w", text=st, font=("Helvetica", 10, "bold"), fill=col)
            cv.create_text(W * 0.78, y, anchor="w", text=lab, font=("Helvetica", 9), fill=INK2)
        self.readout("tested = a frozen laboratory test passed · derived = follows from tested or adopted pieces · "
                     "adopted = taken from known physics · open")


CHAPTERS = [("Overview", Intro), ("1  The world grows", Grow), ("2  The order of events", Order),
            ("3  Gravity as a response", Gravity), ("4  Matter as a defect", Defect),
            ("5  The field equation", FieldEq), ("6  What is known", Status)]


class App:
    CW, CH = 820, 620

    def __init__(self, start=0):
        import tkinter as tk
        from tkinter import ttk
        self.tk = tk
        self.root = tk.Tk()
        self.root.title("Causal Physics — an interactive demonstration")
        self.root.configure(bg=BG)
        style = ttk.Style()
        try:
            style.theme_use("clam")
        except tk.TclError:
            pass
        left = tk.Frame(self.root, bg=PANEL, padx=8, pady=8)
        left.grid(row=0, column=0, rowspan=2, sticky="ns")
        tk.Label(left, text="Causal Physics", font=("Helvetica", 15, "bold"), bg=PANEL, fg=INK).pack(anchor="w")
        tk.Label(left, text="events · order · gravity", font=("Helvetica", 9), bg=PANEL, fg=INK2).pack(anchor="w",
                                                                                                    pady=(0, 10))
        self.chap_buttons = []
        for k, (name, _) in enumerate(CHAPTERS):
            b = tk.Button(left, text=name, anchor="w", relief="flat", bg=PANEL, fg=INK, activebackground=GRID,
                          font=("Helvetica", 11), command=lambda k=k: self.show(k))
            b.pack(fill="x", pady=1)
            self.chap_buttons.append(b)
        tk.Label(left, text="Controls", font=("Helvetica", 11, "bold"), bg=PANEL, fg=INK).pack(anchor="w", pady=(16, 2))
        self.ctrl = tk.Frame(left, bg=PANEL)
        self.ctrl.pack(fill="x")
        self.canvas = tk.Canvas(self.root, width=self.CW, height=self.CH, bg=BG, highlightthickness=0)
        self.canvas.grid(row=0, column=1, sticky="nsew", padx=6, pady=6)
        self.canvas.bind("<Button-1>", lambda e: self.current.click(e.x, e.y))
        self.readout = tk.StringVar()
        tk.Label(self.root, textvariable=self.readout, bg=BG, fg=INK, font=("Helvetica", 10), anchor="w",
                 justify="left", wraplength=self.CW).grid(row=1, column=1, sticky="we", padx=8, pady=(0, 8))
        right = tk.Frame(self.root, bg=BG)
        right.grid(row=0, column=2, rowspan=2, sticky="ns", padx=(0, 8), pady=6)
        self.title = tk.Label(right, text="", font=("Helvetica", 14, "bold"), bg=BG, fg=INK, anchor="w")
        self.title.pack(anchor="w")
        nb = ttk.Notebook(right)
        nb.pack(fill="both", expand=True)
        self.txt_plain = tk.Text(nb, width=46, wrap="word", font=("Helvetica", 11), bg=BG, fg=INK, relief="flat",
                                 padx=8, pady=8)
        self.txt_phys = tk.Text(nb, width=46, wrap="word", font=("Helvetica", 10), bg=BG, fg=INK, relief="flat",
                                padx=8, pady=8)
        nb.add(self.txt_plain, text="In plain words")
        nb.add(self.txt_phys, text="For physicists")
        self.chapters = [None] * len(CHAPTERS)
        self.current = None
        self.show(start)

    def show(self, k):
        if self.current is not None:
            self.current.stop()
        for i, b in enumerate(self.chap_buttons):
            b.configure(bg=GRID if i == k else PANEL, font=("Helvetica", 11, "bold" if i == k else "normal"))
        if self.chapters[k] is None:
            self.chapters[k] = CHAPTERS[k][1](self)
        ch = self.chapters[k]
        self.current = ch
        for w in self.ctrl.winfo_children():
            w.destroy()
        ch.controls(self.ctrl)
        title, plain, phys = TEXTS[ch.key]
        self.title.configure(text=title)
        for box, s in ((self.txt_plain, plain), (self.txt_phys, phys)):
            box.configure(state="normal")
            box.delete("1.0", "end")
            box.insert("1.0", s)
            box.configure(state="disabled")
        ch.draw()

    def run(self):
        self.root.mainloop()


# ===================================================================================================== selftest

def selftest():
    ok = True
    rng = np.random.default_rng(3)

    def check(name, cond, info):
        nonlocal ok
        ok &= bool(cond)
        print(f"{'PASS' if cond else 'FAIL'}  {name}: {info}")

    t, x, _ = sprinkle_diamond(40.0, 10.0, rng)
    R = relation(t, x)
    d, f = mm_dimension(R)
    check("dimension 1+1", abs(d - 2.0) < 0.15, f"d = {d:.3f} (N {t.size})")
    t3, x3, y3 = sprinkle_diamond(3.0, 10.0, rng, 3)
    d3, _ = mm_dimension(relation(t3, x3, y3))
    check("dimension 2+1", abs(d3 - 3.0) < 0.3, f"d = {d3:.3f} (N {t3.size})")
    a = int(np.argmin(t))
    b = int(np.argmax(t))
    ch = longest_chain(R, a, b, t)
    check("chain is a chain", all(R[i, j] for i, j in zip(ch[:-1], ch[1:])), f"length {len(ch) - 1}")
    psi, c, rc, gam = responses(np.array([-0.1]), "throughput")
    check("throughput: rho_c c = 1", abs(rc[0] * c[0] - 1.0) < 1e-12 and gam == 1.0, f"γ {gam}")
    tg, xg = ball_geodesic(0.08, 20.0)
    check("ball geodesic ends at 0", abs(xg[-1]) < 1e-6, f"apex {xg.max():.3f} (Newton {0.08 * 400 / 8:.3f}; fast throw)")
    # exact kink: with the rule's continuum limit the free wave minus the integral reproduces the kink; check phi_0
    tt = np.array([0.0, 0.0])
    xx = np.array([-1.0, 1.0])
    check("free wave at t = 0 = kink", np.allclose(kink_free(tt, xx, 0.3), 4 * np.arctan(np.exp(xx / math.sqrt(1 - 0.09)))),
          "d'Alembert initial data")
    r = KinkRun(100.0, 16.0, 6.0, 0.0, rng)
    r.step(10 ** 9)
    tr = r.track()
    drift = fit_velocity(tr, 1.0, 6.0)
    check("kink at rest stays", abs(drift) < 0.05 and abs(tr[-1, 1]) < 0.5, f"velocity {drift:+.4f}, centre {tr[-1, 1]:+.3f}")
    # the Fenwick dominance sum against a brute-force causal sum on a small order
    r2 = KinkRun(20.0, 6.0, 3.0, 0.3, rng)
    r2.step(10 ** 9)
    i = int(np.argmax(r2.t))
    past = (r2.t[i] - r2.t > np.abs(r2.x[i] - r2.x))
    brute = r2.phi0[i] - 0.5 / r2.rho * np.sin(r2.phi[past]).sum()
    check("Fenwick sum = causal sum", abs(brute - r2.phi[i]) < 1e-9, f"diff {abs(brute - r2.phi[i]):.1e}")
    resp = ov3_response(0.95245, 0.95241, 0.93697, 0.95447, 0.94949)
    check("OV3 at k 0.06", abs(resp["dust"]["Phi/GR"] - 0.9403) < 2e-4 and abs(resp["dust"]["gamma"] - 1.0330) < 2e-4
          and abs(resp["wall/dust"] + 1.0447) < 2e-4,
          f"Φ/GR {resp['dust']['Phi/GR']:.4f}, γ {resp['dust']['gamma']:.4f}, wall/dust {resp['wall/dust']:+.4f}")
    g0 = ov3_response(1, 1, 1, 1, 1)
    check("OV3 local limit = GR", abs(g0["dust"]["gamma"] - 1) < 1e-12 and abs(g0["wall/dust"] + 1) < 1e-12, "γ 1, wall −1")
    print("selftest", "PASSED" if ok else "FAILED")
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser(description="Causal Physics: an interactive demonstration")
    ap.add_argument("--chapter", type=int, default=0, help="0 overview, 1-6 the chapters")
    ap.add_argument("--selftest", action="store_true")
    q = ap.parse_args()
    if q.selftest:
        return selftest()
    App(max(0, min(q.chapter, len(CHAPTERS) - 1))).run()
    return 0


if __name__ == "__main__":
    sys.exit(main())
