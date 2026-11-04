"""
core.operators.bd_expect_curved
===============================

The EXPECTED smeared Benincasa-Dowker operator of a causal set sprinkled
into the conformally flat static geometry g = Pi^2 eta, and the
self-consistent field equation of the author's choice (a), 18 Oct 2026:

    B_g Pi = kappa rho_s          (the operator of the causal set in the
                                   geometry the field itself creates, acting
                                   directly on Pi)

Why this is computable by quadrature
------------------------------------
A conformal factor does not change the causal order: the causal set of
g = Pi^2 eta has the order of flat space. Only the density changes: the
sprinkling is Poisson with density rho per unit PROPER volume, i.e.
rho Pi(y)^D per unit flat volume. The number of elements in the interval
[y, x] is therefore Poisson with mean rho V_g(x, y), where

    V_g(x, y) = integral over the flat causal diamond of Pi(z)^D dV_flat(z)

is the proper volume of the interval. The expectation of the smeared
operator is then (compare core.operators.bd_expect, the flat case)

    E[B_g f](x) = -a_d (eps rho)^(2/d) [ f(x) + rho eps
                  int_{J^-(x)} f(y) gk(eps rho V_g(x, y)) Pi(y)^D dV_flat(y) ]

with the same kernel gk(w) = sum_i C_i w^i / i! e^(-w) and the FUNDAMENTAL
density rho in the prefactor (the discreteness scale is set per proper
volume). With Pi = 1 everything reduces exactly to the flat operator.

The proper volume of the diamond (the only new ingredient)
----------------------------------------------------------
V_g = V_flat * Phi, Phi = <Pi^D> averaged over the flat diamond. A static
field needs only the SPATIAL projection of the diamond, and that is exact
and simple: between y and x (spatial distance L, time difference dt) the
diamond projects onto the prolate spheroid d1 + d2 <= dt with foci at the
two end points, with density l(z) = dt - (d1 + d2) (the length of the time
interval above the point z). diamond_phi_spheroid integrates Pi^D over it
in prolate spheroidal coordinates (Gauss-Legendre across, Gauss-Legendre or
periodic trapezoid along and around). This holds for every diamond, from
near-null slivers along the light cone (thin spheroids) to timelike ones
(near-spherical). The calibration test K.11 compares Phi with a Monte Carlo
integral over the diamond sampled in its rest frame.

Panels in s. With Pi < 1 near a source, Phi < 1 and the kernel at a given
flat s is evaluated at a SMALLER w = w_flat Phi; the w panels are therefore
extended (doubling beyond 80) until w_flat Phi_min >= 80, with Phi_min =
min(Pi)^D a lower bound of Phi.

Continuum candidates (local limits, for comparison only)
--------------------------------------------------------
Averaged causal-set d'Alembertians tend, in curved spacetime, to
box_g - R/2 in every dimension (Belenchia, Benincasa, Liberati, Marin,
Marino, Ormani 2016, "Universal behaviour of generalized causal set
d'Alembertians in curved spacetime"). For g = Pi^2 eta and a static Pi
(signature -+++, box static = +lap, as in this laboratory):

    box_g Pi   = Pi^-2 lap Pi + (D - 2) Pi^-3 |grad Pi|^2
    -R/2 Pi    = (D - 1) Pi^-2 lap Pi + (D - 1)(D - 4)/2 Pi^-3 |grad Pi|^2

so that the local form of (a) is

    alpha Pi^-2 lap Pi + c Pi^-3 |grad Pi|^2 = kappa rho_s,
    alpha = D, c = (D - 2) + (D - 1)(D - 4)/2      (c = 0 for D = 3, 2 for D = 4)

Every equation of this form is equivalent to lap Pi^p = (p / alpha) kappa
Pi^(1+p) rho_s with p = 1 + c / alpha: outside the source Pi^p is harmonic.
The candidates used by F1.4 and F1.5:
    "a"      (a) with -R/2:             alpha = D, c = c_D
    "box_g"  box_g alone (no -R/2):     alpha = 1, c = D - 2
    "flat"   no back-reaction:          lap Pi = kappa rho_s (linear)

Functions
---------
    curved_radial_rows(D, rho, eps, t_x, r_rows, R, nodes, Pi_nodes, ..., geom)
            rows of E[B_g] on the radial hat basis (Pi_nodes on `nodes`);
            geom = 'linear' or 'spline' (C^2 geometry, needed in 3+1, K.11c)
    pi_geometry(nodes, Pi_nodes, geom)
            the interpolated geometry Pi(r)
    diamond_phi_spheroid(Pi_of_r, D, r_x, RS, CO, dt, tau, n_mu, n_nu, n_phi)
            Phi = <Pi^D> over flat diamonds (exact spatial projection, vectorised)
    continuum_radial(model, D, R, source, kappa, n_grid)
            local candidate solution Pi(r) on a fine grid
    solve_selfconsistent(...)
            Pi solving E[B_g(Pi)] Pi = rhs with Pi(R) = 1 (Anderson-accelerated
            fixed point; the rows are rebuilt from the current Pi each step)
"""

from __future__ import annotations

import math
import time
from typing import Callable, Dict, Optional, Sequence

import numpy as np

from .bd import BD_NORM, hat_weights
from .bd_expect import _directions, _wall_distance, _S_GEOM, _W_EDGES, g_kernel
from ..order.sprinkle import interval_volume


# =============================================================================
# geometry of the past of x and the diamond average
# =============================================================================

def _s_edges_extended(D: int, phi_min: float) -> np.ndarray:
    """Panel edges in s = w_flat^(2/d): the flat edges, extended by doublings until w_flat * phi_min >= 80."""
    w = list(_W_EDGES)
    while w[-1] * phi_min < _W_EDGES[-1]:
        w.append(w[-1] * 2.0)
    top = np.array(w) ** (2.0 / D)
    geo = [e for e in _S_GEOM if e < top[0]]
    return np.array([0.0] + geo + list(top))


def diamond_phi_spheroid(Pi_of_r: Callable[[np.ndarray], np.ndarray], D: int, r_x: float, RS: np.ndarray,
                         CO: np.ndarray, dt: np.ndarray, tau: np.ndarray, n_mu: int = 4, n_nu: int = 12,
                         n_phi: int = 4) -> np.ndarray:
    """
    Phi = <Pi^D> over the flat diamond, using the EXACT spatial projection of
    the diamond (a static field needs nothing else). Between y and x
    (spatial distance L = RS, time difference dt) the diamond projects onto
    the prolate spheroid d1 + d2 <= dt (d1, d2 the distances to the two end
    points) with the density l(z) = dt - (d1 + d2), the length of the time
    interval above z. Prolate spheroidal coordinates with foci at the ends,
    c = L / 2, sigma = (d1 + d2) / L in [1, dt / L], te = (d1 - d2)/L in [-1, 1]:
        n = 3:  dV = c^3 (sigma^2 - te^2) dsigma dte dphi     (Gauss-Legendre, trapezoid in phi)
        n = 2:  sigma = cosh(mu), te = cos(nu), dV = c^2 (cosh^2 mu - cos^2 nu) dmu dnu,
                nu in [0, 2 pi) (both sides; trapezoid, periodic); sigma = 1 + (smax - 1) x^2
                with Gauss-Legendre in x removes the 1/sqrt(sigma^2 - 1) of dmu
    position: along the axis xi = c (1 + sigma te) from y, across rho = c sqrt((sigma^2 - 1)(1 - te^2)).
    """
    n = D - 1
    L = np.broadcast_to(RS, np.broadcast(RS, dt).shape)
    c = 0.5 * L
    smax = np.maximum(dt / np.maximum(L, 1e-300), 1.0)
    SN = np.sqrt(np.maximum(1.0 - CO * CO, 0.0))
    xg, wg = np.polynomial.legendre.leggauss(n_mu)
    xg, wg = 0.5 * (xg + 1.0), 0.5 * wg
    num = np.zeros(L.shape)
    den = np.zeros(L.shape)

    def add(xi, rperp_e1, rperp2, weight):
        alpha = L - xi                                  # offset from x along n (towards y)
        r2 = r_x * r_x + 2.0 * r_x * (alpha * CO + rperp_e1) + alpha * alpha + rperp2
        return weight * Pi_of_r(np.sqrt(np.maximum(r2, 0.0))) ** D, weight

    if n == 3:
        tg, wt = np.polynomial.legendre.leggauss(n_nu)
        phis = 2.0 * np.pi * (np.arange(n_phi) + 0.5) / n_phi
        for x, wx in zip(xg, wg):
            sig = 1.0 + (smax - 1.0) * x
            jac = (smax - 1.0) * wx
            ell = dt - L * sig
            for te, wte in zip(tg, wt):
                xi = c * (1.0 + sig * te)
                rp = c * np.sqrt(np.maximum((sig * sig - 1.0) * (1.0 - te * te), 0.0))
                wbase = jac * wte * ell * c ** 3 * (sig * sig - te * te) / n_phi
                for ph in phis:
                    a, b = add(xi, rp * SN * np.cos(ph), rp * rp, wbase)
                    num += a
                    den += b
    elif n == 2:
        # sigma = cosh(mu) = 1 + (smax - 1) x^2 makes dmu = 2 sqrt(smax - 1) dx / sqrt(sigma + 1) smooth in x
        nus = 2.0 * np.pi * (np.arange(n_nu) + 0.5) / n_nu
        for x, wx in zip(xg, wg):
            ch = 1.0 + (smax - 1.0) * x * x
            sh = np.sqrt(np.maximum(ch * ch - 1.0, 0.0))
            jac = wx * 2.0 * np.sqrt(smax - 1.0) / np.sqrt(ch + 1.0)
            ell = dt - L * ch
            for nu in nus:
                xi = c * (1.0 + ch * np.cos(nu))
                rp = c * sh * np.sin(nu)                 # signed: both sides of the axis
                wbase = jac * ell * c * c * (ch * ch - np.cos(nu) ** 2) * (2.0 * np.pi / n_nu)
                a, b = add(xi, rp * SN, rp * rp, wbase)
                num += a
                den += b
    else:
        raise ValueError("D must be 3 or 4")
    return num / np.where(den > 0, den, 1.0)


def _pi_interp(nodes: np.ndarray, Pi_nodes: np.ndarray) -> Callable[[np.ndarray], np.ndarray]:
    """Piecewise-linear Pi(r) on the nodes, Pi = 1 beyond the last node (vacuum outside the region)."""
    nodes = np.asarray(nodes, float)
    Pi_nodes = np.asarray(Pi_nodes, float)
    return lambda r: np.interp(r, nodes, Pi_nodes, right=1.0)


def _pi_spline(nodes: np.ndarray, Pi_nodes: np.ndarray) -> Callable[[np.ndarray], np.ndarray]:
    """
    C^2 cubic spline Pi(r) through the nodal values (K.11c): Pi'(0) = 0
    (a smooth radial field), not-a-knot at the outer node, continued as the
    same cubic beyond it (only the thin edges of diamonds next to the wall
    reach there). Linear in the nodal values; a constant field is
    reproduced exactly (so Pi = 1 still gives the flat operator).
    Why C^2: the curved operator responds to the CURVATURE of the geometry
    (its local limit carries -R/2, i.e. lap Pi): long, nearly null diamonds
    are thin needles whose proper volume changes with their thickness in
    proportion to the transverse second derivative of Pi^D. A piecewise-
    linear geometry has that second derivative as delta functions at the
    nodes; in 3+1 the far rows then carry a spurious residue larger than the
    source that converges only slowly in the quadrature (K.11b far_d4).
    A C^1 interpolant (PCHIP, jumps of the second derivative) is not enough
    (development runs of 19 Oct, recorded in K.11c).
    """
    from scipy.interpolate import CubicSpline
    nodes = np.asarray(nodes, float)
    Pi_nodes = np.asarray(Pi_nodes, float)
    cs = CubicSpline(nodes, Pi_nodes, bc_type=((1, 0.0), "not-a-knot"), extrapolate=True)
    return lambda r: cs(np.asarray(r, float))


def pi_geometry(nodes: np.ndarray, Pi_nodes: np.ndarray, geom: str = "linear") -> Callable[[np.ndarray], np.ndarray]:
    """The geometry Pi(r) the rows are built in: 'linear' (K.11, K.11b, F1.4, F1.5) or 'spline' (C^2, K.11c)."""
    if geom == "linear":
        return _pi_interp(nodes, Pi_nodes)
    if geom == "spline":
        return _pi_spline(nodes, Pi_nodes)
    raise ValueError(f"unknown geometry interpolation {geom!r}")


# =============================================================================
# rows of the expected operator in g = Pi^2 eta
# =============================================================================

def curved_radial_rows(D: int, rho: float, eps: float, t_x: float, r_rows: Sequence[float], R: float,
                       nodes: np.ndarray, Pi_nodes: np.ndarray, n_rho: int = 160, n_s: int = 8,
                       n_ang: int = 16, n_mu: int = 4, n_nu: int = 12, n_phi: int = 4,
                       geom: str = "linear", geom_nodes: Optional[np.ndarray] = None,
                       geom_values: Optional[np.ndarray] = None) -> np.ndarray:
    """
    Rows A of E[B_g] on the radial hat basis, E[B_g f](r_i) = sum_k A[i, k] f_k
    for a static radial f piecewise linear on `nodes`, in the geometry
    g = Pi^2 eta given by Pi_nodes (same nodes). geom: how the geometry is
    interpolated between the nodes, 'linear' (Pi = 1 beyond nodes[-1] >= R;
    the instrument of K.11, K.11b, F1.4, F1.5) or 'spline' (C^2 cubic
    spline, K.11c). The field ACTED ON stays in the hat basis in both cases
    (pointwise values, as in the flat operator of K.4); only the geometry,
    whose curvature the operator responds to, needs to be smooth.
    geom_nodes, geom_values: optional separate description of the geometry
    (e.g. a continuum solution on a fine grid, the reference of K.11c);
    default: the geometry is the field itself (nodes, Pi_nodes).
    With Pi_nodes = 1 this equals core.operators.bd_expect.expected_radial_rows
    with the same quadrature parameters for either geom (core check
    `curved_rows_flat_limit`).
    """
    n = D - 1
    nodes = np.asarray(nodes, float)
    Pi_nodes = np.asarray(Pi_nodes, float)
    gn = nodes if geom_nodes is None else np.asarray(geom_nodes, float)
    gv = Pi_nodes if geom_values is None else np.asarray(geom_values, float)
    Pi = pi_geometry(gn, gv, geom)
    kappa_D = interval_volume(D, 1.0)
    lam = eps * rho * kappa_D                    # w_flat = lam tau^d
    l2 = lam ** (2.0 / D)
    pref = -BD_NORM[D] * (eps * rho) ** (2.0 / D)
    # lower bound of Phi for the panel extension (a spline may dip slightly below the smallest nodal value)
    pi_lo = float(np.min(gv)) if geom == "linear" else \
        float(np.min(Pi(np.linspace(0.0, float(gn[-1]), 20 * gn.size))))
    phi_min = max(pi_lo, 1e-6) ** D
    edges = _s_edges_extended(D, min(1.0, phi_min))
    cos, wa = _directions(n, n_ang)
    z, wz = np.polynomial.legendre.leggauss(n_rho)
    z, wz = 0.5 * (z + 1.0), 0.5 * wz
    u, wu = np.polynomial.legendre.leggauss(n_s)
    u, wu = 0.5 * (u + 1.0), 0.5 * wu
    K = len(nodes)
    A = np.zeros((len(r_rows), K))
    for i, r_x in enumerate(r_rows):
        r_x = float(r_x)
        row = np.zeros(K + 1)
        for a in range(cos.size):                # one direction at a time (memory)
            CO = cos[a]
            rho_lim = min(t_x, float(_wall_distance(r_x, np.array([CO]), R)[0]))
            rs = rho_lim * z ** 2                # (P,)
            wr = wz * 2.0 * rho_lim * z
            s_max = l2 * np.maximum(t_x * t_x - rs * rs, 0.0)
            S, WS = [], []
            for lo, hi_e in zip(edges[:-1], edges[1:]):
                hi = np.minimum(hi_e, s_max)[:, None]
                span = np.maximum(hi - lo, 0.0)
                S.append(lo + span * u)
                WS.append(span * wu)
            s = np.concatenate(S, axis=1)        # (P, Q)
            ws = np.concatenate(WS, axis=1)
            RS = rs[:, None]
            tau2 = s / l2
            tau = np.sqrt(tau2)
            dt = np.sqrt(RS * RS + tau2)
            r_y = np.sqrt(np.maximum(r_x * r_x + RS * RS + 2.0 * r_x * RS * CO, 0.0))
            radial = RS ** (n - 1) if n > 1 else np.ones_like(RS)
            meas = wa[a] * wr[:, None] * ws * radial / (l2 * 2.0 * dt)
            phi = diamond_phi_spheroid(Pi, D, r_x, RS, CO, dt, tau, n_mu, n_nu, n_phi)
            Wt = meas * g_kernel((s ** (D / 2.0)) * phi, D) * Pi(r_y) ** D
            rr = np.broadcast_to(np.minimum(r_y, nodes[-1]), Wt.shape).ravel()
            k, w0, w1 = hat_weights(rr, nodes)
            Wf = Wt.ravel()
            row += np.bincount(k, weights=Wf * w0, minlength=K + 1)[:K + 1]
            row += np.bincount(k + 1, weights=Wf * w1, minlength=K + 1)[:K + 1]
        kx, a0, a1 = hat_weights(np.array([min(r_x, nodes[-1])]), nodes)
        own = np.zeros(K)
        own[kx[0]] += a0[0]
        if kx[0] + 1 < K:
            own[kx[0] + 1] += a1[0]
        A[i] = pref * (own + rho * eps * row[:K])
    return A


def rows_task(D, rho, eps, t_x, r_rows, R, nodes, Pi_nodes, quad, geom="linear"):
    """Module-level wrapper for worker processes: curved rows for a chunk of row radii."""
    return curved_radial_rows(D, rho, eps, t_x, r_rows, R, nodes, Pi_nodes, n_rho=int(quad[0]), n_s=int(quad[1]),
                              n_ang=int(quad[2]), n_mu=int(quad[3]), n_nu=int(quad[4]), n_phi=int(quad[5]),
                              geom=geom)


# =============================================================================
# continuum candidates (local limits)
# =============================================================================

def model_constants(model: str, D: int):
    """(alpha, c, p) of the local equation alpha Pi^-2 lap Pi + c Pi^-3 |grad Pi|^2 = kappa rho_s; 'flat' -> None."""
    if model == "a":
        alpha, c = float(D), (D - 2) + 0.5 * (D - 1) * (D - 4)
    elif model == "box_g":
        alpha, c = 1.0, float(D - 2)
    elif model == "flat":
        return None
    else:
        raise ValueError(model)
    return alpha, float(c), 1.0 + float(c) / alpha


def omega_n(n: int) -> float:
    """Surface area of the unit sphere in n dimensions."""
    return 2.0 * math.pi ** (n / 2.0) / math.gamma(n / 2.0)


def _poisson_radial(n: int, r: np.ndarray, f: np.ndarray) -> np.ndarray:
    """Psi with lap Psi = f (radial, n dims), Psi'(0) = 0, Psi(R) = 0, by two trapezoid integrations (Gauss law)."""
    dens = f * omega_n(n) * r ** (n - 1)
    Q = np.concatenate([[0.0], np.cumsum(0.5 * (dens[1:] + dens[:-1]) * np.diff(r))])
    du = np.zeros_like(r)
    du[1:] = Q[1:] / (omega_n(n) * r[1:] ** (n - 1))
    U = np.concatenate([[0.0], np.cumsum(0.5 * (du[1:] + du[:-1]) * np.diff(r))])
    return U - U[-1]


def continuum_radial(model: str, D: int, R: float, source: Callable[[np.ndarray], np.ndarray], kappa: float = 1.0,
                     n_grid: int = 40001, tol: float = 1e-12, max_iter: int = 2000) -> Dict[str, np.ndarray]:
    """
    Local candidate solution Pi(r) with Pi(R) = 1 on a uniform grid:
      flat:  lap Pi = kappa rho_s;
      a, box_g:  lap Psi = (p / alpha) kappa Pi^(1+p) rho_s, Psi = Pi^p (relaxed fixed point).
    Raises if Pi reaches 0 (no static solution of that strength).
    """
    n = D - 1
    r = np.linspace(0.0, R, n_grid)
    rs = source(r)
    if model == "flat":
        return {"r": r, "Pi": 1.0 + _poisson_radial(n, r, kappa * rs), "iterations": 1}
    alpha, c, p = model_constants(model, D)
    Pi = np.ones_like(r)
    it_total = 0
    # continuation in the source strength (t -> 1): each stage starts from the previous solution, so a strong
    # source is reached through solutions that exist; a stage that drives Pi^p to 0 means no static solution
    for t in (0.125, 0.25, 0.5, 0.75, 1.0):
        for it in range(1, max_iter + 1):
            psi = 1.0 + _poisson_radial(n, r, t * (p / alpha) * kappa * Pi ** (1.0 + p) * rs)
            if np.min(psi) <= 0.0:
                raise ValueError(f"continuum {model}: Pi^p reaches 0 (source too strong)")
            new = psi ** (1.0 / p)
            step = float(np.max(np.abs(new - Pi)))
            Pi = 0.5 * Pi + 0.5 * new
            it_total += 1
            if step < tol:
                break
    return {"r": r, "Pi": Pi, "iterations": it_total}


# =============================================================================
# the self-consistent causal-set field of (a)
# =============================================================================

def solve_selfconsistent(D: int, R: float, nodes: np.ndarray, rows: np.ndarray, rhs: np.ndarray, quad: Sequence[int],
                         parallel, workers: int, progress: Callable[[str], None], label: str = "",
                         rho: float = 20.0, eps: float = 0.05, tol: float = 1e-7, max_iter: int = 40,
                         chunk: int = 16, flat: bool = False, depth: int = 4,
                         Pi0: Optional[np.ndarray] = None, geom: str = "linear") -> Dict[str, object]:
    """
    Nodal Pi (Pi(R) = 1) solving E[B_g(Pi)] Pi = rhs at the row radii (least
    squares), rows at t_x = 3 R. Lengths are in units of l_k when eps rho = 1.

    Vacuum subtraction: the equation is solved as
        E[B_g(Pi)] Pi - E[B_flat] 1 = rhs,
    i.e. the flat operator applied to the vacuum Pi = 1 is subtracted. In the
    exact theory E[B_flat] 1 = 0 in the bulk, so this changes nothing there;
    it removes (i) the quadrature error of the row sums, common to every
    geometry, which would otherwise act as a fake source of the size of a
    weak field, and (ii) the cut-off past of rows next to the wall, which
    would otherwise demand a source to keep the vacuum Pi = 1. The flat
    operator acting on u = Pi - 1 (F1.3) is exactly this difference form.

    Iteration (preconditioned residual, Anderson-accelerated). Each step
    builds the rows in the geometry of the current Pi and evaluates the
    residual F(Pi) = E[B_g(Pi)] Pi - E[B_flat] 1 - rhs; the correction is
    delta = -M^+ F with the FIXED preconditioner M = D * E[B_flat] (the
    linearization of (a) in the weak field: the local limit gives D lap u,
    K.11 local). The fixed point (F = 0) does not depend on M, M only sets
    the speed. A first version froze the geometry and solved the linear
    equation for the field in it (Picard); that is unstable at moderate
    fields, because the frozen -R/2 term of the old geometry acts on the new
    field like a potential of the wrong sign (F1.4 2+1, 18 Oct: Pi driven
    negative). Steps that would push Pi below 0.05 are halved and the
    Anderson history is reset. flat=True: one solve with the flat rows (no
    back-reaction), for reference. Pi0: starting field (e.g. the solution for
    a weaker source, scaled), default the vacuum Pi = 1. geom: interpolation of
    the geometry between the nodes ('linear' or 'spline', see
    curved_radial_rows; 3+1 needs 'spline', K.11c). Convergence: the
    largest correction relative to max |Pi - 1| is below tol (and at least two
    steps were taken).
    """
    nodes = np.asarray(nodes, float)
    rows = np.asarray(rows, float)
    chunks = [rows[i:i + chunk] for i in range(0, rows.size, chunk)]
    t0 = time.time()

    def build(Pi_nodes, it):
        parts = parallel(rows_task, [(D, rho, eps, 3.0 * R, ch, R, nodes, Pi_nodes, tuple(quad), geom) for ch in chunks],
                         workers, lambda s: None,
                         lambda i, r_, d, n_: "")
        A = np.vstack(parts)
        progress(f"{label} iteration {it}: rows built ({time.time()-t0:.0f}s)")
        return A

    Pi = np.ones(len(nodes))
    A_vac = build(Pi, 0)
    e_vac = A_vac @ Pi                           # E[B_flat] 1 at the rows (vacuum subtraction)
    if flat:
        b = rhs + e_vac - A_vac[:, -1]
        sol, *_ = np.linalg.lstsq(A_vac[:, :-1], b, rcond=None)
        Pi_new = np.r_[sol, 1.0]
        res = float(np.sqrt(np.mean((A_vac @ Pi_new - e_vac - rhs) ** 2)))
        return {"Pi": Pi_new, "iterations": 1, "converged": True, "residual": res, "steps": []}
    Minv = np.linalg.pinv(D * A_vac[:, :-1])      # fixed preconditioner (free nodes; Pi(R) = 1 is fixed)
    if Pi0 is not None:
        Pi = np.asarray(Pi0, float).copy()
        Pi[-1] = 1.0
    X, G = [], []                                # Anderson history: iterates and corrections g
    steps, resids = [], []
    converged = False
    A = A_vac
    for it in range(1, max_iter + 1):
        A = build(Pi, it)
        F = A @ Pi - e_vac - rhs
        g = np.r_[-(Minv @ F), 0.0]
        # RELATIVE step: the preconditioner is the predicted weak-field linearization of (a), so its first step
        # already lands on the prediction; an absolute tolerance would stop there for a weak source and bias
        # the result towards (a). The step is measured against the size of the field itself.
        scale = max(float(np.max(np.abs(Pi - 1.0))), 1e-12)
        step = float(np.max(np.abs(g))) / scale
        steps.append(step)
        resids.append(float(np.sqrt(np.mean(F ** 2))))
        if step < tol and it > 1:
            converged = True
            break
        X.append(Pi.copy())
        G.append(g.copy())
        X, G = X[-(depth + 1):], G[-(depth + 1):]
        if len(G) > 1:
            dG = np.array([G[j + 1] - G[j] for j in range(len(G) - 1)]).T
            dX = np.array([X[j + 1] - X[j] for j in range(len(X) - 1)]).T
            gam, *_ = np.linalg.lstsq(dG, g, rcond=None)
            Pi_next = Pi + g - (dX + dG) @ gam
        else:
            Pi_next = Pi + g
        Pi_next[-1] = 1.0
        h = 1.0
        while np.min(Pi_next) <= 0.05 and h > 1e-3:   # safeguard: shorter plain step, history reset
            h *= 0.5
            Pi_next = Pi + h * g
            X, G = [], []
        Pi = Pi_next
    # residual of the returned field in the last geometry built (at convergence the step is below tol)
    res = float(np.sqrt(np.mean((A @ Pi - e_vac - rhs) ** 2)))
    return {"Pi": Pi, "iterations": len(steps), "converged": converged, "residual": res, "steps": steps,
            "residuals": resids}
