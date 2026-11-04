"""
core.operators.bd_expect
========================

The EXPECTED smeared Benincasa-Dowker operator on a static field, computed
by quadrature (no sampling, no noise). This is F1.6 (E45b) variant (a), chosen by
the author on 15 Oct 2026: "compare the discrete solution with the
solution of the nonlocal equation with the same eps; test 1/r only where
r >> l_k".

For a Poisson sprinkling of density rho, the expectation of the smeared
operator at an element x is

    E[B_eps f](x) = -a_d (eps rho)^(2/d) [ f(x) + rho eps
                    int_{J^-(x) in region} f(y) g(eps rho V(x, y)) dV(y) ]

with g(w) = sum_i C_i w^i / i! e^(-w), where C_i are the BD coefficients
(BD_COEFFS). g is the Poisson average of the smeared kernel: the
probability that the interval [y, x] holds n elements is Poisson with mean
rho V, and sum_n K_eps(n) P(n) = eps g(eps rho V). It tends to the
d'Alembertian when f varies slowly on the nonlocality scale
l_k = (eps rho)^(-1/d). Below that scale it is a nonlocal operator with a
computable response, and its expectation is what one sprinkling
estimates, with noise.

Coordinates of the past of x = (t_x, x_s):
- the spatial offset rho_s = |y_s - x_s| and its direction;
- s = (eps rho kappa)^(2/d) tau^2, where tau is the proper time
  (tau^2 = dt^2 - rho_s^2), so that w = eps rho V = s^(d/2).

In these variables the volume element is regular on the light cone:

    dV = rho_s^(n-1) d rho_s dOmega * ds * (eps rho kappa)^(-2/d) / (2 dt),
    dt = sqrt(rho_s^2 + tau^2).

The s integral of g(s^(d/2)) vanishes for the BD coefficients. That is the
cancellation that keeps the operator finite: far contributions cancel ray
by ray. The difference form f(y) - f(x) together with the separate kernel
integral I_one keeps the field-dependent part free of that cancellation.

Region: t_y >= 0 and |y_s| <= R. No elements lie outside, exactly as for
the sprinkled operator. Both boundaries are integrated exactly:
- the wall |y_s| = R gives the upper limit of rho_s for each direction;
- the bottom t_y = 0 gives the upper limit of s for each rho_s.
Quadrature is Gauss-Legendre in z with rho_s = rho_lim z^2, which
concentrates nodes at the apex, and composite Gauss-Legendre panels in s
(geometric near s = 0, where 1/dt is peaked for small rho_s).

Validated by the calibration test K.4 and the core test `bd_expectation`:
- against the independent light-cone integral in 1+1 (0.6395 / 0.9822);
- the local limit (ratio -> 1 as l_k / w -> 0) in 1+1, 2+1 and 3+1.

Functions:
    expected_apply(D, rho, eps, t_x, r_x, R, f_radial, f_x)
                         E[B f](x) for a static radial f
    expected_apply_general(D, rho, eps, t_x, r_x, R, f, f_x)
                         the same for f(|y_s - x_s|, |y_s|, t_y), e.g. centred at x
    expected_radial_rows(D, rho, eps, t_x, r_rows, R, nodes)
                         rows of E[B] on the radial hat basis (as bd_static_rows)
"""

from __future__ import annotations

import math
from typing import Callable, Dict, Sequence

import numpy as np

from .bd import BD_COEFFS, BD_NORM
from ..order.sprinkle import interval_volume

# panel edges in w = eps rho V (the kernel is e^-w suppressed beyond ~60)
_W_EDGES = (0.05, 0.25, 1.0, 2.5, 5.0, 10.0, 20.0, 40.0, 80.0)
# geometric panels in s below the first edge: near the apex the volume
# element carries 1/dt = 1/sqrt(rho_s^2 + tau^2), peaked at s ~ lam^(2/d) rho_s^2,
# which is tiny for small rho_s; panels growing by 10 resolve it for every rho_s
_S_GEOM = tuple(10.0 ** k for k in range(-14, 0))


def _s_edges(D: int) -> np.ndarray:
    """Panel edges in s = w^(2/d): 0, geometric panels near 0, then the w edges."""
    top = np.array(_W_EDGES) ** (2.0 / D)
    geo = [e for e in _S_GEOM if e < top[0]]
    return np.array([0.0] + geo + list(top))


def g_kernel(w: np.ndarray, d: int) -> np.ndarray:
    """Poisson average of the smeared kernel per unit eps: sum_i C_i w^i/i! e^-w."""
    s = np.zeros_like(w)
    for i, c in enumerate(BD_COEFFS[d]):
        s = s + c * w ** i / math.factorial(i)
    return s * np.exp(-w)


def _directions(n: int, n_ang: int):
    """Spatial directions as cos(angle to the x-axis) and their solid-angle weights."""
    if n == 1:
        return np.array([1.0, -1.0]), np.array([1.0, 1.0])
    if n == 2:
        # phi in [0, pi], doubled (the integrand is symmetric in phi -> -phi)
        phi = (np.arange(n_ang) + 0.5) * math.pi / n_ang
        return np.cos(phi), np.full(n_ang, 2.0 * math.pi / n_ang)
    c, wc = np.polynomial.legendre.leggauss(n_ang)          # c = cos(theta)
    return c, 2.0 * math.pi * wc


def _wall_distance(r_x: float, cos: np.ndarray, R: float) -> np.ndarray:
    """Distance from a point at radius r_x along a direction (cos to its radius vector) to |y| = R."""
    return -r_x * cos + np.sqrt(np.maximum(R * R - r_x * r_x * (1.0 - cos * cos), 0.0))


def _nodes(D, rho, eps, t_x, r_x, R, n_rho, n_s, n_ang, with_offset=False):
    """
    All quadrature nodes for the past of x. Returns (W, r_y, t_y), where W
    already holds the measure and the kernel g, so that
    rho eps sum(W f(y)) approximates rho eps int f g dV. Shape:
    (directions, rho_s nodes, s nodes). with_offset=True also returns
    rho_s = |y_s - x_s| at every node (for fields centred at x).
    """
    n = D - 1
    kappa = interval_volume(D, 1.0)
    lam = eps * rho * kappa                                  # w = lam tau^d
    cos, wa = _directions(n, n_ang)
    rho_lim = np.minimum(t_x, _wall_distance(r_x, cos, R))  # (A,)
    z, wz = np.polynomial.legendre.leggauss(n_rho)
    z = 0.5 * (z + 1.0)
    wz = 0.5 * wz
    rs = rho_lim[:, None] * z[None, :] ** 2                  # (A, P)
    wr = wz[None, :] * 2.0 * rho_lim[:, None] * z[None, :]
    # s = lam^(2/d) tau^2 ; the bottom t_y >= 0 means dt <= t_x, s <= lam^(2/d) (t_x^2 - rho_s^2)
    l2 = lam ** (2.0 / D)
    s_max = l2 * np.maximum(t_x * t_x - rs * rs, 0.0)        # (A, P)
    u, wu = np.polynomial.legendre.leggauss(n_s)
    u = 0.5 * (u + 1.0)
    wu = 0.5 * wu
    edges = _s_edges(D)                                      # panel edges in s
    S, WS = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        hi = np.minimum(b, s_max)[..., None]                 # (A, P, 1)
        span = np.maximum(hi - a, 0.0)
        S.append(a + span * u)
        WS.append(span * wu)
    s = np.concatenate(S, axis=-1)                           # (A, P, Q)
    ws = np.concatenate(WS, axis=-1)
    RS = rs[..., None]
    tau2 = s / l2
    dt = np.sqrt(RS * RS + tau2)
    t_y = t_x - dt
    CO = cos[:, None, None]
    r_y = np.sqrt(np.maximum(r_x * r_x + RS * RS + 2.0 * r_x * RS * CO, 0.0))
    radial = RS ** (n - 1) if n > 1 else np.ones_like(RS)
    meas = radial / (l2 * 2.0 * dt)
    W = wa[:, None, None] * wr[..., None] * ws * meas * g_kernel(s ** (D / 2.0), D)
    if with_offset:
        return W, r_y, t_y, np.broadcast_to(RS, W.shape)
    return W, r_y, t_y


def _past_integral(D, rho, eps, t_x, r_x, R, values_of, n_rho=300, n_s=12, n_ang=32):
    """
    Returns (I_diff, I_one), where
        I_diff = rho eps int (f(y) - f(x)) g dV   (the field-dependent part)
        I_one  = rho eps int g dV                 (the kernel alone)
    over the past of x = (t_x, r_x along axis 1) inside the region. Then
    rho eps int f g dV = I_diff + f(x) I_one. values_of(radius, t_y) gives f(y).
    """
    W, r_y, t_y = _nodes(D, rho, eps, t_x, r_x, R, n_rho, n_s, n_ang)
    I_one = rho * eps * float(W.sum())
    fx = float(values_of(np.array([r_x]), np.array([t_x]))[0])
    I_diff = rho * eps * float((W * (values_of(np.broadcast_to(r_y, W.shape), t_y) - fx)).sum())
    return I_diff, I_one


def expected_apply(D: int, rho: float, eps: float, t_x: float, r_x: float, R: float,
                   f_radial: Callable[[np.ndarray], np.ndarray], f_x: float, **kw) -> float:
    """E[B f](x) for a static radial function f(|y_s|); f_x = f at x."""
    pref = -BD_NORM[D] * (eps * rho) ** (2.0 / D)
    I_diff, I_one = _past_integral(D, rho, eps, t_x, r_x, R, lambda r, t: f_radial(r), **kw)
    return pref * (f_x * (1.0 + I_one) + I_diff)


def expected_apply_general(D: int, rho: float, eps: float, t_x: float, r_x: float, R: float,
                           f: Callable[[np.ndarray, np.ndarray, np.ndarray], np.ndarray], f_x: float,
                           n_rho=300, n_s=12, n_ang=32) -> float:
    """
    E[B f](x) for a static field given as f(rho_s, r_y, t_y): rho_s = |y_s - x_s|,
    r_y = |y_s| (radius from the region's axis), t_y. Example: the Gaussian
    centred at x is f = lambda d, r, t: exp(-d^2 / 2 w^2).
    """
    W, r_y, t_y, rs = _nodes(D, rho, eps, t_x, r_x, R, n_rho, n_s, n_ang, with_offset=True)
    pref = -BD_NORM[D] * (eps * rho) ** (2.0 / D)
    I_one = rho * eps * float(W.sum())
    I_diff = rho * eps * float((W * (f(rs, r_y, t_y) - f_x)).sum())
    return pref * (f_x * (1.0 + I_one) + I_diff)


def expected_radial_rows(D: int, rho: float, eps: float, t_x: float, r_rows: Sequence[float],
                         R: float, nodes: np.ndarray, n_rho=300, n_s=12, n_ang=32) -> np.ndarray:
    """
    Rows of the expected operator on the radial hat basis: E[B u](r_i) =
    sum_k A[i, k] u_k, for u piecewise linear on `nodes`. nodes[-1] must be
    >= R, because the region holds no elements beyond R.

    The constant u = 1 gives the row sum pref * (1 + I_one). It is not
    zero near the bottom and the wall, where the past is cut off.
    """
    from .bd import hat_weights
    pref = -BD_NORM[D] * (eps * rho) ** (2.0 / D)
    K = len(nodes)
    A = np.zeros((len(r_rows), K))
    for i, r_x in enumerate(r_rows):
        W, r_y, _ = _nodes(D, rho, eps, t_x, float(r_x), R, n_rho, n_s, n_ang)
        rr = np.minimum(np.broadcast_to(r_y, W.shape).ravel(), nodes[-1])
        k, w0, w1 = hat_weights(rr, nodes)
        Wf = W.ravel()
        row = np.bincount(k, weights=Wf * w0, minlength=K + 1)[:K]
        row = row + np.bincount(k + 1, weights=Wf * w1, minlength=K + 1)[:K]
        kx, a0, a1 = hat_weights(np.array([min(float(r_x), nodes[-1])]), nodes)
        own = np.zeros(K)
        own[kx[0]] += a0[0]
        if kx[0] + 1 < K:
            own[kx[0] + 1] += a1[0]
        A[i] = pref * (own + rho * eps * row)
    return A


# =============================================================================
# The nonlocal field equation in units of l_k (F1.6 (E45b), deterministic shape test)
# =============================================================================

def shape_nodes(R: float, a: float, dr: float, ratio: float) -> np.ndarray:
    """
    Radial nodes in units of l_k: uniform with step dr up to max(4, 2a)
    (the source and the nonlocality scale), then geometric with `ratio`
    up to R.
    """
    inner = np.arange(0.0, max(4.0, 2.0 * a) + 1e-9, dr)
    geo = [float(inner[-1])]
    while geo[-1] * ratio < R:
        geo.append(geo[-1] * ratio)
    return np.unique(np.r_[inner, geo[1:], R])


def expected_green_profile(D: int, R: float, a: float, dr: float, ratio: float, rows_per_element: int,
                           q: float = 1.0, **kw) -> Dict[str, np.ndarray]:
    """
    Solve E[B] u = q 1(r < a) for a static radial u with u(R) = 0, in a
    cylinder of radius R, all lengths in units of l_k = (eps rho)^(-1/d)
    (the expected operator depends on eps and rho only through l_k, so
    one computation per dimension covers every eps and N). Rows are placed
    rows_per_element times per element of the node grid (at the centres
    of equal sub-intervals) at t_x = 3 R, so the bottom of the region never
    enters (the past is cut only by the wall, as in a tall cylinder).
    Least squares over the rows; nodal values u (u at R = 0).
    """
    eps = 0.05
    rho = 1.0 / eps                                        # eps rho = 1 -> l_k = 1
    nodes = shape_nodes(R, a, dr, ratio)
    fr = (np.arange(rows_per_element) + 0.5) / rows_per_element
    rows = np.sort((nodes[:-1, None] + fr[None, :] * np.diff(nodes)[:, None]).ravel())
    A = expected_radial_rows(D, rho, eps, 3.0 * R, rows, R, nodes, **kw)
    rhs = q * (rows < a)
    sol, *_ = np.linalg.lstsq(A[:, :-1], rhs, rcond=None)
    u = np.r_[sol, 0.0]
    return {"nodes": nodes, "u": u, "rows": rows, "residual": A @ u - rhs,
            "cond": float(np.linalg.cond(A[:, :-1]))}
