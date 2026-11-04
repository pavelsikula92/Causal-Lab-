"""
core.order.horizon
==================

Horizon molecules of a causal set: the entropy instrument of the Jacobson
route (author's decision of 22 Oct 2026: "try to derive gravity as
thermodynamics of the deeper layer"; docs/derivations/jacobson_route.md).

Definition (Barton, Counsell, Dowker, Gould, Jubb, Taylor, "Horizon
molecules in causal set theory", arXiv:1909.08620, the 1-molecule)
-------------------------------------------------------------------------
A Cauchy surface Sigma and a causal horizon H split the spacetime into
four regions; C^-_- are the elements in the past of Sigma and OUTSIDE the
horizon (they can still signal to the observer whose horizon H is), C^-_+
those in the past of Sigma and INSIDE. A horizon molecule is a pair
p- < p+ with
    p- in C^-_-,  p+ in C^-_+,
    p+ is the ONLY element of C^- (the past of Sigma) in the future of p-.
(It is a link.) The paper shows that the expected number tends to
a^(d) rho^((d-2)/d) Area(J), J = H cap Sigma, a^(d) a constant of the
dimension only: entropy proportional to area, counted in discreteness units.

This module: the planar case used for calibration,
    Sigma = {t = 0},   H = the null plane t = x  (outside: x - t > 0),
in a flat or conformally flat geometry Pi^2 eta (the order is flat, the
density is rho Pi^d per coordinate volume). In flat space the expected
number per unit area is EXACTLY a^(d) rho^((d-2)/d) for every rho (the
integral scales), with
    <H> / A = rho^2 int_{p in C^-_-} V_+(p) e^{-rho V(p)} dV_p ,
V(p) the volume of the future cone of p below Sigma, V_+(p) its part inside
H (Poisson: exactly one element there, and it is inside). With rho = 1,
h = -t_p and xi = x_p - t_p (dt dx = dh dxi):
    V(h) = b_{d-1} h^d / d,
    V_+(h, xi) = int_{xi/2}^{h} B_{d-1}(s, s - xi) ds  (0 if h < xi/2),
B_n(R, c) the volume of the n-ball of radius R with x_1 <= c. Values
(module function, quadrature): a^(2) = 1/3, a^(3) = 0.218853,
a^(4) = 0.173205 (= sqrt(3)/10 to 9 digits).

Functions
---------
    molecule_constant(d)            a^(d) by quadrature
    count_molecules(X, h_max, ...)  the molecules of a sprinkled point set (planar Sigma and H)
"""

from __future__ import annotations

import math
from typing import Dict

import numpy as np


def _ball(n: int, R: float = 1.0) -> float:
    return math.pi ** (n / 2.0) / math.gamma(n / 2.0 + 1.0) * R ** n


def _ball_cut(n: int, R: float, c: float) -> float:
    """Volume of the n-ball of radius R with x_1 <= c (n = 1, 2, 3)."""
    c = min(max(c, -R), R)
    if n == 1:
        return c + R
    if n == 2:
        return R * R * (math.pi - math.acos(c / R)) + c * math.sqrt(max(R * R - c * c, 0.0))
    if n == 3:
        return math.pi * (R + c) ** 2 * (2.0 * R - c) / 3.0
    raise ValueError("d must be 2, 3 or 4")


def cone_volume(d: int, h: float) -> float:
    """Volume of the future cone of height h in d spacetime dimensions (rho = 1)."""
    return _ball(d - 1) * h ** d / d


def inside_volume(d: int, h: float, xi: float) -> float:
    """V_+(h, xi): the part of the cone of height h, apex at null distance xi outside H, that lies inside H."""
    from scipy import integrate
    if h <= 0.5 * xi:
        return 0.0
    return integrate.quad(lambda s: _ball_cut(d - 1, s, s - xi), 0.5 * xi, h, epsabs=1e-13, epsrel=1e-12)[0]


def molecule_constant(d: int) -> float:
    """a^(d) = int_0^inf dh e^{-V(h)} int_0^{2h} V_+(h, xi) dxi (module notes), by adaptive quadrature."""
    from scipy import integrate
    h_top = (40.0 * d / _ball(d - 1)) ** (1.0 / d)              # e^{-V} < e^{-40} beyond
    f = lambda h: math.exp(-cone_volume(d, h)) * integrate.quad(lambda xi: inside_volume(d, h, xi), 0.0, 2.0 * h,
                                                                 epsabs=1e-13, epsrel=1e-12)[0]
    return integrate.quad(f, 0.0, h_top, limit=400, epsabs=1e-13, epsrel=1e-11)[0]


def h_cut(d: int, rho_min: float, tail: float = 30.0) -> float:
    """Height below Sigma beyond which a molecule's lower element has probability < e^-tail (density >= rho_min)."""
    return (tail * d / (_ball(d - 1) * rho_min)) ** (1.0 / d)


def count_molecules(X: np.ndarray, h_max: float, sel_lo, sel_hi) -> Dict[str, np.ndarray]:
    """
    Horizon molecules of the point set X (rows (t, x, y, ...), all with t < 0, i.e. the past of Sigma = {t = 0})
    for the horizon H = {t = x}. The lower element p- must have h = -t <= h_max and transverse coordinates
    (columns 2..) in [sel_lo, sel_hi) (the counting window; the sprinkled region must extend at least h_max
    beyond it in every direction, so that the whole future cone of p- below Sigma is sprinkled). Returns the
    indices of p- and p+ and the transverse coordinates of p-.
    The future cone of p below Sigma has height -t_p <= h_max, so every point of it is within Euclidean
    distance sqrt(2) h_max of p: a k-d tree query of that radius, then the exact cone test.
    """
    from scipy.spatial import cKDTree
    X = np.asarray(X, float)
    D = X.shape[1]
    if np.any(X[:, 0] >= 0.0):
        raise ValueError("all points must lie in the past of Sigma (t < 0)")
    t, x = X[:, 0], X[:, 1]
    outside = (x - t) > 0.0
    cand = outside & (-t <= h_max)
    if D > 2:
        tr = X[:, 2:]
        cand &= np.all((tr >= np.asarray(sel_lo, float)) & (tr < np.asarray(sel_hi, float)), axis=1)
    ci = np.nonzero(cand)[0]
    tree = cKDTree(X)
    radii = math.sqrt(2.0) * (-t[ci]) * (1.0 + 1e-12)
    lists = tree.query_ball_point(X[ci], radii)
    lo_idx, hi_idx = [], []
    for p, nb in zip(ci, lists):
        nb = np.asarray(nb, dtype=np.int64)
        nb = nb[nb != p]
        if nb.size == 0:
            continue
        dt = t[nb] - t[p]
        dx2 = ((X[nb, 1:] - X[p, 1:]) ** 2).sum(axis=1)
        fut = nb[(dt > 0.0) & (dt * dt >= dx2)]                  # the future of p (all points are below Sigma)
        if fut.size == 1 and (x[fut[0]] - t[fut[0]]) <= 0.0:     # exactly one, and it is inside H
            lo_idx.append(p)
            hi_idx.append(int(fut[0]))
    lo_idx = np.asarray(lo_idx, dtype=np.int64)
    return {"lower": lo_idx, "upper": np.asarray(hi_idx, dtype=np.int64),
            "transverse": X[lo_idx, 2:] if D > 2 else np.zeros((lo_idx.size, 0)), "candidates": int(ci.size)}
