"""
core.chains.geodesics
=====================

Continuum predictions for the motion of a body between two fixed events.

The action of Causal Physics 2.0 for a body is S = -nu_B W with
W = sum Pi dtau; in the continuum W = int Pi sqrt(1 - |v|^2) dt (Pi the
clock rate, c = 1). Its maximum between (0, x0) and (T, x0) is the
geodesic of the conformal metric Pi^2 eta.

For a radial field Pi(|x - m|) and both end events at the same place x0,
the maximising curve stays on the straight line through m and x0 (by
symmetry), so the problem is one-dimensional along that line and
`geodesic_1p1` solves it exactly (port of the 2.16 routine, unchanged
mathematics): the Lagrangian has no explicit t, so
E = omega / sqrt(1 - xdot^2) is conserved and

    T/2 = int_0^{y_top} dy / sqrt(1 - f(y)^2 / f(y_top)^2)

with f(y) = omega(x0 + s y) increasing along the path direction s.
"""

from __future__ import annotations

import math
from typing import Any, Callable, Dict, Sequence

import numpy as np


def geodesic_1p1(omega: Callable[[float], float], T: float, x0: float,
                 fracs: Sequence[float] = (0.25, 0.5), nodes: int = 96) -> Dict[str, Any]:
    """
    Maximal int omega(x) sqrt(1 - xdot^2) dt between (0, x0) and (T, x0).
    Returns the displacement at the time fractions `fracs` (signed, along
    +x), the turning point x_top and the Newtonian estimate
    omega'(x0)/omega(x0) * T^2 / 8. omega must increase monotonically
    along the path direction on [0, T/2] (checked).
    """
    zs, ws = np.polynomial.legendre.leggauss(nodes)
    zs, ws = 0.5 * (zs + 1.0), 0.5 * ws
    h = 1e-6
    slope = (omega(x0 + h) - omega(x0 - h)) / (2 * h)
    newton = slope / omega(x0) * T * T / 8.0
    if abs(slope) < 1e-12:
        return {"x_top": 0.0, "disp": {f: 0.0 for f in fracs}, "newton": 0.0, "E": omega(x0)}
    s = 1.0 if slope > 0 else -1.0
    f = lambda y: float(omega(x0 + s * y))
    ys = np.linspace(0, T / 2, 401)
    fy = [f(y) for y in ys]
    if any(b <= a for a, b in zip(fy, fy[1:])):
        raise ValueError("omega is not increasing along the path direction on [0, T/2]")

    def g(y_top: float, r: float) -> float:
        # int_{y_top - r}^{y_top} dy / sqrt(1 - f^2/E^2), E = f(y_top); y = y_top - r z^2
        E = f(y_top)
        tot = 0.0
        for z, wz in zip(zs, ws):
            y = y_top - r * z * z
            den = 1.0 - (f(y) / E) ** 2
            tot += wz * 2.0 * r * z / math.sqrt(max(den, 1e-300))
        return tot

    lo, hi = 0.0, T / 2
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if g(mid, mid) < T / 2:
            lo = mid
        else:
            hi = mid
    y_top = 0.5 * (lo + hi)
    disp = {}
    for fr in fracs:
        t_target = fr * T if fr <= 0.5 else (1 - fr) * T
        a, b = 0.0, y_top
        for _ in range(60):
            m = 0.5 * (a + b)
            if T / 2 - g(y_top, y_top - m) < t_target:
                a = m
            else:
                b = m
        disp[fr] = s * 0.5 * (a + b)
    return {"x_top": s * y_top, "disp": disp, "newton": newton, "E": f(y_top),
            "half_time_check": g(y_top, y_top)}


def radial_pi(eps: float, soft: float) -> Callable[[np.ndarray], np.ndarray]:
    """Pi(r) = 1 - eps / sqrt(r^2 + soft^2): a softened point source (kernel units)."""
    def Pi(r):
        return 1.0 - eps / np.sqrt(np.asarray(r, float) ** 2 + soft * soft)
    return Pi


def radial_prediction(eps: float, soft: float, b: float, T: float) -> Dict[str, Any]:
    """
    Mid-time displacement AWAY from the source for end events at distance
    b from it, separated by T in time, in the field radial_pi(eps, soft).
    """
    Pi = radial_pi(eps, soft)
    geo = geodesic_1p1(lambda y: float(Pi(b + y)), T, 0.0, fracs=(0.25, 0.5))
    return {"mid": geo["disp"][0.5], "quarter": geo["disp"][0.25], "newton": geo["newton"],
            "x_top": geo["x_top"]}


def radial_prediction_k(eps: float, soft: float, b: float, T: float, k: float) -> Dict[str, Any]:
    """
    Mid-time displacement of the continuum geodesic of the metric Pi^(2k) eta
    (the continuum limit of the heaviest chains of sum Pi^k on a flat order;
    k = 1 is radial_prediction). Used by F2.4 (E46b) for the competing functional H2.
    """
    Pi = radial_pi(eps, soft)
    geo = geodesic_1p1(lambda y: float(Pi(b + y)) ** k, T, 0.0, fracs=(0.5,))
    return {"mid": geo["disp"][0.5], "k": k}
