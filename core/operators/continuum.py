"""
core.operators.continuum
========================

Continuum solutions of the field equation A5' in kernel form,

    lap Pi = kappa rho_s,     Pi -> 1 (Pi = 1 on the boundary |x| = R),

for a static source on the axis, in n = D - 1 spatial dimensions. Used as
the reference of F1.6 (E45b) (discrete operators) and as the input profile of
K.3 (E45c) (clock-rate readback).

    ball_source_solution   uniform source density q inside |x| < a, the
                           exact radial solution with Pi(R) = 1
    softened_point         Pi = 1 - eps * (G(sqrt(r^2 + s^2)) - G(sqrt(R^2 + s^2)))
                           with the Green-function shape G = 1/r (n = 3),
                           -log r (n = 2), -r (n = 1): a smooth stand-in
                           for a point source
"""

from __future__ import annotations

import math
from typing import Callable

import numpy as np


def green_shape(n: int) -> Callable[[np.ndarray], np.ndarray]:
    """Radial shape of the Green function of the Laplacian in n dimensions (up to a constant)."""
    if n == 3:
        return lambda r: 1.0 / r
    if n == 2:
        return lambda r: -np.log(r)
    if n == 1:
        return lambda r: -r
    return lambda r: r ** (2.0 - n)


def softened_point(n: int, eps: float, soft: float, R: float) -> Callable[[np.ndarray], np.ndarray]:
    """Pi(r) = 1 - eps (G(sqrt(r^2+s^2)) - G(sqrt(R^2+s^2))); Pi(R) = 1, Pi < 1 inside."""
    G = green_shape(n)
    gR = float(G(math.sqrt(R * R + soft * soft)))

    def Pi(r):
        r = np.asarray(r, float)
        return 1.0 - eps * (G(np.sqrt(r * r + soft * soft)) - gR)
    return Pi


def ball_source_solution(n: int, q: float, a: float, R: float) -> Callable[[np.ndarray], np.ndarray]:
    """
    u = Pi - 1 for lap u = q inside |x| < a (0 outside), u(R) = 0, regular at 0.
    Total source Q = q * ball(n, a). Outside: u = -(Q / (omega_n)) * (G(r) - G(R)) * c_n
    with the Green function of the Laplacian; inside: matched quadratic.
    Written out per n (n = 1, 2, 3).
    """
    if n == 3:
        Q = q * 4.0 / 3.0 * math.pi * a ** 3
        def u(r):
            r = np.asarray(r, float)
            out = -Q / (4 * math.pi) * (1.0 / np.maximum(r, 1e-300) - 1.0 / R)
            ins = r < a
            # inside: u = q r^2 / 6 + C, continuous at a
            ua = -Q / (4 * math.pi) * (1.0 / a - 1.0 / R)
            out = np.where(ins, ua + q * (r * r - a * a) / 6.0, out)
            return out
        return u
    if n == 2:
        Q = q * math.pi * a * a
        def u(r):
            r = np.asarray(r, float)
            out = Q / (2 * math.pi) * np.log(np.maximum(r, 1e-300) / R)
            ua = Q / (2 * math.pi) * math.log(a / R)
            return np.where(r < a, ua + q * (r * r - a * a) / 4.0, out)
        return u
    if n == 1:
        Q = q * 2 * a
        def u(r):
            r = np.asarray(r, float)
            out = Q / 2.0 * (r - R)
            ua = Q / 2.0 * (a - R)
            return np.where(r < a, ua + q * (r * r - a * a) / 2.0, out)
        return u
    raise ValueError("n must be 1, 2 or 3")
