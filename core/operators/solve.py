"""
core.operators.solve
====================

Solving the static field equation with a discrete Laplacian, and reading
the radial profile.

    solve_radial_ls(A, rhs, ...)     least squares for the radial-basis
                                     (BD-static) rows, outer node fixed to 0
    fit_profile(r, u, window)        u = A (r^g - 1)/g + B (Box-Cox family;
                                     g = 0 is the logarithm): best g by a scan,
                                     A and B by linear least squares
"""

from __future__ import annotations

from typing import Any, Dict, Sequence

import numpy as np


def solve_radial_ls(A: np.ndarray, rhs: np.ndarray) -> np.ndarray:
    """Nodal values u_0..u_{K-1} (u_K = 0 at the boundary) from rows A u = rhs."""
    sol, *_ = np.linalg.lstsq(A[:, :-1], rhs, rcond=None)
    return np.r_[sol, 0.0]


def boxcox(r: np.ndarray, g: float) -> np.ndarray:
    return np.log(r) if abs(g) < 1e-9 else (r ** g - 1.0) / g


def fit_profile(r: Sequence[float], u: Sequence[float], window: Sequence[float],
                g_grid: Sequence[float] = tuple(np.round(np.arange(-3.0, 2.0001, 0.01), 4))) -> Dict[str, Any]:
    """Best Box-Cox exponent g of u(r) inside the window (see module notes)."""
    r, u = np.asarray(r, float), np.asarray(u, float)
    m = (r >= window[0]) & (r <= window[1]) & np.isfinite(u)
    r, u = r[m], u[m]
    if r.size < 4:
        return {"g": float("nan"), "A": float("nan"), "B": float("nan"), "rss": float("nan"), "n": int(r.size)}
    best = None
    for g in g_grid:
        M = np.column_stack([boxcox(r, g), np.ones_like(r)])
        coef, *_ = np.linalg.lstsq(M, u, rcond=None)
        rss = float(((M @ coef - u) ** 2).sum())
        if best is None or rss < best[0]:
            best = (rss, g, coef)
    return {"g": float(best[1]), "A": float(best[2][0]), "B": float(best[2][1]), "rss": best[0], "n": int(r.size)}


def binned(r: np.ndarray, u: np.ndarray, edges: Sequence[float]):
    """Mean r and mean u per radial bin (bins with no node are nan)."""
    rc, uc = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        m = (r >= a) & (r < b)
        rc.append(float(r[m].mean()) if m.any() else float("nan"))
        uc.append(float(u[m].mean()) if m.any() else float("nan"))
    return np.array(rc), np.array(uc)
