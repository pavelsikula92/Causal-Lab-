"""
core.chains.readings
====================

Readings of the maximal-chain functionals between two fixed events, as
used by the motion tests (F2.x; in GTS Lab 3.x E46 and E46b) and by the
mirror calibration K.7.

Setup: base points P uniform in the causal interval I(p, q), p = (0, 0),
q = (T, 0); marks u uniform in [0, 1). A field Pi_m(x) = 1 - eps /
sqrt(|x - m|^2 + soft^2) of a source at m enters in one of two
representations:

  conformal  points kept where u < Pi^D (density proportional to Pi^D, flat relation)
             L  uniform average over all longest chains
             S  uniform average over all maximal chains with the fewest links
             U  uniform average over ALL maximal chains (link paths)
  weighted   all points, weights Pi (or Pi^2)
             H1 uniform average over chains maximising sum Pi
             H2 the same for sum Pi^2
             G  greedy walk along links to the largest Pi (diagnostic)

A reading d is the displacement at mid-time t = T/2 along the direction e
("away from the source" when the source sits at -b e). Every functional is
a family of chains with weights summing to 1 (weight check).
"""

from __future__ import annotations

from typing import Any, Dict, Sequence

import numpy as np

from ..order.sprinkle import Interval
from ..order.links import links_csr
from .chains import heaviest_between, shortest_between, allpaths_between, greedy_between
from .geodesics import radial_pi

CONFORMAL = ("L", "S", "U")
WEIGHTED = ("H1", "H2", "G")


def with_ends(Y: np.ndarray, D: int, T: float) -> np.ndarray:
    """Add the end events p = (0, 0..) and q = (T, 0..) and sort by time (p first, q last)."""
    X = np.vstack([np.zeros(D), Y, np.r_[T, np.zeros(D - 1)]])
    return X[np.argsort(X[:, 0], kind="stable")]


def read_conformal(P: np.ndarray, u: np.ndarray, D: int, N0: int, T: float, eps: float, soft: float,
                   b: float, sides: Sequence[int] = (1, -1)) -> Dict[str, Any]:
    """L, S, U on the conformal sprinkling, source at -sgn b e for each sgn in sides."""
    e = np.zeros(D - 1)
    e[0] = 1.0
    Pi = radial_pi(eps, soft)
    rho_base = N0 / Interval(D, T).volume()
    out = {k: {} for k in CONFORMAL}
    wchk = []
    for sgn in sides:
        m = -sgn * b * e
        pi = Pi(np.sqrt(((P[:, 1:] - m) ** 2).sum(axis=1)))
        Xc = with_ends(P[u < pi ** D], D, T)
        q = len(Xc) - 1
        hl = heaviest_between(Xc, np.ones(len(Xc)), 0, q, sgn * e)
        ptr, idx = links_csr(Xc, rho_min=rho_base * float(min(pi.min(), 1.0)) ** D)
        hs = shortest_between(Xc, ptr, idx, 0, q, sgn * e)
        hu = allpaths_between(Xc, ptr, idx, 0, q, sgn * e)
        out["L"][sgn], out["S"][sgn], out["U"][sgn] = hl["disp"], hs["disp"], hu["disp"]
        wchk += [hl.get("weight_check", 1.0), hs.get("weight_check", 1.0), hu.get("weight_check", 1.0)]
    return {"d": out, "weight_worst": max(abs(c - 1.0) for c in wchk)}


def read_weighted(P: np.ndarray, D: int, N0: int, T: float, eps: float, soft: float, b: float,
                  sides: Sequence[int] = (1, -1)) -> Dict[str, Any]:
    """H1, H2 and G on the flat sprinkling with weights Pi and Pi^2."""
    e = np.zeros(D - 1)
    e[0] = 1.0
    Pi = radial_pi(eps, soft)
    Xw = with_ends(P, D, T)
    q = len(Xw) - 1
    ptr, idx = links_csr(Xw, rho_min=N0 / Interval(D, T).volume())
    out = {k: {} for k in WEIGHTED}
    wchk = []
    for sgn in sides:
        m = -sgn * b * e
        w = Pi(np.sqrt(((Xw[:, 1:] - m) ** 2).sum(axis=1)))
        h1 = heaviest_between(Xw, w, 0, q, sgn * e)
        h2 = heaviest_between(Xw, w * w, 0, q, sgn * e)
        g = greedy_between(Xw, ptr, idx, w, 0, q, sgn * e)
        out["H1"][sgn], out["H2"][sgn], out["G"][sgn] = h1["disp"], h2["disp"], g["disp"]
        wchk += [h1.get("weight_check", 1.0), h2.get("weight_check", 1.0)]
    return {"d": out, "weight_worst": max(abs(c - 1.0) for c in wchk)}


def mirror_residuals(D: int, N0: int, seed: int, T: float) -> Dict[str, Any]:
    """
    eps = 0, one-sided readings of every functional on P and on its mirror
    MP (M: x1 -> -x1). The identity d(P) + d(MP) = 0 holds exactly when the
    reading is equivariant: links use squared differences only, the
    dynamic programmes depend only on the relation and the index order.
    Returns |d(P) + d(MP)| per functional, d(P) itself and the weight check.
    """
    rng = np.random.default_rng(seed)
    P = Interval(D, T).propose(rng, N0)
    u = rng.random(N0)
    PM = P.copy()
    PM[:, 1] = -PM[:, 1]
    res, val, ww = {}, {}, 0.0
    for reader in ("conformal", "weighted"):
        if reader == "conformal":
            r1 = read_conformal(P, u, D, N0, T, 0.0, 0.2, 0.5, sides=(1,))
            r2 = read_conformal(PM, u, D, N0, T, 0.0, 0.2, 0.5, sides=(1,))
        else:
            r1 = read_weighted(P, D, N0, T, 0.0, 0.2, 0.5, sides=(1,))
            r2 = read_weighted(PM, D, N0, T, 0.0, 0.2, 0.5, sides=(1,))
        ww = max(ww, r1["weight_worst"], r2["weight_worst"])
        for k in r1["d"]:
            a, c = r1["d"][k][1], r2["d"][k][1]
            res[k] = abs(a + c) if (a == a and c == c) else float("nan")
            val[k] = a
    return {"residual": res, "d": val, "weight_worst": ww}
