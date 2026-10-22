"""
core.order.sprinkle
===================

Sprinklings WITH coordinates: the fixed backgrounds of phase A.

A sprinkling here is N points drawn uniformly (or with a prescribed
density) in a region of D-dimensional Minkowski coordinates, sorted by the
time coordinate. The row order is therefore a linear extension of every
causal relation used in the lab (flat, conformal, H1), because each of
them only relates an earlier time to a later one.

Regions
-------
    Diamond(D)              the 2.x causal diamond between (0, 0..) and
                            (1, 0..); `diamond_2x` reproduces the 2.x
                            point set of `sprinkled_cone` exactly
    Cylinder(D, t0, t1, R)  t in [t0, t1], |x| <= R (a ball of radius R in
                            the D - 1 spatial directions): the region for
                            static fields around a source on the axis
    Box(D, lo, hi)          an axis-aligned box (screens and source slabs
                            of the cone instrument, K.2 (E47))
    Interval(D, p, q)       the Alexandrov interval between two timelike
                            separated points (F2.2 (E46): chains between fixed
                            end events)

Densities
---------
`sprinkle(region, n, seed, density=None)` draws exactly n points (a
binomial process, as the 2.x lab did). With a density function f(X) <= 1
(relative to its maximum) the points are thinned by acceptance f, which
gives a sprinkling of density proportional to f: the conformal
sprinklings of a metric Pi^2 eta use f = (Pi / Pi_max)^D.

Seeds
-----
`derive_seed(base, *keys)` maps a base seed and any labels (experiment,
size, replica, arm) to an independent 63-bit seed by SHA-256, so seeds do
not depend on the order in which runs are made, and a replication with
another base seed K changes every stream at once.
"""

from __future__ import annotations

import hashlib
import math
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, Optional, Sequence

import numpy as np


# =============================================================================
# Seeds
# =============================================================================

def derive_seed(base: int, *keys: Any) -> int:
    """
    Independent seed for one stream: SHA-256 of the base seed and the keys
    (their repr), first 8 bytes, top bit cleared. Stable across platforms
    and Python versions (repr of ints, floats and strings is stable).
    """
    h = hashlib.sha256(repr((int(base),) + tuple(keys)).encode("utf-8")).digest()
    return int.from_bytes(h[:8], "little") & ((1 << 63) - 1)


# =============================================================================
# Regions
# =============================================================================

def ball_volume(n: int, R: float = 1.0) -> float:
    """Volume of a ball of radius R in n dimensions (n = 0 gives 1)."""
    return math.pi ** (n / 2.0) / math.gamma(n / 2.0 + 1.0) * R ** n


def interval_volume(D: int, tau: float) -> float:
    """
    Volume of the Alexandrov interval of proper time tau in D-dimensional
    Minkowski space: two cones of height tau/2 on a ball of radius tau/2,
    V = 2 * ball(D-1, tau/2) * (tau/2) / D.
    (D = 2: tau^2/2; D = 3: pi tau^3/12; D = 4: pi tau^4/24.)
    """
    return 2.0 * ball_volume(D - 1, tau / 2.0) * (tau / 2.0) / D


class Region:
    """A region of D-dimensional coordinates (t, x_1, ..., x_{D-1})."""
    D: int

    def volume(self) -> float:
        raise NotImplementedError

    def propose(self, rng: np.random.Generator, m: int) -> np.ndarray:
        """m candidate points, uniform in the region (rejection inside)."""
        raise NotImplementedError

    def describe(self) -> Dict[str, Any]:
        raise NotImplementedError


@dataclass
class Cylinder(Region):
    """t in [t0, t1], spatial ball |x| <= R."""
    D: int
    t0: float = 0.0
    t1: float = 1.0
    R: float = 0.5

    def volume(self) -> float:
        return (self.t1 - self.t0) * ball_volume(self.D - 1, self.R)

    def propose(self, rng, m):
        t = rng.uniform(self.t0, self.t1, size=m)
        n = self.D - 1
        if n == 0:
            return t[:, None]
        # uniform in the ball by rejection from the cube (efficiency >= 52 % for n <= 3)
        out = []
        need = m
        while need > 0:
            x = rng.uniform(-self.R, self.R, size=(2 * need + 16, n))
            x = x[(x ** 2).sum(axis=1) <= self.R ** 2]
            out.append(x[:need])
            need -= len(out[-1])
        x = np.concatenate(out)[:m]
        return np.column_stack([t, x])

    def describe(self):
        return {"region": "cylinder", "D": self.D, "t0": self.t0, "t1": self.t1, "R": self.R}


@dataclass
class Box(Region):
    """Axis-aligned box: lo[k] <= X[:, k] <= hi[k] for k = 0..D-1 (k = 0 is t)."""
    D: int
    lo: Sequence[float] = ()
    hi: Sequence[float] = ()

    def volume(self) -> float:
        return float(np.prod(np.asarray(self.hi, float) - np.asarray(self.lo, float)))

    def propose(self, rng, m):
        lo, hi = np.asarray(self.lo, float), np.asarray(self.hi, float)
        return lo + (hi - lo) * rng.random((m, self.D))

    def describe(self):
        return {"region": "box", "D": self.D, "lo": list(self.lo), "hi": list(self.hi)}


@dataclass
class Interval(Region):
    """
    The Alexandrov interval between p = (0, x0) and q = (T, x0): points with
    |x - x0| <= min(t, T - t). (Both end points share the spatial position.)
    """
    D: int
    T: float = 1.0
    x0: Sequence[float] = ()

    def _x0(self):
        return np.zeros(self.D - 1) if len(self.x0) == 0 else np.asarray(self.x0, float)

    def volume(self) -> float:
        return interval_volume(self.D, self.T)

    def propose(self, rng, m):
        x0 = self._x0()
        out = []
        need = m
        while need > 0:
            k = 3 * need + 32
            t = rng.uniform(0.0, self.T, size=k)
            x = rng.uniform(-self.T / 2, self.T / 2, size=(k, self.D - 1))
            ok = np.sqrt((x ** 2).sum(axis=1)) <= np.minimum(t, self.T - t)
            pts = np.column_stack([t[ok], x[ok] + x0])
            out.append(pts[:need])
            need -= len(out[-1])
        return np.concatenate(out)[:m]

    def describe(self):
        return {"region": "interval", "D": self.D, "T": self.T, "x0": list(self._x0())}


@dataclass
class Diamond(Region):
    """The 2.x unit diamond between (0, 0..) and (1, 0..)."""
    D: int

    def volume(self) -> float:
        return interval_volume(self.D, 1.0)

    def propose(self, rng, m):
        return Interval(self.D, 1.0).propose(rng, m)

    def describe(self):
        return {"region": "diamond", "D": self.D}


# =============================================================================
# Sprinklings
# =============================================================================

@dataclass
class Sprinkling:
    """
    Points sorted by time, the region, and the NOMINAL density rho: points
    per unit of PROPER volume (for a conformal sprinkling of Pi^2 eta the
    proper volume element is Pi^D dV, and rho = N / integral Pi^D dV).
    ell = rho^(-1/D) is the discreteness length.
    """
    X: np.ndarray
    D: int
    region: Region
    rho: float
    label: str = ""
    meta: Dict[str, Any] = field(default_factory=dict)

    @property
    def n(self) -> int:
        return len(self.X)

    @property
    def ell(self) -> float:
        return self.rho ** (-1.0 / self.D)

    def r(self, center: Optional[Sequence[float]] = None) -> np.ndarray:
        """Spatial distance of every point from `center` (default: the axis x = 0)."""
        S = self.X[:, 1:]
        if center is not None:
            S = S - np.asarray(center, float)
        return np.sqrt((S ** 2).sum(axis=1))


def sprinkle(region: Region, n: int, seed: int,
             density: Optional[Callable[[np.ndarray], np.ndarray]] = None,
             label: str = "") -> Sprinkling:
    """
    Exactly n points in `region`, uniform or with relative density
    `density(X)` in [0, 1] (thinning). Sorted by time. For a thinned
    sprinkling the nominal density is n / integral(density) dV, where the
    integral is estimated from the acceptance rate of all proposals
    (relative accuracy ~ 1/sqrt(proposals)); `meta["acceptance"]` keeps it.
    """
    rng = np.random.default_rng(int(seed))
    D = region.D
    if density is None:
        X = region.propose(rng, int(n))
        rho = n / region.volume()
        acc = 1.0
    else:
        kept = []
        got, proposed, accepted = 0, 0, 0
        while got < n:
            m = max(1024, 2 * (n - got))
            P = region.propose(rng, m)
            u = rng.random(m)
            ok = u < np.clip(density(P), 0.0, 1.0)
            proposed += m
            accepted += int(ok.sum())
            kept.append(P[ok])
            got += int(ok.sum())
        X = np.concatenate(kept)[:n]
        acc = accepted / proposed
        rho = n / (region.volume() * acc)
    X = X[np.argsort(X[:, 0], kind="stable")]
    return Sprinkling(X=X, D=D, region=region, rho=float(rho), label=label,
                      meta={"seed": int(seed), "acceptance": acc, **region.describe()})


def diamond_2x(D: int, n: int, seed: int) -> Sprinkling:
    """
    The exact point set of the 2.x `sprinkled_cone(D, n, seed)` (same
    RNG calls in the same order), INCLUDING the genesis at row 0, so that
    3.0 instruments can be checked element by element against the archive
    (checked in GTS Lab 3.x against the 2.x archive).
    """
    rng = np.random.default_rng(int(seed))
    pts = []
    while len(pts) < n:
        m = max(64, 4 * (n - len(pts)))
        t = rng.random(m)
        x = rng.uniform(-0.5, 0.5, size=(m, max(1, D - 1))) if D > 1 else None
        for i in range(m):
            rad = min(t[i], 1.0 - t[i])
            if D == 1 or float(np.sqrt((x[i] ** 2).sum())) <= rad:
                pts.append((t[i],) + (tuple(x[i]) if D > 1 else ()))
                if len(pts) >= n:
                    break
    pts.sort(key=lambda z: z[0])
    X = np.array([(0.0,) + (0.0,) * (D - 1)] + pts)
    return Sprinkling(X=X, D=D, region=Diamond(D), rho=n / interval_volume(D, 1.0),
                      label=f"diamond_2x {D}", meta={"seed": int(seed), "genesis": True})


def sprinkle_poisson(region: Region, rho_max: float, rng: np.random.Generator,
                     density: Optional[Callable[[np.ndarray], np.ndarray]] = None) -> np.ndarray:
    """
    A true Poisson process in `region`: Poisson(rho_max * volume) proposals,
    thinned by density(X) in [0, 1] if given. Returns time-sorted points.
    Restricting a sprinkling to a subregion gives exactly this process, so
    small regions (e.g. the diamonds of the clock-rate reading) can be
    sprinkled on their own.
    """
    m = int(rng.poisson(rho_max * region.volume()))
    X = region.propose(rng, m) if m else np.zeros((0, region.D))
    if density is not None and m:
        X = X[rng.random(m) < np.clip(density(X), 0.0, 1.0)]
    return X[np.argsort(X[:, 0], kind="stable")]
