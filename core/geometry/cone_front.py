"""
core.geometry.cone_front
========================

The cone-opening instrument (K.2 (E47) calibration, F7.1 (E50) use).

Input: screen events Q (coordinates of the screen's clocks and rulers --
the observer's frame, not the bulk) and a boolean per event: is it in the
causal future of the source set? Nothing else of the order is used.

Front time per screen cell: within each cell (a bin of the screen
coordinate y), the front lies between the latest event NOT in the future
that precedes the earliest event IN the future, and that earliest event;
the estimate is the midpoint. A delay profile is the difference of two
front profiles read on the SAME screen events (e.g. a curved relation
against the flat one), which cancels the discreteness of the estimate.

Deflection: the front's tilt, alpha(y) = -d(delay)/dy, by a local linear
fit over `half` cells on each side (a plane front crossing a lens at
impact parameter y is turned by the gradient of its delay).
"""

from __future__ import annotations

from typing import Dict, Sequence

import numpy as np


def front_profile(Q: np.ndarray, future: np.ndarray, y: np.ndarray, edges: Sequence[float]) -> np.ndarray:
    """Front time in every cell [edges[k], edges[k+1]) of the screen coordinate y (nan if unreadable)."""
    out = np.full(len(edges) - 1, np.nan)
    t = Q[:, 0]
    for k in range(len(edges) - 1):
        m = (y >= edges[k]) & (y < edges[k + 1])
        if not m.any():
            continue
        tf = t[m & future]
        if tf.size == 0:
            continue
        first = tf.min()
        tn = t[m & ~future & (t < first)]
        if tn.size == 0:
            continue
        out[k] = 0.5 * (first + tn.max())
    return out


def local_tilt(yc: np.ndarray, delay: np.ndarray, half: int = 2) -> np.ndarray:
    """-d(delay)/dy by a linear fit over 2*half+1 cells around each cell."""
    out = np.full(len(yc), np.nan)
    for k in range(len(yc)):
        a, b = max(0, k - half), min(len(yc), k + half + 1)
        yy, dd = yc[a:b], delay[a:b]
        ok = np.isfinite(dd)
        if ok.sum() >= 3:
            out[k] = -np.polyfit(yy[ok], dd[ok], 1)[0]
    return out
