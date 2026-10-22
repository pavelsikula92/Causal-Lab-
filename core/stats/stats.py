"""
core.stats.stats
================

Small statistics shared by the tests.

    mean_se(xs)                  mean and standard error (nan-safe)
    boot_se(values, B, seed)     mean, bootstrap SE (resampling whole units) and 95 % percentiles
    equivalence(m, se, margin, n_se)
                                 EQUIVALENT / EXCESS / INCONCLUSIVE against a margin around 0
    slope_loglog(x, y)           least-squares slope of log y against log x
    fmt(x, nd)                   a number for a table ('-' for nan)
"""

from __future__ import annotations

import math
import statistics
from typing import Sequence, Tuple

import numpy as np


def mean_se(xs: Sequence[float]) -> Tuple[float, float]:
    """Mean and standard error (nan-safe; SE nan for fewer than 2 values)."""
    v = [float(x) for x in xs if x is not None and not (isinstance(x, float) and math.isnan(x))]
    if not v:
        return float("nan"), float("nan")
    m = statistics.mean(v)
    se = statistics.stdev(v) / math.sqrt(len(v)) if len(v) > 1 else float("nan")
    return m, se


def boot_se(values, B: int, seed: int):
    """Mean, bootstrap SE (sd of B resampled means of whole units) and the 2.5 / 97.5 % percentiles."""
    v = np.asarray([x for x in values if np.isfinite(x)], float)
    if v.size < 2:
        return (float(v.mean()) if v.size else float("nan")), float("nan"), (float("nan"), float("nan"))
    rng = np.random.default_rng(seed)
    m = v[rng.integers(0, v.size, size=(B, v.size))].mean(axis=1)
    return float(v.mean()), float(m.std(ddof=1)), (float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5)))


def equivalence(mean: float, se: float, margin: float, n_se: float) -> str:
    """EQUIVALENT if |mean| + n_se se <= margin; EXCESS if |mean| - n_se se > margin; else INCONCLUSIVE."""
    if not (np.isfinite(mean) and np.isfinite(se)):
        return "INCONCLUSIVE"
    if abs(mean) + n_se * se <= margin:
        return "EQUIVALENT"
    if abs(mean) - n_se * se > margin:
        return "EXCESS"
    return "INCONCLUSIVE"


def slope_loglog(x, y) -> float:
    """Least-squares slope of log y against log x over the positive, finite pairs (nan if < 2)."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    if m.sum() < 2:
        return float("nan")
    return float(np.polyfit(np.log(x[m]), np.log(y[m]), 1)[0])


def fmt(x, nd: int = 3) -> str:
    """A number for a table, '-' for nan/None."""
    if x is None:
        return "-"
    try:
        xf = float(x)
    except (TypeError, ValueError):
        return str(x)
    if math.isnan(xf) or math.isinf(xf):
        return "-"
    return f"{xf:.{nd}f}"
