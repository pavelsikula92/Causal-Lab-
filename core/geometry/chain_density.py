"""
core.geometry.chain_density
===========================

The clock rate Pi read from the order (F2.1 (E44b) in 1+1, K.3 (E45c) in 2+1 and 3+1).

In a sprinkling of uniform PROPER density rho, the longest chain of an
Alexandrov interval of proper time tau has length L ~ m_D (rho V)^(1/D),
V = kappa_D tau^D (Brightwell-Gregory; Bollobas-Winkler), i.e. L is
proportional to tau. A small diamond at rest of coordinate height H at a
place with clock rate Pi has proper time Pi H, so L(x) is proportional to
Pi(x): the chain density is the clock. The proportionality constant (and
its finite-size correction) cancels in the ratio L(x) / mean L, which is
what is compared with Pi(x) / mean Pi.
"""

from __future__ import annotations

import numpy as np

from ..order.relation import chain_dp


def longest_chain_length(X: np.ndarray) -> int:
    """Number of elements of the longest chain of the time-sorted points X (0 if empty)."""
    if len(X) == 0:
        return 0
    return int(chain_dp(X)["W"].max())


def readback(L_by_site, Pi_by_site):
    """
    Compare the normalised readings L/mean(L) with Pi/mean(Pi) over the
    sites: Pearson r and the largest relative error (F2.1 (E44b) statistic).
    """
    L = np.asarray(L_by_site, float)
    P = np.asarray(Pi_by_site, float)
    a, b = L / L.mean(), P / P.mean()
    if np.std(a) == 0 or np.std(b) == 0:
        r = float("nan")
    else:
        r = float(np.corrcoef(a, b)[0, 1])
    err = float(np.max(np.abs(a - b) / b))
    return {"r": r, "max_rel_err": err, "pi_order": a.tolist(), "pi_true": b.tolist()}
