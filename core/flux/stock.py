"""
core.flux.stock
===============

Capacity-limited routing in which the STOCK is a carried state (F3.3b).

In core.flux.capacity (F3.1, F3.2, F3.3) the excess of an event over its
capacity is stored as a local increment dQ_v and never acts again: the
model has no memory. Here the stock moves forward along the causal links
and is offered again downstream (author's decision, 21 Oct 2026: the
smallest extension that adds memory, because axiom 3 already contains the
stock; no inheritance of Pi, no rigidity I, no feedback on Pi).

Rule at event v (one pass over a linear extension of the causal order):
    F~_v  = s_v + sum_u Phi_{u->v}            offered throughput (source + flux of the parents)
    Q~_v  = sum_u Psi_{u->v}                  stock arriving from the parents
    A_v   = min(F~_v + Q~_v, C)               accepted: the new flux first, the stock fills what is left
    Q_v   = F~_v + Q~_v - A_v                 stock held by v
    Phi_{v->k} = A_v / d+_v,  Psi_{v->k} = Q_v / d+_v   (the reference uniform law for both)
Events that forward nothing (the outflow boundary |x| > R_abs and the
maximal events) pass neither flux nor stock: both leave the region.

Balance (axiom 3, exact, every event):
    s_v + sum_u Phi_uv + sum_u Psi_uv = sum_k Phi_vk + sum_k Psi_vk   (forwarding events)
so over the region  sum_v s_v = sum over non-forwarding events of (A_v + Q_v).

Order of evaluation and foliation
---------------------------------
The state of v depends only on its parents' outputs, so ANY linear
extension of the causal order gives exactly the same result: the rule
carries no preferred foliation (no global time step "all parents, then
all children"). `order` lets a test evaluate in another linear extension
(for example the time of a boosted frame) and compare; the routine checks
that every parent precedes its child in the order given.

Saturation indicator: chi_v = 1 if the event holds stock (Q_v > 0).
With the stock discarded instead of carried (keep_stock=False) this is
exactly route_capacity (F3.2): A_v = min(F~_v, C), dQ_v = F~_v - A_v.
"""

from __future__ import annotations

from typing import Dict, Optional

import numpy as np


def route_stock(ptr: np.ndarray, idx: np.ndarray, s: np.ndarray, C: float, absorb: Optional[np.ndarray],
                keep_stock: bool = True, order: Optional[np.ndarray] = None) -> Dict[str, np.ndarray]:
    """Routing with carried stock (module notes). Returns offered, stock_in, accepted, stock, Phi, Psi, forwards, leaving."""
    N = len(ptr) - 1
    s = np.asarray(s, float)
    dplus = np.bincount(idx, minlength=N)
    fwd = dplus > 0
    if absorb is not None:
        fwd &= ~np.asarray(absorb, bool)
    share = np.where(fwd, 1.0 / np.maximum(dplus, 1), 0.0)
    if order is None:
        order = np.arange(N)
    else:
        order = np.asarray(order)
        pos = np.empty(N, dtype=np.int64)
        pos[order] = np.arange(N)
        v_of = np.repeat(np.arange(N), np.diff(ptr))
        if np.any(pos[idx] >= pos[v_of]):
            raise ValueError("order is not a linear extension of the causal order (a parent after its child)")
    offered = s.copy()
    stock_in = np.zeros(N)
    acc = np.zeros(N)
    Q = np.zeros(N)
    Phi = np.zeros(idx.size)
    Psi = np.zeros(idx.size)
    for v in order:
        a, b = ptr[v], ptr[v + 1]
        if b > a:
            par = idx[a:b]
            Phi[a:b] = acc[par] * share[par]
            offered[v] += Phi[a:b].sum()
            if keep_stock:
                Psi[a:b] = Q[par] * share[par]
                stock_in[v] = Psi[a:b].sum()
        acc[v] = min(offered[v] + stock_in[v], C)
        Q[v] = offered[v] + stock_in[v] - acc[v]
    leaving = np.where(fwd, 0.0, acc + (Q if keep_stock else 0.0))
    return {"offered": offered, "stock_in": stock_in, "accepted": acc, "stock": Q, "Phi": Phi, "Psi": Psi,
            "forwards": fwd, "leaving": leaving}


def stock_balance(s: np.ndarray, r: Dict[str, np.ndarray], keep_stock: bool = True) -> float:
    """Relative residual of the regional balance: sum s = leaving (+ local stock if it is not carried)."""
    tot = float(np.sum(s))
    out = float(r["leaving"].sum()) + (0.0 if keep_stock else float(r["stock"].sum()))
    return abs(tot - out) / max(tot, 1e-300)
