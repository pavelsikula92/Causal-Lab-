"""
core.flux.routing
=================

Flux on the links of an order (the flux balance of Causal Physics 2.0):

    sum_in Phi + s_v = sum_out Phi + dQ_v          (Axiom: flux balance)

with the F1.8 (E41-A) UNIFORM routing law and no storage (dQ_v = 0 except at
maximal elements, whose throughput leaves the region):

    F_v = s_v + sum_{u -> v} Phi_uv,     Phi_uv = keep * F_u / d_u^+

d_u^+ is the number of future links (children) of u. keep = 1 is the
conserving law. keep < 1 is the F1.1 (E45a) NULL MODEL: every element discards
a fraction 1 - keep of what it forwards, which breaks continuity.

Because routing is linear in s, the flux of a source alone ("net
current" j of the source) is obtained by routing the source vector alone;
the vacuum flux of minimal elements is a separate vector (F1.1 (E45a) uses it as
the instrument check: by isotropy its net flux through a sphere is 0).

Crossings of a timelike cylinder
--------------------------------
`side_crossings` places each link on the straight segment between its
end events and finds where the segment crosses the cylinder |x| = r (0, 1
or 2 crossings), with the crossing time and direction. The net outward
flux through the side of the cylinder during [t1, t2] is

    J(r) = sum over outward crossings with t in [t1, t2] of Phi
         - sum over inward crossings with t in [t1, t2] of Phi

`node_set_balance` is the exact discrete Gauss identity for any node set
B: sum_B s = (flux out of B) - (flux into B) + (throughput of maximal
elements in B) - (what the null discards inside B). It holds to machine
precision for every routing (core test `routing_identities`).
"""

from __future__ import annotations

from typing import Dict, Optional, Sequence, Tuple

import numpy as np


def children_counts(ptr: np.ndarray, idx: np.ndarray, N: int) -> np.ndarray:
    """Number of future links of every element."""
    return np.bincount(idx, minlength=N).astype(np.int64)


def route_uniform(ptr: np.ndarray, idx: np.ndarray, s: np.ndarray, keep: float = 1.0,
                  absorb: Optional[np.ndarray] = None) -> Tuple[np.ndarray, np.ndarray]:
    """
    Throughput F (length N) and the flux carried by each link in CSR
    order (Phi, aligned with idx). Elements are processed in time order,
    which is a linear extension, so every parent is final when used.
    `absorb` (mask): elements that forward nothing -- an outflow boundary
    standing in for the unbounded space outside the region (their
    throughput leaves, like that of maximal elements).
    """
    N = len(ptr) - 1
    dplus = children_counts(ptr, idx, N)
    fwd = dplus > 0
    if absorb is not None:
        fwd &= ~np.asarray(absorb, bool)
    out = np.where(fwd, float(keep) / np.maximum(dplus, 1), 0.0)   # share per child link
    F = np.asarray(s, float).copy()
    for v in range(N):
        a, b = ptr[v], ptr[v + 1]
        if b > a:
            par = idx[a:b]
            F[v] += float((F[par] * out[par]).sum())
    Phi = np.empty(idx.size)
    for v in range(N):
        a, b = ptr[v], ptr[v + 1]
        if b > a:
            par = idx[a:b]
            Phi[a:b] = F[par] * out[par]
    return F, Phi


def link_endpoints(ptr: np.ndarray, idx: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """(u, v) arrays of all links u -> v in CSR order."""
    N = len(ptr) - 1
    v = np.repeat(np.arange(N), np.diff(ptr))
    return idx.copy(), v


def side_crossings(X: np.ndarray, u: np.ndarray, v: np.ndarray, r: float,
                   center: Optional[Sequence[float]] = None, chunk: int = 2_000_000
                   ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Crossings of the straight link segments u -> v with the cylinder
    |x - center| = r. Returns (link index, crossing time, sign) with sign
    +1 outward, -1 inward. A segment through the cylinder twice gives two
    crossings (in and out).
    """
    li, tc, sg = [], [], []
    c = np.zeros(X.shape[1] - 1) if center is None else np.asarray(center, float)
    for c0 in range(0, u.size, chunk):
        uu, vv = u[c0:c0 + chunk], v[c0:c0 + chunk]
        a = X[uu, 1:] - c
        d = X[vv, 1:] - X[uu, 1:]
        A = (d * d).sum(axis=1)
        B = 2.0 * (a * d).sum(axis=1)
        C = (a * a).sum(axis=1) - r * r
        disc = B * B - 4.0 * A * C
        ok = (A > 0) & (disc > 0)
        sq = np.sqrt(np.where(ok, disc, 0.0))
        for sgn in (-1.0, 1.0):
            with np.errstate(divide="ignore", invalid="ignore"):
                s = (-B + sgn * sq) / (2.0 * np.where(A > 0, A, 1.0))
            hit = ok & (s >= 0.0) & (s < 1.0)
            k = np.nonzero(hit)[0]
            if k.size:
                t = X[uu[k], 0] + s[k] * (X[vv[k], 0] - X[uu[k], 0])
                # d|x|^2/ds = B + 2 A s: positive when moving outward
                dirn = np.sign(B[k] + 2.0 * A[k] * s[k])
                li.append(c0 + k)
                tc.append(t)
                sg.append(dirn)
    if not li:
        z = np.zeros(0)
        return z.astype(np.int64), z, z
    return np.concatenate(li), np.concatenate(tc), np.concatenate(sg)


def cap_flux(X: np.ndarray, u: np.ndarray, v: np.ndarray, Phi: np.ndarray, t: float, r: float,
             center=None) -> float:
    """Flux on links crossing the slice at time t inside |x| < r (upward)."""
    tu, tv = X[u, 0], X[v, 0]
    m = (tu < t) & (tv >= t)
    if not m.any():
        return 0.0
    c = np.zeros(X.shape[1] - 1) if center is None else np.asarray(center, float)
    f = (t - tu[m]) / (tv[m] - tu[m])
    pos = X[u[m], 1:] + f[:, None] * (X[v[m], 1:] - X[u[m], 1:]) - c
    ins = (pos ** 2).sum(axis=1) < r * r
    return float(Phi[m][ins].sum())


def stock_profile(X: np.ndarray, u: np.ndarray, v: np.ndarray, Phi: np.ndarray,
                  times, radii, center=None) -> np.ndarray:
    """
    Stock in transit: flux on links crossing the slice at time t inside
    |x| < r, for every (t, r) at once (rows = times, columns = radii).
    For a unit impulse this is P_inside(t; r), the probability that flux
    emitted at the source is still inside radius r at time t.

    Use (F1.1 (E45a), settling): in a linear, time-invariant routing a source
    switched on at t = 0 has net outward flux J(t; r) = S (1 - P(t; r)),
    so the window mean of P is the settling deficit of J/S expected in
    that window. Each row equals cap_flux(X, u, v, Phi, t, r).
    """
    tu, tv = X[u, 0], X[v, 0]
    c = np.zeros(X.shape[1] - 1) if center is None else np.asarray(center, float)
    radii = np.asarray(radii, float)
    out = np.zeros((len(times), radii.size))
    for k, t in enumerate(times):
        m = (tu < t) & (tv >= t)
        if not m.any():
            continue
        f = (t - tu[m]) / (tv[m] - tu[m])
        pos = X[u[m], 1:] + f[:, None] * (X[v[m], 1:] - X[u[m], 1:]) - c
        d = np.sqrt((pos ** 2).sum(axis=1))
        ph = Phi[m]
        order = np.argsort(d)
        cum = np.r_[0.0, np.cumsum(ph[order])]
        # flux strictly inside each radius
        out[k] = cum[np.searchsorted(d[order], radii, side="left")]
    return out


def net_side_flux(X: np.ndarray, u: np.ndarray, v: np.ndarray, Phi: np.ndarray, r: float,
                  t1: float, t2: float, center=None) -> Dict[str, float]:
    """Outward, inward and net flux through the side |x| = r during [t1, t2]."""
    li, tc, sg = side_crossings(X, u, v, r, center)
    w = (tc >= t1) & (tc <= t2)
    out = float(Phi[li[w & (sg > 0)]].sum())
    inn = float(Phi[li[w & (sg < 0)]].sum())
    return {"out": out, "in": inn, "net": out - inn, "n_cross": int(w.sum())}


def node_set_balance(ptr: np.ndarray, idx: np.ndarray, s: np.ndarray, F: np.ndarray,
                     Phi: np.ndarray, inside: np.ndarray, keep: float = 1.0,
                     absorb: Optional[np.ndarray] = None) -> Dict[str, float]:
    """
    The exact discrete Gauss identity for the node set `inside` (mask):
        sum_B s  =  out(B) - in(B) + sinks(B) + discarded(B)
    out(B)/in(B): flux on links leaving / entering B; sinks(B): throughput
    of elements of B with no children; discarded(B): (1 - keep) F of the
    elements of B with children. Returns both sides and the residual.
    """
    N = len(ptr) - 1
    u, v = link_endpoints(ptr, idx)
    inside = np.asarray(inside, bool)
    out_flux = float(Phi[inside[u] & ~inside[v]].sum())
    in_flux = float(Phi[~inside[u] & inside[v]].sum())
    dplus = children_counts(ptr, idx, N)
    sink = dplus == 0
    if absorb is not None:
        sink = sink | np.asarray(absorb, bool)
    sinks = float(F[inside & sink].sum())
    disc = float(((1.0 - keep) * F)[inside & ~sink].sum())
    lhs = float(np.asarray(s, float)[inside].sum())
    rhs = out_flux - in_flux + sinks + disc
    return {"sources": lhs, "out": out_flux, "in": in_flux, "sinks": sinks,
            "discarded": disc, "residual": lhs - rhs}
