"""
core.flux.feedback
==================

Capacity-limited routing WITH feedback of the load on the transport
(F3.3). It extends core.flux.capacity.route_capacity (F3.1, F3.2: no
feedback) by the response of axiom 5,

    sigma_v = F_v / C_v (capped at 1),     Pi_v = f(sigma_v) = 1 - (1 - Pi_sat) sigma_v^m ,

and ONE channel through which Pi acts back on the transport (the author
decided on 21 Oct 2026 to test the channels separately, never several at
once):

  'off'       no feedback: exactly route_capacity (Pi is only recorded);
  'capacity'  slow clocks process less: C_v = C Pi_v (A);
  'routing'   flux prefers low impedance: an event forwards its accepted flux to
              its children with weights Pi_k / I_k, normalised over its
              children (B; rigidity I_k = 1 until F0.5 defines it, so the
              weights are Pi_k / sum Pi);
  'source'    internal events run at the local clock rate: s_v Pi_v (C).

Flux balance (axiom 3) holds exactly in every mode:
    sum_v s_eff_v = sum_v dQ_v + sum_v outflow_v .

Locality and the order of evaluation
------------------------------------
'off', 'capacity' and 'source' are ONE pass in causal (time) order: the
state of v depends only on its causal past (the flux its parents sent) and
on itself. Pi_v then solves a local fixed-point equation, because the load
depends on Pi_v:
    capacity: Pi = f(min(F~ / (C Pi), 1)),  F~ = s + inflow;
    source:   Pi = f(min((s Pi + inflow) / C, 1)).
For 'source' the solution is unique (the right side decreases in Pi). For
'capacity' the right side increases in Pi and there can be two solutions
(an unloaded and a saturated branch: the positive feedback); the rule
takes the LARGEST, i.e. the state reached by relaxing from the vacuum
Pi = 1 (monotone iteration from Pi = 1 converges to it). An event jumps to
the saturated branch only when the unloaded branch no longer exists.

'routing' cannot be a causal pass: when a parent forwards its flux, the
state of its children is not yet determined (they are in its future), so
a rule that weights children by their Pi is not a function of the causal
past. It is computed as the STATIONARY self-consistent state: routing
passes repeated with the children's Pi from the previous pass (damped,
Pi <- (1 - damp) Pi + damp Pi_new) until the largest change of Pi is
below tol. Starting from Pi = 1 the first pass is exactly 'off'. Whether
it converges is reported (instrument). This is a statement about
stationary states, not about an update rule; tests say so.

With Pi_sat = 1 (f = 1) every mode reproduces route_capacity exactly
(core check `feedback_reduces_to_capacity`).
"""

from __future__ import annotations

from typing import Dict, Optional

import numpy as np

MODES = ("off", "capacity", "routing", "source")


def response(sigma, Pi_sat: float, m: float):
    """Pi = 1 - (1 - Pi_sat) min(sigma, 1)^m (axiom 5 response; Pi = Pi_sat at full load)."""
    s = np.minimum(np.asarray(sigma, float), 1.0)
    return 1.0 - (1.0 - Pi_sat) * s ** m


def _largest_root(h, lo: float, hi: float = 1.0, n_scan: int = 48, n_bisect: int = 50) -> float:
    """Largest x in [lo, hi] with h(x) = 0, given h(hi) >= 0 and h(lo) <= 0 (scan down, then bisect)."""
    xs = np.linspace(hi, lo, n_scan + 1)
    prev = hi
    for x in xs[1:]:
        if h(x) < 0.0:
            a, b = x, prev                      # h(a) < 0 <= h(b): the largest root lies in [a, b]
            for _ in range(n_bisect):
                c = 0.5 * (a + b)
                if h(c) < 0.0:
                    a = c
                else:
                    b = c
            return b
        prev = x
    return lo                                   # h >= 0 down to lo: lo itself is the root (h(lo) = 0)


def _pi_capacity(F: float, C: float, Pi_sat: float, m: float) -> float:
    """Largest solution of Pi = f(min(F / (C Pi), 1)) (relaxation from Pi = 1)."""
    if F <= 0.0 or Pi_sat >= 1.0:
        return 1.0
    k = 1.0 - Pi_sat
    h = lambda x: x - (1.0 - k * min(F / (C * x), 1.0) ** m)
    return _largest_root(h, Pi_sat)


def _pi_source(inflow: float, s: float, C: float, Pi_sat: float, m: float) -> float:
    """Unique solution of Pi = f(min((s Pi + inflow) / C, 1)) (the right side decreases in Pi)."""
    if Pi_sat >= 1.0 or (s <= 0.0 and inflow <= 0.0):
        return 1.0
    k = 1.0 - Pi_sat
    h = lambda x: x - (1.0 - k * min((s * x + inflow) / C, 1.0) ** m)   # increasing in x
    a, b = Pi_sat, 1.0
    if h(a) >= 0.0:
        return a
    for _ in range(60):
        c = 0.5 * (a + b)
        if h(c) < 0.0:
            a = c
        else:
            b = c
    return b


def _pass(ptr, idx, s, C, fwd, wl, mode, Pi_sat, m):
    """One routing pass in time order with per-link weights wl (fraction of the parent's accepted flux)."""
    N = len(ptr) - 1
    acc = np.zeros(N)
    offered = np.zeros(N)
    s_eff = np.asarray(s, float).copy()
    cap = np.full(N, float(C))
    Pi = np.ones(N)
    Phi = np.zeros(idx.size)
    send = np.where(fwd, 1.0, 0.0)               # absorbed / maximal events forward nothing
    for v in range(N):
        a, b = ptr[v], ptr[v + 1]
        inflow = 0.0
        if b > a:
            par = idx[a:b]
            ph = acc[par] * send[par] * wl[a:b]
            Phi[a:b] = ph
            inflow = float(ph.sum())
        sv = float(s_eff[v])
        if mode == "capacity":
            p = _pi_capacity(sv + inflow, C, Pi_sat, m)
            cap[v] = C * p
        elif mode == "source":
            p = _pi_source(inflow, sv, C, Pi_sat, m)
            s_eff[v] = sv * p
        else:                                    # off / routing: capacity C, Pi recorded from the load
            p = float(response((sv + inflow) / C, Pi_sat, m))
        Pi[v] = p
        offered[v] = s_eff[v] + inflow
        acc[v] = min(offered[v], cap[v])
    return acc, offered, s_eff, cap, Pi, Phi


def route_feedback(ptr: np.ndarray, idx: np.ndarray, s: np.ndarray, C: float, absorb: Optional[np.ndarray],
                   mode: str, Pi_sat: float, m: float, I: Optional[np.ndarray] = None, tol: float = 1e-6,
                   max_iter: int = 200, damp: float = 0.5) -> Dict[str, np.ndarray]:
    """
    Capacity routing with one feedback channel (module notes). Returns the
    fields of route_capacity (offered, accepted, stored, outflow, Phi,
    forwards) plus Pi, the effective source s_eff, the capacity per event,
    and for 'routing' the number of passes and whether it converged.
    """
    if mode not in MODES:
        raise ValueError(f"mode must be one of {MODES}")
    N = len(ptr) - 1
    dplus = np.bincount(idx, minlength=N)
    fwd = dplus > 0
    if absorb is not None:
        fwd &= ~np.asarray(absorb, bool)
    v_of = np.repeat(np.arange(N), np.diff(ptr))              # child of every link (idx holds the parent)
    uniform = 1.0 / np.maximum(dplus[idx], 1)
    iters, converged, last_change = 1, True, 0.0
    if mode != "routing":
        acc, offered, s_eff, cap, Pi, Phi = _pass(ptr, idx, s, C, fwd, uniform, mode, Pi_sat, m)
    else:
        Ik = np.ones(N) if I is None else np.asarray(I, float)
        Pi_w = np.ones(N)                                     # the children's Pi used for the weights
        converged = False
        for it in range(1, max_iter + 1):
            q = Pi_w / Ik
            norm = np.bincount(idx, weights=q[v_of], minlength=N)
            wl = q[v_of] / np.where(norm[idx] > 0, norm[idx], 1.0)
            acc, offered, s_eff, cap, Pi, Phi = _pass(ptr, idx, s, C, fwd, wl, "routing", Pi_sat, m)
            last_change = float(np.max(np.abs(Pi - Pi_w)))
            iters = it
            if last_change < tol:
                converged = True
                break
            Pi_w = (1.0 - damp) * Pi_w + damp * Pi
    dQ = offered - acc
    outflow = np.where(fwd, 0.0, acc)
    return {"offered": offered, "accepted": acc, "stored": dQ, "outflow": outflow, "Phi": Phi, "forwards": fwd,
            "Pi": Pi, "s_eff": s_eff, "capacity": cap, "iterations": iters, "converged": converged,
            "last_change": last_change}


def balance_residual_fb(r: Dict[str, np.ndarray]) -> float:
    """|sum s_eff - sum dQ - sum outflow| / sum s_eff (axiom 3 over the whole region)."""
    tot = float(r["s_eff"].sum())
    return abs(tot - float(r["stored"].sum()) - float(r["outflow"].sum())) / max(tot, 1e-300)
