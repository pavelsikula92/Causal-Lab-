"""
core.flux.pulse
===============

Capacity-limited routing with an INTENSIVE, copied state: the pulse (F3.3c).

F3.3b showed that a carried stock is real memory but dilutes: it is split
over the children like the flux (sum out = in) and spreads over the future
cone within a step or two. The author's hypothesis (21 Oct 2026): the
state of a persistent structure is not divided along the links but COPIED,
like a frequency or a phase; every child receives the same value however
many children there are. It is a state of the event, not a conserved
quantity, so axiom 3 (which balances the flux) does not apply to it; the
flux itself keeps its exact balance.

Rule at event v, one pass over any linear extension of the causal order
(the state of v depends only on its parents, so the result does not depend
on the order; tests check a boosted-frame order bit for bit):

  inherited pulse  P_in(v) = combine over the parents u -> v of P(u):
      'max'   the largest P(u) (a dominant excitation);
      'main'  P(u*) of the main parent u* = argmax_u Phi(u -> v), the parent
              that sends v the most flux (ties: the lowest index; no flux
              in, no main parent, P_in = 0): the pulse follows the carrying flux;
      'none'  P_in = 0 (no pulse at all: exactly route_capacity, F3.2).
  effect   the pulse acts as internal events: s_eff(v) = s(v) + P_in(v);
  load     offered F~ = s_eff + inflow; accepted A = min(F~, C); dQ = F~ - A
           (local, not carried: the pulse is the only memory, F3.2 otherwise);
  survival and birth   the pulse is alive at v if dQ > 0 (dq_crit = 0, F3.3c)
           or dQ >= dq_crit (dq_crit > 0, F3.3d: the materialization
           threshold). Where it is alive, v keeps the inherited pulse or
           creates one of value omega if it inherited none; elsewhere P = 0
           (the pulse dies at once below the threshold: no decay parameter).
           So the pulse lives only where the excess load is high enough, and
           its value never changes along the links (omega or 0): it is not
           diluted.

Feedback of the clock rate (F3.3f, optional; channel = 'none' is exactly
the rule above, bit for bit): the load sets Pi_v by the response of axiom 5,
Pi = f(sigma) = 1 - (1 - Pi_sat) min(sigma, 1)^m (core.flux.feedback), and
ONE channel acts back, as in F3.3:
  'capacity'  C_v = C Pi_v: the largest solution of Pi = f(min(F~ / (C Pi), 1)),
              F~ = s + P_in + inflow (slow clocks process less: positive feedback);
  'source'    the internal events, the pulse included, run at the clock rate:
              s_eff = (s + P_in) Pi, Pi the unique solution of
              Pi = f(min(((s + P_in) Pi + inflow) / C, 1)) (negative feedback).
Pi_v depends only on v's own load, so the pass stays order-free.

Narrowed future cones (F3.3h, optional; narrow_power = 0 is everything above,
bit for bit): the author's hypothesis that a slow clock narrows the causal
cone of the event's children (axiom 7 in its smallest, local form). An event
v forwards its accepted flux only to the children k inside its narrowed cone
    |x_k - x_v| <= Pi_v^p (t_k - t_v),    p = narrow_power,
uniformly among them (F7: full light bending needs c_eff = c Pi^2, so the
laboratory uses p = 2). Pi_v is v's own clock rate, known when v forwards
(with channel 'none' it is the response Pi = f(min(F~ / C, 1)) of axiom 5,
recorded without acting on the capacity). With Pi_v = 1 the narrowed cone
is the whole cone (every link is causal), so nothing changes. If no linked
child lies inside (a sprinkling has few nearly timelike links), v forwards
to its single most timelike child (smallest |dx| / dt), the linked child
closest to the narrowed cone: the flux stays in the medium (axiom 3 holds
locally; letting it leave would be a sink inside the region); such events
are counted ('narrow_fallback'). The cone is drawn in the coordinates of the sprinkling:
this is a frame-dependent stand-in for a metric that narrows the cones
itself, not the rule of the theory (tests say so). The pass stays causal
and order-free (the allowed children of v depend only on v).

Balance (axiom 3 for the flux, exact): sum_v s_eff(v) = sum_v dQ_v +
sum over non-forwarding events of A_v. The pulse's own injection is part of
s_eff (internal events are sources).
"""

from __future__ import annotations

from typing import Dict, Optional

import numpy as np

from .feedback import _pi_capacity, _pi_source, response

RULES = ("none", "max", "main")
CHANNELS = ("none", "capacity", "source")


def route_pulse(ptr: np.ndarray, idx: np.ndarray, s: np.ndarray, C: float, absorb: Optional[np.ndarray],
                rule: str, omega: float, order: Optional[np.ndarray] = None, dq_crit: float = 0.0,
                channel: str = "none", Pi_sat: float = 1.0, m: float = 2.0, narrow_power: float = 0.0,
                X: Optional[np.ndarray] = None) -> Dict[str, np.ndarray]:
    """Routing with the copied pulse (module notes). Returns offered, accepted, stored, Phi, P, P_in, s_eff, leaving,
    Pi and capacity (Pi = 1, capacity = C without a feedback channel), and with narrowed cones narrow_fallback (events
    with no linked child inside the narrowed cone, forwarding to their most timelike child) and n_allowed (children
    inside, -1 where not evaluated)."""
    narrow = float(narrow_power) > 0.0
    if narrow and X is None:
        raise ValueError("narrowed cones need the coordinates X of the events")
    if rule not in RULES:
        raise ValueError(f"rule must be one of {RULES}")
    if channel not in CHANNELS:
        raise ValueError(f"channel must be one of {CHANNELS}")
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
    offered = np.zeros(N)
    acc = np.zeros(N)
    s_eff = s.copy()
    P = np.zeros(N)
    P_in = np.zeros(N)
    Phi = np.zeros(idx.size)
    Pi = np.ones(N)
    cap = np.full(N, float(C))
    link_w = share[idx].copy()                        # weight of each link: the parent's share (uniform without narrowing)
    fallback = np.zeros(N, dtype=bool)
    n_allowed = np.full(N, -1, dtype=np.int64)
    if narrow:
        X = np.asarray(X, float)
        v_link = np.repeat(np.arange(N), np.diff(ptr))  # child of every link (idx holds the parent)
        out = np.argsort(idx, kind="stable")           # links grouped by parent
        optr = np.concatenate([[0], np.cumsum(np.bincount(idx, minlength=N))])
        link_w[:] = 0.0                                # set when the parent is evaluated
    for v in order:
        a, b = ptr[v], ptr[v + 1]
        inflow, pin = 0.0, 0.0
        if b > a:
            par = idx[a:b]
            ph = acc[par] * link_w[a:b]
            Phi[a:b] = ph
            inflow = float(ph.sum())
            if rule == "max":
                pin = float(P[par].max())
            elif rule == "main" and inflow > 0.0:
                j = int(np.argmax(ph))                      # first maximum: the lowest position among equal flux
                cand = par[ph == ph[j]]
                pin = float(P[int(cand.min())])            # ties broken by the lowest event index (order-free)
        P_in[v] = pin
        s_eff[v] = s[v] + pin
        if channel == "capacity":
            Pi[v] = _pi_capacity(s_eff[v] + inflow, C, Pi_sat, m)
            cap[v] = C * Pi[v]
        elif channel == "source":
            Pi[v] = _pi_source(inflow, s_eff[v], C, Pi_sat, m)
            s_eff[v] = s_eff[v] * Pi[v]
        offered[v] = s_eff[v] + inflow
        if narrow and channel == "none":
            Pi[v] = float(response(offered[v] / C, Pi_sat, m))
        acc[v] = min(offered[v], cap[v])
        if narrow and fwd[v]:
            lk = out[optr[v]:optr[v + 1]]               # v's out-links
            kids = v_link[lk]
            dt = X[kids, 0] - X[v, 0]
            dx = np.sqrt(((X[kids, 1:] - X[v, 1:]) ** 2).sum(axis=1))
            ok = dx <= (Pi[v] ** narrow_power) * dt
            n = int(ok.sum())
            n_allowed[v] = n
            if n:
                link_w[lk[ok]] = 1.0 / n
            else:                                     # nothing inside: the most timelike child (lowest index on ties)
                j = int(np.lexsort((kids, dx / dt))[0])
                link_w[lk[j]] = 1.0
                fallback[v] = True
        dq = offered[v] - acc[v]
        alive = (dq > 0.0) if dq_crit <= 0.0 else (dq >= dq_crit)
        if alive and rule != "none":
            P[v] = pin if pin > 0.0 else float(omega)
    dQ = offered - acc
    leaving = np.where(fwd, 0.0, acc)
    return {"offered": offered, "accepted": acc, "stored": dQ, "Phi": Phi, "P": P, "P_in": P_in, "s_eff": s_eff,
            "forwards": fwd, "leaving": leaving, "Pi": Pi, "capacity": cap, "narrow_fallback": fallback, "n_allowed": n_allowed}


def pulse_balance(r: Dict[str, np.ndarray]) -> float:
    """Relative residual of sum s_eff = sum dQ + sum leaving (axiom 3 for the flux)."""
    tot = float(r["s_eff"].sum())
    return abs(tot - float(r["stored"].sum()) - float(r["leaving"].sum())) / max(tot, 1e-300)
