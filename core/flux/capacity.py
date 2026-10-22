"""
core.flux.capacity
==================

Capacity-limited routing on a FIXED order (layer 2, step c; register
R139-R141, F3.1 (E55a)/F3.2 (E55b)). No feedback: the routing does not depend on the
load (that is F3.3 (E56)).

Semantics (decided by the author on 17-18 Oct 2026)
---------------------------------------------------
Flux balance of axiom 3 at every event v:

    sum_u Phi_uv + s_v = sum_w Phi_vw + Delta Q_v

Routing runs ONCE in causal (time) order. The offered flux of v is
F~_v = s_v + sum_u Phi_uv; v accepts A_v = min(F~_v, C) and forwards it
uniformly to its children (sigma0_uv = 1 / d_u^+, the F1.8 (E41-A) law); the
rest is stored, Delta Q_v = F~_v - A_v. The stored amount is a LOCAL
increment at the event: it is not transported to later events (there is
no evolution rule for Q yet). Events beyond the outflow radius, and
maximal events, forward nothing: what they accept leaves the region.
Hence the exact balance

    sum_v s_v = sum_v Delta Q_v + sum_v outflow_v                 (balance)

Proportional acceptance: when F~_v > C every component of the offered
flux (background, probe) is accepted with the same ratio
alpha_v = A_v / F~_v. So a small marked probe injected at v is
transmitted by v with ratio alpha_v and travels downstream with the
same ratios. Its escape fraction (the share that leaves the region
instead of being stored) is, for every v at once, the backward pass

    e_v = alpha_v * (1                          if v forwards nothing,
                     sum_children e_c / d_v^+   otherwise)

which is the measured permeability Pi_v^meas = T_out / T_in of a probe
injected at v. probe_check verifies it against an explicit finite probe
routed forward with the nonlinear min (the probe must be small enough
not to change the background).

Topology control (rewired null)
-------------------------------
rewire_links: degree-preserving double-edge swaps that respect causal
order: (a -> b), (c -> d) become (a -> d), (c -> b) only if t_a < t_d and
t_c < t_b (all edges still point forward in time, so no cycle), the two
new edges do not exist yet, and a != c, b != d. Afterwards the edges
that are not covering relations of the new order (implied by a longer
path) are found exactly by descendant bitsets (transitive reduction) and
removed; the removed fraction and the change of the degree distribution
are reported, because removal changes the degrees.
"""

from __future__ import annotations

from typing import Any, Dict, Optional, Tuple

import numpy as np


def children_csr(ptr: np.ndarray, idx: np.ndarray, N: int):
    """Children CSR (cptr, cidx) from the parents CSR (ptr, idx)."""
    v_of = np.repeat(np.arange(N), np.diff(ptr))
    order = np.argsort(idx, kind="stable")
    cptr = np.zeros(N + 1, dtype=np.int64)
    np.add.at(cptr, idx + 1, 1)
    return np.cumsum(cptr), v_of[order]


def route_capacity(ptr: np.ndarray, idx: np.ndarray, s: np.ndarray, C: float,
                   absorb: Optional[np.ndarray] = None) -> Dict[str, np.ndarray]:
    """
    One pass in time order (module notes). Returns offered F~, accepted A,
    stored dQ, outflow, alpha = A / F~ (1 where nothing is offered), and
    the flux per link Phi (aligned with idx). C = inf gives the F1.1 (E45a) routing.
    """
    N = len(ptr) - 1
    s = np.asarray(s, float)
    dplus = np.bincount(idx, minlength=N)
    fwd = dplus > 0
    if absorb is not None:
        fwd &= ~np.asarray(absorb, bool)
    share = np.where(fwd, 1.0 / np.maximum(dplus, 1), 0.0)
    offered = s.copy()
    acc = np.zeros(N)
    Phi = np.zeros(idx.size)
    for v in range(N):
        a, b = ptr[v], ptr[v + 1]
        if b > a:
            par = idx[a:b]
            Phi[a:b] = acc[par] * share[par]
            offered[v] += Phi[a:b].sum()
        acc[v] = min(offered[v], C)
    dQ = offered - acc
    outflow = np.where(fwd, 0.0, acc)
    alpha = np.where(offered > 0, acc / np.where(offered > 0, offered, 1.0), 1.0)
    return {"offered": offered, "accepted": acc, "stored": dQ, "outflow": outflow,
            "alpha": alpha, "Phi": Phi, "forwards": fwd}


def escape_fraction(ptr: np.ndarray, idx: np.ndarray, alpha: np.ndarray, forwards: np.ndarray) -> np.ndarray:
    """Measured permeability e_v of a small probe injected at v (module notes), one backward pass."""
    N = len(ptr) - 1
    cptr, cidx = children_csr(ptr, idx, N)
    e = np.zeros(N)
    for v in range(N - 1, -1, -1):
        if not forwards[v]:
            e[v] = alpha[v]
        else:
            ch = cidx[cptr[v]:cptr[v + 1]]
            e[v] = alpha[v] * e[ch].mean()
    return e


def balance_residual(s: np.ndarray, r: Dict[str, np.ndarray]) -> float:
    """|sum s - sum dQ - sum outflow| / max(1, sum s)."""
    tot = float(np.sum(s))
    return abs(tot - float(r["stored"].sum()) - float(r["outflow"].sum())) / max(1.0, tot)


def route_marked(ptr: np.ndarray, idx: np.ndarray, alpha: np.ndarray, forwards: np.ndarray,
                 m: np.ndarray) -> Dict[str, np.ndarray]:
    """
    A MARKED component m (e.g. a probe) routed through the same pass with
    the acceptance ratios alpha of the joint routing (proportional
    acceptance): in_v = m_v + inflow, out_v = alpha_v in_v, forwarded
    uniformly like the background. Linear in m.
    """
    N = len(ptr) - 1
    dplus = np.bincount(idx, minlength=N)
    share = np.where(forwards, 1.0 / np.maximum(dplus, 1), 0.0)
    tin = np.asarray(m, float).copy()
    tout = np.zeros(N)
    for v in range(N):
        a, b = ptr[v], ptr[v + 1]
        if b > a:
            par = idx[a:b]
            tin[v] += float(np.dot(tout[par], share[par]))
        tout[v] = alpha[v] * tin[v]
    return {"T_in": tin, "T_out": tout, "stored": tin - tout,
            "escaped": float(np.where(forwards, 0.0, tout).sum())}


def probe_permeability(ptr: np.ndarray, idx: np.ndarray, s: np.ndarray, C: float, absorb,
                       delta_rel: float) -> Dict[str, np.ndarray]:
    """
    Measured local permeability Pi_v^meas = T_out / T_in of a marked probe
    (register R142; author's semantics 5). The probe adds delta = delta_rel
    * mean(s) at EVERY event (so every event has T_in > 0); the joint
    routing (background + probe) fixes the acceptance ratios, the marked
    component is routed with them. NOT computed from a saturation formula.
    Returns Pi^meas per event, the joint routing and the marked routing.
    """
    s = np.asarray(s, float)
    delta = delta_rel * max(float(np.mean(s)), 1e-300)
    m = np.full(s.size, delta)
    joint = route_capacity(ptr, idx, s + m, C, absorb)
    mk = route_marked(ptr, idx, joint["alpha"], joint["forwards"], m)
    pi = mk["T_out"] / mk["T_in"]
    return {"pi": pi, "joint": joint, "marked": mk, "delta": delta}


def probe_check(ptr, idx, s, C, absorb, delta_rel: float, size_factor: float) -> Dict[str, float]:
    """
    Probe instrument (card): (i) marked balance sum m = sum stored + escaped
    to round-off; (ii) size sensitivity max_v |Pi(delta) - Pi(size_factor *
    delta)|; (iii) the probe changes the background little: max_v |alpha
    without probe - Pi(delta)|.
    """
    p1 = probe_permeability(ptr, idx, s, C, absorb, delta_rel)
    p2 = probe_permeability(ptr, idx, s, C, absorb, delta_rel * size_factor)
    base = route_capacity(ptr, idx, s, C, absorb)
    mk = p1["marked"]
    tot = p1["delta"] * len(s)
    bal = abs(tot - float(mk["stored"].sum()) - mk["escaped"]) / tot
    return {"marked_balance": bal,
            "size_sensitivity": float(np.max(np.abs(p1["pi"] - p2["pi"]))),
            "background_shift": float(np.max(np.abs(base["alpha"] - p1["pi"]))),
            "pi": p1["pi"]}


def _descendant_bits(ptr: np.ndarray, idx: np.ndarray, N: int) -> np.ndarray:
    """Descendant sets as packed bit rows (N x ceil(N/8) uint8), reverse time order."""
    cptr, cidx = children_csr(ptr, idx, N)
    nb = (N + 7) // 8
    bits = np.zeros((N, nb), dtype=np.uint8)
    own = np.zeros((N, nb), dtype=np.uint8)
    ar = np.arange(N)
    own[ar, ar >> 3] = (1 << (ar & 7)).astype(np.uint8)
    for v in range(N - 1, -1, -1):
        ch = cidx[cptr[v]:cptr[v + 1]]
        if ch.size:
            bits[v] = np.bitwise_or.reduce(bits[ch] | own[ch], axis=0)
    return bits


def redundant_edges(ptr: np.ndarray, idx: np.ndarray, N: int) -> np.ndarray:
    """
    Boolean per link (aligned with idx): True if the edge u -> v is implied
    by a longer path (not a covering relation). u -> v is redundant iff v is
    a descendant of another child of u.
    """
    cptr, cidx = children_csr(ptr, idx, N)
    bits = _descendant_bits(ptr, idx, N)
    nb = bits.shape[1]
    # R_u = OR of descendants of the children of u (a child itself is not its own descendant)
    R = np.zeros((N, nb), dtype=np.uint8)
    for u in range(N):
        ch = cidx[cptr[u]:cptr[u + 1]]
        if ch.size:
            R[u] = np.bitwise_or.reduce(bits[ch], axis=0)
    v_of = np.repeat(np.arange(N), np.diff(ptr))
    u_of = idx
    return ((R[u_of, v_of >> 3] >> (v_of & 7)) & 1).astype(bool)


def rewire_links(ptr: np.ndarray, idx: np.ndarray, t: np.ndarray, swaps_per_edge: float,
                 rng: np.random.Generator, reduce: bool = True,
                 full_check: bool = True) -> Dict[str, Any]:
    """
    Fast topology control (module notes): time-forward, degree-preserving
    swaps WITHOUT the covering check during the swaps. reduce=True then
    removes the edges that are not covering relations (transitive
    reduction; this changes the degrees, measured: 75-78 % of the edges
    removed in 2+1, L1 change of the degree distribution ~1.7 of 2).
    reduce=False keeps them: a DAG with exactly the original in/out
    degrees and time order, whose edges are in general NOT covering
    relations of its order (the "DAG" variant of the null; the redundant
    fraction is reported). full_check=False (only with reduce=False) skips
    the exact N x N bitset check, which needs ~N^2/4 bytes; use
    redundant_sample instead. Returns the new parents CSR and diagnostics.
    """
    N = len(ptr) - 1
    v_of = np.repeat(np.arange(N), np.diff(ptr)).astype(np.int64)
    u_of = idx.astype(np.int64).copy()
    E = u_of.size
    key = set((u_of * N + v_of).tolist())
    accepted = 0
    n_try = int(swaps_per_edge * E)
    picks = rng.integers(0, E, size=(n_try, 2))
    for i, j in picks:
        if i == j:
            continue
        a, b, c, d = u_of[i], v_of[i], u_of[j], v_of[j]
        if a == c or b == d:
            continue
        if not (t[a] < t[d] and t[c] < t[b]):
            continue
        k1, k2 = a * N + d, c * N + b
        if k1 in key or k2 in key:
            continue
        key.discard(a * N + b)
        key.discard(c * N + d)
        key.add(k1)
        key.add(k2)
        v_of[i], v_of[j] = d, b
        accepted += 1
    # rebuild parents CSR
    order = np.lexsort((u_of, v_of))
    u_s, v_s = u_of[order], v_of[order]
    ptr2 = np.zeros(N + 1, dtype=np.int64)
    np.add.at(ptr2, v_s + 1, 1)
    ptr2 = np.cumsum(ptr2)
    if full_check or reduce:
        red = redundant_edges(ptr2, u_s, N)
    else:
        red = None
    keep = ~red if reduce else np.ones(u_s.size, bool)
    v_k, u_k = v_s[keep], u_s[keep]
    ptr3 = np.zeros(N + 1, dtype=np.int64)
    np.add.at(ptr3, v_k + 1, 1)
    ptr3 = np.cumsum(ptr3)

    def l1(a1, a2):
        m = max(a1.max(initial=0), a2.max(initial=0)) + 1
        h1 = np.bincount(a1, minlength=m) / max(a1.size, 1)
        h2 = np.bincount(a2, minlength=m) / max(a2.size, 1)
        return float(np.abs(h1 - h2).sum())

    din0, dout0 = np.diff(ptr), np.bincount(idx, minlength=N)
    din1, dout1 = np.diff(ptr3), np.bincount(u_k, minlength=N)
    return {"ptr": ptr3, "idx": u_k, "accepted_swaps": accepted, "tried": n_try,
            "redundant_fraction": (float(red.mean()) if red.size else 0.0) if red is not None else None,
            "reduced": reduce,
            "original_edges_kept": len(set((idx.astype(np.int64) * N + np.repeat(np.arange(N), np.diff(ptr))).tolist())
                                       & set((u_k.astype(np.int64) * N + v_k).tolist())) / max(idx.size, 1),
            "deg_in_L1": l1(din0, din1), "deg_out_L1": l1(dout0, dout1)}


def rewire_strict(ptr: np.ndarray, idx: np.ndarray, t: np.ndarray, n_trials: int,
                  rng: np.random.Generator, X: Optional[np.ndarray] = None,
                  local: Optional[float] = None) -> Dict[str, Any]:
    """
    Topology null in its strict form (author's semantics 6): double-edge
    swaps (a -> b), (c -> d) -> (a -> d), (c -> b) accepted only if, in the
    NEW graph,
      * every edge still points forward in time (t_a < t_d, t_c < t_b):
        no cycle, the time order is kept;
      * no duplicate edge, a != c, b != d;
      * the time span of each new edge is <= span_max, the largest span of
        the original links (keeps the check below local; reported);
      * every edge is a covering relation: neither new edge is implied by
        a longer path, and no existing edge becomes implied by a path
        through a new edge.
    In- and out-degrees are preserved exactly. Removing edges cannot make
    an edge redundant, so only paths through the new edges are checked:
    x0 -> y0 is redundant iff y0 is reachable from another child of x0
    (all intermediate times lie in (t_x0, t_y0)); an existing edge x -> y
    becomes redundant iff x is an ancestor-or-self of x0 and y a
    descendant-or-self of y0, (x, y) != (x0, y0); since t_y - t_x <=
    span_max, only ancestors with t >= t_y0 - span_max and descendants
    with t <= t_x0 + span_max are visited. Exactness of the result is
    verified afterwards by redundant_edges (must be 0).
    Pure Python; the cost per trial grows with the size of those time
    windows, so this null is practical only for N of a few thousand.

    Proposals: uniform pairs of edges (local=None), or LOCAL pairs
    (local = a length, X = coordinates): the second edge (c -> d) is drawn
    among the edges whose head d lies within `local` of the head b of the
    first edge (Euclidean distance of the embedding coordinates, time
    included). The acceptance rule is the same; locality only changes
    which swaps are proposed.
    """
    N = len(ptr) - 1
    t = np.asarray(t, float)
    v_of = np.repeat(np.arange(N), np.diff(ptr)).astype(np.int64)
    u_of = idx.astype(np.int64)
    span_max = float(np.max(t[v_of] - t[u_of])) if u_of.size else 0.0
    ch = [set() for _ in range(N)]
    pa = [set() for _ in range(N)]
    E = list(zip(u_of.tolist(), v_of.tolist()))
    for u, v in E:
        ch[u].add(v)
        pa[v].add(u)
    tl = t.tolist()

    def reach(start, nbr, ok):
        seen = {start}
        st = [start]
        while st:
            x = st.pop()
            for y in nbr[x]:
                if y not in seen and ok(y):
                    seen.add(y)
                    st.append(y)
        return seen

    def edge_ok(x0, y0):
        ty0, tx0 = tl[y0], tl[x0]
        # (1) x0 -> y0 implied by a longer path?
        st = [c for c in ch[x0] if c != y0 and tl[c] < ty0]
        seen = set(st)
        while st:
            x = st.pop()
            for y in ch[x]:
                if y == y0:
                    return False
                if y not in seen and tl[y] < ty0:
                    seen.add(y)
                    st.append(y)
        # (2) an existing edge implied through x0 -> y0?
        Y = reach(y0, ch, lambda z: tl[z] <= tx0 + span_max)
        X = reach(x0, pa, lambda z: tl[z] >= ty0 - span_max)
        for x in X:
            for y in ch[x]:
                if y in Y and not (x == x0 and y == y0):
                    return False
        return True

    accepted = time_valid = 0
    if local is None:
        picks = rng.integers(0, len(E), size=(n_trials, 2))
    else:
        from scipy.spatial import cKDTree
        tree = cKDTree(np.asarray(X, float))
        heads = np.asarray(v_of)
        edges_by_head = [[] for _ in range(N)]
        for k_, (_, v_) in enumerate(E):
            edges_by_head[v_].append(k_)
        first = rng.integers(0, len(E), size=n_trials)
        picks = []
        for i in first:
            near = tree.query_ball_point(X[heads[i]], local)
            cand = [k_ for w_ in near for k_ in edges_by_head[w_]]
            picks.append((i, cand[int(rng.integers(0, len(cand)))]))
    for i, j in picks:
        if i == j:
            continue
        (a, b), (c, d) = E[i], E[j]
        if a == c or b == d or not (tl[a] < tl[d] and tl[c] < tl[b]):
            continue
        if tl[d] - tl[a] > span_max or tl[b] - tl[c] > span_max:
            continue
        if d in ch[a] or b in ch[c]:
            continue
        time_valid += 1
        ch[a].discard(b); pa[b].discard(a); ch[c].discard(d); pa[d].discard(c)
        ch[a].add(d); pa[d].add(a); ch[c].add(b); pa[b].add(c)
        if edge_ok(a, d) and edge_ok(c, b):
            E[i], E[j] = (a, d), (c, b)
            accepted += 1
        else:
            ch[a].discard(d); pa[d].discard(a); ch[c].discard(b); pa[b].discard(c)
            ch[a].add(b); pa[b].add(a); ch[c].add(d); pa[d].add(c)
    uu = np.array([e[0] for e in E], dtype=np.int64)
    vv = np.array([e[1] for e in E], dtype=np.int64)
    order = np.lexsort((uu, vv))
    uu, vv = uu[order], vv[order]
    ptr2 = np.zeros(N + 1, dtype=np.int64)
    np.add.at(ptr2, vv + 1, 1)
    ptr2 = np.cumsum(ptr2)
    orig = set((u_of * N + v_of).tolist())
    kept = len(orig & set((uu * N + vv).tolist())) / max(len(orig), 1)
    return {"ptr": ptr2, "idx": uu, "accepted_swaps": accepted, "tried": int(n_trials),
            "time_valid": time_valid, "span_max": span_max, "original_edges_kept": kept,
            "redundant_after": int(redundant_edges(ptr2, uu, N).sum()),
            "deg_in_same": bool(np.array_equal(np.diff(ptr2), np.diff(ptr))),
            "deg_out_same": bool(np.array_equal(np.bincount(uu, minlength=N), np.bincount(idx, minlength=N)))}


def redundant_sample(ptr: np.ndarray, idx: np.ndarray, t: np.ndarray, n_sample: int,
                     rng: np.random.Generator) -> Dict[str, Any]:
    """
    Memory-light estimate for large graphs: for n_sample random edges u -> v,
    is v reachable from another child of u (intermediate times in (t_u, t_v))?
    Exact per sampled edge (time-bounded search). Returns the redundant
    fraction and the time spans of the redundant and non-redundant edges.
    """
    N = len(ptr) - 1
    cptr, cidx = children_csr(ptr, idx, N)
    v_of = np.repeat(np.arange(N), np.diff(ptr))
    E = idx.size
    pick = rng.choice(E, size=min(n_sample, E), replace=False)
    tl = np.asarray(t, float)
    red = np.zeros(pick.size, bool)
    for k, e in enumerate(pick):
        u, v = int(idx[e]), int(v_of[e])
        tv = tl[v]
        st = [int(c) for c in cidx[cptr[u]:cptr[u + 1]] if c != v and tl[c] < tv]
        seen = set(st)
        found = False
        while st and not found:
            x = st.pop()
            for y in cidx[cptr[x]:cptr[x + 1]]:
                y = int(y)
                if y == v:
                    found = True
                    break
                if y not in seen and tl[y] < tv:
                    seen.add(y)
                    st.append(y)
        red[k] = found
    span = tl[v_of[pick]] - tl[idx[pick]]
    return {"fraction": float(red.mean()), "n": int(pick.size),
            "span_redundant": span[red], "span_covering": span[~red]}


def is_time_dag(ptr: np.ndarray, idx: np.ndarray, t: np.ndarray, N: int) -> Dict[str, bool]:
    """Every edge points forward in time (hence a DAG) and no edge is duplicated."""
    v_of = np.repeat(np.arange(N), np.diff(ptr)).astype(np.int64)
    fwd = bool(np.all(np.asarray(t)[idx] < np.asarray(t)[v_of]))
    key = idx.astype(np.int64) * N + v_of
    return {"forward_in_time": fwd, "no_duplicates": bool(np.unique(key).size == key.size)}
