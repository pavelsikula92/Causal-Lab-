"""
F5.1g -- Inertial coherence of a distributed causal flow: does a flow whose events take their memory from the common
flow-past of their flow parents keep its velocity and hold together in a flat 3+1 order, without depending on the
region in which the order is read?
====================================================================================================================

Why this test exists (9 Nov 2026; card docs/F5_1g_card_draft.md, revision 2, approved by the author)
---------------------------------------------------------------------------------------------------
The fall rule of F5.1c/d/e is a proxy (it reads the future and depends on the wall). A flow of single events shared
equally among the children with no deficit against one memory event disperses and depends on the region (F5.1f-a,
INSTRUMENT_FAILED, descriptive M5 identity 0.10): one link carries no proper time, a single chain cannot hold a
direction. The author's decisions: the free flow first and gravity later; the memory is distributed (every event
reconstructs it from its local causal past, no global state of the bundle); the memory of an event is the common past
of all the parents that bring it flow, taken along the flow (the past in the whole order has an infinite maximal
boundary along the past light cones and would violate M5 by construction).

The rule (only the order, counts and the flow)
----------------------------------------------
  flow parents P(v) = {p : Phi_pv > 0}; flow-past J_Phi(v) = v and every event from which flow reached v;
  common flow-past C(v) = the intersection of J_Phi(p) over p in P(v); memory M(v) = its maximal elements (x with no
  y in C(v), x < y). With a single parent, M(v) = {p}.
  deficit of a child c of v: D(c) = sum over x in M(v) of [L(x, c) - L(x, v) - 1] >= 0;
  admissible children A(v) = {c : D(c) = 0} -- zero deficit is a condition of continuing; if A(v) is empty the flow
  stays in v (Delta Q, recorded);
  the inflow of v is shared equally among A(v); streams that meet add. The rule is deterministic.
Launch (the initial condition, not part of the rule): a source tube of radius r_s around the world line
x = x_0 + v (t - t_0) d for t in [t_0, t_0 + T_L]; every tube event receives s = 1 (no parent is created); the
common past of a tube event also contains the tube events in its causal past (the history of the hand); an event of
the first layer of the tube (nothing in its past) has an empty memory and keeps its inflow. After T_L the flow is free.

Branches (same sprinklings and sources)
---------------------------------------
  G0  memory = the flow parents only (depth one link): negative control (the dispersal of F5.1f);
  G1  the carried direction of F3.6 (coordinates; all flux to the one child of smallest rapidity relative to the
      carried u, u = flux-weighted sum of the parents' u, the source adds s u_src): positive control of the readouts,
      not a candidate rule (M3);
  G2  M(v) = maximal elements of the common flow-past (the hypothesis); computed in region B and, for M5, in region A;
  G3  G2 with the memory borrowed from another flow event of the same time slab (chosen by a hash; elements not in the
      causal past of v dropped): does the causal content of the memory matter, or only the presence of a memory?

Computation (not part of the rule)
----------------------------------
L(x, .) by the relaxation of LocalLinks inside the future cone of x within the box of the region, up to a time
window above the current event (cached per x and extended when needed); a child counts only if it lies in the box and
within t(v) + dt_w, and only if its diamond I(x, c) lies inside the box for every x in M(v) (then L is exact). In this
flat order the relation x < y equals the Minkowski relation of the sprinkled points (the order is the transitive
closure of its links), which is used to test maximality and the causal past of tube events. An event is boundary-
affected if a child of it in the time window was dropped by the box (outside it, or a diamond leaving it); a flow
event is tainted if it or an event of its flow-past is boundary-affected.

Regions: B is the box of the order with the window dt_w = 6 (F5.1-diag4: the admissible sets saturate there); A is a
smaller box with dt_A = 4.5. A is boundary-affected only through its box; a difference caused by its shorter time window
is not marked, so it shows up as a mismatch on the interior (M5a then tests the window as well). M5a counts only if the
interior holds at least min_interior = 10 % of the flow events of B; otherwise it cannot be tested and the verdict
cannot be PASSED.

Readouts (flow-weighted, after the release t_r = t_0 + T_L, per slab of one unit of time)
-----------------------------------------------------------------------------------------
centre along the throw, s in [0, S]; Delta = slope of the centre - v; width = RMS across the throw, growth = slope of
the width over s in [S/2, S]; along-RMS; held flow; memory diagnostics (|M|, age = mean L(x, v), number of parents,
eta_common = |C| / |union of the parents' flow-pasts|, share of events with D_min > 0, admissible children).

Criteria (frozen with the approval of this card; SE by bootstrap over sprinklings; checked in this order)
-------------------------------------------------------------------------------------------------------
  INSTRUMENT_FAILED  G1: |Delta(G1, v)| + 2 SE >= 0.05 for v = 0 or 0.3, or growth(G1) + 2 SE >= 0.2.
  FAILED: M5         M5a: on the events untainted in both regions the inflow differs (relative 1e-9) in any
                     sprinkling, or an event untainted in one region is missing in the other (PASSED further needs the
                     interior to hold >= 10 % of the flow events of B); M5b: |Delta_A - Delta_B|
                     + 2 SE >= 0.05 or |w_A - w_B| / w_B + 2 SE >= 0.1 for v = 0 or 0.3 (paired, G2).
  FAILED: inertia    Delta(G2, 0.3) + 2 SE < -0.05, or |Delta(G2, 0)| - 2 SE > 0.05.
  FAILED: coherence  growth(G2) - 2 SE >= 0.2, or growth(G0) - growth(G2) - 2 SE <= 0 (paired) for both v.
  FAILED: memory     growth(G3) - growth(G2) - 2 SE <= 0 (paired) for both v -- this hypothesis about the causal
          content    content of the memory fails, not memory as a mechanism.
  PASSED             M5a and M5b hold, |Delta(G2, v)| + 2 SE < 0.05 for both v, growth(G2) + 2 SE < 0.2 and G2 more
                     coherent than G0 and than G3 by 2 SE (paired) for both v.
  UNDECIDED          otherwise.
The numbers of the launch and of the regions (r_s, T_L, S, box sizes, dt_w) are those of the params below; the card
allowed fixing them after development rows in 2+1; the author chose to run the test as it is (9 Nov), so they are
frozen here.

Reading declared now
--------------------
PASSED: a flow whose memory is its own common past, read by the order alone, keeps its velocity and holds together; it
is the candidate for F5.1h (the same rule in a homogeneous gradient of the order). FAILED: M5: the common flow-past
still lets the flow depend on the region. FAILED: inertia or coherence: the common flow-past is too old or too narrow
to hold a direction (the memory diagnostics say which). FAILED: memory content: what holds the flow is the presence of
a memory, not its causal content. Even PASSED gives a candidate mechanism of inertia, not gravity.
"""

from __future__ import annotations

import heapq
import math
import time

import numpy as np

from core.order.sprinkle import derive_seed
from core.stats.parallel import parallel_map
from core.stats.stats import boot_se
from core.testing import Result
from phases.F5_gravity.F5_window import Window, diamond_margin, id_hash, _SWEEP
from phases.F7_light.F7_01c_throughput_cones_3d import CapacityOrder3

TEST = {
    "id": "F5.1g",
    "phase": "F5",
    "title": "Inertial coherence of a distributed causal flow (F5.1g): does a flow whose memory is the common flow-past "
             "of its parents keep its velocity and hold together in a flat 3+1 order?",
    "theory_part": "F5 gravity and equivalence; ontology (matter = influence that persists; memory = inertia; "
                   "superposition = coexisting admissible continuations), Axiom of continuity, Axiom of the action, "
                   "M2-M5",
    "question": "Does a flow whose events take their memory from the maximal elements of the common flow-past of their "
                "flow parents, and pass their flow equally to the children with zero deficit against every memory "
                "event, keep the velocity of its centre and its coherence in a flat 3+1 order, without depending on the "
                "region in which the order is read?",
    "method": "flat 3+1 order, density 1; launch: a source tube of radius r_s moving with v = 0 or 0.3 along +-x for T_L, "
              "then free for S; branches G0 (memory = parents), G1 (carried direction of F3.6, positive control), G2 "
              "(common flow-past, the hypothesis; regions B and A), G3 (G2 with memories borrowed within a time slab); "
              "flow-weighted centre and width per slab; 8 sprinklings, 2 sources per v and sprinkling",
    "prediction": "G2 keeps the velocity of its centre (Delta = 0), grows in width slower than 0.2 and slower than G0 "
                  "and G3, and its flow does not depend on the region",
    "criterion": "INSTRUMENT_FAILED if G1 fails |Delta| + 2 SE < 0.05 (both v) or growth + 2 SE < 0.2. FAILED: M5 if "
                 "M5a (identical inflow on untainted events) or M5b (|Delta_A - Delta_B| + 2 SE < 0.05, |w_A - w_B| / w_B "
                 "+ 2 SE < 0.1) fails. FAILED: inertia if Delta(G2, 0.3) + 2 SE < -0.05 or |Delta(G2, 0)| - 2 SE > 0.05. "
                 "FAILED: coherence if growth(G2) - 2 SE >= 0.2 or G2 not narrower-growing than G0 by 2 SE (both v). "
                 "FAILED: memory content if G3 not wider-growing than G2 by 2 SE (both v). PASSED if all hold with "
                 "|Delta(G2, v)| + 2 SE < 0.05 and growth(G2) + 2 SE < 0.2. UNDECIDED otherwise.",
    "null_model": "G0 (memory of depth one link) and G3 (memory without its causal content)",
    "level": "mechanism",
    "meters": [],
    "depends_on": ["F5.1f", "F5.1-diag4", "F3.6"],
    "params": {"T": 24.0, "half_x": 9.0, "half_yz": 7.0, "half_x_A": 7.5, "half_yz_A": 5.5, "dt_w": 6.0, "dt_A": 4.5,
               "min_interior": 0.1,
               "width": 4.0, "r_s": 1.5, "t0": 1.0, "T_L": 5.0, "S": 10.0, "slab": 1.0, "y_off": 1.5, "n_seeds": 8,
               "n_se": 2.0, "tol": 0.05, "tol_w": 0.1, "growth_max": 0.2, "rel_ident": 1e-9, "boot_B": 4000,
               "max_flow_events": 60000},
    "parts": {"main": {}},
    "quick": {"n_seeds": 2, "boot_B": 200, "T": 12.0, "half_x": 5.0, "half_yz": 4.0, "half_x_A": 4.0,
              "half_yz_A": 3.5, "dt_w": 2.5, "dt_A": 2.0, "r_s": 1.0, "T_L": 2.5, "S": 4.0, "y_off": 1.0,
              "max_flow_events": 4000},
    "estimate": "unknown before the first full run; per sprinkling 4 sources x 5 flows (G0, G1, G2 in B and A, G3) of "
                "up to max_flow_events events; ~1-3 GB per process",
}

V_THROW = (0.0, 0.3)


# ============================================================================================== longest chains

class LCache:
    """L(x, .) for memory events x: longest chains from x inside the future cone of x within the box, up to a top
    time; extended (recomputed with a higher top) when an event above the top is asked for."""

    def __init__(self, order, half, dt_w):
        self.order, self.half, self.dt_w = order, half, float(dt_w)
        self.data = {}

    def ensure(self, x, t_need):
        if x in self.data and self.data[x][2] >= t_need:
            return
        X = self.order.X
        Tt = X[:, 0]
        top = t_need + 2.0
        tx, xx = float(X[x, 0]), X[x, 1:]
        i0 = int(np.searchsorted(Tt, tx, "left"))
        i1 = int(np.searchsorted(Tt, top, "right"))
        idx = np.arange(i0, i1)
        P = X[idx]
        keep = (np.sqrt(((P[:, 1:] - xx) ** 2).sum(axis=1)) <= (P[:, 0] - tx) + 1e-12) & \
               np.all(np.abs(P[:, 1:]) <= self.half[None, :], axis=1)
        win = Window(self.order, idx[keep])
        jx = win.loc(x)
        Lx = _SWEEP(win.n, jx, win.f_u, win.f_w, win.f_u >= jx) if jx >= 0 else np.full(win.n, -1)
        self.data[x] = (win.ids, Lx, top)

    def L(self, x, ev):
        """L(x, e) for the global events ev (-1 if e is not in the computed future of x)."""
        ids, Lx, _ = self.data[x]
        ev = np.atleast_1d(np.asarray(ev, dtype=np.int64))
        j = np.minimum(np.searchsorted(ids, ev), max(ids.size - 1, 0))
        ok = (ids.size > 0) & (ids[j] == ev) if ids.size else np.zeros(ev.size, bool)
        out = np.full(ev.size, -1, dtype=np.int64)
        out[ok] = Lx[j[ok]]
        return out


# ============================================================================================== bitsets

def bits_of(b):
    """Indices of the set bits of a Python int."""
    if b == 0:
        return np.zeros(0, dtype=np.int64)
    raw = b.to_bytes((b.bit_length() + 7) // 8, "little")
    return np.flatnonzero(np.unpackbits(np.frombuffer(raw, dtype=np.uint8), bitorder="little")).astype(np.int64)


def precedes(X, a, b):
    """a < b in the flat order (the Minkowski relation of the sprinkled points; module notes)."""
    d = X[b] - X[a]
    return (d[..., 0] > 0) & (np.sqrt((d[..., 1:] ** 2).sum(axis=-1)) <= d[..., 0] + 1e-12)


# ============================================================================================== one flow

def run_flow(order, p, src, branch, half, rec_events, dt_w=None):
    """Propagate the flow of one source by one branch in one region. Returns per-event inflow (dict global id ->
    float), the tainted set and the diagnostics. src: dict with x0, d, v, t0, T_L, r_s, u_src."""
    X = order.X
    t_end = src["t0"] + src["T_L"] + float(p["S"]) + 1.0
    dt_w = float(p["dt_w"]) if dt_w is None else float(dt_w)
    lc = LCache(order, half, dt_w) if branch in ("G0", "G2", "G3") else None
    # the source tube: events within r_s of the moving centre during the launch
    Tt = X[:, 0]
    i0 = int(np.searchsorted(Tt, src["t0"], "left"))
    i1 = int(np.searchsorted(Tt, src["t0"] + src["T_L"], "right"))
    cand = np.arange(i0, i1)
    ctr = src["x0"][None, :] + (X[cand, 0:1] - src["t0"]) * (src["v"] * src["d"])[None, :]
    tube = cand[np.sqrt(((X[cand, 1:] - ctr) ** 2).sum(axis=1)) <= src["r_s"]]
    tube_set = set(int(e) for e in tube)
    inflow, parents, held = {}, {}, {}
    for e in tube:
        inflow[int(e)] = inflow.get(int(e), 0.0) + 1.0
    heap = [int(e) for e in tube]
    heapq.heapify(heap)
    queued = set(heap)
    kidx, J, fkids, u_of = {}, [], [], {}              # flow-event index, flow-past bitsets, flow children
    glob = []                                           # flow-event index -> global id
    aff = set()                                         # boundary-affected events
    diag = {"held": 0.0, "dmin_pos": 0, "events": 0, "M": [], "age": [], "npar": [], "eta": [], "adm": []}
    true_mem, mem_t = [], []                            # G3: the true memories of the processed events and their times
    u_src = src["u_src"]
    while heap:
        v = heapq.heappop(heap)
        if X[v, 0] > t_end:
            break
        if len(kidx) >= int(p["max_flow_events"]):
            diag["truncated"] = True
            break
        phi = inflow.get(v, 0.0)
        if phi <= 0.0:
            continue
        diag["events"] += 1
        par = parents.get(v, {})
        k = len(J)
        kidx[v] = k
        glob.append(v)
        fkids.append([])
        # ---------------------------------------------------------------- the flow-past and the common flow-past
        tube_past = 0
        if v in tube_set:
            tp = [w for w in tube_set if w != v and X[w, 0] < X[v, 0]]
            if tp:
                tp = np.array(tp)
                tp = tp[precedes(X, tp, np.full(tp.size, v))]
                for w in tp:
                    if int(w) in kidx:
                        tube_past |= 1 << kidx[int(w)]
        if par:
            C = -1
            U = 0
            for q in par:
                C &= J[kidx[q]]
                U |= J[kidx[q]]
            C |= tube_past
            U |= tube_past
        else:
            C, U = tube_past, tube_past
        J.append((1 << k) | U)
        for q in par:
            fkids[kidx[q]].append(k)
        # memory
        if branch == "G0":
            mem = list(par.keys()) if par else []
            if not par and tube_past:
                mem = [glob[int(i)] for i in maximal_local(X, C, fkids, glob)]
        elif branch in ("G2", "G3"):
            mem = [glob[int(i)] for i in maximal_local(X, C, fkids, glob)]
            if branch == "G3":
                # the memory of another flow event of the same time slab (chosen by a hash), its elements outside the
                # causal past of v dropped
                lo = int(np.searchsorted(np.array(mem_t), float(X[v, 0]) - float(p["slab"]), "left")) if mem_t else 0
                ncand = len(mem_t) - lo
                own = mem
                if ncand > 0:
                    pick = lo + int(id_hash(np.array([v], dtype=np.uint64))[0] % np.uint64(ncand))
                    mem = [x for x in true_mem[pick] if bool(precedes(X, np.array([x]), np.array([v]))[0])]
                mem_t.append(float(X[v, 0]))
                true_mem.append(own)
        else:
            mem = None
        # children in the region
        kids = order.cidx[order.cptr[v]:order.cptr[v + 1]]
        kids = kids[X[kids, 0] <= X[v, 0] + dt_w]
        in_box = np.all(np.abs(X[kids, 1:]) <= half[None, :], axis=1)
        affected = bool((~in_box).any())
        kids = kids[in_box]
        if branch == "G1":
            # carried direction: u_v = flux-weighted parents' u + s u_src, all flux to the child of least rapidity
            u = sum((w_ * u_of[q] for q, w_ in par.items()), np.zeros(4)) + (u_src if v in tube_set else 0.0)
            nrm = u[0] ** 2 - float(u[1:] @ u[1:])
            if nrm <= 0 or kids.size == 0:
                held[v] = phi
                diag["held"] += phi
                continue
            u = u / math.sqrt(nrm)
            u_of[v] = u
            b = X[kids] - X[v]
            tau = np.sqrt(np.maximum(b[:, 0] ** 2 - (b[:, 1:] ** 2).sum(axis=1), 1e-300))
            ch = (u[0] * b[:, 0] - b[:, 1:] @ u[1:]) / tau
            best = kids[ch == ch.min()]
            adm = np.array([best[np.argmin(id_hash(best.astype(np.uint64)))]])
        else:
            if affected:
                aff.add(v)
            if not mem:
                held[v] = phi
                diag["held"] += phi
                continue
            tx = float(X[v, 0]) + dt_w
            Ls = []
            for x in mem:
                lc.ensure(x, tx)
                margin = diamond_margin(X[x, 1:], float(X[x, 0]), X[kids], half) if kids.size else np.zeros(0)
                if kids.size and (margin < 0).any():
                    affected = True
                    kids = kids[margin >= 0]
                Ls.append(x)
            if affected:
                aff.add(v)
            if kids.size == 0:
                held[v] = phi
                diag["held"] += phi
                continue
            D = np.zeros(kids.size, dtype=np.int64)
            ages = []
            for x in mem:
                Lv = int(lc.L(x, [v])[0])
                Lc = lc.L(x, kids)
                D += np.where(Lc >= 0, Lc - Lv - 1, 10 ** 6)
                ages.append(Lv)
            dmin = int(D.min())
            diag["M"].append(len(mem))
            diag["age"].append(float(np.mean(ages)))
            diag["npar"].append(len(par))
            diag["eta"].append(len(bits_of(C)) / max(len(bits_of(U)), 1))
            if dmin > 0:
                diag["dmin_pos"] += 1
                held[v] = phi
                diag["held"] += phi
                continue
            adm = kids[D == 0]
            diag["adm"].append(int(adm.size))
        if affected:
            aff.add(v)
        share = phi / adm.size
        for c in adm:
            c = int(c)
            inflow[c] = inflow.get(c, 0.0) + share
            parents.setdefault(c, {})
            parents[c][v] = parents[c].get(v, 0.0) + share
            if c not in queued:
                heapq.heappush(heap, c)
                queued.add(c)
    if rec_events:
        diag["inflow"] = {int(e): float(f) for e, f in inflow.items() if e in kidx}
        diag["parents"] = {int(e): [int(q) for q in parents.get(e, {})] for e in kidx}
        diag["aff"] = sorted(int(e) for e in aff)
    return inflow, kidx, diag


def taint(parents, order_ev, bad):
    """Events whose flow-past contains an event of bad (both regions' boundary-affected events)."""
    t = set()
    for e in order_ev:                                  # time order: parents come first
        if e in bad or any(q in t for q in parents.get(e, [])):
            t.add(e)
    return t


def maximal_local(X, C, fkids, inv):
    """Maximal elements of the set C (bitset over flow-event indices): no flow child in C, then no element of C later
    in the order (the flat relation)."""
    idx = bits_of(C)
    if idx.size == 0:
        return idx
    cand = np.array([i for i in idx if not any((C >> c) & 1 for c in fkids[i])], dtype=np.int64)
    if cand.size <= 1:
        return cand
    g = np.array([inv[int(i)] for i in cand])
    A = precedes(X, g[:, None], g[None, :])             # A[a, b]: cand a < cand b
    return cand[~A.any(axis=1)]


# ============================================================================================== readouts

def readout(X, inflow, events, src, p):
    """Flow-weighted centre along the throw, RMS across and along, per slab after the release."""
    t_r = src["t0"] + src["T_L"]
    S, slab = float(p["S"]), float(p["slab"])
    ev = np.array([e for e in events if X[e, 0] >= t_r - 0.5 * slab], dtype=np.int64)
    if ev.size == 0:
        return None
    w = np.array([inflow[int(e)] for e in ev])
    s = X[ev, 0] - t_r
    rel = X[ev, 1:] - src["x0"][None, :]
    along = rel @ src["d"]
    across = rel - along[:, None] * src["d"][None, :]
    out = {"s": [], "c": [], "wa": [], "wl": [], "mass": []}
    for k in range(int(S) + 1):
        m_ = np.abs(s - k * slab) < 0.5 * slab
        if w[m_].sum() <= 0:
            continue
        ww = w[m_] / w[m_].sum()
        c = float((ww * along[m_]).sum())
        ca = (ww[:, None] * across[m_]).sum(axis=0)
        out["s"].append(k * slab)
        out["c"].append(c)
        out["wa"].append(float(math.sqrt((ww * ((across[m_] - ca) ** 2).sum(axis=1)).sum())))
        out["wl"].append(float(math.sqrt((ww * (along[m_] - c) ** 2).sum())))
        out["mass"].append(float(w[m_].sum()))
    if len(out["s"]) < 3:
        return None
    s_ = np.array(out["s"])
    out["v_centre"] = float(np.polyfit(s_, out["c"], 1)[0])
    late = s_ >= S / 2.0
    out["growth"] = float(np.polyfit(s_[late], np.array(out["wa"])[late], 1)[0]) if late.sum() >= 2 else float("nan")
    out["w_end"] = out["wa"][-1]
    return out


# ============================================================================================== one sprinkling

def seed_task(p, seed):
    t_start = time.time()
    rng = np.random.default_rng(seed)
    T, hx, hyz = float(p["T"]), float(p["half_x"]), float(p["half_yz"])
    n = rng.poisson(T * (2 * hx) * (2 * hyz) ** 2)
    X = np.column_stack([rng.uniform(0, T, n), rng.uniform(-hx, hx, n),
                         rng.uniform(-hyz, hyz, n), rng.uniform(-hyz, hyz, n)])
    X = X[np.argsort(X[:, 0], kind="stable")]
    order = CapacityOrder3(X, 0.0, float(p["width"]))
    hB = np.array([hx, hyz, hyz])
    hA = np.array([p["half_x_A"], p["half_yz_A"], p["half_yz_A"]], float)
    out = {"seed": int(seed), "N": int(len(X)), "flows": []}
    for v in V_THROW:
        for sgn in (1.0, -1.0):
            span = float(p["T_L"]) + float(p["S"])
            if v > 0.0:
                d = np.array([sgn, 0.0, 0.0])
                x0 = -0.5 * v * span * d                # the walk symmetric in the box
            else:
                # at rest the two sources must be different flows (the same flow read along +x and -x would cancel
                # exactly): both are read along +x, placed at y = +-y_off
                d = np.array([1.0, 0.0, 0.0])
                x0 = np.array([0.0, sgn * float(p["y_off"]), 0.0])
            g = 1.0 / math.sqrt(1.0 - v * v)
            src = {"x0": x0, "d": d, "v": v, "t0": float(p["t0"]), "T_L": float(p["T_L"]), "r_s": float(p["r_s"]),
                   "u_src": np.concatenate([[g], g * v * d])}
            rowf = {"v": v, "sgn": sgn}
            for br, half, tag in (("G0", hB, "G0"), ("G1", hB, "G1"), ("G2", hB, "G2"), ("G2", hA, "G2A"),
                                  ("G3", hB, "G3")):
                inflow, kidx, dg = run_flow(order, p, src, br, half, rec_events=tag in ("G2", "G2A"),
                                            dt_w=p["dt_A"] if tag == "G2A" else p["dt_w"])
                ro = readout(X, inflow, list(kidx.keys()), src, p)
                rowf[tag] = {"ro": ro, "events": dg["events"], "held": dg["held"], "dmin_pos": dg["dmin_pos"],
                             "truncated": dg.get("truncated", False),
                             "M": float(np.mean(dg["M"])) if dg["M"] else float("nan"),
                             "age": float(np.mean(dg["age"])) if dg["age"] else float("nan"),
                             "npar": float(np.mean(dg["npar"])) if dg["npar"] else float("nan"),
                             "eta": float(np.mean(dg["eta"])) if dg["eta"] else float("nan"),
                             "adm": float(np.mean(dg["adm"])) if dg["adm"] else float("nan")}
                if tag in ("G2", "G2A"):
                    rowf[tag]["inflow"] = dg["inflow"]
                    rowf[tag]["parents"] = dg["parents"]
                    rowf[tag]["aff"] = dg["aff"]
            # M5a: identical inflow on the events untainted in both regions (taint: the flow-past contains an event
            # that is boundary-affected in either region)
            fB, fA = rowf["G2"].pop("inflow"), rowf["G2A"].pop("inflow")
            pB, pA = rowf["G2"].pop("parents"), rowf["G2A"].pop("parents")
            bad = set(rowf["G2"].pop("aff")) | set(rowf["G2A"].pop("aff"))
            tB, tA = taint(pB, sorted(pB), bad), taint(pA, sorted(pA), bad)
            mism, n_int = 0, 0
            for e in set(fB) | set(fA):
                if e in tB or e in tA:
                    continue
                if e in fB and e in fA:
                    n_int += 1
                    if abs(fB[e] - fA[e]) > float(p["rel_ident"]) * max(abs(fB[e]), 1e-300):
                        mism += 1
                else:
                    mism += 1
            rowf["m5a"] = {"interior": n_int, "mismatch": mism, "n_B": len(fB)}
            out["flows"].append(rowf)
    out["seconds"] = time.time() - t_start
    return out


# ============================================================================================== the test

def _nm(x):
    """Mean ignoring nan; nan for an empty or all-nan list."""
    x = np.asarray(x, float)
    return float(np.nanmean(x)) if np.isfinite(x).any() else float("nan")


def run(p, ctx):
    t0 = time.time()
    n, B, nse = int(p["n_seeds"]), int(p["boot_B"]), float(p["n_se"])
    tol, tol_w, gmax = float(p["tol"]), float(p["tol_w"]), float(p["growth_max"])
    jobs = [(p, derive_seed(ctx.seed, "F5.1g", ctx.part, k)) for k in range(n)]
    rows = parallel_map(seed_task, jobs, ctx.workers, ctx.progress,
                        lambda i, r_, dn, nn: f"[F5.1g] {dn}/{nn}: N {r_['N']}, {len(r_['flows'])} flows "
                                              f"({time.time() - t0:.0f}s)")

    def val(f, tag, key):
        ro = f[tag]["ro"]
        if ro is None:
            return float("nan")
        if key == "delta":
            return ro["v_centre"] - f["v"]
        return ro[key]

    def stat(tag, key, v, other=None, rel=False):
        per = []
        for r_ in rows:
            vals = []
            for f in r_["flows"]:
                if f["v"] != v:
                    continue
                a = val(f, tag, key)
                if other is not None:
                    b = val(f, other, key)
                    a = (a - b) / b if rel else a - b
                vals.append(a)
            vals = [x for x in vals if np.isfinite(x)]
            if vals:
                per.append(float(np.mean(vals)))
        if len(per) < 2:
            return (float(np.mean(per)) if per else float("nan"), float("nan"))
        return tuple(boot_se(per, B, derive_seed(ctx.seed, "F5.1g-boot", tag, key, v, str(other)))[:2])

    lines = [f"flat 3+1 order, density 1, box t [0, {p['T']:g}], |x| <= {p['half_x']:g}, |y|, |z| <= {p['half_yz']:g} "
             f"(~{np.mean([r_['N'] for r_ in rows]):.0f} events); source r_s {p['r_s']:g}, launch {p['T_L']:g}, free "
             f"{p['S']:g}; window dt_w {p['dt_w']:g}; region A |x| <= {p['half_x_A']:g}, |y|, |z| <= {p['half_yz_A']:g}"]
    S_ = {}
    for v in V_THROW:
        for tag in ("G0", "G1", "G2", "G2A", "G3"):
            S_[(tag, v)] = {"d": stat(tag, "delta", v), "g": stat(tag, "growth", v), "w": stat(tag, "w_end", v)}
            F_ = [f[tag] for r_ in rows for f in r_["flows"] if f["v"] == v]
            lines.append(f"  v {v:.1f} {tag:4s}: Delta {S_[(tag, v)]['d'][0]:+.3f} +- {S_[(tag, v)]['d'][1]:.3f}, "
                         f"growth {S_[(tag, v)]['g'][0]:+.3f} +- {S_[(tag, v)]['g'][1]:.3f}, width at the end "
                         f"{S_[(tag, v)]['w'][0]:.2f}; events {np.mean([f_['events'] for f_ in F_]):.0f}, held "
                         f"{np.mean([f_['held'] for f_ in F_]):.1f}, D_min>0 {np.mean([f_['dmin_pos'] for f_ in F_]):.1f}, "
                         f"|M| {_nm([f_['M'] for f_ in F_]):.1f}, age {_nm([f_['age'] for f_ in F_]):.1f}, "
                         f"parents {_nm([f_['npar'] for f_ in F_]):.1f}, eta {_nm([f_['eta'] for f_ in F_]):.2f}, "
                         f"admissible {_nm([f_['adm'] for f_ in F_]):.1f}, truncated "
                         f"{np.mean([f_['truncated'] for f_ in F_]):.2f}")
    m5a = [f["m5a"] for r_ in rows for f in r_["flows"]]
    m5a_ident = all(x["mismatch"] == 0 for x in m5a)
    frac_int = sum(x["interior"] for x in m5a) / max(sum(x["n_B"] for x in m5a), 1)
    m5a_testable = frac_int >= float(p["min_interior"])
    m5a_ok = m5a_ident and m5a_testable
    lines.append(f"  M5a: events untainted in both regions {sum(x['interior'] for x in m5a)} ({frac_int:.2f} of the flow "
                 f"events of B; needed >= {p['min_interior']:g}), mismatches {sum(x['mismatch'] for x in m5a)} -> "
                 f"{'identical' if m5a_ident else 'DIFFERENT'}{'' if m5a_testable else ' (too few to test)'}")
    m5b_ok = True
    for v in V_THROW:
        dd = stat("G2", "delta", v, "G2A")
        dw = stat("G2A", "w_end", v, "G2", rel=True)
        ok = np.isfinite(dd[1]) and abs(dd[0]) + nse * dd[1] < tol and np.isfinite(dw[1]) and \
            abs(dw[0]) + nse * dw[1] < tol_w
        m5b_ok &= bool(ok)
        lines.append(f"  M5b v {v:.1f}: Delta_B - Delta_A {dd[0]:+.3f} +- {dd[1]:.3f}; (w_A - w_B) / w_B {dw[0]:+.3f} "
                     f"+- {dw[1]:.3f} -> {'ok' if ok else 'violated'}")
    inst = all(np.isfinite(S_[("G1", v)]["d"][1]) and abs(S_[("G1", v)]["d"][0]) + nse * S_[("G1", v)]["d"][1] < tol
               and S_[("G1", v)]["g"][0] + nse * S_[("G1", v)]["g"][1] < gmax for v in V_THROW)
    d0, d3 = S_[("G2", 0.0)]["d"], S_[("G2", 0.3)]["d"]
    lost = d3[0] + nse * d3[1] < -tol or abs(d0[0]) - nse * d0[1] > tol
    kept = all(abs(S_[("G2", v)]["d"][0]) + nse * S_[("G2", v)]["d"][1] < tol for v in V_THROW)
    pair0 = {v: stat("G0", "growth", v, "G2") for v in V_THROW}
    pair3 = {v: stat("G3", "growth", v, "G2") for v in V_THROW}
    for v in V_THROW:
        lines.append(f"  v {v:.1f}: growth G0 - G2 {pair0[v][0]:+.3f} +- {pair0[v][1]:.3f}; G3 - G2 {pair3[v][0]:+.3f} "
                     f"+- {pair3[v][1]:.3f}")
    coh_fail = any(S_[("G2", v)]["g"][0] - nse * S_[("G2", v)]["g"][1] >= gmax for v in V_THROW) or \
        all(pair0[v][0] - nse * pair0[v][1] <= 0 for v in V_THROW)
    coh_ok = all(S_[("G2", v)]["g"][0] + nse * S_[("G2", v)]["g"][1] < gmax and pair0[v][0] - nse * pair0[v][1] > 0
                 for v in V_THROW)
    mem_fail = all(pair3[v][0] - nse * pair3[v][1] <= 0 for v in V_THROW)
    mem_ok = all(pair3[v][0] - nse * pair3[v][1] > 0 for v in V_THROW)
    lines.append(f"instrument (G1): {'OK' if inst else 'FAILS'}; M5a {m5a_ok}; M5b {m5b_ok}; kept {kept}; lost {lost}; "
                 f"coherence ok {coh_ok} / fails {coh_fail}; memory content ok {mem_ok} / fails {mem_fail}")
    if not inst:
        status = "INSTRUMENT_FAILED"
    elif not m5a_ident or not m5b_ok:
        status = "FAILED"
        lines.append("FAILED: M5")
    elif lost:
        status = "FAILED"
        lines.append("FAILED: inertia")
    elif coh_fail:
        status = "FAILED"
        lines.append("FAILED: coherence")
    elif mem_fail:
        status = "FAILED"
        lines.append("FAILED: memory content")
    elif kept and coh_ok and mem_ok and m5a_ok:
        status = "PASSED"
    else:
        status = "UNDECIDED"
    return Result(status, lines, {"rows": rows})
