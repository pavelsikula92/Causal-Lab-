"""
core.geometry.order_read
========================

Is a causal order geometric, and of which dimension? The E38-B reading of
the 2.x archive (gate G2 of the E43 certificate), ported into Causality
Lab so that it no longer needs the archive. Used by K.1 (calibration in
3+1; history E54) and later by F0 and F6 (orders grown by a dynamics).

Two code paths read the same quantities:

  * general orders (`Poset`): ancestor and descendant sets as Python
    integers used as bitsets, the archive's own representation. Exact for
    any order (sprinklings, percolation, null models, grown orders), but
    memory and time grow like N^2 / 8 bytes: used up to ~16 000 elements.
  * sprinklings with coordinates (`read_coords`): the same readings from
    the coordinates (the relation is decided on the fly), so 3+1 sizes up
    to 128 000 are possible. Ported from GTS Lab 3.x
    (gts_lab.instruments.geometry_read).
The two paths agree element by element on the same interval and the same
element (core check `order_read_paths`).

Readings (identical for every order)
------------------------------------
  interval   (E35) u uniform, w uniform among the descendants of u,
             S = I(u, w) + {u, w}, kept if 40 <= |S| <= 4000; ordering
             fraction r of S and longest chain of S. D_MM inverts the MEAN
             r (Myrheim-Meyer); D_chain = 1 / slope of ln(chain) against
             ln(|S|); D_MM of the smaller and of the larger half (scale
             drift); D_sd = spread of the per-interval dimensions
  height     longest chain of the whole order
  homog      (E38-B) D_MM of intervals grouped by the height third of u
  action     smeared Benincasa-Dowker dS_d at eps = 0.05 for d = 2, 3, 4,
             mean over sampled elements of the upper half of the linear
             extension: dS = 1 + sum_n N_n eps f_d(n, eps), N_n = ancestors
             with n elements strictly between (a sprinkling reads ~0 at
             its own dimension)

Decision (unchanged from the archive; thresholds fixed in 2.x)
--------------------------------------------------------------
  classify (E35)   MANIFOLD_LIKE(D) with D = (D_MM + D_chain) / 2 unless:
                   D_sd > 0.40, |D_MM - D_chain| > 0.60, scale drift
                   > 0.60, size drift > 0.50 (over sizes), no intervals
  decide (E38-B)   NOT_MANIFOLD_LIKE / CHAIN (d* = 1) / INHOMOGENEOUS
                   (homogeneity spread > 0.6) / ACTION_MISMATCH (|dS_d*| > 1)
                   / EXPLAINED_BY_NULL (a matched null MANIFOLD_LIKE with
                   |D_null - D| <= 0.6) / GEOMETRIC(D, d*), d* the integer
                   nearest to D

Matched null models (archive kernel.null_models, unchanged)
-----------------------------------------------------------
  tp       originary transitive percolation, same N and relation density
  seq      same number of seeds at every step (for a sprinkling: its link
           counts), chosen uniformly among earlier elements
  layered  same layer profile, same mean number of links per element
  percolation  originary transitive percolation p = 0.01 (negative control)

Port: the archive functions are reproduced statement by statement with
the same random streams, so for the same order and seed they return the
same numbers; this was checked against the archive on 19 Oct 2026 (the
archive is not part of Causality Lab; the check is recorded in K.1).
"""

from __future__ import annotations

import math
import random
import statistics
from functools import lru_cache
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

from ..order.relation import prec_matrix, prec_vector, succ_vector, chain_dp, tau2_matrix
from ..order.sprinkle import interval_volume
from ..operators.bd import BD_COEFFS, kernel_cutoff_n
from ..operators.bd import smeared_kernel as smeared_kernel_vec

ACTION_EPS = 0.05
# E35 classifier and E38-B decision thresholds (archive constants)
SD_MAX = 0.40
AGREE_MAX = 0.60
SCALE_MAX = 0.60
SIZE_DRIFT_MAX = 0.50
D_TOL = 0.60
HOMOG_MAX = SCALE_MAX
ACTION_TOL = 1.0


# =============================================================================
# Myrheim-Meyer
# =============================================================================

def mm_ordering_fraction(D: float) -> float:
    """Expected fraction of related pairs in a D-dimensional Minkowski interval: r(D) = 2 f(D),
    f(D) = Gamma(D+1) Gamma(D/2) / (4 Gamma(3D/2))."""
    if D <= 0:
        return float("nan")
    return 2.0 * math.gamma(D + 1.0) * math.gamma(D / 2.0) / (4.0 * math.gamma(1.5 * D))


def mm_dimension(r: float, *, D_lo: float = 1.0, D_hi: float = 12.0, tol: float = 1e-10) -> float:
    """Invert r(D) by bisection (r falls with D); r >= r(1) -> 1 (a chain); r <= r(12) -> nan."""
    if not (0.0 < r <= 1.0):
        return float("nan")
    if r >= mm_ordering_fraction(D_lo):
        return float(D_lo)
    if r <= mm_ordering_fraction(D_hi):
        return float("nan")
    lo, hi = D_lo, D_hi
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if mm_ordering_fraction(mid) > r:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# =============================================================================
# General orders: bitsets (archive kernel.causal_growth.Poset)
# =============================================================================

def _bits(x: int):
    """Indices of the set bits of x, ascending."""
    while x:
        low = x & -x
        yield low.bit_length() - 1
        x ^= low


def _popcount(x: int) -> int:
    return bin(x).count("1")


class Poset:
    """Elements 0..n-1 in a linear extension; parents[j] are the LINKS below j. Bitsets on demand."""

    def __init__(self, parents: List[Tuple[int, ...]], label: str = ""):
        self.parents = [tuple(p) for p in parents]
        self.n = len(parents)
        self.label = label
        self._anc: Optional[List[int]] = None
        self._desc: Optional[List[int]] = None

    def anc(self) -> List[int]:
        if self._anc is None:
            a = [0] * self.n
            for j, ps in enumerate(self.parents):
                m = 0
                for p in ps:
                    m |= a[p] | (1 << p)
                a[j] = m
            self._anc = a
        return self._anc

    def desc(self) -> List[int]:
        if self._desc is None:
            d = [0] * self.n
            for j in range(self.n - 1, -1, -1):
                for p in self.parents[j]:
                    d[p] |= d[j] | (1 << j)
            self._desc = d
        return self._desc

    def heights(self) -> List[int]:
        h = [0] * self.n
        for j, ps in enumerate(self.parents):
            if ps:
                h[j] = 1 + max(h[p] for p in ps)
        return h


def links_from_relation(anc: List[int]) -> List[Tuple[int, ...]]:
    """Transitive reduction: i is a link below j iff i is a maximal element of anc[j]."""
    out: List[Tuple[int, ...]] = []
    for j, A in enumerate(anc):
        links = []
        covered = 0
        rest = A
        while rest:
            i = rest.bit_length() - 1
            rest ^= 1 << i
            if covered >> i & 1:
                continue
            links.append(i)
            covered |= anc[i]
            rest &= ~covered
        out.append(tuple(sorted(links)))
    return out


def poset_from_points(X: np.ndarray, label: str = "") -> Poset:
    """The order of time-sorted points (row 0 first) as a Poset; same relation as the archive's sprinkled_cone."""
    T = X[:, 0]
    S = X[:, 1:]
    N = len(X)
    anc = [0] * N
    for j in range(1, N):
        dt = T[j] - T[:j]
        ds2 = ((S[:j] - S[j]) ** 2).sum(axis=1) if S.shape[1] else np.zeros(j)
        rel = np.nonzero((dt > 0) & (dt * dt >= ds2))[0]
        A = 0
        for i in rel:
            A |= 1 << int(i)
        anc[j] = A
    P = Poset(links_from_relation(anc), label=label)
    P._anc = anc
    return P


# ---- readings of a Poset -----------------------------------------------------

def interval_stats_poset(P: Poset, *, n_intervals: int = 60, min_size: int = 40, max_size: int = 4000,
                         n_pairs: int = 6000, seed: int = 0) -> Dict[str, Any]:
    """The E35 interval statistics of a general order (archive interval_stats, called with 40 / 4000 by E35 measure)."""
    anc, desc = P.anc(), P.desc()
    rng = random.Random(int(seed))
    rs, sizes, chains, dims = [], [], [], []
    tried = 0
    nodes = list(range(P.n))
    while len(rs) < int(n_intervals) and tried < int(n_pairs):
        tried += 1
        u = rng.choice(nodes)
        if not desc[u]:
            continue
        d = desc[u]
        k = rng.randrange(_popcount(d))
        w = None
        for i, b in enumerate(_bits(d)):
            if i == k:
                w = b
                break
        S = (desc[u] & anc[w]) | (1 << u) | (1 << w)
        m = _popcount(S)
        if not (int(min_size) <= m <= int(max_size)):
            continue
        rel = 0
        for x in _bits(S):
            rel += _popcount(anc[x] & S)
        r = rel / (m * (m - 1) / 2.0)
        hh: Dict[int, int] = {}
        best = 1
        for x in _bits(S):
            hx = 1 + max((hh[p] for p in P.parents[x] if p in hh), default=0)
            hh[x] = hx
            best = max(best, hx)
        rs.append(r)
        sizes.append(m)
        chains.append(best)
        dd = mm_dimension(r)
        if not math.isnan(dd):
            dims.append(dd)
    return _interval_summary(rs, sizes, chains, dims)


def _interval_summary(rs, sizes, chains, dims) -> Dict[str, Any]:
    """D_MM, scale halves, D_sd, D_chain from the sampled intervals (archive code, shared by both paths)."""
    nan = float("nan")
    out: Dict[str, Any] = {"n_used": len(rs), "D_MM": nan, "D_sd": nan, "D_chain": nan, "chain_slope": nan,
                           "size_mean": nan, "chain_mean": nan, "D_MM_small": nan, "D_MM_large": nan, "size_split": nan}
    if not rs:
        return out
    out["D_MM"] = mm_dimension(sum(rs) / len(rs))
    if len(rs) >= 8:
        order = sorted(range(len(rs)), key=lambda i: sizes[i])
        half = len(order) // 2
        lo = [rs[i] for i in order[:half]]
        hi = [rs[i] for i in order[half:]]
        out["D_MM_small"] = mm_dimension(sum(lo) / len(lo))
        out["D_MM_large"] = mm_dimension(sum(hi) / len(hi))
        out["size_split"] = sizes[order[half]]
    if len(dims) > 1:
        mu = sum(dims) / len(dims)
        out["D_sd"] = math.sqrt(sum((x - mu) ** 2 for x in dims) / (len(dims) - 1))
    out["size_mean"] = sum(sizes) / len(sizes)
    out["chain_mean"] = sum(chains) / len(chains)
    xs = [math.log(s) for s in sizes]
    ys = [math.log(c) for c in chains]
    if len(set(sizes)) > 2:
        mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
        sxx = sum((a - mx) ** 2 for a in xs)
        if sxx > 0:
            sl = sum((a - mx) * (b - my) for a, b in zip(xs, ys)) / sxx
            out["chain_slope"] = sl
            out["D_chain"] = 1.0 / sl if sl > 1e-6 else float("inf")
    return out


def homogeneity_poset(P: Poset, n_intervals: int = 90, min_size: int = 20, seed: int = 0) -> Dict[str, Any]:
    """E38-B homogeneity of a general order: D_MM by the height third of the interval's bottom."""
    anc, desc = P.anc(), P.desc()
    h = P.heights()
    H = max(h) if h else 0
    rng = random.Random(int(seed) ^ 0x4040)
    groups: Dict[int, List[float]] = {0: [], 1: [], 2: []}
    tried = 0
    got = 0
    while got < int(n_intervals) and tried < 40 * int(n_intervals):
        tried += 1
        u = rng.randrange(P.n)
        d = desc[u]
        if not d:
            continue
        k = rng.randrange(_popcount(d))
        w = None
        for i, b in enumerate(_bits(d)):
            if i == k:
                w = b
                break
        S = (desc[u] & anc[w]) | (1 << u) | (1 << w)
        m = _popcount(S)
        if m < int(min_size) or m > 4000:
            continue
        rel = sum(_popcount(anc[x] & S) for x in _bits(S))
        groups[min(2, int(3 * h[u] / H)) if H else 0].append(rel / (m * (m - 1) / 2.0))
        got += 1
    return _homog_summary(groups)


def _homog_summary(groups) -> Dict[str, Any]:
    dims = {g: (mm_dimension(sum(v) / len(v)) if len(v) >= 8 else float("nan")) for g, v in groups.items()}
    ok = [x for x in dims.values() if not math.isnan(x)]
    return {"D_by_third": [dims[0], dims[1], dims[2]], "n_by_third": [len(groups[k]) for k in (0, 1, 2)],
            "spread": (max(ok) - min(ok)) if len(ok) >= 2 else float("nan")}


def interval_counts_poset(anc: List[int], desc: List[int], x: int) -> Dict[int, int]:
    """N_n(x): ancestors y of x with exactly n elements strictly between, |anc(x) & desc(y)| (exact)."""
    A = anc[x]
    out: Dict[int, int] = {}
    for y in _bits(A):
        k = _popcount(A & desc[y])
        out[k] = out.get(k, 0) + 1
    return out


def smeared_kernel_scalar(n: int, d: int, eps: float) -> float:
    """eps * f_d(n, eps) for one n (archive cone_growth.smeared_kernel)."""
    C = BD_COEFFS[int(d)]
    if eps >= 1.0:
        return C[n] if n < len(C) else 0.0
    q = 1.0 - eps
    r = eps / q
    s = 0.0
    for i, c in enumerate(C):
        if i > n:
            break
        s += c * math.comb(n, i) * r ** i
    return eps * (q ** n) * s


def action_dS(counts: Dict[int, int], d: int, eps: float = ACTION_EPS) -> float:
    """dS_eps(x) = 1 + sum_n N_n(x) eps f_d(n, eps)."""
    return 1.0 + sum(cnt * smeared_kernel_scalar(n, d, eps) for n, cnt in counts.items())


def relation_density(P: Poset) -> float:
    n = P.n
    if n < 2:
        return float("nan")
    return sum(_popcount(a) for a in P.anc()) / (n * (n - 1) / 2.0)


def layer_profile(P: Poset) -> List[int]:
    h = P.heights()
    if not h:
        return []
    prof = [0] * (max(h) + 1)
    for x in h:
        prof[x] += 1
    return prof


def mean_links(P: Poset) -> float:
    return sum(len(p) for p in P.parents) / P.n if P.n else float("nan")


def read_poset(P: Poset, seed: int, n_intervals: int, n_action: int) -> Dict[str, Any]:
    """The E38-B `read` of a general order (the keys used by classify / subject / decide)."""
    I = interval_stats_poset(P, n_intervals=n_intervals, min_size=40, max_size=4000, seed=seed)
    h = P.heights()
    anc, desc = P.anc(), P.desc()
    upper = list(range(P.n // 2, P.n))
    rng = random.Random(int(seed))
    xs = rng.sample(upper, min(int(n_action), len(upper))) if upper else []
    counts = [interval_counts_poset(anc, desc, x) for x in xs]
    dS = {d: (sum(action_dS(c, d) for c in counts) / len(counts)) if counts else float("nan") for d in (2, 3, 4)}
    return {"n": P.n, "height": (max(h) + 1) if h else 0, "D_MM": I["D_MM"], "D_sd": I["D_sd"], "D_chain": I["D_chain"],
            "D_MM_small": I["D_MM_small"], "D_MM_large": I["D_MM_large"], "size_mean": I["size_mean"],
            "n_intervals": I["n_used"], "homog": homogeneity_poset(P, seed=seed), "dS": dS,
            "mean_links": mean_links(P), "rho": relation_density(P)}


# =============================================================================
# Matched null models (archive kernel.null_models, unchanged)
# =============================================================================

def _closure_poset(seed_lists: List[List[int]], label: str) -> Poset:
    anc: List[int] = []
    for seeds in seed_lists:
        A = 0
        for s in seeds:
            A |= anc[s] | (1 << s)
        anc.append(A)
    P = Poset(links_from_relation(anc), label=label)
    P._anc = anc
    return P


def transitive_percolation_np(n: int, p: float, seed: int = 0, originary: bool = True) -> Poset:
    """Originary transitive percolation (each earlier element directly below with probability p; >= 1 by rejection)."""
    rng = np.random.default_rng(int(seed))
    lists: List[List[int]] = [[]]
    for j in range(1, int(n)):
        while True:
            idx = np.nonzero(rng.random(j) < p)[0]
            if len(idx) or not originary:
                break
        lists.append([int(i) for i in idx])
    return _closure_poset(lists, f"percolation p={p:.3g}")


def tp_matched_density(target: Poset, seed: int = 0, rel_tol: float = 0.02, max_iter: int = 30):
    """Originary percolation at the target's N with p bisected (log p) to the target's relation density."""
    n = target.n
    goal = relation_density(target)
    lo, hi = math.log(1e-6), math.log(0.9)
    best = None
    for _ in range(int(max_iter)):
        mid = 0.5 * (lo + hi)
        P = transitive_percolation_np(n, math.exp(mid), seed=seed)
        rho = relation_density(P)
        err = abs(rho - goal) / goal if goal > 0 else float("inf")
        if best is None or err < best[2]:
            best = (math.exp(mid), P, err)
        if err <= rel_tol:
            break
        if rho < goal:
            lo = mid
        else:
            hi = mid
    p, P, err = best
    return P, {"null": "tp_matched_density", "p": p, "target_rho": goal, "achieved_rho": relation_density(P),
               "rel_error": err, "matched": err <= rel_tol}


def seed_sequence_null(seed_counts: Sequence[int], seed: int = 0):
    """Element x draws seed_counts[x] distinct earlier elements uniformly; its past is their down-closure."""
    rng = random.Random(int(seed))
    lists: List[List[int]] = [[]]
    for x in range(1, len(seed_counts)):
        k = max(1, min(int(seed_counts[x]), x))
        lists.append(rng.sample(range(x), k))
    P = _closure_poset(lists, "seed-sequence null")
    return P, {"null": "seed_sequence", "n": P.n, "matched": True}


def _layered(profile: Sequence[int], q: float, seed: int) -> Poset:
    rng = np.random.default_rng(int(seed))
    lists: List[List[int]] = []
    start: List[int] = []
    pos = 0
    for size in profile:
        start.append(pos)
        pos += size
    for h, size in enumerate(profile):
        for _ in range(size):
            if h == 0:
                lists.append([])
                continue
            prev0, prev1 = start[h - 1], start[h]
            seeds = [int(rng.integers(prev0, prev1))]
            lower = start[h]
            m = int(rng.binomial(lower, q)) if q > 0 else 0
            if m:
                seeds.extend(int(v) for v in rng.choice(lower, size=min(m, lower), replace=False))
            lists.append(seeds)
    return _closure_poset(lists, f"layered null q={q:.3g}")


def layered_null(target: Poset, seed: int = 0, rel_tol: float = 0.05):
    """Same layer profile, extra-seed rate q tuned to the target's mean link count (archive layered_null)."""
    prof = layer_profile(target)
    goal = mean_links(target)
    cache: Dict[float, Tuple[Poset, float]] = {}

    def ev(q):
        if q not in cache:
            P = _layered(prof, q, seed)
            cache[q] = (P, mean_links(P))
        return cache[q]

    def err_of(q):
        return abs(ev(q)[1] - goal) / goal if goal > 0 else float("inf")

    coarse = [0.0] + [10.0 ** e for e in range(-6, 1)]
    for q in coarse:
        ev(q)
    best_q = min(coarse, key=err_of)
    lo_q = hi_q = None
    for a, b in zip(coarse, coarse[1:]):
        if (ev(a)[1] - goal) * (ev(b)[1] - goal) <= 0:
            lo_q, hi_q = a, b
            break
    if lo_q is not None and hi_q is not None and hi_q > 0:
        lo_l = math.log(max(lo_q, 1e-7))
        hi_l = math.log(hi_q)
        up = ev(hi_q)[1] >= goal
        for _ in range(8):
            mid = math.exp(0.5 * (lo_l + hi_l))
            if err_of(mid) < err_of(best_q):
                best_q = mid
            if err_of(mid) <= 0.25 * rel_tol:
                break
            if (ev(mid)[1] < goal) == up:
                lo_l = math.log(mid)
            else:
                hi_l = math.log(mid)
    q, P, err = best_q, ev(best_q)[0], err_of(best_q)
    return P, {"null": "layered", "q": q, "target_links": goal, "achieved_links": mean_links(P), "rel_error": err,
               "profile_equal": layer_profile(P) == prof, "matched": err <= rel_tol and layer_profile(P) == prof}


def matched_nulls(P: Poset, seed: int) -> Dict[str, Poset]:
    """The three matched nulls of an order, with the archive's seed offsets (E54: s + 1, s + 2, s + 3)."""
    counts = [0] + [len(P.parents[x]) for x in range(1, P.n)]
    return {"tp": tp_matched_density(P, seed=seed + 1)[0], "seq": seed_sequence_null(counts, seed=seed + 2)[0],
            "layered": layered_null(P, seed=seed + 3)[0]}


# =============================================================================
# Sprinklings from coordinates (GTS Lab 3.x instruments.geometry_read)
# =============================================================================

def interval_reading(X: np.ndarray, u: int, w: int) -> Dict[str, Any]:
    """|S|, related pairs inside S and the longest chain of S = I(u, w) + {u, w}, from coordinates."""
    lo, hi = u + 1, w
    B = X[lo:hi]
    inside = lo + np.nonzero(prec_vector(X[u], B) & succ_vector(X[w], B))[0] if hi > lo else np.zeros(0, int)
    S = np.r_[u, inside, w].astype(np.int64)
    XS = X[S]
    return {"m": int(S.size), "rel": int(prec_matrix(XS, XS).sum()), "chain": int(chain_dp(XS)["W"].max())}


def _sample_pairs(X, rng, n_target, min_size, max_size, n_pairs):
    N = len(X)
    got, tried = 0, 0
    while got < n_target and tried < n_pairs:
        tried += 1
        u = rng.randrange(N)
        fut = u + 1 + np.nonzero(prec_vector(X[u], X[u + 1:]))[0]
        if fut.size == 0:
            continue
        w = int(fut[rng.randrange(fut.size)])
        B = X[u + 1:w]
        cnt = int((prec_vector(X[u], B) & succ_vector(X[w], B)).sum()) + 2
        if not (min_size <= cnt <= max_size):
            continue
        got += 1
        yield u, w, interval_reading(X, u, w)


def interval_stats_coords(X: np.ndarray, *, n_intervals: int = 60, min_size: int = 40, max_size: int = 4000,
                          n_pairs: int = 6000, seed: int = 0) -> Dict[str, Any]:
    rng = random.Random(int(seed))
    rs, sizes, chains, dims = [], [], [], []
    for u, w, rd in _sample_pairs(X, rng, n_intervals, min_size, max_size, n_pairs):
        m = rd["m"]
        r = rd["rel"] / (m * (m - 1) / 2.0)
        rs.append(r)
        sizes.append(m)
        chains.append(rd["chain"])
        dd = mm_dimension(r)
        if not math.isnan(dd):
            dims.append(dd)
    return _interval_summary(rs, sizes, chains, dims)


def homogeneity_coords(X: np.ndarray, heights: np.ndarray, *, n_intervals: int = 90, min_size: int = 20,
                       seed: int = 0) -> Dict[str, Any]:
    rng = random.Random(int(seed) ^ 0x4040)
    H = int(heights.max()) if heights.size else 0
    groups: Dict[int, List[float]] = {0: [], 1: [], 2: []}
    for u, w, rd in _sample_pairs(X, rng, n_intervals, min_size, 4000, 40 * n_intervals):
        m = rd["m"]
        groups[min(2, int(3 * heights[u] / H)) if H else 0].append(rd["rel"] / (m * (m - 1) / 2.0))
    return _homog_summary(groups)


@lru_cache(maxsize=None)
def _ncut(eps: float) -> int:
    return max(kernel_cutoff_n(eps, dd) for dd in (2, 3, 4))


def action_counts_coords(X: np.ndarray, x: int, rho: float, eps: float = ACTION_EPS) -> Dict[int, int]:
    """N_n(x) for every n up to the kernel cut (ancestors beyond it weigh < 1e-13 and are skipped)."""
    D = X.shape[1]
    ncut = _ncut(eps)
    lam = 1.5 * ncut + 50.0
    tc2 = (lam / (rho * interval_volume(D, 1.0))) ** (2.0 / D)
    dt, t2 = tau2_matrix(X[:x], X[x:x + 1])
    cand = np.nonzero((dt[:, 0] > 0) & (t2[:, 0] >= 0) & (t2[:, 0] <= tc2))[0]
    out: Dict[int, int] = {}
    if cand.size == 0:
        return out
    Xc = X[cand]
    nb = np.zeros(cand.size, dtype=np.int64)
    for c0 in range(0, cand.size, 1024):
        nb[c0:c0 + 1024] = prec_matrix(Xc[c0:c0 + 1024], Xc).sum(axis=1)
    vals, cnts = np.unique(nb[nb <= ncut], return_counts=True)
    return {int(v): int(c) for v, c in zip(vals, cnts)}


def read_coords(X: np.ndarray, rho: float, seed: int, n_intervals: int, n_action: int) -> Dict[str, Any]:
    """The E38-B `read` of a sprinkling from coordinates (same keys as read_poset; mean_links / rho not computed)."""
    heights = (chain_dp(X)["W"] - 1.0).astype(np.int64)
    I = interval_stats_coords(X, n_intervals=n_intervals, min_size=40, max_size=4000, seed=seed)
    rng = random.Random(int(seed))
    upper = list(range(len(X) // 2, len(X)))
    xs = rng.sample(upper, min(int(n_action), len(upper)))
    counts = [action_counts_coords(X, x, rho) for x in xs]
    dS = {d: float(np.mean([action_dS(c, d) for c in counts])) if counts else float("nan") for d in (2, 3, 4)}
    return {"n": len(X), "height": int(heights.max()) + 1, "D_MM": I["D_MM"], "D_sd": I["D_sd"], "D_chain": I["D_chain"],
            "D_MM_small": I["D_MM_small"], "D_MM_large": I["D_MM_large"], "size_mean": I["size_mean"],
            "n_intervals": I["n_used"], "homog": homogeneity_coords(X, heights, seed=seed), "dS": dS,
            "mean_links": float("nan"), "rho": float("nan")}


# =============================================================================
# Decision (archive E35 classify, E40 subject, E38-B decide)
# =============================================================================

def _m(xs):
    xs = [x for x in xs if x is not None and not (isinstance(x, float) and (math.isnan(x) or math.isinf(x)))]
    return statistics.mean(xs) if xs else float("nan")


def classify(rows_by_size: Dict[int, List[Dict[str, Any]]]) -> Tuple[str, Dict[str, Any]]:
    """MANIFOLD_LIKE(D) or NOT_MANIFOLD_LIKE(reasons) from a size series (E35)."""
    sizes = sorted(rows_by_size)
    big = rows_by_size[sizes[-1]]

    def m(k, rows):
        ok = [r[k] for r in rows if r[k] is not None and not math.isnan(r[k]) and not math.isinf(r[k])]
        return statistics.mean(ok) if ok else float("nan")
    D_MM, D_ch, sd = m("D_MM", big), m("D_chain", big), m("D_sd", big)
    small, large = m("D_MM_small", big), m("D_MM_large", big)
    per_size = [x for x in (m("D_MM", rows_by_size[n]) for n in sizes) if not math.isnan(x)]
    drift = (max(per_size) - min(per_size)) if len(per_size) >= 2 else float("nan")
    Hs = [statistics.mean(r["height"] for r in rows_by_size[n]) for n in sizes]
    Ns = [statistics.mean(r["n"] for r in rows_by_size[n]) for n in sizes]
    aH = ((math.log(Hs[-1]) - math.log(Hs[0])) / (math.log(Ns[-1]) - math.log(Ns[0]))
          if len(sizes) >= 2 and Hs[0] > 0 and Ns[-1] > Ns[0] else float("nan"))
    nums = {"D_MM": D_MM, "D_chain": D_ch, "D_sd": sd, "D_MM_small": small, "D_MM_large": large, "size_drift": drift,
            "height_exponent": aH, "D_height": (1.0 / aH) if aH and not math.isnan(aH) and aH > 0 else float("nan"),
            "D": (D_MM + D_ch) / 2.0 if not (math.isnan(D_MM) or math.isnan(D_ch)) else float("nan"),
            "heights": Hs, "sizes": Ns}
    fails = []
    if math.isnan(D_MM) or math.isnan(D_ch):
        fails.append("no usable intervals" if math.isnan(D_MM) else "chain estimator infinite")
    else:
        if math.isnan(sd) or sd > SD_MAX:
            fails.append(f"spread {sd:.2f} > {SD_MAX}")
        if abs(D_MM - D_ch) > AGREE_MAX:
            fails.append(f"D_MM {D_MM:.2f} vs D_chain {D_ch:.2f}")
        if math.isnan(small) or math.isnan(large) or abs(large - small) > SCALE_MAX:
            fails.append(f"scale drift {small:.2f} -> {large:.2f}")
        if not math.isnan(drift) and drift > SIZE_DRIFT_MAX:
            fails.append(f"size drift {drift:.2f}")
    if fails:
        return "NOT_MANIFOLD_LIKE(" + "; ".join(fails) + ")", nums
    return f"MANIFOLD_LIKE(D={nums['D']:.2f})", nums


def subject(rows_by_size: Dict[int, List[Dict[str, Any]]]) -> Dict[str, Any]:
    """The E38-B subject from readings (E40): label, numbers, homogeneity spread and dS of the largest size."""
    sizes = sorted(rows_by_size)
    lbl, nums = classify(rows_by_size)
    b = rows_by_size[sizes[-1]]
    return {"label": lbl, "numbers": nums, "homog_spread": _m([r["homog"]["spread"] for r in b]),
            "dS": {dd: _m([r["dS"][dd] for r in b]) for dd in (2, 3, 4)},
            "links": _m([r["mean_links"] for r in b]), "rho": _m([r["rho"] for r in b])}


def decide(subj: Dict[str, Any], nulls: Dict[str, Dict[str, Any]]) -> Tuple[str, Optional[int]]:
    """The E38-B decision (module notes)."""
    lbl, D = subj["label"], subj["numbers"].get("D", float("nan"))
    if not lbl.startswith("MANIFOLD_LIKE") or math.isnan(D):
        return "NOT_MANIFOLD_LIKE", None
    dstar = min((1, 2, 3, 4), key=lambda d: abs(d - D))
    if dstar == 1:
        return "CHAIN", None
    if math.isnan(subj["homog_spread"]) or subj["homog_spread"] > HOMOG_MAX:
        return f"INHOMOGENEOUS(spread {subj['homog_spread']:.2f})", dstar
    if math.isnan(subj["dS"][dstar]) or abs(subj["dS"][dstar]) > ACTION_TOL:
        return f"ACTION_MISMATCH(dS_{dstar} {subj['dS'][dstar]:.2f})", dstar
    hit = [k for k, v in nulls.items() if v["label"].startswith("MANIFOLD_LIKE")
           and not math.isnan(v["numbers"].get("D", float("nan"))) and abs(v["numbers"]["D"] - D) <= D_TOL]
    if hit:
        return "EXPLAINED_BY_NULL(" + ",".join(hit) + ")", dstar
    return f"GEOMETRIC(D={D:.2f}, d*={dstar})", dstar
