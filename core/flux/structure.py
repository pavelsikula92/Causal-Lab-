"""
core.flux.structure
===================

Statistics of the arrangement of saturated events (F3.2, F3.3), ported
statement by statement from GTS Lab 3.x E55b (fields of the analysis
events of a cylinder sprinkling).

The analysis events are |x| < R_abs, t in [t1, t2). chi_v = 1 for a
saturated event (dQ_v > 0), else 0.

cells(X, chi, mask, p)
    per time slab of width `slab` in [t1, t2): {cell key: (saturated
    fraction, count)} on square cells of side `cell` in the (x1, x2) plane
rs_successive(seq, min_count)
    mean Spearman r_s between the cell values of successive slabs over the
    cells present in both with >= min_count events (nan if undefined)
moran(cells_of_slab, min_count)
    Moran's I of the cell values, rook neighbours (nan if undefined)
localization(X, r, t, chi, mask, p, rng)
    E_loc = Moran's I minus its mean under n_perm permutations of chi
    among the analysis events of the same radial bin (n_rbins equal bins of
    |x| in [0, R_abs)) and slab: keeps the radial and temporal profile,
    destroys local clustering; averaged over slabs. One value per sprinkling.
pairs_persistence(cells_list, min_count)
    Delta_pers per disjoint pair (2j, 2j+1): within = mean of the two
    sprinklings' own r_s; between = slabs taken alternately from the two
    sprinklings (keeps the common envelope, breaks sprinkling-specific
    structure); Delta = within - between. One value per pair.
"""

from __future__ import annotations

import warnings
from typing import Any, Dict, List, Tuple

import numpy as np

Cells = List[Dict[int, Tuple[float, int]]]


def _slab_edges(p: Dict[str, Any]) -> np.ndarray:
    return np.arange(float(p["t1"]), float(p["t2"]) + 1e-9, float(p["slab"]))


def cells(X: np.ndarray, chi: np.ndarray, mask: np.ndarray, p: Dict[str, Any]) -> Cells:
    """Per slab: {cell key: (saturated fraction, count)} for the analysis events (mask)."""
    t = X[:, 0]
    cell, Rabs = float(p["cell"]), float(p["R_abs"])
    g = np.arange(-Rabs, Rabs + cell, cell)
    key = np.digitize(X[:, 1], g) * 1000 + np.digitize(X[:, 2], g)
    slabs = _slab_edges(p)
    out = []
    for k in range(len(slabs) - 1):
        m = mask & (t >= slabs[k]) & (t < slabs[k + 1])
        u, inv = np.unique(key[m], return_inverse=True)
        n = np.bincount(inv, minlength=u.size)
        f = np.bincount(inv, weights=chi[m], minlength=u.size) / np.maximum(n, 1)
        out.append({int(a): (float(b), int(c)) for a, b, c in zip(u, f, n)})
    return out


def rs_successive(seq: Cells, min_count: int) -> float:
    """Mean Spearman r_s between successive slabs over the common cells (nan if undefined)."""
    from scipy.stats import spearmanr
    rs = []
    for k in range(len(seq) - 1):
        common = [c for c in seq[k] if c in seq[k + 1] and seq[k][c][1] >= min_count and seq[k + 1][c][1] >= min_count]
        if len(common) > 3:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")            # a constant slab gives nan
                v = spearmanr([seq[k][c][0] for c in common], [seq[k + 1][c][0] for c in common])[0]
            rs.append(v)
    rs = [x for x in rs if np.isfinite(x)]
    return float(np.mean(rs)) if rs else float("nan")


def moran(fr: Dict[int, Tuple[float, int]], min_count: int) -> float:
    """Moran's I of the cell fractions with rook neighbours (cells with >= min_count events)."""
    vals = {c: v[0] for c, v in fr.items() if v[1] >= min_count}
    if len(vals) < 4:
        return float("nan")
    keys = list(vals)
    z = np.array([vals[c] for c in keys])
    z = z - z.mean()
    den = float((z * z).sum())
    if den == 0.0:
        return float("nan")
    pos = {c: i for i, c in enumerate(keys)}
    num, W = 0.0, 0
    for c, i in pos.items():
        for d in (1, -1, 1000, -1000):
            j = pos.get(c + d)
            if j is not None:
                num += z[i] * z[j]
                W += 1
    return float(len(keys) / W * num / den) if W else float("nan")


def localization(X: np.ndarray, r: np.ndarray, chi: np.ndarray, mask: np.ndarray, p: Dict[str, Any],
                 rng: np.random.Generator) -> Dict[str, float]:
    """E_loc = Moran's I minus its mean under radial-bin-and-slab permutations, averaged over slabs."""
    t = X[:, 0]
    mc, nb = int(p["min_count"]), int(p["n_rbins"])
    rb = np.minimum((r / float(p["R_abs"]) * nb).astype(int), nb - 1)
    slabs = _slab_edges(p)
    sl = np.clip(np.searchsorted(slabs, t, side="right") - 1, 0, len(slabs) - 2)
    groups = [np.where(mask & (sl == k) & (rb == b))[0] for k in range(len(slabs) - 1) for b in range(nb)]
    obs = cells(X, chi, mask, p)
    I_obs = [moran(f, mc) for f in obs]
    I_perm = np.zeros((int(p["n_perm"]), len(obs)))
    for q in range(int(p["n_perm"])):
        c2 = chi.copy()
        for g in groups:
            if g.size > 1:
                c2[g] = chi[rng.permutation(g)]
        I_perm[q] = [moran(f, mc) for f in cells(X, c2, mask, p)]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ex = np.array(I_obs) - np.nanmean(I_perm, axis=0)
        I_mean = float(np.nanmean(I_obs)) if np.isfinite(I_obs).any() else float("nan")
    ex = ex[np.isfinite(ex)]
    return {"E_loc": float(ex.mean()) if ex.size else float("nan"), "I_obs": I_mean}


def pairs_persistence(cells_list: List[Cells], min_count: int) -> List[Dict[str, float]]:
    """Delta_pers per disjoint pair of sprinklings (module notes); nan where undefined."""
    out = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for j in range(len(cells_list) // 2):
            A, B = cells_list[2 * j], cells_list[2 * j + 1]
            wA, wB = rs_successive(A, min_count), rs_successive(B, min_count)
            mix = [A[k] if k % 2 == 0 else B[k] for k in range(len(A))]
            mix2 = [B[k] if k % 2 == 0 else A[k] for k in range(len(A))]
            within = float(np.nanmean([wA, wB]))
            between = float(np.nanmean([rs_successive(mix, min_count), rs_successive(mix2, min_count)]))
            out.append({"within": within, "between": between, "delta": within - between})
    return out
