"""
C1b: can the rows tell the exact mean from the proper-time kernel? (15 Nov; the author's decision of 08:04)
=======================================================================================================
A projection statistic that uses each pair as one sample, so the covariance between bins (the same sprinkling enters
every bin) is carried by the spread over pairs instead of being assumed away:

    s_k = sum_j w_j (D_kj - P1_kj) / sum_j w_j dP_kj,    dP = P0 - P1,  w_j = dP_j / sigma_j^2,

D_kj the curved - flat difference of pair k in bin j, P0 / P1 the exact mean / proper-time kernel averaged over the
events of that pair and bin (minus the flat exact mean over the flat events). E s = 1 if the exact mean holds, 0 if the
proper-time kernel does; SE(s) = sd_k(s_k) / sqrt(n). sigma_j: the spread over pairs of D - P0 per bin on the a 0.05
development rows (the noise does not depend on a to first order).

The difference of two chi^2 values is not used: both models are compared with the same data, so the two chi^2 are
strongly correlated, and their difference has no simple null distribution.

    python -m diag.c1b_discrim
"""
from __future__ import annotations

import json
import math
import os

import numpy as np

TB, RB, NT, NR = 0.5, 0.5, 12, 9          # coarse bins of compare(): t [0, 6], r [0, 4.5]
T_MIN = 1.0


def load_pairs(fn=os.path.join("diag", "out", "c1b_pairs_rho60_a0.05.json")):
    return json.load(open(fn))


def coarse(F, C):
    """Event-weighted coarse-bin sums of a fine-grid field F (fine cells t 0.05, r 0.1) with counts C."""
    it = np.arange(C.shape[0])[:, None] * 0.05
    ir = np.arange(C.shape[1])[None, :] * 0.1
    ib = np.minimum((it / TB + 1e-9).astype(int) + 0 * ir.astype(int), NT)
    jb = (ir / RB + 1e-9).astype(int) + 0 * it.astype(int)
    S = np.zeros((NT + 1, NR + 1))
    np.add.at(S, (ib.ravel(), np.minimum(jb.ravel(), NR)), (F[: C.shape[0], : C.shape[1]] * C).ravel())
    return S[:NT, :NR]


def per_pair(pairs, f0, f1, ff):
    """D, P0, P1 per pair and coarse bin (n, NT, NR)."""
    D, P0, P1 = [], [], []
    for p in pairs:
        Cc, Sc = np.array(p["curved"]["C"], float), np.array(p["curved"]["S"])
        Cf, Sf = np.array(p["flat"]["C"], float), np.array(p["flat"]["S"])
        Nc, Nf = coarse(np.ones_like(Cc), Cc), coarse(np.ones_like(Cf), Cf)
        with np.errstate(invalid="ignore", divide="ignore"):
            yc, yf = coarse(np.ones_like(Cc), Sc) / Nc, coarse(np.ones_like(Cf), Sf) / Nf
            pf = coarse(ff, Cf) / Nf
            D.append(yc - yf)
            P0.append(coarse(f0, Cc) / Nc - pf)
            P1.append(coarse(f1, Cc) / Nc - pf)
    return np.array(D), np.array(P0), np.array(P1)


def sel_mask(D, region):
    tc = (np.arange(NT) + 0.5) * TB
    rc = (np.arange(NR) + 0.5) * RB
    m = np.all(np.isfinite(D), 0) & (tc[:, None] >= T_MIN)
    if region == "inside":
        m &= rc[None, :] < 1.5
    elif region == "vacuum":
        m &= rc[None, :] > 1.5
    return m


def stat(D, P0, P1, sigma, m):
    dP = P0 - P1
    w = np.where(m, dP.mean(0) / sigma ** 2, 0.0)
    s = (w * np.nan_to_num(D - P1)).sum((1, 2)) / (w * np.nan_to_num(dP)).sum((1, 2))
    return s.mean(), s.std(ddof=1) / math.sqrt(len(s)), s


def main():
    pairs = load_pairs()
    Z = np.load(os.path.join("diag", "out", "c1b_scan_rho60.npz"))
    ff = Z["flat"]
    D, P0, P1 = per_pair(pairs, Z["f0_0.05"], Z["f1_0.05"], ff)
    resid = D - P0
    sigma = np.nanstd(resid, 0, ddof=1)
    sigma = np.where(np.isfinite(sigma) & (sigma > 0), sigma, np.inf)
    lines = ["C1b discrimination, rho 60; per-pair projection statistic s (1 = exact mean, 0 = proper-time kernel)"]
    lines.append("measured on the a 0.05 development rows (8 pairs):")
    for reg in ("all", "inside", "vacuum"):
        m = sel_mask(D, reg)
        mu, se, _ = stat(D, P0, P1, sigma, m)
        lines.append(f"  {reg:7s}: s = {mu:+.2f} +- {se:.2f}  (separation of the two models {1 / se:.1f} SE)")
    lines.append("expected separation at other a (noise of the a 0.05 rows, the predicted dP of a; 8 and 16 pairs):")
    for a in (0.025, 0.05, 0.10, 0.15):
        if f"f0_{a}" not in Z:
            continue
        _, Q0, Q1 = per_pair(pairs, Z[f"f0_{a}"], Z[f"f1_{a}"], ff)
        dP = Q0 - Q1
        row = []
        for reg in ("all", "inside", "vacuum"):
            m = sel_mask(D, reg)
            w = np.where(m, dP.mean(0) / sigma ** 2, 0.0)
            num = (w * np.nan_to_num(resid)).sum((1, 2))
            den = (w * np.nan_to_num(dP)).sum((1, 2))
            se = (num / den).std(ddof=1) / math.sqrt(len(num))
            row.append(f"{reg} {1 / se:.1f} / {1 / (se / math.sqrt(2)):.1f}")
        mx = np.nanmax(np.abs(np.where(sel_mask(D, 'all'), dP.mean(0), np.nan)))
        lines.append(f"  a {a:<5}: max |dP| {mx:.3f}; separation in SE: " + "; ".join(row))
    txt = "\n".join(lines)
    print(txt)
    open(os.path.join("diag", "out", "c1b_discrim_rho60.txt"), "w").write(txt + "\n")


if __name__ == "__main__":
    main()
