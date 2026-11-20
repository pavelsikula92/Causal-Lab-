"""OV1 2D development summary (decides nothing): python -m diag.ov1_2d_summ [eps ...]. Mean +- SE over sprinklings.
Seeds >= CAL (100000) are the calibration set of the phase bands (the author, 19 Nov: calibration seeds separate from
the validation seeds); everything else is read from the validation seeds (< CAL)."""
import glob
import json
import os
import sys

import numpy as np

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "ov1_2d")
CAL = 100000


def ms(x):
    x = np.asarray(x, float)
    return x.mean(), (x.std(ddof=1) / np.sqrt(x.size) if x.size > 1 else float("nan")), (x.std(ddof=1) if x.size > 1 else float("nan"))


def main():
    eps_list = [float(e) for e in sys.argv[1:]] or sorted({json.load(open(f))["eps"] for f in glob.glob(f"{OUT}/ov1_*.json")})
    for eps in eps_list:
        ALL = [json.load(open(f)) for f in sorted(glob.glob(f"{OUT}/ov1_{eps:g}_*.json"))]
        R = [r for r in ALL if r["seed"] < CAL]
        RC = [r for r in ALL if r["seed"] >= CAL]
        if not R:
            continue
        levs = sorted(R[0]["levels"], key=int)
        print(f"eps {eps:g}: n_loc {R[0]['n_loc']}, n_cut {R[0]['n_cut']}, {len(R)} sprinklings, "
              f"{np.mean([r['sec'] for r in R]):.0f} s each")
        print("  vacuum (L1): <S> per level, sd per sprinkling")
        for L in levs:
            S = [r["levels"][L]["S"] for r in R]
            N = [r["levels"][L]["N"] for r in R]
            m, se, sd = ms(S)
            print(f"    T/w {L:>2}: N {np.mean(N):8.0f}  <S> {m:+9.3f} +- {se:7.3f}  (sd {sd:7.3f})")
        for a, b in zip(levs[:-1], levs[1:]):
            d = [r["levels"][b]["S"] - r["levels"][a]["S"] for r in R]
            m, se, sd = ms(d)
            print(f"    S({b}) - S({a}): {m:+9.3f} +- {se:7.3f}  (sd {sd:7.3f})")
        c = []
        for r in R:
            N = np.array([r["levels"][L]["N"] for L in levs], float)
            S = np.array([r["levels"][L]["S"] for L in levs])
            c.append(np.polyfit(N, S, 1)[0])
        m, se, sd = ms(c)
        Nmax = np.mean([r["levels"][levs[-1]]["N"] for r in R])
        print(f"    slope c = dS/dN: {m:+.3e} +- {se:.3e}   (c N_max {m * Nmax:+.3f} +- {se * Nmax:.3f})")
        src = RC if len(RC) >= 20 else R
        print(f"  phase bands (flat) from the {'calibration' if src is RC else 'VALIDATION (no calibration set yet)'} "
              f"seeds (n {len(src)}): mean, sd, [0.5 %, 99.5 %]; validation outside the band")
        for L in levs:
            for q in ("chi", "d_mm", "h"):
                x = np.array([r["levels"][L][q] for r in src], float)
                lo, hi = np.percentile(x, [0.5, 99.5])
                y = np.array([r["levels"][L][q] for r in R], float)
                out = np.mean((y < lo) | (y > hi)) if src is RC else float("nan")
                print(f"    T/w {L:>2} {q:4}: {x.mean():.5f} sd {x.std(ddof=1):.5f} [{lo:.5f}, {hi:.5f}]  "
                      f"validation outside {out:.3f}")
        print("  balance (L2): interior elements; per-sprinkling means, then mean +- SE over sprinklings")
        for L in levs:
            RB = [r for r in R if "balance" in r["levels"][L]]
            if len(RB) < 2:
                continue
            for k in ("U", "V"):
                B = [r["levels"][L]["balance"][k] for r in RB]
                Dm = ms([b["D_mean"] for b in B])
                Dt = ms([b["D_tslope"] for b in B])
                D2 = ms([b["D2_mean"] for b in B])
                print(f"    T/w {L:>2} {k}: n_inner {np.mean([b['n_inner'] for b in B]):.0f}  "
                      f"<D S> {Dm[0]:+.5f} +- {Dm[1]:.5f} (per-element sd {np.mean([b['D_sd'] for b in B]):.4f})  "
                      f"d<D S>/dt {Dt[0]:+.5f} +- {Dt[1]:.5f}   <D2 S> {D2[0]:+.5f} +- {D2[1]:.5f} "
                      f"(sd {np.mean([b['D2_sd'] for b in B]):.4f})")
                nb_ = len(B[0]["profile"])
                prof = [ms([b["profile"][j][0] for b in B if b["profile"][j][0] is not None]) for j in range(nb_)]
                print("      profile of <D S> in time from the past tip (8 bins, all elements): "
                      + " ".join(f"{p[0]:+.4f}({p[1]:.4f})" for p in prof))


if __name__ == "__main__":
    main()
