"""Compact summary of the C0 stability runs: python -m diag.c0_stab_summary [file.json] (default: the newest)."""
import glob
import json
import math
import sys

import numpy as np

fn = sys.argv[1] if len(sys.argv) > 1 else sorted(glob.glob("diag/out/c0_stab_2*.json"))[-1]
R = json.load(open(fn))
print(fn)
for kn in sorted({r["kernel"] for r in R}):
    S = [r for r in R if r["kernel"] == kn]
    print(f"=== {kn} (n {len(S)}, pairs/event {np.mean([r['pairs'] / r['N'] for r in S]):.0f})")
    for md in S[0]["modes"]:
        n = min(len(r["modes"][md]) for r in S)
        T = np.array([w["t"] for w in S[0]["modes"][md][:n]])
        sd = np.mean([[w["sd"] for w in r["modes"][md][:n]] for r in S], 0)
        mn = np.mean([[w["mean"] for w in r["modes"][md][:n]] for r in S], 0)
        idx = [i for i, t in enumerate(T) if abs(t - round(t) - 0.125) < 1e-6 and t < 7.0]
        ls = np.log(sd)
        g = [(ls[i + 4] - ls[i]) for i in range(4, len(ls) - 4, 4) if T[i + 4] < 7.0]
        print(f"  {md:7s} sd  " + " ".join(f"{sd[i]:.2f}" for i in idx) + "   mean " + " ".join(f"{mn[i]:+.2f}" for i in idx)
              + "   dln(sd)/dt " + " ".join(f"{v:+.2f}" for v in g))
