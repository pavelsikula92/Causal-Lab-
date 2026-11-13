"""Summary of the B3 development rows (diag/out/b3_dev_cache or a saved b3_dev_*.json). Decides nothing."""
from __future__ import annotations

import glob
import json
import math
import sys

import numpy as np


def ms(a):
    a = np.asarray([v for v in a if np.isfinite(v)], float)
    if a.size < 2:
        return f"{a.mean() if a.size else float('nan'):+.4f} (n {a.size})"
    return f"{a.mean():+.4f} +- {a.std(ddof=1) / math.sqrt(a.size):.4f} (sd {a.std(ddof=1):.4f}, n {a.size})"


def main():
    if len(sys.argv) > 1:
        R = json.load(open(sys.argv[1]))["results"]
    else:
        R = [json.load(open(f)) for f in glob.glob("diag/out/b3_dev_cache/*.json")]
    rho = sorted({r["rho"] for r in R})
    for rh in rho:
        for mu in sorted({r["muA"] for r in R}, reverse=True):
            K = {r["k"]: r for r in R if r["rho"] == rh and r["muA"] == mu and r["kappa"] > 0}
            Z = {r["k"]: r for r in R if r["rho"] == rh and r["muA"] == mu and r["kappa"] == 0}
            if not K:
                continue
            ks = sorted(K)
            print(f"== rho {rh} muA {mu} kappa {K[ks[0]]['kappa']}: N {np.mean([K[k]['N'] for k in ks]):.0f} "
                  f"(cand {np.mean([K[k]['Ncand'] for k in ks]):.0f}), {np.mean([K[k]['sec'] for k in ks]):.0f} s; "
                  f"min L {min(K[k]['Lmin_seen'] for k in ks):+.3f}, overflow {sum(K[k]['over'] for k in ks)}")
            late = [i for i, w in enumerate(K[ks[0]]["rows"]) if w["t"] >= 5.0]
            rL = [sum(K[k]["rows"][i]["LL"] for i in late) / sum(K[k]["rows"][i]["LcLc"] for i in late) for k in ks]
            print(f"   r_L (t >= 5, |x| < 6, band) {ms(rL)}")
            for i in late:
                print(f"      t {K[ks[0]]['rows'][i]['t']}: r_L {ms([K[k]['rows'][i]['LL'] / K[k]['rows'][i]['LcLc'] for k in ks])}")
            b = np.array([[bb[:2] for bb in K[k]["rows"][late[-1]]["bins"]] for k in ks])
            print("   L by x (last read slab), discrete: " + " ".join(f"{v:+.3f}" for v in np.nanmean(b[:, :, 0], 0)))
            print("                          continuum: " + " ".join(f"{v:+.3f}" for v in np.nanmean(b[:, :, 1], 0)))
            xc = [K[k]["rows"][late[-1]]["xc"] - K[k]["rows"][0]["xc"] for k in ks]
            xz = [Z[k]["rows"][late[-1]]["xc"] - Z[k]["rows"][0]["xc"] for k in ks if k in Z]
            print(f"   kink drift (last - first) kappa {ms(xc)}; kappa 0 {ms(xz)}; continuum "
                  f"{K[ks[0]]['rows'][late[-1]]['xc_ex'] - K[ks[0]]['rows'][0]['xc_ex']:+.4f}")
            wr = [K[k]["rows"][late[-1]]["w"] / K[k]["rows"][late[-1]]["w_ex"] - 1 for k in ks]
            wz = [Z[k]["rows"][late[-1]]["w"] / Z[k]["rows"][late[-1]]["w_ex"] - 1 for k in ks if k in Z]
            print(f"   width last slab w/w_ex - 1: kappa {ms(wr)} (w_ex {K[ks[0]]['rows'][late[-1]]['w_ex']:.3f}); kappa 0 {ms(wz)}")
            print(f"   density realised / expected: {ms([K[k]['dens'][0] / K[k]['dens'][1] - 1 for k in ks])}")
            for what, key, ckey in (("clocks (chain length)", "clen", "tau_cont"), ("cones (reach in y)", "reach", "Pi2_mid")):
                rat, cr = [], []
                for k in ks:
                    if k not in Z:
                        continue
                    pk, pz = K[k]["probe"], Z[k]["probe"]
                    if not all(n in pk and n in pz for n in ("near", "far")):
                        continue
                    a = np.mean(pk["near"][key]) / np.mean(pk["far"][key])
                    z = np.mean(pz["near"][key]) / np.mean(pz["far"][key])
                    rat.append(a / z)
                    cr.append(np.mean(pk["near"][ckey]) / np.mean(pk["far"][ckey]))
                print(f"   {what}: (near/far)_kappa / (near/far)_0 {ms(rat)}; continuum {np.mean(cr):.4f}")


if __name__ == "__main__":
    main()
