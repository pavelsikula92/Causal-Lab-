"""
OV3-lin-a from the frozen records (development; decisions 41, 44, 46): the static linear response of the order to a
dust slab and to a BPS wall, per k, in the young-field reading (n = kT = 4) and for the bare ensemble action.
==================================================================================================================
Inputs, read from the result files (nothing typed by hand):
  F6.8g  (json): K_hat_c, K_hat_ab, K_hat_bb at n = 4, at k = 0.06 and 0.12 (data["ov3"]);
  F6.10b (json + txt): the source form factors F_omega, F_nu at the same k (json: full precision; txt: the variants at
         age T = 80 and grid 1.5, four decimals).
The solver is diag/ov3_lin.py (local limit = linearised GR, verified with sympy). Every output is a ratio to
linearised GR, so alpha = -1/(8 pi G) drops out.
The systematic uncertainty of the source (decision 41: the wall's lapse source -4 F_nu + 3 F_omega must carry the
covariance of the two form factors) is taken directly: the whole response is recomputed with the form factors of each
F6.10b variant (age, grid), so their correlation enters automatically. The stiffness side varies by <= 3e-5 between
levels and platforms (F6.8g I3, decision 46) and is not propagated.

    python -m diag.ov3_lin_a --f68g <F6.8g json> --f610b <F6.10b json>      (the txt next to the json is read too)
"""
import argparse
import json
import re
from pathlib import Path

from diag.ov3_lin import response


def f610b_forms(jpath):
    d = json.load(open(jpath))["data"]["ratios"]
    main = {s: {float(k): v["ratio"] for k, v in d[s].items()} for s in ("omega", "nu")}
    txt = Path(jpath).with_suffix(".txt").read_text()
    var = {"T80": {"omega": {}, "nu": {}}, "grid": {"omega": {}, "nu": {}}}
    for s, k, r, ra, rg in re.findall(r"(omega|nu)\s+k ([0-9.]+)\s*:.*?ratio ([0-9.]+); T 80: ([0-9.]+); grid 1\.5: "
                                      r"([0-9.]+)", txt):
        var["T80"][s][float(k)] = float(ra)
        var["grid"][s][float(k)] = float(rg)
        var.setdefault("main4", {"omega": {}, "nu": {}})[s][float(k)] = float(r)
    return main, var


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--f68g", required=True)
    ap.add_argument("--f610b", required=True)
    q = ap.parse_args()
    g = json.load(open(q.f68g))
    assert not g.get("quick"), "a quick record"
    ov3 = {float(k): v for k, v in g["data"]["ov3"].items()}
    Fm, var = f610b_forms(q.f610b)
    print(f"F6.8g {Path(q.f68g).name} ({g['status']}, sha {g['test_sha256'][:8]}); "
          f"F6.10b {Path(q.f610b).name}")
    print("k     | K       K_ab    K_bb    | F_om    F_nu    | dust Phi/GR  Psi/GR  gamma  | wall Phi/GR  Psi/Phi  "
          "wall/dust | wall source -4F_nu+3F_om (GR -1)")
    rows = {}
    for k in sorted(ov3):
        K, Kab, Kbb = ov3[k]["K_c"], ov3[k]["K_ab"], ov3[k]["K_bb"]
        Fw, Fn = Fm["omega"][k], Fm["nu"][k]
        r = response(K, Kab, Kbb, Fw, Fn)
        d, w = r["dust"], r["wall"]
        src = -4 * Fn + 3 * Fw
        sysd = {}
        for vname in ("T80", "grid"):
            fw4, fn4 = var[vname]["omega"][k], var[vname]["nu"][k]
            fw0, fn0 = var["main4"]["omega"][k], var["main4"]["nu"][k]
            rv = response(K, Kab, Kbb, Fw + (fw4 - fw0), Fn + (fn4 - fn0))
            sysd[vname] = {"dust_Phi": rv["dust"]["Phi/GR"] - d["Phi/GR"], "gamma": rv["dust"]["Psi/Phi"] - d["Psi/Phi"],
                           "wall_Phi": rv["wall"]["Phi/GR"] - w["Phi/GR"],
                           "src": (-4 * (fn4 - fn0) + 3 * (fw4 - fw0))}
        sd = {key: sum(abs(sysd[v][key]) for v in sysd) for key in ("dust_Phi", "gamma", "wall_Phi", "src")}
        rows[k] = {"dust": d, "wall": w, "src": src, "sys": sd}
        print(f"{k:<5g} | {K:.5f} {Kab:.5f} {Kbb:.5f} | {Fw:.5f} {Fn:.5f} | {d['Phi/GR']:.4f}±{sd['dust_Phi']:.4f}  "
              f"{d['Psi/GR']:.4f}  {d['Psi/Phi']:.4f}±{sd['gamma']:.4f} | {w['Phi/GR']:.4f}±{sd['wall_Phi']:.4f}  "
              f"{w['Psi/Phi']:+.4f}  {w['Phi'] / d['Phi']:+.4f} | {src:+.4f}±{sd['src']:.4f}")
    print("(± = |change| summed over the F6.10b variants age T 80 and grid 1.5, the correlated systematic of the source)")
    return rows


if __name__ == "__main__":
    main()
