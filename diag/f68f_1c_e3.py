"""
F6.8f-dev, step 1c: the b -> 0 check of the series against F6.8e3 (development, decides nothing)
================================================================================================

docs/derivations/f68f_dev.md, step 1c. At large n the F6.8e3 bulk density per direction must be

    bulk(p) = 2 mu_3 k^5 n^5 C2(p) / l2^4,   C2(p) = lim_{eps -> 0} < s4hat A0^-2 >_theta / eps^2,

with no free constant (corrections O(1/n^2), O(k^2)). Reads the records of diag.f68f_1c (two small eps, C2 by
Richardson) and of diag.f68f_layer --bulk-only (the same p nodes) given on the command line, and compares them.

Run from the repository root:
    python -m diag.f68f_1c_e3 --series diag/out/f68f_1c_<time>.json --bulk diag/out/f68f_layer_<time>.json
Output diag/out/f68f_1c_e3_<UTC time>.json and .txt.
"""
from __future__ import annotations

import argparse
import json
import os
import time

import numpy as np

from diag import f68f_common as C

L2 = 0.3618006272791338
MU3 = 0.5


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--series", required=True)
    ap.add_argument("--bulk", required=True)
    a = ap.parse_args()
    t0 = time.time()
    S = json.load(open(a.series))
    Bk = json.load(open(a.bulk))
    pts = sorted(S["data"]["points"], key=lambda d: d["eps"])
    e1, e2 = pts[0]["eps"], pts[1]["eps"]
    c1 = np.array(pts[0]["profile"]) / e1 ** 2
    c2 = np.array(pts[1]["profile"]) / e2 ** 2
    C2 = (4 * c1 - c2) / 3 if abs(e2 - 2 * e1) < 1e-12 else c1
    row = Bk["data"][0]
    k, n = float(Bk["args"]["k"]), float(row["n"])
    p, w, bulk = np.array(row["p"]), np.array(row["w"]), np.array(row["bulk"])
    ps = np.array(pts[0]["p"])
    if p.size != ps.size or np.max(np.abs(p - ps)) > 1e-12:
        raise SystemExit("the p nodes of the two records differ")
    pred = 2 * MU3 * k ** 5 * n ** 5 * C2 / L2 ** 4
    lines = [f"F6.8f-dev step 1c: b -> 0 check against F6.8e3 at k = {k:g}, n = {n:g} (C2 from eps {e1:g}, {e2:g}); "
             "development, decides nothing", "",
             "     p        F6.8e3 bulk        series prediction     ratio"]
    for x, b, q in zip(p, bulk, pred):
        lines.append(f"  {x:7.3f}  {b:+.8e}    {q:+.8e}    {b / q:.6f}")
    tot = float((w * bulk).sum()) / float((w * pred).sum())
    lines += ["", f"  p-integrated: F6.8e3 / series = {tot:.6f} (relative difference {tot - 1:+.2e})",
              f"  C2 from eps {e1:g} alone against Richardson: largest relative difference "
              f"{float(np.max(np.abs(c1 / C2 - 1))):.1e}",
              f"  records: {a.series} (code {S['code_md5']}), {a.bulk} (code {Bk['code_md5']})"]
    stem = C.save("f68f_1c_e3", vars(a), C.code_md5([os.path.abspath(__file__)]),
                  {"p": p, "bulk": bulk, "pred": pred, "ratio_total": tot}, lines, t0)
    print("\n".join(lines))
    print("saved", stem + ".json")


if __name__ == "__main__":
    main()
