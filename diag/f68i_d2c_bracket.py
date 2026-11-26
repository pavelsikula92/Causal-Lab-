"""
F6.8i-dev, check for f68i_d2c revision 2.2 (development, decides nothing; docs/derivations/f68i_dev.md section 6.5).

At given (eps, p_c) computes the needle A(q), N'(q) of f68i_d2c twice: with the plain window of
f68f_1b_exact.phase_A and with the bracketing guide rays psi +- BRACKET W(s). Where the plain window works, the two
must agree (the bracket only moves the window edges outward, where Phi > q^2/2 contributes nothing).
Usage: python -m diag.f68i_d2c_bracket --pairs 30:5.2023443,20:4.3813227
"""
import argparse
import json
import os
import time

from diag import f68i_d2c as D

C = D.C


def job(pair):
    eps, pc = pair
    out = {"eps": eps, "p_c": pc}
    for key, br in (("plain", 0.0), ("bracket", D.BRACKET)):
        t = time.time()
        try:
            out[key] = D.needle(eps, pc, br)
        except RuntimeError as e:
            out[key] = {"error": str(e)}
        out[key]["seconds"] = time.time() - t
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", type=str, required=True)
    ap.add_argument("--workers", type=int, default=0)
    a = ap.parse_args()
    pairs = [tuple(float(v) for v in x.split(":")) for x in a.pairs.split(",")]
    t0 = time.time()
    code = C.code_md5([os.path.abspath(__file__), D.__file__, D.EX.__file__, D.B.__file__])
    if a.workers > 1:
        from multiprocessing import Pool
        with Pool(min(a.workers, len(pairs))) as pool:
            res = pool.map(job, pairs)
    else:
        res = [job(p) for p in pairs]
    L = ["F6.8i-dev D2c revision 2.2 check: plain window vs bracketed window of the needle "
         "(development, decides nothing)", ""]
    for r in res:
        L.append(f"eps {r['eps']:g}  p_c {r['p_c']:.7f}")
        for key in ("plain", "bracket"):
            n = r[key]
            if "error" in n:
                L.append(f"  {key:8s} refused: {n['error']}")
            else:
                L.append(f"  {key:8s} N' " + ", ".join(f"{v:.6f}" for v in n["N"]) +
                         "  slopes " + ", ".join(f"{v:.4f}" for v in n["slopes"]) + f"  ({n['seconds']:.0f} s)")
        if "error" not in r["plain"] and "error" not in r["bracket"]:
            d = max(abs(x / y - 1) for x, y in zip(r["bracket"]["N"], r["plain"]["N"]))
            L.append(f"  max relative difference of N': {d:.2e}")
    L.append(f"total {time.time() - t0:.0f} s")
    print("\n".join(L))
    stem = C.save("f68i_d2c_bracket", vars(a), code, res, L, t0)
    print(f"saved {stem}.json")


if __name__ == "__main__":
    main()
