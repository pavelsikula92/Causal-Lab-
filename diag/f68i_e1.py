"""F6.8i-dev E1 (development, decides nothing): the smooth part of <s4hat A0^-2> beyond the first caustic, coarse grid.
Usage: python diag/f68i_e1.py <out.json>  (16 phases; eps 2-40; p 0.5, 2, 4, 8, 16). docs/derivations/f68i_dev.md section 3."""
import sys, time, json, numpy as np
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.dirname(__import__("os").path.abspath(__file__))))
from diag import f68f_1c as S
from multiprocessing import Pool
EPS = [2, 4, 8, 12, 16, 24, 40]
PS = [0.5, 2, 4, 8, 16]
def run(a):
    eps, p = a
    t0 = time.time()
    try:
        r = S.point(eps, p, 16)
        sh = r["s4hat"]; A0 = r["A"][:, 0]
        return {"eps": eps, "p": p, "phi": r["phi"], "A0": r["A0"], "nray": r["nray"], "minK": r["minK"], "lin": r["lin"],
                "spread": float(np.std(sh * A0 ** -2.0)), "t": time.time() - t0}
    except Exception as e:
        return {"eps": eps, "p": p, "err": repr(e)}
if __name__ == "__main__":
    jobs = [(e, p) for e in EPS for p in PS]
    with Pool(2) as pool:
        out = pool.map(run, jobs, chunksize=1)
    json.dump(out, open(sys.argv[1], "w"), indent=1)
    for r in out:
        print(r)
