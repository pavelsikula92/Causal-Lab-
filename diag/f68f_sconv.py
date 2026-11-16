"""F6.8f-dev step 0': s quadrature of the bulk density (F6.8e3 instrument, full s range) at k = 0.12, n = 24, single
directions p; N_S = 16 (as F6.8e3), 24, 32, 48. Development. Output diag/out/f68f_sconv.txt"""
import math, os, sys, time
import numpy as np
_h = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, _h); sys.path.insert(0, os.path.dirname(_h))
import f68f_layer as F
from multiprocessing import Pool
k, n = 0.12, 24.0


def bulk(args):
    p, ns = args
    T = n / k
    s, ws = F._s_nodes(np.array([1e12]), ns, True)
    s, ws = s[0], ws[0]
    keep = ws > 0
    s, ws = s[keep], ws[keep]
    Fb = F.integrand(k, p / n, T, s)
    return p, ns, float((ws * T * T / (2.0 * F._L2 * np.sqrt(T * T + s / F._L2)) * Fb).sum()), int(s.size)


if __name__ == "__main__":
    t0 = time.time()
    jobs = [(p, ns) for p in (0.5, 3.0, 10.0) for ns in (16, 24, 32, 48)]
    with Pool(2) as pool:
        res = pool.map(bulk, jobs)
    lines = ["F6.8f-dev: s quadrature of the bulk pair density at r = T (k = 0.12, n = 24); development", ""]
    for p, ns, v, m in res:
        lines.append(f"  p {p:5.1f}  N_S {ns:3d} ({m} nodes): {v:+.10e}")
    lines.append(f"  cost {time.time() - t0:.0f} s")
    txt = "\n".join(lines)
    print(txt)
    open(os.path.join(_h, "out", "f68f_sconv.txt"), "w").write(txt + "\n")
