"""F6.8f-dev step 0': p (direction) quadrature of the F6.8e3 bulk slope at k = 0.12, n = 24, 48 nodes (development)."""
import sys, os, math, time, numpy as np
_h = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, _h); sys.path.insert(0, os.path.dirname(_h))
import f68f_layer as F
from multiprocessing import Pool
k = 0.12
def bulk(args):
    n, p = args
    T = n / k; CO = p / n
    s, ws = F._s_nodes(np.array([1e12]), F.N_S, True); s, ws = s[0], ws[0]; keep = ws > 0; s, ws = s[keep], ws[keep]
    Fb = F.integrand(k, CO, T, s)
    return float((ws * T * T / (2.0 * F._L2 * np.sqrt(T * T + s / F._L2)) * Fb).sum())
if __name__ == "__main__":
    for n in (24.0,):
        for panels in ([[0, n/4, 12], [n/4, n/2, 12], [n/2, 3*n/4, 12], [3*n/4, n, 12]],):
            pn, wp = F._p_nodes(panels)
            t = time.time()
            with Pool(2) as pool: vals = np.array(pool.map(bulk, [(n, float(p)) for p in pn]))
            rho = 2 * math.pi / n * 2 * float((wp * vals).sum())
            print(n, len(pn), -F._A4 * rho / k ** 5, f"{time.time()-t:.0f}s", flush=True)
