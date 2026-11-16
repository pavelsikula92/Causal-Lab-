"""F6.8f-dev step 1a: quadrature levels of the exact volume, rows 2 and 4, fixed tau (development). printed; saved by hand to diag/out/f68f_1a_levels.txt"""
import sys, math, numpy as np
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
import f68f_1a as A
from multiprocessing import Pool
rows = [A.ROWS[1], A.ROWS[3]]
levs = [(24, 16, 24, None), (32, 24, 32, None)]
ths = [2 * math.pi * j / 8 for j in range(8)]
jobs = []
for r in rows:
    for lev in levs:
        jobs.append(r + (0.0, 0.0, "tau", lev))
        for b in (0.01, 0.02):
            for th in ths:
                jobs.append(r + (b, th, "tau", lev))
if __name__ == "__main__":
    with Pool(2) as p:
        res = p.map(A.task, jobs, chunksize=1)
    for r in rows:
        for lev in levs:
            V0 = [V for (a, V, T) in res if a[:4] == r and a[7] == lev and a[4] == 0][0]
            R = {}
            for b in (0.01, 0.02):
                Vs = np.array([V for (a, V, T) in res if a[:4] == r and a[7] == lev and a[4] == b])
                R[b] = 2 * (Vs.mean() / V0 - 1) / b ** 2
            g = R[0.01] - (R[0.02] - R[0.01]) / 3
            print(r, lev[:3], "R_tau(0.01) %+.8e R_tau(0.02) %+.8e extrapolated %+.8e" % (R[0.01], R[0.02], g), flush=True)
