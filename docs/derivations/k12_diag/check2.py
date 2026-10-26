"""B[r^2 - R^2] must be 4 far from the wall (quadratic: nonlocal corrections vanish); B[Gaussian] vs symbol."""
import numpy as np, math
from lin import rows_exact
R = 400.0
nodes = np.unique(np.r_[np.arange(0, R, 0.5), R])
rows = np.array([5.0, 100.0, 300.0, 380.0, 395.0])
A = rows_exact(nodes, rows, R)
print("B[r^2-R^2] at", rows, "=", A @ (nodes**2 - R*R))
