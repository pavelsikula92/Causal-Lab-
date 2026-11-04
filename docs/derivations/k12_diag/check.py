"""Validate rows_exact: apply to a smooth Gaussian u = exp(-r^2/2w^2) on a fine grid; compare with the
symbol prediction (Hankel transform with B(k)) and with lap u."""
import numpy as np, math
from scipy import integrate, special
from lin import rows_exact, hat
from kern import symbol_ratio
w = 8.0; R = 120.0
nodes = np.unique(np.r_[np.arange(0, R, 0.25), R])
rows = np.array([0.5, 4.0, 8.0, 12.0, 20.0])
u = np.exp(-nodes**2/(2*w*w))
A = rows_exact(nodes, rows, R)
Bu = A @ u
lap = lambda r: np.exp(-r**2/(2*w*w))*(r*r/w**4 - 2/w**2)
# symbol prediction: u_hat(k) = 2 pi w^2 exp(-k^2 w^2/2); B u(r) = int k dk/2pi J0(kr) (-k^2 r(k)) u_hat
ks = np.linspace(1e-4, 1.2, 3000); rk = np.array([symbol_ratio(k) for k in ks])
uh = 2*math.pi*w*w*np.exp(-ks**2*w*w/2)
for r, b in zip(rows, Bu):
    pred = np.trapezoid(ks/(2*math.pi)*special.j0(ks*r)*(-ks**2*rk)*uh, ks)
    print(f"r {r}: rows {b:.6e}  symbol {pred:.6e}  lap {lap(r):.6e}")
