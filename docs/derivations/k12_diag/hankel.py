"""Exact expected flat 2+1 BD operator, NO grid and NO wall (infinite plane), linear static response
3 B[u] = s to the Gaussian source: flux r u'(r) relative to the continuum (Gauss law)."""
import numpy as np, math
from scipy import special
from kern import symbol_ratio
ks = np.r_[np.geomspace(1e-6, 1e-2, 400), np.linspace(1.0001e-2, 3.0, 6000)]
rk = np.array([symbol_ratio(k) for k in ks])
for sig in (8, 16, 32, 64, 128):
    sh = np.exp(-ks**2*sig**2/2)              # source shape (charge normalised out)
    out = []
    for x in np.array([1.25, 2.5, 3.0, 3.75, 6.25, 9.0])*sig:
        # flux r u' = r * int k dk/2pi * (s_hat/(3 k^2 r(k))) * k J1(kr); ratio to r(k) = 1
        num = np.trapezoid(sh*special.j1(ks*x)/rk, ks); den = np.trapezoid(sh*special.j1(ks*x), ks)
        out.append(num/den)
    print(f"sigma {sig}: flux ratio at r/sigma = 1.25, 2.5, 3, 3.75, 6.25, 9: " + " ".join(f"{v:.4f}" for v in out))
