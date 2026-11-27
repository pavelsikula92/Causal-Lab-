"""N'(q) at the cusp (p_c, theta_c = -p_c/2) by the exact needle of f68i_d2d (per-decade windows), q -> 0."""
import math, sys, numpy as np
sys.path.insert(0, "/home/claude/gh/Causal-Lab-")
import diag.f68i_d2d as D
from scipy.integrate import quad
from scipy.optimize import curve_fit
BQ = 2 * quad(lambda x: (1 - x ** 4) ** 1.5, 0, 1, epsabs=1e-15, epsrel=1e-15)[0]
qs = np.array([0.0625, 0.044, 0.031, 0.0221, 0.0156, 0.011, 0.0078, 0.0055, 0.0039, 0.00276])
for eps, pc in ((float(a), float(b)) for a, b in (x.split(":") for x in sys.argv[1:])):
    cs = qs ** 2 / 2
    N = D.Needle(eps, pc, (-pc / 2) % (2 * math.pi), float(cs.max()))
    Np = N.A(cs) * cs ** 0.25 / (4 ** 0.25 * BQ)
    p, _ = curve_fit(lambda q, N0, a, b: N0 + a * q ** b, qs, Np, p0=(Np[-1], 0.05, 1.0), maxfev=20000)
    print(f"eps {eps:g} p_c {pc}: N'(q) " + " ".join(f"{x:.6f}" for x in Np) + f" | fit N0 {p[0]:.6f} a {p[1]:.4f} b {p[2]:.3f}; "
          f"N'(0.0625)/N0 {Np[0]/p[0]:.5f}, Psi factor (N'(0.0625)/N0)^(9/7) {(Np[0]/p[0])**(9/7):.5f}; nray {N.nray}", flush=True)
