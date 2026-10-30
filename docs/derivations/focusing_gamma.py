"""
Focusing of light by matter fixes the spatial part of the metric (23 Oct 2026, docs/derivations/spatial_part.md).

Weak static field, signature (-, +, +, +): g = diag(-(1 + 2 phi), (1 - 2 psi) delta_ij), phi (clocks) and psi
(space) independent unknown functions. The linearised Ricci tensor gives for a null direction k = (1, n), |n| = 1:
    R_kk = lap phi + lap psi + n^i n^j d_i d_j (psi - phi).
The focusing principle (Raychaudhuri / Clausius on every local horizon): R_kk = 8 pi G T_kk for every null k. For
matter at rest T_kk = rho (k^0)^2 = rho, the same in every direction. The script prints R_kk minus its value for
psi = phi and checks that the direction-dependent part is exactly n^i n^j d_i d_j (psi - phi).
Usage: python docs/derivations/focusing_gamma.py
"""
import sympy as sp

t, x, y, z, th, ph = sp.symbols('t x y z theta varphi', real=True)
X = [t, x, y, z]
phi = sp.Function('phi')(x, y, z)                       # temporal potential (clock rate)
psi = sp.Function('psi')(x, y, z)                       # spatial potential (rulers, width of the cones)
eta = sp.diag(-1, 1, 1, 1)
h = sp.diag(-2 * phi, -2 * psi, -2 * psi, -2 * psi)
ei = eta.inv()
hup = ei * h
tr = sum(hup[i, i] for i in range(4))
box = lambda F: sum(ei[i, i] * sp.diff(F, X[i], 2) for i in range(4))
R = sp.Matrix(4, 4, lambda m, n: sp.Rational(1, 2) * (
    sum(sp.diff(hup[s, n], X[s], X[m]) + sp.diff(hup[s, m], X[s], X[n]) for s in range(4))
    - box(h[m, n]) - sp.diff(tr, X[m], X[n])))
nv = [sp.sin(th) * sp.cos(ph), sp.sin(th) * sp.sin(ph), sp.cos(th)]
k = sp.Matrix([1] + nv)
Rkk = (k.T * R * k)[0]
lap = lambda F: sum(sp.diff(F, c, 2) for c in (x, y, z))
hess = lambda F: sum(nv[i] * nv[j] * sp.diff(F, (x, y, z)[i], (x, y, z)[j]) for i in range(3) for j in range(3))
expected = lap(phi) + lap(psi) + hess(psi - phi)
print("R_kk - [lap phi + lap psi + n n : dd (psi - phi)] =", sp.simplify(sp.expand(Rkk - expected)))
print("Focusing equal in every direction for matter at rest requires the traceless Hessian of psi - phi to vanish;")
print("with psi, phi -> 0 far away this gives psi = phi (PPN gamma = 1), and then 2 lap phi = 8 pi G rho (Newton).")
