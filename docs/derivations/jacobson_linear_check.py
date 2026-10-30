"""
Linear check of Jacobson's condition for the weak field of H1 (23 Oct 2026, docs/derivations/jacobson_route.md, section 6).

For local Rindler horizons, the Clausius relation dS = dQ/T with S proportional to the area is equivalent to
    R_ab k^a k^b = 8 pi G T_ab k^a k^b     for every null k
(Jacobson 1995). Here the weak field of H1 (gamma = 1) is written as
    h_mn = -2 U (eta_mn + 2 n_m n_n),
with n a constant unit timelike vector: the frame in which H1 is written. The source is dust T = rho u u moving
with velocity v; U is static in the source frame and 4 pi G rho = (source-frame Laplacian of U).
The script evaluates R_kk - 8 pi G T_kk (linearised Ricci tensor) at a point for four null directions, twice:
  n = laboratory time  (H1 as a scalar Pi in a fixed frame)          -> nonzero: the condition fails;
  n = u of the source  (H1 written with the source's own 4-velocity) -> zero: the condition holds.
Usage: python docs/derivations/jacobson_linear_check.py
"""
import sympy as sp

t, x, y, z, a = sp.symbols('t x y z a', real=True)
v = sp.Rational(3, 5)                                   # source speed (exact rational, any |v| < 1)
X = [t, x, y, z]
eta = sp.diag(-1, 1, 1, 1)
gm = 1 / sp.sqrt(1 - v ** 2)
u_up = sp.Matrix([gm, gm * v, 0, 0])                    # 4-velocity of the source
f = sp.exp(-(a ** 2 + 2 * y ** 2 + 3 * z ** 2)) * (1 + a * y)   # arbitrary smooth profile in the source frame
lapf = sp.diff(f, a, 2) + sp.diff(f, y, 2) + sp.diff(f, z, 2)    # = 4 pi G rho (source frame)
xi = gm * (x - v * t)                                   # comoving coordinate along the motion
U, lapU = f.subs(a, xi), lapf.subs(a, xi)


def ricci_lin(h):
    """Linearised Ricci tensor R_mn = 1/2 (d_s d_m h^s_n + d_s d_n h^s_m - box h_mn - d_m d_n h)."""
    ei = eta.inv()
    hup = ei * h
    tr = sum(hup[i, i] for i in range(4))
    box = lambda F: sum(ei[i, i] * sp.diff(F, X[i], 2) for i in range(4))
    return sp.Matrix(4, 4, lambda m, n: sp.Rational(1, 2) * (
        sum(sp.diff(hup[s, n], X[s], X[m]) + sp.diff(hup[s, m], X[s], X[n]) for s in range(4))
        - box(h[m, n]) - sp.diff(tr, X[m], X[n])))


def h_of(n_up):
    """Weak field of H1 written in the frame n: h_mn = -2 U (eta_mn + 2 n_m n_n)."""
    n_dn = eta * n_up
    return sp.Matrix(4, 4, lambda i, j: -2 * U * (eta[i, j] + 2 * n_dn[i] * n_dn[j]))


ks = [sp.Matrix([1, 1, 0, 0]), sp.Matrix([1, -1, 0, 0]), sp.Matrix([1, 0, 1, 0]),
      sp.Matrix([1, sp.Rational(3, 5), sp.Rational(4, 5), 0])]
pt = {t: sp.Rational(1, 3), x: sp.Rational(1, 2), y: sp.Rational(1, 5), z: sp.Rational(-1, 4)}
for name, n_up in [("H1 in the laboratory frame (n = laboratory time)", sp.Matrix([1, 0, 0, 0])),
                   ("H1 covariant (n = u of the source)", u_up)]:
    R = ricci_lin(h_of(n_up))
    out = []
    for k in ks:
        uk = (u_up.T * eta * k)[0]
        out.append(float(sp.N(((k.T * R * k)[0] - 2 * lapU * uk ** 2).subs(pt), 15)))   # R_kk - 8 pi G rho (u.k)^2
    print(f"{name}: R_kk - 8 pi G T_kk at a point = " + ", ".join(f"{o:+.2e}" for o in out))
