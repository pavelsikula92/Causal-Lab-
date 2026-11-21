"""
OV3 section 10 (development, decides nothing): the PPN parameter beta from the local cubic static action
=========================================================================================================

Static isotropic metric  g = -e^{2 omega} dt^2 + e^{2 psi} dx^2   (psi = omega + nu in the paper's variables
g = e^{2 omega} [-dt^2 + e^{2 nu} dx^2]). The most general local static action with two derivatives, to cubic order:

    L = q1 w'^2 + q2 w'p' + q3 p'^2                                   (quadratic; EH: 0, 4, 2)
      + c1 w w'^2 + c2 w w'p' + c3 w p'^2 + c4 p w'^2 + c5 p w'p' + c6 p p'^2    (cubic; every two-derivative
                                                                                    cubic term reduces to these)

(w = omega, p = psi; in 3D, ' -> grad, products -> dot products). The point mass: first order w1 = -m/r, p1 = +m/r
(omega + psi has no source and is regular, so p1 = -w1). Second order, outside the source: w2 = a m^2/r^2,
p2 = b m^2/r^2. With -g_00 = e^{2w} = 1 - 2U + 2 beta U^2 and U = m/r, beta = 1 + a.

Checks printed:
  1. the general formula beta(c1..c6) with the EH quadratic part;
  2. EH (c = 0, 4, 2, 0, 4, 2, from L = e^{w + p} [4 w'p' + 2 p'^2]) gives beta = 1 and b = -1/4 (Schwarzschild in
     isotropic coordinates: psi = 2 ln(1 + m/2r) = m/r - m^2/(4 r^2));
  3. beta depends on the cubic coefficients only through kappa3 = L3(w = -U, p = U) / (U U'^2); beta - 1 = -kappa3/8;
  4. along the direction (w, p) = (-U, U), i.e. (omega, nu) = (-U, 2U), the full EH action is -2 U'^2 at every order;
  5. the planar two-mode normalisation: u = a cos kz + c cos 2kz, <u u'^2> = (3/4) a^2 c k^2 + ..., <u'^2> =
     (a^2 + 4 c^2) k^2 / 2; and the pure conformal check: EH cubic of the conformal sector is 12 u u'^2.

    python -m diag.ov3_beta_sym
"""
import sympy as sp


def main():
    r, m = sp.symbols("r m", positive=True)
    c = sp.symbols("c1:7")
    a, b, eps = sp.symbols("a b eps")
    w = sp.Function("w")(r)
    p = sp.Function("p")(r)
    wd, pd = w.diff(r), p.diff(r)
    L2 = 4 * wd * pd + 2 * pd ** 2
    L3 = c[0] * w * wd ** 2 + c[1] * w * wd * pd + c[2] * w * pd ** 2 + c[3] * p * wd ** 2 + c[4] * p * wd * pd \
        + c[5] * p * pd ** 2
    lag = r ** 2 * (L2 + L3)

    def el(f):
        return sp.diff(lag, f) - sp.diff(sp.diff(lag, f.diff(r)), r)

    sub = {w: -eps * m / r + eps ** 2 * a * m ** 2 / r ** 2, p: eps * m / r + eps ** 2 * b * m ** 2 / r ** 2}
    eqs = [sp.simplify(sp.expand(el(f).subs(sub).doit()).coeff(eps, 2) * r ** 2 / m ** 2) for f in (w, p)]
    sol = sp.solve(eqs, [a, b], dict=True)[0]
    beta = sp.simplify(1 + sol[a])
    print("1. second-order equations (outside the source):", eqs)
    print("   a =", sol[a], ";  b =", sol[b])
    print("   beta =", beta)
    eh = dict(zip(c, (0, 4, 2, 0, 4, 2)))
    print("2. EH: beta =", beta.subs(eh), "; b =", sol[b].subs(eh), "(Schwarzschild isotropic: -1/4)")
    U = sp.Function("U")(r)
    k3 = sp.simplify(L3.subs({w: -U, p: U}).doit() / (U * U.diff(r) ** 2))
    print("3. kappa3 = L3(w=-U, p=U)/(U U'^2) =", k3, ";  1 - kappa3/8 - beta =", sp.simplify(1 - k3 / 8 - beta))
    W, P = sp.symbols("W P")
    full = sp.exp(W + P) * (4 * wd * pd + 2 * pd ** 2)
    along = sp.simplify(full.subs({W: -U, P: U}).subs({w: -U, p: U}).doit())
    print("4. EH along (omega, psi) = (-U, U):", along)
    z, k, A, C = sp.symbols("z k A C", positive=True)
    u = A * sp.cos(k * z) + C * sp.cos(2 * k * z)
    per = 2 * sp.pi / k
    avg = lambda f: sp.simplify(sp.integrate(sp.expand(sp.expand_trig(f)), (z, 0, per)) / per)
    cub = sp.expand(avg(u * u.diff(z) ** 2))
    quad = sp.expand(avg(u.diff(z) ** 2))
    print("5. <u u'^2> =", cub, ";  <u'^2> =", quad)
    print("   single mode (C = 0): <u u'^2> =", cub.subs(C, 0), "(the exact null: a half-period shift maps A -> -A)")


if __name__ == "__main__":
    main()
