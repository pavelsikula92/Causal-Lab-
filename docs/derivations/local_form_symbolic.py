"""
Symbolic derivation check: the local form of the causal-set field equation
==========================================================================

What this script verifies (sympy, exact algebra, no numerics)
-------------------------------------------------------------
The averaged causal-set (Benincasa-Dowker) operator tends locally to
box_g - R/2 (Belenchia et al. 2016). Acting on a static field Pi(x) in the
static metric family

    ds^2 = -Pi^2 dt^2 + Pi^(-2 gamma) dx^2          (signature - + + +)

(gamma = -1: the conformal geometry Pi^2 eta of equation (a);
 gamma = +1: hypothesis H1; in the weak field gamma is the PPN gamma)

the script computes box_g Pi and the Ricci scalar R from the metric
(Christoffel symbols -> Riemann -> Ricci) and checks, in 2+1 (n = 2) and
3+1 (n = 3) spatial dimensions:

  C0  convention: R of the conformal metric Pi^2 eta equals the textbook
      formula (Wald, Appendix D) for g~ = Omega^2 g with Omega = Pi,
          R~ = -2(D-1) Pi^-3 lap Pi - (D-1)(D-4) Pi^-4 |grad Pi|^2 ;
  C1  conformal case (gamma = -1):
          (box_g - R/2) Pi = l1 Pi^-2 lap Pi + l2 Pi^-3 |grad Pi|^2,
          l1 = D,  l2 = (D-2) + (D-1)(D-4)/2   (2+1: l2 = 0, 3+1: l2 = 2);
  C2  the equivalent power form: with p = 1 + l2/l1,
          lap(Pi^p) = (p/l1) Pi^(1+p) [l1 Pi^-2 lap Pi + l2 Pi^-3 |grad Pi|^2]
      identically, so  l1 Pi^-2 lap Pi + l2 Pi^-3 |grad Pi|^2 = kappa rho_s
      <=>  lap(Pi^p) = (p/l1) kappa Pi^(1+p) rho_s;
  C3  the gamma family (gamma symbolic):
          (box_g - R/2) Pi = Pi^(2 gamma) [ l1(gamma) lap Pi + l2(gamma) Pi^-1 |grad Pi|^2 ],
          l1(gamma) = 2 - gamma (n - 1),
          l2(gamma) = 1 + gamma      (2+1),   1 + gamma^2   (3+1);
  C4  weak field: with Pi = 1 + e u, the term linear in e is l1(gamma) lap u
      (the weak-field coupling of B_g Pi = kappa rho_s is kappa / l1);
  C5  3+1, gamma = 1 (H1): (box_g - R/2) Pi = 2 Pi |grad Pi|^2 exactly (no
      lap Pi term), so outside the source grad Pi = 0.

Every check prints OK or FAILS; the script exits with status 1 if any fails.
Run:  python docs/derivations/local_form_symbolic.py      (needs: pip install sympy)
Runtime: a few seconds.
"""

from __future__ import annotations

import sys

import sympy as sp


def metric_objects(n: int, gamma):
    """Coordinates, the field Pi(x), the metric diag(-Pi^2, Pi^(-2 gamma) ...), its inverse and sqrt(-det g)."""
    t = sp.Symbol("t")
    xs = sp.symbols(f"x0:{n}")
    Pi = sp.Function("Pi")(*xs)                      # static field: depends on space only
    coords = [t] + list(xs)
    g = sp.diag(*([-Pi ** 2] + [Pi ** (-2 * gamma)] * n))
    return coords, xs, Pi, g, g.inv(), sp.sqrt(-g.det())


def ricci_scalar(coords, g, ginv):
    """R = g^{bc} R_{bc} for a DIAGONAL metric, R_{bc} = d_a Gamma^a_bc - d_c Gamma^a_ba + Gamma^a_ad Gamma^d_bc - Gamma^a_cd Gamma^d_ba."""
    D = len(coords)
    Gam = [[[sum(ginv[a, d] * (sp.diff(g[d, b], coords[c]) + sp.diff(g[d, c], coords[b])
                                - sp.diff(g[b, c], coords[d])) for d in range(D)) / 2
             for c in range(D)] for b in range(D)] for a in range(D)]

    def ric(b, c):
        return sum(sp.diff(Gam[a][b][c], coords[a]) - sp.diff(Gam[a][b][a], coords[c])
                   + sum(Gam[a][a][d] * Gam[d][b][c] - Gam[a][c][d] * Gam[d][b][a] for d in range(D))
                   for a in range(D))
    return sum(ginv[b, b] * ric(b, b) for b in range(D))      # diagonal metric: only g^{bb} R_bb


def box(coords, ginv, sqrtg, f):
    """box_g f = (1/sqrt(-g)) d_mu (sqrt(-g) g^{mu nu} d_nu f) for a diagonal metric."""
    return sum(sp.diff(sqrtg * ginv[i, i] * sp.diff(f, coords[i]), coords[i]) for i in range(len(coords))) / sqrtg


def lap_grad2(xs, f):
    """Flat spatial Laplacian and squared gradient of f."""
    return sum(sp.diff(f, x, 2) for x in xs), sum(sp.diff(f, x) ** 2 for x in xs)


def is_zero(expr) -> bool:
    """Exact zero test (simplify after expanding powers)."""
    return sp.simplify(sp.powsimp(sp.expand(expr), force=True)) == 0


def main() -> int:
    results = []
    gam = sp.Symbol("gamma")
    for n in (2, 3):
        D = n + 1
        tag = f"{n}+1"
        # ---- conformal geometry Pi^2 eta (gamma = -1) -------------------------------------------------------------
        coords, xs, Pi, g, ginv, sqrtg = metric_objects(n, -1)
        R = ricci_scalar(coords, g, ginv)
        L, G2 = lap_grad2(xs, Pi)
        wald = -2 * (D - 1) * Pi ** -3 * L - (D - 1) * (D - 4) * Pi ** -4 * G2
        results.append((f"C0 {tag} Ricci scalar of Pi^2 eta = Wald conformal formula", is_zero(R - wald)))
        op = box(coords, ginv, sqrtg, Pi) - R / 2 * Pi
        l1, l2 = D, (D - 2) + sp.Rational((D - 1) * (D - 4), 2)
        results.append((f"C1 {tag} (box_g - R/2) Pi = {l1} Pi^-2 lap Pi + {l2} Pi^-3 |grad Pi|^2",
                        is_zero(op - (l1 * Pi ** -2 * L + l2 * Pi ** -3 * G2))))
        p = 1 + sp.Rational(l2, l1)
        Lp, _ = lap_grad2(xs, Pi ** p)
        results.append((f"C2 {tag} lap(Pi^p) = (p/l1) Pi^(1+p) [local form], p = {p}",
                        is_zero(Lp - p / l1 * Pi ** (1 + p) * (l1 * Pi ** -2 * L + l2 * Pi ** -3 * G2))))
        # ---- the gamma family, gamma symbolic ----------------------------------------------------------------------
        coords, xs, Pi, g, ginv, sqrtg = metric_objects(n, gam)
        R = ricci_scalar(coords, g, ginv)
        op = box(coords, ginv, sqrtg, Pi) - R / 2 * Pi
        L, G2 = lap_grad2(xs, Pi)
        l1g = 2 - gam * (n - 1)
        l2g = 1 + gam if n == 2 else 1 + gam ** 2
        results.append((f"C3 {tag} (box_g - R/2) Pi = Pi^(2g) [({l1g}) lap Pi + ({l2g}) Pi^-1 |grad Pi|^2]",
                        is_zero(op - Pi ** (2 * gam) * (l1g * L + l2g * Pi ** -1 * G2))))
        # weak field: Pi = 1 + e u, linear coefficient
        e = sp.Symbol("e")
        u = sp.Function("u")(*xs)
        weak = op.subs(Pi, 1 + e * u).doit()
        lin = sp.diff(weak, e).subs(e, 0)
        Lu, _ = lap_grad2(xs, u)
        results.append((f"C4 {tag} weak field: linear term = ({l1g}) lap u", is_zero(lin - l1g * Lu)))
        if n == 3:
            op1 = sp.simplify(op.subs(gam, 1))
            results.append(("C5 3+1 gamma = 1: (box_g - R/2) Pi = 2 Pi |grad Pi|^2 (no lap Pi term)",
                            is_zero(op1 - 2 * Pi * G2)))
    for name, ok in results:
        print(f"{'OK   ' if ok else 'FAILS'}  {name}")
    print(f"{sum(ok for _, ok in results)}/{len(results)} checks passed")
    return 0 if all(ok for _, ok in results) else 1


if __name__ == "__main__":
    sys.exit(main())
