"""Exact static kernel of the expected 2+1 BD operator (eps*rho = 1), its symbol, and checks."""
import numpy as np, math
from scipy import integrate, special
LAM = math.pi/12.0                      # w = eps rho V, V = pi tau^3/12
C = (-1.0, 27/8, -9/4)
PREF = -(math.pi/(3*math.sqrt(2)))**(2/3)/math.gamma(5/3)
def g(w): return (C[0] + C[1]*w + C[2]*w*w/2)*np.exp(-w)
TAU_MAX = (80/LAM)**(1/3)
def K(rho):
    """static kernel K(rho) = int_0^inf tau g(lam tau^3)/sqrt(tau^2+rho^2) dtau"""
    f = lambda t: t*g(LAM*t**3)/math.sqrt(t*t+rho*rho)
    return integrate.quad(f, 0, TAU_MAX, limit=400, epsabs=1e-15, epsrel=1e-12, points=[1,2,3,4])[0]
def Mm(m): return integrate.quad(lambda t: t**(1+m)*g(LAM*t**3), 0, TAU_MAX, limit=400, epsabs=1e-15)[0]
def symbol_ratio(k):
    """B(k)/(-k^2): 1 for the local operator"""
    I = integrate.quad(lambda t: t*g(LAM*t**3)*math.exp(-k*t), 0, TAU_MAX, limit=400, epsabs=1e-16, epsrel=1e-13)[0]
    return PREF*(1 + 2*math.pi/k*I)/(-k*k)
if __name__ == "__main__":
    print("M0..M4", [Mm(m) for m in range(5)], "1/(2pi)", 1/(2*math.pi))
    for k in [0.01,0.03,0.1,0.2,0.3,0.5,1.0]:
        print(f"k {k}: B/(-k^2) = {symbol_ratio(k):.5f}")
    for r in [0,0.5,1,2,4,8,16,32]:
        print(f"K({r}) = {K(r):.4e}")
