"""OV2 sec. 14, L1 check (development, decides nothing): the expected crossing-pair sum of sec. 11 (F-a as first written)
vanishes; option --calib computes the kinetic action per unit time instead (expect -4m).

L1 check: expected crossing-pair sum X = eps^2 sum_{pairs crossing t0} f(n) (dphi)^2 for a static kink,
against the kink energy H = int (1/2 phi'^2 + V) dz = 8 m. Light-cone pair coordinates p = du, q = dv (>0),
dt = (p+q)/2, dz = (q-p)/2; density: N per unit (u,v) area (proper 2N), w = eps N p q, <f> = g(w).
Crossing measure for fixed separation: the first point lies in a strip of time width dt; du dv = 2 dt dz, so
X = eps^2 N^2 int dp dq g(w) * 2 * dt * Psi(dz),  Psi(D) = int dz (phi(z+D)-phi(z))^2."""
import math, sys
import numpy as np

def phi(z, m): return 4*np.arctan(np.exp(m*z))

def Psi_tab(m, Dmax, n=4001):
    D = np.linspace(-Dmax, Dmax, n)
    z = np.linspace(-40/m - Dmax, 40/m + Dmax, 200001); dz = z[1]-z[0]
    out = np.array([np.sum((phi(z+d, m)-phi(z, m))**2)*dz for d in D])
    return D, out

def X_of(m, eps, N, L):
    D, Ps = Psi_tab(m, L/2)
    g = lambda w: np.exp(-w)*(1-2*w+0.5*w*w)
    # variables: t = ln p in [ln pmin, ln L], w in [0, min(60, eps N p L)]
    x, wx = np.polynomial.legendre.leggauss(64)
    tot = 0.0
    pmin = 1e-12/(eps*N)
    edges = np.linspace(math.log(pmin), math.log(L), 161)
    for a, b in zip(edges[:-1], edges[1:]):
        t = 0.5*(a+b) + 0.5*(b-a)*x; wt = 0.5*(b-a)*wx
        for ti, wti in zip(t, wt):
            p = math.exp(ti)
            wmax = min(60.0, eps*N*p*L)
            # w panels
            we = np.linspace(0, wmax, 25)
            for c, d in zip(we[:-1], we[1:]):
                w = 0.5*(c+d)+0.5*(d-c)*x; ww = 0.5*(d-c)*wx
                q = w/(eps*N*p)
                dq = ww/(eps*N*p)
                dt = 0.5*(p+q); dzz = 0.5*(q-p)
                ps = np.interp(dzz, D, Ps)
                tot += wti*p*np.sum(g(w)*2*(1.0 if CALIB else dt)*ps*dq)
    return eps**2*N**2*tot

CALIB = "--calib" in sys.argv
args = [a for a in sys.argv[1:] if not a.startswith("--")]
m = float(args[0]) if args else 9.0
for N in [16640, 66560, 266240]:
    X = X_of(m, 0.12, N, 1.0)
    if CALIB:
        print(f"m {m} N {N}: S_kin per unit t = {-X:+.4f} (expect {-4*m:+.1f})", flush=True)
    else:
        print(f"m {m} N {N}: X = {X:+.4f}   H = 8m = {8*m:.1f}   X/H = {X/(8*m):+.4f}", flush=True)
