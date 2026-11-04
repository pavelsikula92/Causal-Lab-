"""
Linear (weak-field) K.12 problem with the EXACT static kernel of the expected flat 2+1 BD operator:
    3 B[u](r_i) = s(r_i),  u(R) = 0 (u = 0 beyond the wall; exact in flat geometry thanks to the vacuum subtraction),
u piecewise linear on nodes, collocation at 2 rows per element, least squares (the K.12 scheme).
Compared with the continuum 3 lap u = s with u(R) = 0.
"""
import numpy as np, math, sys
from scipy import integrate
from kern import K, PREF, symbol_ratio

# --- tabulated static kernel (lengths in units with eps*rho = 1, as in the lab)
RHO_TAB = np.r_[np.linspace(0, 12, 2401), np.geomspace(12.05, 4000, 400)]
K_TAB = np.array([K(r) for r in RHO_TAB[RHO_TAB <= 60]])
far = RHO_TAB[RHO_TAB > 60]
A4 = K_TAB[-1]*60**5                                     # tail K ~ A4 / rho^5 (checked: M4 term)
K_TAB = np.r_[K_TAB, A4/far**5]
def Kf(rho): return np.interp(rho, RHO_TAB, K_TAB, right=0.0)

def Kring(r, rp):
    """int_0^{2pi} K(|x-y|) dphi for |x| = r, |y| = rp (vectorised over rp); phi grid resolves the kernel scale
    (arc spacing <= 0.05 l_k for rho < 30) at every radius"""
    a = min(math.pi, 30.0/max(r, 1e-9))
    ph = np.unique(np.r_[np.linspace(0, a, 1201), np.linspace(0, math.pi, 801)])
    w = np.r_[0.5*(ph[1]-ph[0]), 0.5*(ph[2:]-ph[:-2]), 0.5*(ph[-1]-ph[-2])]*2.0
    rp = np.atleast_1d(rp)[:, None]
    rho = np.sqrt(np.maximum(r*r + rp*rp - 2*r*rp*np.cos(ph)[None, :], 0.0))
    return (Kf(rho)*w[None, :]).sum(axis=1)

def rp_grid(r, R, Rout):
    """integration grid in r': dense near r (kernel scale ~1), coarser away; up to Rout"""
    pts = np.unique(np.r_[np.linspace(max(0, r-8), min(Rout, r+8), 1601),
                          np.linspace(0, Rout, 4001), np.geomspace(1e-3, Rout, 600), [0.0, R, Rout]])
    return pts[pts <= Rout]

def hat(nodes, x):
    """matrix H[j, k] = hat_k(x_j)"""
    H = np.zeros((len(x), len(nodes)))
    for k in range(len(nodes)):
        e = np.zeros(len(nodes)); e[k] = 1
        H[:, k] = np.interp(x, nodes, e, right=0.0)
    return H

def rows_exact(nodes, rows, R):
    A = np.zeros((len(rows), len(nodes)))
    for i, r in enumerate(rows):
        rp = rp_grid(r, R, R)
        kr = Kring(r, rp)*rp
        H = hat(nodes, rp); h0 = hat(nodes, np.array([r]))[0]
        w = np.gradient(rp)  # trapezoid-like weights
        w = np.r_[0.5*(rp[1]-rp[0]), 0.5*(rp[2:]-rp[:-2]), 0.5*(rp[-1]-rp[-2])]
        inner = (w*kr) @ (H - h0[None, :])
        # outside the disc: u = 0 there, contributes -u(x) int_out K
        ro = np.r_[np.linspace(R, R+60, 3001), np.geomspace(R+60.01, R+4000, 800)]
        ko = Kring(r, ro)*ro
        Kout = np.trapezoid(ko, ro)
        A[i] = PREF*(inner - h0*Kout)
    return A

def rows_local(nodes, rows):
    """collocation of the LOCAL operator lap on the piecewise-linear field at interior points: f'/r (f''=0 there)"""
    A = np.zeros((len(rows), len(nodes)))
    for i, r in enumerate(rows):
        k = np.searchsorted(nodes, r) - 1
        h = nodes[k+1]-nodes[k]
        A[i, k] = -1/(h*r); A[i, k+1] = 1/(h*r)
    return A

def shape_nodes(R, a, dr, ratio):
    inner = np.arange(0.0, max(4.0, 2.0*a)+1e-9, dr)
    geo = [float(inner[-1])]
    while geo[-1]*ratio < R: geo.append(geo[-1]*ratio)
    return np.unique(np.r_[inner, geo[1:], R])

def continuum(sig, q, R, n=200001):
    r = np.linspace(0, R, n)
    s = q*np.exp(-r**2/(2*sig**2))
    Q = integrate.cumulative_trapezoid(2*math.pi*r*s, r, initial=0)
    du = np.zeros_like(r); du[1:] = Q[1:]/(2*math.pi*r[1:])/3.0
    u = integrate.cumulative_trapezoid(du, r, initial=0); u -= u[-1]
    return r, u

def solve(A, nodes, rows, sig, q):
    s = q*np.exp(-rows**2/(2*sig**2))
    sol, *_ = np.linalg.lstsq(3*A[:, :-1], s, rcond=None)
    return np.r_[sol, 0.0]

def flux(r, u, x, h):
    f = lambda y: np.interp(y, r, u)
    return x*(f(x+h)-f(x-h))/(2*h)

if __name__ == "__main__":
    sig = float(sys.argv[1]); mode = sys.argv[2]       # mode: k12 | fine
    q = 0.064/sig**2; R = 12.5*sig
    if mode == "k12":
        nodes = shape_nodes(R, 2*sig, sig/4, 1.12)
    else:
        step = float(sys.argv[3])
        nodes = np.unique(np.r_[np.arange(0, R, step), R])
    per = 2; fr = (np.arange(per)+0.5)/per
    rows = np.sort((nodes[:-1, None] + fr[None, :]*np.diff(nodes)[:, None]).ravel())
    rc, uc = continuum(sig, q, R)
    out = {}
    for name, A in (("exact", rows_exact(nodes, rows, R)), ("local", rows_local(nodes, rows))):
        u = solve(A, nodes, rows, sig, q)
        xs = np.array([1.25, 2.5, 3.75, 3.0])*sig
        h = sig/16
        # field of the scheme is piecewise linear: evaluate its flux on element midpoints nearest to xs via interp
        fl = flux(nodes, u, xs, h); flc = flux(rc, uc, xs, h)
        xf = np.array([0.25, 0.5, 0.75, 0.9])*R
        ff = flux(nodes, u, xf, h); ffc = flux(rc, uc, xf, h)
        print(f"sigma {sig:g} {mode} nodes {len(nodes)} op {name}: u(0) {u[0]:.5f} cont {uc[0]:.5f} ratio {u[0]/uc[0]:.4f}; "
              f"slope ratio 1.25/2.5/3.75 sig: " + " ".join(f"{a/b:.4f}" for a, b in zip(fl[:3], flc[:3]))
              + f"; flux(3sig) ratio {fl[3]/flc[3]:.4f}; flux ratio at R/4,R/2,3R/4,0.9R: " + " ".join(f"{a/b:.4f}" for a, b in zip(ff, ffc)))
