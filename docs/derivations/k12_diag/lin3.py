"""
Linear-order model of the K.12 scheme: the curved operator E[B_g(Pi)] Pi, Pi = 1 + u, splits at first order into
  E[B_flat] u_PL          (the field ACTED ON, piecewise linear on the nodes; local limit 1 x lap u)
+ geometry response      (the operator acting on 1 in the geometry u_spline, C^2 spline; local limit 2 x lap u).
Variants for the geometry part: 'lap' (exactly local on the spline) or 'B' (the flat nonlocal kernel on the spline,
a stand-in for the true nonlocal geometry response). Fields are measured on the spline of the nodal values, as K.12.
"""
import numpy as np, math, sys
from scipy.interpolate import CubicSpline
from lin import rows_exact, rows_local, shape_nodes, continuum, flux, Kring, rp_grid, PREF

def spline_basis(nodes):
    """cubic splines S_k (Pi'(0) = 0, not-a-knot at R) through the unit vectors"""
    return [CubicSpline(nodes, np.eye(len(nodes))[k], bc_type=((1, 0.0), "not-a-knot")) for k in range(len(nodes))]

def rows_geom(nodes, rows, R, kind):
    S = spline_basis(nodes)
    A = np.zeros((len(rows), len(nodes)))
    for i, r in enumerate(rows):
        if kind == "lap":
            A[i] = [s(r, 2) + s(r, 1)/r for s in S]
            continue
        rp = rp_grid(r, R, R)
        kr = Kring(r, rp)*rp
        w = np.r_[0.5*(rp[1]-rp[0]), 0.5*(rp[2:]-rp[:-2]), 0.5*(rp[-1]-rp[-2])]
        Sv = np.array([s(rp) for s in S]).T; S0 = np.array([s(r) for s in S])
        ro = np.r_[np.linspace(R, R+60, 3001), np.geomspace(R+60.01, R+4000, 800)]
        Kout = np.trapezoid(Kring(r, ro)*ro, ro)
        A[i] = PREF*((w*kr) @ (Sv - S0[None, :]) - S0*Kout)
    return A

def metrics(nodes, u, sig, R, q, tag):
    rc, uc = continuum(sig, q, R)
    f = CubicSpline(nodes, u, bc_type=((1, 0.0), "not-a-knot"))
    fs = lambda x: f(x)
    h = sig/16
    xs = np.array([1.25, 2.5, 3.75, 3.0])*sig
    fl = xs*(fs(xs+h)-fs(xs-h))/(2*h); flc = flux(rc, uc, xs, h)
    xf = np.array([0.25, 0.5, 0.75, 0.9])*R
    ff = xf*(fs(xf+h)-fs(xf-h))/(2*h); ffc = flux(rc, uc, xf, h)
    print(f"  {tag}: u(0)/cont {u[0]/uc[0]:.4f}; slope ratio 1.25/2.5/3.75 sig: " + " ".join(f"{a/b:.4f}" for a, b in zip(fl[:3], flc[:3]))
          + f"; flux(3sig)/cont {fl[3]/flc[3]:.4f}; flux/cont at R/4,R/2,3R/4,0.9R: " + " ".join(f"{a/b:.3f}" for a, b in zip(ff, ffc))
          + f"; Gauss deficit 1-flux(R/2)/flux(3sig) {1-ff[1]/fl[3]:+.3f}", flush=True)

if __name__ == "__main__":
    # usage: lin3.py sigma R_over_sigma mode [step]   (mode k12 = the K.12 grid, fine = uniform step in l_k)
    sig = float(sys.argv[1]); RoS = float(sys.argv[2]); mode = sys.argv[3]
    q = 0.064/sig**2; R = RoS*sig
    if mode == "k12":
        nodes = shape_nodes(R, 2*sig, sig/4, 1.12)
    else:
        nodes = np.unique(np.r_[np.arange(0, R, float(sys.argv[4])), R])
    fr = np.array([0.25, 0.75])
    rows = np.sort((nodes[:-1, None] + fr[None, :]*np.diff(nodes)[:, None]).ravel())
    s = q*np.exp(-rows**2/(2*sig**2))
    import time; T0 = time.time()
    print(f"sigma {sig:g} R {R:g} ({mode}, {len(nodes)} nodes, h/l_k near source {nodes[1]-nodes[0]:g}, largest {np.max(np.diff(nodes)):.1f})")
    APL = rows_exact(nodes, rows, R); ALOC = rows_local(nodes, rows)
    GL = rows_geom(nodes, rows, R, "lap"); GB = rows_geom(nodes, rows, R, "B")
    for tag, A in (("field PL nonlocal + geometry nonlocal (model of K.12)", APL + 2*GB),
                   ("field PL nonlocal + geometry local", APL + 2*GL),
                   ("field PL local + geometry local (l_k -> 0 limit of the scheme)", ALOC + 2*GL),
                   ("all on the spline, nonlocal (no PL)", 3*GB),
                   ("all on the spline, local (consistent collocation)", 3*GL)):
        sol, *_ = np.linalg.lstsq(A[:, :-1], s, rcond=None)
        metrics(nodes, np.r_[sol, 0.0], sig, R, q, tag)
    print(f'  time {time.time()-T0:.0f}s', flush=True)
