"""First-order response of the interval volume to the theory metric -Pi^2 dt^2 + Pi^-2 dx^2, N = Pi^-2 = 1 + 2 a g(x).
eps/a = dV/V0 - d(tau_g^4)/tau0^4 : the departure of the volume from the flat formula in the geodesic proper time."""
import numpy as np
R = 1.5
def g_ball(x):
    r = np.sqrt((x ** 2).sum(-1))
    return np.where(r < R, 0.5 * (3 - (r / R) ** 2), R / np.maximum(r, 1e-12))
gl, wl = np.polynomial.legendre.leggauss(40)
def line(g, p, q):                      # int_p^q g dl, p, q arrays (..., 3)
    L = np.sqrt(((q - p) ** 2).sum(-1))
    out = 0.0
    for x, w in zip(gl, wl):
        out = out + 0.5 * w * g(p + (q - p) * (0.5 * (x + 1)))
    return out * L
def eps(g, u, v, dt, ns=40, ne=40, nphi=48):
    l = np.linalg.norm(v - u); e = (v - u) / l; m = 0.5 * (u + v)
    a1 = np.cross(e, [1.0, 0.3, 0.2]); a1 /= np.linalg.norm(a1); a2 = np.cross(e, a1)
    xs, ws = np.polynomial.legendre.leggauss(ns); xe, we = np.polynomial.legendre.leggauss(ne)
    s = l + (dt - l) * 0.5 * (xs + 1); wsv = 0.5 * (dt - l) * ws
    ph = 2 * np.pi * (np.arange(nphi) + 0.5) / nphi
    S, E, P = np.meshgrid(s, xe, ph, indexing="ij")
    W = wsv[:, None, None] * we[None, :, None] * (2 * np.pi / nphi) * (S ** 2 - l * l * E ** 2) / 8
    zeta = S * E / 2; rp = 0.5 * np.sqrt(np.maximum((S ** 2 - l * l) * (1 - E ** 2), 0))
    w = m + zeta[..., None] * e + rp[..., None] * (np.cos(P)[..., None] * a1 + np.sin(P)[..., None] * a2)
    V0 = (W * (dt - S)).sum()
    A = (W * (dt - S) * 2 * g(w)).sum()
    B = (W * 2 * (line(g, u + 0 * w, w) + line(g, w, v + 0 * w))).sum()
    tau2 = dt * dt - l * l
    gbar = line(g, u, v) / l
    return V0 / (np.pi * tau2 ** 2 / 24), (A - B) / V0, 4 * (dt * dt + l * l) * gbar / tau2
def show(lab, g, u, v, tau):
    u, v = np.array(u, float), np.array(v, float); l = np.linalg.norm(v - u); dt = np.sqrt(l * l + tau * tau)
    c, dv, dtau = eps(g, u, v, dt)
    print(f"{lab:34s} l {l:.2f} tau {tau:.2f}: V0 check {c:.6f}; dV/V0 {dv:+10.4f}; -dtau4 {dtau:+10.4f}; eps/a {dv + dtau:+.5f}")
print("controls (flat to first order; eps must vanish)")
show("constant g = 1", lambda x: 1.0 + 0 * x[..., 0], [3, 0, -2], [3, 0, 2], 0.7)
show("uniform field g = z", lambda x: x[..., 2], [3, 0, -2], [3, 0, 2], 0.7)
show("uniform field g = x", lambda x: x[..., 0], [3, 0, -2], [3, 0, 2], 0.7)
print("vacuum (the whole diamond outside the ball R 1.5)")
for tau in (0.4, 0.8, 1.5):
    show("tangential, b 3", g_ball, [3, 0, -2], [3, 0, 2], tau)
for tau in (0.4, 0.8):
    show("tangential, b 2.5, long", g_ball, [2.5, 0, -3], [2.5, 0, 3], tau)
    show("radial, r 2.2 to 5.2", g_ball, [0, 0, 2.2], [0, 0, 5.2], tau)
    show("oblique", g_ball, [2.5, 0.5, -1], [3.5, -1, 2.5], tau)
print("matter (the segment crosses the ball)")
for tau in (0.4, 0.8):
    show("through the centre", g_ball, [0, 0, -3], [0, 0, 3], tau)
    show("b 0.75", g_ball, [0.75, 0, -3], [0.75, 0, 3], tau)
    show("inside, centre, short", g_ball, [0, 0, -0.5], [0, 0, 0.5], tau)
    show("grazing b 1.6 (vacuum)", g_ball, [1.6, 0, -3], [1.6, 0, 3], tau if tau < 0.2 else 0.15)
