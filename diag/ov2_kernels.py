"""Numba kernels of OV2 in 2D (diag/ov2_2d.py); a module never run as __main__ (see diag/ov1_kernels.py)."""
import numba as nb
import numpy as np

from diag.ov1_kernels import small_future


@nb.njit(cache=True)
def matter_sums(vr, iv, ncut, f, phi):
    """Over all related pairs y < x (x later) with n(x, y) <= ncut, for the elements in u-rank order:
    D = sum f(n) (phi_x - phi_y)^2           (pairwise, Dirichlet form)
    P = sum f(n) phi_x phi_y                  (the cross term of the plain form)
    R = sum f(n) (phi_y - phi_x) phi_x        (row-difference form, x the later element)
    F = per element: sum over its past of f(n) (for B1)."""
    N = vr.size
    tree = np.zeros(N + 1, np.int64)
    undo = np.empty(N, np.int64)
    oy = np.empty(N, np.int64)
    on = np.empty(N, np.int64)
    D = 0.0
    P = 0.0
    R = 0.0
    Fpast = np.zeros(N)
    for x in range(N):                       # x earlier, y later
        m = small_future(x, vr, iv, ncut, tree, undo, oy, on)
        for k in range(m):
            y = oy[k]
            w = f[on[k]]
            d = phi[y] - phi[x]
            D += w * d * d
            P += w * phi[x] * phi[y]
            R += w * (phi[x] - phi[y]) * phi[y]
            Fpast[y] += w
    return D, P, R, Fpast
