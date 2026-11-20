"""Numba kernels of OV1 in 2D (diag/ov1_2d*.py), in a module that is never run as __main__.

Numba's on-disk cache records the module of each compiled function. A kernel defined in a script that runs a process
pool is cached by the workers under the module name __mp_main__ (Windows spawn), and a later normal import of the
same script cannot load that cache entry (ModuleNotFoundError: __mp_main__, 19 Nov). Kernels imported from here are
always cached under diag.ov1_kernels. The code is copied unchanged from diag/ov1_2d.py, ov1_2d_mcmc.py and
ov1_2d_mcmc2.py, where it was checked against brute force.
"""
import math

import numba as nb
import numpy as np

@nb.njit(cache=True, inline="always")
def _add(tree, i, d):
    i += 1
    while i < tree.size:
        tree[i] += d
        i += i & (-i)


@nb.njit(cache=True, inline="always")
def _pre(tree, i):          # sum over [0, i]
    s = 0
    i += 1
    while i > 0:
        s += tree[i]
        i -= i & (-i)
    return s


@nb.njit(cache=True)
def small_future(x, vr, iv, ncut, tree, undo, oy, on):
    """All y > x with n(x, y) <= ncut; exact. Writes (y, n) into oy, on; returns the count. tree must be zero."""
    N = vr.size
    vx = vr[x]
    m = 0
    # case A: dv >= du (in ranks), sweep in u
    nu = 0
    cs = 0
    d = 1
    while True:
        i = x + d
        if i >= N:
            break
        if d >= 2:                       # grow the square (x, x+d) x (vx, vx+d)
            j = x + d - 1
            if vx < vr[j] < vx + d:
                cs += 1
            r = vx + d - 1
            if r > vx and r < N:
                jj = iv[r]
                if x < jj < x + d - 1:
                    cs += 1
        if cs > ncut:
            break
        if vr[i] > vx and vr[i] - vx >= d:
            n = _pre(tree, vr[i] - 1) - _pre(tree, vx)
            if n <= ncut:
                oy[m] = i
                on[m] = n
                m += 1
        if vr[i] > vx:
            _add(tree, vr[i], 1)
            undo[nu] = vr[i]
            nu += 1
        d += 1
    for k in range(nu):
        _add(tree, undo[k], -1)
    # case B: du > dv, sweep in v
    nu = 0
    cs = 0
    d = 1
    while True:
        r = vx + d
        if r >= N:
            break
        if d >= 2:
            j = x + d - 1
            if j < N and vx < vr[j] < vx + d - 1:
                cs += 1
            jj = iv[vx + d - 1]
            if x < jj < x + d:
                cs += 1
        if cs > ncut:
            break
        j = iv[r]
        if j > x and j - x > d:
            n = _pre(tree, j - 1) - _pre(tree, x)
            if n <= ncut:
                oy[m] = j
                on[m] = n
                m += 1
        if j > x:
            _add(tree, j, 1)
            undo[nu] = j
            nu += 1
        d += 1
    for k in range(nu):
        _add(tree, undo[k], -1)
    return m


@nb.njit(cache=True)
def abundances(vr, iv, ncut, lev_of, nlev):
    """N_n (n <= ncut) per nested level; a pair belongs to every level that contains its lower element."""
    N = vr.size
    H = np.zeros((nlev, ncut + 1), np.int64)
    tree = np.zeros(N + 1, np.int64)
    undo = np.empty(N, np.int64)
    oy = np.empty(N, np.int64)
    on = np.empty(N, np.int64)
    for x in range(N):
        m = small_future(x, vr, iv, ncut, tree, undo, oy, on)
        for k in range(m):
            for L in range(lev_of[x] + 1):
                H[L, on[k]] += 1
    return H


@nb.njit(cache=True)
def rel_count_height(vr):
    """Number of relations and the longest chain (elements)."""
    N = vr.size
    tree = np.zeros(N + 1, np.int64)
    R = 0
    for i in range(N):
        R += _pre(tree, vr[i] - 1) if vr[i] > 0 else 0
        _add(tree, vr[i], 1)
    tails = np.empty(N, np.int64)
    L = 0
    for i in range(N):
        a, b = 0, L
        while a < b:
            c = (a + b) // 2
            if tails[c] < vr[i]:
                a = c + 1
            else:
                b = c
        tails[a] = vr[i]
        if a == L:
            L += 1
    return R, L


@nb.njit(cache=True)
def _past_lists(x, vrR, ivR, ncut, tree, undo, oy, on):
    N = vrR.size
    m = small_future(N - 1 - x, vrR, ivR, ncut, tree, undo, oy, on)
    for k in range(m):
        oy[k] = N - 1 - oy[k]
    return m


@nb.njit(cache=True)
def _dF(p, q, vr, f, fy, fn, mf, py, pn, mp, gy, gn, mg, qy, qn, mq):
    """Change of sum_pairs f(n) when the U-adjacent pair (p, q = p+1) is swapped.
    (fy, fn, mf) future of p, (py, pn, mp) past of p, (gy, gn, mg) future of q, (qy, qn, mq) past of q."""
    nf = f.size - 2
    if vr[q] > vr[p]:                    # link p < q removed
        dF = -f[0]
        for k in range(mf):
            y = fy[k]
            if y != q and vr[y] > vr[q]:
                n = fn[k]
                dF += f[n - 1] - f[n]
        for k in range(mq):
            x = qy[k]
            if x != p and vr[x] < vr[p]:
                n = qn[k]
                dF += f[n - 1] - f[n]
    else:                                # link q < p added
        dF = f[0]
        for k in range(mg):
            y = gy[k]
            if vr[y] > vr[p]:
                n = gn[k]
                if n + 1 <= nf:
                    dF += f[n + 1] - f[n]
        for k in range(mp):
            x = py[k]
            if vr[x] < vr[q]:
                n = pn[k]
                if n + 1 <= nf:
                    dF += f[n + 1] - f[n]
    return dF


@nb.njit(cache=True)
def balance(vr, iv, vrR, ivR, f, ncut, eps):
    """Delta S of the swap of every U-adjacent pair (p, p+1); returns dS[p] for p = 0..N-2."""
    N = vr.size
    tree = np.zeros(N + 1, np.int64)
    undo = np.empty(N, np.int64)
    A = [np.empty(N, np.int64) for _ in range(8)]
    fy, fn, py, pn, gy, gn, qy, qn = A[0], A[1], A[2], A[3], A[4], A[5], A[6], A[7]
    dS = np.empty(N - 1)
    mf = small_future(0, vr, iv, ncut, tree, undo, fy, fn)
    mp = _past_lists(0, vrR, ivR, ncut, tree, undo, py, pn)
    for p in range(N - 1):
        q = p + 1
        mg = small_future(q, vr, iv, ncut, tree, undo, gy, gn)
        mq = _past_lists(q, vrR, ivR, ncut, tree, undo, qy, qn)
        dS[p] = -8.0 * eps * eps * _dF(p, q, vr, f, fy, fn, mf, py, pn, mp, gy, gn, mg, qy, qn, mq)
        fy, gy = gy, fy
        fn, gn = gn, fn
        py, qy = qy, py
        pn, qn = qn, pn
        mf, mp = mg, mq
    return dS


@nb.njit(cache=True)
def _apply_U(vr, iv, vrR, ivR, p):
    N = vr.size
    a, b = vr[p], vr[p + 1]
    vr[p], vr[p + 1] = b, a
    iv[b], iv[a] = p, p + 1
    for k in (N - 1 - p, N - 2 - p):
        vrR[k] = N - 1 - vr[N - 1 - k]
        ivR[vrR[k]] = k


@nb.njit(cache=True)
def _apply_V(vr, iv, vrR, ivR, r):
    N = vr.size
    i, j = iv[r], iv[r + 1]
    vr[i], vr[j] = r + 1, r
    iv[r], iv[r + 1] = j, i
    for e in (i, j):
        k = N - 1 - e
        vrR[k] = N - 1 - vr[e]
        ivR[vrR[k]] = k


@nb.njit(cache=True)
def _dS_U(p, vr, iv, vrR, ivR, f, ncut, eps, tree, undo, L):
    q = p + 1
    mf = small_future(p, vr, iv, ncut, tree, undo, L[0], L[1])
    mp = _past_lists(p, vrR, ivR, ncut, tree, undo, L[2], L[3])
    mg = small_future(q, vr, iv, ncut, tree, undo, L[4], L[5])
    mq = _past_lists(q, vrR, ivR, ncut, tree, undo, L[6], L[7])
    return -8.0 * eps * eps * _dF(p, q, vr, f, L[0], L[1], mf, L[2], L[3], mp, L[4], L[5], mg, L[6], L[7], mq)


@nb.njit(cache=True)
def _seed(s):
    np.random.seed(s)


@nb.njit(cache=True)
def run_sweeps(A, B, S, R, beta, f, ncut, eps, n_prop):
    """n_prop proposals; A = (vr, iv, vrR, ivR) natural, B = the transposed representation (iv, vr, ...).
    Returns S, R (relations), accepted."""
    vr, iv, vrR, ivR = A
    tv, tiv, tvR, tivR = B
    N = vr.size
    tree = np.zeros(N + 1, np.int64)
    undo = np.empty(N, np.int64)
    L = [np.empty(N, np.int64) for _ in range(8)]
    acc = 0
    for _ in range(n_prop):
        useU = np.random.random() < 0.5
        p = np.random.randint(0, N - 1)
        if useU:
            dS = _dS_U(p, vr, iv, vrR, ivR, f, ncut, eps, tree, undo, L)
            rel = vr[p] < vr[p + 1]
        else:
            dS = _dS_U(p, tv, tiv, tvR, tivR, f, ncut, eps, tree, undo, L)
            rel = tv[p] < tv[p + 1]
        if dS <= 0.0 or np.random.random() < math.exp(-beta * dS):
            if useU:
                _apply_U(vr, iv, vrR, ivR, p)
                _apply_V(tv, tiv, tvR, tivR, p)
            else:
                _apply_U(tv, tiv, tvR, tivR, p)
                _apply_V(vr, iv, vrR, ivR, p)
            S += dS
            R += -1 if rel else 1
            acc += 1
    return S, R, acc


@nb.njit(cache=True)
def _init(pu, pv):
    N = pu.size
    C = np.zeros((N, N), np.int32)
    for x in range(N):
        for y in range(N):
            if pu[x] < pu[y] and pv[x] < pv[y]:
                c = 0
                for z in range(N):
                    if pu[x] < pu[z] < pu[y] and pv[x] < pv[z] < pv[y]:
                        c += 1
                C[x, y] = c
    return C


@nb.njit(cache=True)
def _F(pu, pv, C, f):
    N = pu.size
    s = 0.0
    R = 0
    L0 = 0
    for x in range(N):
        for y in range(N):
            if pu[x] < pu[y] and pv[x] < pv[y]:
                s += f[C[x, y]]
                R += 1
                if C[x, y] == 0:
                    L0 += 1
    return s, R, L0


@nb.njit(cache=True)
def _move(P, Q, eP, p, C, f, accept_beta, eps, u01):
    """Proposal in the order P (positions P, the other order Q, eP[k] = element at position k of P): swap the
    elements at positions p, p+1. Applies it if accepted; returns (accepted, dS, dR)."""
    a, b = eP[p], eP[p + 1]
    N = P.size
    if Q[a] < Q[b]:                      # link a < b removed
        dF = -f[0]
        for y in range(N):
            if y != b and P[y] > P[b] and Q[y] > Q[b] and Q[y] > Q[a]:
                n = C[a, y]
                dF += f[n - 1] - f[n]
        for x in range(N):
            if x != a and P[x] < P[a] and Q[x] < Q[a] and Q[x] < Q[b]:
                n = C[x, b]
                dF += f[n - 1] - f[n]
        dS = -8.0 * eps * eps * dF
        if dS <= 0.0 or u01 < math.exp(-accept_beta * dS):
            for y in range(N):
                if y != b and P[y] > P[b] and Q[y] > Q[b] and Q[y] > Q[a]:
                    C[a, y] -= 1
            for x in range(N):
                if x != a and P[x] < P[a] and Q[x] < Q[a] and Q[x] < Q[b]:
                    C[x, b] -= 1
            C[a, b] = 0
            P[a], P[b] = P[b], P[a]
            eP[p], eP[p + 1] = b, a
            return True, dS, -1
        return False, dS, 0
    else:                                # link b < a created
        dF = f[0]
        for y in range(N):
            if P[y] > P[b] and Q[y] > Q[a] and y != a:
                n = C[b, y]
                dF += f[n + 1] - f[n]
        for x in range(N):
            if P[x] < P[a] and Q[x] < Q[b] and x != b:
                n = C[x, a]
                dF += f[n + 1] - f[n]
        dS = -8.0 * eps * eps * dF
        if dS <= 0.0 or u01 < math.exp(-accept_beta * dS):
            for y in range(N):
                if P[y] > P[b] and Q[y] > Q[a] and y != a:
                    C[b, y] += 1
            for x in range(N):
                if P[x] < P[a] and Q[x] < Q[b] and x != b:
                    C[x, a] += 1
            C[b, a] = 0
            P[a], P[b] = P[b], P[a]
            eP[p], eP[p + 1] = b, a
            return True, dS, 1
        return False, dS, 0


@nb.njit(cache=True)
def sweeps(pu, pv, eu, ev, C, f, beta, eps, n_prop, S, R):
    acc = 0
    N = pu.size
    for _ in range(n_prop):
        useU = np.random.random() < 0.5
        p = np.random.randint(0, N - 1)
        u01 = np.random.random()
        if useU:
            ok, dS, dR = _move(pu, pv, eu, p, C, f, beta, eps, u01)
        else:
            ok, dS, dR = _move(pv, pu, ev, p, C, f, beta, eps, u01)
        if ok:
            acc += 1
            S += dS
            R += dR
    return S, R, acc


@nb.njit(cache=True)
def height(pu, pv):
    N = pu.size
    seq = np.empty(N, np.int64)
    for e in range(N):
        seq[pu[e]] = pv[e]
    tails = np.empty(N, np.int64)
    L = 0
    for i in range(N):
        a, b = 0, L
        while a < b:
            c = (a + b) // 2
            if tails[c] < seq[i]:
                a = c + 1
            else:
                b = c
        tails[a] = seq[i]
        if a == L:
            L += 1
    return L - 1


@nb.njit(cache=True)
def _rel_count(pu, pv):
    N = pu.size
    R = 0
    for x in range(N):
        for y in range(N):
            if pu[x] < pu[y] and pv[x] < pv[y]:
                R += 1
    return R

