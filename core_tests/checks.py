"""
core_tests.checks
=================

Fast checks of the core and the instruments (python run.py core). Each
check returns (ok, message). They take seconds; they do not replace the
calibration tests of phase K, which measure instruments at the sizes the
tests use.
"""

from __future__ import annotations

import math
import sys
import tempfile
import time
from pathlib import Path
from typing import Callable, Dict, List, Tuple

import numpy as np

CHECKS: Dict[str, Callable[[], Tuple[bool, str]]] = {}


def check(fn):
    CHECKS[fn.__name__] = fn
    return fn


@check
def links_bruteforce():
    """Shell-cut links equal brute-force links (uniform and conformal sprinklings, 1+1 to 3+1)."""
    from core.order.sprinkle import sprinkle, Cylinder
    from core.order.links import links_csr, links_bruteforce as brute
    bad = []
    for D in (2, 3, 4):
        S = sprinkle(Cylinder(D, 0, 1, 0.5), 800, seed=D)
        a, b = links_csr(S.X, S.rho), brute(S.X)
        if not (np.array_equal(a[0], b[0]) and np.array_equal(a[1], b[1])):
            bad.append(f"uniform D={D}")
        dens = lambda P: (1.0 - 0.3 / np.sqrt((P[:, 1:] ** 2).sum(axis=1) + 0.1 ** 2) * 0.1) ** D
        S2 = sprinkle(Cylinder(D, 0, 1, 0.5), 800, seed=10 + D, density=dens)
        a, b = links_csr(S2.X, S2.rho * 0.5), brute(S2.X)
        if not (np.array_equal(a[0], b[0]) and np.array_equal(a[1], b[1])):
            bad.append(f"conformal D={D}")
    return (not bad), ("all equal" if not bad else "; ".join(bad))


@check
def longest_chains():
    """chain_dp equals a brute-force longest chain; heaviest-chain counts equal enumeration."""
    from core.order.sprinkle import Interval
    from core.order.relation import chain_dp, count_heaviest, prec_matrix
    rng = np.random.default_rng(3)
    bad = []
    for D in (2, 3, 4):
        X = Interval(D, 1.0).propose(rng, 60)
        X = X[np.argsort(X[:, 0])]
        R = prec_matrix(X, X)
        L = np.ones(len(X))
        for j in range(len(X)):
            par = np.nonzero(R[:j, j])[0]
            if par.size:
                L[j] = 1 + L[par].max()
        if not np.array_equal(L, chain_dp(X, block=7)["W"]):
            bad.append(f"longest D={D}")
    X = Interval(3, 1.0).propose(rng, 14)
    X = X[np.argsort(X[:, 0])]
    R = prec_matrix(X, X)
    w = np.ones(len(X))
    W = chain_dp(X, w)["W"]
    lc = count_heaviest(X, w, W)
    best = int(W.max())
    brute = 0
    for k in range(len(X)):
        def rec(path):
            nonlocal brute
            last = path[-1]
            if len(path) == best and last == k:
                brute += 1
            for nx in np.nonzero(R[last])[0]:
                if nx <= k:
                    rec(path + [nx])
        for s in range(len(X)):
            rec([s])
    tot = sum(math.exp(lc[k]) for k in range(len(X)) if W[k] == best)
    if abs(tot - brute) > 1e-6:
        bad.append(f"counts {tot} vs {brute}")
    return (not bad), ("all equal" if not bad else "; ".join(bad))


@check
def routing_identities():
    """Flux continuity: node-set balance and closed-surface identity to machine precision."""
    from core.order.sprinkle import sprinkle, Cylinder
    from core.order.links import links_csr
    from core.flux.routing import route_uniform, link_endpoints, node_set_balance, net_side_flux, cap_flux
    S = sprinkle(Cylinder(3, 0, 2, 0.5), 2000, seed=7)
    ptr, idx = links_csr(S.X, S.rho)
    r = S.r()
    src = (r < 0.08).astype(float)
    worst = 0.0
    for keep, absorb in ((1.0, None), (0.9, r > 0.4), (1.0, r > 0.4)):
        F, Phi = route_uniform(ptr, idx, src, keep=keep, absorb=absorb)
        u, v = link_endpoints(ptr, idx)
        for rr in (0.1, 0.25):
            ins = (r < rr) & (S.X[:, 0] >= 0.5) & (S.X[:, 0] <= 1.5)
            bal = node_set_balance(ptr, idx, src, F, Phi, ins, keep=keep, absorb=absorb)
            side = net_side_flux(S.X, u, v, Phi, rr, 0.5, 1.5)["net"]
            caps = cap_flux(S.X, u, v, Phi, 1.5, rr) - cap_flux(S.X, u, v, Phi, 0.5, rr)
            worst = max(worst, abs(bal["residual"]), abs(bal["sources"] - (side + caps + bal["sinks"] + bal["discarded"])))
    return worst < 1e-9, f"worst residual {worst:.1e}"


@check
def stock_profile_equals_caps():
    """stock_profile (all radii at once) equals cap_flux at every (t, r)."""
    from core.order.sprinkle import sprinkle, Cylinder
    from core.order.links import links_csr
    from core.flux.routing import route_uniform, link_endpoints, cap_flux, stock_profile
    S = sprinkle(Cylinder(3, 0, 2, 0.5), 2000, seed=7)
    ptr, idx = links_csr(S.X, S.rho)
    r = S.r()
    imp = ((r < 0.08) & (S.X[:, 0] < 0.25)).astype(float)
    F, Phi = route_uniform(ptr, idx, imp, absorb=r > 0.4)
    u, v = link_endpoints(ptr, idx)
    times, radii = [0.5, 1.0, 1.5], [0.1, 0.2, 0.3]
    P = stock_profile(S.X, u, v, Phi, times, radii)
    ref = np.array([[cap_flux(S.X, u, v, Phi, t, rr) for rr in radii] for t in times])
    err = float(np.abs(P - ref).max())
    return err < 1e-12, f"max |difference| {err:.1e}"


@check
def bd_constants_and_null_rows():
    """BD normalisations; the kernel at eps = 1 is the coefficient list; permuted-weight null rows keep the row sum on f = 1."""
    from core.operators.bd import smeared_kernel, BD_NORM, BD_COEFFS, bd_static_rows, stratified_rows
    from core.order.sprinkle import sprinkle, Cylinder
    ok_n = abs(BD_NORM[2] - 2.0) < 1e-15 and abs(BD_NORM[4] - 4 / math.sqrt(6)) < 1e-15 and abs(BD_NORM[3] - 0.90668) < 1e-4
    ok_k = all(np.allclose(smeared_kernel(np.arange(len(C)), d, 1.0), C) for d, C in BD_COEFFS.items())
    S = sprinkle(Cylinder(3, 0, 1, 0.5), 1500, seed=4)
    nodes = np.linspace(0, 0.5, 11)
    rows = stratified_rows(S.X, S.r(), 0.8, 1.0, nodes, 0.45, 1, np.random.default_rng(0))
    B = bd_static_rows(S.X, S.rho, 0.2, rows, nodes, 0.5, null_rng=np.random.default_rng(1))
    one = np.ones(len(nodes))
    ok_r = np.allclose(B["A"] @ one, B["A_null"] @ one, rtol=1e-12, atol=1e-12)
    return ok_n and ok_k and ok_r, f"a_3 = {BD_NORM[3]:.5f}; kernel at eps = 1 ok {ok_k}; null rows keep f = 1 sums {ok_r}"


@check
def bd_expectation():
    """Expected BD operator: 1+1 light cone, local limit, hat rows vs direct quadrature (as K.4 quadrature)."""
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from phases.K_meridla.K_04_bd_expected import _quadrature, TEST
    ok, lines, data = _quadrature(TEST["params"])
    return ok, f"light cone {data['lightcone']:.4f} vs {data['quadrature']:.4f}; local " + \
        ", ".join(f"{D - 1}+1 {v:.4f}" for D, v in data["local_limit"].items())


@check
def curved_rows_flat_limit():
    """Curved expected operator (g = Pi^2 eta): Pi = 1 gives the flat rows; constant Pi = c gives Phi = c^D exactly."""
    import numpy as np
    from core.operators.bd_expect import expected_radial_rows
    from core.operators.bd_expect_curved import curved_radial_rows, diamond_phi_spheroid
    nodes = np.linspace(0.0, 12.0, 25)
    worst = 0.0
    for D in (3, 4):
        q = (40, 6, 8, 3, 6, 3)
        A0 = expected_radial_rows(D, 20.0, 0.05, 36.0, [1.0, 5.0], 12.0, nodes, n_rho=q[0], n_s=q[1], n_ang=q[2])
        for geom in ("linear", "spline"):          # a constant field is exact in both interpolations of the geometry
            A1 = curved_radial_rows(D, 20.0, 0.05, 36.0, [1.0, 5.0], 12.0, nodes, np.ones(nodes.size), *q, geom=geom)
            worst = max(worst, float(np.max(np.abs(A1 - A0)) / np.max(np.abs(A0))))
    phi = diamond_phi_spheroid(lambda r: np.full(np.shape(r), 0.8), 4, 3.0, np.array([[2.0, 7.0]]), -0.4,
                               np.array([[2.5, 7.1]]), np.array([[1.5, 1.19]]))
    dphi = float(np.max(np.abs(phi - 0.8 ** 4)))
    return worst < 1e-10 and dphi < 1e-12, f"rows rel. difference {worst:.1e}; constant-Pi diamond error {dphi:.1e}"


@check
def feedback_reduces_to_capacity():
    """F3.3: with Pi_sat = 1 every feedback channel is exactly the capacity routing; the balance of axiom 3 holds."""
    import numpy as np
    from core.order.sprinkle import sprinkle, Cylinder
    from core.order.links import links_csr
    from core.flux.capacity import route_capacity
    from core.flux.feedback import route_feedback, balance_residual_fb, MODES
    S = sprinkle(Cylinder(3, 0.0, 3.0, 0.5), 1500, 5)
    ptr, idx = links_csr(S.X, S.rho)
    absorb = S.r() > 0.4
    s = np.ones(S.n)
    C = float(np.quantile(route_capacity(ptr, idx, s, np.inf, absorb)["offered"], 0.95))
    ref = route_capacity(ptr, idx, s, C, absorb)
    same = all(np.array_equal(route_feedback(ptr, idx, s, C, absorb, md, 1.0, 2.0)[k], ref[k])
               for md in MODES for k in ("offered", "accepted", "stored", "Phi"))
    bal = max(balance_residual_fb(route_feedback(ptr, idx, s, C, absorb, md, 0.5, 2.0)) for md in MODES)
    return same and bal < 1e-12, f"Pi_sat = 1 identical in all modes: {same}; balance at Pi_sat = 0.5 worst {bal:.1e}"


@check
def stock_routing():
    """F3.3b: carried stock balances exactly, discarding it is route_capacity, any linear extension gives the same result."""
    import numpy as np
    from core.order.sprinkle import sprinkle, Cylinder
    from core.order.links import links_csr
    from core.flux.capacity import route_capacity
    from core.flux.stock import route_stock, stock_balance
    S = sprinkle(Cylinder(3, 0.0, 3.0, 0.5), 1500, 6)
    ptr, idx = links_csr(S.X, S.rho)
    absorb = S.r() > 0.4
    s = np.ones(S.n)
    C = float(np.quantile(route_capacity(ptr, idx, s, np.inf, absorb)["offered"], 0.95))
    ref = route_capacity(ptr, idx, s, C, absorb)
    off = route_stock(ptr, idx, s, C, absorb, keep_stock=False)
    same = np.array_equal(off["accepted"], ref["accepted"]) and np.array_equal(off["stock"], ref["stored"])
    on = route_stock(ptr, idx, s, C, absorb)
    tb = (S.X[:, 0] - 0.6 * S.X[:, 1]) / 0.8
    boost = route_stock(ptr, idx, s, C, absorb, order=np.argsort(tb, kind="stable"))
    ident = all(np.array_equal(on[k], boost[k]) for k in ("accepted", "stock", "Phi", "Psi"))
    bal = stock_balance(s, on)
    return same and ident and bal < 1e-12, f"discarded = capacity routing {same}; boosted order identical {ident}; balance {bal:.1e}"


@check
def pulse_routing():
    """F3.3c: no pulse is route_capacity; the copied pulse keeps the flux balance; any linear extension gives the same result."""
    import numpy as np
    from core.order.sprinkle import sprinkle, Cylinder
    from core.order.links import links_csr
    from core.flux.capacity import route_capacity
    from core.flux.pulse import route_pulse, pulse_balance
    S = sprinkle(Cylinder(3, 0.0, 3.0, 0.5), 1500, 7)
    ptr, idx = links_csr(S.X, S.rho)
    absorb = S.r() > 0.4
    s = np.ones(S.n)
    C = float(np.quantile(route_capacity(ptr, idx, s, np.inf, absorb)["offered"], 0.9))
    ref = route_capacity(ptr, idx, s, C, absorb)
    none_ = route_pulse(ptr, idx, s, C, absorb, "none", 1.0)
    same = np.array_equal(none_["accepted"], ref["accepted"]) and np.array_equal(none_["stored"], ref["stored"])
    tb = (S.X[:, 0] - 0.6 * S.X[:, 1]) / 0.8
    ok_order, bal = True, 0.0
    for rule in ("max", "main"):
        a = route_pulse(ptr, idx, s, C, absorb, rule, 1.0)
        b = route_pulse(ptr, idx, s, C, absorb, rule, 1.0, order=np.argsort(tb, kind="stable"))
        ok_order &= all(np.array_equal(a[k], b[k]) for k in ("accepted", "stored", "P"))
        bal = max(bal, pulse_balance(a))
    return same and ok_order and bal < 1e-12, f"no pulse = capacity routing {same}; boosted order identical {ok_order}; balance {bal:.1e}"


@check
def pulse_with_feedback():
    """F3.3f: pulse 'none' + a Pi channel = route_feedback of that channel; with the pulse: balance and order-free."""
    import numpy as np
    from core.order.sprinkle import sprinkle, Cylinder
    from core.order.links import links_csr
    from core.flux.capacity import route_capacity
    from core.flux.feedback import route_feedback
    from core.flux.pulse import route_pulse, pulse_balance
    S = sprinkle(Cylinder(3, 0.0, 3.0, 0.5), 1500, 11)
    ptr, idx = links_csr(S.X, S.rho)
    absorb = S.r() > 0.4
    s = np.ones(S.n)
    C = float(np.quantile(route_capacity(ptr, idx, s, np.inf, absorb)["offered"], 0.9))
    tb = (S.X[:, 0] - 0.6 * S.X[:, 1]) / 0.8
    same, ok_order, bal = True, True, 0.0
    for ch in ("capacity", "source"):
        a = route_pulse(ptr, idx, s, C, absorb, "none", 1.0, channel=ch, Pi_sat=0.5, m=2.0)
        f = route_feedback(ptr, idx, s, C, absorb, ch, 0.5, 2.0)
        same &= np.allclose(a["accepted"], f["accepted"], rtol=0, atol=1e-12) and np.allclose(a["Pi"], f["Pi"], rtol=0, atol=1e-12)
        b = route_pulse(ptr, idx, s, C, absorb, "main", 1.0, dq_crit=0.2, channel=ch, Pi_sat=0.5, m=2.0)
        c = route_pulse(ptr, idx, s, C, absorb, "main", 1.0, order=np.argsort(tb, kind="stable"), dq_crit=0.2,
                        channel=ch, Pi_sat=0.5, m=2.0)
        ok_order &= all(np.array_equal(b[k], c[k]) for k in ("accepted", "stored", "P", "Pi"))
        bal = max(bal, pulse_balance(b))
    return same and ok_order and bal < 1e-12, f"= route_feedback {same}; boosted order identical {ok_order}; balance {bal:.1e}"


@check
def order_read_paths():
    """Geometry reading (K.1): coordinate and bitset paths agree on intervals and action counts; Myrheim-Meyer inverts."""
    import random as _r
    import numpy as np
    from core.order.sprinkle import diamond_2x
    from core.geometry import order_read as OR
    worst = 0
    for D in (3, 4):
        S = diamond_2x(D, 500, 7)
        P = OR.poset_from_points(S.X)
        anc, desc = P.anc(), P.desc()
        rng = _r.Random(3)
        for _ in range(10):
            u = rng.randrange(P.n // 3)
            d = list(OR._bits(desc[u]))
            if not d:
                continue
            w = d[rng.randrange(len(d))]
            rc = OR.interval_reading(S.X, u, w)
            Sb = (desc[u] & anc[w]) | (1 << u) | (1 << w)
            rel = sum(OR._popcount(anc[x] & Sb) for x in OR._bits(Sb))
            worst = max(worst, abs(rc["m"] - OR._popcount(Sb)) + abs(rc["rel"] - rel))
        x = P.n - 3
        cc, cb = OR.action_counts_coords(S.X, x, S.rho), OR.interval_counts_poset(anc, desc, x)
        worst = max(worst, sum(abs(cc.get(k, 0) - cb.get(k, 0)) for k in range(OR._ncut(0.05) + 1)))
    mm = max(abs(OR.mm_dimension(OR.mm_ordering_fraction(Dd)) - Dd) for Dd in (2.0, 3.0, 4.0, 5.5))
    return worst == 0 and mm < 1e-8, f"worst path difference {worst}; Myrheim-Meyer roundtrip {mm:.1e}"


@check
def continuum_ball_solution():
    """The continuum ball-source solution solves lap u = q 1(r < a) with u(R) = 0 (finite differences)."""
    from core.operators.continuum import ball_source_solution
    worst = 0.0
    for n in (1, 2, 3):
        u = ball_source_solution(n, 1.0, 0.1, 0.5)
        r = np.linspace(0.02, 0.48, 4001)
        h = r[1] - r[0]
        ur = u(r)
        lap = (ur[2:] - 2 * ur[1:-1] + ur[:-2]) / h ** 2 + (n - 1) / r[1:-1] * (ur[2:] - ur[:-2]) / (2 * h)
        target = (r[1:-1] < 0.1).astype(float)
        m = np.abs(r[1:-1] - 0.1) > 3 * h
        worst = max(worst, float(np.max(np.abs(lap[m] - target[m]))), abs(float(u(np.array([0.5]))[0])))
    return worst < 1e-4, f"worst residual {worst:.1e}"


@check
def chain_functionals_bruteforce():
    """S and U equal brute-force enumeration of all maximal chains; weights sum to 1; greedy mirror identity."""
    from core.order.sprinkle import Interval
    from core.order.links import links_csr
    from core.chains.chains import allpaths_between, shortest_between, greedy_between, _mid_position
    from core.chains.readings import with_ends
    worst, wsum = 0.0, 0.0
    for sd in range(3):
        rng = np.random.default_rng(sd)
        P = Interval(3, 1.0).propose(rng, 40)
        X = with_ends(P, 3, 1.0)
        N = len(X)
        rmin = 40 / Interval(3, 1.0).volume() * 0.5
        ptr, idx = links_csr(X, rho_min=rmin)
        ch = [[] for _ in range(N)]
        for v in range(N):
            for u in idx[ptr[v]:ptr[v + 1]]:
                ch[u].append(v)
        paths, stack = [], [[0]]
        while stack:
            pth = stack.pop()
            if pth[-1] == N - 1:
                paths.append(pth)
                continue
            for c in ch[pth[-1]]:
                stack.append(pth + [c])
        vals = []
        for pth in paths:
            for a, b in zip(pth[:-1], pth[1:]):
                if X[a, 0] <= 0.5 < X[b, 0]:
                    vals.append(_mid_position(X[a], X[b], 0.5)[0])
                    break
        L = min(len(p_) for p_ in paths)
        u = allpaths_between(X, ptr, idx, 0, N - 1, [1.0, 0.0])
        s_ = shortest_between(X, ptr, idx, 0, N - 1, [1.0, 0.0])
        worst = max(worst, abs(np.mean(vals) - u["disp"]),
                    abs(np.mean([v for v, p_ in zip(vals, paths) if len(p_) == L]) - s_["disp"]))
        wsum = max(wsum, abs(u["weight_check"] - 1.0))
        PM = P.copy()
        PM[:, 1] = -PM[:, 1]
        XM = with_ends(PM, 3, 1.0)
        pM, iM = links_csr(XM, rho_min=rmin)
        g1 = greedy_between(X, ptr, idx, np.ones(N), 0, N - 1, [1.0, 0.0])["disp"]
        g2 = greedy_between(XM, pM, iM, np.ones(N), 0, N - 1, [1.0, 0.0])["disp"]
        worst = max(worst, abs(g1 + g2))
    return worst < 1e-12 and wsum < 1e-9, f"worst |difference| {worst:.1e}; weights sum to 1 within {wsum:.1e}"


@check
def mirror_identity_small():
    """d(P) + d(MP) = 0 exactly for all six functionals at eps = 0 (2+1, 3+1, small sets)."""
    from core.chains.readings import mirror_residuals
    worst = 0.0
    for D in (3, 4):
        for k in range(2):
            r = mirror_residuals(D, 300, 100 * D + k, 1.0)
            worst = max(worst, max(r["residual"].values()))
    return worst <= 1e-12, f"worst {worst:.1e}"


@check
def capacity_routing():
    """C = inf equals uniform routing; capacity balance exact; marked probe balance, size and background checks."""
    from core.order.sprinkle import sprinkle, Cylinder
    from core.order.links import links_csr
    from core.flux.routing import route_uniform
    from core.flux.capacity import route_capacity, balance_residual, probe_check
    S = sprinkle(Cylinder(3, 0, 2, 0.5), 2000, seed=7)
    ptr, idx = links_csr(S.X, S.rho)
    absorb = S.r() > 0.4
    s = np.ones(S.n)
    F, Phi = route_uniform(ptr, idx, s, absorb=absorb)
    ri = route_capacity(ptr, idx, s, np.inf, absorb)
    e_inf = float(np.abs(ri["offered"] - F).max() + np.abs(ri["Phi"] - Phi).max())
    C = float(np.quantile(ri["offered"], 0.6))
    rc = route_capacity(ptr, idx, s, C, absorb)
    bal = balance_residual(s, rc)
    pc = probe_check(ptr, idx, s, C, absorb, 1e-6, 10.0)
    ok = e_inf < 1e-12 and bal < 1e-12 and pc["marked_balance"] < 1e-9 and pc["size_sensitivity"] < 1e-4
    return ok, f"C=inf diff {e_inf:.1e}; balance {bal:.1e}; probe balance {pc['marked_balance']:.1e}, size {pc['size_sensitivity']:.1e}"


@check
def covering_check_and_strict_null():
    """redundant_edges equals brute force; links have none; strict rewiring keeps degrees and covering relations."""
    from core.order.sprinkle import sprinkle, Cylinder
    from core.order.links import links_csr
    from core.flux.capacity import redundant_edges, rewire_strict
    rng = np.random.default_rng(3)
    n = 40
    A = np.triu(rng.random((n, n)) < 0.15, 1)
    vv, uu = np.nonzero(A.T)
    ptr = np.cumsum(np.r_[0, np.bincount(vv, minlength=n)]).astype(np.int64)
    red = redundant_edges(ptr, uu, n)
    M = A.astype(int)
    P = M.copy()
    for _ in range(n):
        P = ((P + P @ M) > 0).astype(int)
    S = sprinkle(Cylinder(3, 0, 2, 0.5), 500, seed=5)
    p2, i2 = links_csr(S.X, S.rho)
    o = rewire_strict(p2, i2, S.X[:, 0], 4000, np.random.default_rng(1), X=S.X, local=S.ell)
    ok = (np.array_equal(red, ((M @ P) > 0)[uu, vv]) and int(redundant_edges(p2, i2, S.n).sum()) == 0
          and o["redundant_after"] == 0 and o["deg_in_same"] and o["deg_out_same"])
    return ok, f"brute force equal; strict null {o['accepted_swaps']} swaps, redundant after {o['redundant_after']}"


@check
def statistics():
    """equivalence rule, bootstrap SE of a mean, log-log slope."""
    from core.stats.stats import equivalence, boot_se, slope_loglog
    eq = (equivalence(0.01, 0.01, 0.05, 2.0), equivalence(0.2, 0.01, 0.05, 2.0), equivalence(0.04, 0.02, 0.05, 2.0))
    v = np.random.default_rng(0).normal(0, 1, 400)
    m, se, _ = boot_se(v, 2000, 1)
    sl = slope_loglog([1, 2, 4, 8], [1, 0.5, 0.25, 0.125])
    ok = eq == ("EQUIVALENT", "EXCESS", "INCONCLUSIVE") and abs(se - 1 / math.sqrt(400)) < 0.01 and abs(sl + 1) < 1e-12
    return ok, f"equivalence {eq}; bootstrap SE {se:.4f} (theory 0.05); slope {sl:.3f}"


def _square_plus_seed(x, seed):
    """Top-level task for the parallel check (pure function of its arguments)."""
    return float(np.random.default_rng(seed).normal() + x * x)


@check
def parallel_identical():
    """parallel_map on 2 processes returns exactly the sequential results, in input order."""
    from core.stats.parallel import parallel_map
    args = [(k, 1000 + k) for k in range(8)]
    a = parallel_map(_square_plus_seed, args, 1, lambda m: None, lambda *x: "")
    b = parallel_map(_square_plus_seed, args, 2, lambda m: None, lambda *x: "")
    return a == b, f"{len(a)} tasks identical: {a == b}"


@check
def test_headers():
    """Every test file in phases/ has a valid header, a run() and quick parameters; ids are unique."""
    from core.registry import discover
    tests = discover()
    missing = [t for t, e in tests.items() if not e["test"].get("quick")]
    return not missing, f"{len(tests)} tests valid" + (f"; no quick parameters: {missing}" if missing else "")


@check
def settable_overrides():
    """--set: anything in a quick run; in a full run only the parameters the part declares settable."""
    from core.testing import check_overrides, merged_params, validate
    test = {"params": {"a": 1, "m": 0.05}, "parts": {"s": {"b": 2}, "c": {}}, "quick": {"a": 0},
            "settable": {"s": ["b"]}}
    full_ok = merged_params(test, "s", False, {"b": 5})
    quick_any = merged_params(test, "c", True, {"m": 0.1})
    refused = [check_overrides(test, "s", False, {"m": 0.1}) is not None,
               check_overrides(test, "c", False, {"b": 5}) is not None]
    try:
        merged_params(test, "s", False, {"m": 0.1})
        raised = False
    except ValueError:
        raised = True
    bad_header = validate(dict(test, settable={"nope": ["b"]}))
    ok = (full_ok == {"a": 1, "m": 0.05, "b": 5} and quick_any["m"] == 0.1 and quick_any["a"] == 0 and all(refused)
          and raised and any("unknown part" in b for b in bad_header))
    return ok, f"full settable applied, other keys refused {refused}, raised {raised}, header check {bool(bad_header)}"


@check
def results_roundtrip():
    """A result is saved with provenance (test and core hashes, versions) and read back; nothing is overwritten."""
    from core import provenance
    from core.registry import discover
    old = provenance.RESULTS
    try:
        with tempfile.TemporaryDirectory() as d:
            provenance.RESULTS = Path(d)
            entry = next(iter(discover().values()))
            t0 = time.time()
            f1 = provenance.save("X.1", entry["path"], "p", 0, False, {"a": 1, "b": float("nan")}, "PASSED", ["line"],
                                 {"v": np.arange(3)}, t0)
            f2 = provenance.save("X.1", entry["path"], "p", 0, False, {"a": 1}, "PASSED", ["line"], {}, t0)
            # the command line passes the hashes taken at the START of a run (core.cli.cmd_run); a stale
            # provenance.py without this argument lost a finished F3.3e result on 21 Oct 2026
            h = {"test_sha256": "0" * 64, "core_sha256": "1" * 64}
            f3 = provenance.save("X.1", entry["path"], "p", 0, False, {"a": 1}, "PASSED", ["line"], {}, t0, hashes=h)
            recs = provenance.load_all("X.1")
            first = next(r for r in recs if r["_file"].endswith(f1.name))      # file order follows the hash prefix
            third = next(r for r in recs if r["_file"].endswith(f3.name))
            ok = (len({f1, f2, f3}) == 3 and len(recs) == 3 and first["data"]["v"] == [0, 1, 2]
                  and first["params"]["b"] is None and len(first["test_sha256"]) == 64 and len(first["core_sha256"]) == 64
                  and third["test_sha256"] == h["test_sha256"] and third["files_changed_during_run"] is True)
    finally:
        provenance.RESULTS = old
    return ok, f"three records (one with start hashes), distinct files, provenance present: {ok}"


def main(argv: List[str]) -> int:
    names = [a for a in argv if a in CHECKS] or list(CHECKS)
    t_all = time.time()
    n_ok = 0
    for name in names:
        t0 = time.time()
        try:
            ok, msg = CHECKS[name]()
        except Exception as exc:                      # a crashing check is a failing check
            ok, msg = False, f"{type(exc).__name__}: {exc}"
        n_ok += bool(ok)
        print(f"{'PASS' if ok else 'FAIL'}  {name:30s} {msg}  ({time.time() - t0:.1f}s)", flush=True)
    print(f"\n{n_ok}/{len(names)} core checks passed in {time.time() - t_all:.0f}s")
    return 0 if n_ok == len(names) else 1
