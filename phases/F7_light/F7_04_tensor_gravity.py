"""
F7.4 -- Tensor gravity from the momentum of flux: the field of a bundle, and the bending of fast probes (3+1)
===========================================================================================================

Why this test exists (23 Oct 2026, the author's decision "G2")
------------------------------------------------------------
Two results of the same day point here:
  * the bending rule of the F3 series (towards slower clocks, F3.8-F3.10)
    moves along the geodesics of the conformal geometry Pi^2 eta (PPN
    gamma = -1): a probe moving with speed v across the field is bent with
    the factor (1 - v^2) of Newton, light not at all (docs/findings.md,
    derived consequences). General relativity gives (1 + v^2), for light
    twice Newton, as measured;
  * Jacobson's condition (Clausius on local horizons, equivalent to
    R_kk = 8 pi G T_kk for every null k) holds in the linear order exactly
    when the weak field is written with the source's own 4-velocity:
    h_ab = -2U (2 u_a u_b - eta_ab) (signature + - - -), i.e. the field is
    a TENSOR sourced by the momentum of matter, not a scalar sourced by its
    amount (docs/derivations/jacobson_route.md, section 6;
    docs/derivations/jacobson_linear_check.py). F3.6 (replicated) supplies
    that momentum: flux carries its 4-velocity u as a state.
This test builds the tensor field from the flux of a bundle on a 3+1
causal set and lets fast probes move in it, next to the old scalar rule.

What is derived, what is tested, what is input
----------------------------------------------
  derived (analytic, before this test): with T = sum F u u as the source and
    the linearised Einstein equation (trace reversal), the spatial part of
    the metric equals the temporal part (gamma = 1) and a probe crossing the
    field with speed v is bent with the factor (1 + v^2); the scalar rule
    gives (1 - v^2). These follow algebraically from the structure;
  tested here: that the discrete construction realises them - the 3+1 field
    of the discrete flux has the Newtonian 1/r profile with the coefficient
    of the continuum Green function (part field), the discrete transport of
    fast probes on a causal set gives the two factors (parts probe_tensor and
    probe_scalar), and everything is frame independent (frame checks);
  input: the kernel (the retarded Green function of the wave operator in
    3+1), the coupling lambda, the trace reversal (the linearised Einstein
    equation, which Jacobson's argument supplies but the lab has not
    derived from the order), and geodesic motion of the probe (equivalence
    principle, built into the rule as in F3.10).
So this is a verification of a construction plus the numerical check of the
derived scalar factor (1 - v^2); it is not a derivation of gamma = 1 from
the causal order alone.

Model (locked before any full run; the quick runs changed it, see 'split' and s0)
-----------------------------
3+1, Poisson sprinkling rho = 1 (lengths in l), signature (+, -, -, -).
Links in 3+1 are expensive (about 100 s for 36 000 events), so the source
and each probe live in their own narrow tube, sprinkled separately (a
Poisson process restricted to a region is the Poisson process of that
region). The tubes are separate causal sets: they share only the field
(where they overlap in coordinates, nothing else is exchanged).
  source  a bundle at rest at the spatial origin: tube of radius R_s = 6 for
          t in [0, T]; source events within radius 1.5 and t < 4 emit F = 1
          each with u = (1, 0, 0, 0); transport by the carried rule of F3.6
          (the F3.10 code with zero coupling: the source does not respond to
          any field, it is an external source).
  field   of the flux events w of the source (flux F_w, carried 4-velocity
          u_w), at any point x:
              U(x)        = lambda sum_w F_w g(sigma),
              hbar_ab(x)  = -4 lambda sum_w F_w u_a(w) u_b(w) g(sigma),
              h_ab        = hbar_ab - (1/2) eta_ab hbar,   hbar = eta^ab hbar_ab,
          sigma = (x - w).(x - w) the squared interval, w in the causal past
          of x. g is a smeared version of the retarded Green function
          theta(t) delta(sigma) / (2 pi) of the wave operator in 3+1:
              g(sigma) = sigma exp(-sigma / s0) / (2 pi s0^2)   (sigma > 0),
          Lorentz invariant, causal (only timelike past), zero on the light
          cone (so its first derivatives have no delta term), normalised as
          the delta function: integral over sigma = 1 / (2 pi). For a static
          source with flux mu per unit time, outside it, the continuum limit
          of this construction is
              U_c(r) = (lambda mu / 4 pi) int_0^inf sigma exp(-sigma/s0) / (s0^2 sqrt(r^2 + sigma)) d sigma,
          which tends to the Newtonian lambda mu / (4 pi r) for r^2 >> s0
          (the smearing in retarded time is about s0 / r). s0 = 3. (Changed in
          the first quick run, before any full run: with s0 = 6 the quick run
          gave C = 0.81 against the bare 1/r, which is this smearing factor;
          the prediction is now U_c, and s0 smaller.) For a source at rest,
          h_00 = h_xx = h_yy = h_zz = -2U: gamma = 1 by the trace reversal.
          Derivatives analytic: d_a g = g'(sigma) 2 (x - w)_a.
  split   (decided after the first quick run of probe_tensor, before any full
          run) in 3+1 the retarded Green function lives ON the light cone, so
          the gradient of the discrete sum is dominated by the shot noise of
          the source events crossing the smeared cone: about lambda mu / s0 at
          every distance, while the true gradient falls as lambda mu / (4 pi
          r^2). The tensor rule (which also uses time derivatives, weighted by
          gamma^2) then breaks down (NaN in the quick run). In 2+1 (F3.8) the
          kernel fills the cone and is smooth. Therefore the two halves are
          tested separately: part field checks that the discrete flux with
          carried momentum gives on average the continuum tensor field (and
          reports the gradient noise); the probe parts move probes on a causal
          set in the continuum tensor field of the same kind of source (a
          static point mass at the origin, U = K / r, h_ab = -2U (2 u_a u_b -
          eta_ab) with u at rest, K set by U = U_target at r = b), so that they
          test the transport rules, not the noise of the 3+1 kernel.
  probes  a narrow bundle with velocity v along x passing the source at the
          impact parameter b = 6 (path x from -8 to +8, closest approach at
          x = 0): a tube of radius 2.5 sprinkled in the probe's rest frame
          and boosted (Poisson is Lorentz invariant), source events within
          rest-frame radius 1.5 and rest-frame time < 2 emit F = 1 with u =
          gamma (1, v, 0, 0). The probe starts at t_s = 16 (the field along
          the path is static by then). v = 0.6, 0.75, 0.9. The probe does not
          act on the field (test body).
  rules   each probe event takes the carried direction u (F3.6), chooses the
          child of least rapidity relative to u, and sends to it u + delta u
          with the field evaluated at the event and tau the proper time of
          the chosen link:
            newton   delta u^i = gamma^2 d_i U tau   (coordinate acceleration
                     d^2 x/dt^2 = grad U; a lab-frame REFERENCE, not covariant,
                     only to normalise the other two on the same path and field);
            tensor   delta u^a = -Gamma^a_bc u^b u^c tau, Gamma from h_ab
                     (geodesic of eta + h, first order);
            scalar   the F3.8 rule with ln Pi = -U: bending in the rest frame
                     of u towards -grad ln Pi by the rapidity |grad_perp ln Pi| tau.
          After the kick u is normalised with eta (first order).
lambda: set once from the calibration stream so that U = 0.01 at r = b
(weak field; deflections of a few 0.01 rad). Not tuned on any result.

Readout and criteria
--------------------
field   (calibration of the discrete 3+1 kernel) U at points at distances
        r = 10, 12, 15, 18, 22 in random directions (outside the source bundle,
        whose rms radius grows to about 4-5 in 3+1: seen in the quick runs,
        before any full run), times 30-40, against the continuum limit U_c(r)
        of the construction, mu the flux per unit time of the source bundle
        measured over the retarded times the field points see (it fluctuates
        by ~20 % between time slabs; seen in the quick runs). Exponent k of
        U against r (the same fit of U_c gives k_c, about 1; alternatives
        k_c -+ 1) and coefficient C = mean U / U_c (pred 1, alternative a
        factor 2 off). PASSED if |k - k_c| + 2 SE <= 1/4 and |C - 1| + 2 SE
        <= 1/4; FAILED if |k - k_c| - 2 SE > 1/4 or |C - 1| - 2 SE > 1/4. The
        1/r law itself is the kernel's (input); the part calibrates that the
        discrete sum over the actual flux realises it.
        Frame check: hbar recomputed from boosted coordinates and directions
        (beta' = 0.6) equals the boost-transformed tensor (relative 1e-9),
        and the ratio hbar'_0x / hbar'_00 in the frame where the source moves
        is reported against u'_x / u'_0 (gravitomagnetism; a scalar field
        has 0 there). Reported: gamma_field = mean(h_ii) / h_00 (1 by
        construction), the rms radius of the source bundle, and the relative
        noise of the discrete gradient |dU/dr - dU_c/dr| / |dU_c/dr| at the
        field points (the reason for the split).
probe_tensor  f_t(v) = (mean delta u^y, tensor) / (mean delta u^y, newton)
        at the end of the path (flux-weighted u^y of the carriers in the last
        two time units; u^y = 0 at the start), the same sprinklings, field and
        probe sources for both rules (paired bootstrap over sprinklings).
        Predicted 1 + v^2 (1.36, 1.5625, 1.81); the scalar alternative 1 - v^2.
        PASSED if |f_t - (1 + v^2)| + 2 SE <= v^2 / 2 for every v (a quarter
        of the distance 2 v^2 between the two predictions); FAILED if
        |f_t - (1 + v^2)| - 2 SE > v^2 / 2 for some v.
probe_scalar  the same for f_s(v) against 1 - v^2 (the derived consequence
        of 23 Oct about the F3 rule), the same tolerance.
Instrument (all parts): the transport of the tensor and scalar rules in the
frame boosted by beta' = 0.6 chooses identically (a difference starting at a
numerical tie, margin < 1e-9, is reported as a tie); the probe flux at the
end of the path is at least 0.8 of the flux after the emission (the bundle
stays inside its tube); in part field the rms radius of the source bundle is
below half the smallest field distance (the field points lie outside the
source, where the potential of the bundle is that of a point). Otherwise
INSTRUMENT_FAILED. The source bundle cannot leave its tube (there are no
events outside); its rms radius is reported. Replication with --seed 1.
The light limit (v -> 1: tensor 2, scalar 0) is reported by extrapolation,
without a verdict.

Seen before this test was written: F3.6-F3.10, F3.10b quick runs, K.13a,
F1.7, the derivations of 23 Oct. Nothing of F7.4 has been run.
"""

from __future__ import annotations

import math
import time

import numpy as np

from core.order.links import links_csr
from core.order.sprinkle import Cylinder, derive_seed, sprinkle_poisson
from core.stats.parallel import parallel_map
from core.stats.stats import boot_se
from core.testing import Result
from scipy.integrate import quad
from phases.F3_structure.F3_10_two_bundles import TIE, boost_points, children_csr, transport_species

TEST = {
    "id": "F7.4",
    "phase": "F7",
    "title": "Tensor gravity from the momentum of the flow: the field of a bundle and the bending of fast probes "
             "(3+1)",
    "theory_part": "gravity as a tensor field excited by the momentum of the flow (Jacobson's condition, F3.6); "
                   "bending of light γ = 1 against the scalar rule of the F3 series",
    "question": "Does a tensor field from a discrete flow with carried momentum in 3+1 give a Newtonian profile with "
                "the correct coefficient, bend fast probes with the factor 1 + v² (light 2×), and the old scalar "
                "rule with 1 − v²?",
    "method": "source = a bundle at rest (rule F3.6, no response), field h̄_ab = −4λ Σ F u_a u_b g(σ) with a smeared "
              "retarded Green function in 3+1, h = h̄ − ½ η h̄; probes at speeds 0.5 / 0.7 / 0.9 in narrow tubes, "
              "rules newton (reference), tensor (geodesic of η + h) and scalar (F3.8, ln Π = −U) on the same data",
    "prediction": "field: k = k_c (≈ 1), C = 1; probe_tensor: f = 1 + v²; probe_scalar: f = 1 − v²",
    "criterion": "field: PASSED |k − k_c| + 2 SE <= 1/4 and |C − 1| + 2 SE <= 1/4, FAILED − 2 SE > 1/4. probe_tensor "
                 "/ probe_scalar: PASSED |f − prediction| + 2 SE <= v²/2 for every v, FAILED − 2 SE > v²/2 for some "
                 "v. INSTRUMENT: frame check β' = 0.6, probe flow at the end >= 0.8, field points outside the source "
                 "bundle. Otherwise UNDECIDED. Replication --seed 1.",
    "null_model": "the newton rule (coordinate acceleration ∇U) on the same path and field; the scalar rule of F3.8",
    "level": "mechanism",
    "meters": ["K.6"],
    "depends_on": ["F3.6", "F3.10"],
    "params": {"D": 4, "T": 44.0, "R_s": 6.0, "r_emit": 1.5, "t_emit": 4.0, "s0": 3.0, "b": 6.0, "x_half": 8.0,
               "t_start": 16.0, "R_p": 2.5, "r_emit_p": 1.5, "t_emit_p": 2.0, "speeds": [0.6, 0.75, 0.9],
               "U_target": 0.01, "field_r": [10.0, 12.0, 15.0, 18.0, 22.0], "field_t": [30.0, 34.0, 38.0, 40.0],
               "field_dirs": 8, "end_window": 2.0, "beta_check": 0.6, "L_g": 2.0, "min_parents": 4, "tau_c": 1.0,
               "n_cal": 2, "n_seeds": 16, "n_se": 2.0, "boot_B": 4000, "flux_min": 0.8},
    "parts": {"field": {"mode": "field"}, "probe_tensor": {"mode": "probe", "rule": "tensor"},
              "probe_scalar": {"mode": "probe", "rule": "scalar"}},
    "quick": {"T": 30.0, "R_s": 5.0, "n_seeds": 2, "n_cal": 1, "boot_B": 200, "t_start": 12.0, "x_half": 5.0,
              "b": 5.0, "speeds": [0.6, 0.9], "field_r": [8.0, 10.0, 12.0], "field_t": [22.0, 26.0], "field_dirs": 4},
    "estimate": "field ~10 min, probe_tensor and probe_scalar ~10–15 min each (10 processes; source tube ~40 000 "
                "events, links ~2 min; probes ~1 000–3 000 events)",
}

ETA = np.diag([1.0, -1.0, -1.0, -1.0])        # metric of the flat background, signature (+, -, -, -)


# ============================================================================================== small geometry

def mink4(a, b):
    """Minkowski product (+, -, -, -) along the last axis."""
    return a[..., 0] * b[..., 0] - (a[..., 1:] * b[..., 1:]).sum(axis=-1)


def boost_matrix(beta):
    """Matrix B of the coordinate boost x' = B x along the first spatial axis (the frame moving with +beta)."""
    g = 1.0 / math.sqrt(1.0 - beta * beta)
    B = np.eye(4)
    B[0, 0] = B[1, 1] = g
    B[0, 1] = B[1, 0] = -g * beta
    return B


def tube_rest(p, seed, R, dur, t0=-1.0):
    """A cylinder of radius R and time [t0, t0 + dur] at rest at the origin, links and children."""
    rng = np.random.default_rng(seed)
    X = sprinkle_poisson(Cylinder(4, t0, t0 + dur, R), 1.0, rng)
    ptr, idx = links_csr(X, 1.0)
    child = np.repeat(np.arange(len(X)), np.diff(ptr))
    assert bool((idx < child).all()), "processing order is not a linear extension of the causal order"
    cptr, cidx = children_csr(ptr, idx, len(X))
    return X, ptr, idx, cptr, cidx


# ============================================================================================== the field

class TensorField:
    """
    The field of the flux events of a source (module notes): U = lambda sum F g(sigma), hbar_ab = -4 lambda sum F
    u_a u_b g(sigma), h_ab = hbar_ab - 1/2 eta_ab hbar; analytic first derivatives. P: event coordinates (n, 4); F:
    their flux; Uu: their carried 4-velocities (upper index).
    """

    def __init__(self, P, F, Uu, lam, s0):
        self.P, self.F, self.lam, self.s0 = P, F, float(lam), float(s0)
        self.ul = Uu * np.array([1.0, -1.0, -1.0, -1.0])           # lower-index 4-velocities u_a
        self.uu = np.einsum("na,nb->nab", self.ul, self.ul)        # u_a u_b per source event

    def _kernel(self, x):
        """g and g' of every source event at the point x, with d = x - w (zero outside the timelike past)."""
        d = x[None, :] - self.P
        sig = mink4(d, d)
        ok = (d[:, 0] > 0.0) & (sig > 0.0)
        e = np.where(ok, np.exp(-np.where(ok, sig, 0.0) / self.s0), 0.0)
        norm = 1.0 / (2.0 * math.pi * self.s0 ** 2)
        g = norm * np.where(ok, sig, 0.0) * e                   # sigma exp(-sigma/s0) / (2 pi s0^2)
        gp = norm * (1.0 - np.where(ok, sig, 0.0) / self.s0) * e   # dg/dsigma
        return d, g, gp

    def U(self, x):
        _, g, _ = self._kernel(x)
        return self.lam * float((self.F * g).sum())

    def hbar(self, x):
        _, g, _ = self._kernel(x)
        return -4.0 * self.lam * np.einsum("n,nab->ab", self.F * g, self.uu)

    def derivs(self, x):
        """U, dU/dx^a (4), h_ab (4, 4) and dh_bc/dx^a as an array D[a, b, c]."""
        d, g, gp = self._kernel(x)
        dl = d * np.array([1.0, -1.0, -1.0, -1.0])                 # d_a (lower index)
        dsig = 2.0 * dl                                            # d sigma / d x^a = 2 d_a
        U = self.lam * float((self.F * g).sum())
        dU = self.lam * (self.F * gp) @ dsig                       # dU/dx^a
        hb = -4.0 * self.lam * np.einsum("n,nab->ab", self.F * g, self.uu)
        dhb = -4.0 * self.lam * np.einsum("n,na,nbc->abc", self.F * gp, dsig, self.uu)
        tr = np.einsum("ab,ab->", np.linalg.inv(ETA), hb)          # hbar = eta^ab hbar_ab
        dtr = np.einsum("bc,abc->a", np.linalg.inv(ETA), dhb)
        h = hb - 0.5 * ETA * tr
        dh = dhb - 0.5 * ETA[None, :, :] * dtr[:, None, None]
        return U, dU, h, dh


class PointField:
    """
    Continuum tensor field of a static point source at the spatial origin (probe parts, module notes): U = K / r,
    h_ab = -2U (2 u_a u_b - eta_ab) = -2U delta_ab with u at rest. With a boost matrix B the field is evaluated in
    the boosted coordinates x' = B x (scalar U, lower-index tensors transformed with B^-1).
    """

    def __init__(self, K, B=None):
        self.K = float(K)
        self.B = np.eye(4) if B is None else B
        self.Bi = np.linalg.inv(self.B)

    def derivs(self, xp):
        x = self.Bi @ xp                                         # coordinates in the rest frame of the source
        r = math.sqrt(float((x[1:] ** 2).sum()))
        U = self.K / r
        dU = np.zeros(4)
        dU[1:] = -self.K * x[1:] / r ** 3                        # dU/dx^i in the rest frame (static: d_0 U = 0)
        h = -2.0 * U * np.eye(4)                                 # h_ab = -2U for a = b, 0 otherwise
        dh = -2.0 * dU[:, None, None] * np.eye(4)[None, :, :]    # d_a h_bc = -2 d_a U delta_bc
        # to the boosted coordinates: d'_a = (B^-1)^alpha_a d_alpha, lower indices with B^-1 as well
        Bi = self.Bi
        dU_p = Bi.T @ dU
        h_p = Bi.T @ h @ Bi
        dh_p = np.einsum("Aa,Bb,Cc,ABC->abc", Bi, Bi, Bi, dh)
        return U, dU_p, h_p, dh_p


def kick(rule, u, tau, dU, dh):
    """
    The velocity change of one link for the covariant rules (module notes). u: upper-index 4-velocity; dU = dU/dx^a
    and dh[a, b, c] = d h_bc / dx^a at the event; tau: proper time of the chosen link. Returns the new direction,
    normalised with eta.
    """
    if rule == "tensor":
        # geodesic of eta + h to first order: du^m/dtau = -Gamma^m_ab u^a u^b with
        # Gamma_c,ab = 1/2 (d_a h_cb + d_b h_ca - d_c h_ab), so a_c = -(1/2)(2 u^a u^b d_a h_cb - u^a u^b d_c h_ab)
        t1 = np.einsum("a,b,acb->c", u, u, dh)                  # u^a u^b d_a h_cb
        t2 = np.einsum("a,b,cab->c", u, u, dh)                  # u^a u^b d_c h_ab
        a_low = -0.5 * (2.0 * t1 - t2)
        new = u + (np.linalg.inv(ETA) @ a_low) * tau            # raise the index with eta^-1
    elif rule == "scalar":
        # the F3.8 rule with ln Pi = -U: in the rest frame of u, bend towards decreasing ln Pi (slower clocks)
        # by the rapidity |grad_perp ln Pi| tau
        G = -dU                                                 # G_a = d ln Pi / dx^a
        Gp = G - (G @ u) * (ETA @ u)                            # projection orthogonal to u: G_a - (G.u) u_a
        Gp_up = np.linalg.inv(ETA) @ Gp
        n2 = -float(Gp @ Gp_up)                                 # |grad_perp ln Pi|^2 (spacelike, so G.G < 0)
        if n2 <= 0.0:
            return u
        gn = math.sqrt(n2)
        # the spatial gradient vector with the positive spatial metric -eta is -Gp_up; bending goes against it
        e = Gp_up / gn
        eta_k = gn * tau
        new = math.cosh(eta_k) * u + math.sinh(eta_k) * e
    else:
        raise ValueError(rule)
    return new / math.sqrt(max(float(new @ ETA @ new), 1e-300))


def newton_kick(u, tau, dU):
    """Lab-frame reference: d^2 x^i / dt^2 = dU/dx^i, i.e. delta u^i = gamma^2 (dU/dx^i) tau; then normalised."""
    g0 = u[0]
    new = u.copy()
    new[1:] = u[1:] + g0 ** 2 * dU[1:] * tau
    new[0] = math.sqrt(1.0 + float((new[1:] ** 2).sum()))
    return new


# ============================================================================================== transport of a probe

def probe_transport(X, ptr, idx, cptr, cidx, src, u_src, rule, field):
    """
    One causal pass of a probe bundle (carried rule of F3.6) responding to the external field with the given rule.
    Returns flux, directions, decisions and decision margins per event.
    """
    N = len(X)
    inflow = np.zeros(N)
    Usum = np.zeros((N, 4))
    flux = np.zeros(N)
    dirs = np.zeros((N, 4))
    dec = np.full(N, -1, dtype=np.int64)
    gap = np.full(N, np.inf)
    for v in range(N):                                           # increasing index = linear extension of the order
        tot = inflow[v] + src[v]
        if tot <= 0.0:
            continue
        P = Usum[v] + src[v] * u_src
        u = P / math.sqrt(max(float(mink4(P, P)), 1e-300))       # carried direction (F3.6)
        flux[v], dirs[v] = tot, u
        a, b = cptr[v], cptr[v + 1]
        if b == a:
            continue                                             # no child: the flux ends here
        ks = cidx[a:b]
        B = X[ks] - X[v]
        tau = np.sqrt(np.maximum(mink4(B, B), 1e-300))
        ch = mink4(u[None, :], B) / tau                          # cosh of the rapidity of each child relative to u
        j = int(np.argmin(ch))
        k = int(ks[j])
        dec[v] = k
        if ch.size > 1:
            gap[v] = float(np.partition(ch, 1)[1] - ch[j])
        U, dU, h, dh = field.derivs(X[v])                       # the external field at the event (its causal past)
        if rule == "newton":
            send = newton_kick(u, float(tau[j]), dU)
        else:
            send = kick(rule, u, float(tau[j]), dU, dh)
        inflow[k] += tot
        Usum[k] += tot * send
    return {"flux": flux, "dir": dirs, "dec": dec, "gap": gap}


def first_difference(ra, rb, map_dir):
    """'identical' / 'tie' / 'differs' for two probe runs (choices, flux bit for bit, directions within 1e-9)."""
    dif = np.nonzero(ra["dec"] != rb["dec"])[0]
    if dif.size:
        v = int(dif[0])
        return "tie" if min(ra["gap"][v], rb["gap"][v]) < TIE else "differs"
    if not np.array_equal(ra["flux"], rb["flux"]):
        return "differs"
    return "identical" if float(np.max(np.abs(ra["dir"] - map_dir(rb["dir"])))) <= 1e-9 else "differs"


# ============================================================================================== source

def source_run(p, seed):
    """The source bundle at rest (carried rule, no response), its flux events and the measured flux per unit time."""
    T = float(p["T"])
    X, ptr, idx, cptr, cidx = tube_rest(p, seed, float(p["R_s"]), T + 1.0, -1.0)
    r2 = (X[:, 1:] ** 2).sum(axis=1)
    src = ((X[:, 0] >= 0.0) & (X[:, 0] < float(p["t_emit"])) & (r2 <= float(p["r_emit"]) ** 2)).astype(float)
    u0 = np.array([1.0, 0.0, 0.0, 0.0])
    r = transport_species(X, ptr, idx, cptr, cidx, [src], [[0]], u0, 0.0, p)
    fl, dirs, dec, gap = r["flux"][0], r["dir"][0], r["dec"][0], r["gap"][0]
    car = fl > 0.0
    w = (X[:, 0] >= 0.4 * T) & (X[:, 0] <= 0.9 * T) & car
    mu = float(fl[w].sum() / (0.5 * T))                          # flux per unit time in the middle of the run
    width = float(np.sqrt((fl[w] * r2[w]).sum() / fl[w].sum()))  # flux-weighted rms radius of the bundle
    return {"X": X, "ptr": ptr, "idx": idx, "cptr": cptr, "cidx": cidx, "src": src, "flux": fl, "dir": dirs,
            "dec": dec, "gap": gap, "mu": mu, "width": width}


def make_field(s, lam, p, B=None):
    """TensorField of the source run s; with B (a boost matrix) in the boosted coordinates."""
    car = s["flux"] > 0.0
    P, Uu = s["X"][car], s["dir"][car]
    if B is not None:
        P, Uu = P @ B.T, Uu @ B.T
    return TensorField(P, s["flux"][car], Uu, lam, float(p["s0"]))


def calibration_task(p, seed):
    """lambda such that U = U_target at r = b: the field of a source with lambda = 1 at r = b, at the probe times."""
    s = source_run(p, seed)
    f1 = make_field(s, 1.0, p)
    rng = np.random.default_rng(seed + 1)
    vals = []
    for t in np.linspace(float(p["t_start"]), float(p["T"]) - 4.0, 6):
        for _ in range(6):
            n = rng.normal(size=3)
            n /= np.linalg.norm(n)
            vals.append(f1.U(np.concatenate([[t], float(p["b"]) * n])))
    return float(np.mean(vals))


def U_continuum(r, lam, mu, s0):
    """Continuum limit of the construction for a static source outside it (module notes): tends to lam mu / (4 pi r)."""
    val, _ = quad(lambda sg: sg * math.exp(-sg / s0) / (s0 * s0 * math.sqrt(r * r + sg)), 0.0, math.inf)
    return lam * mu / (4.0 * math.pi) * val


# ============================================================================================== tasks

def field_task(p, seed, lam):
    """Part field: U at the sample points, the frame check of the tensor, gamma_field, the source width."""
    s = source_run(p, seed)
    fld = make_field(s, lam, p)
    rng = np.random.default_rng(seed + 7)
    rows = []
    pts = []
    for t in p["field_t"]:
        for r in p["field_r"]:
            for _ in range(int(p["field_dirs"])):
                n = rng.normal(size=3)
                n /= np.linalg.norm(n)
                x = np.concatenate([[float(t)], float(r) * n])
                pts.append(x)
                rows.append((float(r), fld.U(x)))
    # frame check on a few points: hbar from boosted coordinates against the transformed tensor
    bc = float(p["beta_check"])
    B = boost_matrix(bc)
    Binv = np.linalg.inv(B)
    fb = make_field(s, lam, p, B)
    dev, gm_meas = 0.0, []
    for x in pts[:: max(1, len(pts) // 12)]:
        hb = fld.hbar(x)
        hb_b = fb.hbar(B @ x)
        hb_t = Binv.T @ hb @ Binv                                # lower-index transformation
        dev = max(dev, float(np.max(np.abs(hb_b - hb_t)) / max(np.max(np.abs(hb)), 1e-300)))
        gm_meas.append(hb_b[0, 1] / hb_b[0, 0])
    # noise of the discrete gradient: radial derivative of the discrete U against the continuum dU_c/dr
    s0 = float(p["s0"])
    grad_err = []
    for x in pts[:: max(1, len(pts) // 24)]:
        r = float(np.linalg.norm(x[1:]))
        _, dU, _, _ = fld.derivs(x)
        d_disc = float(dU[1:] @ (x[1:] / r))                     # dU/dr of the discrete field
        d_cont = (U_continuum(r * 1.001, lam, 1.0, s0) - U_continuum(r * 0.999, lam, 1.0, s0)) / (0.002 * r)
        grad_err.append((d_disc, d_cont))
    ub = B @ np.array([1.0, 0.0, 0.0, 0.0])
    ub_low = ub * np.array([1.0, -1.0, -1.0, -1.0])
    _, _, h, _ = fld.derivs(pts[len(pts) // 2])
    # flux per unit time over the retarded times the field points see (t - r): the flux per unit time of a bundle
    # fluctuates by ~20 % between time slabs (the event rate of its few chains), so it is measured where it acts
    t_lo = min(float(t) for t in p["field_t"]) - max(float(r) for r in p["field_r"]) - 1.0
    t_hi = max(float(t) for t in p["field_t"]) - min(float(r) for r in p["field_r"]) + 1.0
    m = (s["flux"] > 0.0) & (s["X"][:, 0] >= t_lo) & (s["X"][:, 0] <= t_hi)
    mu_ret = float(s["flux"][m].sum() / (t_hi - t_lo))
    grad_rel = [abs(dd - dc * mu_ret) / abs(dc * mu_ret) for dd, dc in grad_err]
    return {"rows": rows, "mu": mu_ret, "width": s["width"], "frame_dev": dev, "grad_rel": float(np.median(grad_rel)),
            "gm_meas": float(np.mean(gm_meas)), "gm_pred": float(ub_low[1] / ub_low[0]),
            "gamma_field": float(np.mean([h[i, i] for i in (1, 2, 3)]) / h[0, 0])}


def probe_task(p, seed, K, rule):
    """Part probe: for every speed, newton and the rule on the same probe sprinkling and field; frame check."""
    fld = PointField(K)
    bc = float(p["beta_check"])
    B = boost_matrix(bc)
    fld_b = PointField(K, B)
    out = {"per_v": {}}
    for iv, v in enumerate(p["speeds"]):
        v = float(v)
        gam = 1.0 / math.sqrt(1.0 - v * v)
        x_half, b, ts = float(p["x_half"]), float(p["b"]), float(p["t_start"])
        dur_lab = 2.0 * x_half / v
        # probe tube in its rest frame, then boosted to the lab and moved to the start of the path
        Xr, ptr, idx, cptr, cidx = tube_rest(p, derive_seed(seed, "probe", iv), float(p["R_p"]),
                                             (dur_lab + 2.0) / gam + 1.0, -1.0)
        srcp = ((Xr[:, 0] >= 0.0) & (Xr[:, 0] < float(p["t_emit_p"]))
                & ((Xr[:, 1:] ** 2).sum(axis=1) <= float(p["r_emit_p"]) ** 2)).astype(float)
        X = boost_points(Xr, -v) + np.array([ts, -x_half, b, 0.0])   # lab coordinates (order and links unchanged)
        u_src = np.array([gam, gam * v, 0.0, 0.0])
        rn = probe_transport(X, ptr, idx, cptr, cidx, srcp, u_src, "newton", fld)
        rr = probe_transport(X, ptr, idx, cptr, cidx, srcp, u_src, rule, fld)
        # readout: flux-weighted u^y at the end of the path (lab time window), and the flux that reached it
        t_end = ts + dur_lab
        w_end = (X[:, 0] >= t_end - float(p["end_window"])) & (X[:, 0] <= t_end)
        t_b = ts + gam * float(p["t_emit_p"]) + 1.0               # first window after the emission has ended
        w_beg = (X[:, 0] >= t_b) & (X[:, 0] <= t_b + float(p["end_window"]))

        def uy(r_):
            f_ = r_["flux"][w_end]
            return float((f_ * r_["dir"][w_end, 2]).sum() / f_.sum()) if f_.sum() > 0 else float("nan")
        f_end, f_beg = float(rr["flux"][w_end].sum()), float(rr["flux"][w_beg].sum())
        # frame check: the same probe (coordinates and directions boosted by beta') in the boosted field
        rb = probe_transport(X @ B.T, ptr, idx, cptr, cidx, srcp, B @ u_src, rule, fld_b)
        frame = first_difference(rr, rb, lambda D_: D_ @ np.linalg.inv(B).T)
        out["per_v"][v] = {"uy_newton": uy(rn), "uy_rule": uy(rr), "flux_ratio": f_end / f_beg if f_beg > 0 else 0.0,
                           "frame": frame, "N": len(X)}
    return out


# ============================================================================================== run

def run(p, ctx):
    mode = str(p["mode"])
    t0 = time.time()
    n, Bn, nse = int(p["n_seeds"]), int(p["boot_B"]), float(p["n_se"])
    seeds = [derive_seed(ctx.seed, "F7.4", mode, k) for k in range(n)]
    if mode == "field":
        cal = parallel_map(calibration_task, [(p, derive_seed(ctx.seed, "F7.4-cal", k)) for k in range(int(p["n_cal"]))],
                           ctx.workers, ctx.progress,
                           lambda i, r_, dd, nn: f"[F7.4] calibration {dd}/{nn} ({time.time()-t0:.0f}s)")
        lam = float(p["U_target"]) / float(np.mean(cal))
        lines = [f"3+1, rho = 1; source bundle at rest (tube R_s = {p['R_s']:g}, emission radius {p['r_emit']:g}, t < "
                 f"{p['t_emit']:g}), field hbar_ab = -4 lambda sum F u_a u_b g(sigma), g = sigma exp(-sigma/s0)/(2 pi s0^2), "
                 f"s0 = {p['s0']:g}; lambda = {lam:.4e} (U = {p['U_target']} at r = b = {p['b']:g})"]
        res = parallel_map(field_task, [(p, sd, lam) for sd in seeds], ctx.workers, ctx.progress,
                           lambda i, r_, dd, nn: f"[F7.4] field sprinkling {dd}/{nn} ({time.time()-t0:.0f}s)")
        rs = np.array(p["field_r"], float)
        # per sprinkling: mean U at each r, the exponent and the coefficient against lambda mu / (4 pi r)
        s0 = float(p["s0"])
        shape = np.array([U_continuum(r0, 1.0, 1.0, s0) for r0 in rs])     # U_c(r) per unit lambda mu
        k_c = float(-np.polyfit(np.log(rs), np.log(shape), 1)[0])          # exponent of the continuum limit (about 1)
        ks, Cs = [], []
        for r_ in res:
            Ur = np.array([np.mean([u for rr_, u in r_["rows"] if rr_ == r0]) for r0 in rs])
            ks.append(float(-np.polyfit(np.log(rs), np.log(Ur), 1)[0]) if np.all(Ur > 0) else float("nan"))
            Cs.append(float(np.mean(Ur / (lam * r_["mu"] * shape))))       # measured / continuum, per sprinkling
        k_m, k_se = boot_se(ks, Bn, derive_seed(ctx.seed, "F7.4-k"))[:2]
        C_m, C_se = boot_se(Cs, Bn, derive_seed(ctx.seed, "F7.4-C"))[:2]
        lines.append(f"  U against r = {', '.join(f'{x:g}' for x in rs)}: exponent k = {k_m:.3f} +- {k_se:.3f} (continuum "
                     f"limit of the construction k_c = {k_c:.3f}; bare Newton 1), coefficient C = U / U_c = {C_m:.3f} +- "
                     f"{C_se:.3f} (pred 1); U_c / Newton at these r: "
                     + ", ".join(f"{sh * 4 * math.pi * r0:.3f}" for sh, r0 in zip(shape, rs))
                     + f"; flux per unit time mu = {np.mean([r_['mu'] for r_ in res]):.1f}")
        lines.append(f"  gamma_field = mean h_ii / h_00 = {np.mean([r_['gamma_field'] for r_ in res]):.4f} (1 by construction; "
                     f"a scalar field gives 0); gravitomagnetic ratio hbar'_0x / hbar'_00 in the frame beta' = "
                     f"{p['beta_check']}: {np.mean([r_['gm_meas'] for r_ in res]):+.4f} against u'_x / u'_0 = "
                     f"{res[0]['gm_pred']:+.4f} (a scalar field: 0)")
        inst = {"frame check of the tensor (relative 1e-9)": max(r_["frame_dev"] for r_ in res) <= 1e-9,
                "field points outside the source (rms radius < r_min / 2)":
                    max(r_["width"] for r_ in res) < 0.5 * min(float(x) for x in p["field_r"])}
        lines.append(f"  largest frame deviation {max(r_['frame_dev'] for r_ in res):.1e}; source width "
                     f"{np.mean([r_['width'] for r_ in res]):.2f} (max {max(r_['width'] for r_ in res):.2f}); median relative "
                     f"noise of the discrete gradient dU/dr at the field points "
                     f"{np.median([r_['grad_rel'] for r_ in res]):.2f} (the values themselves: C above)")
        status = "UNDECIDED"
        if abs(k_m - k_c) + nse * k_se <= 0.25 and abs(C_m - 1) + nse * C_se <= 0.25:
            status = "PASSED"
        elif abs(k_m - k_c) - nse * k_se > 0.25 or abs(C_m - 1) - nse * C_se > 0.25:
            status = "FAILED"
        data = {"lambda": lam, "k": [k_m, k_se], "k_c": k_c, "C": [C_m, C_se], "instrument": inst}
    else:
        rule = str(p["rule"])
        K = float(p["U_target"]) * float(p["b"])                 # U = U_target at the impact parameter
        lam = K
        lines = [f"3+1, rho = 1; probes on causal-set tubes (radius {p['R_p']:g}) in the continuum tensor field of a static "
                 f"point source, U = K / r, K = {K:g} (U = {p['U_target']} at b = {p['b']:g}), h_ab = -2U delta_ab; path x "
                 f"from -{p['x_half']:g} to +{p['x_half']:g}, start t = {p['t_start']:g}"]
        res = parallel_map(probe_task, [(p, sd, K, rule) for sd in seeds], ctx.workers, ctx.progress,
                           lambda i, r_, dd, nn: f"[F7.4] probe ({rule}) sprinkling {dd}/{nn} ({time.time()-t0:.0f}s)")
        pred = (lambda v: 1.0 + v * v) if rule == "tensor" else (lambda v: 1.0 - v * v)
        rng = np.random.default_rng(derive_seed(ctx.seed, "F7.4-boot", rule))
        ok_all, fail_any, fs = True, False, {}
        for v in [float(x) for x in p["speeds"]]:
            un = np.array([r_["per_v"][v]["uy_newton"] for r_ in res])
            ur = np.array([r_["per_v"][v]["uy_rule"] for r_ in res])
            f_m = float(ur.mean() / un.mean())
            bs = []
            for _ in range(Bn):
                j = rng.integers(0, n, n)
                bs.append(ur[j].mean() / un[j].mean())
            f_se = float(np.std(bs, ddof=1))
            tol = 0.5 * v * v
            fs[v] = (f_m, f_se)
            ok_all &= abs(f_m - pred(v)) + nse * f_se <= tol
            fail_any |= abs(f_m - pred(v)) - nse * f_se > tol
            lines.append(f"  v = {v:g}: deflection newton {un.mean():+.5f} (u^y at the end), {rule} {ur.mean():+.5f}; "
                         f"f = {f_m:.4f} +- {f_se:.4f}, predicted {pred(v):.3f} (tolerance {tol:.3f}); probe flux at the end "
                         f"{np.mean([r_['per_v'][v]['flux_ratio'] for r_ in res]):.3f} of the start")
        vs = np.array(sorted(fs))
        fit = np.polyfit(vs ** 2, [fs[v][0] for v in vs], 1)
        lines.append(f"  light limit (no verdict): linear fit of f in v^2 gives f(1) = {np.polyval(fit, 1.0):.3f} "
                     f"(tensor: 2, scalar: 0; Newton: 1)")
        frames = [r_["per_v"][v]["frame"] for r_ in res for v in r_["per_v"]]
        inst = {f"frame check beta' = {p['beta_check']} ({rule})": all(f_ in ("identical", "tie") for f_ in frames),
                f"probe flux at the end >= {p['flux_min']}": all(r_["per_v"][v]["flux_ratio"] >= float(p["flux_min"])
                                                                for r_ in res for v in r_["per_v"])}
        lines.append(f"  numerical ties in the frame check: {sum(f_ == 'tie' for f_ in frames)} of {len(frames)}")
        status = "PASSED" if ok_all else ("FAILED" if fail_any else "UNDECIDED")
        data = {"lambda": lam, "f": {str(k_): list(v_) for k_, v_ in fs.items()}, "instrument": inst}
    if not all(inst.values()):
        status = "INSTRUMENT_FAILED"
    lines.append("instrument: " + "; ".join(f"{k_}: {'OK' if v_ else 'FAILS'}" for k_, v_ in inst.items()))
    return Result(status, lines, data)
