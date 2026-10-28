# Causal Gravity: an interactive simulator of the gravity law of causal physics
# ===========================================================================
# A stand-alone interactive application; it does not change the code of the causality laboratory.
# Masses, geometry and boundary conditions are inputs of the simulator.
#
# Run:
#   python demo/causal_gravity.py                       # opens the default scene
#   python demo/causal_gravity.py --selftest            # runs the numerical self-tests (no window)
#   python demo/causal_gravity.py --scene 2             # opens the scene with the given index
#   python demo/causal_gravity.py --scene 2 --autorun   # and starts the simulation at once
# Requirements: Python 3.9+, tkinter and numpy.
#
# THE GRAVITY LAW (state of the theory, 28 Oct 2026)
# --------------------------------------------------
# Causal gravity is the gravity law of the laboratory after decision 13 (route A):
#   * the principle of conserved causal throughput (where the light cones narrow, the levels of the cone saturate)
#     gives the PPN parameter gamma = 1 in 3+1 (F7.1c);
#   * in the linear order the principle is the same thing as coupling to the full energy-momentum tensor:
#     hbar_ab = 4U u_a u_b for dust (docs/derivations/tenzorova_vazba.md); the substrate realises this coupling
#     with the rule "every flow narrows the cone in its own rest frame, capacity deficits add" (F7.4a PASSED:
#     counter-streaming beams gravitate as Gamma^2 (1 + v^2), the cones tilt along the flow);
#   * the nonlinear completion is general relativity (decision 13), to be derived on the substrate later
#     (Jacobson's route).
# So in every scene of this simulator causal gravity and Einstein's general relativity give the same numbers;
# that is the expected result, not a coincidence. What the simulator adds is the comparison with the frameworks
# that fail (Newton, the scalar stage of the theory, the fixed order of the F3 rule) and the measurement of each
# classical test in each framework.
#
# WHAT IS COMPUTED (approximations of the simulator)
# --------------------------------------------------
#   * Field: U = sum_j M_j / sqrt(|x - x_j|^2 + eps^2) (G = c = 1), linear superposition, quasi-static (no
#     retardation). The metric of causal gravity and of GR is the Schwarzschild metric in isotropic coordinates
#     with this U: ds^2 = -A(U) dt^2 + B(U) dx^2, A = ((1 - U/2)/(1 + U/2))^2, B = (1 + U/2)^4 (exact for one
#     static mass; for several masses correct to the first post-Newtonian order).
#   * Motion: geodesics of the static metric, every body as a test body in the field of the others:
#       dv/dt = -(A'/2B) grad U - v (v.grad U) B'/B + |v|^2 (B'/2B) grad U + v (v.grad U) A'/A.
#   * Frame dragging: spinning bodies (spin J, geometric units m^2) carry g_0i = -2 (J x r)_i / r^3. Bodies in the
#     plane feel a = v x (curl g_0i); gyroscopes precess with Omega = (gamma + 1/2) v x grad U (geodetic) plus
#     ((1 + gamma)/2) [3 (J.n) n - J] / r^3 (frame dragging). Mass currents of moving bodies are not included.
#   * Gravitational radiation: the 2.5PN radiation-reaction acceleration of a binary (harmonic gauge),
#     a = (8/5) eta M^2 / r^3 [rdot n (3 v^2 + 17 M / 3r) - v (v^2 + 3 M / r)], times the radiation factor of the
#     framework (GR and causal 1, scalar frameworks 1/6 for circular orbits, Newton 0).
#   * Moving matter as a source (F7.4a): straight counter-streaming beams act with the active line density
#     lambda * f(v), f = Gamma^2 (1 + v^2) for the tensor coupling, 1 for the scalar stage and Newton.
#   * Strong field: outside the black-hole scene, a body that enters |U| > 0.2 of another source merges with it
#     (the superposition of point fields fails there); photons are absorbed.
#
# COMPARISON FRAMEWORKS (the same scene integrated in parallel, dashed trails)
# ---------------------------------------------------------------------------
#   GR       general relativity (Einstein): Schwarzschild metric, frame dragging, quadrupole radiation;
#   N        Newton (light as a corpuscle at speed c): no clock effects, half the light bending, no precession;
#   H1       scalar stage of the theory: the field Pi in a fixed frame with the exponential metric H1
#            (A = e^-2U, B = e^2U): light bending and perihelion correct, no frame dragging, scalar radiation
#            (1/6 of the quadrupole power for circular orbits), the black hole has a throat and no horizon;
#   F3       fixed order (transport rule F3, Pi^2 eta, gamma = -1): no light bending, probes bend with 1 - v^2;
#   PPN      parametrised post-Newtonian metric with the beta, gamma sliders.
#
# SCENES
# ------
#  0 Perihelion precession                 11 Black hole (strong field)
#  1 Solar system (exaggerated)            12 Alpha Centauri A/B
#  2 Light bending                         13 Sirius A/B (white dwarf)
#  3 Shapiro delay                         14 Solar neighbourhood to 12 ly (3D)
#  4 Binary star and planet                15 Free scene
#  5 Clocks in a field                     16 Fast probes and light (1 + gamma v^2)
#  6 Solar system (real, J2000)            17 Frame dragging (exaggerated)
#  7 Earth surface: apple and clocks       18 Gravity Probe B (gyroscope)
#  8 Earth: ISS and GPS                    19 Binary pulsar PSR B1913+16
#  9 UAP: scripted manoeuvre (hypothesis)  20 Inspiral (strong radiation)
# 10 UAP: keyboard control                 21 What gravitates: counter-streaming beams (F7.4a)
#
# CONTROLS
# --------
# Bodies, photons and clocks can be added or deleted. A click places the object, a drag sets direction and
# speed. Masses of new bodies are entered in scientific notation with a unit (kg, Earth or solar masses in the SI
# scenes, model units G = c = 1 in the geometric scenes), e.g. 5.972e24 kg. The pace is simulated time per second
# of real time. The 3D view of planar scenes is only a projection of the plane z = 0.
#
# CLASSICAL TESTS
# ---------------
# The button "Run classical tests" integrates, for every framework, light bending at the Sun, probe bending,
# Mercury's perihelion, the Shapiro delay, the Pound-Rebka redshift, Gravity Probe B, the binary-pulsar period
# decay, the active mass of moving matter and the black-hole shadow, and prints a table with the observations.
#
# SELF-TESTS (--selftest) check the implementation of the formulas, not the validity of the theory.
# ===========================================================================

from __future__ import annotations

import math
import sys
import time
from collections import deque

import numpy as np

APP_NAME = "Causal Gravity"
APP_VERSION = "2.0 (28 Oct 2026)"

EPS = 1e-3  # numerical softening of point sources [length units]
C_SI = 299792458.0  # speed of light [m/s]; in the SI scenes time is stored internally as ct [m]
G_C2 = 7.42591e-28  # G/c^2 [m/kg]: converts kilograms to the geometric length unit
M_EARTH = 4.43502e-3  # G M_Earth / c^2 [m]
M_EARTH_KG = 5.9722e24  # mass of the Earth [kg]
R_EARTH = 6.371e6  # mean radius of the Earth [m]
G_SURF = M_EARTH / R_EARTH**2 * C_SI**2  # surface gravity [m/s^2]
# Spin of the Earth in geometric units: G I omega / c^3 with I = 8.034e37 kg m^2, omega = 7.2921e-5 rad/s.
J_EARTH = 6.6743e-11 * 8.034e37 * 7.2921159e-5 / C_SI**3  # [m^2], ~0.0145
GPB_A = 7.0274e6  # semi-major axis of the orbit of Gravity Probe B [m] (altitude ~642 km)
AU = 1.495978707e11  # astronomical unit [m]
YEAR = 365.25 * 86400.0  # Julian year [s]
M_SUN = 1476.625  # G M_sun / c^2 [m]
M_SUN_KG = 1.98847e30  # mass of the Sun [kg]
R_SUN = 6.957e8  # nominal solar radius [m]
ARCSEC = math.pi / (180.0 * 3600.0)  # arcsecond [rad]
MAS = ARCSEC / 1000.0  # milliarcsecond [rad]
LY = 9.4607304725808e15  # light year [m]
TAU_NULL = 0.05 * C_SI  # "null field" hypothesis: velocity response time 0.05 s, stored in metres
L_MAX = 0.2  # weak-field limit; for |U| > L_MAX the superposition of point fields is not trusted

# Units offered for the mass of a new body: label -> geometric mass per unit [m].
MASS_UNITS_SI = {"kg": G_C2, "Earth masses (M⊕)": M_EARTH, "solar masses (M☉)": M_SUN}
MASS_UNITS_GEO = {"model units (G = c = 1)": 1.0}


def fmt_len(m):
    """A length in metres in a convenient display unit."""
    a = abs(m)
    if a >= 0.1 * LY:
        return f"{m / LY:.4g} ly"
    if a >= 0.01 * AU:
        return f"{m / AU:.4g} AU"
    if a >= 1000.0:
        return f"{m / 1000.0:.4g} km"
    return f"{m:.4g} m"


def fmt_time(sec):
    """A time in seconds as seconds, minutes, hours, days or years."""
    a = abs(sec)
    if a < 120:
        return f"{sec:.2f} s"
    if a < 7200:
        return f"{sec / 60:.2f} min"
    if a < 2 * 86400:
        return f"{sec / 3600:.2f} h"
    if a < 2 * YEAR:
        return f"{sec / 86400:.2f} d"
    return f"{sec / YEAR:.3f} yr"


def fmt_angle(rad):
    """An angle in degrees, or in arcseconds / milliarcseconds when it is small."""
    a = abs(rad)
    if a >= math.radians(0.01):
        return f"{math.degrees(rad):.4f}°"
    if a >= ARCSEC * 0.01:
        return f"{rad / ARCSEC:.4f}″"
    return f"{rad / MAS:.4f} mas"


def fmt_sci(x, digits=4):
    """Scientific notation used for every mass in the simulator, e.g. 5.9720e+24."""
    return f"{x:.{digits}e}"


def parse_sci(text):
    """
    A positive finite number in scientific notation from user text. Accepts 5.972e24, 5,972E+24, 5.972x10^24,
    5.972*10^24, 5.972·10^24 and 5.972×10^24 (spaces ignored). Returns None for anything else.
    """
    s = str(text).strip().replace(" ", "").replace(",", ".")
    for sep in ("×10^", "x10^", "X10^", "*10^", "·10^"):
        s = s.replace(sep, "e")
    s = s.replace("−", "-")
    try:
        v = float(s)
    except ValueError:
        return None
    return v if math.isfinite(v) and v > 0 else None


# =============================================================================
# Physics: bodies, gravity laws (frameworks)
# =============================================================================


class Body:
    """
    A body of a scene. Kinds:
      "mass"   a massive body (a source if mass > 0); probe = True marks a fast test probe, gyro a gyroscope;
      "photon" a massless particle;
      "clock"  a static clock (fixed);
      "craft"  a body with a scripted motion (UAP scene; its proper acceleration is the thrust it would need);
      "pilot"  a body steered from the keyboard;
      "rod"    a straight pair of counter-streaming beams along x (F7.4a): rest line density lam [m per m], half
               length half [m], beam speed vb [c]; light = True makes it a pair of light beams of energy density lam.
    """

    def __init__(self, name, x, y, vx=0.0, vy=0.0, mass=0.0, kind="mass", fixed=False, color="#ffffff",
                 radius=0.0):
        self.name, self.kind, self.color = name, kind, color
        self.radius = float(radius)  # a positive radius marks a solid surface; other bodies land on it
        self.z = 0.0  # third coordinate, used by the real 3D star map
        self.thrust = (0.0, 0.0)  # steered acceleration [m/s^2] (pilot)
        self.info = ""  # extra description (star map, binaries)
        self.landed_t = None  # simulation time of landing (Earth scenes)
        self.script = None  # scripted motion (UAP scene), see CraftScript
        self.max_g = 0.0  # craft/pilot: largest proper acceleration so far [m/s^2]
        self.now_g = 0.0
        self.mass = float(mass) if kind == "mass" else 0.0
        self.fixed = bool(fixed) or kind in ("clock", "rod")
        self.x, self.y, self.vx, self.vy = float(x), float(y), float(vx), float(vy)
        self.tau = 0.0  # accumulated proper time
        self.tau_def = 0.0  # accumulated t - tau; keeps precision in very weak fields
        self.t_birth = None  # time the body entered the scene (set by the application)
        self.trail = deque(maxlen=4000)
        self.alive = True
        self.J = (0.0, 0.0, 0.0)  # spin angular momentum G J / c^3 [m^2] (frame dragging)
        self.probe = False  # fast test probe: measured by the deflection tracker
        self.gyro = None  # gyroscope state, see Sim.update_gyros
        self.ground = False  # clock on the ground (reference clock of the Earth scenes)
        self.captured = False  # photon or body captured by a black hole
        self.lam, self.half, self.vb, self.light = 0.0, 0.0, 0.0, False  # rod parameters


def make_gyro(body, s0=(1.0, 0.0)):
    """Turns a body into a gyroscope whose spin points along the unit vector s0 in the plane (guide star)."""
    n = math.hypot(*s0)
    body.gyro = {"s0": (s0[0] / n, s0[1] / n), "geo": 0.0, "fd": 0.0, "fd_z": 0.0, "t": 0.0,
                 "turn": 0.0, "ang": None, "snaps": []}  # snaps: (t, geo, fd, fd_z) after every whole orbit
    return body


# ---------------------------------------------------------------- scalar stage: exact vacuum profile of H1

_BH = {}


def bh_profile():
    """
    Exterior vacuum solution of the scalar-stage field equation B_g ln Pi = 0 with the metric H1 in 3+1 for a
    spherical source M = 1:  e^(2L) [(1 - L) Laplace L + L |grad L|^2] = 0. With s = ln r: dL/ds = P,
    dP/ds = -P - L P^2 / (1 - L), P = r dL/dr. Integrated by RK4 from r = 1e5 to r = 1e-4 (ds = 1e-3) from
    L = -1/r + 1/(6 r^3). The closed form is Ei(1 - L) = Ei(1) + e M / r (docs/derivations/gravitacni_rovnice.md).
    Returns (s, L, P, info) with the photon sphere (r dL/dr = 1/2), the critical impact parameter b = r e^(-2L) and
    the throat (r dL/dr = 1). Decision 13: these strong-field properties are NOT predictions of the theory; they
    describe the scalar stage only.
    """
    if _BH:
        return _BH["s"], _BH["L"], _BH["P"], _BH["info"]

    def rhs(L, P):
        return P, -P - L * P * P / (1.0 - L)

    s0, s1, h = math.log(1e5), math.log(1e-4), -1e-3
    r0 = 1e5
    L, P = -1.0 / r0 + 1.0 / (6 * r0**3), 1.0 / r0 - 1.0 / (2 * r0**3)
    ss, LL, PP = [s0], [L], [P]
    s_ = s0
    while s_ > s1:
        k1 = rhs(L, P)
        k2 = rhs(L + 0.5 * h * k1[0], P + 0.5 * h * k1[1])
        k3 = rhs(L + 0.5 * h * k2[0], P + 0.5 * h * k2[1])
        k4 = rhs(L + h * k3[0], P + h * k3[1])
        L += h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        P += h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        s_ += h
        ss.append(s_), LL.append(L), PP.append(P)
    s_arr, L_arr, P_arr = np.array(ss[::-1]), np.array(LL[::-1]), np.array(PP[::-1])  # increasing s

    def crossing(level):
        k = np.where((P_arr[:-1] - level) * (P_arr[1:] - level) <= 0)[0][-1]  # outermost crossing
        f = (level - P_arr[k]) / (P_arr[k + 1] - P_arr[k])
        return math.exp(s_arr[k] + f * (s_arr[k + 1] - s_arr[k])), L_arr[k] + f * (L_arr[k + 1] - L_arr[k])

    r_ph, L_ph = crossing(0.5)
    r_th, L_th = crossing(1.0)
    info = {"r_ph": r_ph, "R_ph": r_ph * math.exp(-L_ph), "b_crit": r_ph * math.exp(-2 * L_ph),
            "r_throat": r_th, "R_throat": r_th * math.exp(-L_th)}
    _BH.update(s=s_arr, L=L_arr, P=P_arr, info=info)
    return s_arr, L_arr, P_arr, info


def bh_L(r, M=1.0):
    """L and dL/dr of the scalar-stage profile at radii r for mass M; logarithmic tail below 1e-4 M."""
    s_arr, L_arr, P_arr, _ = bh_profile()
    r = np.maximum(np.asarray(r, float), 1e-300)
    sx = np.log(r / M)
    L = np.interp(sx, s_arr, L_arr)
    P = np.interp(sx, s_arr, P_arr)
    lo = sx < s_arr[0]
    L = np.where(lo, L_arr[0] + P_arr[0] * (sx - s_arr[0]), L)
    hi = sx > s_arr[-1]
    L = np.where(hi, -M / r, L)
    P = np.where(hi, M / r, P)
    return L, P / r


# ---------------------------------------------------------------- frameworks

def lorentz(v):
    """Lorentz factor Gamma_v = 1 / sqrt(1 - v^2) (not to be confused with the PPN gamma)."""
    return 1.0 / math.sqrt(max(1.0 - v * v, 1e-300))


class Law:
    """
    A gravity law (framework). U = -L = sum M / r is the Newtonian potential (G = c = 1).
    kind "metric": geodesics of the static metric ds^2 = -A(U) dt^2 + B(U) dx^2; light speed sqrt(A/B), static clocks
      dtau/dt = sqrt(A); A - 1 is computed separately (Am1) to keep clock differences of 1e-15.
    kind "newton": dv/dt = grad U for bodies and photons, clocks feel nothing.
    beta, gamma, eps: post-Newtonian parameters of g_tt = 1 - 2U + 2 beta U^2 and g_ij = 1 + 2 gamma U + (3/2) eps U^2;
    gm: factor of frame dragging (0 = none); rad: factor of the quadrupole radiation reaction; src(v): active mass of
    counter-streaming matter at speed v relative to rest (F7.4a); light_src: active mass of a light beam per unit
    energy; bh: strong-field description shown in the readouts.
    """

    def __init__(self, key, name, color, kind="metric", A=None, dA=None, B=None, dB=None, Am1=None,
                 beta=1.0, gamma=1.0, eps=1.0, gm=0.0, rad=0.0, src=None, light_src=0.0, short=None):
        self.key, self.name, self.color, self.kind = key, name, color, kind
        self.A, self.dA, self.B, self.dB, self.Am1 = A, dA, B, dB, Am1
        self.beta, self.gamma, self.eps = beta, gamma, eps
        self.gm, self.rad, self.light_src = gm, rad, light_src
        self.src = src if src is not None else (lambda v: 1.0)
        self.short = short or key

    def static_rate_m1(self, U):
        """dtau/dt - 1 of a static clock at potential U (numerically stable for |U| ~ 1e-10)."""
        if self.kind == "newton":
            return 0.0 * np.asarray(U, float)
        return np.expm1(0.5 * np.log1p(np.maximum(self.Am1(U), -1.0 + 1e-300)))

    def deflection2(self):
        """Coefficient of the second-order light bending, pi (8 + 8 gamma - 4 beta + 3 eps) / 4 (times (M/b)^2)."""
        if self.kind == "newton":
            return 0.0
        return math.pi * (8 + 8 * self.gamma - 4 * self.beta + 3 * self.eps) / 4.0

    def perihelion_factor(self):
        """Perihelion advance relative to GR, (2 + 2 gamma - beta) / 3; Newton 0."""
        return 0.0 if self.kind == "newton" else (2 + 2 * self.gamma - self.beta) / 3.0


def _schwarzschild_iso():
    """Schwarzschild metric in isotropic coordinates: A = ((1 - U/2)/(1 + U/2))^2, B = (1 + U/2)^4."""
    return dict(A=lambda U: ((1 - U / 2) / (1 + U / 2)) ** 2,
                dA=lambda U: -2 * (1 - U / 2) / (1 + U / 2) ** 3,
                B=lambda U: (1 + U / 2) ** 4,
                dB=lambda U: 2 * (1 + U / 2) ** 3,
                Am1=lambda U: -2 * U / (1 + U / 2) ** 2)


def _tensor_src(v):
    """Active mass of two counter-streaming beams relative to rest, T_00 + T_ii: Gamma_v^2 (1 + v^2) (F7.4a)."""
    return lorentz(v) ** 2 * (1.0 + v * v)


def _law_cg():
    """Causal gravity: tensor coupling (F7.4a) with general relativity as the nonlinear completion (decision 13)."""
    return Law("CG", "causal gravity (tensor coupling)", "#ffffff", beta=1.0, gamma=1.0, eps=1.0, gm=1.0, rad=1.0,
               src=_tensor_src, light_src=2.0, short="Causal", **_schwarzschild_iso())


def _law_gr():
    """General relativity (Einstein): Schwarzschild metric (exact for one mass, 1PN-correct when superposed)."""
    return Law("GR", "general relativity (Einstein)", "#ff7ad9", beta=1.0, gamma=1.0, eps=1.0, gm=1.0, rad=1.0,
               src=_tensor_src, light_src=2.0, short="Einstein", **_schwarzschild_iso())


def _law_newton():
    """Newton: light as a corpuscle at speed c, clocks unaffected, no radiation, rest mass gravitates."""
    return Law("N", "Newton (light as a corpuscle)", "#7aff8a", kind="newton", beta=0.0, gamma=0.0, eps=0.0,
               gm=0.0, rad=0.0, light_src=1.0, short="Newton")


def _law_h1():
    """Scalar stage: field Pi in a fixed frame, metric H1, A = e^(2L) = e^(-2U), B = e^(-2L) = e^(2U)."""
    return Law("H1", "scalar stage (Π in a fixed frame, H1)", "#c3a6ff", beta=1.0, gamma=1.0, eps=4.0 / 3.0,
               gm=0.0, rad=1.0 / 6.0, light_src=0.0, short="Scalar H1",
               A=lambda U: np.exp(-2 * U), dA=lambda U: -2 * np.exp(-2 * U),
               B=lambda U: np.exp(2 * U), dB=lambda U: 2 * np.exp(2 * U), Am1=lambda U: np.expm1(-2 * U))


def _law_f3():
    """Fixed order: the transport rule F3 follows the geodesics of Pi^2 eta (gamma = -1), A = B = e^(-2U)."""
    return Law("F3", "fixed order (F3, γ = −1)", "#ffb37a", beta=1.0, gamma=-1.0, eps=4.0 / 3.0, gm=0.0,
               rad=1.0 / 6.0, light_src=0.0, short="F3",
               A=lambda U: np.exp(-2 * U), dA=lambda U: -2 * np.exp(-2 * U),
               B=lambda U: np.exp(-2 * U), dB=lambda U: -2 * np.exp(-2 * U), Am1=lambda U: np.expm1(-2 * U))


def _law_ppn(beta, gamma):
    """Parametrised post-Newtonian metric to first order: A = 1 - 2U + 2 beta U^2, B = 1 + 2 gamma U; frame
    dragging (1 + gamma)/2 of GR; radiation as GR (the PPN metric does not fix it); moving matter
    Gamma^2 (1 + gamma v^2)."""
    return Law("PPN", f"PPN β = {beta:g}, γ = {gamma:g}", "#7ad0ff", beta=beta, gamma=gamma, eps=0.0,
               gm=(1.0 + gamma) / 2.0, rad=1.0, light_src=1.0 + gamma, short="PPN",
               src=lambda v: lorentz(v) ** 2 * (1.0 + gamma * v * v),
               A=lambda U: 1 - 2 * U + 2 * beta * U * U, dA=lambda U: -2 + 4 * beta * U,
               B=lambda U: 1 + 2 * gamma * U, dB=lambda U: 2 * gamma + 0 * U,
               Am1=lambda U: -2 * U + 2 * beta * U * U)


MAIN_LAW = _law_cg


# =============================================================================
# Physics: field, geodesic acceleration, frame dragging, radiation and the integrator
# =============================================================================


class Sim:
    """A system of bodies and its integrator without graphics, shared by the application and the self-tests."""

    def __init__(self, law=None):
        self.law = law if law is not None else MAIN_LAW()
        self.profile = None  # source with the exact scalar-stage profile (black-hole scene, framework H1 only)
        self.no_guard = False  # True switches off merging/absorption in the strong field (black-hole scene)
        self.capture = None  # (centre body, radius): photons and bodies closer than radius are captured
        self.radiation = False  # True applies the radiation reaction between massive moving bodies
        self.bodies = []
        self.t = 0.0
        self.eta = 0.02  # safety factor of the time step
        self.dt_max = 0.5
        self.events = []  # records of merging and absorption
        self.rad_E = {}  # radiated energy per pair of bodies (names) [geometric units]
        self._extra = None
        self._vt = []
        self._spins = []
        self._spin_ids = []

    def remove(self, body):
        """Removes a body. Removing the profile source (scalar-stage black hole) also removes its exact field."""
        body.alive = False
        if body in self.bodies:
            self.bodies.remove(body)
        if body is self.profile or getattr(body, "is_profile", False):
            self.profile = None
            self.no_guard = False
            self.capture = None

    # ---- state arrays
    def arrays(self):
        live = [b for b in self.bodies if b.alive and b.kind != "rod"]
        pos = np.array([[b.x, b.y] for b in live], float).reshape(-1, 2)
        vel = np.array([[b.vx, b.vy] for b in live], float).reshape(-1, 2)
        # the field of the profile source is added by ext_field, so it is not counted among the point masses
        mass = np.array([0.0 if b is self.profile else b.mass for b in live], float)
        photon = np.array([b.kind == "photon" for b in live], bool)
        fixed = np.array([b.fixed or b.kind == "craft" for b in live], bool)  # scripted bodies follow the script
        return live, pos, vel, mass, photon, fixed

    def rods(self):
        return [b for b in self.bodies if b.alive and b.kind == "rod"]

    def spins(self):
        return [b for b in self.bodies if b.alive and any(b.J)]

    # ---- fields
    def ext_field(self, points):
        """Field of the profile source: the exact exterior solution of the scalar-stage equation (framework H1)."""
        if self.profile is None:
            return 0.0, 0.0
        c = self.profile
        d = points - np.array([c.x, c.y])
        r = np.sqrt(np.sum(d * d, axis=1))
        L, dL = bh_L(r, c.mass)
        with np.errstate(over="ignore", invalid="ignore"):
            g = np.nan_to_num((dL / np.maximum(r, 1e-300))[:, None] * d)  # zero at the centre (the source is fixed)
        return L, g

    def rod_factor(self, rod):
        """Active line density of a rod in the current framework (F7.4a)."""
        return rod.lam * (self.law.light_src if rod.light else self.law.src(rod.vb))

    def rod_field(self, points):
        """L = -U and grad L of all rods (uniform segments along x) at points (P, 2)."""
        L = np.zeros(len(points))
        g = np.zeros((len(points), 2))
        for rod in self.rods():
            lam = self.rod_factor(rod)
            if lam == 0.0:
                continue
            dx = points[:, 0] - rod.x
            dy = points[:, 1] - rod.y
            s2 = dy * dy + EPS * EPS
            s = np.sqrt(s2)
            ap, am = rod.half - dx, rod.half + dx
            rp, rm = np.sqrt(ap * ap + s2), np.sqrt(am * am + s2)
            U = lam * (np.arcsinh(ap / s) + np.arcsinh(am / s))
            dUdx = lam * (-1.0 / rp + 1.0 / rm)
            dUdy = -lam * dy / s2 * (ap / rp + am / rm)
            L -= U
            g[:, 0] -= dUdx
            g[:, 1] -= dUdy
        return L, g

    @staticmethod
    def field_at(points, src_pos, src_mass, exclude=None):
        """
        L and grad L at points (P, 2) from point sources. exclude[p] is the index of a source left out for point p
        (its own field), or -1; exclude="self" leaves out the diagonal when the points are the sources.
        """
        d = points[:, None, :] - src_pos[None, :, :]  # (P, S, 2)
        s2 = np.sum(d * d, axis=2) + EPS * EPS
        s = np.sqrt(s2)
        w = np.broadcast_to(src_mass[None, :], s.shape).copy()
        if isinstance(exclude, str) and exclude == "self":
            np.fill_diagonal(w, 0.0)
        elif exclude is not None:
            ex = np.asarray([-1 if j is None else j for j in exclude], int)
            p = np.flatnonzero(ex >= 0)
            w[p, ex[p]] = 0.0
        L = -np.sum(w / s, axis=1)
        g = np.sum((w / (s2 * s))[:, :, None] * d, axis=1)  # grad of -M/s is +M d / s^3
        return L, g

    def total_field(self, points, pos, mass, exclude=None):
        """Point masses, the profile source and the rods."""
        if len(mass):
            L, g = self.field_at(points, pos, mass, exclude=exclude)
        else:
            L, g = np.zeros(len(points)), np.zeros((len(points), 2))
        if self.profile is not None:
            Le, ge = self.ext_field(points)
            L, g = L + Le, g + ge
        if self.rods():
            Lr, gr = self.rod_field(points)
            L, g = L + Lr, g + gr
        return L, g

    def L_points(self, pts):
        """L at points from all sources."""
        live, pos, vel, mass, photon, fixed = self.arrays()
        return self.total_field(pts, pos, mass)[0]

    def L_at(self, x, y, skip=None):
        """L at one point from all sources; skip is a body left out."""
        src = [b for b in self.bodies if b.alive and b.mass > 0 and b is not skip and b is not self.profile
               and b.kind != "rod"]
        pt = np.array([[x, y]])
        L = 0.0
        if src:
            L += float(self.field_at(pt, np.array([[b.x, b.y] for b in src]), np.array([b.mass for b in src]))[0][0])
        if self.profile is not None:
            L += float(self.ext_field(pt)[0][0])
        if self.rods():
            L += float(self.rod_field(pt)[0][0])
        return L

    def gm_bz(self, points, skip=None):
        """
        z-component of curl g_0i from spinning bodies at in-plane points, times the frame-dragging factor of the
        framework. g_0i = -2 (J x d)_i / r^3; in the plane only J_z enters: B_z = 2 J_z / r^3. A body moving in the
        plane feels a = v x B = (v_y B_z, -v_x B_z).
        """
        bz = np.zeros(len(points))
        if not self.law.gm:
            return bz
        for k, sb in enumerate(self.spins()):
            if sb.J[2] == 0.0:
                continue
            d = points - np.array([sb.x, sb.y])
            r2 = np.sum(d * d, axis=1) + EPS * EPS
            contrib = 2.0 * sb.J[2] / (r2 * np.sqrt(r2))
            if skip is not None:
                contrib = np.where(np.asarray(skip) == id(sb), 0.0, contrib)
            bz += contrib
        return self.law.gm * bz

    def extra_acc(self, live):
        """External acceleration of a steered body: SI m/s^2 converted to the internal unit 1/length (divided by c^2)."""
        return np.array([([b.thrust[0] / C_SI**2, b.thrust[1] / C_SI**2] if b.kind == "pilot" else [0.0, 0.0])
                         for b in live], float).reshape(-1, 2)

    @staticmethod
    def geodesic_accel(vel, L, gL):
        """Geodesic acceleration of H1 in closed form: dv/dt = -e^(4L) grad L - |v|^2 grad L + 4 v (v.grad L)."""
        v2 = np.sum(vel * vel, axis=1)[:, None]
        vg = np.sum(vel * gL, axis=1)[:, None]
        return -np.exp(4.0 * L)[:, None] * gL - v2 * gL + 4.0 * vel * vg

    def law_accel(self, vel, L, gL):
        """Acceleration of the current framework: H1 in closed form, the others from the static metric."""
        law = self.law
        if law.key == "H1":
            return self.geodesic_accel(vel, L, gL)
        U, gU = -L, -gL
        if law.kind == "newton":
            return gU.copy()
        A, dA, B, dB = law.A(U), law.dA(U), law.B(U), law.dB(U)
        A = np.maximum(A, 1e-300)
        v2 = np.sum(vel * vel, axis=1)
        vg = np.sum(vel * gU, axis=1)
        return (-dA / (2 * B) + v2 * dB / (2 * B))[:, None] * gU + ((dA / A - dB / B) * vg)[:, None] * vel

    def radiation_pairs(self, live, mass, fixed):
        """Index pairs of moving massive bodies that radiate (only when the scene switches radiation on)."""
        if not self.radiation or self.law.rad == 0.0:
            return []
        idx = [i for i, b in enumerate(live) if b.kind == "mass" and mass[i] > 0 and not fixed[i] and not b.probe]
        return [(i, j) for k, i in enumerate(idx) for j in idx[k + 1:]]

    def radiation_rel(self, x, v, m1, m2):
        """2.5PN radiation-reaction acceleration of the relative orbit x = x1 - x2 (harmonic gauge), times the
        radiation factor of the framework."""
        M = m1 + m2
        eta = m1 * m2 / (M * M)
        r = math.hypot(x[0], x[1])
        n = x / r
        rd = float(n @ v)
        v2 = float(v @ v)
        return self.law.rad * 1.6 * eta * M * M / r**3 * (rd * (3 * v2 + 17 * M / (3 * r)) * n - (v2 + 3 * M / r) * v)

    def deriv(self, pos, vel, mass, fixed, pairs=()):
        """Time derivatives of the joint state; source positions follow every stage of RK4."""
        L, gL = self.total_field(pos, pos, mass, exclude="self")
        acc = self.law_accel(vel, L, gL)
        if self.law.gm and self._spins:
            bz = self.gm_bz(pos, skip=self._spin_ids)
            acc = acc + np.stack([vel[:, 1] * bz, -vel[:, 0] * bz], axis=1)
        for i, j in pairs:
            a_rel = self.radiation_rel(pos[i] - pos[j], vel[i] - vel[j], mass[i], mass[j])
            M = mass[i] + mass[j]
            acc[i] += mass[j] / M * a_rel
            acc[j] -= mass[i] / M * a_rel
        if self._extra is not None:
            acc = acc + self._extra
        for i, vt in self._vt:
            acc[i] = (vt - vel[i]) / TAU_NULL
        acc[fixed] = 0.0
        dpos = vel.copy()
        dpos[fixed] = 0.0
        return dpos, acc

    def choose_dt(self, pos, vel, mass, photon):
        """
        Adaptive step: a small fraction of the shortest orbital/crossing time of the bodies that MOVE (static clocks,
        landed and fixed bodies do not limit it). With nothing moving the field is static and one long step is exact:
        returns inf (the caller takes the rest of the frame), or dt_max if a scripted craft is present. Near a solid
        surface the step is also limited by the time to impact, so a fall ends in ~one RK4 step on the surface.
        """
        live = [b for b in self.bodies if b.alive and b.kind != "rod"]
        moving = np.array([not (b.fixed or b.kind == "craft") for b in live], bool)
        if not np.any(moving):
            return self.dt_max if any(b.kind == "craft" for b in live) else math.inf
        speed = np.sqrt(np.sum(vel * vel, axis=1)) + 1e-12
        if self.profile is not None:  # time scales near the profile source
            c = self.profile
            d = np.sqrt((pos[:, 0] - c.x) ** 2 + (pos[:, 1] - c.y) ** 2) + 1e-12
            tb = np.minimum(d / speed, np.where(photon, np.inf, np.sqrt(d**3 / c.mass)))
            dt = float(np.clip(self.eta * np.min(tb[moving]), 1e-6, self.dt_max))
        else:
            # sources whose field is negligible (probes, apples, satellites: < 1e-9 of the heaviest) do not limit it
            src = mass > 1e-9 * (float(mass.max()) if len(mass) else 0.0)
            if not np.any(src) or len(mass) < 2:
                dt = self.dt_max
            else:
                d = pos[:, None, :] - pos[None, src, :]
                sdist = np.sqrt(np.sum(d * d, axis=2) + EPS * EPS)
                idx = np.where(src)[0]
                sdist[idx, np.arange(len(idx))] = np.inf  # a source does not limit its own step
                t_cross = np.min(sdist, axis=1) / speed
                t_orb = np.min(np.sqrt(sdist**3 / (mass[src][None, :] + 1e-30)), axis=1)
                t_body = np.where(photon, t_cross, np.minimum(t_cross, t_orb))
                dt = float(np.clip(self.eta * np.min(t_body[moving]), 1e-6, self.dt_max))
        for rod in self.rods():  # time scales near the beams: distance / sqrt(2 lambda)
            lam = max(abs(self.rod_factor(rod)), 1e-30)
            dx = np.maximum(np.abs(pos[:, 0] - rod.x) - rod.half, 0.0)
            s = np.hypot(dx, pos[:, 1] - rod.y) + 1e-12
            t_rod = np.minimum(s / np.sqrt(2.0 * lam), s / speed)
            dt = min(dt, float(np.clip(self.eta * np.min(t_rod[moving]), 1e-6, self.dt_max)))
        # time to impact on a solid surface; only close to the surface
        for o in (b for b in live if b.radius > 0):
            for i, b in enumerate(live):
                if not moving[i] or photon[i] or b.radius > 0:
                    continue
                rx, ry = pos[i, 0] - o.x, pos[i, 1] - o.y
                r = math.hypot(rx, ry)
                h = r - o.radius
                if r <= 0.0 or h <= 0.0 or h > 0.05 * o.radius:
                    continue
                vr = (vel[i, 0] * rx + vel[i, 1] * ry) / r
                vt2 = max(float(vel[i] @ vel[i]) - vr * vr, 0.0)
                g = max(o.mass / (r * r) - vt2 / r, 0.0)  # gravity minus the centripetal part of an orbit
                u = max(-vr, 0.0)  # inward radial speed
                if g > 0:
                    t_imp = (-u + math.sqrt(u * u + 2.0 * g * h)) / g
                else:
                    t_imp = h / u if u > 0 else math.inf
                dt = min(dt, t_imp * (1.0 + 1e-7) + 1e-9)  # the step ends just past the surface
        return dt

    def step(self, dt=None):
        """One RK4 step of the living bodies; photons are then put back on the null cone |v| = sqrt(A/B)."""
        live, pos, vel, mass, photon, fixed = self.arrays()
        if not live:
            return 0.0
        self._extra = self.extra_acc(live) if any(b.kind == "pilot" for b in live) else None
        # "null field" hypothesis: the velocity of a steered body relaxes to the target within TAU_NULL
        self._vt = [(i, np.array(b.vtarget, float)) for i, b in enumerate(live)
                    if b.kind == "pilot" and getattr(b, "vtarget", None) is not None]
        self._spins = [b for b in self.spins()] if self.law.gm else []
        self._spin_ids = [id(b) for b in live]  # a spinning body does not drag itself
        pairs = self.radiation_pairs(live, mass, fixed)
        if dt is None:
            dt = self.choose_dt(pos, vel, mass, photon)
            if not math.isfinite(dt):
                dt = self.dt_max
        k1p, k1v = self.deriv(pos, vel, mass, fixed, pairs)
        k2p, k2v = self.deriv(pos + 0.5 * dt * k1p, vel + 0.5 * dt * k1v, mass, fixed, pairs)
        k3p, k3v = self.deriv(pos + 0.5 * dt * k2p, vel + 0.5 * dt * k2v, mass, fixed, pairs)
        k4p, k4v = self.deriv(pos + dt * k3p, vel + dt * k3v, mass, fixed, pairs)
        pos_n = pos + dt / 6.0 * (k1p + 2 * k2p + 2 * k3p + k4p)
        vel_n = vel + dt / 6.0 * (k1v + 2 * k2v + 2 * k3v + k4v)
        L, _ = self.total_field(pos_n, pos_n, mass, exclude="self")
        law = self.law
        for i, b in enumerate(live):
            if b.fixed:
                v2 = 0.0
            elif b.kind == "craft":
                v2 = b.vx * b.vx + b.vy * b.vy
            else:
                b.x, b.y = pos_n[i]
                if b.kind == "photon" and law.kind == "metric":  # keep the null condition |v| = sqrt(A/B)
                    sp = math.hypot(*vel_n[i])
                    U = -float(L[i])
                    vel_n[i] *= math.sqrt(max(float(law.A(U)), 0.0) / float(law.B(U))) / max(sp, 1e-300)
                b.vx, b.vy = vel_n[i]
                v2 = float(vel_n[i] @ vel_n[i])
            if b.kind == "photon":
                continue
            if law.kind == "newton":
                b.tau += dt  # Newton: time flows alike everywhere
            elif law.key != "H1":
                U = -max(float(L[i]), -20.0)
                delta = float(law.Am1(U)) - float(law.B(U)) * v2  # A - B v^2 - 1
                d = -math.expm1(0.5 * math.log1p(max(delta, -1.0 + 1e-16)))
                b.tau_def += d * dt
                b.tau += (1.0 - d) * dt
            else:
                # dtau/dt = sqrt(e^(2L) - e^(-2L)|v|^2); the deficit 1 - dtau/dt in a stable logarithmic form keeps
                # clock differences of 1e-15 in the Earth field (|L| ~ 1e-10)
                Li = max(float(L[i]), -20.0)
                q = min(math.exp(-4.0 * Li) * v2, 1.0 - 1e-16)
                d = -math.expm1(Li + 0.5 * math.log1p(-q))
                b.tau_def += d * dt
                b.tau += (1.0 - d) * dt
        for i, j in pairs:  # radiated energy: minus the work of the reaction force on the relative orbit
            a_rel = self.radiation_rel(pos_n[i] - pos_n[j], vel_n[i] - vel_n[j], mass[i], mass[j])
            mu = mass[i] * mass[j] / (mass[i] + mass[j])
            key = (live[i].name, live[j].name)
            self.rad_E[key] = self.rad_E.get(key, 0.0) - mu * float(a_rel @ (vel_n[i] - vel_n[j])) * dt
        self.t += dt
        self.update_gyros(live, pos_n, vel_n, mass, dt)
        for i, b in enumerate(live):
            if b.kind == "craft" and b.script is not None:
                self.move_craft(b, pos_n, mass, i)
        self.apply_capture()
        self.land_on_solids()
        if not self.no_guard:
            self.strong_field_guard()
        return dt

    def update_gyros(self, live, pos, vel, mass, dt):
        """
        Gyroscope precession accumulated over the step. Geodetic: Omega = (gamma + 1/2) v x grad U (out of the plane,
        rotates the spin within the plane). Frame dragging: Omega = gm [3 (J.n) n - J] / r^3 (GR: gm = 1); its
        component along s0 x z tilts the spin out of the plane. Newton: no precession.
        """
        law = self.law
        for i, b in enumerate(live):
            g = b.gyro
            if g is None or not b.alive:
                continue
            g["t"] += dt
            self._gyro_orbit(b, g)
            if law.kind == "newton":
                continue
            gU = -self.field_at(pos[i:i + 1], pos, mass, exclude=[i])[1][0]
            vx, vy = vel[i]
            g["geo"] += (law.gamma + 0.5) * (vx * gU[1] - vy * gU[0]) * dt
            if not law.gm:
                continue
            sx, sy = g["s0"]
            axis = np.array([sy, -sx, 0.0])  # s0 x z
            for sb in self.spins():
                if sb is b:
                    continue
                d = np.array([b.x - sb.x, b.y - sb.y, 0.0])
                r = float(np.linalg.norm(d))
                if r <= 0.0:
                    continue
                n = d / r
                J = np.array(sb.J, float)
                om = law.gm * (3.0 * float(J @ n) * n - J) / r**3
                g["fd"] += float(om @ axis) * dt  # tilt out of the plane
                g["fd_z"] += float(om[2]) * dt  # rotation within the plane

    def _gyro_orbit(self, b, g):
        """Counts whole orbits of a gyroscope around the heaviest other body and stores the accumulated precession
        after each of them (orbit averages are exact only over whole orbits)."""
        c = max((o for o in self.bodies if o.alive and o.mass > 0 and o is not b and o.kind == "mass"),
                key=lambda o: o.mass, default=None)
        if c is None:
            return
        ang = math.atan2(b.y - c.y, b.x - c.x)
        if g["ang"] is not None:
            g["turn"] += math.atan2(math.sin(ang - g["ang"]), math.cos(ang - g["ang"]))
            if abs(g["turn"]) >= 2 * math.pi * (len(g["snaps"]) + 1):
                g["snaps"].append((g["t"], g["geo"], g["fd"], g["fd_z"]))
        g["ang"] = ang

    def apply_capture(self):
        """Black-hole scene: photons and bodies inside the capture radius stop (infinite redshift region)."""
        if self.capture is None:
            return
        c, rc = self.capture
        if not c.alive:
            self.capture = None
            return
        for b in self.bodies:
            if b.alive and not b.fixed and b is not c and b.kind in ("photon", "mass"):
                if math.hypot(b.x - c.x, b.y - c.y) < rc:
                    b.fixed, b.captured = True, True
                    b.vx = b.vy = 0.0

    def move_craft(self, b, pos, mass, i):
        """Scripted motion in SI; the proper acceleration a_kin - a_grav is what a thrust would have to supply."""
        x, y, vx, vy, ax, ay = b.script.state(self.t / C_SI)
        b.x, b.y, b.vx, b.vy = x, y, vx / C_SI, vy / C_SI
        _, g = self.field_at(np.array([[x, y]]), pos, mass, exclude=[i])
        agx, agy = -g[0, 0] * C_SI**2, -g[0, 1] * C_SI**2  # Newtonian gravitational acceleration [m/s^2]
        b.now_g = math.hypot(ax - agx, ay - agy)
        b.max_g = max(b.max_g, b.now_g)

    def land_on_solids(self):
        """A body inside a solid stops on its surface; a photon is absorbed. Rods are thin solids for the probes."""
        solids = [o for o in self.bodies if o.alive and o.radius > 0]
        rods = self.rods()
        for b in self.bodies:
            if not b.alive or b.fixed or b.radius > 0 or b.kind in ("craft", "rod"):
                continue
            for rod in rods:
                if abs(b.x - rod.x) <= rod.half and abs(b.y - rod.y) <= 0.02 * rod.half:
                    if b.kind == "photon":
                        b.alive = False
                    else:
                        b.vx = b.vy = 0.0
                        b.fixed = True
                        b.landed_t = self.t
            if b.fixed or not b.alive:
                continue
            for o in solids:
                dx, dy = b.x - o.x, b.y - o.y
                r = math.hypot(dx, dy)
                if r < o.radius and b.kind == "pilot":  # a steered body stops on the surface and may take off
                    ux, uy = dx / max(r, 1e-12), dy / max(r, 1e-12)
                    b.x, b.y = o.x + ux * o.radius, o.y + uy * o.radius
                    vr = b.vx * ux + b.vy * uy
                    if vr < 0:
                        b.landing_speed = -vr * C_SI
                        b.vx, b.vy = b.vx - vr * ux, b.vy - vr * uy
                        b.vx *= 0.98  # some ground friction
                        b.vy *= 0.98
                    continue
                if r < o.radius:
                    if r < 1e-12 * max(o.radius, 1.0):  # exactly at the centre the radial direction is undefined
                        dx, dy, r = 0.0, 1.0, 1.0
                    if b.kind == "photon":
                        b.alive = False
                    else:
                        b.x, b.y = o.x + dx / r * o.radius, o.y + dy / r * o.radius
                        b.vx = b.vy = 0.0
                        b.fixed = True
                        b.landed_t = self.t

    def strong_field_guard(self):
        """
        If the other sources make |U| > L_MAX at a body, the superposition of point fields is not trusted: a photon
        is absorbed and a massive body merges with the source that dominates its field (Newtonian mass and momentum
        conservation). Each event is written to self.events.
        """
        src = [b for b in self.bodies if b.alive and b.mass > 0 and b.kind == "mass"]
        cand = [b for b in self.bodies if b.alive and not b.fixed and b.kind in ("mass", "photon")]
        if not src or not cand:
            return
        sp = np.array([[o.x, o.y] for o in src])
        sm = np.array([o.mass for o in src])
        pts = np.array([[b.x, b.y] for b in cand])
        ex = [next((k for k, o in enumerate(src) if o is b), -1) for b in cand]
        L, _ = self.field_at(pts, sp, sm, exclude=ex)
        for k in np.flatnonzero(L < -L_MAX):
            b = cand[k]
            if not b.alive:  # merged earlier in this pass; its mass is already counted
                continue
            others = [o for o in src if o is not b and o.alive]
            if not others:
                continue
            dom = max(others, key=lambda o: o.mass / math.sqrt((b.x - o.x) ** 2 + (b.y - o.y) ** 2 + EPS * EPS))
            if b.kind == "photon":
                b.alive = False
                self.events.append(f"t = {self.t:.4g}: {b.name} absorbed by {dom.name} (|U| > {L_MAX:g}, "
                                   "outside the weak field)")
            else:
                big, small = (dom, b) if dom.mass >= b.mass else (b, dom)
                if small.fixed and not big.fixed:
                    big, small = small, big
                M = big.mass + small.mass
                if not big.fixed:
                    big.x = (big.mass * big.x + small.mass * small.x) / M
                    big.y = (big.mass * big.y + small.mass * small.y) / M
                    big.vx = (big.mass * big.vx + small.mass * small.vx) / M
                    big.vy = (big.mass * big.vy + small.mass * small.vy) / M
                big.mass = M
                small.alive = False
                self.events.append(f"t = {self.t:.4g}: {small.name} merged with {big.name} (|U| > {L_MAX:g}, "
                                   "outside the weak field)")

    def advance(self, T, max_steps=5000, on_step=None):
        """Advances the simulation by T or up to the step limit; calls on_step(sim) after each step."""
        t_end = self.t + T
        steps = 0
        while self.t < t_end and steps < max_steps:
            live, pos, vel, mass, photon, _ = self.arrays()
            dt = min(self.choose_dt(pos, vel, mass, photon), t_end - self.t)
            self.step(max(dt, 1e-9))
            steps += 1
            if on_step:
                on_step(self)
        return steps

    def pull_split(self, body, main):
        """|grad U| at a body from a chosen main source and from all other sources (vectorised)."""
        src = [o for o in self.bodies if o.alive and o.mass > 0 and o is not body and o.kind != "rod"]
        if not src:
            return 0.0, 0.0
        d = np.array([[body.x - o.x, body.y - o.y] for o in src])
        m = np.array([o.mass for o in src])
        s2 = np.sum(d * d, axis=1) + EPS * EPS
        g = (m / (s2 * np.sqrt(s2)))[:, None] * d
        is_main = np.array([o is main for o in src])
        gm = g[is_main].sum(axis=0) if is_main.any() else np.zeros(2)
        go = g[~is_main].sum(axis=0) if (~is_main).any() else np.zeros(2)
        return float(np.hypot(*gm)), float(np.hypot(*go))


# =============================================================================
# Measurements: perihelion advance, deflection of light and probes, Shapiro delay, binary radiation
# =============================================================================


def heaviest_other(sim, body):
    others = [b for b in sim.bodies if b.alive and b.mass > 0 and b is not body and b.kind == "mass"]
    return max(others, key=lambda b: b.mass) if others else None


class PerihelionTracker:
    """Measures the perihelion direction of a body around a central source and its advance."""

    def __init__(self, body, center):
        self.body, self.center = body, center
        self.r_hist = deque(maxlen=3)
        self.peri = []  # (time, angle, distance)
        self.r_max = 0.0
        self.r_min = math.inf
        self.other_max = 0.0  # largest pull of the other bodies relative to the central source
        self.lrl = []  # (time, direction of the Laplace-Runge-Lenz vector at perihelion)
        self.sense = 1.0
        self._n = 0

    def lrl_angle(self):
        """Direction of the osculating Newtonian eccentricity vector e = (v^2/mu - 1/r) r - (r.v/mu) v. Sampled at
        each perihelion it removes the periodic wobble of the osculating orbit; its secular drift is the advance."""
        b, c = self.body, self.center
        rx, ry = b.x - c.x, b.y - c.y
        vx, vy = b.vx - c.vx, b.vy - c.vy
        mu = c.mass + b.mass
        r = math.hypot(rx, ry)
        v2, rv = vx * vx + vy * vy, rx * vx + ry * vy
        ex = (v2 / mu - 1.0 / r) * rx - rv / mu * vx
        ey = (v2 / mu - 1.0 / r) * ry - rv / mu * vy
        return math.atan2(ey, ex)

    def rate(self):
        """Secular advance [rad per unit time] from a linear fit of the Laplace-Runge-Lenz angles."""
        if len(self.lrl) < 3:
            return None
        t = np.array([x[0] for x in self.lrl])
        a = np.unwrap([x[1] for x in self.lrl])
        return float(np.polyfit(t, a, 1)[0]) * self.sense

    def update(self, sim):
        b, c = self.body, self.center
        dx, dy = b.x - c.x, b.y - c.y
        r = math.hypot(dx, dy)
        self.r_hist.append((r, math.atan2(dy, dx), sim.t))
        self.r_max, self.r_min = max(self.r_max, r), min(self.r_min, r)
        self._n += 1
        if self._n % 8 == 1:  # the disturbance check is slow-varying; every 8th step is enough
            gm, go = sim.pull_split(b, c)
            if gm > 0:
                self.other_max = max(self.other_max, go / gm)
        lz = dx * (b.vy - c.vy) - dy * (b.vx - c.vx)
        self.sense = 1.0 if lz >= 0 else -1.0  # +1 counter-clockwise, -1 clockwise
        if len(self.r_hist) == 3:
            (r0, a0, t0), (r1, a1, t1), (r2, a2, t2) = self.r_hist
            near_min = r1 - self.r_min <= 0.3 * (self.r_max - self.r_min) + 1e-12  # ignore small wobbles
            if r1 < r0 and r1 <= r2 and near_min:  # the minimum refined by a parabola
                den = r0 - 2 * r1 + r2
                u = 0.5 * (r0 - r2) / den if den > 0 else 0.0
                da_m = math.atan2(math.sin(a1 - a0), math.cos(a1 - a0))
                da_p = math.atan2(math.sin(a2 - a1), math.cos(a2 - a1))
                ang = a1 + (u * (da_p if u > 0 else da_m))
                if not self.peri or t1 - self.peri[-1][0] > 1e-9:
                    self.peri.append((t1, ang, r1))
                    self.lrl.append((sim.t, self.lrl_angle()))

    def result(self):
        """(advance per orbit, GR prediction 6 pi M / (a (1 - e^2)), a, e, orbits); None before two full orbits or for
        a nearly circular orbit (e < 0.05)."""
        if len(self.peri) < 3 or self.r_max <= 0:
            return None
        a = 0.5 * (self.r_min + self.r_max)
        e = (self.r_max - self.r_min) / (self.r_max + self.r_min)
        if e < 0.05:
            return None
        steps = [math.atan2(math.sin(p2[1] - p1[1]), math.cos(p2[1] - p1[1]))
                 for p1, p2 in zip(self.peri[:-1], self.peri[1:])]
        adv = float(np.median(steps)) * self.sense
        M = self.center.mass + self.body.mass
        pred = 6.0 * math.pi * M / (a * (1.0 - e * e))
        return adv, pred, a, e, len(steps)


class DeflectionTracker:
    """Deflection of a photon or a fast probe and, for photons, the Shapiro delay, against first-order predictions
    for the straight path through the same end points."""

    def __init__(self, particle, mass_body, sim):
        self.ph, self.m = particle, mass_body
        sp = math.hypot(particle.vx, particle.vy)
        self.v0 = sp
        self.d0 = (particle.vx / sp, particle.vy / sp)
        rx, ry = particle.x - mass_body.x, particle.y - mass_body.y
        self.x1 = rx * self.d0[0] + ry * self.d0[1]  # position along the initial line
        self.b = abs(rx * self.d0[1] - ry * self.d0[0])  # impact parameter
        self.t0 = sim.t
        self.path = 0.0  # Euclidean length of the actual path
        self.last = (particle.x, particle.y)
        self.int_main, self.int_other = 0.0, 0.0  # pull of the main source and of the others along the path
        self.turn = 0.0  # accumulated signed turn of the direction of motion
        self.last_dir = math.atan2(particle.vy, particle.vx)
        self.photon = particle.kind == "photon"
        self._n = 0

    def update(self, sim):
        x, y = self.ph.x, self.ph.y
        dl = math.hypot(x - self.last[0], y - self.last[1])
        self.path += dl
        self.last = (x, y)
        self._n += 1
        if self._n % 4 == 1:
            gm, go = sim.pull_split(self.ph, self.m)
            self.int_main += gm * dl * 4
            self.int_other += go * dl * 4
        ang = math.atan2(self.ph.vy, self.ph.vx)
        # the angle is summed step by step, so a photon that circles a black hole can turn by more than 180 degrees
        self.turn += math.atan2(math.sin(ang - self.last_dir), math.cos(ang - self.last_dir))
        self.last_dir = ang

    @staticmethod
    def _x_plus_r(x, b):
        """x + sqrt(x^2 + b^2) without cancellation for x < 0 (then b^2 / (r - x))."""
        r = math.hypot(x, b)
        return x + r if x >= 0 else b * b / (r - x)

    def state(self, sim):
        ph, m, d0 = self.ph, self.m, self.d0
        law = sim.law
        sp = math.hypot(ph.vx, ph.vy)
        cosang = max(-1.0, min(1.0, (ph.vx * d0[0] + ph.vy * d0[1]) / max(sp, 1e-300)))
        defl = abs(self.turn) if self.turn != 0.0 else math.acos(cosang)
        rx, ry = ph.x - m.x, ph.y - m.y
        x2 = rx * d0[0] + ry * d0[1]
        b, x1, M = self.b, self.x1, m.mass
        other = self.int_other / self.int_main if self.int_main > 0 else 0.0
        out = {"defl": defl, "b": b, "other": other, "photon": self.photon, "v0": self.v0,
               "delay": (sim.t - self.t0) - self.path if self.photon else None}
        # aimed at the source, or so close that the guard absorbs it: no weak-field prediction
        if b <= (EPS if sim.no_guard else max(EPS, M / L_MAX)):
            out.update(pred=None, pred1=None, newton=None, delay_pred=None, direct=True)
            return out
        geom = (x2 / math.hypot(x2, b) - x1 / math.hypot(x1, b)) * M / b
        g = law.gamma
        if self.photon:
            pred1 = (1.0 + g) * geom
            out.update(pred1=pred1, pred=pred1 + law.deflection2() * (M / b) ** 2, newton=geom,
                       delay_pred=(1.0 + g) * M * math.log(self._x_plus_r(x2, b) / self._x_plus_r(x1, b)),
                       direct=False)
        else:
            v2 = self.v0 * self.v0
            newton = geom / v2
            out.update(newton=newton, pred1=newton * (1.0 + g * v2), pred=newton * (1.0 + g * v2),
                       delay_pred=None, direct=False)
        return out


def make_trackers(sim, orbits, old_p=None, old_l=None):
    """Perihelion and deflection trackers. Existing trackers are kept by body identity, so adding an object does not
    erase the measurements so far."""
    old_p, old_l = old_p or {}, old_l or {}
    trackers, light = [], []
    for b in sim.bodies:
        if not b.alive or b.kind == "rod":
            continue
        if b.kind == "photon" or (b.kind == "mass" and b.probe):
            if id(b) in old_l:
                light.append(old_l[id(b)])
            else:
                m = heaviest_other(sim, b)
                if m is not None:
                    light.append(DeflectionTracker(b, m, sim))
            continue
        if b.kind == "mass" and not b.fixed and orbits:
            c = heaviest_other(sim, b)
            if id(b) in old_p and old_p[id(b)].center is c:
                trackers.append(old_p[id(b)])
            elif c is not None and (c.mass > 10 * b.mass or (orbits == "binary" and c.mass >= b.mass)):
                trackers.append(PerihelionTracker(b, c))
    return trackers, light


# ---------------------------------------------------------------- binary radiation

def chirp_mass(m1, m2):
    return (m1 * m2) ** 0.6 / (m1 + m2) ** 0.2


def peters_f(e):
    """Eccentricity enhancement of the quadrupole radiation, f(e) = (1 + 73/24 e^2 + 37/96 e^4) / (1 - e^2)^(7/2)."""
    return (1 + 73 / 24 * e * e + 37 / 96 * e**4) / (1 - e * e) ** 3.5


def peters_pdot(m1, m2, P, e):
    """Quadrupole (GR) period derivative dP/dt = -(192 pi/5) (2 pi Mc / P)^(5/3) f(e) (G = c = 1, dimensionless)."""
    return -192.0 * math.pi / 5.0 * (2.0 * math.pi * chirp_mass(m1, m2) / P) ** (5.0 / 3.0) * peters_f(e)


def peters_lum(m1, m2, a, e):
    """Orbit-averaged quadrupole luminosity (32/5) eta^2 M^5 / a^5 f(e) (G = c = 1)."""
    M = m1 + m2
    eta = m1 * m2 / (M * M)
    return 32.0 / 5.0 * eta**2 * M**5 / a**5 * peters_f(e)


def peters_merge_time(m1, m2, a):
    """Time to coalescence of a circular binary, 5 a^4 / (256 mu M^2) (G = c = 1)."""
    M = m1 + m2
    return 5.0 * a**4 / (256.0 * (m1 * m2 / M) * M * M)


class BinaryTracker:
    """
    Radiation of a binary: periastron passages of the relative orbit and the energy radiated between them (the work
    of the reaction force, Sim.rad_E). dP/dt = (3/2) P L / E with E = -mu M / (2a) (Keplerian), L the luminosity
    averaged over whole orbits.
    """

    def __init__(self, b1, b2, sim, start_at_periastron=False):
        self.b1, self.b2 = b1, b2
        self.r_hist = deque(maxlen=3)
        self.peri = []  # (time, radiated energy)
        self.r_min, self.r_max = math.inf, 0.0
        if start_at_periastron:
            self.peri.append((sim.t, self.energy(sim)))

    def energy(self, sim):
        return sim.rad_E.get((self.b1.name, self.b2.name), 0.0) + sim.rad_E.get((self.b2.name, self.b1.name), 0.0)

    def sep(self):
        return math.hypot(self.b1.x - self.b2.x, self.b1.y - self.b2.y)

    def update(self, sim):
        r = self.sep()
        self.r_min, self.r_max = min(self.r_min, r), max(self.r_max, r)
        self.r_hist.append((r, sim.t, self.energy(sim)))
        if len(self.r_hist) == 3:
            (r0, t0, e0), (r1, t1, e1), (r2, t2, e2) = self.r_hist
            if r1 < r0 and r1 <= r2 and (not self.peri or t1 - self.peri[-1][0] > 1e-9):
                self.peri.append((t1, e1))

    def result(self):
        """(period, luminosity, dP/dt, a, e) from whole orbits, or None before one full orbit."""
        if len(self.peri) < 2 or self.r_max <= 0:
            return None
        (ta, ea), (tb, eb) = self.peri[0], self.peri[-1]
        n = len(self.peri) - 1
        P = (tb - ta) / n
        lum = (eb - ea) / (tb - ta)
        m1, m2 = self.b1.mass, self.b2.mass
        M = m1 + m2
        a = (M * (P / (2.0 * math.pi)) ** 2) ** (1.0 / 3.0)
        e = (self.r_max - self.r_min) / (self.r_max + self.r_min)
        E = -m1 * m2 / M * M / (2.0 * a)
        return P, lum, 1.5 * P * lum / E, a, e


# =============================================================================
# Scenes (presets)
# =============================================================================


def kepler_peri_velocity(M, a, e):
    """Newtonian speed at perihelion of an orbit with semi-major axis a and eccentricity e."""
    return math.sqrt(M * (1 + e) / (a * (1 - e)))


# Mean elements J2000: mass [kg], a [AU], e, longitude of perihelion and mean longitude [deg].
PLANETS = [
    ("Mercury", 3.3011e23, 0.387098, 0.205630, 77.456, 252.251, "#c9a37a"),
    ("Venus", 4.8675e24, 0.723332, 0.006772, 131.533, 181.980, "#ffe0a0"),
    ("Earth+Moon", 6.0458e24, 1.000001, 0.016709, 102.947, 100.464, "#7ac8ff"),
    ("Mars", 6.4171e23, 1.523679, 0.093400, 336.041, 355.453, "#ff7a5a"),
    ("Jupiter", 1.89819e27, 5.20260, 0.048900, 14.753, 34.404, "#ffcf8a"),
    ("Saturn", 5.6834e26, 9.55490, 0.056500, 92.432, 49.944, "#ffe7b0"),
    ("Uranus", 8.6813e25, 19.2184, 0.046381, 170.964, 313.232, "#a0f0ff"),
    ("Neptune", 1.02413e26, 30.1104, 0.009456, 44.971, 304.880, "#7a9aff"),
]

# PSR B1913+16 (Weisberg & Huang 2016): masses [M_sun], orbital period [d], eccentricity.
HULSE_TAYLOR = {"m_p": 1.4398, "m_c": 1.3886, "Pb_d": 0.322997448918, "e": 0.6171334,
                "pdot_obs": -2.398e-12, "pdot_obs_err": 0.004e-12, "omdot_obs": 4.226598}


def kepler_state(mu, a, e, varpi, lam):
    """Position and velocity on a Kepler ellipse (G = c = 1, velocity in units of c). mu [m], a [m], e, longitude of
    perihelion varpi and mean longitude lam [rad]."""
    M = (lam - varpi) % (2 * math.pi)
    E = M if e < 0.8 else math.pi
    for _ in range(50):  # Newton iterations of Kepler's equation
        E -= (E - e * math.sin(E) - M) / (1 - e * math.cos(E))
    n = math.sqrt(mu / a**3)
    x, y = a * (math.cos(E) - e), a * math.sqrt(1 - e * e) * math.sin(E)
    k = a * n / (1 - e * math.cos(E))
    vx, vy = -k * math.sin(E), k * math.sqrt(1 - e * e) * math.cos(E)
    c, s_ = math.cos(varpi), math.sin(varpi)
    return x * c - y * s_, x * s_ + y * c, vx * c - vy * s_, vx * s_ + vy * c


def solar_system(planets=True, only=None):
    """The Sun and planets at J2000 in 2D; the Sun moves so that the total momentum is zero."""
    sun = Body("Sun", 0, 0, mass=M_SUN, color="#ffd27a")
    B = [sun]
    if planets or only:
        for nm, kg, a, e, vp, lam, col in PLANETS:
            if only is not None and nm not in only:
                continue
            x, y, vx, vy = kepler_state(M_SUN, a * AU, e, math.radians(vp), math.radians(lam))
            B.append(Body(nm, x, y, vx, vy, mass=kg * G_C2, color=col))
    px = sum(b.mass * b.vx for b in B[1:])
    py = sum(b.mass * b.vy for b in B[1:])
    sun.vx, sun.vy = -px / sun.mass, -py / sun.mass
    return B


def binary_bodies(n1, m1, c1, n2, m2, c2, a, e):
    """Two bodies on a Kepler orbit about their centre of mass, starting at periastron (G = c = 1)."""
    x, y, vx, vy = kepler_state(m1 + m2, a, e, 0.0, 0.0)
    f1, f2 = -m2 / (m1 + m2), m1 / (m1 + m2)
    return (Body(n1, f1 * x, f1 * y, f1 * vx, f1 * vy, mass=m1, color=c1),
            Body(n2, f2 * x, f2 * y, f2 * vx, f2 * vy, mass=m2, color=c2))


def circular_speed(m_src, d_src, r_orbit, law=None):
    """
    Speed of a circular orbit of a test body in the static metric of the framework (default: causal gravity): the
    source of mass m_src is at distance d_src, the orbit has radius r_orbit about the centre of rotation (for a binary
    the centre of mass). From the geodesic equation with v.grad U = 0:
        v^2 (1/r_orbit - B' g / 2B) = -A' g / 2B,   U = m_src / d_src,  g = m_src / d_src^2.
    """
    law = law or MAIN_LAW()
    U = m_src / d_src
    g = m_src / d_src**2
    if law.kind == "newton":
        return math.sqrt(g * r_orbit)
    dA, B, dB = law.dA(U), law.B(U), law.dB(U)
    return math.sqrt((-dA * g / (2 * B)) / (1.0 / r_orbit - dB * g / (2 * B)))


def iso_from_areal(R, M=1.0):
    """Isotropic radius of the Schwarzschild metric for the areal radius R: R = r (1 + M/2r)^2."""
    return 0.5 * ((R - M) + math.sqrt(max(R * R - 2 * M * R, 0.0)))


class CraftScript:
    """
    Scripted motion in SI: an initial state and phases given by duration [s] and constant kinematic acceleration
    ax, ay [m/s^2]. After the last phase the body keeps its last velocity; the prepared scripts end at rest.
    """

    def __init__(self, x0, y0, vx0, vy0, phases):
        self.phases = []
        t, x, y, vx, vy = 0.0, float(x0), float(y0), float(vx0), float(vy0)
        for dur, ax, ay in phases:
            self.phases.append((t, t + dur, x, y, vx, vy, ax, ay))
            x, y = x + vx * dur + 0.5 * ax * dur * dur, y + vy * dur + 0.5 * ay * dur * dur
            vx, vy = vx + ax * dur, vy + ay * dur
            t += dur
        self.end = (t, x, y, vx, vy)

    def state(self, t):
        """(x, y, vx, vy, ax, ay) at time t [s]."""
        for t0, t1, x, y, vx, vy, ax, ay in self.phases:
            if t < t1:
                d = max(t - t0, 0.0)
                return x + vx * d + 0.5 * ax * d * d, y + vy * d + 0.5 * ay * d * d, vx + ax * d, vy + ay * d, ax, ay
        t0, x, y, vx, vy = self.end
        d = t - t0
        return x + vx * d, y + vy * d, vx, vy, 0.0, 0.0

    def phase_name(self, t, names):
        for k, (t0, t1, *_rest) in enumerate(self.phases):
            if t < t1:
                return names[k] if k < len(names) else ""
        return names[-1] if names else ""


SCENES = [
    "Perihelion precession",
    "Solar system (exaggerated)",
    "Light bending",
    "Shapiro delay",
    "Binary star and planet",
    "Clocks in a field",
    "Solar system (real, J2000)",
    "Earth surface: apple and clocks",
    "Earth: ISS and GPS",
    "UAP: scripted manoeuvre (hypothesis)",
    "UAP: keyboard control",
    "Black hole (strong field)",
    "α Centauri A/B",
    "Sirius A/B (white dwarf)",
    "Solar neighbourhood to 12 ly (3D)",
    "Free scene",
    "Fast probes and light (1 + γv²)",
    "Frame dragging (exaggerated)",
    "Gravity Probe B (gyroscope)",
    "Binary pulsar PSR B1913+16",
    "Inspiral (strong radiation)",
    "What gravitates: counter-streaming beams (F7.4a)",
]

ROD_SPEEDS = (0.0, 0.3, 0.6, 0.9)  # beam speeds of the F7.4a scene (plus one pair of light beams)


def scene(name):
    """
    Bodies of a preset scene, its view (centre x, y, half-width) and meta:
      si         True: lengths in metres, time shown in seconds; False: geometric units G = c = 1
      dt_max     largest integration step;  eta  safety factor of the step
      speed      default log10 pace (simulated time per second of real time)
      display    fixed colour scale of the field ("log", lo, hi) or deviation from a reference ("rel", L0, range)
      click      "orbit": a click makes a circular orbit; "rest": a click releases the object at rest
      mass       default mass of a new body (value, unit label); None: 1e-6 of the reference mass in model units
      flags      orbits, bh, map3d, pilot, radiation, pulsar, inspiral, probes, dragging, gpb, rods
    """
    B = []
    meta = {"si": False, "dt_max": 0.5, "speed": math.log10(30.0), "display": ("log", -4.0, -1.0), "click": "orbit",
            "mass": None}
    if name == "Perihelion precession":
        M = 0.01
        B.append(Body("star", 0, 0, mass=M, color="#ffd27a"))
        a, e = 1.0, 0.5
        # Newtonian perihelion speed; the 1PN effects change the orbit slightly, a and e are measured afterwards
        B.append(Body("planet", a * (1 - e), 0, 0, kepler_peri_velocity(M, a, e), mass=1e-7, color="#7ac8ff"))
        view = (-0.6, 0.0, 3.0)
    elif name == "Solar system (exaggerated)":
        M = 0.005
        B.append(Body("star", 0, 0, mass=M, color="#ffd27a"))
        for nm, a, e, m, col in (("A", 0.5, 0.3, 2e-7, "#c9a37a"), ("B", 0.9, 0.05, 1e-6, "#7affb0"),
                                 ("C", 1.5, 0.1, 5e-6, "#7ac8ff"), ("D", 2.6, 0.2, 2e-5, "#ff9a7a")):
            B.append(Body(nm, a * (1 - e), 0, 0, kepler_peri_velocity(M, a, e), mass=m, color=col))
        view = (0.0, 0.0, 3.2)
    elif name == "Light bending":
        B.append(Body("mass", 0, 0, mass=0.01, fixed=True, color="#ffd27a"))
        for bb in (0.25, 0.4, 0.7, 1.0, 1.6):
            B.append(Body(f"photon b={bb:g}", -6.0, bb, 1.0, 0.0, kind="photon", color="#fff27a"))
        view = (0.0, 0.0, 6.5)
    elif name == "Shapiro delay":
        B.append(Body("mass", 0, 0, mass=0.01, fixed=True, color="#ffd27a"))
        B.append(Body("photon b=0.3", -20.0, 0.3, 1.0, 0.0, kind="photon", color="#fff27a"))
        B.append(Body("photon b=3", -20.0, 3.0, 1.0, 0.0, kind="photon", color="#ffb37a"))
        view = (0.0, 0.0, 21.0)
    elif name == "Binary star and planet":
        M1 = M2 = 0.005
        sep = 0.6
        v = 0.5 * math.sqrt((M1 + M2) / sep)
        B.append(Body("star 1", -sep / 2, 0, 0, -v, mass=M1, color="#ffd27a"))
        B.append(Body("star 2", sep / 2, 0, 0, v, mass=M2, color="#ff9a5a"))
        rp = 2.4
        B.append(Body("planet", rp, 0, 0, math.sqrt((M1 + M2) / rp), mass=1e-6, color="#7ac8ff"))
        view = (0.0, 0.0, 3.2)
    elif name == "Clocks in a field":
        M = 0.02
        B.append(Body("mass", 0, 0, mass=M, fixed=True, color="#ffd27a"))
        for r in (0.3, 1.0, 4.0):
            B.append(Body(f"clock r={r:g}", r, 0, kind="clock", color="#b0b0ff"))
        a = 1.0
        B.append(Body("orbiting clock r=1", 0, a, -math.sqrt(M / a), 0, mass=1e-9, color="#7affb0"))
        view = (0.0, 0.0, 4.5)
    elif name == "Earth surface: apple and clocks":
        R = R_EARTH
        earth = Body("Earth", 0, 0, mass=M_EARTH, fixed=True, color="#3f6b3a", radius=R)
        earth.J = (0.0, 0.0, J_EARTH)
        B.append(earth)
        B.append(Body("apple (10 m)", -12.0, R + 10.0, mass=0.2 * G_C2, color="#ff5a4a"))
        v0, ang = 14.0, math.radians(50)  # thrown ball: 14 m/s at 50 degrees from 1.5 m
        B.append(Body("ball (14 m/s, 50°)", -30.0, R + 1.5, v0 * math.cos(ang) / C_SI, v0 * math.sin(ang) / C_SI,
                      mass=0.4 * G_C2, color="#ffd27a"))
        g = Body("clock on the ground", 20.0, R + 0.01, kind="clock", color="#b0b0ff")
        g.ground = True
        B.append(g)
        B.append(Body("clock on a tower (30 m)", 20.0, R + 30.0, kind="clock", color="#b0b0ff"))
        view = (0.0, R + 15.0, 40.0)
        meta.update(si=True, dt_max=10.0 * C_SI, speed=0.0, display=("rel", -M_EARTH / R, M_EARTH / R**2 * 60.0),
                    click="rest", mass=(0.2, "kg"))
    elif name == "Earth: ISS and GPS":
        R = R_EARTH
        earth = Body("Earth", 0, 0, mass=M_EARTH, fixed=True, color="#3f6b3a", radius=R)
        earth.J = (0.0, 0.0, J_EARTH)
        B.append(earth)
        g = Body("clock on the ground", 0, R + 0.01, kind="clock", color="#b0b0ff")
        g.ground = True
        B.append(g)
        for nm, r, col in (("ISS", 6.771e6, "#ffffff"), ("GPS", 2.6560e7, "#ffd27a")):
            B.append(Body(nm, r, 0, 0, math.sqrt(M_EARTH / r), mass=1e3 * G_C2, color=col))
        view = (0.0, 0.0, 3.0e7)
        meta.update(si=True, dt_max=10.0 * C_SI, speed=math.log10(2000.0), display=("log", -10.2, -9.1),
                    click="orbit", mass=(1e3, "kg"))
    elif name == "Solar system (real, J2000)":
        B = solar_system()
        view = (0.0, 0.0, 1.8 * AU)
        meta.update(si=True, dt_max=2.0 * 86400 * C_SI, speed=math.log10(30 * 86400.0), eta=0.004,
                    display=("log", -9.7, -7.0), click="orbit", orbits=True, mass=(1.0, "Earth masses (M⊕)"))
    elif name == "UAP: scripted manoeuvre (hypothesis)":
        R = R_EARTH
        B.append(Body("Earth", 0, 0, mass=M_EARTH, fixed=True, color="#3f6b3a", radius=R))
        x0, y0 = -3000.0, R + 1000.0
        craft = Body("UAP", x0, y0, kind="craft", color="#9affff")
        # illustrative manoeuvre: hover, fast flight, abrupt stop, hover, vertical climb and stop
        craft.script = CraftScript(x0, y0, 0.0, 0.0, [(2.0, 0, 0), (1.0, 1000.0, 0), (2.0, 0, 0), (0.2, -5000.0, 0),
                                                       (1.5, 0, 0), (1.0, 0, 1500.0), (1.0, 0, -1500.0),
                                                       (30.0, 0, 0)])
        craft.phase_names = ["hovering", "accelerating at 1000 m/s²", "flying at 1000 m/s", "stopping within 0.2 s",
                             "hovering", "climbing at 1500 m/s²", "braking", "hovering"]
        B.append(craft)
        B.append(Body("stone (free fall from 1000 m)", x0, y0, mass=1.0 * G_C2, color="#ff9a5a"))
        view = (0.0, R + 1200.0, 4000.0)
        meta.update(si=True, dt_max=0.01 * C_SI, speed=0.0, display=("rel", -M_EARTH / R, M_EARTH / R**2 * 3000.0),
                    click="rest", mass=(1.0, "kg"))
    elif name == "Black hole (strong field)":
        # G = c = 1 and M = 1, lengths in GM/c^2. Causal gravity = GR here (decision 13): Schwarzschild metric in
        # isotropic coordinates. The framework H1 (scalar stage) uses the exact profile of its own field equation.
        bh = Body("black hole (M = 1)", 0, 0, mass=1.0, fixed=True, color="#202020")
        bh.is_profile = True
        B.append(bh)
        for bb in (3.5, 4.0, 4.5, 5.0, 5.15, 5.21, 5.23, 5.26, 5.5, 6.0, 7.0, 8.0):
            # the conserved impact parameter is the initial height times the refractive index n = sqrt(B/A) of the
            # Schwarzschild metric, n = (1 + U/2)^3 / (1 - U/2); hence y0 = b / n
            y0 = bb
            for _ in range(6):
                U = 1.0 / math.hypot(40.0, y0)
                y0 = bb * (1 - U / 2) / (1 + U / 2) ** 3
            B.append(Body(f"photon b={bb:g}", -40.0, y0, 1.0, 0.0, kind="photon", color="#fff27a"))
        for R_ar in (3.0, 4.0, 6.0, 10.0):  # static clocks at areal radii R
            B.append(Body(f"clock R={R_ar:g}M", -iso_from_areal(R_ar), 0.0, kind="clock", color="#b0b0ff"))
        r0 = 12.0
        B.append(Body("particle (periapsis 12 M)", r0, 0, 0, 0.85 * math.sqrt(1.0 / r0), mass=1e-9, color="#7ac8ff"))
        view = (0.0, 0.0, 14.0)
        meta.update(dt_max=1.0, eta=0.01, display=("log", -2.0, 0.7), bh=True, speed=math.log10(10**(-0.2) / 0.03))
    elif name == "UAP: keyboard control":
        R = R_EARTH
        B.append(Body("Earth", 0, 0, mass=M_EARTH, fixed=True, color="#3f6b3a", radius=R))
        pl = Body("UAP (pilot)", 0.0, R + 500.0, kind="pilot", color="#9affff")
        pl.hover = True
        B.append(pl)
        B.append(Body("stone (free fall from 500 m)", 60.0, R + 500.0, mass=1.0 * G_C2, color="#ff9a5a"))
        view = (0.0, R + 500.0, 1500.0)
        meta.update(si=True, dt_max=0.005 * C_SI, speed=0.0, display=("rel", -M_EARTH / R, M_EARTH / R**2 * 3000.0),
                    click="rest", pilot=True, mass=(1.0, "kg"))
    elif name in ("α Centauri A/B", "Sirius A/B (white dwarf)"):
        if name.startswith("α"):  # masses [M_sun], a [AU], e; period 79.91 yr
            (n1, m1, c1), (n2, m2, c2), a_au, e = ("α Cen A", 1.0788, "#fff4d0"), ("α Cen B", 0.9092, "#ffd29a"), \
                23.33, 0.5208
        else:  # Sirius: period 50.13 yr
            (n1, m1, c1), (n2, m2, c2), a_au, e = ("Sirius A", 2.063, "#cfe0ff"), ("Sirius B", 1.018, "#ffffff"), \
                19.8, 0.5923
        B.extend(binary_bodies(n1, m1 * M_SUN, c1, n2, m2 * M_SUN, c2, a_au * AU, e))
        view = (0.0, 0.0, 1.3 * a_au * AU)
        meta.update(si=True, dt_max=20 * 86400 * C_SI, speed=math.log10(2 * YEAR), eta=0.004,
                    display=("log", -9.5, -6.5), click="orbit", orbits="binary", mass=(1.0, "solar masses (M☉)"))
    elif name == "Solar neighbourhood to 12 ly (3D)":
        # Equatorial positions J2000 (OpenStax Astronomy, appendix I) and approximate masses [M_sun].
        STARS = [
            ("Sun", 0.0, 0.0, 0.0, "G2 V", 1.0),
            ("Proxima Centauri", 14.48, -62.67, 4.2, "M5.5 V", 0.12),
            ("α Centauri A", 14.65, -60.83, 4.4, "G2 V", 1.08),
            ("α Centauri B", 14.65, -60.83, 4.4, "K1 V", 0.91),
            ("Barnard's Star", 17.95, 4.70, 6.0, "M4 V", 0.14),
            ("Wolf 359", 10.93, 7.00, 7.8, "M6 V", 0.11),
            ("Lalande 21185", 11.05, 35.97, 8.3, "M2 V", 0.39),
            ("Sirius A", 6.75, -16.70, 8.6, "A1 V", 2.06),
            ("Sirius B", 6.75, -16.72, 8.6, "DA2 (white dwarf)", 1.02),
            ("Luyten 726-8 A", 1.65, -17.95, 8.7, "M5.5 V", 0.12),
            ("Ross 154", 18.82, -23.83, 9.7, "M3.5 V", 0.17),
            ("Ross 248", 23.68, 44.17, 10.3, "M5.5 V", 0.14),
            ("ε Eridani", 3.53, -9.45, 10.5, "K2 V", 0.82),
            ("Lacaille 9352", 23.08, -35.85, 10.7, "M0.5 V", 0.48),
            ("Ross 128", 11.78, 0.80, 10.9, "M4 V", 0.17),
            ("EZ Aquarii A", 22.63, -15.28, 11.3, "M5 V", 0.11),
            ("Procyon A", 7.65, 5.22, 11.4, "F5 IV", 1.50),
            ("Procyon B", 7.65, 5.22, 11.4, "white dwarf", 0.60),
            ("61 Cygni A", 21.10, 38.73, 11.4, "K5 V", 0.70),
            ("61 Cygni B", 21.10, 38.73, 11.4, "K7 V", 0.63),
            ("Struve 2398 A", 18.70, 59.62, 11.5, "M3 V", 0.33),
            ("Groombridge 34 A", 0.30, 44.02, 11.6, "M1.5 V", 0.38),
            ("ε Indi A", 22.05, -56.77, 11.8, "K5 V", 0.76),
            ("τ Ceti", 1.73, -15.93, 11.9, "G8.5 V", 0.78),
        ]
        # Proper motion mu [arcsec/yr], position angle [deg, from north through east] and radial velocity [km/s]:
        # mu and the position angle from the RECONS table of the 100 nearest systems (GJ 15 A from Kervella et al.
        # 2019); radial velocities from Kervella et al. 2019, A&A 623, A72, table 1; None = unknown (sky motion only).
        # Companions without own data share the motion of the primary. A missing mu means no motion data.
        MOTION = {
            "Proxima Centauri": (3.853, 281.5, None), "α Centauri A": (3.710, 277.5, None),
            "α Centauri B": (3.724, 284.8, None), "Barnard's Star": (10.358, 355.6, -110.40),
            "Wolf 359": (4.696, 234.6, None), "Lalande 21185": (4.802, 186.9, None),
            "Sirius A": (1.339, 204.1, None), "Sirius B": (1.339, 204.1, None),
            "Luyten 726-8 A": (None, None, None), "Ross 154": (0.666, 106.8, -10.41),
            "Ross 248": (1.617, 177.0, None), "ε Eridani": (0.977, 271.1, 16.20),
            "Lacaille 9352": (6.896, 78.9, None), "Ross 128": (1.361, 153.6, -30.99),
            "EZ Aquarii A": (3.254, 46.6, None), "Procyon A": (1.259, 214.7, None),
            "Procyon B": (1.259, 214.7, None), "61 Cygni A": (5.281, 51.9, -65.73),
            "61 Cygni B": (5.172, 52.6, -64.04), "Struve 2398 A": (None, None, None),
            "Groombridge 34 A": (2.919, 81.9, 11.81), "ε Indi A": (4.704, 122.7, -40.10),
            "τ Ceti": (1.922, 296.4, -16.75),
        }
        colors = {"O": "#9bb0ff", "B": "#aabfff", "A": "#cad7ff", "F": "#f8f7ff", "G": "#fff4ea", "K": "#ffd2a1",
                  "M": "#ffb07a", "D": "#ffffff", "w": "#ffffff"}
        for nm, ra_h, dec_d, d_ly, sp, msun in STARS:
            a, dlt, d = math.radians(15 * ra_h), math.radians(dec_d), d_ly * LY
            st = Body(nm, d * math.cos(dlt) * math.cos(a), d * math.cos(dlt) * math.sin(a), mass=msun * M_SUN,
                      fixed=True, color=colors.get(sp[0], "#ffffff"))
            st.z = d * math.sin(dlt)
            st.info = f"{sp}, {msun:g} M☉, {d_ly:g} ly"
            mu, pa, rv = MOTION.get(nm, (None, None, None))
            st.v3 = (0.0, 0.0, 0.0)  # space velocity in units of c relative to the Sun at rest
            st.motion_note = "Sun (origin)" if nm == "Sun" else "no motion data"
            if mu is not None:
                vt = 4.74047 * mu * d_ly / 3.26156 * 1000.0  # tangential speed [m/s]
                rr = (math.cos(dlt) * math.cos(a), math.cos(dlt) * math.sin(a), math.sin(dlt))
                ee = (-math.sin(a), math.cos(a), 0.0)  # east
                nn = (-math.sin(dlt) * math.cos(a), -math.sin(dlt) * math.sin(a), math.cos(dlt))  # north
                p_ = math.radians(pa)
                vr = (rv or 0.0) * 1000.0
                st.v3 = tuple((vr * rr[k] + vt * (math.sin(p_) * ee[k] + math.cos(p_) * nn[k])) / C_SI
                              for k in range(3))
                st.motion_note = (f"v_t {vt / 1000:.0f} km/s, v_r {rv:+.1f} km/s" if rv is not None
                                  else f"v_t {vt / 1000:.0f} km/s, v_r unknown (sky motion only)")
            st.trail3 = deque(maxlen=4000)
            B.append(st)
        view = (0.0, 0.0, 13.5 * LY)
        meta.update(si=True, dt_max=YEAR * C_SI, speed=math.log10(1000 * YEAR), display=("log", -14.0, -12.0),
                    map3d=True, mass=(0.1, "solar masses (M☉)"))
    elif name == "Fast probes and light (1 + γv²)":
        M = 0.005
        B.append(Body("mass", 0, 0, mass=M, fixed=True, color="#ffd27a"))
        cols = ("#7affb0", "#7ad0ff", "#ff9a7a", "#ff7ad9")
        for v, col in zip((0.3, 0.5, 0.7, 0.9), cols):
            p = Body(f"probe v={v:g}c", -30.0, 2.0, v, 0.0, mass=1e-12, color=col)
            p.probe = True
            B.append(p)
        B.append(Body("photon b=2", -30.0, 2.0, 1.0, 0.0, kind="photon", color="#fff27a"))
        view = (0.0, 0.0, 32.0)
        meta.update(dt_max=0.5, eta=0.01, speed=math.log10(6.0), probes=True, display=("log", -4.0, -1.5))
    elif name == "Frame dragging (exaggerated)":
        M = 0.01
        star = Body("spinning star (J = M²)", 0, 0, mass=M, fixed=True, color="#ffd27a")
        star.J = (0.0, 0.0, M * M)  # extremal spin, exaggerated on purpose
        B.append(star)
        a, e = 0.6, 0.5
        vp = kepler_peri_velocity(M, a, e)
        B.append(Body("prograde planet", a * (1 - e), 0, 0, vp, mass=1e-9, color="#7affb0"))
        B.append(Body("retrograde planet", -a * (1 - e), 0, 0, vp, mass=1e-9, color="#ff9a7a"))
        rg = 1.4
        B.append(make_gyro(Body("gyroscope (circular, r = 1.4)", 0, rg, -circular_speed(M, rg, rg), 0, mass=1e-9,
                                color="#9ad9ff"), s0=(1.0, 0.0)))
        view = (0.0, 0.0, 2.0)
        meta.update(dt_max=0.5, eta=0.005, speed=math.log10(60.0), dragging=True)
    elif name == "Gravity Probe B (gyroscope)":
        R = R_EARTH
        earth = Body("Earth (spin axis ↑)", 0, 0, mass=M_EARTH, fixed=True, color="#3f6b3a", radius=R)
        earth.J = (0.0, J_EARTH, 0.0)  # spin axis in the plane: the orbit is polar, as GP-B's
        B.append(earth)
        r = GPB_A
        B.append(make_gyro(Body("Gravity Probe B", r, 0.0, 0.0, math.sqrt(M_EARTH / r), mass=1e-24,
                                color="#ffffff"), s0=(1.0, 0.0)))
        view = (0.0, 0.0, 9.0e6)
        meta.update(si=True, dt_max=30.0 * C_SI, speed=math.log10(600.0), display=("log", -10.2, -9.1),
                    click="orbit", gpb=True, mass=(1e3, "kg"))
    elif name == "Binary pulsar PSR B1913+16":
        ht = HULSE_TAYLOR
        m1, m2 = ht["m_p"] * M_SUN, ht["m_c"] * M_SUN
        P = ht["Pb_d"] * 86400.0 * C_SI
        a = ((m1 + m2) * (P / (2 * math.pi)) ** 2) ** (1.0 / 3.0)
        B.extend(binary_bodies("pulsar", m1, "#9ad9ff", "companion", m2, "#ffd29a", a, ht["e"]))
        view = (0.0, 0.0, 1.4 * a)
        meta.update(si=True, dt_max=60.0 * C_SI, speed=math.log10(3600.0), eta=0.004, display=("log", -6.5, -4.0),
                    click="orbit", orbits="binary", radiation=True, pulsar=True, mass=(1.4, "solar masses (M☉)"))
    elif name == "Inspiral (strong radiation)":
        a = 20.0
        v = circular_speed(0.5, a, 0.5 * a)  # each body: a test body in the field of the other, circular in 1PN
        B.append(Body("black hole A", -0.5 * a, 0, 0, -v, mass=0.5, color="#c0c0c0"))
        B.append(Body("black hole B", 0.5 * a, 0, 0, v, mass=0.5, color="#9ad9ff"))
        view = (0.0, 0.0, 16.0)
        meta.update(dt_max=1.0, eta=0.01, speed=math.log10(300.0), display=("log", -2.5, -0.7), radiation=True,
                    inspiral=True, a0=a)
    elif name == "What gravitates: counter-streaming beams (F7.4a)":
        lam, half = 1.5e-4, 3.0  # weak enough that the nonlinear metric factor (1 - 2U) stays ~1 %
        cols = ("#b0b0b0", "#7ad0ff", "#7affb0", "#ff9a7a", "#fff27a")
        for k, v in enumerate(ROD_SPEEDS + (None,)):
            x = 15.0 * (k - 2)  # side by side along x: the field of a neighbour is mostly horizontal and ~1 % weak
            rod = Body(f"beams v = {v:g}c" if v is not None else "light beams", x, 0.0, kind="rod",
                       color=cols[k])
            rod.lam, rod.half = lam, half
            rod.vb, rod.light = (v, False) if v is not None else (1.0, True)
            B.append(rod)
            pr = Body(f"probe {k + 1}", x, 1.0, mass=1e-14, color=cols[k])
            pr.rod = rod
            B.append(pr)
        view = (0.0, 0.5, 34.0)
        meta.update(dt_max=0.5, eta=0.01, speed=math.log10(6.0), display=("log", -4.0, -2.0), rods=True,
                    click="rest", mass=(1e-12, "model units (G = c = 1)"))
    else:  # free scene: a fixed source keeps added orbits stable
        B.append(Body("mass (fixed)", 0, 0, mass=0.01, fixed=True, color="#ffd27a"))
        view = (0.0, 0.0, 3.0)
    return B, view, meta


def category_si(kg):
    """Name of an object by its mass in kilograms."""
    for lim, name in ((1.6e29, "star"), (2.5e28, "brown dwarf"), (1.9e26, "giant planet"), (1e23, "planet"),
                      (1e19, "moon"), (1e15, "asteroid"), (1e8, "boulder / building"), (1e3, "car / aircraft"),
                      (5.0, "stone"), (0.05, "apple")):
        if kg >= lim:
            return name
    return "pebble"


def category_relative(r):
    """Name of an object by its mass relative to the reference (heaviest) mass of a geometric scene, read as a star."""
    for lim, name in ((0.08, "star"), (0.013, "brown dwarf"), (1e-4, "giant planet"), (3e-7, "planet"),
                      (1e-10, "moon")):
        if r >= lim:
            return name
    return "asteroid"


# =============================================================================
# Experiments: Mercury and the classical tests of gravity in every framework
# =============================================================================


def mercury_experiment(laws, years, planets):
    """
    Generator: perihelion advance of Mercury in arcseconds per century under each law over `years` years, with or
    without the other planets (J2000, 2D). Yields (law name, fraction done) while running and returns {law key: rate}.
    The rate is the slope of the Laplace-Runge-Lenz direction at the perihelia (PerihelionTracker). Expected
    (observation and GR): ~43''/century from relativity alone; the planets add several hundred arcseconds of
    Newtonian perturbation (the real 3D value is ~532''/century; this 2D model differs).
    """
    out = {}
    for law in laws:
        sim = Sim(law)
        sim.bodies = solar_system(planets=planets, only=None if planets else ["Mercury"])
        sim.eta, sim.dt_max = 0.004, 2.0 * 86400 * C_SI
        sun, mer = sim.bodies[0], next(b for b in sim.bodies if b.name == "Mercury")
        tr = PerihelionTracker(mer, sun)
        T = years * YEAR * C_SI
        k = 0
        while sim.t < T:
            live, pos, vel, mass, photon, _ = sim.arrays()
            sim.step(max(min(sim.choose_dt(pos, vel, mass, photon), T - sim.t), 1e-9))
            tr.update(sim)
            k += 1
            if k % 200 == 0:
                yield law.name, sim.t / T
        r = tr.rate()
        out[law.key] = None if r is None else r * 100 * YEAR * C_SI / ARCSEC
    return out


def _march(sim, until, on_step=None, chunk=250, max_steps=2_000_000):
    """Steps sim until until(sim) is true; yields every `chunk` steps so a GUI stays responsive."""
    n = 0
    while not until(sim) and n < max_steps:
        sim.step()
        if on_step:
            on_step(sim)
        n += 1
        if n % chunk == 0:
            yield None


def ct_light(law):
    """Light grazing the Sun (b = R_sun) from -2000 R_sun to +2000 R_sun; deflection [arcsec]."""
    sim = Sim(law)
    sun = Body("Sun", 0, 0, mass=M_SUN, fixed=True)
    X = 2000.0 * R_SUN
    ph = Body("photon", -X, R_SUN, 1.0, 0.0, kind="photon")
    sim.bodies = [sun, ph]
    sim.eta, sim.dt_max = 0.01, 1e10
    tr = DeflectionTracker(ph, sun, sim)
    yield from _march(sim, lambda s: ph.x >= X or not ph.alive, tr.update)
    return tr.state(sim)["defl"] / ARCSEC


def ct_probe(law):
    """A probe at v = 0.5c past M = 1e-4 at b = 1 (path -400 .. +400): deflection divided by Newton's."""
    sim = Sim(law)
    m = Body("mass", 0, 0, mass=1e-4, fixed=True)
    p = Body("probe", -400.0, 1.0, 0.5, 0.0, mass=1e-15)
    p.probe = True
    sim.bodies = [m, p]
    sim.eta, sim.dt_max = 0.005, 5.0
    tr = DeflectionTracker(p, m, sim)
    yield from _march(sim, lambda s: p.x >= 400.0, tr.update)
    st = tr.state(sim)
    return st["defl"] / st["newton"]


def ct_mercury(law, years=2.0):
    """Mercury without the planets: relativistic perihelion advance [arcsec/century]."""
    gen = mercury_experiment([law], years, False)
    try:
        while True:
            next(gen)
            yield None
    except StopIteration as stop:
        return stop.value.get(law.key)


def ct_shapiro(law):
    """Excess one-way travel time of a signal Earth (1 AU) -> Cassini (8.43 AU) passing the Sun at 1.6 R_sun [us]."""
    sim = Sim(law)
    sun = Body("Sun", 0, 0, mass=M_SUN, fixed=True)
    ph = Body("signal", -AU, 1.6 * R_SUN, 1.0, 0.0, kind="photon")
    sim.bodies = [sun, ph]
    sim.eta, sim.dt_max = 0.01, 2e10
    tr = DeflectionTracker(ph, sun, sim)
    x2 = 8.43 * AU
    yield from _march(sim, lambda s: ph.x >= x2, tr.update)
    # the last step may overshoot x2; the excess time is measured along the actual path, so it does not matter
    return tr.state(sim)["delay"] / C_SI * 1e6


def ct_redshift(law):
    """Pound-Rebka: fractional rate difference of clocks at 0 and 22.5 m above the ground [1e-15]."""
    sim = Sim(law)
    sim.bodies = [Body("Earth", 0, 0, mass=M_EARTH, fixed=True, radius=R_EARTH),
                  Body("bottom", 0, R_EARTH + 0.01, kind="clock"), Body("top", 0, R_EARTH + 22.51, kind="clock")]
    sim.dt_max = 0.1 * C_SI
    sim.advance(1.0 * C_SI, max_steps=100)
    yield None
    return (sim.bodies[1].tau_def - sim.bodies[2].tau_def) / sim.t * 1e15


def ct_gpb(law):
    """Gravity Probe B: gyroscope on a polar circular orbit of radius 7027.4 km, one full orbit. Returns (geodetic,
    frame-dragging) drift [mas/yr]."""
    sim = Sim(law)
    earth = Body("Earth", 0, 0, mass=M_EARTH, fixed=True, radius=R_EARTH)
    earth.J = (0.0, J_EARTH, 0.0)
    r = GPB_A
    gy = make_gyro(Body("GP-B", r, 0.0, 0.0, math.sqrt(M_EARTH / r), mass=1e-24), s0=(1.0, 0.0))
    sim.bodies = [earth, gy]
    sim.eta, sim.dt_max = 0.01, 20.0 * C_SI
    T = 2 * math.pi * math.sqrt(r**3 / M_EARTH)
    yield from _march(sim, lambda s: s.t >= T)
    g = gy.gyro
    k = YEAR * C_SI / MAS / g["t"]
    return g["geo"] * k, abs(g["fd"]) * k


def ct_pulsar(law):
    """PSR B1913+16 with radiation reaction, two orbits: period derivative dP/dt [dimensionless]."""
    ht = HULSE_TAYLOR
    sim = Sim(law)
    m1, m2 = ht["m_p"] * M_SUN, ht["m_c"] * M_SUN
    P = ht["Pb_d"] * 86400.0 * C_SI
    a = ((m1 + m2) * (P / (2 * math.pi)) ** 2) ** (1.0 / 3.0)
    sim.bodies = list(binary_bodies("pulsar", m1, "#fff", "companion", m2, "#fff", a, ht["e"]))
    sim.radiation, sim.eta, sim.dt_max = True, 0.005, 60.0 * C_SI
    tr = BinaryTracker(sim.bodies[0], sim.bodies[1], sim, start_at_periastron=True)
    yield from _march(sim, lambda s: s.t >= 2.02 * P, tr.update)
    res = tr.result()
    return None if res is None else res[2]


def ct_beams(law, light=False):
    """Active mass of two counter-streaming beams at 0.6c (or of light beams per unit energy) relative to the same
    beams at rest, read from the field of the rods."""
    sim = Sim(law)
    pt = np.array([[0.0, 1.0]])
    rest = Body("rest", 0, 0, kind="rod")
    rest.lam, rest.half, rest.vb = 1e-4, 3.0, 0.0
    mov = Body("moving", 0, 0, kind="rod")
    mov.lam, mov.half, mov.vb, mov.light = 1e-4, 3.0, (1.0 if light else 0.6), light
    sim.bodies = [rest]
    g0 = sim.rod_field(pt)[1][0, 1]
    sim.bodies = [mov]
    g1 = sim.rod_field(pt)[1][0, 1]
    yield None
    return g1 / g0


def metric_b_crit(law, M=1.0):
    """Critical impact parameter b = min r sqrt(B/A) of the static metric of one mass (U = M/r), or None when there is
    no photon sphere (no interior minimum). The scalar stage uses its exact profile."""
    if law.kind == "newton":
        return None
    if law.key == "H1":
        return bh_profile()[3]["b_crit"] * M
    r = np.geomspace(0.05 * M, 50.0 * M, 20001)
    U = M / r
    A = law.A(U)
    ok = A > 0
    if not np.any(ok):
        return None
    f = np.where(ok, r * np.sqrt(np.maximum(law.B(U), 0.0) / np.maximum(A, 1e-300)), np.inf)
    k = int(np.argmin(f))
    if k == 0 or k == len(f) - 1 or not np.isfinite(f[k]) or not ok[k - 1]:
        return None
    return float(f[k])


def ct_shadow(law):
    yield None
    return metric_b_crit(law)


def classical_test_rows():
    """The rows of the classical-tests table: (key, label, unit, scale, observation)."""
    ht = HULSE_TAYLOR
    return [
        ("light", "light bending at the solar limb", "″", 1.0, "1.7512 (VLBI: γ − 1 = (−0.8 ± 1.2)×10⁻⁴)"),
        ("probe", "bending of a probe at 0.5c ÷ Newton", "×", 1.0, "1 + γv² = 1.25 for γ = 1"),
        ("mercury", "Mercury perihelion, relativistic part", "″/cy", 1.0, "42.98 ± 0.04"),
        ("shapiro", "Shapiro delay Earth–Cassini at 1.6 R☉", "μs", 1.0, "GR × (1 ± 2×10⁻⁵) (Cassini 2003)"),
        ("redshift", "Pound–Rebka redshift over 22.5 m", "10⁻¹⁵", 1.0, "2.458 × (0.999 ± 0.008)"),
        ("gpb_geo", "Gravity Probe B geodetic", "mas/yr", 1.0, "6601.8 ± 18.3 (GR 6606.1)"),
        ("gpb_fd", "Gravity Probe B frame dragging", "mas/yr", 1.0, "37.2 ± 7.2 (GR 39.2)"),
        ("pulsar", "PSR B1913+16 dP/dt", "10⁻¹²", 1e12, f"{ht['pdot_obs'] * 1e12:.3f} ± {ht['pdot_obs_err'] * 1e12:.3f}"),
        ("beams", "active mass, beams at 0.6c ÷ rest", "×", 1.0, "lab F7.4a: 2.401 ± 0.053 (exact rule 2.361)"),
        ("lightbeam", "active mass, light ÷ same energy at rest", "×", 1.0, "GR: 2 (Tolman–Ehrenfest–Podolsky)"),
        ("shadow", "black-hole shadow b_crit / M", "M", 1.0, "EHT: consistent with GR (~10 %)"),
    ]


def classical_tests(laws):
    """
    Generator over all classical tests and all frameworks. Yields (text, fraction done) and returns
    {test key: {law key: value or None}}.
    """
    tests = [("light", ct_light), ("probe", ct_probe), ("mercury", ct_mercury), ("shapiro", ct_shapiro),
             ("redshift", ct_redshift), ("gpb", ct_gpb), ("pulsar", ct_pulsar), ("beams", ct_beams),
             ("lightbeam", lambda law: ct_beams(law, light=True)), ("shadow", ct_shadow)]
    res = {key: {} for key, *_ in classical_test_rows()}
    n_total = len(tests) * len(laws)
    done = 0
    for key, fn in tests:
        for law in laws:
            gen = fn(law)
            try:
                while True:
                    next(gen)
                    yield f"{key}: {law.short}", done / n_total
            except StopIteration as stop:
                val = stop.value
            except (ValueError, FloatingPointError, ZeroDivisionError, OverflowError):
                val = None
            if key == "gpb":
                res["gpb_geo"][law.key], res["gpb_fd"][law.key] = val if val is not None else (None, None)
            else:
                res[key][law.key] = val
            done += 1
    return res


def classical_tests_table(res, laws):
    """Readout lines of the classical-tests table."""
    cols = [law.short for law in laws]
    w = 11
    head = f"{'test':44s}" + "".join(f"{c[:w]:>{w + 1}s}" for c in cols) + "   observed"
    lines = ["CLASSICAL TESTS: every framework integrated by this simulator (same code, same scenes)", head,
             "-" * len(head)]
    for key, label, unit, scale, obs in classical_test_rows():
        cells = []
        for law in laws:
            v = res.get(key, {}).get(law.key)
            if v is None:
                cells.append("none" if key == "shadow" else "–")
            else:
                v = v * scale
                cells.append("0" if v == 0 else (f"{v:.4g}" if abs(v) >= 1e-3 else f"{v:.2e}"))
        lines.append(f"{(label + ' [' + unit + ']'):44s}" + "".join(f"{c:>{w + 1}s}" for c in cells) + f"   {obs}")
    lines += [
        "Reading: causal gravity and Einstein agree in every row (decision 13: the same metric, tensor coupling and "
        "radiation); differences could only appear at the discreteness scale of the order.",
        "Newton fails light (half), perihelion, Shapiro, clocks, gyroscopes and radiation; the scalar stage H1 passes "
        "light, perihelion, Shapiro and clocks but fails frame dragging (0), the binary pulsar (1/6 for circular "
        "orbits) and the active mass of moving matter; the fixed order F3 (γ = −1) fails light bending completely.",
        "Simplifications: 2D, test-body dynamics, finite paths (light: ±2000 R☉), GP-B without the Earth's "
        "oblateness (simplified GR values ~6604 and ~40.8 mas/yr against the mission's 6606.1 and 39.2), scalar radiation as 1/6 of the quadrupole.",
    ]
    return lines


# =============================================================================
# State of the theory, about text, prediction board
# =============================================================================

STATUS = [
    ("Newton's law, Gauss flux of a source",
     "derived from the discrete operator and tested in the lab (F1.3, K.12); in the simulator INPUT (U = ΣM/r with "
     "given masses), not derived"),
    ("Coupling κ ↔ G", "input (calibration), not derived"),
    ("Motion along the geodesics of the field", "tested in 2+1 (F2.7: motion exactly on the geodesic of the causal-set "
     "field)"),
    ("Discrete operator → continuum equation", "tested: K.11, K.13a; K.12 without screening (r(0) = 1)"),
    ("γ = 1 (light bending, Shapiro delay)",
     "derived from the throughput principle (decision 12); lab F7.1c in 3+1: 1.000 / 0.972 and 0.948 / 0.979 "
     "(formally undecided only by the statistics of the control); no longer an input"),
    ("Tensor coupling to T_ab (frame dragging, pressure and motion gravitate)",
     "tested: F7.4a PASSED (beams at 0.6c: R = 2.401 ± 0.053 against 2.361; anisotropy K = 1.221 ± 0.012 against "
     "1.192; the cones tilt along the flow)"),
    ("Nonlinear completion = general relativity (strong field, horizon)",
     "decision 13 (route A): input until derived on the substrate (Jacobson's route, steps 2 and 4)"),
    ("β = 1 (perihelion)", "continuum (ln Π, F1.7); in GR form by decision 13; discrete F1.9a: instrument failed in the "
     "quadrature check, descriptively c → continuum (β ≈ 1 at σ = 16)"),
    ("Gravitational radiation (quadrupole)",
     "follows from the tensor coupling at linear order; the simulator uses the 2.5PN reaction force of GR; not tested "
     "in the lab"),
    ("Mass as point particles", "input of the simulator; the origin of matter from the dynamics is open"),
    ("Superposition, quasi-static field", "linear superposition of point fields, no retardation (slow sources)"),
    ("Many-body dynamics", "every body a test body in the field of the others; gravitomagnetism of mass currents and "
     "the EIH terms are not modelled (spins only)"),
    ("a₀, flat rotation curves", "not in the model (not derived)"),
    ("Earth: fall, clocks, GPS", "weak field: the same as GR and the observations (Pound–Rebka, GPS); nothing new"),
    ("UAP 'null field'", "HYPOTHESIS of the author, neither in the theory nor in the lab; the scene only shows what it "
     "would have to explain"),
    ("Scalar stage: no horizon, throat 2.341 M, shadow +0.91 %",
     "NOT a prediction of the theory (decision 13); shown only as the framework H1"),
]

ABOUT = """Causal gravity (decision 13, 28 Oct 2026)

Throughput principle: where the light cones narrow (c_eff ↓)
the levels of the cone saturate (more events), ρ_c · c = const.
In 3+1 this gives γ = 1 (F7.1c) and, in the linear order, the
coupling of the field to the full energy-momentum tensor:
   h̄_ab = 4U u_a u_b  (dust),
realised on the substrate by "every flow narrows the cone in
its own rest frame, capacity deficits add" (F7.4a PASSED).
The nonlinear completion is general relativity.

Metric used (G = c = 1, U = ΣM/r):
   ds² = −A dt² + B dx²,  A = ((1−U/2)/(1+U/2))²,  B = (1+U/2)⁴
Spinning bodies: g_0i = −2 (J × r)_i / r³ (frame dragging).
Radiation: 2.5PN reaction force of a binary (quadrupole).
Moving matter as a source: T_00 + T_ii, beams Γ²(1 + v²).

Motion (geodesics of the static metric):
   dv/dt = −(A′/2B)∇U − v(v·∇U)B′/B + v²(B′/2B)∇U + v(v·∇U)A′/A
         + v × (∇ × g_0i)
Proper time: dτ = √(A − B v²) dt.

Measured against (1PN, test body, one mass):
   perihelion (2 + 2γ − β)/3 · 6πM/(a(1−e²)),
   light (1+γ)(M/b)(sin θ₂ − sin θ₁) + π(8+8γ−4β+3ε)/4 (M/b)²,
   probe (1 + γv²) × Newton,  Shapiro (1+γ) M ln((x₂+r₂)/(x₁+r₁)),
   gyroscope (γ + ½) v × ∇U and (1+γ)/2 [3(J·n)n − J]/r³,
   binary dP/dt = −(192π/5)(2πM_c/P)^(5/3) f(e).

Honest reading: in this simulator causal gravity and Einstein
coincide by construction; the content is that the substrate
rule produces the tensor coupling (F7.4a) and that Newton, the
scalar stage and the fixed order fail the listed tests."""

# Risky claims of the theory and the tests that measure them. "reading" is the interpretation in docs/nalezy.md
# with its date; the status of the tests is loaded from the saved records.
PREDICTIONS = [
    ("The far field of a source is Newtonian (Gauss flux)", ["F1.3"], "a different exponent or strength",
     "holds: 1/r field with the flux of the charge (F1.3)"),
    ("Discrete operator → □_g − R/2 (also in H1)", ["K.11", "K.11c", "K.13a"], "no convergence to the continuum form",
     "holds at first order (K.11, K.13a, 0.8–2.4 % at σ = 32)"),
    ("The operator does not screen the field (2+1)", ["K.12"], "a floor independent of σ after the instrument",
     "24–25 Oct: screening excluded analytically (r(0) = 1); nonlinear remainder ~0.06 open"),
    ("Motion along the geodesics of the causal-set field", ["F2.7"], "a path off the geodesic",
     "motion exactly on the geodesic (2.027 against 2.023)"),
    ("γ = 1: light bending from the substrate (throughput principle)", ["K.16c", "F7.1a", "F7.1b", "F7.1c"],
     "γ ≠ 1 (Cassini: |γ − 1| < 2.3·10⁻⁵)",
     "3+1: γ ≈ 1 in two runs (1.000 / 0.972, 0.948 / 0.979); cones only 1/3 as predicted; formally undecided by the "
     "statistics of the control (27 Oct)"),
    ("Tensor coupling: moving matter gravitates as T_00 + T_ii, cones tilt with the flow", ["F7.4a"],
     "active mass Γ², Γ or 1 instead of Γ²(1 + v²); isotropic cones",
     "PASSED 28 Oct: R = 2.401 ± 0.053 (exact rule 2.361, energy-only 1.5625 is 16 SE away), K = 1.221 ± 0.012"),
    ("β = 1 in the discrete field", ["F1.9a"], "c away from 0 by more than 1/4 (β = 1/2 or 3/2)",
     "F1.9a INSTRUMENT FAILED (quadrature σ = 4); descriptively c(σ = 16) = −0.002 and +0.055 against +0.008 and "
     "+0.058"),
    ("Matter arises from the dynamics of the causal order", [], "no stable structure", "not measured: the largest gap"),
    ("a₀ = cH₀/2π from the dynamics (flat rotation curves)", ["F6.7"], "another scale or evolution of a₀(z)",
     "not derived (analysis F6.7, 23 Oct)"),
    ("Gravitational waves: tensor modes, speed c, quadrupole formula", [],
     "scalar or dipole radiation (Hulse–Taylor pulsar, LIGO)",
     "follows from the tensor coupling (F7.4a) at linear order; the dynamics of the order is not yet derived"),
    ("Strong field: horizon, photon sphere, shadow", [], "shadow against EHT (M87*, Sgr A*)",
     "general relativity by decision 13 (horizon 2M, photon sphere 3M, b_crit 5.196 M); derivation pending "
     "(Jacobson's route); the scalar-stage throat is withdrawn as a prediction"),
]


# =============================================================================
# Bridge to the laboratory: records, prediction board and the continuum reference
# =============================================================================


def default_results_dir():
    """results/ of the laboratory next to demo/ (D:\\Projekty\\causality_lab\\results when the simulator is in demo/)."""
    import os

    return os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results"))


def load_records(folder, include_quick=False):
    """Saved test records (JSON), newest first. Damaged files are skipped."""
    import json
    import os

    recs = []
    if not os.path.isdir(folder):
        return recs
    for tid in sorted(os.listdir(folder)):
        sub = os.path.join(folder, tid)
        if not os.path.isdir(sub):
            continue
        paths = [os.path.join(sub, f) for f in os.listdir(sub) if f.endswith(".json")]
        if include_quick and os.path.isdir(os.path.join(sub, "quick")):
            paths += [os.path.join(sub, "quick", f) for f in os.listdir(os.path.join(sub, "quick"))
                      if f.endswith(".json")]
        for pth in paths:
            try:
                with open(pth, encoding="utf-8") as fh:
                    rec = json.load(fh)
            except (OSError, ValueError):
                continue
            if not isinstance(rec, dict):
                continue
            rec["_path"] = pth
            recs.append(rec)
    recs.sort(key=lambda r: str(r.get("time_utc") or ""), reverse=True)
    return recs


def continuum_a_2p1(R, sigma, q, n_grid=20001):
    """Continuum solution of equation (a) in 2+1: Laplace Pi = (1/3) Pi^2 rho_s with Pi(R) = 1 and kappa = 1, by
    relaxation with a gradually switched-on source. An independent reference for the K.12/F2.7 profiles."""
    r = np.linspace(0.0, R, n_grid)
    rs = q * np.exp(-(r**2) / (2 * sigma**2))

    def poisson(f):  # Laplace Psi = f in radial symmetry; Psi'(0) = 0 and Psi(R) = 0
        dens = f * 2 * math.pi * r
        Q = np.concatenate([[0.0], np.cumsum(0.5 * (dens[1:] + dens[:-1]) * np.diff(r))])
        du = np.zeros_like(r)
        du[1:] = Q[1:] / (2 * math.pi * r[1:])
        U = np.concatenate([[0.0], np.cumsum(0.5 * (du[1:] + du[:-1]) * np.diff(r))])
        return U - U[-1]

    Pi = np.ones_like(r)
    for t in (0.125, 0.25, 0.5, 0.75, 1.0):
        for _ in range(3000):
            new = 1.0 + poisson(t * Pi**2 * rs / 3.0)
            step = float(np.max(np.abs(new - Pi)))
            Pi = 0.5 * Pi + 0.5 * new
            if step < 1e-11:
                break
    return r, Pi


# =============================================================================
# Self-test of the simulator code without a window
# =============================================================================


def _finish(gen):
    """Runs a test generator to the end and returns its value."""
    try:
        while True:
            next(gen)
    except StopIteration as stop:
        return stop.value


def selftest():
    """Checks the numerical implementation of the formulas; it is not a proof of the physical theory."""
    ok_all = True

    def report(name, val, ref, tol):
        nonlocal ok_all
        rel = abs(val - ref) / abs(ref)
        ok = rel <= tol
        ok_all &= ok
        print(f"  {'OK' if ok else 'FAIL':5s} {name}: {val:.6g}; expected {ref:.6g} (relative deviation {rel:.2e}, "
              f"tolerance {tol:g})")

    print(f"Self-tests of causal_gravity.py ({APP_NAME} {APP_VERSION}; units G = c = 1)")
    # 1. Newtonian limit: period of a nearly circular orbit around a small mass.
    s = Sim()
    M = 1e-5
    s.bodies = [Body("c", 0, 0, mass=M, fixed=True), Body("p", 1.0, 0, 0, math.sqrt(M), mass=0.0)]
    crossings, prev_y = [], 0.0
    p = s.bodies[1]

    def on(sm):
        nonlocal prev_y
        if prev_y < 0 <= p.y and p.x > 0:
            crossings.append(sm.t)
        prev_y = p.y

    s.advance(2.2 * 2 * math.pi / math.sqrt(M), max_steps=200000, on_step=on)
    if len(crossings) >= 2:
        report("Kepler period (M = 1e-5, r = 1)", crossings[1] - crossings[0], 2 * math.pi / math.sqrt(M), 1e-3)
    else:
        print("  FAIL  Kepler period: the orbit did not close")
        ok_all = False
    # 2. Perihelion advance against 6 pi M / (a (1 - e^2)); higher orders are small.
    s = Sim()
    M, a, e = 0.002, 1.0, 0.5
    s.bodies = [Body("c", 0, 0, mass=M, fixed=True), Body("p", a * (1 - e), 0, 0, kepler_peri_velocity(M, a, e))]
    tr = PerihelionTracker(s.bodies[1], s.bodies[0])
    s.advance(5.3 * 2 * math.pi * math.sqrt(a**3 / M), max_steps=400000, on_step=tr.update)
    res = tr.result()
    if res:
        report(f"perihelion advance per orbit ({res[4]} orbits)", res[0], res[1], 0.03)
    else:
        print("  FAIL  perihelion advance: no perihelion found")
        ok_all = False
    # 3.-4. Light bending and Shapiro delay against first order + the second order of the Schwarzschild metric.
    s = Sim()
    s.bodies = [Body("m", 0, 0, mass=1e-3, fixed=True), Body("ph", -50.0, 0.2, 1.0, 0.0, kind="photon")]
    lt = DeflectionTracker(s.bodies[1], s.bodies[0], s)
    s.advance(100.0, max_steps=400000, on_step=lt.update)
    st = lt.state(s)
    report("light bending (M = 1e-3, b = 0.2)", st["defl"], st["pred"], 2e-3)
    report("Shapiro delay", st["delay"], st["delay_pred"], 0.03)
    # 5. Null condition of the photon and the rate of a static clock.
    ph = s.bodies[1]
    U = -s.L_at(ph.x, ph.y)
    law = s.law
    report("photon speed |v| = sqrt(A/B)", math.hypot(ph.vx, ph.vy), math.sqrt(law.A(U) / law.B(U)), 1e-9)
    s = Sim()
    s.bodies = [Body("m", 0, 0, mass=0.02, fixed=True), Body("clk", 0.5, 0, kind="clock")]
    s.advance(10.0)
    U = 0.02 / math.hypot(0.5, EPS)
    report("static clock rate dtau/dt = sqrt(A)", s.bodies[1].tau / s.t, (1 - U / 2) / (1 + U / 2), 1e-9)
    # 6. Earth: fall time from 10 m against sqrt(2h/g).
    s = Sim()
    s.dt_max = 0.002 * C_SI
    h = 10.0
    s.bodies = [Body("E", 0, 0, mass=M_EARTH, fixed=True, radius=R_EARTH), Body("a", 0, R_EARTH + h, mass=1e-28)]
    while s.bodies[1].landed_t is None and s.t < 5 * C_SI:
        s.step()
    g = M_EARTH / (R_EARTH + h) ** 2 * C_SI**2
    report("fall time from 10 m [s]", s.bodies[1].landed_t / C_SI, math.sqrt(2 * h / g), 2e-3)
    # 7. GPS: rate difference of an orbiting clock and a ground clock per day.
    s = Sim()
    s.dt_max = 10.0 * C_SI
    r = 2.6560e7
    s.bodies = [Body("E", 0, 0, mass=M_EARTH, fixed=True, radius=R_EARTH), Body("g", 0, R_EARTH + 0.01, kind="clock"),
                Body("GPS", r, 0, 0, math.sqrt(M_EARTH / r), mass=1e-24)]
    s.advance(86400 * C_SI, max_steps=100000)
    d_meas = (s.bodies[1].tau_def - s.bodies[2].tau_def) / s.t * 86400e6
    d_pred = (M_EARTH / (R_EARTH + 0.01) - M_EARTH / r - 0.5 * M_EARTH / r) * 86400e6
    report("GPS minus ground clock [us/day]", d_meas, d_pred, 1e-2)
    # 8. Earth: clock 30 m higher, the gravitational redshift.
    report("Pound-Rebka redshift over 22.5 m [1e-15]", _finish(ct_redshift(MAIN_LAW())),
           M_EARTH / (R_EARTH + 0.01) * 1e15 - M_EARTH / (R_EARTH + 22.51) * 1e15, 1e-3)
    # 9. Mercury: relativistic perihelion advance in every framework without the planets.
    res = _finish(mercury_experiment([_law_newton(), _law_cg(), _law_gr(), _law_h1(), _law_ppn(0.5, 1.0)], 3.0,
                                     False))
    report("Mercury, Newton without planets [''/century]", abs(res["N"]) + 1.0, 1.0, 0.05)
    report("Mercury, causal gravity [''/century]", res["CG"], 42.98, 0.01)
    report("Mercury, GR [''/century]", res["GR"], 42.98, 0.01)
    report("Mercury, scalar stage H1 [''/century]", res["H1"], 42.98, 0.01)
    report("Mercury, PPN beta = 1/2 [''/century]", res["PPN"], 42.98 * 3.5 / 3.0, 0.01)
    # 10. Black hole: b = 5.23 escapes in causal gravity (Schwarzschild) and is captured in the scalar stage.
    B0, _, m = scene("Black hole (strong field)")
    fate = {}
    for law in (_law_cg(), _law_h1()):
        B, _, m = scene("Black hole (strong field)")
        s = Sim(law)
        s.bodies = [B[0]] + [b for b in B if b.name in ("photon b=5.23", "photon b=5.15")]
        s.dt_max, s.eta, s.no_guard = m["dt_max"], m["eta"], True
        if law.key == "H1":
            s.profile = B[0]
        s.capture = (B[0], 1.0)
        while s.t < 250 and any(not b.captured and math.hypot(b.x, b.y) < 45 for b in s.bodies[1:]):
            s.step()
        fate[law.key] = {b.name: b.captured for b in s.bodies[1:]}
    report("black hole: b = 5.23 escapes (causal) and is captured (scalar stage); b = 5.15 captured (causal)",
           1.0 if (not fate["CG"]["photon b=5.23"] and fate["H1"]["photon b=5.23"] and fate["CG"]["photon b=5.15"])
           else 0.0, 1.0, 1e-9)
    report("black hole: critical impact parameter, causal = Schwarzschild [M]", metric_b_crit(_law_cg()),
           3 * math.sqrt(3), 1e-4)
    report("black hole: critical impact parameter, scalar stage [M]", bh_profile()[3]["b_crit"], 5.2435, 1e-3)
    # 11. Lense-Thirring: prograde minus retrograde perihelion advance around a spinning mass (equatorial orbits).
    M, a, e, J = 0.001, 0.3, 0.3, 1e-6
    adv = {}
    for sense in (1, -1):
        s = Sim()
        c = Body("c", 0, 0, mass=M, fixed=True)
        c.J = (0.0, 0.0, J)
        pl = Body("p", sense * a * (1 - e), 0, 0, kepler_peri_velocity(M, a, e), mass=0.0)
        s.bodies = [c, pl]
        s.eta = 0.01
        tr = PerihelionTracker(pl, c)
        s.advance(4.3 * 2 * math.pi * math.sqrt(a**3 / M), max_steps=400000, on_step=tr.update)
        r_ = tr.result()
        adv[sense] = r_[0] if r_ else float("nan")
    report("Lense-Thirring: prograde minus retrograde advance per orbit [rad]", adv[1] - adv[-1],
           -16 * math.pi * J / (math.sqrt(M * a**3) * (1 - e * e) ** 1.5), 0.05)
    # 12. Gravity Probe B: one polar orbit; geodetic and frame-dragging drift against the analytic averages.
    geo, fd = _finish(ct_gpb(MAIN_LAW()))
    r = GPB_A
    k = YEAR * C_SI / MAS
    report("Gravity Probe B geodetic [mas/yr]", geo, 1.5 * M_EARTH**1.5 / r**2.5 * k, 1e-2)
    report("Gravity Probe B frame dragging [mas/yr]", fd, J_EARTH / (2 * r**3) * k, 1e-2)
    # 13. Radiation reaction: radiated energy over whole orbits against the quadrupole formula (Peters).
    s = Sim()
    a, e = 200.0, 0.5
    s.bodies = list(binary_bodies("A", 0.5, "#fff", "B", 0.5, "#fff", a, e))
    s.radiation, s.eta = True, 0.005
    bt = BinaryTracker(s.bodies[0], s.bodies[1], s, start_at_periastron=True)
    P = 2 * math.pi * math.sqrt(a**3)
    s.advance(2.02 * P, max_steps=2_000_000, on_step=bt.update)
    P_m, lum, pdot, a_m, e_m = bt.result()
    # the 1PN dynamics turns the Newtonian initial state into a slightly different orbit: compare with the measured
    # period and eccentricity
    report(f"quadrupole luminosity (measured e = {e_m:.3f})", lum, peters_lum(0.5, 0.5, a_m, e_m), 2e-2)
    report("period derivative dP/dt", pdot, peters_pdot(0.5, 0.5, P_m, e_m), 2e-2)
    # 14. Fast probes: bending relative to Newton, 1 + gamma v^2.
    report("probe at 0.5c, causal gravity (1 + v^2)", _finish(ct_probe(MAIN_LAW())), 1.25, 1e-2)
    report("probe at 0.5c, fixed order F3 (1 - v^2)", _finish(ct_probe(_law_f3())), 0.75, 1e-2)
    # 15. Counter-streaming beams (F7.4a): active mass relative to rest.
    report("beams at 0.6c, causal gravity Gamma^2 (1 + v^2)", _finish(ct_beams(MAIN_LAW())), 2.125, 1e-9)
    report("beams at 0.6c, scalar stage (rest count)", _finish(ct_beams(_law_h1())), 1.0, 1e-9)
    print("RESULT:", "all checks passed" if ok_all else "some checks FAILED")
    return 0 if ok_all else 1


# =============================================================================
# Application and its Tkinter user interface
# =============================================================================

BG, PANEL, TXT, DIM, ACCENT = "#14161c", "#1c1f27", "#d8dce6", "#9aa0ad", "#9ad9ff"


def run_app(start_scene=None, autorun=False):
    import tkinter as tk
    from tkinter import ttk

    root = tk.Tk()
    root.title(f"{APP_NAME} – interactive simulator of the gravity of causal physics ({APP_VERSION})")
    root.geometry("1400x880")
    root.configure(bg=BG)

    state = {
        "sim": Sim(), "view": (0.0, 0.0, 3.0), "running": False, "trackers": [], "light": [], "binaries": [],
        "ghosts": [], "exp": None, "exp_result": None, "tests": None, "tests_result": None, "tests_laws": None,
        "tests_status": "", "hover": None, "drag": None, "pan": None, "frame": 0, "field_img": None,
        "scene": SCENES[0], "meta": {"si": False, "display": ("log", -4.0, -1.0), "click": "orbit"},
        "name_seq": 0, "mass_reference": None, "view3d": False, "readout_hold": False, "readout_text": "",
        "readout_pending": "", "mass_unit_prev": None, "achieved": None,
    }

    # ---------- layout
    left = tk.Frame(root, bg=PANEL, width=340)
    left.pack(side="left", fill="y")
    left.pack_propagate(False)
    right = tk.Frame(root, bg=BG)
    right.pack(side="left", fill="both", expand=True)
    canvas = tk.Canvas(right, bg="#0b0d12", highlightthickness=0)
    canvas.pack(fill="both", expand=True)
    out_frame = tk.Frame(right, bg="#10131a")
    out_frame.pack(fill="x")
    out_toolbar = tk.Frame(out_frame, bg=PANEL)
    out_toolbar.pack(side="top", fill="x")
    out_status = tk.Label(out_toolbar, text="Measurement readout", bg=PANEL, fg=DIM, anchor="w")
    out_status.pack(side="left", fill="x", expand=True, padx=6)
    out_text_frame = tk.Frame(out_frame, bg="#10131a")
    out_text_frame.pack(side="top", fill="both", expand=True)
    out_sb = tk.Scrollbar(out_text_frame, orient="vertical")
    out_xsb = tk.Scrollbar(out_text_frame, orient="horizontal")
    out_sb.pack(side="right", fill="y")
    out_xsb.pack(side="bottom", fill="x")
    out = tk.Text(out_text_frame, height=12, bg="#10131a", fg=TXT, font=("Consolas", 10), relief="flat", wrap="none",
                  yscrollcommand=out_sb.set, xscrollcommand=out_xsb.set, undo=False, exportselection=True)
    out.pack(side="left", fill="both", expand=True)
    out_sb.config(command=out.yview)
    out_xsb.config(command=out.xview)
    err_lbl = tk.Label(right, text="", bg=BG, fg="#ff6b6b", anchor="w", justify="left", wraplength=1000)
    err_lbl.pack(fill="x")
    seen_errors = set()

    def report_error(where, exc):
        """Shows an error under the readout and keeps the simulation running; the traceback is printed only once."""
        import traceback

        err_lbl.config(text=f"Error ({where}): {type(exc).__name__}: {exc}  (the simulation goes on)")
        key = (where, type(exc).__name__, str(exc))
        if key not in seen_errors:
            seen_errors.add(key)
            traceback.print_exc()

    def safe(where, fn):
        """Wraps an event handler so that an error does not stop the whole GUI."""

        def wrapped(*a, **kw):
            try:
                return fn(*a, **kw)
            except Exception as exc:  # noqa: BLE001 -- shown to the user, the app goes on
                report_error(where, exc)

        return wrapped

    def _apply_readout(text):
        """Updates the readout without rewriting unchanged content."""
        text = str(text)
        if text == state.get("readout_text", ""):
            return
        y_top = out.yview()[0]
        at_bottom = out.yview()[1] >= 0.98
        out.configure(state="normal")
        out.delete("1.0", "end")
        out.insert("1.0", text)
        out.configure(state="disabled")
        state["readout_text"] = text
        if at_bottom:
            out.see("end")
        else:
            out.yview_moveto(y_top)

    def set_readout(text):
        """Updates the data, but keeps the content stable while text is being selected for copying."""
        text = str(text)
        state["readout_pending"] = text
        if state.get("readout_hold"):
            return
        _apply_readout(text)

    def copy_text_widget(widget, only_selection=False, status=None):
        """Copies the selected text or the whole content of a text widget to the clipboard."""
        try:
            no_selection = False
            if only_selection:
                try:
                    content = widget.get("sel.first", "sel.last")
                except tk.TclError:
                    content = widget.get("1.0", "end-1c")
                    no_selection = True
            else:
                content = widget.get("1.0", "end-1c")
            widget.winfo_toplevel().clipboard_clear()
            widget.winfo_toplevel().clipboard_append(content)
            widget.winfo_toplevel().update_idletasks()
            if status is not None:
                if no_selection:
                    status.config(text=f"No text selected; the whole readout was copied ({len(content)} characters).")
                elif content:
                    status.config(text=f"Copied {len(content)} characters to the clipboard.")
        except Exception as exc:  # noqa: BLE001
            if status is not None:
                status.config(text=f"Copying failed: {exc}")
        return "break"

    def select_all_text(widget):
        widget.tag_add("sel", "1.0", "end-1c")
        widget.mark_set("insert", "1.0")
        widget.see("insert")
        return "break"

    def bind_text_copy(widget, status=None):
        """Common copy shortcuts and a context menu."""
        for seq in ("<Control-a>", "<Control-A>"):
            widget.bind(seq, lambda _e, w=widget: select_all_text(w))
        for seq in ("<Control-c>", "<Control-C>", "<Control-Insert>"):
            widget.bind(seq, lambda _e, w=widget, st=status: copy_text_widget(w, True, st))
        menu = tk.Menu(widget, tearoff=0)
        menu.add_command(label="Copy selection", command=lambda w=widget, st=status: copy_text_widget(w, True, st))
        menu.add_command(label="Copy all", command=lambda w=widget, st=status: copy_text_widget(w, False, st))
        menu.add_command(label="Select all", command=lambda w=widget: select_all_text(w))

        def show_menu(event, m=menu):
            try:
                m.tk_popup(event.x_root, event.y_root)
            finally:
                m.grab_release()
            return "break"

        widget.bind("<Button-3>", show_menu)
        widget.bind("<Button-2>", show_menu)

    def _release_readout():
        """After leaving the readout the live data return (not while selecting or copying)."""
        try:
            focus = root.focus_get()
        except (tk.TclError, KeyError):
            focus = None
        if focus is out:
            return
        state["readout_hold"] = False
        if state.get("readout_pending") != state.get("readout_text"):
            _apply_readout(state.get("readout_pending", ""))

    out.bind("<FocusIn>", lambda _e=None: state.update(readout_hold=True))
    out.bind("<FocusOut>", lambda _e=None: root.after(350, _release_readout))  # lets "Copy selection" finish first
    bind_text_copy(out, out_status)
    for txt_, cmd in (("Copy all", lambda: copy_text_widget(out, False, out_status)),
                      ("Copy selection", lambda: copy_text_widget(out, True, out_status)),
                      ("Select all", lambda: select_all_text(out)),
                      ("Refresh", lambda: (state.update(readout_hold=False), readouts()))):
        tk.Button(out_toolbar, text=txt_, command=cmd).pack(side="left", padx=2, pady=2)

    nb = ttk.Notebook(left)
    nb.pack(fill="both", expand=True, padx=6, pady=6)

    def scroll_tab(title):
        """A tab that scrolls vertically even in a small window."""
        wrap = tk.Frame(nb, bg=PANEL)
        nb.add(wrap, text=title)
        area = tk.Canvas(wrap, bg=PANEL, highlightthickness=0, bd=0)
        bar = ttk.Scrollbar(wrap, orient="vertical", command=area.yview)
        content = tk.Frame(area, bg=PANEL)
        win = area.create_window((0, 0), window=content, anchor="nw")
        area.configure(yscrollcommand=bar.set)
        area.pack(side="left", fill="both", expand=True)
        bar.pack(side="right", fill="y")
        content.bind("<Configure>", lambda _e, c=area: c.configure(scrollregion=c.bbox("all")))
        area.bind("<Configure>", lambda e, c=area, w=win: c.itemconfigure(w, width=e.width))

        def wheel(event, c=area, o=wrap):
            widget = event.widget
            while widget is not None:
                if widget is o:
                    delta = getattr(event, "delta", 0)
                    if delta:
                        c.yview_scroll(-int(delta / 120) or (-1 if delta > 0 else 1), "units")
                    elif getattr(event, "num", None) == 4:
                        c.yview_scroll(-1, "units")
                    elif getattr(event, "num", None) == 5:
                        c.yview_scroll(1, "units")
                    return "break"
                widget = getattr(widget, "master", None)
            return None

        for seq in ("<MouseWheel>", "<Button-4>", "<Button-5>"):
            wrap.bind_all(seq, wheel, add="+")
        return wrap, content

    _, tab_c = scroll_tab("Controls")
    _, tab_g = scroll_tab("Comparison")
    _, tab_s = scroll_tab("Theory status")
    _, tab_a = scroll_tab("About")

    def label(parent, text, **kw):
        kw.setdefault("fg", TXT)
        return tk.Label(parent, text=text, bg=PANEL, anchor="w", justify="left", **kw)

    def check(parent, text, var, command=None, **kw):
        return tk.Checkbutton(parent, text=text, variable=var, command=command, bg=PANEL, fg=TXT,
                              selectcolor="#2b2f3a", activebackground=PANEL, anchor="w", justify="left", **kw)

    # ---------- Controls tab
    label(tab_c, "Scene").pack(fill="x", padx=8, pady=(8, 0))
    scene_var = tk.StringVar(value=SCENES[0])
    ttk.Combobox(tab_c, textvariable=scene_var, values=SCENES, state="readonly", height=24).pack(fill="x", padx=8)
    btns = tk.Frame(tab_c, bg=PANEL)
    btns.pack(fill="x", padx=8, pady=8)
    run_btn = tk.Button(btns, text="Run", width=9)
    run_btn.pack(side="left")
    tk.Button(btns, text="Step", width=6, command=lambda: tick(single=True)).pack(side="left", padx=4)
    tk.Button(btns, text="Reset", width=7, command=lambda: load_scene(scene_var.get())).pack(side="left")

    speed_lbl = label(tab_c, "Pace (log₁₀ of simulated time per second of real time)", wraplength=300)
    speed_lbl.pack(fill="x", padx=8)
    speed_var = tk.DoubleVar(value=0.0)
    speed_scale = tk.Scale(tab_c, from_=-3.0, to=12.5, resolution=0.05, orient="horizontal", variable=speed_var,
                           bg=PANEL, fg=TXT, highlightthickness=0)
    speed_scale.pack(fill="x", padx=8)
    speed_val = tk.Label(tab_c, text="", bg=PANEL, fg=ACCENT, anchor="w", justify="left", wraplength=300)
    speed_val.pack(fill="x", padx=8)

    def speed_value_label(*_):
        """The same pace control in every scene, independent of the frame rate."""
        pace = 10.0 ** speed_var.get()
        if state["meta"].get("si"):
            txt = f"1 s of real time = {fmt_time(pace)} simulated"
            if state.get("achieved") is not None:
                txt += f"  (achieved {fmt_time(state['achieved'])}/s)"
        else:
            txt = f"{pace:.4g} model time units per second of real time"
            if state.get("achieved") is not None:
                txt += f"  (achieved {state['achieved']:.4g}/s)"
        speed_val.config(text=txt)

    speed_var.trace_add("write", speed_value_label)

    trails_var, field_var, names_var = tk.BooleanVar(value=True), tk.BooleanVar(value=True), tk.BooleanVar(value=True)
    for txt_, var in (("Trails", trails_var), ("Field U (fixed scale of the scene, legend bottom left)", field_var),
                      ("Labels", names_var)):
        check(tab_c, txt_, var).pack(fill="x", padx=8)
    projection_var = tk.BooleanVar(value=False)
    projection_check = check(tab_c, "3D view of a planar scene (projection only)", projection_var,
                             command=lambda: change_projection(), wraplength=290)
    projection_check.pack(fill="x", padx=8, pady=(4, 0))

    label(tab_c, "Mouse: wheel = zoom, right button = pan").pack(fill="x", padx=8, pady=(10, 0))
    click_lbl = label(tab_c, "", wraplength=300)
    click_lbl.pack(fill="x", padx=8)
    tool_var = tk.StringVar(value="none")
    tool_widgets = {}
    for txt_, val in (("none", "none"), ("add body", "mass"), ("add photon", "photon"), ("add clock", "clock"),
                      ("delete (click on an object)", "delete")):
        rb = tk.Radiobutton(tab_c, text=txt_, value=val, variable=tool_var, bg=PANEL, fg=TXT, selectcolor="#2b2f3a",
                            activebackground=PANEL, anchor="w")
        rb.pack(fill="x", padx=16)
        tool_widgets[val] = rb

    # unified mass input: a number in scientific notation and a unit
    label(tab_c, "Mass of a new body (scientific notation, e.g. 5.972e24)", wraplength=300).pack(fill="x", padx=8,
                                                                                               pady=(6, 0))
    mass_row = tk.Frame(tab_c, bg=PANEL)
    mass_row.pack(fill="x", padx=8)
    mass_var = tk.StringVar(value="1.0000e-06")
    mass_entry = tk.Entry(mass_row, textvariable=mass_var, width=14)
    mass_entry.pack(side="left", fill="x", expand=True)
    unit_var = tk.StringVar(value=list(MASS_UNITS_GEO)[0])
    unit_box = ttk.Combobox(mass_row, textvariable=unit_var, values=list(MASS_UNITS_GEO), state="readonly", width=20)
    unit_box.pack(side="left", padx=(4, 0))
    cat_lbl = tk.Label(tab_c, text="", bg=PANEL, fg=ACCENT, anchor="w", justify="left", wraplength=300)
    cat_lbl.pack(fill="x", padx=8)

    pilot_frame = tk.Frame(tab_c, bg=PANEL)
    pilot_frame.pack(fill="x", padx=8, pady=(6, 0))
    tk.Label(pilot_frame, text="UAP control: arrows/WASD, Shift = ×10, space toggles hovering (1 g). The slider sets "
             "the thrust [log₁₀ g] or the speed in the 'null field' hypothesis [log₁₀(m/s) − 2]:", bg=PANEL, fg=TXT,
             anchor="w", justify="left", wraplength=300).pack(fill="x")
    thrust_var = tk.DoubleVar(value=0.5)
    tk.Scale(pilot_frame, from_=-1.0, to=3.0, resolution=0.05, orient="horizontal", variable=thrust_var, bg=PANEL,
             fg=TXT, highlightthickness=0).pack(fill="x")
    nullf_var, follow_var = tk.BooleanVar(value=False), tk.BooleanVar(value=True)
    check(pilot_frame, "'Null field' hypothesis: velocity control without inertia, crew at 0 g", nullf_var,
          wraplength=290).pack(fill="x")
    check(pilot_frame, "Camera follows the UAP", follow_var).pack(fill="x")
    pilot_frame.pack_forget()

    # form for adding an object to the real 3D star map
    map3d_frame = tk.LabelFrame(tab_c, text="New object in the 3D map", bg=PANEL, fg=TXT, bd=1, padx=6, pady=4)
    map_name_var = tk.StringVar(value="New body")
    map_xyz_vars = [tk.StringVar(value="0") for _ in range(3)]
    map_v_vars = [tk.StringVar(value="0") for _ in range(3)]
    tk.Label(map3d_frame, text="Name", bg=PANEL, fg=TXT, anchor="w").grid(row=0, column=0, columnspan=3, sticky="w")
    tk.Entry(map3d_frame, textvariable=map_name_var).grid(row=1, column=0, columnspan=3, sticky="ew", pady=(0, 4))
    for j, ax in enumerate(("X", "Y", "Z")):
        tk.Label(map3d_frame, text=f"{ax} [ly]", bg=PANEL, fg=TXT).grid(row=2, column=j, sticky="w")
        tk.Entry(map3d_frame, textvariable=map_xyz_vars[j], width=8).grid(row=3, column=j, sticky="ew", padx=(0, 3))
    for j, ax in enumerate(("Vx", "Vy", "Vz")):
        tk.Label(map3d_frame, text=f"{ax} [km/s]", bg=PANEL, fg=TXT).grid(row=4, column=j, sticky="w")
        tk.Entry(map3d_frame, textvariable=map_v_vars[j], width=8).grid(row=5, column=j, sticky="ew", padx=(0, 3))
    for j in range(3):
        map3d_frame.columnconfigure(j, weight=1)
    tk.Button(map3d_frame, text="Add object (mass from the field above)", command=lambda: add_map3d_object()).grid(
        row=6, column=0, columnspan=3, sticky="ew", pady=(5, 0))
    map3d_frame.pack_forget()
    keys = set()

    def mass_units():
        return MASS_UNITS_SI if state["meta"].get("si") else MASS_UNITS_GEO

    def parse_mass():
        """Mass of a new body in geometric units [m] from the entry and the unit; None for invalid input."""
        v = parse_sci(mass_var.get())
        f = mass_units().get(unit_var.get())
        return None if v is None or f is None else v * f

    def reference_mass(sim=None):
        """Mass of the heaviest source when the scene was loaded (geometric units)."""
        sim = sim or state["sim"]
        if sim is state.get("sim") and state.get("mass_reference") is not None:
            return float(state["mass_reference"])
        vals = [b.mass for b in sim.bodies if b.alive and b.mass > 0 and b is not sim.profile]
        return max(vals, default=(M_SUN if state.get("meta", {}).get("map3d") else 1.0))

    def update_mass_info(*_):
        m = parse_mass()
        meta = state.get("meta", {})
        if m is None:
            text = "Invalid mass: enter a positive number in scientific notation, e.g. 5.972e24."
        elif meta.get("map3d"):
            text = f"= {fmt_sci(m / M_SUN)} M☉ = {fmt_sci(m / G_C2)} kg"
        elif meta.get("si"):
            kg = m / G_C2
            text = (f"= {fmt_sci(kg)} kg = {fmt_sci(kg / M_EARTH_KG)} M⊕ = {fmt_sci(kg / M_SUN_KG)} M☉ "
                    f"→ {category_si(kg)}")
        else:
            r = m / reference_mass()
            text = f"= {fmt_sci(m)} model units = {fmt_sci(r)} × heaviest body → {category_relative(r)}"
        cat_lbl.config(text=text)

    def on_unit_change(*_):
        """Changing the unit keeps the physical mass and rewrites the number."""
        prev = state.get("mass_unit_prev")
        units = mass_units()
        new = unit_var.get()
        v = parse_sci(mass_var.get())
        if prev in units and new in units and v is not None and prev != new:
            mass_var.set(fmt_sci(v * units[prev] / units[new]))
        state["mass_unit_prev"] = new
        update_mass_info()

    def normalise_mass(_e=None):
        v = parse_sci(mass_var.get())
        if v is not None:
            mass_var.set(fmt_sci(v))

    mass_var.trace_add("write", update_mass_info)
    unit_var.trace_add("write", on_unit_change)
    mass_entry.bind("<FocusOut>", normalise_mass)
    mass_entry.bind("<Return>", normalise_mass)
    note_lbl = label(tab_c, "", wraplength=300)
    note_lbl.pack(fill="x", padx=8)

    # ---------- Theory status and About tabs
    for i, (k, v) in enumerate(STATUS):
        warn = v.startswith(("INPUT", "HYPOTHESIS", "NOT")) or "input" in v or "not derived" in v or "not tested" in v
        label(tab_s, k, font=("Segoe UI", 9, "bold"), wraplength=300).pack(fill="x", padx=8, pady=(6 if i else 10, 0))
        tk.Label(tab_s, text=v, bg=PANEL, fg="#ffb37a" if warn else "#9ad99a", anchor="w", justify="left",
                 wraplength=290).pack(fill="x", padx=16)
    tk.Label(tab_a, text=ABOUT, bg=PANEL, fg=TXT, anchor="nw", justify="left", font=("Consolas", 9)).pack(
        fill="both", padx=8, pady=8)

    # ---------- Comparison tab: frameworks, PPN sliders, Mercury, classical tests
    label(tab_g, "Frameworks: the same scene integrated in parallel in another framework (dashed trails; the main "
          "scene is causal gravity). Ticking restarts the scene.", wraplength=300).pack(fill="x", padx=8, pady=(8, 4))
    ghost_vars = {}
    for key, txt_ in (("GR", "general relativity (Einstein)"), ("N", "Newton (light as a corpuscle)"),
                      ("H1", "scalar stage (Π in a fixed frame, H1)"), ("F3", "fixed order (F3, γ = −1)"),
                      ("PPN", "PPN with the β, γ sliders")):
        ghost_vars[key] = tk.BooleanVar(value=False)
        check(tab_g, txt_, ghost_vars[key]).pack(fill="x", padx=8)
        ghost_vars[key].trace_add("write", lambda *a: safe("frameworks", load_scene)(scene_var.get()))
    beta_var, gamma_var = tk.DoubleVar(value=1.0), tk.DoubleVar(value=1.0)
    for txt_, var, lo, hi in (("β (nonlinearity, GR = 1)", beta_var, 0.0, 2.0),
                              ("γ (space curvature, GR = 1)", gamma_var, -1.0, 2.0)):
        label(tab_g, txt_).pack(fill="x", padx=8)
        tk.Scale(tab_g, from_=lo, to=hi, resolution=0.05, orient="horizontal", variable=var, bg=PANEL, fg=TXT,
                 highlightthickness=0).pack(fill="x", padx=8)
    tk.Button(tab_g, text="Apply β, γ (the scene restarts)", command=lambda: load_scene(scene_var.get())).pack(
        fill="x", padx=8, pady=4)

    label(tab_g, "\nClassical tests: light bending, probe bending, Mercury, Shapiro, Pound–Rebka, Gravity Probe B, "
          "binary pulsar, moving matter and black-hole shadow, measured in causal gravity, Einstein, Newton, the "
          "scalar stage, F3 and PPN (sliders). Runs in the background (~1–2 min); the table goes to the readout.",
          wraplength=300).pack(fill="x", padx=8)
    tk.Button(tab_g, text="Run classical tests", command=lambda: start_tests()).pack(fill="x", padx=8, pady=4)
    tests_lbl = tk.Label(tab_g, text="", bg=PANEL, fg=ACCENT, anchor="w", justify="left", wraplength=300)
    tests_lbl.pack(fill="x", padx=8)

    label(tab_g, "\nMercury experiment: perihelion advance in ″/century (Newton, causal, Einstein, scalar stage and "
          "PPN from the sliders); runs in the background.", wraplength=300).pack(fill="x", padx=8)
    exp_frame = tk.Frame(tab_g, bg=PANEL)
    exp_frame.pack(fill="x", padx=8)
    tk.Label(exp_frame, text="years:", bg=PANEL, fg=TXT).pack(side="left")
    years_var = tk.StringVar(value="10")
    tk.Entry(exp_frame, textvariable=years_var, width=5).pack(side="left", padx=4)
    planets_var = tk.BooleanVar(value=False)
    tk.Checkbutton(exp_frame, text="with planets", variable=planets_var, bg=PANEL, fg=TXT, selectcolor="#2b2f3a",
                   activebackground=PANEL).pack(side="left")
    tk.Button(tab_g, text="Run Mercury experiment", command=lambda: start_experiment()).pack(fill="x", padx=8, pady=4)
    exp_lbl = tk.Label(tab_g, text="", bg=PANEL, fg=ACCENT, anchor="w", justify="left", wraplength=300)
    exp_lbl.pack(fill="x", padx=8)
    label(tab_g, "\nBridge to the laboratory: test records, the prediction board, field profiles against the "
          "continuum solution.", wraplength=300).pack(fill="x", padx=8)
    tk.Button(tab_g, text="Open the laboratory and the board", command=lambda: safe("laboratory", open_lab)()).pack(
        fill="x", padx=8, pady=4)

    def open_lab():
        """A window with the laboratory records, the prediction board and a viewer of results."""
        win = tk.Toplevel(root)
        win.title("Laboratory: prediction board and records")
        win.geometry("1250x820")
        win.configure(bg=BG)
        top = tk.Frame(win, bg=BG)
        top.pack(fill="x", padx=6, pady=6)
        tk.Label(top, text="results folder:", bg=BG, fg=TXT).pack(side="left")
        path_var = tk.StringVar(value=default_results_dir())
        tk.Entry(top, textvariable=path_var, width=70).pack(side="left", padx=4)
        quick_var = tk.BooleanVar(value=False)
        tk.Checkbutton(top, text="include quick runs", variable=quick_var, bg=BG, fg=TXT, selectcolor="#2b2f3a",
                       activebackground=BG).pack(side="left")
        lnb = ttk.Notebook(win)
        lnb.pack(fill="both", expand=True, padx=6, pady=6)
        f_board, f_rec = tk.Frame(lnb, bg=BG), tk.Frame(lnb, bg=BG)
        lnb.add(f_board, text="Prediction board")
        lnb.add(f_rec, text="Records and field profiles")
        board_tools = tk.Frame(f_board, bg=BG)
        board_tools.pack(side="top", fill="x")
        bsb = tk.Scrollbar(f_board)
        bsb.pack(side="right", fill="y")
        board = tk.Text(f_board, bg="#10131a", fg=TXT, font=("Segoe UI", 10), wrap="word", yscrollcommand=bsb.set)
        bsb.config(command=board.yview)
        board_status = tk.Label(board_tools, text="", bg=BG, fg=DIM, anchor="w")
        board_status.pack(side="right", padx=6)
        for txt_, cmd in (("Copy all", lambda: copy_text_widget(board, False, board_status)),
                          ("Copy selection", lambda: copy_text_widget(board, True, board_status)),
                          ("Select all", lambda: select_all_text(board))):
            tk.Button(board_tools, text=txt_, command=cmd).pack(side="left", padx=3, pady=3)
        board.pack(fill="both", expand=True)
        bind_text_copy(board, board_status)
        for tag, col in (("h", "#ffffff"), ("ok", "#9ad99a"), ("bad", "#ff8a7a"), ("inst", "#ffb37a"),
                         ("none", DIM), ("txt", TXT)):
            board.tag_config(tag, foreground=col)
        board.tag_config("h", font=("Segoe UI", 11, "bold"))
        left_f = tk.Frame(f_rec, bg=BG)
        left_f.pack(side="left", fill="y")
        lb = tk.Listbox(left_f, width=58, bg="#10131a", fg=TXT, font=("Consolas", 9), selectbackground="#2b4a6a")
        lsb = tk.Scrollbar(left_f, command=lb.yview)
        lb.config(yscrollcommand=lsb.set)
        lsb.pack(side="right", fill="y")
        lb.pack(side="left", fill="y")
        right_f = tk.Frame(f_rec, bg=BG)
        right_f.pack(side="left", fill="both", expand=True)
        plot = tk.Canvas(right_f, bg="#0b0d12", height=330, highlightthickness=0)
        plot.pack(fill="x")
        record_tools = tk.Frame(right_f, bg=BG)
        record_tools.pack(side="top", fill="x")
        record_status = tk.Label(record_tools, text="", bg=BG, fg=DIM, anchor="w")
        record_status.pack(side="right", padx=6)
        tsb = tk.Scrollbar(right_f)
        tsb.pack(side="right", fill="y")
        txt = tk.Text(right_f, bg="#10131a", fg=TXT, font=("Consolas", 9), wrap="word", yscrollcommand=tsb.set)
        txt.pack(fill="both", expand=True)
        tsb.config(command=txt.yview)
        for txt_, cmd in (("Copy all", lambda: copy_text_widget(txt, False, record_status)),
                          ("Copy selection", lambda: copy_text_widget(txt, True, record_status)),
                          ("Select all", lambda: select_all_text(txt))):
            tk.Button(record_tools, text=txt_, command=cmd).pack(side="left", padx=3, pady=3)
        bind_text_copy(txt, record_status)
        lab = {"recs": []}

        def status_tag(st):
            return {"PASSED": "ok", "FAILED": "bad", "INSTRUMENT_FAILED": "inst"}.get(st, "none")

        def fill_board():
            board.delete("1.0", "end")
            recs = lab["recs"]
            board.insert("end", f"Prediction board (records: {len(recs)} in {path_var.get()})\nFor every claim: what "
                         "would refute it, the latest records of its tests (automatic) and the reading from "
                         "docs/nalezy.md.\n\n", "txt")
            for claim, tests, falsifier, reading in PREDICTIONS:
                board.insert("end", claim + "\n", "h")
                board.insert("end", f"   refuted by: {falsifier}\n", "txt")
                found = False
                for tid in tests:
                    latest = {}
                    for rec in recs:  # newest record per part
                        if rec.get("id") == tid and not rec.get("quick"):
                            latest.setdefault(rec.get("part") or "", rec)
                    for part, rec in sorted(latest.items()):
                        found = True
                        st = str(rec.get("status", "?"))
                        board.insert("end", f"   {tid} {part or ''}: ", "txt")
                        board.insert("end", st, status_tag(st))
                        board.insert("end", f"  ({rec.get('time_utc', '')})\n", "none")
                if not found:
                    board.insert("end", "   records: none (not measured)\n", "none")
                board.insert("end", f"   reading: {reading}\n\n", "txt")

        def fill_list():
            lb.delete(0, "end")
            for rec in lab["recs"]:
                q = " quick" if rec.get("quick") else ""
                lb.insert("end", f"{str(rec.get('id') or '?'):7s} {str(rec.get('part') or ''):10s} "
                                 f"{str(rec.get('status') or '?'):17s} {rec.get('time_utc') or ''}{q}")

        def load():
            lab["recs"] = load_records(path_var.get(), quick_var.get())
            fill_board()
            fill_list()

        def show(_ev=None):
            sel = lb.curselection()
            if not sel:
                return
            rec = lab["recs"][sel[0]]
            txt.delete("1.0", "end")
            txt.insert("end", f"{rec.get('id')} {rec.get('part') or ''}: {rec.get('status')}  ({rec.get('_path')})\n\n")
            txt.insert("end", "\n".join(str(x) for x in (rec.get("lines") or [])))
            plot.delete("all")
            d, prm = rec.get("data") or {}, rec.get("params") or {}
            if isinstance(d, dict) and d.get("nodes") and d.get("Pi_cs") and int(prm.get("D", 0)) == 3:
                draw_profile(rec)
            elif rec.get("id") == "F1.9a" and isinstance(d, dict) and d.get("table"):
                draw_f19a(rec)
            else:
                plot.create_text(10, 10, anchor="nw", fill=DIM, text="No plot for this record (plots: the field (a) "
                                 "in 2+1 from K.12/F2.7 and the F1.9a table).")

        def axes(xmin, xmax, ymin, ymax, xl, yl, xticks=None):
            w, h = max(plot.winfo_width(), 400), 330
            X = lambda x: 60 + (x - xmin) / (xmax - xmin) * (w - 80)  # noqa: E731
            Y = lambda y: h - 40 - (y - ymin) / (ymax - ymin) * (h - 70)  # noqa: E731
            plot.create_line(60, h - 40, w - 20, h - 40, fill=DIM)
            plot.create_line(60, 30, 60, h - 40, fill=DIM)
            for k in range(6):
                xv = xmin + (xmax - xmin) * k / 5
                yv = ymin + (ymax - ymin) * k / 5
                if xticks is None:
                    plot.create_text(X(xv), h - 28, text=f"{xv:.3g}", fill=DIM, font=("Segoe UI", 8))
                plot.create_text(52, Y(yv), text=f"{yv:.3g}", fill=DIM, anchor="e", font=("Segoe UI", 8))
            for xv, lab_ in xticks or []:
                plot.create_text(X(xv), h - 28, text=lab_, fill=DIM, font=("Segoe UI", 8))
            plot.create_text(w / 2, h - 12, text=xl, fill="#c8ccd6", font=("Segoe UI", 9))
            plot.create_text(10, 14, text=yl, fill="#c8ccd6", anchor="w", font=("Segoe UI", 9))
            return X, Y

        def draw_profile(rec):
            d, prm = rec["data"], rec["params"]
            nodes, pic = np.array(d["nodes"], float), np.array(d["Pi_cs"], float)
            R, sig, q = float(prm["R"]), float(prm["sigma"]), float(prm["q"])
            r, pa = continuum_a_2p1(R, sig, q)
            X, Y = axes(0, R, min(pic.min(), pa.min()) - 0.01, 1.005, "r [ℓ_k]", "Π(r)")
            pts = []
            for rr, pp in zip(r[::50], pa[::50]):
                pts += [X(rr), Y(pp)]
            plot.create_line(*pts, fill="#7ad0ff", width=2)
            pts = []
            for rr, pp in zip(nodes, pic):
                pts += [X(rr), Y(pp)]
                plot.create_oval(X(rr) - 2, Y(pp) - 2, X(rr) + 2, Y(pp) + 2, fill="#ffd27a", outline="")
            if len(pts) >= 4:
                plot.create_line(*pts, fill="#ffd27a")
            w = max(plot.winfo_width(), 400)
            plot.create_text(w - 20, 250, anchor="ne", fill="#ffd27a", text="● field of the causal set (nodes)")
            plot.create_text(w - 20, 266, anchor="ne", fill="#7ad0ff",
                             text=f"— continuum solution (a), same wall, σ = {sig:g}, q = {q:g}, R = {R:g}")

        def draw_f19a(rec):
            tab = rec["data"]["table"]
            rows = [t for t in tab if t.get("gamma") in ("h1", "conf")]
            if not rows:
                return
            vals = [v for t in rows for v in t["c_disc"] + t["c_cont"]]
            ns = len(rows[0]["c_disc"])
            sig = (rec.get("params") or {}).get("sigmas", list(range(ns)))
            X, Y = axes(-0.15, max(ns - 1, 1) + 0.15, min(vals + [-0.5]) - 0.1, max(vals + [2.0]) + 0.1,
                        "σ / ℓ_k (both rows r/σ)", "c = dB_L[L] / |∇L|²",
                        xticks=[(i, f"σ = {sig[i]:g}") for i in range(min(ns, len(sig)))])
            cols = {"h1": "#ffffff", "conf": "#ffb37a"}
            for t in rows:
                pts = []
                for i, v in enumerate(t["c_disc"]):
                    pts += [X(i), Y(v)]
                    plot.create_oval(X(i) - 3, Y(v) - 3, X(i) + 3, Y(v) + 3, fill=cols[t["gamma"]], outline="")
                if len(pts) >= 4:
                    plot.create_line(*pts, fill=cols[t["gamma"]])
                plot.create_line(X(0), Y(t["c_cont"][0]), X(max(ns - 1, 1)), Y(t["c_cont"][0]),
                                 fill=cols[t["gamma"]], dash=(3, 3))
            plot.create_text(70, 40, anchor="nw", fill="#ffffff", text="● H1 (prediction c = 0, dashed: continuum)")
            plot.create_text(70, 56, anchor="nw", fill="#ffb37a", text="● control γ = −1 (c = 2)")

        tk.Button(top, text="Load", command=safe("laboratory", load)).pack(side="left", padx=6)
        lb.bind("<<ListboxSelect>>", safe("laboratory", show))
        safe("laboratory", load)()

    def selected_laws():
        laws = {"GR": _law_gr, "N": _law_newton, "H1": _law_h1, "F3": _law_f3}
        out_ = [laws[k]() for k in ("GR", "N", "H1", "F3") if ghost_vars[k].get()]
        if ghost_vars["PPN"].get():
            out_.append(_law_ppn(beta_var.get(), gamma_var.get()))
        return out_

    def start_experiment():
        try:
            years = max(1.0, min(200.0, float(years_var.get().replace(",", "."))))
        except ValueError:
            years = 10.0
        laws = [_law_newton(), _law_cg(), _law_gr(), _law_h1(), _law_ppn(beta_var.get(), gamma_var.get())]
        state["exp"] = mercury_experiment(laws, years, planets_var.get())
        state["exp_meta"] = (years, planets_var.get(), beta_var.get(), gamma_var.get(), laws)
        state["exp_result"] = None
        exp_lbl.config(text="running …")

    def start_tests():
        laws = [_law_cg(), _law_gr(), _law_newton(), _law_h1(), _law_f3(), _law_ppn(beta_var.get(), gamma_var.get())]
        state["tests"] = classical_tests(laws)
        state["tests_laws"] = laws
        state["tests_result"] = None
        tests_lbl.config(text="running …")

    def background_work(budget):
        """Advances the Mercury experiment and the classical tests within a time budget [s] per frame."""
        t2 = time.perf_counter()
        for key, lbl, done_text in (("exp", exp_lbl, "done, results in the readout below"),
                                    ("tests", tests_lbl, "done, the table is at the end of the readout")):
            if state[key] is None:
                continue
            try:
                while time.perf_counter() - t2 < budget:
                    text, frac = next(state[key])
                    lbl.config(text=f"running: {text}, {100 * frac:.0f} %")
            except StopIteration as stop:
                state[key + "_result"] = stop.value
                state[key] = None
                lbl.config(text=done_text)
                readouts()
            break  # one job per frame keeps the window responsive

    # ---------- coordinates
    def W2S(x, y):
        """Scene to screen, optionally through the tilted 3D projection."""
        cx, cy, hw = state["view"]
        w, h = max(canvas.winfo_width(), 10), max(canvas.winfo_height(), 10)
        sc = 0.5 * w / hw
        if state.get("view3d") and not state["meta"].get("map3d"):
            yaw, pitch = state.get("rot", [0.25, 0.9])
            dx, dy = x - cx, y - cy
            xr = math.cos(yaw) * dx - math.sin(yaw) * dy
            yr = math.sin(yaw) * dx + math.cos(yaw) * dy
            return w / 2 + xr * sc, h / 2 - math.sin(pitch) * yr * sc
        return w / 2 + (x - cx) * sc, h / 2 - (y - cy) * sc

    def S2W(sx, sy):
        """Screen to scene, including the 3D projection of the plane z = 0."""
        cx, cy, hw = state["view"]
        w, h = max(canvas.winfo_width(), 10), max(canvas.winfo_height(), 10)
        sc = 0.5 * w / hw
        xr = (sx - w / 2) / sc
        if state.get("view3d") and not state["meta"].get("map3d"):
            yaw, pitch = state.get("rot", [0.25, 0.9])
            yr = (h / 2 - sy) / (sc * max(abs(math.sin(pitch)), 0.15))
            dx = math.cos(yaw) * xr + math.sin(yaw) * yr
            dy = -math.sin(yaw) * xr + math.cos(yaw) * yr
            return cx + dx, cy + dy
        return cx + xr, cy - (sy - h / 2) / sc

    # ---------- scene management
    def orbits_tracked():
        return (not state["meta"]["si"]) or state["meta"].get("orbits", False)

    def binary_trackers(sim):
        """Radiation trackers of the scenes with radiation (the two heaviest moving bodies)."""
        if not sim.radiation:
            return []
        movers = sorted([b for b in sim.bodies if b.alive and b.kind == "mass" and not b.fixed and b.mass > 0],
                        key=lambda b: -b.mass)
        if len(movers) < 2:
            return []
        return [BinaryTracker(movers[0], movers[1], sim, start_at_periastron=bool(state["meta"].get("pulsar")))]

    def setup_trackers(fresh=False):
        """Measurement trackers of the main scene and the frameworks; the running measurements are kept."""
        sim = state["sim"]
        old_p = {id(t.body): t for t in state["trackers"] if t.body.alive and t.center.alive}
        old_l = {id(t.ph): t for t in state["light"] if t.ph.alive and t.m.alive}
        state["trackers"], state["light"] = make_trackers(sim, orbits_tracked(), old_p, old_l)
        if fresh:
            state["binaries"] = binary_trackers(sim)
        for g in state["ghosts"]:
            op = {id(t.body): t for t in g["trackers"] if t.body.alive and t.center.alive}
            ol = {id(t.ph): t for t in g["light"] if t.ph.alive and t.m.alive}
            g["trackers"], g["light"] = make_trackers(g["sim"], orbits_tracked(), op, ol)
            if fresh:
                g["binaries"] = binary_trackers(g["sim"])

    def configure_sim(sim, meta, bodies):
        """Common settings of the main simulation and of the frameworks."""
        sim.bodies = bodies
        sim.dt_max, sim.eta = meta["dt_max"], meta.get("eta", 0.02)
        sim.radiation = bool(meta.get("radiation"))
        if meta.get("bh"):  # strong field: no guard; the scalar stage uses its exact profile
            c = next(b for b in sim.bodies if getattr(b, "is_profile", False))
            sim.no_guard = True
            sim.capture = (c, 1.0)  # inside the photon sphere of both (isotropic r = 1.866 and 1.798)
            if sim.law.key == "H1":
                sim.profile = c
        for b in sim.bodies:
            b.t_birth, b.start = 0.0, (b.x, b.y, b.vx, b.vy)
        return sim

    def build_ghosts(name, meta):
        """Independent copies of the scene for every selected framework."""
        state["ghosts"] = []
        if meta.get("map3d"):  # the star map has no gravitational dynamics to compare
            return
        for law in selected_laws():
            gs = configure_sim(Sim(law), meta, [b for b in scene(name)[0] if b.kind != "pilot"])
            state["ghosts"].append({"sim": gs, "trackers": [], "light": [], "binaries": [], "law": law})

    def load_scene(name):
        bodies, state["view"], meta = scene(name)
        sim = configure_sim(Sim(MAIN_LAW()), meta, bodies)
        state["sim"], state["scene"], state["meta"] = sim, name, meta
        state["mass_reference"] = max((b.mass for b in sim.bodies if b.alive and b.mass > 0 and b is not sim.profile),
                                      default=(M_SUN if meta.get("map3d") else 1.0))
        tool_var.set("none")
        state["rot"] = [0.6, 0.5]  # rotation of the 3D map: azimuth and elevation [rad]
        build_ghosts(name, meta)
        speed_scale.config(from_=-3.0, to=12.5)
        # mass entry: units and default of the scene
        units = mass_units()
        unit_box.config(values=list(units))
        if meta.get("mass"):
            val, unit = meta["mass"]
        else:
            val, unit = 1e-6 * state["mass_reference"], list(units)[0]
        state["mass_unit_prev"] = None
        unit_var.set(unit)
        state["mass_unit_prev"] = unit
        mass_var.set(fmt_sci(val))
        if meta.get("map3d"):
            note_lbl.config(text="\nReal 3D star map. Positions in ly, velocities in km/s. With the tool 'add body' a "
                            "click puts an object into the plane of view; the form also takes the depth Z and a space "
                            "velocity.")
            click_lbl.config(text="The left button rotates the map. Choose 'add body' or 'delete' to edit objects.")
            map3d_frame.pack(fill="x", padx=8, pady=(6, 0))
            projection_check.pack_forget()
            for key in ("photon", "clock"):
                tool_widgets[key].pack_forget()
            if not tool_widgets["mass"].winfo_manager():
                tool_widgets["mass"].pack(fill="x", padx=16, before=tool_widgets["delete"])
            state["view3d"] = False
            projection_var.set(False)
        else:
            map3d_frame.pack_forget()
            if not projection_check.winfo_manager():
                projection_check.pack(fill="x", padx=8, pady=(4, 0), before=click_lbl)
            for key in ("photon", "clock"):
                if not tool_widgets[key].winfo_manager():
                    tool_widgets[key].pack(fill="x", padx=16, before=tool_widgets["delete"])
            if meta.get("si"):
                note_lbl.config(text="\nSI inputs; the computation uses internal geometric units.")
                click_lbl.config(text=("Click: release at rest; drag: throw (80 px = 20 m/s)." if meta["click"] == "rest"
                                       else "Click: circular orbit; drag sets direction and multiple of the circular "
                                            "speed (80 px = 1×)."))
            else:
                note_lbl.config(text="\nGeometric units G = c = 1 (lengths and masses in the same model unit).")
                click_lbl.config(text=("Click: release at rest; drag: velocity (80 px = 0.01 c)."
                                       if meta["click"] == "rest" else
                                       "Click: circular orbit; drag sets direction and speed (80 px = 1× circular)."))
            state["view3d"] = bool(projection_var.get())
        if meta.get("pilot"):
            pilot_frame.pack(fill="x", padx=8, pady=(6, 0))
        else:
            pilot_frame.pack_forget()
        speed_var.set(meta["speed"])
        out.xview_moveto(0.0)
        out.yview_moveto(0.0)
        state["trackers"], state["light"], state["binaries"], state["n_events"] = [], [], [], 0
        state["achieved"] = None
        speed_value_label()
        update_mass_info()
        err_lbl.config(text="")
        setup_trackers(fresh=True)
        canvas.focus_set()  # the canvas, not the scene list, takes the keys
        draw(force_field=True)

    def change_projection():
        state["view3d"] = bool(projection_var.get()) and not state["meta"].get("map3d")
        state["rot"] = [0.25, 0.95] if state["view3d"] else [0.0, 0.0]
        draw(force_field=True)

    scene_var.trace_add("write", lambda *a: load_scene(scene_var.get()))

    def toggle_run():
        state["running"] = not state["running"]
        run_btn.config(text="Pause" if state["running"] else "Run")

    run_btn.config(command=toggle_run)

    # ---------- image of the field with a fixed colour scale
    def field_image():
        sim = state["sim"]
        if not any(b.alive and (b.mass > 0 or b.kind == "rod") for b in sim.bodies):
            return None
        w, h = canvas.winfo_width(), canvas.winfo_height()
        k = 6
        gw, gh = max(w // k, 2), max(h // k, 2)
        if state.get("view3d") and not state["meta"].get("map3d"):
            sx, sy = np.meshgrid(np.arange(gw) * k + k / 2, np.arange(gh) * k + k / 2)
            cx, cy, hw = state["view"]
            scale = 0.5 * max(w, 10) / hw
            yaw, pitch = state.get("rot", [0.25, 0.95])
            xr = (sx - w / 2) / scale
            yr = (h / 2 - sy) / (scale * max(abs(math.sin(pitch)), 0.15))
            X = cx + math.cos(yaw) * xr + math.sin(yaw) * yr
            Y = cy - math.sin(yaw) * xr + math.cos(yaw) * yr
        else:
            xs = np.array([S2W(i * k + k / 2, h / 2)[0] for i in range(gw)])
            ys = np.array([S2W(w / 2, j * k + k / 2)[1] for j in range(gh)])
            X, Y = np.meshgrid(xs, ys)
        pts = np.stack([X.ravel(), Y.ravel()], axis=1)
        L = sim.L_points(pts)
        disp = state["meta"]["display"]
        if disp[0] == "log":  # fixed logarithmic scale of |U|
            lo, hi = disp[1], disp[2]
            t = np.clip((np.log10(np.maximum(-L, 1e-300)) - lo) / (hi - lo), 0.0, 1.0).reshape(gh, gw)
        else:  # "rel": deviation of L from a reference value
            t = np.clip(1.0 - (L - disp[1]) / disp[2], 0.0, 1.0).reshape(gh, gw)
        # colour scale: dark blue -> violet -> orange -> light yellow; fixed range per scene
        r = np.clip(40 + 230 * t**0.9, 0, 255)
        g = np.clip(18 + 200 * t**2.2, 0, 255)
        bl = np.clip(60 + 120 * np.sin(np.pi * np.minimum(t * 1.2, 1.0)) - 40 * t, 0, 255)
        rgb = np.stack([r, g, bl], axis=2).astype(np.uint8)
        dim = ((t > 0) | (disp[0] == "rel")).astype(np.uint8)[..., None]  # "log": U ~ 0 far away is dark
        rgb = (rgb * dim + np.array([11, 13, 18], np.uint8) * (1 - dim)).astype(np.uint8)
        ppm = b"P6 %d %d 255\n" % (gw, gh) + rgb.tobytes()
        img = tk.PhotoImage(data=ppm, format="PPM")
        return img.zoom(k)

    # ---------- drawing
    def draw_map3d():
        """The real 3D star map, orthographic projection after a rotation (azimuth, elevation)."""
        sim = state["sim"]
        w, h = max(canvas.winfo_width(), 10), max(canvas.winfo_height(), 10)
        cx, cy, hw = state["view"]
        sc = 0.5 * w / hw
        yaw, pitch = state["rot"]
        cyw, syw, cp, sp = math.cos(yaw), math.sin(yaw), math.cos(pitch), math.sin(pitch)

        def proj(x, y, z):
            x1, y1 = cyw * x - syw * y, syw * x + cyw * y  # rotation about the celestial pole
            y2, z2 = cp * y1 - sp * z, sp * y1 + cp * z  # elevation
            return w / 2 + x1 * sc, h / 2 - z2 * sc, y2  # the depth y2 sets the drawing order

        for rr in (4, 8, 12):  # circles in the equatorial plane
            pts = []
            for k in range(73):
                a = 2 * math.pi * k / 72
                X, Y, _ = proj(rr * LY * math.cos(a), rr * LY * math.sin(a), 0.0)
                pts += [X, Y]
            canvas.create_line(*pts, fill="#33415a", dash=(2, 3))
            X, Y, _ = proj(rr * LY, 0, 0)
            canvas.create_text(X + 4, Y, text=f"{rr} ly", fill="#5a6a88", anchor="w", font=("Segoe UI", 8))
        for x, y, z, lab_ in ((13 * LY, 0, 0, "RA 0h"), (0, 13 * LY, 0, "RA 6h"), (0, 0, 13 * LY, "north")):
            X0, Y0, _ = proj(0, 0, 0)
            X1, Y1, _ = proj(x, y, z)
            canvas.create_line(X0, Y0, X1, Y1, fill="#46506a")
            canvas.create_text(X1, Y1, text=lab_, fill="#7a86a0", font=("Segoe UI", 8))
        if trails_var.get():  # paths of the stars in straight space motion
            for b in sim.bodies:
                tr = getattr(b, "trail3", None)
                if tr and len(tr) > 1:
                    pts = []
                    for x, y, z in list(tr)[::2] + [(b.x, b.y, b.z)]:
                        X, Y, _ = proj(x, y, z)
                        pts += [X, Y]
                    canvas.create_line(*pts, fill=b.color, width=1)
        items, projected = [], {}
        for b in sim.bodies:
            X, Y, dep = proj(b.x, b.y, b.z)
            items.append((dep, X, Y, b))
            projected[id(b)] = (X, Y, b)
        state["map3d_points"] = projected
        for dep, X, Y, b in sorted(items, key=lambda t: -t[0]):  # the most distant stars first
            rad = max(2.0, min(9.0, 3.0 + 3.0 * math.log10(max(b.mass / M_SUN, 0.05) / 0.05)))
            X0, Y0, _ = proj(b.x, b.y, 0.0)
            canvas.create_line(X, Y, X0, Y0, fill="#2a3346")  # drop line to the equatorial plane
            canvas.create_oval(X - rad, Y - rad, X + rad, Y + rad, fill=b.color, outline="")
            if names_var.get():
                canvas.create_text(X + rad + 3, Y - 3, text=f"{b.name} ({b.info.split(', ')[-1]})", fill="#c8ccd6",
                                   anchor="w", font=("Segoe UI", 8))
        canvas.create_text(10, 10, anchor="nw", fill=DIM, font=("Segoe UI", 9),
                           text=(f"{state['scene']}   t = {fmt_time(sim.t / C_SI)} from J2000   (tool: "
                                 f"{tool_var.get()}, left button = "
                                 f"{'rotate' if tool_var.get() == 'none' else 'edit object'}, wheel = zoom; the stars "
                                 "move in straight lines)"))
        scale_bar()
        readouts()

    def draw_circle(cx_, cy_, r, **kw):
        x0, y0 = W2S(cx_ - r, cy_ + r)
        x1, y1 = W2S(cx_ + r, cy_ - r)
        canvas.create_oval(x0, y0, x1, y1, **kw)

    def draw_rod(rod):
        """Two counter-streaming beams: dashes moving in opposite directions at a speed that grows with v."""
        xa, ya = W2S(rod.x - rod.half, rod.y)
        xb, yb = W2S(rod.x + rod.half, rod.y)
        moving = rod.vb > 0
        for sgn in (1, -1):
            off = int(state["frame"] * rod.vb * 3.0 * sgn) % 14 if moving else 0
            canvas.create_line(xa, ya - 3 * sgn, xb, yb - 3 * sgn, fill=rod.color, width=2,
                               dash=(8, 6) if moving else (), dashoffset=off)
        if names_var.get():
            canvas.create_text(xb + 8, yb, text=rod.name, fill="#c8ccd6", anchor="w", font=("Segoe UI", 8))

    def draw_spin(b):
        """Spin axis of a rotating body: a dashed axis if it lies in the plane, a circular arrow if it points at us."""
        Jx, Jy, Jz = b.J
        r = max(b.radius, 0.08 * state["view"][2])
        if math.hypot(Jx, Jy) > 0:
            n = math.hypot(Jx, Jy)
            x0, y0 = W2S(b.x - 1.25 * r * Jx / n, b.y - 1.25 * r * Jy / n)
            x1, y1 = W2S(b.x + 1.25 * r * Jx / n, b.y + 1.25 * r * Jy / n)
            canvas.create_line(x0, y0, x1, y1, fill="#9ad9ff", dash=(5, 3), arrow="last")
        if Jz != 0:
            sx, sy = W2S(b.x, b.y)
            rr = 16
            canvas.create_arc(sx - rr, sy - rr, sx + rr, sy + rr, start=30, extent=300 if Jz > 0 else -300,
                              style="arc", outline="#9ad9ff", width=2)

    def draw_gyro(b):
        """A gyroscope: its spin arrow in the plane (rotated by the measured in-plane precession)."""
        g = b.gyro
        sx, sy = W2S(b.x, b.y)
        ang = math.atan2(g["s0"][1], g["s0"][0]) + g["geo"] + g["fd_z"]
        canvas.create_line(sx, sy, sx + 22 * math.cos(ang), sy - 22 * math.sin(ang), fill="#ffffff", width=2,
                           arrow="last")

    def draw_bh_overlay(sim):
        """Horizon, photon sphere and critical impact parameter of causal gravity (= Schwarzschild); with the
        framework H1 also its photon sphere and throat."""
        c = next((b for b in sim.bodies if getattr(b, "is_profile", False) and b.alive), None)
        if c is None:
            return
        draw_circle(c.x, c.y, 0.5 * c.mass, outline="#ffffff", width=2)  # horizon (isotropic r = M/2)
        draw_circle(c.x, c.y, iso_from_areal(3.0, c.mass), outline="#ffffff", dash=(4, 3))  # photon sphere R = 3M
        for bb in (3 * math.sqrt(3) * c.mass, -3 * math.sqrt(3) * c.mass):
            xa, ya = W2S(-40, bb)
            xb, yb = W2S(-6, bb)
            canvas.create_line(xa, ya, xb, yb, fill="#ffffff", dash=(1, 4))
        h1 = next((g for g in state["ghosts"] if g["law"].key == "H1"), None)
        if h1 is not None:
            info = bh_profile()[3]
            col = h1["law"].color
            draw_circle(c.x, c.y, info["r_ph"] * c.mass, outline=col, dash=(2, 3))
            draw_circle(c.x, c.y, info["r_throat"] * c.mass, outline=col, dash=(1, 2))
            xa, ya = W2S(-40, info["b_crit"])
            xb, yb = W2S(-6, info["b_crit"])
            canvas.create_line(xa, ya, xb, yb, fill=col, dash=(1, 4))

    def draw(force_field=False, refresh_field=True):
        canvas.delete("all")
        sim = state["sim"]
        if state["meta"].get("map3d"):
            draw_map3d()
            return
        if field_var.get() and (force_field or (refresh_field and state["frame"] % 4 == 0)
                                or state["field_img"] is None):
            try:
                state["field_img"] = field_image()
            except tk.TclError:
                state["field_img"] = None
        if field_var.get() and state["field_img"] is not None:
            canvas.create_image(0, 0, image=state["field_img"], anchor="nw")
        for b in sim.bodies:  # solids (e.g. the Earth) first
            if b.alive and b.radius > 0:
                draw_solid(b)
        for b in sim.bodies:
            if b.alive and b.kind == "rod":
                draw_rod(b)
            if b.alive and any(b.J):
                draw_spin(b)
        if state["meta"].get("bh"):
            draw_bh_overlay(sim)
        for g in state["ghosts"]:  # frameworks: thin dashed trails, hollow markers
            col = g["law"].color
            for b in g["sim"].bodies:
                if not b.alive or b.radius > 0 or b.fixed or b.kind == "rod":
                    continue
                if trails_var.get() and len(b.trail) > 1:
                    pts = []
                    for x, y in list(b.trail)[::2]:
                        pts.extend(W2S(x, y))
                    if len(pts) >= 4:
                        canvas.create_line(*pts, fill=col, width=1, dash=(4, 3))
                sx, sy = W2S(b.x, b.y)
                canvas.create_oval(sx - 4, sy - 4, sx + 4, sy + 4, outline=col, width=1)
        for b in sim.bodies:
            if not b.alive or b.radius > 0 or b.kind == "rod":
                continue
            if trails_var.get() and len(b.trail) > 1:
                pts = []
                for x, y in list(b.trail)[::2]:
                    pts.extend(W2S(x, y))
                if len(pts) >= 4:
                    canvas.create_line(*pts, fill=b.color, width=1)
            sx, sy = W2S(b.x, b.y)
            if b.kind == "mass" and b.probe:
                canvas.create_polygon(sx, sy - 5, sx + 5, sy, sx, sy + 5, sx - 5, sy, fill=b.color, outline="")
            elif b.kind == "mass":
                rad = max(2.5, min(12.0, 3.0 + 1.1 * math.log10(max(b.mass, 1e-9) / 1e-9)))
                canvas.create_oval(sx - rad, sy - rad, sx + rad, sy + rad, fill=b.color, outline="")
            elif b.kind == "photon":
                canvas.create_oval(sx - 2, sy - 2, sx + 2, sy + 2, fill=b.color, outline="")
            elif b.kind in ("craft", "pilot"):  # a saucer-shaped marker
                canvas.create_oval(sx - 12, sy - 4, sx + 12, sy + 4, fill=b.color, outline="#ffffff")
                canvas.create_oval(sx - 5, sy - 8, sx + 5, sy - 1, fill="#d0ffff", outline="")
            else:
                canvas.create_rectangle(sx - 4, sy - 4, sx + 4, sy + 4, outline=b.color, width=2)
            if b.gyro is not None:
                draw_gyro(b)
            if names_var.get():
                canvas.create_text(sx + 8, sy - 8, text=b.name, fill="#c8ccd6", anchor="w", font=("Segoe UI", 8))
        if state["drag"] is not None:
            (x0, y0), (x1, y1) = state["drag"]
            canvas.create_line(x0, y0, x1, y1, fill="#ffffff", arrow="last", dash=(3, 2))
            if tool_var.get() == "mass" and math.hypot(x1 - x0, y1 - y0) >= 6:
                k = math.hypot(x1 - x0, y1 - y0) / 80.0
                if state["meta"]["click"] == "rest":
                    txt = f"v = {20 * k:.1f} m/s" if state["meta"]["si"] else f"v = {0.01 * k:.4f} c"
                else:
                    txt = f"v = {k:.2f} × circular speed"
                canvas.create_text(x1 + 10, y1, anchor="w", fill="#ffffff", font=("Segoe UI", 9), text=txt)
        if state["meta"]["si"]:
            title = f"{state['scene']}   t = {fmt_time(sim.t / C_SI)}   (SI)"
        else:
            title = f"{state['scene']}   t = {sim.t:.2f}   (G = c = 1)"
        canvas.create_text(10, 10, anchor="nw", fill=DIM, font=("Segoe UI", 9), text=title)
        scale_bar()
        legend()
        w = max(canvas.winfo_width(), 10)
        y = 10
        canvas.create_text(w - 10, y, anchor="ne", fill="#ffffff", font=("Segoe UI", 9),
                           text="— causal gravity (solid trails)" + ("; frameworks:" if state["ghosts"] else ""))
        for g in state["ghosts"]:
            y += 15
            canvas.create_text(w - 10, y, anchor="ne", fill=g["law"].color, font=("Segoe UI", 9),
                               text=f"- - {g['law'].name}")
        if state["hover"] is not None:
            hover_text()
        readouts()

    def hover_text():
        """Values under the cursor: the potential U and the rate of a static clock (per day in SI)."""
        sim = state["sim"]
        hx, hy = state["hover"]
        wx, wy = S2W(hx, hy)
        U = -sim.L_at(wx, wy)
        rate_m1 = float(sim.law.static_rate_m1(min(U, 20.0)))
        h = max(canvas.winfo_height(), 10)
        if state["meta"]["si"]:
            src = [b for b in sim.bodies if b.alive and b.mass > 0]
            c = max(src, key=lambda b: b.mass) if src else None
            dist = f"{fmt_len(math.hypot(wx - c.x, wy - c.y))} from {c.name} | " if c else ""
            txt = f"{dist}U = {U:.4e} | static clock: {rate_m1 * 86400e6:+.3f} μs/day against infinity"
        else:
            txt = f"x = {wx:.3f}, y = {wy:.3f} | U = {U:.4e} | static clock rate dτ/dt = {1 + rate_m1:.6f}"
        canvas.create_text(10, h - 28, anchor="w", text=txt, fill="#ffffff", font=("Segoe UI", 9))

    def draw_solid(b):
        """A solid as a disc if it fits the view, otherwise the visible part of its surface."""
        cx, cy, hw = state["view"]
        w = max(canvas.winfo_width(), 10)
        sc = 0.5 * w / hw
        if b.radius * sc < 5000:
            draw_circle(b.x, b.y, b.radius, fill=b.color, outline="#7fb07a")
        else:
            phi0 = math.atan2(cy - b.y, cx - b.x)  # direction from the centre to the observer
            dphi = 4.0 * hw / b.radius
            pts = []
            for k in range(81):
                ph = phi0 - dphi + 2 * dphi * k / 80
                pts.extend(W2S(b.x + b.radius * math.cos(ph), b.y + b.radius * math.sin(ph)))
            deep = b.radius - 10.0 * hw
            for ph in (phi0 + dphi, phi0 - dphi):
                pts.extend(W2S(b.x + deep * math.cos(ph), b.y + deep * math.sin(ph)))
            canvas.create_polygon(*pts, fill=b.color, outline="#7fb07a")
        if names_var.get() and b.radius * sc < 5000:
            sx, sy = W2S(b.x, b.y)
            canvas.create_text(sx, sy, text=b.name, fill="#e0f0e0", font=("Segoe UI", 9))

    def scale_bar():
        """A scale bar of 1, 2 or 5 × 10^k, about a fifth of the view."""
        cx, cy, hw = state["view"]
        w, h = max(canvas.winfo_width(), 10), max(canvas.winfo_height(), 10)
        sc = 0.5 * w / hw
        target = 0.4 * hw
        unit = 1.0
        if state["meta"]["si"]:
            unit = LY if target >= 0.1 * LY else (AU if target >= 0.01 * AU else 1.0)
        p10 = 10 ** math.floor(math.log10(target / unit))
        length = max(m * p10 for m in (1, 2, 5) if m * p10 <= target / unit) * unit
        px = length * sc
        x0, y0 = w - 30 - px, h - 22
        canvas.create_line(x0, y0, x0 + px, y0, fill="#ffffff", width=2)
        for xx in (x0, x0 + px):
            canvas.create_line(xx, y0 - 5, xx, y0 + 5, fill="#ffffff", width=2)
        if state["meta"]["si"]:
            txt = (f"{length / LY:g} ly" if unit == LY else (f"{length / AU:g} AU" if unit == AU else
                                                              (f"{length / 1000:g} km" if length >= 1000
                                                               else f"{length:g} m")))
        else:
            txt = f"{length:g} length units"
        canvas.create_text(x0 + px / 2, y0 - 10, text=txt, fill="#ffffff", font=("Segoe UI", 9))

    def legend():
        if not field_var.get():
            return
        h = max(canvas.winfo_height(), 10)
        disp = state["meta"]["display"]
        if disp[0] == "log":
            txt = f"colour: log10 U from {disp[1]:g} (dark) to {disp[2]:g} (light), U = ΣM/r (G = c = 1)"
        else:
            txt = (f"colour: clock slowdown against height; light = surface, dark = faster by {disp[2]:.2e} "
                   f"(ΔU = gh/c² over {disp[2] * C_SI**2 / G_SURF:.0f} m)")
        canvas.create_text(10, h - 12, anchor="w", text=txt, fill="#c8ccd6", font=("Segoe UI", 8))

    # ---------- readouts
    def readouts_si():
        """Earth scenes: falls, throws, clock differences and the scripted UAP manoeuvre."""
        sim = state["sim"]
        earth = next((b for b in sim.bodies if b.alive and b.radius > 0), None)
        lines = [f"Earth: M = GM/c² = 4.435 mm, R = 6371 km, g = {G_SURF:.3f} m/s². Weak field (U ≈ 7·10⁻¹⁰): causal "
                 "gravity = GR, Newton suffices on the surface; the simulator shows the same as known physics."]
        if earth is None:
            return lines
        for b in sim.bodies:
            if not b.alive or b.kind != "mass" or b.radius > 0 or not hasattr(b, "start"):
                continue
            x0, y0, vx0, vy0 = b.start
            r0 = math.hypot(x0 - earth.x, y0 - earth.y)
            if r0 <= 0.0:
                continue
            h0 = max(r0 - earth.radius, 0.0)
            g = earth.mass / r0**2 * C_SI**2
            ux, uy = (x0 - earth.x) / r0, (y0 - earth.y) / r0  # local vertical
            vv, vh = (vx0 * ux + vy0 * uy) * C_SI, (vx0 * uy - vy0 * ux) * C_SI
            if h0 > 2e5:  # an orbiting body: no free-fall readout
                continue
            t_pred = (vv + math.sqrt(vv * vv + 2 * g * h0)) / g  # throw in a uniform field from height h0
            if b.landed_t is not None:
                t_meas = (b.landed_t - (b.t_birth or 0.0)) / C_SI
                dist = math.hypot(b.x - x0, b.y - y0)
                lines.append(f"LANDED {b.name}: from {h0:.1f} m in {t_meas:.3f} s (Newton {t_pred:.3f} s), distance "
                             f"{dist:.2f} m (Newton {abs(vh) * t_pred:.2f} m)")
            else:
                lines.append(f"FLIGHT {b.name}: height {math.hypot(b.x - earth.x, b.y - earth.y) - earth.radius:.2f} "
                             f"m, speed {math.hypot(b.vx, b.vy) * C_SI:.2f} m/s")
        ground = next((b for b in sim.bodies if b.alive and b.kind == "clock" and b.ground), None)
        if ground is None:  # without a ground clock every clock is compared with a clock at infinity
            lines += clock_lines()
        else:
            Lg = sim.L_at(ground.x, ground.y)
            age_g = sim.t - (ground.t_birth or 0.0)
            def_g = ground.tau_def / age_g if age_g > 0 else float("nan")  # 1 minus the ground clock rate
            for b in sim.bodies:
                if not b.alive or b is ground or b.kind not in ("clock", "mass") or b.radius > 0:
                    continue
                if b.kind == "mass" and math.hypot(b.x - earth.x, b.y - earth.y) - earth.radius < 2e5:
                    continue
                age = sim.t - (b.t_birth or 0.0)
                if age <= 0:
                    continue
                d_meas = def_g - b.tau_def / age  # rate(b) - rate(ground)
                d_pred = (sim.L_at(b.x, b.y, skip=b) - Lg) - 0.5 * (b.vx * b.vx + b.vy * b.vy)  # 1PN: -ΔU - v²/2
                f, unit = (86400e6, "μs/day") if abs(d_pred) * 86400e6 >= 0.1 else (86400e9, "ns/day")
                lines.append(f"CLOCK {b.name} against the ground clock: {d_meas * f:+.4f} {unit} (prediction −ΔU − v²/2: "
                             f"{d_pred * f:+.4f} {unit}; relative {d_meas:+.3e})")
        for b in sim.bodies:
            if b.alive and b.kind == "craft" and b.script is not None:
                ph = b.script.phase_name(sim.t / C_SI, getattr(b, "phase_names", []))
                lines.append(f"UAP: {ph}; speed {math.hypot(b.vx, b.vy) * C_SI:.0f} m/s; proper acceleration now "
                             f"{b.now_g / 9.81:.1f} g, at most {b.max_g / 9.81:.0f} g (what a thrust would have to "
                             "supply; hovering = 1 g, free fall = 0 g)")
                lines.append("   The author's 'null field' hypothesis (NOT part of the theory or the lab): the craft "
                             "would change the field around itself and fall freely in it, the crew would feel 0 g. "
                             f"Needed gradient |∇U| = a/c² ≈ {b.max_g / C_SI ** 2:.1e} 1/m. The field equation knows "
                             "only positive, attracting sources; a source steered by the craft and momentum "
                             "conservation without reaction mass are not solved in it.")
        return lines

    def clock_lines():
        """Every clock: measured rate against a clock at infinity and the static prediction of the framework."""
        sim = state["sim"]
        si = state["meta"]["si"]
        lines = []
        for b in sim.bodies:
            if not (b.alive and b.kind == "clock"):
                continue
            age = sim.t - (b.t_birth or 0.0)
            if age <= 0:
                lines.append(f"CLOCK {b.name}: waiting for time to pass …")
                continue
            d = b.tau_def / age  # measured 1 - dtau/dt
            pred = -float(sim.law.static_rate_m1(min(-sim.L_at(b.x, b.y), 20.0)))  # static 1 - sqrt(A)
            if si:
                f, unit = (86400e6, "μs/day") if abs(pred) * 86400e6 >= 0.1 else (86400e9, "ns/day")
                lines.append(f"CLOCK {b.name}: against infinity {-d * f:+.4f} {unit} (static prediction "
                             f"{-pred * f:+.4f} {unit}; relative {-d:+.3e})")
            else:
                lines.append(f"CLOCK {b.name}: dτ/dt = {1 - d:.6f} (static prediction √A = {1 - pred:.6f})")
        return lines

    def experiment_lines():
        """Results of the Mercury experiment against the observation and the PPN formula."""
        r = state.get("exp_result")
        if not r:
            return []
        years, planets, beta, gamma, laws = state["exp_meta"]
        names = {law.key: law.name for law in laws}
        lines = [f"MERCURY EXPERIMENT ({years:g} years, {'with planets' if planets else 'Sun + Mercury only'}), "
                 "advance ″/century:"]
        base = r.get("N")
        for k in ("N", "CG", "GR", "H1", "PPN"):
            if r.get(k) is None:
                lines.append(f"   {names.get(k, k)}: too few perihelia")
                continue
            extra = "" if k == "N" or base is None else f"   (difference from Newton {r[k] - base:+.2f}″)"
            lines.append(f"   {names[k]}: {r[k]:+.2f}″{extra}")
        lines.append(f"   observed relativistic part 42.98 ± 0.04″/century (GR 42.98″); PPN (2+2γ−β)/3 · 42.98 = "
                     f"{(2 + 2 * gamma - beta) / 3 * 42.98:.2f}″ for β = {beta:g}, γ = {gamma:g}")
        if planets:
            lines.append("   the Newtonian part from the planets differs in this 2D system from the real ~532″ (no "
                         "inclinations); only the difference from Newton is compared")
        return lines

    def tests_lines():
        if state.get("tests") is not None:
            return ["", "CLASSICAL TESTS running in the background … (see the Comparison tab)"]
        if state.get("tests_result"):
            return [""] + classical_tests_table(state["tests_result"], state["tests_laws"])
        return []

    def twin(g, b):
        return next((x for x in g["sim"].bodies if x.name == b.name), None)

    def gyro_rates(b):
        """(geodetic, in-plane frame dragging, tilt frame dragging, whole orbits) rates [rad per time unit] of a
        gyroscope; averaged over the whole orbits completed so far (running average before the first one)."""
        g = b.gyro
        if g is None or g["t"] <= 0:
            return None
        if g["snaps"]:
            t, geo, fd, fd_z = g["snaps"][-1]
            return geo / t, fd_z / t, fd / t, len(g["snaps"])
        return g["geo"] / g["t"], g["fd_z"] / g["t"], g["fd"] / g["t"], 0

    def comparison_lines():
        """The same measured quantity in the main scene and in every framework."""
        if not state["ghosts"]:
            return []
        sim = state["sim"]
        si = state["meta"]["si"]
        lines = []
        for tr in [t for t in state["trackers"] if t.body.alive]:
            vals = [("Causal", tr)] + [(g["law"].short, next((t for t in g["trackers"] if t.body.name == tr.body.name),
                                                              None)) for g in state["ghosts"]]
            parts = []
            for key, t in vals:
                if t is None:
                    continue
                if si:
                    rr = t.rate()
                    parts.append(f"{key} {rr * 100 * YEAR * C_SI / ARCSEC:+.2f}″/cy" if rr is not None else f"{key} …")
                else:
                    res = t.result()
                    parts.append(f"{key} {math.degrees(res[0]):+.3f}°/orbit" if res else f"{key} …")
            lines.append(f"FRAMEWORKS precession {tr.body.name}: " + " | ".join(parts))
        for lt in [t for t in state["light"] if t.ph.alive]:
            st = lt.state(sim)
            def fmt(s_, body):
                if body.captured:
                    return "captured"
                return fmt_angle(s_["defl"]) + (f" ({s_['defl'] / s_['newton']:.3f}× N)"
                                                 if not s_["photon"] and s_.get("newton") else "")

            parts = [f"Causal {fmt(st, lt.ph)}"]
            for g in state["ghosts"]:
                t = next((x for x in g["light"] if x.ph.name == lt.ph.name), None)
                if t is not None:
                    parts.append(f"{g['law'].short} {fmt(t.state(g['sim']), t.ph)}")
            lines.append(f"FRAMEWORKS bending {lt.ph.name}: " + " | ".join(parts))
        for b in [x for x in sim.bodies if x.alive and x.kind == "clock"]:
            parts = []
            for key, sm, bb in [("Causal", sim, b)] + [(g["law"].short, g["sim"], twin(g, b)) for g in state["ghosts"]]:
                if bb is None:
                    continue
                age = sm.t - (bb.t_birth or 0.0)
                if age <= 0:
                    continue
                d = bb.tau_def / age
                parts.append(f"{key} {-d * 86400e6:+.4f} μs/day" if si else f"{key} {1 - d:.6f}")
            if parts:
                lines.append(f"FRAMEWORKS clock {b.name} ({'against infinity' if si else 'dτ/dt'}): " + " | ".join(parts)
                             + ("  (Newton: clocks do not feel the field)" if any(g["law"].key == "N"
                                                                                  for g in state["ghosts"]) else ""))
        for b in [x for x in sim.bodies if x.alive and x.kind == "mass" and x.landed_t is not None
                  and getattr(x, "rod", None) is None]:
            parts = [f"Causal {(b.landed_t - (b.t_birth or 0.0)) / (C_SI if si else 1.0):.6f}"]
            for g in state["ghosts"]:
                bb = twin(g, b)
                if bb is not None and bb.landed_t is not None:
                    parts.append(f"{g['law'].short} {(bb.landed_t - (bb.t_birth or 0.0)) / (C_SI if si else 1.0):.6f}")
            lines.append(f"FRAMEWORKS landing {b.name} [{'s' if si else 'time'}]: " + " | ".join(parts))
        for b in [x for x in sim.bodies if x.alive and x.gyro is not None]:
            parts = []
            for key, bb in [("Causal", b)] + [(g["law"].short, twin(g, b)) for g in state["ghosts"]]:
                rates = gyro_rates(bb) if bb is not None else None
                if rates is None:
                    continue
                if si:
                    k = YEAR * C_SI / MAS
                    parts.append(f"{key} geo {rates[0] * k:.1f}, drag {abs(rates[2]) * k:.2f} mas/yr")
                else:
                    parts.append(f"{key} geo {rates[0]:.3e}, drag {rates[1]:+.3e} rad/time")
            lines.append(f"FRAMEWORKS gyroscope {b.name}: " + " | ".join(parts))
        if state["binaries"]:
            parts = []
            for key, trs in [("Causal", state["binaries"])] + [(g["law"].short, g["binaries"]) for g in state["ghosts"]]:
                res = trs[0].result() if trs else None
                if state["meta"].get("inspiral"):
                    sep = trs[0].sep() if trs and trs[0].b1.alive and trs[0].b2.alive else None
                    parts.append(f"{key} separation {sep:.3f}" if sep is not None else f"{key} merged")
                else:
                    parts.append(f"{key} dP/dt {res[2]:.4e}" if res else f"{key} …")
            lines.append("FRAMEWORKS radiation: " + " | ".join(parts))
        for b in [x for x in sim.bodies if x.alive and getattr(x, "rod", None) is not None]:
            parts = []
            for key, sm in [("Causal", sim)] + [(g["law"].short, g["sim"]) for g in state["ghosts"]]:
                bb = b if sm is sim else next((x for x in sm.bodies if x.name == b.name), None)
                ref = next((x for x in sm.bodies if getattr(x, "rod", None) is not None and x.rod.vb == 0.0
                            and not x.rod.light), None)
                if bb is not None and getattr(bb, "a_meas", None) is None and sm.t - (bb.t_birth or 0.0) > 400:
                    parts.append(f"{key} 0 (no fall)")
                    continue
                if bb is None or ref is None or getattr(bb, "a_meas", None) is None or \
                        getattr(ref, "a_meas", None) is None:
                    parts.append(f"{key} …")
                    continue
                parts.append(f"{key} {bb.a_meas / ref.a_meas:.3f}")
            lines.append(f"FRAMEWORKS pull of {b.rod.name} ÷ beams at rest: " + " | ".join(parts))
        return lines

    def readouts_solar():
        """Real solar system: perihelion advance of the planets in arcseconds per century."""
        lines = ["Solar system J2000 (2D, inclinations neglected), real masses; weak field (U ≈ 10⁻⁸): causal gravity "
                 "= GR, it differs from Newton by the 1PN effects. Advance from the Laplace–Runge–Lenz vector at the "
                 "perihelia."]
        for tr in [t for t in state["trackers"] if t.body.alive]:
            rr = tr.rate()
            M = tr.center.mass + tr.body.mass
            a = 0.5 * (tr.r_min + tr.r_max) if tr.r_max > 0 else 0.0
            e = (tr.r_max - tr.r_min) / (tr.r_max + tr.r_min) if tr.r_max > 0 else 0.0
            if rr is None or a <= 0:
                lines.append(f"PRECESSION {tr.body.name}: waiting for three perihelion passages …")
                continue
            period = 2 * math.pi * math.sqrt(a**3 / M)
            gr = 6 * math.pi * M / (a * (1 - e * e)) * (100 * YEAR * C_SI / period) / ARCSEC
            lines.append(f"PRECESSION {tr.body.name}: {rr * 100 * YEAR * C_SI / ARCSEC:+.2f}″/century in total "
                         f"({len(tr.lrl)} perihelia; relativity alone would give {gr:.2f}″, the rest are planetary "
                         "perturbations)")
        return lines + clock_lines()

    def readouts_bh():
        """Black hole: causal gravity (= Schwarzschild) and, with the framework H1, the scalar stage."""
        sim = state["sim"]
        c = next((b for b in sim.bodies if getattr(b, "is_profile", False) and b.alive), None)
        if c is None:
            return ["BLACK HOLE DELETED: its field is gone; photons and particles move in the field of the remaining "
                    "bodies (Reset restores the scene)."]
        info = bh_profile()[3]
        lines = [
            "BLACK HOLE (M = 1, lengths in GM/c²). Causal gravity = general relativity in the strong field (decision "
            "13): Schwarzschild metric. Not yet derived on the substrate (Jacobson's route).",
            "   horizon: R = 2M (white circle, isotropic r = M/2) | photon sphere: R = 3M (white dashed) | critical "
            f"impact parameter (shadow): b = 3√3 M = {3 * math.sqrt(3):.4f} M (dotted lines)",
        ]
        h1 = next((g for g in state["ghosts"] if g["law"].key == "H1"), None)
        if h1 is not None:
            lines.append(f"   scalar stage H1 (NOT a prediction, decision 13): no horizon, throat R = "
                         f"{info['R_throat']:.3f} M, photon sphere R = {info['R_ph']:.3f} M, b = {info['b_crit']:.4f} M "
                         f"(+{100 * (info['b_crit'] / (3 * math.sqrt(3)) - 1):.2f} %)")
        for b in sim.bodies:
            if b.kind != "photon":
                continue
            st = "captured" if b.captured else ("absorbed" if not b.alive else "flying / escaped")
            lt = next((t for t in state["light"] if t.ph is b), None)
            dfl = f", bending {math.degrees(lt.state(sim)['defl']):.1f}°" if lt and not b.captured else ""
            line = f"PHOTON {b.name}: causal {st}{dfl}"
            for g in state["ghosts"]:
                gb = twin(g, b)
                if gb is not None:
                    line += f" | {g['law'].short} {'captured' if gb.captured else 'flying / escaped'}"
            lines.append(line)
        for b in sim.bodies:
            if b.kind != "clock":
                continue
            r_iso = math.hypot(b.x - c.x, b.y - c.y)
            R_ar = r_iso * (1 + c.mass / (2 * r_iso)) ** 2
            age = sim.t - (b.t_birth or 0.0)
            meas = 1 - b.tau_def / age if age > 0 else float("nan")
            pred = math.sqrt(max(1 - 2 * c.mass / R_ar, 0.0))
            z = 1 / pred - 1 if pred > 0 else float("inf")
            lines.append(f"CLOCK {b.name}: rate {meas:.4f} (Schwarzschild √(1−2M/R) = {pred:.4f}); redshift of light "
                         f"from here z = {z:.3f}")
        if h1 is None:
            lines.append("Tip: tick the scalar stage H1 in the Comparison tab: photons b = 5.21 and 5.23 escape here and "
                         "are captured in the scalar stage (its shadow is 0.91 % larger).")
        return lines

    def readouts_pilot():
        sim = state["sim"]
        earth = next((b for b in sim.bodies if b.radius > 0), None)
        lines = ["UAP CONTROL: arrows/WASD thrust, Shift ×10, space hovering. Real Earth gravity (SI); thrust = the "
                 "proper acceleration the craft would have to supply."]
        for b in sim.bodies:
            if b.kind != "pilot":
                continue
            h_ = math.hypot(b.x - earth.x, b.y - earth.y) - earth.radius if earth else 0.0
            v = math.hypot(b.vx, b.vy) * C_SI
            g_now, g_max = b.now_g / 9.81, b.max_g / 9.81
            limit = ("safe" if g_max < 9 else "only for pilots in g-suits (≤ 9 g)" if g_max < 15 else
                     "briefly survivable (≲ 45 g)" if g_max < 45 else "LETHAL for the crew (> 45 g)" if g_max < 100
                     else "LETHAL, even the machine would struggle (> 100 g)")
            lines.append(f"UAP: height {fmt_len(h_)}, speed {v:.1f} m/s ({v * 3.6:.0f} km/h), hovering "
                         f"{'ON' if getattr(b, 'hover', False) else 'off'}")
            if nullf_var.get():
                lines.append(f"   null-field mode (HYPOTHESIS): the keys set a speed of {getattr(b, 'v_cmd', 0.0):.0f} "
                             f"m/s, reached within {TAU_NULL / C_SI:.2f} s whatever the gravity (instant stops and "
                             "turns); the crew feels 0 g.")
                lines.append(f"   a conventional craft would need {g_now:.0f} g of thrust now, at most {g_max:.0f} g "
                             f"({limit}); a field of the craft would need |∇U| = a/c² ≈ {b.now_g / C_SI ** 2:.2e} 1/m, "
                             f"max {b.max_g / C_SI ** 2:.2e} 1/m. Neither the theory nor the lab has a source or "
                             "mechanism for it.")
            else:
                lines.append(f"   load on the crew now {g_now:.1f} g, at most {g_max:.1f} g → {limit}")
            if getattr(b, "landing_speed", None):
                lines.append(f"   last touchdown at {b.landing_speed:.1f} m/s" + (" (crash!)" if b.landing_speed > 10
                                                                                 else ""))
        return lines

    def readouts_binary():
        lines = [f"{state['scene']}: real masses and orbit (2D). Weak field: causal gravity = GR; the periastron "
                 "advance is shown against 6πM/(a(1−e²)) per orbit."]
        for tr in state["trackers"]:
            rr = tr.rate()
            M = tr.center.mass + tr.body.mass
            a = 0.5 * (tr.r_min + tr.r_max) if tr.r_max > 0 else 0
            e = (tr.r_max - tr.r_min) / (tr.r_max + tr.r_min) if tr.r_max > 0 else 0
            if a > 0:
                period = 2 * math.pi * math.sqrt(a**3 / M)
                gr = 6 * math.pi * M / (a * (1 - e * e)) * (YEAR * C_SI / period)  # rad per year
                meas = rr * YEAR * C_SI if rr is not None else None
                if state["meta"].get("pulsar"):
                    txt = (f"{math.degrees(meas):.4f}°/yr" if meas is not None else "… (waiting for 3 periastra)")
                    lines.append(f"PERIASTRON {tr.body.name}: period {fmt_time(period / C_SI)}, e = {e:.4f}; advance "
                                 f"{txt} | GR {math.degrees(gr):.4f}°/yr | observed "
                                 f"{HULSE_TAYLOR['omdot_obs']:.4f}°/yr")
                else:
                    txt = (f"{meas * 100 / ARCSEC:.4f}″/century" if meas is not None else "… (waiting for 3 periastra)")
                    lines.append(f"PERIASTRON {tr.body.name} around {tr.center.name}: period {fmt_time(period / C_SI)}, "
                                 f"a = {fmt_len(a)}, e = {e:.3f}; advance {txt} | 1PN prediction "
                                 f"{gr * 100 / ARCSEC:.4f}″/century")
        if state["meta"].get("pulsar"):
            ht = HULSE_TAYLOR
            res = state["binaries"][0].result() if state["binaries"] else None
            if res is None:
                lines.append("RADIATION: waiting for a full orbit (the radiated energy is the work of the reaction "
                             "force between periastra) …")
            else:
                P, lum, pdot, a, e = res
                gr = peters_pdot(state["binaries"][0].b1.mass, state["binaries"][0].b2.mass, P, e)
                lines.append(f"RADIATION: dP/dt = {pdot:.4e} measured over whole orbits | quadrupole formula "
                             f"{gr:.4e} | observed {ht['pdot_obs']:.3e} ± {ht['pdot_obs_err']:.0e} (intrinsic, "
                             f"Weisberg & Huang 2016) | luminosity {lum * C_SI**5 / 6.6743e-11:.3e} W")
                lines.append("   scalar frameworks (H1, F3): 1/6 of the quadrupole for circular orbits (monopole and "
                             "dipole radiation of eccentric orbits not modelled); Newton: no radiation.")
        if state["scene"].startswith("Sirius"):
            M_B, R_B = 1.018 * M_SUN, 0.0081 * R_SUN
            z = 1 / math.sqrt(1 - 2 * M_B / R_B) - 1
            lines.append(f"SIRIUS B (white dwarf): gravitational redshift from the surface c·z = {z * C_SI / 1000:.1f} "
                         "km/s for M = 1.018 M☉, R = 0.0081 R☉ (causal gravity = GR here); observed 80.0 ± 1.6 km/s "
                         "(Joyce et al. 2018, HST). The value depends mainly on the radius of the star.")
        return lines + clock_lines()

    def readouts_map():
        sim = state["sim"]
        sun = sim.bodies[0]
        yrs = sim.t / C_SI / YEAR
        lines = [f"SOLAR NEIGHBOURHOOD, {yrs:,.0f} years after J2000 (initial positions J2000). The stars move in "
                 "straight lines with their real space velocity relative to the Sun (proper motion RECONS, radial "
                 "velocity Kervella et al. 2019 where known); mutual gravity and the Galaxy are ignored, a good "
                 "approximation over thousands of years. The field of the stars is tiny here (U ~ 10⁻¹³).",
                 "Editing: add an object with the form or by clicking with the tool 'add body'; the tool 'delete' "
                 "removes the nearest marked star. New objects move with the given space velocity."]
        rows = []
        for b in sim.bodies[1:]:
            r = (b.x - sun.x, b.y - sun.y, b.z - sun.z)
            d = math.sqrt(sum(c * c for c in r))
            v = getattr(b, "v3", (0.0, 0.0, 0.0))
            v2 = sum(c * c for c in v)
            extra = ""
            if v2 > 0 and "unknown" not in b.motion_note:  # closest approach on the straight path
                t_min = -sum(r[k] * v[k] for k in range(3)) / v2
                d_min = math.sqrt(sum((r[k] + v[k] * t_min) ** 2 for k in range(3)))
                when = (sim.t + t_min) / C_SI / YEAR
                extra = f"; closest {d_min / LY:.2f} ly in the year {2000 + when:,.0f}"
            rows.append((d, b, extra))
        for d, b, extra in sorted(rows, key=lambda t: t[0]):
            lines.append(f"   {b.name:20s} {d / LY:6.2f} ly   {b.info.split(', ')[0]:22s} {b.motion_note}{extra}")
        return lines

    def readouts_probes():
        sim = state["sim"]
        lines = ["FAST PROBES AND LIGHT: bending relative to Newton's value for the same speed and path. Prediction "
                 "(1PN): 1 + γv² × Newton; causal gravity, Einstein and the scalar stage γ = 1, the fixed order F3 "
                 "γ = −1, Newton 1. Light (v = 1): twice Newton's corpuscle in causal gravity."]
        for lt in [t for t in state["light"] if t.ph.alive]:
            st = lt.state(sim)
            if st.get("direct"):
                lines.append(f"{lt.ph.name}: aimed at the source, no weak-field prediction")
                continue
            v2 = lt.v0 ** 2
            ratio = st["defl"] / st["newton"] if st.get("newton") else float("nan")
            lines.append(f"{lt.ph.name} (b = {st['b']:.2f}): bending {fmt_angle(st['defl'])} = {ratio:.4f} × "
                         f"Newton | predicted: causal 1 + v² = {1 + v2:.4f}, Newton 1, F3 1 − v² = {1 - v2:.4f}")
        return lines

    def readouts_dragging():
        sim = state["sim"]
        star = next((b for b in sim.bodies if b.alive and any(b.J)), None)
        lines = ["FRAME DRAGGING (exaggerated: J = M², an extremal spin). The spinning star drags the light cones: "
                 "prograde orbits precess less, retrograde more (Lense–Thirring), a gyroscope precesses "
                 "geodetically and is dragged. Frameworks without g_0i (Newton, scalar stage, F3) show no "
                 "difference."]
        if star is None:
            return lines + ["The spinning star was deleted."]
        M, J = star.mass, star.J[2]
        adv = {}
        for tr in [t for t in state["trackers"] if t.body.alive]:
            res = tr.result()
            if not res:
                lines.append(f"PRECESSION {tr.body.name}: waiting for two orbits …")
                continue
            a_, pred_s, a, e, n = res
            lt_ = -8 * math.pi * J / (math.sqrt(M * a**3) * (1 - e * e) ** 1.5) * tr.sense
            adv[tr.sense] = (a_, pred_s, lt_)
            lines.append(f"PRECESSION {tr.body.name}: {math.degrees(a_):+.3f}°/orbit ({n} orbits, a = {a:.3f}, "
                         f"e = {e:.3f}) | Schwarzschild {math.degrees(pred_s):+.3f}° + Lense–Thirring "
                         f"{math.degrees(lt_):+.3f}° = {math.degrees(pred_s + lt_):+.3f}°")
        if 1.0 in adv and -1.0 in adv:
            lines.append(f"   prograde − retrograde: {math.degrees(adv[1.0][0] - adv[-1.0][0]):+.3f}°/orbit (prediction "
                         f"{math.degrees(adv[1.0][2] - adv[-1.0][2]):+.3f}°: the pure frame-dragging part)")
        for b in [x for x in sim.bodies if x.alive and x.gyro is not None]:
            rates = gyro_rates(b)
            if rates is None:
                continue
            r = math.hypot(b.x - star.x, b.y - star.y)
            P = 2 * math.pi * math.sqrt(r**3 / M)
            lines.append(f"GYROSCOPE {b.name}: in-plane precession per orbit: geodetic {math.degrees(rates[0] * P):+.3f}° "
                         f"(prediction (γ + ½)·2πM/r = {math.degrees(3 * math.pi * M / r):+.3f}°), frame dragging "
                         f"{math.degrees(rates[1] * P):+.4f}° (prediction −2πJ/√(Mr³) = "
                         f"{math.degrees(-2 * math.pi * J / math.sqrt(M * r**3)):+.4f}°)")
        return lines

    def readouts_gpb():
        sim = state["sim"]
        earth = next((b for b in sim.bodies if b.alive and b.radius > 0), None)
        lines = ["GRAVITY PROBE B (2004–2005): a gyroscope on a polar orbit (the Earth's spin axis lies in the "
                 "orbital plane, dashed arrow), spin pointing at a guide star in the plane. Geodetic precession turns "
                 "the spin within the plane; frame dragging by the rotating Earth tilts it out of the plane."]
        for b in [x for x in sim.bodies if x.alive and x.gyro is not None]:
            rates = gyro_rates(b)
            if rates is None or earth is None:
                continue
            k = YEAR * C_SI / MAS
            r = GPB_A if state["scene"].startswith("Gravity Probe B") else math.hypot(b.x - earth.x, b.y - earth.y)
            P = 2 * math.pi * math.sqrt(r**3 / earth.mass)
            geo_p = 1.5 * earth.mass**1.5 / r**2.5 * k
            fd_p = abs(earth.J[1]) / (2 * r**3) * k
            avg = (f"averaged over {rates[3]} whole orbit{'s' if rates[3] > 1 else ''}" if rates[3]
                   else f"running average over {b.gyro['t'] / P:.2f} orbit (exact only after a whole orbit)")
            lines.append(f"GYROSCOPE {b.name}: geodetic {rates[0] * k:.1f} mas/yr (GR, this orbit {geo_p:.1f}); frame "
                         f"dragging {abs(rates[2]) * k:.2f} mas/yr (GR, orbit average {fd_p:.2f}); {avg}")
            lines.append("   GP-B mission (Everitt et al. 2011): predicted 6606.1 and 39.2 mas/yr, measured 6601.8 ± "
                         "18.3 and 37.2 ± 7.2 mas/yr. Newton: 0 and 0; scalar stage: 6606 and 0 (no g_0i); F3: −2202 "
                         "and 0.")
        return lines

    def readouts_inspiral():
        sim = state["sim"]
        meta = state["meta"]
        lines = ["INSPIRAL: two black holes of 0.5 each, initial separation "
                 f"{meta.get('a0', 0):g} (G = c = 1), started on a circular orbit of the causal metric. Radiation reaction (2.5PN) "
                 "shrinks the orbit; below separation "
                 f"{0.5 / L_MAX:.1f} the weak-field guard merges them. Scalar frameworks radiate 1/6, Newton nothing."]
        trs = state["binaries"]
        if not trs:
            return lines + ["The binary is gone (merged or deleted)."]
        tr = trs[0]
        m1, m2 = tr.b1.mass, tr.b2.mass
        a0 = meta.get("a0", 20.0)
        M, mu = m1 + m2, m1 * m2 / (m1 + m2)
        x = a0**4 - 256.0 / 5.0 * mu * M * M * sim.t
        a_p = x ** 0.25 if x > 0 else 0.0
        T = peters_merge_time(m1, m2, a0)
        if tr.b1.alive and tr.b2.alive:
            lines.append(f"t = {sim.t:.1f}: separation {tr.sep():.3f} | quadrupole (Peters, circular) {a_p:.3f} | "
                         f"radiated energy {tr.energy(sim):.4e} | Peters time to coalescence {T:.0f}")
        else:
            lines.append(f"MERGED (see the events below); Peters time to coalescence from {a0:g}: {T:.0f}")
        return lines + sim.events[-3:]

    def measure_rod_probes(sim):
        """Fall acceleration of the probes released above the beams: a = 2 |Δy| / t² once |Δy| ≥ 0.1."""
        for b in sim.bodies:
            rod = getattr(b, "rod", None)
            if rod is None or not b.alive or getattr(b, "a_meas", None) is not None or not hasattr(b, "start"):
                continue
            dy = abs(b.y - b.start[1])
            t = sim.t - (b.t_birth or 0.0)
            if (dy >= 0.1 or b.landed_t is not None) and t > 0:
                b.a_meas = 2.0 * dy / (t * t)

    def readouts_rods():
        sim = state["sim"]
        lines = ["WHAT GRAVITATES (F7.4a): each rod is a pair of beams streaming in opposite directions with the same "
                 "REST density (the light pair: the same energy density). Momentum cancels, energy and stress remain. "
                 "A probe released 1 unit above each rod falls with an acceleration proportional to the active "
                 "mass. Candidates: tensor Γ²(1 + v²) (causal, Einstein; lab F7.4a PASSED: 2.401 ± 0.053 at 0.6c), "
                 "energy only Γ², lab-frame density Γ, rest count 1 (Newton, scalar stage)."]
        ref = next((b for b in sim.bodies if getattr(b, "rod", None) is not None and b.rod.vb == 0.0
                    and not b.rod.light), None)
        for b in sim.bodies:
            rod = getattr(b, "rod", None)
            if rod is None or not b.alive:
                continue
            if rod.light:
                cand = "causal/Einstein 2, Newton 1, scalar stage 0 (per unit energy)"
            else:
                G = lorentz(rod.vb)
                cand = (f"tensor {G * G * (1 + rod.vb ** 2):.3f}, energy {G * G:.3f}, lab density {G:.3f}, "
                        "rest count 1")
            meas = (b.a_meas / ref.a_meas) if (ref is not None and getattr(b, "a_meas", None) is not None
                                               and getattr(ref, "a_meas", None) is not None) else None
            lines.append(f"{rod.name}: pull ÷ beams at rest " + (f"{meas:.3f}" if meas is not None else "… (falling)")
                         + f" | {cand}")
        return lines

    def readouts():
        sim = state["sim"]
        meta = state["meta"]
        if meta.get("map3d"):
            set_readout("\n".join(readouts_map()))
            return
        if meta.get("bh"):
            body = readouts_bh() + comparison_lines()
        elif meta.get("pilot"):
            note = (["Frameworks: the steered UAP is not computed in them (the keys steer only the main scene)."]
                    if state["ghosts"] else [])
            body = readouts_pilot() + readouts_si()[1:] + note + comparison_lines()
        elif meta.get("orbits") == "binary":
            body = readouts_binary() + comparison_lines()
        elif meta.get("inspiral"):
            body = readouts_inspiral() + comparison_lines()
        elif meta.get("gpb"):
            body = readouts_gpb() + comparison_lines()
        elif meta.get("rods"):
            body = readouts_rods() + comparison_lines()
        elif meta.get("probes"):
            body = readouts_probes() + comparison_lines()
        elif meta.get("dragging"):
            body = readouts_dragging() + comparison_lines()
        elif meta["si"]:
            body = (readouts_solar() if meta.get("orbits") else readouts_si()) + comparison_lines() + sim.events[-3:]
        else:
            body = readouts_default()
        set_readout("\n".join(body + experiment_lines() + tests_lines()))

    def readouts_default():
        sim = state["sim"]
        law = sim.law
        lines = []
        for tr in [t for t in state["trackers"] if t.body.alive and t.center.alive]:
            r = tr.result()
            if r:
                adv, pred, a, e, n = r
                line = (f"PRECESSION {tr.body.name}: measured {math.degrees(adv):+.3f}°/orbit ({n} orbits, a = {a:.3f}, "
                        f"e = {e:.3f}) | causal (β = γ = 1): {math.degrees(pred * law.perihelion_factor()):.3f}° | "
                        f"β = 1/2: {math.degrees(pred * 7 / 6):.3f}° | β = 3/2: {math.degrees(pred * 5 / 6):.3f}°")
                # the 1PN effect is only ~3M/p of the Newtonian pull, so weak pulls of other bodies can matter
                rel_1pn = 3.0 * (tr.center.mass + tr.body.mass) / (a * (1.0 - e * e))
                if tr.other_max > 1.0:
                    line += (f"  !! DISTURBED: another body pulls harder than {tr.center.name} (up to "
                             f"{tr.other_max:.0f}×); not an orbit around one mass, the prediction does not apply")
                elif tr.other_max > 0.1 * rel_1pn:
                    line += (f"  !! DISTURBED: the other bodies pull up to {100 * tr.other_max:.2f} % of "
                             f"{tr.center.name} (the 1PN effect is only {100 * rel_1pn:.2f} %); Newtonian "
                             "perturbations may dominate")
                lines.append(line)
            else:
                lines.append(f"PRECESSION {tr.body.name}: waiting for two orbits (nearly circular orbits, e < 0.05, are "
                             "not reported) …")
        for lt in [t for t in state["light"] if t.ph.alive and t.m.alive]:
            st = lt.state(sim)
            if st.get("direct"):
                lines.append(f"LIGHT {lt.ph.name}: aimed at {lt.m.name} (b = {st['b']:.3g}, closer than "
                             f"{lt.m.mass / L_MAX:.3g} where |U| > {L_MAX:g}: it will be absorbed); no prediction")
                continue
            if not st["photon"]:
                lines.append(f"PROBE {lt.ph.name}: bending {fmt_angle(st['defl'])} = "
                             f"{st['defl'] / st['newton']:.4f} × Newton (causal 1 + v² = {1 + lt.v0 ** 2:.4f})")
                continue
            lines.append(f"LIGHT b={st['b']:.2f}: bending {math.degrees(st['defl']):.3f}° | causal γ = 1: "
                         f"{math.degrees(st['pred']):.3f}° (1st order {math.degrees(st['pred1']):.3f}°) | Newton "
                         f"(corpuscle): {math.degrees(st['newton']):.3f}° | F3 γ = −1: 0° | Shapiro {st['delay']:+.4f} "
                         f"against {st['delay_pred']:+.4f}"
                         + (f"  !! DISTURBED: other bodies {100 * st['other']:.0f} % of the pull of {lt.m.name}; the "
                            "prediction counts only with it" if st["other"] > 0.02 else ""))
        for b in sim.bodies:
            if b.alive and (b.kind == "clock" or (b.kind == "mass" and b.mass < 1e-3 and not b.fixed and not b.probe)):
                U = -sim.L_at(b.x, b.y, skip=b)
                age = sim.t - (b.t_birth or 0.0)  # proper time since the body entered the scene
                rate = b.tau / age if age > 0 else float("nan")
                lines.append(f"CLOCK {b.name}: τ = {b.tau:.3f}, τ/t since insertion = {rate:.6f} | static clock here "
                             f"now: √A = {1 + float(law.static_rate_m1(min(U, 20.0))):.6f}")
        if sim.events:
            lines.extend(sim.events[-4:])
        head = ["Predictions are 1PN for a test body at ONE mass (light with the second order of the metric); "
                "differences of order M/r are higher orders, not errors. Disturbance by other bodies is flagged."]
        lines += comparison_lines()
        return head + (lines if lines else ["No measurements in this scene (add a body or a photon)."])

    # ---------- time step
    def pilot_control(sim):
        """
        The pilot from the keys. Normal mode: thrust = direction × 10^slider g (Shift ×10) plus hovering (1 g up);
        the craft has inertia. 'Null field' mode (the author's HYPOTHESIS): the keys command a VELOCITY
        (10^(slider+2) m/s, Shift ×10; no key = stop), reached within TAU_NULL whatever the gravity -- instant stops
        and turns, crew at 0 g; shown is what a conventional thrust (or a field gradient of the craft) would need.
        """
        for b in sim.bodies:
            if b.kind != "pilot" or not b.alive:
                continue
            dx = (1 if ({"Right", "d", "D"} & keys) else 0) - (1 if ({"Left", "a", "A"} & keys) else 0)
            dy = (1 if ({"Up", "w", "W"} & keys) else 0) - (1 if ({"Down", "s", "S"} & keys) else 0)
            shift = 10.0 if ({"Shift_L", "Shift_R"} & keys) else 1.0
            if nullf_var.get():
                V = 10 ** (thrust_var.get() + 2.0) * shift  # m/s
                n = math.hypot(dx, dy)
                vt = (V * dx / n / C_SI, V * dy / n / C_SI) if n else (0.0, 0.0)
                b.vtarget, b.thrust = vt, (0.0, 0.0)
                gx = gy = 0.0  # Newtonian gravitational acceleration at the craft [m/s^2]
                for o in sim.bodies:
                    if o.alive and o.mass > 0 and o is not b:
                        ddx, ddy = b.x - o.x, b.y - o.y
                        r3 = (ddx * ddx + ddy * ddy) ** 1.5
                        gx -= o.mass * ddx / r3 * C_SI**2
                        gy -= o.mass * ddy / r3 * C_SI**2
                ax = (vt[0] - b.vx) * C_SI**2 / TAU_NULL  # kinematic acceleration [m/s^2]
                ay = (vt[1] - b.vy) * C_SI**2 / TAU_NULL
                b.v_cmd = V if n else 0.0
                b.now_g = math.hypot(ax - gx, ay - gy)  # what a thrust would have to supply
                b.max_g = max(b.max_g, b.now_g)
                continue
            b.vtarget = None
            mag = 9.81 * 10 ** thrust_var.get() * shift
            n = math.hypot(dx, dy)
            ax, ay = (mag * dx / n, mag * dy / n) if n else (0.0, 0.0)
            if getattr(b, "hover", False):  # compensation of the local gravity (Newtonian) [m/s^2]
                g = sim.pull_split(b, None)[1] * C_SI**2
                ex = next((o for o in sim.bodies if o.radius > 0), None)
                if ex is not None:
                    r = math.hypot(b.x - ex.x, b.y - ex.y)
                    ax += g * (b.x - ex.x) / r
                    ay += g * (b.y - ex.y) / r
            b.thrust = (ax, ay)
            b.now_g = math.hypot(ax, ay)
            b.max_g = max(b.max_g, b.now_g)

    def tick(single=False):
        sim = state["sim"]
        now = time.perf_counter()
        real_dt = min(max(now - state.get("last_tick", now - 0.03), 0.001), 0.1)
        state["last_tick"] = now
        # the same control in every scene: simulated time per second of real time; SI scenes store time as ct
        unit_ct = C_SI if state["meta"].get("si") else 1.0
        T = 10.0 ** speed_var.get() * real_dt * unit_ct
        if state["meta"].get("map3d"):  # star map: every star moves with its space velocity
            hw = state["view"][2]
            for b in sim.bodies:
                vx, vy, vz = getattr(b, "v3", (0.0, 0.0, 0.0))
                if vx or vy or vz:
                    b.x, b.y, b.z = b.x + vx * T, b.y + vy * T, b.z + vz * T
                    tr = b.trail3
                    if not tr or math.dist(tr[-1], (b.x, b.y, b.z)) > 0.004 * hw:
                        tr.append((b.x, b.y, b.z))
            sim.t += T
            state["achieved"] = T / C_SI / real_dt
            speed_value_label()
            draw(refresh_field=False)
            return
        if state["meta"].get("pilot"):
            pilot_control(sim)
        t_start = sim.t
        # frameworks that lag behind (heavier to integrate) get the frame: the main scene waits for them
        lagging = any(g["sim"].t < sim.t - 1e-9 * max(abs(sim.t), 1.0) for g in state["ghosts"])

        def on(sm, trackers, light, binaries):
            for tr in trackers:
                tr.update(sm)
            for lt in light:
                lt.update(sm)
            for bt in binaries:
                if bt.b1.alive and bt.b2.alive:
                    bt.update(sm)

        t0 = time.perf_counter()
        if single:
            sim.step()
            on(sim, state["trackers"], state["light"], state["binaries"])
        elif lagging:
            state["achieved"] = 0.0
        else:
            t_end = sim.t + T
            while sim.t < t_end and time.perf_counter() - t0 < 0.025:  # keeps the window responsive
                live, pos, vel, mass, photon, _ = sim.arrays()
                dt = min(sim.choose_dt(pos, vel, mass, photon), t_end - sim.t)
                sim.step(max(dt, 1e-9))
                on(sim, state["trackers"], state["light"], state["binaries"])
            state["achieved"] = (sim.t - t_start) / unit_ct / real_dt  # may be lower than asked
            speed_value_label()
        if state["meta"].get("pilot") and follow_var.get():
            pl = next((b for b in sim.bodies if b.kind == "pilot" and b.alive), None)
            if pl is not None:
                state["view"] = (pl.x, pl.y, state["view"][2])
        # frameworks: the same scene under other laws, advanced to the same time with their own steps
        t1 = time.perf_counter()
        behind = [g for g in state["ghosts"] if g["sim"].t < sim.t]
        while behind and time.perf_counter() - t1 < 0.06:  # round robin: no framework starves the others
            for g in behind:
                gs = g["sim"]
                live, pos, vel, mass, photon, _ = gs.arrays()
                dt = min(gs.choose_dt(pos, vel, mass, photon), sim.t - gs.t)
                gs.step(max(dt, 1e-9))
                on(gs, g["trackers"], g["light"], g["binaries"])
            behind = [g for g in behind if g["sim"].t < sim.t]
        if state["meta"].get("rods"):
            for sm in [sim] + [g["sim"] for g in state["ghosts"]]:
                measure_rod_probes(sm)
        background_work(0.015)
        if sim.events and state.get("n_events", 0) != len(sim.events):
            state["n_events"] = len(sim.events)
            setup_trackers()  # a body merged or a photon was absorbed
        for sm in [sim] + [g["sim"] for g in state["ghosts"]]:
            for b in sm.bodies:
                if b.alive and not b.fixed and b.kind != "rod":
                    if not b.trail or math.hypot(b.x - b.trail[-1][0], b.y - b.trail[-1][1]) > state["view"][2] * 0.003:
                        b.trail.append((b.x, b.y))
        state["frame"] += 1
        draw()

    def loop():
        try:
            if state["running"]:
                tick()
            elif state["exp"] is not None or state["tests"] is not None:  # background jobs run also when paused
                background_work(0.025)
        except Exception as exc:  # noqa: BLE001 -- an error of one frame must not stop the animation
            report_error("frame", exc)
        finally:
            root.after(30, loop)

    # ---------- 3D map tools and mouse
    def map3d_point_from_click(sx, sy):
        """A click in the 3D map as a point in the plane of view (depth y2 = 0 of the current projection)."""
        w, h = max(canvas.winfo_width(), 10), max(canvas.winfo_height(), 10)
        hw = state["view"][2]
        scale = 0.5 * w / hw
        yaw, pitch = state["rot"]
        x1 = (sx - w / 2) / scale
        z2 = (h / 2 - sy) / scale
        y1 = math.sin(pitch) * z2
        z = math.cos(pitch) * z2
        x = math.cos(yaw) * x1 + math.sin(yaw) * y1
        y = -math.sin(yaw) * x1 + math.cos(yaw) * y1
        return x, y, z

    def delete_nearest_3d(sx, sy):
        """Removes the nearest point of the map within 16 pixels of the click."""
        best, dist = None, 16.0
        for px, py, b in state.get("map3d_points", {}).values():
            d = math.hypot(px - sx, py - sy)
            if d < dist:
                best, dist = b, d
        if best is not None:
            if best.name == "Sun":
                report_error("3D map", ValueError("The Sun is the origin of the coordinates and cannot be removed."))
                return
            state["sim"].remove(best)
            update_mass_info()
            draw(refresh_field=False)

    def on_press(ev):
        if state["meta"].get("map3d"):
            tool = tool_var.get()
            if tool == "delete":
                delete_nearest_3d(ev.x, ev.y)
            elif tool == "mass":
                state["map3d_add_click"] = (ev.x, ev.y)
            else:
                state["rot_drag"] = (ev.x, ev.y, list(state["rot"]))
            return
        if state.get("view3d") and tool_var.get() == "none":
            state["rot_drag"] = (ev.x, ev.y, list(state.get("rot", [0.25, 0.95])))
            return
        if tool_var.get() == "delete":  # removes the nearest object within 15 px
            sim = state["sim"]
            best, bd = None, 15.0
            for b in sim.bodies:
                if b.alive:
                    sx, sy = W2S(b.x, b.y)
                    d = math.hypot(sx - ev.x, sy - ev.y)
                    if d < bd:
                        best, bd = b, d
            if best is not None:
                sim.remove(best)  # removing the black hole also removes its exact field
                for g in state["ghosts"]:  # the same unique name in every framework
                    for gb in [b for b in g["sim"].bodies if b.name == best.name]:
                        g["sim"].remove(gb)
                setup_trackers()
                draw(force_field=True)
            return
        if tool_var.get() != "none":
            state["drag"] = ((ev.x, ev.y), (ev.x, ev.y))

    def on_motion(ev):
        if (state["meta"].get("map3d") or state.get("view3d")) and state.get("rot_drag"):
            x0, y0, (yw, pt) = state["rot_drag"]
            state["rot"] = [yw + (ev.x - x0) * 0.01, max(-1.45, min(1.45, pt + (ev.y - y0) * 0.01))]
            draw(refresh_field=False)
            return
        if state["drag"] is not None:
            state["drag"] = (state["drag"][0], (ev.x, ev.y))
            if not state["running"]:
                draw()

    def on_release(ev):
        if state["meta"].get("map3d"):
            start = state.get("map3d_add_click")
            state["map3d_add_click"] = None
            if start is not None:
                x3, y3, z3 = map3d_point_from_click(ev.x, ev.y)
                for var, val in zip(map_xyz_vars, (x3, y3, z3)):
                    var.set(f"{val / LY:.8g}")
                add_map3d_object(pos_m=(x3, y3, z3))
            state["rot_drag"] = None
            return
        if state.get("view3d") and state.get("rot_drag") is not None:
            state["rot_drag"] = None
            return
        if state["drag"] is None:
            return
        (x0, y0), (x1, y1) = state["drag"]
        state["drag"] = None
        sim = state["sim"]
        wx, wy = S2W(x0, y0)
        ex, ey = S2W(x1, y1)
        dvx, dvy = ex - wx, ey - wy
        tool = tool_var.get()
        n = len(sim.bodies)
        state["name_seq"] += 1
        seq = state["name_seq"]  # sequence numbers are not reused after deleting, names stay unique
        if tool in ("mass", "clock"):
            for o in [b for b in sim.bodies if b.alive and b.radius > 0]:  # a click inside the Earth moves to the surface
                rx, ry = wx - o.x, wy - o.y
                r = math.hypot(rx, ry)
                if r < o.radius:
                    ux, uy = (rx / r, ry / r) if r > 0 else (0.0, 1.0)
                    wx, wy = o.x + ux * o.radius * (1 + 1e-9), o.y + uy * o.radius * (1 + 1e-9)
            c = next((b for b in sim.bodies if getattr(b, "is_profile", False) and b.alive), None)
            if c is not None:  # not exactly at the centre of the black hole
                rx, ry = wx - c.x, wy - c.y
                r = math.hypot(rx, ry)
                r_min = 0.05 * c.mass
                if r < r_min:
                    ux, uy = (rx / r, ry / r) if r > 0 else (1.0, 0.0)
                    wx, wy = c.x + ux * r_min, c.y + uy * r_min
        if tool == "mass":
            m_len = parse_mass()
            if m_len is None:
                report_error("mass input", ValueError("The mass must be a positive number in scientific notation, "
                                                      "e.g. 5.972e24."))
                state["name_seq"] -= 1
                return
            meta = state["meta"]
            if meta["si"]:
                name = f"{category_si(m_len / G_C2)} {seq}"
            else:
                name = f"{category_relative(m_len / reference_mass(sim))} {seq}"
            if meta["click"] == "rest":  # a click releases the body; a drag gives it a velocity
                dpx, dpy = x1 - x0, y1 - y0
                dl = max(math.hypot(dpx, dpy), 1e-9)
                k = dl / 80.0 * (20.0 / C_SI if meta["si"] else 0.01)
                vx, vy = (k * dpx / dl, -k * dpy / dl) if dl >= 6 else (0.0, 0.0)
            else:
                vx, vy, _ = new_body_velocity(wx, wy, m_len, x0, y0, x1, y1)
            sim.bodies.append(Body(name, wx, wy, vx, vy, mass=m_len, color="#9ad9ff"))
            sim.bodies[-1].start = (wx, wy, vx, vy)
        elif tool == "photon":
            sp = math.hypot(dvx, dvy)
            if sp == 0:
                dvx, dvy, sp = 1.0, 0.0, 1.0
            sim.bodies.append(Body(f"photon {seq}", wx, wy, dvx / sp, dvy / sp, kind="photon", color="#fff27a"))
        elif tool == "clock":
            sim.bodies.append(Body(f"clock {seq}", wx, wy, kind="clock", color="#b0b0ff"))
            sim.bodies[-1].start = (wx, wy, 0.0, 0.0)
        else:
            state["name_seq"] -= 1
            return
        sim.bodies[-1].t_birth = sim.t
        if len(sim.bodies) > n:  # copy the new object into every framework
            nb_ = sim.bodies[-1]
            for g in state["ghosts"]:
                cp = Body(nb_.name, nb_.x, nb_.y, nb_.vx, nb_.vy, mass=nb_.mass, kind=nb_.kind, color=nb_.color)
                cp.t_birth, cp.start = g["sim"].t, getattr(nb_, "start", (nb_.x, nb_.y, nb_.vx, nb_.vy))
                g["sim"].bodies.append(cp)
        update_mass_info()
        setup_trackers()
        draw(force_field=True)

    def add_map3d_object(pos_m=None):
        """Adds an object to the 3D map; coordinates in ly, velocities in km/s, the mass from the mass field."""
        if not state.get("meta", {}).get("map3d"):
            return
        try:
            if pos_m is None:
                xyz = tuple(float(v.get().strip().replace(",", ".")) * LY for v in map_xyz_vars)
            else:
                xyz = tuple(float(v) for v in pos_m)
            vel = tuple(float(v.get().strip().replace(",", ".")) * 1000.0 / C_SI for v in map_v_vars)
            if not all(math.isfinite(v) for v in (*xyz, *vel)):
                raise ValueError("Coordinates and velocities must be finite numbers.")
        except ValueError as exc:
            report_error("3D map object", exc)
            return
        m_len = parse_mass()
        if m_len is None:
            report_error("mass input", ValueError("The mass must be a positive number in scientific notation."))
            return
        name = map_name_var.get().strip() or f"New body {state['name_seq'] + 1}"
        if any(b.alive and b.name == name for b in state["sim"].bodies):
            report_error("3D map object", ValueError(f"An object named '{name}' already exists."))
            return
        b = Body(name, xyz[0], xyz[1], mass=m_len, fixed=True, color="#9ad9ff")
        b.z = xyz[2]
        b.v3 = vel
        b.trail3 = deque(maxlen=4000)
        dist = math.sqrt(sum(v * v for v in xyz)) / LY
        b.info = f"user object, {fmt_sci(m_len / M_SUN, 3)} M☉, {dist:.4g} ly"
        v_kms = tuple(v * C_SI / 1000.0 for v in vel)
        b.motion_note = f"v = ({v_kms[0]:.4g}, {v_kms[1]:.4g}, {v_kms[2]:.4g}) km/s"
        b.t_birth = state["sim"].t
        b.start = (b.x, b.y, 0.0, 0.0)
        state["sim"].bodies.append(b)
        state["name_seq"] += 1
        map_name_var.set(f"New body {state['name_seq'] + 1}")
        update_mass_info()
        draw(refresh_field=False)

    def new_body_velocity(wx, wy, m, x0, y0, x1, y1):
        """
        Velocity of a new body: a click without dragging gives a circular orbit (counter-clockwise) around the heaviest
        source; a drag sets direction and multiple of the circular speed (80 px = 1×). Without a source a tenth of
        the drag in scene coordinates. Returns (vx, vy, multiple of the circular speed or None).
        """
        sim = state["sim"]
        src = [b for b in sim.bodies if b.alive and b.mass > 0]
        dpx, dpy = x1 - x0, y1 - y0
        dlen = math.hypot(dpx, dpy)
        if not src:
            ex, ey = S2W(x1, y1)
            return 0.1 * (ex - wx), 0.1 * (ey - wy), None
        c = max(src, key=lambda b: b.mass)
        rx, ry = wx - c.x, wy - c.y
        r = max(math.hypot(rx, ry), 1e-6)
        vc = math.sqrt((c.mass + m) / r)
        if dlen < 6:  # a click alone makes a circular orbit
            return c.vx - vc * ry / r, c.vy + vc * rx / r, 1.0
        k = dlen / 80.0
        return c.vx + vc * k * dpx / dlen, c.vy - vc * k * dpy / dlen, k  # the screen y axis points down

    def on_wheel(ev):
        delta = ev.delta if getattr(ev, "delta", 0) else (120 if ev.num == 4 else -120)
        cx, cy, hw = state["view"]
        wx, wy = S2W(ev.x, ev.y)
        f = 0.85 if delta > 0 else 1 / 0.85
        state["view"] = (wx + (cx - wx) * f, wy + (cy - wy) * f, hw * f)
        draw(force_field=True)

    def on_pan_press(ev):
        state["pan"] = (ev.x, ev.y, state["view"])

    def on_pan_motion(ev):
        if state["pan"] is None:
            return
        x0, y0, (cx, cy, hw) = state["pan"]
        w = max(canvas.winfo_width(), 10)
        sc = 0.5 * w / hw
        state["view"] = (cx - (ev.x - x0) / sc, cy + (ev.y - y0) / sc, hw)
        draw(force_field=True)

    canvas.bind("<ButtonPress-1>", safe("mouse", lambda ev: (canvas.focus_set(), on_press(ev))))
    canvas.bind("<B1-Motion>", safe("mouse", on_motion))
    canvas.bind("<ButtonRelease-1>", safe("mouse", on_release))
    for seq_ in ("<MouseWheel>", "<Button-4>", "<Button-5>"):
        canvas.bind(seq_, safe("mouse", on_wheel))
    canvas.bind("<ButtonPress-3>", safe("mouse", on_pan_press))
    canvas.bind("<B3-Motion>", safe("mouse", on_pan_motion))

    CONTROL_KEYS = {"Up", "Down", "Left", "Right", "w", "a", "s", "d", "W", "A", "S", "D", "space"}

    def on_key(ev):
        k = ev.keysym
        if not state["meta"].get("pilot"):
            return None
        if k in CONTROL_KEYS and not state["running"]:  # the control scene starts with the first control key
            toggle_run()
        if k == "space":
            if "space" not in keys:  # a held key must not toggle again
                for b in state["sim"].bodies:
                    if b.kind == "pilot":
                        b.hover = not getattr(b, "hover", False)
            keys.add("space")
            return "break"  # the key must not also press a focused button
        keys.add(k)
        return "break" if k in CONTROL_KEYS else None  # arrows must not move a focused slider or list

    def on_key_up(ev):
        keys.discard(ev.keysym)
        return "break" if state["meta"].get("pilot") and ev.keysym in CONTROL_KEYS else None

    # the canvas has the focus (after loading a scene or when the mouse enters); window bindings catch the rest
    canvas.bind("<KeyPress>", safe("keyboard", on_key))
    canvas.bind("<KeyRelease>", safe("keyboard", on_key_up))
    root.bind("<KeyPress>", safe("keyboard", lambda ev: on_key(ev) if ev.widget is not canvas else None))
    root.bind("<KeyRelease>", safe("keyboard", lambda ev: on_key_up(ev) if ev.widget is not canvas else None))
    root.bind("<FocusOut>", lambda ev: keys.clear())  # no key may stay 'pressed' when the window loses focus
    canvas.bind("<Enter>", lambda ev: canvas.focus_set())

    def on_hover(ev):
        state["hover"] = (ev.x, ev.y)
        if not state["running"]:
            draw(refresh_field=False)

    canvas.bind("<Motion>", safe("mouse", on_hover))
    canvas.bind("<Leave>", lambda e: state.update(hover=None))
    canvas.bind("<Configure>", safe("window", lambda e: draw(force_field=True)))

    if start_scene is not None:
        scene_var.set(SCENES[start_scene])  # the variable trace calls load_scene
    else:
        load_scene(SCENES[0])
    if autorun:
        toggle_run()
    loop()
    root.mainloop()


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        raise SystemExit(selftest())
    # options: --scene N opens the scene with index N, --autorun starts the simulation at once
    sc = None
    if "--scene" in sys.argv:
        try:
            sc = int(sys.argv[sys.argv.index("--scene") + 1]) % len(SCENES)
        except (IndexError, ValueError):
            sc = None
    run_app(sc, "--autorun" in sys.argv)
