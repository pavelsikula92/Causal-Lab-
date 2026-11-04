# Gravitational equations: the present formulation of the theory

Written on 27 Oct 2026 after the replication of F7.1c. It is an analytic summary, not the result of a new test. The
text summarises what the equation says, what of it is derived, what the laboratory has tested and what is still an
input. The decisions it rests on: 1 (field equation, fourth round), 6 (source) and 12 (gravity changes the causal
order), see `docs/author_decisions.md`.

Conventions: signature (−, +, +, +), c = 1, in 3+1. Π is the clock rate (dτ = Π dt along a static world line) and
L = ln Π. Near matter L < 0. M denotes GM/c² (a length). Δ and ∇ are the flat coordinate operators of isotropic
coordinates.

## 1. The equation

**Discretely (causal set):**

    B_g ln Π = κ ρ_s

B_g is the Benincasa–Dowker operator on an order with the geometry g. ρ_s is the density of saturated events
(decision 6: S_v = 1 when ΔQ_v ≥ ΔQ_crit).

**The metric (H1) that the equation assumes:**

    ds² = −Π² dt² + Π⁻² dx² = −e^{2L} dt² + e^{−2L} dx²

**Continuum limit:** B_g → □_g − R/2. In the H1 metric exactly:

    R = 2 e^{2L} (ΔL − |∇L|²)
    (□_g − R/2) L = e^{2L} [ (1 − L) ΔL + L |∇L|² ] = κ ρ_s

**Weak field:** ΔL = κρ_s, so L = −U (the Newtonian potential). The calibration is κρ_s ↔ 4πGρ.

**Vacuum:** (1 − L) ΔL + L |∇L|² = 0.

## 2. Where the metric comes from: the principle of conserved causal throughput

In detail in `why_gamma1.md` (sections 4 and 7). In short:

- The metric is read from the order by three instruments: clocks (number of events in a diamond) → φ; the boundary of
  the relation (speed of the cone) → φ + ψ; the coordinate density → φ − d_s ψ.
- Principle: the density of events times the speed of the cone is constant (ρ_c · c = const; book 5.3, 23, 26.1). In a
  weak field it follows that γ = 2/(d_s − 1), so γ = 1 in 3+1.
- With the clocks Π, H1 is fixed exactly: A·B = 1, c_eff = Π², ρ_c = Π⁻².
- A6 (S = Σ Π Δτ) gives the effective metric Π² g_u. H1 arises exactly when g_u has the cone speed Π².
- Schwarzschild satisfies A·B = 1 only to first order: A·B = (1 − M²/4r²)².
- Two alternatives give another γ: a fixed order gives γ = −1 (light does not bend) and the narrowing of the cones
  alone γ = 1/d_s.

## 3. Checks the equation has passed (analytically, 27 Oct)

All points were recomputed in the container (sympy and scipy).

**(a) Curvature and the form of the equation.** R and the left-hand side (□_g − R/2)L were computed from the metric
(Christoffel symbols → Ricci tensor) for a general L(r). The difference from the formulas of section 1 is identically 0.

**(b) A first integral and a closed vacuum solution (new).**

- With spherical symmetry and w = r²L′ the vacuum equation reads (1 − L) w′ + L L′ w = 0, i.e.
  d ln w / dL = 1 − 1/(1 − L).
- An exact first integral follows:

      r² L′ = M e^{L} (1 − L)

  M is the asymptotic charge (r²L′ → M for r → ∞).
- The second integration leads to the exponential integral:

      Ei(1 − L) = Ei(1) + e M / r

- A numerical integration of the ODE from r = 200 M agrees with the closed form to 5·10⁻⁸.
- Expansion: L = −M/r + M³/(6r³) − M⁴/(12r⁴) + … The coefficient of M²/r² is zero; this is a = 0 from before.

**(c) PPN.**

- g_tt = −e^{2L} = −(1 − 2M/r + 2M²/r² − …), so β = 1.
- g_rr = e^{−2L} = 1 + 2M/r + …, so γ = 1.
- γ now follows from the principle of throughput (section 2), not only from the ansatz.
- β = 1 + a, where a is the coefficient of M²/r² in L. The vacuum equation forces a = 0.
- β is, however, not an independent prediction, because f = ln Π was selected by the perihelion (decision 1, fourth
  round). Within the family (Π^k − 1)/k only ln Π gives β = 1.

**(d) The second post-Newtonian order (a prediction, not yet tested).**

| isotropic, x = M/r | our equation | general relativity |
|---|---|---|
| −g_tt | 1 − 2x + 2x² − **1**·x³ | 1 − 2x + 2x² − **3/2**·x³ |
| g_rr | 1 + 2x + **2**·x² + **1**·x³ | 1 + 2x + **3/2**·x² + **1/2**·x³ |

The first deviation is in g_rr at second order. It is exactly the failure of A·B = 1 for Schwarzschild. It would be
measurable at second order of the bending of light or of the Shapiro delay; no measurement that precise exists yet.

## 4. Strong field (the vacuum solution of a point mass)

> **Not valid as a prediction of the theory (decision 13, 28 Oct, path A):** the principle of throughput is a law of
> linear order and the non-linear completion is general relativity. The throat, the missing horizon, the +0.91 % of
> the shadow and the 2PN coefficients of point 3d are therefore not predictions of the theory. The section stays as a
> record of the properties of the scalar equation outside its validity.

The values follow from the closed solution of point 3b.

| quantity | our equation | general relativity |
|---|---|---|
| horizon | none: g_tt = e^{2L} > 0 for every r > 0 | r_areal = 2M |
| photon sphere | rL′ = 1/2: r_iso = 1.798 M, L = −0.535, r_areal = 2M(1 − L) = **3.070 M** | 3 M |
| critical impact parameter | b_crit = min r e^{−2L} = **5.2435 M** | 3√3 M = 5.1962 M (+0.91 %) |
| throat | rL′ = 1: r_iso = 0.612 M, L = −1.341, g_tt = 0.068, r_areal = M(1 − L) = **2.341 M** | – |

**Reading.**

- **Shadow.** It is 0.91 % larger than in general relativity. The present precision of shadow images (EHT, errors of
  order 10 %) cannot resolve it, but it is a concrete number.
- **Throat.**
  - Below 2.341 M the area of the spheres grows again.
  - g_tt goes to 0 only for r → 0, and at an infinite proper distance (∫ e^{−L} dr diverges).
  - A purely vacuum point mass is therefore a throat as for the exponential metric, not a black hole with a horizon.
  - For comparison: the exponential metric L = −M/r has a throat at e·M = 2.718 M, a photon sphere at 3.297 M and
    b_crit 5.437 M. Our solution lies between it and general relativity.
- **Physical relevance of the throat.**
  - It depends on whether a body with a surface below the photon sphere or below the throat can exist statically at
    all.
  - F1.7 N found a bound on static solutions for a Gaussian source (ln Π: u_lin(0) between −0.3 and −0.5).
  - The conversion to the compactness GM/(Rc²) has, however, not been done and the source has no pressure. A
    comparison with neutron stars and with collapse is therefore not yet possible.
- **Limits.** All points of this section are the continuum limit. The laboratory has no discrete solution in a strong
  field.

## 5. The side of the source (open)

> **Added 27 Oct (`tensor_coupling.md`):** at linear order the source is the T_ab of the flow, not the scalar ρ_s.
> In 3+1 the principle of throughput is the same as a pure tensor coupling (h̄_ab = 4U u_a u_b). A scalar as a
> complete dynamical law is excluded: the quadrupole radiation of a circular binary is 1/6 of the value of general
> relativity and frame dragging is missing. The points below describe the scalar formulation.

- ρ_s is a scalar, the density of saturated events. The equation is scalar within the fixed ansatz of the H1 metric,
  it is not tensorial.
- Light, pressure and motion do not carry a source until it is defined how they contribute to the saturation.
  - This does not affect the bending of light by matter (γ).
  - It does affect whether light itself gravitates.
  - It affects the active mass of hot and pressurised sources. F1.7 N: the non-linear field contributes to the
    exterior mass (1.02–1.51 × the linear value).
- K.15 has an instrument of the carried flow T^ab. The coupling T^ab → ρ_s, however, is not defined.
- κ ↔ G is a calibration (κρ_s ↔ 4πGρ), not a derivation.
- **Dynamics and waves are missing.** The equation is static so far.
  - Scalar theories of gravity usually radiate as a monopole or a dipole. Observations of binary pulsars (the decay of
    the period agrees with the quadrupole radiation of general relativity) are therefore the main future test of the
    equation.
  - Without a dynamical formulation nothing can be said about them.

## 6. What is derived, what tested and what is an input

| | item | state |
|---|---|---|
| **input** | the principle of conserved throughput ρ_c · c = const | book 26.1 as a condition of the Einstein regime; why the substrate keeps it is open |
| **input** | clocks Π, action S = Σ Π Δτ (A6) | axiom |
| **input** | the BD operator (□ − R/2 in the continuum) | taken from the theory of causal sets; K.11, K.13a |
| **input** | f = ln Π | selected by the perihelion, so β = 1 is an input |
| **input** | ρ_s = ρ★ S_v (saturated events), κ ↔ G | decision 6; calibration |
| **derived** | H1 from the principle and the clocks; γ = 2/(d_s − 1), γ = 1 in 3+1 | `why_gamma1.md` |
| **derived** | R and the continuum form of the equation in H1 | sympy, identity (3a) |
| **derived** | first integral r²L′ = M e^L(1 − L), solution Ei(1 − L) = Ei(1) + eM/r | 3b, verified numerically |
| **derived** | a = 0 (β = 1 for the given choice ln Π), 2PN coefficients | 3c, 3d |
| **derived** | no horizon, photon sphere 3.070 M, b_crit 5.2435 M, throat 2.341 M | section 4 |
| **tested** | instrument of γ from the order: 2+1 (K.16b) and 3+1 (K.16c) | PASSED |
| **tested** | the narrowing of the cones alone gives γ = 1/d_s | F7.1a 2+1: 0.47 / 0.43 (FAILED as predicted); F7.1c 3+1: 0.349 / 0.328 and 0.318 / 0.310 (prediction 1/3) |
| **tested** | throughput gives γ = 2 (2+1) | F7.1b 2× UNDECIDED: 1.81 / 1.94 and 1.77 / 1.96; a systematic shift of A downwards and of B upwards |
| **tested** | throughput gives γ = 1 (3+1) | F7.1c 2× UNDECIDED: 1.000 / 0.972 and 0.948 / 0.979; path A passes in both runs, path B fails because of the noise of the comparison with the control (SE ~0.1) |
| **tested** | the BD operator in the H1 geometry: cancellation of the linear term for f = Π | K.13a |
| **tested** | discrete response outside the source c ≈ 0 (β = 1) | F1.9a INSTRUMENT_FAILED (quadrature σ = 4); descriptively c agrees with the continuum value and the deviation falls as (ℓ_k/σ)² |
| **tested** | bound on static strong-field solutions | F1.7 N (continuum, Gaussian source) |
| **tested** | coupling to T_ab: counter-streaming flows, anisotropy and dragging of the cones (decision 13) | F7.4a PASSED: R = 2.401 ± 0.053 against 2.361 (energy 1.5625), K = 1.221 ± 0.012 against 1.192 |

Summary: γ = 1 has moved from an ansatz to a principle, which the laboratory confirms on the raw values in 2+1 and
3+1. A formal verdict is missing only because of the statistical weakness of the comparison. β = 1 remains an input.
Everything in the strong field is a derived prediction of the continuum limit, without a discrete test.

## 7. What would move the theory forward most

In the order of how much it decides about the truth:

1. **Dynamics.** A time-dependent form of the equation and the question of radiation (monopole and dipole against
   binary pulsars). This is where the scalar formulation can fall. It is the most important unresolved point.
2. **The side of the source.** How light, pressure and motion contribute to ρ_s (through the T^ab of K.15). Without
   this, compact bodies and cosmology cannot be compared.
3. **Why the substrate keeps the throughput.** A derivation of the principle from the axioms instead of taking it
   from book 26.1.
4. **A formal verdict of F7.1c.** Only if the author wants it: a new test with statistics fixed in advance (more
   sprinklings for the control and for the candidate). The criteria of F7.1c do not change (R132).
