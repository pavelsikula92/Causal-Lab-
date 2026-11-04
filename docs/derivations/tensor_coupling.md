# What the field carries: a scalar of saturation, or the tensor of the flow

Written on 27 Oct 2026 (fundamental research, at the author's request). It is an analytic analysis, not a result of
the laboratory. The computations were recomputed in the container (sympy). It follows `gravity_equations.md`,
`why_gamma1.md` and `jacobson_route.md` (section 6).

Conventions: signature (−, +, +, +), c = 1, g = η + h. A static weak field has h_00 = 2U and h_ij = 2γU δ_ij, U > 0 is
the Newtonian potential (U = −L). The trace-reversed field is h̄_ab = h_ab − ½ η_ab h. u is the 4-velocity of the
source.

## 1. The question

The equation B_g ln Π = κρ_s is scalar. It has one field Π, one source ρ_s and the H1 metric written in one frame. Is
it a complete dynamical law of gravity, or only a static section of something larger? Three observations decide it,
which a scalar in a fixed frame cannot meet:

- **Frame dragging** (g_0i from moving or rotating matter). Gravity Probe B measured −37.2 ± 7.2 mas/yr against −39.2
  in general relativity. A scalar Π in a fixed frame has g_0i = 0. The laboratory has already shown this linearly:
  `jacobson_linear_check.py`, the condition does not hold for a moving source.
- **The period decay of binary pulsars.** The Hulse–Taylor pulsar agrees with the quadrupole formula of general
  relativity to ~0.2 %, the double pulsar J0737−3039 to ~0.01 %.
- **The polarisation of gravitational waves.** Measurements with three detectors (GW170814 and others) strongly prefer
  purely tensor polarisations over purely scalar ones.

## 2. What a scalar would give as a complete law

A scalar field normalised so that it gives Newton statically (ΔL = 4πGρ) radiates as a quadrupole with the power

    P_scalar = (G/60) [ (tr I⃛)² + 2 I⃛_ij I⃛_ij ]

I_ij = ∫ ρ x_i x_j. General relativity has P = (G/5) Q⃛_ij Q⃛_ij with a traceless Q. For a circular binary
P_scalar / P_GR = **1/6** (sympy, exactly).

To this comes monopole and dipole radiation, depending on how ρ_s depends on the velocity and on the binding energy of
the bodies. For compact bodies (neutron stars) the dipole is usually even stronger than the quadrupole.

Caveat: the equation of the laboratory has no fixed action of the field, so the normalisation of the energy flux is
taken from the standard scalar field. A factor of order 1/6, however, cannot be corrected by the normalisation that
gives Newton. Only another form of the field can correct it.

**Conclusion:** the scalar equation as a complete dynamical law is excluded by observation. It holds as a static
section in the rest frame of the source.

## 3. The main result: γ is the ratio of the tensor and the trace coupling

A static weak field with the parameter γ in the rest frame of the source, d_s spatial dimensions:

    h̄_00 = U (1 + γ d_s),     h̄_ij = U (1 + 2γ − γ d_s) δ_ij

In 3+1 (d_s = 3):

    h̄_ab = U [ (1 − γ) η_ab + 2(1 + γ) u_a u_b ]

- **γ = −1 (fixed order, A6 on a sprinkling):** h̄_ab = 2U η_ab. The field couples only to the trace of the
  energy–momentum tensor (spin 0, Nordström's theory). Light does not bend. Light itself does not gravitate (the
  electromagnetic field has zero trace).
- **γ = 1 (principle of throughput):** h̄_ab = 4U u_a u_b. The field couples to the whole T_ab of dust without a trace
  part. This is exactly linearised general relativity (□h̄_ab = −16πG T_ab, T = ρ u u).
- **γ = 1/3 (narrower cones only):** a mixture of the two.

The principle of throughput √−g · c = 1 gives γ = 2/(d_s − 1). The condition of a pure tensor coupling (h̄_ij = 0 for
dust) gives γ = 1/(d_s − 2). The two agree **only for d_s = 3**.

Three consequences:

1. **The result of F7.1c has a physical content that goes beyond the reading of γ.** A substrate that keeps its
   throughput couples in 3+1 at linear order to the source in the same way as general relativity: to T_ab, not to its
   trace.
   - The principle of throughput is the same as h̄_ij = 0: the deformation of the order lies only along the stream
     (u_a u_b) and across the stream it is zero.
2. **F7.1b (2+1) tests the principle of throughput, not general relativity.** General relativity in 2+1 has no
   Newtonian attraction (jacobson_route, section 6). That the principle and general relativity meet exactly in 3+1 is
   either a coincidence or a hint of why the world has three spatial dimensions. It is not derived.
3. **The question of the side of the source is answered at linear order.** The source is not the scalar ρ_s but the
   T_ab of the flow. It follows that:
   - light gravitates;
   - pressure gravitates as ρ + 3p (the active mass of hot and pressurised matter);
   - the motion of a source creates g_0i (frame dragging).

   In this reading saturated events are T_ab u^a u^b in the rest frame of the flow. A localised flow carries its
   4-velocity as a state (F3.6), K.15 measures T^ab, and the merge rule F = |P| conserves the 4-momentum exactly. The
   conservation ∂_a T^ab = 0 is exactly the condition without which a tensor coupling is not consistent.

## 4. The tensor reading at linear order: what it predicts

The field is a sum over the flows, h̄_ab(x) = Σ 4 U_s u_a u_b (retarded), i.e. the tensor generalisation of the
Σ F/(2πτ) of F3.8. The predictions are then those of linearised general relativity:

- frame dragging with the value of general relativity;
- the quadrupole radiation formula exactly (when ∂T = 0);
- two transverse polarisations + and × propagating at c;
- the bending of light with γ = 1 and the gravity of light and of pressure.

The book (36.2) says the same in words: in the equilibrium regime (plateau) the dynamics of the impedance reproduces
the tensor wave forms of general relativity.

## 5. The non-linear order: a fork

The so-called bootstrap (Gupta 1954, Kraichnan 1955, Feynman, Weinberg 1965, Deser 1970) states: a massless spin-2
field that couples consistently to the total conserved T_ab including the energy of the field itself leads uniquely
to general relativity. The assumptions are Lorentz invariance, locality and at most two derivatives. In 4D this is
supported also by Lovelock's theorem.

General relativity, however, does not satisfy the principle of throughput in the exact form A·B = 1. Schwarzschild has
A·B = (1 − M²/4r²)², so the deviation is at second order. Hence a fork:

**A. Throughput is a plateau law (linear order) and the non-linear completion is general relativity.**
- This corresponds to the Jacobson route (decision 9) and to book 36.2.
- The scalar equation B_g ln Π = κρ_s stays correct to 1PN (β = 1 agrees).
- Its strong-field and 2PN predictions (throat at 2.341 M, no horizon, shadow +0.91 %, g_rr with 2x²) would, however,
  not be predictions of the theory. They would be artefacts of extending a weak-field principle beyond its validity.
- The theory would then have to derive general relativity (Clausius on a causal set), not choose it.

**B. Throughput holds exactly, also non-linearly.**
- The theory then differs from general relativity from 2PN on and must violate some assumption of the bootstrap: the
  non-locality of the BD operator at the scale ℓ, or a preferred frame.
- At the same time it must show that radiation and frame dragging nevertheless agree with observation.
- It is a markedly riskier route and today there is no mechanism for it.

**The author chose A (decision 13, 28 Oct).** The laboratory recommended the same: A rests on finished results (K.14,
K.15, F3.6, the linear Jacobson check, F7.1c), does not collide with the uniqueness theorems and agrees with the book.
The strong-field section of `gravity_equations.md` is therefore not valid as a prediction of the theory.

## 6. Further fundamental steps (proposals, nothing running)

1. **A discrete test of the tensor coupling (F7.4a, the author called it F7.2): do counter-streaming flows gravitate
   according to T_ab?** **PASSED 28 Oct:** R = 2.401 ± 0.053 against 2.361 (tensor) and 1.5625 (energy), anisotropy
   1.221 ± 0.012 against 1.192; a single stream tilts the cone (dragging). The rule "own frame of the flow + addition
   of the capacity deficits" gives a coupling to T_ab on a causal set. Details in `docs/findings.md`.
   - Two counter-streaming bundles (configuration of K.15) with velocities ±v have, in the rest frame of the sum, the
     active mass T_00 + T_ii.
   - Three candidates give different predictions:
     - tensor coupling: 2ργ²(1 + v²);
     - rest density of saturated events (a Lorentz scalar): 2ρ;
     - density in the laboratory frame: 2ργ.
   - The instrument is K.16c (φ from the clocks).
   - Assumption: the substrate rule must decide how saturations with different u compose. Exactly this is where the
     candidates differ, and that is the fundamental question.
2. **The Jacobson route, step 4 non-linearly:** the Clausius relation on a causal set with a flux across a horizon. It
   decides between A and B: if the substrate keeps the Clausius relation, general relativity and option A follow.
3. **Why the substrate keeps the throughput:** in the form h̄_ij = 0 it is the question why saturation deforms the
   order only along the stream. A candidate is the conservation of the flow (axiom IV): a flow without transverse
   stress (dust) does not transfer momentum across.
