# The Jacobson route: gravity as an equation of state of the causal structure

Written on 22 Oct 2026 at the author's request ("let's try it"). **Status: a plan and an analysis of the assumptions;
nothing of it is yet a result of the laboratory.** The instrument of the first step is K.14.

## 1. Jacobson's argument (1995)

T. Jacobson, *Thermodynamics of Spacetime: The Einstein Equation of State*, PRL 75, 1260 (1995). Units c = k_B = 1.

Assumptions:

1. **The entropy of a horizon is proportional to its area.** Every local causal horizon (the Rindler horizon of an
   accelerated observer around an arbitrary point, in an arbitrary direction) carries the entropy S = η·A with a
   **universal** constant η.
2. **The Clausius relation** δQ = T δS for the heat that flows through the horizon, in local equilibrium.
3. **The Unruh temperature** T = ħκ/2π, where κ is the acceleration of the observer. It comes from quantum field
   theory.
4. **The geometry of null congruences** (the Raychaudhuri equation). The area of a cut of the horizon changes
   according to how the curvature focuses it: δA = −∫ λ R_ab k^a k^b dλ dA.
5. **The energy flux** across the horizon: δQ = ∫ λ T_ab k^a k^b dλ dA (T_ab is the energy–momentum tensor of matter).

Result: requirement 2 for all points and all null directions k requires R_ab k^a k^b = (2π/ħη) T_ab k^a k^b. From
this follow the **Einstein equations**

    R_ab − ½ R g_ab + Λ g_ab = 8πG T_ab,     G = 1/(4ħη),

with the cosmological constant Λ as an integration constant (not fixed).

The Einstein equations are not a fundamental law here. They are equations of state, much as the equation of an ideal
gas follows from the statistics of molecules that we do not know individually.

## 2. What follows for the theory

If the assumptions held for the substrate of Causal Physics, gravity would be **derived**, not put in:

- no field equation that we would choose (B_g f(Π), the choice of f, ln Π against Π − 1);
- no ansatz H1: γ = 1, the Shapiro delay, the perihelion advance and gravitational waves would follow from general
  relativity;
- both halves of gravity (the clocks and space) would arise together, so the Nordström problem (γ = −1 for a scalar)
  does not occur.

**The bridge to SI (decision 8).** If the entropy is the number of horizon molecules times the entropy s₀ of one
molecule, then η = s₀ a^(4)/ℓ² and the condition G = 1/(4ħη) gives

    ℓ² = 4 s₀ a^(4) ℓ_P² = (2√3/5) s₀ ℓ_P²,    ℓ ≈ 1.075 √s₀ ℓ_P

(a^(4) = √3/10, ℓ_P the Planck length). The length of the discreteness is no longer a free parameter; it is tied to G
up to s₀. It becomes a prediction only when ℓ or s₀ is measured independently (e.g. from Sorkin's estimate of the
fluctuations Λ ∝ 1/√V, or from the non-locality of the operator). Until then it is a consistency relation.

**The connection with slow clocks.** In a static gravitational field the Tolman–Ehrenfest relation holds: the local
temperature T_loc = T_∞/N, N = √(−g₀₀). In the metric of the laboratory N = Π. Where the clocks run slower it is
locally warmer. The starting intuition of the theory (gravity = slowing of clocks) has a precise meaning in the
thermodynamic reading: Π is the redshift factor of the temperature.

## 3. Which assumptions the laboratory can supply and which not

| Assumption | What in a causal set | State |
|---|---|---|
| 1. S ∝ A, universal η | the number of horizon molecules (Dou–Sorkin; 1-molecules Barton et al. 2019) | **K.14** (flat and conformally flat space); further curved cones |
| 2. Clausius relation | a requirement on the dynamics of the order in local equilibrium | input (a principle), not yet testable: the laboratory has no dynamics |
| 3. Unruh temperature | a quantum field on a causal set (Sorkin–Johnston) | input from quantum field theory |
| 4. Raychaudhuri | holds in the continuum limit of a sprinkling | taken over |
| 5. T_ab of matter | the flow of the laboratory **and momentum**, which the laboratory does not have (audit, section 3) | **missing**: without momentum matter cannot be put into the equation |

The most sensitive points are two. The universality of η (assumption 1) is what can fail numerically: if the density
of molecules per unit proper area depended on the curvature already at leading order, the argument would fall.
Assumption 5 is what is missing conceptually: the matter of the theory must have an energy–momentum tensor, i.e. also
momentum.

## 4. Steps

1. **K.14:** number of molecules = a^(d) × area in flat space (2, 3, 4 dimensions, two densities) and proportionality
   to the **proper** area in conformally flat space (alternatives: coordinate area, volume density).
2. **Curved cones:** a horizon in a geometry in which the cones change too, not only the density (weak field of
   Schwarzschild / H1). It needs a sprinkling with a non-conformal causal relation. For a start the first order
   suffices, as in K.13a.
3. **Matter with momentum:** define the energy flux across a horizon from the theory's own matter (flow + momentum).
   Here the Jacobson route meets the F3 series.
4. **Clausius on a causal set:** in a geometry with a flux of matter across a horizon, compare the change of the number
   of molecules with δQ/T. For a geometry that already solves the Einstein equations this is automatic if step 1 holds.
   It becomes non-trivial only with the matter of step 3 and with a dynamics of the order.

## 5. Caveats

- Jacobson's argument assumes local equilibrium. Away from it, it gives corrections (viscosity, higher curvature
  terms; Eling, Guedens, Jacobson 2006).
- Identifying the entropy with the number of molecules is a hypothesis. Another option is the entanglement entropy of
  a quantum field on a causal set (Sorkin, Yazdi), which needs a truncation of the spectrum.
- The argument says what gravity must look like in equilibrium. It does not say which microscopic dynamics enforces the
  Clausius relation. The growth of the order (audit, section 5) remains open.
- Λ remains undetermined. That is, however, in keeping with causal sets, where Λ appears as a fluctuation (Sorkin).

## 6. State after 23 Oct 2026

- **Step 1 (K.14):** passed; an instrument, not a result of the theory.
- **Step 3 (matter with momentum):** the condition that the audit listed as missing is met. F3.6 (replicated with
  seeds 0 and 1) showed that a localised flow carries its 4-velocity u as a state. The energy–momentum tensor of the
  flow is therefore definable: T^ab = Σ F u^a u^b per unit volume. Limitation: the model conserves the scalar flow, not
  the four-momentum (an excess at merging of 10⁻⁵ to 5·10⁻⁴ of the flow, F3.10 unit), so ∂_a T^ab = 0 holds only
  approximately.
- **Step 4 at linear order (derived, `jacobson_linear_check.py`):** for local Rindler horizons the Clausius relation is
  equivalent to the condition R_ab k^a k^b = 8πG T_ab k^a k^b for every light-like direction k. The weak field H1
  written in the frame n is h_ab = −2U(η_ab + 2 n_a n_b).
  - n = the frame of the laboratory (a scalar Π in a fixed frame), the source moves: the condition **does not hold**
    (differences of order unity in four light-like directions). A scalar field does not know in which direction the
    matter moves.
  - n = the own 4-velocity of the source u: the condition **holds exactly** (zero in all directions). This is exactly
    linearised general relativity for dust.
  - For many sources with different u the contributions add: the field is tensorial, h̄_ab ∝ Σ F u_a u_b G_ret (a
    generalisation of the retarded field Σ F/(2πτ) of F3.8). Its transverse traceless part is gravitational waves with
    the polarisations + and ×. The cones fixed by the whole metric η + h then deform asymmetrically, not only narrow
    (c_eff = cΠ²).
  - Consequence: Jacobson's condition requires that the field carries the own 4-velocity of the sources, i.e. the
    momentum of matter (F3.6). A scalar Π is only a static special case (the 00 component for a source at rest).
- **Step 2 (a horizon with curved cones):** feasible with the tools of K.13a (optical time, volumes with curved cones).
  In H1 the density is √−g = Π^(1 − d_s γ) and the proper area of a cut of the horizon Π^(−(d_s − 1)γ); the ratio is
  Π^(1 − γ), so 1 at γ = 1. Point density and proper area are not distinguished in H1; the distinction must come from
  the shape of the cones. A deviation from proportionality to the proper area is expected only at order
  (ℓ/L_curvature)²: local flatness.
- **Open:** Clausius at non-linear order and with a dynamics of the order that responds to matter (F6), and the Unruh
  temperature (an input from quantum field theory). In 2+1 general relativity has neither local gravitational degrees
  of freedom nor Newtonian attraction; the attraction in F3.10 (2+1) is therefore a property of the scalar model, not of
  general relativity. A test of waves and of the tensor field must be in 3+1.
