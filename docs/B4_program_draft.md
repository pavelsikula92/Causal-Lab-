# B4: a localized defect as a candidate for matter (programme draft, 13 Nov 2026, 22:20; NOT a card, decides nothing)

Order of the programme (the author, 13 Nov, 22:17): **C0 → B4-O in 3+1 → identity and motion → B4-S → gravitational
feedback.** Nothing from this plan is added to C0 or, retroactively, to earlier experiments.

## Principle

First a mechanism that explains why a localized configuration survives; only then whether it has the properties of
matter. The main criterion is not the lifetime but **the robustness of the identity under perturbations**. An object
that survives only a perfectly spherical initial profile is not a convincing physical object.

The author's heuristic "an extreme provokes the opposite extreme, which forms a new stability" is made measurable as
an attractor: a periodic or quasi-periodic orbit in configuration space to which the system returns after small
perturbations (Floquet multipliers below 1; different initial profiles converge to the same (E, ω)).

## Two parallel hypotheses

**B4-O, dynamic localization (minimal change of the model).**
- A real scalar field; an oscillating localized configuration (an oscillon).
- Measured: the lifetime, the radiation and the stability against deformations.
- No exactly conserved charge is assumed.

**B4-S, a genuinely stable soliton (a stronger hypothesis).**
- A multi-component field or an additional internal structure.
- A possible topological or Noether invariant.
- Stability against collapse and against radiation.
- The charge must follow from the dynamics of the substrate and be conserved through evolution, collisions and changes
  of density. It is never assigned to make the object survive.

## Measurements (B4-O)

1. **Lifetime is not stability.** Measured together:
   - the total energy;
   - the energy inside the core (in the core's own frame when it moves);
   - the flux of energy across the core's boundary;
   - the frequency spectrum of the radiation;
   - the radius and shape of the core in time.

   These separate a stable object, a metastable oscillon and a lump that slowly radiates and decays.
2. **A family of potentials, not sine-Gordon alone.** The family keeps the minimum, the mass m and the leading nonlinear
   coefficient, and varies only the higher terms. The rule takes J = V′(φ).
   - Why this is decisive: small-amplitude oscillons in 3+1 are unstable (the 3D nonlinear Schrödinger limit
     collapses), so the long-lived ones have large amplitude and their lifetime is set exactly by the higher terms.
   - The literature reports much longer lifetimes for some classes of potentials than for sine-Gordon (Salmi and
     Hindmarsh, Phys. Rev. D).
   - If the oscillon disappears under a small change of the higher terms, it is a special numerical configuration,
     not a robust mechanism.
3. **Perturbations of the identity:**
   - a changed initial profile;
   - an asymmetric deformation;
   - a hit by a weak wave;
   - a collision of two objects;
   - a change of the substrate density;
   - a change of the numerical resolution.
4. **Identity.** The pair (E, ω), energy and internal frequency ω < m, is the invariant tag of the object. Tracking the
   maximum is not enough: after a collision we must be able to say whether it is still the same object.

## Identity and motion

The tests:
- motion without an external force;
- conservation of energy and momentum;
- Lorentz contraction in the corresponding continuum;
- elastic and inelastic collisions;
- transfer of energy into wave modes;
- decay channels.

The velocity must follow from the dynamics of the field (boosted initial data). It is never put into the transport
rule, unlike the earlier F3.13c.
- The internal frequency is the object's clock: a moving oscillon should show ω/γ in the laboratory. This is a test of
  Lorentz invariance of the substrate itself, since the sprinkling has no rest frame.
- The domain (the pyramid) does have a rest frame, so the object must stay far from its boundary.

Properties of matter checked here:

| property | requirement |
|---|---|
| inertia | the object resists a change of its state of motion |
| universality | the same object has the same inertial response however it was formed |
| identity | after a collision or deformation we can tell whether it is the same object |

## Gravitational feedback (last)

Only after the oscillon exists in a homogeneous substrate:
- the feedback through the source ρ + Σp is switched on (decision 24);
- then it is measured whether its presence changes the local clocks, the causal cones and the flow of the
  surroundings.

The red shift of the oscillon's internal clock continues the clock measurements of F5.2d. This is the real test of
"matter is a localized defect of the substrate", not a visualization of a dense region.

## Practical gates (from what C0 has shown)

- **Time and cost.** Oscillon lifetimes in 3+1 are tens to thousands of time units; the substrate reaches T ≈ 6 at
  ρ 30 (700 s per sprinkling). B4-O therefore has two parts:
  - the mechanism (the potential family, long times, the energy balance) is studied in the continuum on a lattice;
  - the substrate is compared in short windows: the first periods, the early radiation rate, the dependence on ρ.
- **Noise.** In the present 3+1 rule the noise accumulates roughly linearly in time, and by t ≈ 5 it is comparable
  with the signal even for the quietest kernel. C0 must deliver the growth exponent of the noise. A computational
  formulation to examine within C0, not a change of the physics: the 3+1 Huygens property (the massless Green function
  lives on the cone) allows the field at an event to be composed from data on a recent slice (the Kirchhoff form)
  instead of the whole past, so the noise would not sum over the whole history.
- **Energy on the substrate.** The rule φ = φ0 − G*J is not derived from an action on the causal set, so the energy is
  read by interpolation to the continuum (coordinates as a diagnostic only). The substrate noise is a channel of its
  own, so the energy balance needs paired sprinklings, with and without the object (as κ against κ = 0 in F5.2d).
