# The field of the order and the throughput: from φ to resistance

Written on 10 Nov 2026 after F3.14a and F5.2a2, at the author's request. It incorporates his four remarks of 15:59.
**Status: a working note, not a result.** Every statement is marked:
- **tested:** a laboratory result, with the test named;
- **derived:** follows from the text of the theory;
- **choice:** a decision that is not derived;
- **open:** neither tested nor derived.

Nothing here changes the axioms until the author decides it.

## 0. Where we are

| statement | status |
|---|---|
| A defect of a field φ on all events of a 1+1 order, evolved by the retarded Green function of the order, has inertia, with M_in = M_rest within 2 % | tested: F3.14a (seed 0 PASSED; replication failed only on the drift tolerance; physical conditions met twice) |
| In a 1+1 order whose only change is the number of events, the same defect under the same rule falls as its continuum body, and masses 4, 8, 16 fall alike within ~2 % | tested: F5.2a2 (seed 0 PASSED; replication running) |
| The defect compresses the order in proportion to its energy (active mass) | open (F5.2b) |

## 1. What φ is for

In the tests φ was an input: one number per event. The dynamics computes it from the causal past:

  φ(x) = φ₀(x) − (1/2ρ) Σ_{y ≺ x} V′(φ(y))   (1+1).

It did four things that nothing in the laboratory did without it.

| role | without φ | with φ |
|---|---|---|
| inertia | a flow carried its four-velocity (a coordinate rule, F3.6) or spread out like radiation (F5.1g) | conservative propagation of a state (F3.14a) |
| a sharp boundary of matter | a cluster in noise has no boundary | a winding of φ has an integer charge |
| energy and momentum | an axiom (flows carry u^a) | in the continuum limit, quantities of the field (section 3); on the order, not yet defined (section 5) |
| an identity under motion | none | charge plus profile in the defect's own frame (F3.14a, section 7 of the card) |

What φ **is** physically, the tests do not say. They test the dynamical properties of a field, not its meaning.

## 2. φ as the phase of the clocks: a working ontological hypothesis

**Hypothesis (choice, not tested).**
- Π_v is the rate of the clocks of event v (Axiom of local response).
- φ_v is their **phase**.
- The vacuum is the state in which all events are in phase with the substrate.
- A kink is a topological slip of the phase by 2π, a "dislocation in time" that cannot be undone locally.

In this reading matter and gravity are both properties of the clocks of the order: matter as their phase, gravity as
their rate.

**The relative phase (the author's remark 1).** The potential V(φ) = 1 − cos φ depends on the absolute phase and
singles out φ = 0 without a reason. Only the phase relative to the substrate should be physical:

  ψ_v = φ_v − θ_v,   V(ψ) = 1 − cos ψ,

where θ_v is the reference phase of the order at v.
- **θ_v = const (what the tests used, θ = 0):** the choice is a convention, the origin of the phase. The model has no
  physical preference. The vacuum is any θ + 2πn.
- **θ_v varies with the order:** for example, θ_v advances with the proper time of the substrate. Then the phase is
  coupled to the clock rate Π, and gravity would act on matter also through ψ, not only through the counting.
  - θ_v must be a Lorentz-invariant function of the order. A phase that advances with proper time needs a reference
    event or a reference rhythm, so it is not one automatically. **Open.**

**Statement for the text (the author's formulation):** the interpretation of φ as the phase of the clocks of an event
is a working ontological hypothesis. Experiments such as F3.14a can test the dynamical properties of the field; they
cannot by themselves confirm its meaning as a clock.

## 3. Three currents that must not be identified (the author's remark 2)

| quantity | definition (continuum, 1+1) | what it carries | conserved because |
|---|---|---|---|
| **j^μ_top**, topological current | j^μ = (1/2π) ε^{μν} ∂_ν φ | the topological charge; its spatial integral is the number of kinks | ∂_μ j^μ = 0 identically (symmetry of the second derivatives), whatever the dynamics |
| **T^{μν}**, energy and momentum of the field | T⁰⁰ = ½φ_t² + ½φ_x² + V(φ), T⁰¹ = −φ_t φ_x (usual sign convention) | energy and momentum | ∂_μ T^{μν} = 0 by the field equation and the symmetry of translations |
| **F**, causal flow | Axiom of causal continuity: Σ_in F + s_v = Σ_out F + ΔQ_v | "influence" along the relations of the order | by the axiom |

**Hypothesis (open):** F can be derived from the dynamics of φ. It is **not** shown that F equals any of the two
currents above. The structure of conservation is similar in all three, but that is not an argument.
- j_top would make F a count of defects;
- T^{0μ} would make F energy;
- F may also be neither.

What decides it: F must satisfy the axiom at every event of the order, as a discrete statement, while j_top and T are
so far quantities of the continuum limit.

## 4. The state of the conservative dynamics (the author's remark 4)

A second-order dynamics in the continuum needs the data φ|_Σ₀ and ∂_tφ|_Σ₀. The value φ(x) at one event is not a
complete state.

**What the construction of F3.14a does (derived from its definition).**
- The rule is the integral form of the field equation (Duhamel, or Yang–Feldman):
  φ(x) = φ₀(x) − (1/2ρ) Σ_{y ≺ x, y after Σ₀} V′(φ(y)).
- **The initial data enter in full through φ₀.** In 1+1, φ₀ = ½[g(x − t) + g(x + t)] + ½∫_{x−t}^{x+t} h is the free
  wave of the data g = φ|_Σ₀ and h = ∂_tφ|_Σ₀. The data on Σ₀ are given as continuum functions at the coordinates of the
  events (the launch; coordinates declared).
- **After Σ₀, the state is the whole causal past, not the value at an event.** φ(x) depends on all values of φ in the
  past of x, weighted by the Green function. The information that a time derivative carries in the continuum is held in
  the **pattern** of past values: two events with the same φ but different pasts continue differently.
- Consistency with the continuum is tested: dispersion within 0.02 % (I2), velocity, contraction and response to force
  of the kink (F3.14a), and a deterministic quadrature of the same rule agreeing with the continuum solver (F5.2a,
  diagnosis).

**What is open.**
- **Intrinsic initial data.** Σ₀ and the functions g, h are given in coordinates. On a causal set the data would be the
  values of φ on a region of finite thickness (a thickened antichain, as in the initial-value problem of causal set
  field theory) rather than a slice. How a pair (φ, ∂_tφ) is encoded in such a region, and how the free wave φ₀ is
  computed from it intrinsically, is not constructed.
- **An eternal order.** Without an initial slice the massless 1+1 sum over the infinite past diverges (an infrared
  property of 1+1). A massive propagator or the higher dimension changes this. Not constructed.

## 5. The load of an event (the author's remark 3)

**Proposed form:**

  σ_v = 𝓔_v[φ] / C_v,

where 𝓔_v is a **not yet defined** local energy contribution of the field at event v, defined directly on the order.
The continuum density T⁰⁰ is a quantity of the continuum limit. It is not yet a discrete object.

**What the definition of 𝓔_v must decide:**
1. **What is counted.** Three different quantities are possible:
   - (a) the energy stored in the configuration of the field at v;
   - (b) the energy carried through v;
   - (c) the contribution of v to the total energy of the field.

   In the continuum, (a) and (c) correspond to the density T⁰⁰ and (b) to the flux T^{0i}. They have different
   dimensions and integrate differently.
2. **In which frame.** An energy needs a time direction. A frame-free order offers two candidates:
   - the rest frame of the defect itself, T_ab u^a u^b with u its four-velocity, which is itself a readout;
   - an invariant scalar of the field, such as the Lagrangian density or the trace of T.

   The choice changes what loads the order.
3. **Units and normalisation.** T⁰⁰ is energy per proper volume, so its contribution per event is T⁰⁰/ρ. The capacity
   C must be in the same units for σ to be dimensionless. The Axiom of finite capacity gives C as a bound on the flow F,
   so this is tied to section 3 (what F is).
4. **Composition.** The contributions must add over events (an extensive quantity), and the sum over the events of a
   region must reproduce the energy of the field in it in the continuum limit.

**A natural discrete candidate (choice, to be examined).** The field already has a discrete dynamics, so it can have a
discrete action, a sum over events of the kinetic term of the Green function (or of the BD operator) and of V. The
contribution of each event to that action is a local, extensive, frame-free quantity. Whether it is the right "load",
and how it relates to T⁰⁰, is open.

**Development result (10 Nov; findings, "Load candidates"; decides nothing).** In 1+1 the trace of T^{μν} of a scalar
field is exactly 2V(φ), because the kinetic term is conformal in two dimensions. The candidate

  𝓔_v = 2V(ψ_v)/ρ   (ψ = φ − θ, the relative phase)

is local, frame-free, periodic in the phase and zero in the vacuum. Summed over the world tube of a kink it gives M·Δτ,
measured within 1–2 % at rest, at v = 0.3 and 0.6, for M = 8 and 16, and in the nonuniform order of F5.2a2. This is the
action of a flow, S = −Σ m Δτ: **the load of the order summed along the history of matter is its action.** On radiation
it gives the trace, not the energy.

The BD-based estimators of the kinetic term failed per event: the noise of probe P-a again.

**Two hints for later, both interpretations:**
- In the clock picture the load of an event is its mismatch of phase with the substrate, 1 − cos ψ.
- The cosine potential is bounded: the load is largest at anti-phase (ψ = π). A bounded load is what the Axiom of
  finite capacity asks for. This is a hint for H5-D, not a derivation.

**Limits:**
- the identity trace = 2V is special to 1+1; in 3+1 the trace contains −(∂φ)², and a discrete gradient estimator
  better than the BD sums is needed;
- the trace is the scalar source; matter in motion and pressure need the tensor.

## 6. From load to response: what the theory already says

The author's reading: gravity is the throughput limited by capacity, and the field shows itself as resistance.

The text of the theory (ontology "Load", "Throughput", "The field"; Axioms 5–9; decision 10) gives this chain:
1. **Load** (Axiom of finite capacity): the flow through an event is bounded by C.
2. **Response** (Axiom of local response): every event has a clock rate Π_v ≤ 1, a stiffness I_v ≥ 1 and an impedance
   Z_v = I_v/Π_v. The law of response is not part of the axiom (decision 10).
3. **Throughput** (derived from the action of the order at the linear order, Part VI): where clocks slow down and cones
   narrow, events are added so that ρ_c c_eff = const. In 3+1 this is the linearised coupling of general relativity
   to T_ab, and γ = 1.
4. **Resistance:** around matter the order offers a larger impedance to passing influence.

**What the laboratory has already said about the step from load to response** (decision 10):
- **R1, response from the load alone** (Π_v = f(σ_v)): gives a field **without range**. It heats a bundle and does not
  bind it (F3.7). (tested)
- **R2, response from a field equation with a source** (B_g ln Π = κρ_s, or the linear response of the BD action): gives
  a field **with range**, which attracts (F3.8, F3.11). (tested)
- **A5′:** R1 and R2 are compatible if the load acts only in saturation. (choice)

**Consequence for the chain with φ.** The capacity bound alone does not produce gravity at a distance. The long-range
field comes from the action of the order (R2), sourced by the energy of matter. The capacity enters in saturation,
where the field is strong. The chain is therefore:

  φ → 𝓔_v (section 5) → source of the response of the order (R2, the action of the order) → Π, cones, density
  (throughput) → the impedance that a probe feels (resistance); capacity in saturation (A5′).

## 7. The field as resistance, and the three masses

- **Inertial mass** (F3.14a, tested in 1+1): M_in = E of the configuration of φ.
- **Passive gravitational mass** (F5.2a2, tested in 1+1, seed 0): defects fall alike whatever their energy.
- **Active gravitational mass** (open): if the source of the response is 𝓔_v and its sum over the defect is E, the
  defect compresses the order in proportion to E. All three masses would then be the same number, from one rule and one
  field.

**The open link of the axioms, now concrete.** The Causal Resistance Functional of the first version,
R = Σ Z_v Δτ, equals the action of a flow, S = −Σ m Δτ, only if Z ∝ m (Axiom of the action of a flow, status). With φ
the question becomes measurable: does the response of the order to a defect of energy E produce an impedance whose
integral along a probe's history is proportional to E?

## 8. How the order adds events causally (open)

The response must be causal: a new event may depend only on its past.
- **Self-consistent iteration over a whole region** (F6.3b): the present depends on the future. It is a **proxy**
  (axiom audit). It is excluded as a dynamics.
- **Sequential growth:** the density and the cones of the next region follow from the load in its past. A growth law
  of the order is open (Axiom of the dynamics of the order). Such a rule must not define layers by coordinate time, or
  it becomes frame-dependent.
- **Sum over orders weighted by S_BD + S_φ:** the candidate of the Axiom of the dynamics of the order. In the classical
  limit the stationary order is the linear response of general relativity to T_ab of φ (derived at the linear order,
  Part VI).
- **An intermediate test without a growth law:** prescribe the response that the theory derives (the linear response
  to T_ab of the kink), sprinkle the order with it as in F7.1, and check whether a second defect falls towards the
  first in proportion to the first one's energy. This tests the chain of sections 5–7 without solving growth. It is not
  a dynamics of the order.

## 9. Dimension

In 1+1 the active mass cannot be tested.
- The Einstein tensor vanishes identically.
- The throughput formula γ = 2/(d_s − 1) diverges.
- The conformal factor does not change the density per proper volume in the way the principle requires.

Section 8 is therefore meaningful from 2+1 on (F7.1b: γ ≈ 2 in 2+1). The author postponed the change of dimension;
sections 2–7 can be worked out now.

## 10. Decisions for the author

1. **φ as the phase of the clocks:** adopt as a working ontological hypothesis (the formulation of section 2)? And θ:
   a constant (a convention) or a function of the order (to be constructed)?
2. **The three currents:** keep j_top, T^{μν} and F distinct in the text, with the hypothesis "F is derivable from the
   dynamics of φ" stated as open?
3. **𝓔_v:** which of (a) stored, (b) carried, (c) contribution; in which frame; and whether to examine the candidate of
   section 5 (the contribution of an event to a discrete action of the field).
4. **The response: DECIDED (10 Nov, decision 15).** The chain of section 6 holds: R2 gives the range, and the capacity
   acts in saturation by A5′. The author: "F3.7 showed clearly which way not to go."
5. **The initial data** (section 4): an intrinsic construction of initial data on a region of the order as a separate
   task, before the field enters the axioms.

## 11. Summary table

| step | status |
|---|---|
| φ on all events, propagated by the Green function of the order | tested (F3.14a, F5.2a2); the field is an input |
| φ as the phase of the clocks | choice (working hypothesis) |
| relative phase ψ = φ − θ | choice; θ open |
| j_top, T^{μν}, F distinct; F derivable from φ | the distinction derived; the derivation of F open |
| the state of the dynamics = the causal past + the initial data in φ₀ | derived from the construction; intrinsic initial data open |
| 𝓔_v, the discrete energy contribution | open |
| response with range from the action of the order (R2); capacity in saturation | tested (F3.7, F3.8, F3.11) and derived (Part VI); A5′ is a choice |
| active mass = E | open (needs sections 5 and 8, from 2+1 on) |
| R = Σ Z Δτ against S = −Σ m Δτ | open; now a measurable question |
