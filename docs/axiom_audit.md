# Audit of the axioms and rules

Written on 22 Oct 2026 at the author's request ("this is the priority, the other tests are secondary"). Sources:
*Physics of Causality* (Part IV, Axioms; chapters 5, 22, 23, 26), the docstrings of the core and of the tests of the
laboratory, `findings.md`, `author_decisions.md`. The audit neither rejects nor approves anything. For every rule it
determines what it rests on, and marks what is an **input** of the theory and what an **output**.

## Criteria

Every rule is judged by five questions. A rule that does not pass (a)–(d) is a **proxy**: an aid that may show a
direction, but whose results are not physics.

| | Question | Why |
|---|---|---|
| (a) | Does it use only the order and the numbers of elements (possibly a state of events derived from them), not coordinates? | Coordinates are not available in the physics of the substrate; order and number are ("order + number = geometry"). |
| (b) | Is it causal: does the state of an event depend only on its past? | A rule that needs the future is not a dynamics. |
| (c) | Is it independent of the reference frame and of the order of computation? | Lorentz invariance is verified to tens of orders of magnitude; the order of computation is not physical. |
| (d) | Is it finite and well defined without the wall of the region and without a regulator? | A result that holds only thanks to the edge of the simulation is a property of the simulation. |
| (e) | Is it an input (postulate, tuned parameter) or an output (it followed without being put in)? | A test of an output that was put in as an input proves nothing (tautology). |

## 0. The choice of the substrate: made in the tests, but stated nowhere

*Correction after the author's remark (22 Oct): the laboratory has already excluded the bound K₀. The first version
of this section overlooked it.*

This is the most important finding of the audit. **The text of the theory and the laboratory rest on two different
substrates.**

**The substrate of the text of the theory** (*Physics of Causality*, Part IV and chapter 26):
- a finite degree deg(v) ≤ K₀ (section 15);
- locality by path length in a graph, dimension from the number of events in a ball of graph distance N(r) ∝ r^D
  (sections 15, 22);
- time as a number of updates τ (section 20.1);
- special relativity derived from clocks "at rest relative to the substrate", the Lorentz transformation as a
  "reparametrisation of coordinates, not a fundamental symmetry" (sections 26.2, 26.5).

Such a substrate has a **rest frame** (a position close to Lorentz's ether theory).

**The substrate of the laboratory** is a Poisson sprinkling of Minkowski spacetime, i.e. a causal set. It has no rest
frame, it is statistically Lorentz invariant, and every event has infinitely many links along the light cone. They
are finite only thanks to the wall of the region.

Both cannot be had at once: a Lorentz-invariant discrete structure cannot select a finite set of neighbours locally
(Bombelli, Henson, Sorkin, "Discreteness without symmetry breaking: a theorem", 2006/2009). A finite degree and the
invariance of the substrate exclude each other. The options are:

- **A. A Lorentz-invariant substrate (causal set):** relativity for free, non-locality unavoidable, a capacity may
  limit the flow, not the number of links.
- **B. A substrate with a rest frame (a graph of finite degree, updates):** Lorentz invariance must be derived for
  the observers from the substrate with extreme precision. A small violation at the Planck scale is, moreover,
  amplified quantum mechanically (Collins, Perez, Sudarsky, Urrutia, Vucetich, PRL 93, 191301, 2004).

**The laboratory has in fact already ended up in A:**
- F0.3 (E42b): the degree of an event grows with N, the bound K₀ does not hold (listed under "Does not hold");
- F1.8 (E40, E41-A): the capacity limits the flow, not the number of links (supported, not robust against the seed);
- F6.1 and F6.2: local growth rules of type B did not produce a manifold (0 of 7) nor a three-dimensional space.

The tests did not support B in its simplest forms, and the laboratory moved to sprinkling. That is choice A.

Two things are missing:
1. **To state A as the author's decision** and rewrite the text of the theory. The remnants of B are K₀ in the
   axioms, locality and dimension from graph distance, time as a number of updates and the derivation of the Lorentz
   factor from the rest frame of the substrate (section 26.5). In A this part is unnecessary, because invariance is a
   property of the substrate.
2. **To get rid of the remnants of B in the rules of the laboratory.** Narrowed cones in coordinates and "the most
   timelike descendant" are intuitions of a rest frame carried over into A. The quick runs of F0.6 show exactly this
   conflict. The uniform division of the flow among the descendants is an intuition of a finite degree: in A an event
   has infinitely many descendants.

A warning for A: a growth dynamics is hard there too. Classical sequential growth (Rideout–Sorkin) generally produces
orders that do not resemble manifolds. Finding a growth that gives a manifold is an open problem of the whole field.
F6.1 is its local analogue.

## 1. The axioms of the theory (*Physics of Causality*, Part IV)

| Axiom / notion | (a) | (b) | (c) | (d) | (e) | Note |
|---|---|---|---|---|---|---|
| I. Primacy of events (events without properties, everything relational) | ✓ | ✓ | ✓ | ✓ | input | a clean foundation, compatible with A and B |
| II. Partial order, no loops | ✓ | ✓ | ✓ | ✓ | input | the same in the theory of causal sets |
| Flow Φ(v) = Σ_u w(u,v) Φ(u) over the past | ✓ if w from the order | ✓ | ? | ✗ in A | input | the recursion is causal; the weight w(u,v) is not fixed; in A, v has infinitely many parents and children |
| Impulse ρ, tension χ_ρ = ln(ρ/ρ_vac), χ_K = ln(ρ/K₀) | ? | ? | ? | ? | input | ρ "drives the rewiring", but an operational definition is missing; the laboratory does not implement it |
| III. Local update operator Û (ontogenesis) | ✓ | ✓ | **?** | ✓ | input | it is growth, which fits; but if "time = sequence of updates" is global, (c) is violated. A necessary condition of independence of the order of growth (general covariance in the sense of Rideout–Sorkin) |
| IV. Causal continuity (balance of the flow at every event) | ✓ | ✓ with a reservoir | ✓ | ✓ | input | normalising the outgoing weights needs the number of children (the future); the reservoir ΔΦ solves it (decision 2) |
| Finite capacity deg(v) ≤ K₀ | ✓ | ✓ | **✗ in A** | ✓ | input | **refuted in the laboratory** (F0.3: the degree grows with N); the capacity limits the flow (F1.8); a remnant of B in the text of the theory (section 0) |
| Locality = path length in a graph | ✓ | ✓ | **✗ in A** | **✗ in A** | input | in a causal set every event is one link away from events arbitrarily far along the light cone |
| Π, I, Z = I/Π | ? | ? | ? | ? | input | only verbal definitions ("throughput", "resistance to rewiring"); operationally so far only Π (laboratory) and the planned I (F0.5) |
| D_eff from N(r) ∝ r^D (graph balls) | ✓ | ✓ | **✗ in A** | **✗ in A** | output (claim) | in A the graph balls are infinite; in causal sets covariant estimates of dimension exist (Myrheim–Meyer, from the number of relations in an interval) |
| Gravity: a ∝ −Φ(r), Φ ∝ 1/r² from the flow through shells (section 23) | ✗ | ? | ✗ | ? | **input** | "minimal assumption": the proportionality of the acceleration to the flow is a postulate and the shells r are coordinate shells; Gauss's law is correct, the link to the acceleration is put in |
| Motion: minimum of R = Σ Z Δτ (section 5.8) | ✓ if Z from the order | ✓ | ? | ? | input | momentum is missing ("persistence of motion is not conservation of momentum"), so inertia has no carrier |
| SR: γ = 1/√(1 − v²/c²) by a "linear" division of the throughput (section 26.5) | ✗ | ✓ | **✗** | ✓ | **input** | presupposes a rest frame of the substrate; a "linear division" would give 1 − v/c, the square root is put in (a Pythagorean division) |

## 2. The rules of the laboratory

| Rule (where) | (a) | (b) | (c) | (d) | (e) | State |
|---|---|---|---|---|---|---|
| Sprinkling, causal relation, links (`core.order`) | ✓ | ✓ | ✓ | ✗ (wall) | input | substrate A; links exist only thanks to the wall (F6.3 links at the wall) |
| The BD operator and its expectation value (K.4, K.11, K.13a) | ✓ | ✓ | ✓ | ✓ (smearing ε) | input | **firm ground**: covariant and finite, the non-locality tamed by the construction |
| Uniform division of the flow among the descendants (`route_capacity`, F3) | ✓ | **✗** | ✓ | **✗** | input | needs the future number of children; in an infinite spacetime there are infinitely many → **proxy** |
| Capacity C = Q₀.₉₉ of the flow under unit driving | ✓ | ✓ | ✓ | ✗ | parameter | the calibration depends on the density and the region |
| Response Π = f(σ) (axiom 5; from 27 Oct the candidate law R1, decision 10) | ✓ | ✓ | ✓ | ✓ | input | local; the shape of f, Π_sat, m are working parameters (decision 5); axiom 5 is from 27 Oct only a general response |
| Copied pulse, threshold ΔQ_crit, capacity channel (F3.3c–h) | ✓ | ✓ | ✓ | inherits the division of the flow | input + parameters | the pulse as a source contradicts decision 4 |
| Narrowed cones \|Δx\| ≤ Π²Δt, "the most timelike descendant" (F3.3h) | **✗** | ✓ | **✗** | inherits | input (H1) | **proxy**; F0.6 fixed FAILED (23 Oct): in a moving frame extinction and a drift of the centroid +0.34 |
| Local relation of the dynamic order (F6.3b) | **✗** | ✓ | **✗** | ✓ | input (H1) | **proxy** (F0.6 dynamic, quick runs) |
| Self-consistent iteration of the order and the flow over the whole region (F6.3b) | ✓ | **✗** | ✓ | ✗ | a choice of solution | the present depends on the future → **proxy** |
| Action W = Σ Π on maximal chains (F2) | ✓ | ✓ | ✓ | ✓ | input | **firm ground** (chains are covariant); a difference of ~10 % against the geodesic open |
| Geometry H1 (−Π²dt² + Π⁻²dx²) and γ = 1 | – | – | ✓ | ✓ | **input** | γ = 1 is put in, a test of γ is tautological (decision 1) |
| Field equation B_g f(Π) = ρ_s, f = ln Π | ✓ | ? | ✓ | ✓ | input | the static radial solver is an instrument of the continuum, not a dynamics |
| Bridge to SI, G = κc³ℓ_d²/(4πλ₁ħ) | – | – | – | – | **input** | a conditional bridge (decision 8) |
| Readouts (radii, shells, centroids in coordinates) | ✗ | – | – | – | – | **admissible**: an instrument may use coordinates if it declares them; a rule of the dynamics may not |

## 3. Summary

**On firm ground** (passes (a)–(d)): the causal order, the BD operator and its limit □_g − R/2 (K.4, K.11, K.13a),
the action on chains, the local clock response Π = f(σ).

**Proxies** (the results are not physics until a rule that passes replaces them): the uniform division of the flow
among the descendants, narrowed cones in coordinates, the local relation of the dynamic order, the self-consistent
iteration. The whole F3 series and F6.3b rest on them.

**Inputs that are, or could be, presented as outputs:** γ = 1 (H1), a ∝ Φ (gravity), the Lorentz factor (section
26.5), G.

**Outputs:** none non-trivial so far. The closest are the geodesic behaviour of the action in 2+1 (F2) and the
convergence of the BD operator. Both, however, are properties of the known construction of causal sets, not new
predictions of the theory.

**Missing degrees of freedom:**
- **The momentum of matter.** Without it there is neither inertia nor a covariant narrowing of the cones.
- **The second gravitational degree of freedom.** Clocks and space must be independent, otherwise Nordström's scalar
  gravity repeats itself (the laboratory found it again: γ = −1 for Π²η).
- **The quantum side.**

## 4. The author's intuition about the dilution of causality (22 Oct)

A record of the idea: light from far away and causality in time get diluted. A tree next to me has more detail than a
tree on the horizon, and memories and monuments of the same age are subject to entropy.

A candidate principle, not a rule (R132). Physics says three things about it:

1. **Geometric dilution** (the tree on the horizon): the same flow is spread into a larger solid angle, ∝ 1/r². The
   theory already has this (section 23), and it is the same thing that gives Gauss's law.
2. **Entropic dilution** (memories, monuments): loss of detail through many intermediate interactions. In physics
   this is a property of a coarse description. The fundamental dynamics (unitary quantum theory) does not lose
   information.
3. **Sharpness is not lost in vacuum.** The light of distant galaxies and quasars arrives sharp down to its angular
   size. The principle is therefore **testable**:
   - in a substrate **B** with finite channels (the laboratory abandoned it through its tests, section 0) even "pure"
     causality across the vacuum would be diluted, and the images of very distant point sources would blur with the
     causal distance;
   - observations of sharp images of distant quasars and gamma-ray bursts already constrain such models of
     "spacetime foam" (e.g. Perlman et al., ApJ 805, 10, 2015);
   - in a substrate **A** no dilution happens along the light cone (distant events on the cone are causally "close"),
     dilution happens only through the volume of the interval, i.e. through intermediate events.

   The intuition thus itself points to the decision of section 0 and gives it an observable consequence.

## 5. Proposed order

1. State choice A, which the tests have in fact made (section 0), as the author's decision. If the author wanted B,
   it must be explicit and with F0.3, F6.1 and F6.2 in mind.
2. Rewrite the axioms of section 1 that did not pass as rules from order and number (locality, dimension, time,
   capacity of the flow instead of a degree) and drop the derivation of the Lorentz factor from the rest frame.
3. Dynamics as growth of the order, causal and independent of the order of growth. For A the mathematical basis is
   ready: classical sequential growth (Rideout, Sorkin, PRD 61, 024002, 2000), whose conditions (general covariance,
   Bell causality) are identical to the conditions of decision 2.
4. Only then matter (with momentum), the existence gate and the test of the universality of free fall.

The running and finished tests (K.12, K.13a, F1.7, F3.3h, F6.3b, F0.6) stay valid according to their criteria. What
changes is only what may be read from them as physics.
