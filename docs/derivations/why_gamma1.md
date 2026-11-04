# Why γ = 1: an analysis of principles and the candidate "conservation of causal throughput"

27 Oct 2026, after F7.1a (at the author's request: "refresh the foundations of physics and the principles of nature
and of society, then an analysis"). Analytic work, no test. The algebra was checked with sympy (outputs in the text).

## 1. What exactly needs to be explained

Static weak field, isotropic coordinates, d_s spatial dimensions:

    ds² = −e^{2φ} dt² + e^{−2ψ} dx²,    γ = ψ/φ.

In a causal set the metric is order + number (the instrument K.16b reads both components blindly):

| quantity of the substrate | weak-field value |
|---|---|
| clocks (number of elements of a diamond at rest) | φ |
| speed of light c (boundary of the relation) | φ + ψ |
| density of elements per coordinate volume ρ_c = √−g | φ − d_s ψ |

What the laboratory knows so far:
- a fixed order (A6 on a sprinkling, the rule of F3): c = 1, ρ_c = 1, so γ = −1 (the Nordström class);
- a capacity that only removes links (F7.1a): ρ_c = 1, c < 1, so γ = 1/d_s (2+1: measured 0.47 / 0.43);
- general relativity needs in 3+1 at the same time c = e^{2φ} and ρ_c = e^{−2φ}.

One principle is sought from which both the narrowing of the cones and the densification of the events follow in
the right ratio, without two independently set functions.

## 2. What the foundations of physics say about γ

- **The equivalence principle alone does not fix γ.** Einstein in 1911 obtained from it only the slowing of clocks
  and half the bending of light (γ = 0). The full bending (γ = 1) requires a curvature of space, i.e. a field
  equation.
- **The type of coupling decides.** Scalar gravity coupled to the trace of the energy–momentum tensor (Nordström) has
  γ = −1. Light has zero trace, so it neither feels nor exerts gravity. A field coupled to the whole energy–momentum
  tensor (spin 2) gives γ = 1: every energy gravitates, including light and moving matter. This is the universality
  of the coupling.
- **In the continuum limit γ = 1 is the same as "the vacuum does not focus light".** The linearised Ricci curvature in
  a light-like direction k = (1, n) is R_kk = ∇²φ + ∇²ψ + nⁱnʲ∂ᵢ∂ⱼ(ψ − φ) (`focusing_gamma.py`). Requiring R_kk = 0 in
  vacuum for all directions n forces ψ − φ to have no traceless Hessian, with the boundary condition ψ = φ. This is
  the focusing principle (`spatial_part.md`), i.e. the null–null part of the Einstein equations. As a postulate it
  only renames γ; the question why the substrate would do it remains.
- **The radial version is known.** In Schwarzschild coordinates g_tt g_rr = −1 exactly when radial light rays carry
  no energy flux (T^t_t = T^r_r; Jacobson, "When is g_tt g_rr = −1?", 2007; quoted from memory, to be checked).

## 3. Principles of nature and of society as a source of candidates

I went through principles that recur in nature and in society and asked what each of them would force for the pair
(c, ρ_c).

| principle | example | what it would force | γ |
|---|---|---|---|
| least action | optics, mechanics | a dynamics, not the ratio c : ρ_c | – |
| equivalence | free fall | only the clocks | 0 |
| scalar coupling (density only) | Nordström | ρ_c, c = 1 | −1 |
| limited capacity (removal only) | F7.1a | c < 1, ρ_c = 1 | 1/d_s |
| Le Chatelier (a system compensates a change) | chemistry, market | a compensation, but not its size | – |
| **conservation of flow through a constriction** | traffic q = ρv, continuity, Kirchhoff, money MV = PT | ρ_c · c = const. | **2/(d_s − 1)** |
| isotropy and universality of the coupling | spin 2 | ψ = φ (continuum) | 1 |

The row with the conservation of flow is the only one that gives a concrete ratio without a free function. The rest
of this analysis is about it.

## 4. Candidate: conservation of causal throughput

Intuition from traffic: in a jam, vehicles slow down and pack together so that the throughput of the road q = ρ·v
stays the same. The continuity equation says the same, as do Kirchhoff's law and the quantity theory of money (when
the velocity of circulation falls, there must be more money for the same number of transactions to flow).

**Principle:** where the substrate slows the propagation of causes (narrower cones), events become denser so that the
number of events that the front of light crosses per unit time and unit area does not change:

    q = ρ_c · c = const.  (= 1 far from matter).

**Consequence (sympy):** ln(ρ_c c) = 2φ − (d_s − 1)ψ, so

    ψ = 2φ / (d_s − 1),     γ = 2 / (d_s − 1).

- **d_s = 3 (our world): γ = 1.**
- d_s = 2: γ = 2. In 2+1 general relativity has no Newtonian limit, so no "correct" γ in 2+1 exists. The principle
  thus gives γ = 1 exactly in the dimension in which gravity exists as in our world.

**An even stronger property:** the H1 metric (−Π²dt² + Π⁻²dx²) satisfies ρ_c · c = 1 in 3+1 **exactly, to all
orders** (sympy). Schwarzschild in isotropic coordinates only to first order: ρ_c c = (1 − M²/4r²)². H1 is therefore
exactly "clocks Π + conservation of causal throughput":

    A = Π² (clocks),  ρ_c · c = 1 (throughput)  ⇒  A·B = 1  ⇒  B = Π⁻², c = Π², ρ_c = Π⁻².

Two independent choices of the ansatz H1 (the temporal and the spatial part) are replaced by one physical condition.
The prediction of the capacity rule of F7.1a thereby turns into a requirement: the narrowing of the cones is not a
free function c(σ) but c = Π², and the density must grow to Π⁻². The relation to F1.9a and β = 1 remains: the
canonical equation B_g ln Π = κρ_s was derived precisely in H1.

## 5. Caveats (so that no more is read than holds)

1. **It is a postulate, not a derivation.** It replaces two choices by one and gives it a meaning (continuity of the
   causal traffic), but why the substrate conserves the throughput is not yet explained.
2. **It depends on the choice of coordinates.** The throughput refers to the Killing time (clocks at infinity) and to
   isotropic spatial coordinates. Per unit proper area of the front it comes out as ρ_c c / B^{(d_s − 1)/2} = e^{2φ},
   a different quantity, which does not fix γ. Isotropic coordinates are those in which the substrate "looks the same
   in all directions", which is natural for an isotropic substrate. A covariant formulation, however, is missing. The
   nearest candidate is that the area of the front is also measured by a number of events (horizon molecules, K.14).
   Then the conservation of throughput would be the counting version of "the vacuum does not focus" (section 2) and
   the two routes would meet.
3. **It needs the creation of events.** Densification to Π⁻² means that matter adds events, not only removes links.
   Without a growth dynamics (axiom 7) the principle remains a description of the target state, not a mechanism.
4. **A numerical coincidence is possible.** γ = 1 exactly for d_s = 3 is tempting. A test that can refute the
   principle, not only confirm it, is therefore necessary.

## 6. What follows for the laboratory (proposals, nothing running)

1. **F7.1b, a discrete realisation of the principle (2+1):** the capacity narrowing of the cones of F7.1a and in
   addition events added up to ρ_c = 1/c. The prediction in 2+1 is γ = 2, not 1. The test checks that a consistent
   metric (clocks and cones) with the right γ comes out of a discrete order with added events. It is a check of the
   principle and of the instrument outside its calibrated range (K.16b knew γ only from −1 to 1). In 2+1 it is cheap.
2. **K.16 in 3+1** (density → φ − 3ψ) and F7.1b in 3+1 with the prediction γ = 1. More expensive, only if point 1
   passes.
3. **A theoretical question for the author:** is the conservation of causal throughput acceptable as a principle of
   the substrate? If so, the next step is to look for a dynamics (growth) that satisfies it by itself, and for a
   covariant formulation through the numbers of events on the front.

## 7. The ontology of *Physics of Causality* and the principle of throughput (27 Oct, at the author's request)

The author's formulation: **"throughput according to saturation; the more saturated a level of the light cone, the
more materialised it is".** I went through the text of the book. The principle is there, scattered over the chapters;
it was only not joined into a formula.

**What the book already contains:**
- **5.3, 7.7:** c is "the maximal stable throughput that avoids congestion". c_eff ≤ c falls when the impedance grows
  or the permeability falls, because fewer admissible channels remain.
- **26.1:** the Einstein limit requires that the propagation admits a "uniform maximal throughput".
- **23:** gravity is "a spatial redistribution of update density induced by gradients of impedance in a substrate
  with a conserved flow".
- **27.5:** rigid structures push configurational freedom into their surroundings and increase the branching there.
- **4.6:** in extreme saturation the effective connectivity grows (D_eff → 4).
- **6.5:** Π is a decreasing function of the proximity to saturation. **19.1:** matter is Π → 0, Z → ∞.
- **12:** the tension χ_ρ = ln(ρ/ρ_vac), an activity relative to the vacuum.

**Formalisation:** a level (layer) of the light cone has the saturation S = ρ_c/ρ_vac (events per coordinate volume
relative to the vacuum). Its throughput is q = S · c_eff. By 26.1 and 5.3 the throughput in the regime of the
Einstein limit is uniform and maximal, so q = const. Hence c_eff = 1/S: a more saturated layer is narrower. Full
materialisation is S → ∞, c_eff → 0, Π → 0, in agreement with 19.1. In 3+1 with the clocks Π (26.6: the clock rate
∝ Π/I) it gives S = Π⁻² and c_eff = Π², i.e. exactly H1 and γ = 1 (section 4). If ρ/ρ_vac of section 12 is read as a
density of events (an identification, not a derivation), χ_ρ = ln S = −2 ln Π.

**A contradiction in the text that the author must resolve:** section 26.6 claims that gravitational time dilation
does not affect the causal order ("the partial order defined by causal reachability remains invariant across the
substrate"). For γ = 1 the opposite holds: the narrowing of the cones is a change of the order. A fixed order gives
γ = −1 (A6, F3, measured by the laboratory). Sections 5.3, 5.6 and 7.7 at the same time admit a change of
reachability (c_eff, reshaping of paths around a defect). Proposed amendment: causality (no loops, continuity of the
flow) stays, but the order changes near matter. The cones narrow and the layers saturate.

**What does not follow from the ontology:** why the throughput stays constant when the saturation changes. Section
26.1 states it as a condition of the regime, not as a law. The Lorentz factor of 26.5, moreover, puts in a
Pythagorean division (axiom audit). The principle of throughput is therefore in the ontology as a description of the
Einstein regime. The laboratory can test it (F7.1b); it cannot yet derive it from a dynamics.
