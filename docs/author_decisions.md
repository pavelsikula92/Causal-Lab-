# The author's decisions on open questions (22 Oct 2026)

A record of the decisions and of their direct consequences for the laboratory. Decisions are not changed
retroactively according to results (R132).

## 1. The field equation and the bending of light in 3+1

**Decision:** γ = 1 remains the physical target; the field equation will be reworked. The present equation (a) is
not rescued by additional tuning, and the term −R/2 is neither added nor removed just for the sake of a number.

- H1 is split into two things: a **law of propagation** c_eff = cΠ² (a working postulate for the cones, F6.3b) and
  an **effective gravitational metric**, from which the temporal behaviour, the spatial curvature and the motion of
  light follow. It cannot be claimed that H1 follows from the conformal metric Π²η (a conformal metric has flat
  cones).
- Weak-field target: ds² = −(1 + 2φ/c²)c²dt² + (1 − 2ψ/c²)dx², ψ = φ (γ = 1).
- **Why (a) gives a zero response in 3+1 at γ = 1 (derived 22 Oct, sympy, weak field Π = 1 + εu):** in the metric
  −Π²dt² + Π^(−2γ)dx² with d_s spatial dimensions (D = d_s + 1 in total)
  - □_g Π → +1·∇²u (independent of γ and d_s),
  - −(R/2) Π → (1 − γ(d_s − 1))·∇²u,
  - sum λ₁ = 2 − γ(d_s − 1).

  d_s = 3, γ = 1 ⇒ λ₁ = 0: −R/2 gives exactly −∇²u and cancels □_g, the Newtonian response disappears (an exact
  algebraic cancellation of the linear term, not a numerical artefact). d_s = 2, γ = 1 ⇒ λ₁ = 1: the contribution
  of −R/2 is zero, which is why (a) with H1 works in 2+1. The term −R/2 is not arbitrary: it is the continuum limit
  of the Benincasa–Dowker operator (□ − R/2, K.11), so removing it would mean abandoning the BD operator.
- **Decision (second round, 22 Oct): the BD operator stays; what changes is what it acts on: B_g f(Π) = ρ_s.**
  At linear order □_g f(Π) → f′(1)∇²u and −(R/2) f(Π) → f(1)(1 − γ(d_s − 1))∇²u, so in 3+1 at γ = 1
  λ₁ = f′(1) − f(1): f = Π gives 0, f = Π^k gives k − 1, and every function with f(1) = 0, f′(1) = 1 (ln Π, Π − 1,
  (Π^k − 1)/k) gives 1. The linear response therefore does not separate them.
  - primary candidate f = ln Π (motive: clock rates compose multiplicatively, ln turns them into an additive
    quantity; a motive, not a derivation – it must be shown that the definition of Π in the theory has such a
    composition); a candidate, not a new axiom;
  - a mandatory linear control f = Π − 1; non-linear controls (Π^k − 1)/k for several k;
  - the difference is sought beyond the linear order: B_g f(Π) = f′(Π) □_g Π + f″(Π)(∇Π)² − (R/2) f(Π),
    (∇Π)² = g^{μν}∂_μΠ ∂_νΠ; ln Π carries an extra −(∇Π)²/Π², Π − 1 nothing.
  - Benchmark (a static spherically symmetric source in 3+1, the same ρ_s(r), boundary conditions and weak-field
    normalisation for all candidates): (1) weak field – the same linear response; (2) non-linear solution – a regular
    profile returning to Π → 1; (3) the second potential ψ leads to γ = 1; (4) an acceptable solution when the size
    and strength of the source change, without additional tuning of f. ψ and γ must follow from the whole model, not
    from a choice for an individual test.
  - Note of the laboratory on point (3): for f with f(1) = 0 the term −R/2 drops out entirely at linear order (R is
    O(ε) and is multiplied by f(1) = 0), so in a weak field the equation is □ f(Π) = ρ_s whatever the geometry; the
    difference between BD and □_g shows only non-linearly. The second potential: if the metric is given by the ansatz
    H1 (−Π²dt² + Π⁻²dx²), then ψ = φ and γ = 1 are an input, not a result, and test (3) would be tautological. A
    non-trivial test of γ requires reading both potentials independently – the temporal one from clocks (Π, lengths
    of chains), the spatial one from spatial distances or from the width of the cones in the dynamic order (F6.3b) –
    and only then comparing them.
  - **Confirmation by the author (23 Oct, after K.13a):** H1 (γ = 1) stays and the field equation changes, i.e. the
    route of the second round. K.13a verified the reason with its own instrument: in the H1 geometry in 3+1 the
    coefficient of ∇²u in B_g 1 goes to −1, the linear term of equation (a) cancels. The laboratory's proposal of the
    same day, B_g Π − Π·B_g 1 = κρ_s, is withdrawn: in the continuum limit it is □_g Π, i.e. the removal of the term
    −R/2, which this decision excludes. At linear order it coincides with the control f = Π − 1 and differs only at
    second order. The candidates f are separated by F1.7 N, which was waiting for K.13a.

- **Decision (fourth round, 23 Oct, after F1.7 N): f = ln Π is the primary, canonical scalarisation of the field
  equation.**
  - Equation: B_g ln Π = κρ_s. The BD operator stays whole including −R/2; in the continuum limit and in H1 (γ = 1,
    3+1), with L = ln Π: e^(2L) [ (1 − L) ∇²L + L |∇L|² ] = κρ_s; in a weak field ∇²L = κρ_s (Newtonian potential
    L = −U).
  - Reason: in the family (Π^k − 1)/k we have β = 1 − k/2, and only ln Π (k → 0) gives β = 1, i.e. the perihelion
    advance of general relativity (derived 23 Oct, docs/findings.md, Derived consequences). In addition the motive of
    composing clock ratios.
  - Note of the laboratory on the status (so that later no more is read into it than holds):
    - the term −R/2 was not removed: it cancels only at linear order, because it multiplies f(Π) and f(1) = 0;
    - the choice is SELECTED BY THE DATA (perihelion), not a prediction: β = 1 is now an input, not a confirmation of
      the theory; the observation fixes only f″(1) = −1, "the only one" holds within the family (Π^k − 1)/k;
    - γ = 1 remains an input of the ansatz H1, not a result;
    - the testable predictions of the equation come later: the discrete solution on a causal set, the bound on the
      existence of static strong-field solutions (F1.7 N: ln Π between u_lin(0) = −0.3 and −0.5), the 2PN order and
      the narrowing of the cones in the order (F7.1).
- **Decision (fourth round, 23 Oct, after F1.7 N): the canonical field equation is B_g ln Π = κρ_s** (in the
  continuum limit (□_g − R/2) ln Π; the BD operator stays including the term −R/2, only the argument changes).
  - Reason: in 3+1 with H1 the equation B_g f(Π) = κρ_s has the PPN parameter β = (1 − f″(1))/2, so in the family
    (Π^k − 1)/k β = 1 − k/2; the perihelion (|β − 1| of order 10⁻⁴) requires f″(1) = −1, which within the family only
    ln Π meets (sympy, findings).
  - Limits: γ = 1 is an input of the ansatz H1, not a result; the observation fixes only f(1), f′(1), f″(1); the
    higher derivatives of ln Π (higher PN orders, strong field, the bound on static solutions from F1.7 N) remain a
    hypothesis; κ and the conversion to G are still to be derived.
  - Note of the laboratory on the author's notation: the form "B_g ln Π − Π·B_g 1 = κρ" (from the withdrawn
    proposal) would add a term +(R/2)Π, which is of first order: it doubles the linear coupling (coefficient 2
    instead of 1) and gives β = 5/4, a perihelion advance excluded by observation. The form recorded is therefore
    B_g ln Π = κρ_s, on which the derivation of β = 1 rests. The term −R/2 does not affect β (for f(1) = 0 it is of
    third order in vacuum; □_g ln Π without it also gives β = 1); it shows inside the source and at higher orders.

- **Decision (third round, 22 Oct): K.13a stays the main discrete test, K.13b in its original form is dropped, the
  note on ln Π is a separate analytic justification.**
  - Reason: the BD operator is linear in the field. For f(1) = 0, E_g[f(Π)] = ε f′(1) B_0[u] + O(ε²), and the
    difference ln Π − (Π − 1) at second order is −½ ε² B_0[u²] (the flat operator; ∇²(u²) only in the continuum).
    The response of the operator therefore does not separate f. Only a principle and the whole non-linear equation
    (F1.7 N) separate it.
  - K.13a measures the operator itself: E_g[B_g w] = B_0[w] + ε δB_u[w] + O(ε²) for fixed fields w independent of
    ε. The field w = 1 gives −R/2 = −ε[γ(d_s − 1) − 1]∇²u, a Gaussian w gives the correction of □_g, and changing
    ℓ_k/σ the convergence; the systematics of the finite scale are given separately. Mandatory checks: Π = 1 (flat
    operator) and γ = −1 (K.11).
  - Principle for f: the composition of clock ratios (docs/derivations/ln_pi_potential.md). It is a motivation, not
    a derivation. Uniqueness would require a new axiom of additivity of the field variable. The action and the
    conservation of the flow are two independent control routes.
  - Note of the laboratory: in H1 the first order of the operator is linear in γ, δB = (1 − γd_s) P + (1 + γ) Q. P
    is the density part (flat cones, instrument K.11), Q is the part from the change of the speed of the cones.
    K.13a measures both parts separately.

## 2. The flow in a dynamic order

**Decision:** a global normalisation by the number of future descendants is not a basic physical rule; in F6.3b it
stays only as a reference computation. Physical candidate: a local reservoir model (stock Q_v; when a link is
created it receives a share of the available stock, the rest stays), with the balance Σ Φ_in + s_v = Σ Φ_out + ΔQ_v
in every step.

Conditions: (1) no knowledge of future descendants, (2) balance in every step, (3) results independent of the order
of the computational updates.

**A tension that must be resolved before the rule is chosen (note of the laboratory):** the descendants of an event
arise in an order that is a linear extension of the order, and a different extension (e.g. another inertial frame)
sorts them differently. A rule of the type "every new link receives a share α of the remaining stock" gives the
descendants shares that depend on this order, and so violates condition 3. A fixed share of the original stock per
link is independent of the order only until the stock is exhausted; after that it matters which descendants came
first. Conditions 1–3 together therefore restrict the rule strongly; the independence of the order must be measured
(bitwise agreement in different linear extensions, as the laboratory has done since F3.3b).

## 3. The background flow

**Decision:** no universal background flow is introduced. The basic test is at s_bg = 0. Regimes: s_bg = 0 (a
stand-alone droplet), s_bg > 0 (a droplet maintained by the surrounding flow), a scan of s_bg (a critical value or a
robust region). A droplet that exists only at s_bg > 0 is a "structure maintained by an external flow", not a
self-bound object. A cosmological reading of the background only after numerical causes are excluded, and only from
a further rule. (F6.3b A has s_bg = 0 in its grid.)

## 4. The copied pulse

**Decision:** the pulse is not an axiom, it is a hypothesis of the laboratory; no physical meaning (frequency, phase)
is attributed to it yet. The following are distinguished: the state of the pulse p_v (a carried tag), the flow Φ
(balanced), the stock Q_v, the clock response Π_v. Copying an informational state must not copy a conserved quantity.

**Consequence for the models so far (note of the laboratory):** in F3.3c–F3.3h and F6.3/F6.3b the pulse acts as an
internal source, s_eff = s + P_in, so every event with a pulse adds ω to the flow. The balance in the code holds only
because the pulse is counted as a source; physically it is a flow created by copying a state, which this decision
excludes. The self-sustaining phase of F3.3f after the source is switched off, and the droplet of F6.3b, are driven
precisely by this injection. The results of these tests hold for the model "pulse as a source"; the model according
to decision 4 (the pulse changes the response, e.g. Π or the capacity, but does not add flow) is new and not yet
tested.

## 5. The parameters Π_sat, m, ΔQ_crit

**Decision:** working parameters, not constants. F3.3e (scan of ΔQ_crit, selection by a score fixed in advance,
confirmation on new seeds) is done: the optimum 0.2 FAILED on new data. F3.4: sensitivity to Π_sat and m including
the interaction with the threshold; a combination is not selected and confirmed on the same data. A stable region is
sought, not a sharp peak; a variant that does not hold on new data is rejected or reformulated.

## 6. Saturated matter as a source of gravity

**Decision:** the simplest mapping S_v = 1 when ΔQ_v ≥ ΔQ_crit, otherwise 0; ρ_s(v) = ρ★ S_v (ρ★ a single
normalisation, not derived; the first test follows the shape of the field, not the absolute strength). Loop: load →
S → ρ_s → Π → transport → new load, with a local balance after every step, a fixed order of updates and a
convergence condition. On failure, first determine whether the mapping, the field equation or the feedback fails.

Note of the laboratory: the field solver (K.11c, F2.7, K.12) is static and radial; for a droplet on the axis the
radial approximation is possible (the source averaged over sprinklings and in time), a general non-static source
requires a new solver. Decision 1 holds: until the field equation in 3+1 is decided, the loop is tested in 2+1.

## 7. Growth of the order

**Decision:** yes, but not in the first version. F6.3b: dynamics of the relations (fixed points). Next phase:
dynamics of the birth of events with a local rule; a new event must not arise in the past of existing ones;
predictions must not depend on an arbitrary order of birth (template: sequential growth models of causal sets).
Growth of N is not by itself an expansion: the spatial geometry must be reconstructed and a change of physical
distances shown.

## 8. The length ℓ_d and the bridge to SI

**Decision:** G is not a predicted constant. G = κc³ℓ_d²/(4πλ₁ħ) remains a conditional bridge; ℓ_d = ℓ_k is not
introduced without a derivation. Routes: derive ℓ_d from a mechanism (e.g. a stable defect), or determine it from one
observation and then predict other quantities without further tuning (G is then not a prediction).

## 9. Gravity as the thermodynamics of the causal structure (22 Oct)

**Decision:** try the Jacobson route (Jacobson 1995): entropy of local horizons proportional to the area + the
Clausius relation + the Unruh temperature ⇒ the Einstein equations as an equation of state. The causal set itself
supplies the entropy as the number of horizon molecules. Analysis of the assumptions: `docs/derivations/jacobson_route.md`;
first instrument K.14.

Notes of the laboratory:
- The choice of the substrate (axiom audit, section 0) is not yet recorded formally. The author made it conditional
  on finding the essence of gravity. The Jacobson route rests on a Lorentz-invariant substrate (local Rindler
  horizons in every direction), so its success would support choice A.
- The route of the field equation B_g f(Π) (decision 1, K.13a, F1.7) stays valid as a comparison. If the Jacobson
  route passes, the field equation is not chosen but derived, and γ = 1 stops being an input.
- The missing piece is matter with an energy–momentum tensor (assumption 5). Without momentum the matter of the
  theory cannot be put into the Clausius relation.

## 10. Axiom 5: a general response and candidate response laws (27 Oct)

**Decision:** axiom 5 is split into a general principle of response and separately testable response laws. The
original relation Π_v = f(σ_v) thereby changes from an axiom into a candidate closure of the model.

- **Axiom 5 (local response):** every event v has a local response variable Π_v, 0 < Π_v ≤ 1. The rigidity I_v ≥ 1
  measures the resistance to a change of the order, the impedance is Z_v = I_v/Π_v. The axiom does not fix the
  functional dependence of Π_v on the load, the internal causal structure or other causal variables. An
  independently testable response law must fix it.
- **Candidate response laws:**
  - **R1, response from the load:** Π_v = f(σ_v), σ_v = F_v/C, f decreasing;
  - **R2, response from the field:** Π_v is fixed by a field equation with a source from internal causal events
    (today B_g ln Π = κρ_s, decision 1); the dependence on the load is limited to the saturation regime (A5′).
  - R1 and R2 are alternative closures of the model until their compatibility is shown. Neither follows from the
    general axiom.

Notes of the laboratory:
- What the tests have already said about both laws: R1 alone gives a field without range, which heats the bundle and
  does not bind it (F3.7). R2 gives a field with range, which attracts the bundle (F3.8, F3.11). They are compatible
  if R1 acts only in saturation, as A5′ formulates it.
- The conflict of equation (a) with H1 in 3+1 (decision 1) was caused neither by A5 nor by A6. It was caused by the
  form of the field equation, i.e. an additional assumption (λ₁ = 2 − γ(d_s − 1) = 0), and it was resolved by the
  equation B_g ln Π = κρ_s.
- The axiom audit (section 2) lists "Response Π = f(σ) (axiom 5)" as firm ground. After this decision it is firm
  ground as the local rule R1, not as an axiom.
- A related open matter: the content of axiom 6 depends on the definition of Δτ (decision A1 of the variational
  audit), see `docs/findings.md`, Derived consequences, 27 Oct.

## 11. The instrument for F6.4b: manifold phase and transition (27 Oct, after F6.4a and F6.4c)

**Decision:** F6.4b (matter in the action) may build on the instrument of 2D orders with the BD action, which
reproduces the manifold phase and the position of the transition, but not the height of the crystalline phase. This
is a decision after a failure, not a change of the criteria: F6.4a keeps FAILED and F6.4c UNDECIDED.

Reason: phase I (χ, height, action, d_MM) and the transition (β ≈ 2.4) agree with the paper to hundredths; the exact
equilibrium and β = 0 against direct sampling passed. The discrepancy is confined to the metastable crystalline
phase. It has at least two long-lived states (F6.4c), and its height 4.18 in the paper cannot be reproduced by any
preparation. Matter is put into the ensemble in the manifold phase.

Notes of the laboratory:
- Measuring γ (φ from clocks, ψ from the width of the cones) needs a transverse direction, i.e. at least 2+1. Before
  F6.4b there will therefore be a calibration of 3D orders (literature to be checked against the full text, as for
  F6.4a).
- The crystalline phase is not used in F6.4b. If something needs it later, the caveat of F6.4c applies to it.

## 12. Gravity changes the causal order (27 Oct, according to the results of the laboratory)

**Decision:** section 26.6 of *Physics of Causality* ("gravitational time dilation does not change the causal
order") is amended according to the results of the laboratory. Causality stays: no closed loops, conservation of the
flow (axiom IV). But the causal order changes near matter: the cones narrow (c_eff falls) and the levels of the light
cone saturate (more events per coordinate volume). The author's formulation holds: throughput according to
saturation; the more saturated a level of the cone, the more materialised.

Reasons (laboratory):
- a fixed order gives γ = −1, light does not bend (A6 on a sprinkling, the F3 series; derived 27 Oct);
- narrowed cones alone, without added events, give γ = 1/d_s (F7.1a, 2+1: 0.47 / 0.43);
- conservation of the throughput q = S·c_eff gives γ = 2/(d_s − 1), so γ = 1 in 3+1, and with the clocks Π exactly
  H1 (`docs/derivations/why_gamma1.md`).

Notes of the laboratory:
- A constant throughput is in the book a condition of the Einstein regime (26.1), not a derived law. Why the
  substrate keeps it remains open.
- The discrete realisation of the principle is tested by F7.1b (2+1, prediction γ = 2).


## 13. Throughput is a plateau law, the non-linear completion is general relativity (28 Oct, path A)

**Decision:** fork A of `docs/derivations/tensor_coupling.md` holds.
- The principle of conserved throughput is a law of linear order (plateau). In 3+1 it is the same as a coupling to
  the whole energy–momentum tensor: h̄_ab = 4U u_a u_b.
- The non-linear completion is general relativity. It is to be derived by the Jacobson route (decision 9), not
  chosen.

Consequences:
- The scalar equation B_g ln Π = κρ_s stays valid as a static section to 1PN (β = 1). It is not a complete dynamical
  law: the scalar radiates 1/6 of the quadrupole power of general relativity and has no frame dragging.
- The predictions of the scalar equation beyond 1PN are not predictions of the theory: the throat at 2.341 M, the
  missing horizon, the shadow +0.91 % and the 2PN coefficients of `docs/derivations/gravity_equations.md`.
- The source is the T_ab of the flow (K.15), not the scalar ρ_s. Decision 6 (S_v as the source) holds as the rest
  component T_ab u^a u^b.
- What remains to be checked on the substrate is how the saturations of flows with different 4-velocities compose
  (F7.4a). The non-linear order belongs to the Jacobson route (steps 2 and 4).

Note of the laboratory: the test that the author called F7.2 has the ID F7.4a in the registry. F7.2 is already
assigned to the Shapiro delay, and F7.4a follows the postponed F7.4 (tensor gravity from the momentum of the flow).
Unlike F7.4 it puts in neither the trace reversal nor a Green function. The tensor structure can arise only from the
response rule of the order.
