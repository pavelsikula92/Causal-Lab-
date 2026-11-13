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


## 14. H5: inertia of a locked defect from a conservative dynamics of the order (10 Nov)

**Decision:** a locked defect must move as a whole through its own internal binding. No carrier flow enters the base
hypothesis. The card is `docs/H5_card_draft.md`.

| point | decision |
|---|---|
| D1 | the field lives on all events of the order, valued in an internal target space |
| D2 | the dynamics is conservative and of second order (relaxation cannot supply inertia) |
| D2′ | in 1+1 it is implemented through the retarded Green function of the order (½ × the causal relation) |
| D3 | 1+1 first, with the sine-Gordon kink as the analytic reference |
| D4 / H5-D | the potential is to be derived from the Axiom of finite capacity of the flow through events, last and separately |
| BD operator | stays in the laboratory as a separate instrument; it is not the dynamical base of H5 |
| 3+1 | postponed (Johnston's propagator is supported on links; the light-cone problem returns) |
| K₀ | removed from H5 and from its interpretation (a superseded version of the model, not a parameter of the theory) |

Conditions of the author:
- the instrument first (H5-0): stability, dispersion, convergence, independence of the order of computation;
- the rest mass and the inertial mass are measured independently;
- if the capacity axiom does not give the sine-Gordon potential, that is a result about the admissible localised
  configurations, not a failure of the instrument.

Note of the laboratory: the BD operator in retarded form amplified the Poisson noise in the development probes (a
constant field diverged within two units of time). The Green function of the order was stable.


## 15. The chain from the field to gravity: the response comes from the action of the order, capacity in saturation (10 Nov)

**Decision:** the author accepts the chain of `docs/derivations/field_and_throughput.md`, section 6:

  φ → energy (the local contribution 𝓔_v, not yet defined on the order) → the source of the response of the order
  (the action of the order, R2) → clocks, cones, density (throughput) → the resistance that a probe feels.

The capacity acts in a strong field (saturation, A5′).

Reason (the author's): F3.7 showed clearly which way not to go. A response from the load alone (R1) gives a field
without range, which heats a bundle and does not bind it.

Notes of the laboratory:
- The other decisions of the note (φ as the phase of the clocks, the distinction of the three currents, the definition
  of 𝓔_v, the intrinsic initial data) stay open.
- The active gravitational mass (F5.2b) needs 𝓔_v and a dimension from 2+1 on; in 1+1 the throughput is degenerate.


## 16. Direction B (2+1): the four decisions (11 Nov)

**Decisions:**

| | decision |
|---|---|
| D1 | the proper time τ in the Green function of the order (1/2πτ in 2+1) is estimated from the **longest chain**, τ̂ = L_max / (α₃ ρ^{1/3}). α₃ is calibrated before freezing, also at several levels of thinning (so that a change of density is not taken for a change of proper time). The interval count N ≈ (π/12) ρ τ³ is a diagnostic only. A coordinate stage of the instrument is allowed for development (form of G, calibration of α₃, comparison of the estimators, near-cone and finite-size biases); coordinate τ must not enter a frozen test. |
| D2 | develop a discrete gradient. Order of work: ∇φ from the causal neighbourhood → T_ab → its components against the continuum → its trace → only then the gradient in the load. Signature and normalisation of the Lagrangian are fixed before freezing; T = −½(∂φ)² + 3V is the target expression of the chosen convention, not yet a verified discrete identity. |
| D3 | the local rest frame of matter as the timelike eigenvector of T^a_b (T^a_b u^b = λ u^a, u·u = −1, future-pointing) is accepted as a **candidate**, needing its own instrument and its own card. No unique timelike eigenvector → no cut and the event is marked unclassifiable (no fallback to coordinates). The direction of a link relative to u must be read relationally; coordinates only for development. The strength of the cut must be a function of the load or of Π fixed in advance, not tuned to the bending of light. It must be stated whether only the propagation changes or the relation ≺ itself, and how transitivity is kept. Thinning alone does not narrow the cones. In 2+1, γ = 2 is a coefficient of the theory, not the PPN parameter of general relativity. |
| D4 | a domain wall (the kink extended along y) is the first 2+1 object. Energy per unit length; the dependence on the box size in y is followed. The 1/r law belongs to a localised defect (B4). |

**Stages:**
- B0 (Green function, intrinsic τ, stability, dispersion, convergence): approved.
- B1 (inertia and stability of the wall): approved.
- B2 (response in a prescribed order, the intrinsic rest frame and the narrowing of cones): approved conditionally.
- B3 (causal growth with a dynamic change of the relation, active mass): approved as the next phase.
- B4 (localised defect, radial response): a separate research goal.

**Two boundaries:**
- B2 (a prescribed order, motion in its response) and B3 (matter changes the relation during causal growth) are
  different experiments. The two changes must not be switched on together.
- Success in B2 is not a confirmation of the active mass. It shows a consistent mechanism by which the stress tensor
  fixes the rest frame and the rest frame affects propagation. The active mass is tested in B3.


## 17. B0: the calibration of the Green weight of the longest chain (11 Nov, 18:45)

**Decisions:**

| | decision |
|---|---|
| condition | the weight of a pair with longest chain L is calibrated by **E[w(L) \| τ] = G(τ)** (the expected contribution at a given proper time equals the continuum Green function), not by E[G(τ) \| L]. The aim is to keep the expected Green integral of a source, not to estimate G at a known L. |
| scaling | w(L; ρ) = ρ^{1/3} w̃(L), s = ρ^{1/3} τ; checked on independent sprinklings at several densities. |
| method | two moment conditions on the residual r(s) = E[w̃(L) \| s] − 1/(2πs): ∫₀^{s_max} s r ds = 0 and ∫₀^{s_max} s² r ds = 0. They cancel the leading error for a smooth source with local pair measure (a s − b s² + O(s³ ln s)) ds. s_max is set at the calibration and then frozen; it is never changed according to which source comes out better. |
| form | w̃(1) = a (links separately), w̃(L) = c / (L + β) for L ≥ 2; three parameters with a > 0, c > 0, β > −2 (a positive, monotonically falling tail). The scratch values a ≈ 0.257, c ≈ 0.200, β ≈ 0.74 are a starting point only, not frozen constants. |
| free table | only a fallback, opened as a new calibration hypothesis if the parametric form fails the independent validation; never a switch because the parametric form gives an unwelcome result. |
| validity | sources smooth on the scale ℓ ≫ ρ^{−1/3}, on the base Poisson order (homogeneous, unthinned, causally convex region). |

**What is exact and what is not.** The calibration is not "exact by construction". Exact (within the numerical
uncertainty) are the two enforced moment conditions. The accuracy of the whole Green integral for the declared class of
sources is a consequence to be tested. The two moments do not make r(s) = 0 pointwise, so the local error is reported
next to the moments (the scratch value E[w | s]/G ≈ 1.36 at s = 1–1.5 must not hide behind a zero total bias).

**Why the calibration does not depend on the source (on the base order).** The past of a source event y and the open
interval I(y, x) are disjoint; Poisson increments in disjoint regions are independent; a causally convex region does not
cut the interval. So P(L | τ) does not depend on J(y) when J(y) is fixed by the past of y.

**Limit.** This does not hold automatically for a dynamically thinned order (B3): there the realisation of the end
point x is decided from its past, and conditioning on x being realised can change the distribution of the events in
I(y, x). The B0 calibration is a calibration of the Green function on the base Poisson order; its transfer to the
dynamic order needs its own check.

**Checks required before B0 is closed:**

| check | requirement |
|---|---|
| density | w ∝ ρ^{1/3} on independent sprinklings |
| links L = 1 | bias and uncertainty reported separately |
| moment 1, moment 2 | zero within the calibration uncertainty |
| smooth source | several widths and orientations not used in the fit |
| non-smooth source | reported as a diagnostic outside the declared class, not as an automatic failure |
| dynamic thinning | whether the conditional P(L \| τ) stays usable (for B3) |


## 18. The aim of the classical programme: gravity and the dimensions of spacetime; the defect is a test body (11 Nov, 20:17)

**Decision.** The classical programme asks one question: **how the causal order carries gravity, and what changes with
the dimension** (1+1 → 2+1 → 3+1). What a particle is belongs to the quantum part, which is built on the classical base
afterwards.

**Consequences:**
- The kink (and any later localised defect) is a **test body**: the simplest localised energy that keeps its shape on
  the order and can be weighed, launched and let fall. Results with it are statements about the order and gravity, not
  about the nature of particles or matter.
- The field φ and its potential are chosen as an instrument. That matter does not arise from the order by itself
  (F3.2, F3.3, F3.13d) is not a gap of the classical programme.
- In 2+1 and 3+1 a plain scalar has no stable localised soliton (Derrick). The localised test body for the 1/r law (B4)
  is chosen pragmatically, as an instrument, and recorded as such, so that it does not later turn into a claim about
  particles.
- What each dimension adds: 1+1 counting alone (inertia, universal fall, active mass); 2+1 the change of the relation
  (bending of light, γ of the theory) and a force falling with distance; 3+1 Huygens (the Green function on the cone
  only) and the comparison with general relativity, including γ = 1.


## 19. D2: how to break the circle gradient ↔ frame (11 Nov, 20:23)

The development rows of D2 (findings, Undecided, 11 Nov evening) showed that a compact neighbourhood of an event needs a
frame, while the frame of D3 is to come from T_ab, which needs the gradient; and that the pointwise T_ab of one scalar
field has a degenerate timelike eigenspace wherever the gradient is spacelike (every wall and kink).

**Decisions (the laboratory's recommendation, accepted):**
- **The load (B-b) without components:** only the scalar (d phi)^2 = 1/2 box(phi^2) + phi V'(phi) is used (signature
  +, −, −; box phi = −V'(phi) is the field rule itself), so the load needs no frame.
- **The rest frame (D3) from coarse-grained matter:** the timelike eigenvector of T_ab summed over a region (the momentum
  of the defect), not the pointwise T_ab. For localised matter the sum converges over a Lorentz invariant region because
  the source is compact. For a wall the boost along the wall stays free (physical).
- The self-consistent iteration (frame → compact interval → gradient → frame) is kept as a fallback.


## 20. B0 is closed; the intrinsic instrument of 2+1 is variant B (11 Nov, 23:17)

**Decision.** Stage B0 (the Green function of the order in 2+1 as a field evolution) is closed as a development stage.
The intrinsic instrument for the frozen tests of direction B is **variant B**: the longest chain from the order (forward
dynamic programme over the links) with the weight of decision 17, w(1) = a ρ^{1/3}, w(L ≥ 2) = c ρ^{1/3}/(L + β),
a = 0.2733, c = 0.2131, β = 1.203 (diag/out/b0_weight_calib_20261111T174914.json). Coordinate τ stays a development
reference only; variant A (one α₃) is not used (its ω bias of 1–3 % does not shrink with ρ).

**What B0 established (findings, Undecided, B0 development rows):**
- algorithm: links identical to core.order.links; L from sources equal to an independent interval chain on 373 / 373
  pairs; forward kernel equal to the backward DP to 3e-15; the field independent of the linear extension to 3e-15;
- calibration: the moments of r(s) zero on an independent sample; density scaling and boost independence of
  E[w] / E[G] to ~0.01; Green integral 0.996–1.006 at ρ 200 for sources not used in the fit;
- field (ρ 25 → 50): ω/ω_c 1.004–1.007 → 1.000–1.002; vacuum λ 0.0068 → 0.0038; the wall moves at the launched
  velocity;
- the growth of the vacuum seen with coordinate τ came from its cap (λ follows the cap).

**Open items carried into B1 (reported there, not hidden):**
- the noise of single events is larger under B (per-event rms in the wall core 0.30 against 0.20 at ρ 25); a
  distortion of the profile is not resolved;
- single waves decay by 1.2–3.0 % over T = 6 at ρ 50 (shrinking with ρ);
- the dynamic thinning (B3) needs its own check of P(L | τ) (decision 17);
- cost: the intrinsic field is O(N · E_H), ~7 min per run at ρ 50, T 6.


## 21. B1 is frozen as F3.15 with the proposals of the card, revision 1 (12 Nov, 11:01)

**Decision.** The author left the open choices of the B1 card to Claude ("dodej test a uvidíme", 12 Nov, 11:01). They are
fixed as proposed in docs/B1_card_draft.md, revision 1, without change after the development rows; the card is frozen as
revision 2 and implemented as **F3.15** (phases/F3_structure/F3_15_wall_inertia_2p1.py). Seed 0 decides; a replication
uses seed 1.

**The choices:**
- deciding density ρ 30, T 8, pyramid X = T + 8, Y = T + 2.5, band |y| < 2; **32 sprinklings** × 4 arms (rest, v 0.3,
  v 0.6, kick f 0.05 for t ∈ [1, 3));
- the slab rule: the last slab is not read in any arm (geometric; at v 0.6 the plateau window leaves the pyramid);
- tolerances (rule of 10 Nov, derivations in the card): rest Q 0.01, drift 0.15, w 0.05, E 0.05; velocity 0.03,
  contraction 0.05, width 0.05, energy 0.05 (each v); r_M 0.10;
- instrument: I1 vacuum |λ_eff| + 2 SE < 0.02 and I2 dispersion per k < 0.015 (4 sprinklings, the B0 waves); I3 the rest
  wall at ρ 60 (8 sprinklings) |w/w_ex − 1| + 2 SE < 0.05; I4 kernel, linear extension, region < 1e-12;
- controls on 8 sprinklings: telegraph Γ = 1, 2 (launched and kick; instrument, as F3.14b), V = 0 readout check
  (instrument), causality-blind (a PASSED condition, as F3.14b). The blind sum keeps the weights of each source and
  randomises only the receiving events (uniform over all later events): in 2+1 the weight depends on the pair, so the
  unit-weight form of F3.14b does not transfer;
- the y-box dependence is not part of B1 (reported later on a subsample, if wanted).

**One change after the smoke run (quick, ρ 10, decides nothing), before the deciding run:** the telegraph damping is
Γ = 1, 2 instead of the 2, 5 of F3.14b. At Γ = 5 the wall dissolved (Q 0.83 → 0.3, width 1.4 → 8) instead of
stopping; the most likely cause is that the source of ψ = e^{Γt/2}φ grows as e^{Γt/2} and so weights the region near the
cone, where the chain weight is pointwise least accurate (decision 17: E[w | s]/G from 0.42 to 1.46). At Γ = 2 the wall
stayed intact. Recorded as a limit of the intrinsic instrument for fast-growing sources (relevant for B3).

**Why the author's choices could be delegated:** none of them was taken after seeing a deciding run; the tolerances follow
the rule of 10 Nov from the development rows, and the development rows decide nothing.


## 22. B2 is frozen as F5.2c: walls in a prescribed curved order of the theory (12 Nov, 21:00)

**Decision.** The author asked for B2 to be prepared and delivered as a test ("ano prosím, dodej pak test pustím ho u
sebe", 12 Nov, 17:01), as for B1 (decision 21). The choices below were made by Claude after the development rows and
before any deciding run; card docs/B2_card_draft.md, revision 2; test **F5.2c**
(phases/F5_gravity/F5_02c_walls_curved_2p1.py). Seed 0 decides; a replication uses seed 1.

**What B2 became, and why (the continuum, before any discrete run):** in the metric of the theory
ds² = −Π² dt² + Π⁻² (dx² + dy²) a wall uniform along y has the passive weight E − σ (energy minus tension;
docs/derivations/walls_in_curved_order.md). The kink (σ = E) does not fall at all; it is an exact static solution at
every position. "The wall falls as the continuum" (decision 16) therefore became three conditions: the kink stays, the
breather wall (ω 0.7, σ/E 0.22) falls as the continuum, and a launched kink follows the cones of the order.

**The choices:**
- the order: an exact sprinkling of the metric of the theory with Π² = (x + ℓ)/ℓ, ℓ 24 (closed causal relation: the
  optical metric is the hyperbolic plane), density ρ 30 per proper volume, domain of dependence of −13 ≤ x ≤ 22,
  |y| ≤ 12.5, T 8. Not the link-cutting rule of F7.1b: that rule reads speeds in coordinates (a proxy in the audit) and
  is itself undecided (γ 1.85); a prescribed order in B2 should be the order of the theory itself. The rule-produced
  order is B3;
- the field rule of F3.15 unchanged (variant B, frozen constants); φ₀ and the exact reference from the y-independent
  continuum solver in the same metric;
- arms: kink at rest, kink launched at proper velocity 0.3, breather ω 0.7 (the same events); the breather on a flat
  order (ℓ = ∞, −16 ≤ x ≤ 22) as the control;
- readouts as B1 (band |y| < 2, slab rule); the breather centre from a fit of its profile (the energy centroid of B1
  is too noisy for an object that oscillates once in the run: per-sprinkling sd 0.40 against 0.13; reported);
- 32 sprinklings; the flat control on 8;
- tolerances by the rule of 10 Nov from the development rows (card, section 4);
- the light bending γ is not part of B2: in an exact sprinkling it is built into the relation.


## 23. B3a is frozen as F5.2d: the active side of the growth in 2+1 (13 Nov, 08:30)

**The author's choices (13 Nov, 07:35 and 07:50):** B3 goes by B3a, the cones of new events are cut in the rest frame of
the static source (the coordinate frame for a wall at rest, decision 19; a declared proxy as F7.1b; the intrinsic
reading of the frame is B3b); and B3 is the active side only. The fall of a test body moves to B4: walls do not fall in
the order of the theory (B2, F5.2c), and in the continuum the fall of a breather in a deepening potential is swamped by
its redshifted phase (diag/b3_continuum.py).

**What the development found, before freezing (diag/b3_continuum.py, diag/b3_growth.py):**
- the rule with the flat free wave differs from the theory by 1.6–2.5× if gravity acts from t = 0; with gravity
  switched on from t 3 to 4 the free wave is constant where the order is curved, and the rule agrees with the theory
  within 1–2 %;
- the load identity of decision 19 with f = φ² amplified the calibration error of the weight (f = 4π² on a plateau,
  the potential became asymmetric and the kink drifted); with f = φ(2π − φ)/2, zero on both vacua, it does not;
- the cone at the parent alone made long near-null links too wide where the potential deepens in time (r_L 1.6 at
  ρ 10); the cone from the binned potential at the midpoint of the pair gives r_L 1.00.

**Decision.** Card docs/B3_card_draft.md, revision 2; test **F5.2d** (phases/F5_gravity/F5_02d_growth_active_2p1.py,
md5 5a1028555a2b5c2fefced0705b5e0057; helpers diag/b3_growth.py 6bada302b8418ec1d4c0763d86f5df38,
diag/b3_continuum.py 3890f3a8630a7b918058ec9dc2d491f1). ρ 30, κ 0.01, T 9, kink sources μ 1 and 0.25, 32 sprinklings,
each paired with κ = 0; tolerances by the rule of 10 Nov from the development rows (card, section 4). Seed 0 decides.
