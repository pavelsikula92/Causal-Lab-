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


## 24. The active source of gravity is ρ + Σp in every dimension (13 Nov, 15:37)

**Decision (the author's):** the source of the response of the order is

    ρ_active = T_ab u^a u^b + h^ab T_ab = ρ + Σ_i p_i,   h^ab = g^ab + u^a u^b,

in the local rest frame u of the source, **in every dimension**. It is the rule of the substrate tested in F7.4a ("every
flow narrows the cone in its own rest frame, the capacity deficits add"), taken as a general rule; it is not adjusted to
the dimension-dependent weak limit of general relativity.

**Why:** consistency of the theory. F7.4a tests the addition of capacity deficits in the local frame of a flow; if that
mechanism is the rule of the substrate, it must not be replaced in 2+1 by another combination only because general
relativity has a different weak limit there.

**The boundary of interpretation:**
- In 3+1, ρ + Σp is the source of the Newtonian potential of general relativity (½(ρ + 3p) in R_00); the agreement
  with GR is specific to 3+1.
- In D spacetime dimensions general relativity gives ((D − 3)ρ + Σp)/(D − 2); in 2+1 only Σp (dust does not gravitate,
  point masses make conical defects). In 2+1 the rule (a) is therefore a prediction of the causal substrate of its
  own, not a reproduction of Einstein's gravity.
- A scalar ρ_active put into one scalar equation does not settle the question: for moving sources, anisotropic
  stresses and a general spacetime the response must be tensorial (the cones narrowed and tilted, as in F7.4a). That
  is open (B3b, C3).

**Consequences for the laboratory (recorded now, before any new run):**
- **The trace of T as the load (decisions 15 and 19, direction B) is the scalar branch (Nordström type).** For static
  dust it coincides with ρ + Σp; for anything with stress or motion it does not:

  | source | trace (direction B so far) | ρ + Σp (decision 24) |
  |---|---|---|
  | static kink in 1+1 | E | E |
  | kink in 1+1 moving at v | E/γ | γE(1 + v²) |
  | static wall in 2+1 | 2E | 0 |
  | wall in 3+1 | 3E | −E (repels) |
  | light | 0 | 2E |

- **F5.2b and F5.2d keep their verdicts (R132)**: they passed their own criteria, which tested the response to the
  trace. They are **not** read as a confirmation of the tensor source. F5.2b's static arms are the same under both
  sources; its launched arm (A2) and all of F5.2d tested the scalar branch only.
- In 2+1 under (a) a static kink wall does not gravitate (ρ + Σp = 0). An object uniform along y with energy E and
  tension σ has the active mass E − σ, the same combination as its passive weight in the order of the theory
  (docs/derivations/walls_in_curved_order.md); derived, not tested. A test of the active side in 2+1 under (a) needs a
  source with ρ + Σp ≠ 0 (a breather wall, a launched wall, light), and the load of the growth rule must be computed
  as ρ + Σp in the rest frame of the source.
- The identity of decision 19 (the load without the components of T) gives only the trace; ρ + Σp needs the
  components in the rest frame, i.e. the frame of the source (decision 19: from the coarse-grained T_ab of the defect).

## 25. C0: the apex term is removed by the local counterterm; MV2+ct and MV8+ct go on, MVA is a reference (13 Nov, 22:49)

**What the development found (docs/findings.md, "C0 development rows"; diag/c0_apex.py, diag/c0_minvar.py v2):** the
deficit at the apex of the cone and the mean dynamics of the linear 3+1 rule are the exact continuum expectation of the
layer kernel at finite density (predicted before the run, confirmed at ρ 30 and 60). Their integrated part, the apex
term δ = (a/8√π) Σ w_n H_n ∝ ρ^−1/2, shifts the mass. It can be removed from the continuum alone in two ways: the local
counterterm φ(v) − δ J(v) = φ0(v) − c Σ w J(u) (ct), or the constraint Σ w H_n = 0 in the kernel (MVA). At equal
systematics the constraint was about twice as noisy.

**Decision.**
- The apex term is removed by the local counterterm (+ct), with δ computed from (A2) for the kernel and the density, no
  measured number.
- The candidates of the next development phase are **MV2+ct** (the accurate one) and **MV8+ct** (the quiet one).
- MVA (MVA8, MVA16) stays as a comparison reference, not as a candidate branch.
- **The card C0 is frozen only after** the nonlinear response (J = μ sin φ with the counterterm), the higher modes
  (k ≥ 2) and the long-time noise (its growth exponent) are measured.

## 26. C0: MV8+ct is the kernel, MV2+ct the reference; the nonlinear response against two predictions (14 Nov, 05:03)

**What the development found (docs/findings.md, "C0 development rows after decision 25"):**
- The linear modes k = 0–3 follow the predicted finite-density expectation of the kernel.
- The nonlinear modes deviate from it as the substrate noise renormalizes the source, μ sin φ → μ e^(−s²/2) sin φ.
- MV8+ct is about 1.6× quieter than MV2+ct and softens the potential less: −0.64 against −0.94 at A 2, t 6.75, ρ 30.
  Its linear systematics against the continuum are larger (0.16 against 0.07 at ρ 30), predicted, and falling with ρ.

**Decision.**
- **MV8+ct is the main kernel of C0; MV2+ct is a reference check.** For C0 the nonlinear behaviour weighs more than
  the accuracy of the linear limit (B4 follows).
- Predictable systematics is not zero systematics. The card keeps both separately: the deviation from the continuum,
  and the agreement with the predicted expectation.
- The nonlinear response is compared with two predictions of different status:
  - **A. The noise-free continuum**, J = μ sin φ: the reference limit without fluctuations of the substrate. The
    difference measures the total finite-density and fluctuation deviation from the classical continuum.
  - **B. The fluctuation-corrected mean response**, J_eff = μ e^(−s²(t,ρ)/2) sin φ with the measured s. Agreement tests
    whether the fluctuations explain the measured softening. It holds for zero-mean Gaussian fluctuations of the phase
    and the corresponding averaging, so the card checks that condition directly: the Gaussianity of the fluctuations,
    their correlation with φ, and ⟨sin φ⟩ against e^(−s²/2) sin⟨φ⟩ event by event.
- The evaluation reports separately:
  - the nonlinear deviation (measurement against A);
  - the renormalization (measurement against B);
  - the residual (the part B does not explain);
  - the density scaling of s and of the renormalization;
  - the time scaling of s and of the renormalization.
- **No tolerance against the continuum may grow with the measured noise** (no rubber ruler). The tolerances are fixed
  in advance from the development rows; the noise-free deviation, the independently computed correction and the
  residual are reported side by side.
- s² ∝ t² ρ^−1/2 is a tested approximation, not a confirmed law (s grows as t^0.85–1.2, so s² as t^1.7–2.4).
- For B4: the large amplitude of an oscillon means a larger change of its effective potential, which does not imply a
  longer life (frequency, profile, radius and radiation channels change). B4 is not tuned on assumed robustness to
  noise; first the effective potential as a function of amplitude, time and density must be known.

## 27. Card C0 revision 1 is frozen as F3.16 (14 Nov, 06:29)

The author started the deciding run of card docs/C0_card_draft.md revision 1 (commit 3488bde). The model, the
readouts, the criteria and the tolerances are those of the card; nothing is changed after the run (R132).

**Test F3.16**: phases/F3_structure/F3_16_green_3p1.py, md5 3de531220105c4c332e23a88b56899fd.

**Helpers:**
- diag/c0_card.py 1f2fc1071acc0e480a6f34a14a5af08f
- diag/c0_apex.py 12c93bbaff2d1ee384848da128c90d4d
- diag/c0_ct.py 8da25cd17f5df94e725441696efda5c2
- diag/c0_ct_renorm.py 9686e935fbd7d924f62f2ca877d305e9
- diag/c0_minvar.py 5b31241ec688e2eceec4bd326a11f63f
- diag/c0_links4.py 0f160d3dc53ef91b6f6c80624dba4a5b

The md5 sums are of the files with LF line ends, as in the repository; a Windows checkout with CRLF has the same
content.

**Tolerances**: diag/c0_card_tol.py (md5 00fa61213a3f889de5b0ae8fdcca1371) from the development rows
diag/out/c0_card_dev_20261114T052221.json, seeds 13000–13003.

**The deciding run**: base seed 0; ρ 30 with 16 sprinklings, ρ 60 with 4; MV8+ct decides, MV2+ct is reported.

**Risk declared before the run**: the fall of the renormalization at A 1 with density (development +0.020 ± 0.013)
may end UNDECIDED.

## 28. C0 revision b: statistical criteria; card F3.16b frozen (14 Nov, 15:49)

**What led here.**
- F3.16 ended INSTRUMENT_FAILED by design errors of its criteria (docs/findings.md).
- The late deviation was a fluctuation: the mean of the linear rule obeys (A3) exactly
  (docs/derivations/c0_mean_exactness.md, diag/c0_late.py).
- A bootstrap of the new development rows (diag/out/c0_card_b_dev_20261114T120302.json, seeds 19000+, at the card's
  numbers of sprinklings) showed that the tolerance 1.25 (|dev| + 2 SE) on about 45 near-zero quantities passes data
  identical to the development in 1 % of cases.

**Decision (the author confirmed the proposal of 14 Nov, ~15:40, and started the run).**
- **Instrument criteria** (predictions exact in expectation: ε of 10 shells against (A1), the normalisation, (A3) of
  the linear rows in every slab) become a significance test: |d| / √(SE² + 0.002²) < t_crit. t_crit is Student t with
  n − 1 degrees of freedom at family-wise α 1 % over all items (Bonferroni).
- **Physics tolerances** keep the rule of 10 Nov (the development value + 2 SE + a margin), with a statistical margin:
  the statistic is the maximum over slabs of |d_j| + 2 SE_j, and the tolerance is its development value + z SE (the
  bootstrap SE of the statistic), with z one-sided at family-wise 5 % over the criteria.
- **Sign criteria** as before; the SE of the renormalization comes from a jackknife over the sprinklings.
- In the bootstrap of the development this form passes 89 %, UNDECIDED 10 %, INSTRUMENT_FAILED 1 % (optimistic, the
  same data).
- This form is the default for later cards with many near-zero criteria, unless the author decides otherwise.

**Frozen.** Card docs/C0_card_b.md.
- **Test F3.16b**: phases/F3_structure/F3_16b_green_3p1_rev.py, md5 a58572c4c3cacc7ef602def7316e1362.
- **Helpers**: diag/c0_card_b.py 2148fa9ace3716a0e613cae488746b0f, diag/c0_fast.py
  1dc98997a7e414a7d7eeb0b7a8b6a5da (bit-identical to diag/c0_card.evolve_card), and the helpers of decision 27.
- **Run**: ρ 30 with 16 sprinklings and ρ 60 with 12; base seed 0 decides.

## 29. Card C1 is frozen as F5.3 (14 Nov, 19:00)

**The author's instruction (14 Nov, 18:43):** "Dost dobrý, na chvilku se vzdálím tak jak bude ten test rovnou ho můžeš
pustit." ("Good, I'll step away for a while; as soon as the test is ready, you can start it straight away.") The run
was started by Claude through the computer-use link (VS Code, Terminal → Run Task "Lab: F5.3").

**What led here.** docs/findings.md, "C1 development rows": on a conformally flat planar well the C0 rule follows its
exact mean (the curved volume along each pair), which is not a local ξR field. The author agreed on 18:39 to fix this
as a property of the substrate with a card (option 1 of 18:3x).

**Frozen.** Card docs/C1_card.md.
- **Test F5.3**: phases/F5_gravity/F5_03_field_curved_order_3p1.py, md5 4ff14c8a75c4dba678b999a669d287bd.
- **Helpers**: diag/c1_conf.py 7f05b01eefdbef0ce2899846fe8ae4eb, diag/c1_apex.py e12a8e24b2983f522abf4ae3b53d31e6,
  diag/c0_fast.py 1dc98997a7e414a7d7eeb0b7a8b6a5da, and the helpers of decisions 27 and 28.
- **Criteria**: the form of decision 28 (significance tests against predictions computed by frozen code; one sign
  criterion).
- **Run**: ρ 30 and 60, 8 pairs each, base seed 0 decides.

## 30. Card F6.8e is frozen; the run is made by Claude in its own workspace (15 Nov, 10:23)

**The author's instruction (15 Nov, 10:23):** "WSL neznám, spusť ho u sebe." ("I don't know WSL, run it on your side.")
Given after the card draft and the test were delivered (10:2x) with the note that the test needs an 80-bit long double,
which numpy on Windows does not have.

**What led here.** The author's order of work (08:30): B (F6.8e, F6.8d), then C, then A, then D. docs/findings.md, "F6.8e
development": the optical sector of F6.8c is not stable beyond r ~ 100 for a known reason; a stable instrument shows the
bulk density growing as r⁵, i.e. an optical age term ∝ (kT)⁵.

**Frozen.** Card docs/F6_8e_card_draft.md.
- **Test F6.8e**: phases/F6_order_dynamics/F6_08e_bd_action_optical_asymptotic.py, md5 4314e2193b0eb0f8fabb989c19529f82
  (LF), with the frozen helpers of F6.8a, F6.8b2 and F6.8c that it imports.
- **Criteria**: as in the card, section 5; the judged pair n 1 536 → 3 072 was not computed in the development.
- **Run**: full run, deterministic (no seed), in Claude's Linux workspace (Python 3.13, numpy 2.5.3, long double eps
  1.1e-19); the output is copied into results/F6.8e unchanged.

- **Revision 1 (15 Nov, 13:2x), the author's choice "Revize před během":** the first full run in Claude's workspace was
  stopped at 175 of 277 tasks by the workspace's two-hour limit per process; nothing of it was printed or saved. Before
  any result, the profile is computed at one wavenumber (k = 0.15; k = 0.25 only in the scaling check I5) with 10 + 6
  nodes in p. Status on α (n 1 536 → 3 072) at k = 0.15. Test md5 9178baf0348e46dd8b74fe758cc9edde (LF); card section 7.

## 31. F6.8e2 and F6.8d are frozen; both run by Claude in its workspace (15 Nov, 21:2x)

**The author's choice (15 Nov, 21:2x):** "Zmrazit a pustit obě" ("freeze and run both"), after the instruction of 21:14
("ano") to prepare the revision of F6.8e and to fix I5 in F6.8d.

**What led here.** F6.8e ended INSTRUMENT_FAILED on I5 (the k⁵ normalisation; docs/findings.md). F6.8d was
reformulated by the author's choice of 13:0x and his remarks of 14:19 (an action density, not a pair count; the
conformal-difference control; Δ diagnostic only; separation of the data).

**Frozen.**
- **F6.8e2**: card docs/F6_8e2_card_draft.md; test phases/F6_order_dynamics/F6_08e2_bd_action_optical_asymptotic_rev.py,
  md5 507378756bb4692d5270b2a10ef08b79 (LF). Judged: α at k = 0.1, n 3 072 → 6 144. The F6.8e record of α (seen) is
  not a confirmation.
- **F6.8d**: card docs/F6_8d_card_draft.md; test phases/F6_order_dynamics/F6_08d_bd_action_mixed_age.py, md5
  9442f6b3b1482afd9422a93f1e1efe13 (LF). Judged: α_ab and α_ab − α_c at k = 0.15, n 1 536 → 3 072. Its form of I5 was
  taken over from F6.8e2 after the F6.8e result was known; no number of F6.8d changed because of the F6.8e exponent.
- **Run**: deterministic, in Claude's Linux workspace (long double), F6.8e2 first; outputs copied into results unchanged.

## 32. F6.8e3, the second revision of F6.8e, is frozen and run by Claude (15 Nov, 22:4x)

**The author's instruction (15 Nov, 22:38):** "dobrá ale test by se měl opravit a spustit poté znovu" ("fine, but the
test should be fixed and then run again"), after F6.8e2 ended INSTRUMENT_FAILED on I3, I5 and I2.

**Frozen.** Card docs/F6_8e3_card_draft.md; test phases/F6_order_dynamics/F6_08e3_bd_action_optical_asymptotic_rev2.py,
md5 b25a87203b72285904e2fe8f25424965 (LF). The instrument checks are stated
in the exponent (levels and the two wavenumbers within 0.01 at the judged pair); the conformal control at n = 384;
k = 0.13 and 0.21; judged α at k = 0.13, n 3 072 → 6 144. The exponent is not blind (F6.8e, F6.8e2 recorded 5.000).

**Run**: deterministic, in Claude's Linux workspace, after F6.8d; output copied into results unchanged.


## 33. F6.8d is not closed as a physical failure (15 Nov, 23:03)

**The author's decision:** F6.8d is not closed as a physical failure. The verdict INSTRUMENT_FAILED stays unchanged.
I2 failed on numerical precision at the largest age (two quadrature levels 1.7e-3 against 1e-3 at n = 3 072), not on
the physical hypothesis. The exponents stay undecided. The values up to n = 1 536 (α_ab 3.0000, α_ab − α_c 1.9999) are
supporting development data, not a substitute for the deciding check.

## 34. K.18 and F6.8f are frozen (17 Nov, 08:22)

**The author's decision (17 Nov, 08:22):**
1. Brute force: the tolerance 1e-3 applies to the extrapolated limit kappa -> infinity, not to the raw finite b.
2. The extrapolation uses kappa = 10, 20, 40 sqrt(eps), i.e. b = 0.01, 0.0025, 0.000625.
3. F6.8f: A at gamma_T(2 -> 4) <= 2.5, B at gamma_T(2 -> 4) >= 5 or a confirmed sign change of M. A also requires a
   nonnegative exponent that decreases over the three intervals.
4. C: a result for young / weak fields, without a decision on the static limit of an old field.
5. R and N are recorded, with the separate test |R(4) - 1| <= 0.05.
Then: "Poté zmraž karty" ("then freeze the cards").

**What led here.** The order of 16 Nov 22:49 (K calibration, K_eff assembly, card F6.8f); docs/derivations/f68f_dev.md,
steps 1b-2. The F6.8f criteria were proposed on 16 Nov 23:40 UTC, before any K_eff development row was seen. The
author decided them after the development rows (R(0.5) = 1.014, R(3) = 1.69 at reduced resolution) had been reported.

**Operational forms written by Claude.** "Confirmed sign change": R < -1e-3 with both t quadratures (G and G - 1). The
separate null test is reported next to the reading and does not enter the status. Under the nonnegative exponent a
K_eff that peaks inside eps <= 4 and then falls reads C.

**Frozen.** md5 with LF line endings.
- **Card K.18:** docs/K18_card_draft.md 60a3fd338f9424171fc0920fbf54b7ab.
- **Test K.18:** phases/K_instruments/K_18_needle_volume.py cb76f2c286753db1934bbfd6ef1ad10d.
- **Card F6.8f:** docs/F6_8f_card_draft.md 3ca658995cd056ca1e0320d8318f2e8d.
- **Test F6.8f:** phases/F6_order_dynamics/F6_08f_bd_action_optical_finite_amplitude.py
  60a229d79b94882a185efa1dde22bb75.
- **Helpers:**
  - diag/f68f_1b.py 516c661d1ec0d67093fbca16b9b9d8d7
  - diag/f68f_1c.py a13c26cc6c497bc33eea6fae903eae1c
  - diag/f68f_1b_exact.py 691e97f189ca2c44b4172962c502b700
  - diag/f68f_1d.py 8d23537b1010d0e32b14dc4444e551a8
  - diag/f68f_layer.py 13611e43d9d51bd98c3ea50817d3e82b
  - diag/f68f_keff.py 17b8452e729f431ff2c6973976b9ec98
  - diag/f68f_common.py ed25abf95e9a56ac197ca3c6c9c3417a
  - and the frozen tests F6.8e3, F6.8a, F6.8b2 they import.
- **Quick runs (mechanics only):** K.18 PASSED; F6.8f INSTRUMENT_FAILED at quick resolution (8 phases), as its
  checks should give.

**Run:**
- K.18 first. Part e3 needs an 80-bit long double, so it runs under WSL.
- F6.8f after K.18 has PASSED.
- Full runs decide; base seed 0 (deterministic).

## 35. K.18: the brute-force window of diag/f68f_1d.py corrected, K.18 re-frozen and run again (17 Nov)

**The author's decision (17 Nov):** correct the instrument and repeat the test.

**What led here.** The full K.18 run (Claude's workspace, diag/run_resumable.py, code as frozen in decision 34) stopped
in part brute without a record: 9 of the 48 brute-force jobs raised "X window outside the tables" (p = 4, 8 at
kappa = 10, 20 sqrt(eps); diag/out/K18_brute_missing.txt). The quick run has only the point eps 4, p 2 and did not
reach them. The Delta z windows of the tables were built on l = rho, while a query point has
rho_x^2 = (sa l - CO X)^2 + Y^2; at large CO = p / kappa the queries overshoot the window by 1-8 % at rho -> 0 and
rho -> 1 (docs/derivations/f68f_dev.md, section 8, "Window coverage at large p/kappa").

**The correction (Claude).** At every rho node the window is the union of the previous band and the exact Delta z range
of the queries that use the node (X grids of the l nodes, Y^2 <= 8 sigma^2 l (1 - l), T_xy), widened by 5 %.
Criteria, parameters, levels and points of K.18 are unchanged; the test file and the card are unchanged.

**Verification** (diag/out/f68f_1d_cover_validation.txt):
- all 48 brute-force jobs run;
- the 39 jobs the previous version computed change by <= 2.5e-9;
- flat space (eps 1e-9, kappa 20, q 1/4): A - 1 = 3.1710e-6, 3.1715e-6, 3.1761e-6 at p = 0, 2, 8 against 3.1764e-6
  with exact distances; the previous version fails at p = 8;
- K.18 quick run PASSED (results/K.18/quick/20261117T141241Z_main_seed0_f7513da2.json, gitignored).

**Seen before the correction.** The raw volumes A of the brute-force jobs the previous version computed, and of the 9
jobs with the corrected windows. Their comparison with the q^2 series and the extrapolation were not computed (part
brute assembles them only at the end of the run).

**Frozen.** md5 with LF line endings. Everything as in decision 34 except
- diag/f68f_1d.py 037bcf568241655c4bf2370fdfc9a56e (was 8d23537b1010d0e32b14dc4444e551a8).

**Run:** the full K.18 run is repeated with the corrected instrument; the cache of diag/run_resumable.py is keyed by the
md5 of the code, so nothing from the stopped run is reused. F6.8f after K.18 has PASSED (F6.8f does not import
diag/f68f_1d.py).

**Result (17 Nov).** results/K.18/20261117T182826Z_main_seed0_f7513da2.json (code md5 as frozen above; Claude's
workspace, diag/run_resumable.py in three passes, the workspace stops a process after 2 h): **FAILED**.
- Instrument parts pass: single (one ray per phase); jacobi delta_J 7.2e-14, series A_0 vs plane wave 3.8e-12; series vs
  exact needle 3.4e-8, window failures 0.
- qconv 1.8e-6 (<= 1e-3) and e3 +5.5e-4 (<= 5e-3) pass. brute: largest |A_ext - A_series| 8.9e-7 (<= 1e-3).
- The O(1 / kappa^2) form fails at 2 counted ratios, both at eps 4, p 8, theta 4, first ratio (kappa 10 -> 20):
  3.133 (q 1/4) and 3.103 (q 1/8); level differences -2.43e-4, -7.8e-5, -2.1e-5. The second ratios (3.725, 3.857) have
  numerators below 1e-4 and are not counted. All other counted ratios lie in 3.958-4.054. At eps 2, p 2, theta 4 the
  ratios are 4.86-5.50 with numerators below 1e-4 (not counted).
- The point eps 4, p 8 was not run with three levels in development (step 1d: eps 4 at p 2 and 4).
- By this decision F6.8f counts only after K.18 has PASSED; with this record the condition is not met.

## 36. K.18b: K.18 revised with four brute-force levels; the F6.8f data usable for development (17 Nov)

**The author's decision (17 Nov):**
- The F6.8f run already started (author's machine) runs to the end. Its data may be used for further research as
  development (decides nothing). Its record is not counted while no calibration of its instrument has PASSED
  (decisions 34, 35).
- K.18 is revised as K.18b with the criterion proposed by Claude, frozen before the result of F6.8f is known. If K.18b
  PASSES, the F6.8f record counts; F6.8f does not import diag/f68f_1d.py, the only code of the failed part.

**What led here.** K.18 FAILED only in the form check of part brute (decision 35, Result). A fit of the three level
differences to a/kappa^2 + c/kappa^4 + e/kappa^6 (kappa = kT, b = eps/kappa^2): |c| ~ 1e-5 to 2e-4 at every point,
while a varies by two orders and changes sign in (p, theta); at eps 4, p 8, theta 4, a = -3.4e-4 and c/a = -0.39
(q 1/4), -0.16 (q 1/8). The O(b^2) term is not negligible at kappa = 10 sqrt(eps) where the O(b) coefficient is small;
the absolute cut 1e-4 of the ratio check did not account for it. The three-level extrapolation removes both terms.

**K.18b** (phases/K_instruments/K_18b_needle_volume_four_levels.py). Only part brute differs from K.18:
- levels kappa = 10, 20, 40, 80 sqrt(eps) at the same 8 points and q = 1/4, 1/8;
- |A_ext(20, 40, 80) - A_series| <= 1e-3 (the author's criterion) and |A_ext(20, 40, 80) - A_ext(10, 20, 40)| <= 1e-5
  [proposal; the O(b^3) remainder e/64 of the coarser extrapolation is <= 9e-7 at every K.18 point];
- the ratio check of K.18 is removed; level differences and ratios are reported, not judged.
All other parts, parameters, points and the verdict logic are as in K.18.

**Seen before freezing.** The K.18 record and the fit above; the value of one fourth-level job (eps 2, p 0, theta 2,
q 1/4; cost 10 s); all 16 fourth-level jobs run to completion (no window or coverage error), values not read; quick run
PASSED mechanically (results/K.18b/quick/20261117T185257Z_main_seed0_0a03f218.json, point eps 4, p 2, theta 2, not a
judged point).

**Frozen.** md5 with LF line endings: phases/K_instruments/K_18b_needle_volume_four_levels.py
159785ec55f67294515405332325ce21; the diag files and F6.8e3 rev2 as in decision 35.

**Run:** Claude's workspace (diag/run_resumable.py; part e3 needs an 80-bit long double).

**Result (18 Nov).** K.18b full run PASSED: results/K.18b/20261117T230319Z_main_seed0_0a03f218.{json,txt} (test file md5
159785ec55f67294515405332325ce21, LF; 1671 s in the final pass, resumed from the cache over four passes after two
workspace restarts). Parts: single 1 ray per phase; jacobi 7.24e-14; series A0 vs plane wave 3.77e-12; series vs exact
needle 3.36e-8, 0 window failures; qconv 1.80e-6; e3 +5.51e-4; brute |A_ext(20, 40, 80) - A_series| <= 1.07e-6 (<= 1e-3),
|A_ext(20, 40, 80) - A_ext(10, 20, 40)| <= 1.44e-6 (<= 1e-5). At eps 4, p 8, theta 4 the level ratios are 3.13, 3.73, 3.61
(q 1/4) and 3.10, 3.86, 3.52 (q 1/8), reported, not judged.

With K.18b PASSED, the F6.8f record counts (this decision): results/F6.8f/20261117T205949Z_main_seed0_c3454e9d.{json,txt},
FAILED, reading B: R(eps) = 1.0034, 1.0564, 1.2482, 2.5178 at eps = 0.25, 1, 2, 4; gamma_T = 5.07, 5.48, 7.02 (>= 5);
separate null test |R(4) - 1| = 1.518 FAILED (not in the status); instrument checks within tolerance.

## 37. F5.2c2: F5.2c revised with the flat control on 32 sprinklings; replications of F5.2c and F5.2d (18 Nov)

**The author's decision (18 Nov):**
- F5.2c is revised as F5.2c2 with one change, the flat control on n_ctrl = 32 sprinklings instead of 8, and the author
  runs it after F5.1e seed 1 (seed 0 decides, seed 1 replicates).
- F5.2d is replicated unchanged with seed 1 (test sha256 decef7cc).
- Order of the work on gravity in 3+1 (direction C): passive fall of walls of the C0 field (C2), then the active side
  (a 3+1 analogue of F5.2d), then two bodies (a 3+1 analogue of F5.2b). Calibration of the dynamics of the order
  (F6.4a FAILED, F6.4c UNDECIDED) separately.

**What led here.** F5.2c seed 1 (16 Nov) was INSTRUMENT_FAILED on the flat control alone: D_0 = -0.152 +- 0.091,
|D_0| + 2 SE = 0.335 against 0.25, 1.7 SE from zero; all physics conditions were met (D_K -0.007 +- 0.021, dv_L
+0.0038 +- 0.0035, r_B 0.943 +- 0.039). The control ran on 8 sprinklings, the curved arms on 32. With a per-sprinkling
sd of D_0 of about 0.17-0.26 (seeds 0 and 1), 2 SE at n = 8 is 0.12-0.18 of the tolerance 0.25: an instrument with no
drift fails the check in roughly one run of four or five. A design error of the card (too little power in the
control), not a property of the instrument. At n = 32, 2 SE is about 0.08.

**F5.2c2** (phases/F5_gravity/F5_02c2_walls_curved_2p1_rev.py): differs from F5.2c only in n_ctrl 32, the id, the
seed labels "F5.2c2" (new sprinklings; no data of F5.2c reused) and the cache folder. Model, arms, readouts, tolerances
and verdict logic unchanged. Quick run (rho 12, 2 + 1 sprinklings) ran mechanically in Claude's workspace (decides
nothing). The verdicts of F5.2c stay: seed 0 PASSED, seed 1 INSTRUMENT_FAILED (R132).

**Frozen.** md5 with LF line endings: phases/F5_gravity/F5_02c2_walls_curved_2p1_rev.py
82aa4955fccc6c1ce43721361a8357bf. Run: python run.py F5.2c2 --workers 8 (author's machine; ~32 core-hours).

**Result (18 Nov).** F5.2d seed 1 (unchanged, decef7cc): PASSED, r_L 0.965 +- 0.006 (load 16) and 0.946 +- 0.013
(load 8), u +0.020 +- 0.008; the deficit of the potential in the last slab is replicated (findings, F5.2d).
F5.2c2 (19 Nov; author's machine, test sha256 d651784d, frozen md5 82aa4955 unchanged, 29 361 s; results/F5.2c2/20261118T234947Z_all_seed0_d651784d.json): **PASSED**.
Instrument OK, flat control D_0 -0.026 +- 0.035 (|D_0| + 2 SE 0.097 / 0.25); D_K +0.020 +- 0.019, dv_L +0.0043 +- 0.0035,
slope - 0.3 +0.0174 +- 0.0035, r_B 0.952 +- 0.041 (findings, F5.2c). The seed-0 result of F5.2c is replicated.

## 38. OV1 (variations of the order) and F6.9a, the vacuum control in 2D (19-20 Nov)

**The author's decisions (19 Nov):**
- Route A of the roadmap: the question is the rule by which matter selects the admissible, stable order (M -> H_M ->
  Pi_M = P[H_M]), not a field equation for Pi. Steps OV1 (admissible variations of H and the variational rule of
  S_BD, no matter), OV2 (the flow F_v, then S_m[H, F]), OV3 (finite-amplitude dynamics with a source). The names
  OV1-OV3 replace A1-A3 of the first draft; the open decisions A1, A2 of variational_audit.md and A2-B are unrelated.
- OV1 (docs/derivations/order_variation_OV1.md): level L3 (the individual action and the ensemble mean, separate
  verdicts); stationarity S-b (centred and second differences; S-a rejected, S-c a secondary diagnostic); eps_0 = 0.12
  with 0.06 and 0.24 in 2D (the F6.4a/c value; not transferred to other dimensions; F6.8 has only eps rho = 1);
  nested Alexandrov diamonds with a common future tip, T/w = 8, 16, 32; w and n_loc from 99 % of the absolute kernel
  weight; family I = all reversible local rewrites keeping the partial order (in 2D: adjacent transpositions of 2D
  orders); family M = the frozen basis of F6.8 (3+1); dimensions 2D -> 2+1 -> 3+1; the flat-vacuum control c(T) = 0
  is blocking.
- Phase: bands from flat sprinklings at the actual N, calibration seeds (>= 100000) disjoint from validation seeds;
  the escape is probed by an F6.4-type sampler with the local moves of family I, labelled an exploration of the
  configuration space (not OV3); size scan N = 50, 100, 200, 400 at eps 0.12; stage-B rule: the stage-A neighbours
  where the fraction outside the flat h band crosses 0.5 (the rule "maximum of var(S)" failed at N = 100 and stays on
  record).

**Development (decides nothing; OV1 sec. 10, 13).** Validation on the author's machine, c/(4 eps) per element:
(+2.2 +- 5.4)e-4, (-1.7 +- 2.5)e-4, (+0.9 +- 1.7)e-4 at eps 0.24 / 0.12 / 0.06 (512 / 512 / 256 sprinklings); the
constant <S> at T/w = 8 is 3.1 +- 3.1, 6.0 +- 1.5, 2.8 +- 1.0; balance <D S> zero; the calibration bands carry over
(0.4-2.0 % of validation sprinklings outside, nominal 1 %). Escape point of the manifold-like phase under local moves
(eps 0.12, 3500 sweeps): beta_e = 2.42, 1.26, 0.97, 0.83 at N = 50, 100, 200, 400; the escape is chain by chain
(first-order-like metastability); the crystalline side is not equilibrated.

**F6.9a (20 Nov).** ID F6.9a; tolerance and numbers of sprinklings as proposed; the constant term is reported, not
judged. Card docs/OV1_vacuum_2d_card_draft.md, revision 2 (I3 at most 5 % outside the bands, I4 on all sprinklings,
seeds by derive_seed). Criterion per eps: PASSED if |c/(4 eps)| + 2 SE <= 1e-3, FAILED if |c/(4 eps)| - 2 SE > 1e-3,
else UNDECIDED; 5400 / 1200 / 280 sprinklings at eps 0.24 / 0.12 / 0.06. Quick run (Claude's container) mechanical
OK, I1 to 1.1e-10.

**Frozen.** md5 with LF line endings: phases/F6_order_dynamics/F6_09a_bd_vacuum_2d_diamonds.py
9a1da81c85b9a6dcbaf6e6dad131e8fd. Run: python run.py F6.9a --workers 8 (author's machine, started 20 Nov ~10:20;
6885 tasks, ~40 core-hours).

**Result (20 Nov).** F6.9a seed 0 (author's machine, test sha256 7929a6b8, frozen md5 9a1da81c unchanged, 7855 s;
results/F6.9a/20261120T113217Z_main_seed0_7929a6b8.json): **PASSED**. c/(4 eps) = (+0.5 +- 1.7)e-4, (-0.7 +- 1.6)e-4, (+0.6 +- 1.5)e-4 at eps
0.24 / 0.12 / 0.06 (|x| + 2 SE <= 4.0e-4 against 1e-3); instrument I1-I4 OK (findings, F6.9a).

**Result (21 Nov).** F6.9a seed 1 (author's machine, 12 735 s; results/F6.9a/20261121T013343Z_main_seed1_98240566.json): **PASSED (replicated)**.
c/(4 eps) = (+0.15 +- 1.70)e-4, (-0.50 +- 1.63)e-4, (-0.61 +- 1.64)e-4 at eps 0.24 / 0.12 / 0.06 (|x| + 2 SE <= 3.9e-4
against 1e-3); instrument I1-I4 OK. The file had become CRLF on the author's machine (raw sha256 98240566); its
LF-normalised sha256 is 7929a6b8 and its LF md5 9a1da81c, unchanged (findings, F6.9a).

## 39. OV2: the choices of the draft specification (20 Nov)

The author's decisions on docs/derivations/order_variation_OV2.md, section 8:
1. Both candidates of the matter action: S_m^phi (the action of the field C0 on the order) and S_m^tau (the proper
   time of a flow counted by the order), with a comparison.
2. The kinetic operator in difference form (minimal coupling; the plain form couples R phi^2 and breaks the 2 pi
   periodicity), checked by T0.
3. S_m as the coupling functional at fixed phi (T_ab is defined off shell); the conservation under the retarded C0
   dynamics is the separate test T3.
4. Order of the tests: T0 -> T1 -> F1/F2 in 2D, then T2 and T3 in 3+1.
5. Branch F (the flow F_v) is pursued in that order, after T1.

## 40. OV2: the cards T1 and T2 frozen (21 Nov, 09:17)

The author froze the cards of docs/OV2_coupling_cards_draft.md (revision 3) with Claude's recommendations of the
morning summary:
1. IDs F6.10a (T1, 2D, conformal) and F6.10b (T2, 3+1, the omega and nu sectors); the tolerances and the fresh
   parameters of the drafts are accepted.
2. T2 is one card for both sectors.
3. The sprinkled cross-check of T1 (I4) is not required and is not part of the test.

Test files (md5, LF): phases/F6_order_dynamics/F6_10a_matter_conformal_2d.py a53d3a3188eefef7173b80f383615223 and
phases/F6_order_dynamics/F6_10b_matter_metric_3p1.py bdae39b0e674bc190f0fb64eb4229254. Only the header text changed from
the drafts (syntax trees without string constants identical); quick runs PASSED (F6.10a 21 Nov 08:18 UTC; F6.10b 01:12
UTC on the draft code).
Disclosure (card T1): the grid check computed R at the fresh parameters before freezing (-0.014 and -0.005); the
criterion had been written before, and the quadrature is deterministic.
Still open: the F ontology (F-a or F-b), the OV3 decisions (docs/derivations/order_variation_OV3.md, section 7), T3, and
the sha alias of F6.9a after its CRLF conversion.

**Result (21 Nov).** F6.10b (author's machine, test sha256 14801b85 = the frozen file bdae39b0, 2955 s;
results/F6.10b/20261121T091000Z_main_seed0_14801b85.json): **PASSED**. Ratio to the continuum at k = 0: omega 0.9965, nu 1.0060;
at k l_k = 0.03 / 0.06 / 0.12: omega 0.985 / 0.955 / 0.855, nu 0.990 / 0.950 / 0.818; instrument I1-I4 OK (findings,
F6.10b). F6.10a ran on the draft file (sha256 514b1ecf; the frozen file was delivered during the run, flag
files_changed_during_run): PASSED with R = -0.0139 / -0.0047, the same code; it is re-run on the frozen file for the
record.
F6.10a re-run on the frozen file (test sha256 f1fd6835 = a53d3a31, unchanged during the run, 183 s;
results/F6.10a/20261121T095829Z_main_seed0_f1fd6835.json): **PASSED**, R = -0.0139 (N 66 560) and -0.0047 (N 266 240),
instrument I1-I3 OK (findings, F6.10a).

## 41. The F ontology, the T3 pilot, OV3-lin-a, K_c(k), the wall card and the selection principle (21 Nov, 11:01)

1. **The F ontology.**
   - F-a, the energy and momentum derived from S_m (the corresponding stress-energy current), is the main physical
     hypothesis.
   - F-b, the topological charge, is an independent control of the measuring procedure.
   - The third option (e.g. the throughput, the number of causal steps) is not chosen now. It stays an alternative
     hypothesis in case F-a fails the conservation law. It is not a candidate for F until it has its own definition of
     a balance: a throughput measures the capacity to transmit, not automatically a conserved extensive quantity.
   - **Two separate questions** (to be kept apart in the documentation):
     1. Are the energy and the momentum of F-a conserved between cuts, including the contributions of the lateral
        boundaries? This is the test of F-a.
     2. Can F-a be decomposed into non-negative contributions on the causal links so that the local balance
        In_v + S_v = Out_v + ΔQ_v holds? This decides whether F-a can be directly the flow of the axiom of continuity.
     - If (1) passes and (2) fails, the source T_ab has not failed, but the axiom F needs another or a more general
       expression.
   - F-a is naturally nonlocal: E(A) = eps^2 Σ over pairs crossing A of f(n)(Δφ)^2 + E_pot(A) includes the pairs that
     cross the antichain.
   - F-b checks the choice of the cuts, the antichains, the orientation and the topological readout, but not the
     implementation of the weights f(n)(Δφ)^2. For those, small exact synthetic tests on hand-built posets are added.
2. **T3: development pilot approved.**
   - About 1e4 events, rho = 100-200, a few tens of independent sprinklings.
   - Diagnostic: E_x = ∂S_m/∂φ_x and B_ξ = Σ_{x in slab} E_x (£_ξ φ)_x with ξ = ∂_t. B_ξ is compared with the change
     of the energy between two antichain cuts, including the flux through the lateral boundary of the region.
   - The exact sign and normalisation of the identity are derived from the definition of S_m used and frozen before
     the evaluation.
   - The identity assumes the translation symmetry of the action. For a nonlocal pairwise action the nonlocal and
     boundary contributions must be included. On a single sprinkling translation invariance is not exact, so the
     scatter between sprinklings is part of the result, not a numerical nuisance.
   - The pilot determines:
     - the distribution of B_ξ over independent sprinklings;
     - the difference between the smooth reference kink and the noisy C0 field;
     - the dependence of the residual on the density and on the size of the region;
     - the number of sprinklings needed for a preset precision, n ≈ (1.96 s_B / δ)^2, with δ the half-width of the
       95 % interval. The independent statistical unit is the whole sprinkling, not a pair or an event.
   - The 1 % target is an estimated required precision, not automatically a criterion of physical success.
     - Hundreds of seeds are still a realistic validation path. Millions mean another formulation is needed.
     - Smoothing on l_k would have to be a separately defined variant of the test, not a later correction of the
       result.
   - F-b is evaluated in the pilot on the same cuts as the topological control, with a separate readout.
3. **OV3-lin-a: computed now in the young-field reading.** The output is a conditional linear response, not the
   propagator of an old static body. The card and the findings must state:
   - that the local part of the stiffness is used;
   - that validity refers to a young field, or a field weak enough for the age term to be negligible;
   - that the problem of finite amplitude and age of F6.8f stays open;
   - that the result by itself closes neither C.1 nor the 1/r law.
   - It stays a result for the bare ensemble action until the entropic correction to the second variation is
     evaluated (item 6).
4. **The source of K_c(k): option A.**
   - A new F6.8c-type computation for k l_k = 0.02 to 0.12, with the local part extracted by the same procedure as for
     K in F6.8a.
   - F6.8c2's C(0) = 2.002 at the preregistered kT = 8 contains the age contribution; it is not a separate measurement
     of the local K_c(k).
   - The card reports separately: (1) the total measured optical stiffness; (2) the age component; (3) the local
     part; (4) the sensitivity of this separation to the chosen ranges of k and T.
   - The covariance relation C = 2 + 4 φ_0 + φ_2 is an independent check. φ_2 must actually be measured and may not
     be computed from the same C it is to check: ideally a separate computation or a separate instrumental procedure.
5. **The repulsion of the BPS wall: a separate card, approved.**
   - For the static BPS wall: rho = rho_E, T^a_a = -3 rho_E, Σp = -2 rho_E; the lapse source is
     2σ - τ = -4 rho_E + 3 rho_E = -rho_E, opposite in sign to dust. In the local limit Ψ/Φ = -1 for the wall and +1
     for dust. At k l_k = 0.12 the form factors give -4(0.818) + 3(0.855) = -0.707, about -0.71 rho_E.
   - Three separately preregistered evaluations:
     - the sign of the resulting response;
     - the size of the deviation from -rho_E;
     - the ratio of the response of the wall to that of dust at the same k.
   - The amplification of the error when two large contributions are subtracted must be computed with their
     covariance: σ_src^2 = 16 σ_ν^2 + 9 σ_ω^2 - 24 Cov(F_ν, F_ω), not by multiplying one error by four.
6. **The selection principle (a revision of the author's earlier choice).**
   - The stationary branch of an individual action is not enough. OV1 found locally reachable configurations with an
     action lower by three orders, and the manifold-like phase is held by the statistical weight of the orders. There
     is no reason to assume that the physical order is a stationary point of the individual action.
   - For OV3-lin-a: option 1, the saddle of the ensemble mean action <S>, conditioned on the manifold-like phase.
   - For nonlinear OV3, the working principle is option 2, an explicit balance of action and entropy. An effective
     functional of a macrostate M is introduced, Γ_β(M) = -(1/β) ln Σ_{H in M} e^{-β S_tot[H]}. Its exact form and
     measure belong to a separate specification. The formula shows the intended structure: the weight depends on the
     action and on the number of microscopic orders that realise the same macrostate.
   - Option 3, a stochastic dynamics weighted by the action, is not chosen now. It would need derived transition
     probabilities and a time dynamics of the configuration space.
   - Entropy may change the second variation. Whether the entropic term changes the effective stiffness already at
     quadratic order must be measured. OV3-lin-a therefore stays a result for the bare ensemble action until this
     correction is evaluated. If it is significant, the result cannot be presented as the physical linear propagator
     without the correction.
   - This avoids two wrong extremes: minimising the individual action, which OV1 refuted, and introducing a complete
     stochastic dynamics that is not yet defined.

## 42. T3: the readout of the C0 field approved; OV2 section 14.2-14.4 fixed for stage 2 (21 Nov, 11:42)

- **Background.** Section 14 corrected section 11: the crossing-pair sum E(A) has zero expectation, and the energy of
  F-a is the Noether charge H = -Q of the windowed identity <B_chi> = <R_chi> (exact in expectation). The decision-41
  formula E(A) is replaced by this H.
- **The readout, approved by the author.** For a C0 field, (L_xi phi)_x = d_t Phi_fit(x), where Phi_fit is the kink
  fitted to the C0 run (centre and velocity, as in F3.14b). The ramp terms are evaluated on the C0 values.
- **What is fixed.** Claude takes this approval as fixing OV2 section 14.2-14.4 for stage 2 of the pilot: the
  identity, the sign <B> = -(H(t2) - H(t1)), the normalisation H = -Q and the readout. The author may narrow this.
- **Status.** Stage 2 stays development (decision 41). The 1 % target is an estimated required precision, not a
  criterion.

## 43. The card F6.8g (the local optical stiffness K_c(k)) approved as drafted (21 Nov, 12:15)

On docs/F6_8g_card_draft.md, revision 0:
1. **The ID F6.8g and the rows: approved.** k = 0.02, 0.03, 0.045, 0.06, 0.08, 0.10, 0.12, at the age ladder
   T in {20, 28, 40, 57, 80, 113, 160} with n = kT <= 12. T = 10 and 14 are reported only.
2. **The model and the criterion: approved.**
   - The model of the age component: K_hat_s = K_loc,s(k) - k^2 F_s(kT), with F_s(n) = a_1 n + a_2 n^2 (J = 2) and
     F_s(0) = 0, fitted on the rows with n <= 8.
   - The criterion of the sensitivity spread: <= 0.005 at every k <= 0.1.
3. **F6.8h (the TT shear, an independent phi_2) is not developed now.**

The test file is implemented from the card. As the laboratory standard requires, every judged point is exercised at
full resolution before freezing. The card is frozen by the author.

## 44. F6.8g revised: the young-field stiffness read at fixed n = kT = 4 for all sectors (21 Nov, 13:44)

**What the quick run and the probe of the draft of decision 43 found:**
- The k = 0 null contains a young-age term that does not depend on k. It falls as T^-4: J_bb = 8.2e-5, 2.1e-5,
  5.1e-6 and 1.3e-6 at T = 20, 28, 40 and 57. In units of K_hat_bb it grows as 1/k^2, so a row is clean only for
  n T >~ 100. The limit n -> 0 is out of numerical reach.
- The conformal sector at small n differs from F6.8a's local part by about 1 k^2, and the offset scales with n
  (F_c = 0.92 and 0.97 at n = 1.2 for k = 0.06 and 0.03). F6.8a's local part is the intercept of the asymptotic linear
  age law (n >~ 8), not the limit n -> 0. The local part is therefore not unique at order k^2. Instrument I1 of the
  draft compared two different definitions; the design error was Claude's.

**The author's choice (option B):**
- The young-field stiffness is read at fixed n_y = kT = 4 for all sectors (the field of age 4/k).
- Rows at n = 3 and 6 are added for the n-dependence.
- Instrument I1 becomes the known answer at fixed n: the conformal intercept K0_c = 1 within 1e-3, as I5 of F6.8c2.
- The separation into a local and an age part, its sensitivity spread and the criterion 0.005 of decision 43 lapse
  with the reading they belonged to.
- The mixed sector at the same rows is included if the code of F6.8b2 allows it.
- OV3-lin-a reads all inputs at n = 4. The source form factors F_omega and F_nu of F6.10b do not depend on the age
  (T 80 against 120: below 0.1 %).

## 45. Test IDs: the laboratory's numbers are kept (21 Nov, 15:07)

The "next tests" table of the published paper (paper/chapters/07_predictions.tex, l. 237-243) reserved F6.8g for the
optical age term of old fields and F6.9 for the age term in the full static propagator. The laboratory has since used
F6.9a for the OV1 vacuum control (decision 38) and F6.8g for the young-field stiffness card (decisions 43, 44). The
author keeps the laboratory's numbers. The planned tests of the paper's table are renumbered in the next version of the
paper. Whether that version includes F6.8g and OV3-lin-a depends on how today's runs end.

## 46. F6.8g revision 2 frozen: judged rows k = 0.045 - 0.12 (21 Nov, 15:48)

**What led here.** Two development runs of revision 1 (decision 44) were made before freezing, and they are not
verdicts: Claude's workspace (results/F6.8g/20261121T143903Z) and the author's machine (results/F6.8g/20261121T143802Z),
both with sha 6bfb3471. Instruments I1-I4 passed in both. Two findings:
- **The old optical rows are not numerically robust at small k.** The two platforms (Linux, numpy 2.5; Windows,
  numpy 2.2) differ in K_hat_bb by 2.7e-4 at (0.03, T 133), 1.4e-4 at (0.045, T 133) and 5e-5 at (0.06, T 100), but by
  <= 6e-6 at T <= 89. J_bb is a small remainder of a large cancellation, and it is divided by k^2.
- **The row (0.03, T 133) lies about 3.5e-3 below the n-dependence of every other k.** This gave a fit residual of 1e-3
  at n = 4, against 5e-5 at n = 3 and 6. A level-2 probe of this row was stopped unfinished, so its cause is not
  settled.

**The author's choice:** "Pojďme cestou různých k. Ano zmraž to." The judged rows are k = 0.045, 0.06, 0.08, 0.10 and
0.12 at n = 4 (T = 89 ... 33); k = 0.03 is dropped. I3 is checked at (0.045, n 4) and (0.12, n 4). The inputs of OV3-lin-a
are taken at k = 0.06 and 0.12.

**Exercised before freezing.**
- The judged rows, from the development runs: K0_bb(n 4) = 0.99944 (Claude) and 0.99943 (author), fit residual
  4.7e-5; K0_c = 0.99993; K0_ab = 0.99995.
- The new I3 point at level 2: K_hat_bb 0.961251 against 0.961281 at level 1 (3.0e-5).
- These numbers were seen before freezing. The criterion (|K0_bb - 1| <= 0.01) was set in decision 44, before them.

**Frozen.** md5 with LF line endings:
- test phases/F6_order_dynamics/F6_08g_bd_action_optical_local.py: d594c74386a690e4e978fb4bb65debad;
- card docs/F6_8g_card_draft.md (rev. 2): 00c780d57b04e7259e4aebb5a76ebbd8.

**Result (21 Nov).** F6.8g, run by the author on the frozen file (sha256 86ec1e05, LF md5 d594c743, unchanged during the
run), 8 workers, 4091 s; results/F6.8g/20261121T183904Z_main_seed0_86ec1e05.json: **PASSED**.
- K0_bb(n 4) = 0.99943, fit residual 4.9e-5. I1 K0_c = 0.99993; I2 OK; I3 2.1e-5 and 2.5e-6; I4 <= 1.2e-3.
- A first start on 16 workers stopped at task 28/30 on a memory error, without a record. It was restarted unchanged.
- OV3-lin-a from this record and F6.10b is recorded in findings as a derived consequence (the author, 21 Nov 17:23):
  - dust: Phi/Phi_GR 0.940 / 0.809 and gamma 1.033 / 1.122 at k l_k 0.06 / 0.12;
  - the wall repels: wall/dust -1.045 / -1.062.

## 47. The next version of the paper: planned-test numbers F6.8i and F6.11; the DOI later; new figures (21 Nov, 20:56)

The paper's new version (chapters 00-08, A, B; commit bc8eb3a on claude/c0-minvar) renumbers the planned tests that
the previous version had called F6.8g and F6.9 (decision 45). Claude proposed the numbers; the author: "ano".
- **F6.8i**: the optical age term of old fields, eps = b (kT)^2 >> pi^2 (the previous version's F6.8g).
- **F6.11**: the age term in the full static propagator (the previous version's F6.9).
- F6.8h stays reserved for the transverse traceless shear (phi_2 measured independently, decision 43).

The DOI in paper/main.tex (author block, Data availability) is left as it is; the author sets it when the new Zenodo
version is created ("to upravím potom").

The author asked for something new that describes or visualises exactly what the programme does ("hodilo by se mít
něco nového a přesně popisující či vizualizující co děláme"). Claude added three figures to paper/figures/make_figures.py,
every number read from the records:
- fig_defect (Part II): inertia (F3.14b, F3.15), swerves (F3.14c), fall under the growth rule (F5.2b);
- fig_chain (Part V): a diagram of the chain of the field equation, with the tests and the status of each step;
- fig_fieldeq (Part V): stiffness at kT = 4 (F6.8g), source form factors (F6.10b), and the response computed from them.

## 48. β: the reduction, the perturbative route, stage β1 first, F6.12; the results folder by phase (21 Nov, 21:43)

**β (OV3 section 10; card docs/F6_12_card_draft.md, draft 0).** The author: "1. ANo 2. Poruchová cesta 3. Ano 4. ANo".
1. The reduction is accepted: β − 1 = −κ₃/8, with κ₃ the cubic static action of the order along the exponential metric
   (ω, ν) = (−U, 2U), where the Einstein-Hilbert action has no cubic term.
2. Route (a), perturbative: interval volumes to third order (a new Fermat term T₃ for a stratified medium) and the
   Mecke/kernel expansion to third order; no finite differences in amplitude.
3. Stage β1 first: the pure conformal sector at cubic order (EH: 12 ω(∇ω)², i.e. 9A²Ck² for u = A cos kz + C cos 2kz),
   as development; later instrument I1 of the card.
4. The ID F6.12 and the proposed criterion (PASSED |β₀ − 1| ≤ 0.01; FAILED > 0.03). The card stays a draft: it is
   completed and frozen only after the code exists and every judged point has been exercised.

**The results folder.** The author: "Oprav složku results - ať obsahuje podsložky F1-K s popisem pro daný objective
testů v podsložce - vše v eng".
- The 98 test folders were moved unchanged into `results/<phase folder>/<ID>/`; the phase folders are those of
  `phases/` (K_instruments, F0_substrate … F8_data). Each has a README.md with the objective of the phase and the
  question of every test (from docs/test_list.md). `results/README.md` describes the layout. history_3x/ and
  sha_aliases.json stay at the top. Script: diag/reorg_results.py (moves only, deletes nothing, never overwrites).
- `core/provenance.py` writes to the new layout and still reads the old flat one (`results_dir`, `_record_dirs`).
  This is a change of the core, so records from now on carry a new core hash. No gate is affected that was not
  already closed: F1.9a's K.13a gate compares core hashes, and the core had changed since K.13a (958b2e63) before.
- The frozen test F1.9a keeps its cache at `results/F1.9a/cache/`; a run recreates that folder.
- The F3.14c replication that was running during the move started with the old core in memory. Its record will
  appear in `results/F3.14c/`; `load_all` reads it there, and it can be moved into `results/F3_structure/F3.14c/`.
- docs/findings.md: the paths `results/<ID>` were rewritten to the new layout. The decisions, cards and derivations
  are historical and keep the old paths (`results/README.md` gives the mapping).

## 49. The quantum layer opened as a section of its own (28 Nov, 10:29)

**The author** (after Claude's assessment that the laboratory is still entirely classical and that H3 and H4 point to
the need for a quantum level): "Dobře, než doběhne test připrav mi nové testy - první pro kvantové. Bude to chtít
vlastní results, findings, kartu,... vlastní sekci v labině pro odlišení."

- The quantum layer is opened on 28 Nov, before all classical pillars have verdicts. This replaces the first rule of
  docs/quantum/README.md section 6 (9 Nov: "nothing starts before the classical pillars have their verdicts").
- The danger named then stays, so the first quantum tests are calibrations of quantum instruments on a bare Poisson
  sprinkling of flat spacetime. They depend on no open classical pillar.
- **Its own section, kept apart from the classical one** (docs/quantum/README.md section 7):
  - IDs Q.n;
  - planned tests in docs/quantum/planned_tests.csv;
  - cards in docs/quantum/Q<n>_card_draft.md;
  - findings in docs/quantum/findings.md;
  - test files in phases/Q_quantum/;
  - records in results/Q_quantum/<ID>/.
- The phase Q is added to core/registry.py when the first quantum test file is written. That is a change of the core,
  so records from then on carry a new core hash, as with decision 48.
- Every quantum test has a null model (the classical, non-negative flow, or the same construction with non-negative
  weights). A quantum effect counts only where it differs from the null.
- **The first card, Q.1 (draft 0):** Johnston's sum-over-chains propagator on a 1+1 sprinkling. The mean ½J₀(mτ) is
  exact at every density. The question is the self-averaging of a single order. The card is proposed; the author
  freezes it.

## 50. Q.1 frozen (28 Nov, 10:47)

**The author on draft 0** (10:39): the card is approved for completion and freezing, but not unchanged. Four points:
1. Define exactly what is held fixed when N changes: region, density, mass, construction, the definition of sd and
   the reference.
2. Distinguish one realisation from the ensemble: several independent realisations per N, the distribution of the
   deviations, the median and the spread of the slope.
3. Complete the decision bands: INCONCLUSIVE for the gap, not a widening of the bands after the results; a verdict
   for slopes below −0.65 (a deviation from the expected scaling, not a failure).
4. Keep the calibration (the exact mean, the non-negative null) apart from the physics verdict, and keep the
   development as a diagnostic outside the result, including the error growing with m (2.1 / 6.0 / 7.2 %).
Q.1 tests the numerical realisation of the propagator on a causal order, not quantum mechanics as a whole; that is a
correctly limited first step.

**Draft 1** (docs/quantum/Q1_card_draft.md) answers the four points:
- Fixed region and mass; the density is the only scaled quantity.
- M1 = sd of the error at the fixed pair; M2 = median over realisations of the error of one order over its own
  elements. Bootstrap intervals of both slopes.
- Complete bands: in band / fail / deviation / inconclusive. INCONCLUSIVE is recorded as UNDECIDED with a sub-label,
  because the statuses of the core are fixed.
- I1–I3 as instruments; the development in its own section.

**The author** (10:47): "Schváleno".
- The card is frozen.
- Test file phases/Q_quantum/Q_01_chain_propagator_2d.py (sha256 d58fd228…, LF md5 a596dfaa).
- The phase Q is added to core/registry.py (the phase list, the file-name pattern; LF md5 3c3376b9). This is a change
  of the core: records from now on carry a new core hash. No closed gate is affected; the classical tests are
  unchanged.
- Smoke run (--quick, 40 realisations, 7 s): I1 6e-16, I2 max |z| 1.95, I3 +0.22 SE; the bands are not read from it.
- The core checks pass (25/25).

## 51. Q.1c frozen (28 Nov, 14:41)

**Route (11:29).** The author: "C, poté A - bez něho bychom A museli opakovat". The complex (Feynman) propagator is
calibrated before the interference test.

**Construction (11:31).** The author: the Feynman propagator must be defined purely from G_ret, G_adv of the given
order, through the Pauli–Jordan operator and the SJ spectral decomposition of its positive part, without any
external intervention. The card follows this: no truncation of the spectrum, constant amplitudes, and the continuum
used only as the reference.

**Reference (11:43).** The author: the offset against the Minkowski vacuum is only reported. The reference is the
continuum SJ state of the same diamond, computed independently of any causal set.

**Draft 0** (14:0x) was approved by the author at 14:24 ("Shváleno"). The smoke run of the test file then found four
points, corrected in draft 1:
1. 8 probe groups instead of 5. No inserted point may lie inside or on the interval of a timelike pair of its group;
   otherwise E K_R is not exactly ½J₀ and I2 fails.
2. |Δ| from the SVD, because the eigendecomposition of ΔᵀΔ reached only 8e-9 against I1's 1e-10.
3. m²/ρ is at most 0.20, not 0.10.
4. I2's threshold.

**The author (14:41):** "Schvaluji všechny čtyři změny. U bodu 4 rozhoduji: I2 bude mít práh 4,0 SE."
- The card (docs/quantum/Q1c_card_draft.md draft 1) is frozen.
- Test file phases/Q_quantum/Q_01c_feynman_vacuum_2d.py (sha256 2defb344…, LF md5 e751fb95).
- Smoke run (--quick): I1 5.8e-15, I2 max |z| 3.43, I3 OK, I4 χ²/12 = 443. It decides nothing.
- The full run is on the author's PC with the workers of his choice: about 3 GB per worker at N_exp = 8000, about 10 h
  on 10 workers. An interrupted run continues from its cache (results/Q_quantum/Q.1c/cache).
