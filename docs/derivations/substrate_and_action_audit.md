# Audit of A and B: the quantities of the substrate and the action

Written on the night of 15–16 Nov 2026, while the author was away, as agreed at 22:4x: "audit A + B as a working
document, read only, deciding nothing". **Status: an audit. It changes no axiom, no card and no verdict.** Its output
is the list of decisions in section 5. The algebra marked "derived here" is short and written out in the text; none of
it is tested.

**Sources read:**
- paper: 015_ontology (transport, matter, gravity), 02_foundations (substrate, flows, action, response), 03_units
  (symbols), 04_axioms;
- docs: axiom_audit.md, author_decisions.md (decisions 1, 2, 6, 7, 10, 12–15, 18, 19, 24), findings.md (Derived
  consequences of 24 and 27 Oct, F3.12a2);
- derivations: variational_audit, field_and_throughput, rule_from_action, why_gamma1, ln_pi_potential,
  gravity_equations;
- planned_tests.csv (F0.5, F3.5, F5.2d), roadmap.md (programme A–F of 15 Nov).

## 1. Map of the quantities

| Quantity | Where it is defined | What it is now | Operational? | Conflict or gap |
|---|---|---|---|---|
| event | Ax. Primacy of events; 02 §1 | element of a causal set | yes (sprinkling) | how an event **arises** is undefined (growth open, decision 7) |
| link | 03 symbols; 02 §1 | u ≺ v with nothing between | yes | not an independent object: the links of a new event are fixed by its past (section 2, A.2) |
| K₀ | book §15; 02 §1.2; decision 14 | a bound on the **degree** of an event | tested: F0.3 (E42b) asked whether the stiffness I can be read from the degree; it cannot (the degree grows with N, clustering 0); K.10: E42b not robust to the seed | the programme of 22:32 lists K₀ as "capacity of an event"; the paper calls the capacity **C**, a bound on the **flow** (section 2, A.1) |
| C | Ax. Finite capacity; 03 | the largest flow an event passes on | only as a calibration (Q₀.₉₉ of the flow under unit driving) | depends on the density and the region (axiom audit §2); supported by F1.8, not robust to the seed |
| Π | Ax. Local response; 03 ("permeability, read as the clock rate") | a state of every event, 0 < Π ≤ 1 | as a **readout** of the order: clocks (K.3, K.16b/c) | state of the event or readout of the geometry? (section 3, B.2); the programme calls it "local throughput", the paper's throughput is q = ρ_c c_eff, a different quantity |
| Φ (flow) | book ("Φ(v) = Σ w Φ(u)"); paper uses **F** | flow carried by an event, with u^a (Ax. Momentum) | in tests (F1.1, F3.6, K.15) | the symbol Φ is also the temporal potential (rule_from_action, F6.8 sectors), φ the clock potential (02 §2) and φ the kink field (F3.14a) |
| I_v | Ax. Local response (I ≥ 1); decision 10 | "resistance of the order to reconfiguration" | **no** | one candidate tested and failed: I from the degree (F0.3); two open, neither adopted: F0.5 (density of diamonds in the world tube) and the author's rerouting cost (F3.12a2) |
| Z_v | Ax. Local response | I_v/Π_v | no (I is not defined) | — |
| Δτ | Ax. Action of a flow; 02 §4 | proper time between successive events of a flow = longest chain | yes on a given order | in R it carries an event index (Δτ_v); in the book also "time as a number of updates" (remnant B, axiom audit §0) |
| R | 02 §4; Ax. Action status; book 5.8 | Σ_v Z_v Δτ_v, "minimum" (book) | no | its relation to S_m and to S_BD is open (section 3) |
| S_m | Ax. Action of a flow | −Σ m Δτ | yes (F2 on maximal chains) | two forms in the documents: −Σ m Δτ (paper) and Σ Π Δτ (gravity_equations §2, findings 27 Oct); see B.3 |
| S_BD | Ax. Dynamics of the order; 02 §BD | Benincasa–Dowker action of the order | expectation computed (F6.8a–e) | static stiffness depends on the age T of the field (F6.8a, b2, c2, e) |
| ρ_c, c_eff, q | 02 §5; 03 | coordinate density, cone speed, throughput q = ρ_c c_eff | readouts (K.16b/c) | the throughput principle is a plateau law (decision 13), derived at linear order from S_BD + S_m |

## 2. A: the dynamics of the substrate

### A.1 K₀ and C are not the same quantity in the documents

- **K₀** was the bound on the degree, deg(v) ≤ K₀ (book §15). The laboratory tested it in F0.3 (E42b, "Can the
  stiffness I be read from the degree of an event?"): the degree grows with N and the clustering is 0, so neither a
  bound K₀ nor a reading of I from the degree holds (K.10: E42b is not robust to the seed, so the numbers are weak, the
  growth of the degree with N is not), the paper explains why it cannot hold on a Lorentz-invariant substrate (02 §1.2, Bombelli–Henson–Sorkin), and
  decision 14 removed K₀ from H5 "as a superseded version of the model, not a parameter of the theory".
- **C** is the capacity in the paper: a bound on the flow through an event (Ax. Finite capacity). Its only realisation
  is a calibration from the flow itself (axiom audit §2), so it is not yet an intrinsic quantity of the order.
- **The programme of 22:32** lists K₀ as "capacity of an event". If this means C, the two names should become one. If
  it means a bound on the number of relations, it contradicts F0.3 and Ax. Lorentz-invariant substrate.

### A.2 How events and links arise

- On a causal set a link is not created independently. A new event is defined by its past (a down-set of the order);
  all its relations and links follow from that choice. "How new links arise" is therefore the same question as "which
  past a new event receives".
- The framework that answers it with conditions equal to decision 2 is classical sequential growth (Rideout–Sorkin:
  general covariance, Bell causality). Generic growth of that kind does not produce manifold-like orders; a growth that
  does is an open problem of the field (axiom audit §0).
- What the laboratory has instead:
  - prescribed orders, sprinkled into a metric the theory gives (F7.1b/c, F5.2c, F5.3);
  - one growth rule, F5.2d (2+1, decision 23): the density and the cones of new events follow the potential produced
    by matter. It is a declared proxy (cones cut in the coordinate rest frame of a static source) and the scalar branch
    (load = trace of T, decision 24).
- **Not available:** an evolution law for the order that uses only the order and counts (M3), the past only (M4), and
  is independent of the linear extension (M4).

### A.3 What limits creation and transfer

- Transfer: the capacity bounds the flow, not the number of relations (F1.8, not robust to the seed; F3.12a2: a cut of
  capacity changes the routing more than the removal of the same paths).
- Creation: nothing in the axioms limits the creation of events. In the response of the order events **are added**
  near matter (Ax. Dynamics, decision 12), ρ_c = Π⁻² in 3+1. Whether capacity limits creation (A5′: in saturation
  only) is a choice (decision 15), not derived.

### A.4 The flow and its conservation when the substrate changes

- Ax. Causal continuity is tested on a fixed sprinkling (F1.1, Gauss).
- On a growing order the balance is the reservoir model of decision 2. The laboratory noted at the time that its three
  conditions (no future, balance, independence of the order of computation) restrict the rule strongly: "a share of
  the remaining stock" violates condition 3, and a fixed share holds only until the stock runs out. No rule satisfying
  all three is written down.
- The division of the flow among descendants is a proxy (it needs the number of future descendants; in an infinite
  order there are infinitely many). Ax. Momentum solved localisation by passing the whole flow to one successor
  (02, "Contact"). Which successor, without coordinates, is defined by the four-velocity carried as a state; F3.6 shows it is
  frame independent in the tested boost.

### A.5 Two models of matter side by side

- **Flows** (F, u^a, T^ab of flows, merging F = |P|): Ax. Continuity, Ax. Momentum, K.15, F4.1, F7.4a.
- **A field on all events** (φ, retarded Green function, kink): F3.14a, F5.2a2, F5.2c/d, F3.15, F3.16b, F5.3;
  decision 18 calls the defect a test body.
- field_and_throughput §3: the flow F, the topological current j_top and T^{μν} of φ must not be identified; that F
  follows from φ is open. The capacity C is defined on F, while the load in the field tests is built from φ (decision
  19: (dφ)², decision 24: ρ + Σp).
- So "how the causal flow is defined" (programme A) currently has two answers, used in different tests.

### A.6 Π: a state of the event, or a readout of the order?

This is the central point of A and B together (developed in B.2):
- Ax. Local response gives every event a state Π_v.
- Ax. Dynamics and Part VI let the order change: the cones narrow and events are added. On such an order the clocks are
  read from the order (K.3, K.16b/c), and Π is the lapse of the geometry, not an extra state.
- ln_pi_potential §3 and §4.2: Π is a clock ratio relative to a reference (Π_∞ = 1, static observers, Killing time).
  For a non-static field the reference is undefined, so Π as a state of an event is not covariant (M4) without a
  further construction.

## 3. B: the action and the variational principle

### B.1 R is a functional of a history, not of the order

- Δτ is defined between **successive events of a flow** (Ax. Action of a flow). Δτ_v of an event that carries no flow
  is undefined. R = Σ_v Z_v Δτ_v is therefore a sum along a history in a given order, as variational_audit §2.3
  already concluded.
- R contains no mass and no label of the history. Two different bodies on the same history have the same R.

### B.2 Two pictures of Π, and what the cancellation Z Δτ = I Δt means in each

The documents use two pictures without saying which one holds.

**P1. The order is a sprinkling of the gravitational metric itself** (H1 in 3+1: −Π²dt² + Π⁻²dx²). Used by F7.1c,
F5.3, K.16c, decision 13 and Part VI. Π is read from the order (clocks); along a static history Δτ = Π Δt, so
Z Δτ = I Δt.

- **Derived here:** R = ∫ Z dτ_g has the extremals of the conformal metric Z² g. With Z = I/Π and g = H1:
  Z² g = I² (−dt² + Π⁻⁴ dx²). For I = 1 the time component is −1: a slowly moving body feels **no Newtonian
  acceleration**. The factor 1/Π in Z removes exactly the clock slowing that produces the fall.
- So in P1 "what carries the gravitational influence of Π" has a definite answer: the order (through Δτ, i.e. S_m),
  not R. R moves bodies along the H1 geodesics only if Z is constant along the history (B.4).

**P2. The order is a sprinkling of a metric g_u with narrowed cones but flat clocks** (g_u = −dt² + Π⁻⁴dx²), and Π is
an extra state of the events. This is the reading of findings, Derived consequences 27 Oct, and of gravity_equations
§2 ("A6, S = Σ Π Δτ, gives the effective metric Π² g_u; H1 arises when g_u has the cone speed Π²").

- Along a static history Δτ = Δt, so Z Δτ = (I/Π) Δt: **no cancellation**.
- **Derived here:** R = ∫ (I/Π) dτ_u has the extremals of (I/Π)² g_u = I² (−Π⁻² dt² + Π⁻⁶ dx²). For I = 1 the clocks
  run **faster** near matter (Π⁻² > 1): the static force points **away** from the source.
- **Derived here:** the coordinate density of a sprinkling of g_u in 3+1 is √−g_u = Π⁻⁶, not the Π⁻² that the
  throughput principle requires of the order (ρ_c c_eff = 1 needs ρ_c = Π⁻² when c_eff = Π²). P2 therefore meets the
  throughput principle only in the effective metric of motion, not in the count of events.
- In P2 the BD action of the order sees g_u only. Π is invisible to it and needs its own equation: this is the route
  of the field equation B_g ln Π = κρ_s (decision 1), with B_g on the order of geometry g.

**Summary of B.2.** The cancellation Z Δτ = I Δt is not a choice of the definition of Δτ once Δτ is the proper time
counted by the order: it holds in P1 and fails in P2. In neither picture does R with I = 1 give the Newtonian fall.

### B.3 Two forms of the action of a flow in the documents

| Form | Where | Consistent with |
|---|---|---|
| S_m = −Σ m Δτ | paper (Ax. Action of a flow, 02 §4) | P1: on an order of geometry H1, geodesics of H1 |
| S = −ħω_B Σ Π Δτ | findings 27 Oct; gravity_equations §2 and §6 ("input: clocks Π, action S = Σ Π Δτ (A6)") | P2: Π² g_u = H1 |

**Derived here:** using the second form on the order of P1 counts Π twice: Π² · H1 = −Π⁴ dt² + dx², i.e. twice the
Newtonian acceleration. The two forms are equivalent only if the picture is fixed with them. gravity_equations §6
should say which one it uses.

### B.4 R against S_m

- The paper says R coincides with S_m "only if Z_v is proportional to the mass" (02 §4, Ax. Action status). More
  precisely: the extremals coincide only if Z is **the same at every event of the history** (Z_v = λ m). A Z that
  varies along the history changes the metric of motion from g to Z² g (B.2).
- **Derived here:** in P1, Z_v = λ m means I_v = λ m Π_v. With Ax. Local response (I ≥ 1, Π ≤ 1) this fails wherever
  Π_v < 1/(λ m), in particular in matter, where the ontology has Π → 0 (book 19.1). If Z ∝ m is wanted, either I ≥ 1
  or "Π → 0 in matter" has to go.
- **The extremum.** Book 5.8 states motion as the **minimum** of R. For timelike histories in a Lorentzian order proper
  time has no minimum: histories close to null have arbitrarily small Δτ, and on a causal set the shortest chain
  between two related events is the relation itself. The laboratory and the ontology use the **maximum** (F2: action on maximal chains; ontology
  "Falling": "the history along which the most events happen"). A minimum principle needs a non-Lorentzian time
  (time as a count of updates, remnant B), or Z must diverge near the cone. As written, "minimum of R" and "maximum of
  proper time" are different principles.

### B.5 R against the expected BD action

- S_BD is a functional of the order: all events, no Π, no I, no Z. S_m and R are functionals of a history in a given
  order. They are different kinds of object; R cannot be the action of the order.
- The derivation that exists (rule_from_action §2–5, Part VI): S_total = S_BD + S_m. Its stationary point at linear
  order is linearised general relativity coupled to T_ab (γ = 1 in 3+1). F6.8a/b2/c2 compute the coefficients of the
  expected S_BD (6, 8, 2 against Einstein–Hilbert). R does not occur in this chain.
- **Where R could fit** (options, none derived):
  1. **R is S_m.** Then Z is a constant of the history (λ m), I = λ m Π, and the conflict of B.4 has to be resolved.
  2. **R is a functional of the resistance of the substrate**, not an action: the cost of rerouting a flow (the
     author's hypothesis of 31 Oct, I_v = C(H → H′ | Φ, Q, K₀); F3.12a2 measured that capacity dominates). It then
     belongs to A (saturation, A5′), not to the equations of motion.
  3. **I_v is the stiffness of the order**, i.e. the second variation of S_BD. This is what F6.8 measures: a
     non-local quadratic form with a local part K₀ = 1.00005 (F6.8a; this K₀ is the local stiffness ratio of F6.8, not
     the degree bound) and age terms. It is not a per-event scalar and not ≥ 1 in general (it changes sign below
     ~13 ℓ_k, and with the age in the optical sector, F6.8e development).
  4. **R is a stability or Lyapunov functional** (variational_audit §1). Then it needs a dynamics to be monotone along,
     and that is A.

### B.6 What equations of motion exist now, and their conditions

- **Matter:** geodesics of the geometry of the order, from S_m (P1). Tested on fixed orders (F2.2, F2.4, F2.5); the
  fall of a localised body from the order: 2+1 passed (F5.1d), 3+1 failed (F5.1e); in 1+1 with the field dynamics
  F5.2a2 and F5.3 (C1, development of the well).
- **The order:** the linear response from S_BD + S_m, derived in the continuum limit, with the coefficients of S_BD
  computed. Conditions of validity found by the laboratory:
  - linear order (decision 13: throughput is a plateau law);
  - static fields;
  - the age of the field: the expected static stiffness depends on how long the field has existed (conformal and mixed
    sectors linearly, the optical sector faster: F6.8c2, F6.8e development ∝ (kT)⁵, F6.8d running). The weak-field
    limit is clean only where these terms are small (k³ ℓ_k² T ≪ 1 in the conformal sector). Which age-term form holds
    in the optical sector is not yet decided (F6.8e3, F6.8d).
- **No classical evolution equation of the order** is derived from S_BD. On a causal set S_BD weights a sum over
  orders. A classical "stationary order" is a continuum statement. This is why C.1 waits for B: the propagator must be
  the second variation of a fixed action, including or explicitly excluding the age terms.

## 4. Stand-ins used in tests now

| Stand-in | Where | Stands in for |
|---|---|---|
| order sprinkled into a prescribed metric | F7.1b/c, F5.2c, F5.3, C1b development | the dynamics of the order (A) |
| cones cut in the coordinate rest frame of a static source | F5.2d | an intrinsic frame (B3b) |
| load = trace of T | F5.2b, F5.2d | the source ρ + Σp (decision 24) |
| C = Q₀.₉₉ of the flow under unit driving | F1.8, F3 | an intrinsic capacity |
| uniform division of the flow among descendants | F3 series | a local transport rule (proxy) |
| kink with a chosen potential | F3.14a ff. | matter (decision 18: a test body, by design) |
| H1 from the throughput principle | gravity_equations, F7.1c | the response of the order beyond linear order (decision 13: GR completion) |

These are declared in their cards. The table only collects them, because programme A asks that none of them be read
as the dynamics of the substrate.

## 5. Decisions for the author

Each item lists the options found in the documents and what each one implies. The audit recommends none.

**A1. K₀.** F0.3 tested K₀ as the degree (and I read from it) and it did not hold. Is K₀ in the programme the capacity
C of the paper (a bound on the flow)?
- (a) Yes: one name for one quantity; decision 14 and 02 §1.2 keep their content (the degree bound stays refuted).
- (b) No, K₀ is a bound on relations: it contradicts F0.3 and Ax. Lorentz-invariant substrate, and the choice of
  substrate A would have to be reopened.

**A2. Π.** Is Π a state of an event (P2) or the clock rate read from the order (P1)? (B.2)
- (a) P1: Π is a readout. The order carries gravity; the action is S_BD + S_m; Ax. Local response becomes a statement
  about the order. The throughput principle holds in the count of events (ρ_c = Π⁻²).
- (b) P2: Π is a field on an order with flat clocks. It needs its own field equation (B_g ln Π = κρ_s, decision 1) and
  a reference for "clock rate" in non-static fields; the order then has density Π⁻⁶, not Π⁻².
- Also: should the programme's "Π = local throughput" become "permeability" (03), since the throughput of the paper is
  q = ρ_c c_eff, conserved, not Π?

**A3. The flow.** Which is the matter whose flow is conserved: flows F with u^a, the field φ, or a relation between
them (F derived from φ, field_and_throughput §3)? Capacity is defined on F; the tests of 1+1–2+1 gravity use φ.

**A4. Notation.** Φ is the flow (book, programme), the temporal potential (rule_from_action, F6.8), φ the clock
potential (02) and φ the kink field (F3.14a). One symbol per quantity?

**A5. Growth.** Is the evolution of the order to be a sequential growth (Rideout–Sorkin conditions = decision 2), a sum
over orders weighted by S_BD + S_m with a classical limit, or the laboratory's F5.2d-type rule made intrinsic? The
first two are open problems of the field; the third is a proxy now.

**A6. Conservation on a growing order.** Which reservoir rule satisfies the three conditions of decision 2 at once, or
which condition is relaxed?

**B1. Δτ.** Is Δτ the proper time counted by the order (longest chain, Ax. Action), or a count of updates (book §20.1,
remnant B)? With the first, Z Δτ = I Δt holds in P1 and the minimum principle of book 5.8 is ill-posed (B.4).

**B2. The extremum.** Motion as the maximum of proper time (F2, ontology "Falling") or as the minimum of R (book 5.8)?
They are different principles (B.4).

**B3. What R is.** One of: (1) the action of a flow (then Z = λ m along a history, and I ≥ 1 or Π → 0 in matter has to
go); (2) a resistance functional of rerouting (A, saturation); (3) I as the stiffness of the order, i.e. the second
variation of S_BD; (4) a stability or Lyapunov functional. Or R is retired as superseded, with S_m and S_BD kept.

**B4. I_v.** The degree is excluded as a reading of I (F0.3). If I_v stays a local quantity: F0.5 (density of diamonds in the world tube, "independent of Π") or the
rerouting cost of F3.12a2? Note for F0.5: in P1 every count in the world tube depends on the geometry, and the geometry
is set by Π (ρ_c = Π⁻²); "independent of Π" may not be possible for a count of the order.

**B5. The form of the action of a flow** in gravity_equations §2 and §6: −Σ m Δτ (P1) or Σ Π Δτ (P2)? Mixing them
doubles the Newtonian acceleration (B.3).

**B6. The action from which C.1 is derived.** The expected S_BD with its age terms, or its local part only, with the
age terms as a stated condition of validity? This fixes what "the asymptotic static propagator" is.

## 6. What this audit does not do

- It derives no new equation of motion and runs no test.
- The "derived here" lines are one-line conformal-factor algebra (extremals of ∫ Ω dτ_g are the geodesics of Ω² g).
  They are statements about the documents' own definitions, not about nature.
- Literature statements (Bombelli–Henson–Sorkin, Rideout–Sorkin) are quoted as the axiom audit quotes them and have
  not been re-checked tonight.

## 7. The author's working positions on A2 and B6, and what follows from them (16 Nov)

**The author (16 Nov, 09:48), working positions, not frozen:**
- **A2:** Π is a readout of the order, Π = 𝒫[H], not an independent state. A state variant would be accepted only if a
  test showed that the same local order can carry different Π according to its own history.
- **B6:** the full expected BD action is the basic object, and its local part is a controlled approximation. The full
  action needs a precise definition of the static configuration, of the boundary of the past, and of how its age
  changes under a transformation.
- **C.1:** the response of the order, δH from S⁽²⁾_BD, turned into δΠ by 𝒫. It is not a propagator chosen for Π
  first.
- **Request:** write out the consequences of A2-B for the definition of the action, then of B6-full for its second
  variation, and see whether the combination closes.

This section does that. It decides nothing.

### 7.1 What follows from A2-B (Π = 𝒫[H])

1. **One kind of degree of freedom.** The action is varied over H (and over matter) only. Π is not varied. Ax. Local
   response turns from "every event has a state Π_v" into "Π is read from the order near v". B.2 then holds in
   picture P1: the order is laid out in the gravitational metric and the clocks are read from it.

2. **Π is not a local readout, only a ratio of two clocks.**
   - A region of uniform Π is flat spacetime in rescaled coordinates (in H1 both t and x are rescaled by constants). Its
     order is a sprinkling of flat spacetime, identical in distribution to the order at Π = 1: no count inside the
     region distinguishes Π = 0.9 from Π = 1.
   - What the order fixes is the ratio Π(x₁)/Π(x₂) between two histories joined by signals, the redshift. Π itself is
     fixed only with a normalisation (Π_∞ = 1 for a static field, ln_pi_potential §4.2).
   - In a non-static field the lapse depends on the slicing, so a "Π at an event" needs a chosen family of observers.
   - **Consequence for the proposed discriminating test.** "The same local order with different Π" holds trivially for
     a uniform difference of Π, so that test would wrongly favour the state variant. The test has to be stated on
     ratios between two places (or on the gradient of ln Π), where A2-B makes a definite claim.

3. **Π has no field equation of its own.**
   - In general relativity, for a static metric −N²dt² + h_ij dx^i dx^j the 00 equation is the Whittaker equation
     D²N = 4πG N (ρ + Σp), with D² the Laplacian of h.
   - The source ρ + Σp is the one of decision 24. The scalar law B_g ln Π = κ ρ_s (decision 1) becomes the weak-field
     limit of this component, not an extra law.
   - Π fixes only one component. The cones and γ come from the spatial equations. Under A2-B, C.1 must therefore give
     both potentials (the lapse and the spatial part), and γ is their ratio. The readout 𝒫 of Π alone does not close
     the response.

4. **R.** In P1, R with I = 1 gives no Newtonian fall (B.2). R is the action of motion only if Z is constant along a
   history (B.4). A2-B therefore points B3 towards options 2–4 (a resistance functional, the stiffness of the order, a
   stability functional) or towards retiring R. This is a consequence, not a decision.

5. **"Π → 0 in matter"** (book 19.1). With a GR-like response, the lapse inside a static star stays finite; it
   vanishes only at a horizon. Under A2-B the ontology's statement has to become "smaller", or it applies only at a
   horizon.

### 7.2 What follows from B6-full for the second variation

1. **S⁽²⁾ depends on the history of the source.** The expected action of a static field counts intervals reaching
   back over the age T of the field, so the static quadratic form is K̂(k, T), not K̂(k). The sectors so far:
   - conformal: linear in kT (c_T = 0.1035; the paper's window, Sec. "window");
   - mixed: (kT)³ (development; F6.8d undecided, decision 33);
   - optical: (kT)⁵ (F6.8e3 PASSED, c5 = 3.05e-7).

2. **The window with the optical sector, if the second order is taken literally.**
   - Relative size δ_opt = c5 (kℓ_k)² (kT)⁵. With k = 1/L it is of order one at L*_opt = (c5 ℓ_k² T⁵)^(1/7). The
     conformal form L* = (c_T ℓ_k² T)^(1/3) is Eq. (Lstar) of the paper.
   - **Derived here, arithmetic only.** The same four readings of T as the paper's table:

     | reading of T | L* conformal, ℓ_P | L*_opt, ℓ_P | L* conformal, 10⁻¹⁹ m | L*_opt, 10⁻¹⁹ m |
     |---|---|---|---|---|
     | age of the universe | 1.5e-15 m | 6.0e7 m | 51 μm | 2.0e12 m |
     | age of Earth's field | 1.1e-15 m | 2.7e7 m | 35 μm | 8.9e11 m |
     | one year | 6.4e-19 m | 3.5 m | 21 nm | 1.1e5 m |
     | one day | 8.9e-20 m | 5.1 cm | 3 nm | 1.7 km |

   - The conformal columns reproduce the paper's table, which is a check of the arithmetic.
   - Taken literally, the optical age term would change gravity at order one at lengths where gravity is measured and
     agrees with general relativity: laboratory scales, satellites, the solar system. This holds for every reading of
     T and for any ℓ_k ≥ ℓ_P.

3. **But the second order does not apply there.**
   - The expansion in the amplitude b of the field fails near b²(kT)⁴ ~ 1 (card F6.8e, Limits).
   - At the L*_opt of the ℓ_P column above, with b of the order of Earth's surface potential (7e-10), b²(kT)⁴ is 3e40 to
     1e55. The second order would need b < (L/T)², e.g. b < 4e-37 for a mode of Earth's field at 2.7e7 m. No real
     static field meets this.
   - So the T⁵ law of F6.8e3 holds in its computed range (young or extremely weak fields). It says nothing yet about the
     response of a real, old static field, and neither the size nor the sign of that response is known.

4. **What this means for the age.** The age that enters is the time over which the source has existed as a static
   source: the proper time along its worldline since it formed. That singles out the rest frame of the source, so a
   static field with an age is not a boost-invariant configuration. Whether the departure from B/8 = A/6 comes from
   this (physics of a source with a beginning) or from the sharp cut is still open. Point 3 does not decide it.

### 7.3 Does A2-B + B6-full close?

- **A2-B on its own closes:** one kind of degree of freedom (H); Π read as a ratio of clocks; the Π equation is the
  lapse component of the response; the source is the one of decision 24. Its open items are the ones the author named
  (uniqueness and independence of 𝒫 from the instrument), together with point 7.1.2 (Π is defined only as a ratio, or
  with a normalisation).
- **B6-full does not close as "second variation → propagator" for physical fields.** S⁽²⁾ of the full expected action
  is valid only for b²(kT)⁴ ≪ 1. Real static fields are outside this range by many orders of magnitude. There, the
  literal second-order result would contradict measured gravity for every reading of T.
- **B6-local with the age terms as a stated condition is not a way out either.** The condition is violated by the same
  factors (δ_opt ≫ 1), so the local propagator would be used where the full action says it does not hold.

**Ways forward found in the documents (options; none chosen):**
1. **A finite-b computation of the optical sector:** the response resummed in b(kT)², for an old static field. This is
   a new computation with the stable instrument of F6.8e, and the most direct one.
2. **The reading of T.** If the relevant past of a dynamic order is not a sprinkling of the static metric back to a
   sharp onset (the past is itself produced by the response, A), the long thin intervals that carry the (kT)⁵ growth
   may not exist in the order. This needs the dynamics of A, which does not exist yet.
3. **The expected S_BD is not the action of the order** for old fields, and another action has to be found. This would
   reopen Part VI.
4. **C.1 only for young fields** (kT small), declared as such. This is well-defined but does not give the far field of
   a real body (F8.1, P1).

**For C.1 now:** under A2-B + B6-full the static propagator of a real body is not available from the second variation.
C.1 can be derived now only in form 4, or for the local part with the explicit statement that its condition fails for
real fields. Option 1 is the computation that would tell whether the full action has a usable static limit at all.

## 8. The author's working conclusions after section 7 (16 Nov, 09:59; not frozen)

- **A2-B accepted as a working conclusion.** The only dynamical object is H; Π = 𝒫[H] is not varied and has no
  independent equation. It is a relative temporal response of the order: a constant value is not locally measurable,
  and the physical content is carried by ∇ ln Π, or by the ratio of clocks at two places.
- **For C.1:** no propagator of a separate scalar field Π. The response of H has to give at least two geometric
  outputs, the temporal one (read as the relative rate of clocks) and the spatial one (the cones, the spatial metric).
  These are two outputs of the same substrate, not automatically two new degrees of freedom. The Whittaker equation is
  a static reference condition, not C.1 itself; the temporal part does not settle the spatial response or γ.
- **B6-full: the second variation is not enough for an old field.** The deciding parameter is ε = b(kT)². Once
  ε² ≫ 1 the second variation describes only the neighbourhood of b = 0. From the T⁵ term neither the size nor the
  sign of the correction for an old static field can be inferred, and no condition for the local approximation can be
  derived from it.
- **New development branch F6.8f-dev: the optical response at finite amplitude.**
  - Aim: whether the expected BD action gives a usable static response at finite b and large T. Not a test of gravity.
  - Quantity: E_opt(k, T, b) = ∂E[S_BD]/∂b and K_eff(k, T, b) = [E_opt(b) − E_opt(0)]/b, with the projection on the
    optical sector exactly as in F6.8e3. For b → 0 it must reproduce F6.8e3 (a check, not a result).
  - Grid (development, not criteria): |ε| = 0, 0.1, 0.25, 0.5, 1, 2, 4, each at several combinations of b, k and T, to
    see whether the physics is governed by b(kT)² or the scaling holds only in the perturbative limit.
  - Checks: the b → 0 limit; zero age; k → 0; both signs of b; the conformal and mixed sectors not broken; convergence
    with the stable instrument of F6.8e.
  - Output: the curves K_eff(k, T, b), their sign, convergence, saturation and T-dependence.
  - Readings: (A) the response converges at large T: on to the full static C.1, B6-full stays a candidate; (B) it grows
    without bound or changes sign: the expected BD action gives no usable old configuration, options 2 or 3 of §7.3
    open, a local propagator does not solve it; (C) a static response exists only for young fields: option 4 is
    legitimate as a limited result, and C.1 must not be presented as the propagator of a real old body.
  - None of this changes F6.8e3 or F6.8d.
- **B3:** R is not used as a basis of gravitational dynamics. It can stay only if it is shown what it measures and with
  which discriminating experiment; otherwise retiring it is more reasonable.
- **Order of work:** F6.8f-dev; the static limit; C.1; B3.

**Notes by Claude on the design (16 Nov):**
1. **Both signs of b, phase-averaged.** For ν = b cos kz the change b → −b is the phase shift z → z + π/k. The
   phase-averaged action is therefore exactly even in b, and check 4 is met identically. It tests something only per
   phase, before the average (as a code check).
2. **The core of the F6.8e instrument is perturbative.** The shifted coordinates z = w + bδ use the first-order
   deviation of the ray. At finite b this core must be replaced; the long double and the tip clustering can stay.
   - The medium is stratified, n = e^{b cos kz}. Snell's invariant q = n sin θ (θ measured from the z-axis) is
     conserved along a ray. The travel time and the transverse displacement between two points are therefore 1D
     quadratures in z, with turning points where n(z) = q.
   - The first arrival is the minimum over the ray branches, so caustics are included exactly, with no expansion in
     b. Its derivative in b is ∫ cos(kz) n dl along the arriving ray (Fermat's principle), so E_opt can be computed
     directly as well as by differences.
3. **Where the nonlinearity starts (derived here, paraxial, one line).** A ray perpendicular to k obeys
   z'' = ∂_z ln n = −bk sin kz. About the layer of largest index (kz₀ = 0 for b > 0) neighbouring rays oscillate with
   ω = k√b (a guided ray), so rays leaving a point refocus at L_f = π/(k√b): the first caustic lies at b(kL)² = π² ≈ 10.
   Rays along k do not focus in this way.
   - The τ-independent part ∝ b²(kL)⁴ = ε² is the onset of this focusing.
   - The grid should therefore reach beyond the first caustic, to ε ≈ 20–40, not only to 4.
4. **Precision.** At second order the BD moments cancel the τ⁰, τ², τ⁴ parts of V(τ) and leave ~10⁻⁹ of them. At finite
   b the same structural cancellation is expected, so the exact travel time must be accurate to ~10⁻¹³ relative. This
   is the main risk of the branch.
5. **Status of F6.8d:** INSTRUMENT_FAILED, with the exponents undecided (decision 33), not UNDECIDED.
