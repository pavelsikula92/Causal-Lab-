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
