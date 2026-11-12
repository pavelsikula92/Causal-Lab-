# Direction B: the field, the growth rule and gravity in 2+1 (preparation)

Written on 11 Nov 2026 at the author's request, while F5.2b is being replicated. **Status: a plan with decision points.
Nothing here is tested.** Labels: *carries over* (works in 1+1 and transfers unchanged), *new* (needed in 2+1), *open*
(a decision of the author or an unsolved problem).

## 0. What 1+1 gave and what it cannot give

In 1+1 one field and one causal growth rule gave all three:
- the inertial mass (F3.14b, replicated);
- the passive mass (F5.2a3, replicated);
- the active mass (F5.2b, seed 0).

1+1 cannot give three things:
- **gravity at a distance:** the response of a static source grows linearly, with no 1/r;
- **a change of the causal relation:** every 1+1 metric is conformally flat, so the counts alone carry the geometry;
- **a check of light bending:** there is no transverse direction.

## 1. What 2+1 is for, and an honest warning

- **The goal:** a bound defect that changes the order around it so that
  - another defect falls towards it with a force that falls off with distance;
  - **the cones narrow**, so the relation itself changes, not only the density.
- **The warning:** in 2+1 general relativity has **no** Newtonian attraction (a point mass makes a conical deficit, and
  there are no local forces). The scalar response of the theory (R2, decision 15) does give one: in 2+1 the static
  solution of ∇² ln Π = κq is ln Π ∝ ln r, a force 1/r. A 2+1 test therefore checks **the response of the theory**,
  not general relativity. The comparison with general relativity belongs to 3+1, where the theory adopts it as its
  completion (decision 13).

## 2. What carries over unchanged

| piece | 1+1 | 2+1 |
|---|---|---|
| the field rule | φ(x) = φ₀(x) − Σ_{y≺x} G(x, y) V′(φ(y))/ρ | the same form, with the 2+1 Green function (B-a) |
| the growth rule by realisation | probability Π² from the realised past | density: the same principle; the exponent follows from the throughput (B-c) |
| marks fixed at the sprinkling, causality, order independence | yes | yes |
| the load as the trace of T | 2V (exact, the kinetic term is conformal) | −½(∂φ)² + 3V: the kinetic part no longer drops out (B-b) |

## 3. What is new: four problems, each with options

### B-a. The Green function of the order in 2+1 (the instrument)

In the continuum the retarded Green function of ∂_t² − ∇² is G = 1/(2πτ) inside the cone, where τ is the proper time.
It is no longer a constant, so it is no longer ½ × the relation.

**Intrinsic options:**
1. **τ from the interval count n_xy** (the Alexandrov volume in 2+1 is πτ³/12). Use the weight E[1/τ | n], which is
   finite also for links (n = 0): w(n) ∝ ρ^{1/3} Γ(n + 2/3)/Γ(n + 1).
2. **τ from the longest chain L(x, y)** (Brightwell–Gregory). It is coarser for short intervals, but cheaper.

**Problems:**
- **The light cone:** the number of links to the past is unbounded. With an initial slice the sum is finite, and the
  error near the cone shrinks like ρ^{−1/3} (an estimate, to be measured). This is criterion (d) of the audit again.
- **The cost:** the field needs the weight of every pair in the past. With interval counts that is O(N²) counts of
  O(N). For a kink resolved in 2+1 (ρ ~ 200 per unit³, N ~ 3·10⁵) this is out of reach. With coordinate τ it is
  feasible (~5·10¹⁰ pair evaluations), but that is a **declared proxy** for the instrument stage only.
- **Noise:** the shot-noise formula (docs/derivations/swerves_of_a_defect.md) predicts the noise in advance. It has to
  be recomputed for G ∝ 1/τ.

**Recommendation for stage B0:**
- develop the instrument with coordinate τ (declared, exactly as the readouts use coordinates);
- in parallel, validate the intrinsic weight (option 1 or 2) on small boxes against it;
- the frozen tests use only the intrinsic version.

**Decision for the author:** option 1 or 2, and whether the coordinate stage is allowed as an instrument.

### B-b. The load in 2+1

The trace of the energy-momentum tensor of a scalar field in 2+1 is −½(∂φ)² + 3V (signature −++). The kinetic term no
longer drops out. The development rows of 10 Nov showed that the BD estimators of the kinetic term are unusable per
event.

**Options:**
1. **A discrete gradient** from the Green function itself (the same weights): local fits of φ on the past of an event.
   To be developed.
2. **For a static defect (dust), its mass from the virial relation**, as in 1+1. This holds only on solitons, not for
   radiation.
3. **The load from the potential alone** (V), declared as an approximation for bound matter at rest.

**Recommendation:** first measure option 1 on the defect (as in "Load candidates", 10 Nov), and only then decide.

### B-c. The change of the relation (the core of B)

**Why density alone is not enough in 2+1:**
- The order of the theory has the metric −Π²dt² + Π^{−2}dx² (H1, throughput, decision 12). Its cones narrow
  (c = Π²) and its density per coordinate volume grows (√g = Π^{−1}).
- A density change alone, at Minkowski cones, gives a conformally flat metric Π²(−dt² + dx²).
- Light is not deflected in a conformally flat metric (null geodesics are conformally invariant). That is the wrong
  value of γ: it is excluded by observation in 3+1, and it is not the theory.

**So the relation has to change.**

**What the laboratory already has** (F7.1b, 2+1, replicated as UNDECIDED, γ ≈ 1.85 ± 0.1 against the prediction 2):
- the events have coordinate density e^σ;
- a link is kept only if its speed ≤ e^{−σ} at the parent, **in the rest frame of the source**;
- the relation is reachability through the kept links.

This realises "narrower cones and more events" (S·c = const) and gives a metric order.

**The problem:** the speed is read in coordinates, and the audit classifies a coordinate rule for narrowed cones as a
proxy (F3.3h, failed the frame test F0.6).

**An intrinsic candidate (new, open):** "the rest frame of the source" is the rest frame **of the matter that loads
the event**. The field φ itself defines it: its four-velocity u from T_ab, or the local direction of the flow of its
load.
- A link is kept if its rapidity relative to that u is below the limit set by Π (computed from the realised past).
- Where there is no matter (vacuum) there is no preferred frame, and no link is cut.
- This is consistent with the theory (the throughput is a property of loaded events) and does not need coordinates.
- It does need u read from the order. In 1+1 that was a readout; here it would become part of the dynamics. **To be
  constructed.**

**Decision for the author:** whether the relation is to change through kept links (the F7.1b route) with the frame of
the matter that loads the event, or by another mechanism.

### B-d. The defect in 2+1

**Derrick's theorem in 2 spatial dimensions:** a scalar field with only a gradient term and a potential has no stable
localised static soliton. The gradient energy is scale invariant, the potential energy pushes the soliton to shrink.

**Options:**
1. **A wall** (the kink of 1+1 extended along y, periodic in y). Stable, with a known profile. It carries the whole
   machinery of 2+1 (Green function, relation change, load). It has no 1/r: the field of a wall is that of a plane.
2. **A baby Skyrmion** (an O(3) field with a Skyrme term and a potential). Localised and stable, but a new field type
   and a higher term.
3. **A vortex of a complex field.** Without a gauge field its energy diverges logarithmically.

**Recommendation:** walls first (stages B0–B2), then a decision on a localised defect for 1/r (B3).

**Decision for the author:** whether walls are acceptable as the first 2+1 object, and later which localised defect.

## 4. Stages

| stage | what | criterion of success | depends on |
|---|---|---|---|
| **B0** | the 2+1 Green function of the order as a field evolution (as H5-0): stability, dispersion ω(k), convergence, order independence; cost; intrinsic weight against the coordinate one | the instrument conditions of H5-0 transferred | B-a  **closed 11 Nov (decision 20): variant B, the chain weight of decision 17** |
| **B1** | a wall in 2+1 on a flat order: inertia (as F3.14b) | M_in = M_rest, contraction | B0 |
| **B2** | a wall in a prescribed curved order of the theory (density e^σ and cones e^{−σ}, the order of F7.1b): fall, and the light bending γ read by a null probe | fall as the continuum; γ = 2 in 2+1 (throughput) | B1, B-c |
| **B3** | the growth rule with a changing relation: a wall changes density and cones from its realised past; a second wall falls | active mass = energy, universality, γ | B2, B-b, B-c |
| **B4** | a localised defect and the 1/r law of the response R2 | force ∝ M_A / r | B3, B-d |

**Costs:**
- B0–B2 are hours on 6 processes with coordinate τ.
- The intrinsic weight for full-size runs needs an algorithm for the weights of all pairs better than O(N³): a
  computational problem of its own.

## 5. Decisions for the author (summary)

1. **B-a:** the τ estimator for the Green function (interval count or longest chain), and the coordinate stage as an
   instrument.
2. **B-b:** develop the discrete gradient for the load (recommended), or use the virial or potential approximation.
3. **B-c:** the change of the relation through kept links in the frame of the matter that loads the event (to be
   constructed), or another mechanism.
4. **B-d:** walls as the first 2+1 object.
5. The 2+1 test checks the response of the theory (R2), not general relativity, which has no Newtonian force in 2+1.
