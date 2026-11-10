# Card H5: inertia of a locked causal defect

Draft, revision 2, for the author's approval.

- Revision 1: 10 Nov 2026.
- Revision 2: 10 Nov 2026, after the author's decisions D1–D4, his three technical conditions and the order of the
  experiment he set.

A new hypothesis after F5.1g (FAILED: M5) and F3.13c (UNDECIDED). Their verdicts and criteria do not change (R132).
No test code is written before the author approves this revision. The probes of section 3 are development rows: they
decide nothing.

## 0. What F3.13c did (the author's question, answered from the code)

**No event is created, and the source term is exactly zero after launch.** The events are the sprinkling. A source
only adds the flow s_v = 1 for 0 ≤ t_s < 2 and r_s < 5. After t_s = 2, s_v = 0 exactly.

**The defect moved as dust.** Each carrying event sends its whole flow to its child of least rapidity relative to its
own u (the F3.6 rule, which uses coordinates). Each world line therefore kept its own velocity. The orientation n rode
on top: the alignment changes n and never u or the route. Binding and motion were decoupled. With the alignment
switched off, the defect would have moved in the same way.

**F3.13c did not test internal binding.**

Three consequences, recorded in revision 1 and accepted by the author:
- without flow, the orientation of F3.13 has no carrier;
- relaxational alignment gives no inertia, because defects of an overdamped dynamics move with v ∝ force and stop
  when the force stops;
- in an order every pattern is rebuilt ahead and left behind, so motion is told from deformation by **shape**, not by
  the identity of events.

## 1. Decisions of the author (10 Nov)

| | decision | reason (the author's) |
|---|---|---|
| D1 | the field lives on **all events** | a field of the order must not need a flow as its carrier, otherwise H5 only moves the problem to the definition of the flow |
| D2 | **conservative, second-order** dynamics | relaxation cannot supply inertia; using it would build a limitation into the experiment |
| D3 | **1+1 first** | the sine-Gordon kink gives an analytically known reference for stability, velocity, contraction and response to a force |
| D4 | **finite capacity as a tested interpretation** of the potential | the potential may represent a local capacity bound, but the identification is not derived |

## 2. The author's three technical conditions

### C1. A wave operator is not yet a well-defined evolution

Before any kink, H5-0 must give all of the following:
- the discrete operator in the given dimension;
- the explicit evolution equation and the initial data;
- the stability of the linear modes of the vacuum;
- the dispersion relation ω(k) against the continuum;
- convergence with the density of the sprinkling;
- independence of the solution from the order of computation and from the division into layers.

The last point is a check of the implementation, not a proof of causality: causality is in the construction.

### C2. Mass is measured two independent ways

The two masses are M_in = F / a and M_rest = E_kink. The test reads M_in / M_rest → 1. Neither may be computed from
the other through an assumed equation of motion. The definition of F must be written out (section 4).

### C3. The potential is not the capacity until it is derived

The base model has a potential with explicit parameters. H5-D, the capacity interpretation, comes last and must derive
the potential from the capacity bound without tuning to results. If that fails, the inertia result can still stand; the
interpretation cannot.

## 3. Development probes of the evolution (10 Nov; decide nothing; 1+1, density ρ, one to three seeds)

Scripts are in the scratchpad. They become `diag/h5_probe_1d.py` when H5-0 is coded.

### P-a. The BD operator in retarded form is unstable

The smeared 2D BD operator of the lab (core.operators.bd) was solved for φ(x) at each event:
φ(x) = −[Σ_{y≺x} K_ε(n_xy) φ(y)] · 2ερ / (2ερ + m²).

**On a constant field the normalisation is right.** The sum Σ K over the past has mean −1.0 ± 0.1. Per event it
fluctuates with sd 1.1 (ρ = 100, ε = 0.2), 0.9 (ρ = 400, ε = 0.1) and 0.5 (ρ = 400, ε = 0.05), while Σ|K| ≈ 7–11.

**The evolution of a constant field blows up.** At ρ = 100 and ε = 0.2 the rms error was 1.7 one unit after the
start, 18 after two units and 97 after 3.5 units. Each event's value is a sum with large, cancelling weights over noisy
past values, so the Poisson noise is amplified at every step.

**Reading:** the differential form of the operator is fine in expectation, but as an evolution it amplifies noise.
This is exactly the risk C1 names.

### P-b. The order's retarded Green function is stable

The continuum retarded Green function of ∂_t² − ∂_x² in 1+1 is ½ inside the past light cone. On a causal set it is ½
times the causal relation itself (Johnston 2008; check before citing). The field equation
φ_tt − φ_xx + V′(φ) = J then becomes an explicit, retarded rule that uses only the order and counts:

  φ(x) = φ_0(x) − (1 / 2ρ) Σ_{y ≺ x, y after the initial slice} [ V′(φ(y)) − J(y) ]

φ_0 is the free wave of the initial data (d'Alembert's formula, section 5). Each event is computed once, from its
past. No weight is divided by a small coefficient.

With V = 1 − cos φ, ρ = 100, a box |x| ≤ 20, t ∈ [0, 8] (or 10), and readout only inside the domain of dependence of
the initial slice:

| arm | result |
|---|---|
| kink at rest | Q = 1.00 in every slab; rms error against the exact kink 0.007 → 0.06 over 8 units |
| kink with v = 0.5 | centre 0.26, 0.75, 1.26, 1.78, 2.31, 2.82, 3.27, 3.71 against the exact 0.25 … 3.75; Q = 1.00 |
| kink with v = 0.5, ρ = 400 | centre 0.26 … 3.95 against 3.75; error not smaller than at ρ = 100 (one seed: noise or bias, open) |
| kick f = 0.08 for 1 ≤ t < 3, paired with f = 0 on the same sprinklings, three seeds | displacement at t = 8.5 is 0.67, 0.73 and 0.64 against 0.82 for M = 8. This is only indicative, because the readout is not calibrated |
| kink at rest, three seeds | drifted +0.20, +0.33 and +0.22 by t = 8.5, all in the same direction: noise or bias, H5-0 must resolve it |

**Reading:** with the order's own Green function, the kink holds its charge, keeps its launch velocity within about
0.05 over 7 units, and answers a force in the right direction with roughly the right size. The interaction sum over
the causal past is what holds it together. Without it (V = 0) the initial data split into two half-kinks running at
±c (d'Alembert's formula).

### Consequence for D2: a new decision for the author

| | |
|---|---|
| **D2′ (proposed)** | Implement the conservative dynamics through the **retarded Green function of the order**, not through the BD operator. In the continuum the Green function is the inverse of the wave operator, so the physics is the same. In 1+1 it is the causal relation itself. P-a is recorded in H5-0 as a documented failure of the instrument. |
| The cost in higher dimensions | Johnston's massless propagator in 3+1 is proportional to the **link** matrix. The links are unbounded along the light cone (the problem of F5.1, criterion (d)), so the difficulty returns there. 2+1 has no equally simple form. This belongs to the later stages, not to H5 in 1+1. |
| Criterion (d) in 1+1 | The sum is finite because the initial slice bounds the past. The slice is part of the declared launch (as in any initial-value problem). For an eternal static kink the 2D massless sum over the infinite past would diverge. This is an infrared property of 1+1, recorded as a limit. |

## 4. Definitions of the two masses (condition C2)

**Equation (continuum, units m = 1):**

  φ_tt − φ_xx + μ sin φ = −f χ(t),  where χ(t) = 1 for t₁ ≤ t < t₂ and 0 otherwise; μ = 1 in the base.

**Rest mass, from the profile at rest:**
- M_rest = E = ∫ [½ φ_t² + ½ φ_x² + μ(1 − cos φ)] dx, read on slabs of the kink at rest (H5-A), from the field
  binned in x (readout with declared coordinates);
- continuum value: for φ = 4 arctan e^{√μ x}, ½φ_x² = μ(1 − cos φ) = 2μ sech²(√μ x), so E = 8√μ = 8.

**Force, from the input and the measured charge (not from the kink's equation of motion):**
- field momentum P = −∫ φ_t φ_x dx;
- from the field equation alone,
  dP/dt = −∫ (φ_tt φ_x + φ_t φ_tx) dx = −∫ (φ_xx − μ sin φ − fχ) φ_x dx − ∫ ∂_x(½ φ_t²) dx = 2π Q f χ,
  because the boundary terms vanish (φ_x → 0 and cos φ equal at both ends) and ∫ φ_x dx = 2πQ;
- so F := 2π Q f, where f is the input and Q the measured integer charge.

This is an identity of the field equation (the momentum balance), not a law of motion of the kink. On the order it
holds only approximately: if the discrete dynamics does not conserve momentum, M_in ≠ M_rest. That is part of what
the test measures.

**Inertial mass, from the trajectory:**
- M_in = F τ_F / Δ(γv), where τ_F = t₂ − t₁ and Δ(γv) comes from the slope of the centre after the kick against the
  slope before it, paired with the unkicked run on the same sprinklings;
- velocities stay small (Δ(γv) ≈ 0.08 for f = 0.05, τ_F = 2), so the relativistic correction is small but included;
- continuum value: Δ(γv) = 2πf τ_F / 8.

**Independence:**
- M_rest uses only the profile at rest;
- M_in uses only the input f, the integer Q and the kinematics of the centre;
- neither number enters the other.

As a diagnostic (decides nothing), P is also read directly, to see whether the discrete momentum balance holds.

**Side effect of the force:** during the window the vacuum shifts to sin φ₀ = −f/μ. After the switch it oscillates
with amplitude ~f and frequency √μ, uniformly. This is reported, and it is kept small (f ≤ 0.05).

## 5. The model (base, before H5-D)

- **Substrate:** 1+1 Poisson sprinkling with density ρ, box |x| ≤ X, t ∈ [0, T]. Lengths are in units of the kink
  width (μ = 1).
- **Field:** φ on all events (D1). The internal target space is an angle: the vacuum is φ ∈ 2πℤ and the charge
  Q = [φ(+) − φ(−)] / 2π.
- **Dynamics (D2′):**
  φ(x) = φ_0(x) − (1/2ρ) Σ_{y ≺ x, t_y ≥ 0} [μ sin φ(y) + f χ(t_y)].
  - The relation y ≺ x is the order.
  - The cut t_y ≥ 0 is the initial slice of the launch.
  - Events are computed in any linear extension.
- **Initial data (launch):** φ_0(x) = ½[g(x − t) + g(x + t)] + ½ ∫_{x−t}^{x+t} h, with g and h the profile and the
  time derivative of the kink with velocity v at t = 0 (exact continuum functions evaluated at the coordinates of the
  event).
  - This is the free wave of the initial data on the slice, the only place where coordinates enter the dynamics. It is
    the initial condition, not a source.
  - After t = 0 there is no external term except fχ in H5-C, which is exactly zero outside [t₁, t₂).
- **Potential:** μ(1 − cos φ), with μ explicit (C3). The null V = 0 is a localisation control of the base model, not
  the capacity test.

## 6. Order of the experiment (the author's)

### H5-0. Instrumentation

These conditions must hold before any kink is read.
- **I1 linear stability of the vacuum:** initial data of small random smooth modes (amplitude 0.01, k ≤ 4). The rms
  amplitude at the end divided by that at the start must be ≤ 1.2; in the continuum it is bounded, with ω² = k² + 1.
- **I2 dispersion:** plane-wave packets with k ∈ {0.5, 1, 2, 4}. The phase evolution gives ω(k), and
  |ω_meas / √(k² + 1) − 1| must stay below the tolerance. The group velocity of the packet is reported.
- **I3 convergence:** ρ ∈ {100, 200, 400}. The rms error of the kink at rest against the exact solution must not grow
  with ρ. The expected trend is ∝ ρ^(−1/2), reported. The open point of P-b (no improvement from 100 to 400) is
  resolved here.
- **I4 order and layers:**
  - two linear extensions (sorted by t, and by the boosted time t + 0.5x) give the same field (relative 1e-12);
  - so does computing in time layers against all at once;
  - two box sizes agree on the events of their common domain of dependence (relative 1e-12).
- **I5 record:** P-a (BD in retarded form) is reported as the documented failure of the differential form.

### H5-A. Static kink

Over the run: Q = 1 in every slab; drift of the centre |Δx_c| + 2 SE < δ_x; width w constant within δ_w; E constant
within δ_E; and M_rest = E is read here. The common drift seen in P-b is resolved here.

### H5-B. Free motion

v = 0.3 and 0.6, launched on the slice, then nothing. Required:
- |v_end − v_launch| + 2 SE < δ_v;
- w·γ = w₀ within δ_w (Lorentz contraction);
- E = γ M_rest within δ_E;
- the identity holds.

### H5-C. Impulse

f = 0.05 for τ_F = 2, then exactly zero; paired with the unkicked run on the same sprinklings. Read:
- M_in = 2πQf τ_F / Δ(γv), against M_rest from H5-A;
- the free motion after the kick (velocity kept).

### H5-E. Relaxational null

The same potential, initial data and kick, with an overdamped (telegraph) dynamics φ_tt + Γφ_t − φ_xx + V′ = −fχ,
Γ ∈ {2, 5}. Its retarded Green function in 1+1 is ½ e^{−Γt/2} I₀((Γ/2)√(t² − x²)) inside the cone.
- It uses the coordinate time of a medium frame. That is allowed for a null model and declared, like G1 in F5.1g.
- Expected: the launched kink slows and stops; under the kick, v follows f and decays after the switch.
- The instrument must separate H5-E from H5-B and H5-C by more than 2 SE.

### Controls of the base model

- **V = 0:** no localisation; the half-kinks run off at ±c.
- **Relation-blind sum:** each event sums over a random set of earlier events of the same size as its causal past.
  This destroys the causal structure, and the kink must not propagate as a wave.
- **Continuum reference:** finite differences of the same equation.

### H5-D. Capacity interpretation (last, separate)

Derive μ (or the form of V) from the capacity bound without tuning, then check whether H5-A to H5-C still hold with it.

**A point for the author:** K₀ (a finite degree) was the bound of the first version, and the laboratory has excluded it
(audit, section 0). The present axiom is the **finite capacity of the flow through an event**. H5-D should start from
that axiom unless the author decides otherwise.

## 7. Operational identity (readouts; coordinates allowed and declared)

Per slab of width 1, inside the domain of dependence |x| < X − t − 2:
- Q from the plateaus at both ends of the readable part;
- x_c, the centroid of the binned increments of φ (the topological density);
- w and the fit residual from a fit of the kink profile;
- E in a window co-moving with x_c, and the energy flux across its edges;
- v as the slope of x_c over slabs.

**Identity [D(t₁)] = [D(t₂)]:**
- Q is the same integer;
- w·γ is constant within δ_w;
- the residual does not grow;
- E in the co-moving window is constant within δ_E.

**Motion** is x_c moving with the identity kept. **Deformation** is a failure of the identity, whatever x_c does.

## 8. Criteria

The shape is fixed now. The numbers are fixed after the H5-0 development rows, before any full run.

- **INSTRUMENT_FAILED:** any of I1–I4 fails, or H5-A does not keep Q and w, or H5-E is not separated from H5-B and H5-C.
- **PASSED:**
  - H5-B: |v_end − v_launch| + 2 SE < δ_v, contraction within δ_w, identity kept;
  - H5-C: |M_in / M_rest − 1| + 2 SE < δ_M;
  - the V = 0 and relation-blind controls lose the identity.
- **FAILED:**
  - inertia: H5-B loses more than δ_v of its velocity (|Δv| − 2 SE > δ_v);
  - deformation: the identity fails while moving;
  - mass: |M_in / M_rest − 1| − 2 SE > δ_M.
- **UNDECIDED** otherwise.

Proposed numbers (to confirm after H5-0):
- tolerances δ_v = 0.03, δ_w = 0.1, δ_E = 0.1, δ_x = 0.1, δ_M = 0.1;
- ρ ∈ {100, 200, 400}, X = 20, T = 10, 16 sprinklings.

## 9. Reading declared now

**PASSED (1+1):** the order's own retarded Green function, a sum over the causal past, gives a defect with inertia.
Its inertial mass equals the energy of its configuration, and an overdamped dynamics does not do this. The inertia
of matter is then a property of a conservative dynamics of the order, not of a flow that carries the defect. The next
steps:
- 2+1 and 3+1, with the links problem of the propagator and Derrick's theorem;
- F5.1i: the gravitational mass of the same defect.

**FAILED:** the discreteness breaks inertia or the equality of the masses at the tested scale. The diagnostics (which
condition, the dependence on ρ) say whether this vanishes for w ≫ ℓ.

**Limits:**
- 1+1 has no transverse directions;
- the potential is chosen, not derived;
- the internal target space is a new assumption (D1);
- the launch uses coordinates (the free wave of the initial data);
- even PASSED says nothing about gravity.

## 10. New against the axioms (to be recorded in the audit when approved)

- **The field φ on every event:** a new state of the order. **Input.**
- **The retarded sum over the causal past:** passes (a), (b) and (c). For (d), the initial slice is the declared
  regulator in 1+1.
- **The flow of the Axiom of continuity as the conserved current of this dynamics:** an interpretation, for the
  author.
- **The potential as finite capacity:** H5-D, tested last.
