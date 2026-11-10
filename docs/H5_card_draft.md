# Card H5: inertia of a locked causal defect

Revision 4, 10 Nov 2026: the criteria of section 8 with the author's four additions, for freezing. Test ID F3.14a.

- Revision 1: 10 Nov 2026.
- Revision 2: 10 Nov 2026, after the author's decisions D1–D4, his three technical conditions and the order of the
  experiment he set.
- Revision 3: after the decisions of 10 Nov, 11:07 (section 1).
- Revision 4: after the author's review of the numbers (10 Nov, 11:52): I3 defined; the H5-E interval and metric
  fixed; V = 0 measured, not predicted; blind separation defined; general rules (section 8.0).

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

**Decisions of 10 Nov, 11:07:**

| point | decision |
|---|---|
| D2′ | the retarded Green function of the order is the base of the conservative dynamics of H5 in 1+1 |
| BD operator | stays in the laboratory as a separate instrument for other tests; it is not the dynamical base of this test |
| 3+1 | postponed: Johnston's propagator is supported on links, and the light-cone problem returns |
| H5-D | the potential is to be derived from the Axiom of finite capacity of the flow through events |
| K₀ | removed from H5 and from its physical interpretation: it belongs to a superseded version of the model and is not a parameter of the theory |

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
the potential from the bound on the flow through events without tuning to results. If that fails, the inertia result
can still stand; the interpretation cannot.

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

### P-c. Development rows of H5-0 to H5-C (10 Nov; `diag/h5_instrument.py`; 16 sprinklings per density)

The field is sprinkled in the domain of dependence of the initial slice (|x| ≤ X − t), so every computed value is
exact for its order and no wall enters the dynamics. Fenwick engine; 66 s in the container. Full numbers are in
findings, "H5-0 development rows".

| | ρ = 100 | ρ = 200 | ρ = 400 |
|---|---|---|---|
| I4 order / engine / region (largest relative difference) | ≤ 1.4e-15 | | |
| I1 vacuum rms last/first (exact) | 1.041 (1.033) | 0.988 (0.982) | 0.971 (0.967) |
| I1 relative error, last slab | 0.150 | 0.107 | 0.074 |
| I2 \|ω/ω_cont − 1\|, every k | ≤ 0.001 | ≤ 0.0004 | ≤ 0.0002 |
| I3 rms error of the rest kink | 0.121 | 0.116 | 0.046 |
| H5-A drift of the centre over 9 units | +0.11 ± 0.09 | +0.07 ± 0.09 | −0.03 ± 0.04 |
| H5-A E (same readout on the exact field ≈ 8.01) | 8.093 | 8.058 | 8.026 |
| H5-B v (launch 0.3) | 0.318 ± 0.014 | 0.302 ± 0.009 | 0.302 ± 0.003 |
| H5-B v (launch 0.6) | 0.597 ± 0.007 | 0.597 ± 0.006 | 0.599 ± 0.003 |
| H5-B w·γ | 0.995–1.003 | 0.999–1.004 | 0.999–1.001 |
| H5-C M_in = 2πQfτ / Δ(γv) | 8.62 | 8.34 | 8.13 |
| H5-C M_in / E | 1.065 | 1.035 | 1.013 |

The common drift of the rest kink seen in probe P-b was noise: it falls with the density. The SE of M_in / E is about
0.035–0.06.

**A limit of the content (the lesson of F3.13c), stated before the test.** The sprinkling is Lorentz invariant and the
rule is covariant, so a moving kink is, in distribution, a boosted resting kink: once H5-A holds, H5-B is largely
fixed by symmetry. M_in = E follows from the momentum balance of the field, which holds on average for a rule
invariant in distribution. What is **not** fixed by symmetry:
- that the discrete sum is stable at all (BD in retarded form is not);
- the size and the sign of the discreteness corrections (M_in / E − 1 falls roughly as ρ^(−1) in these rows);
- the separation from a relaxational dynamics (H5-E).

A PASSED of H5 in 1+1 therefore says: *the order has a conservative dynamics of its own, and its defect has inertia
with M_in = E up to discreteness corrections of a measured size*. It does not say that inertia is explained beyond
Lorentz invariance. That content lies in F5.1i (gravitational mass of the same defect) and in 3+1.

### P-d. Development rows of the controls (10 Nov; 16 sprinklings, ρ = 100 and 400; `--controls`)

| arm (launch v = 0.3 unless kick) | p = γv over t ∈ [5, 10), ρ = 400 | conservative | paired difference | localisation L (last/first) |
|---|---|---|---|---|
| telegraph Γ = 2, launched | +0.0009 | +0.3205 | −0.3195 ± 0.0079 | 0.714 |
| telegraph Γ = 5, launched | +0.0046 | +0.3205 | −0.3159 ± 0.0125 | 0.880 |
| telegraph Γ = 2, kick | +0.0008 | +0.0736 | −0.0727 ± 0.0070 | 1.009 |
| telegraph Γ = 5, kick | +0.0053 | +0.0736 | −0.0683 ± 0.0116 | 1.065 |
| V = 0 | (no track) | +0.3205 | — | 0.000 (conservative 1.003) |
| causality-blind sum | (no track) | +0.3205 | — | 0.376 ± 0.004 (conservative 1.003) |

- **Telegraph:** the damped kink stops within the first units (continuum: dP/dt = −ΓP exactly). The separation from the
  conservative arms is 10–40 SE at 16 sprinklings.
- **V = 0:** the sum term is identically zero, so the field is exactly the free wave of the initial data (d'Alembert).
  This control therefore tests the **readout**, which must detect a field that is not localised. It does not test the
  order. The fitted width grows from 1.0 to 150 in 9 units.
- **Causality-blind sum:** the fitted width grows from 1.0 to 80, and the energy in the window first rises (the
  non-local mixing pumps energy in) and then leaves.
- In both controls the fit "succeeds" with an absurd width. Fit failure is therefore not a usable criterion; the
  localisation L is.

The frozen criteria are in section 8.


### Consequence for D2 (accepted by the author on 10 Nov as D2′)

| | |
|---|---|
| **D2′ (accepted)** | Implement the conservative dynamics through the **retarded Green function of the order**, not through the BD operator. In the continuum the Green function is the inverse of the wave operator, so the physics is the same. In 1+1 it is the causal relation itself. P-a is recorded in H5-0 as a documented failure of the retarded differential form. The BD operator stays in the laboratory for other tests. |
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

- **V = 0:** a check of the readout. Without the potential the sum vanishes and the field is the free wave of the initial
  data; the loss of localisation is measured (section 8.2), not predicted.
- **Relation-blind sum:** each event sums over a random set of earlier events of the same size as its causal past.
  This destroys the causal structure, and the kink must not propagate as a wave.
- **Continuum reference:** finite differences of the same equation.

### H5-D. The effective potential from the finite capacity of the flow (last, separate)

**Statement.** The potential of the field is neither a freely chosen function nor a function of a finite degree of an
event. Its form is to be **derived** from the bound on the flow through events in the causal order.

**This is a task, not yet a derivation: a concrete mechanism has to be shown.** One route, to be worked out:
- the bound F_v ≤ C_v (C_v the admissible capacity of the event) enters a variational problem through a Lagrange
  multiplier λ_v ≥ 0, which acts as a local effective potential;
- the dependence of that potential on φ must follow from the definition of the flow (how φ enters F_v, e.g. through
  the conserved current of the dynamics) and of the capacity, and must not be chosen so that sine-Gordon comes out.

**A hard bound need not give a smooth potential.** It may give only:
- a one-sided condition (a wall);
- a non-linear response;
- an effective potential only after coarse-graining.

The computation decides which.

**Order of work:**
1. Define the flow of the field from the conserved current of the H5 dynamics (Axiom of continuity), and the capacity
   from the axiom. No tuning to the result of H5-A to H5-C.
2. Derive the effective potential: analytically where possible, otherwise by coarse-graining a discrete computation.
3. **Comparison, not fitting:** set the derived potential beside the analytic sine-Gordon potential. Is it periodic
   (does it admit a topological charge at all)? Is it bounded? What is its curvature at the vacuum?
4. If it admits localised configurations, repeat H5-A to H5-C with it.

**Reading.** If the sine-Gordon potential does not come out of the capacity axiom, this is **not** a failure of the
instrument. It is a result about which kind of localised configuration the theory actually admits. The inertia result
of the base model (H5-A to H5-C) does not depend on H5-D.

H5-D gets its own criteria in a separate card after H5-0 to H5-C, because its form depends on the mechanism found in
step 2.

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

## 8. Frozen criteria (revision 4, for freezing by the author)

Test ID **F3.14a**: *Inertia of a localised defect of a field of the causal order*.
- Branch: a defect of a field of the order (separate from F3.13, which tested a topological defect of an orientation
  carried by a flow).
- Hypothesis: H5 of `docs/cosmogenesis.md`.

### 8.0 General rules (the author's methodological note, binding)

- **Sprinklings:** 32 new sprinklings **per density**, ρ ∈ {100, 200, 400} (the same number at each density), seed
  label "F3.14a". Every sprinkling is used by every arm at its density, so the arms are paired.
- **Deciding density** ρ = 400. The values at ρ = 100 and 200 are reported and enter only I3.
- **Controls** (H5-E, V = 0, blind) run at ρ = 400 only, on the same 32 sprinklings.
- **Statistics:** SE = standard deviation over sprinklings / √32. Paired quantities are differences on the same
  sprinkling. The SE of a ratio of means comes from a bootstrap over sprinklings (2000 resamples, fixed seed).
- **Every condition must hold on its own.** Conditions are not averaged into a score, a failed condition is not
  balanced by another, and nothing is chosen after the run among densities, Γ, velocities or readouts. Only the
  quantities defined here decide.
- **Slabs:** width 1, centres t = 0.5 … 9.5 (10 slabs). All events of the trapezoid are exact for their order, so
  there is no excluded region.

### 8.1 Readouts (fixed)

- **Kink fit per slab:** φ = 2πn₀ + 4 arctan exp((x − x_c − v_s(t − t_m))/w), on the events with |x − x_g| < 4, with
  x_g = v_launch · t_m (0 for rest and kick).
- **Q per slab:** (mean φ on x_g + 5 … x_g + 7 − mean φ on x_g − 7 … x_g − 5) / 2π.
- **Energy per slab:** E = ∫ [φ_t²/2 + φ_x²/2 + μ(1 − cos φ)] dx over |x − x_c| < 5, from local regressions in cells
  0.5 × 0.5, step 0.25.
- **Localisation per slab:** E_loc = the same energy over |x − x_g| < 3, with the model's own μ (μ = 0 for V = 0). It
  does not use the fit.
- **Velocity:** the slope of x_c over the slabs named in each condition. p = γv.
- **Exact readouts:** the same readouts applied to the exact continuum solution at the same events: Q_ex, w_ex, E_ex.

### 8.2 Instrument (INSTRUMENT_FAILED if any condition fails)

- **I1 vacuum stability.** Per sprinkling, R = rms(φ, last slab) / rms(φ, first slab) and R_c, the same ratio of the
  exact linear solution at the same events. δ₁ = R / R_c − 1, the relative deviation of the amplitude ratio from the
  continuum. Required: |mean δ₁| + 2 SE < 0.1, and the mean of rms(φ − φ_ex) / rms(φ_ex) in the last slab < 0.2.
  Initial data: 8 modes, k uniform in (0.25, 4], total amplitude 0.01.
- **I2 dispersion.** For **each** k ∈ {0.5, 1, 2, 4} separately: ω_meas is the slope of the unwrapped phase of the
  projection on cos kx, sin kx, slabs of 0.25 over [0, 10]. δ₂ = ω_meas / √(k² + 1) − 1, the relative error of the
  frequency. Required: |mean δ₂| + 2 SE < 0.005.
- **I3 convergence.** Per sprinkling of the rest arm, ε(ρ) = mean over the 10 slabs of |w − w_ex| / w_ex, with w and
  w_ex from the same fit. Required: ε̄(100) − ε̄(400) > 2 √(SE₁₀₀² + SE₄₀₀²) (independent sprinklings). ε̄(200) is
  reported, with whether the trend is monotonic.
- **I4 order and layers.** Each comparison uses the **same sprinkling, the same initial data and the same parameters**
  (μ = 1, f = 0.05, t ∈ [1, 3), kink v = 0.3). It applies to 4 sprinklings at ρ = 100, X = 12, T = 6:
  - linear extension t against t + 0.5x;
  - the Fenwick engine against the brute-force block engine;
  - the trapezoid X against X − 2, on their common events.

  Required: largest relative difference < 1e-12 in each comparison.
- **H5-A static kink** (ρ = 400; all four conditions):
  - |Q̄ − Q̄_ex| + 2 SE < 0.02;
  - drift Δx = x_c(slab 10) − x_c(slab 1), with |mean| + 2 SE < 0.1;
  - |w̄ / w̄_ex − 1| + 2 SE < 0.1;
  - |Ē − Ē_ex| / 8 + 2 SE < 0.05.

  Here ¯ means the mean over the slabs of one sprinkling, and SE is over sprinklings. M_rest := Ē of this arm.
- **H5-E relaxational null** (telegraph φ_tt + Γφ_t − φ_xx + sin φ = −fχ, same initial data and kick; implemented
  through ψ = e^{Γt/2}φ with the same sum over the causal past and the coordinate time of a medium frame, declared).
  - For **each** Γ ∈ {2, 5} and **each** of the two arms (launched v = 0.3; kick), p_Γ = γv over the slabs t ∈ [5, 10)
    (centres 5.5 … 9.5), paired with p_cons, the same quantity of the conservative arm (H5-B v = 0.3, or H5-C kick).
  - Required for each of the four cases: (p_cons − p_Γ) − 2 SE_Δp > 0, **and** the minimal drop p̄_Γ ≤ 0.5 p̄_cons.
  - Continuum expectation: dP/dt = −ΓP, so p_Γ ≈ 0 by t = 5.
  - If either requirement fails, the instrument cannot tell inertia from relaxation, and no conclusion about inertia is
    drawn.
- **V = 0 readout check** (launched v = 0.3). Here the sum term vanishes and φ is the free wave of the initial data.
  The control checks that the readout recognises a field that is not localised. Per sprinkling,
  L = E_loc(slab 10) / E_loc(slab 1), paired with the conservative v = 0.3 arm: D_L = L_cons − L_V0. Required:
  D_L − 2 SE > 0.2.

### 8.3 Inertia (ρ = 400)

**H5-B**, for v = 0.3 and v = 0.6 **separately**. All conditions must hold for both.
- velocity: v_meas = slope of x_c over all 10 slabs;
- contraction: C = w̄(v) · γ(v) / w̄(rest), paired on the same sprinkling;
- energy: e_B = (Ē(v) − Ē_ex(v)) / (8γ);
- identity: |Q̄ − Q̄_ex| + 2 SE < 0.02.

**H5-C** (f = 0.05, t ∈ [1, 3), τ = 2):
- Δp = γv(kick) − γv(rest), with the velocities as slopes over the slabs t ∈ [4, 10) and paired;
- M_in = 2π Q̄_kick f τ / mean(Δp);
- r_M = M_in / M_rest, with SE by bootstrap over the pairs.

Reported only (decides nothing): the free motion after the kick (γv over [4, 7) against [7, 10)), the field momentum
P, and r_M at ρ = 100 and 200 (the trend).

**Causality-blind control** (launched v = 0.3).
- Each event sums over m_x events drawn uniformly, with replacement, from all earlier events. m_x is the number of
  events of its causal past, counted by the same engine with unit weights; earlier events cover the same time range
  0 ≤ t_y < t_x as the causal past.
- The normalisation 1/(2ρ), φ₀ and the parameters are the same. Only the index set of the sum differs.
- Separation: D_L = L_cons − L_blind (paired), required D_L − 2 SE > 0.2. This is one metric, fixed now. The width and
  the velocity are reported.

### 8.4 Verdict (checked in this order)

1. **INSTRUMENT_FAILED:** any condition of 8.2 fails.
2. **FAILED**, if any one of these holds:
   - H5-B at v = 0.3 or 0.6: |mean v_meas − v| − 2 SE > 0.03 (inertia), or |mean C − 1| − 2 SE > 0.1 (contraction), or
     |mean e_B| − 2 SE > 0.05 (energy);
   - H5-C: |r_M − 1| − 2 SE > 0.1 (inertial and rest mass differ).
3. **PASSED:** all of the following hold:
   - H5-B at both velocities: |mean v_meas − v| + 2 SE < 0.03, |mean C − 1| + 2 SE < 0.1, |mean e_B| + 2 SE < 0.05,
     and the identity condition;
   - H5-C: |r_M − 1| + 2 SE < 0.1;
   - the blind control separated.
4. **UNDECIDED:** everything else.

If the blind control is **not** separated while H5-B and H5-C pass, the verdict is UNDECIDED, with the reading that the
sum and not the causal structure carries the kink.

The tolerance 0.1 of H5-C is the resolving power of this experiment, fixed in advance. It is not a claim about the
precision of the theory.


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
- in the base model the potential is chosen, not derived (H5-D derives it from the capacity of the flow);
- the internal target space is a new assumption (D1);
- the launch uses coordinates (the free wave of the initial data);
- even PASSED says nothing about gravity.

## 10. New against the axioms (to be recorded in the audit when approved)

- **The field φ on every event:** a new state of the order. **Input.**
- **The retarded sum over the causal past:** passes (a), (b) and (c). For (d), the initial slice is the declared
  regulator in 1+1.
- **The flow of the Axiom of continuity as the conserved current of this dynamics:** an interpretation, for the
  author.
- **The potential from the finite capacity of the flow:** H5-D, derived and compared last; K₀ plays no part.
