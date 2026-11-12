# Card F5.2a: free fall of a field defect in a 1+1 order of nonuniform density

**FROZEN** by the author on 10 Nov 2026, 13:53 (revision 2). Test file `phases/F5_gravity/F5_02a_field_defect_fall.py`;
helpers `phases/F5_gravity/field1d_curved.py`. No number or condition of section 5b changes after this point (R132).

**Result of the revision: F5.2a2 PASSED** on 10 Nov 2026 (seed 0, run by the author, 463 s; findings, "Holds"). Not yet
replicated.

**Result and implementation error (10 Nov, after the run).**
- F5.2a ran (run by the author, 128 sprinklings at ρ = 800, 326 s) and came out UNDECIDED: everything held except the
  fall of M 16, a = 1.030 ± 0.014.
- An implementation error was then found: the sprinkling placed the points in the whole box [0, T] × [−X, X] and
  dropped those outside the trapezoid. The density was therefore 0.808 ρ Ω², while the rule normalises with ρ. Every
  arm ran with an effective coupling 0.808 μ: the kinks were 11 % wider than their continuum references.
- It showed as the launched arm: D = +0.72 against +0.42 in the continuum, unchanged from ρ = 400 to 6400. A
  deterministic quadrature of the same rule agreed with the continuum.
- **F5.2a gives no conclusion** (findings, "Instrument failures").
- **The revision F5.2a2** (`phases/F5_gravity/F5_02a2_field_defect_fall_fixed.py`) uses the corrected sprinkling
  `sprinkle_curved2` and new seeds. Every number and condition of section 5b is unchanged.
- Development check of the correction (16 sprinklings, ρ = 800, decides nothing): D / D_cont = 0.996 ± 0.005 (M 4),
  0.994 ± 0.016 (M 8), 1.022 ± 0.040 (M 16); the launched arm 0.95 ± 0.09.
- The development rows of section 4b were taken with the same error.

**Note after freezing (decides nothing; observed during the smoke test, 10 Nov; explained by the implementation error
above).** The reported launched arm starts
from a translated, non-boosted profile, a state far from equilibrium. Its discrete trajectory converges only slowly to
the continuum one, even in a flat order: at t = 9.5 it reaches +3.23 at ρ = 400 and +2.92 at ρ = 1600, against +2.73
in the continuum. The short-wavelength part of the shed radiation is where the order differs from the continuum.
The arm stays reported only, as frozen. The resting kinks of the decision are near equilibrium and fell within 2 % in
the development rows.

- Revision 1 (10 Nov): approved by the author at 12:30. He accepted the ID F5.2a, the profile, the three masses and the
  counting null as an instrument condition.
- Development rows decide nothing.

## 0. ID and place

**Proposed ID: F5.2a.** F5.1 to F5.1g are the series about the fall and inertia of flows and histories. Most of their
rules are proxies (F5.1c/d/e), and F5.1g failed on M5. The fall of a bound defect is a new series:
- **F5.2a:** passive gravitational mass. Does a defect fall, and do all defects fall alike?
- **F5.2b (later):** active gravitational mass. Does a defect compress the order in proportion to its energy? This
  needs a rule for the back-reaction of the field on the order, a decision of the theory, and is meaningful only in
  3+1, because in 1+1 the Einstein tensor vanishes identically.

The planned row F5.1i ("equivalence of inertia and gravitational response") is carried over to the F5.2 series.
Hypothesis H5 of `docs/cosmogenesis.md` (the gravitational half). The instrument is that of F3.14a. Its physical
conditions were met in both runs; its replication failed only on the drift tolerance (findings).

## 1. The idea: in 1+1, gravity is the number of events alone

Every 1+1 metric is conformally flat: ds² = Ω²(t, x)(dt² − dx²). The conformal factor does not change the light cones,
so **the causal order is the same as in flat spacetime; only the number of events changes**. A Poisson sprinkling of
constant proper density ρ has the coordinate density ρΩ². "Order + number = geometry" holds here literally: the
geometry, and with it gravity, is carried only by how many events lie where.

**The rule of F3.14a is not changed at all.**

  φ(x) = φ_0(x) − (1/2ρ) Σ_{y ≺ x, t_y ≥ 0} V′(φ(y))

In the continuum the sum (1/ρ) Σ approximates ∫ V′ √−g d²x = ∫ V′ Ω² dt dx. The massless wave operator is conformally
invariant in 1+1, so the Green function ½ × the causal relation stays exact. The same rule, on an order of
nonuniform density, is therefore the field equation in the curved spacetime. **The rule knows nothing about the metric.
If the kink falls, it falls because of the counting.**

**A sign specific to 1+1, stated in advance.** In the conformal gauge the density of events (∝ Ω²) and the rate of
clocks (dτ = Ω dt) are the same function. A body falls towards lower Φ = ln Ω, i.e. towards the **sparser** order and
the slower clocks. In 3+1 the relation of density and clock rate is different (in the weak field √−g ≈ 1 − 2φ: denser
near a mass, with slower clocks). The 1+1 sign is therefore not a test of the throughput principle.

## 2. The hypothesis and the continuum reference

**F5.2a.** In a 1+1 order of nonuniform event density, a defect of the field of F3.14a:
- falls along the geodesic of the metric that the density defines;
- falls with the same acceleration whatever its mass (universality of free fall, the weak equivalence principle);
- does not fall when the density is uniform, or when the sum is made blind to the number of events.

**Continuum reference.** For ds² = e^{2Φ(x)}(dt² − dx²) the geodesic obeys d²x/dt² = −Φ′(x)(1 − (dx/dt)²). For a kink
of potential μ(1 − cos φ) the mass is M = 8√μ and the proper width is 1/√μ. All masses follow the same geodesic, up to
tidal corrections of order (w/L)², which are negligible here.

## 3. The model (proposed)

- **Profile:** Φ(x) = gL tanh(x/L), with g = 0.05 and L = 10.
  - Near x = 0 the field is about uniform, with acceleration ≈ g.
  - The curvature is R ∝ Φ″ ≠ 0, so the spacetime is not merely Rindler. In a uniform field Ω = e^{gx} the
    spacetime is flat (Rindler), and the test would only repeat F3.14a.
  - The coordinate density varies by a factor e^{±2gL} = e^{±1} across the box.
- **Sprinkling:** a Poisson process of coordinate density ρΩ²(x) in the trapezoid |x| ≤ X − t, with X = 26 and T = 10.
- **Masses:** μ ∈ {0.25, 1, 4}, so M = 4, 8, 16 and proper widths 2, 1, 0.5. All three are on the **same sprinkling**
  (paired).
- **Launch:** the kink at rest at x₀ = 0, with its profile in proper distance s = ∫Ω dx. Its initial data enter through
  the free wave φ_0, as in F3.14a. In curved spacetime a resting kink is not a solution: it starts to fall, which is
  the point.
- **Arms:**
  - curved (the three masses);
  - **flat control** (g = 0, the same masses): no fall;
  - **counting null:** the same curved sprinkling, but every term of the sum divided by the local coordinate density
    Ω²(y), i.e. the sum is made to measure coordinate area instead of counting events. The dynamics becomes flat and
    the kink must not fall. This uses coordinates (declared; a null model). It shows that the fall comes from the
    number of events;
  - reported only: a kink launched at v = 0.3 against the field (inertia and gravity together).
- **Readout:** the kink fit of F3.14a in coordinates, slab centres 0.5 … 9.5. The displacement D = x_c(9.5) − x_c(0.5)
  and the acceleration from a fit to the geodesic family (one free parameter, the scale of g) are read.

## 4. Development probe (10 Nov; decides nothing; ρ = 200, three sprinklings, g = 0.05, L = 10)

Displacement of the centre at t = 9.5 (geodesic −2.17):

| | seed 0 | seed 1 | seed 2 |
|---|---|---|---|
| μ = 1 (M = 8) | −1.77 | −2.03 | −2.33 |
| μ = 4 (M = 16) | −2.10 | −1.94 | −2.36 |
| flat (g = 0), μ = 1 | +0.26 | −0.56 | |

**Reading (development):** the kink falls, in the direction and of the size of the geodesic. The rule was not changed;
only the number of events per coordinate area was. The wandering seen in F3.14a (sd ~0.3 at ρ = 200 over 9 units) is
the main noise. The universality of free fall at the level of a few % therefore needs many sprinklings or a paired
design. The counting null was not run in the probe.

## 4b. Development rows (10 Nov; `diag/f52a_fall.py`; 16 sprinklings per density; decide nothing)

**References.** The geodesic from rest gives D_geo = −2.164 between slab centres 0.5 and 9.5.

The continuum field equation (finite differences, the same initial data, the same centre readout) gives D_cont:

| μ | M | D_cont | D_cont / D_geo |
|---|---|---|---|
| 0.25 | 4 | −2.098 | 0.970 |
| 1 | 8 | −2.145 | 0.991 |
| 4 | 16 | −2.158 | 0.997 |

**The continuum itself does not fall exactly along the geodesic.** The kink is an extended body, and its initial
profile is not a static solution (it breathes and radiates). The wide kink (μ = 0.25, w = 2) falls 3 % less than a
point. The raw ratio D_i / D_j therefore contains a continuum effect of up to 3 % that is not a violation of the
equivalence principle. **The fair reference for each body is its own continuum fall D_cont.**

**I4 on the curved order:** order, engine and region agree to ≤ 2.2e-15.

| | ρ = 200 | ρ = 400 | ρ = 800 |
|---|---|---|---|
| D / D_cont, μ = 0.25 | 1.001 ± 0.008 | 1.025 ± 0.009 | 1.016 ± 0.007 |
| D / D_cont, μ = 1 | 0.939 ± 0.016 | 1.002 ± 0.017 | 0.986 ± 0.016 |
| D / D_cont, μ = 4 | 0.903 ± 0.074 | 0.968 ± 0.043 | 0.994 ± 0.033 |
| sd of D per sprinkling, μ = 0.25 / 1 / 4 | 0.07 / 0.13 / 0.64 | 0.07 / 0.15 / 0.37 | 0.06 / 0.14 / 0.28 |
| flat control D, μ = 0.25 / 1 / 4 | −0.01 / −0.01 / −0.00 | +0.01 / +0.05 / −0.01 | +0.01 / +0.04 / +0.13 ± 0.12 |
| counting null D (μ = 1) | +0.09 ± 0.05 | −0.05 ± 0.04 | +0.02 ± 0.04 |
| raw η (M 4 vs 16) | +0.20 ± 0.10 | +0.06 ± 0.05 | +0.007 ± 0.029 |

- **The kink falls, and it falls because of the counting.** In the same order with the counting removed (the null) it
  does not fall: D = +0.02 ± 0.04 against −2.12.
- **At ρ = 800 every mass falls as its continuum body within 2 %.**
- **At low density the discreteness slows the fall:** a = 0.94 (μ = 1) and 0.90 (μ = 4) at ρ = 200. The effect is
  larger for the narrower kink, and it is gone by ρ = 400–800.
- **The narrow, heavy kink (μ = 4, w = 0.5) wanders most** (sd 0.28–0.64). It sets the number of sprinklings.
- **Launched at v = 0.3 against the field:** D = +0.68 to +0.75, against +0.54 for the geodesic. The initial data of
  this arm are a translated, not a boosted, profile, so the geodesic is not its reference. The test computes the
  continuum reference for the same initial data. The arm is reported only.

## 5b. Frozen criteria (revision 2, frozen 10 Nov 2026, 13:53)

Changes against revision 1, for approval:
- the reference of each body is its own **continuum fall D_cont** (the same initial data and readout), not the
  geodesic, because the continuum body differs from the geodesic by up to 3 % (section 4b); the geodesic ratio is
  reported;
- I1 and I2 of F3.14a are not repeated on the curved order (there is no exact solution there). The comparison with
  the continuum field equation is the direct check of the dynamics on the curved order, and I1 and I2 hold on the flat
  order (F3.14a, passed in both runs).

**General rules (as F3.14a):**
- deciding density ρ = 800, with **128 new sprinklings**, seed label "F5.2a";
- all arms on the same sprinkling of each seed (the flat arms on the flat sprinkling of the same seed);
- SE = sd / √128, paired where stated;
- every condition holds on its own; nothing is averaged and nothing is chosen after the run;
- ρ = 400 (32 sprinklings) is reported for the trend and decides nothing.

**Instrument (INSTRUMENT_FAILED if any fails):**
- **I4 on the curved order:** < 1e-12 (order, engine, region).
- **Flat control:** for each mass, |mean D| + 2 SE < 0.2 (about a tenth of the fall). Expected 2 SE ≈ 0.08 for μ = 4.
- **Counting null:** |mean D| + 2 SE < 0.2. If the kink still falls with the counting removed, the fall does not come
  from the number of events.
- **The fall is resolved:** for each mass, |mean D| − 2 SE > 1.0.

**Decision:**
- **Fall:** a_i = D_i / D_cont,i for each mass, with |a_i − 1| + 2 SE < **0.05**. Expected 2 SE: 0.005 / 0.012 / 0.033;
  largest development offset 0.016.
- **Equivalence (WEP):** η_ij = a_i / a_j − 1, paired, for the pairs (4, 8), (8, 16) and (4, 16), with |η_ij| + 2 SE <
  **0.06**. Expected 2 SE ≈ 0.02 for (4, 16); the development offset is up to 0.022.
- **PASSED:** all six conditions hold.
- **FAILED:** |a_i − 1| − 2 SE > 0.05 or |η_ij| − 2 SE > 0.06 for any mass or pair.
- **UNDECIDED** otherwise.

**Why these numbers:**
- **δ_a = 0.05** is above the worst expected 2 SE (0.033) plus the largest development offset (0.016).
- **δ_η = 0.06** is above 2 SE (0.02) plus the offset (0.022), with a margin of about 0.02.
- These tolerances are the resolving power of this experiment, not a claim about the theory (MICROSCOPE: η ~ 10⁻¹⁵).
- A discreteness effect smaller than 5 % cannot be seen here, and the test does not claim it.

**Reported:**
- the raw η (without the continuum correction);
- D / D_geo;
- the launched arm against its continuum reference;
- ρ = 400;
- the spread of D (the wandering) for each mass.

**Cost:** about 400 000 events per arm at ρ = 800, 8 arms × 128 sprinklings. Roughly 20–30 min with 8 processes.

## 5. Criteria (shape; the numbers after development rows, before any full run)

Learned from F3.14a: every tolerance must be set **above the measured resolving power**, from the spread in the
development rows with a margin, not from an estimate.

**Instrument (INSTRUMENT_FAILED if any fails):**
- the instrument conditions of F3.14a that do not depend on the drift tolerance: I1, I2 and I4 on the curved order;
- **flat control:** |mean D| + 2 SE < δ_flat for each mass;
- **counting null:** |mean D| + 2 SE < δ_flat. If the kink still falls with the counting removed, the fall does not
  come from the number of events.

**Decision (ρ = 400, N_s sprinklings, paired):**
- **geodesic:** a = D / D_geo for each mass. PASSED if |a − 1| + 2 SE < δ_geo for all three; FAILED if
  |a − 1| − 2 SE > δ_geo for any.
- **universality (WEP):** η_ij = D_i / D_j − 1, paired on the same sprinkling, for the pairs (4, 8), (8, 16) and
  (4, 16). PASSED if |η| + 2 SE < δ_η for all pairs; FAILED if |η| − 2 SE > δ_η for any.
- **UNDECIDED** otherwise.
- Every condition holds on its own: nothing is averaged into a score and nothing is chosen after the run.

**Tolerances (to be set from development rows):** δ_flat, δ_geo and δ_η, N_s, and possibly g (a stronger field
improves the ratio of signal to wandering).

## 6. Reading declared now

- **PASSED:** in an order whose only change is the number of events, a bound defect falls along the geodesic of the
  geometry that the number defines, and all defects fall alike. This is the weak equivalence principle derived from
  counting, not assumed. With F3.14a, inertia and passive gravitational mass come from the same rule.
- **FAILED (geodesic):** the discrete counting does not give the continuum fall at this scale.
- **FAILED (WEP):** defects of different mass or size fall differently. This would be a violation of the equivalence
  principle by discreteness, a physical prediction (MICROSCOPE bounds η ~ 10⁻¹⁵).

**Limits:**
- 1+1, with the opposite sign of density and clock rate compared with 3+1;
- the geometry is prescribed (no back-reaction: F5.2b);
- the potential is chosen (H5-D);
- the launch and the readout use coordinates, and so does the counting null (declared).

## 7. Decisions for the author

1. ID **F5.2a** (and F5.2b for the active mass), instead of F5.1i.
2. The profile Φ = gL tanh(x/L), g = 0.05, L = 10, or another one (a stronger g, or a localised well).
3. The three masses μ = 0.25, 1, 4, on the same sprinkling.
4. The counting null as an instrument condition.
5. The tolerances from development rows (as for F3.14a), with the rule that every tolerance sits above the measured
   resolving power.

## Revision F5.2a3 (frozen by the author on 10 Nov 2026, 19:09): the sample only

The replication of F5.2a2 came out UNDECIDED on the two WEP pairs with M = 16. Sixteen development sprinklings had
underestimated the wandering of the narrow kink (design error; findings).

**Changed, only the sample:** **512** sprinklings at ρ = 800 instead of 128, seed label "F5.2a3".

**Why this number.** The sd of η(4, 16) is 0.17–0.22 per sprinkling. With a margin of 1.5, 2 SE ≈ 0.03 < 0.06.

**Unchanged:** the model, the corrected sprinkling, the readouts, the references and every criterion of section 5b.

**Test file:** `phases/F5_gravity/F5_02a3_field_defect_fall_n512.py`. Smoke run OK (10 Nov).
