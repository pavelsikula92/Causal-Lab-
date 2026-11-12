# Card F5.2b: active gravitational mass from a causal growth rule of the order (1+1, scalar response)

Revision 1, 11 Nov 2026. Run by the author unchanged on 11 Nov, which froze it (R132).

**Result:** PASSED (seed 0, 5065 s): r_B = 1.029 / 1.017 / 1.009 for M_A = 4 / 8 / 16; u = −0.0007 ± 0.0081;
r_sep = 1.014 ± 0.008; counted proper times near matter 0.62 of the coordinate ones. Replication pending.

Originally a draft for freezing by the author. Direction A was approved on 11 Nov at 11:59, with the
author's corrections. Test file `phases/F5_gravity/F5_02b_growth_active_mass.py`; helpers
`phases/F5_gravity/growth1d.py`; smoke run OK. Development rows: findings, Discoveries, "Growth probe" (11 Nov).

## 1. Hypothesis

A causal growth rule of the order, together with the unchanged rule of the field of F3.14a, gives an active
gravitational mass equal to the energy of the source. In 1+1 the response is scalar.

**Growth rule (candidate for the Axiom of the dynamics of the order).**
- Every candidate event x of a Poisson sprinkling of proper density ρ carries an independent mark u_x, fixed at the
  sprinkling.
- x is realised iff u_x < Π(x)² = e^{2L(x)}, where

  L(x) = −(κ/2ρ) Σ_{realised y ≺ x} q(y),   q = 2V_A(φ_A) + 2V_B(φ_B).

  q is the load, the trace of T per event (findings, "Load candidates"). The factor ½ is the normalisation of the
  retarded Green function of 1+1.
- **The author's correction:** the realised events have coordinate density ρΠ², so the sum carries the factor Π² itself.

**Continuum law.** (∂_t² − ∂_x²) L = −κΠ² q, i.e. □_g L = −κq in ds² = Π²(dt² − dx²). For the fields,
(∂_t² − ∂_x²) φ_i = −Π² V_i′(φ_i).

**1+1 (stated in advance).**
- The metric is conformally flat, so a Poisson process of density ρΠ² **is** a sprinkling of it.
- The causal relation stays that of Minkowski space. The geometry lies in the counts (proper times), not in the cones.
- In 2+1 the relation itself will have to change.

**The coupling to the test defect is not tuned.** Both fields use the unchanged rule on the realised events. That
counting carries the fall is F5.2a3, replicated.

## 2. Model (frozen with the card)

- ρ = 800, κ = 0.003, T = 20, trapezoid |x| ≤ T + 16 − t.
- Two separate fields, A (μ_A) and B (μ_B), interacting only through the order. Kinks at rest at −3 (A) and +3 (B).
- Combinations (M_A, M_B): (4, 8), (8, 8), (16, 8), (8, 4).
- **Boosted launch (A2):** (8, 8) with both kinks at v = 0.3.
- **Pairing:** every run is paired with κ = 0 on the same candidates and marks. All combinations of one seed use the
  same candidates and marks.
- **Continuum reference:** the coupled equations by finite differences, with the same initial data (boosted for A2).
- **256 new sprinklings**, seed label "F5.2b".

## 3. Readout (frozen)

- **Kink centre:** the kink fit of F3.14a per slab of width 1, for each field.
- **Displacements:**
  - d_B = [x_B(19.5) − x_B(0.5)]_κ − [the same]_{κ=0}; towards A is negative. d_A likewise.
  - For A2, the change of the separation x_B − x_A.
- **Ratios:** r = mean(d)/d_cont, a ratio of means (lesson of F5.2a3). The SE is by bootstrap over sprinklings.
- **Universality:** u = r_B(M_B 4) − r_B(M_B 8) at M_A = 8, paired by seed.
- **Density:** realised fraction / mean Π² per bin of x, for t > 2T/3.
- **Intervals:** 60 realised timelike pairs near A on 8 sprinklings. The count inside / ρ·(proper area ∫Π²) and
  / ρ·(coordinate area).

## 4. Criteria (proposed; 2 SE; every condition on its own; nothing averaged; nothing chosen after the run)

| level | condition | tolerance | development value (32 sprinklings) |
|---|---|---|---|
| instrument A0 | a change of later marks leaves earlier events unchanged (4 sprinklings) | exact | yes |
| instrument A1 | linear extension t + 0.5x: same realised set and fields | 1e-10 | 3.7e-13 |
| instrument A3a | \|realised/Π² − 1\| + 2 SE | < 0.01 | 1.000 ± 0.002 |
| instrument A3b | \|count/proper area − 1\| + 2 SE | < 0.03 | 1.004–1.006 |
| instrument, null | \|mean B displacement at κ = 0\| + 2 SE, every combination | < 0.1 | (wandering only) |
| **level 1, growth** | count/coordinate area + 2 SE | < 0.9 | 0.63 |
| **level 2, active** | \|r_B − 1\| + 2 SE for M_A = 4, 8, 16 | < 0.15 | 1.051, 1.029, 1.027 (SE 0.036 / 0.023 / 0.014) |
| **level 3, universality** | \|u\| + 2 SE | < 0.06 | −0.012 |
| **A2, covariance** | \|r_sep − 1\| + 2 SE (boosted launch) | < 0.15 | (smoke only) |

**FAILED** if any level 2, level 3 or A2 condition has |·| − 2 SE above its tolerance.
**PASSED** if all conditions of levels 1–3 and A2 hold. **UNDECIDED** otherwise.

**Why these numbers (the rule of 10 Nov: expected 2 SE + the largest development offset + a margin).**
- **r:** at 256 sprinklings 2 SE ≈ 0.01–0.026, and the development offset is up to 0.05 (the unexplained excess of
  2–5 %), so 0.15.
- **u:** 2 SE ≈ 0.02 and the offset 0.012, so 0.06.

**Reported:**
- r_A (A towards B) and the third law M_A d_A / (M_B |d_B|);
- L along x against the continuum;
- the realised fraction.

**Cost:** about 2 h with 6 processes.

## 5. Reading declared now

**PASSED.** A causal growth rule uses only the realised past, the counts and independent marks. Under it, matter thins
its future order (slower counted clocks), and a test defect falls towards the source in proportion to the source's
energy, independently of its own mass, as box_g L = −κq predicts, also for a moving system. With F3.14b and F5.2a3, the
inertial, passive and active mass come from one field and one rule (in 1+1). The rule is a candidate for the Axiom of
the dynamics of the order.

**FAILED.** The discrete growth does not reproduce the continuum response, or it breaks universality.

**Levels (the author's).**
1. Growth: the source changes the probability of realising future events.
2. Active response: the change of the order moves a test defect.
3. Universal gravity: the active mass equals the energy of the source, and the acceleration does not depend on the test
   mass. Only level 3 closes the test of active mass.

**Limits:**
- 1+1: scalar response R2, no 1/r, and the potential of a static source grows in space and time;
- the causal relation does not change (conformal flatness), as it will have to in 2+1;
- κ is chosen (G is not derived);
- the potential is chosen;
- the excess of 2–5 % over the continuum is unexplained.
