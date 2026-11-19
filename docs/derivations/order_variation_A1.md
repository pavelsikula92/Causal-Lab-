# A1: admissible variations of the order and the variational rule (specification)

19 Nov 2026. Status: a specification, no runs, decides nothing. It fixes what has to be decided before any
variational computation on the order. Nothing here is frozen; the open choices are listed in section 12 and are the
author's.

**Direction (the author, 19 Nov).**
- The question is not a field equation for Π. It is the rule by which the presence of matter selects, among the
  possible causal orders, the order that is admissible and stable.
- Under A2-B the only geometric degree of freedom is the order H; Π = 𝒫[H] is a derived readout with Π_∞ = 1. A
  matter-induced change of Π must arise through a change of H: M → H_M → 𝒫[H_M] = Π_M.
- The work is split into three steps:
  - **A1** (this document): which changes of H are admissible, and the variational rule of S_BD on them, without matter;
  - **A2**: the action of matter S_m[H, F], after the flow F_v is defined operationally (a separate branch F);
  - **A3**: finite-amplitude dynamics with a source, its stability and the derived Π_M.

**Naming.** In this document A1–A3 are the three steps above. They are not the open decisions A1 (Δτ_v) and A2 (I_v)
of docs/derivations/variational_audit.md, not the working position A2-B (substrate_and_action_audit.md §7), and not
part A of the roadmap, of which they are a sub-programme. A renaming is the author's choice.

## 1. Scope

**A1 defines:**
- the configuration space of the order for a variational principle (section 2);
- two families of admissible changes δH, kept separate (section 3);
- the boundary conditions of a variation (section 4);
- the level of the action, a single order or the ensemble mean (section 5);
- the role and the fixing of the nonlocality ε (section 6);
- stationarity and the second variation in each family (section 7);
- three controls that any implementation of A1 must pass before A2 and A3 use it (sections 8–10);
- their verdict rules (section 11).

**A1 does not:**
- couple matter, and does not define the flow F_v (A2 and branch F);
- solve for a sourced order or derive Π_M, 1/r, G or C.1 (A3);
- evaluate Π. Π = 𝒫[H] is not part of the configuration and appears in A1 nowhere;
- treat the stochastic route, in which the action weights histories or transitions. That route includes classical
  sequential growth (Rideout and Sorkin 1999), which rests on its own assumptions and is not derived from S_BD. It is
  recorded as an alternative and is not mixed with the variational route here;
- revise any verdict. F6.8f stays FAILED, F6.8d undecided, F6.4a FAILED, F6.4c UNDECIDED (R132).

## 2. Configuration space

1. **The order.** H = (C, ≺) with C a finite set and ≺ irreflexive and transitive, hence acyclic. Local finiteness
   (every interval finite) holds automatically for finite C.
2. **Labels.** The action and every criterion are functions of the order up to isomorphism (discrete general
   covariance). Embedding coordinates, where a family uses them, are bookkeeping, not degrees of freedom.
3. **Domain.** H is the order of a region Ω, split into an interior Ω° and a collar Ω_c (section 4). Variations act
   in Ω° only.
4. **Density.** The density of elements per unit spacetime volume, ρ, is fixed and fixes the units. In family M it is
   the Poisson intensity; in family I the change of the number of elements is a recorded part of each move.
5. **Admissibility.** A configuration is admissible if (i) it is a partial order and (ii) it lies in the manifold-like
   class of section 10. Condition (ii) holds by construction in family M. In family I it is a constraint that is
   checked, never assumed.

## 3. Variation families

The two families answer different questions and must not merge. For each, the table states what it changes and what
it keeps. Results are recorded per family. A disagreement between them is a result and is not averaged.

| | family M (metric-parametrised) | family I (intrinsic moves) |
|---|---|---|
| question | is the reference order stationary and stable against geometric deformations? | against arbitrary changes of the order? |
| variable | δg on Ω°, in a fixed finite basis | elementary moves of the order in Ω° |
| changes | cones (relation) and density (conformal factor) | relations; number of elements |
| keeps | partial order; manifold-likeness; Lorentz invariance in distribution | partial order |
| imports | a continuum metric and chart | nothing beyond the order |
| role | reference level | the check that the result does not depend on the continuum family |

### 3.1 Family M

- **Definition.** H_g is a sprinkling of (Ω, g) at density ρ, and a variation is g → g + h with h supported in Ω°.
- **Two realisations, recorded separately:**
  - **M-a:** a realisation with the sprinklings coupled across g, by a common Poisson seed and thinning. K.17 measured
    how little the coupling helps. In the laboratory box the correlation of S₀ and S₊ was 0.55–0.85, because thinning
    changes the elements, and the noise of the second difference was 10–10⁴ times the limit of feasibility.
  - **M-b:** the expectation over sprinklings by quadrature, with no realisation at all. This is the F6.8 route,
    calibrated by K.4.
- **Basis of h.** The static sectors of F6.8 (conformal ω, mixed, optical ν). A list of localised bumps or plane-wave
  packets with compact support in Ω° is fixed before any computation. Pure-gauge directions h = L_ξ g are identified
  and kept as a null check (section 7).
- **What this family can and cannot show.** It inserts the continuum. Stationarity of ⟨S_BD⟩ under δg is, in the local
  limit, the vacuum Einstein–Hilbert condition plus nonlocal and age terms. Family M therefore tests whether a candidate
  order is a stationary point of the discrete action. It cannot test whether the order selects geometry.

### 3.2 Family I

- **Elementary moves:**
  - **I-1, remove a link.** A link is a covering pair u ⋖ v. Removing it keeps transitivity, because no element lies
    between u and v, so no pair of relations implies u ≺ v.
  - **I-2, add a relation u ≺ v** with v ⊀ u, followed by transitive closure: every a ⪯ u becomes related to every
    b ⪰ v. The closure can add many relations, and their number is recorded with the move.
  - **I-3, insert an element** with a past set P (a down-set) and a future set F (an up-set), such that every element
    of P precedes every element of F, followed by closure.
  - **I-4, remove an element.** The relations among the other elements are kept; transitivity holds because the
    relations through the removed element are already present.
- **Locality of a move.** A move is local if all elements it touches lie in one interval of the reference order with
  at most n_loc elements. n_loc is fixed in advance and ties to ε (section 6).
- **Inverse pairs.** I-1/I-2 restricted to links and I-3/I-4 are inverse classes. Section 7 uses them.
- **Known behaviour in the laboratory.** F6.4a reproduced, in 2D, the manifold-like phase at β = 0.1 and the
  transition at β ≈ 2.4 of Surya (2012). The setup was orders of N = 50 weighted by exp(−βS) with ε = 0.12; F6.4c
  confirmed the instrument. At large β the ensemble is crystalline, not manifold-like (published height 4.18; the
  laboratory measured 3.70–3.89, F6.4a FAILED and F6.4c UNDECIDED on the height). Driving S_BD down over unrestricted
  orders therefore leaves the manifold-like class in 2D. This is the reason for section 10.

## 4. Boundary conditions

1. **Compact support.** Every variation is supported in Ω°. The collar Ω_c is held fixed, which is a Dirichlet
   condition on the order.
2. **Width of the collar.** It must contain every element that can carry non-negligible kernel weight into an interval
   that reaches Ω°. With the smeared kernel this is a width of several ℓ_k. It is fixed before a run and checked by
   doubling it (section 11, instrument).
3. **Boundary terms.** For a flat causal interval the mean BDG action tends to the volume of the codimension-2 joint of
   the two light-cone boundaries; the bulk contribution vanishes (Buck, Dowker, Jubb, Surya, arXiv:1502.05388; worked
   out for d = 2, 3, 4). Contributions of other boundaries and corners are studied separately (arXiv:2007.13206,
   2501.00139; not re-checked here).
   - With compact support and a sufficient collar, such terms do not change under a variation. Their size is reported
     nevertheless, because the controls of sections 8 and 9 compare small differences.
4. **Past boundary and age.** Ω has a past boundary at age T, as in F6.8a–f. T is an explicit parameter of every A1
   computation and every result is reported as a function of T. The age terms are the reason: conformal linear in kT,
   mixed (kT)³ (F6.8d undecided), optical (kT)⁵ (F6.8e3), and the expansion fails at ε_b² = b²(kT)⁴ ≳ 1 (F6.8f, decision
   36).
5. **Future boundary.** A variational principle on a slab [0, T] fixes the future collar as well. A retarded growth
   dynamics would not. This difference between the two routes is recorded and not resolved in A1.
6. **Normalisation Π_∞ = 1.** Not used in A1. It belongs to A3, where Π is read.

## 5. Action level

- **S_ε[H], a single order.** The smeared Benincasa–Dowker(–Glaser) action in d = 4 with nonlocality ε and
  ℓ_k = (ερ)^(−1/4). Sign and normalisation are those of F6.8, where K0 = 1 corresponds to the Einstein–Hilbert value.
  It is defined in both families.
- **⟨S_ε⟩[g], the ensemble mean.** The expectation over sprinklings of g. It exists in family M only and is computed
  by quadrature (M-b).
- **They are different principles.** Stationarity of ⟨S⟩ means that the order is a typical sprinkling of a stationary
  g: a mean-field statement. Stationarity of S is a statement about one order. In the laboratory box (K.17) the
  fluctuation of a second difference of S was 28–1045 against an expected 2.73 for ε = 0.03–0.3, so on orders of that
  size a stationarity of S under a small δg is not measurable and may not be meaningful.
- **Options (section 12):**
  - **L1:** the principle on ⟨S⟩ (family M only), which imports the continuum;
  - **L2:** the principle on S (family I, and M-a), with a statistical stationarity criterion (section 7);
  - **L3:** both, with L1 as the reference and L2 as the check, and the results recorded separately.

## 6. The nonlocality parameter ε

- **Meaning.** The smearing scale of the kernel. ℓ_k = (ερ)^(−1/4) is the length below which the action is not local.
- **What ε controls:**
  - **The fluctuations of S.** K.17: the sd of the second difference fell from 1045 to 28 between ε = 0.3 and 0.03.
  - **Where the local Einstein–Hilbert form fails.** F6.8a: K_loc changes sign below a wavelength of about 13 ℓ_k.
  - **The size of the age terms.** For example, c_T k³ ℓ_k² T in the conformal sector.
- **Fixing ε.** ε is a parameter of the action, not of a result.
  - It is fixed before any A1 computation, and the same value serves all three controls.
  - Every result is also reported at two more values, to show its dependence.
  - No value is chosen after seeing which one gives stationarity.
  - The physical value of ℓ_k is open (the paper's window) and is not decided here.
- **Range in family I.** On small orders, ℓ_k must be well above the discreteness and well below Ω. This limits the
  usable N, and through it the dimension (section 10).

## 7. Stationarity and the second variation

1. **Stationary point, not a minimum.**
   - The static Einstein–Hilbert form in the F6.8 variables, 6ω² + 8ων + 2ν², is indefinite: the matrix [[6, 4],
     [4, 2]] has determinant −4.
   - A1 therefore defines stationarity, and judges the second variation by its sign structure against the
     Einstein–Hilbert form, never by positivity.
   - "Stable" in the sense of A3 (dynamical stability of a solution) is a different notion and is not defined here.
2. **Family M (L1).**
   - **First variation:** δ⟨S⟩[g; h] = 0 for every h of the fixed basis, within the quadrature error calibrated by K.4.
   - **Second variation:** the matrix δ²⟨S⟩[h_i, h_j] on the basis.
   - **Gauge:** along pure-gauge h = L_ξ g the first and second variations must vanish (an instrument null). This holds
     exactly, because ⟨S⟩ is defined on the order.
3. **Family I (L2).** S is a function on a discrete set, so stationarity has to be defined. Three candidates, and the
   choice is the author's (section 12):
   - **S-a, local extremum:** every local move raises S, or every local move lowers it. This is expected to fail for
     an indefinite action, and a failure says nothing about geometry.
   - **S-b, balance:** for each inverse pair of move classes, the mean ΔS over the class and the mean over its inverse
     sum to zero within errors, and ΔS has no systematic sign within each class at the configuration.
   - **S-c, coarse-grained:** moves are grouped into composite deformations that approximate the basis of family M,
     and ΔS of these composites is treated as in family M. This bridges the two families, but by construction it
     imports family M's continuum.
4. **Discrete second variation.** For commuting local moves m₁, m₂:
   Δ²S(m₁, m₂) = S(m₁m₂) − S(m₁) − S(m₂) + S.
   Its matrix on a fixed set of move pairs is the family-I counterpart of the Hessian.

## 8. Control 1: the flat vacuum

- **Claim tested.** Without matter, the order of a flat slab Ω is a stationary point in each family and at each level.
- **Family M, L1.**
  - By Lorentz invariance, the only first-order response of ⟨S⟩ around flat space in the interior is δ⟨S⟩ = c ∫ tr h.
  - Flat space is stationary if and only if this cosmological coefficient c vanishes. In F6.8a, β = 0 was checked
    by the null k = 0 (J0/J ~ 1e−5).
  - The control repeats this with compact support and the A1 boundary conditions, as a function of the age T and of
    the distance from the past boundary.
  - This is a real check: a c ≠ 0 near the past boundary would be an age term at first order.
- **Family I, L2.** Whether a flat sprinkling is stationary under intrinsic moves in the chosen sense (S-a/b/c).
  - This control can fail legitimately.
  - In 2D the crystalline orders have a lower action than the flat ones (section 3.2).
  - A failure is a result about family I, not about family M.

## 9. Control 2: consistency with F6.8

- **Family M, L1.** The second variation around flat at the A1 boundary conditions must reproduce:
  - in the local regime, the static Einstein–Hilbert coefficients 6, 8, 2 of F6.8a, F6.8b2 and F6.8c2;
  - at the same parameters, the known age terms (c_T = 0.1035 conformal; the optical (kT)⁵ with c5 = 3.05e−7).
- **Why only a check.** F6.8 used phase-averaged plane waves with a sharp past cut, while A1 uses compact support with
  a fixed collar.
  - Agreement is expected only where the two boundary conditions coincide in the limit (kT small, support far from the
    collar), and that mapping is specified before the run.
  - This tests the A1 implementation. It is not a new result.
- **Family I.** No requirement. A coarse-grained (S-c) reproduction of the coefficients would be welcome, but K.17
  makes it infeasible in a laboratory box at present, and it is not demanded.

## 10. Control 3: escape from the manifold-like class

**Purpose.** A configuration that a procedure reports as stationary or as having a lower action, but that is not
manifold-like, can never count as a gravitational solution, whatever its action. This is a verdict of its own,
reported next to every result of family I (and of any optimiser), never merged with the others.

**Instruments available in the laboratory:**
- **K.1 (E38-B reading):** interval dimensions, homogeneity, BD action and paired null models. It is calibrated in
  3+1, where it reads a sprinkling as GEOMETRIC with d* = 4 from N* = 32 000 (replicated). Below N* it has no
  calibration, and a 3+1 family-I computation on smaller orders would need its own calibration first.
- **Measures of F6.4a/c:** the Myrheim–Meyer dimension d_MM, the height h and the fraction of relations χ, calibrated
  in 2D at N = 50.

**Criteria (to be frozen before any run).** A configuration is manifold-like if all of the following hold:
1. K.1 reads it as GEOMETRIC with d* = 4 (3+1, N ≥ N*), or the 2D counterpart at the 2D calibration;
2. d_MM is within the null band of flat sprinklings of the same N in each of a fixed set of sub-intervals;
3. the abundance of intervals N_n (n ≤ n_max) is within the null band;
4. the height and the longest chain are within the null band.

Each null band is set from flat sprinklings before the run, at a stated quantile.

**Verdict:** NO_ESCAPE, ESCAPE or UNDECIDED. An ESCAPE voids every stationarity or lower-action claim of that
configuration. It does not turn the control into FAILED.

**Recommendation.** Run family I in 2D first. The instruments and the published phase structure exist there (F6.4a/c),
and the 3+1 instrument needs N ≥ 32 000.

## 11. Verdict rules

**Separate cards.** Each of the controls of sections 8, 9 and 10 becomes a separate card in each family. A1 as such
has no single verdict.

**Verdict of a control card** (lab standard, |x| + 2 SE against the tolerance):

| verdict | when |
|---|---|
| PASSED | the criterion is met and the instruments are OK |
| FAILED | the criterion is violated beyond the tolerance and the instruments are OK |
| INSTRUMENT_FAILED | an instrument check fails |
| UNDECIDED | the errors do not allow a decision |

Instrument checks for INSTRUMENT_FAILED:
- quadrature against K.4;
- the pure-gauge null of section 7;
- the flat-sprinkling null of the manifold-like criteria;
- invariance under doubling the collar;
- the ε and T dependence reported as required in sections 4 and 6.

**Phase verdict (section 10).** Separate. An ESCAPE excludes the configuration from any claim.

**Power of the controls.** Every null or control is sized before the run so that 2 SE is at most a third of its
tolerance. Its false-alarm rate (5 % for a 2 SE rule) is stated on the card. Two earlier runs failed on
underpowered controls: F5.2c seed 1 and the C2 diagnostic cl14.

**Other rules:**
- **Families.** Verdicts per family are separate, and a disagreement between the families is recorded as a result.
- **R132.** No criterion changes after the run. No earlier verdict (F6.8d, F6.8f, F6.4a, F6.4c) is reread in the
  light of A1.

## 12. Open decisions (the author's)

1. **The level:** L1, L2 or L3 (section 5).
2. **Stationarity in family I:** S-a, S-b or S-c (section 7.3).
3. **ε:** the reference value of ε (or ℓ_k in units of ρ^(−1/4)) and the two comparison values (section 6).
4. **Ω:** its shape, size and age T, and the width of the collar (section 4).
5. **Family I moves:** the move set and n_loc (section 3.2).
6. **Order of dimensions:** family I in 2D first, then 3+1 (section 10).
7. **The basis of family M:** the fixed list of h, and the mapping to the F6.8 boundary conditions (section 9).
8. **Renaming:** whether A1–A3 are renamed to avoid the clash with variational_audit.md (section "Naming").

## References

- Benincasa and Dowker, Phys. Rev. Lett. 104, 181301 (2010): the BD action. Laboratory use: rule_from_action.md §5,
  §8; F6.8.
- Surya, Class. Quantum Grav. 29, 132001 (2012): phases of 2D orders with the BD action (calibrated in F6.4a/c).
- Buck, Dowker, Jubb and Surya, arXiv:1502.05388: boundary terms of the mean action of a flat interval.
- Rideout and Sorkin, Phys. Rev. D 61, 024002 (2000): classical sequential growth (the stochastic alternative, not
  used).
- Laboratory: docs/derivations/substrate_and_action_audit.md §7–8 (A2-B, B6-full); variational_audit.md;
  f68f_dev.md; findings F6.8a, F6.8b2, F6.8c2, F6.8e3, F6.8f, F6.4a, F6.4c, K.1, K.4, K.17.
