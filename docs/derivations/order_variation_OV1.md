# OV1: admissible variations of the order and the variational rule (specification)

19 Nov 2026. Status: a specification with the author's decisions of 19 Nov (10:42), plus a development script for 2D
that decides nothing (section 13). No criterion is frozen. A card follows only after the development numbers, by a
decision of the author.

**Direction (the author, 19 Nov).**
- The question is not a field equation for Π. It is the rule by which the presence of matter selects, among the
  possible causal orders, the order that is admissible and stable.
- Under A2-B the only geometric degree of freedom is the order H; Π = 𝒫[H] is a derived readout with Π_∞ = 1. A
  matter-induced change of Π must arise through a change of H: M → H_M → 𝒫[H_M] = Π_M.
- **Steps:**
  - **OV1** (this document): the admissible changes of H and the variational rule of S_BD on them, without matter;
  - **OV2**: the operational definition of the flow F_v (an independent branch), then the action of matter S_m[H, F];
  - **OV3**: finite-amplitude dynamics with a source, its stability, and the derived Π_M.

**Naming.** OV1–OV3 replace the labels A1–A3 of the first draft (file order_variation_A1.md, renamed). The open
decisions A1 (Δτ_v) and A2 (I_v) of docs/derivations/variational_audit.md and the working position A2-B are unchanged
and unrelated.

**Hierarchy (the author).** The action of a single order, S[H], defines the microscopic object. The ensemble mean
⟨S⟩_μ = ∫ S[H] dμ(H) defines its statistical response. The two levels are never merged into one result.

## 1. Scope

**OV1 defines:**
- the configuration space of the order (section 2);
- two families of admissible changes, M and I (section 3);
- the region, the age and the collar (section 4);
- the two levels of the action (section 5);
- ε (section 6);
- stationarity (section 7);
- three controls: the flat vacuum, which blocks; the consistency with F6.8; and the phase (sections 8–10);
- the verdict rules (section 11);
- the order of dimensions (section 12).

**OV1 does not:**
- couple matter or define F_v (OV2);
- solve for a sourced order, or derive Π_M, 1/r, G or C.1 (OV3);
- evaluate Π;
- treat the stochastic route (sequential growth; Rideout and Sorkin 2000). That route rests on its own assumptions
  and is not derived from S_BD;
- reread any verdict (R132). F6.8f stays FAILED, F6.8d undecided, F6.4a FAILED and F6.4c UNDECIDED.

## 2. Configuration space

1. H = (C, ≺): C is a finite set, and ≺ is irreflexive and transitive, hence acyclic. Every interval is finite.
2. The action and every criterion depend on the order only up to isomorphism. Relabellings and isomorphic
   representations are not changes of the order.
3. H is the order of a region Ω (section 4), with a fixed collar Ω_c and an interior Ω°.
4. The density ρ fixes the units. In family I the number of elements is fixed (the author: "the number of events and
   the fixed boundary events are kept"). The flat-vacuum ensemble uses Poisson sprinklings, so N varies between
   realisations.
5. A configuration is admissible if it is a partial order in the manifold-like class (section 10). In family M this
   holds by construction; in family I it is checked.

## 3. Variation families

The two families answer different questions. Each states what it changes and what it keeps. They are recorded
separately, and a disagreement between them is a result.

| | family M | family I |
|---|---|---|
| variable | δg in the frozen basis of F6.8 | reversible local rewrites of the order |
| changes | cones and density | relations |
| keeps | partial order; manifold-likeness; Lorentz invariance in distribution | partial order; N; collar |
| imports | a continuum metric and chart | nothing beyond the order |

### 3.1 Family M (the author: the frozen basis of F6.8)

- **Basis.** The basis of F6.8 exactly:
  - the static sectors ω and ν;
  - the same normalisations and the same gauge;
  - the same plane waves k, phase averaging and past cut;
  - every sector F6.8 judged, the optical one included.

  No mode is dropped after a result has been seen.
- **The whole quadratic form is checked**, not three numbers along the axes:
  Q_EH(ω, ν) = 6ω² + 8ων + 2ν², that is the matrix [[6, 4], [4, 2]] with determinant −4. This covers the diagonal
  terms, the mixed term and the sign structure.
- **Two realisations, recorded separately:**
  - **M-a:** sprinklings coupled across g. K.17 showed that the coupling cancels little noise: the correlation
    between sprinklings was 0.55–0.85, and the noise of the second difference was 10–10⁴ times the limit of
    feasibility.
  - **M-b:** the expectation by quadrature, the F6.8 route calibrated by K.4.
- **Family M is a 3+1 matter.** In 2D every metric is locally conformally flat, and a conformal factor does not change
  the order, only the density. The bulk Einstein–Hilbert action in 2D is moreover topological (Gauss–Bonnet), so its
  second variation in the bulk is zero. A 2D family M would therefore test only nonlocal and boundary terms, and it is
  not a reproduction of the 3+1 stiffness of F6.8. The direct comparison with F6.8 belongs to 3+1 (section 9).

### 3.2 Family I

**The author's rule.** All reversible local rewrites that keep the partial order, with these properties:
- the number of elements is kept, and so are the boundary elements;
- the support lies inside a convex interval;
- acyclicity and transitivity hold after each rewrite;
- the inverse rewrite is included whenever it exists;
- pure relabellings are excluded.

The set is generated by a complete rule, not picked by hand.

**Realisation in 2D (proposed here, for the author to confirm).** The class is the 2D orders of F6.4a/c: an order is
the intersection of two total orders U and V on the N elements, a sprinkling of a diamond in light-cone coordinates.
- **Elementary rewrite:** the transposition of two elements adjacent in U, or adjacent in V.
- **What it does.** It changes exactly one relation, between two elements whose interval is empty: a link is created
  or removed. It keeps N and the class of 2D orders, and it is an involution, so it is its own inverse.
- **Completeness.** Adjacent transpositions generate every permutation of U and V, hence all 2D orders on N elements.
  A rewrite whose relation change has a non-empty interval is a product of elementary ones.
- **Where n_loc enters.** The elementary rewrite needs no range: its relation change is always a link. n_loc (section
  4) sets the collar and the region from which a change of the action is collected; it is derived from the kernel, as
  decided.
- **What is left out.** Rewrites that leave the class of 2D orders, that is general posets of higher order dimension,
  are not included. In that larger space the phase instruments of F6.4 are not calibrated. Including them is a later
  step, to be decided.

**The symmetric pair (S-b).** For an element a and the order k ∈ {U, V}:
- H₊ swaps a with its k-successor and H₋ swaps a with its k-predecessor, a rank shift of a by ±1;
- the first and the last element of each order have only a one-sided move. They are recorded as one-sided
  constraints and not treated as symmetric directions (two per order and level).

## 4. Region, age and collar

1. **Ω.** A symmetric Alexandrov diamond Ω_T = J⁺(p_T) ∩ J⁻(q). The age is varied by moving the past tip p_T at a
   fixed future tip q, which gives nested diamonds; inside them the same interior readout region can be kept.
2. **Ages.** T/w_∂ ∈ {8, 16, 32}.
3. **w_∂ and n_loc from the kernel (the author).**
   - n_loc(ε) is the smallest n such that Σ_{m≤n} |f(m, ε)| ≥ 0.99 Σ_m |f(m, ε)|.
   - w_∂ is the proper time of an interval with n_loc expected elements. In 2D, with the units of section 13,
     w_∂² N_exp = n_loc, so a diamond of age T holds N_exp = (T/w_∂)² n_loc elements.
   - The same definition holds for every ε, and neither is changed after a result (section 11).

   | ε | n_loc | N_exp at T/w = 8, 16, 32 |
   |---|---|---|
   | 0.06 | 135 | 8 640, 34 560, 138 240 |
   | 0.12 | 65 | 4 160, 16 640, 66 560 |
   | 0.24 | 30 | 1 920, 7 680, 30 720 |

4. **Boundary terms.** No term proportional to the volume of the boundary is added. For a flat interval the mean BDG
   action tends to the volume of the codimension-2 joint of the two light-cone boundaries (Buck, Dowker, Jubb, Surya,
   arXiv:1502.05388). In 2D that joint is the two side corners of the diamond. w_∂ only separates the collar from the
   interior for the diagnostics; it does not replace the boundary term.
5. **Collar.** In family I the elements of Ω_c are fixed. For diagnostics, the interior of a 2D diamond is the set of
   points whose own diamond of proper time w_∂ lies inside Ω: in light-cone coordinates, u and v at least w_∂ from
   each edge.
6. **Normalisation Π_∞ = 1.** Not used in OV1.

## 5. Levels of the action (L3)

- **Individual action.** S: H ↦ S[H] is the smeared BD(G) action of one order, in the F6.4a convention in 2D and in the
  F6.8 convention in 3+1.
- **Ensemble.** ⟨S⟩_μ = ∫ S[H] dμ(H). The measure μ is specified and frozen for every variation. In OV1, μ is the
  Poisson sprinkling of the region at density ρ. If a variation changed μ, the change of μ would be counted separately
  (the author's condition).
- **Verdicts, separate:**
  - **L1:** stationarity of the ensemble mean, the vanishing of the cosmological coefficient c(T), and, in 3+1, the
    reproduction of Q_EH.
  - **L2:** the distribution of the individual realisations: the spread and any systematic direction of the
    individual centred differences.
  - **L3 needs both.** It does **not** require each order to be exactly stationary, since that would read Poisson noise
    as a failure of the variational principle.

## 6. ε

**Decision (the author, 19 Nov, 11:30).** For the initial 2D OV1: ε₀ = 0.12, ε₋ = 0.06, ε₊ = 0.24.
- 0.12 is the 2D reference value of F6.4a/c (Surya 2012). It is not a derived constant of the theory.
- The factor two on each side tests whether stationarity and the phase behaviour depend on the coefficient.
- The earlier "ε from F6.8" is withdrawn. F6.8 works in the units ερ = 1 (ℓ_k = 1), not with a dimensionless
  per-element coefficient, and the two are not interchangeable.

**The role of ε, defined (the author's condition).** ε is the dimensionless smearing coefficient of the kernel of the
action, per element.
- **Where it enters.** It enters the functional only:
  - in S = 4ε (N − 2ε Σ_n N_n f(n, ε)) in 2D, through the prefactor and the weights f(n, ε) of the interval
    abundances;
  - in the definitions derived from that kernel: n_loc(ε) and w_∂(ε) (section 4), hence the sizes N_exp = (T/w)²
    n_loc and the collar.
- **The local move does not depend on ε.** An adjacent transposition is defined by the order alone; ε enters only
  ΔS of the move and, in the phase sampler, the acceptance through S.
- **ε is not a density, a length or a sampling parameter.** The density sets the units: N_exp is the expected number
  of elements. The length ℓ_k = (ερ)^(−1/d) is derived from ε and ρ and is not an independent input. The sampling
  parameters (seeds, sweeps, β) are separate.
- **No automatic transfer.** The value is not carried over to 2+1 or 3+1 until the same normalisation is verified
  there.

## 7. Stationarity (S-b) and the second variation

1. **Stationarity, not a minimum.**
   - The primary definition is S-b: δS = 0, with the signature of δ²S equal to the signature of δ²S_EH in the
     corresponding family.
   - S-a (a minimum) is rejected as a criterion.
   - S-c (coarse graining) is a secondary diagnostic only. It never rescues a failed S-b.
2. **Discrete differences (the author).** For a direction δ with both perturbations H₊ and H₋ admissible:
   - D_δ S = (S[H₊] − S[H₋])/2;
   - D²_δ S = S[H₊] + S[H₋] − 2S[H].

   A direction without an admissible inverse is a one-sided constraint, recorded as such.
3. **Family M (L1).** δ⟨S⟩[h] = 0 on the basis, and the matrix δ²⟨S⟩[h_i, h_j] has the signature of Q_EH. Along pure
   gauge both must vanish; this is an instrument null.
4. **Family I (L2 within L3).** Over the interior elements, the distribution of D S has zero mean and no systematic
   drift with position or age. D² S is reported.
   - In 2D the bulk Einstein–Hilbert action is topological, so there is no continuum signature to compare with. In 2D
     the signature item therefore reduces to reporting D² S, and it decides nothing.

## 8. Control 1: the flat vacuum (blocking)

- **Statement.** Lorentz invariance allows the linear response δ⟨S⟩ = c(T) ∫_Ω tr h dV. The requirement is
  **c(T) = 0** at every pre-registered age, including the contribution of the past boundary.
- **Readout in 2D.**
  - c(T) is the coefficient of the volume N in ⟨S⟩ over the nested diamonds T/w = 8, 16, 32, read from the slope of
    S against N fitted on each sprinkling.
  - The paired differences S(16) − S(8) and S(32) − S(16) are reported. With the volumes growing by a factor 4 per
    step, a volume term would show as differences growing by a factor 4.
  - The continuum expectation is ⟨S⟩ → a constant, the joint term (section 4.4).
- **L2 part.** The individual orders fluctuate. The distribution of D S over interior elements must show no nonzero
  systematic drift.
- **If c(T) ≠ 0 because of the past boundary**, it is not subtracted as an artefact. It is a result: the construction
  of the action and the boundary fails the vacuum control, and the origin of the term is explained first.
- **Blocking.** OV2 and OV3 do not use a construction that fails this control.

## 9. Control 2: consistency with F6.8 (3+1)

- **Family M-b.** Around flat space, with the frozen basis of F6.8 and its boundary conditions, the method must
  reproduce:
  - the whole form Q_EH, namely the coefficients 6, 8 and 2 of F6.8a, F6.8b2 and F6.8c2 and the signature;
  - at the same parameters, the age terms: c_T = 0.1035, and the optical (kT)⁵ with c5 = 3.05e−7.
- **Purpose.** This tests the implementation, not new physics.
- **Not part of this control:** a 2D result is never presented as this reproduction.

## 10. Control 3: the phase

**Two separate outputs (the author):**
- **A. Stationarity inside the manifold-like phase:** whether the flat order is stationary (sections 7 and 8).
- **B. Stability of the phase:** whether changes of the order lead to one of the non-manifold phases known from
  F6.4a/c.

**Diagnostics.** The F6.4 diagnostics, with no new criterion: the ordering fraction χ, d_MM from χ, the height h (links
of the longest chain) and the interval abundances N_n.

**Classification.** Each run reports one of four outcomes:
- stationary inside the manifold-like phase;
- stationary outside it;
- unstable, or an escape from it;
- undecidable by the instrument.

**Bands at the actual N (the author, 19 Nov, 11:30).**
- The reference band of the manifold-like phase is set from flat sprinklings at the same N, the same geometry of the
  region and the same boundary conditions as the later test. It is never extrapolated from the crystalline reference
  at N = 50.
- For each N there is a separate calibration set of flat sprinklings, with seeds disjoint from the validation seeds
  of OV1:
  - calibration: seeds ≥ 100 000;
  - validation: seeds < 100 000.

  diag/ov1_2d_summ.py takes the bands from the calibration seeds and reports the fraction of validation sprinklings
  outside them.
- **Limit of the reading.** The flat band fixes only the manifold-like reference. A departure from it does not by
  itself say whether the structure is crystalline or another non-manifold structure. That classification stays with
  the F6.4 diagnostics, which are reported separately: χ, d_MM, h, N₀ and the layered profile.

**The escape: an F6.4-type sampler with local moves (the author, 19 Nov, 11:30), not OV3.** OV3 is not yet defined,
and using it here would make an unfinished theory its own control instrument. The sampler (diag/ov1_2d_mcmc.py):
- proposals: the adjacent transpositions of family I (U or V with probability 1/2, a position uniformly);
- acceptance and sampling: those of F6.4a, unchanged:
  - Metropolis with exp(−βS);
  - one sweep = N(N − 1)/2 proposals;
  - one sample per sweep after the thermalisation;
  - per chain: S, var S, χ, d_MM, h, N₀, the acceptance and τ_int;
- the proposal is symmetric and the moves connect all 2D orders, so the stationary distribution is the Boltzmann
  distribution of F6.4a. Only the mixing differs;
- the result is labelled an exploration of the configuration space, not a physical time evolution.

**P0, the instrument.** Does the local sampler reproduce the laboratory's F6.4a measurement at N = 50, ε = 0.12 on the
F6.4a β grid?
- Checks:
  - the exact Boltzmann distribution at N = 4 (TV 0.0079 ≤ 0.01);
  - β = 0 against 20 000 directly sampled orders (|ΔS| 0.016, |Δχ| 0.0019, |Δh| 0.000);
  - the local ΔS accumulated over 10 000 proposals against a full recomputation (4e−13).
- **Results** (N = 50, ε = 0.12, 4 chains × (500 + 3000) sweeps; diag/out/ov1_2d_mcmc):
  - **Phase I, β = 0.1.** S 3.828 ± 0.011, χ 0.4980, h 10.230. The laboratory's F6.4a measured 3.835, 0.4981 and
    10.237; Surya published 3.846, 0.499 and 10.232.
  - **Transition.** The maximum of var(S) is at β = 2.4 (var 103), as in F6.4a (2.4); the published range is 2.2–2.5.
  - **Phase II, β = 3.5.** S −41.27 ± 0.07, χ 0.566, h 3.94. F6.4a measured −41.07, 0.593 and 3.70; F6.4c's
    preparations gave h 3.71–3.89; Surya published −41.37, 0.579 and 4.18.
    - The four chains differ: h 3.73–4.05 and χ 0.546–0.591. At β = 4.0 and 4.5 χ moves non-monotonically
      (0.603, 0.598).
    - The crystalline phase has metastable structures that neither sampler equilibrates within a chain, so the two
      samplers need not agree on h and χ there. This matches the bimodality found in F6.4c.
  - **Reading (decides nothing).**
    - The local sampler reproduces the manifold-like phase and the location of the transition of F6.4a to within
      their errors.
    - In phase II it lands within the spread of the known metastable states, not on one value.
    - The local moves are adequate for the escape at N = 50.

**P1, the escape: a size scan (the author, 19 Nov, 12:34: option 1, the more reliable one).** Local sampler at
N = 50, 100, 200, 400, ε = 0.12, run by diag/ov1_2d_mcmc2.py.
- **Engine.** Dense: the matrix C of the interval counts, so ΔS costs O(N).
  - It agrees with brute force to 3e−13 over 6 000 proposals at N = 12, 40 and 90.
  - It is statistically consistent with the sparse engine of P0 at N = 50 (β 0.1, 2.0).
  - Cost per sweep on one core: 0.02 s at N = 100, 0.07 s at 200 and 0.6 s at 400.
- **Sampling.** That of F6.4a: 4 chains per β, 500 + 3 000 sweeps, one sweep = N(N − 1)/2 proposals.
- **β grid, fixed before the run.** The transition moves with N (at N = 200 the chain is crystalline already at
  β = 1), so the F6.4a grid is not reused.
  - **Stage A**, the same for every N: 0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.4, 2.0, 2.8, 4.0.
  - **Stage B:** 9 equally spaced interior points between the stage-A neighbours of the stage-A maximum of var(S).
    This rule is fixed now and not chosen afterwards.
- **Calibration at each N.** 20 000 uniform 2D orders of the same N, drawn directly, seed 100 000 + N, disjoint from
  the chain seeds; the bands are [0.5 %, 99.5 %].

  | N | χ | d_MM | h |
  |---|---|---|---|
  | 50 | [0.376, 0.625] | [1.70, 2.37] | [7, 14] |
  | 100 | [0.414, 0.587] | [1.78, 2.25] | [12, 20] |
  | 200 | [0.440, 0.561] | [1.85, 2.17] | [19, 29] |
  | 400 | [0.457, 0.545] | [1.89, 2.12] | [30, 41] |

  - At N = 50 the crystalline χ (0.55–0.60) lies inside the flat χ band. There χ alone does not classify; h does
    (crystalline 3.7–4.2 against the band [7, 14]).
  - This is the author's point: the flat band fixes the manifold-like reference only, and the separate F6.4
    diagnostics classify the departure.
- **Reported per chain:** the means, var(S), the acceptance, τ_int of S and of χ, and the fraction of samples outside
  each band.
- **Results at N = 50 and 100** (Claude's container, 19 Nov; diag/out/ov1_2d_mcmc2/report_50_100.txt). "Out h" is
  the fraction of samples outside the flat h band; χ and d_MM leave their bands only at N = 100.

  | N | β | S | χ | h | out h |
  |---|---|---|---|---|---|
  | 50 | 2.2 | −2.4 | 0.493 | 8.3 | 0.05 |
  | 50 | 2.4 | −14.6 | 0.525 | 6.9 | 0.42 |
  | 50 | 2.6 | −29.4 | 0.555 | 5.4 | 0.89 |
  | 100 | 1.0 | −3.4 | 0.478 | 13.9 | 0.03 |
  | 100 | 1.4 | −181.7 | 0.579 | 5.9 | 1.00 |

  - **N = 50.** The manifold-like phase holds up to β ≈ 2.2; the escape (out h crossing 0.5) lies between β 2.4 and
    2.6; var(S) peaks at 2.4 (104).
  - **N = 100.** The escape lies between stage-A points 1.0 and 1.4: S jumps from −3.4 to −182, h from 13.9 to 5.9.
  - **In the crystalline phase** at N = 100 the chains do not equilibrate: τ_int of S is up to 150, the acceptance
    0.1–0.2, and the chains differ, so the phase-II values there are not reliable. The manifold side and the location
    of the escape are (τ_int ≤ 3 for β ≤ 1).
  - **The stage-B rule failed at N = 100.** The stage-A maximum of var(S) fell at β = 4.0 (163, a metastable
    crystalline chain), not at the transition, so stage B was run on [2.8, 4.0] instead of [1.0, 1.4]. The rule as
    fixed is kept on record; a replacement rule is the author's decision (proposed: stage B between the stage-A
    neighbours where out h crosses 0.5, the quantity of the escape itself). At N = 50 the old rule gave [2.0, 4.0]
    (maximum of var(S) at 2.8), which happened to contain the escape; the proposed rule gives [2.0, 2.8].
- **Stage A at N = 200 and 400** (the author's machine, 19 Nov, 4 chains × (500 + 3 000) sweeps, about 1 590 s per
  chain at N = 400):

  | N | β | S | χ | h | out h | τ_int(S) |
  |---|---|---|---|---|---|---|
  | 200 | 0.7 | −11.5 | 0.486 | 21.5 | 0.03 | 0.9 |
  | 200 | 1.0 | −911 | 0.638 | 6.6 | 1.00 | 255 |
  | 400 | 0.5 | −24.9 | 0.491 | 32.3 | 0.08 | 0.7 |
  | 400 | 0.7 | −41.2 | 0.489 | 31.4 | 0.18 | 0.8 |
  | 400 | 1.0 | −3 005 | 0.643 | 6.9 | 1.00 | 548 |

  - **Escape interval.** At both sizes the escape lies between β 0.7 and 1.0. The sequence of escape intervals is
    2.4–2.6 (N = 50), 1.0–1.4 (100), 0.7–1.0 (200) and 0.7–1.0 (400). It is not 1/N: between 200 and 400 it does not
    move at the resolution of stage A.
  - **Below the escape** the chains mix well (τ_int ≤ 1.4) and all three diagnostics lie inside their flat bands. At
    N = 400, β = 0.7, h has started to fall (31.4 against the flat 34.7; 18 % of the samples are outside the band).
  - **Above it** the chains do not equilibrate: τ_int of S is 200–550, the acceptance 0.08–0.3, and ⟨S⟩ is not
    monotonic in β (at N = 400, −3 005 at β 1.0 against −1 459 at 4.0). The crystalline values are not reliable.
  - **What the escape measures.** The chains start from flat orders and move by local moves. The observed escape is
    therefore where the manifold-like phase stops being stable under local moves within 3 500 sweeps. It is an upper
    bound for an equilibrium transition, not the transition itself; hysteresis from the crystalline side is not
    measured.
  - **Decision (the author, 19 Nov, 17:24).** Stage B uses the rule `outh`: 9 interior points between the stage-A
    neighbours where the fraction of samples outside the flat h band crosses 0.5. The failed rule `varS` stays on
    record and its stage-B chains at N = 50 and 100 are reported, not used.
  - **Stage B (`outh`) at N = 50 and 100** (Claude's container; diag/out/ov1_2d_mcmc2/report_50_100_outh.txt). β_e is
    the point where out h crosses 0.5, by linear interpolation between neighbouring β.

    | N | β | out h | | β | out h | β_e |
    |---|---|---|---|---|---|---|
    | 50 | 2.40 | 0.42 | | 2.48 | 0.71 | 2.42 |
    | 100 | 1.24 | 0.27 | | 1.28 | 0.78 | 1.26 |

    - **N = 50: a smooth crossover.** S, χ and h change continuously over β 2.2–2.6, var(S) peaks at 2.40 (104), and
      τ_int is at most 34.
    - **N = 100: abrupt.** Up to β = 1.20 everything is inside the bands (out h ≤ 0.06, τ_int ≤ 2.3). Between 1.24 and
      1.32 the chains split, with τ_int of S 160–320, var(S) up to 3 600 and ⟨S⟩ jumping from −6 to −139. That is
      coexistence, or metastability between the two phases.
    - N = 50 lies at about twice β_e of N = 100.
  - **Stage B (`outh`) at N = 200 and 400** (the author's machine, 19 Nov, finished 21:2x;
    diag/out/ov1_2d_mcmc2/report_200_400.txt, md5 69d4e231). β_e is the out-h crossing; "escaped" counts the chains
    with out h ≥ 0.8, out of 4.

    | N | β | out h | escaped | β | out h | escaped | β_e |
    |---|---|---|---|---|---|---|---|
    | 200 | 0.94 | 0.12 | 0/4 | 0.97 | 0.56 | 2/4 | 0.97 |
    | 400 | 0.82 | 0.47 | 1/4 | 0.85 | 0.57 | 2/4 | 0.83 |

    - **Escape per chain.** At both sizes the escape happens chain by chain, not as a gradual shift: a chain is either
      still manifold-like or has jumped, to S ≈ −900 at N = 200 and −2 100 to −3 400 at N = 400. At N = 400 all four
      chains have escaped at β = 0.88. This is first-order-like metastability under local moves.
    - **Before the escape at N = 400** h is already depressed: 31.0 against the flat 34.7 at β = 0.79, with 20–25 % of
      the samples below the band, while χ and d_MM stay at the flat values.

  **Escape point over the sizes** (local moves, ε = 0.12, 3 500 sweeps):

  | N | 50 | 100 | 200 | 400 |
  |---|---|---|---|---|
  | β_e | 2.42 | 1.26 | 0.97 | 0.83 |

  - β_e falls more slowly than 1/N. The local exponents of β_e ∝ N^(−a) are 0.94, 0.38 and 0.22 for the three
    doublings.
  - **The form β_e = β_∞ + A/N** fitted to N = 100 and 200 gives β_∞ = 0.68, A = 58, and predicts 0.83 at N = 400, as
    measured. N = 50 does not follow it (predicted 1.84, measured 2.42).
    - This is one degree of freedom, so it is consistent with a finite limit near 0.7 but does not establish one.
    - A further size (N = 800) would test it: predicted β_e ≈ 0.75, at about 8× the cost of N = 400 per sweep.
  - **Reading (development, decides nothing).**
    - Within the measured range the BD action at ε = 0.12, used as a weight exp(−βS), keeps the flat order in the
      manifold-like phase only for β < β_e(N). Above it, local moves reach crystalline orders with a much lower
      action (ΔS of order −10³ at N = 400).
    - The flat order is therefore not a minimum of S_BD under family I, which is why S-a is rejected (section 7), and
      it is stationary only in the sense of S-b.
    - **For OV3:** any dynamics driven down S_BD must be shown not to run into these orders. The escape point is a
      stability limit under local moves within the sweeps used, an upper bound for an equilibrium transition.
  - Before the decision: Under the proposed rule both brackets are
    [0.7, 1.0]; under the old rule N = 200 would get [2.8, 4.0] (maximum of var(S) at 4.0), a second failure.
- **Who runs what.** N = 50 and 100 run in Claude's container. N = 200 and 400 run on the author's machine:
  stage A plus B at N = 400 is about 6 h on 8 workers, at N = 200 under 1 h.

## 11. Verdict rules

- **Separate verdicts.** Each control in each family and each level is a separate card.

  | verdict | when |
  |---|---|
  | PASSED | |x| + 2 SE within the tolerance, and the instruments OK |
  | FAILED | outside the tolerance, with the instruments OK |
  | INSTRUMENT_FAILED | an instrument check failed |
  | UNDECIDED | the errors do not allow a decision |
- **Instrument checks:**
  - brute-force S and ΔS on small orders, which the development script runs (`--check`);
  - quadrature against K.4;
  - the pure-gauge null;
  - the flat null of the phase bands;
  - invariance under doubling the collar.
- **Phase.** The phase classification (section 10) is reported next to the verdict and never merged with it.
- **Locked before the first run:** n_loc, w_∂, the collar, the move set and the phase bands. None is changed after a
  result (the author).
- **Power.** Every null is sized so that 2 SE is at most a third of its tolerance, and its false-alarm rate is stated on
  the card. The controls of F5.2c seed 1 and of cl14 were underpowered.
- **R132.**

## 12. Order of dimensions

1. **2D:** family I, the phase bands, the flat vacuum and the dependence on T, using the instruments of F6.4a/c.
2. **2+1:** the same frozen logic, with no retuning of the move set.
3. **3+1:** family M against F6.8, and family I only for N ≥ 32 000 (K.1). Smaller orders serve development and
   checks of the implementation, not a verdict.

## 13. Development in 2D (diag/ov1_2d.py; decides nothing)

**Model.**
- The diamond is [0, 1]² in light-cone coordinates (u, v), sprinkled by a Poisson process of mean N_exp. x ≺ y iff
  u_x < u_y and v_x < v_y.
- An interval of proper time τ (τ² = Δu Δv) holds N_exp τ² elements on average.
- The action is S = 4ε (N − 2ε Σ_n N_n f(n, ε)), the F6.4a convention.

**Abundances.** The N_n are counted exactly for every pair with |f(n)| ≥ 1e−14 (n ≤ n_cut), by sweeps with an exact
stopping rule: the square at the corner of a rectangle bounds the count of every later rectangle.

**Nested levels.** The diamonds T/w = 32, 16, 8 are the top-corner squares of side 1, 1/2 and 1/4 of one sprinkling,
so they share the future tip q. A pair belongs to every level that contains its lower element.

**Moves.** ΔS of each adjacent transposition is computed locally, from the small-interval futures and pasts of the
two elements.

**Instrument checks (`--check`).** S and every ΔS are compared with brute force, the interval counts by matrix
product:
- small orders, N = 102–294, all swaps;
- larger orders, N = 1 804 and 2 549, where n_cut < N and the stopping rule is active, 12 swaps each.

The largest difference was 1e−10 for S and 5.5e−13 for ΔS.

**Cost.** On one core of the container, one sprinkling with all three levels and the balance on T/w = 8 and 16 takes
about 9 s at ε = 0.24, 40 s at 0.12 and 180 s at 0.06.

**Development rows** (ε = 0.12, 16 sprinklings; diag/out/ov1_2d; mean ± SE, with the sd per sprinkling).

*Vacuum.*

| T/w | N | ⟨S⟩ | sd |
|---|---|---|---|
| 8 | 4 167 | +4.3 ± 8.5 | 34 |
| 16 | 16 650 | −18.2 ± 17.1 | 69 |
| 32 | 66 606 | −17.4 ± 29.8 | 119 |

- The sd of S grows as √N (×2 per level).
- The paired differences are S(16) − S(8) = −22.5 ± 15.5 and S(32) − S(16) = +0.8 ± 25.8.
- The slope is c = dS/dN = (−2.5 ± 4.8)·10⁻⁴, that is c N_max = −16 ± 32. Relative to the bare volume term 4ε per
  element this is c/(4ε) = (−5 ± 10)·10⁻⁴: the N term cancels to 10⁻³, consistent with c = 0.
- With 16 sprinklings the constant ⟨S⟩ (the joint term, section 4.4) is not resolved.

*Phase nulls (flat).* χ = 0.5001 / 0.5004 / 0.4999 and d_MM = 1.9997 / 1.9989 / 2.0003 at T/w = 8 / 16 / 32. The
sd of d_MM is 0.015 / 0.0065 / 0.0032, falling as N^(−1/2). h = 121 / 248 / 503.

*Balance, interior elements.*

| T/w | ⟨D S⟩, U | ⟨D S⟩, V | sd per element |
|---|---|---|---|
| 8 | +0.0002 ± 0.0003 | −0.0000 ± 0.0002 | 0.18 |
| 16 | +0.00013 ± 0.00008 | +0.00000 ± 0.00007 | 0.18 |

- The slope of ⟨D S⟩ in time is zero within ≤ 2e−3.
- The profile from the past to the future tip is flat within the errors; the end bins are at most 1.4σ.
- ⟨D² S⟩ (reported) is +0.0012 ± 0.0030 / −0.0057 ± 0.0027 (U / V) at 8 and −0.0017 ± 0.0009 / −0.0004 ± 0.0012
  at 16, against an sd of 0.38.

**Validation rows on the author's machine, ε = 0.24** (512 sprinklings, 19 Nov; diag/out/ov1_2d on the author's
machine; mean ± SE, with the sd per sprinkling).

*Vacuum.*

| T/w | N | ⟨S⟩ | sd |
|---|---|---|---|
| 8 | 1 921 | +3.1 ± 3.1 | 70 |
| 16 | 7 678 | +2.8 ± 6.9 | 155 |
| 32 | 30 722 | +8.7 ± 14.9 | 337 |

- The paired differences are S(16) − S(8) = −0.3 ± 6.1 and S(32) − S(16) = +5.9 ± 13.1.
- The slope is c = (+2.1 ± 5.1)·10⁻⁴, that is c N_max = +6.5 ± 15.8 and c/(4ε) = (+2.2 ± 5.4)·10⁻⁴.
- ⟨S⟩ is flat at about 3 over a 16-fold volume, consistent with a constant (joint) term and no volume term.

*Balance, interior elements.*

| T/w | ⟨D S⟩, U | ⟨D S⟩, V | sd per element |
|---|---|---|---|
| 8 | −0.00043 ± 0.00022 | +0.00014 ± 0.00025 | 0.77 |
| 16 | −0.00007 ± 0.00007 | +0.00001 ± 0.00007 | 0.76 |

- The time slopes are zero within 3e−3.
- ⟨D² S⟩ is within 1.2σ of 0, against an sd of 1.7.

*Phase.* d_MM = 2.000 at every level (sd 0.020 / 0.010 / 0.005). The calibration seeds (≥ 100 000) are not run yet.

**Validation rows on the author's machine, ε = 0.12** (512 sprinklings, 19 Nov, finished 21:5x).

*Vacuum.*

| T/w | N | ⟨S⟩ | sd |
|---|---|---|---|
| 8 | 4 164 | +6.0 ± 1.5 | 34 |
| 16 | 16 639 | +6.2 ± 3.2 | 73 |
| 32 | 66 572 | +1.4 ± 7.5 | 169 |

- The paired differences are S(16) − S(8) = +0.2 ± 2.9 and S(32) − S(16) = −4.8 ± 6.6.
- The slope is c = (−0.8 ± 1.2)·10⁻⁴, that is c N_max = −5.3 ± 7.9 and c/(4ε) = (−1.7 ± 2.5)·10⁻⁴.
- There is no volume term: the volume cancels per element to a few 10⁻⁴.
- The constant is resolved at T/w = 8: ⟨S⟩ = 6.0 ± 1.5 at ε 0.12, against 3.1 ± 3.1 at ε 0.24. Whether it scales as
  1/ε waits for ε 0.06.

*Balance, interior elements.*

| T/w | ⟨D S⟩, U | ⟨D S⟩, V | sd per element |
|---|---|---|---|
| 8 | −0.00004 ± 0.00003 | +0.00002 ± 0.00003 | 0.18 |
| 16 | −0.00001 ± 0.00001 | +0.00001 ± 0.00001 | 0.18 |

- The time slopes are zero within 3e−4, and the profiles from the past to the future tip are flat within 1e−4.
- ⟨D² S⟩ (reported) is −0.0009 ± 0.0005 / −0.0005 ± 0.0005 (U / V) at 8 and −0.00043 ± 0.00023 / −0.00011 ± 0.00022
  at 16, against an sd of 0.38. The values are slightly negative, below 2σ each.

**Validation rows on the author's machine, ε = 0.06** (256 sprinklings, 19 Nov, finished 23:1x; 164 s each).

*Vacuum.*

| T/w | N | ⟨S⟩ | sd |
|---|---|---|---|
| 8 | 8 642 | +2.8 ± 1.0 | 16 |
| 16 | 34 551 | +0.4 ± 2.3 | 37 |
| 32 | 138 228 | +4.5 ± 5.2 | 84 |

- The paired differences are S(16) − S(8) = −2.4 ± 2.1 and S(32) − S(16) = +4.1 ± 4.9.
- The slope is c = (+2.1 ± 4.1)·10⁻⁵, that is c N_max = +2.9 ± 5.7 and c/(4ε) = (+0.9 ± 1.7)·10⁻⁴.

*Balance.*
- ⟨D S⟩ = +0.00001 ± 0.00001 at T/w 8 and 0.00000 ± 0.00000 at 16, in both orders, against a per-element sd of 0.045;
  the profiles are flat.
- ⟨D² S⟩ = +0.00012 ± 0.00012 / +0.00011 ± 0.00012 at 8 and +0.00006 ± 0.00005 / +0.00002 ± 0.00005 at 16, now
  slightly positive. The small negative values at ε 0.12 are therefore not a systematic sign.

**The three ε together (vacuum, validation seeds).**

| ε | c/(4ε) | ⟨S⟩ at T/w = 8 | sd of S at T/w = 8 |
|---|---|---|---|
| 0.24 | (+2.2 ± 5.4)·10⁻⁴ | +3.1 ± 3.1 | 70 |
| 0.12 | (−1.7 ± 2.5)·10⁻⁴ | +6.0 ± 1.5 | 34 |
| 0.06 | (+0.9 ± 1.7)·10⁻⁴ | +2.8 ± 1.0 | 16 |

- **No volume term at any ε.** The per-element cancellation is a few 10⁻⁴, and the slopes are consistent with zero
  within 1σ.
- **The constant does not scale as 1/ε.** At ε 0.06, 1/ε scaling would predict ≈ 12, and 2.8 ± 1.0 is measured. The
  three values agree at 1.8σ (ε 0.12 against 0.06); their weighted mean is 3.7 ± 0.8.
- **The sd of S scales as ε√N.** It halves with ε at fixed T/w and grows ×2 per level. The constant is resolved best
  at small ε.

**Reading (development, decides nothing).**
- A zero mean of D S over a time-symmetric interior is guaranteed by the time reversal of the flat diamond, so in flat
  space it checks the implementation, not physics. The content of L2 lies in the profile near the tips, and later in
  curved or sourced orders.
- The content of the vacuum control is c(T). The fluctuation of S grows as √N, so a tolerance on c N_max of the order
  of the joint term needs thousands of sprinklings. A tolerance on c/(4ε), the cancellation of the volume term per
  element, needs hundreds.
- Which scale the card uses is the author's decision, to be taken after the run of section 14.

## 14. Development run on the author's machine (decides nothing)

```
python -m diag.ov1_2d --check
python -m diag.ov1_2d --eps 0.24 --seeds 0-511 --levels 8,16,32 --drift 8,16 --workers 8
python -m diag.ov1_2d --eps 0.12 --seeds 0-511 --levels 8,16,32 --drift 8,16 --workers 8
python -m diag.ov1_2d --eps 0.06 --seeds 0-255 --levels 8,16,32 --drift 8,16 --workers 8
python -m diag.ov1_2d --eps 0.24 --seeds 100000-100511 --levels 8,16,32 --drift "" --workers 8
python -m diag.ov1_2d --eps 0.12 --seeds 100000-100511 --levels 8,16,32 --drift "" --workers 8
python -m diag.ov1_2d --eps 0.06 --seeds 100000-100255 --levels 8,16,32 --drift "" --workers 8
python -m diag.ov1_2d_summ > diag/out/ov1_2d/summary.txt
```

- **Calibration of the phase bands (section 10):** seeds ≥ 100 000, without the balance part (--drift ""), so they are
  cheaper.
- **Phase sampler P0** (section 10):
  ```
  python -m diag.ov1_2d_mcmc --check
  python -m diag.ov1_2d_mcmc --scan --N 50 --eps 0.12 --chains 4 --therm 500 --meas 3000 --workers 8
  ```
  This takes about 0.01 s per sweep and chain, about 6 min on 8 workers.

- **Estimated time on 8 workers**, scaled from one core of the container: about 10 min, 45 min and 1.5–2 h.
- **Resumable:** an existing JSON is skipped.
- **Expected precision at 512 sprinklings (ε = 0.12):** c N_max ± 6 and c/(4ε) ± 2·10⁻⁴; ⟨D S⟩ at T/w = 16 to about
  ±2·10⁻⁵.

## References

- Benincasa and Dowker, Phys. Rev. Lett. 104, 181301 (2010).
- Surya, Class. Quantum Grav. 29, 132001 (2012) (calibrated in F6.4a/c).
- Buck, Dowker, Jubb and Surya, arXiv:1502.05388.
- Rideout and Sorkin, Phys. Rev. D 61, 024002 (2000).
- Laboratory: substrate_and_action_audit.md §7–8; variational_audit.md; f68f_dev.md; findings F6.8a, b2, c2, e3, f;
  F6.4a, F6.4c; K.1, K.4, K.17.
