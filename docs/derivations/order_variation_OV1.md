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

- **Reference values.** ε₀ is the reference value, with the comparison values ε₀/2 and 2ε₀, all inside the admissible
  interval 0 < ε < 1 of the smeared kernel. If one fell outside it, both factors would shrink equally on a log scale.
  ε is not chosen per dimension.
- **Open point for the author: what is ε_F6.8?**
  - F6.8 works in the units ερ = 1 (ℓ_k = 1), and its expectation depends on ε and ρ only through ερ. F6.8 therefore
    has no per-element ε: it fixes the unit of length, not a dimensionless number.
  - A discrete realisation (family I, L2) needs a per-element ε, the ratio of discreteness to nonlocality.
  - **Provisional value for 2D:** ε₀ = 0.12, the value of F6.4a/c (Surya 2012), with 0.06 and 0.24 for comparison.
    This is the existing calibrated construction in 2D.
  - In 3+1, ε₀ must be chosen before the first 3+1 run. Only ℓ_k = 1 carries over from F6.8.
  - **To be confirmed** before the 2D card is frozen.

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

**Open point for the author: the size.** The decided sizes (T/w = 8, 16, 32 with w from 99 % of the kernel) give
N = 4 160 – 66 560 at ε = 0.12. F6.4a/c calibrated its phases at N = 50 only, so:
- the F6.4 diagnostics are defined at every N, but the crystalline reference values (χ 0.579, h 4.18) exist at N = 50
  only;
- the manifold-like band at each N is set from flat sprinklings at that N (development, section 13), before any card;
- a crystalline reference at large N would need new MCMC runs, which are not made here.

**What OV1 can and cannot test.** In OV1 no procedure moves an order further than one rewrite from a flat sprinkling.
- In OV1 output B is therefore a calibration of the classifier: the bands at each N, and the check that H₊ and H₋ stay
  inside them.
- The escape proper needs a procedure that moves through configuration space: a dynamics (OV3), or an ensemble
  weighted by the action with the local rewrite set of section 3.2 (an F6.4-type run with local moves).
- Which of the two is used is the author's choice.

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

**Cost.** At ε = 0.12, one sprinkling with all three levels and the balance on T/w = 8 takes about 30 s on one core.

**Results.** [to be filled from the development runs]

## References

- Benincasa and Dowker, Phys. Rev. Lett. 104, 181301 (2010).
- Surya, Class. Quantum Grav. 29, 132001 (2012) (calibrated in F6.4a/c).
- Buck, Dowker, Jubb and Surya, arXiv:1502.05388.
- Rideout and Sorkin, Phys. Rev. D 61, 024002 (2000).
- Laboratory: substrate_and_action_audit.md §7–8; variational_audit.md; f68f_dev.md; findings F6.8a, b2, c2, e3, f;
  F6.4a, F6.4c; K.1, K.4, K.17.
