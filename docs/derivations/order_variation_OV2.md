# OV2: matter on the order — the flow F_v and the action S_m[H, ·] (specification, draft)

20 Nov 2026. Status: a draft specification. No runs, no frozen criteria; the open choices are listed in section 8 and
are the author's. It follows OV1 (docs/derivations/order_variation_OV1.md; the vacuum control F6.9a is running).

## 1. Scope

OV2 has two branches, which stay separate until section 6 joins them (the author, 19 Nov):

- **F.** The operational definition of the flow F_v on the order, by the axiom of causal continuity:
  Σ_in F + s_v = Σ_out F + ΔQ_v.
- **S.** An action of matter S_m[H, ·] that depends on the order only through the order. Matter must act on the order
  because its own action depends on the order, not through an inserted source ρ_s (the author, 19 Nov).

**OV2 does not:**
- solve for a sourced order or derive Π_M (OV3);
- choose the selection principle of OV3 (stationarity of one order, or the ensemble level; see section 7);
- rename or merge the three currents of field_and_throughput.md §3 (j_top, T^{μν}, F), which must not be identified
  because they behave alike.

## 2. What any S_m must give: the known answer

In the continuum the source of gravity is T_ab, the response of the matter action to the metric at fixed matter:

    T_ab = −(2/√|g|) δS_m/δg^{ab}   (sign convention to be fixed on the card).

This is a definition at fixed matter (off shell). It does not need the matter to obey its equation of motion; the
conservation ∇_a T^{ab} = 0 does.

**Requirement.** In family M of OV1, the first variation of the ensemble mean ⟨S_m⟩ under δg, at fixed matter
configuration, must reproduce ½ ∫ √|g| T^{ab} δg_ab of the continuum matter, to a stated accuracy and in a stated
range of ℓ_k. This is the OV2 counterpart of OV1's F6.8 consistency (OV1 §9). It tests the coupling before any
response of the order is computed.

## 3. Candidates for S_m

### 3.1 S_m^φ: the action of the field C0 on the order

    S_m^φ[H, φ] = (1/ρ) Σ_v [ ½ φ_v (B̃φ)_v − V(φ_v) ],   V = μ (1 − cos φ).

- **The kinetic operator, B̃φ = Bφ − φ·(B1).** The smeared BD operator B is written in difference form,
  (B̃φ)(x) ∝ Σ_{y≺x} f(n_xy) (φ_y − φ_x). On a curved order B1 ≠ 0, and the continuum mean of B1 is proportional to
  the scalar curvature R; its sum is the BD action itself. The plain form ½ φ Bφ therefore contains a coupling
  ∝ R φ².
  - For an angle-valued field (sine-Gordon, φ ~ φ + 2π) this breaks the periodicity: the action of the same physical
    configuration would depend on the winding representative.
  - The difference form removes the φ² R term and respects φ → φ + 2π. It is the minimal coupling.
  - Which form holds is checked, not assumed: test T0 (section 5).
- **Symmetric operator.** The action needs a symmetric operator B̄ = (B̃ + B̃ᵀ)/2. The laboratory's dynamics C0 is
  retarded, and the retarded rule is not the Euler–Lagrange equation of any action (a known property of causal-set
  field theory).
  - **Resolution proposed for OV2.** S_m^φ is used as the coupling functional at fixed φ. Its first variation in H
    defines the source whatever dynamics produced φ (section 2: the definition is off shell).
  - The question "does the retarded C0 conserve the T_ab read from S_m^φ" is separate (test T3). It is not assumed.

### 3.2 S_m^τ: the proper time of a flow counted by the order (A6)

    S_m^τ[H, γ] = −Σ_flows m τ_H(γ),

where τ_H is the discrete proper time along the history γ of the flow: the longest chain, or the volume of the
diamond, in ρ^(−1/d) units (the clocks of K.16).
- **Advantage.** It is the axiom of the action of a flow itself, and it gives ½ T_ab directly at linear order
  (rule_from_action.md §2).
- **Cost.** It needs the history γ of a flow on the order, the world line or world sheet of a defect. For a wall the
  world sheet of the kink is a readout of φ (the zero of φ − π), so S_m^τ is not independent of S_m^φ for walls.
- **Comparison.** For a wall, both must give the same T_ab in the continuum: the energy of the kink equals its tension.

### 3.3 Branch F: the flow F_v

- **Natural candidate.** The C0 rule distributes influence along relations, φ(v) − δJ(v) = φ₀(v) − c Σ_{u<v} w_n J(u),
  which suggests the flow on a relation u → v as F_{u→v} = c w_{n(u,v)} J(u).
- **Tests:**
  - **F1, conservation:** does it satisfy the axiom Σ_in F + s_v = Σ_out F + ΔQ_v at every event? Identically, to
    numerical precision, or not at all.
  - **F2, which current:** what is its coarse-grained continuum limit, j_top, T^{0μ}, or neither?
- **The outcome decides whether S_m can be written through F (the author's S_m[H, F]).** If F is neither T^{0μ} nor a
  function of it, S_m[H, F] cannot source gravity at linear order. That would be a result about the axiom of causal
  continuity, not a failure of OV2.

## 4. Order of dimensions

As OV1: 2D, then 2+1, then 3+1.

**In 2D only the trace of T_ab can be tested through the order.** Every 2D metric is locally conformal to flat, and a
conformal factor does not change the order, only the density. For g → e^{2σ} g:
- δS_m = ∫ √|g| σ T^a_a;
- for sine-Gordon T^a_a = −2V (the kinetic part is conformally invariant in 2D).

The 2D test therefore checks the measure (the density of elements) and the conformal invariance of the kinetic term
in expectation. The full T_ab, in particular the kinetic part through the cones, needs ν perturbations (2+1 or 3+1,
the basis of F6.8).

## 5. Tests (development first, cards later)

| | dimension | question | known answer |
|---|---|---|---|
| T0 | 2D | Does S_m^φ in difference form keep φ → φ + 2π (the same configuration, two winding representatives) exactly, and the plain form not? | difference form: invariant; plain form: shifts by ∝ Σ B1 φ² |
| T1 | 2D | Conformal δg = 2σ g with compact support, a static kink at fixed φ: is δ⟨S_m^φ⟩ = ∫ σ T^a_a √|g| = −2 ∫ σ V √|g|? | the continuum trace |
| T2 | 3+1 (family M-b) | ω and ν perturbations of F6.8 at fixed φ (a wall): is δ⟨S_m^φ⟩ = ½ ∫ T^{ab} δg_ab, with T_ab of the wall? | the continuum T_ab, kinetic and potential parts |
| T3 | 2D then 3+1 | Does the retarded C0 dynamics conserve the T_ab of S_m^φ (∂_a T^{ab} = 0 coarse-grained), within the C2 precision? | 0 in the continuum; the C1 offset is the known nonlocal departure |
| F1 | 2D | Conservation of F_{u→v} at every event | exact or not |
| F2 | 2D | Coarse-grained limit of F against j_top and T^{0μ} on the same walls | which, if any |

The vacuum of OV1 (F6.9a) must pass before T1 and T2 are read. The instrument items reuse the OV1 machinery: exact
interval counts, brute force on small orders and nested diamonds.

## 6. Joining the branches

After T0–T3 and F1–F2:
- **A matter action S_m[H, ·] with a verified coupling** (it reproduces T_ab at linear order). This is the input of
  OV3.
- **A statement about F:** whether the axiom of causal continuity is realised by a current that sources gravity.

OV3 then combines S_BD and S_m under the principle chosen there.

## 7. A decision OV3 will need, raised now

OV1's phase sampler showed that the flat order is not a minimum of S_BD under local moves (family I): orders with an
action lower by 10³ exist, and the manifold-like phase is stabilised by the number of orders (entropy), below β_e(N).

A selection rule "stationary point of S_BD + S_m on one order" may therefore find solutions outside the manifold-like
phase. The alternatives are:
- the ensemble level of OV1 (L1): a saddle of ⟨S⟩ over a frozen measure μ;
- an explicit balance of action and entropy;
- a stochastic dynamics weighted by the action.

S_m enters each of them in the same way, so OV2 does not depend on the choice. OV3 does.

## 8. Decisions (the author, 20 Nov; decision 39)

1. Both candidates, S_m^φ and S_m^τ, with a comparison.
2. The kinetic operator in difference form, checked by T0.
3. S_m is the coupling functional at fixed φ; T3 is a separate test.
4. Order: T0 → T1 → F1/F2 in 2D, then T2 and T3 in 3+1.
5. Branch F is pursued in that order, after T1.

**Refinement found in the implementation (for the author to confirm).** "Difference form" has two realisations:
- the row form K = B − diag(B1);
- the pairwise (Dirichlet) form S_kin = −ε² Σ_pairs f(n) (φ_x − φ_y)².

Their quadratic forms differ by ε² Σ_v φ_v² (F_future(v) − F_past(v)), a time-asymmetric diagonal that does not
vanish near the past and future boundaries. Only the pairwise form is exactly invariant under φ → φ + 2πk; T0 confirms
it.

## 9. Development in 2D (diag/ov2_2d.py, diag/ov2_kernels.py; decides nothing)

**Model.**
- The light-cone square, ds² = −du dv, proper density ρ = 2N_exp, ε = 0.12.
- Fixed static kink φ = 4 arctan(e^{mx}), m = 12, V = m²(1 − cos φ).
- Conformal bump σ = a exp(−r²/(2s²)), a = ±0.2, s = 0.12. The configurations +a, 0 and −a are thinned from one
  master sprinkling, so they are coupled.
- Continuum: S_kin = ∫ φ_u φ_v du dv = −42.45 and S_pot = −42.45 (the virial equality of the kink).

| quantity | N_exp 4 160 (T/w 8), 256 sprinklings | N_exp 16 640 (T/w 16), 64 sprinklings | continuum |
|---|---|---|---|
| T0 shift 2π, Dirichlet | 2e−15 ± 5e−15 | 4e−14 ± 2e−14 | 0 exactly |
| T0 shift 2π, row form | −5.7 ± 4.1 | +21 ± 16 | — |
| T0 shift 2π, plain | −37 ± 42 | −532 ± 200 | — |
| flat S_kin (Dirichlet) | −38.0 ± 1.1 | −36.0 ± 4.7 | −42.45 |
| flat S_pot | −42.49 ± 0.06 | −42.46 ± 0.05 | −42.45 |
| T1 δS_pot (+a / −a) | −3.545 ± 0.020 / +2.783 ± 0.017 | −3.533 ± 0.019 / +2.803 ± 0.017 | −3.537 / +2.781 |
| T1 δS_kin Dirichlet (+a / −a) | +0.52 ± 0.43 / −0.65 ± 0.40 | −1.49 ± 1.78 / −1.66 ± 1.63 | 0 / 0 (conformal invariance) |
| first order (S(+a) − S(−a))/2 | −2.58 ± 0.30 | −3.08 ± 1.19 | −3.12 (= ∫ σ T^a_a √g) |

**Update (20 Nov, 14:2x): 2 048 sprinklings at T/w 8 and 320 at T/w 16.**

| quantity | T/w 8 (2 048) | T/w 16 (320) | continuum |
|---|---|---|---|
| flat S_kin (Dirichlet) | −37.65 ± 0.40 | −43.07 ± 2.08 | −42.45 |
| δS_kin (+a / −a) | −0.04 ± 0.15 / −0.43 ± 0.15 | +1.06 ± 0.80 / −0.96 ± 0.75 | 0 / 0 |
| first order of S_kin, (δ₊ − δ₋)/2 | +0.20 ± 0.11 | +1.01 ± 0.55 | 0 |
| second order of S_kin, δ₊ + δ₋ | −0.47 ± 0.21 | +0.10 ± 1.1 | 0 |
| first order of S_m | −2.954 ± 0.107 | −2.16 ± 0.56 | −3.122 |
| δS_pot (+a / −a) | −3.530 ± 0.007 / +2.772 ± 0.006 | −3.527 ± 0.008 / +2.801 ± 0.008 | −3.537 / +2.781 |

- **The flat S_kin is 11 % low at T/w 8 and agrees at T/w 16.** That is a finite-nonlocality effect that disappears
  with density.
- **The first-order conformal response of S_kin fell from +0.58 ± 0.30 (256 sprinklings) to +0.20 ± 0.11 (2 048).**
  It is 1.9σ, about 6 % of δS_m. Not resolved; the earlier value was partly a fluctuation.
- **On the author's machine, a = 0.1 at T/w 16** (1 024 sprinklings, 20 Nov).
  - δS_kin (+a / −a): +0.15 ± 0.30 / +0.16 ± 0.30.
  - The first order of S_kin is −0.01 ± 0.21; the first order of S_m is −1.573 ± 0.219 against the continuum −1.561.
  - δS_pot: −1.6568 ± 0.0033 / +1.4736 ± 0.0030 against −1.6598 / +1.4717.
  - The flat S_kin is −40.49 ± 1.11.
  - At T/w 16 no first-order conformal response of S_kin is seen. The T/w 8 hint (+0.20 ± 0.11 at a = 0.2) is
    consistent with a finite-ℓ_k effect.
- **T/w 32** (a = 0.2, 1 024 sprinklings, the author's machine, 20 Nov, 100 s each):
  - δS_kin (+a / −a): −0.02 ± 0.87 / −1.35 ± 0.90, so the first order of S_kin is +0.66 ± 0.62;
  - the first order of S_m is −2.50 ± 0.61 against −3.12;
  - δS_pot: −3.5369 ± 0.0023 / +2.7798 ± 0.0021 against −3.5371 / +2.7805;
  - the flat S_kin is −39.99 ± 2.28.

**First order of S_kin over the densities (development).**

| T/w | a | sprinklings | first order of S_kin |
|---|---|---|---|
| 8 | 0.2 | 2 048 | +0.20 ± 0.11 |
| 16 | 0.2 | 320 | +1.01 ± 0.55 |
| 16 | 0.1 | 1 024 | −0.01 ± 0.21 |
| 32 | 0.2 | 1 024 | +0.66 ± 0.62 |

- **No row departs from 0 by more than 2σ, and there is no trend with density.** The conformal invariance of the
  pairwise kinetic term holds within the errors.
- **The measure part is exact to 1e−4** at T/w 32.
- **The precision of the kinetic part is limited by fluctuations.** The sd per sprinkling of S_kin is 18, 35 and 73 at
  T/w 8, 16 and 32: it grows as √N, like the BD action. The ensemble mean converges, while a single order fluctuates
  (the L1/L2 hierarchy of OV1 §5). Sprinklings alone therefore cannot give a tight T1 at high density.

**Proposal for the T1 card: the ensemble level by quadrature (family M-b in 2D).**
- For n ~ Poisson(ρV), E[f(n, ε)] is exactly the continuum smeared kernel g(x) = e^{−x}(1 − 2x + x²/2), x = ερV.
- So ⟨S_kin⟩[σ] = −ε² ρ² ∫∫_{y≺x} √g_x √g_y g(ερ V_g(x, y)) (φ_x − φ_y)² is a deterministic integral over pairs, with
  V_g the interval volume in the conformal metric. ⟨S_pot⟩ is −∫ √g V.
- This gives T1 without noise at any ℓ_k, like F6.8 for S_BD.
- The sprinkled rows above stay the L2 check (the individual orders fluctuate, with no systematic direction).

**Reading.**
- **T0.** Only the pairwise form is invariant, exactly, as derived. The plain form fails by ∝ Σ B1 φ², and the row
  form by its asymmetric diagonal.
- **The measure part of T1 holds to 0.6 %.** δS_pot = −½ ∫ (e^{2σ} − 1) V du dv.
- **The kinetic part is the open item.**
  - At T/w 8 the flat S_kin is 10 % below the continuum (ratio 0.895 ± 0.027). This is plausible as the finite
    ℓ_k m = 0.38 of the kink; F6.8a found a non-analytic (kℓ)² ln correction.
  - The first-order conformal response of S_kin is +0.58 ± 0.30 (2σ) at T/w 8 and unresolved at T/w 16.
  - The noise of S_kin grows with N (sd 18 → 37 per sprinkling), like the BD action itself.
- **Before a T1 card:**
  - the convergence of the flat S_kin to −42.45 with ℓ_k m (a density scan);
  - the conformal invariance of S_kin to a stated precision, which needs about 10³ sprinklings per density.

## 10. T1 at the ensemble level by quadrature (diag/ov2_quad2d.py; development, decides nothing)

**Method.**
- ⟨S_kin⟩ = −ε² ∫∫_{x≺y} λ_x λ_y g(εΛ) (φ_x − φ_y)², with g(z) = e^{−z}(1 − 2z + z²/2), by the Mecke formula. This is
  exact for Poisson sprinklings.
- Variables for each x: t = ln Δu, and w = εN Δu Δv.
- The outer nodes cluster at the future edges.
- The first order in a is the analytic derivative under the integral at a = 0, d/da = ∫∫ N² [2(σ̂_x + σ̂_y) g +
  g′ εN · 2∫_I σ̂] (Δφ)². The rectangle integrals of the Gaussian are exact (erf), so the derivative is not a difference
  of two large numbers.
- **An error of the first version, corrected.** A fixed lower cutoff Δu ≥ 10⁻⁷ omitted the near-null pairs (links that
  cross the kink, with g = 1). The omitted part is ε² N² Δu_min ∫(1 − v)(Δφ)² and grows as N², which gave the false
  growth +259 at N_exp 266 240. The cutoff is now 10⁻¹⁰/(εN).

**Results** (a = 0.2, s = 0.12, m = 12; converged to 3–4 digits between the two finer grids).

| T/w | N_exp | flat ⟨S_kin⟩ | sprinkled | first order of ⟨S_kin⟩ | its share of ∫σT^a_a√g = −3.122 |
|---|---|---|---|---|---|
| 8 | 4 160 | −37.91 | −37.65 ± 0.40 | +0.270 | 8.6 % |
| 16 | 16 640 | −40.71 | −40.49 ± 1.11 | +0.124 | 4.0 % |
| 32 | 66 560 | −41.85 | −39.99 ± 2.28 | +0.048 | 1.5 % |
| 64 | 266 240 | −42.26 | — | +0.016 | 0.5 % |
| 128 | 1 064 960 | −42.40 | — | +0.004 to +0.005 | 0.15 % |

The continuum flat value is −42.45, and the continuum first order of S_kin is 0.

- **The flat kinetic action converges to the continuum.** The quadrature agrees with the sprinkled means at T/w 8, 16
  and 32 (0.6σ, 0.2σ, 0.8σ), which checks both instruments.
- **The conformal response of S_kin vanishes in the continuum limit,** falling by a factor 2.2–3 per doubling of T/w.
  At finite ℓ_k it is a positive nonlocal correction: 8.6 % of the response at T/w 8, below 2 % from T/w 32.
- **The sprinkled first orders agree with it within their errors.** At T/w 8, +0.20 ± 0.11 against +0.27; at T/w 16,
  a = 0.1, −0.01 ± 0.21 against +0.06; at T/w 32, +0.66 ± 0.62 against +0.05.
- **T1 in 2D therefore holds at the ensemble level.** The pairwise matter action responds to a conformal change of
  the metric as the continuum does: the measure part exactly, and the kinetic part conformally invariant up to a
  correction that vanishes with ℓ_k. The response is ∫ σ T^a_a √g with T^a_a = −2V.

**Proposal for a T1 card (the author's decision).** A deterministic quadrature on fresh parameters not used in the
development (another kink width, bump width and position), at T/w 32 and 64.
- **Criterion:** the first order of ⟨S_m⟩ equals ∫ σ T^a_a √g within 3 %, and the kinetic share is below 3 %.
- **Instruments:** grid convergence (two refinements), and the flat ⟨S_kin⟩ against the continuum.
- **Report:** the sprinkled L2 rows.

## 11. F1 for the candidate of section 3.3: analytic (decides nothing until the author reads it)

**The rule and the candidate.** The retarded rule is φ(v) = φ₀(v) − Σ_{u≺v} K(u, v) J(u):
- in 1+1, K = 1/(2ρ) on every relation;
- in 3+1 (C0), K = c w_n on the relations with n ≤ 8;
- J = V′(φ).

The candidate flow is F_{u→v} = K(u, v) J(u).

**The two sides of the balance.**
- Inflow, exact by the rule: Σ_in F(v) = φ₀(v) − φ(v).
- Outflow: Σ_out F(v) = J(v) κ(v), with κ(v) = Σ_{w≻v} K(v, w).
- In 1+1, κ(v) = |J⁺(v) ∩ Ω| / (2ρ), on average half the volume of the future of v inside the region.
- In 3+1, κ(v) = c Σ_n w_n N⁺_n(v), where N⁺_n(v) counts the future elements at interval size n ≤ 8. That number grows
  with the extent of the region along the light cone, through the long thin intervals.

**Consequence.**
- The outflow of an event is set by the future boundary of the region, not by anything at the event. No local
  storage ΔQ_v and source s_v satisfy Σ_in F + s_v = Σ_out F + ΔQ_v: the would-be ΔQ_v depends on where the region
  ends.
- The influence of the retarded rule is therefore a propagated potential, a superposition of Green functions, not a
  conserved flow.
- **F1 for this candidate: the axiom is not satisfied at all.** This is a statement about the candidate, not about the
  axiom.

**What a test of the axiom needs instead.**
- **Slice conservation.** For two slices A < A′ of the order (maximal antichains, or time cuts), the flux through A
  plus the sources between them equals the flux through A′ plus what is stored between them. This is the discrete
  form of ∂_μ j^μ = s, and a flow that satisfies the axiom passes it.
- **Two candidates with a conservation law** (field_and_throughput.md §3 forbids identifying them by fiat):
  - **F-a, energy and momentum from S_m.** The slice energy is E(A) = ε² Σ_{pairs crossing A} f(n)(Δφ)² plus the
    potential of the elements next to A, with the normalisation of the continuum ∫ T^{00}. It is tested on dynamic C0
    configurations (a moving kink, a kink–antikink pair): is E(A) the same on every slice within the C2 precision,
    and does it equal ∫ T^{00} dx? This is the current that can source gravity (it is the T_ab of T1 and T2).
  - **F-b, the topological current.** The winding of φ across a slice. It is conserved by the boundary values of φ,
    and it counts defects; it does not source gravity at linear order.
- **Which of them is the F of the axiom of causal continuity** (or neither, if the axiom is about something else) is a
  decision about the ontology, the author's.

## 12. T2 in 3+1, stage 1: the uniform conformal sector (diag/ov2_t2.py, diag/ov2_kernels4.py; development)

**Model.**
- A static wall φ(z) = 4 arctan(e^{mz}) on a flat 3+1 order of age T, in units ερ = 1 (ℓ_k = 1), as F6.8.
- Pairwise kinetic action by the Mecke formula, per unit time and area:
  s_kin = (a₄/4) ∫dz ∫_{J⁻(x), age T} g(ερV)(φ_y − φ_x)² dV, with a₄ = 4/√6 and g(w) = e^{−w}(−1 + 9w − 8w² + (4/3)w³).
  - Normalisation: its local limit is −½ ∫ φ′² dz = −4m, because ρε ∫ g (φ_y − φ_x)² dV → −2(∂φ)²/(a₄ (ερ)^{1/2}).
- A uniform conformal factor e^b multiplies the density and every interval volume by e^{4b}. Then d s_kin/db =
  (a₄/4) ∫∫ [8g + 4w g′] (Δφ)² dV, and it must tend to the kinetic part of ∫ T^a_a, −∫ φ′² dz = −8m. In 4D
  T^a_a = −(∂φ)² − 4V; the potential part −16m is the measure, exact as in 2D.

**Quadrature.**
- Variables (ρ_s, μ, s) for the past of x, with s = √κ τ².
- A(ρ_s) = ∫ g / (2√κ dt) ds on panels clipped exactly at the age cut, and H = ∫ (Δφ)² dμ in the variable z_y.
- **The age cut needs its own panels in ρ_s.** Near ρ_s = T the clipped kernel gives a large integrand on a width of
  about 12/T. Without those panels the result jumped by 10³. With them it converges, and the leading age term cancels
  (∫ s g(s²) ds = 0).

**Results** (converged in the radial grid and in T):

| wall width 1/m | s_kin / continuum | d s_kin/db / continuum |
|---|---|---|
| 2.5 ℓ_k | 0.943 | 0.982 |
| 5 ℓ_k | 0.979 | 0.995 |
| 10 ℓ_k | 0.993 | 0.998 |

- **Both tend to 1 as ℓ_k m → 0.** The age T changes them by below 10⁻³ once T ≥ 8/m.
- **Widths of 20 and 40 ℓ_k at T ≥ 16/m are not converged** with the present grids (the outer and radial panels).
  These rows are not used.

**Reading.** In 3+1 the pairwise matter action responds to a uniform conformal change of the metric as the continuum
does. This includes the kinetic part, which in 4D is not conformally invariant: −∫ φ′² against −8m, to 0.2 % at a
width of 10 ℓ_k. The age of the order gives no growing term for the wall.

**Stage 2, plane-wave ω = b cos kz** (21 Nov; wall width 1/m = 10 ℓ_k).
- **Method.** The diamond average of cos kz is boosted into the rest frame of each interval: ω′ = kΔz/τ,
  q′ = k√(1 + Δz²/τ²), h = τ/2 (port of diamond_ft4 of F6.8a). At k = 0 the kernel reproduces stage 1 exactly
  (−0.79864).
- **Results** for the kinetic part, against the continuum −∫ cos(kz) φ′² dz = −4πk/sinh(πk/2m). Converged in the
  radial grid (40 → 60: identical); in T within 1 % (T 40 / 80 / 160: 0.881 / 0.890 / 0.889), with no growth with
  age.

  | k ℓ_k | d s_kin/db | continuum | ratio | F6.8a K_loc(k) |
  |---|---|---|---|---|
  | 0.02 | −0.7808 | −0.7870 | 0.992 | 0.990 (k = 0.025) |
  | 0.05 | −0.6981 | −0.7233 | 0.965 | 0.968 |
  | 0.1 | −0.4862 | −0.5461 | 0.890 | 0.899 |
  | 0.2 | −0.1457 | −0.2176 | 0.670 | 0.703 |

- **Reading (development).**
  - The response of the matter's kinetic part to a conformal mode of wavenumber k is the continuum one times a
    nonlocal form factor that tends to 1 as kℓ_k → 0.
  - The form factor is close to the conformal stiffness K_loc(k) of the BD action itself (F6.8a): within 0.3–3 % for
    kℓ_k ≤ 0.1, and 5 % at 0.2.
  - If the source and the stiffness carried the same form factor, it would cancel in the linear response of OV3.
    This is an observation to test, not a result.
  - The potential part (the measure) is local and carries no form factor.

**Stage 3, the cone sector ν = b cos kz** (21 Nov, night; diag/ov2_t2nu.py; wall width 10 ℓ_k).
- **Method.** g = −dt² + e^{2bh} dx². The first-order cones and the interval volume are those of F6.8b2: the optical
  length of the straight segment, the past parametrised by s at a fixed cone, the measure 1/(2 l₂ dt) and the age
  cut depending on b, and dV/db from the prolate quadrature spheroid_pieces.
  - The wall enters through Φ(Δ) = ∫ dz e^{ikz} (φ(z + Δ) − φ(z))².
  - The boundary term of the age cut is integrated in S = s_max as F6.8b2 does.
- **Target:** the kinetic part of ½ ∫ T^{ab} δg_ab, which is −½ ∫ cos(kz) φ′² dz = −2πk / sinh(πk/2m).
- **Results** (T = 80):

  | k ℓ_k | d s_kin/db (ν) | continuum | ratio ν | ratio ω (stage 2) |
  |---|---|---|---|---|
  | 0 | −0.40181 | −0.40000 | 1.0045 | 0.998 |
  | 0.02 | −0.39202 | −0.39350 | 0.996 | 0.992 |
  | 0.05 | −0.34745 | −0.36166 | 0.961 | 0.965 |
  | 0.1 | −0.23554 | −0.27303 | 0.863 | 0.890 |
  | 0.2 | −0.06242 | −0.10881 | 0.574 | 0.670 |

  - Age: T = 40 / 80 / 160 give 0.855 / 0.863 / 0.861 at k = 0.1.
  - Grid: mult 1 against 1.5 is identical to 5 digits.
  - The boundary term of the age cut is ≤ 1.4·10⁻³ and falls with T.
- **Reading (development).**
  - The response of the wall's kinetic action to a change of the cones, which is a change of the order itself,
    reproduces the spatial part of T_ab in the local limit: 0.45 % at k = 0 and 0.4 % at kℓ_k = 0.02.
  - It carries a nonlocal form factor that falls somewhat faster than that of the conformal sector at larger k.
  - Together with stage 1 and stage 2 and the exact measure, **T2 holds at the ensemble level**: the pairwise matter
    action responds to the static metric perturbations of both sectors as ½ ∫ T^{ab} δg_ab, up to form factors
    → 1 as kℓ_k → 0.
  - The draft card is docs/OV2_coupling_cards_draft.md (T2, proposed F6.10b).

**Next stages:**
- **Stage 2 (done above):** plane-wave ω = b cos kz, using the diamond transform of F6.8a (diamond_ft4, boosted to the frame of each
  interval). It must give −∫ cos(kz)(φ′² + 4V) dz = −12πk / sinh(πk/2m).
- **Stage 3:** the cone sector ν with the first-order cones of F6.8b2. Expected: −∫ ν(½φ′² + 3V) dz for the static
  wall.

