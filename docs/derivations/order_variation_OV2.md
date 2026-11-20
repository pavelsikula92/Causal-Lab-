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
- **T/w 32** (a = 0.2, 1 024 sprinklings) runs on the author's machine.

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
