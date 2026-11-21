# OV3: the response of the order to matter — the linear step (draft, development)

21 Nov 2026 (night; the author asked Claude to continue and prepare tests). Status: a draft. It decides nothing, no
criterion is frozen, and the open decisions in section 7 are the author's.
- It joins OV1 (the variations of the order; the vacuum control F6.9a PASSED) and OV2 (the coupling of matter: T1
  and T2 at the ensemble level, development).
- It derives what the linear static response of the order to a source must be when both ingredients are taken at
  their measured values. It is the step "M → H_M" of the author's direction (19 Nov), restricted to linear order,
  family M and the ensemble level (L1).

## 1. The linear static response in family M (L1)

**Variables.** The static metric of the order, g = e^{2ω}(−dt² + e^{2ν} dx²) (the variables of F6.8: ω conformal, ν
optical). For a static source varying along z, mode by mode in k, the second variation of the expected BD action
(F6.8) and the first variation of the expected matter action (OV2 T2) give the stationarity conditions (L1)

    α k² [ 6 K(k) ω + 4 K_ab(k) ν ] = F_ω(k) τ_k,       α k² [ 4 K_ab(k) ω + 2 K_c(k) ν ] = F_ν(k) σ_k,

with:
- τ = T^a_a and σ = Σ_i T_ii, the sources (OV2 §2 and §12);
- K, K_ab and K_c, the form factors of the BD stiffness, from F6.8a, b2 and c2 (all → 1 as kℓ_k → 0);
- F_ω and F_ν, the form factors of the matter coupling, from OV2 §12 (the kinetic part; the potential part has
  F = 1);
- α the overall stiffness, set by the normalisation of the BD action and ℓ_k, which fixes G. Its laboratory value is
  not derived here.

**The local limit** (all form factors 1; the matrix [[6, 4], [4, 2]], determinant −4; solved and checked with sympy)
gives

    ω = (2σ − τ)/(2 α k²) = (ρ + Σp)/(2 α k²),      ν = (2τ − 3σ)/(2 α k²),     Ψ = −(ω + ν) = ρ/(2 α k²),

where ρ = T_tt and Σp = Σ_i T_ii. With α = −1/(8πG) these are the static linearised Einstein equations:
- the lapse, ∇²Φ = 4πG(ρ + Σp), the weak-field Whittaker equation of substrate_and_action_audit.md §7.1;
- the spatial part, ∇²Ψ = 4πG ρ.

This is the algebra of the EH form, not a new result. What is new is that both ingredients now come from the order:
- **the stiffness.** F6.8a, b2 and c2: 6, 8 and 2 from counting intervals;
- **the source.** OV2 T2: the response of the pairwise matter action to ω and to ν (a change of the cones)
  reproduces T^a_a and Σ T_ii in the local limit (k = 0: 0.998 and 1.0045, wall width 10 ℓ_k).

## 2. A known answer for walls: repulsion

For the static sine-Gordon wall (BPS: ½φ′² = V):
- ρ = T_tt = ρ_E(z), and Σp = (½φ′² − V) − 2(½φ′² + V) = −2ρ_E;
- the tension equals the energy, the passive weight of F5.2c is E − σ = 0, and the kink stays at rest.

So the source of the lapse is ρ + Σp = −ρ_E < 0: **a domain wall is gravitationally repulsive** in linearised GR (the
Vilenkin and Ipser–Sikivie result for walls, quoted from the literature, not re-derived here beyond the algebra
above). The spatial potential Ψ is sourced by +ρ_E.

This makes the wall a sharp known-answer case for OV3:
- **sign.** The response of the order to a wall must give the lapse of opposite sign to that of a dust source of the
  same energy;
- **ratio.** The ratio of the two potentials of the order, Ψ/Φ = ρ/(ρ + Σp), is −1 for the wall and +1 for dust
  (γ = 1).

## 3. What is known about the form factors (development numbers)

| kℓ_k | K (F6.8a, K_loc) | F_ω (OV2 stage 2) | F_ν (OV2 stage 3) |
|---|---|---|---|
| 0.02–0.025 | 0.990 | 0.992 | 0.996 |
| 0.05 | 0.968 | 0.965 | 0.961 |
| 0.1 | 0.899 | 0.890 | 0.863 |
| 0.2 | 0.703 | 0.670 | 0.574 |

- **K_ab.** F6.8b2: K_ab − K = Δ(k), with |Δ| ≤ 0.01 for kℓ_k ≤ 0.1 (PASSED).
- **K_c, the optical sector.** F6.8c2 gives 𝒞(0) = 2.002. Its k-dependence at fixed kT, together with the age term
  (kT)⁵ of F6.8e3, is known only in the ranges of those tests.
- **The coupling's form factors are close to the stiffness's** for kℓ_k ≤ 0.05. If they were equal sector by sector,
  they would cancel in ω and ν for any source, and the response of the order would be local even at finite ℓ_k. They
  are not equal at kℓ_k = 0.1–0.2 (F_ν falls faster), so the response there departs from GR at order (kℓ_k)². It is
  computable from the table once K_c(k) is known.

## 4. What blocks the static response of a real (old) source

- **The age of the field.**
  - K(k, T) contains the age terms: c_T k³ T in the conformal sector, (kT)³ in the mixed one (F6.8d undecided) and
    (kT)⁵ in the optical one (F6.8e3).
  - F6.8f FAILED its finite-amplitude reading for old fields (b²(kT)⁴ ≫ 1).
  - The table above is the local part; using it is the "young field" reading (option 4 of
    substrate_and_action_audit.md §7.3), and it must be declared as such.
  - The matter side has no growing age term (OV2 §12: the wall converges in T within 1 %).
- **The selection principle.**
  - Stationarity of ⟨S⟩ (L1) is used here.
  - OV1's phase sampler showed that S_BD alone does not keep the flat order: the manifold-like phase is held by
    entropy below β_e(N).
  - Whether the physical selection is the L1 saddle, an action-and-entropy balance or a stochastic dynamics is the
    author's decision (OV2 §7). At linear order around the flat order all three may agree; beyond it they need not.
- **The readout 𝒫.** Π = 𝒫[H] is read as the ratio of clocks (A2-B). At linear order its lapse part is e^{ω} at
  fixed coordinates, but 𝒫 as an operation on the order (the clocks of K.16) is still to be fixed. The response
  above is in the metric variables of family M, not yet in Π.

## 5. Proposed development steps (decide nothing)

1. **OV3-lin-a (cheap, deterministic).** Solve the 2 × 2 system per k with the measured form factors for kℓ_k ≤ 0.1,
   for a dust slab and for the wall. Report Φ(z) and Ψ(z), their ratio, and the departure from linearised GR. This
   needs K_c(k) at the same k: a run of the F6.8c2 machinery at kT fixed, or the covariance relation
   C = 2 + 4φ₀ + φ₂ with φ₂ measured.
2. **OV3-lin-b (the sign test).** The wall against a dust slab of the same energy: the order's lapse response must be
   repulsive for the wall. As a card it is a sharp, sign-level known answer.
3. **OV3-lin-c (L2, individual orders).** Whether single sprinklings perturbed along the response direction show the
   expected change of the action (S-b). The noise of S_kin (OV2 §9) makes this costly, so it is optional.

## 6. Relation to the roadmap

- This is C.1 at linear order for young or weak fields, derived from the order (A2-B: one degree of freedom, H; Π a
  readout), with the source from a matter action that depends only on the order.
- It is not the static propagator of a real old body. The age problem (F6.8f) stands.

## 7. Open decisions (the author's)

1. Whether OV3-lin-a is computed now with the young-field reading declared, or only after the age question is
   settled.
2. The source of K_c(k): a new F6.8c-type run at the needed k, or the covariance relation.
3. Whether the wall-repulsion sign test becomes a card.
4. The selection principle (OV2 §7), before any step beyond linear order.

## 8. Decisions (the author, 21 Nov 11:01; decision 41)

1. **OV3-lin-a is computed now, in the young-field reading.**
   - The output is a conditional linear response, not the propagator of an old static body.
   - The card and the findings state:
     - that the local part of the stiffness is used;
     - that validity refers to a young field, or to a field weak enough for the age term to be negligible;
     - that the finite-amplitude and age problem of F6.8f stays open;
     - that the result by itself closes neither C.1 nor the 1/r law.
   - **The bare action.** It is a result for the bare ensemble action ⟨S⟩ until the entropic correction to the second
     variation (item 4) is evaluated. If that correction is significant, the result is not the physical linear
     propagator without it.
2. **K_c(k): option A.** A new F6.8c-type computation at kℓ_k = 0.02–0.12, the local part extracted as for K in F6.8a.
   - The card reports separately:
     - the total measured optical stiffness;
     - the age component;
     - the local part;
     - the sensitivity of the separation to the chosen ranges of k and T.
   - The covariance relation C = 2 + 4φ₀ + φ₂ is an independent check. φ₂ is measured separately, never computed from
     the C it checks.
3. **The BPS wall: a separate card.**
   - Three preregistered evaluations:
     - the sign of the response;
     - the deviation from −ρ_E;
     - the ratio of the wall to dust at the same k.
   - **The error of the source.** The lapse source of the wall is −4F_ν + 3F_ω (in units of ρ_E, at equal stiffness
     factors). Its uncertainty must include the covariance of the two form factors:
     σ² = 16σ_ν² + 9σ_ω² − 24 Cov(F_ν, F_ω).
     - In the deterministic quadratures of F6.10b the uncertainties are systematic (age, grid), not statistical.
     - The cleanest way to get the covariance is to compute the combination itself in one quadrature, so that its age
       and grid changes are measured directly. This is proposed for the card.
4. **The selection principle.**
   - The stationary point of an individual action is rejected: OV1 found locally reachable orders with an action lower
     by 10³, and the manifold-like phase is held by the statistical weight.
   - **Linear OV3** uses the saddle of the ensemble mean ⟨S⟩, conditioned on the manifold-like phase (option 1).
   - **Nonlinear OV3** uses the explicit balance of action and entropy (option 2), through an effective functional of a
     macrostate, Γ_β(M) = −(1/β) ln Σ_{H ∈ M} e^{−β S_tot[H]}. Its exact form and measure belong to a separate
     specification.
   - A stochastic dynamics (option 3) is not chosen. It would need transition probabilities and a time dynamics of
     the configuration space.
   - **Measured before any physical reading:** whether the entropic term changes the effective stiffness already at
     quadratic order.

## 9. The field equation of the order, by order in the field (21 Nov, evening; development, decides nothing)

The field equation of the theory is the stationarity of the effective functional of the order (decision 41):

    δΓ[g] / δg_ab = 0,     Γ = ⟨S_BD⟩ + ⟨S_m⟩  (+ the entropic term of Γ_β beyond linear order).

Nothing of GR is put in. ⟨S_BD⟩ counts intervals; ⟨S_m⟩ is the pairwise matter action, which depends only on the order
(OV2). Below, "order n" means order n in the field amplitude b of the equation (order n + 1 of the action).

### 9.1 First order: established for young fields

**The gravity side.** By covariance (docs/derivations/path_to_predictions.md §2) the quadratic part of ⟨S_BD⟩ is
∫√−g [R + R F₀(□) R + C F₂(□) C]. Its three static sectors are

    A = 6 + 9φ₀,   B = 8 + 12φ₀,   C = 2 + 4φ₀ + φ₂,   with φ₀ = 4k²F₀ and φ₂ = (4/3)k²F₂.

F6.8g (PASSED, n = kT = 4) measured all three on the same rows. φ₀ from the conformal sector and from the mixed sector
agree to 1e-4: this is a real check of covariance (one form factor, two sectors).

| kℓ_k | φ₀ (n = 3 / 4 / 6) | φ₂ (n = 3 / 4 / 6), from C |
|---|---|---|
| 0.045 | −0.0196 / −0.0194 / −0.0188 | +0.0030 / +0.0000 / −0.0010 |
| 0.06 | −0.0322 / −0.0317 / −0.0307 | +0.0052 / +0.0007 / −0.0027 |
| 0.08 | −0.0523 / −0.0515 / −0.0497 | +0.0092 / +0.0013 / −0.0050 |
| 0.10 | −0.0759 / −0.0746 / −0.0718 | +0.0142 / +0.0019 / −0.0081 |
| 0.12 | −0.1024 / −0.1005 / −0.0965 | +0.0203 / +0.0025 / −0.0118 |

- **The spin-0 channel carries the departure from EH.** At n = 4, φ₀ = −k²ℓ_k²(2.63 ln(1/kℓ_k) + 1.41) (max residual
  1.2e-5), i.e. F₀ = −ℓ_k²(0.66 ln(1/kℓ_k) + 0.35). It is nearly independent of the age (±2 % between n = 3 and 6).
- **The spin-2 (Weyl) channel is small at young ages but not robust.** |φ₂| ≤ 0.02 at kℓ_k ≤ 0.12. It changes sign
  between n = 4 and 6, so φ₂ ≈ 0 at n = 4 is a crossing, not a property. φ₂ is computed from C, the same number it
  would check, so it is not an independent measurement (decision 43); the TT-shear sector (F6.8h) would measure it.

**The matter side.** The first variation of ⟨S_m⟩ is T_ab with two form factors (F6.10b, PASSED): F_ω on the trace and
F_ν on Σ T_ii. Both → 1 as kℓ_k → 0; at kℓ_k = 0.06 / 0.12 they are 0.954 / 0.855 (F_ω) and 0.949 / 0.818 (F_ν). They
do not depend on the age (T 80 against 120: < 0.1 %).

**The linear field equation** (young field, bare ensemble action; F₂ dropped at n = 4):

    G_ab⁽¹⁾ + 2 (η_ab □ − ∂_a ∂_b) F₀(□) R⁽¹⁾  =  8πG T_ab^eff,

where T_ab^eff carries F_ω and F_ν, and α = −1/(8πG) is the overall stiffness (G is not derived).
- As kℓ_k → 0 it is the linearised Einstein equation, with γ = 1 and the source ρ + Σp for the lapse.
- At finite kℓ_k it gives OV3-lin-a (findings, derived consequence, 21 Nov): Φ/Φ_GR = 0.940 / 0.809 and γ(k) = 1.033 /
  1.122 at kℓ_k = 0.06 / 0.12; a BPS wall repels (wall/dust −1.045 / −1.062).
- It belongs to the known class of nonlocal quadratic gravity. What is new is that the form factor F₀ and the source
  form factors come from counting intervals, not from a choice.

### 9.2 Second order (the cubic action): not computed

- **What it decides in GR:** the PPN parameter β (the perihelion of Mercury), the gravitational self-energy, and the
  equivalence principle for bound systems.
- **The continuum expectation.** The mean BD action on a curved manifold tends to the EH action in the continuum limit
  (Benincasa and Dowker 2010; the paper relies on this). In the local limit the full nonlinear equation should
  therefore be G_ab = 8πG T_ab, so β = 1. The laboratory has not checked this beyond first order.
- **What the order adds: the expansion parameter is ε = b (kT)², not b.**
  - F6.8f (FAILED, reading B) measured the optical sector at finite amplitude. The response relative to second order
    is R = 1.003, 1.056, 1.25, 2.52 at ε = 0.25, 1, 2, 4, i.e. R − 1 ≈ 0.055 ε² at small ε. The expansion fails at
    ε ~ 1–4 even for a weak field.
  - Caveat: F6.8f measured in the asymptotic (old-field) regime, not at n = 4.
  - For the young field of 9.1 (n = 4), ε = 16 b. In the solar system (b ~ 1e-6) this is negligible.
  - For old fields ε grows as (kT)². The age problem and the nonlinearity problem of the order are one problem.
- **The next computation:** the cubic order of ⟨S_BD⟩ and ⟨S_m⟩ for a static source at fixed n = 4. It needs the
  method of F6.8g extended to third-order interval volumes. It gives β, to be compared with the perihelion of Mercury
  (planned F8.7).

### 9.3 Higher orders and the strong field

Only indirect: the finite-amplitude section of the paper (the first conjugate point at ε = π², where the field focuses
light within intervals). The selection principle beyond linear order is the balance of action and entropy, Γ_β
(decision 41). Whether its entropic term cures the age growth is open; it is the main open question of the field
equation.

### 9.4 Status, in one table

| Part of the field equation | Status |
|---|---|
| Linear, static, young field, kℓ_k → 0 | Einstein, γ = 1 (F6.8a/b2/c2, F6.8g, F6.10a/b) |
| Linear, static, young field, finite kℓ_k | spin-0 form factor F₀ measured; source form factors measured; Φ, γ(k) and the wall computed (OV3-lin-a) |
| Linear, time-dependent and tensor (TT) sectors | not measured (F6.8h not developed) |
| Second order (β) | not computed; the next derivation |
| Old fields (any order) | open: the age growth (kT)⁵ (F6.8e3), no saturation at finite amplitude (F6.8f) |
| Entropic correction to the stiffness | not measured (decision 41) |
| The value of G | not derived |

## 10. The second order: what decides β, and how the order can be asked (21 Nov, 22:xx; development, decides nothing)

The author chose this direction (21 Nov, 21:32: "Druhý řád polní rovnice (β) při kT=4 … jdeme tímto směrem").
Script: diag/ov3_beta_sym.py (sympy; every statement marked [sympy] below is printed by it).

### 10.1 β needs only one cubic number

Static isotropic metric, written as g = −e^{2ω} dt² + e^{2ψ} dx² (the paper's ν is ψ − ω). The most general **local**
static action with two derivatives has, to cubic order,

    L = 4 ∇ω·∇ψ + 2 (∇ψ)²                                    (the quadratic part at its EH values, 6, 8, 2 in ω, ν)
      + c₁ ω(∇ω)² + c₂ ω ∇ω·∇ψ + c₃ ω(∇ψ)² + c₄ ψ(∇ω)² + c₅ ψ ∇ω·∇ψ + c₆ ψ(∇ψ)²,

since every two-derivative cubic term reduces to these by integration by parts. For a point mass the first order is
ω₁ = −U, ψ₁ = +U with U = m/r (ω + ψ has no source and is regular). Outside the source the second order is
ω₂ = a m²/r², ψ₂ = b m²/r², and with −g₀₀ = e^{2ω} = 1 − 2U + 2βU²:

    β = 1 + (c₁ − c₂ + c₃ − c₄ + c₅ − c₆)/8   [sympy]

- **Check.** EH, L = e^{ω+ψ}[4∇ω·∇ψ + 2(∇ψ)²], has c = (0, 4, 2, 0, 4, 2). It gives β = 1 and b = −1/4, which is
  Schwarzschild in isotropic coordinates, ψ = 2 ln(1 + m/2r) = m/r − m²/4r² [sympy].
- **The reduction.** The combination is exactly the cubic action evaluated along the first-order solution:
  κ₃ ≡ L₃(ω = −U, ψ = U)/(U(∇U)²) = −c₁ + c₂ − c₃ + c₄ − c₅ + c₆, so **β − 1 = −κ₃/8** [sympy]. In general
  normalisation, β − 1 = (C₃/Q₂)/4, where Q₂ and C₃ are the coefficients of (∇U)² and U(∇U)² along that direction.
- **The direction is the exponential metric.** (ω, ν) = (−U, 2U), i.e. g = −e^{−2U} dt² + e^{2U} dx². Along it the full
  EH action is −2(∇U)² at **every** order [sympy], because the prefactor e^{ω+ψ} = 1.

So the question "is β = 1?" becomes a null test: **does the cubic part of the order's action vanish along the
exponential metric?** No other cubic coefficient enters. In particular the six coefficients need not be separated.

### 10.2 What the order must supply

The quantity is E(A, C) = ⟨S_BD⟩ per unit four-volume for the static planar metric

    g = −e^{−2u} dt² + e^{2u} dx²,    u(z) = A cos kz + C cos 2kz,

at fixed age n = kT = 4 (the young field of F6.8g). The averages are ⟨u(u′)²⟩ = (3/4)A²Ck² and
⟨(u′)²⟩ = (A² + 4C²)k²/2 [sympy], so

    E = −K̂_dir (A² + 4C²) k² + κ₃(k) (3/4) A² C k² + O(4),
    β(k) − 1 = −κ₃(k) / (8 K̂_dir(k)),

with EH values K̂_dir = 1 and κ₃ = 0.
- **The quadratic part is known.** K̂_dir = (16 K̂_ab − 6 K̂_c − 8 K̂_bb)/2 from the three sectors of F6.8g. At
  kℓ_k = 0.06 / 0.12 it is 1.014 / 1.044 (computed from the F6.8g rows; above 1, because the optical sector falls
  fastest).
- **The causal structure** of g depends only on ν = 2u (the cones). The density of events is e^{4ω+3ν} = e^{2u}.
- **An exact null.** A single mode (C = 0) has no cubic term: a half-period shift maps A → −A. The A³ coefficient of
  the order must therefore vanish identically. This is a free instrument check.
- **Collinear configuration.** The planar modes (k, k, −2k) probe the cubic vertex at one configuration of momenta.
  In the local limit (kℓ_k → 0) that is all β needs. At finite kℓ_k, κ₃(k) is the collinear form factor, not the β of
  a point mass at that scale.

### 10.3 How to compute it: two routes

Both need what F6.8 never needed: interval volumes to **third** order in a perturbation that bends the cones, with
two modes.
- **(a) Perturbative.** Extend the Fermat expansion of the optical sector (F6.8c: T₂) to T₃, and the Mecke/kernel
  expansion of the expected action (F6.8a/b2) to third order. Exact in amplitude order by order, so there are no
  finite-difference cancellations. This is a new derivation: T₃ for a stratified medium, the kernel g''' terms, and
  the age cut at third order. It has an exact check: at C = 0 every A³ term must cancel.
- **(b) Finite amplitude.** Use brute-force interval volumes in the two-mode medium (generalise the needle volumes of
  K.18b / diag/f68f_1d.py from one mode to the profile u(z)), then a third finite difference in (A, C). At n = 4,
  ε = b(kT)² = 32A, so A ≲ 1e-3 keeps the quartic terms below the cubic. A third difference at A ~ 1e-3 needs
  volumes to ~1e-10 relative. K.18b reached 1.1e-6, so (b) needs a precision this code has not shown.

**Recommendation: route (a), in three stages.**
- **β1, calibration (cheap).** The pure conformal sector at cubic order: the order is unchanged and only the density
  e^{4ω} changes. EH gives 12 ω(∇ω)², i.e. 9A²Ck² for the two modes. This is F6.8a extended to third order.
  It checks the Mecke and kernel expansion and the age cut at third order, and gives no β.
- **β2, the cones at third order (the derivation).** T₃ of rays in a stratified medium. Check it against brute-force
  volumes at small amplitude (K.18b code, one and two modes).
- **β3, the test.** κ₃(k) along the exponential metric at n = 4, kℓ_k = 0.045–0.12, fitted as in F6.8g; β(k → 0).

### 10.4 What a β test can and cannot say (for the author)

- **The young-field local limit is expected to give β = 1.** ⟨S_BD⟩ tends to the EH action in the continuum limit of a
  fixed geometry (Benincasa–Dowker 2010), and k → 0 at fixed n is such a limit. The test checks this at cubic order,
  as F6.8g did at quadratic order, and measures the finite-kℓ_k departure.
- **For Mercury it changes nothing directly.** The departure is of order (kℓ_k)² ln(1/kℓ_k), with ℓ_k ≤ 1e-19 m against
  solar-system scales.
- **The β that Mercury actually measures belongs to an old field:** the Sun's field has existed for Gyr. There
  ε = b(kT)² is huge, and F6.8f shows that the expansion in ε does not saturate. **The decisive β question is
  therefore the age problem (F6.8i), not the young cubic vertex.** The young β test is still the necessary first
  step: without it, the cubic order of the order's action is unchecked.

### 10.5 Proposed card (draft, not frozen): F6.12

See docs/F6_12_card_draft.md.

### 10.6 Stage β1, first results (21 Nov, 22:4x; development, decides nothing)

diag/ov3_beta1.py: the conformal sector at cubic order, u = A cos kz + C cos 2kz, F6.8a's quadrature, kernel and diamond
transform imported unchanged; level 1 (mult 1), Claude's container.
- **The code reproduces the quadratic order.** K2_A at (0.12, n 4) = 0.849174, against F6.8a's plane_wave_row 0.849175
  and F6.8g's K_c 0.84918; at (0.06, n 4) 0.952443 against 0.95245.
- **The phase-average nulls hold.** The A³ and C³ coefficients are ≤ 1e-11.
- **K3_c (the A²C coefficient relative to EH, 4.5 k²) at n = 4:**

  | kℓ_k | 0.03 | 0.045 | 0.06 | 0.08 | 0.10 | 0.12 |
  |---|---|---|---|---|---|---|
  | K3_c | 1.00480 | 1.01073 | 1.01888 | 1.03305 | 1.05073 | 1.07160 |

  - The fit K0 + k²(c₁ ln(1/k) + c₀) gives K0 = 0.99952 over k = 0.03–0.12 (max residual 9e-5) and 0.99926 over
    0.045–0.12 (6e-5).
  - **The conformal cubic of the discrete action tends to the Einstein–Hilbert value.** At finite kℓ_k it lies
    above it, with the opposite sign of the quadratic departure.
- **The age dependence is not monotonic.** K3_c at n = 3 / 4 / 6 is 1.0240 / 1.0189 / 1.0251 (k 0.06) and
  1.0922 / 1.0716 / 1.0959 (k 0.12). The quadratic K2_A is monotonic (0.9517 / 0.9524 / 0.9539 at k 0.06). The cubic
  involves the mode 2k, whose own age is 2n. Whether this is physics or quadrature is to be checked at level 2 (running).
