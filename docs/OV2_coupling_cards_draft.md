# Cards OV2-T1 (2D) and OV2-T2 (3+1): the response of the matter action to the metric (FROZEN)

**Revision 3, FROZEN by the author, 21 Nov 2026 09:17 (decision 40).** F6.10a (T1) and F6.10b (T2); T2 is one card
for both sectors; the sprinkled cross-check of T1 (I4) is not required. Tolerances and fresh parameters as below.
Test files (md5): phases/F6_order_dynamics/F6_10a_matter_conformal_2d.py a53d3a3188eefef7173b80f383615223,
phases/F6_order_dynamics/F6_10b_matter_metric_3p1.py bdae39b0e674bc190f0fb64eb4229254.
- Revision 2 (21 Nov, night): drafts prepared by Claude at the author's request.
- Rev 2: both cards implemented as draft test files (quick runs only):
  phases/F6_order_dynamics/F6_10a_matter_conformal_2d.py and phases/F6_order_dynamics/F6_10b_matter_metric_3p1.py.
  The grid levels of T1 are fixed from a convergence check; the Δ of T2 is defined; the sprinkled cross-check of T1
  (I4) is not in the file, pending decision 3.
- Specification: docs/derivations/order_variation_OV2.md, sections 2, 3.1, 4, 9, 10 and 12, and decision 39.
- Development: diag/ov2_2d.py, diag/ov2_quad2d.py, diag/ov2_t2.py, diag/ov2_t2nu.py.

Both cards are deterministic quadratures at the ensemble level (L1 of OV1 §5). The tolerances follow from the
development rows and from the convergence of the grids. Every card uses fresh model parameters that the development
did not use: a different wall width, bump width and position, and wavenumbers.

## Common ground

- **The matter action** (OV2 §3.1, decision 39):
  - S_m = S_kin + S_pot, with the pairwise (Dirichlet) kinetic action S_kin = −c_d Σ_pairs f(n, ε)(φ_x − φ_y)²
    (c_2 = ε², c_4 = a₄ (ερ)^{1/2} ε/(4ρ)) and S_pot = −(1/ρ) Σ V(φ_v).
  - It is evaluated at a fixed field configuration, a static sine-Gordon kink (wall); T_ab is defined off shell.
- **The expectation** over Poisson sprinklings, by the Mecke formula. The kernel average is exactly g(ερV), so
  ⟨S_kin⟩ is a deterministic double integral over pairs.
- **What is judged.** The first-order response of ⟨S_m⟩ to a metric perturbation against the continuum
  ½ ∫ T^{ab} δg_ab √|g|.
- **What is not judged:**
  - the potential (measure) part, which is exact by construction and serves as a check of the normalisation only;
  - the absolute value of ⟨S_kin⟩, which is an instrument item.

## T1 (2D, conformal; proposed F6.10a)

**Hypothesis.** For g → e^{2σ} g in 2D, the first-order response of ⟨S_m⟩ is ∫ σ T^a_a √|g| with T^a_a = −2V. The
kinetic part is conformally invariant up to a nonlocal correction that is small once the kink width exceeds ℓ_k.

**Model.**
- The light-cone square, ε = 0.12, N_exp = 66 560 and 266 240 (T/w = 32 and 64 of OV1).
- Kink m = 9 (development: 12), centred on the diagonal.
- Bump σ = a exp(−((u − u₀)² + (v − v₀)²)/(2s²)) with s = 0.09 and (u₀, v₀) = (0.45, 0.55) (development: 0.12 and
  (0.5, 0.5)).
- The first order is read from the analytic derivative under the integral (diag/ov2_quad2d.py, kin_quad_d1),
  generalised to (u₀, v₀).
- **Grid levels** (n_outer, n_t, n_w) = (192, 128, 96) and (256, 160, 128), Δ_quad the change between them. Chosen
  from a convergence check at the fresh parameters (R at levels 64/96/128/192/256 of n_outer: +0.551, −0.006, −0.0139,
  −0.0139, −0.0139 at N 66 560; −8.49, −0.341, −0.0048, −0.0047, −0.0047 at N 266 240). The coarse levels are far off
  at the higher density, so the judged pair is one level above the first converged one.
- **Disclosure.** The grid check computed R at the fresh parameters before the card was frozen. The quadrature is
  deterministic and the 3 % criterion was written in revision 1 before this, so nothing was tuned on it; the numbers
  are R = −0.014 (N 66 560) and −0.005 (N 266 240), I1 0.9966.

**Criterion** (per density).
- **PASSED** if |δ₁⟨S_m⟩/∫σ T^a_a √g − 1| + Δ_quad ≤ 0.03 and |δ₁⟨S_kin⟩/∫σ T^a_a √g| + Δ_quad ≤ 0.03.
  - Δ_quad is the change between the two finer grids.
  - In development the kinetic share was 1.5 % at T/w 32 and 0.5 % at 64, with m = 12.
- **FAILED** if either quantity minus Δ_quad exceeds 0.06.
- **UNDECIDED** otherwise.

**Instrument** (INSTRUMENT_FAILED if any fails):
- **I1:** the flat ⟨S_kin⟩ against the continuum ∫ φ_u φ_v du dv within 1 % at T/w 64. In development −42.26
  against −42.45.
- **I2:** grid convergence, the two finer grids within 1 % of the signal.
- **I3:** at a = 0 the derivative of the potential part equals −∫ σ̂ V du dv to 1e−6 (normalisation).
- **I4 (not required, decision 40; not part of the test):** the cross-check against a sprinkled estimate at T/w 32. The sprinkled mean of the first order (1 024
  sprinklings, as in development) must lie within 3 SE of the quadrature.

## T2 (3+1, the conformal and the cone sector; proposed F6.10b)

**Hypothesis.** For the static 3+1 metric g = e^{2ω}(−dt² + e^{2ν} dx²) with ω = a cos kz and ν = b cos kz, the
first-order response of ⟨S_kin⟩ of a wall equals the continuum response times a nonlocal form factor that tends to 1
as kℓ_k → 0.
- ω: −∫ cos(kz) φ′² dz.
- ν: −½ ∫ cos(kz) φ′² dz.
- In 4D the kinetic part is not conformally invariant, so ω tests it through the density and the volumes, and ν
  through the order (the cones).

**Model.**
- ερ = 1; wall widths 1/m = 8 ℓ_k (development: 10).
- k ℓ_k = 0.03, 0.06, 0.12 (development: 0.02, 0.05, 0.1, 0.2).
- Ages T = 120 (judged) and 80 (I1); grid multiplier 1 (judged) and 1.5 (I2).
- Δ = |ratio(T 120) − ratio(T 80)| + |ratio(grid 1.5) − ratio(grid 1)|.
- The quadrature code is copied verbatim into the test file from diag/ov2_kernels4.py, diag/ov2_t2.py and
  diag/ov2_t2nu.py, so the frozen test does not depend on the development files.
- Development row seen at the fresh width while timing (k = 0.06, T = 80): ω 0.954, ν 0.949.
- ω by the boosted diamond average (diag/ov2_t2.py); ν by the first-order cones and the prolate volume quadrature of
  F6.8b2 (diag/ov2_t2nu.py).

**Criterion.**
- **The local limit:**
  - **PASSED** if, at k → 0 (the k = 0 rows), both sectors give |ratio − 1| + Δ ≤ 0.02, and if the ratios at the
    three k do not increase with k and |ratio(k) − 1| + Δ ≤ 0.15 at kℓ_k ≤ 0.06.
  - **FAILED** if a k = 0 ratio is off by more than 0.05 + Δ, or a ratio at kℓ_k ≤ 0.06 by more than 0.25 + Δ.
  - **UNDECIDED** otherwise.
- **Reported, not judged:** the form factors against K_loc(k) of F6.8a (ω) and the mixed K_loc,ab of F6.8b2 (ν).
  Whether source and stiffness share their form factor is an observation for OV3.

**Instrument:**
- **I1:** age convergence, rows at T = 80 and 120 within 1 %.
- **I2:** grid convergence (mult 1 and 1.5) within 1 %.
- **I3:** the ω sector at k = 0 reproduces the uniform conformal value of stage 1 to 1e−5.
- **I4:** the geometry of F6.8b2 (its I1 check: dV/db of the test diamond and the uniform null) passes as there.

## Decisions of the author (decision 40, 21 Nov 09:17)

1. The IDs F6.10a and F6.10b, the tolerances and the fresh parameters above: accepted.
2. T2 is one card (both sectors).
3. The sprinkled cross-check (T1 I4) is not required.
