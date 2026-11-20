# Cards OV2-T1 (2D) and OV2-T2 (3+1): the response of the matter action to the metric (drafts)

**Revision 1, draft, 21 Nov 2026 (night; the author asked Claude to continue and prepare tests). Not frozen.** IDs are
the author's choice; proposed F6.10a (T1) and F6.10b (T2).
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
- **I4:** the cross-check against a sprinkled estimate at T/w 32. The sprinkled mean of the first order (1 024
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
- Ages T = 80 and 120.
- ω by the boosted diamond average (diag/ov2_t2.py); ν by the first-order cones and the prolate volume quadrature of
  F6.8b2 (diag/ov2_t2nu.py).

**Criterion.**
- **The local limit:**
  - **PASSED** if, at k → 0 (the k = 0 rows), both sectors give |ratio − 1| + Δ ≤ 0.02, and if the ratios at the
    three k decrease monotonically with |ratio(k) − 1| ≤ 0.15 at kℓ_k ≤ 0.06.
  - **FAILED** if a k = 0 ratio is off by more than 0.05, or the ratios at kℓ_k ≤ 0.06 are off by more than 0.25.
  - **UNDECIDED** otherwise.
- **Reported, not judged:** the form factors against K_loc(k) of F6.8a (ω) and the mixed K_loc,ab of F6.8b2 (ν).
  Whether source and stiffness share their form factor is an observation for OV3.

**Instrument:**
- **I1:** age convergence, rows at T = 80 and 120 within 1 %.
- **I2:** grid convergence (mult 1 and 1.5) within 1 %.
- **I3:** the ω sector at k = 0 reproduces the uniform conformal value of stage 1 to 1e−5.
- **I4:** the geometry of F6.8b2 (its I1 check: dV/db of the test diamond and the uniform null) passes as there.

## Decisions for the author

1. The IDs, the tolerances and the fresh parameters above.
2. Whether T2 is one card (both sectors) or two.
3. Whether the sprinkled cross-check (T1 I4) is required or only reported. It costs about 1 h on 8 workers.
