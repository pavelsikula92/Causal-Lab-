# The path to predictions: the static propagator of the order

Written on 28 Oct 2026 after F6.8a. The computations were checked in sympy; the script is reproduced below. It follows
`rule_from_action.md`, section 8.

## 1. What a prediction of the theory is

At linear order and long wavelengths the theory is general relativity. This follows from decision 13, from the
derivation of the substrate rule and from F6.8a. Its own prediction, distinguishable from general relativity, can
therefore come only from the discreteness of the order. Concretely: the shape of the deviation of the static field of a
point mass from Newton and from γ = 1 at scales close to ℓ_k, and the age term.

Such a prediction has **a fixed shape and a single free parameter ℓ_k**:

- the shape is fixed by the coefficients of the BD kernel, not by tuning;
- experiments (torsion balances below 1 mm, γ from Shapiro/Cassini) then either bound ℓ_k or see the deviation.

## 2. General structure (covariance)

The mean BD action is a functional of the geometry only, so it is diffeomorphism invariant. Its quadratic part around
flat space must therefore have the form

    S₂ = ∫√−g [ R + R F₀(□) R + C_abcd F₂(□) C^abcd ]

The basis {R², C²} is complete, because R_ab² is converted into these two through the Gauss–Bonnet term.

Static plane waves ω = a cos kx (conformal factor) and ν = b cos kx (optical index) give (sympy):

    R⁽¹⁾ = 2k²(3a + 2b) cos kx,    C⁽¹⁾·C⁽¹⁾ = (4/3) b² k⁴ cos² kx

The Weyl tensor thus gives only ν; the conformal sector is Weyl-free. Averaged over a period, in units of k²/2:

    A (ω²) = 6 + 9φ₀,   B (ων) = 8 + 12φ₀,   C (ν²) = 2 + 4φ₀ + φ₂,
    φ₀ = 4k²F₀,   φ₂ = (4/3)k²F₂.

Einstein–Hilbert (φ = 0) gives 6, 8, 2, checked by a direct computation of √−g R to second order.

## 3. The field of a point mass

Matter couples only to ω = Φ. Solving the quadratic form with a source gives (sympy):

    Φ/Φ_GR = (4/3) S₂(k) − (1/3) S₀(k)
    Ψ/Ψ_GR = (2/3) S₂(k) + (1/3) S₀(k)
    γ(k)   = (2S₂ + S₀)/(4S₂ − S₀)
    S₀ = 1/(1 + 3φ₀/2),   S₂ = 1/(1 − 3φ₂/2)

This is the standard decomposition into spin 2 and spin 0. **F6.8a measured the whole spin-0 channel:**
1 + 3φ₀/2 = A/6 = K̂, so

    S₀(k) = 1/K̂(k, T) = 1/[K_loc(k) − 0.10 k³T],   1 − K_loc ≈ (kℓ)²(4.0 ln(1/kℓ) + 0.9)

The spin-0 propagator has a static pole where K̂ = 0:

- without the age (for T → 0) at k ≈ 0.49/ℓ_k;
- with the age term the pole moves to longer wavelengths, k ≈ (10/T)^{1/3} ℓ_k^{-2/3}, i.e. to the length scale
  (0.1 ℓ_k² T)^{1/3}.

For a Planckian ℓ_k and T = the age of the universe this is ~10⁻¹⁵ m. What the pole means physically (a standing wave,
an instability of the static sector, or only the limit of validity) is open.

## 4. What is missing: the spin-2 form factor φ₂(k)

φ₂ appears only in the ν² sector (Weyl ∝ ν). We therefore need the second order in the optical perturbation, i.e.
F6.8c:

- the optical length to second order: index ½∫h² and the bending of the ray −½∫∫G ∇⊥h·∇⊥h (Fermat, Green function with
  fixed ends);
- the volume of an interval to second order: a surface term on the boundary of the spheroid δ(ℓ₀) (δ_yz + δ_zx)².

Without φ₂ neither Φ(r) nor γ(r) can be closed.

**F6.8b (the mixed sector ων)** needs only first-order cones. Covariance predicts **B/8 = K̂(k) for every k**: the same
function as in F6.8a, without a free parameter. It is a sharp test of the assumption that the discrete action has a
covariant quadratic form, including the question whether the age term also behaves covariantly. F6.8b at the same time
builds the geometry of the cones (delays along segments inside diamonds), which F6.8c extends to second order.

## 5. Order

1. **F6.8b:** the covariance test B/8 = K̂ and the construction of first-order optics.
   **Done (F6.8b2, PASSED 29 Oct):** K0_ab = 0.9979, c_mix = 7.983. Beyond the leading order Δ ≈ 0.46 k² and
   c_T,ab ≈ 0.195 against 0.110, which a covariant action cannot give (suspected age cut, F6.8d). The field of a point
   mass is therefore to be composed directly from A, B, C: γ = B/(2C) − 1, Φ/Φ_GR = −2/(A − B²/4C); the form factors
   φ₀, φ₂ hold only at leading order.
2. **F6.8c:** φ₂(k), the second order of the optics. It is the critical path to a prediction.
3. **Prediction:** Φ(r), Ψ(r) and γ(r) of a point mass from S₀ and S₂ (Fourier inversion), with the shape in units of
   r/ℓ_k and with the age T.
4. **Comparison:**
   - torsion balances (deviations from 1/r below 1 mm) bound ℓ_k;
   - γ(r) and Cassini;
   - a novelty search for the age term and the form factors of the BD action.

## 6. The sympy script (check of section 2)

Metric g = e^{2ω} diag(−1, e^{2ν}, e^{2ν}, e^{2ν}) with ω = εa cos kx and ν = εb cos kx. Procedure:

- √−g R to ε², averaged over a period: (k²/2)(6a² + 8ab + 2b²);
- R to ε: 2k²(3a + 2b) cos kx;
- Weyl to ε, contraction with η: (4/3) b² k⁴ cos² kx.

The decomposition into S₀, S₂ follows from A_eff = A − B²/(4C) = (2 + 3φ₀)(3φ₂ − 2)/(2 + 4φ₀ + φ₂) and from
ν = −Bω/(2C).
