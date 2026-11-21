# Card F6.12 (draft, NOT frozen): β of the order, the cubic static action along the exponential metric at n = kT = 4

**Status: draft 0, 21 Nov 2026, for the author's decisions.** The ID F6.12 is a proposal. There is no code yet. The card
cannot be frozen until the code exists and every judged point has been exercised at full resolution.

Derivation: docs/derivations/order_variation_OV3.md, section 10; diag/ov3_beta_sym.py (sympy).

## Why

The paper's field equation of the order is computed to first order only (Part V, "The field equation, by order in the
field": second order (β) not computed). β is the first number of the second order.

**One number is needed (OV3 10.1, sympy).** For the most general local static action with two derivatives, the
point-mass solution to second order gives β − 1 = −κ₃/8, where κ₃ is the cubic action along the first-order solution:
- that solution is (ω, ν) = (−U, 2U), the exponential metric g = −e^{−2U} dt² + e^{2U} dx²;
- along it the EH action is −2(∇U)² at every order, so EH has κ₃ = 0 and β = 1.

The test is therefore a **null test of the order's cubic action** along one direction. The six cubic coefficients need
not be separated.

## What is computed

E(A, C): the expected smeared BD action per unit four-volume for the static planar metric

    g = −e^{−2u} dt² + e^{2u} dx²,   u(z) = A cos kz + C cos 2kz,

in the young-field reading of F6.8g: a past cut at age T = n/k with n = 4. Units ερ = 1, ℓ_k = 1.

- **Quadratic part:** E₂ = −K̂_dir (A² + 4C²) k².
- **Cubic part:** E₃ = κ₃(k) (3/4) A²C k², i.e. κ₃ = (2/3) ∂³E/∂A²∂C / k² at A = C = 0.
- **The response:** β(k) − 1 = −κ₃(k) / (8 K̂_dir(k)).
- **Cross-check of the quadratic part:** K̂_dir = (16 K̂_ab − 6 K̂_c − 8 K̂_bb)/2 from the F6.8g rows (1.014 / 1.044 at
  k 0.06 / 0.12). Recomputing it is a check of the code.

**Method: route (a) of OV3 10.3, perturbative and exact order by order.**
- the Mecke formula with the exact expected kernel g_k, expanded to third order;
- the interval volumes to third order in the cone perturbation ν = 2u, a new Fermat term T₃ for a stratified medium;
- the density e^{4ω+3ν} = e^{2u}, expanded to third order;
- the age cut at third order.

No finite differences in amplitude.

**Rows (judged):** k = 0.045, 0.06, 0.08, 0.10, 0.12 at n = 4 (T = 89 … 33), as in F6.8g.
**Reported:** the same at n = 3 and 6.

## Reading (to be stated before the run)

At fixed n, fit κ₃(k) = κ₃0 − k² (d₁ ln(1/k) + d₀) over the five k (the fit form of F6.8g). Then β₀ = 1 − κ₃0/8, since
K̂_dir → 1.

## Instruments (INSTRUMENT_FAILED if any fails)

- **I1, the conformal cubic (stage β1 of OV3 10.3; the known answer).** The same code with ν = 0 (ω = u: the order is
  unchanged, only the density changes) must give the EH cubic 12 ω(∇ω)², i.e. 9 A²C k², at k → 0. The conformal cubic
  intercept must be 1 within a tolerance to be fixed after development (proposal 2e-3).
- **I2, the third-order volumes (stage β2).** T₃ against brute-force interval volumes (the K.18b code generalised to
  the two-mode profile), at small amplitude, for a set of intervals, including near-null ones. Tolerance to be fixed
  after development.
- **I3, the exact null.** A single mode (C = 0) has no cubic term, because a half-period shift maps A → −A. The A³
  coefficient computed by the full code must vanish to rounding (proposal: below 1e-9 of the A²C coefficient).
- **I4, the level.** The rows (0.045, n 4) and (0.12, n 4) at level 2 against level 1, |Δκ₃| ≤ 2e-3 (as F6.8g I3).
- **I5, the quadratic part.** K̂_dir from this code against the F6.8g rows, within 1e-4.

## Criterion (proposal)

- **PASSED** if |β₀ − 1| ≤ 0.01 (|κ₃0| ≤ 0.08);
- **FAILED** if |β₀ − 1| > 0.03;
- **UNDECIDED** otherwise.

## Reported, not judged

- κ₃(k) and β(k) for every row. This is the collinear form factor of the cubic vertex at that k, not the β of a point
  mass at that scale (OV3 10.2).
- The fits at n = 3, 4, 6.
- K̂_dir(k).

## What this test cannot say (OV3 10.4)

- **The young-field local limit is expected to give β = 1.** The BD continuum theorem plus covariance predict it. The
  test checks this at cubic order, as F6.8g did at quadratic order.
- **It does not give the β of the Sun's field.** That field is old: ε = b(kT)² ≫ 1, where F6.8f shows no saturation.
  That is the age problem (F6.8i).

## Cost

- **Unknown until stage β2 exists.** β1 is an extension of F6.8a, about one day of development. β2 (T₃ and its
  brute-force check) is the main work, several days.
- One row costs a few times an F6.8g row; F6.8g took 4091 s on 8 workers.

## Decisions for the author

1. The reduction: β from the single cubic number κ₃ along the exponential metric (OV3 10.1). Accept?
2. Route (a), perturbative (recommended), or route (b), finite amplitude with brute-force volumes (needs a precision
   of ~1e-10 relative, not shown by the existing code).
3. The order of work: β1 (conformal cubic) first, as a development step and later as I1; then β2; then the card is
   completed and frozen.
4. The ID F6.12, and the criterion |β₀ − 1| ≤ 0.01 / > 0.03.
