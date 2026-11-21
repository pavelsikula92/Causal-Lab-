# Card F6.8g (proposed ID): the local optical stiffness K_c(k) at kℓ_k = 0.02–0.12 (draft)

**Revision 0, draft, 21 Nov 2026 (Claude; decision 41, item 4: option A). Not frozen.** The ID is the author's choice.

## Why

OV3-lin-a (docs/derivations/order_variation_OV3.md, sections 1, 5 and 8) solves the static linear response of the
order per k. It needs the local part of the optical (ν²) stiffness K_c(k) = C(k)/2 at the same k as:
- the conformal K_loc(k) (F6.8a);
- the mixed K_loc,ab(k) (F6.8b2);
- the source form factors F_ω, F_ν (F6.10b: k = 0.03, 0.06, 0.12).

F6.8c2 measured 𝒞(0) = 2.002 at the preregistered n = kT = 8. That value contains the age contribution, so it is not a
measurement of the local K_c(k) (decision 41). F6.8e3 showed that the optical age term grows as (kT)⁵ at large n. At
small n its form is not fixed by the data: F6.8c2 found s n + q n² with coefficients that depend on the range.

## What is computed

- The code is unchanged: phases/F6_order_dynamics/F6_08c_bd_action_optical.py, optical_row(k, T).
- Each row gives the phase-averaged second-order densities of both the optical sector K̂_bb(k, T) and the conformal
  sector K̂_c(k, T) on the same nodes.
- **Rows.** k = 0.02, 0.03, 0.045, 0.06, 0.08, 0.10, 0.12, at the ages of the ladder
  T ∈ {20, 28, 40, 57, 80, 113, 160} with n = kT ≤ 12. Rows at T = 10 and 14 are reported, not used in the fit.
  - T ≤ 160 stays below the ages where F6.8c lost its scaling.
  - The young end, n ≈ 0.4–1, is reached at small k.
- **Why T ≥ 20** (timing run, 21 Nov, k = 0.06). The optical row at T = 10 is dominated by the boundary term of the
  age cut: K̂_bb = 1.600, boundary part −0.0090 against J_bb 0.0035. At T = 20 the row gives K̂_bb = 0.988. The
  conformal row is 0.947 at T = 10 and 0.953 at T = 20, against F6.8a's local 0.956.
  - Below T ≈ 20 a young-age term that depends on T alone is not small. The model k² F(kT) does not describe it.
  - I1 checks on the conformal sector that the fitted range is free of it.

## The separation (the reading, stated before the run)

**The model** for each sector s (bb and c), fitted jointly over all rows:

    K̂_s(k, T) = K_loc,s(k) − k² F_s(kT),     F_s(n) = Σ_{j=1}^{J} a_{s,j} n^j,   F_s(0) = 0.

- K_loc,s(k) is free per k.
- The age component F_s is one function of n = kT for all k. This is the scaling that the conformal and mixed data
  obey to < 1e-4 (F6.8c2).
- F_s(0) = 0 defines the age component: a field of zero age has none.
- The main reading uses J = 2 and all rows with n ≤ 8.

**Reported separately for every k** (decision 41, item 4):
1. the total measured optical stiffness K̂_bb(k, T), every row;
2. the age component k² F_bb(kT), every row;
3. the local part K_loc,bb(k);
4. the sensitivity of the separation: the spread of K_loc,bb(k) over the variants
   - n ≤ 4 / 8 / 12;
   - T ≥ 20 / 28;
   - J = 2 / 3;
   - and K_loc taken directly as the youngest row (n ≤ 1 where it exists).

## Instruments (INSTRUMENT_FAILED if any fails)

- **I1, the method on the known sector.** The same separation applied to the conformal K̂_c of the same rows must
  reproduce F6.8a's K_loc(k) (its fit K0 − k²(c1 ln(1/k) + c0), K0 = 1.00005, c1 = 3.98, c0 = 0.91) within 1e-3 at
  every k ≤ 0.1. This validates the separation on a sector whose local part is known.
- **I2, the geometry of F6.8c** as there:
  - ⟨T2⟩ against the numerical Fermat ray ≤ 2e-3;
  - d²V/db² against numerical roots ≤ 2e-3;
  - the uniform null ≤ 1e-9.
- **I3, the level.** Level 2 against level 1 on three check rows ((0.03, 28), (0.06, 57), (0.12, 40)) ≤ 2e-3.
- **I4, the k = 0 null** at T = 20, 40, 160 in units of K̂_bb at the smallest k ≤ 2e-3.

## Criterion

This card is a measurement: its numbers are inputs of OV3-lin-a. The known answer judged is the long-wavelength limit
of the local part. Fit K_loc,bb(k) = K0_bb − k²(c1 ln(1/k) + c0) over k ≤ 0.1:
- **PASSED** if |K0_bb − 1| ≤ 0.01 (C → 2, as F6.8c2) and the sensitivity spread of item 4 is ≤ 0.005 at every
  k ≤ 0.1;
- **FAILED** if |K0_bb − 1| > 0.03;
- **UNDECIDED** otherwise, including a larger spread: the separation is then not stable enough for OV3-lin-a.

## The covariance check (decision 41: φ₂ measured independently)

- **The relation.** In the covariant form S₂ = ∫√−g [R + R F₀(□) R + C F₂(□) C] (path_to_predictions.md §2):
  - A = 6 + 9φ₀;
  - B = 8 + 12φ₀;
  - C = 2 + 4φ₀ + φ₂.
- **Why this card cannot test it.** Here φ₂ would have to be computed from the C it is meant to check.
- **The independent route** is a sector where only the Weyl part enters: a static transverse-traceless shear
  h_xy = b cos(kz). Its R⁽¹⁾ = 0, its form factor is 1 + c φ₂ alone, and its EH coefficient is computed from √−g R to
  second order.
- **What it needs.** New second-order geometry: anisotropic cones and the volumes of sheared intervals. It is a
  separate card (proposed F6.8h) after a development step. Its feasibility is not yet known.
- **Until then.** This card reports φ₂^(C) = C_loc − 2 − 4φ₀ only as a number to be checked later, never as a check.

## Cost

Timing run (21 Nov, Claude's workspace, two loaded cores): optical_row costs 136 s, 154 s and 171 s at (0.06, 10),
(0.06, 20) and (0.12, 40), about 70–90 s on a free core. The ladder above has about 40 rows, plus the check rows at
level 2 (several times a row). That is about 1–2 CPU-h, or about 10–15 min wall with 16 processes.

## Decisions for the author

1. The ID and the rows (k and the age ladder).
2. The model of the age component (J = 2, n ≤ 8) and the criterion of the sensitivity spread (0.005).
3. Whether F6.8h (the TT shear, φ₂) is developed now.
