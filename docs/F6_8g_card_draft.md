# Card F6.8g: the young-field stiffness of the order at fixed n = kT = 4, kℓ_k = 0.045–0.12 (FROZEN)

**Revision 2, FROZEN by the author, 21 Nov 2026 15:48 (decision 46).** Test file
phases/F6_order_dynamics/F6_08g_bd_action_optical_local.py, md5 d594c74386a690e4e978fb4bb65debad. Before freezing, every
judged point was exercised at full resolution: the rows k ≥ 0.045 at n = 4 on both machines, and the I3 point
(0.045, n = 4) at level 2 (3.0e-5).
- Revision 2 drops k = 0.03. Two development runs of revision 1 (Claude's workspace and the author's machine, not
  verdicts) showed that the old optical rows are not numerically robust at small k:
  - between the two platforms, K̂_bb differs by 2.7e-4 at (0.03, T 133), 1.4e-4 at (0.045, T 133) and 5e-5 at
    (0.06, T 100), but by ≤ 6e-6 at T ≤ 89;
  - the row (0.03, T 133) lies about 3.5e-3 below the n-dependence of every other k.
  On k ≥ 0.045 both runs give K0_bb(n = 4) = 0.9994 with a fit residual of 4.7e-5.
- Revision 0 (decision 43) separated a local part from an age part. Its quick run and a probe showed that this
  separation is not defined at order k². The probe numbers are in decision 44.
- Revision 1 reads the stiffness at a fixed age in units of the wavelength.

## Why

OV3-lin-a (docs/derivations/order_variation_OV3.md, sections 1, 5 and 8) solves the static linear response of the
order per k in the young-field reading (decision 41). It needs, at the same k and the same age:
- the conformal stiffness K;
- the mixed stiffness K_ab;
- the optical stiffness K_bb = C/2.

**The young field is defined as the field of age T = n_y/k with n_y = 4.** The source form factors F_ω and F_ν of
F6.10b do not depend on the age (T 80 against 120: below 0.1 %).

**Two facts set the design (decision 44):**
- **The k = 0 null.** It contains a young-age term independent of k, J_bb ≈ 13/T⁴. A row is clean only for
  nT ≳ 100, so n → 0 cannot be reached.
- **No unique local part.** The conformal sector at small n differs from F6.8a's asymptotic-intercept local part by
  about 1·k², growing with n. A local part is not defined at order k².

## What is computed

The code is unchanged:
- phases/F6_order_dynamics/F6_08c_bd_action_optical.py, optical_row(k, T): K̂_bb and K̂_c on the same nodes;
- phases/F6_order_dynamics/F6_08b2_bd_action_mixed_rev.py, mixed_row(k, T): K̂_ab and K̂_c.

The units are ερ = 1 and ℓ_k = 1.

**Rows:**
- **Judged:** k = 0.045, 0.06, 0.08, 0.10, 0.12 at n = 4, so T = 4/k = 89, 67, 50, 40, 33. Optical and mixed
  rows.
- **Reported:** optical rows at n = 3 (T = 100 … 25) and n = 6 (T ≤ 160; the row (0.045, n = 6) has T = 133 and is reported only), for the n-dependence.

## Reading (stated before the run)

For each sector s at fixed n, the fit (the method of F6.8c2) is

    K̂_s(k, n/k) = K0_s − k² (c1 ln(1/k) + c0)

over the five k.

## Instruments (INSTRUMENT_FAILED if any fails)

- **I1, the known answer at fixed n.** The conformal intercept at n = 4 is K0_c = 1 within 1e-3 (F6.8c2: 0.99976 at
  n = 6, 8 and 10).
- **I2, the geometry of F6.8c:**
  - ⟨T2⟩ against the Fermat ray ≤ 2e-3;
  - d²V/db² against numerical roots ≤ 2e-3;
  - the uniform null ≤ 1e-9.
  - Also the geometry check of F6.8b2: dV/db and the uniform null, with the tolerances of F6.8b2.
- **I3, the level.** The optical rows (0.045, n = 4) and (0.12, n = 4) at level 2 against level 1, |ΔK̂_bb| ≤ 2e-3.
- **I4, the null at the judged ages.** The k = 0 optical row at T = 4/k for each judged k: a₄|J_bb|/k² ≤ 2e-3. The
  predicted largest value is about 1.2e-3 at k = 0.12.

## Criterion

From the fit of the optical sector at n = 4:
- **PASSED** if |K0_bb − 1| ≤ 0.01 (C → 2 at long wavelength, at this age as at n = 8 in F6.8c2);
- **FAILED** if |K0_bb − 1| > 0.03;
- **UNDECIDED** otherwise.

## Reported, not judged

- K̂_bb, K̂_c and K̂_ab of every row.
- The fits at n = 3, 4 and 6 for the optical and conformal sectors, and at n = 4 for the mixed sector (F6.8b2:
  K0_ab = 0.9979 at its ages).
- **The inputs of OV3-lin-a** at k = 0.06 and 0.12: K̂_c, K̂_ab and K̂_bb at n = 4.
- φ₂^(C) = 2K̂_bb − 2 − 4φ₀ with φ₀ = (2/3)(K̂_c − 1). It is a number to be checked by an independent φ₂ later
  (decision 43: F6.8h is not developed now), never a check here.

## Cost

- About 14 optical rows of 70–170 s each, two level-2 rows, five null rows, and five mixed rows. F6.8b2 rows take
  2–30 min, depending on the age.
- In total about 2.5 CPU-h. With 16 processes the wall time is set by the longest single rows: about 1 h.

## Decisions for the author

1. Freeze this revision. The rows, the instruments and the criterion follow decision 44.
