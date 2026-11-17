# Card F6.8e3: second revision of F6.8e (the instrument checks stated in the exponent)

**15 Nov 2026.** Prepared and frozen on the author's instruction of 22:38 ("test by se měl opravit a spustit poté
znovu"; decision 32).

**Test**: phases/F6_order_dynamics/F6_08e3_bd_action_optical_asymptotic_rev2.py.

## 1. Why a second revision

F6.8e and F6.8e2 both ended INSTRUMENT_FAILED on checks whose tolerances or ranges were set wrongly. The judged
quantity α was never the failing item, and both runs recorded α = 5.000 at the largest ages. The four failed checks
were:

| Run | Check | What went wrong |
|---|---|---|
| F6.8e | I5 | the density ratio of two k was tested against (k₁/k₂)⁵, but the normalisation carries a 0.2 % k-dependence |
| F6.8e2 | I3 | the conformal control was placed at r = 61 440, beyond the range of its double-precision path |
| F6.8e2 | I5 | a ratio drift of 1.7 × 10⁻³ (an exponent difference of 8 × 10⁻⁴) against 10⁻³ |
| F6.8e2 | I2 | a level difference of 1.2 × 10⁻³ at the smallest age, which does not matter for an exponent |

**The common fault.** The checks were stated in raw densities, with tolerances unrelated to the judged quantity.

## 2. Changes (nothing else changes)

- **I2, levels.** The exponent from the p = 0 densities at the judged pair n 3 072 → 6 144, at two quadrature levels.
  Required: |Δα₀| ≤ 0.01.
- **I3, conformal control.** Moved to n = 384, r = 2 954, where the conformal path is exact. Tolerance 3 % as before.
- **I5, the kT-shape independent of k.** The exponent from the p = 0 densities at the judged pair at k = 0.13 and 0.21.
  Required: |Δα₀| ≤ 0.01. The ratio against (k₁/k₂)⁵ is reported only.
- **Wavenumbers.** k = 0.13 (profile) and 0.21. The largest r is 47 262, where the optical levels agree.
  - k = 0.07 was checked before this card without printing values: its levels differ by 7.5 × 10⁻³ at r = 87 771, so
    it was not used.
- **Ages.** n = 192 … 6 144. The judged pair is n 3 072 → 6 144 at k = 0.13.
- **Status:** PASSED |α − 5| ≤ 0.1; FAILED > 0.5; UNDECIDED otherwise.

The check tolerance of 0.01 is ten times tighter than the status tolerance.

## 3. What the verdict means

The exponent is not blind: F6.8e and F6.8e2 recorded 5.000. This run gives the verdict with a corrected instrument and
is not new evidence beyond them.

Reading and limits are as in F6.8e (docs/F6_8e_card_draft.md, sections 6 and 7):
- bulk of the past only;
- second order in b;
- whether the age dependence is physical or an artefact of the cut is open.
