# Card F6.8e2 (draft): revision of F6.8e after INSTRUMENT_FAILED on I5

**15 Nov 2026. Draft for the author's decision; not frozen.** Prepared on the author's instruction of 21:14 ("ano").

**Test**: phases/F6_order_dynamics/F6_08e2_bd_action_optical_asymptotic_rev.py.

## 1. What happened in F6.8e

F6.8e (decision 30, revision 1; results/F6.8e/20261115T134513Z) ended INSTRUMENT_FAILED on I5 alone.
- **The failed check.** The ratio of the optical densities at k = 0.15 and 0.25 was 7.793 × 10⁻² against
  (0.15/0.25)⁵ = 7.776 × 10⁻², a difference of 2.2 × 10⁻³ (n = 192) and 2.3 × 10⁻³ (n = 1 536). The tolerance was
  10⁻³, set from the development pair 0.12 / 0.2 (8 × 10⁻⁴).
- **What the difference is.** It grows with k (about 0.6 k⁴) and does not depend on n. It is a k-dependent
  normalisation of k² F(kT), not a change of its dependence on kT.
- **The design fault.** I5 tested the normalisation, while the question is the shape in kT.

**Seen, not a verdict (R132):** α = 4.937 … 5.000 for n = 96 … 3 072 at k = 0.15. This revision is therefore not blind
to the exponent. Its decision rests on rows not computed before it.

## 2. Changes (nothing else changes)

- **I5 tests what the question needs:** that the kT-dependence does not depend on k.
  - R(n) = I(k₁, n/k₁) / I(k₂, n/k₂) at p = 0 must be the same at two ages: |R(384)/R(3 072) − 1| ≤ 10⁻³. The F6.8e
    rows give 4 × 10⁻⁵.
  - The deviation of R from (k₁/k₂)⁵ is reported, not judged.
- **Wavenumbers:** k = 0.1 (profile) and 0.18 (scaling partner), used neither in the development nor in F6.8e.
- **Ages:** n = 192 … 6 144. The judged pair is n = 3 072 → 6 144, beyond the range of F6.8e.

## 3. Criteria (as F6.8e otherwise)

**Instrument:**
- I0: long double;
- I1: stable V'' against F6.8c, ≤ 10⁻⁹;
- I2: two levels, ≤ 10⁻³;
- I3: conformal c_T within 3 %;
- I4: tail ≤ 10⁻³;
- I5: as above.

**Status** (α between n = 3 072 and 6 144 at k = 0.1):

| Status | Condition |
|---|---|
| PASSED | \|α − 5\| ≤ 0.1 |
| FAILED | \|α − 5\| > 0.5 |
| UNDECIDED | otherwise |

**Run:** in Claude's Linux workspace, about 80 min. Deterministic.

**Reading and limits:** as F6.8e, sections 6 and 7.
