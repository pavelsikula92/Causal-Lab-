# Card OV1-V (2D): no volume term in the smeared BD action of flat causal diamonds (the vacuum control of OV1)

**Revision 2, 20 Nov 2026, for freezing.** ID F6.9a (the author, 20 Nov 10:10). Test:
phases/F6_order_dynamics/F6_09a_bd_vacuum_2d_diamonds.py (self-contained; kernels copied from diag/ov1_kernels.py).
Revision 2 differs from revision 1 only in the instrument items I3 and I4 and the seeds (section 8). Specification:
docs/derivations/order_variation_OV1.md, sections 4, 5, 6, 8 and 13. Development: diag/ov1_2d.py, diag/ov1_2d_summ.py.
The tolerances and the numbers of sprinklings below follow from the development and validation rows of OV1 §13, by
the power rule of OV1 §11.

## 1. Hypothesis

The vacuum control of OV1 (§8) is blocking. Lorentz invariance allows the linear response δ⟨S⟩ = c(T) ∫ tr h dV, and
the flat order is stationary only if c(T) = 0.

In 2D the readout is the volume term of the ensemble mean. The smeared BD action of a flat causal diamond has no term
proportional to the number of elements: ⟨S⟩ does not grow with the volume over the nested diamonds T/w = 8, 16, 32
(volumes 1 : 4 : 16) at ε = 0.06, 0.12 and 0.24. The continuum expects a constant, the joint term
(arXiv:1502.05388). The constant is reported, not judged.

## 2. Model (as OV1 §13, unchanged)

- **Order.** A Poisson sprinkling of the light-cone square [0, 1]² of mean N_exp = 1 024 n_loc(ε), with x ≺ y iff
  u_x < u_y and v_x < v_y.
- **Nested diamonds.** The top-corner squares of side 1, 1/2 and 1/4, with a common future tip, give T/w = 32, 16 and
  8. Here n_loc(ε) is the smallest n with 99 % of Σ|f(m, ε)| and w² N_exp = n_loc (OV1 §4).
- **Action.** S = 4ε (N − 2ε Σ_n N_n f(n, ε)) (F6.4a convention). The N_n are counted exactly for n ≤ n_cut, where
  |f| ≥ 1e−14.
- **ε.** 0.06, 0.12 and 0.24 (the author, 19 Nov), with ε defined as in OV1 §6.

## 3. Readout

- Per sprinkling, c is the least-squares slope of S against N over the three nested levels.
- The readout per ε is c/(4ε), the volume term per element in units of the bare term 4ε: the mean ± SE over
  sprinklings.
- **Reported, not judged:**
  - ⟨S⟩ per level and the paired differences S(16) − S(8) and S(32) − S(16);
  - the constant (⟨S⟩ at T/w = 8);
  - the phase diagnostics χ, d_MM and h per level, with the fraction outside the calibration bands.

## 4. Criteria (proposed; frozen only by the author's decision)

- **Tolerance.** |c/(4ε)| + 2 SE ≤ 1·10⁻³ at each ε: the volume term cancels per element to 0.1 %.
  - The validation rows give |c/(4ε)| ≤ 2.2·10⁻⁴. The tolerance was not set from them; it follows from the power
    rule.
- **Power** (2 SE ≤ tol/3, OV1 §11). SE per sprinkling, from the validation rows: 1.22·10⁻² at ε 0.24, 5.7·10⁻³ at
  0.12, 2.7·10⁻³ at 0.06.

  | ε | sprinklings | expected SE | cost per sprinkling | core-hours |
  |---|---|---|---|---|
  | 0.24 | 5 400 | 1.66·10⁻⁴ | 8 s | 12 |
  | 0.12 | 1 200 | 1.65·10⁻⁴ | 30 s | 10 |
  | 0.06 | 280 | 1.65·10⁻⁴ | 164 s | 13 |

  In total about 35 core-hours, or 4.5 h on 8 workers. With the balance of I4 on all sprinklings it is about 40
  core-hours, 5 h.
- **New sprinklings.** Seeds derive_seed(base, "F6.9a", ε, k), the laboratory's SHA-256 streams. They are independent
  of the development, validation and calibration seeds.
- **Verdict per ε:**
  - PASSED if |c/(4ε)| + 2 SE ≤ 1·10⁻³;
  - FAILED if |c/(4ε)| − 2 SE > 1·10⁻³;
  - UNDECIDED otherwise.

  The card is PASSED only if all three ε pass, and FAILED if any fails.
- **Instrument (INSTRUMENT_FAILED if any fails):**
  - **I1:** `python -m diag.ov1_2d --check` passes: S and every ΔS against brute force to 1e−8, including two orders
    where n_cut < N.
  - **I2:** the mean number of elements per level is within 3 SE of N_exp s².
  - **I3:** the phase diagnostics of the run lie within the calibration bands (written into the test from
    diag/out/ov1_2d/summary.txt, md5 91afb053), with at most 5 % of the sprinklings outside at any level and ε
    (nominal 1 %).
  - **I4:** the time-reversal null. ⟨D S⟩ over the interior is zero within |x| + 2 SE ≤ 1·10⁻³ at T/w = 8, over all
    sprinklings. It is guaranteed by the symmetry of the flat diamond, so it checks the code, not physics.

## 5. Decisions for the author

1. The ID.
2. The tolerance (1·10⁻³ per element) and the numbers of sprinklings above.
3. Whether the constant is only reported (as proposed) or gets a criterion of its own. A criterion would need its
   continuum value in the F6.4a normalisation, derived first, and the value is not yet derived here.

## 6. Reading declared now

- **PASSED:** the BD action of the flat 2D order carries no volume term at the three ε. OV1's blocking vacuum control
  holds in 2D, and 2+1 follows (OV1 §12).
- **FAILED at some ε:** the construction (action with these boundaries) fails the vacuum control there. The origin of
  the term is explained before anything else, and it is not subtracted as an artefact (OV1 §8).
- **What a PASS does not show:**
  - that the flat order is stable against changes of the order; the phase sampler P1 shows it is not above β_e(N);
  - anything in 2+1 or 3+1;
  - the 3+1 consistency with F6.8 (OV1 §9).

## 7. Development and validation rows (decide nothing)

OV1 §13: development at ε = 0.12 (16 sprinklings); validation on the author's machine, 19–20 Nov.

| ε | sprinklings | c/(4ε) | ⟨S⟩ at T/w = 8 |
|---|---|---|---|
| 0.24 | 512 | (+2.2 ± 5.4)·10⁻⁴ | +3.1 ± 3.1 |
| 0.12 | 512 | (−1.7 ± 2.5)·10⁻⁴ | +6.0 ± 1.5 |
| 0.06 | 256 | (+0.9 ± 1.7)·10⁻⁴ | +2.8 ± 1.0 |

Calibration bands (512 / 512 / 256 seeds ≥ 100 000): 0.4–2.0 % of the validation sprinklings fall outside, against a
nominal 1 %.

## 8. Revision 2 (20 Nov 2026): instrument settings before freezing

- **I3, 3 % → 5 %.**
  - The bands at ε 0.06 come from 256 calibration seeds, and their realised rate outside was 0.4–2.0 % in the
    validation rows.
  - At 280 sprinklings and a true rate of 1.6 %, a 3 % threshold fails by chance with a probability of about 4 % per
    item, and there are six independent items per ε. That is an underpowered control of the kind that failed F5.2c
    seed 1.
  - At 5 % the false alarm is below 1·10⁻³ per item. An escape from the phase would move every sprinkling out (P1).
- **I4, 256 sprinklings → all.** With 256 sprinklings at ε 0.24, 2 SE would be 6·10⁻⁴, against the rule 2 SE ≤ tol/3 =
  3.3·10⁻⁴. With all sprinklings 2 SE ≤ 1.4·10⁻⁴ at every ε.
- **Seeds.** derive_seed instead of 200 000 + k, the laboratory convention.
- **Unchanged:** the hypothesis, the readout, the criterion (tolerance 1·10⁻³, 2 SE), the numbers of sprinklings and
  the reported constant (the author, 20 Nov 10:10).
- **Quick run** (Claude's container, `python run.py F6.9a --quick --workers 2`, decides nothing):
  - it ran mechanically in 150 s;
  - I1: S and ΔS against brute force to 1.1e−10, including N = 2 523 and 1 823 with n_cut < N;
  - I2: OK;
  - I3 and I4 fail at 2–6 sprinklings, as expected for a smoke run.
