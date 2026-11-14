# Card C0, revision b (F3.16b): the 3+1 Green function with statistical criteria

**14 Nov 2026, ready for freezing, NOT frozen.** It follows F3.16 (INSTRUMENT_FAILED by design errors of its criteria;
docs/findings.md) and the author's plan of 14 Nov, 08:52.

**Test F3.16b**: phases/F3_structure/F3_16b_green_3p1_rev.py, md5 a58572c4c3cacc7ef602def7316e1362.

**Helpers:**
- diag/c0_card_b.py 2148fa9ace3716a0e613cae488746b0f
- diag/c0_fast.py 1dc98997a7e414a7d7eeb0b7a8b6a5da
- otherwise those of F3.16 (decision 27), unchanged.

The md5 sums are of the LF files.

## 1–3. Hypothesis, model, readouts

As card revision 1 (docs/C0_card_draft.md). The kernel is MV8+ct; MV2+ct is not reported in this revision. The rule
runs through diag/c0_fast.evolve_fast, which is bit-identical to diag/c0_card.evolve_card (checked on 25 136 and
77 330 events) and about 2× faster. ρ 30 with 16 sprinklings and ρ 60 with 12; seeds derive_seed(base, "F3.16b", ρ, k),
base 0 decides.

## 4. Criteria

### Why the form changed

A bootstrap of the development rows (diag/out/c0_card_b_dev_20261114T120302.json, seeds 19000+) resampled sets of the
card's size from sprinklings that behave exactly like the development. It found:

- **Revision 1 form**, tolerance 1.25 (|dev| + 2 SE) on about 45 quantities whose true value is near zero: PASSED in
  1 % of 100 resamples, INSTRUMENT_FAILED in 98 %. A near-zero quantity passes only when |d| < 0.5 SE.
- **The form below**: PASSED 89 %, UNDECIDED 10 %, INSTRUMENT_FAILED 1 % of 100 resamples. The tolerances were
  derived from the same data, so this rate is optimistic.

### Instrument

These predictions are exact in expectation (docs/derivations/c0_mean_exactness.md):
- ε(r) of the 10 shells against (A1);
- the normalisation;
- (A3) of the linear rows k 0–3, in every slab.

At both densities every item must have |d| / √(SE² + 0.002²) < t_crit:
- t_crit is Student t with n − 1 degrees of freedom at family-wise α 1 % over all 110 items (Bonferroni):
  5.29 at ρ 30, 5.99 at ρ 60;
- the floor 0.002 covers the evaluation of the predictions, which differ by about 1e-4 when dt is halved, and the
  first slabs, where the scatter between sprinklings is about 1e-4 and the predictions are off by up to 6e-4;
- power: a systematic error of a linear mode of about 0.03–0.04 in a slab is detected.

### Physics, tolerance criteria (MV8+ct)

The statistic is M = max over slabs (and over amplitudes for the conditions of B) of |d_j| + 2 SE_j.
- PASSED side: M < tolerance.
- FAILED side: max_j (|d_j| − 2 SE_j) > tolerance.
- The growth exponent p: its mean + 2 SE against the tolerance.
- tolerance = M_dev + z SE_M, rounded up, where SE_M is the bootstrap SE of the statistic M itself over the
  development sprinklings and z = 2.79 (one-sided, family-wise 5 % over the 19 criteria). This is the rule of 10 Nov
  (the development value + 2 SE + a margin) with the margin z SE_M of the statistic in place of the 25 %.

| criterion | ρ 30 | ρ 60 | M_dev (ρ 30 / ρ 60) |
|---|---|---|---|
| residual against B, A 0.5 | 0.055 | 0.045 | 0.029 / 0.025 |
| residual against B, A 1 | 0.04 | 0.04 | 0.020 / 0.023 |
| residual against B, A 2 | 0.065 | 0.04 | 0.038 / 0.022 |
| \|direct check\| | 0.014 | 0.011 | 0.012 / 0.008 |
| \|skewness\|, t ≥ 1.5 | 0.41 | 0.31 | 0.32 / 0.29 |
| \|excess kurtosis\|, t ≥ 1.5 | 0.6 | 0.615 | 0.52 / 0.57 |
| growth exponent p | 1.2 | 1.15 | 1.17 / 1.11 |
| linear (A), k 0–3 | — | 0.08 | — / 0.068 |
| nonlinear (A) against the pendulum, A 0.5 / 1 / 2 | — | 0.075 / 0.1 / 0.25 | — / 0.064 / 0.088 / 0.216 |
| nonlinear (A) against sine-Gordon, A 1 cos x | — | 0.09 | — / 0.072 |

Derivation: .

### Physics, sign criteria

Difference − 2 SE > 0 PASSED, + 2 SE < 0 FAILED:
- **convergence**: the largest |amp-dev| at ρ 30 minus that at ρ 60, k 0–3;
- **the noise falls**: 1 − sd(60)/sd(30) at t 4.75, every mode;
- **the renormalization falls**: renormalization(30) − renormalization(60) at A 1 and 2, with the SE by jackknife over
  the sprinklings (B depends on the run only through the pooled sd).

Verdict order: INSTRUMENT_FAILED → FAILED → PASSED → UNDECIDED.

## 5. Reading declared now

As card revision 1. **PASSED**: a kernel reading only the counts in intervals of the order, with a local term fixed by
the continuum, carries the 3+1 retarded Green function as predicted, and the linear and nonlinear field within fixed
tolerances of the continuum. The nonlinear softening is explained by the fluctuations, and errors and noise fall with
the density. **FAILED**: a criterion fails beyond its SE. **Limits**:
- flat order only;
- T ≤ 6;
- μ 1;
- one kernel family.
