# Q.1c — the vacuum of one 1+1 sprinkled order: Johnston's Feynman propagator (card, draft 0, 28 Nov 2026)

**Status: draft 0. Nothing here is frozen.** The author freezes the card. It is the second test of the quantum section,
chosen by the author before the interference test (route "C, then A", 28 Nov 11:29). Development: diag/q1c_dev.py,
diag/q1c_sj.py, diag/q1c_dev2.py.

## Question in plain words

- Q.1 showed that one sprinkled order carries the **propagation** of a quantum field: the retarded propagator.
- Q.1c asks whether one order also carries its **vacuum**, the ground state. In quantum field theory, even an empty
  region has a correlation of the field between two points; that correlation is the imaginary part of the Feynman
  propagator.
- The vacuum here is built **from the order alone**, with no external choice (the Sorkin–Johnston state, SJ). The
  question is whether it agrees with the vacuum of continuous spacetime in the **same region**, and how the
  difference falls as the density grows.

## Construction (only from the matrices of the order)

Johnston, arXiv:0909.0944:
1. K_R = aC(1 − abC)⁻¹, a = 1/2, b = −m²/ρ: the retarded propagator from the causal matrix C, as in Q.1.
2. K_A = K_Rᵀ.
3. Δ = K_R − K_A; iΔ is Hermitian.
4. Q = Σ_{λ>0} λ u u†, the positive part of iΔ, over the full spectrum. **No eigenvalue is dropped.**
5. K_F = K_R + iQ.

Then Re K_F = (K_R + K_A)/2, whose mean is exact (Q.1). The content of the card is **Im K_F = Re Q**: a non-linear
(spectral) function of the order, whose mean is not exact. Because Δ is real antisymmetric, Q = (|Δ| + iΔ)/2 with
|Δ| = (ΔᵀΔ)^{1/2}, so Re Q = |Δ|/2. This identity is checked as an instrument (I1).

## Reference: the continuum SJ state of the same diamond

- **Computed independently of any causal set** (diag/q1c_sj.py, record diag/out/q1c_sj_20261128T105641Z):
  - the exact continuum Pauli–Jordan kernel (½J₀(mτ) inside the light cone, half the jump on it);
  - product trapezoid quadrature on nested grids n = 24, 48, 96, and its positive part;
  - extrapolated in n with the order of convergence estimated from the three grids.
- **Uncertainty:** the size of the last extrapolation step, ≤ 0.006 in units of the amplitude ¼, mostly ≤ 0.002.
- **Not the reference: the Minkowski vacuum** (−¼Y₀ timelike, K₀/2π spacelike). The vacuum of a finite region differs
  from it by a near-constant offset that grows as mL falls (continuum SJ against Minkowski, in units of ¼):
  +0.035 … +0.039 at m = 5, +0.003 … +0.005 at m = 10, ≈ 0 ± 0.005 at m = 20. That offset is reported, not judged
  (the author, 11:43).

## Model

- **Region:** the light-cone square u, v ∈ [0, 1], as in Q.1. The mass is fixed per point: m ∈ {5, 10, 20} (mL).
- **The only scaled quantity:** the density, N_exp ∈ {1000, 2000, 4000, 8000}; m²/ρ = 0.10 at most.
- **Probe pairs:** 12 fixed pairs in the centre (u, v = ½ ± k/24, k ≤ 8), added as fixed elements to the sprinkling:
  - 4 timelike pairs with Δu = Δv;
  - 4 timelike pairs with Δu = 2Δv;
  - 4 spacelike pairs;
  - mτ or mσ from 0.29 to 6.7. Null-separated pairs are excluded: the function diverges on the light cone.
- **Probes in general position.** No two inserted points share u or v. The pairs come in 5 groups and each group gets
  its own sprinklings.
  - Reason (development, records diag/out/q1c_dev2_20261128T115328Z and 124921Z): a point exactly on the light cones
    of a pair, which has probability zero in a sprinkling, moved Re Q of that pair by −0.32 (in units of ¼).
- **Realisations:** R = 48 per (m, N_exp, group). Seeds derive_seed(base, "Q.1c", m, N_exp, group, r).

## Metrics

For pair j at (m, N_exp), realisation r:
- d_rj = (Re Q_r(p_j, q_j) − W_SJ(p_j, q_j)) / (1/4).

Two metrics:
- **Systematic deviation (bias):** b_j = mean_r d_rj, with SE_j² = var_r(d_rj)/R + (reference uncertainty_j)². The
  test statistic at each N_exp is χ² = Σ_j (b_j/SE_j)² over the 12 pairs (12 degrees of freedom).
- **Fluctuation of one order:** S = median over pairs of sd_r(d_rj), the typical error of one order at one pair. Its
  slope s against ln N_exp over the four densities has a bootstrap 95 % interval [lo, hi] (realisations resampled
  within each pair and density, B = 2000).

## Criterion (proposed)

Per m:
- **bias:**
  - consistent if χ² at N_exp = 8000 has p ≥ 0.01;
  - deviant if p < 0.001;
  - unclear otherwise.
- **fluctuation:**
  - falls if hi < 0;
  - grows if lo > 0;
  - unclear otherwise.

The verdict:
- **INSTRUMENT_FAILED** if any instrument fails. No statement is made then.
- **FAILED** if the bias is deviant at some m (a systematic difference from the continuum vacuum remains at the
  highest density), or the fluctuation grows at some m.
- **PASSED** if the bias is consistent and the fluctuation falls at all three m.
- **UNDECIDED** otherwise, with the sub-label "inconclusive". The bands are not widened after the run.

No rate is predicted for the fluctuation, unlike Q.1. The spectral construction need not self-average at the rate
ρ^{−1/2}, so the card judges only that it falls and reports the rate.

## Instruments (INSTRUMENT_FAILED if any fails)

- **I1 (the identity):** on 10 orders of N = 50–300, Re Q from |Δ|/2 against the explicit complex eigendecomposition
  of iΔ (positive part), max |difference| ≤ 1e-10; and Im Q = Δ/2, ≤ 1e-12.
- **I2 (the retarded part):** at every probe pair (timelike), the mean of K_R(p, q) equals ½J₀(mτ) within 3.5 SE
  (Q.1's exact mean, here at the inserted pairs).
- **I3 (general position):** no two inserted points share u or v (a check of the construction).
- **I4 (the null, non-negative weights):** the same construction with b = +m²/ρ at m = 10, N_exp = 2000. Its "vacuum"
  must differ from W_SJ: χ² p < 1e-6. If it agrees, the agreement of the quantum case would not be specific to the
  oscillating amplitude.

## Reported, not judged

- b_j, sd_j per pair and density.
- The slope of the rms bias against ln N_exp.
- The fluctuation slope s per m.
- The offset of the continuum SJ state against the Minkowski vacuum.
- The rows at N_exp ≤ 4000 for comparison with the development.

## Reading declared now

- **PASSED:**
  - One sprinkled 1+1 order carries the vacuum of the same region: the imaginary part of Johnston's Feynman propagator,
    built only from the order. No systematic difference from the continuum SJ state remains at N_exp = 8000 (within
    the statistics and the reference), and the error of one order falls with the density.
  - Together with Q.1, the order carries both propagation and the ground state. Those are the two ingredients of the
    interference test (A), so the quantum instrument is ready for it.
  - It says nothing about 3+1 (links; open), interference, or the Born rule.
- **FAILED:** a systematic difference remains, or one order does not settle. Test A would then need another vacuum
  (or an ensemble), and that is decided before it is designed.

## Development diagnostics (not a result of Q.1c)

Record diag/out/q1c_dev2_20261128T124921Z: probes in general position, 24 realisations, N_exp 1000 / 2000 / 4000.

| m | rms bias over the 12 pairs | χ²/12 at 4000 | median sd of one order | sd slope (rough) |
|---|---|---|---|---|
| 5 | 0.018 / 0.010 / 0.009 | 1.47 | 0.052 / 0.039 / 0.033 | ≈ −0.33 |
| 10 | 0.025 / 0.014 / 0.008 | 0.94 | 0.063 / 0.046 / 0.034 | ≈ −0.45 |
| 20 | 0.061 / 0.031 / 0.014 | 1.40 | 0.118 / 0.083 / 0.061 | ≈ −0.48 |

- The bias falls with the density at all m.
- The fluctuation of one order at one pair is larger than in Q.1 (3–6 % of ¼ at N_exp 4000) and falls more slowly at
  m = 5.
- The first version, with all probes inserted together, is in record 115328Z. It shows the null-coincidence artifact.

**Cost:** dense matrices of N ≈ N_exp + 6. At N_exp = 8000, about 2 GB per worker and a few minutes per realisation.
The full run is several hours on the author's PC; the workers are to be chosen by its memory.
