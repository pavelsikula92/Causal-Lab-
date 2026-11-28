# Q.1 — the sum-over-chains propagator on a 1+1 sprinkling (card, draft 1, 28 Nov 2026)

**Status: FROZEN** by the author, 28 Nov 2026, 10:47 ("Schváleno"; decision 50). The test file is phases/Q_quantum/Q_01_chain_propagator_2d.py (sha256 d58fd228…, LF md5 a596dfaa). Its smoke run (--quick, 40 realisations) decides nothing. From now on, the metric, the range and the bands do not change; a change would be a new revision with a new ID. It is the first test of the quantum section
(decision 49; docs/quantum/README.md section 7).

Draft 1 answers the author's four points on draft 0 (28 Nov, 10:39):
1. what is held fixed when N changes;
2. one realisation against the ensemble;
3. complete decision bands;
4. calibration kept apart from the physics verdict, and the development kept out of the result.

## Question

A particle on a causal set propagates from x to y as a sum over chains x = z₀ ≺ z₁ ≺ … ≺ z_n = y. Each relation
carries an amplitude a and each intermediate element an amplitude b (S. Johnston, arXiv:0806.3083, section 3; 1+1:
the causal matrix C):

    K = aC + a²b C² + a³b² C³ + … = aC (1 − abC)⁻¹,   a = 1/2,  b = −m²/ρ.

For two fixed points at proper time τ and a Poisson sprinkling of density ρ, the expected number of chains with k
intermediate elements is (ρV)^k/(k!)² with V = τ²/2. So E K(x, y) = ½J₀(mτ) **exactly at every density**: the
continuum retarded propagator of a free scalar of mass m in 1+1 (Johnston, eq. 3.23).

> **Does a single sprinkled order carry Johnston's propagator, with its error falling at the expected rate ρ^{−1/2},
> or does the propagator appear only in the average over many orders?**

Q.1 tests the numerical realisation of a propagator on a causal order. It does not test quantum mechanics as a whole.
Whether interference differs from the classical non-negative flow is the question of Q.2.

## 1. What is held fixed, and what scales

- **Region, fixed:** the light-cone square u, v ∈ [0, 1] (u = t + x, v = t − x; dt dx = du dv/2). The two fixed elements
  are x = (0, 0) and y = (1, 1): τ_xy = 1, V = 1/2. The same region at every density.
- **Mass, fixed:** m (in units of 1/τ_xy) ∈ {2, 5, 10}. The reference is ½J₀(m) = +0.11195, −0.08880, −0.12297.
- **The only varied quantity:** the density ρ = 2 N_exp, with N_exp ∈ {1000, 4000, 16000} the expected number of
  interior elements. "Scaling with N" means scaling with ρ in a fixed region at fixed m. The discreteness parameter
  m²/ρ falls with it (0.050 at most).
- **Construction, fixed:**
  - Each realisation is an independent Poisson sprinkling: N ~ Poisson(N_exp), positions uniform in (u, v). There is no
    conditioning on N.
  - Relation: z ≺ w iff u_z < u_w and v_z < v_w.
  - K(x, z) for every element z from the forward recursion k_j = a(C_xj + b Σ_{i<j} k_i C_ij), in the order of u.
    K(x, y) is the last entry.
- **Realisations:** R = 2000 per (m, N_exp). Seeds derive_seed(base, "Q.1", m, N_exp, r), independent of the
  development seeds.

## 2. One realisation and the ensemble: the metrics

For realisation r at a given (m, N_exp):

- **e_r = K_r(x, y) − ½J₀(m):** the error of one order at the fixed pair.
- **E_r = rms over the elements z of r with τ_xz ≥ 1/2 of K_r(x, z) − ½J₀(mτ_xz):** the error of one order over its
  own elements. This is a property of the single realisation, not of the ensemble.

Two judged metrics per (m, N_exp):

- **M1 = sd(e_r)** over the R realisations (ddof 1): the spread between orders at the fixed pair. Since E e_r = 0
  exactly, this is also the rms error of one realisation there.
- **M2 = median(E_r)** over the R realisations: the typical error of one realisation over its own elements.

The slope of each metric:
- s = the least-squares slope of ln M against ln N_exp over the three densities.
- Its uncertainty comes from a bootstrap. The R realisations are resampled with replacement, independently at each
  density, B = 2000 times, giving the median and the 2.5–97.5 % interval [lo, hi].

**Reported, not judged:** median |e_r| and q90 |e_r| with their slopes, and the full distributions of e_r and E_r
(stored in the record).

## 3. Decision bands (complete)

For each of the six pairs (m, metric), with [lo, hi] the bootstrap interval of its slope:

| band | condition |
|---|---|
| in band | −0.65 ≤ lo and hi ≤ −0.40 |
| fail | lo > −0.25 (the whole interval above −0.25: no self-averaging) |
| deviation | hi < −0.65 (faster than the expected ρ^{−1/2}) |
| inconclusive | anything else (the interval reaches into −0.40 … −0.25, or straddles −0.65 or −0.40) |

The verdict:
- **INSTRUMENT_FAILED** if any instrument of section 4 fails. No physics statement is made then.
- **FAILED** if any pair is "fail".
- **PASSED** if all six pairs are "in band".
- **UNDECIDED** otherwise. The record names the sub-case:
  - "UNDECIDED (deviation)" if any pair is "deviation" and none fails: not a failure, a departure from the expected
    scaling, explained first;
  - "UNDECIDED (inconclusive)" otherwise.

The bands are not widened after the run. A run that lands in "inconclusive" stays inconclusive. A new question would
need a new card (Q.1b, …).

*Note on the name:* the core of the laboratory knows the statuses PASSED, FAILED, UNDECIDED, INSTRUMENT_FAILED and
RECORDED (core/testing.py). The author's INCONCLUSIVE is recorded as UNDECIDED with the sub-label above. A separate
status would be a change of the core.

## 4. Calibration (instruments; kept apart from the verdict)

INSTRUMENT_FAILED if any fails:

- **I1 (the recursion):** K(x, ·) from the recursion against the dense solve aC(1 − abC)⁻¹ on 20 orders of N = 50–300,
  max |difference| ≤ 1e-12.
- **I2 (the exact mean):** at each of the nine points, |mean_r e_r| ≤ 3.5 SE (SE = M1/√R). This is a check of the
  implementation against the exact expectation, not evidence that one order carries the propagator.
- **I3 (the null, non-negative weights):** the same sum with b = +m²/ρ, at m = 5, N_exp = 4000, R = 2000. The mean is
  ½I₀(5) = 13.620 (positive, no zeros); |mean − ½I₀(5)| ≤ 3.5 SE. This checks the sign and the character of the
  solution: with non-negative weights there is no oscillation.

## 5. Reading declared now

- **PASSED:**
  - One sprinkled 1+1 order carries Johnston's propagator. Its error falls as ρ^{−1/2} at fixed mτ, both at a fixed
    pair and over the order's own elements.
  - The instrument is usable on the order that is realised, not only on average.
  - This is a calibration of the quantum section's first instrument. It has no physical conclusion of its own.
  - It says nothing about 3+1 (Q.1b, links, where the mean is not exact), interference (Q.2), the Born rule (Q.3), or
    masses near the discreteness scale (m²/ρ → 1).
- **FAILED:** a single order does not carry the propagator; it appears only in the average. Q.2–Q.4 would then need an
  ensemble or another instrument, decided before they are designed.
- **UNDECIDED:** reported with its sub-case; no conclusion.

## 6. Development diagnostics (not a result of Q.1)

These are the records diag/out/q1_dev_20261128T093204Z and q1_dev_20261128T094112Z: 400 realisations,
diag/q1_dev.py, other seeds. They show that the design is feasible; they are not part of the test.

| m | M1 slope (bootstrap 95 %) | M2 slope (bootstrap 95 %) | q90 slope | M1 at N_exp 16000 | M1 / \|½J₀\| | M2 at 16000 |
|---|---|---|---|---|---|---|
| 2 | −0.500 [−0.535, −0.469] | −0.502 [−0.531, −0.470] | −0.503 | 0.00237 | 2.1 % | 0.00171 |
| 5 | −0.530 [−0.566, −0.494] | −0.517 [−0.537, −0.494] | −0.556 | 0.00533 | 6.0 % | 0.00415 |
| 10 | −0.499 [−0.536, −0.463] | −0.513 [−0.529, −0.499] | −0.508 | 0.00882 | 7.2 % | 0.00924 |

- I2 in development: |z| ≤ 1.43 at all points.
- From N_exp = 250, the adjacent slopes are −0.47 … −0.68 (pre-asymptotic at m = 10).
- **The relative error grows with m** at fixed density (2.1 % → 6.0 % → 7.2 % at N_exp 16000). This matters for the
  range of validity: the instrument is controlled by m²/ρ, and at a given density heavier fields are carried less
  precisely by one order.
- With R = 2000 the bootstrap intervals are expected to be about √5 narrower (±0.015).
- Cost: a full run is a few minutes.
