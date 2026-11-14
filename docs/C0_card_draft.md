# Card C0: the retarded Green function of a 3+1 order and the field it carries

**Revision 1 (tolerances from the development rows), 14 Nov 2026, ready for freezing, NOT frozen.** Decisions 25 and
26. Test **F3.16** (phases/F3_structure/F3_16_green_3p1.py, md5 3de531220105c4c332e23a88b56899fd). Helpers:
- diag/c0_card.py 1f2fc1071acc0e480a6f34a14a5af08f
- diag/c0_apex.py 12c93bbaff2d1ee384848da128c90d4d
- diag/c0_ct.py 8da25cd17f5df94e725441696efda5c2
- diag/c0_ct_renorm.py 9686e935fbd7d924f62f2ca877d305e9
- diag/c0_minvar.py 5b31241ec688e2eceec4bd326a11f63f
- diag/c0_links4.py 0f160d3dc53ef91b6f6c80624dba4a5b

The tolerances follow the rule of 10 Nov and are derived by diag/c0_card_tol.py (md5
00fa61213a3f889de5b0ae8fdcca1371) from the development rows diag/out/c0_card_dev_20261114T052221.json (seeds 13000–
13003 at ρ 30 and 60). The rule: tolerance = 1.25 × (development value + 2 SE at the card's n), rounded up.

## 1. Hypothesis

On a Poisson sprinkling of flat 3+1 spacetime, the order alone (the number n of events in the interval of a pair u < v)
with a fixed layer kernel w(n) and a local counterterm δ J(v), both derived from the continuum Green function and not
from any measurement, carries the retarded Green function of the wave operator, G = δ(t − r)/(4πr):

    φ(v) − δ J(v) = φ0(v) − c Σ_{u<v} w(n) J(u),   c = 1/(2π √6 √ρ),   δ = (a/8√π) Σ_n w_n H_n,   a = √(24/(πρ)),

and through it the linear (J = μφ) and nonlinear (J = μ sin φ) field of the continuum, with deviations that

- equal the finite-density expectation of the kernel ((A1)–(A3) of diag/c0_apex.py, no free parameter), and
- fall with the density;

and with noise that accumulates without a growing mode.

## 2. Model

- Kernel: **MV8+ct** (decision 26). MV2+ct is computed on the same pairs as a reference check and decides nothing.
- Enumeration: pairs with ρV ≤ 22, layers n ≤ 16.
- The counterterm is solved per event (Newton; 1 − δμ cos φ > 0, unique).
- Domain: pyramids X 7, T 6; densities ρ 30 (16 sprinklings) and ρ 60 (4); seeds derive_seed(base, "F3.16", ρ, k);
  base 0 decides (`python run.py F3.16`).
- Coordinates only in the readouts.

## 3. Readouts

- **ε(r)**: the kernel's sum per unit cone measure r dr, by shells r < 2.5 (sample events at t ∈ [3.5, 5.5)), against
  (A1).
- **Normalisation**: the constant-source sum against its expectation at each event's own t, S/E[S](t_v) − 1, t ≥ 3.
- **Modes**: the amplitude deviation of the mode from the exact solution, by projection on cos(kx), in slabs
  |t − b| < 0.25 of the inner events.
  - Linear k = 0, 1, 2, 3 (data A cos(kx) cos(kt), μ 1).
  - Nonlinear homogeneous A = 0.5, 1, 2 (exact: the pendulum).
  - Nonlinear A = 1 cos x (exact: 1+1 sine-Gordon).
- **Comparisons (decision 26), each reported separately**:
  - **A**, the noise-free continuum: amp-dev itself. Linear and nonlinear.
  - **(A3)**, the finite-density expectation of the kernel without noise, predicted.
  - **B**, the fluctuation-corrected mean response J_eff = μ e^(−s²/2) sin f, with s the measured sd of the same row,
    pooled over the sprinklings of the density. Nonlinear homogeneous rows.
  - **Residual** = measurement − B.
  - **Renormalization** = B − (A3).
- **Conditions of B**, for the nonlinear homogeneous rows per slab:
  - the skewness and excess kurtosis of the event noise η;
  - the direct check ⟨sin φ⟩ − e^(−s²/2) ⟨sin(φ − η)⟩, which is zero for zero-mean Gaussian noise independent of φ.
- **Noise**: the sd of the remainder by time; its growth exponent p (sd ∝ t^p, t 1.5–5.5); sd(ρ 60)/sd(ρ 30) per slab.
  The scaling s² ∝ t² ρ^−1/2 is a tested approximation, not a law (decision 26).

## 4. Criteria (revision 1; MV8+ct decides)

No tolerance grows with the measured noise (decision 26). Every condition on its own: |·| + 2 SE < tolerance on the
PASSED side, |·| − 2 SE > tolerance on the FAILED side. Sign criteria: difference − 2 SE > 0 (PASSED), + 2 SE < 0
(FAILED). Verdict order: INSTRUMENT_FAILED → FAILED → PASSED → UNDECIDED.

| level | condition | ρ 30 | ρ 60 | development (ρ 30 / ρ 60) |
|---|---|---|---|---|
| instrument | max_r \|ε(r) − (A1)\|, r < 2.5 | 0.013 | 0.008 | 0.0080 / 0.0038 |
| instrument | \|normalisation − 1\| (t ≥ 3) | 0.011 | 0.005 | 0.0060 / 0.0006 |
| instrument | max \|amp-dev − (A3)\|, linear k 0–3 and nonlinear A 0.5 | 0.065 | 0.045 | 0.014 (lin1, SE 0.035) / 0.014 |
| linear (A) | max \|amp-dev\|, k 0–3 | — | 0.09 | — / 0.062 |
| nonlinear (A) | max \|amp-dev\| against the pendulum, A 0.5 / 1 / 2 | — | 0.085 / 0.105 / 0.285 | — / 0.064 / 0.079 / 0.205 |
| nonlinear (A) | max \|amp-dev\| against sine-Gordon, A 1 cos x | — | 0.09 | — / 0.061 |
| nonlinear (B) | residual max \|meas − B\|, A 0.5 / 1 / 2 | 0.02 / 0.025 / 0.065 | 0.015 / 0.025 / 0.04 | 0.007 / 0.008 / 0.028; 0.005 / 0.004 / 0.011 |
| conditions of B | \|direct check\| | 0.017 | 0.012 | 0.011 / 0.0065 |
| conditions of B | \|skewness\|, t ≥ 1.5 | 0.34 | 0.385 | 0.25 / 0.30 |
| conditions of B | \|excess kurtosis\|, t ≥ 1.5 | 0.69 | 0.715 | 0.51 / 0.53 |
| noise | growth exponent p, every mode | 1.45 | 1.45 | 1.14 / 1.11 |
| convergence | max \|amp-dev\|(30) − (60), k 0–3 | sign | | +0.024 … +0.045 |
| noise | 1 − sd(60)/sd(30) at t 4.75, every mode | sign | | +0.08 … +0.18 |
| scaling | renormalization(30) − (60), A 1 and 2 | sign | | +0.020 ± 0.013, +0.060 ± 0.023 |

The instrument level is the agreement with the predictions of the kernel itself ((A1), (A3)). The physics levels
are the deviation from the continuum (A), the explanation of the nonlinear softening by the fluctuations (B and its
conditions), and the scaling with density and time.

Risk declared now: the renormalization of A 1 falls by +0.020 ± 0.013 in the development rows; at the card's n it may
end UNDECIDED.

## 5. Reading declared now

**PASSED:** a kernel that reads only the order (the counts in intervals), with a local term fixed by the continuum,
carries the retarded Green function of 3+1 and the linear and nonlinear field on it. Its errors are the predicted
finite-density effects and fall with the density. It is the instrument for C1 onwards (curved 3+1 orders, B4).

**FAILED:** the kernel does not carry the continuum field within tolerance, its errors do not fall with the density,
or the noise has a growing mode.

**Limits:**
- flat order only (no gravity);
- T ≤ 6 (the noise accumulates; long-lived objects need another formulation; the development rows reach T 8);
- the counterterm is a finite-density correction derived in the continuum;
- J local, μ 1;
- one kernel family.

## 6. Development rows (decide nothing)

- **ρ 30 and ρ 60, X 7, T 6**: diag/c0_minvar.py v1 and v2. See docs/findings.md, "C0 development rows", and
  docs/c0_apex_predictions.md.
- **ρ 20 and ρ 30, X 9.5, T 8**: diag/c0_ct.py, run diag/out/c0_ct_20261114T010642.json, predictions
  docs/c0_ct_predictions.md; check diag/c0_ct_renorm.py.
- **The card's configuration**: diag/c0_card.py, run diag/out/c0_card_dev_20261114T052221.json, seeds 13000–13003,
  ρ 30 and 60 (the tolerances above).

## 7. Settled before freezing

- **Normalisation**: it agrees with (A3) within the scatter between sprinklings (per-sprinkling sd 0.003–0.005 at
  ρ 30, shared by all kernels).
- **Noise readout**: it removes the expectation at each event's own t; the field modes are read event by event.
- **Nonlinear**: compared with A and B separately (decision 26). The conditions of B are measured directly; in the
  development rows the direct check is ≤ 0.011, and the skewness and kurtosis are small after t 1.5.
- **Kernel**: MV8+ct (decision 26).
