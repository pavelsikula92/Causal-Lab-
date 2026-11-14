# Card C0: the retarded Green function of a 3+1 order and the field it carries

**Revision 0.1 (structure after decision 26), 14 Nov 2026, NOT frozen.** Decisions 25 and 26. The tolerances follow
the rule of 10 Nov (2 SE at the card's number of sprinklings + the largest development offset + a margin) and are filled
from the development rows of diag/c0_card.py, which runs the instrument and every readout in the card's configuration
(dev seeds 13000+; the card runs on fresh seeds). Every number marked [ ] is open. Helpers: diag/c0_apex.py (the
continuum expectation), diag/c0_minvar.py (kernels), diag/c0_ct_renorm.py (J_eff), diag/c0_links4.py (sprinkling).

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
- Domain: pyramids X 7, T 6; densities ρ 30 and 60; [16 and 4] sprinklings; seeds [ ].
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

## 4. Criteria (to be fixed from the development rows of diag/c0_card.py; structure only)

No tolerance against the continuum grows with the measured noise (decision 26): every tolerance below is a fixed
number set before the run.

| level | condition | tolerance | derivation |
|---|---|---|---|
| instrument | max_r \|ε(r) − (A1)\|, r < 2.5, each ρ | [ ] | [ ] |
| instrument | \|normalisation − 1\|, each ρ | [ ] | [ ] (per-sprinkling sd 0.003 at ρ 30) |
| instrument | max \|amp-dev − (A3)\|, linear k = 0–3 and nonlinear A 0.5, each ρ | [ ] | [ ] |
| linear (A) | max \|amp-dev\|, k = 0–3, ρ 60 | [ ] | [ ] |
| convergence | max \|amp-dev\|(ρ 60) < max \|amp-dev\|(ρ 30), k = 0–3 | − 2 SE | [ ] |
| nonlinear (A) | max \|amp-dev\| against the pendulum, A 0.5, 1, 2, ρ 60 | [ ] | [ ] |
| nonlinear (A) | max \|amp-dev\| against sine-Gordon, A 1 cos x, ρ 60 | [ ] | [ ] |
| nonlinear (B) | residual max \|meas − B\|, A 0.5, 1, 2, each ρ | [ ] | [ ] |
| conditions of B | \|skewness\|, \|excess kurtosis\|, \|direct check\| | [ ] | [ ] |
| noise | growth exponent p, every mode | p < [ ] | [ ] |
| noise | sd(ρ 60)/sd(ρ 30) | [ ] (ρ^−1/4 = 0.84) | [ ] |
| scaling | renormalization max \|B − (A3)\| (ρ 60) < (ρ 30), A 1, 2 | − 2 SE | [ ] |

Every condition on its own: |·| + 2 SE on the PASSED side, |·| − 2 SE on the FAILED side. Verdict order:
INSTRUMENT_FAILED → FAILED → PASSED → UNDECIDED.

## 5. Reading declared now

**PASSED:** a kernel that reads only the order (the counts in intervals), with a local term fixed by the continuum,
carries the retarded Green function of 3+1 and the linear and nonlinear field on it. Its errors are the predicted
finite-density effects and fall with the density. It is the instrument for C1 onwards (curved 3+1 orders, B4).

**FAILED:** the kernel does not carry the continuum field within tolerance, its errors do not fall with the density,
or the noise has a growing mode.

**Limits:**
- flat order only (no gravity);
- T ≤ [ ] (the noise accumulates; long-lived objects need another formulation);
- the counterterm is a finite-density correction derived in the continuum;
- J local, μ 1;
- one kernel family.

## 6. Development rows (decide nothing)

- **ρ 30 and ρ 60, X 7, T 6**: diag/c0_minvar.py v1 and v2. See docs/findings.md, "C0 development rows", and
  docs/c0_apex_predictions.md.
- **ρ 20 and ρ 30, X 9.5, T 8**: diag/c0_ct.py, [the night run of 13/14 Nov], predictions docs/c0_ct_predictions.md.

## 7. Open before freezing

- **Normalisation: settled (14 Nov).** It agrees with (A3) within the scatter between sprinklings. The per-sprinkling sd
  of 0.003 at ρ 30 is shared by all kernels and enters the SE of the criterion. The earlier +0.5 % was seed 0 (+0.010)
  and its neighbours.
- **The noise readout must remove the expectation at each event's own t** (S − E[S](t_v)). Without it, the deterministic
  spread of t²/2 in a slab of width 0.5 (var 0.22–0.54) poses as a kernel-independent floor. After removal, var =
  2.6 var_diag + 0.19 at ρ 30.
- **For the field modes the same holds**: the deviation from the exact solution is already taken event by event at the
  event's own t, so their sd is not affected.
- Tolerances from the development rows of diag/c0_card.py.
- The number of sprinklings from the cost (ρ 60 is about 6× ρ 30).
