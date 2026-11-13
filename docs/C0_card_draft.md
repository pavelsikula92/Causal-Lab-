# Card C0: the retarded Green function of a 3+1 order and the field it carries

**Revision 0 (skeleton), 14 Nov 2026, NOT frozen.** Decision 25: the card is frozen only after the nonlinear response,
the higher modes and the long-time noise have been measured (diag/c0_ct.py, the run of 13/14 Nov). The tolerances
follow the rule of 10 Nov (2 SE at the card's number of sprinklings + the largest development offset + a margin) and
are filled from those development rows. Every number marked [ ] is open. Helpers: diag/c0_apex.py (the continuum
expectation), diag/c0_minvar.py (kernels), diag/c0_ct.py (counterterm, modes, noise), diag/c0_links4.py (sprinkling).

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

- Kernel: [MV2+ct or MV8+ct, the author's choice after the night run]. The other candidate and MVA8 are computed on the
  same pairs as references and decide nothing.
- Enumeration: pairs with ρV ≤ 22, layers n ≤ 16.
- The counterterm is solved per event (Newton; 1 − δμ cos φ > 0, unique).
- Domain: pyramids X [ ], T [ ]; densities ρ [30 and 60]; [ ] sprinklings per density; seed [ ] decides.
- Coordinates only in the readouts.

## 3. Readouts

- **ε(r)**: the kernel's sum per unit cone measure r dr, resolved by the distance r of the pair (sample events in a late
  time window), against (A1).
- **Normalisation**: the constant-source sum against t²/2 and against (A3).
- **Modes**: the amplitude deviation of the mode from the exact solution, by projection on cos(kx), in slabs of time.
  - Linear k = 0, 1, 2, 3 (data A cos(kx) cos(kt), μ 1).
  - Nonlinear homogeneous A = 0.5, 1, 2 (against the pendulum) and A = 1 cos x (against 1+1 sine-Gordon).
- **Noise**: the sd of the remainder by time, with the deterministic spread of the slab removed (see 7); its growth
  exponent p (sd ∝ t^p); its density scaling.

## 4. Criteria (to be fixed from the development rows; structure only)

| level | condition | tolerance | derivation |
|---|---|---|---|
| instrument | \|ε(r) − (A1)\| per shell, r < 2.5, each ρ | [ ] | [ ] |
| instrument | \|normalisation − (A3)\|, each ρ | [ ] | [ ] (per-sprinkling sd 0.003 at ρ 30, see 7) |
| linear | \|amp-dev − (A3) prediction\|, k = 0–3, t ≤ [ ] | [ ] | [ ] |
| linear | max \|amp-dev\| (against the continuum), k = 0–3, t ≤ [ ], ρ 60 | [ ] | [ ] |
| convergence | max \|amp-dev\|(ρ 60) < max \|amp-dev\|(ρ 30), k = 0–3 | − 2 SE | [ ] |
| nonlinear | max \|amp-dev\| against the pendulum, A = 0.5, 1, 2, t ≤ [ ] | [ ] | [ ] |
| nonlinear | max \|amp-dev\| against sine-Gordon, A = 1 cos x | [ ] | [ ] |
| noise | growth exponent p of the sd, every mode | p < [ ] (no growing mode) | [ ] |
| noise | sd(ρ 60)/sd(ρ 30) at fixed t | [ ] (ρ^−1/4 = 0.84) | [ ] |
| noise | sd/A independent of A for the nonlinear modes at small A | [ ] | [ ] |

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
- The choice of the kernel (MV2+ct or MV8+ct).
- The number of sprinklings and the time limit T from the cost (ρ 60 is about 6× ρ 30).
