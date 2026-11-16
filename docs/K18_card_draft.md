# Card K.18 (draft, not frozen): the needle volume of the optical sector at finite amplitude

**16 Nov 2026, night.** Prepared on the author's order of 22:49 (point 1). Not frozen: freezing is the author's
decision, with the md5 of the test file and of the helper files.

**Test:** phases/K_instruments/K_18_needle_volume.py.
**Helpers (frozen with the card by md5):** diag/f68f_1b.py, diag/f68f_1c.py, diag/f68f_1b_exact.py,
diag/f68f_common.py. Precedent: the frozen card of B3 listed its diag helpers with md5.

## 1. What the instrument measures

The volume of the causal interval of a near-cone pair (a needle) in the optical sector n = exp(b cos(k·w + θ)), relative
to flat space at the same σ: A(ε, p, θ; q). It is computed in the scaling limit κ = kT → ∞, b → 0 at fixed ε = bκ² and
p = k∥T, at q = kσ. F6.8f uses it for the age term of the expected BD action at finite ε. Derivations are in
docs/derivations/f68f_dev.md, steps 1b and 1c.

Three ways of computing it:
- **Plane wave, A₀:** from the Jacobi fields of the first-arrival ray (step 1b).
- **The q² series to q⁶:** from the Taylor series (order 8) of the first-arrival actions (step 1c).
- **The exact needle:** lower envelopes of the ray families, no expansion (step 1b). Used as the independent check.

## 2. Range

- ε = 0.1, 0.25, 0.5, 1, 2, 4: the author's grid. This is below the first caustic, so there is one ray at every
  phase.
- p on the panels of F6.8e3, [0, 24].
- 128 phases.
- q = 1/2, 1/4, 1/8.

## 3. Parts and criteria

| Part | Criterion | Source | Failing gives |
|---|---|---|---|
| single | one ray at every phase of the range | needed by the series; the author's Δε_c ≤ 0.01 has no object below the first caustic | INSTRUMENT_FAILED |
| jacobi | δ_J = max \|A₀(N) − A₀(2N)\| ≤ 1e-4; series A₀ against plane-wave A₀ ≤ 1e-9 | author (δ_J); proposal (A₀) | INSTRUMENT_FAILED |
| series | \|A_exact − A_series\| ≤ 1e-6 at q = 1/2, 1/4 on 12 phases; at most 2 window failures of the exact needle | proposal | INSTRUMENT_FAILED |
| qconv | δ_A(1/4 → 1/8) ≤ 1e-3 at every (ε, p) of the range | author | FAILED |
| e3 | b → 0 against the F6.8e3 instrument (an independent code): p-integrated bulk density at k = 0.12, n = 384 against 2μ₃k⁵n⁵C₂(p)/l₂⁴, within 5e-3 | proposal | FAILED |

**Verdict:** INSTRUMENT_FAILED if single, jacobi or series fails; PASSED if qconv and e3 pass; FAILED otherwise.

## 4. What is missing, and blocks the freeze

The author's criterion **"agreement with brute force at moderate ε (2, 4) ≤ 1e-3"**, with the two levels κ = 10√ε and
20√ε (b = 0.01 and 0.0025), needs the exact interval volume at finite κ and b in the real medium. That instrument is
not built. A draft was removed untested on 16 Nov.

The chord-centred brute force of step 1a cannot be used here. The ray bends by about ε/κ, while the needle has width
q/(2κ); at ε = 2, q = 1/4 that is 16 times its width, so the chord lies outside the interval.

**The author decides one of two things:**
- (a) part e3 (an independent code, but only at b → 0) replaces the brute force for this card, and the finite-κ check
  becomes a separate card;
- (b) the card waits for the finite-κ instrument.

## 5. Development rows already seen (decide nothing)

- Series against the exact needle, after the two fixes of the exact needle (window coverage, quintic Hermite):
  - three phases at ε = 2 and 4: ≤ 5e-9 (q = 1/2), ≤ 2e-8 (1/4), ≤ 1.2e-7 (1/8);
  - before the fixes: 3e-7 at q = 1/8, or a window failure.
- A₀ against the plane wave of step 1b: 1e-12. A₁ against the O(q²) expansion: 3e-8 relative.
- The finite-q expansion of 16 Nov (order q² only) gave δ_A ≤ 2e-5 away from caustics at every p up to ε = 40.
- The plane-wave numerical error δ_J of step 1b: ≤ 1e-10 on its grid.
