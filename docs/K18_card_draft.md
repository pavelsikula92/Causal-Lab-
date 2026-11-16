# Card K.18 (draft, not frozen): the needle volume of the optical sector at finite amplitude

**16 Nov 2026, night.** Prepared on the author's order of 22:49 (point 1). Not frozen: freezing is the author's
decision, with the md5 of the test file and of the helper files.

**Test:** phases/K_instruments/K_18_needle_volume.py.
**Helpers (frozen with the card by md5):** diag/f68f_1b.py, diag/f68f_1c.py, diag/f68f_1b_exact.py, diag/f68f_1d.py,
diag/f68f_layer.py, diag/f68f_common.py. Precedent: the frozen card of B3 listed its diag helpers with md5.

## 1. What the instrument measures

The volume of the causal interval of a near-cone pair (a needle) in the optical sector n = exp(b cos(k·w + θ)), relative
to flat space at the same σ: A(ε, p, θ; q). It is computed in the scaling limit κ = kT → ∞, b → 0 at fixed ε = bκ² and
p = k∥T, at q = kσ. F6.8f uses it for the age term of the expected BD action at finite ε. Derivations are in
docs/derivations/f68f_dev.md, steps 1b and 1c.

Four ways of computing it:
- **Plane wave, A₀:** from the Jacobi fields of the first-arrival ray (step 1b).
- **The q² series to q⁶:** from the Taylor series (order 8) of the first-arrival actions (step 1c). This is what
  F6.8f uses.
- **The exact needle:** lower envelopes of the ray families of the scaling limit, no expansion (step 1b). It checks
  the series.
- **The brute force at finite κ and b:** exact rays of the real 3+1 medium, lower envelopes, the interval volume
  integrated directly; no expansion in b, q or 1/κ (step 1d). It checks the scaling limit itself.

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
| e3 | b → 0 against the F6.8e3 instrument (an independent code): p-integrated bulk density at k = 0.12, n = 1536 against 4μ₃k⁵n⁵C₂(p)/l₂⁴, within 5e-3 | proposal | FAILED |
| brute | finite κ = 10√ε and 20√ε (b = 0.01, 0.0025) at the 8 series points with ε = 2, 4, q = 1/4, 1/8: the 1/κ² extrapolation A_ext = (4A₂₀ − A₁₀)/3 within 1e-3 of the series; the ratio of the level differences (A₁₀ − A_s)/(A₂₀ − A_s) in [3.5, 4.5] where \|A₁₀ − A_s\| > 1e-4 (below that the instrument's ~1e-6 error moves the ratio: 5.01 at a difference of 1.5e-5 in development) | author (1e-3, levels); proposal (extrapolation, ratio) | FAILED |

**Verdict:** INSTRUMENT_FAILED if single, jacobi or series fails; PASSED if qconv, e3 and brute pass; FAILED otherwise.

**Why brute compares the extrapolation.** The finite-κ volume differs from the scaling limit by a physical O(b) =
O(1/κ²) correction, not by an instrument error. The raw differences in development reach 3.7e-3 at κ = 10√ε (ε = 4,
p = 4) and 9.4e-4 at 20√ε. The configuration comparable with the series is κ → ∞. The author decides whether his 1e-3
applies to the extrapolation (as proposed) or to the raw level 20√ε.

## 4. The brute force (formerly missing; built 16–17 Nov night, diag/f68f_1d.py)

- **The geometry.** The medium is stratified along k, so the optical distance depends only on the separation along k
  (Δz), the horizontal distance ρ and the start phase. Rays are integrated as graphs z(ρ) with the exact Lagrangian
  n(z)√(1 + z′²). At each ρ node the first arrival τ(Δz) is the lower envelope over the family.
- **The volume.** V = ∫dl dX dY n³ h₊, h = Δt − T_x − T_y, with Δt fixed by σ. The Y root at each (l, X), Gauss in Y,
  a fine grid in X, Gauss in l.
- **A fault found and fixed before any use.** The Y root used plain regula falsi. At isolated (l, X) points it stalled
  with one end fixed (a = 0), which gave Y*² = root/2 (ratio 1/√2 exactly). In flat space that gave A = 0.99939
  instead of 1.0000032. Now it uses the Illinois variant with an explicit convergence error. The flat test agrees with
  exact distances to 4e-9.
- **The development rows (step 1d)** are 10 points at ε = 2 and 4, p 0–4, q = 1/8, 1/4, 1/2:
  - |A_ext − A_series| ≤ 3.7e-6;
  - the ratio of the level differences is 3.97–4.07 where the difference is above 1e-4 (the 1/κ² form), and 5.01, 3.75 at differences of 1.5e-5, 3e-5 (instrument noise ~1e-6);
  - cost 16–40 s per volume.
  - The card's 8 points (other phases) were not computed.

**Still for the author:** the freeze itself (the md5 of the test file and the six helpers).

## 5. Development rows already seen (decide nothing)

- e3 (b → 0 against F6.8e3, k = 0.12):
  - +2.9e-3 at n = 384, −5.5e-4 at n = 1536;
  - per direction, a part c/n² (+3.7e-3 at n = 384) and a remainder of about −8e-4 that does not vanish with n (the
    extrapolation n → ∞ gives 0.99922). It is smooth in p (+7e-4 at p = 0.16 to −2.9e-3 at p = 12.4). This could be
    O(k²) or the F6.8e3 quadrature; not separated.
  - Hence n = 1536 in the card: the 1/n² part is 2e-4 there.

- Series against the exact needle, after the two fixes of the exact needle (window coverage, quintic Hermite):
  - three phases at ε = 2 and 4: ≤ 5e-9 (q = 1/2), ≤ 2e-8 (1/4), ≤ 1.2e-7 (1/8);
  - before the fixes: 3e-7 at q = 1/8, or a window failure.
- A₀ against the plane wave of step 1b: 1e-12. A₁ against the O(q²) expansion: 3e-8 relative.
- The finite-q expansion of 16 Nov (order q² only) gave δ_A ≤ 2e-5 away from caustics at every p up to ε = 40.
- The plane-wave numerical error δ_J of step 1b: ≤ 1e-10 on its grid.
