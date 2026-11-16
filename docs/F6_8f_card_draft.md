# Card F6.8f (draft, not frozen): the optical age term of the expected BD action at finite amplitude

**16 Nov 2026, night.** Prepared on the author's order of 22:49 (point 3). Not frozen: the criteria below are
proposals for the author's decision.

**Test:** phases/F6_order_dynamics/F6_08f_bd_action_optical_finite_amplitude.py.
**Instrument:** K.18 (draft).

## 1. Question

Does the second-order extrapolation describe the optical age term of the expected BD action up to ε = b(kT)² = 4, and
how does the full term depart from it?

## 2. What the development found (docs/derivations/f68f_dev.md, steps 0–1c)

- **The age derivative at finite b, exactly.** dE[S]/dT at finite b is an integral over the events of the past of x at
  time distance T (first arrival ≤ T).
- **The scaling limit.** In the limit κ → ∞, b → 0 at fixed ε, p only near-cone pairs count. Per direction the density
  is (T/2l₂) ∫ds ⟨g₀(s²A)⟩_θ, with A the needle volume.
- **The BD moments select the q⁶ order.** μ₀ = μ₁ = μ₂ = 0, so the survivor is dE/dn ∝ k⁵ Φ(ε), with
  Φ(ε) = ∫dp ⟨ŝ₄A₀^{−2}⟩_θ. At small ε, Φ = C₂ε², which is the (kT)⁵ law of F6.8e3.
- **The static limit of an old field exists if and only if ε^{−1/2}Φ(ε) is integrable at infinity.**

## 3. What this card can and cannot decide

- **It can:** on ε ≤ 4 (the calibrated range of K.18), decide whether the second order still describes the transition
  region.
- **It cannot:** decide the readings A, B, C of the author's message of 09:59, which are statements about Φ at ε → ∞.
  That needs the multi-ray range: K.18b, with the cut loci, and a card F6.8g up to ε ~ 40 and beyond.

## 4. Quantities and proposed criteria

- **Definitions:**
  - N(ε) = Φ(ε)/(C₂ε²);
  - C₂ from ε = 0.1 and 0.2 by Richardson;
  - judged ε = 0.25, 1, 2, 4.
- **Development rows (seen, decide nothing):** only ε = 0.1, 0.2, 0.5, 3. The judged rows are not seen. The quick run
  of the test also uses ε = 0.5 and 3 only.
- **INSTRUMENT_FAILED:**
  - K.18 not PASSED;
  - b → 0 against F6.8e3 beyond 5e-3;
  - more than one ray at some phase;
  - linear term of Φ above 1e-8.
- **Hypothesis H2, the second order holds on ε ≤ 4:**
  - PASSED if |N(ε) − 1| ≤ 0.05 at every judged ε;
  - FAILED if |N(4) − 1| > 0.2 or Φ changes sign;
  - UNDECIDED otherwise.
- **Recorded, no verdict:** the local exponent a(ε) = d ln|Φ|/d ln ε between judged neighbours. This indicates the
  static limit: a < 2 is weaker than the second order, and a < −1/2 at large ε would be needed for reading A.
- **Level:** mechanism.

## 5. The author's checks of 09:59

1. b → 0 against F6.8e3: part e3 of K.18 and the C₂ comparison.
2. Zero age: Φ(ε → 0) → 0.
3. k = 0: A = 1, so Φ = 0.
4. Both signs of b: Φ is even in ε. This is exact: A(−ε, θ) = A(ε, θ + π).
5. Conformal and mixed sectors: not computed here; unchanged.
6. Stable instrument: K.18.
