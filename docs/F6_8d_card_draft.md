# Card F6.8d (draft): the age term of the mixed sector

**15 Nov 2026. Draft for the author's decision; not frozen.** The question was reformulated by the author's choice
"Přeformulovat" (13:0x), after the preparation docs/derivations/f68d_sectors_beyond_leading_order.md.

**Test**: phases/F6_order_dynamics/F6_08d_bd_action_mixed_age.py.

## 1. Question

Does the age term of the mixed (ω × ν) sector of the expected discrete BD action grow as (kT)³? If it does, its
departure from the conformal sector at large ages comes from the age, the long pairs focused by the optical index,
and is not a linear term with a non-covariant coefficient.

## 2. What led here

- **F6.8b2's two departures.** Δ = K_loc,ab − K_loc ≈ 0.46 k², and c_T,ab ≈ 0.195 against c_T ≈ 0.110 (ratio 1.78),
  both from linear fits over three ages at kT = 8–16.
- **Covariance rules out the ratio.** Any covariant quadratic form R + R F₀ R + C F₂ C, with any F₀ depending on the
  age, gives B/8 = A/6. So the ratio 1.78 is not of covariant form.
- **Not the explicit boundary term.** The boundary term of the b-dependent cut is ≤ 2.4 × 10⁻³ and carries neither
  departure (record of F6.8b2, refitted).
- **Development rows.**
  - Mixed density integrated over angle, ∝ r²: local exponents 1.50, 1.87, 1.96 at r = 200, 400, 800, k = 0.2.
  - Mixed density at fixed p = 0, in long double, ∝ r³: exponents 2.87, 2.97, 2.99, 3.00 at r = 400 … 3 200, stable
    between levels to 10⁻⁵.
  - Conformal density on the same nodes: constant.
- **Reading.** The mixed age term grows as (kT)³. With F6.8e (optical age term ∝ (kT)⁵) each factor of ν brings (kL)².

## 3. Prediction (before the deciding computation)

**Mechanism.** The first-order focusing of a long interval by the optical index grows as (kL)². It survives the BD
moments in the τ⁶ part.

**Consequences:**
- I_ab(r, p) ∝ k⁵ (kr)³ ψ(p);
- ρ_ab ∝ r²;
- the local exponent α_ab → 3.

## 4. What is computed

- **The density.** The bulk mixed pair density of F6.8b2, full s range, with the first-order pieces in long
  double. The long-double form computes σ − 1 from τ² instead of from dt/L − 1, which loses digits.
- **The grid.** k = 0.15 at n = 96 … 3 072, with p nodes 10 + 6 over [0, 24] as in F6.8e revision 1.
- **What is new.** The development reached n = 640, at k = 0.2 and p = 0 only. The judged pair n = 1 536 → 3 072 was
  not computed before this card.

## 5. Criteria

**Instrument** (as F6.8e):
- I0: long double;
- I1: pieces against F6.8b2, ≤ 10⁻⁹;
- I2: two levels, ≤ 10⁻³;
- I3: conformal c_T within 3 %;
- I4: tail ≤ 10⁻³;
- I5: k⁵ scaling, ≤ 10⁻³.

**Status** (α_ab between n = 1 536 and 3 072 at k = 0.15):

| Status | Condition |
|---|---|
| PASSED | \|α_ab − 3\| ≤ 0.1 |
| FAILED | \|α_ab − 3\| > 0.5, in particular a linear term (α_ab ≈ 1) |
| UNDECIDED | otherwise |

**Reported, not judged:**
- α_ab at every pair of n;
- c3;
- the ratio of the mixed to the conformal age slope at each n (F6.8b2 read 1.78 at n = 8–16).

**Level:** chain. Deterministic.

**Runs** in Claude's Linux workspace, as F6.8e (about 45 min).

## 6. What it does not decide

Whether the T-independent part (Δ) is covariant is not decided here. At small n it cannot be separated from the age
part, because the cut is within a wavelength of x: in the development of F6.8b2, Δ = 0.006 at k = 0.1, kT = 2. It is
reported only.

The test covers the bulk of the past only: the explicit boundary term was small at kT ≤ 16 and is not computed at large
n.
