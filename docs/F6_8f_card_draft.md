# Card F6.8f: the static response K_eff of the optical sector at finite amplitude

**Frozen 17 Nov 2026 by the author's decision of 08:22** (docs/author_decisions.md, decision 34), with the md5 of the
test file and of the helpers listed there. The file keeps its draft name.

Prepared on the author's order of 16 Nov 22:49 (point 3), for the goal of 09:59. The criteria were proposed at
23:40 UTC on 16 Nov, **before any K_eff row of the development was seen** (the development run of step 2 had finished 4
of 84 jobs and printed no values). The author decided them on 17 Nov at 08:22, after the development rows of step 2
(R(0.5) = 1.014, R(3) = 1.69 at reduced resolution) had been reported to him.

**Test:** phases/F6_order_dynamics/F6_08f_bd_action_optical_finite_amplitude.py.
**Instrument:** K.18 (frozen together with this card). **Assembly:** diag/f68f_keff.py (step 2 of docs/derivations/f68f_dev.md).

## 1. Question

Does K_eff = (E_opt(b) − E_opt(0))/b, with E_opt = ∂E[S_BD]/∂b of the optical sector, keep the second-order growth
(kT)⁵, bend towards saturation, or change sign on ε = b(kT)² ≤ 4?

## 2. From the development (steps 1c and 2; decides nothing)

- **The age density.** In the scaling limit (κ = kT → ∞, b → 0 at fixed ε, p) it is dE/dn = c Φ(bn²), with
  Φ(ε) = ∫dp ⟨ŝ₄A₀^{−2}⟩_θ from the needle volumes of K.18. The constant cC₂ is the b² coefficient of F6.8e3.
  Checked: the p-profile and the constant agree with F6.8e3 to 5.5e-4 at n = 1536.
- **K_eff at fixed k and age n,** with M(ε) = ∫₀^ε √u Φ′(u) du:
  - E_opt = (c/2b^{3/2}) M(ε);
  - K_eff = (c/2b^{5/2}) M(ε);
  - K₂ = lim_{b→0} K_eff = (2/5)cC₂n⁵.
- **One curve.** All K_eff curves of (b, k, T) collapse onto R(ε) = K_eff/K₂ = M(ε)/((4/5)C₂ε^{5/2}).
  - Null model (the second-order extrapolation from F6.8e3): R = 1.
  - At fixed b: K_eff(T) = K₂(T)·R(bk²T²).
  - Local exponent: γ_T = d ln|K_eff|/d ln T = 2 d ln|M|/d ln ε. It is 5 for the second order and → 0 for a static
    limit.
- **Computation.** M = √εΦ(ε) − ∫₀^{√ε} Φ(t²) dt, by Gauss–Legendre on panels in t = √u.
  - Selftest of the assembly: Φ = u², u³ exact to 2e-16.

## 3. What this card can and cannot decide

- **The calibrated range.** K.18 is calibrated for ε ≤ 4: one ray per phase; the first caustic is at ε = π² for p = 0.
- **It can:** read the K_eff curve on that range: whether it keeps the second-order growth, bends, or changes sign.
- **It cannot:** decide whether an old field (ε → ∞) has a finite static response. That is a property of Φ at large ε.
  It needs the multi-ray range: K.18 extended with the cut loci, and a card F6.8g. **Reading A here is an indication
  on ε ≤ 4, not a proof of the static limit.**

## 4. Quantities and criteria (frozen)

**Grid.**
- Judged ε = 0.25, 1, 2, 4: the ends of the t panels [0, ½], [½, 1], [1, √2], [√2, 2], with G = 4 Gauss nodes each.
- C₂ from ε = 0.1, 0.2 by Richardson.
- Development used only:
  - ε = 0.1, 0.2, 0.5, 3;
  - the nodes of panels ending at 0.5 and 3, at reduced resolution (16 phases, 4 directions).
  - The judged rows at full resolution are not seen.
- γ_T on the judged intervals 0.25 → 1, 1 → 2, 2 → 4: γ = 2 ln(M_b/M_a)/ln(ε_b/ε_a).

**INSTRUMENT_FAILED** if any of:
- K.18 not PASSED (separate record);
- more than one ray at some phase;
- stationarity of the ray (linear term of the action sum) above 1e-8;
- |R(G) − R(G−1)| above 1e-3 at a judged ε;
- Φ with 128 against 64 phases above 1e-3 (relative) at a judged ε.

**Readings (the author's A, B, C of 09:59; thresholds by his decision of 17 Nov 08:22), checked in this order:**

| Reading | Status | Condition |
|---|---|---|
| B: K_eff grows or changes sign | FAILED | γ_T(2 → 4) ≥ 5 (no slowdown against the second order), or a confirmed sign change of M: at some judged ε both R(G) and R(G−1) are below −1e-3 (K_eff opposite to the second order beyond the quadrature tolerance) |
| A: the K_eff curve bends towards saturation | PASSED | not B; M keeps the sign of C₂ at every judged ε; γ_T(2 → 4) ≤ 2.5; the exponent is nonnegative and decreases over the three intervals: γ_T(0.25 → 1) > γ_T(1 → 2) > γ_T(2 → 4) ≥ 0 |
| C: young / weak fields | UNDECIDED | otherwise: a result for young / weak fields, without a decision on the static limit of an old field (the author's wording) |

- **"Confirmed" sign change** is Claude's operational form of the author's word: the sign of R must hold with both
  t quadratures (G and G − 1) and exceed their tolerance 1e-3.
- **A consequence of the nonnegative exponent.** A K_eff that peaks inside ε ≤ 4 and then falls (γ_T(2 → 4) < 0) reads
  C, not A. Checked on analytic Φ before the freeze: Φ = −ε²/(1+ε²)^{1.2} gives γ_T = 3.76, 1.23, −0.56 and reading C;
  Φ = −ε²/(1+ε²)^{0.8} gives 4.18, 2.62, 1.82 and reading A; Φ = −ε² gives 5, 5, 5 and reading B; a sign change gives B.

**Recorded with every status:**
- R(ε) and N(ε) = Φ/(C₂ε²) at the judged ε (null model R = N = 1).
- **The separate null test |R(4) − 1| ≤ 0.05** (whether the second order still holds at ε = 4). It is reported as its
  own PASSED / FAILED next to the reading and does not enter the status.

**Level:** mechanism.

## 5. The author's checks of 09:59

1. **b → 0:** R, N → 1 as ε → 0. F6.8e3 itself is compared in K.18 part e3.
2. **Zero age:** M ~ ε^{5/2} → 0.
3. **k = 0:** Φ = 0, computed at ε = 0 at every direction.
4. **Both signs of b:** exact. n(−b, θ) = n(b, θ + π), and an even number of phases maps the grid onto itself. So
   K_eff(−b) = K_eff(b).
5. **Conformal and mixed sectors:** no code path here; unchanged.
6. **Stable instrument:** K.18, and the G and phase checks above.

## 6. Cost

- Full run:
  - 22 + 12 values of ε × 16 directions × 128 phases of the series, ~115 s each;
  - the phase check at 64 phases.
  - About 9 h with 2 processes, 2–3 h with 8.
  - No 80-bit long double is needed (float64 only); not yet run on Windows.
- Quick run: development ε only (0.5, 3), G = 3, 8 phases, 3 directions.
