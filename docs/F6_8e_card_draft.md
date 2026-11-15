# Card F6.8e: the asymptotic age term of the optical sector

**15 Nov 2026. Frozen by decision 30 (10:23); revision 1 before any result (13:2x, the author's choice "Revize před
během"), see section 7.** The order of work is the author's (08:30: "B, then C,
then A, then D"; B = F6.8e and F6.8d).

**Test**: phases/F6_order_dynamics/F6_08e_bd_action_optical_asymptotic.py (md5 in the decision that freezes it).

## 1. Question

How does the optical (ν²) coefficient of the expected discrete BD action depend on the age T of a static field as
kT → ∞?

The conformal and mixed sectors weaken linearly with T (c_T = 0.110, c_T,ab = 0.195). The optical sector weakened
faster than linearly in the computed range (F6.8c2: k = 0.12, D(n) = 13.9, 14.3, 15.0, 18.0, 28.0 at n = kT = 6, 8, 10,
16, 24), and its form is open (paper, Sec. "The age term of the optical sector" and P2). It decides:
- the far-field law of a point mass (F8.1);
- the window of l_k (Part VIII);
- γ(k, T) away from k → 0.

## 2. What the development found (decides nothing)

**The numerical loss of F6.8c at large ages is real, and its cause is known.**
- For a long thin interval the second-order relative volume V''/V0 has a part that does not depend on τ and grows as
  (kL)⁴: −25.4, −348, −5 343 at L = 100, 200, 400 (k = 0.12).
- The BD moments remove it exactly from the s-integral, together with the τ² and τ⁴ parts. What survives, the τ⁶ part,
  is about 10⁻⁹ of it.
- F6.8c computes V'' as a bulk integral plus a surface integral, which cancel to O(L⁴/τ²). Its per-pair quadrature
  error (10⁻⁹–10⁻¹⁰) is the size of the survivor.
- The s-integrated density at r = 200 and 400 changed by factors of 2–4 between quadrature levels.
- This is the "numerical loss at T ≳ 170" in the record of F6.8c. The verdicts of F6.8c and F6.8c2 stand (R132). The
  rows of F6.8c2 had T ≤ 167 and passed their level check (≤ 6.7 × 10⁻⁴).

**A stable instrument.**
- (a) Coordinates that follow the bent ray, z = w + b δ(l_w), where δ is the first-order transverse deviation of the
  ray. The Jacobian is 1, so the b² coefficient is the same number, but the bulk–surface cancellation is gone.
- (b) Extended precision (80-bit long double) for the pointwise cancellation inside l''.
- (c) Axial nodes clustered at the tips of the diamond.
- Against F6.8c where F6.8c is converged: 6 × 10⁻¹³ (F6.8c's own test diamond), 5 × 10⁻¹³ and 2 × 10⁻¹² (k = 0.12,
  L = 5, 20).
- The s-integrated density is stable between levels to 10⁻¹¹ at r = 100 and 5 × 10⁻⁴ at r = 400.

**The growth.** Bulk pair density I(r, p) at the phase advance p = 0, k = 0.12:

| r | 50 | 100 | 200 | 400 | 800 | 1 600 | 3 200 | 6 400 |
|---|---|---|---|---|---|---|---|---|
| I(r, 0) | −7.40e−7 | −9.06e−7 | −6.17e−6 | −1.10e−4 | −2.95e−3 | −9.01e−2 | −2.85 | −90.9 |
| local exponent | | 0.29 | 2.77 | 4.16 | 4.74 | 4.93 | 4.98 | 4.995 |

- The conformal density on the same nodes goes as r¹ (the known linear age term). Integrated over p in the same way
  it returns c_T = 0.105 (F6.8a 0.1035).
- The p-profile has a fixed shape: width about 10, symmetric in p, below 10⁻³ of its peak at |p| = 24.
- I(k, r) = k⁵ Ψ(kr) to 10⁻³ (k = 0.12 against 0.2 at equal kr): the age part has the form k² F(kT), as in the other
  sectors.
- **Agreement with the full rows of F6.8c2.** At n = 24 (k = 0.12) the bulk density gives an age slope of 1.19 per unit
  n, in units of k². The full rows give 1.24 between n = 16 and 24, so the layer and boundary terms of the cut do not
  change the slope there.

## 3. Prediction (before the deciding computation)

- **Mechanism.** The τ-independent part grows as (kL)⁴, the second-order focusing of the cross-section of a long
  interval by the optical index. The surviving τ⁶ part carries the same factor, because the transverse structure at a
  fixed phase advance does not depend on L.
- **Consequences in sequence:**
  - I(r, p) ∝ k⁵ (kr)⁵ ψ(p);
  - the bulk density ρ(r) ∝ r⁴;
  - K̂_bb(k, T) = K0 − k²[c1 ln(1/k) + c0 + c5 (kT)⁵] for kT → ∞.
- **Exponent.** The local exponent α(n) = 1 + d ln ρ / d ln n of the age term tends to 5.
- **Excluded by a pass:** a linear age term (α = 1), as in the conformal and mixed sectors, or saturation (α < 1).

## 4. What is computed

- Two wavenumbers, k = 0.15 and 0.25, at n = kT = 96, 192, 384, 768, 1 536 and 3 072, with r = n/k.
- The bulk pair density of F6.8c (_pair_terms, full s range), with β and γ from the stable instrument.
- A profile in p on Gauss nodes over [0, 12] (12 nodes) and [12, 24] (8 nodes), doubled by symmetry.
- From these: ρ_bb(r) = (2π/kr) ∫ dp I and F'_bb(n) = −a₄ ρ_bb / k⁵.

**What is new and what is not.** Because I(k, r) = k⁵Ψ(kr), a new wavenumber gives the same exponent at the same n, so
the two k are a check of the code path, not new information about α. The new information is the range: the development
reached n = 768 (k = 0.12, p = 0 only), and the judged pair, n = 1 536 → 3 072, was not computed before this card. Only
the instrument at n = 3 072 was checked: two quadrature levels agree to 5 × 10⁻⁷, with the values not printed.

## 5. Criteria (fixed before the run)

**Instrument** (INSTRUMENT_FAILED if any fails):
- I0: extended precision available (long double eps ≤ 10⁻¹⁸). Windows numpy has none, so the test must run under
  Linux or WSL.
- I1: the stable per-pair V''/V0 against F6.8c at three diamonds, ≤ 10⁻⁹ relative.
- I2: two quadrature levels at p = 0, at the smallest and largest n of each k, ≤ 10⁻³ relative.
- I3: the same profile integration on the conformal sector gives c_T within 3 % of 0.1035.
- I4: the tail |I(24)|/|I(0)| ≤ 10⁻³ at the largest n.
- I5: k⁵ scaling between the two k at n = 192 and 1 536, ≤ 10⁻³.

**Status** (α between n = 1 536 and 3 072, at both k):
- PASSED: |α − 5| ≤ 0.1. The optical age term grows as (kT)⁵.
- FAILED: |α − 5| > 0.5 at either k.
- UNDECIDED: otherwise.

**Reported, not judged:**
- α at every pair of successive n;
- c5;
- the age n* ≈ (1/(k²c5))^(1/5) at which the optical stiffness K̂_bb changes sign (development: c5 ≈ 3.2 × 10⁻⁷, n* ≈ 42
  at k = 0.15);
- the slope at k = 0.12, n = 24 against F6.8c2.

**Level:** chain (as F6.8a–c2). Deterministic: no seed, no replication by seed; the I1 and I2 checks take its place.

**Cost:** about 270 evaluations of 30–100 s each, about 1 h with 8 processes, needing Linux or WSL.

## 6. Reading declared now

**PASSED.** The b² coefficient of the expected BD action in the optical sector grows with the age of the static field as
T⁵, not linearly. The static quadratic form of the order has no limit for old fields in this sector.
- The age term of P2, read in the conformal form, does not describe the optical sector.
- γ(k, T) and Φ(k, T) of the propagator (Eq. γ = B/2C − 1) are dominated by it beyond n*. For the window, caveat (i) of
  the paper applies in its strong form.

**FAILED.** The growth is not a fifth power: the scaling argument or the instrument is wrong.

**Limits:**
- Second order in b only. The expansion in b breaks down near b²(kT)⁴ ~ 1, and the response at finite b is not
  computed.
- Bulk of the past only. The layer and boundary terms of the age cut are not computed at large n (at n = 24 they do not
  change the slope).
- Whether the cut itself shapes the sectors is F6.8d.
- One field shape (a plane wave, phase-averaged), static.

## 7. Revision 1 (15 Nov, 13:2x, before any result)

The first full run, in Claude's workspace, was stopped at 175 of 277 tasks by the workspace's two-hour limit per
process. Nothing of it was printed or saved; the result table is written only at the end. To fit the limit:
- the profile is computed at one wavenumber, k = 0.15. The k⁵ scaling is exact, so a second wavenumber adds no
  information about α; k = 0.25 stays only in the scaling check I5;
- the p nodes are 10 + 6 instead of 12 + 8.

**Status (revised):** α between n = 1 536 and 3 072 at k = 0.15. PASSED |α − 5| ≤ 0.1; FAILED |α − 5| > 0.5;
UNDECIDED otherwise. The instrument conditions are unchanged.

Test md5 9178baf0348e46dd8b74fe758cc9edde (LF).
