# F6.8f-dev: the optical response at finite amplitude (development log)

Opened 16 Nov 2026 on the author's plan (docs/derivations/substrate_and_action_audit.md, section 8; roadmap, working
order 1; "Ano prosím", 10:22, for step 0). **Development: it decides nothing, it is not a test and it has no card.**
None of it changes F6.8e3 or F6.8d.

**Question of the branch.** Does the expected BD action give a usable static response at finite amplitude b and large
age T? Quantity: E_opt(k, T, b) = ∂E[S_BD]/∂b and K_eff = [E_opt(b) − E_opt(0)]/b, with the optical projection of
F6.8e3. The second variation holds only for ε² = b²(kT)⁴ ≪ 1.

**Steps:**
- **Step 0:** the exact first-arrival time at finite b.
- **Step 1:** V(b)/V0 of one long thin interval against ε at fixed kL. The b² coefficient must give V'' of F6.8e3.
- **Step 2:** the density with the BD kernel at finite b, and the precision of the survivor.
- **Then:** the grid in ε (to ε ≈ 20–40, beyond the first caustic at π²) over several (b, k, T).

## Step 0: the exact first arrival (16 Nov)

**Code:** diag/f68f_dev.py. **Output:** diag/out/f68f_step0.txt.

**Method.**
- Medium n = e^{b cos(kz + θ)}: the optical sector of F6.8c, g = −dt² + n² dx².
- The medium is stratified, so Snell's invariant q = n sin θ is conserved along a ray. The travel time and the
  transverse distance are 1D integrals in z: t = qX + ∫√(n² − q²) dz and X = ∫ q/√(n² − q²) dz.
- **Ray families:**
  - monotone rays for q < e^{−|b|};
  - bouncing rays inside the allowed layer around a maximum of n, with their phase unfolded (start up or down, end up
    or down, number of half-periods).
- The first arrival is the smallest time over all families that reach the point. Caustics and multiple arrivals are
  included exactly, with no expansion in b.
- dt/db is taken along the arriving ray (Fermat): ∫ cos(kz + θ) n dl.
- **Precision devices:**
  - turning-point integrals in the variable u = m − a cos ψ, which removes the end singularities;
  - cos u − cos a written as a product of sines;
  - n² − q² through expm1;
  - the families parametrised by the half-width a of the layer and by λ = −ln(q e^{|b|}), not by q.
  - Double precision; long double is not needed at this step.

**Checks:**

| Check | What | Result |
|---|---|---|
| T1 | b = 10⁻⁹ against the flat distance | differences 1.3e-9, 4.0e-9, 1.7e-9, all first order in b |
| T2 | second-order coefficient, phase-averaged (16 phases, b = ±10⁻⁴), against the closed form t2_mean of F6.8c | 7.5e-10, 1.8e-7, 2.6e-8, 1.0e-7 (d = 3, 5, 2, 40), at the level of the finite difference |
| T3 | b = 0.3 and 0.6, against the polyline minimisation of F6.8c (n = 200 and 400, Richardson) | 6.8e-12, 3.2e-11 |
| T4 | the first caustic from the layer of largest index (b = 10⁻², k = 1) | one arrival (the axis ray) for X < π/(k√b); three from X = 1.01 π/(k√b), the bounced ones earlier than the axis ray (−1.3e-4, −1.1e-2, −0.18 at X/X_c = 1.01, 1.1, 1.5) |
| T5 | dt/db along the ray against a central difference in b | 3.1e-10, 3.0e-10 |

**Reading:**
- The instrument for the causal relation at finite b exists. It reproduces the second-order optics of F6.8c and the
  finite-b ray minimisation of F6.8c to their accuracy.
- The first caustic sits where the paraxial estimate puts it, b(kL)² = π², for rays along the layer.
- What is open for step 1: the volume of an interval needs the first arrival at every point of the diamond, with a
  precision that survives the BD moments (~10⁻¹³ relative at second order). Whether that needs long double, and how the
  thin-diamond geometry is best sampled at ε ≳ 1, is the work of step 1.

## Step 1: design note before computing (16 Nov, for the author)

Three points came up while setting up step 1. They change how the step should be done, so they are written down
before any number.

1. **At fixed flat τ the per-pair volume is not meaningful for long intervals at finite b.**
   - The optical distance between x and y shifts at first order by b T1 with T1 ~ L. This moves the causal threshold
     of the pair: τ² changes by about 2L b T1, far more than τ² ~ 1 itself. At fixed coordinates a long pair can leave
     the causal relation before ε reaches 1.
   - In the s-integral at a fixed spatial separation, ds = dΔt² = dσ² with σ² = Δt² − T_xy(b)². The shift of the
     threshold therefore drops out exactly, apart from the age cut at the upper end.
   - For g = −dt² + n² dx² (ultrastatic), σ is exactly the proper time of the pair along its geodesic.
   - **So the per-pair object at finite b is Ṽ(σ) = V at Δt = √(T_xy² + σ²).** The b² check against γ of F6.8e3 at
     fixed τ is meaningful only for short intervals (L = 2, 5, 20, as I1 of F6.8e).
   - **Correction (same day, 0′ below):** F6.8c already parametrises the past this way. In its header: "s = l2 (dt² −
     T(x,y)²) (fixed cone)", with dt' and dt'' in l', l'' and the measure m = 1/(2 l2 dt) giving (ln m)' and (ln m)''.
     The γ of F6.8c/e/e3 is therefore already at fixed σ, and the b² check against it is direct at any L. What is new
     at finite b is only that T_xy, dt and m must be taken exactly, not to second order.

2. **The structure of Ṽ/V0 and what survives the BD moments.**
   - Ṽ/V0 = A + B σ² + C σ⁴ + D σ⁶ + …. The moments remove A, B and C; the survivor is D (and higher).
   - At second order A − 1 ∝ ε² is the "τ-independent part ∝ (kL)⁴" of F6.8e.

3. **A route that avoids the 10⁻¹³ requirement (an argument, to be checked).**
   - In the rest frame of a long pair (boost L/σ) the curvature components along the pair grow by (L/σ)². The pair then
     sees curvature ~ b k² L²/σ², so ε = b (kL)² is its curvature on the scale of its own size.
   - For L → ∞ at fixed ε and fixed phase advance p, the geometry around the central ray tends to its Penrose limit: a
     plane wave whose profile along the ray carries the wave of n.
   - The causal structure of a plane wave is exact through linear Jacobi equations, with conjugate points (caustics)
     where they are. So A(ε, p) follows from linear ODEs, exactly in ε.
   - The σ-dependence (B, C, D) comes from the transverse structure of the wave on the scale σ, i.e. in powers of
     (kσ)². D can be computed as a perturbation in (kσ)² around the exact plane wave. That gives the survivor directly,
     without subtracting large numbers.
   - **Check of the route:** at small ε it must give the survivor of F6.8e3 (I ∝ k⁵ (kr)⁵).

**Proposed order of step 1:**
- **1a:** the exact volume by brute force for short intervals at fixed τ against γ of F6.8e3, to check the volume code.
- **1b:** A(ε, p) for long intervals in the plane-wave limit up to ε ≈ 40, against brute force at moderate ε. This shows
  where the nonlinearity starts and what the caustic at π² does to the volume.
- **1c:** D at finite ε, by the (kσ)² expansion, checked against F6.8e3 at small ε.

## Step 0′: σ, the Jacobian and the BD moments with the age cut (16 Nov, the author's order of 10:34)

The author asked for two things before 1b: whether σ is invariant, and an analytic check of what the BD moments remove
when the age cut is finite.

### σ and the change of variable

**Derived here: σ is the Lorentzian distance of the pair.**
- In the optical sector, g = −dt² + h, with h = n² dx² (ultrastatic, g_tt = −1).
- Take a timelike curve from x to y. Its spatial path has h-length ℓ and it takes the time Δt. Its proper time is
  ∫√(1 − v²) dt, with v the h-speed. By concavity this is ≤ Δt √(1 − (ℓ/Δt)²), with equality at constant speed.
- The supremum over curves takes the shortest spatial path, ℓ = T_xy, the h-distance (the first arrival, the global
  minimum, so caustics are included). Hence sup τ = √(Δt² − T_xy²) = σ.
- The Lorentzian distance (the supremum of proper time over causal curves) is a coordinate-free quantity. Its expression
  through t and T_xy uses the Killing time t of the static field.
- This holds in the optical sector, where Ω = 1. With a conformal factor Ω ≠ 1, σ stays a valid integration variable
  (T_xy is the optical distance, which the light cones fix), but it is no longer the proper time.

**The Jacobian.**
- s = l2 σ² at fixed spatial end points gives ds = 2 l2 Δt dΔt. The pair measure dΔt becomes m ds with
  m = 1/(2 l2 Δt), where Δt = √(T_xy(b)² + s/l2).
- m depends on b through T_xy. This is part of the integrand (the terms (ln m)' and (ln m)'' of F6.8c), not an artefact
  of the parametrisation. At finite b it must be taken exactly.

### The BD moments (computed)

M_j[g] = ∫₀^∞ s^j g(s²) ds for the three kernel functions of the expansion (kernel_derivs of F6.8a, BD coefficients −1,
9, −16, 8; w = s², s = l2 τ²):

| j (power τ^2j) | 0 | 1 | 2 | 3 | 4 |
|---|---|---|---|---|---|
| g0 = g_k | −1.6e-16 | −3.2e-16 | −4.3e-16 | +0.500 | +2.659 |
| g1 = w g_k' | +4.9e-16 | +8.9e-16 | +2.1e-15 | −1.000 | −6.647 |
| g2 = w² g_k'' | −2.0e-15 | −3.9e-15 | −9.9e-15 | +3.000 | +23.26 |

- The parts τ⁰, τ², τ⁴ of each coefficient vanish to rounding, for all three kernel functions; τ⁶ is the first that
  survives.
- With m ≈ 1/(2 l2 L)(1 − s/(2 l2 L²) + …) for a long pair, the s-dependence of m moves a coefficient to higher
  moments with factors 1/L².

### The age cut, exactly

The past of x is cut at Δt ≤ T. With the flat cut (b = 0 in the cut), s ≤ l2 (T² − r²), i.e. r ≤ R(s) = √(T² − s/l2).
Exchanging the order of the integrals, per direction:

  J(T) = ∫₀^∞ ds ∫₀^{R(s)} r² dr m F(s, r),     dJ/dT = ∫₀^{l2 T²} ds (T/R) R² m F(s, R(s)),

with F = Σ g_i(s²) P_i(s, r) the second-order integrand. This is exact.

- **F6.8e/e3 computed only the leading term,** the bulk density at r = T over the full s range: ∫ ds T² m F(s, T).
  The rest is the layer of pairs with r close to T: expand in R − T = −s/(2 l2 T) − s²/(8 l2² T³) − ….
- **The layer term of order n** carries (R − T)ⁿ ∂ⁿ_r [r m F]. It adds n powers of s, with a factor T⁻ⁿ from (R − T)ⁿ
  and about T⁻ⁿ from ∂ⁿ_r (the coefficients are powers of r).

**Answers to the author's four questions:**

| Question | Answer |
|---|---|
| Which terms cancel exactly? | In the bulk at full s range, A, B and C (τ⁰, τ², τ⁴) for all three kernel functions, to rounding. Truncation at s_max is exponentially small (g ∝ e^{−s²}) except in the layer r ≈ T. |
| Which cancel only as T → ∞? | In the layer: C first contributes at n = 1 (relative ~ T⁻²), B at n = 2 (~ T⁻⁴), A at n = 3 (~ T⁻⁶), each times its own size. |
| Does the age edge produce a contribution from A − 1? | Yes, but only at n = 3: about (A − 1) M₃ s³/T⁶ against the bulk survivor D M₃. A − 1 is ~10⁹ D at r = 400 (F6.8e card), and T⁶ = 4e15 at T = 400, so the ratio is ~10⁻⁷ there and falls further with T. |
| Does that contribution give the optical term of F6.8e3? | No. The F6.8e3 term is the n = 0 bulk density, i.e. the τ⁶ survivor D (and higher), with A, B, C removed by the moments. The edge corrections are a separate, decreasing piece. |

**What is not settled analytically:**
1. **The b-dependence of the cut.** With b in the cut, R(s) is defined through T_xy(b). The region in r moves with b,
   which is the explicit boundary term of F6.8c. It is ≤ 2.4e-3 and falls as T⁻⁴ in the k = 0 null of F6.8c, but it is
   not computed at large n.
2. **The size of the C and B layer terms at moderate n.** At n = 24 (k = 0.12) the record shows bulk against full rows:
   1.19 (F6.8e development) or 1.109 (F6.8e3 report) against 1.244 (F6.8c2 full rows). That is a 4–11 % difference,
   so the layer is not negligible there.
   - Proposed check: compute dJ/dT exactly at second order with the stable instrument at r = R(s), at n = 24, 96 and
     384. It should approach the bulk value as 1/n² if the analysis is right, and its value at n = 24 should explain the
     4–11 %.

**For 1b and 1c:**
- A(ε, p) does not create the response. It is removed by the moments and enters the edge only at n = 3.
- 1b is still useful, as the map of where the nonlinearity and the caustics start in the cross-section, and as the
  zeroth order of 1c. But the response is in D (and higher).
- **1c and 1d carry the result.** The quantity to compare with F6.8e3 is the full s-integrated bulk density (all g_i
  with their moments) at small ε, as the author asked: not D or A alone.

### Step 0′, numerical: the layer of the flat age cut (16 Nov; diag/f68f_layer.py, diag/out/f68f_layer.txt)

dJ/dT computed exactly at second order (flat cut) with the stable instrument of F6.8e3 (lev1, the same s and p nodes as
its pair_density), against its bulk term alone. k = 0.12; age slope −a₄ρ/k⁵ in units of k² per unit n. Cost 3 972 s on
2 processes. Development, decides nothing.

| n | bulk (r = T, as F6.8e3) | exact with layer | layer / bulk | (layer / bulk) · n² |
|---|---|---|---|---|
| 24 | 1.109022 | 1.489463 | +0.343 | 198 |
| 96 | 137.4478 | 142.4668 | +0.0365 | 337 |
| 384 | 33 320.06 | 33 399.55 | +0.00239 | 352 |

- **The 1/n² law holds asymptotically.** (layer/bulk)·n² changes by 4.5 % between n = 96 and 384 and approaches a
  constant c ≈ 3.5e2. At n = 24 the higher layer terms are still large (198 instead of ~350). So the layer is a relative
  correction ≈ c/n²: 0.24 % at n = 384, below 1e-4 from n ≈ 2 000.
- **Effect on the exponent of F6.8e3.** Between n and 2n the layer shifts the local exponent by
  log₂[(1 + c/4n²)/(1 + c/n²)] ≈ −1.08 c/n²: −1.0e-2 at n = 192, −1.6e-4 at n = 1 536, −4e-5 at n = 3 072. The exponent
  5 is the large-n limit and is unchanged; the reported values α = 4.9841 … 4.9999 would move down by these amounts.
- **The prediction for n = 24 failed.** I expected the layer to explain the 4–11 % gap to the F6.8c2 full rows (1.244).
  It moves the value from 1.109 to 1.489, past 1.244, so the layer alone does not explain the gap.
- **The comparison at n = 24 is not like for like.** The F6.8c2 number is a secant: (K̂_bb(n = 16) − K̂_bb(n = 24))/(8 k²)
  = (0.741976 − 0.598664)/(8 · 0.0144) = 1.244, and it contains the explicit b-dependent cut. The values here are the
  derivative at n = 24 with the flat cut.
- **The secant with the layer, first version: wrong at n = 16 and 20.** I used the p panels of F6.8e3, [0, 24], for every
  n. p = n cos(k, pair) runs only over [0, n], so for n = 16 and 20 the sum included directions with cos > 1 (kperp
  clipped to 0). It gave 0.649 and 0.984 (bulk 0.478 and 0.708) and a Simpson secant of 1.012. These numbers are void;
  n = 24, 96 and 384 are not affected. Rerun with panels ending at p = n (diag/f68f_layer.py, p_panels).
- **The explicit boundary term is already on record and is small.** F6.8c2 stores it per row: bnd = −1.3706e-5
  (n = 16) and −1.8363e-5 (n = 24) in J_bb. In the secant that is a₄ · 4.656e-6/(8 k⁴) = **0.0046** of the 1.244.
  The interior of F6.8c2 (bulk and layer, flat cut) gives 1.2395. So the boundary does not carry a 19 % gap. The
  comparison is between the interior of F6.8c2 and the bulk plus layer here, after the rerun.
- **The rerun with p ≤ n** (diag/out/f68f_layer_secant.txt; n = 24 from the cache):

  | n | bulk | exact with layer |
  |---|---|---|
  | 16 | 0.459039 | 0.664136 |
  | 20 | 0.707912 | 0.985576 |
  | 24 | 1.109022 | 1.489463 |

  The Simpson secant is 0.733 (bulk) and **1.016 with the layer**, against 1.2395 for the interior of F6.8c2. The
  wrong p range had barely mattered (1.012 → 1.016).
- **A gap of 0.224 (18 %) remains between the two implementations.** It is not identified and is kept as such, not
  given to the boundary. Next: nodes and quadrature of both, starting with the p (direction) quadrature of the
  F6.8e3 instrument at n = 16 and 24. F6.8c2 has 36 Gauss nodes per half of cos ∈ [−1, 1] at n = 24; the F6.8e3
  panels have 16. The F6.8e development row (1.19 at n = 24 against 1.109) already showed that the n = 24 value
  depends on the quadrature.
- **For F6.8e3 nothing changes:** the layer falls as 1/n² and the exponent is a large-n statement. **For the static
  limit it matters:** at finite age the layer and probably the boundary are not small, so a finite-age response needs
  the full cut, not the bulk.

## Step 1a: the exact volume of short intervals at finite b (16 Nov; diag/f68f_1a.py, diag/out/f68f_1a.txt)

**Instrument.** V(b, θ) = ∫ d³w n³ [dt − T(x,w) − T(w,y)]₊ with T the exact first arrival:
- each segment is minimised by Fermat (chord deviation u(1−u) Σ c_j T_j(2u−1), 14 terms, Newton);
- against the exact Snell first arrival of step 0 the agreement is 0 to 5e-16 (six random segments, b up to 0.16);
- prolate coordinates around the chord, the root of the boundary on every ray (safeguarded secant on a bracket);
- the θ average on 8 phases; flat V0 reproduced to 1e-15.

Fixed σ (dt = √(T(x,y)² + τ²), as F6.8c/e3) and fixed τ (dt = √(L² + τ²)). R(b) = 2(⟨V⟩ − V0)/(b² V0) → γ; γ
extrapolated from b = 0.01 and 0.02 (R = γ + c b²). Development, decides nothing.

| Row (k, CO, L, τ) | γ fixed σ: exact / F6.8e3 | γ fixed τ: exact / second order with dt′ = dt″ = 0 |
|---|---|---|
| 0.7, 0.4, 2, 1.2 (the test diamond of F6.8c) | −0.0159487 / −0.0159484 | 19.26881 / 19.26882 |
| 1, 0, 2, 1 | −0.089168 / −0.089202 | 48.98888 / 48.98885 |
| 1, 1, 2, 1 | 0.1213981 / 0.1213981 | 28.85486 / 28.85487 |
| 1, 0.7, 4, 1.5 | 0.19616 / 0.19600 | 74.22884 / 74.22868 |

β (fixed σ and fixed τ) agrees to 6 digits in every row.

**The two rows with a visible difference** (3e-5 and 1.6e-4 absolute, the same at fixed σ and fixed τ, so in the
main part of V″, not in the parametrisation) are the quadrature of the exact volume:
- the second-order instrument is constant to 9 digits between its levels (24,16,16), (32,24,24) and (48,32,32);
- the exact volume moves towards it with its levels (fixed τ, diag/f68f_1a_levels.py):
  - row 2: 48.988877 → 48.988855 → 48.988849, against 48.988846;
  - row 4: 74.228837 → 74.228725 → 74.228698, against 74.228682.

The difference falls by a factor 3–4 per level.

**Readings:**
- **The second-order instrument of F6.8e3 is confirmed** against a computation with no expansion in b and no closed
  forms of the optics, on four short intervals including k ⊥ the pair, k ∥ the pair and kL = 4.
- **Fixed σ is the right variable.** At fixed τ, γ is 100–1 200 times larger: the shift of the cone, which fixed σ
  removes. |δ(τ²)|/τ² is already 5–10 % at b = 0.01.
- **Finite amplitude in short intervals.** At fixed σ, R(b)/γ − 1 is:
  - −0.2, −0.8 and −3.3 % at b = 0.04, 0.08 and 0.16 (row 1, ε = b(kL)² up to 0.31);
  - +0.3, +0.9 and +3.3 % in row 4 (ε up to 2.56);
  - below 3e-4 in row 3.

  The departure grows like b², not like ε. ε reaches only 2.6 here, so this says nothing about old fields (ε ≫ 1).
  That is 1b and 1d.
- **At fixed τ and b = 0.16 some phases give no interval** (T(x,y) ≥ dt, the pair is not timelike): R_τ is
  undefined there, outside |δ(τ²)/τ²| ≪ 1.

**Instrument faults found and fixed in development:**
- results were written only at the end, so a run cut off at the time limit lost everything; now a cache per job;
- at b = 0.16 the first bracket dt/L lay below σ = 1 when T(x,y) < L; the bracket now starts above 1.

## Step 1b: development protocol, fixed before the computation (the author, 16 Nov 17:11)

A development protocol, not a verdict on a physical hypothesis.

**Grid:**
- ε = 0, 1, 2, 4, 8, 16, 40;
- denser near the expected caustic: ε = 9, 9.5, 9.8, 9.87, 9.95, 10, 10.5, 12. This is denser sampling, not an
  assumption that the caustic lies at π²: its position is determined by the computation;
- for every point, L/σ = 4, 8, 16.

**Normalisation:** A is the normalised volume, A(0, p) = 1; A_scale(p) = |A(0, p)|.

**Criteria:**

| Quantity | Working criterion |
|---|---|
| Change between L/σ = 8 and 16, away from caustics | δ_A = \|A₁₆ − A₈\| / max(\|A₁₆\|, \|A₈\|, 1) ≤ 1e-3 |
| Numerical error of the Jacobi solution (carried into A, from independent changes of tolerance and step) | δ_J ≤ 1e-4 |
| Agreement with brute force at moderate ε (only where the brute force reaches that accuracy; otherwise both uncertainties are reported and the check is not decisive) | \|A_PW − A_brute\| / max(\|A_PW\|, \|A_brute\|, 1) ≤ 1e-3 |
| First conjugate point (from the Jacobi determinants, not from a drop of the volume) | \|ε_c⁽¹⁶⁾ − ε_c⁽⁸⁾\| ≤ 0.01 |

**Near the caustic:**
- no smooth relative convergence of A is required;
- the position of the conjugate point and the causal structure are judged instead;
- a change in the number of branches is recorded separately;
- the volume is evaluated on both sides, without smoothing.

**Exception:** only ε = 40. If δ_A ≤ 1e-3 is not met there, the result stays exploratory; it is kept, but not
called an exact limit value.

**Output per point:**
- A for L/σ = 4, 8, 16;
- δ_A (8 → 16) and δ_J;
- the position of the conjugate point and a caustic flag.

The error of the Penrose approximation and the numerical error of the Jacobi solution are reported separately.

**The corrected secant** (step 0′) is judged on its own against the interior of F6.8c2, 1.2395. The record of the
first run is not rewritten.

## Step 1b: working criteria (the author, 16 Nov 17:11–17:20; a development protocol, frozen before computing, not a verdict)

**Grid:**
- ε = 0, 1, 2, 4, 8, 16, 40;
- denser ε = 9, 9.5, 9.8, 9.87, 9.95, 10, 10.5, 12. This samples the neighbourhood of the expected caustic; it does
  not assume it is at π². Its position is computed;
- for every point, L/σ = 4, 8, 16.

**Normalisation:** A = V/V0, A(0, p) = 1, A_scale = max(|A_16|, |A_8|, 1).

| Quantity | Working criterion |
|---|---|
| change of A between L/σ = 8 and 16, outside caustics | \|A₁₆ − A₈\| / max(\|A₁₆\|, \|A₈\|, 1) ≤ 1e-3 |
| numerical error of the Jacobi solution, in A | ≤ 1e-4 (independent changes of tolerance and step) |
| agreement with brute force at moderate ε | ≤ 1e-3 in the same normalisation, where the brute force itself reaches it; otherwise both uncertainties are given and the check is not decisive |
| first conjugate point | \|ε_c(16) − ε_c(8)\| ≤ 0.01, from the Jacobi fields, not from a drop of the volume; the change of the number of branches recorded separately |

**Caustic and exception rules:**
- Near a caustic the relative criterion for A is not used; the position of the conjugate point and the causal
  structure are.
- The volume is evaluated on both sides without smoothing.
- ε = 40 only: if the convergence criterion fails there, the result is exploratory, not a limit value.

**Output per point:** A, the convergence estimate, the error of the Penrose approximation and the numerical error of
the Jacobi solution (reported separately), the position of the conjugate point and a caustic flag.

### Step 1b: the plane-wave limit, derivation (16 Nov; diag/f68f_1b.py)

**The scaling limit.**
- κ = kL → ∞ at fixed ε = bκ² and p = k∥L (so cos(k, pair) = p/κ → 0 and k⊥ → k).
- σ/L → 0 with kσ → 0: the needle of the interval is thin compared with 1/k. With s = l/L ∈ [0, 1] and ψ = kδ (the
  transverse displacement of the ray along k⊥, in radians of phase), the optical length is
  T = L + (L/κ²) S[ψ] + O(b²), with

  S[ψ] = ∫₀¹ [ψ′²/2 + ε cos(ps + θ + ψ)] ds,   ψ(0) = ψ(1) = 0.

- Stationary rays solve ψ″ = −ε sin(ps + θ + ψ) (a driven pendulum). The first arrival is the global minimum of S.
  - The ray is not straight at finite ε: kδ = ψ is O(ε). So the profile of the plane wave is taken along the bent ray,
    not along the chord.
  - Every candidate for the minimum has |ψ′(0)| ≤ ε + 2√ε, since S[0] ≤ ε and S ≥ ½∫ψ′² − ε.

**The plane wave.**
- The spacetime is ultrastatic, so the null Jacobi equation is the Jacobi equation of the optical metric along the
  ray.
- In the limit it is J″ = −ε cos(ps + θ + ψ(s)) J in the plane of k, and J″ = 0 across it (the b² terms of the
  conformal curvature vanish as b → 0).
- In Brinkmann coordinates around the ray, x = 0 and y = (U, V, 0). A point w = (u, v, X) is:
  - after x iff v ≥ ½ Xᵀ M₁(u) X;
  - before y iff v ≤ V − ½ Xᵀ M₂(u) X.

  Here M₁ = E₁′E₁⁻¹ with E₁(0) = 0, E₁′(0) = 1, and M₂ = −E₂′E₂⁻¹ with E₂(U) = 0, E₂′(U) = −1.
- Since ∫ d²X [V − ½XᵀKX]₊ = πV²/√det K, the volume is V = πV² ∫ du (det(M₁ + M₂))^(−1/2).
- In the plane of k, M₁ + M₂ = W/(E₁E₂), and the Wronskian W = E₁′E₂ − E₁E₂′ is constant and equals E₁(1). Across
  the plane, 1/s + 1/(1 − s).
- With σ² = 2UV and V0 = πσ⁴/24:

  **A = V/V0 = 6 ∫₀¹ √( s(1 − s) E₁(s) E₂(s) / E₁(1) ) ds**   (flat: E₁ = s, E₂ = 1 − s, A = 1).

**Checks and meaning.**
- At p = 0, θ = 0 (straight ray, E₁ = sin(√ε s)/√ε) the RK value equals the quadrature of the closed form to 4e-15
  (ε = 5: A = 1.28324).
- **E₁(1) = ∂ψ(1)/∂ψ′(0)** is the slope of the shooting map. E₁(1) = 0 is a conjugate point of y, a fold of the family
  of rays from x, where two rays are born or die.
- **The branch count is computed independently of the Jacobi fields** (roots of the shooting map, with every fold
  interval resolved). Positivity of E₁ on (0, 1] and E₂ on [0, 1) for the first-arrival ray is checked, not assumed:
  a minimising ray has no interior conjugate point.
- **A diverges only where the first-arrival ray itself has E₁(1) = 0.** Where the first arrival changes branch (the
  cut locus), A jumps.

**A point for the author: the convergence parameter of the Penrose limit is not L/σ alone.**
- The corrections are in (kσ)² (transverse width against the wavelength) and in b = ε/κ² (the ray's own b² terms).
- L/σ = κ/(kσ), so the plane-wave limit needs L/σ ≫ κ ≫ √ε. At ε = 40, κ ≥ 20 is needed (b ≤ 0.1), and then
  L/σ = 4, 8, 16 means kσ = 5, 2.5, 1.25, far from the limit.
- **Proposal** (to be decided before the finite-σ run): fix κ so that b ≤ 0.01 (κ = √(100ε), at least 10). Take
  kσ = 1/2, 1/4, 1/8, i.e. δ = (kσ)² = 1/4, 1/16, 1/64. Apply the criteria of the author to kσ = 1/4 → 1/8 in place
  of L/σ = 8 → 16, and report the b-dependence separately (two κ).
  In my message I wrote "κ ≳ 20 at ε = 40"; that was wrong (κ = 20 gives b = 0.1). κ ≥ 10√ε = 63.2.

**Revised working criteria (the author, 16 Nov 17:42; replace the L/σ form above; a development protocol, not a
verdict).** κ = kL, q = kσ, ε = bκ², L/σ = κ/q. The limit needs κ → ∞ (b = ε/κ² → 0) and q → 0 independently.

| Level | κ | b |
|---|---|---|
| 1 | 10√ε | 0.01 |
| 2 | 20√ε | 0.0025 |
| 3, only if 1 and 2 disagree | 40√ε | 0.000625 |

Two levels do not prove asymptotic convergence in b.

| Check | Criterion |
|---|---|
| volume, q = 1/4 → 1/8 | δ_A = \|Ā(1/8) − Ā(1/4)\| / max(1, \|Ā(1/8)\|, \|Ā(1/4)\|) ≤ 1e-3; Ā the mean over the same 128 phases θ_j = 2πj/128, compared pairwise phase by phase |
| numerical error of the Jacobi fields and of the volume integral | δ_J ≤ 1e-4, by a change of step and a tighter tolerance; the quadrature error reported separately and included |
| brute force at moderate ε | δ_brute ≤ 1e-3 in the same normalisation, at ε = 2 and 4, outside caustics, only between comparable configurations (finite-κ, q corrections included or extrapolated on both sides) |
| conjugate point | \|ε_c(1/8) − ε_c(1/4)\| ≤ 0.01 |

- **Record with every mean:** the distribution of the pairwise differences over θ, the number of switches of the
  minimising branch, and the range of A.
- **Two different events, stored separately:**
  - ε_c, the conjugate point of a continued stationary branch (zero of E₁(1));
  - ε_cut, where the first arrival starts to be realised by another branch.

  For each ε_c, record whether it lies on the minimising branch or on another one. The minimising ray has no interior
  conjugate point, so ε_c is not required on it. Near the cut locus the smooth criterion for A is not used: the
  position of the switch and the volume on both sides are.
- **Grid:** ε = 0, 1, 2, 4, 8, 9, 9.5, 9.8, 9.87, 9.95, 10, 10.5, 12, 16, 40; p = 0, 1, 2, 4, 8; q = 1/2, 1/4, 1/8.
- **ε = 40:** exploratory if the q-convergence fails. The κ check is reported separately there, not waived.


### Step 1b: results in the plane-wave limit q → 0, b → 0 (diag/out/f68f_1b.txt, 2 463 s)

**Numerics.** Mean of A over 128 phases. The numerical error of the Jacobi solution, |A(N) − A(2N)| (N = 4 000), is
≤ 1e-10 everywhere except at ε = 9.87, p = 0 (1.9e-8, next to the divergence below). No phase has E₁ or E₂ ≤ 0
inside the interval on the first-arrival ray.

**The onset of a second ray** (bisection in ε on the 128 phases):

| p | 0 | 1 | 2 | 4 | 8 |
|---|---|---|---|---|---|
| ε of the onset | 9.8699 ± 0.0024 | 10.692 | 12.146 | 18.716 | none below 40 |

At p = 0 this is π² = 9.8696, an output of the computation. The caustic moves to larger ε with p.

**Ā (mean over θ) and its range:**

| ε | p = 0 | p = 1 | p = 2 | p = 4 | p = 8 |
|---|---|---|---|---|---|
| 1 | 0.99950 | 0.99952 | 0.99958 | 0.99976 | 1.00002 |
| 4 | 0.99189 | 0.99225 | 0.99324 | 0.99622 | 1.00034 |
| 8 | 0.96405 | 0.96608 | 0.97132 | 0.98491 | 1.00133 |
| 9.87 | 1.57529 (max 85.5) | 0.94056 | 0.95288 | 0.97696 | 1.00202 |
| 12 | 0.86561 | 0.87276 | 0.89868 | 0.96558 | 1.00295 |
| 16 | 0.79911 | 0.80395 | 0.82034 | 0.93408 | 1.00514 |
| 40 | 0.64255 (0.627–0.684) | 0.64456 | 0.65120 | 0.69196 | 1.02802 (0.730–1.366) |

**Readings (development):**
- **Small ε:** Ā − 1 ∝ ε² (p = 0: −5.0e-4 ε²). The sign changes between p = 4 and p = 8. This is the plane-wave
  second-order coefficient; its link to γ of F6.8e3 in the needle limit is part of 1c.
- **p = 0, θ = 0, ε = π²:** the straight ray is both the minimising ray and conjugate (a pitchfork, the cusp of the
  cut locus). A of that phase diverges (85 on the nearest grid point, 9.87). On both sides A is finite. At p ≥ 1 there
  is no such point on the grid: the folds that appear are generic and the minimising ray avoids them.
- **Large ε (40):**
  - Ā falls to 0.64 at small p, with a narrow range (0.63–0.68) although up to 5 rays exist.
  - At p = 8, Ā rises above 1, still with a single ray.
  - The plane-wave volume does not blow up with ε; it is the multi-ray structure that changes.
- **The column "jumps" is unreliable** (it counts jumps of ψ′(0) > 0.5 between neighbouring phases, which also occur
  with one ray). To be replaced by the switch of the minimising branch (the sign of the action difference of two
  branches), with ε_c and ε_cut stored separately, as the revised criteria ask.
- **Not yet done:** finite q, finite b, and the brute force at ε = 2, 4.
