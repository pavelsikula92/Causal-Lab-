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
- **The p quadrature is not the cause** (diag/f68f_pconv.py, diag/out/f68f_pconv.txt, 2 788 s): the bulk slope at
  n = 24 with 48 p nodes is 1.1090252, against 1.109022 with the 16 nodes of F6.8e3, a difference of 3e-6. Remaining
  candidates, not checked:
  - the s quadrature (16 nodes) and the levels of the volume instrument;
  - the quadrature of F6.8c2 itself (mult = 1, no mult = 2 row at k = 0.12);
  - a difference of definition between the two codes.
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

The corrected secant (step 0′) is judged on its own against the interior of F6.8c2, 1.2395; the record of the
first run is not rewritten.

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

### Step 1b at finite q (b → 0): first results (diag/f68f_1b_q.py, diag/out/f68f_1b_q.txt, 1 266 s)

**Method.** A(q) = 32√2/(πq⁴) ∫ds √(s(1−s)) ∫dψ [q²/2 − Φ]₊^{3/2}:
- the direction out of the plane of k is integrated exactly;
- Φ = G₁ + G₂ − S_min is expanded around every ray to fourth order in the transverse offset (K, c₃, c₄ from the
  variational equations of order 1–3), which gives the factor 1 + c(5c₃²/(2K³) − c₄/K²), c = q²/2 − (S_b − S_min);
- every ray with S_b − S_min < q²/2 adds its own needle.

**Checks:**
- K reproduces A_PW of step 1b to all digits.
- The exact needle (boundary-value problems solved by Newton at every node, no expansion) against the expansion,
  diag/out/f68f_1b_q_check.txt:
  - ε = 4, p = 0, θ = 0.7: −7.9e-7, −5.0e-8 and −3.1e-9 at q = 1/2, 1/4, 1/8;
  - ε = 4, p = 2, θ = 2: −5.1e-9 and below.
  - **The check stopped at the third point** (ε = 12, p = 1, θ = 0.3): the support of the needle was wider than the
    search window (fixed at 1.6 × the paraxial width). The window must be adaptive. The remaining check points were
    not computed.

**Results, δ_A = |Ā(1/8) − Ā(1/4)| / max(1, …) (criterion 1e-3):**
- **Away from caustics δ_A ≤ 2e-5** at every p and ε, also at ε = 40 (1.3e-5 to 3.9e-5). The q-convergence holds
  there. The largest pairwise phase difference is ≤ 3e-4 at ε = 40.
- **Near the caustics the expansion is not valid**, and the numbers there are not results:
  - "corr", the largest relative q² correction at q = 1/4, exceeds 0.1;
  - p = 0, ε = 9.5–10.5: A(1/2) at ε = 9.8 is negative and A diverges at ε = 9.87;
  - smaller in p = 1 near ε = 10 (δ_A up to 1.5e-3, corr 0.07) and p = 4 near ε = 16 (corr 0.018).

  These points need the exact needle with all branches, as the criteria for the caustic neighbourhood say.
- **A grid artefact at p = 0, ε ≥ 10.5.** Ā(q) is about 0.011 above A_PW for every q (ε = 12: 0.8766 against 0.8656),
  not ∝ q².
  - The cause: θ = 0 and π lie exactly on the grid. At p = 0 the two branches born at the pitchfork have exactly equal
    action there, so both needles are counted at every q.
  - The A_PW column (and Ā of step 1b) takes one of them.
  - In the continuum the double count holds only in a window of θ of width ∝ q² around the cut locus. Its weight is
    O(q²), not 1/128.
  - To fix: locate the cut loci in θ and add their contribution analytically, instead of sampling ties.

**Not yet done:** the exact needle near caustics (adaptive window, all branches), the cut-locus contribution, finite b
(κ = 10√ε, 20√ε), the brute force at ε = 2, 4.

### Step 0′: the F6.8c2 code as a derivative; the integrands and the s quadrature (16 Nov, evening)

The 18 % gap of the secant comparison (step 0′ above) was checked at one age instead of over n = 16 … 24:
- **F6.8c2 code as a derivative** (optical_row of F6.8c, mult 1, the F6.8c2 nodes; T = 197.5 and 199.9 at k = 0.12,
  equal node counts; diag/out/f68f_c2fd.txt). The central-difference age slope at n = 23.844 is 2.6219 total, 2.6134
  interior and +0.0085 boundary (k² per unit n). The interior differs from the recorded F6.8c2 row at T = 200 by the
  expected trend: K̂_bb 0.598992 at T = 199.9 against 0.598664 at T = 200.
- **The bulk + layer of this log at n = 24 is 1.489**, about 0.57 of it. The same computation at n = 23.844 was lost in
  a restart of the container.
- **The integrands are identical.** F6.8c `_pair_terms` (bb) and the integrand of diag/f68f_layer.py (F6.8e3
  instrument) agree to all printed digits at four (CO, r) × four s points. The values there are ±2e2 to ±4e2 and the
  result is a small remainder after the cancellation of the BD moments.
- **The quadratures of the F6.8e3 side are converged** at n = 24:
  - in p: 48 nodes against 16 give 1.1090252 against 1.109022 (diag/out/f68f_pconv.txt);
  - in s: N_S = 16, 24, 32, 48 agree to about 1e-5 relative at p = 0.5, 3, 10 (diag/out/f68f_sconv.txt).

**Reading (development, not a verdict on F6.8c2, R132):**
- The difference lies in the assembly or the quadrature of one of the two codes, not in the formulas.
- The candidate is the quadrature of the F6.8c2 code at mult 1. Its only level check (mult 2) was at n = 8.
- The decisive check is the same two ages at mult 2: python -m diag.f68f_c2fd --k 0.12 --T 197.5,199.9 --mult 2. It
  runs in float64, also on Windows. Open.

### Step 1b at finite q: the exact needle by lower envelopes, the caustic band (16 Nov, evening)

**Method (diag/f68f_1b_exact.py).** A(q) as in the finite-q section, with Φ = G₁ + G₂ − S_min computed from the true
first-arrival actions:
- G₁ and G₂ are the lower envelopes over the whole family of rays from x and from y (3 001 slopes, RK4 2 000 steps,
  quintic interpolation in s, cubic Hermite in ψ);
- no ray is chosen in advance;
- the ψ window is widened until Φ > q²/2 at both ends.

**Checks:**
- Against the expansion at ε = 4, p = 2, θ = 2: 1e-8 (q = 1/2) to 3e-7 (q = 1/8).
- minΦ (must not be below −1e-9) is −4.7e-11 to −2.6e-10 everywhere except one point below.

**Results** (Ā over 128 phases; δ_A = |Ā(1/8) − Ā(1/4)| / max(1, …)):

| ε | p | Ā(1/2) | Ā(1/4) | Ā(1/8) | δ_A | largest pairwise difference |
|---|---|---|---|---|---|---|
| 9.5 | 0 | 0.944737 | 0.947005 | 0.948510 | 1.5e-3 | 0.23 |
| 9.8 | 0 | 0.942161 | 0.950940 | 0.961661 | 1.1e-2 | 1.4 |
| 9.87 | 0 | 0.941976 | 0.953870 | 0.972224 | 1.8e-2 | 2.4 |
| 9.95 | 0 | 0.941720 | 0.956950 | 0.982186 | 2.5e-2 | 3.3 |
| 10 | 0 | 0.941389 | 0.957849 | 0.982021 | 2.4e-2 | 3.1 |
| 10.5 | 0 | 0.930312 | 0.927945 | 0.926008 | 1.9e-3 | 0.22 |
| 9.5 | 1 | 0.947645 | 0.947658 | 0.947677 | 1.9e-5 | 3.6e-2 |
| 9.8 | 1 | 0.942759 | 0.942771 | 0.942500 | 2.7e-4 | 1.4e-2 |
| 9.87 | 1 | 0.941514 | 0.941363 | 0.940838 | 5.3e-4 | 2.8e-2 |
| 9.95 | 1 | 0.940009 | 0.939486 | 0.938553 | 9.3e-4 | 8.1e-2 |
| 10 | 1 | 0.939013 | 0.938123 | 0.936897 | 1.2e-3 | 0.12 |
| 10.5 | 1 | 0.925201 | 0.917686 | 0.916659 | 1.0e-3 | 0.10 |
| 12 | 1 | 0.873560 | 0.872949 | 0.872798 | 1.5e-4 | 4.7e-3 |
| 10.5 | 2 | 0.944144 | 0.943353 | 0.943037 | 3.2e-4 | 1.2e-2 |
| 12 | 2 | 0.901029 | 0.899174 | 0.898793 | 3.8e-4 | 1.9e-2 |
| 16 | 4 | 0.934585 | 0.934257 | 0.934122 | 1.3e-4 | 2.9e-2 |

**Records:**
- p = 0: diag/out/f68f_1b_exact.txt and diag/out/f68f_1b_exact_e*_p0.npy, the version before the diag form.
- p = 1, 2, 4: the author's run on his machine, diag/out/f68f_1b_exact_20261116T212730Z.json/.txt (8 processes,
  1 354 s, Windows). The code md5 in the record equals the checked version of diag/f68f_1b_exact.py and diag/f68f_1b.py.

**Readings (development):**
- **p = 0 near the pitchfork (ε = 9.8–10):** no convergence in q. Ā grows as q falls, towards the divergence of the
  phase θ = 0 at ε = π², but stays finite at every finite q. By the criteria for the caustic neighbourhood the smooth
  criterion is not applied there.
- **p = 1, 2, 4:** δ_A ≤ 1e-3 except p = 1 at ε = 10 (1.2e-3) and 10.5 (1.0e-3). Both lie below the onset of the
  second ray at p = 1 (10.69), in its neighbourhood.
- **The p = 0 means contain the grid artefact of the tie at θ = 0** (both needles counted at a grid point that is
  exactly a tie). This is the same artefact as in the expansion. The sampled means at p = 0, ε ≥ 9.87 are therefore
  not the θ-averages of the continuum; the cut-locus correction below is needed.
- **ε = 16, p = 4: minΦ = −3.4e-7**, beyond the −1e-9 self-check: the envelope found an action 3.4e-7 below S_min.
  Hypothesis, not checked: the different discretisation of the family (2 000 steps) and of S_min (4 000 steps) at
  large ε. Until checked, that row is not used.

### Step 1b: the cut loci in θ (16 Nov, evening; partial)

At a cut locus θ_c the first arrival switches between rays a and b. At finite q both needles count while
|S_b − S_a| < q²/2, which adds (A_a + A_b) c / (6πα) to the θ-mean, with c = q²/2 and α = |d(S_b − S_a)/dθ|
(diag/f68f_1b_cut.py). For p = 0 (diag/out/f68f_1b_cut_partial.txt; the run was stopped by a restart of the
container):

| ε | θ_c | α | A_a = A_b | dA(1/2) | dA(1/4) | dA(1/8) |
|---|---|---|---|---|---|---|
| 9.95 | 0 | 3.203 | 6.0526 | 2.51e-2 | 6.27e-3 | 1.57e-3 |
| 10 | 0 | 4.075 | 4.7790 | 1.56e-2 | 3.89e-3 | 9.72e-4 |
| 10.5 | 0 | 8.880 | 2.2853 | 3.41e-3 | 8.53e-4 | 2.13e-4 |
| 12 | 0 | 15.93 | 1.3952 | 1.16e-3 | 2.90e-4 | 7.26e-5 |
| 16 | 0 | 25.68 | 0.9975 | 5.15e-4 | 1.29e-4 | 3.22e-5 |

The only cut locus at p = 0 is the tie θ = 0 of the pitchfork. Its weight on the 128-phase grid (one phase counted
twice, A_b/128) explains the 0.011 offset of the expansion table at ε = 12 (1.395/128 = 0.0109). The remaining
points (ε = 9.87 and 40 at p = 0; p = 1, 2, 4) are being computed with the converted script.

### Scripts and records: the diag form (16 Nov, evening)

**The diag convention** (diag/b2_dev.py, diag/c0_minvar.py), now followed by every f68f script through
diag/f68f_common.py:
- run from the repository root as a module, `python -m diag.<name>`, with argparse options;
- `--workers N` via core.stats.parallel.parallel_map (0 = CAUSALITY_WORKERS or all cores);
- every job cached atomically in diag/out/<name>_cache/, with the md5 of the code in the key;
- the output of a run is diag/out/<name>_<UTC time>.json and .txt, with the md5 of the code; nothing is overwritten;
- scripts that use the 80-bit instrument of F6.8e3 refuse to run without an 80-bit long double, as F6.8e3 does. These
  are f68f_layer and f68f_1a; they need Linux or WSL.

| Script | Computes | Checked after the conversion against the record of 16 Nov |
|---|---|---|
| diag/f68f_dev.py | step 0, checks T1–T5 | output identical to diag/out/f68f_step0.txt |
| diag/f68f_layer.py | step 0′: bulk and layer; options for the p and s quadrature checks (replace diag/f68f_pconv.py and diag/f68f_sconv.py) | one direction at n = 24 (bulk and layer) identical to the cached values of the run; end to end with `--p-list 0.5 --n-s 16 --bulk-only`: −6.0904258135e-06, equal to f68f_sconv.txt; the slope assembly of the full mode is the code of 16 Nov unchanged, not rerun end to end |
| diag/f68f_c2fd.py | the F6.8c2 code as a derivative in the age; `--mult` | runs end to end (the same call of optical_row as before); at k = 0.12, T = 50 it gives K̂_bb 0.8010366682 against 0.8010372305 in the F6.8c2 record of 29 Oct (diag/out/f68f_c2fd_20261116T214446Z.txt) |
| diag/f68f_1a.py | step 1a; `--levels` replaces diag/f68f_1a_levels.py | three volumes identical to the cached values (to 6e-17); the full run is being compared with f68f_1a.txt |
| diag/f68f_1b.py | step 1b, plane wave | ε = 4, 12 at p = 2 and the onset at p = 1 identical to f68f_1b.txt |
| diag/f68f_1b_q.py | step 1b at finite q, expansion (the failed exact check removed) | ε = 4, 12 at p = 2 identical to f68f_1b_q.txt |
| diag/f68f_1b_exact.py | step 1b at finite q, exact envelopes | ε = 9.5, p = 0 on 4 phases identical to the per-phase record; ran on the author's machine (Windows, 8 processes) |
| diag/f68f_1b_cut.py | cut loci | ε = 12, p = 0 being compared with f68f_1b_cut_partial.txt |

**Removed:**
- diag/f68f_pconv.py, diag/f68f_sconv.py: now options of f68f_layer;
- diag/f68f_1a_levels.py: now `f68f_1a --levels`;
- diag/f68f_1b_kappa.py: an untested draft.

The records made before the diag form keep their names and are cited above as they are.

**A reproducibility note.** Today's F6.8c code on Linux (numpy 2.5.3) reproduces the F6.8c2 record of 29 Oct (Windows,
numpy 2.2.6) at k = 0.12, T = 50 to 7e-7 relative in K̂_bb. The cause is not identified. Candidates are the numpy
version or platform, or a change of the code between 29 Oct and 10 Nov, which the git history (from 10 Nov) does not
cover. It is far below the 2.61 / 1.49 difference above and does not explain it.

## Step 1c: the age term of the expected BD action at finite ε (16 Nov, night; diag/f68f_1c.py)

The author's goal of 09:59 is E_opt(k,T,b) = ∂E[S_BD]/∂b and K_eff = (E_opt(b) − E_opt(0))/b at finite ε = b(kT)², and
whether the static response converges with T. Steps 0–1b produced the parts; this step assembles them.

### The age derivative at finite b, exactly
With the age cut at coordinate time T, dE[S]/dT is the integral over the events y of the past of x that lie at time
distance exactly T. Those are the y with first arrival T_opt(x, y) ≤ T. With weight √−g(y) = n(y)³ and the kernel of
the interval volume V(x, y):

  dE/dT ∝ P_x ∫_{T_opt ≤ T} d³y n(y)³ g(V(x, y)),   at Δt = T.

This holds at any b (no expansion). In flat space it is the layer formula of step 0′.

### The scaling limit (κ = kT → ∞, b → 0 at fixed ε and p)
- **Only needles count.** The kernel decays with V, so only pairs near the cone contribute: σ² = T² − T_opt² ≪ T²,
  r ≈ T.
- **Constant factors drop out.** n³, P_x and the measure tend to constants, which the moment μ₀ = 0 removes. Their
  cross terms with the volume are O(b) = O(ε/κ²).
- **Per direction p,** with s = l₂σ² (the variable of F6.8c) and q² = (kσ)² = k²s/l₂:

  ρ(p) = (T/(2l₂)) ∫₀^∞ ds ⟨g₀(s² A(ε, p, θ; q))⟩_θ,

  with A the needle volume of step 1b (relative to flat at the same σ). For ε < π² there is one ray at every phase,
  and A is analytic in q²: A = A₀ + A₁q² + A₂q⁴ + A₃q⁶ + …, A₀ = A_PW.

### The BD moments select the q⁶ order
With u = s (A(s)/A₀)^{1/2}, ∫ds g₀(s²A) = Σ_m d_m A₀^{−(m+1)/2} μ_m, where d_m are the coefficients of ds/du. The
moments of g₀ are (computed):
- μ₀ = μ₁ = μ₂ = 0;
- μ₃ = 1/2, μ₄ = 2.6587, μ₅ = 10.

So every term below s³ vanishes, and the first survivor is

  ρ(p) − ρ_flat(p) = (T/(2l₂)) 4 μ₃ (k²/l₂)³ ⟨ŝ₄ A₀^{−2}⟩_θ + O(k⁸),
  ŝ₄ = −e₄ + 5e₂e₃ − 5e₂³,   e₂ = a₁/2, e₃ = a₂/2 − a₁²/8, e₄ = a₃/2 − a₁a₂/4 + a₁³/16,   a_j = A_j/A₀.

This is the "D survivor" of step 0′ (the τ⁶ part), now at finite ε.

### What it means for the age of the field
- With dΩ = 2π dp/κ and T = κ/k, the age density is dE/dn ∝ k⁵ Φ(ε), with ε = b n², and

  Φ(ε) = ∫ dp ⟨ŝ₄ A₀^{−2}⟩_θ.

- At small ε, Φ = C₂ε², which gives E ∝ b² n⁵. This is the (kT)⁵ law of F6.8e3, and the check below compares it with
  F6.8e3 directly.
- At finite b:

  E(n) ∝ ∫ dn Φ(b n²) = (1/(2√b)) ∫₀^{bn²} ε^{−1/2} Φ(ε) dε.

- **The static limit of an old field exists if and only if ∫^∞ ε^{−1/2} Φ(ε) dε converges,** that is, if Φ decays
  faster than ε^{−1/2}. The author's readings A, B and C become statements about the large-ε behaviour of one function
  of one variable.

### The check against F6.8e3 (the author's check 1, b → 0)
At large n, the F6.8e3 bulk density per direction must be

  bulk(p) = 4 μ₃ k⁵ n⁵ C₂(p) / l₂⁴,   C₂(p) = lim ⟨ŝ₄A₀^{−2}⟩_θ / ε²   (corrected; see below: first written with 2),

with no free constant; corrections are O(1/n²) and O(k²). This compares two independent codes in their common range:
- the F6.8e3 second-order instrument (80-bit, finite n);
- the series of this step (scaling limit, all orders in ε).

### Range of this step, and what it cannot decide
- The q² series needs one ray per phase, so it holds for ε < π² (at p = 0; up to 18.7 at p = 4).
- The author's grid (ε ≤ 4) measures the departure from the second order, N(ε) = Φ/(C₂ε²).
- The asymptotic question (Φ at ε → ∞) needs the multi-ray range. There the cut loci enter (step 1b), each adding
  terms polynomial in c = q²/2 per window, so the same moment selection applies. An extension of the instrument is
  needed and is not done here.

### The other checks of 09:59
- **Check 2, zero age.** At n → 0, ε → 0 and Φ → 0.
- **Check 3, k = 0.** A = 1 at every q, so Φ = 0 identically.
- **Check 4, both signs of b.** A(−ε, θ) = A(ε, θ + π), so the θ-mean and Φ are even in ε. This is exact, not a
  computation.
- **Check 5, conformal and mixed sectors.** They are not touched: only the optical sector is computed.
- **Check 6, stable instrument.** The series is checked against the exact needle (envelopes) and against the O(q²)
  expansion. The linear term of Φ, which must vanish, is reported.

### Step 1c: the series instrument, first checks (diag/out/f68f_1c_check_20261116T220043Z.txt)
- **A₀ equals the plane-wave A of step 1b to 1e-12.**
- **A₁ equals the O(q²) coefficient of diag/f68f_1b_q to ~3e-8 relative.** The two codes use different RK steps.
- **The linear term of Φ, which must vanish, is ≤ 3e-10.**
- **Against the exact needle at ε = 4, p = 2, θ = 2:** −1.6e-8, −6.0e-8 and −3.2e-7 at q = 1/2, 1/4, 1/8 with the
  exact needle as it was. The growth towards small q is the exact needle's own error (below).
- **At the other three check phases the exact needle failed** ("window": the ray families did not cover the window
  after it was widened).

### The exact needle: two faults found by the series, and fixed (diag/f68f_1b_exact.py, 16 Nov, night)
1. **Coverage.** The ψ window had a constant pad (0.05). Near s = 0 and 1 this pad exceeded the reach of the ray
   families, and the window search failed. Now the window is ∝ √(s(1 − s)) (the paraxial shape). The slope range of
   the families is chosen so that they cover the window widened twice at every node. A window outside the families is
   now an explicit error ("coverage"), not a silent widening.
2. **Interpolation.** S(ψ) along a family was interpolated by cubic Hermite (S and dS/dψ), with errors ~dψ⁴. With
   d²S/dψ² = (dP/dv)/(dψ/dv) the interpolation is now quintic Hermite (error ~dψ⁶). Against the series at three
   phases (ε = 4 and 2):
   - |exact − series| ≤ 5e-9 (q = 1/2), ≤ 2e-8 (q = 1/4), ≤ 1.2e-7 (q = 1/8);
   - before the fix: up to 3e-7 at q = 1/8, or the window failure.
- **Consequence for the records made with the old version.** These are the p = 0 band in
  diag/out/f68f_1b_exact.txt/.npy and the author's run diag/out/f68f_1b_exact_20261116T212730Z.*. Against the fixed
  version at ε = 9.5, p = 0 (two phases) they differ by 1e-7 to 4e-7 (q = 1/2) and 2e-6 to 4e-6 (q = 1/8). This is far
  below δ_A (1e-3), so the readings of those tables stand. The minΦ = −3.4e-7 at ε = 16, p = 4 is of the same origin
  (the cubic interpolation) and is not repeated with the fixed version at the phases checked. That row is still not
  used until rerun.

### Step 1c: the b → 0 check at n = 1536 (diag/out/f68f_1c_e3_20261116T233417Z.txt; F6.8e3 bulk diag/out/f68f_layer_20261116T233417Z.json)
The question left open at n = 384 (ratio 1.0029): is the residual c/n² or O(k²)?
- **n = 1536, k = 0.12:** p-integrated F6.8e3/series = 0.999449 (−5.5e-4). Pure c/n² would give +1.8e-4.
- **Per direction, extrapolated as c/n² from 384 and 1536:** the total is 0.99922. The ratio runs smoothly from 1.0007
  (p = 0.16) to 0.9971 (p = 12.4). The directions p > 12 carry little weight, and their F6.8e3 values have lost digits
  at n = 1536 (values quantised at ~1e-3 absolute).
- **Reading.** The residual has two parts:
  - c/n²: +3.7e-3 at n = 384;
  - a remainder of about −8e-4 that does not vanish with n. This is O(k²) or the F6.8e3 quadrature; it is not
    separated (that needs another k; k = 0.06 has lost its digits, k = 0.085 is not run).
- **Consequence.** The agreement of two independent codes in the b → 0 limit is at the 1e-3 level, with no free
  constant. The K.18 draft now uses n = 1536 (the 1/n² part is 2e-4 there).

## Step 1d: the needle volume at finite κ and b in the real medium (17 Nov, night; diag/f68f_1d.py)

The author's 1b protocol asks for agreement with brute force at moderate ε (2, 4), at the levels κ = 10√ε and 20√ε.
The chord-centred brute force of step 1a cannot do it: the ray bends by ~ε/κ, far beyond the needle width q/(2κ).

### The instrument
- **The medium.** n = exp(b cos(k·w + θ)) is stratified along k. The optical distance between two events depends only
  on their separation along k (Δz), the horizontal distance ρ and the start phase.
- **The families.** For the families from x (phase θ) and from y (phase κ cos α + θ), rays are integrated as graphs
  z(ρ) with the exact Lagrangian n(z)√(1 + z′²) (RK4; the excess τ = T − ρ carried separately), with the slope
  derivatives.
- **The tables.** At each ρ node the first arrival τ(Δz) is the lower envelope over the family (quintic Hermite in
  z), tabulated on a window around the ray. It is interpolated by 7-point Lagrange in ρ and cubic Hermite in Δz.
- **The volume.** V = ∫dl dX dY n³ h₊, with h = Δt − T_x − T_y and Δt = √(T_xy² + σ²) (fixed σ), normalised by
  πσ⁴/24. Per (l, X) the root Y* of h = 0, Gauss in Y on [0, Y*], a fine grid in X, Gauss in l.
- **No expansion** in b, q or 1/κ.

### A fault found in the flat test, and fixed
- **The symptom.** In flat space (ε = 1e-9, κ = 20, q = 1/4) the instrument gave A = 0.99939. The same code with exact
  distances in place of the tables gave 1.0000032.
- **Not the tables.** Recorded over all 4297 calls, the tables agree with exact distances to ≤ 1.2e-10 with no
  infinities.
- **The cause.** The Y root used plain regula falsi. At isolated (l, X) points it stalled with one end fixed at
  a = 0. The returned Y*² = (a + b)/2 was then half the root (Y* off by exactly 1/√2 at those points, on about half
  of the slices).
- **The fix.** The Illinois variant, plus an explicit "Y root not converged" error. Flat test after the fix: 1.0000031721
  (tables) against 1.0000031764 (exact distances), 4e-9.

### Against the series (development, decides nothing; diag/out/f68f_1d_20261116T233144Z, …233451Z, …233644Z)
| ε | p | θ | q | A(10√ε) − series | A(20√ε) − series | ratio | extrapolated − series |
|---|---|---|---|---|---|---|---|
| 2 | 1 | 1 | 1/4 | −5.66e-4 | −1.41e-4 | 4.00 | +3.7e-8 |
| 4 | 2 | 2 | 1/4 | −1.19e-3 | −2.97e-4 | 4.00 | −2.7e-7 |
| 2 | 0 | 0.5 | 1/8 | −1.72e-4 | −4.23e-5 | 4.07 | +1.0e-6 |
| 4 | 0 | 3 | 1/8 | −1.48e-5 | −2.96e-6 | 5.01 | +1.0e-6 |
| 4 | 4 | 1 | 1/8 | −3.74e-3 | −9.38e-4 | 3.99 | −3.0e-6 |
| 4 | 1 | 5 | 1/8 | −4.82e-4 | −1.21e-4 | 3.98 | −9.1e-7 |
| 2 | 0 | 0.5 | 1/2 | −2.17e-4 | −5.44e-5 | 4.00 | −6.2e-8 |
| 4 | 0 | 3 | 1/2 | +3.02e-5 | +8.07e-6 | 3.75 | +6.8e-7 |
| 4 | 4 | 1 | 1/2 | −3.71e-3 | −9.30e-4 | 3.99 | −3.7e-6 |
| 4 | 1 | 5 | 1/2 | −5.10e-4 | −1.28e-4 | 3.97 | −1.1e-6 |

- **The finite-κ corrections are O(1/κ²).** Dividing b by 4 divides them by 3.97–4.07 wherever they are above 1e-4. Where they are
  ~1e-5 the instrument's own ~1e-6 error shows in the ratio (5.01, 3.75).
- **After the 1/κ² extrapolation, the brute force agrees with the series to ≤ 3.7e-6** at all ten points. This checks
  the whole chain: rays, Jacobi fields, q² series and the scaling limit, against the real 3+1 medium.
- **For the author's 1e-3.** The raw difference at κ = 10√ε exceeds it at ε = 4, p = 4 (3.7e-3); at 20√ε it is
  9.4e-4. These are physical O(b) corrections, not instrument errors. The K.18 draft compares the extrapolation; the
  author decides.
- **Cost:** 16–40 s per volume.

## Step 2: the assembly of K_eff (17 Nov, night; the author's order of 22:49, point 2; diag/f68f_keff.py)

### From Φ to K_eff
The goal of 09:59: E_opt(k, T, b) = ∂E[S_BD]/∂b and K_eff = (E_opt(b) − E_opt(0))/b. In the scaling limit the age
density per unit age n = kT is dE/dn = c Φ(bn²) (step 1c); c C₂ is the b² coefficient of F6.8e3. At fixed k and n:
- E(b, n) = c ∫₀ⁿ Φ(bm²) dm;
- E_opt = c ∫₀ⁿ m² Φ′(bm²) dm = (c/2b^{3/2}) M(ε), with M(ε) = ∫₀^ε √u Φ′(u) du;
- E_opt(0) = 0, because Φ is even in ε (check 4);
- K_eff = (c/2b^{5/2}) M(ε);
- K₂ = lim_{b→0} K_eff = (2/5) c C₂ n⁵: the (kT)⁵ law of F6.8e3.

So **R(ε) = K_eff/K₂ = M(ε)/((4/5)C₂ε^{5/2})**:
- All K_eff curves of (b, k, T) collapse onto R. "Several (b, k, T) per ε" of 09:59 is automatic in the scaling limit;
  the departures from it are the O(b) corrections of step 1d.
- R = 1 is the null model (the second-order extrapolation from F6.8e3).
- At fixed b: K_eff(T) = K₂(T) R(bk²T²).
- Local exponent: γ_T = d ln|K_eff|/d ln T = 2 d ln|M|/d ln ε. It is 5 for the second order, 0 for a static limit.
- Convergence of K_eff with T needs ∫^∞ √u Φ′ du to converge, i.e. Φ to approach its limit faster than u^{−1/2}: the
  criterion of step 1c, now for K_eff.

### Computation
- **By parts, without Φ′:** M(ε) = √ε Φ(ε) − ½∫₀^ε u^{−1/2}Φ(u) du = √ε Φ(ε) − ∫₀^{√ε} Φ(t²) dt.
- **The t integral** by Gauss–Legendre on panels in t = √u. Φ(t²) is smooth and even in t; in u the integrand would be
  ~u^{3/2} at 0.
- **Φ(ε)** = Σ_p w_p ⟨ŝ₄A₀^{−2}⟩_θ from diag.f68f_1c (its per-job cache is reused).
- **Selftest:** Φ = u² gives R = 1, and Φ = u³ gives M = (6/7)ε^{7/2}, both to 2e-16.

### The proposed card (written before any K_eff row was seen)
The F6.8f card draft (docs/F6_8f_card_draft.md) was rewritten around K_eff at 23:40 UTC, while the development run
below had printed no values:
- judged ε = 0.25, 1, 2, 4 (panel ends);
- readings A / B / C as PASSED / FAILED / UNDECIDED, by the sign of M and γ_T(2 → 4) against 2.5 and 5.

### Step 2: development rows (diag/out/f68f_keff_20261117T001242Z.txt; decides nothing)
Resolution reduced against the card: 4 directions (p panels 0:12:3, 12:24:1), 16 phases. Panels end at ε = 0.5 and 3
(not judged); G = 4.

| ε | Φ(ε) | N = Φ/(C₂ε²) | R = K_eff/K₂ |
|---|---|---|---|
| 0.5 | −1.610e-10 | 1.0125 | 1.0139 |
| 3 | −9.231e-9 | 1.6126 | 1.6888 |

- **γ_T(0.5 → 3) = 5.57.** The second order gives 5. On this range K_eff grows faster than the second-order
  extrapolation and keeps the sign of C₂. There is no bend towards saturation up to ε = 3 at this resolution.
- **Checks.**
  - 1 (b → 0): N(0.2) = 1.0020; R(0.5) = 1.014.
  - 2 (zero age): N at ε = 2.4e-3 is 1.000001.
  - 3 (k = 0): Φ = 0 exactly.
  - 4: exact.
  - 5: no code path.
  - 6 (stable instrument):
    - one ray per phase, stationarity 2.7e-10;
    - R with G = 3: differences 4.2e-5 (ε 0.5), 4.5e-4 (ε 3), below the proposed 1e-3;
    - Φ with 32 against 16 phases: 3e-12 (ε 0.5) and 2.3e-3 (ε 3). 16 phases are not enough at ε = 3. The card uses
      128 against 64.
- **This decides nothing.** The criteria of the card (written before these rows) are not changed by them. The judged
  ε (0.25, 1, 2, 4 at 16 directions, 128 phases) are not computed.
- **What it suggests for the physics question, as an indication only.** Up to ε ≈ 3 the finite-amplitude optical age
  term strengthens the second-order growth instead of saturating it. Whether Φ turns over beyond the first caustic
  (ε > π²) is outside this instrument (F6.8g).
