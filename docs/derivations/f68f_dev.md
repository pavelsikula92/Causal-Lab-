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
