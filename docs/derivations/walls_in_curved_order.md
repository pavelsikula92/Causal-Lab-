# Walls in the 2+1 order of the theory: the passive weight is energy minus tension (12 Nov 2026)

Note for B2 (decision 16). Continuum derivation, checked numerically in `diag/b2_curved.py`; decides nothing.

## 1. The metric and a profile with a closed causal relation

The order of the theory (H1, throughput, decision 12) has the metric

    ds² = −Π² dt² + Π⁻² (dx² + dy²),   √(−g) = Π⁻¹ (2+1).

For B2 the profile Π² = (x + ℓ)/ℓ = H/ℓ (H = x + ℓ) is used: g_tt is linear in x (a uniform-like field), and the metric
is conformal to −dt² + h with the optical metric h = ℓ² (dH² + dy²)/H², the hyperbolic plane of radius ℓ. Hence

    u ≺ v  ⇔  t_v − t_u > d_h(u, v),   cosh(d_h/ℓ) = 1 + ((H_u − H_v)² + (y_u − y_v)²) / (2 H_u H_v),

and the domain of dependence of a coordinate rectangle has closed boundaries (horocycles x = const, geodesics y = const).
A test particle at rest at x = 0 moves as x(t) ≈ −t²/(4ℓ) (Γ^x_tt = Π²/(2ℓ)).

## 2. A wall uniform in y feels no force if its tension equals its energy

For a field independent of y, box_g φ = V′(φ) gives φ_tt = Π³ ∂_x(Π φ_x) − Π² V′(φ). The static energy per unit
coordinate length in y is

    E_K[φ] = ∫ dx √(−g) [½ g^xx φ_x² + V] = ∫ dx [½ Π φ_x² + Π⁻¹ V] = ∫ ds [½ φ_s² + V],   ds = Π⁻¹ dx,

independent of Π. So **the kink in proper distance, φ = 4 arctan e^s, is an exact static solution at every position**
(checked: the continuum centre stays within 1e-4 for t ≤ 8 at ℓ = 24, while a test particle falls 0.44).

The general statement: for an object uniform along y with energy E and tension σ along y (per proper length), the
Killing energy per unit coordinate y is Π · E_y(L), L = Π⁻¹ the proper length per dy, and dE_y/dL = σ. The force is
∝ −∂_x [Π E_y(Π⁻¹)] = −Π′ (E − σ) (per proper length). **The passive weight of a wall is E − σ**: the redshift of the
clocks (Π) and the stretching of the proper length along the wall (Π⁻¹) cancel for a wall whose tension equals its
energy (the kink: T_yy = −ε). Dust (σ = 0) falls as a particle.

The sine-Gordon breather wall of frequency ω has E = 16 √(1 − ω²) and σ = −⟨L⟩ = E − ωJ = 16 (√(1 − ω²) − ω arccos ω)
(J the action variable); ω = 0.7: σ/E = 0.22, so it falls at about 0.78 of the particle rate (the continuum readout of
the energy centroid: −0.30 against −0.34 for 0.78 × particle at t = 6.5; the internal oscillation adds to the readout).

## 3. What this means for B2

- B2 as "the wall falls as the continuum" is sharper than planned: the continuum says the kink does **not** fall and the
  breather does. A discrete order must reproduce both: no spurious force on the kink from the density gradient (the
  coordinate density varies 1.7× across the domain) and the correct fall of an object with little tension.
- The launched kink keeps its proper speed (E_K = 8γ conserved); its coordinate speed is Π² v (kinematics of the cones).
- In 3+1 the same argument gives E_K ∝ Π · Π⁻² for a wall in two directions: a wall is pushed towards larger Π (away
  from matter) in the order of the theory. Not tested; recorded for B4.
