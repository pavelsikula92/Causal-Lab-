# The horizon route to a₀: Deser–Levin, Jacobson and Milgrom

A revised check of the derivation, 24 Oct 2026 (the author's text; a note of the laboratory at the end).
Status: an analytic check on paper, not the result of a simulation or an observational test.

## Summary

The Deser–Levin temperature alone provides a natural cosmological acceleration scale, but it fixes neither the law of
gravity nor the value of the MOND parameter without a further assumption. Milgrom (1999) showed that the temperature
difference T(a) − T(0) has a functional form resembling the modified inertia of MOND. He did not, however, thereby
derive a modified Poisson equation or a change of the gravitational constant.

If it is further assumed that this temperature difference can be put into Jacobson's Clausius argument as an effective
temperature and the result interpreted as G_eff, a simple model with deep-MOND scaling follows, but with a₀ = 2cH_dS
and with an approximately constant additional radial acceleration cH_dS also at high acceleration. Numerically a₀ is
about 9–11 times larger than the usual empirical value. A universal additional acceleration of this type is, moreover,
strongly constrained by planetary dynamics, the exact bound depending on the model tested and on the analysis.

Conclusion: this particular gravitational interpretation does not work as a direct derivation of the observed a₀. This
does not refute Milgrom's original hypothesis of modified inertia nor all possible horizon mechanisms. K.14-dS makes
sense as a separate test of the discrete geometry of horizons, but not as a direct confirmation of this derivation of
a₀, until a measurable link between the horizon statistics and a dynamical equation is defined.

## 1. Jacobson's argument and the choice of the local horizon

In natural units (c = k_B = 1) Jacobson's argument has the schematic form
- energy flux across a local Rindler horizon: δQ ∝ κ ∫ λ T_kk dλ dA,
- change of area from the Raychaudhuri equation: δA ∝ −∫ λ R_kk dλ dA,
- area entropy: δS = η δA,
- the Clausius relation: δQ = T δS.

With the Unruh temperature 2πT/ħ = κ the dependence on the normalisation of the acceleration κ cancels. The result is
the Einstein equation with a relation between Newton's constant and the coefficient of the area entropy, in the
convention used G = 1/(4ħη).

An important limit: in Jacobson's construction κ is a parameter of the local Rindler observer/horizon, not
automatically the gravitational acceleration g created by a source. One therefore cannot simply identify κ = g and
regard a change of the temperature as a change of the gravitational constant. This adds a new physical assumption,
which does not follow from the original Clausius argument.

## 2. The Deser–Levin temperature

For a uniformly accelerated observer in exact de Sitter space: 2πT(a)/ħ = √(a² + H_dS²), in SI
T(g) = ħ/(2π k_B c) · √(g² + c² H_dS²), g the proper acceleration, H_dS the Hubble parameter of exact de Sitter space.
For g = 0 the de Sitter temperature of the background remains.

It is not correct simply to substitute an arbitrary current H(t) of a general FLRW universe into this exact result.
In ΛCDM in general H(t) ≠ c√(Λ/3); a universe with a positive constant Λ approaches the de Sitter relation only
asymptotically in the future. Using H(t) is therefore a further extrapolating hypothesis, not the exact Deser–Levin
result.

## 3. What Milgrom (1999) actually proposed

The difference ΔT = T(a) − T(0), in natural units 2πΔT = √(a² + H_dS²) − H_dS. With â₀ = 2H_dS one can write
a μ̂(a/â₀),

    μ̂(x) = √(1 + (2x)⁻²) − (2x)⁻¹,

for x ≪ 1 μ̂ ≈ x, for x ≫ 1 μ̂ ≈ 1 − 1/(2x). In SI â₀ = 2c H_dS = 2c²√(Λ/3). This is a motivation for modified
inertia, not a separate derivation of modified gravity. Milgrom pointed out explicitly that it is not clear why ΔT in
particular should measure inertia, and that for circular trajectories the interpolating function may be different. In
a more general cosmological space, moreover, the de Sitter argument alone does not decide whether a₀ should refer to H,
to the curvature or to the constant Λ.

## 4. Two gravitational closures (an additional model ansatz: the temperature factor as G_eff/G)

**4.1 Full temperature.** G_eff/G = g/√(g² + c²H_dS²); with g = (G_eff/G) g_N it follows that g_N = √(g² + c²H_dS²),
g² = g_N² − c²H_dS². For g_N < cH_dS no real solution exists and the coupling falls to zero as g → 0. This closure does
not give the deep MOND regime.

**4.2 Temperature difference.** G_eff/G = g/[√(g² + c²H_dS²) − cH_dS], so g_N = cH_dS [√(1 + (g/(cH_dS))²) − 1].
- g ≪ cH_dS: g² ≈ a₀ g_N, a₀ = 2cH_dS;
- g ≫ cH_dS: g ≈ g_N + cH_dS + O((cH_dS)²/g_N): an approximately constant additional radial acceleration also at high
  acceleration. This is a separate problem of the model, not only a mismatch of the value of a₀.

Interpolating function: with x_H = g/(cH_dS), μ_H = [√(1 + x_H²) − 1]/x_H ≈ 1 − 1/x_H; in the standard variable
x = g/a₀ = g/(2cH_dS), μ(x) = √(1 + (2x)⁻²) − (2x)⁻¹ ≈ 1 − 1/(2x). The variables differ by a factor 2.

## 5. A numerical check

H₀ = 70 km s⁻¹ Mpc⁻¹: H₀ ≈ 2.27·10⁻¹⁸ s⁻¹, cH₀ ≈ 6.80·10⁻¹⁰, 2cH₀ ≈ 1.36·10⁻⁹, cH₀/(2π) ≈ 1.08·10⁻¹⁰ m s⁻².
Λ = 1.1·10⁻⁵² m⁻²: c²√(Λ/3) ≈ 5.44·10⁻¹⁰, 2c²√(Λ/3) ≈ 1.09·10⁻⁹, c²√(Λ/3)/(2π) ≈ 8.66·10⁻¹¹ m s⁻².
Against a₀ ≈ 1.2·10⁻¹⁰ m s⁻², 2cH₀ is about 11.3× and 2c²√(Λ/3) about 9.1× larger. cH₀/(2π) is close to the empirical
scale, but the factor 1/(2π) does not arise from the derivation; adding it for the sake of agreement would be a
calibration, not a prediction.
- H1: a₀(t) = 2cH(t), evolving with the expansion; not derived directly for a general FLRW universe from the exact de
  Sitter result.
- H2: a₀ = 2c²√(Λ/3), constant; it arises in the exact de Sitter limit; for today's ΛCDM it must be justified why the
  local dynamics uses this scale and not H(t).

## 6. Planetary tests: what can and cannot be claimed

In variant 4.2 the high-acceleration asymptote is an approximately universal additional radial acceleration cH_dS:
~6.8·10⁻¹⁰ (H_dS ≈ H₀) or ~5.4·10⁻¹⁰ m s⁻² (H_dS = c√(Λ/3)). The bound is not one constant independent of the model:
- Iorio (2011) from published planetary residuals for a Rindler-type additional acceleration: of order 10⁻¹⁵ m s⁻²
  (Mars) and 10⁻¹⁴ m s⁻² (Saturn); he distinguishes direct bounds from residuals from a full reprocessing of the
  ephemerides.
- Carloni, Grumiller and Preis (2011): a bound "without caveats" from the radar delay |a| < 3·10⁻⁹ m s⁻², which by
  itself does not exclude cH₀; according to them the tightest constraints come from perihelion precessions, from which
  it follows that the Pioneer anomaly (~8.7·10⁻¹⁰ m s⁻²) cannot be caused by a universal Rindler acceleration.

The precise claim is narrower: a constant universal radial term of order cH₀ is in strong tension with analyses of
planetary residuals, but the strength of the conclusion depends on whether the dynamical prediction is equivalent to
the Rindler model tested. A final exclusion would require propagating the exact equation of the model into the
orbital dynamics. The objection concerns the gravitational closure 4.2, not automatically modified inertia (different
predictions for orbits and the external field).

## 7. Decision for Causal Physics and K.14-dS

Supported: (1) the de Sitter background gives a natural scale cH_dS = c²√(Λ/3); (2) ΔT ∝ a²/H_dS at small acceleration
can motivate a MOND-like scaling; (3) interpreting ΔT as an effective gravitational coupling gives a₀ = 2cH_dS, not the
observed a₀, and creates an unwanted high-acceleration term.
Not supported: (1) that Jacobson's argument itself requires T(a) − T(0); (2) that the Deser–Levin temperature itself
derives a modified Poisson equation; (3) that 1/(2π) or the observed a₀ follows from this route without a further
parameter; (4) that the failure of the gravitational interpretation refutes every horizon explanation of MOND or the
whole of Causal Physics.
K.14-dS: not as a direct test of a₀ from this derivation (the link from the horizon statistics to a dynamical equation
is missing); yes as a separate geometric test (in advance: the measured quantity, the null hypothesis of area scaling,
the dependence on the density, the limits a → 0 and H → 0); as a test of MOND only with a model from which μ(x) and the
time dependence of a₀ follow without calibration.

## References
1. T. Jacobson, PRL 75, 1260 (1995), arXiv:gr-qc/9504004.
2. S. Deser, O. Levin, Class. Quantum Grav. 14, L163 (1997), arXiv:gr-qc/9706018.
3. M. Milgrom, Phys. Lett. A 253, 273 (1999), arXiv:astro-ph/9805346.
4. L. Iorio, arXiv:1012.0226 (2011).
5. S. Carloni, D. Grumiller, F. Preis, Phys. Rev. D 83, 124024 (2011), arXiv:1103.0274.

## Note of the laboratory (24 Oct)
Checked against the abstracts: Iorio gives 1·10⁻¹⁵ m s⁻² for Mars and 1·10⁻¹⁴ m s⁻² for Saturn (the original wording
of the text "10⁻¹⁵ for Mars and Saturn" corrected). Carloni et al. give a radar bound of 3 nm/s² "without caveats" and
state explicitly that the perihelia give the tightest constraints and exclude a universal Rindler acceleration of the
size of the Pioneer anomaly; cH₀ is of the same order. The analytic steps 3–5 agree with a check in sympy.
