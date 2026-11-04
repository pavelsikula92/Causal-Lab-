# The substrate rule from the action: a derivation instead of an insertion

Written on 28 Oct 2026 at the author's request. It is an analytic analysis; the computations were checked in sympy.
It follows `tensor_coupling.md`, F7.4a and decision 13.

## 1. What is to be derived

The rule that F7.4a verified on a causal set has three parts:

- **(a) own frame:** every flow narrows the cone in its own rest frame;
- **(b) addition:** the capacity deficits of different flows add;
- **(c) ratio:** the narrowing of the cone and the saturation (adding events) are in the ratio that gives γ = 1 (the
  principle of throughput).

In F7.4a the rule was an input. Below it is shown that (a) and (b) follow from axiom A6 and from linearity without a
further assumption about the order. (c) is a single number, and the action of the order itself fixes it.

Conventions: static weak field ds² = −(1 + 2Φ) dt² + (1 − 2Ψ) dx², d spatial dimensions, G = c = 1. Newtonian
Φ = −U. γ = Ψ/Φ.

## 2. Step 1: the matter side (A6) gives the own frame and the addition

A6: the action of a flow is the proper time along its world line, S_m = −Σ m τ.

- On a causal set the proper time is read by counting: the length of the longest chain or the volume of a diamond
  (the clocks of K.16).
- For a small perturbation g = η + h, δτ = −½ h_ab u^a u^b τ, where u is the 4-velocity of the flow.
- The variation of the matter action with respect to h is therefore ½ Σ m u_a u_b δ(world line) = ½ T_ab.

It follows that:

- **(a)** The order "sees" a flow only through u_a u_b. The own frame of the flow is therefore the frame in which its
  proper time is measured, and it does not need to be assumed.
- **(b)** The action is a sum over the flows, so at linear order the contributions add.

For this part nothing needs to be known about the dynamics of the order. A6 and a linear response of the order
suffice.

## 3. Step 2: Lorentz invariance leaves a single number

A linear, local and Lorentz-invariant response h_ab = G ∗ T has two channels: T_ab and η_ab T. For static dust this
leads to a one-parameter family (`tensor_coupling.md`):

    h̄_ab = U [ (1 − γ) η_ab + 2(1 + γ) u_a u_b ]

The structure of the rule of F7.4a is thereby forced up to one number γ:

- γ = −1: a pure trace coupling (Nordström, the fixed order of F3);
- γ = 1: a pure coupling to T_ab.

## 4. Step 3: γ is fixed by the action of the order

A general static quadratic action of the order with a coupling to dust is

    L = c_ΦΦ |∇Φ|² + c_ΦΨ ∇Φ·∇Ψ + c_ΨΨ |∇Ψ|² − κ' ρ Φ

(coupling only through Φ, because dust at rest has only T_00). The variation with respect to Ψ gives
c_ΦΨ ΔΦ + 2 c_ΨΨ ΔΨ = 0, so

    γ = Ψ/Φ = −c_ΦΨ / (2 c_ΨΨ)

independently of c_ΦΦ. The variation with respect to Φ fixes only the amplitude, i.e. the relation κ ↔ G.

**The Einstein–Hilbert action** (sympy, √−g R to second order, d = 2, 3, 4):

    c_ΦΦ = 0,   c_ΦΨ = −2(d − 1),   c_ΨΨ = (d − 1)(d − 2)   ⇒   γ = 1/(d − 2)

In 3+1, γ = 1. In 2+1, c_ΨΨ = 0 and Φ has no source, so there is no Newtonian attraction there.

**The principle of throughput** requires γ = 2/(d − 1), i.e. c_ΨΨ / c_ΦΨ = −(d − 1)/4. This agrees with the
Einstein–Hilbert action only for d = 3.

**Physical reading:**

- **The clocks have no stiffness of their own (c_ΦΦ = 0).** The temporal potential Φ appears only as a coupling, and
  time is determined from space.
- **Throughput is a consequence, not a postulate.** Dust has no stress and the order has the Einstein–Hilbert
  stiffness, so the spatial trace-reversed part does not deform (h̄_ij = 0). In 3+1 this is the same as
  ρ_c·c = const.
- **Density alone as a degree of freedom gives Nordström.** If the order had only the density (the conformal factor)
  as a dynamical degree of freedom, γ = −1 would follow, i.e. exactly the fixed order of the laboratory.

## 5. A candidate action of the order: the Benincasa–Dowker(–Glaser) action

- **What is known.** The expectation value of the BD(G) action on sprinklings goes in the continuum limit to the
  Einstein–Hilbert action plus boundary terms. Benincasa and Dowker showed this for the BD action (2010). For small
  causal diamonds it is examined by arXiv:2007.13192: d = 2 is derived explicitly, d = 3, 4, 5 agree to first order in
  the curvature.
- **The laboratory already uses the same operator:** B_g of K.11 and K.13a.
- **The result of the chain:**

      A6 (proper time of a flow counted by the order) + BD action of the order  ⇒  the rule of F7.4a at linear order (γ = 1)

- **What follows:** the substrate rule is the linear response of an order with the BD stiffness to the action of the
  flows. Part (c), the throughput, is a consequence of the BD action having an Einstein–Hilbert limit.

## 6. What is new and what is not

- **Known (textbook):** in the continuum the chain A6 + Einstein–Hilbert ⇒ linearised general relativity. It is also
  known that the BD action has an Einstein–Hilbert limit.
- **New for the theory:** the rule of the laboratory and the principle of throughput are not separate postulates.
  They are the linear response of an order with the BD stiffness to the action A6. Postulates thus become a
  derivation, but in the continuum limit.
- **New and publishable would be:** to show that the discrete BD action on a causal set itself selects γ = 1, i.e.
  to measure c_ΦΦ, c_ΦΨ, c_ΨΨ directly from the order. The sharpest prediction is **c_ΦΦ = 0**: a static perturbation
  of the clocks alone has no quadratic action. It does not depend on the normalisation, and neither a fixed order nor
  a scalar action meets it.
- **The main obstacle is noise.**
  - The signal is of order a² ∫|∇Φ|² d³x · T ~ 0.1–1 for a = 0.1 in the box of the laboratory (ℓ = 1).
  - The fluctuations of the local BD action grow with the number of elements, an estimated √N · O(1–10) ~ 10²–10³
    for N ~ 5·10⁴.
  - Directly it is therefore several orders of magnitude out of reach. A smeared (non-local) BD action and correlated
    sampling (the same points for different Φ, Ψ) may help. It is an order-of-magnitude estimate and a feasibility
    test must confirm it.
- **Out of reach of the laboratory:** a real dynamics (a path integral with the BD action and the non-manifold phase)
  is an open problem of the theory of causal sets. The chain of section 5 assumes the manifold phase, as F6.4b does.

## 7. Proposals (nothing running)

1. **K.17, a feasibility test** (no longer a prerequisite for F6.8, see section 8): the variance of the smeared BD
   action on a 3+1 sprinkling as a function of N and ε, and the variance of the difference under correlated sampling
   of two nearby metrics. It decides whether measuring c_ΦΦ, c_ΦΨ, c_ΨΨ is possible at all.
2. **F6.8** (now deterministic, section 8): the quadratic form of the discrete BD action in the static sector (Φ, Ψ).
   Predictions: c_ΦΦ = 0 and γ = −c_ΦΨ/(2c_ΨΨ) = 1 in 3+1. Alternatives: a scalar action of the density γ = −1, the
   principle of throughput outside 3+1.
3. **Novelty search:** the coupling of flows through the length of chains (A6) and the response of the BD action to
   matter in a weak field. So far only Sverdlov–Bombelli 2009 (a scalar field, not particles) and the work on the
   limit of the BD(G) action have been found.

## 8. F6.8 deterministically: the expectation value of the action instead of sprinkling (28 Oct 2026)

The route through sprinkling (K.17) runs into noise. The expectation value of the BD action, however, can be computed
directly by quadrature, without noise, so F6.8 does not depend on K.17. The author required that the values must not
be prescribed: they must come out of the kernel and the geometry.

**Sectors.** The static metric −(1+2Φ)dt² + (1−2Ψ)dx² can be written as Ω² · [−dt² + n² dx²], with ω = ln Ω = Φ and
ν = n − 1 = −(Φ + Ψ). The conformal factor does not change the causal order, only the density; ν changes the cones.
The static Einstein–Hilbert form in these variables is 6ω² + 8ων + 2ν² (× |∇·|²), and matter couples only to ω. It
follows that γ = −1 + B/(2C) = 1, with B = 8 and C = 2.

- **F6.8a, conformal sector (coefficient 6):** exact quadrature, implemented.
- **F6.8b, mixed sector (coefficient 8):** first-order cones suffice (optical length of a straight segment, as in
  K.13a). It is an independent check.
- **F6.8c, purely optical sector (coefficient 2):** needs the bending of rays at second order. It is demanding, but it
  suffices if locality and covariance hold (see below).

**Why the conformal sector plus locality suffices.** The mean BD action is exactly diffeomorphism invariant, because
it is defined on the order. A local limit with two derivatives can therefore only be α∫√−g R + β∫√−g.

- β = 0 is checked by the null k = 0.
- α is fixed by F6.8a.
- Locality in the other sectors is checked by F6.8b. **Result of F6.8b2 (revision): PASSED**, c_mix = 7.983 (EH 8),
  the covariant ratio at leading order; beyond it Δ ≈ 0.46 k² and an age term 1.78× larger (in detail in
  `findings.md`).

**Result of F6.8a: PASSED** (full run 28 Oct: K0 = 1.00005, c_conf = 6.0003, c1 = 3.98, c0 = 0.91, c_T = 0.1035; in
detail in `findings.md`). The development numbers (they agree with the full run; ℓ_k = (ερ)^{-1/4} = 1; a plane wave
f = cos kz averaged over the phase; K̂ = ratio to the Einstein–Hilbert value):

- **Instruments:**
  - the analytic second derivative agrees with finite differences of the kernel operator to 3e−4;
  - the closed Fourier transform of the diamond agrees with the quadrature of the kernel to 7e−8;
  - the null k = 0 gives J0/J ~ 1e−5.
- **Local part:** K_loc = 0.990 / 0.968 / 0.899 / 0.703 / 0.470 / −0.015 for k = 0.025 / 0.05 / 0.1 / 0.2 / 0.3 / 0.5.
  - The fit gives **K0 = 1.001**: the Einstein–Hilbert conformal stiffness 6 comes out of the discrete action.
  - The correction for a finite ℓ_k is non-analytic: 1 − K_loc ≈ (kℓ)²(3.9 ln(1/kℓ) + 1.1).
  - For wavelengths below ~13 ℓ_k the stiffness changes sign.
- **The age term (new, non-local):** K̂(k, T) = K_loc(k) − c_T k³T, c_T = 0.10.
  - Linear in T to 0.5 % and proportional to k³ to 2 %.
  - It was predicted before it showed in the numbers. Long, nearly light-like intervals perpendicular to k stay in
    phase with the field along their whole length. The moment conditions of the BD kernel kill their contribution up
    to order τ⁴, but not τ⁶, so every unit length of history adds ~k⁵.
  - The asymptotic 2D integral gives c_T = 0.102 (k → 0), the full quadrature 0.099–0.100.
  - Physically: the static conformal stiffness of the smeared BD action weakens with the time T for which the field
    has existed. The relative size is 0.1 (kℓ_k)³ (T/ℓ_k) = 0.1 k³ ℓ_k² T. For a Planckian ℓ_k it is negligible (the
    Sun at 1 AU after 4.6·10⁹ years ~10⁻⁷⁹).
  - For a smeared action with a mesoscopic ℓ_k it would give constraints. The stiffness at a scale L stays
    gravitational (deviation below 10⁻²) only when ℓ_k ≲ √(0.1 L³ / T). For the field of the Earth (L = 6400 km,
    T = 4.5·10⁹ years) this is ℓ_k ≲ 1 mm. For a laboratory source (L = 1 mm, T ~ a day) this is ℓ_k ≲ 10⁻¹² m.
  - Open: what T is in a real situation (the lifetime of the field, the cosmological past). It is the static
    conformal sector in a region with a sharp past boundary. Novelty has not yet been checked in the literature.
