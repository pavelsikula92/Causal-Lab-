# The spatial part of the metric and momentum: from one principle of how the order responds to matter

A derivation on paper, 23 Oct 2026 (the author's assignment: "we must create the spatial part and the momentum").
Check: `focusing_gamma.py` (sympy). None of it is a result of a test.

## 1. Why they cannot be "added"

Neither the spatial part of the metric nor the momentum of the field are further fields that could be added to the
model. In a causal set the whole geometry is order and number: the cones (order) fix the metric up to a scale, the
number of elements fixes the volume. The spatial part of the metric is therefore **the shape of the cones**, i.e. how
the order deforms around matter. On a fixed sprinkling (F1–F3, F7.4) the order does not deform, which is why the
spatial part is missing there and can only be put in. One principle is missing: **how the order responds to matter**.

## 2. Principle: focusing of light proportional to the energy flux of matter, in every direction

A candidate in the language of Causal Physics (Jacobson 1995 in counting form; the area instrument K.14):

> The number of molecules of every local horizon (the area of a cut of the light cone) decreases in proportion to the
> energy flux of matter that passes through the horizon. The flux in a light-like direction k is T_kk = Σ F (u·k)²
> (matter with a carried 4-velocity u, F3.6).

In the continuum limit this is the Raychaudhuri condition R_kk = 8πG T_kk for every light-like direction k.

## 3. What follows from the principle (weak static field)

The metric −(1 + 2φ)dt² + (1 − 2ψ)dx² with **independent** unknowns φ (clocks) and ψ (rulers, width of the cones).
The linearised Ricci curvature gives for k = (1, n):

    R_kk = ∇²φ + ∇²ψ + nⁱnʲ ∂ᵢ∂ⱼ(ψ − φ)        (sympy, focusing_gamma.py: difference exactly 0)

For matter at rest T_kk = ρ is the same in all directions. The direction-dependent term must therefore vanish: the
traceless Hessian of ψ − φ is zero, and with ψ, φ → 0 far from the source **ψ = φ, so γ = 1**. The rest gives
2∇²φ = 8πGρ, i.e. **Newton's law**. Neither a scalar field equation nor the ansatz H1 is needed for this: both follow.
For moving matter T_kk = ρ(u·k)² is direction dependent and the principle gives the components h_0i
(gravitomagnetism): `jacobson_linear_check.py` shows that it is satisfied exactly by the tensor field Σ F u u, not by a
scalar.

Comparison: the scalar rule of the F3 series (clocks only) has ψ = −φ (geometry Π²η, γ = −1), H1 without the principle
has ψ = φ as an assumption. The focusing principle **enforces** ψ = φ.

## 4. Momentum

- **Matter.** The principle is consistent only if matter conserves its 4-momentum (Bianchi identity: ∇_a T^ab = 0 is a
  condition of integrability of R_ab ∝ T_ab). The present model conserves the scalar flow, not the 4-momentum: when
  directions merge it discards the excess |P| − F (F3.10 unit: 10⁻⁵ to 5·10⁻⁴). The derived correction of the rule:
  **at merging the flow equals |P|** (the rest mass grows in an inelastic merging, as in special relativity), the
  direction P/|P|. Then the 4-momentum of matter is conserved exactly.
- **Field.** In general relativity the momentum of the field is not a local quantity (a pseudotensor). It shows as the
  momentum carried away by waves and as a balance between bodies through the deformation of the geometry. In Causal
  Physics it would be carried by **a deformation of the order that propagates**. Matter then does not bend by a "kick":
  it moves freely (geodesically) along the deformed order, and the change of its momentum is an exchange with the
  order. As long as the order is fixed, the momentum of the field does not exist and can only be put in.

## 5. What follows for the laboratory

1. Matter with 4-momentum: the merge rule F = |P|. A small test: does the localisation of F3.6 survive?
2. An order that responds to matter according to the focusing principle (phase F6): the cones deform so that the number
   of molecules of horizons decreases in proportion to T_kk. Measured **independently** from the order: φ from the
   clocks (lengths of chains), ψ from the width of the cones or from spatial distances. γ = ψ/φ is then measured, not
   put in.
3. The momentum of the field and waves: two bodies and a moving source on a dynamic order (balance with retardation,
   polarisation).

Limits: the principle is a postulate (Jacobson's assumption of local equilibrium and of entropy proportional to area),
not a derivation from an even deeper rule. It is, however, the only postulate about the order from which γ = 1, Newton,
the tensor structure and waves follow. The Unruh temperature remains an input in it.
