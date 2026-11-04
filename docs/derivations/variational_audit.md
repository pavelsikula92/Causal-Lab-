# The variational principle of the substrate: a mathematical audit and open decisions

24 Oct 2026. Status: a mathematical audit and a formulation of open decisions. The computations of sections 2–4 are
complete and checked (integrating factor, Euler–Lagrange equation, second variation). Section 5 is a formulation of a
future test, not a result.

## 1. Purpose

To distinguish three mathematically different objects:

1. **The action S**: a functional that fixes the dynamics or the equations of motion.
2. **A stability functional**: a quantity by which the stability of a configuration against admissible perturbations
   is judged.
3. **A Lyapunov functional**: a quantity whose monotonic evolution along a defined dynamics proves stability or
   convergence.

They must not be confused merely because all of them can be written as functionals.

## 2. The functional of causal resistance on a fixed graph

    R = Σ_v (I_v / Π_v) Δτ_v,     u_v = ln Π_v   ⇒   R[u] = Σ_v I_v e^(−u_v) Δτ_v.

### 2.1 First variation

On a fixed graph, and if neither I_v nor Δτ_v depends on u_v:

    δR[η] = −Σ_v I_v e^(−u_v) η_v Δτ_v.

For I_v > 0, Δτ_v > 0 and perturbations with independent local components η_v the first variation is not identically
zero.

**Conclusion:** R in this form has no interior stationary point with respect to independent u_v. On a bounded
admissible region the optimum may lie on the boundary. Without bounds a finite minimum need not be reached. The
conclusion changes if the admissible perturbations are subject to constraints (a global constraint, a conservation
law). Such constraints must be defined explicitly.

### 2.2 Dependence on the choice of time

If Δτ_v = Π_v Δt_v, then

    R = Σ_v (I_v / Π_v) Π_v Δt_v = Σ_v I_v Δt_v

and the dependence on Π_v vanishes exactly. If Δτ_v does not change under a variation of Π_v, the dependence of R on
the field remains.

**Open decision A1:** is Δτ_v the proper time, the coordinate time or another operationally defined quantity? How does
it change under a perturbation of the field Π?

### 2.3 Operational definition of I_v

The first variation assumes that I_v does not depend on Π_v. If I_v = I_v(Π, Φ, H), the corresponding variational
terms must be included.

**Open decision A2:** how is I_v measured or computed from the local state of the substrate? The definition must be
independent of the desired result, in particular of whether the field shows a stable configuration. (Axiom audit: I
so far has only a verbal definition; operationally it is planned in F0.5.)

Until A1 and A2 are resolved, R is defined as a **functional of a path in a given field** (section 5.8 of the theory),
not as a fully specified functional of the field.

## 3. The existence of a variational formulation does not fix a physical law

The static equation

    A(L) ∇²L + B(L) |∇L|² = f(L, x),     A(L) ≠ 0,

after multiplication by an integrating factor μ(L) turns into the Euler–Lagrange equation of the functional
∫ [½ a(L) |∇L|² + V(L, x)] with a = μA, if

    μ′/μ = (2B − A′) / A

and the right-hand side μf can be integrated into a potential (∂V/∂L = μf). Derivation: the Euler–Lagrange equation
gives a∇²L + ½ a′ |∇L|² = ∂V/∂L, so μB = ½ (μA)′.

For example equation (a) in 2+1,

    ∇²Π = (κ/3) ρ Π²,

is, for a fixed ρ and suitable boundary conditions, the Euler–Lagrange equation of the functional

    S[Π] = ∫_Ω [ ½ |∇Π|² + (κ/9) ρ Π³ ] d^d x.

**Conclusion:** a variational formulation by itself derives neither the exponent nor the spatial dimension nor a
concrete field law. Every choice of the argument f(Π) in the equation B_g f(Π) = κρ_s is variational, including ln Π.
It therefore selects neither β nor ln Π. The physical content must come from independently defined quantities of the
substrate and constraints, not from an additionally chosen free function a(L).

## 4. The static second variation

For the functional of section 3 and a perturbation Π ↦ Π + εη:

    δ²S[η, η] = ∫_Ω [ |∇η|² + (2κ/3) ρ Π η² ] d^d x.

For κρΠ ≥ 0, suitable boundary conditions and a non-zero admissible perturbation the form is non-negative, under
corresponding conditions positive.

This shows a local variational stability of the static functional against perturbations of the field. It does **not**
show the dynamical stability of the causal graph, the stability of the growth process, or a preference for a concrete
effective dimension.

## 5. Phase C: structural stability of a causal set

It is proposed to study separately a discrete action of a causal set of the Benincasa–Dowker type, possibly its
extension by an independently defined coupling to the flow. Research question:

> Under what conditions, if any, does the chosen action prefer manifold-like causal structures over non-manifold
> configurations and lead to a robust effective dimension D_eff ≈ 3?

Turning this into a test requires:
- an exact definition of the action and of the sign convention;
- a fixed ensemble of admissible causal sets;
- a definition of the dynamics or of the statistical ensemble from which the configurations are drawn;
- an independent measurement of manifold-likeness and of the effective dimension;
- a check of the dependence on the size of the system, the parameters, the boundary conditions and the samples;
- a criterion fixed in advance that can refute the hypothesis.

The mere existence of an action, or the suppression of some configurations, does not prove that D_eff = 3 is a stable
or preferred result.

**Status of the literature claim:** the connection with the entropic dominance of Kleitman–Rothschild posets, with the
phase behaviour in lower dimensions and with the suppression of non-manifold configurations is recorded as an
**unverified search hypothesis**. Before use in a scientific text the original papers, the range of dimensions, the
definitions of the action and the exact conditions of the results must be checked.

## 6. Decisions for the next steps

1. Close the definition of Δτ_v and its transformation rule (A1).
2. Fix an operational definition of I_v and determine its dependences (A2).
3. Keep R as a functional of a path until its variation with respect to the field is derived separately.
4. Do not infer structural stability from a positive static second variation.
5. Move the main test of the hypothesis D_eff ≈ 3 to phase C, based on an explicit discrete action and a protocol
   defined in advance.

**Working conclusion:** first define the quantities and the admissible changes. Only then decide whether an action or
a stability criterion selects a physically relevant configuration.
