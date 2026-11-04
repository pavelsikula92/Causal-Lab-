# ln Π as a potential: a motivation from principle for the choice of f (an analytic note)

Written on 22 Oct 2026 at the author's request. **Status: a motivation from principle, not a derivation.** The note
explains why f = ln Π is the primary candidate for the argument of the operator in the equation B_g f(Π) = ρ_s
(decision 1). It also delimits what this argument proves and what not. The choice of f is not derived retroactively
from the result of a numerical test (R132).

## 1. Why the discrete operator does not decide about f

The Benincasa–Dowker operator is, for a fixed geometry, linear in the field on which it acts. For f with f(1) = 0,
f(Π) = f′(1) εu + ½ f″(1) ε²u² + O(ε³), where Π = 1 + εu. Hence:

- **Linear order.** E_g[B_g f(Π)] = ε f′(1) E_0[B_0 u] + O(ε²). The geometry (H1 or another) does not show at this
  order, and all functions with f′(1) = 1 give the same.
- **Second order.** E_g[B_g ln Π] − E_g[B_g (Π − 1)] = −½ ε² E_0[B_0 u²] + O(ε³). The difference is therefore the flat
  operator acting on u². Only in the continuum is it equal to −½ ε² ∇²(u²); on a discrete grid it is a statement about
  convergence. The H1 geometry enters only at order ε³.

A test of the response of the operator (the planned K.13b) therefore does not decide between the candidates. The
difference between them must come from a principle, and is then tested in the whole non-linear equation (F1.7 N). The
discrete test K.13a checks only the operator itself in the H1 geometry on fixed test fields.

## 2. Composition of clock ratios

Two quantities of the laboratory compose differently and must not be confused:

- **Along a chain the proper time adds:** W = Σᵢ Πᵢ. It is an additive sum of local contributions (dτ = Π dt along a
  static world line).
- **Between places the ratios of clock rates multiply.** Let r_AB = Π_A / Π_B be the ratio of the clock rates at the
  static places A and B. When transferred through an intermediate place C, r_AB = r_AC · r_CB. This is an exact
  identity, independent of the path of the signal.

The logarithm turns multiplication into addition: ln r_AB = ln r_AC + ln r_CB. The quantity ln Π is therefore additive
with respect to the composition of ratios; Π − 1 is not:

(xy − 1) − [(x − 1) + (y − 1)] = (x − 1)(y − 1).

The deviation is of second order, i.e. exactly at the order where the candidates begin to differ (section 1).

**Uniqueness under an assumption.** Require that the field variable is a function of the clock ratio relative to a
reference observer (Π_∞ = 1) and is additive under composition: f(xy) = f(x) + f(y) for x, y > 0. For a continuous f
the only solution is f = a ln x (Cauchy's functional equation). The weak-field normalisation f′(1) = 1 gives a = 1. The
candidate ln Π thus **follows from the requirement of additivity**. This requirement itself, however, is a new
assumption, not a theorem of the theory.

## 3. The connection with general relativity

In a static spacetime with a Killing vector ∂_t and lapse N = √(−g₀₀) two things hold:

- The gravitational redshift between static observers is ν_B/ν_A = N_A/N_B, and over several levels it multiplies.
- The four-acceleration of a static observer is **exactly** a_i = ∂_i ln N (covariant components, c = 1). It is an
  identity for every static metric, not a weak-field approximation.

The quantity c² ln N is therefore the exact "acceleration potential" of static observers. In a weak field it gives the
Newtonian potential φ = c² ln N ≈ c²(N − 1). In the H1 metric of the laboratory N = Π, so φ = c² ln Π. The size of the
acceleration measured by the spatial metric Π^(−2γ)dx² is |a| = Π^γ |∇ ln Π|.

The argument holds for static fields and static observers. For a non-static situation (motion of the source, the
dynamics of the order F6) it must be checked again what the reference observer is and whether the clock ratio stays
well defined.

## 4. What this argument does not prove

1. That the operator should act precisely on a potential. The equation B_g f(Π) = ρ_s can be written for any variable.
   Additivity is a natural property of a potential, but the requirement "the BD operator acts on an additive quantity"
   is a new axiom. It must be stated in the theory, not assumed.
2. The normalisation. A constant shift ln Π → ln Π + c corresponds to another choice of the reference clocks. Unlike
   □_g, however, the BD operator does not feel the shift only where R = 0, because B_g 1 = −R/2. The choice Π_∞ = 1 is
   therefore part of the equation, not an arbitrary convention. This point is open.
3. The non-linear behaviour. The difference between ln Π and Π − 1 is −(∇Π)²/Π² in local form
   (B_g f = f′□_gΠ + f″(∇Π)² − (R/2) f). Whether it leads to regular solutions and to the right second potential is
   decided only by the whole equation (F1.7 N, decision 1, point 3), not by this note.

## 5. Two independent routes that can support or refute it

Both routes are so far only formulated; nothing of them is computed.

- **Action.** Is there an action of the causal set (the BD action or an action with a field Π) whose variation gives an
  equation for ln Π and not for another function? Check: the variable in which the kinetic term of the action is
  quadratic.
- **Conservation of the flow.** Can the local flow balance of the laboratory (Φ, Q_v, decisions 2 and 4) be written as
  a conservation law ∇·J = ρ_s with J proportional to ∇ ln Π (Gauss's law for the acceleration potential), or with
  another function?

If both routes lead to ln Π, it is support independent of the numerics. If they lead elsewhere, the primary candidate
changes from principle, not according to the result of a test.
