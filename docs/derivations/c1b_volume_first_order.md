# C1b: the interval volume in the theory metric, to first order

**14 Nov 2026. A derivation and a numerical check; no test, decides nothing.** Script: diag/c1b_volume.py; output:
diag/out/c1b_volume.txt.

## 1. Model

- **Metric**: g = −Π² dt² + Π⁻² dx² (3D space), with N = Π⁻² = 1 + 2 a g(x). The cones tilt: |dx/dt| = Π².
  The density of events per coordinate volume is √−g = Π⁻² = N.
- **Weak field**: with Φ = −a g this is the Newtonian metric −(1 + 2Φ) dt² + (1 − 2Φ) dx². Its Ricci tensor is
  first order in ∇²Φ, so it is vacuum to first order where ∇²g = 0.
- **Source**: a uniform ball of radius R = 1.5: g = (3 − r²/R²)/2 inside, R/r outside (exact first-order vacuum
  outside; only the Weyl, tidal curvature there).
- **Relation**: u ≺ v iff t_v − t_u ≥ τ(x_u, x_v), with τ the straight-ray travel time ∫ N dl. The true cones differ
  by the bending of the rays, which is second order in a. The relation respects t, but is transitive only to
  first order (to be quantified before any card).

## 2. The exact mean still holds

The derivation of docs/derivations/c0_mean_exactness.md uses only that the relation respects a time function: the
interval of (u, v) lies above t_u and everything φ(u) depends on lies below it. So on a Poisson sprinkling of
density ρN the count in the interval is Poisson with mean ρV(u, v) and independent of φ(u), and

    f(v) (1 − δμ) = A − c μ ∫ d⁴u ρ N(u) K(ρ V(u, v)) f(u),    u ≺ v.

The volume is an exact 3D integral (the time extent at each spatial point w is what the two delays leave):

    V(u, v) = ∫ d³w N(w) (Δt − τ(u, w) − τ(w, v))₊ .

In prolate coordinates around the pair (s = |w−u| + |w−v|, η, azimuth) the measure is (s² − l²η²)/8 ds dη dφ, and
the flat case gives π(Δt² − l²)²/24 (checked to 1e-6).

## 3. First order: does the volume see the vacuum curvature?

A local field theory weighs a near-null pair by the geodesic proper time τ_g (to first order
τ_g² = Δt² − l² − 2a(Δt² + l²)⟨g⟩ along the segment). The causal set weighs it by V. They agree where
V = (π/24) τ_g⁴. Define ε = V / ((π/24) τ_g⁴) − 1, to first order in a.

**Argument (withdrawn 15 Nov).** The first version argued that a linear scalar of the Riemann tensor over the diamond
can only contain Ricci(T, T) and R, so that the Weyl tensor drops out. That is wrong: the diamond is not a point, and a
scalar such as ∫ f(t, |x|) E_ij x^i x^j over it (E the electric part of the Weyl tensor in the rest frame) is not
excluded. No general argument is given here. The zeros below are a numerical result for the pairs tested; a possible
reason, not shown, is that g is harmonic outside the ball (mean-value property).

**Numerical check** (per unit a; the two terms cancel from several hundred):

| pair | l | τ | ε / a |
|---|---|---|---|
| constant g (flat) | 4 | 0.7 | 0.00000 |
| uniform field, along and across (flat) | 4 | 0.7 | 0.00000 |
| vacuum, tangential, impact 3 | 4 | 0.4 / 0.8 / 1.5 | 0.00000 |
| vacuum, tangential, impact 2.5 | 6 | 0.4 / 0.8 | 0.00000 |
| vacuum, radial 2.2 → 5.2 | 3 | 0.4 / 0.8 | 0.00000 |
| vacuum, oblique | 3.9 | 0.4 / 0.8 | 0.00000 |
| vacuum, grazing the ball, impact 1.6 | 6 | 0.15 | 0.00000 |
| through the centre of the ball | 6 | 0.4 / 0.8 | +2.59 / +2.56 |
| impact 0.75 | 6 | 0.4 / 0.8 | +2.33 / +2.30 |
| inside the ball, short | 1 | 0.4 / 0.8 | +0.09 / +0.11 |

## 4. Reading

- **Vacuum**: for the pairs tested, to first order, the interval volume is the flat formula in the geodesic proper
  time, also for long pairs. Whether this holds for every vacuum diamond is not shown.
- **Matter**: a pair whose diamond crosses the ball has ε of order 2.5 a, independent of τ: a finite, nonlocal
  weight, as in the conformally flat wells (F5.3, the plateau, the spherical well).
- For this source, to first order, the departure from the proper-time kernel comes from pairs whose diamond crosses
  the ball. That the nonlocality is tied to the Ricci curvature in general, and not to the Weyl curvature, is not
  shown. A field point in vacuum still departs from the local theory through the pairs that cross the ball.

**Limits**: first order in a; the straight-ray relation; the high-density reading of the weights (the finite-ρ mean
needs the solver); one source shape. Second order (Weyl², the bending) is not covered.
