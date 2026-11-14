# C0: the mean of the field rule on a Poisson sprinkling (14 Nov 2026)

Written before the diagnostic run of the late deviation of F3.16 (diag/c0_late.py). It derives what the
"second-order feedback of the noise" can and cannot be, before any form is fitted.

## Setting

The rule of C0, event by event in a sprinkling P of density ρ:

    φ(v) − δ J(φ(v)) = φ0(v) − c Σ_{u ∈ P, u < v} w(n_uv) J(φ(u)),

where:
- n_uv = |P ∩ I(u, v)|, the number of events in the open interval;
- the sum runs over the pairs with ρV ≤ λ and n ≤ 16;
- J = μφ (linear) or μ sin φ.

The domain is a pyramid |x|, |y|, |z| ≤ X − t.

## Three facts

1. **The domain is past-complete.** For u in the pyramid, J⁻(u) lies in the pyramid: |x| ≤ |x_u| + (t_u − t) ≤ X − t.
   Hence φ(u) depends only on P ∩ J⁻(u), and its distribution is that of an unbounded sprinkling. For data
   homogeneous in x, its mean is the same function of t everywhere in the domain, not only in the inner region.
2. **Disjointness.** The interval I(u, v) lies in the future of u, so I(u, v) ∩ J⁻(u) = ∅. Hence n_uv and φ(u) are
   functions of the Poisson process on disjoint sets, and are independent given u ∈ P.
3. **Mecke (Slivnyak).** E Σ_{u ∈ P} f(u, P) = ρ ∫ E f(u, P ∪ {u}) d⁴u. The event v itself does not enter J⁻(u) for
   u < v.

## Consequence

Taking the expectation of the rule at a sprinkled event v:

    E[φ(v)] − δ E[J(φ(v))] = φ0(v) − c ρ ∫ d⁴u K(σ_uv) E[J(φ(u))],   K(σ) = Σ_n w_n Pois(n; ρπσ²/24),

**exactly**, with no term of second order in the noise.

- **Linear J.** E[J(φ)] = μ E[φ], so E[φ] obeys the Volterra equation (A3) of diag/c0_apex.py exactly. The mean of
  the linear modes therefore has no second-order noise feedback at any order, any density and any kernel. A measured
  deviation from (A3) is either statistical or an error in evaluating (A3) (the numerical integration and the slab
  weighting).
- **Nonlinear J.** Everything beyond (A3) lies in the one-point distribution of φ(u): E[μ sin φ(u)]. B replaces it by
  μ e^(−s²/2) sin f (Gaussian, zero mean, variance s²).
  - The error of B at u is exactly D(u) = E[sin φ(u)] − e^(−s²/2) sin f(u). Its slab estimate is the direct check of
    the card: ⟨sin φ⟩ − e^(−s²/2) ⟨sin(φ − η)⟩, per unit μ.
  - The residual R = measurement − B is the response of the rule to −μ D, linearized about B:

        R(t) ≈ −μ ∫_0^t F_0(Δ) D(t − Δ) dΔ  +  (the change of the mean field by R itself, solved together)    (S)

  - The leading non-Gaussian terms are E[cos η] − e^(−s²/2) ≈ κ s⁴/24 and E[sin η] ≈ −γ s³/6 (excess kurtosis κ,
    skewness γ). With s ≤ 0.2 rad at A 0.5 they are below 10⁻⁴.

## Predictions for the diagnostic (fixed now, before any new sprinkling)

- **P1.** For lin0 and lin1, the mean of meas − (A3) over the last two slabs is 0 within 2 SE, at every density, size
  of domain and kernel.
- **P2.** For the nonlinear homogeneous rows, R = meas − B equals (S) evaluated with the measured direct check D, within
  2 SE. For A 0.5 and A 1 the measured D is ≤ 0.002, so R ≈ 0 there.
- **P3.** No dependence of these residuals on the number of events per slab beyond their SE: a larger X puts more
  events into the same slabs, and the residuals must shrink as SE, not stay at −0.02.

If P1 holds and the −0.02 of F3.16 does not reappear, it was a fluctuation of the last slabs (they hold 355–1 911
events and their noise is shared across the modes of one sprinkling). If P1 fails, the evaluation of (A3) or the
readout has an error. That would be an instrument error, not physics.

Run: `python -m diag.c0_late --workers 4` (diag/c0_late.py md5 a74a98110fdd99f5114ef6ffe5026afd, LF line ends). ρ 20 and 30, X 7 and 8.5, T 6,
16 fresh sprinklings per configuration (seeds 15000+ and 16000+). Kernels MV8+ct and MV2+ct, through diag/c0_card.task
unchanged.

**Change of the run, 14 Nov, 10:25, before any result was read.** One sprinkling at ρ 30, X 8.5 took 4 040 s (the
search for candidates is O(N²) in the events, and X 8.5 has 2.2× the events), so the full plan would take about 7 h.
P3 is answered instead by the number of sprinklings: `python -m diag.c0_late --xs 7 --n 32 --workers 4` (ρ 20 and 30,
X 7, 32 fresh sprinklings each, seeds 15000–15031). The 4 finished sprinklings at ρ 30, X 8.5 are read from the cache
with `--xs 8.5 --rhos 30 --n 4`. The predictions P1–P3 are unchanged. diag/c0_late.py md5 97aea6e61ad371dfe5aea1409ae3bc55.

