# K: Instruments and calibration

## Objective

Instruments that read geometry and physics from a causal order, and their calibration against known answers: dimension, proper time from chains, cones and density (the blind instrument), the Benincasa-Dowker operator, horizon molecules, the energy-momentum tensor of flows, and interval volumes. A calibration has no physical conclusion of its own; it makes the physics tests possible, and a physics test that depends on a failed calibration cannot conclude.

## Contents

One folder per test, `K_instruments/<ID>/`, with the records of its full runs: `<UTC time>_<part>_seed<seed>_<test hash 8>.json` (parameters, seed, hashes of the test file and of the core, versions, the printed lines and the data) and the same lines as `.txt`. Reduced runs are in `<ID>/quick/` and do not count. Nothing is overwritten.

The verdicts and their reading are in `docs/findings.md`; the list of tests with their latest states in `docs/test_list.md` (`python run.py report`); the records of one test: `python run.py results <ID>`.

| Test | Question | Level |
|---|---|---|
| K.2 | Does the causal cone work as an instrument of bending? | calibrated: 0.913 and 0.919 of the GR value |
| K.3 | Can Π be read back from the density of chains? | calibrated in the tested range, replicated |
| K.4 | Does the quadrature compute the expected value of the smeared BD operator correctly, including the local limit and convergence? | calibration |
| K.5 | Is the discrete BD operator unbiased? | verified: +0.004 ± 0.008 (clustered by sprinkling) |
| K.6 | Does routing conserve the flow and is the window settled? | verified: residual ~10⁻¹⁵, decay at a rate of ~1 |
| K.7 | Are the chain functionals mirror symmetric? | verified: exactly 0 |
| K.8 | Does the probe measure the throughput T_out/T_in? | verified: balance 0, sensitivity 9·10⁻⁶ |
| K.9 | How faithful are the topology null models? | DAG null 90 % non-covering edges, narrower tail of the flow, the strict null shuffles 18 % of edges – influence of covering relations undecided |
| K.11 | Does the quadrature compute the expected BD operator of a causal set in the geometry Π²η correctly: the proper volume of intervals, convergence of the quadrature and the local limit □_g − R/2? | calibration |
| K.11b | Does the volume-weighted remainder of the curved operator in the far rows (an effective spurious source) converge in the quadrature to the level that F1.4 and F1.5 need? | calibration |
| K.11c | Does the volume-weighted remainder of the far rows of the curved operator in 3+1 converge when the geometry between the nodes is interpolated by a C² cubic spline instead of piecewise linearly? | calibration |
| K.12 | Does the geodesic in the self-consistent field (a) of a causal set in 2+1 approach the geodesic of the continuum solution as the source grows relative to the non-local scale ℓ_k (σ = 8, 16, 32 in a self-similar family), and from which σ is the deviation within the tolerance of a test of the F2.7 type? Does the position of the wall change it? | calibration |
| K.13a | Does the discrete construction (Poisson sprinkling of the geometry −Π²dt² + Π⁻²dx², the layers and smearing of K.4/K.5) reproduce at first order the continuum limit □_g − R/2 on fixed test fields, including the term −R/2 = −∇²u, and converge to it as ℓ_k/σ → 0? | calibration |
| K.14 | Is the number of horizon molecules (Barton et al. 2019) in flat spacetime equal to a constant a^(d) times the area in units of the discreteness, and in a conformally flat spacetime proportional to the PROPER area of the horizon cut? | calibration |
| K.15 | Does T^ab = (1/ΔV₄) Σ F u^a u^b τ_link on populations with a known answer (a beam, counter-streaming beams, an isotropic gas) give the correct energy density, momentum flux, pressure and shear, and is the 4-momentum conserved on slices under the merge rule F = /P/? | calibration |
| K.16 | Does the instrument determine from the events, the relation and the density (without φ, ψ and γ) γ = −1, 0 and +1 of a weak field in 2+1 along two independent paths (cones + density, clocks + cones)? | calibration |
| K.16b | Do a pooled boundary of the relation (instead of the middle of the gap at each anchor) and smaller diamonds remove the systematic shifts of K.16, and does the instrument determine γ = −1, 0, +1 within the tolerance 0.25 along both paths? | calibration |
| K.16c | Does the K.16b instrument carried over to 3+1 (density → φ − 3ψ, diamond πτ⁴/24) determine γ = −1, 1/3 and +1 of a weak field blindly along both paths in the box t ∈ [0, 9], /x/ <= 9? | calibration |
| K.17 | Is the noise of the discrete BD action under correlated sampling small enough to measure the quadratic (Einstein–Hilbert) response to a static field Ψ in the box of the laboratory to 10 % within 1000 sprinklings? | calibration |
| K.18 | Does the needle-volume instrument (plane wave, q^2 series to q^6, exact needle) give the interval volume of the optical sector to the working accuracy in the range eps <= 4 that F6.8f uses? | calibration |
| K.18b | Does the needle-volume instrument (plane wave, q^2 series to q^6, exact needle) give the interval volume of the optical sector to the working accuracy in the range eps <= 4 that F6.8f uses? | calibration |
