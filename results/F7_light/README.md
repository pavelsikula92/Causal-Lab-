# F7: Light and the response of the order

## Objective

How the response of the order fixes the post-Newtonian gamma: fixed order, narrowing cones, conservation of causal throughput; and the tensor coupling (pressure and motion gravitate, cones tilt).

## Contents

One folder per test, `F7_light/<ID>/`, with the records of its full runs: `<UTC time>_<part>_seed<seed>_<test hash 8>.json` (parameters, seed, hashes of the test file and of the core, versions, the printed lines and the data) and the same lines as `.txt`. Reduced runs are in `<ID>/quick/` and do not count. Nothing is overwritten.

The verdicts and their reading are in `docs/findings.md`; the list of tests with their latest states in `docs/test_list.md` (`python run.py report`); the records of one test: `python run.py results <ID>`.

| Test | Question | Level |
|---|---|---|
| F7.1a | When a loaded event keeps only the links with speed <= c(σ) (the cones narrow, the events stay), what γ does the blind instrument read? Prediction from the derivation: 1/2 for any function c(σ), because the density does not change. | mechanism |
| F7.1b | When the cones narrow to c_eff = exp(−σ) and the density of events grows to S = 1/c_eff (throughput S·c_eff conserved), what γ does the blind instrument read? Prediction in 2+1: γ = 2 (in 3+1 it would be 1). | mechanism |
| F7.1c | When in 3+1 the cones narrow to c_eff = exp(−σ) and the density of events grows to exp(σ) (throughput conserved), does the blind instrument K.16c read γ = 1? Alternatives: −1 (fixed order), 1/3 (narrower cones only). | mechanism |
| F7.4a | When each flow narrows the cone in its own rest frame and the capacity deficits add, does the blind instrument read for two counter-streaming flows (±0.6c) a clock field Γ_v²(1 + v²) relative to rest (tensor) and anisotropic cones, and not Γ_v², Γ_v or 1 (scalar rules of the mean frame)? | mechanism |
