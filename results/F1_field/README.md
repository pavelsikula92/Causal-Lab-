# F1: Transport and field

## Objective

How influence is transported along the order and conserved at every event, and what field a source produces: Gauss's law, the far field, superposition and its failure at second order, and the role of capacity.

## Contents

One folder per test, `F1_field/<ID>/`, with the records of its full runs: `<UTC time>_<part>_seed<seed>_<test hash 8>.json` (parameters, seed, hashes of the test file and of the core, versions, the printed lines and the data) and the same lines as `.txt`. Reduced runs are in `<ID>/quick/` and do not count. Nothing is overwritten.

The verdicts and their reading are in `docs/findings.md`; the list of tests with their latest states in `docs/test_list.md` (`python run.py report`); the records of one test: `python run.py results <ID>`.

| Test | Question | Level |
|---|---|---|
| F1.2 | Does the aggregated BD operator of one order applied to the field of a source approach the expected operator, with an error falling as (L/ℓ)^(−p) and the exponent p from an independent calibration? | mechanism |
| F1.2b | Does the aggregated BD operator of one order applied to the field of a source approach the expected operator, with an error falling as (L/ℓ)^(−p) and the exponent p from an independent calibration? | mechanism |
| F1.3 | Does the field of a source solved with the expected non-local BD operator have, for r ≫ l_k, the form of the Green function of the Laplace operator (1/r in 3+1, ln r in 2+1) with an amplitude given by the total source (Gauss's law)? | mechanism |
| F1.4 | Do the fields of two sources add, or does a second-order superposition deficit arise exactly according to the non-linearity of equation (a)? | mechanism |
| F1.5 | Is the field of a weak source proportional to the source with the coupling κ/D (effect of the term −R/2), and does the field strength of a stronger source grow more slowly according to the non-linearity of equation (a)? | mechanism |
| F1.7 | What profile does the field around a source have? | undecided: the load was measured, not Π |
| F1.9a | Does the discrete BD operator in a 1/r field outside the source keep the zero response of H1 (c = 0, β = 1) and the non-zero response of the conformal control (c = 2), and does it converge to them with σ/ℓ_k? | mechanism |
