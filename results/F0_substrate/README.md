# F0: Substrate

## Objective

What the substrate is: events and their order. Whether a substrate with a finite number of neighbours per event (and therefore a rest frame) is viable, or whether a Lorentz-invariant causal set (a Poisson sprinkling) is required; geometry read from order and number.

## Contents

One folder per test, `F0_substrate/<ID>/`, with the records of its full runs: `<UTC time>_<part>_seed<seed>_<test hash 8>.json` (parameters, seed, hashes of the test file and of the core, versions, the printed lines and the data) and the same lines as `.txt`. Reduced runs are in `<ID>/quick/` and do not count. Nothing is overwritten.

The verdicts and their reading are in `docs/findings.md`; the list of tests with their latest states in `docs/test_list.md` (`python run.py report`); the records of one test: `python run.py results <ID>`.

| Test | Question | Level |
|---|---|---|
| F0.6 | Does a rule that draws cones in coordinates give the same physics when it is evaluated in a frame moving at velocity β (the same events, the same order, driving and readout in the rest frame of the source)? Is any difference small compared with what the rule itself does? | mechanism |
