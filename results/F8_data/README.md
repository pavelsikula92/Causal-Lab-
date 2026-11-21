# F8: Confrontation with data

## Objective

Predictions of the theory compared with observation: the static field of a point mass near the non-locality scale, the solar system, and the astrophysical consequences.

## Contents

One folder per test, `F8_data/<ID>/`, with the records of its full runs: `<UTC time>_<part>_seed<seed>_<test hash 8>.json` (parameters, seed, hashes of the test file and of the core, versions, the printed lines and the data) and the same lines as `.txt`. Reduced runs are in `<ID>/quick/` and do not count. Nothing is overwritten.

The verdicts and their reading are in `docs/findings.md`; the list of tests with their latest states in `docs/test_list.md` (`python run.py report`); the records of one test: `python run.py results <ID>`.

| Test | Question | Level |
|---|---|---|
| F8.1 | What deviation from Newton's law and from gamma = 1 does the action of the order predict for the static field of a point mass, with the coefficients of all three sectors computed and nothing fitted to gravity data, and what does the torsion-balance bound on power-law potentials imply for l_k for each reading of the age T? | prediction |
