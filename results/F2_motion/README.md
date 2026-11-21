# F2: Motion in a field

## Objective

How a history moves in a given order: whether maximal (longest) chains follow the geodesics of the effective geometry, in 2+1 and 3+1.

## Contents

One folder per test, `F2_motion/<ID>/`, with the records of its full runs: `<UTC time>_<part>_seed<seed>_<test hash 8>.json` (parameters, seed, hashes of the test file and of the core, versions, the printed lines and the data) and the same lines as `.txt`. Reduced runs are in `<ID>/quick/` and do not count. Nothing is overwritten.

The verdicts and their reading are in `docs/findings.md`; the list of tests with their latest states in `docs/test_list.md` (`python run.py report`); the records of one test: `python run.py results <ID>`.

| Test | Question | Level |
|---|---|---|
| F2.5 | Do the heaviest chains Σ Π in 3+1 follow the geodesic of Π²η and Σ Π² the geodesic of Π⁴η, clearly separated (the form of the action determines the path as in 2+1)? And is the ~10 % excess of Σ Π in 2+1 an effect of finite N? | mechanism |
| F2.7 | Do the heaviest chains Σ Π move, in the field that the causal set itself creates from a source by equation (a), along the geodesic predicted by the continuum chain of the theory, and is this distinguishable from competing field equations and from a zero field? | chain |
