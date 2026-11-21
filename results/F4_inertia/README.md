# F4: Inertia and contact

## Objective

How matter resists and exchanges motion: flows that meet merge with their four-momentum conserved.

## Contents

One folder per test, `F4_inertia/<ID>/`, with the records of its full runs: `<UTC time>_<part>_seed<seed>_<test hash 8>.json` (parameters, seed, hashes of the test file and of the core, versions, the printed lines and the data) and the same lines as `.txt`. Reduced runs are in `<ID>/quick/` and do not count. Nothing is overwritten.

The verdicts and their reading are in `docs/findings.md`; the list of tests with their latest states in `docs/test_list.md` (`python run.py report`); the records of one test: `python run.py results <ID>`.

| Test | Question | Level |
|---|---|---|
| F4.1 | When two flows with carried 4-velocity meet head-on, do their chains merge (contact transfers momentum, the merged flow comes to rest in the CM frame and gains rest mass), or do the flows pass through each other? | mechanism |
