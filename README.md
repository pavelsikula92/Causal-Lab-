# Causality Lab

The laboratory of Causal Physics (Pavel Šikula). It tests, on sprinkled causal orders, whether gravity and inertia are
two manifestations of one dynamics of causal resistance, and whether the geometry of general relativity can be read
from a causal order and its action. The tests are ordered by the phases of the roadmap, and every phase belongs to one
part of the theory (`docs/roadmap.md`). The paper that reports the finished tests is in `paper/`.

## Installation

Python 3.10 or newer, with numpy and scipy (matplotlib only for the figures of the paper):

    pip install numpy scipy matplotlib
    python run.py core          quick checks of the core (about 10 s)

Memory: one process needs up to ~1.5 GB of RAM depending on the test; use at most as many processes as
RAM / 1.5 GB, with a reserve for the system.

## Running

    python run.py list                     phases, tests and their latest status
    python run.py show F1.2                question, prediction, criterion, parameters of a test
    python run.py F1.2                     a test, all its parts
    python run.py F1.2 --part d3           one part only
    python run.py F1.2 --quick             reduced run: checks that the test runs (decides nothing)
    python run.py F1.2 --seed 7            independent replication with another seed
    python run.py F1.2 --workers 6         number of processes (otherwise CAUSALITY_WORKERS or all cores)
    python run.py results F1.2             saved results of a test
    python run.py report                   regenerates docs/test_list.md
    python paper/figures/make_figures.py   draws the figures of the paper from the saved results

## How to read a result

Every test ends with one status: **PASSED** (for an instrument: calibrated in the declared range), **FAILED**,
**UNDECIDED** (the resolution is not sufficient, or the null model meets the criterion as well),
**INSTRUMENT_FAILED** (a check of the instrument failed; no conclusion), **RECORDED** (numbers recorded only).

The level of the conclusion (calibration, mechanism, chain, prediction) says what a result may claim. A result is
saved in `results/<ID>/` as JSON and text, with the hash of the test file and of the core, the parameters, the seed
and the versions of the libraries; nothing is overwritten. Reduced runs go to `results/<ID>/quick/` and do not count.
`results/sha_aliases.json` maps the hashes of test files from before the English translation of 30 Oct 2026 to the
translated files where only text changed (checked by comparing the syntax trees without string constants), so the
records stay current.

## Structure

    core/          the core: order, flow, operators, chains, instruments, statistics, provenance, test framework
    phases/        tests by phase (K_instruments, F0_substrate, F1_field, …); one test = one file
    results/       results; history_3x/ for the results of GTS Lab 2.x and 3.x
    docs/          roadmap, findings, planned tests, author's decisions, derivations, generated list of tests
    paper/         the paper (LaTeX source, figures and the script that draws them)
    Demo/          Causal Gravity, an interactive simulator of the gravity law (stand-alone)
    core_tests/    quick checks of the core

The header of every test (`TEST = {...}`) contains the phase, the part of the theory, the question, the method, the
prediction, the criterion, the null model, the level of the conclusion, the instruments used and the dependencies.
The criterion does not change after the run; a change of the criterion is a new revision of the test with a new ID
(e.g. F1.2b).

The findings, with every negative result, are in `docs/findings.md`; the plan in `docs/planned_tests.csv`. The
results of GTS Lab 2.x and 3.x are converted to the new IDs in `results/history_3x/overview.csv`.

## Licence

Code: MIT (`LICENSE`). Paper: CC BY 4.0.
