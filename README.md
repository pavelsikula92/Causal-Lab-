# Causality Lab

Open numerical laboratory for **Causal Physics**: a programme that takes the causal order of events, not spacetime
geometry, as the primitive of physics, and asks how much of gravitation follows from it. The code samples causal sets
(Poisson sprinklings), reads geometry back from order and number with blind instruments, transports conserved flows
with carried four-velocity, and evaluates the Benincasa–Dowker operator and action of a causal set.

Every test states its question, prediction, null model and numerical criterion before the full run, and ends with one
of four verdicts: passed, failed, undecided or instrument failed. Negative results are kept and a failed test is
never rescued by changing its criteria. The tests are ordered by the phases of the roadmap, and every phase belongs to
one part of the theory (`docs/roadmap.md`). The paper that reports the tests is in `paper/`.

This is an attempt to apply QA methods to scientific research. I am primarily a test engineer, and natural science is
my passion. The work is not peer reviewed; independent replication and criticism are the reason it is open.

## Citation

The paper and the code used for it are archived together on Zenodo:
Pavel Šikula, *Physics of Causality*, https://doi.org/10.5281/zenodo.23853689

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

## Reproducing the main result

The central computation of the paper is the second-order expectation of the smeared Benincasa–Dowker action in
three sectors (deterministic quadrature, no sprinklings). Run in this order; each test reads the recorded results of
the previous ones, and the saved results of the author's runs are in `results/`:

    python run.py F6.8a      conformal sector (EH value 6)        ~30 CPU-min, 4 min recorded with several processes
    python run.py F6.8b2     mixed sector (EH value 8)            ~5-10 CPU-hours, 62 min recorded with many processes
    python run.py F6.8c2     optical sector (EH value 2)          ~1-2 h with 10 processes, 103 min recorded

Compare with the recorded output (`results/<ID>/*.txt`) and with the numbers in the paper (Part V, Sec. 36-37).
Known limits of this result are listed in the paper (Sec. 52) and in `docs/findings.md`; in particular the optical
coefficient depends on how the age dependence is read (fixed `kT`, extrapolation from four wavenumbers).

Some records are kept for honesty and are **not valid results**. F6.8c (optical sector read linearly in the age)
FAILED, and F8.1 (point-mass field) was run on that invalid input, so its recorded numbers carry no information and
are not used in the paper; F8.1 waits for the asymptotic age term of the optical sector (F6.8e).

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
    Demo/          Causal Gravity, an interactive simulator of the gravity law (stand-alone; scenes explained in
                   plain words in Demo/README.md)
    diag/          development diagnostics (decide nothing; not registered tests)
    core_tests/    quick checks of the core

The header of every test (`TEST = {...}`) contains the phase, the part of the theory, the question, the method, the
prediction, the criterion, the null model, the level of the conclusion, the instruments used and the dependencies.
The criterion does not change after the run; a change of the criterion is a new revision of the test with a new ID
(e.g. F1.2b).

The findings, with every negative result, are in `docs/findings.md`; the plan in `docs/planned_tests.csv`. The
results of GTS Lab 2.x and 3.x are converted to the new IDs in `results/history_3x/overview.csv`.

## Licence

Code: MIT (`LICENSE`). Paper: CC BY 4.0.
