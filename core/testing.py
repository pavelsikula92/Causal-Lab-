"""
core.testing
============

The contract between a test file, the runner and the results.

A test is ONE file in phases/<phase folder>/, e.g. phases/F1_pole/F1_02_integral_A5.py.
It defines

    TEST = {
        "id":          "F1.2",               # phase + number; a revision adds a letter: "F1.2b"
        "phase":       "F1",
        "title":       "...",
        "theory_part": "...",                # which part of the theory it tests
        "question":    "...",
        "method":      "...",
        "prediction":  "...",                # and where it comes from
        "criterion":   "...",                # fixed BEFORE the run; never changed after it
        "null_model":  "...",
        "level":       "mechanism",          # calibration | mechanism | chain | prediction
        "meters":      ["K.4", "K.5"],       # calibration tests whose instruments it uses
        "depends_on":  ["F1.1"],
        "old_id":      "E45d",               # designation in GTS Lab 3.x, if any
        "params":      {...},                # common parameters (criterion thresholds included)
        "parts":       {"d3": {...}, ...},   # optional: parameter sets run separately
        "quick":       {...},                # overrides for --quick (a smoke run, decides nothing)
        "settable":    {"search": ["lo"]},   # optional: parameters a FULL run of a part may take from --set
        "estimate":    "...",                # expected running time
    }

    def run(p, ctx) -> Result

`p` is params, updated by the part's parameters, then (quick runs) by
`quick` and any --set overrides, or (full runs) by --set overrides of the
parameters the part declares in `settable` only. `settable` is for
exploratory parts (a scan range chosen at the command line) and for values
chosen from earlier data; it never holds a criterion threshold. The values
used are saved with the result like every other parameter. `ctx` gives the part, the base seed, quick,
the number of worker processes and progress(). The function returns a
Result with one of the statuses below, the lines printed to the terminal
and the data saved as JSON.

Statuses
    PASSED             the criterion is met (for a calibration test: the
                       instrument is calibrated in the declared range)
    FAILED             the criterion is not met
    UNDECIDED          the resolution does not decide, or the null model
                       meets the criterion too (not specific)
    INSTRUMENT_FAILED  an instrument check failed: no statement
    RECORDED           an exploratory part: numbers recorded, no verdict

Rules (docs/roadmapa.md): the criterion lives in the test file and is not
changed after a run (a changed criterion = a new revision with a new id,
the old result stays valid under its own criterion); a verdict counts as
replicated only after a run with another seed (--seed K); the level of the
conclusion limits what the test may claim.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional

STATUSES = ("PASSED", "FAILED", "UNDECIDED", "INSTRUMENT_FAILED", "RECORDED")
STATUS_CS = {"PASSED": "PROŠEL", "FAILED": "NEPROŠEL", "UNDECIDED": "NEROZHODNUTO",
             "INSTRUMENT_FAILED": "MĚŘIDLO SELHALO", "RECORDED": "ZAZNAMENÁNO"}
LEVELS = ("calibration", "mechanism", "chain", "prediction")
LEVEL_CS = {"calibration": "kalibrace", "mechanism": "mechanismus", "chain": "řetězec", "prediction": "predikce"}
REQUIRED = ("id", "phase", "title", "theory_part", "question", "method", "prediction", "criterion",
            "null_model", "level", "meters", "depends_on", "params", "quick", "estimate")


@dataclass
class Result:
    status: str
    lines: List[str] = field(default_factory=list)
    data: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        if self.status not in STATUSES:
            raise ValueError(f"unknown status {self.status}")


@dataclass
class Context:
    part: str
    seed: int
    quick: bool
    workers: int
    progress: Callable[[str], None]


def validate(test: Dict[str, Any], path: str = "") -> List[str]:
    """Problems of a TEST header (empty list = valid)."""
    bad = [f"missing field '{k}'" for k in REQUIRED if k not in test]
    if "level" in test and test["level"] not in LEVELS:
        bad.append(f"level must be one of {LEVELS}")
    if "id" in test and "phase" in test and not str(test["id"]).startswith(str(test["phase"]) + "."):
        bad.append("id must start with the phase and a dot (e.g. F1.2)")
    for k in ("meters", "depends_on"):
        if k in test and not isinstance(test[k], list):
            bad.append(f"'{k}' must be a list")
    if "parts" in test and not isinstance(test["parts"], dict):
        bad.append("'parts' must be a dict")
    if "settable" in test:
        st = test["settable"]
        names = set(test.get("parts") or {}) or {""}
        if not isinstance(st, dict) or not all(isinstance(v, list) for v in st.values()):
            bad.append("'settable' must be a dict part -> list of parameter names")
        else:
            bad += [f"'settable' names an unknown part '{k}'" for k in st if k not in names]
    return [f"{path}: {b}" if path else b for b in bad]


def settable_of(test: Dict[str, Any], part: Optional[str]) -> List[str]:
    """Parameters a full run of this part may take from --set (empty for most tests)."""
    return list((test.get("settable") or {}).get(part or "", []))


def check_overrides(test: Dict[str, Any], part: Optional[str], quick: bool,
                    overrides: Optional[Dict[str, Any]]) -> Optional[str]:
    """None if the --set overrides are allowed, otherwise the reason (quick: anything; full: settable only)."""
    if not overrides or quick:
        return None
    allowed = settable_of(test, part)
    bad = sorted(k for k in overrides if k not in allowed)
    if not bad:
        return None
    where = f"part '{part}'" if part else "this test"
    return (f"--set {', '.join(bad)} is not allowed in a full run of {where} (settable: {allowed or 'none'}); "
            f"other parameters change only with --quick")


def merged_params(test: Dict[str, Any], part: Optional[str], quick: bool,
                  overrides: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """params, then the part's parameters, then (quick) the quick overrides and --set, or (full) the allowed --set."""
    why = check_overrides(test, part, quick, overrides)
    if why:
        raise ValueError(why)
    p = copy.deepcopy(test["params"])
    if part:
        p.update(copy.deepcopy(test.get("parts", {})[part]))
    if quick:
        p.update(copy.deepcopy(test.get("quick", {})))
        qp = test.get("quick_parts", {}).get(part or "", {})
        p.update(copy.deepcopy(qp))
    if overrides:
        p.update(overrides)
    return p


def parts_of(test: Dict[str, Any]) -> List[Optional[str]]:
    """The parts of a test, or [None] for a test without parts."""
    parts = list(test.get("parts", {}) or {})
    return parts if parts else [None]
