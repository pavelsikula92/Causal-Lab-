"""
core.registry
=============

The phases of the roadmap and the discovery of test files.

A test file lives in phases/<folder>/ and is named <ID with _ for .>_<name>.py,
e.g. phases/F1_field/F1_02_integral_A5.py; it must define TEST and run().
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Dict, List, Tuple

from .testing import validate

ROOT = Path(__file__).resolve().parents[1]
PHASES_DIR = ROOT / "phases"

# (phase, folder, name, part of the theory)
PHASES: List[Tuple[str, str, str, str]] = [
    ("K", "K_instruments", "Instruments and calibration", "cross-cutting; instruments without a physical conclusion"),
    ("F0", "F0_substrate", "Substrate", "axioms 1, 2, 8 – events, order, dimension from the order"),
    ("F1", "F1_field", "Transport and field", "axiom 3 (flow), axiom 4 (flow capacity), axiom 5 and A5′"),
    ("F2", "F2_motion", "Motion in a field", "axiom 6 – action and maximal chains"),
    ("F3", "F3_structure", "Formation of structure", "axiom 4 and the non-linear term A5′ (problems 9, 13)"),
    ("F4", "F4_inertia", "Inertia", "the functional R = Σ Z Δτ; mass = internal frequency"),
    ("F5", "F5_gravity", "Gravity and equivalence", "the central claim: one origin of gravity and inertia"),
    ("F6", "F6_order_dynamics", "Dynamics of the order", "axiom 7, D_eff ≈ 3 (problems 2, 5, 8)"),
    ("F7", "F7_light", "Light and waves", "hypotheses H1, H2 (problem 4)"),
    ("F8", "F8_data", "Confrontation with data", "problems 6, 7, 10 and the solar system"),
]

_FILE = re.compile(r"^(K|F\d)_\d{2}[a-z]?\d?_[A-Za-z0-9_]+\.py$")


def _load(path: Path):
    """Import a test as phases.<folder>.<file> (a real module name: worker processes can import it)."""
    import importlib
    return importlib.import_module(f"phases.{path.parent.name}.{path.stem}")


def discover() -> Dict[str, Dict[str, Any]]:
    """{test id: {"path", "module", "test"}} for every valid test file; raises on duplicates or bad headers."""
    out: Dict[str, Dict[str, Any]] = {}
    problems: List[str] = []
    for _, folder, _, _ in PHASES:
        d = PHASES_DIR / folder
        if not d.exists():
            continue
        for path in sorted(d.glob("*.py")):
            if not _FILE.match(path.name):
                continue
            mod = _load(path)
            test = getattr(mod, "TEST", None)
            if test is None or not hasattr(mod, "run"):
                problems.append(f"{path.name}: TEST or run() missing")
                continue
            problems += validate(test, path.name)
            tid = test.get("id")
            if tid in out:
                problems.append(f"{path.name}: duplicate id {tid}")
            out[tid] = {"path": path, "module": mod, "test": test}
    if problems:
        raise SystemExit("invalid test files:\n  " + "\n  ".join(problems))
    return out


def phase_of(tid: str) -> str:
    return tid.split(".")[0]


def sort_key(tid: str):
    """K first, then F0..F8; numbers numerically; revision letters last."""
    ph, rest = tid.split(".", 1)
    m = re.match(r"(\d+)([a-z]?)", rest)
    order = [p[0] for p in PHASES]
    return (order.index(ph) if ph in order else 99, int(m.group(1)) if m else 0, m.group(2) if m else "")
