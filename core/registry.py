"""
core.registry
=============

The phases of the roadmap and the discovery of test files.

A test file lives in phases/<folder>/ and is named <ID with _ for .>_<name>.py,
e.g. phases/F1_pole/F1_02_integral_A5.py; it must define TEST and run().
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Dict, List, Tuple

from .testing import validate

ROOT = Path(__file__).resolve().parents[1]
PHASES_DIR = ROOT / "phases"

# (phase, folder, Czech name, part of the theory)
PHASES: List[Tuple[str, str, str, str]] = [
    ("K", "K_meridla", "Měřidla a kalibrace", "průřezová; měřidla bez fyzikálního závěru"),
    ("F0", "F0_substrat", "Substrát", "axiomy 1, 2, 8 – události, řád, dimenze z řádu"),
    ("F1", "F1_pole", "Transport a pole", "axiom 3 (tok), axiom 4 (kapacita toku), axiom 5 a A5′"),
    ("F2", "F2_pohyb", "Pohyb v poli", "axiom 6 – akce a maximální řetězce"),
    ("F3", "F3_struktura", "Vznik struktury", "axiom 4 a nelineární člen A5′ (problémy 9, 13)"),
    ("F4", "F4_setrvacnost", "Setrvačnost", "funkcionál R = Σ Z Δτ; hmotnost = vnitřní frekvence"),
    ("F5", "F5_gravitace", "Gravitace a ekvivalence", "ústřední tvrzení: jeden původ gravitace a setrvačnosti"),
    ("F6", "F6_dynamika_radu", "Dynamika řádu", "axiom 7, D_eff ≈ 3 (problémy 2, 5, 8)"),
    ("F7", "F7_svetlo", "Světlo a vlny", "hypotézy H1, H2 (problém 4)"),
    ("F8", "F8_data", "Konfrontace s daty", "problémy 6, 7, 10 a sluneční soustava"),
]

_FILE = re.compile(r"^(K|F\d)_\d{2}[a-z]?_[A-Za-z0-9_]+\.py$")


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
