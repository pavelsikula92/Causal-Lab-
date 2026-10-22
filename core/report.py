"""
core.report
===========

python run.py report  ->  docs/seznam_testu.md

One list of all tests, generated so that it never disagrees with the code:

  1. implemented tests: from the TEST headers in phases/ and the latest
     saved result of each part (replication and "older test version" marks);
  2. planned tests: docs/planovane_testy.csv (edited by hand; a row
     disappears from this section as soon as a test file with that id exists);
  3. results of GTS Lab 2.x and 3.x: results/historie_3x/prehled.csv
     (the old designations and their results, assigned to the new ids).
"""

from __future__ import annotations

import csv
import time
from pathlib import Path
from typing import Dict, List

from . import registry, provenance
from .cli import status_line
from .testing import LEVEL_CS

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"


def _csv(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        return []
    with path.open(encoding="utf-8", newline="") as f:
        return [dict(r) for r in csv.DictReader(f, delimiter=";")]


def _cell(s: str) -> str:
    return str(s).replace("|", "/").replace("\n", " ")


def write_report() -> Path:
    tests = registry.discover()
    names = {p[0]: p[2] for p in registry.PHASES}
    out = ["# Seznam testů Causality Lab", "",
           f"Vygenerováno příkazem `python run.py report` ({time.strftime('%Y-%m-%d %H:%M')}). Neupravuj ručně: "
           "hlavičky testů jsou v `phases/`, plán v `docs/planovane_testy.csv`, historie v "
           "`results/historie_3x/prehled.csv`.", "",
           "## Implementované testy", "",
           "| ID | Fáze | Testovaná otázka | Úroveň závěru | Měřidla | Závislosti | Stav (poslední výsledek) |",
           "| --- | --- | --- | --- | --- | --- | --- |"]
    for tid in sorted(tests, key=registry.sort_key):
        t = tests[tid]["test"]
        old = f" ({t['old_id']})" if t.get("old_id") else ""
        out.append(f"| {tid}{old} | {registry.phase_of(tid)} {names.get(registry.phase_of(tid), '')} | {_cell(t['question'])} | "
                   f"{LEVEL_CS[t['level']]} | {', '.join(t['meters']) or '—'} | {', '.join(t['depends_on']) or '—'} | "
                   f"{_cell(status_line(tid, tests[tid]))} |")
    plan = [r for r in _csv(DOCS / "planovane_testy.csv") if r.get("ID") not in tests]
    out += ["", "## Plánované testy", "",
            "| ID | Fáze | Testovaná otázka | Úroveň závěru | Závislosti | Poznámka |", "| --- | --- | --- | --- | --- | --- |"]
    for r in sorted(plan, key=lambda r: registry.sort_key(r["ID"])):
        out.append(f"| {r['ID']} | {r['ID'].split('.')[0]} {names.get(r['ID'].split('.')[0], '')} | {_cell(r['Otázka'])} | "
                   f"{_cell(r['Úroveň'])} | {_cell(r['Závislosti'])} | {_cell(r['Poznámka'])} |")
    hist = _csv(ROOT / "results" / "historie_3x" / "prehled.csv")
    out += ["", "## Výsledky GTS Lab 2.x a 3.x (historie)", "",
            "Testy samotné jsou v archivu; přežívají jejich výsledky. Uložené záznamy 3.x patří do `results/historie_3x/`.", "",
            "| Nové ID | Původní označení | Testovaná otázka | Výsledek | Úroveň závěru |", "| --- | --- | --- | --- | --- |"]
    for r in sorted(hist, key=lambda r: registry.sort_key(r["Nové ID"])):
        out.append(f"| {r['Nové ID']} | {_cell(r['Původní'])} | {_cell(r['Otázka'])} | {_cell(r['Výsledek'])} | {_cell(r['Úroveň'])} |")
    path = DOCS / "seznam_testu.md"
    path.write_text("\n".join(out) + "\n", encoding="utf-8")
    return path
