"""
core.provenance
===============

Every run is saved with what is needed to reproduce it and to see whether
the code changed since:

    test_sha256   hash of the test file
    core_sha256   hash of all files of the core (sorted by path)
    params        the parameters actually used, part, seed, quick
    versions      python, numpy, scipy, platform, laboratory version
    time, duration

Files: results/<ID>/<UTC time>_<part>_<seed>_<test hash 8>.json and a
.txt with the printed lines; quick runs go to results/<ID>/quick/ and are
never counted. Nothing is overwritten.
"""

from __future__ import annotations

import hashlib
import json
import math
import platform
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"


def sha256_file(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def core_sha256() -> str:
    """One hash over all core/*.py files, sorted by relative path."""
    h = hashlib.sha256()
    for p in sorted((ROOT / "core").rglob("*.py")):
        h.update(str(p.relative_to(ROOT)).encode())
        h.update(p.read_bytes())
    return h.hexdigest()


def versions() -> Dict[str, str]:
    import numpy
    import scipy
    from . import __version__
    return {"python": sys.version.split()[0], "numpy": numpy.__version__, "scipy": scipy.__version__,
            "platform": platform.platform(), "laboratory": __version__}


def _clean(x: Any) -> Any:
    """JSON-safe copy: numpy types and arrays to Python, nan/inf to None."""
    try:
        import numpy as np
    except ImportError:                     # pragma: no cover
        np = None
    if isinstance(x, dict):
        return {str(k): _clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_clean(v) for v in x]
    if np is not None:
        if isinstance(x, np.ndarray):
            return _clean(x.tolist())
        if isinstance(x, np.generic):
            x = x.item()
    if isinstance(x, float) and (math.isnan(x) or math.isinf(x)):
        return None
    return x


def save(test_id: str, test_file: Path, part: Optional[str], seed: int, quick: bool, params: Dict[str, Any],
         status: str, lines: List[str], data: Dict[str, Any], t_start: float,
         hashes: Optional[Dict[str, str]] = None) -> Path:
    """
    Write the JSON record and the text lines; returns the JSON path.
    hashes: {"test_sha256", "core_sha256"} taken when the run STARTED (the command line passes them); without
    them the files are hashed now. If a file changed during the run, the record says so (files_changed_during_run).
    """
    now_t, now_c = sha256_file(test_file), core_sha256()
    tsha = hashes["test_sha256"] if hashes else now_t
    csha = hashes["core_sha256"] if hashes else now_c
    changed = bool(hashes) and (now_t != tsha or now_c != csha)
    stamp = time.strftime("%Y%m%dT%H%M%SZ", time.gmtime())
    folder = RESULTS / test_id / ("quick" if quick else "")
    folder.mkdir(parents=True, exist_ok=True)
    base = f"{stamp}_{part or 'all'}_seed{seed}_{tsha[:8]}"
    path = folder / f"{base}.json"
    k = 1
    while path.exists():                    # never overwrite
        path = folder / f"{base}_{k}.json"
        k += 1
    rec = {"id": test_id, "part": part, "seed": seed, "quick": quick, "status": status,
           "test_file": str(Path(test_file).relative_to(ROOT)), "test_sha256": tsha, "core_sha256": csha,
           "files_changed_during_run": changed,
           "versions": versions(), "time_utc": stamp, "duration_s": round(time.time() - t_start, 1),
           "params": _clean(params), "lines": lines, "data": _clean(data)}
    path.write_text(json.dumps(rec, indent=1, ensure_ascii=False), encoding="utf-8")
    path.with_suffix(".txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def load_all(test_id: str, include_quick: bool = False) -> List[Dict[str, Any]]:
    """All saved records of a test, oldest first."""
    folder = RESULTS / test_id
    if not folder.exists():
        return []
    files = sorted(folder.glob("*.json"))
    if include_quick:
        files += sorted((folder / "quick").glob("*.json"))
    out = []
    for f in files:
        try:
            rec = json.loads(f.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue                          # an unreadable file is skipped, never fatal
        rec["_file"] = str(f.relative_to(ROOT)) if f.is_relative_to(ROOT) else str(f)
        out.append(rec)
    return out
