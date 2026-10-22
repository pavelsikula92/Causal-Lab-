"""
Causality Lab -- start a test from a terminal (no installation needed).

    python run.py list
    python run.py F1.2 --part d3
    python run.py F1.2 --quick
    python run.py core
    python run.py help           all commands
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))


def main() -> int:
    from core.cli import main as cli
    return cli(sys.argv[1:])


if __name__ == "__main__":
    raise SystemExit(main())
