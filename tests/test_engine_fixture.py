"""The JS-engine fixture (tests/fixtures/engine_fixture.json) is current: regenerating it from the Python code gives the
same values (numbers compared with rtol 1e-9, since libm can differ in the last digit between platforms).

If this fails after changing an analysis, run: python scripts/make_engine_fixture.py
"""
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import make_engine_fixture  # noqa: E402


def same(a, b, path="$"):
    if isinstance(a, (int, float)) and isinstance(b, (int, float)) and not isinstance(a, bool):
        assert math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-12), f"{path}: {a} != {b}"
    elif isinstance(a, dict):
        assert set(a) == set(b), path
        for k in a:
            same(a[k], b[k], f"{path}.{k}")
    elif isinstance(a, list):
        assert len(a) == len(b), path
        for i, (x, y) in enumerate(zip(a, b)):
            same(x, y, f"{path}[{i}]")
    else:
        assert a == b, f"{path}: {a!r} != {b!r}"


def test_fixture_is_current():
    stored = json.loads(make_engine_fixture.OUT.read_text())
    fresh = json.loads(json.dumps(make_engine_fixture.make(), ensure_ascii=False))
    same(fresh, stored)
