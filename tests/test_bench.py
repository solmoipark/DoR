"""G4: every DoR bench scenario passes (mock)."""

from __future__ import annotations

import pytest

pytest.importorskip("inverse_gems")

from dorgems.pilot.bench import run_bench  # noqa: E402


def test_bench_all_pass(tmp_path, bundles_dir):
    rep = run_bench(out=tmp_path)
    failed = [r for r in rep["results"] if not r["ok"]]
    assert not failed, failed
    assert rep["n"] == 6
