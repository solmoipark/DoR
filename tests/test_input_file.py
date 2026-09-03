"""New-material input template: validation report, grades, and the --input CLI paths (mock)."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

from dorgems.pilot.input_file import load_input_file, validate_report

ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = ROOT / "templates" / "new_material_template.yaml"
EXAMPLE = ROOT / "templates" / "example_ggbs_pohang.yaml"


def test_raw_template_is_rejected_with_clear_errors():
    r = validate_report(TEMPLATE)
    assert not r["ok"] and r["errors"]


def test_example_validates_and_grades():
    r = validate_report(EXAMPLE)
    assert r["ok"], r
    s = r["summary"]
    assert s["slot"] == "slag" and s["can_run"] == {"A_envelope": True, "B_compare": True, "C_infer": True}
    by = {(g["age_d"], g["quantity"], g["unit"]): g for g in r["grades"]}
    assert by[(28.0, "bound_water", "g/100 g binder")]["grade"] == "A"
    assert by[(28.0, "CH_TGA", "wt%")]["grade"] == "D"
    assert by[(28.0, "chem_shrink", "mL/g binder")]["grade"] == "A"
    assert by[(28.0, "DoR_SCM", "%")]["role"] == "validation_only"
    assert s["n_usable_for_likelihood"] == 5
    parts = load_input_file(EXAMPLE)
    assert parts["scm"].name == "GGBS-Pohang-2026" and len(parts["observations"]) == 7


@pytest.mark.parametrize("cmd", ["envelope", "compare", "infer"])
def test_cli_input_paths_mock(tmp_path, bundles_dir, cmd):
    args = [sys.executable, "-m", "dorgems.cli", cmd, "--input", str(EXAMPLE), "--out", str(tmp_path / cmd), "--db", str(tmp_path / "igdb")]
    import os

    env = {**os.environ, "PYTHONIOENCODING": "utf-8"}
    proc = subprocess.run(args, capture_output=True, text=True, encoding="utf-8", errors="replace", timeout=900, cwd=str(ROOT), env=env)
    out = proc.stdout
    payload = json.loads(out[out.index("{"):])
    assert proc.returncode == 0, (proc.stderr[-1500:], out[-800:])
    assert payload.get("ok") or payload.get("summary", {}).get("status") in ("complete",), payload
