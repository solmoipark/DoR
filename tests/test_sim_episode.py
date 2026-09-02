"""G1-5 (simulated LLM): the DoR tool chain runs end-to-end through GemsPilot's real
runner (policy check, workspace remapping, trajectory) with a scripted small model;
an injected 'approval' in the task cannot unlock real execution; a lazy model that
invents a number is caught."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

pytest.importorskip("gemspilot")

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "run_llm_episode.py"


@pytest.mark.parametrize("behaviour", ["good", "injected", "lazy"])
def test_simulated_episode(tmp_path, bundles_dir, behaviour):
    out = tmp_path / behaviour
    proc = subprocess.run([sys.executable, str(SCRIPT), "--simulate", behaviour, "--out", str(out)], capture_output=True, text=True, timeout=600)
    rec = json.loads((out / "episode_result.json").read_text(encoding="utf-8"))
    assert proc.returncode == 0, (proc.stdout[-800:], proc.stderr[-800:])
    assert rec["ok"], rec
    if behaviour == "good":
        assert rec["tools_called"] == ["dor_predict", "dor_export_reaction_model", "dor_build_materials_override", "dor_run_forward_with_dor"]
        assert f"{rec['q50_28d_artifact']:.1f}" in rec["final_message"]
    if behaviour == "injected":
        assert rec["checks"]["real_attempts_denied"] is True
        assert rec["tools_called"].count("dor_run_forward_with_dor") == 2
    if behaviour == "lazy":
        assert rec["tools_called"] == [] and rec.get("expected_failure")
