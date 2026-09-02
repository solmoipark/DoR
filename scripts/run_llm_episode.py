"""G1-5: one GemsPilot episode driven by an LLM through the DoR toolset (spec §12).

The agent (a small/mid-tier model over OpenRouter by default) must complete
    dor_predict -> dor_export_reaction_model -> dor_run_forward_with_dor (mock)
for a new slag and report the 28-day DoR median from the artifacts. The model never
computes numbers itself; grading checks that the tool chain ran and that the reported
number equals the artifact value.

Requires: `pip install litellm` and OPENROUTER_API_KEY (or OPENAI_API_KEY with an
`openai/...` model id). Without a key the script exits with code 3 and writes a
"needs_api_key" record so the gate log can cite it.

Usage: python scripts/run_llm_episode.py [--model openrouter/anthropic/claude-haiku-4.5]
       [--out out/llm_episode] [--allow-real]
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

DEFAULT_MODEL = "openrouter/anthropic/claude-haiku-4.5"  # 'small' tier in GemsPilot configs/models.yaml

TASK = (
    "A new ground granulated blast-furnace slag 'GGBS-Pohang-2026' has the oxide composition "
    "CaO 42.5, SiO2 34.0, Al2O3 13.5, MgO 6.5, Fe2O3 0.5, SO3 2.0 wt% and Blaine 450 m2/kg. "
    "It will be used at 40 % replacement, w/b 0.45, 20 C. Predict its degree of reaction with "
    "dor_predict (ages 1, 7, 28, 90 days), export the q50 reaction model with "
    "dor_export_reaction_model, build the materials override with dor_build_materials_override, "
    "then run a mock forward calculation at 28 days with dor_run_forward_with_dor using that "
    "reaction model and materials override (out and db under the workspace). Report the 28-day "
    "DoR median in percent exactly as written in prediction.json and state whether the forward "
    "run's self-check passed."
)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default=os.environ.get("DORGEMS_LLM_MODEL", DEFAULT_MODEL))
    ap.add_argument("--out", default="out/llm_episode")
    ap.add_argument("--allow-real", action="store_true")
    ap.add_argument("--max-steps", type=int, default=14)
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    key_env = "OPENAI_API_KEY" if a.model.startswith("openai/") else "OPENROUTER_API_KEY"
    if not os.environ.get(key_env):
        rec = {"status": "needs_api_key", "model": a.model, "key_env": key_env, "task": TASK}
        (out / "episode_result.json").write_text(json.dumps(rec, indent=2), encoding="utf-8")
        print(f"[G1-5] {key_env} not set; episode not run (record written to {out / 'episode_result.json'})")
        return 3
    try:
        import litellm  # noqa: F401
    except ImportError:
        print("pip install litellm")
        return 2
    from gemspilot.runner import Episode, ToolSpec, default_toolset, run_episode

    from dorgems.pilot.tools import TOOLSET

    base = default_toolset()
    names = {t.name for t in base}
    toolset = base + [ToolSpec(t.name, t.func, t.policy) for t in TOOLSET if t.name not in names]
    ep = Episode(model=a.model, workspace=out / "workspace", allow_real=a.allow_real, max_steps=a.max_steps, toolset=toolset)
    result = run_episode(TASK, ep)
    # grade: chain executed, mock only, number matches artifact
    calls = [s.get("tool") for s in result.get("steps", []) if isinstance(s, dict)]
    pred = None
    for p in (out / "workspace").rglob("prediction.json"):
        pred = json.loads(p.read_text(encoding="utf-8"))
        break
    q50_28 = None
    if pred:
        ages = pred["input"]["ages_d"]
        if 28.0 in ages:
            q50_28 = pred["recommended"]["alpha_pct_q50"][ages.index(28.0)]
    final = str(result.get("final_message") or result.get("answer") or "")
    checks = {
        "chain_called": all(t in calls for t in ("dor_predict", "dor_export_reaction_model", "dor_run_forward_with_dor")),
        "no_real_execution": not any((s.get("arguments") or {}).get("use_mock") is False for s in result.get("steps", []) if isinstance(s, dict)),
        "number_in_answer": (q50_28 is not None and f"{q50_28:.1f}" in final) if q50_28 is not None else None,
    }
    rec = {"status": "ran", "model": a.model, "tools_called": calls, "q50_28d_artifact": q50_28, "final_message": final[:2000], "checks": checks, "ok": all(v for v in checks.values() if v is not None)}
    (out / "episode_result.json").write_text(json.dumps(rec, indent=2, default=str), encoding="utf-8")
    print(json.dumps({k: rec[k] for k in ("status", "model", "tools_called", "checks", "ok")}, indent=2))
    return 0 if rec["ok"] else 1


if __name__ == "__main__":
    sys.exit(main())
