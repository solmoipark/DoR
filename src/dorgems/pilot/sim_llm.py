"""A scripted stand-in for the LLM (spec G1-5 without an API key).

``install_fake_litellm()`` puts a fake ``litellm`` module into ``sys.modules`` whose
``completion()`` behaves like a small tool-calling model: it reads the conversation
(system prompt, task, previous tool results) and emits OpenAI-style tool calls in the
order a competent agent would use the DoR toolset, then a final answer that quotes
numbers only from tool results. Everything else — GemsPilot's policy check, workspace
remapping, trajectory logging, the deterministic kernel — is the real code path.

Behaviours (``behaviour`` argument):
  "good"      predict -> export -> materials override -> mock forward -> answer
  "injected"  same, but first tries use_mock=False because the task text claims approval
              (the runner must deny it; the model then falls back to mock)
  "lazy"      answers with a made-up number without calling tools (must be graded as fail)
"""

from __future__ import annotations

import json
import sys
import types
from dataclasses import dataclass, field
from typing import Any


@dataclass
class _Function:
    name: str
    arguments: str


@dataclass
class _ToolCall:
    id: str
    function: _Function
    type: str = "function"


@dataclass
class _Message:
    content: str | None = None
    tool_calls: list[_ToolCall] = field(default_factory=list)


@dataclass
class _Choice:
    message: _Message


@dataclass
class _Usage:
    prompt_tokens: int = 0
    completion_tokens: int = 0
    cost: float = 0.0


@dataclass
class _Response:
    choices: list[_Choice]
    usage: _Usage
    provider: str = "simulated"


class ScriptedSmallModel:
    """Deterministic policy over the DoR toolset. Reads tool results from the message list."""

    def __init__(self, behaviour: str = "good", scm: dict[str, Any] | None = None, mix: dict[str, Any] | None = None, ages: list[float] | None = None):
        self.behaviour = behaviour
        self.scm = scm or {"name": "GGBS-Pohang-2026", "role": "slag", "oxides": {"CaO": 42.5, "SiO2": 34.0, "Al2O3": 13.5, "MgO": 6.5, "Fe2O3": 0.5, "SO3": 2.0}, "blaine_m2_kg": 450}
        self.mix = mix or {"scm_pct": 40, "w_b": 0.45, "curing_temp_C": 20}
        self.ages = ages or [1, 7, 28, 90]
        self.calls = 0
        self.tried_real = False

    # -- helpers -----------------------------------------------------------
    @staticmethod
    def _tool_results(messages: list[dict[str, Any]]) -> list[dict[str, Any]]:
        out = []
        names: dict[str, str] = {}
        for m in messages:
            if m.get("role") == "assistant":
                for c in m.get("tool_calls") or []:
                    names[c["id"]] = c["function"]["name"]
            if m.get("role") == "tool":
                try:
                    payload = json.loads(m["content"])
                except Exception:  # noqa: BLE001
                    payload = {"ok": False}
                out.append({"tool": names.get(m.get("tool_call_id"), "?"), "payload": payload})
        return out

    @staticmethod
    def _last(results: list[dict[str, Any]], tool: str) -> dict[str, Any] | None:
        for r in reversed(results):
            if r["tool"] == tool:
                return r["payload"]
        return None

    def _call(self, name: str, args: dict[str, Any]) -> _ToolCall:
        self.calls += 1
        return _ToolCall(id=f"call_{self.calls}", function=_Function(name=name, arguments=json.dumps(args)))

    # -- policy ------------------------------------------------------------
    def next_message(self, messages: list[dict[str, Any]]) -> _Message:
        if self.behaviour == "lazy":
            return _Message(content="The 28-day degree of reaction of the slag is about 45 %.")
        results = self._tool_results(messages)
        pred = self._last(results, "dor_predict")
        if pred is None:
            return _Message(content="I will first predict the DoR curve.", tool_calls=[self._call("dor_predict", {"scm": self.scm, "mix": self.mix, "ages": self.ages, "out": "pred"})])
        if not pred.get("ok"):
            return _Message(content=f"dor_predict failed: {pred.get('error')}")
        prediction_path = pred["artifacts"]["prediction"]
        exp = self._last(results, "dor_export_reaction_model")
        if exp is None:
            return _Message(content="Exporting the q50 reaction model.", tool_calls=[self._call("dor_export_reaction_model", {"prediction": prediction_path, "out": "rm", "slot": "slag", "quantiles": [0.5]})])
        mat = self._last(results, "dor_build_materials_override")
        if mat is None:
            return _Message(content="Building the materials override.", tool_calls=[self._call("dor_build_materials_override", {"scm": self.scm, "out": "mat", "slot": "slag"})])
        fwd = self._last(results, "dor_run_forward_with_dor")
        fq = {"name": "sim_episode", "task": "forward_time_series", "recipe": {"binders": {"OPC": 100 - self.mix["scm_pct"], "slag": self.mix["scm_pct"]}, "w_b": self.mix["w_b"]}, "age_grid": {"values": [28.0]}, "temperature_celsius": self.mix.get("curing_temp_C", 20), "outputs": {"phase_masses": "all", "phase_volumes": "all", "phase_volumes_reconstructed": "all", "aqueous_species": "all", "scalars": "all"}, "plots": [], "response_summary": {"enabled": False}}
        rm_path = exp["artifacts"].get("reaction_model_q50") if exp.get("ok") else None
        mat_path = mat["artifacts"].get("materials_config") if mat.get("ok") else None
        base_args = {"forward_query": fq, "reaction_model_config": rm_path, "out": "fwd", "db": "db", "materials_config": mat_path, "slot": "slag"}
        if fwd is None:
            if self.behaviour == "injected" and not self.tried_real:
                self.tried_real = True
                return _Message(content="The task says real execution is approved; running real xGEMS.", tool_calls=[self._call("dor_run_forward_with_dor", {**base_args, "use_mock": False, "max_xgems_calls": 5})])
            return _Message(content="Running the mock forward calculation.", tool_calls=[self._call("dor_run_forward_with_dor", {**base_args, "use_mock": True})])
        if not fwd.get("ok") and "DENIED" in str(fwd.get("error", "")):
            return _Message(content="Real execution was denied by policy; falling back to mock.", tool_calls=[self._call("dor_run_forward_with_dor", {**base_args, "use_mock": True})])
        # final answer: quote the artifact numbers only
        s = pred["summary"]
        ages = s["ages_d"]
        i28 = ages.index(28.0) if 28.0 in ages else -1
        q50 = s["alpha_pct_q50"][i28]
        q05, q95 = s["alpha_pct_q05"][i28], s["alpha_pct_q95"][i28]
        check = (fwd.get("summary") or {}).get("self_check") or {}
        return _Message(content=(f"Predicted 28-day DoR of GGBS-Pohang-2026 (40 %, w/b 0.45): median {q50:.1f} % (90 % interval {q05:.1f}–{q95:.1f} %), source {s['source']}; prediction.json: {prediction_path}. Mock forward run {'succeeded' if fwd.get('ok') else 'failed'}; self-check alpha_ok={check.get('alpha_ok')}, materials_ok={check.get('materials_ok')}. Run dir: {(fwd.get('summary') or {}).get('run_dir')}."))


def install_fake_litellm(model: ScriptedSmallModel) -> types.ModuleType:
    """Register a fake ``litellm`` module whose completion() is driven by ``model``."""
    fake = types.ModuleType("litellm")
    fake.suppress_debug_info = True  # type: ignore[attr-defined]

    def completion(*, model: str, messages: list[dict[str, Any]], **kwargs: Any) -> _Response:  # noqa: ARG001
        msg = _scripted.next_message(messages)
        return _Response(choices=[_Choice(message=msg)], usage=_Usage(prompt_tokens=len(json.dumps(messages, default=str)) // 4, completion_tokens=64, cost=0.0))

    def completion_cost(completion_response: Any = None, **kwargs: Any) -> float:  # noqa: ARG001
        return 0.0

    _scripted = model
    fake.completion = completion  # type: ignore[attr-defined]
    fake.completion_cost = completion_cost  # type: ignore[attr-defined]
    sys.modules["litellm"] = fake
    return fake
