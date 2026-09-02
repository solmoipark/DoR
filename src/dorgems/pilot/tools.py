"""GemsPilot-compatible tool layer (spec §10).

Every tool returns the ``inverse-gems-tool/1.0`` ToolResult dict. Query
arguments accept dicts, YAML/JSON strings or file paths. Path keywords are
``out``, ``db``, ``session`` so the GemsPilot runner's workspace remapping
applies unchanged. Policies: ``read`` = pure computation; ``mock_ok`` = mock by
default, ``use_mock=False`` only under the host's ``allow_real``. ``real_gated``
is never used (it means "always refuse" in GemsPilot's ``_policy_check``).

Docstrings are read by the LLM: they say when to use a tool and when not to.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

import yaml

from .. import TOOL_CONTRACT, __version__
from ..kinetics.registry import register as _register_kinetics

try:  # GemsPilot present → reuse its ToolSpec and helpers
    from gemspilot.runner import ToolSpec  # type: ignore
    from gemspilot.agent_tools import _load_query_payload, _log_session  # type: ignore
except Exception:  # noqa: BLE001

    @dataclass
    class ToolSpec:  # type: ignore[no-redef]
        name: str
        func: Callable[..., dict[str, Any]]
        policy: str
        description: str = ""

    def _load_query_payload(query: Any) -> dict[str, Any]:  # type: ignore[no-redef]
        if isinstance(query, dict):
            return query
        text = str(query)
        p = Path(text)
        try:
            if p.is_file():
                text = p.read_text(encoding="utf-8")
        except OSError:
            pass
        data = yaml.safe_load(text)
        if not isinstance(data, dict):
            raise ValueError("Query payload must parse to a mapping.")
        return data

    def _log_session(session: str | None, result: dict[str, Any]) -> dict[str, Any]:  # type: ignore[no-redef]
        return result


_register_kinetics()


def tool_result(tool: str, *, ok: bool, summary: dict[str, Any] | None = None, artifacts: dict[str, str] | None = None, warnings: list[str] | None = None, error: str | None = None) -> dict[str, Any]:
    return {"contract": TOOL_CONTRACT, "tool": tool, "ok": ok, "summary": summary or {}, "artifacts": artifacts or {}, "warnings": warnings or [], "error": error}


def _payload_or_list(obj: Any) -> Any:
    if isinstance(obj, (dict, list)):
        return obj
    text = str(obj)
    p = Path(text)
    try:
        if p.is_file():
            text = p.read_text(encoding="utf-8")
    except OSError:
        pass
    return yaml.safe_load(text)


def _lit_db(db_path: str | None) -> Path | None:
    from ..config import literature_db_path

    return Path(db_path) if db_path else literature_db_path(required=False)


def _audit(tool: str, args: dict[str, Any], ok: bool, xgems_calls: int | None = None, run_dir: str | None = None) -> None:
    try:
        from ..inverse.staging import audit

        audit(tool, args, ok, xgems_calls=xgems_calls, run_dir=run_dir)
    except Exception:  # noqa: BLE001
        pass


# ---------------------------------------------------------------------------
# read tools
# ---------------------------------------------------------------------------


def dor_predict(scm: Any, mix: Any, ages: Any = None, out: str = "out/dor_predict", *, ensemble: str | None = None, method_group: str | None = None, seed: int = 0, lit_db: str | None = None, session: str | None = None) -> dict[str, Any]:
    """Predict the SCM degree of reaction (DoR) curve for a NEW SCM from the literature prior.

    Use when no measured DoR exists for this SCM (or a literature prior is wanted).
    If measured DoR values exist, use GemsPilot ``calibrate_scm_kinetics`` instead.
    Inputs: ``scm`` (SCMSpec: name, role, oxides wt%, optional fineness/amorphous_pct
    measured only), ``mix`` (scm_pct, w_b, curing_temp_C), ``ages`` (days).
    Returns prediction.json with q05/q50/q95 curves (latent and observed), OOD
    flags, and literature analogues as evidence. The LLM must not invent numbers:
    quote them from the artifact.
    """
    from ..predict import predict, write_prediction

    try:
        pred = predict(_payload_or_list(scm), _payload_or_list(mix), _payload_or_list(ages) if ages is not None else None, ensemble=ensemble, method_group=method_group, seed=seed, db_path=lit_db)
        path = write_prediction(pred, out)
    except Exception as exc:  # noqa: BLE001
        return _log_session(session, tool_result("dor_predict", ok=False, error=f"{type(exc).__name__}: {exc}"))
    rec = pred["recommended"]
    summary = {
        "prediction_id": pred["id"],
        "role_bayes": pred["role_bayes"],
        "role_gbm": pred["role_gbm"],
        "ages_d": pred["input"]["ages_d"],
        "alpha_pct_q50": [round(v, 1) for v in rec["alpha_pct_q50"]],
        "alpha_pct_q05": [round(v, 1) for v in rec["alpha_pct_q05"]],
        "alpha_pct_q95": [round(v, 1) for v in rec["alpha_pct_q95"]],
        "source": rec["source"],
        "ensemble": pred["ensemble"].get("mode"),
        "ess": pred["ensemble"].get("ess"),
        "ood_flags": pred["ood"]["flags"],
        "sparse_role": pred["ood"]["sparse_role"],
        "n_analogue_mixes": pred["evidence"].get("n_mixes"),
        "n_analogue_papers": pred["evidence"].get("n_papers"),
    }
    return _log_session(session, tool_result("dor_predict", ok=True, summary=summary, artifacts={"prediction": str(path)}, warnings=pred["warnings"]))


def dor_export_reaction_model(prediction: Any, out: str = "out/reaction_models", *, mode: str = "logistic_fit", slot: str | None = None, quantiles: Any = (0.05, 0.5, 0.95), config_id: str | None = None, signature_files: Any = None) -> dict[str, Any]:
    """Convert a prediction.json into InverseGems ``reaction_model_config`` YAML files, one per quantile.

    mode ``logistic_fit`` (default) is CLI-compatible; ``native`` needs the dorgems
    process; ``pin`` is for alpha grids. Alphas are fractions. Pass the returned
    paths as ``reaction_model_config`` to ``dor_run_forward_with_dor`` (or, for
    design queries, ``design_query.reaction_model.config``).
    """
    from ..kinetics.materials_override import slot_for_role
    from ..kinetics.reaction_model import export_reaction_model

    try:
        pred = _payload_or_list(prediction)
        if slot is None:
            slot, _, _ = slot_for_role(pred["input"]["scm"]["role"], pred["input"]["scm"].get("oxides"))
        qs = tuple(float(q) for q in (_payload_or_list(quantiles) or (0.05, 0.5, 0.95)))
        res = export_reaction_model(pred, out, mode=mode, slot=slot, quantiles=qs, config_id=config_id, signature_files=_payload_or_list(signature_files) if signature_files else None)
    except Exception as exc:  # noqa: BLE001
        return tool_result("dor_export_reaction_model", ok=False, error=f"{type(exc).__name__}: {exc}")
    warnings = [w for r in res.values() for w in r["warnings"]]
    return tool_result("dor_export_reaction_model", ok=True, summary={"mode": mode, "slot": slot, "configs": {k: {"id": v["id"], "params": v["params"], "max_abs_dev_pct": (v["fit"] or {}).get("max_abs_dev_pct")} for k, v in res.items()}}, artifacts={f"reaction_model_{k}": v["path"] for k, v in res.items()}, warnings=warnings)


def dor_build_materials_override(scm: Any, out: str = "out/materials", *, slot: str | None = None, alias: str | None = None, cement: Any = None) -> dict[str, Any]:
    """Write a materials.yaml override that puts the new SCM's oxides into an InverseGems slot.

    Needed because InverseGems fixes its SCM names; the user's name survives only as an alias.
    Pass the returned path as ``materials_config`` to the forward tools.
    """
    from ..kinetics.materials_override import build_materials_config

    try:
        from .schemas import coerce_scm

        res = build_materials_config(coerce_scm(_payload_or_list(scm)), out, slot=slot, alias=alias, cement=_payload_or_list(cement) if cement else None)
    except Exception as exc:  # noqa: BLE001
        return tool_result("dor_build_materials_override", ok=False, error=f"{type(exc).__name__}: {exc}")
    return tool_result("dor_build_materials_override", ok=True, summary={k: v for k, v in res.items() if k != "warnings"}, artifacts={"materials_config": res["path"]}, warnings=res["warnings"])


def dor_find_analogues(scm: Any, mix: Any, *, age_days: float | None = None, quantities: Any = None, limit: int = 20, lit_db: str | None = None) -> dict[str, Any]:
    """List the closest literature mixes (same role, similar chemistry/replacement/w/b) with their DoR observations, DOIs and methods.

    Read-only evidence lookup; use to justify or sanity-check a prediction.
    """
    from ..db.analogues import find_analogues
    from ..db.features import scm_input_to_features
    from ..db.reader import open_ro

    try:
        from .schemas import coerce_mix, coerce_scm

        s, m = coerce_scm(_payload_or_list(scm)), coerce_mix(_payload_or_list(mix))
        dbp = _lit_db(lit_db)
        if not dbp or not dbp.is_file():
            raise FileNotFoundError("literature DB not found (set DORGEMS_DB)")
        con = open_ro(dbp)
        try:
            res = find_analogues(con, scm_input_to_features(s, m), s.role, k=int(limit), age_days=age_days, cache_key=str(dbp))
        finally:
            con.close()
    except Exception as exc:  # noqa: BLE001
        return tool_result("dor_find_analogues", ok=False, error=f"{type(exc).__name__}: {exc}")
    return tool_result("dor_find_analogues", ok=True, summary={"n_mixes": res["n_mixes"], "n_papers": res["n_papers"], "flags": res["flags"], "mixes": res["mixes"]}, warnings=res["flags"])


def dor_compare_reaction_models(config_a: str, config_b: str, *, slot: str | None = None, ages: Any = None) -> dict[str, Any]:
    """Tabulate alpha(t) of two reaction_model_config YAMLs (e.g. DoRGems prior vs calibrate_scm_kinetics fit)."""
    from ..kinetics.reaction_model import compare_reaction_models

    try:
        if slot is None:
            raw = yaml.safe_load(Path(config_a).read_text(encoding="utf-8"))
            slot = next(iter((raw.get("scm_reaction") or {}).keys()))
        res = compare_reaction_models(config_a, config_b, slot, _payload_or_list(ages) if ages else None)
    except Exception as exc:  # noqa: BLE001
        return tool_result("dor_compare_reaction_models", ok=False, error=f"{type(exc).__name__}: {exc}")
    return tool_result("dor_compare_reaction_models", ok=True, summary=res)


def dor_model_card() -> dict[str, Any]:
    """Bundle metadata: model versions, training DB hash, LOPO metrics, known limitations."""
    from ..models.bundle import load_bundle

    try:
        b = load_bundle(require_gbm=False)
    except Exception as exc:  # noqa: BLE001
        return tool_result("dor_model_card", ok=False, error=f"{type(exc).__name__}: {exc}")
    card = {
        "dorgems": __version__,
        "bayes_v4": {k: b.bayes.manifest.get(k) for k in ("model", "training_table", "training_db", "training_db_sha256", "n_rows", "n_papers", "convergence", "posterior_kept_draws", "created_at")},
        "bayes_roles": b.bayes.roles,
        "bayes_methods": b.bayes.methods,
        "gbm_v6": ({k: b.gbm.meta.get(k) for k in ("lopo_r2", "lopo_mae_pct", "n_rows", "n_papers", "sigma_point_pct")} if b.gbm else None),
        "ood_roles": {r: v["n_papers"] for r, v in (b.ood.get("roles") or {}).items()},
        "limitations": [
            "Roles other than slag/fly_ash are pooled as 'other' in the Bayesian model (wide intervals).",
            "Measurement-method bias is not identifiable; only noise scales differ by method.",
            "GBM and Bayes share training data: the blend is a heuristic (see docs/gates.md G1-3).",
            "Clinker hydration is Parrot-Killoh only in InverseGems (cannot be pinned).",
        ],
        "provenance": b.provenance,
    }
    return tool_result("dor_model_card", ok=True, summary=card)


def dor_db_lookup(query_name: str, params: Any = None, *, limit: int = 50, lit_db: str | None = None) -> dict[str, Any]:
    """Run one of the named, parameterised read-only literature-DB queries
    (paper, mix, mixes_for_paper, materials_for_paper, dor_observations,
    observations_for_mix, opc_only_reference). Free SQL is not available."""
    from ..db.reader import LiteratureDB, NAMED_QUERIES, run_named_query

    try:
        dbp = _lit_db(lit_db)
        if not dbp or not dbp.is_file():
            raise FileNotFoundError("literature DB not found (set DORGEMS_DB)")
        with LiteratureDB(dbp) as db:
            rows = run_named_query(db, query_name, _payload_or_list(params) if params else {}, limit=limit)
    except Exception as exc:  # noqa: BLE001
        return tool_result("dor_db_lookup", ok=False, error=f"{type(exc).__name__}: {exc}", summary={"allowed_queries": {k: list(v) for k, v in NAMED_QUERIES.items()}})
    return tool_result("dor_db_lookup", ok=True, summary={"query": query_name, "n": len(rows) if isinstance(rows, list) else (1 if rows else 0), "rows": rows})


# ---------------------------------------------------------------------------
# mock_ok tools (xGEMS)
# ---------------------------------------------------------------------------


def dor_run_forward_with_dor(forward_query: Any, reaction_model_config: str, out: str, db: str, *, materials_config: str | None = None, slot: str | None = None, use_mock: bool = True, dat_lst: str | None = None, max_xgems_calls: int | None = None, capture_species: bool = False, session: str | None = None) -> dict[str, Any]:
    """Run an InverseGems forward query with a DoRGems reaction model (and optional materials override).

    use_mock=False triggers real xGEMS and requires host approval and an explicit
    max_xgems_calls (<= 200). The result is self-checked: the kernel's reaction
    degrees must equal the exported curve.
    """
    from ..gems.forward import run_forward

    try:
        fq = _load_query_payload(forward_query)
        res = run_forward(fq, out=out, db=db, reaction_model_config=reaction_model_config, materials_config=materials_config, slot=slot, use_mock=use_mock, dat_lst=dat_lst, max_xgems_calls=max_xgems_calls, capture_species=capture_species)
    except Exception as exc:  # noqa: BLE001
        _audit("dor_run_forward_with_dor", {"out": out}, False)
        return _log_session(session, tool_result("dor_run_forward_with_dor", ok=False, error=f"{type(exc).__name__}: {exc}"))
    tr = res.to_tool_result()
    _audit("dor_run_forward_with_dor", {"out": out, "use_mock": use_mock}, tr["ok"], run_dir=str(res.run_dir))
    return _log_session(session, tr)


def dor_run_envelope(scm: Any, mix: Any, ages: Any, out: str, db: str, *, use_mock: bool = True, dat_lst: str | None = None, max_xgems_calls: int | None = None, ensemble: str | None = None, seed: int = 0, lit_db: str | None = None, session: str | None = None) -> dict[str, Any]:
    """Scenario A in one call: predict DoR → export q05/q50/q95 reaction models →
    materials override → three forward runs → envelope.csv of phase masses,
    porosity and pH. Mock by default; real xGEMS needs use_mock=False, approval
    and max_xgems_calls.
    """
    from ..envelope import run_envelope

    try:
        res = run_envelope(_payload_or_list(scm), _payload_or_list(mix), _payload_or_list(ages) if ages is not None else None, out=out, db=db, use_mock=use_mock, dat_lst=dat_lst, max_xgems_calls=max_xgems_calls, ensemble=ensemble, seed=seed, lit_db=lit_db)
    except Exception as exc:  # noqa: BLE001
        _audit("dor_run_envelope", {"out": out}, False)
        return _log_session(session, tool_result("dor_run_envelope", ok=False, error=f"{type(exc).__name__}: {exc}"))
    _audit("dor_run_envelope", {"out": out, "use_mock": use_mock}, res["ok"], run_dir=res["out"])
    summary = {"status": "complete" if res["ok"] else "incomplete", "task_type": "dor_envelope", "run_dir": res["out"], "answer_available": res["ok"], "missing_outputs": {}, "result_summary": {"slot": res["slot"], "n_rows_envelope": res["n_rows_envelope"], "runs_ok": {k: v["ok"] for k, v in res["runs"].items()}}}
    return _log_session(session, tool_result("dor_run_envelope", ok=res["ok"], summary=summary, artifacts={"prediction": res["prediction"], "envelope": res["envelope"], "summary": res["summary"], "manifest": res["manifest"], "materials_config": res["materials_config"], **{f"reaction_model_{k}": v for k, v in res["reaction_models"].items()}}, warnings=res["warnings"]))


def _lazy(name: str) -> Callable[..., dict[str, Any]]:
    def _f(*a: Any, **k: Any) -> dict[str, Any]:
        from . import tools_b_c

        return getattr(tools_b_c, name)(*a, **k)

    _f.__name__ = name
    _f.__doc__ = f"See dorgems.pilot.tools_b_c.{name}"
    try:
        from . import tools_b_c as _m

        _f.__doc__ = getattr(_m, name).__doc__
    except Exception:  # noqa: BLE001
        pass
    return _f


TOOLSET: list[ToolSpec] = [
    ToolSpec("dor_predict", dor_predict, "read"),
    ToolSpec("dor_export_reaction_model", dor_export_reaction_model, "read"),
    ToolSpec("dor_build_materials_override", dor_build_materials_override, "read"),
    ToolSpec("dor_find_analogues", dor_find_analogues, "read"),
    ToolSpec("dor_compare_reaction_models", dor_compare_reaction_models, "read"),
    ToolSpec("dor_model_card", dor_model_card, "read"),
    ToolSpec("dor_db_lookup", dor_db_lookup, "read"),
    ToolSpec("dor_run_forward_with_dor", dor_run_forward_with_dor, "mock_ok"),
    ToolSpec("dor_run_envelope", dor_run_envelope, "mock_ok"),
    ToolSpec("dor_compare_to_literature", _lazy("dor_compare_to_literature"), "mock_ok"),
    ToolSpec("dor_opc_reference_check", _lazy("dor_opc_reference_check"), "mock_ok"),
    ToolSpec("dor_infer_from_observations", _lazy("dor_infer_from_observations"), "mock_ok"),
    ToolSpec("dor_stage_inferred", _lazy("dor_stage_inferred"), "mock_ok"),
]

TOOLS_BY_NAME: dict[str, ToolSpec] = {t.name: t for t in TOOLSET}


def tool_manifest() -> list[dict[str, Any]]:
    return [{"name": t.name, "policy": t.policy, "doc": (t.func.__doc__ or "").strip().split("\n")[0]} for t in TOOLSET]
