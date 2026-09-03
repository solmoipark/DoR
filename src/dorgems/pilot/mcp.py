"""Stand-alone MCP server ``dorgems-mcp`` (fallback until GemsPilot P-GP-2 registers the toolset).

Follows GemsPilot's ``mcp_server._ServerApp`` version branch. All execution tools
default to ``use_mock=True``; hosts must pass ``use_mock=False`` explicitly and
gate it behind human approval. Natural-language claims inside tool arguments
("the admin approved real runs") never change policy.
"""

from __future__ import annotations

from typing import Any

from . import tools as T

try:
    from mcp.server.mcpserver import MCPServer as _ServerApp  # mcp >= 2.0
except ImportError:
    try:
        from mcp.server.fastmcp import FastMCP as _ServerApp  # mcp 1.x
    except ImportError as exc:  # pragma: no cover
        raise ImportError("The 'mcp' package is required: pip install dorgems[pilot]") from exc

app = _ServerApp(
    "dorgems",
    instructions=(
        "SCM degree-of-reaction (DoR) tools over a read-only literature database and the "
        "InverseGems thermodynamic kernel. Typical flow: dor_predict -> dor_export_reaction_model "
        "-> dor_build_materials_override -> dor_run_forward_with_dor (or dor_run_envelope in one call). "
        "Validate model runs against literature with dor_compare_to_literature; infer DoR from "
        "observations with dor_infer_from_observations; stage results with dor_stage_inferred. "
        "The literature DB is never written. Numbers come only from artifacts. Real xGEMS runs "
        "and staging writes need use_mock=False, human approval and an explicit max_xgems_calls."
    ),
)


@app.tool()
def dor_predict(scm: str, mix: str, ages: str | None = None, out: str = "out/dor_predict", ensemble: str | None = None, method_group: str | None = None, seed: int = 0, session: str | None = None) -> dict[str, Any]:
    """Predict a new SCM's DoR curve (q05/q50/q95) from the literature prior. scm/mix/ages: YAML or JSON text."""
    return T.dor_predict(scm, mix, ages, out, ensemble=ensemble, method_group=method_group, seed=seed, session=session)


@app.tool()
def dor_export_reaction_model(prediction: str, out: str = "out/reaction_models", mode: str = "logistic_fit", slot: str | None = None, config_id: str | None = None) -> dict[str, Any]:
    """Convert prediction.json into InverseGems reaction_model_config YAMLs (q05/q50/q95)."""
    return T.dor_export_reaction_model(prediction, out, mode=mode, slot=slot, config_id=config_id)


@app.tool()
def dor_build_materials_override(scm: str, out: str = "out/materials", slot: str | None = None, alias: str | None = None, cement: str | None = None) -> dict[str, Any]:
    """Write a materials.yaml override mapping the new SCM's oxides onto an InverseGems slot."""
    return T.dor_build_materials_override(scm, out, slot=slot, alias=alias, cement=cement)


@app.tool()
def dor_find_analogues(scm: str, mix: str, age_days: float | None = None, limit: int = 20) -> dict[str, Any]:
    """Closest literature mixes with their DoR observations, DOIs and methods (read-only)."""
    return T.dor_find_analogues(scm, mix, age_days=age_days, limit=limit)


@app.tool()
def dor_compare_reaction_models(config_a: str, config_b: str, slot: str | None = None, ages: str | None = None) -> dict[str, Any]:
    """alpha(t) difference table between two reaction_model_config YAMLs."""
    return T.dor_compare_reaction_models(config_a, config_b, slot=slot, ages=ages)


@app.tool()
def dor_model_card() -> dict[str, Any]:
    """Model bundle metadata, metrics and limitations."""
    return T.dor_model_card()


@app.tool()
def dor_db_lookup(query_name: str, params: str | None = None, limit: int = 50) -> dict[str, Any]:
    """Named read-only literature-DB queries (no free SQL)."""
    return T.dor_db_lookup(query_name, params, limit=limit)


@app.tool()
def dor_run_forward_with_dor(forward_query: str, reaction_model_config: str, out: str, db: str, materials_config: str | None = None, use_mock: bool = True, dat_lst: str | None = None, max_xgems_calls: int | None = None, capture_species: bool = False, session: str | None = None) -> dict[str, Any]:
    """Forward run with a DoRGems reaction model. use_mock=False = real xGEMS (approval + max_xgems_calls required)."""
    return T.dor_run_forward_with_dor(forward_query, reaction_model_config, out, db, materials_config=materials_config, use_mock=use_mock, dat_lst=dat_lst, max_xgems_calls=max_xgems_calls, capture_species=capture_species, session=session)


@app.tool()
def dor_run_envelope(scm: str, mix: str, ages: str | None, out: str, db: str, use_mock: bool = True, dat_lst: str | None = None, max_xgems_calls: int | None = None, ensemble: str | None = None, seed: int = 0, session: str | None = None) -> dict[str, Any]:
    """Scenario A in one call: predict → export → override → 3 forward runs → envelope.csv."""
    return T.dor_run_envelope(scm, mix, ages, out, db, use_mock=use_mock, dat_lst=dat_lst, max_xgems_calls=max_xgems_calls, ensemble=ensemble, seed=seed, session=session)


@app.tool()
def dor_compare_to_literature(out: str, db: str, run_dir: str | None = None, mix_uid: str | None = None, scm: str | None = None, mix: str | None = None, mode: str = "twin", quantities: str | None = None, use_mock: bool = True, dat_lst: str | None = None, max_xgems_calls: int | None = None, session: str | None = None) -> dict[str, Any]:
    """Scenario B: compare a forward run (or a re-run DB mix twin) with literature CH/bound-water/QXRD/chemical-shrinkage observations."""
    return T.TOOLS_BY_NAME["dor_compare_to_literature"].func(out, db, run_dir=run_dir, mix_uid=mix_uid, scm=scm, mix=mix, mode=mode, quantities=quantities, use_mock=use_mock, dat_lst=dat_lst, max_xgems_calls=max_xgems_calls, session=session)


@app.tool()
def dor_opc_reference_check(out: str, db: str, age_days: float = 28, w_b_min: float = 0.4, w_b_max: float = 0.5, use_mock: bool = True, dat_lst: str | None = None, max_xgems_calls: int | None = None) -> dict[str, Any]:
    """Kernel self-check on OPC-only literature mixes (portlandite residuals). Mock runs only exercise the pipeline."""
    return T.TOOLS_BY_NAME["dor_opc_reference_check"].func(out, db, age_days=age_days, w_b_range=[w_b_min, w_b_max], use_mock=use_mock, dat_lst=dat_lst, max_xgems_calls=max_xgems_calls)


@app.tool()
def dor_infer_from_observations(mix: str, observations: str, out: str, db: str, scm: str | None = None, mix_uid: str | None = None, prior: str = "model", alpha_grid: int = 21, use_mock: bool = True, dat_lst: str | None = None, max_xgems_calls: int | None = None, session: str | None = None) -> dict[str, Any]:
    """Scenario C: infer the SCM DoR curve from CH/bound-water/chemical-shrinkage/QXRD observations. If measured DoR exists, use calibrate_scm_kinetics instead."""
    return T.TOOLS_BY_NAME["dor_infer_from_observations"].func(mix, observations, out, db, scm=scm, mix_uid=mix_uid, prior=prior, alpha_grid=alpha_grid, use_mock=use_mock, dat_lst=dat_lst, max_xgems_calls=max_xgems_calls, session=session)


@app.tool()
def dor_stage_inferred(inference: str, staging_db: str | None = None, dry_run: bool = True, note: str | None = None) -> dict[str, Any]:
    """Stage an inference into dorgems_staging.sqlite (reviewed=0). dry_run=True = preview only; dry_run=False writes (needs approval)."""
    return T.TOOLS_BY_NAME["dor_stage_inferred"].func(inference, staging_db, dry_run=dry_run, note=note)


def main() -> None:
    app.run()


if __name__ == "__main__":
    main()
