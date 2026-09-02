"""Scenario B (validation) and C (inverse analysis) tools; staging (spec §10, §11)."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from .tools import _audit, _lit_db, _log_session, _payload_or_list, tool_result


def dor_compare_to_literature(out: str, db: str, *, run_dir: str | None = None, mix_uid: str | None = None, scm: Any = None, mix: Any = None, mode: str = "twin", quantities: Any = None, use_mock: bool = True, dat_lst: str | None = None, max_xgems_calls: int | None = None, lit_db: str | None = None, session: str | None = None, max_mixes: int | None = None) -> dict[str, Any]:
    """Validate the thermodynamic model against literature observations (CH, bound water, chemical shrinkage).

    mode="twin" (default, exact): give ``mix_uid`` (a literature mix) — it is rebuilt as an
    InverseGems recipe, run with its own measured DoR pinned (or the model q50) and
    compared 1:1 with the same mix's observations. Without mix_uid, twin mode runs
    over all literature mixes that have usable observations and a pinned DoR (cap with
    ``max_mixes``). mode="neighbourhood": compare an existing ``run_dir`` with
    analogous literature mixes (looser, one grade lower). Verdicts are deterministic
    thresholds; numbers live in comparison.csv. Mock runs only exercise the pipeline.
    """
    from ..db.reader import LiteratureDB
    from ..validate.twin import twin_compare_mix

    try:
        dbp = _lit_db(lit_db)
        if not dbp or not dbp.is_file():
            raise FileNotFoundError("literature DB not found (set DORGEMS_DB)")
        qs = tuple(_payload_or_list(quantities)) if quantities else ("CH_TGA", "CH_XRD", "bound_water", "chem_shrink")
        if mode == "neighbourhood":
            from ..validate.neighbourhood import neighbourhood_compare

            res = neighbourhood_compare(dbp, run_dir=run_dir, scm=_payload_or_list(scm), mix=_payload_or_list(mix), out=Path(out), quantities=qs)
        elif mix_uid:
            with LiteratureDB(dbp) as db_:
                res = twin_compare_mix(db_, mix_uid, out=Path(out), ig_db=db, use_mock=use_mock, dat_lst=dat_lst, max_xgems_calls=max_xgems_calls, quantities=qs)
        else:
            from ..validate.twin_batch import twin_batch

            res = twin_batch(dbp, out=Path(out), ig_db=db, use_mock=use_mock, dat_lst=dat_lst, max_xgems_calls=max_xgems_calls, quantities=qs, max_mixes=max_mixes)
    except Exception as exc:  # noqa: BLE001
        _audit("dor_compare_to_literature", {"out": out, "mode": mode}, False)
        return _log_session(session, tool_result("dor_compare_to_literature", ok=False, error=f"{type(exc).__name__}: {exc}"))
    _audit("dor_compare_to_literature", {"out": out, "mode": mode, "use_mock": use_mock}, bool(res.get("ok")), run_dir=out)
    return _log_session(session, tool_result("dor_compare_to_literature", ok=bool(res.get("ok")), summary={k: v for k, v in res.items() if k not in ("files", "warnings")}, artifacts=res.get("files") or {}, warnings=res.get("warnings") or [], error=res.get("error")))


def dor_opc_reference_check(out: str, db: str, *, age_days: float = 28, w_b_range: Any = (0.4, 0.5), use_mock: bool = True, dat_lst: str | None = None, max_xgems_calls: int | None = None, lit_db: str | None = None, max_mixes: int | None = None) -> dict[str, Any]:
    """Kernel baseline check on OPC-only literature pastes (spec §8.6): portlandite
    residuals at 28 d for grade-A CH_TGA observations. Independent of the DoR model.
    Gate G2-3: |median r| ≤ 4 g/100 g over ≥ 30 papers. In mock mode the numbers are
    not physical — the tool only proves the pipeline runs.
    """
    from ..db.reader import LiteratureDB
    from ..validate.twin import opc_reference_check

    try:
        dbp = _lit_db(lit_db)
        if not dbp or not dbp.is_file():
            raise FileNotFoundError("literature DB not found (set DORGEMS_DB)")
        wb = tuple(float(x) for x in _payload_or_list(w_b_range))
        with LiteratureDB(dbp) as db_:
            res = opc_reference_check(db_, out=Path(out), ig_db=db, age_days=age_days, w_b_range=wb, use_mock=use_mock, dat_lst=dat_lst, max_xgems_calls=max_xgems_calls, max_mixes=max_mixes)
    except Exception as exc:  # noqa: BLE001
        _audit("dor_opc_reference_check", {"out": out}, False)
        return tool_result("dor_opc_reference_check", ok=False, error=f"{type(exc).__name__}: {exc}")
    _audit("dor_opc_reference_check", {"out": out, "use_mock": use_mock}, True, run_dir=out)
    return tool_result("dor_opc_reference_check", ok=True, summary={k: v for k, v in res.items() if k not in ("files", "warnings")}, artifacts=res["files"], warnings=res["warnings"])


def dor_infer_from_observations(mix: Any, observations: Any, out: str, db: str, *, scm: Any = None, mix_uid: str | None = None, prior: str = "model", alpha_grid: int = 21, use_mock: bool = True, dat_lst: str | None = None, max_xgems_calls: int | None = None, lit_db: str | None = None, seed: int = 0, session: str | None = None) -> dict[str, Any]:
    """Infer the SCM DoR curve from indirect observations (CH_TGA/CH_XRD, bound water,
    chemical shrinkage) with the DoR model as prior (or prior="flat").

    Do NOT use when measured DoR values exist — use GemsPilot ``calibrate_scm_kinetics``.
    Direct DoR observations in the input are kept for validation only. Runs an alpha-grid
    of pinned forward calculations (≤ 21 × n_ages xGEMS calls; real runs need approval
    and max_xgems_calls). Outputs inference.json, reaction model YAMLs and inferred_dor.csv.
    """
    from ..inverse.run import infer_from_observations

    try:
        res = infer_from_observations(_payload_or_list(mix) if mix else {}, _payload_or_list(observations) or [], out=out, ig_db=db, scm=_payload_or_list(scm) if scm else None, mix_uid=mix_uid, lit_db=lit_db or (_lit_db(None) if mix_uid else None), prior=prior, alpha_grid_n=alpha_grid, use_mock=use_mock, dat_lst=dat_lst, max_xgems_calls=max_xgems_calls, seed=seed)
    except Exception as exc:  # noqa: BLE001
        _audit("dor_infer_from_observations", {"out": out}, False)
        return _log_session(session, tool_result("dor_infer_from_observations", ok=False, error=f"{type(exc).__name__}: {exc}"))
    _audit("dor_infer_from_observations", {"out": out, "use_mock": use_mock}, True, xgems_calls=res["summary"].get("xgems_calls"), run_dir=out)
    return _log_session(session, tool_result("dor_infer_from_observations", ok=True, summary={"inference_id": res["inference_id"], "slot": res["slot"], **res["summary"], "validation": res.get("validation")}, artifacts=res["files"], warnings=res["warnings"]))


def dor_stage_inferred(inference: Any, staging_db: str | None = None, *, use_mock: bool = True, note: str | None = None) -> dict[str, Any]:
    """Write an inference into the staging DB table ``inferred_dor`` with reviewed=0.

    use_mock=True (default) is a dry-run that only returns the rows that would be
    written; use_mock=False writes (host approval). The literature DB is never written.
    """
    from ..inverse.staging import stage_inference

    try:
        inf = _payload_or_list(inference)
        if isinstance(inf, str):
            inf = json.loads(Path(inf).read_text(encoding="utf-8"))
        res = stage_inference(inf, path=staging_db, dry_run=use_mock, note=note)
    except Exception as exc:  # noqa: BLE001
        _audit("dor_stage_inferred", {"dry_run": use_mock}, False)
        return tool_result("dor_stage_inferred", ok=False, error=f"{type(exc).__name__}: {exc}")
    _audit("dor_stage_inferred", {"dry_run": use_mock}, True)
    return tool_result("dor_stage_inferred", ok=True, summary=res, warnings=["dry-run: nothing written (use_mock=True)"] if use_mock else [])
