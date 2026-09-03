"""``dorgems`` command line (spec §3.2). Every subcommand is a thin wrapper over
the kernel so the package is fully usable without any agent host."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import yaml


def _load(obj: str | None) -> Any:
    if obj is None:
        return None
    p = Path(obj)
    text = p.read_text(encoding="utf-8") if p.is_file() else obj
    return yaml.safe_load(text)


def _print(obj: Any) -> None:
    print(json.dumps(obj, indent=2, ensure_ascii=False, default=str))


def _from_input(a: argparse.Namespace) -> dict[str, Any] | None:
    """--input <template.yaml> fills scm / mix / ages / observations when given."""
    if not getattr(a, "input", None):
        return None
    from .pilot.input_file import load_input_file

    parts = load_input_file(a.input)
    return {"scm": parts["scm"].model_dump(), "mix": parts["mix"].model_dump(), "ages": parts["ages_d"], "observations": [o.model_dump() for o in parts["observations"]]}


def cmd_validate(a: argparse.Namespace) -> int:
    from .pilot.input_file import validate_report

    r = validate_report(a.input)
    _print(r)
    return 0 if r["ok"] else 1


def cmd_predict(a: argparse.Namespace) -> int:
    from .pilot.tools import dor_predict

    _print(dor_predict(_load(a.scm), _load(a.mix), _load(a.ages), a.out, ensemble=a.ensemble, method_group=a.method_group, seed=a.seed, lit_db=a.lit_db))
    return 0


def cmd_export(a: argparse.Namespace) -> int:
    from .pilot.tools import dor_export_reaction_model

    _print(dor_export_reaction_model(a.prediction, a.out, mode=a.mode, slot=a.slot, config_id=a.config_id))
    return 0


def cmd_materials(a: argparse.Namespace) -> int:
    from .pilot.tools import dor_build_materials_override

    _print(dor_build_materials_override(_load(a.scm), a.out, slot=a.slot, alias=a.alias, cement=_load(a.cement)))
    return 0


def cmd_analogues(a: argparse.Namespace) -> int:
    from .pilot.tools import dor_find_analogues

    _print(dor_find_analogues(_load(a.scm), _load(a.mix), age_days=a.age_days, limit=a.limit, lit_db=a.lit_db))
    return 0


def cmd_run_forward(a: argparse.Namespace) -> int:
    from .pilot.tools import dor_run_forward_with_dor

    r = dor_run_forward_with_dor(a.forward_query, a.reaction_model_config, a.out, a.db, materials_config=a.materials_config, slot=a.slot, use_mock=not a.real, dat_lst=a.dat_lst, max_xgems_calls=a.max_xgems_calls, capture_species=a.capture_species)
    _print(r)
    return 0 if r["ok"] else 1


def cmd_envelope(a: argparse.Namespace) -> int:
    from .pilot.tools import dor_run_envelope

    inp = _from_input(a)
    scm, mix, ages = (inp["scm"], inp["mix"], inp["ages"]) if inp else (_load(a.scm), _load(a.mix), _load(a.ages))
    r = dor_run_envelope(scm, mix, ages, a.out, a.db, use_mock=not a.real, dat_lst=a.dat_lst, max_xgems_calls=a.max_xgems_calls, ensemble=a.ensemble, seed=a.seed, lit_db=a.lit_db)
    _print(r)
    return 0 if r["ok"] else 1


def cmd_compare_models(a: argparse.Namespace) -> int:
    from .pilot.tools import dor_compare_reaction_models

    _print(dor_compare_reaction_models(a.config_a, a.config_b, slot=a.slot, ages=_load(a.ages)))
    return 0


def cmd_model_card(a: argparse.Namespace) -> int:
    from .pilot.tools import dor_model_card

    _print(dor_model_card())
    return 0


def cmd_db_lookup(a: argparse.Namespace) -> int:
    from .pilot.tools import dor_db_lookup

    _print(dor_db_lookup(a.query_name, _load(a.params), limit=a.limit, lit_db=a.lit_db))
    return 0


def cmd_compare(a: argparse.Namespace) -> int:
    from .pilot.tools_b_c import dor_compare_to_literature

    inp = _from_input(a)
    if inp:
        from .validate.user_twin import user_compare

        r = user_compare(inp["scm"], inp["mix"], inp["observations"], out=a.out, ig_db=a.db, use_mock=not a.real, dat_lst=a.dat_lst, max_xgems_calls=a.max_xgems_calls, lit_db=a.lit_db)
        _print(r)
        return 0 if r.get("ok") else 1
    r = dor_compare_to_literature(a.out, a.db, run_dir=a.run_dir, mix_uid=a.mix_uid, scm=_load(a.scm), mix=_load(a.mix), mode=a.mode, quantities=_load(a.quantities), use_mock=not a.real, dat_lst=a.dat_lst, max_xgems_calls=a.max_xgems_calls, lit_db=a.lit_db)
    _print(r)
    return 0 if r["ok"] else 1


def cmd_opc_check(a: argparse.Namespace) -> int:
    from .pilot.tools_b_c import dor_opc_reference_check

    r = dor_opc_reference_check(a.out, a.db, age_days=a.age_days, w_b_range=[a.w_b_min, a.w_b_max], use_mock=not a.real, dat_lst=a.dat_lst, max_xgems_calls=a.max_xgems_calls, lit_db=a.lit_db, max_mixes=a.max_mixes, quantity=a.quantity)
    _print(r)
    return 0 if r["ok"] else 1


def cmd_infer(a: argparse.Namespace) -> int:
    from .pilot.tools_b_c import dor_infer_from_observations

    inp = _from_input(a)
    mix, observations, scm = (inp["mix"], inp["observations"], inp["scm"]) if inp else (_load(a.mix), _load(a.observations), _load(a.scm))
    r = dor_infer_from_observations(mix, observations, a.out, a.db, scm=scm, mix_uid=a.mix_uid, prior=a.prior, alpha_grid=a.alpha_grid, use_mock=not a.real, dat_lst=a.dat_lst, max_xgems_calls=a.max_xgems_calls, lit_db=a.lit_db)
    _print(r)
    return 0 if r["ok"] else 1


def cmd_stage(a: argparse.Namespace) -> int:
    from .pilot.tools_b_c import dor_stage_inferred

    r = dor_stage_inferred(a.inference, a.staging_db, dry_run=not a.write, note=a.note)
    _print(r)
    return 0 if r["ok"] else 1


def cmd_review(a: argparse.Namespace) -> int:
    from .inverse.staging import review_approve, review_list, review_reject

    if a.action == "list":
        _print(review_list(a.staging_db, reviewed=a.reviewed))
    elif a.action == "approve":
        _print(review_approve(a.staging_db, a.inference_id, note=a.note))
    else:
        _print(review_reject(a.staging_db, a.inference_id, note=a.note))
    return 0


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="dorgems", description="SCM degree-of-reaction agent kernel")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("validate", help="check a new-material input file (templates/new_material_template.yaml): schema + observation unit grades")
    p.add_argument("input")
    p.set_defaults(func=cmd_validate)

    p = sub.add_parser("predict", help="scenario A step 1: prediction.json")
    p.add_argument("--scm", required=True, help="YAML/JSON text or file (SCMSpec)")
    p.add_argument("--mix", required=True)
    p.add_argument("--ages", default=None)
    p.add_argument("--out", default="out/dor_predict")
    p.add_argument("--ensemble", default=None, choices=["blend", "bayes", "gbm_anchor_only"])
    p.add_argument("--method-group", default=None)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--lit-db", default=None)
    p.set_defaults(func=cmd_predict)

    p = sub.add_parser("export", help="prediction.json -> reaction_model_config YAMLs")
    p.add_argument("prediction")
    p.add_argument("--out", default="out/reaction_models")
    p.add_argument("--mode", default="logistic_fit", choices=["logistic_fit", "native", "pin"])
    p.add_argument("--slot", default=None)
    p.add_argument("--config-id", default=None)
    p.set_defaults(func=cmd_export)

    p = sub.add_parser("materials", help="materials.yaml override for a new SCM")
    p.add_argument("--scm", required=True)
    p.add_argument("--out", default="out/materials")
    p.add_argument("--slot", default=None)
    p.add_argument("--alias", default=None)
    p.add_argument("--cement", default=None)
    p.set_defaults(func=cmd_materials)

    p = sub.add_parser("analogues", help="closest literature mixes")
    p.add_argument("--scm", required=True)
    p.add_argument("--mix", required=True)
    p.add_argument("--age-days", type=float, default=None)
    p.add_argument("--limit", type=int, default=20)
    p.add_argument("--lit-db", default=None)
    p.set_defaults(func=cmd_analogues)

    p = sub.add_parser("run-forward", help="forward run with a DoRGems reaction model")
    p.add_argument("forward_query")
    p.add_argument("reaction_model_config")
    p.add_argument("--out", required=True)
    p.add_argument("--db", required=True)
    p.add_argument("--materials-config", default=None)
    p.add_argument("--slot", default=None)
    p.add_argument("--real", action="store_true", help="use_mock=False (requires --max-xgems-calls)")
    p.add_argument("--dat-lst", default=None)
    p.add_argument("--max-xgems-calls", type=int, default=None)
    p.add_argument("--capture-species", action="store_true")
    p.set_defaults(func=cmd_run_forward)

    p = sub.add_parser("envelope", help="scenario A end-to-end")
    p.add_argument("--input", default=None, help="new-material input file (replaces --scm/--mix/--ages)")
    p.add_argument("--scm", default=None)
    p.add_argument("--mix", default=None)
    p.add_argument("--ages", default=None)
    p.add_argument("--out", required=True)
    p.add_argument("--db", required=True)
    p.add_argument("--real", action="store_true")
    p.add_argument("--dat-lst", default=None)
    p.add_argument("--max-xgems-calls", type=int, default=None)
    p.add_argument("--ensemble", default=None)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--lit-db", default=None)
    p.set_defaults(func=cmd_envelope)

    p = sub.add_parser("compare-models", help="alpha(t) of two reaction_model_config YAMLs")
    p.add_argument("config_a")
    p.add_argument("config_b")
    p.add_argument("--slot", default=None)
    p.add_argument("--ages", default=None)
    p.set_defaults(func=cmd_compare_models)

    p = sub.add_parser("model-card")
    p.set_defaults(func=cmd_model_card)

    p = sub.add_parser("db-lookup", help="named read-only literature query")
    p.add_argument("query_name")
    p.add_argument("--params", default=None)
    p.add_argument("--limit", type=int, default=50)
    p.add_argument("--lit-db", default=None)
    p.set_defaults(func=cmd_db_lookup)

    p = sub.add_parser("compare", help="scenario B: forward run vs literature observations (or --input: your own material + observations)")
    p.add_argument("--input", default=None, help="new-material input file with observations")
    p.add_argument("--out", required=True)
    p.add_argument("--db", required=True)
    p.add_argument("--run-dir", default=None)
    p.add_argument("--mix-uid", default=None)
    p.add_argument("--scm", default=None)
    p.add_argument("--mix", default=None)
    p.add_argument("--mode", default="twin", choices=["twin", "neighbourhood"])
    p.add_argument("--quantities", default=None)
    p.add_argument("--real", action="store_true")
    p.add_argument("--dat-lst", default=None)
    p.add_argument("--max-xgems-calls", type=int, default=None)
    p.add_argument("--lit-db", default=None)
    p.set_defaults(func=cmd_compare)

    p = sub.add_parser("opc-check", help="scenario B: OPC-only reference check of the kernel")
    p.add_argument("--out", required=True)
    p.add_argument("--db", required=True)
    p.add_argument("--age-days", type=float, default=28)
    p.add_argument("--w-b-min", type=float, default=0.4)
    p.add_argument("--w-b-max", type=float, default=0.5)
    p.add_argument("--max-mixes", type=int, default=None)
    p.add_argument("--quantity", default="CH_TGA", choices=["CH_TGA", "CH_XRD", "bound_water", "chem_shrink"])
    p.add_argument("--real", action="store_true")
    p.add_argument("--dat-lst", default=None)
    p.add_argument("--max-xgems-calls", type=int, default=None)
    p.add_argument("--lit-db", default=None)
    p.set_defaults(func=cmd_opc_check)

    p = sub.add_parser("infer", help="scenario C: observations -> DoR posterior")
    p.add_argument("--input", default=None, help="new-material input file (replaces --scm/--mix/--observations)")
    p.add_argument("--mix", default=None)
    p.add_argument("--observations", default=None)
    p.add_argument("--out", required=True)
    p.add_argument("--db", required=True)
    p.add_argument("--scm", default=None)
    p.add_argument("--mix-uid", default=None)
    p.add_argument("--prior", default="model", choices=["model", "flat"])
    p.add_argument("--alpha-grid", type=int, default=21)
    p.add_argument("--real", action="store_true")
    p.add_argument("--dat-lst", default=None)
    p.add_argument("--max-xgems-calls", type=int, default=None)
    p.add_argument("--lit-db", default=None)
    p.set_defaults(func=cmd_infer)

    p = sub.add_parser("stage", help="stage an inference (dry-run unless --write)")
    p.add_argument("inference")
    p.add_argument("--staging-db", default=None)
    p.add_argument("--write", action="store_true")
    p.add_argument("--note", default=None)
    p.set_defaults(func=cmd_stage)

    p = sub.add_parser("review", help="review staged inferences")
    p.add_argument("action", choices=["list", "approve", "reject"])
    p.add_argument("inference_id", nargs="?")
    p.add_argument("--staging-db", default=None)
    p.add_argument("--reviewed", type=int, default=None)
    p.add_argument("--note", default=None)
    p.set_defaults(func=cmd_review)
    return ap


def main(argv: list[str] | None = None) -> int:
    ap = build_parser()
    a = ap.parse_args(argv)
    return int(a.func(a))


if __name__ == "__main__":
    sys.exit(main())
