"""G3-3: literature mixes with measured DoR AND CH/bound-water observations — hide the DoR,
infer it from the indirect observations with the real kernel, and compare (spec §9.6, §12).

Decision rule (spec G3-3): MAE < GBM LOPO (10.3 %p) or 90 % coverage >= 0.8 -> scenario C
is promoted; otherwise C is demoted to a consistency check.

Resumable: mixes with an existing <out>/<mix>/inference.json are re-read, not re-run
(pass --no-resume to recompute). Candidates without any grade A/B observation for the
likelihood quantities are skipped before any xGEMS call. Partial results are flushed to
g3_3_partial.json after every mix; the final summary is written even if interrupted.

Usage: python scripts/g3_3_real_hidden_dor.py --dat-lst <dat.lst> [--out work/g3_3]
       [--max-mixes N] [--alpha-grid 11] [--no-refine] [--ig-db DIR] [--no-resume]
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

LIK_QUANTITIES = ("CH_TGA", "CH_XRD", "bound_water", "chem_shrink")


def _has_usable_obs(db, mix_uid: str, weights: dict[str, float]) -> bool:
    from dorgems.db.units import harmonize

    mrow = db.mix(mix_uid)
    if mrow is None:
        return False
    for o in db.observations_for_mix(mix_uid, list(LIK_QUANTITIES)):
        if o["age_d"] and o["age_d"] > 0 and weights.get(o["quantity"], 1.0) > 0 and harmonize(o, mrow).usable:
            return True
    return False


def _summarise(rows: list[dict], cands: pd.DataFrame, a: argparse.Namespace, out: Path) -> dict:
    ok = [r for r in rows if r["ok"] and r.get("per_age")]
    all_pts = [dict(p, slot=r["slot"], mix=r["mix_uid"], paper=r["paper_doi"]) for r in ok for p in r["per_age"]]
    summary = {
        "n_candidates": int(len(cands)), "n_tried": len(rows), "n_ok": len(ok), "n_points": len(all_pts),
        "mae_pp": float(np.mean([abs(p["err_pp"]) for p in all_pts])) if all_pts else None,
        "coverage90": float(np.mean([p["covered"] for p in all_pts])) if all_pts else None,
        "bias_pp": float(np.mean([p["err_pp"] for p in all_pts])) if all_pts else None,
        "gbm_lopo_mae_reference_pp": 10.3, "prior": a.prior, "alpha_grid": a.alpha_grid, "refine": not a.no_refine,
        "xgems_calls_total": int(sum(r.get("xgems_calls") or 0 for r in ok)),
    }
    if all_pts:
        df = pd.DataFrame(all_pts)
        summary["by_slot"] = {s: {"n": int(len(g)), "mae_pp": float(g.err_pp.abs().mean()), "bias_pp": float(g.err_pp.mean()), "coverage90": float(g.covered.mean())} for s, g in df.groupby("slot")}
        df.to_csv(out / "g3_3_points.csv", index=False)
    summary["decision"] = "promote" if (summary["mae_pp"] is not None and (summary["mae_pp"] < 10.3 or (summary["coverage90"] or 0) >= 0.8)) else "demote_to_consistency_check"
    (out / "g3_3_real_hidden_dor.json").write_text(json.dumps({**summary, "results": rows}, indent=1, default=str), encoding="utf-8")
    return summary


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dat-lst", required=True)
    ap.add_argument("--out", default=None)
    ap.add_argument("--ig-db", default=None, help="kernel DB directory to reuse (default <out>/igdb)")
    ap.add_argument("--max-mixes", type=int, default=None)
    ap.add_argument("--alpha-grid", type=int, default=11)
    ap.add_argument("--no-refine", action="store_true", help="do not refine the alpha grid inside the prior 90 %% interval")
    ap.add_argument("--prior", default="model", choices=["model", "flat"])
    ap.add_argument("--no-resume", action="store_true")
    a = ap.parse_args()
    from dorgems.config import literature_db_path, modeling_dir
    from dorgems.db.features import build_dor_table
    from dorgems.db.reader import LiteratureDB
    from dorgems.inverse.run import infer_from_observations
    from dorgems.predict import load_defaults
    from dorgems.validate.twin_batch import candidate_mixes

    out = Path(a.out) if a.out else (modeling_dir(required=False) or Path(".")) / "work" / "g3_3"
    out.mkdir(parents=True, exist_ok=True)
    ig_db = Path(a.ig_db) if a.ig_db else out / "igdb"
    lit = literature_db_path()
    weights = {k: float(v) for k, v in (load_defaults().get("likelihood", {}).get("quantity_weights") or {}).items()}
    rows: list[dict] = []
    t0 = time.time()
    with LiteratureDB(lit) as db:
        cands = candidate_mixes(db, LIK_QUANTITIES, min_dor_ages=3, min_common_ages=3)
        table = build_dor_table(db.con)
        usable = [bool(_has_usable_obs(db, m, weights)) for m in cands["mix_uid"]]
    cands = cands.assign(usable=usable)
    print(f"candidates: {len(cands)}, with usable grade A/B observations: {int(cands.usable.sum())}", flush=True)
    todo = cands[cands.usable]
    if a.max_mixes:
        todo = todo.head(int(a.max_mixes))
    for r in cands[~cands.usable].itertuples():
        rows.append({"mix_uid": r.mix_uid, "paper_doi": r.paper_doi, "ok": False, "error": "no usable (grade A/B) observations for the likelihood (pre-filter)"})
    try:
        for k, c in enumerate(todo.itertuples(), 1):
            mix_uid = c.mix_uid
            sub = out / str(mix_uid).replace("/", "_").replace(":", "_")
            rec = {"mix_uid": mix_uid, "paper_doi": c.paper_doi, "ok": False}
            try:
                inf_path = sub / "inference.json"
                if inf_path.is_file() and not a.no_resume:
                    inf = json.loads(inf_path.read_text(encoding="utf-8"))
                    s = {"ages_d": inf["alpha"]["ages_d"], "alpha_q50": inf["alpha"]["q50"], "alpha_q05": inf["alpha"]["q05"], "alpha_q95": inf["alpha"]["q95"], "n_observations_used": inf["n_observations_used"], "ess": inf["ess"], "posterior_method": inf["posterior_method"], "kl": inf.get("prior_vs_posterior_kl"), "xgems_calls": None}
                    slot = inf.get("slot"); warnings = ["resumed from existing inference.json"]
                else:
                    r = infer_from_observations({}, [], out=sub, ig_db=ig_db, mix_uid=mix_uid, lit_db=lit, prior=a.prior, alpha_grid_n=a.alpha_grid, refine=not a.no_refine, use_mock=False, dat_lst=a.dat_lst, max_xgems_calls=200)
                    s = r["summary"]; slot = r["slot"]; warnings = r["warnings"][:8]
                dor = table[table["mix_uid"] == mix_uid].groupby("age_d")["dor_pct"].mean()
                ages_out = np.asarray(s["ages_d"], float)
                per_age = []
                for age, meas in dor.items():
                    la = np.log(float(age))
                    q50 = float(np.interp(la, np.log(ages_out), s["alpha_q50"])); q05 = float(np.interp(la, np.log(ages_out), s["alpha_q05"])); q95 = float(np.interp(la, np.log(ages_out), s["alpha_q95"]))
                    per_age.append({"age_d": float(age), "dor_measured_pct": float(meas), "post_q50_pct": q50 * 100, "post_q05_pct": q05 * 100, "post_q95_pct": q95 * 100, "covered": q05 * 100 <= meas <= q95 * 100, "err_pp": q50 * 100 - float(meas)})
                rec.update({"ok": True, "slot": slot, "n_obs_used": s["n_observations_used"], "ess": s["ess"], "method": s["posterior_method"], "kl": s["kl"], "xgems_calls": s.get("xgems_calls"), "per_age": per_age, "mae_pp": float(np.mean([abs(p["err_pp"]) for p in per_age])) if per_age else None, "coverage": float(np.mean([p["covered"] for p in per_age])) if per_age else None, "warnings": warnings})
            except Exception as exc:  # noqa: BLE001
                rec["error"] = f"{type(exc).__name__}: {exc}"
            rows.append(rec)
            print(f"[{k}/{len(todo)}] {mix_uid}: ok={rec['ok']} mae={rec.get('mae_pp')} cov={rec.get('coverage')} err={rec.get('error')} ({(time.time()-t0)/60:.1f} min)", flush=True)
            (out / "g3_3_partial.json").write_text(json.dumps(rows, indent=1, default=str), encoding="utf-8")
    finally:
        summary = _summarise(rows, cands, a, out)
        print(json.dumps({k: v for k, v in summary.items() if k != "by_slot"}, indent=2))
        print("by_slot:", json.dumps(summary.get("by_slot"), indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
