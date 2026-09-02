"""G3-3: literature mixes with measured DoR AND CH/bound-water observations — hide the DoR,
infer it from the indirect observations with the real kernel, and compare (spec §9.6, §12).

Decision rule (spec G3-3): MAE < GBM LOPO (10.3 %p) or 90 % coverage >= 0.8 -> scenario C
is promoted; otherwise C is demoted to a consistency check.

Usage: python scripts/g3_3_real_hidden_dor.py --dat-lst <dat.lst> [--max-mixes N] [--out work/g3_3]
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


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dat-lst", required=True)
    ap.add_argument("--out", default=None)
    ap.add_argument("--max-mixes", type=int, default=None)
    ap.add_argument("--alpha-grid", type=int, default=21)
    ap.add_argument("--prior", default="model", choices=["model", "flat"])
    a = ap.parse_args()
    from dorgems.config import literature_db_path, modeling_dir
    from dorgems.db.features import build_dor_table
    from dorgems.db.reader import LiteratureDB
    from dorgems.inverse.run import infer_from_observations
    from dorgems.kinetics.curves import stretched_exp
    from dorgems.validate.twin_batch import candidate_mixes

    out = Path(a.out) if a.out else (modeling_dir(required=False) or Path(".")) / "work" / "g3_3"
    out.mkdir(parents=True, exist_ok=True)
    lit = literature_db_path()
    rows = []
    t0 = time.time()
    with LiteratureDB(lit) as db:
        cands = candidate_mixes(db, ("CH_TGA", "CH_XRD", "bound_water", "chem_shrink"), min_dor_ages=3, min_common_ages=3)
        table = build_dor_table(db.con)
    if a.max_mixes:
        cands = cands.head(int(a.max_mixes))
    print(f"candidates: {len(cands)}", flush=True)
    for k, c in cands.iterrows():
        mix_uid = c["mix_uid"]
        sub = out / str(mix_uid).replace("/", "_").replace(":", "_")
        rec = {"mix_uid": mix_uid, "paper_doi": c["paper_doi"], "ok": False}
        try:
            r = infer_from_observations({}, [], out=sub, ig_db=out / "igdb", mix_uid=mix_uid, lit_db=lit, prior=a.prior, alpha_grid_n=a.alpha_grid, use_mock=False, dat_lst=a.dat_lst, max_xgems_calls=200)
            s = r["summary"]
            dor = table[table["mix_uid"] == mix_uid].groupby("age_d")["dor_pct"].mean()
            ages_out = np.asarray(s["ages_d"], float)
            per_age = []
            for age, meas in dor.items():
                # posterior alpha at the measured age from the (a_max, tau) quantiles is not
                # available directly; interpolate the q05/q50/q95 curves in log-age
                la = np.log(float(age))
                q50 = float(np.interp(la, np.log(ages_out), s["alpha_q50"]))
                q05 = float(np.interp(la, np.log(ages_out), s["alpha_q05"]))
                q95 = float(np.interp(la, np.log(ages_out), s["alpha_q95"]))
                per_age.append({"age_d": float(age), "dor_measured_pct": float(meas), "post_q50_pct": q50 * 100, "post_q05_pct": q05 * 100, "post_q95_pct": q95 * 100, "covered": q05 * 100 <= meas <= q95 * 100, "err_pp": q50 * 100 - float(meas)})
            rec.update({"ok": True, "slot": r["slot"], "n_obs_used": s["n_observations_used"], "ess": s["ess"], "method": s["posterior_method"], "kl": s["kl"], "xgems_calls": s["xgems_calls"], "per_age": per_age, "mae_pp": float(np.mean([abs(p["err_pp"]) for p in per_age])) if per_age else None, "coverage": float(np.mean([p["covered"] for p in per_age])) if per_age else None, "warnings": r["warnings"][:8]})
        except Exception as exc:  # noqa: BLE001
            rec["error"] = f"{type(exc).__name__}: {exc}"
        rows.append(rec)
        print(f"[{k+1}/{len(cands)}] {mix_uid}: ok={rec['ok']} mae={rec.get('mae_pp')} cov={rec.get('coverage')} err={rec.get('error')} ({(time.time()-t0)/60:.1f} min)", flush=True)
        (out / "g3_3_partial.json").write_text(json.dumps(rows, indent=1, default=str), encoding="utf-8")
    ok = [r for r in rows if r["ok"] and r.get("per_age")]
    all_pts = [p for r in ok for p in r["per_age"]]
    summary = {
        "n_candidates": int(len(cands)),
        "n_ok": len(ok),
        "n_points": len(all_pts),
        "mae_pp": float(np.mean([abs(p["err_pp"]) for p in all_pts])) if all_pts else None,
        "coverage90": float(np.mean([p["covered"] for p in all_pts])) if all_pts else None,
        "bias_pp": float(np.mean([p["err_pp"] for p in all_pts])) if all_pts else None,
        "gbm_lopo_mae_reference_pp": 10.3,
        "prior": a.prior,
        "xgems_calls_total": int(sum(r.get("xgems_calls") or 0 for r in ok)),
    }
    summary["decision"] = "promote" if (summary["mae_pp"] is not None and (summary["mae_pp"] < 10.3 or (summary["coverage90"] or 0) >= 0.8)) else "demote_to_consistency_check"
    (out / "g3_3_real_hidden_dor.json").write_text(json.dumps({**summary, "results": rows}, indent=1, default=str), encoding="utf-8")
    pd.DataFrame(all_pts).to_csv(out / "g3_3_points.csv", index=False)
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
