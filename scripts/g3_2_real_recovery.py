"""G3-2: synthetic recovery on the REAL kernel (5 cases, spec §9.6 / §12).

For each (role, scm_pct, w_b, a_max*, tau*) case: build the alpha-grid forward map with
real xGEMS, generate CH + bound-water "observations" from alpha*(t) with noise, and check
that alpha(28 d) is inside the posterior 90 % interval and |median error| <= 5 %p.

Usage: python scripts/g3_2_real_recovery.py --dat-lst <dat.lst> [--out work/g3_2]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

CASES = [
    {"role": "slag", "scm_pct": 40, "w_b": 0.45, "a_max": 0.55, "tau": 15.0, "oxides": {"CaO": 41.8, "SiO2": 36.5, "Al2O3": 12.3, "MgO": 7.5, "SO3": 2.0}},
    {"role": "slag", "scm_pct": 60, "w_b": 0.40, "a_max": 0.40, "tau": 25.0, "oxides": {"CaO": 41.8, "SiO2": 36.5, "Al2O3": 12.3, "MgO": 7.5, "SO3": 2.0}},
    {"role": "fly_ash", "scm_pct": 30, "w_b": 0.50, "a_max": 0.30, "tau": 30.0, "oxides": {"CaO": 5.0, "SiO2": 55.0, "Al2O3": 25.0, "Fe2O3": 6.0, "MgO": 1.5}},
    {"role": "fly_ash", "scm_pct": 50, "w_b": 0.45, "a_max": 0.20, "tau": 60.0, "oxides": {"CaO": 5.0, "SiO2": 55.0, "Al2O3": 25.0, "Fe2O3": 6.0, "MgO": 1.5}},
    {"role": "metakaolin", "scm_pct": 20, "w_b": 0.50, "a_max": 0.70, "tau": 5.0, "oxides": {"CaO": 0.2, "SiO2": 54.1, "Al2O3": 43.6, "Fe2O3": 1.1}},
]
AGES = [3.0, 7.0, 28.0, 90.0]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dat-lst", required=True)
    ap.add_argument("--out", default=None)
    ap.add_argument("--noise-frac", type=float, default=0.03, help="noise sd as a fraction of the observable span over alpha")
    ap.add_argument("--n-alpha", type=int, default=21)
    a = ap.parse_args()
    from dorgems.config import modeling_dir
    from dorgems.envelope import build_forward_query
    from dorgems.inverse.alpha_grid import build_forward_map, default_alpha_grid
    from dorgems.inverse.likelihood import Likelihood, build_points
    from dorgems.inverse.posterior import infer, summarise
    from dorgems.kinetics.curves import stretched_exp
    from dorgems.kinetics.materials_override import build_materials_config, slot_for_role
    from dorgems.pilot.schemas import SCMSpec

    out = Path(a.out) if a.out else (modeling_dir(required=False) or Path(".")) / "work" / "g3_2"
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(42)
    results = []
    for i, c in enumerate(CASES):
        case_out = out / f"case{i}"
        scm = SCMSpec(name=f"g32_{c['role']}", role=c["role"], oxides=c["oxides"])
        slot, _, _ = slot_for_role(c["role"], c["oxides"])
        mat = build_materials_config(scm, case_out, slot=slot)
        fq = build_forward_query({"scm_pct": c["scm_pct"], "w_b": c["w_b"]}, slot, AGES, name=f"g32_{i}")
        fmap = build_forward_map(fq, slot=slot, ages=AGES, alphas=default_alpha_grid(a.n_alpha), out=case_out / "fmap", ig_db=out / "igdb", materials_config=mat["path"], use_mock=False, dat_lst=a.dat_lst, max_xgems_calls=200, quantities=("CH_TGA", "bound_water"))
        obs = []
        for j, t in enumerate(AGES):
            alpha_t = float(stretched_exp(np.array([t]), c["a_max"], c["tau"], 0.5)[0])
            for q, key in (("CH_TGA", "CH_g"), ("bound_water", "bound_water_g")):
                row = fmap.table[key][j]
                if not np.isfinite(row).sum() >= 2:
                    continue
                span = float(np.nanmax(row) - np.nanmin(row))
                sig = max(a.noise_frac * span, 0.05)
                obs.append({"obs_uid": f"{q}@{t}", "age_d": t, "quantity": q, "value": float(fmap.value(q, j, alpha_t)) + rng.normal(0, sig), "grade": "A", "uncertainty": sig})
        pts, skipped = build_points(obs, AGES, sigma_obs_default={}, sigma_model={"CH_TGA": 0.5, "bound_water": 0.5})
        lik = Likelihood(fmap, pts, beta=0.5)
        post = infer(lik, prior_a_max=None, prior_tau=None, prior="flat", grid_n=40, rng_seed=i)
        summ = summarise(post, lik, np.array(AGES))
        a28 = float(stretched_exp(np.array([28.0]), c["a_max"], c["tau"], 0.5)[0])
        q05, q50, q95 = summ["alpha"]["q05"][2], summ["alpha"]["q50"][2], summ["alpha"]["q95"][2]
        r = {"case": i, **{k: c[k] for k in ("role", "scm_pct", "w_b", "a_max", "tau")}, "slot": slot, "alpha28_true": a28, "q05": q05, "q50": q50, "q95": q95, "covered": q05 <= a28 <= q95, "abs_err_pp": abs(q50 - a28) * 100, "ess": summ["ess"], "method": summ["posterior_method"], "n_obs": len(pts), "xgems_calls": fmap.meta.get("xgems_calls"), "monotonicity": fmap.monotonicity_report()}
        print(json.dumps({k: v for k, v in r.items() if k != "monotonicity"}), flush=True)
        results.append(r)
    cov = float(np.mean([r["covered"] for r in results]))
    mae = float(np.mean([r["abs_err_pp"] for r in results]))
    summary = {"n": len(results), "coverage90": cov, "mae_pp": mae, "pass": cov >= 0.85 and mae <= 5.0, "results": results}
    (out / "g3_2_real_recovery.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    print("G3-2:", {"coverage90": cov, "mae_pp": mae, "pass": summary["pass"]})
    return 0 if summary["pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
