"""Write InverseGems ``reaction_model_config`` YAML files from a DoRGems prediction (spec §6.1).

Modes
-----
logistic_fit  (default) fit the five-parameter logistic to each quantile curve; CLI-compatible.
native        ``model: dorgems_stretched_exp`` — in-process only (registry.py).
pin           constant alpha (scenario C alpha grid).

All alphas are **fractions**. Provenance goes to a sidecar ``<id>.provenance.json``
(never into the YAML: the whole YAML is hashed into ``reaction_model_signature``).
``signature_files:`` may list a materials override so it enters the signature.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import yaml

from .. import __version__
from .curves import pin_params, stretched_exp
from .fit import default_grid, logistic_fit

MODES = ("logistic_fit", "native", "pin")


def _yaml_dump(data: dict[str, Any]) -> str:
    return yaml.safe_dump(data, sort_keys=False, default_flow_style=False, allow_unicode=True)


def _write(path: Path, data: dict[str, Any], header: str | None = None) -> None:
    text = _yaml_dump(data)
    if header:
        text = "".join(f"# {line}\n" for line in header.splitlines()) + text
    path.write_text(text, encoding="utf-8")


def _quantile_key(q: float) -> str:
    return f"q{int(round(q * 100)):02d}"


def curve_from_prediction(prediction: dict[str, Any], q: float, source: str = "recommended") -> tuple[np.ndarray, np.ndarray]:
    """(ages_d, alpha_fraction) for one quantile from prediction.json."""
    ages = np.asarray(prediction["input"]["ages_d"], float)
    key = _quantile_key(q)
    if source == "recommended":
        rec = prediction["recommended"]
        vals = rec.get(f"alpha_pct_{key}")
        if vals is None:
            raise KeyError(f"recommended has no alpha_pct_{key}")
    else:
        vals = prediction[source]["alpha_pct"]["latent"][key]
    return ages, np.asarray(vals, float) / 100.0


def params_from_prediction(prediction: dict[str, Any], q: float) -> tuple[float, float, float]:
    """(a_max_fraction, tau_d, beta) of the recommended source's parametric quantiles."""
    src = prediction["recommended"]["source"]
    block = prediction[src] if src in prediction else prediction["bayes"]
    key = _quantile_key(q)
    a = float(block["a_max"][key]) / 100.0
    tau = float(block["tau_d"][key])
    beta = float(prediction.get("beta_shape", 0.5))
    return a, tau, beta


def export_reaction_model(
    prediction: dict[str, Any],
    out_dir: str | Path,
    *,
    mode: str = "logistic_fit",
    slot: str | None = None,
    quantiles: tuple[float, ...] = (0.05, 0.5, 0.95),
    config_id: str | None = None,
    signature_files: list[str] | None = None,
    fit_options: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Return {quantile_key: {path, params, fit, warnings}} and write the YAMLs + sidecars."""
    if mode not in MODES:
        raise ValueError(f"mode must be one of {MODES}")
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    slot = slot or prediction.get("slot") or prediction.get("role_gbm") or prediction.get("role_bayes")
    if slot not in ("slag", "fly_ash", "metakaolin", "silica_fume"):
        raise ValueError(f"slot {slot!r} is not an InverseGems SCM slot; map the role first (materials_override.slot_for_role)")
    config_id = config_id or prediction.get("id") or "pred"
    fo = dict(fit_options or {})
    grid = default_grid(fo.get("t_min_d", 0.25), fo.get("t_max_d", 730.0), int(fo.get("n_grid", 48)))
    results: dict[str, Any] = {}
    for q in quantiles:
        key = _quantile_key(q)
        cid = f"dorgems_{config_id}_{key}"
        warnings: list[str] = []
        a_max, tau, beta = params_from_prediction(prediction, q)
        # The parametric quantile curve (stretched-exp with quantile a_max/tau) is the
        # exported object; the pointwise quantile envelope is reported alongside.
        alpha_grid = stretched_exp(grid, a_max, tau, beta)
        fit = None
        if mode == "logistic_fit":
            fit = logistic_fit(
                grid,
                alpha_grid,
                fix_A_zero=bool(fo.get("fix_A_zero", True)),
                report_range_d=tuple(fo.get("report_range_d", (1.0, 365.0))),
                warn_dev_pct=float(fo.get("warn_dev_pct", 2.0)),
                fail_dev_pct=float(fo.get("fail_dev_pct", 5.0)),
            )
            if fit.status == "fail":
                raise RuntimeError(f"{cid}: {fit.message}")
            if fit.status == "warn":
                warnings.append(f"{cid}: {fit.message}")
            scm_block = {slot: {k: round(float(v), 6) for k, v in fit.params.items()}}
            header = None
        elif mode == "native":
            scm_block = {slot: {"model": "dorgems_stretched_exp", "a_max": round(a_max, 6), "tau": round(tau, 6), "beta": beta}}
            header = "model dorgems_stretched_exp is registered only in-process by dorgems.kinetics.registry;\nthe stand-alone inverse-gems CLI cannot load this file. Use the dorgems CLI/pilot."
        else:  # pin
            a28 = float(stretched_exp(np.array([28.0]), a_max, tau, beta)[0])
            scm_block = {slot: pin_params(a28)}
            header = None
            warnings.append("pin mode exports alpha(28 d) as a constant; use pin_reaction_model() for explicit alphas")
        data: dict[str, Any] = {"id": cid, "scm_reaction": scm_block, "availability_modifier": {"enabled": False}}
        if signature_files:
            data["signature_files"] = [str(p) for p in signature_files]
        path = out / f"reaction_parameters.{cid}.yaml"
        _write(path, data, header)
        prov = {
            "dorgems": __version__,
            "config_id": cid,
            "mode": mode,
            "slot": slot,
            "quantile": q,
            "source": prediction["recommended"]["source"],
            "a_max_fraction": a_max,
            "tau_d": tau,
            "beta": beta,
            "fit": fit.to_dict() if fit else None,
            "prediction_provenance": prediction.get("provenance"),
            "availability_modifier_disabled_reason": "replacement-level and CH-availability effects are already in the DoR model; avoid double counting",
            "warnings": warnings,
        }
        (out / f"{cid}.provenance.json").write_text(json.dumps(prov, indent=2), encoding="utf-8")
        results[key] = {"path": str(path), "id": cid, "params": scm_block[slot], "fit": fit.to_dict() if fit else None, "warnings": warnings, "a_max_fraction": a_max, "tau_d": tau}
    return results


def pin_reaction_model(alpha: float, slot: str, out_dir: str | Path, *, config_id: str = "pin", signature_files: list[str] | None = None) -> Path:
    """One YAML pinning alpha for ``slot`` at every age (spec §6.1 mode pin)."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    a = float(np.clip(alpha, 0.0, 1.0))
    cid = f"dorgems_{config_id}_pin_a{a:.3f}"
    data: dict[str, Any] = {"id": cid, "scm_reaction": {slot: pin_params(a)}, "availability_modifier": {"enabled": False}}
    if signature_files:
        data["signature_files"] = [str(p) for p in signature_files]
    path = out / f"reaction_parameters.{cid}.yaml"
    _write(path, data)
    return path


def alpha_from_config(config: str | Path, slot: str, ages: np.ndarray) -> np.ndarray:
    """Evaluate alpha(t) (fraction) of one slot from a reaction_model_config YAML.

    Uses InverseGems' own evaluator when importable (so registered models work);
    otherwise falls back to the local logistic/stretched-exp implementations."""
    raw = yaml.safe_load(Path(config).read_text(encoding="utf-8")) or {}
    block = (raw.get("scm_reaction") or raw.get("scm_parameters") or {}).get(slot)
    if block is None:
        raise KeyError(f"{config} has no scm_reaction.{slot}")
    try:
        from .registry import register

        register()
        from inverse_gems.scm_reaction import scm_alpha

        return np.asarray(scm_alpha(list(map(float, ages)), block), float)
    except ImportError:
        pass
    from .curves import five_param_logistic

    model = block.get("model", "five_param_logistic")
    if model == "five_param_logistic":
        return five_param_logistic(np.asarray(ages, float), block["A"], block["B"], block["C"], block["D"], block.get("G", 1.0))
    if model == "dorgems_stretched_exp":
        return stretched_exp(np.asarray(ages, float), block["a_max"], block["tau"], block.get("beta", 0.5))
    raise ValueError(f"unknown kinetics model {model!r}")


def compare_reaction_models(config_a: str | Path, config_b: str | Path, slot: str, ages: np.ndarray | None = None) -> dict[str, Any]:
    ages = np.asarray(ages if ages is not None else [1, 3, 7, 28, 90, 180, 365], float)
    a = alpha_from_config(config_a, slot, ages)
    b = alpha_from_config(config_b, slot, ages)
    return {
        "slot": slot,
        "ages_d": ages.tolist(),
        "alpha_a": a.tolist(),
        "alpha_b": b.tolist(),
        "diff_pct": (100.0 * (a - b)).tolist(),
        "max_abs_diff_pct": float(np.max(np.abs(100.0 * (a - b)))),
    }
