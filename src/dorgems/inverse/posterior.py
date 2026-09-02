"""Posterior over (a_max, τ) by importance re-weighting of the model prior draws,
SIR refresh when ESS is low, and a grid fallback (spec §9.4)."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .. import INFERENCE_SCHEMA, __version__
from ..kinetics.curves import stretched_exp
from ..models.bayes import weighted_quantile
from .likelihood import Likelihood


def _ess(w: np.ndarray) -> float:
    return float(1.0 / np.sum(w**2))


def _normalise(logw: np.ndarray) -> np.ndarray:
    logw = logw - np.max(logw)
    w = np.exp(logw)
    return w / w.sum()


def _kl(w: np.ndarray) -> float:
    """KL(posterior ‖ prior) for a prior represented by equal-weight draws: Σ w log(S w)."""
    S = len(w)
    nz = w > 0
    return float(np.sum(w[nz] * np.log(S * w[nz])))


def infer(
    lik: Likelihood,
    *,
    prior_a_max: np.ndarray | None,
    prior_tau: np.ndarray | None,
    prior: str = "model",
    ess_min: float = 50.0,
    sir_rounds: int = 1,
    grid_n: int = 40,
    rng_seed: int = 0,
    b_bw_prior: tuple[float, float] | None = None,
    b_bw_points: int = 5,
) -> dict[str, Any]:
    rng = np.random.default_rng(rng_seed)
    offsets_grid = [0.0]
    if b_bw_prior is not None and any(p.quantity == "bound_water" for p in lik.points):
        m, s = b_bw_prior
        offsets_grid = list(np.linspace(m - 2 * s, m + 2 * s, int(b_bw_points))) if s > 0 else [m]
    method = "importance"
    if prior == "model" and prior_a_max is not None:
        a, t = np.asarray(prior_a_max, float), np.asarray(prior_tau, float)
    else:
        prior = "flat"
        S = int(grid_n) ** 2
        a = rng.uniform(0.02, 1.0, S)
        t = np.exp(rng.uniform(np.log(0.5), np.log(3000.0), S))
        method = "importance_flat"

    def logL(a_: np.ndarray, t_: np.ndarray) -> np.ndarray:
        # marginalise the bound-water offset on a small grid (equal prior weight)
        if len(offsets_grid) == 1:
            lik.offsets["bound_water"] = float(offsets_grid[0])
            return lik.loglik(a_, t_)
        acc = np.full(a_.shape[0], -np.inf)
        for b in offsets_grid:
            lik.offsets["bound_water"] = float(b)
            acc = np.logaddexp(acc, lik.loglik(a_, t_))
        return acc - np.log(len(offsets_grid))

    ll = logL(a, t)
    w = _normalise(ll)
    ess = _ess(w)
    rounds = 0
    while ess < ess_min and rounds < sir_rounds:
        rounds += 1
        idx = rng.choice(len(a), size=2 * len(a), replace=True, p=w)
        sd_a = max(np.std(a) * 0.1, 1e-3)
        sd_lt = max(np.std(np.log(t)) * 0.1, 1e-3)
        a = np.clip(a[idx] + rng.normal(0, sd_a, len(idx)), 0.01, 1.0)
        t = np.exp(np.log(t[idx]) + rng.normal(0, sd_lt, len(idx)))
        ll = logL(a, t)
        w = _normalise(ll)
        ess = _ess(w)
        method = "sir"
    if ess < ess_min:
        # grid posterior
        ga = np.linspace(0.02, 1.0, int(grid_n))
        gt = np.exp(np.linspace(np.log(0.5), np.log(3000.0), int(grid_n)))
        A, T = np.meshgrid(ga, gt, indexing="ij")
        a, t = A.ravel(), T.ravel()
        ll = logL(a, t)
        w = _normalise(ll)
        ess = _ess(w)
        method = "grid"
    kl = _kl(w) if method in ("importance",) else None
    return {"a_max": a, "tau": t, "weights": w, "loglik": ll, "ess": ess, "method": method, "prior": prior, "prior_vs_posterior_kl": kl, "offset_grid": offsets_grid}


def summarise(post: dict[str, Any], lik: Likelihood, ages_out: np.ndarray, *, prior_a_max: np.ndarray | None = None, prior_tau: np.ndarray | None = None, quantiles: tuple[float, ...] = (0.05, 0.5, 0.95)) -> dict[str, Any]:
    a, t, w = post["a_max"], post["tau"], post["weights"]
    qs = list(quantiles)
    keys = [f"q{int(round(q * 100)):02d}" for q in qs]
    out: dict[str, Any] = {"a_max": dict(zip(keys, map(float, weighted_quantile(a, qs, w)))), "tau_d": dict(zip(keys, map(float, weighted_quantile(t, qs, w))))}
    alpha = np.stack([stretched_exp(ages_out, ai, ti, lik.beta) for ai, ti in zip(a, t)])  # (S,T)
    arr = np.asarray([weighted_quantile(alpha[:, j], qs, w) for j in range(alpha.shape[1])])  # (T, nq)
    out["alpha"] = {"ages_d": ages_out.tolist(), **{k: [float(x) for x in arr[:, i]] for i, k in enumerate(keys)}, "mean": [float(v) for v in np.average(alpha, axis=0, weights=w)]}
    if prior_a_max is not None:
        pa = np.stack([stretched_exp(ages_out, ai, ti, lik.beta) for ai, ti in zip(prior_a_max, prior_tau)])
        pcols = np.asarray([np.quantile(pa[:, j], qs) for j in range(pa.shape[1])])
        out["prior_alpha"] = {k: [float(x) for x in pcols[:, i]] for i, k in enumerate(keys)}
        out["prior_a_max"] = dict(zip(keys, map(float, np.quantile(prior_a_max, qs))))
        out["prior_tau_d"] = dict(zip(keys, map(float, np.quantile(prior_tau, qs))))
    # posterior predictive per observation + likelihood contribution decomposition
    pred = lik.predicted(a, t)  # (n_points, S)
    _, contrib = lik.loglik(a, t, per_point=True)
    ppc = []
    for k, p in enumerate(lik.points):
        m = pred[k]
        q05, q50, q95 = weighted_quantile(m, [0.05, 0.5, 0.95], w)
        ppc.append({"label": p.label, "quantity": p.quantity, "age_d": p.age_d, "obs": p.y, "pred_q05": float(q05), "pred_q50": float(q50), "pred_q95": float(q95), "residual_q50": float(q50 - p.y), "z": float((q50 - p.y) / np.sqrt(p.var)), "loglik_contribution_mean": float(np.average(contrib[k], weights=w)), "sigma_obs": p.sigma_obs, "sigma_model": p.sigma_model})
    out["ppc"] = ppc
    out["ess"] = post["ess"]
    out["posterior_method"] = post["method"]
    out["prior_vs_posterior_kl"] = post["prior_vs_posterior_kl"]
    out["n_observations_used"] = len(lik.points)
    return out


def write_inference(summary: dict[str, Any], out: Path, *, inference_id: str, header: dict[str, Any]) -> dict[str, str]:
    out.mkdir(parents=True, exist_ok=True)
    payload = {"schema": INFERENCE_SCHEMA, "id": inference_id, "dorgems": __version__, **header, **summary}
    p = out / "inference.json"
    p.write_text(json.dumps(payload, indent=2, ensure_ascii=False, default=_jd), encoding="utf-8")
    pd.DataFrame(summary["ppc"]).to_csv(out / "ppc.csv", index=False)
    return {"inference": str(p), "ppc": str(out / "ppc.csv")}


def _jd(o: Any) -> Any:
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, np.bool_):
        return bool(o)
    return str(o)
