"""GBM-anchored importance re-weighting of the Bayesian posterior (spec §5.4).

    w_s ∝ Π_k N(g_k | alpha_s(t_k), σ_g²),   ESS = (Σw)² / Σw²

If ESS < ``ess_min`` the re-weighting is abandoned and the Bayesian draws are
returned unchanged with a ``model_disagreement`` warning. This is a heuristic:
both models were trained on the same data (see docs/model_card.md); G1-3
decides the default mode empirically.
"""

from __future__ import annotations

from dataclasses import replace
from typing import Any

import numpy as np

from .bayes import CurveDraws, stretched_exp_pct

DEFAULT_ANCHORS = (3.0, 7.0, 28.0, 90.0, 180.0)
ENSEMBLE_MODES = ("blend", "bayes", "gbm_anchor_only")


def importance_weights(draws: CurveDraws, anchors: np.ndarray, gbm_pct: np.ndarray, sigma_g: float, beta: float) -> tuple[np.ndarray, float]:
    alpha_k = stretched_exp_pct(np.asarray(anchors, float), draws.a_max, draws.tau, beta)  # (S,K)
    resid = (np.asarray(gbm_pct, float)[None, :] - alpha_k) / float(sigma_g)
    logw = -0.5 * np.sum(resid**2, axis=1)
    logw -= logw.max()
    w = np.exp(logw)
    w /= w.sum()
    ess = float(1.0 / np.sum(w**2))
    return w, ess


def reweight(
    draws: CurveDraws,
    anchors: np.ndarray | list[float],
    gbm_pct: np.ndarray | list[float],
    *,
    sigma_g: float = 12.0,
    beta: float = 0.5,
    ess_min: float = 100.0,
    rng_seed: int = 0,
) -> tuple[CurveDraws, dict[str, Any]]:
    """Return (draws, info). On success ``draws`` is a resampled (equal-weight) set of
    the same size; ``info`` holds ess, weights-applied flag and warnings."""
    w, ess = importance_weights(draws, np.asarray(anchors, float), np.asarray(gbm_pct, float), sigma_g, beta)
    info: dict[str, Any] = {"ess": ess, "ess_min": ess_min, "sigma_g": sigma_g, "anchors_d": list(map(float, anchors)), "applied": False, "warnings": []}
    if ess < ess_min:
        info["warnings"].append(f"model_disagreement: GBM anchors far from Bayesian prior (ESS={ess:.1f} < {ess_min:g}); ensemble falls back to bayes")
        return draws, info
    rng = np.random.default_rng(rng_seed)
    S = draws.a_max.shape[0]
    idx = rng.choice(S, size=S, replace=True, p=w)
    new = replace(
        draws,
        a_max=draws.a_max[idx],
        tau=draws.tau[idx],
        alpha=draws.alpha[idx],
        sigma_obs=None if draws.sigma_obs is None else draws.sigma_obs[idx],
        weights=None,
        flags=list(draws.flags) + ["ensemble_reweighted"],
    )
    info["applied"] = True
    return new, info


def fit_stretched_exp(anchors: np.ndarray | list[float], values_pct: np.ndarray | list[float], *, beta: float = 0.5) -> dict[str, float]:
    """Diagnostic 'gbm_anchor_only' mode: least-squares (a_max, tau) through the GBM anchors."""
    from scipy.optimize import least_squares

    t = np.asarray(anchors, float)
    y = np.asarray(values_pct, float)

    def resid(p: np.ndarray) -> np.ndarray:
        a, lt = p
        return a * (1.0 - np.exp(-((t / np.exp(lt)) ** beta))) - y

    a0 = max(float(y.max()) * 1.2, 5.0)
    r = least_squares(resid, x0=[a0, np.log(20.0)], bounds=([0.1, np.log(0.05)], [100.0, np.log(5000.0)]))
    a, lt = r.x
    return {"a_max": float(a), "tau_d": float(np.exp(lt)), "beta": beta, "rmse_pct": float(np.sqrt(np.mean(r.fun**2)))}
