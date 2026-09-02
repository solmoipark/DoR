"""Kinetic curve forms, all in **fraction** units (0–1)."""

from __future__ import annotations

import numpy as np


def stretched_exp(t: np.ndarray | float, a_max: float, tau: float, beta: float = 0.5) -> np.ndarray:
    """alpha(t) = a_max · (1 − exp(−(t/τ)^β)), a_max as a fraction."""
    t = np.asarray(t, dtype=float)
    if tau <= 0:
        raise ValueError("tau must be positive")
    return np.clip(a_max * (1.0 - np.exp(-((t / tau) ** beta))), 0.0, 1.0)


def five_param_logistic(t: np.ndarray | float, A: float, B: float, C: float, D: float, G: float = 1.0) -> np.ndarray:
    """InverseGems default (scm_reaction.py:70-74): alpha = D + (A − D)/(1 + (t/C)^B)^G, clipped to [0,1]."""
    t = np.asarray(t, dtype=float)
    if C <= 0:
        raise ValueError("SCM logistic parameter C must be positive.")
    return np.clip(D + (A - D) / (1.0 + (t / C) ** B) ** G, 0.0, 1.0)


def pin_params(alpha: float) -> dict[str, float]:
    """Constant-alpha parameterisation validated in spec §1.2: alpha(t) ≡ alpha for all t."""
    a = float(np.clip(alpha, 0.0, 1.0))
    return {"A": a, "B": 1.0, "C": 1.0, "D": a, "G": 1.0}
