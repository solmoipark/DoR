"""Register the DoRGems stretched-exponential as an InverseGems kinetics model (spec §6.1, mode native).

The registry is process-local (scm_reaction.py:40-63, no entry-point discovery), so
``register()`` must run in the same process **before** ``load_reaction_parameters``.
``dorgems.gems.forward`` and ``dorgems.pilot`` call it on import.
"""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np

MODEL_NAME = "dorgems_stretched_exp"
REQUIRED = ("a_max", "tau", "beta")
_registered = False


def _fn(t: np.ndarray, p: Mapping[str, float]) -> np.ndarray:
    tau = float(p["tau"])
    if tau <= 0:
        raise ValueError("dorgems_stretched_exp parameter tau must be positive.")
    return float(p["a_max"]) * (1.0 - np.exp(-((np.asarray(t, float) / tau) ** float(p["beta"]))))


def register() -> bool:
    """Idempotent. Returns True if InverseGems is importable and the model is registered."""
    global _registered
    if _registered:
        return True
    try:
        from inverse_gems.scm_reaction import _KINETICS_REGISTRY, register_scm_kinetics
    except Exception:  # noqa: BLE001
        return False
    if MODEL_NAME not in _KINETICS_REGISTRY:
        register_scm_kinetics(MODEL_NAME, required=REQUIRED, asymptote_key="a_max")(_fn)
    _registered = True
    return True


def is_registered() -> bool:
    try:
        from inverse_gems.scm_reaction import registered_scm_kinetics_models

        return MODEL_NAME in registered_scm_kinetics_models()
    except Exception:  # noqa: BLE001
        return False
