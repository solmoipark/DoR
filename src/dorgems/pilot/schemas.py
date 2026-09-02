"""Input schemas (spec §4.2). pydantic v2, ``extra="forbid"``.

The kernel accepts these models *or* plain dicts with the same keys; this module
is the single place where user input is validated.
"""

from __future__ import annotations

from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

ROLES = (
    "slag",
    "fly_ash",
    "metakaolin",
    "calcined_clay",
    "silica_fume",
    "limestone",
    "natural_pozzolan",
    "glass_powder",
    "steel_slag",
    "other",
)
OXIDES = ("CaO", "SiO2", "Al2O3", "Fe2O3", "MgO", "SO3", "Na2O", "K2O", "TiO2", "LOI")
QUANTITY_LITERAL = Literal["CH_TGA", "CH_XRD", "bound_water", "QXRD_phase", "chem_shrink", "DoR_SCM", "DoR_clinker"]


class SCMSpec(BaseModel):
    model_config = ConfigDict(extra="forbid")

    name: str
    role: Literal[ROLES]  # type: ignore[valid-type]
    oxides: dict[str, float]
    blaine_m2_kg: float | None = None
    d50_um: float | None = None
    bet_m2_g: float | None = None
    amorphous_pct: float | None = Field(default=None, ge=0, le=100, description="Measured only. Never 100 − Σ crystalline.")
    density_kg_m3: float | None = None

    @field_validator("oxides")
    @classmethod
    def _oxides(cls, v: dict[str, float]) -> dict[str, float]:
        unknown = set(v) - set(OXIDES)
        if unknown:
            raise ValueError(f"unknown oxide keys {sorted(unknown)}; allowed {OXIDES}")
        for k, x in v.items():
            if x is None or x < 0 or x > 100:
                raise ValueError(f"oxide {k} out of range: {x}")
        if "SiO2" not in v:
            raise ValueError("oxides must include SiO2")
        return v

    @model_validator(mode="after")
    def _sum(self) -> "SCMSpec":
        s = sum(x for k, x in self.oxides.items() if k != "LOI")
        if s > 105:
            raise ValueError(f"oxide sum (without LOI) {s:.1f} > 105")
        return self


class MixSpec(BaseModel):
    model_config = ConfigDict(extra="forbid")

    scm_pct: float = Field(ge=0, le=100)
    w_b: float = Field(gt=0, le=2.0)
    curing_temp_C: float = 20.0
    opc_oxides: dict[str, float] | None = None
    other_components: dict[str, float] = Field(default_factory=dict)

    @field_validator("other_components")
    @classmethod
    def _others(cls, v: dict[str, float]) -> dict[str, float]:
        for k, x in v.items():
            if x < 0 or x > 100:
                raise ValueError(f"component {k} out of range: {x}")
        return v

    @model_validator(mode="after")
    def _total(self) -> "MixSpec":
        tot = self.scm_pct + sum(self.other_components.values())
        if tot > 100:
            raise ValueError(f"scm_pct + other components = {tot:.1f} > 100")
        return self


class ObservationSpec(BaseModel):
    model_config = ConfigDict(extra="forbid")

    age_d: float = Field(gt=0)
    quantity: QUANTITY_LITERAL
    value: float
    unit: str
    phase_name: str | None = None
    method: str | None = None
    uncertainty: float | None = Field(default=None, ge=0)

    @model_validator(mode="after")
    def _phase(self) -> "ObservationSpec":
        if self.quantity in ("QXRD_phase", "DoR_clinker") and not self.phase_name:
            raise ValueError(f"{self.quantity} requires phase_name")
        return self


def coerce_scm(obj: Any) -> SCMSpec:
    return obj if isinstance(obj, SCMSpec) else SCMSpec.model_validate(obj)


def coerce_mix(obj: Any) -> MixSpec:
    return obj if isinstance(obj, MixSpec) else MixSpec.model_validate(obj)


def coerce_observations(items: Any) -> list[ObservationSpec]:
    return [o if isinstance(o, ObservationSpec) else ObservationSpec.model_validate(o) for o in (items or [])]
