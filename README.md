# DoRGems

SCM degree-of-reaction (DoR) agent kernel for blended cements. It turns the
literature database (`scm_dor_enriched.db`, 1,350 papers / 44,011 observations)
and the DoR models trained on it into deterministic tools that

* **A.** predict the DoR curve of a *new* SCM with uncertainty and export it as
  InverseGems `reaction_model_config` files (q05 / q50 / q95 envelope);
* **B.** validate xGEMS forward results against literature portlandite, bound
  water and chemical-shrinkage observations;
* **C.** infer a DoR curve from indirect observations (alpha-grid forward map +
  Bayesian re-weighting) and stage the result for review.

The LLM never produces a number: it selects and parameterises tools. The
literature DB is opened read-only; agent-produced values go to a separate
staging DB with `reviewed = 0`.

Specification: [docs/DoRGems_agent_spec_v1.md](docs/DoRGems_agent_spec_v1.md).
Gate log: [docs/gates.md](docs/gates.md). Model card: [docs/model_card.md](docs/model_card.md).

## Install

```bash
pip install -e ../InverseGems          # kernel (BSD-3)
pip install -e ../GemsPilot --no-deps  # agent layer (optional)
pip install -e .[pilot,test]
```

Environment: `DORGEMS_DB` (literature DB), `DORGEMS_STAGING_DB`, `DORGEMS_BUNDLE`,
`DORGEMS_MODELING_DIR`, `INVERSE_GEMS_ROOT`, `DORGEMS_REAL_XGEMS=1` (real-xGEMS tests).
Without `DORGEMS_DB` the DB is discovered under `../DoR of SCMs in blended cements/modeling/`.

## CLI

```bash
dorgems predict   --scm scm.yaml --mix mix.yaml --ages "[1,3,7,28,90,365]" --out out/A
dorgems export    out/A/prediction.json --out out/A/rm
dorgems materials --scm scm.yaml --out out/A
dorgems envelope  --scm scm.yaml --mix mix.yaml --out out/A --db out/igdb          # scenario A, mock
dorgems compare   --mix-uid "<doi>::<mix_id>" --out out/B --db out/igdb            # scenario B twin
dorgems opc-check --out out/B_opc --db out/igdb --max-mixes 20
dorgems infer     --scm scm.yaml --mix mix.yaml --observations obs.yaml --out out/C --db out/igdb
dorgems stage     out/C/inference.json            # dry-run; --write to stage
dorgems review    list | approve <id> | reject <id>
dorgems model-card
```

Real xGEMS runs need `--real --dat-lst <Test-dat.lst> --max-xgems-calls N` (N ≤ 200).

## Agent hosts

* `dorgems-mcp` — stand-alone MCP server with the 13 `dor_*` tools.
* GemsPilot: `dorgems.pilot.tools.TOOLSET` (entry point `gemspilot.toolsets`) once
  P-GP-2 is merged; the tools already follow GemsPilot's `ToolResult` contract and
  `_policy_check` semantics (`read` / `mock_ok`, never `real_gated`).

## Layout

`src/dorgems/{db,models,kinetics,gems,validate,inverse,pilot}` — see the spec §3.2.
`bundles/` holds the frozen inference bundles (`bayes_v4`, `gbm_v6`, `ood_reference.json`);
`scripts/` the one-off export / re-run helpers; `tests/` the gate tests (mock by default).

## Windows notes (this machine)

PyMC runs with `PYTENSOR_FLAGS=mode=NUMBA,cxx=` (no g++) and
`PYTHONPATH=scripts/env_win` (threadpoolctl workaround); the conda env uses
OpenBLAS (`libblas=*=*openblas`) because the MKL build crashed on delay-load.
