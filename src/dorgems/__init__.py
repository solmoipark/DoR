"""DoRGems: SCM degree-of-reaction agent kernel.

The kernel is deterministic and has no LLM or network dependency. The
``dorgems.pilot`` subpackage is the only place that talks to an agent host.
"""

__version__ = "0.1.0"
BUNDLE_SCHEMA_VERSION = 1
PREDICTION_SCHEMA = "dorgems-prediction/1.0"
INFERENCE_SCHEMA = "dorgems-inference/1.0"
TOOL_CONTRACT = "inverse-gems-tool/1.0"
