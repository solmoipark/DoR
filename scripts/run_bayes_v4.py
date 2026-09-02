"""Re-run modeling/bayes_hier_v4.py unmodified, from the modeling directory.

Why a wrapper: the script calls pm.sample(cores=4) at module level. On Windows
multiprocessing uses 'spawn', which re-imports the main module in every child;
a script without an ``if __name__ == "__main__"`` guard would re-execute
itself. Running it through runpy from this guarded file avoids that without
touching the script. Environment: PYTHONPATH=scripts/env_win (threadpoolctl
workaround), PYTENSOR_FLAGS=mode=NUMBA,cxx= (no g++ on this PC).

Usage:  python scripts/run_bayes_v4.py <modeling_dir> [table_csv]
Outputs land in <modeling_dir>/work/ exactly as the script writes them.
"""

from __future__ import annotations

import os
import runpy
import sys
import time
from pathlib import Path


def main() -> None:
    modeling = Path(sys.argv[1]).resolve()
    table = sys.argv[2] if len(sys.argv) > 2 else "dor_scm_final.csv"
    os.chdir(modeling)
    (modeling / "work").mkdir(exist_ok=True)
    sys.argv = ["bayes_hier_v4.py", table]
    t0 = time.time()
    print(f"[run_bayes_v4] cwd={modeling} table={table} start={time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)
    # run_name must NOT be "__main__": runpy would swap sys.modules["__main__"] for the
    # script, and spawned chains would then re-execute the unguarded script.
    runpy.run_path(str(modeling / "bayes_hier_v4.py"), run_name="__dorgems_bayes_v4__")
    print(f"[run_bayes_v4] done in {(time.time() - t0) / 60:.1f} min", flush=True)


if __name__ == "__main__":
    main()
