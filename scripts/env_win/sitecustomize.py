"""Environment-level workaround (Windows, conda-forge BLAS): threadpoolctl raises
OSError 0xc06d007f (delay-load failure) when probing a BLAS/OpenMP DLL, which
aborts pymc.sample(). This makes threadpoolctl tolerant; it does NOT change any
model or script. Activate with PYTHONPATH=<this dir>.
"""
try:
    import threadpoolctl as _t

    _orig = _t.LibController.num_threads.fget

    def _safe_num_threads(self):
        try:
            return _orig(self)
        except OSError:
            return None

    _t.LibController.num_threads = property(_safe_num_threads)

    _orig_set = _t.LibController.set_num_threads

    def _safe_set(self, n):
        try:
            return _orig_set(self, n)
        except OSError:
            return None

    _t.LibController.set_num_threads = _safe_set
except Exception:  # pragma: no cover
    pass
