"""Child-owned BLAS/OpenMP limits and independent Arrow CPU/I/O limits."""

import ctypes
import logging
import os
import sys

LOGGER = logging.getLogger("LDSC.workers")
_LIMITERS = []
_ACCELERATE = None
_CONTROL_ERRORS = []
_INITIALIZED = False
_BUILD_BACKENDS = {}


def _accelerate_state(*, configure=False):
    if sys.platform != "darwin":
        return None
    try:
        library = ctypes.CDLL("/System/Library/Frameworks/Accelerate.framework/Accelerate")
        if configure:
            setter = library.BLASSetThreading
            setter.argtypes, setter.restype = [ctypes.c_uint], ctypes.c_int
            if setter(1) != 0:
                raise RuntimeError("BLASSetThreading rejected single mode")
        getter = library.BLASGetThreading
        getter.argtypes, getter.restype = [], ctypes.c_uint
        observed = {0: "automatic", 1: "single"}.get(getter(), "unknown")
        return {"backend": "Apple Accelerate", "target": "single", "observed": observed,
                "verified": observed == "single"}
    except (AttributeError, OSError, RuntimeError) as error:
        return {"backend": "Apple Accelerate", "target": "single", "observed": None,
                "verified": False, "reason": str(error)}


def initialize_worker(log_level="INFO"):
    """Apply fixed limits only in a managed child, after loading both BLAS users.

    Higher inherited environment limits never override this policy. Unsupported
    numerical controls are diagnostic, not computation failures. Callers report
    runtime verification with ``log_native_thread_state`` at the work boundary.
    """
    global _INITIALIZED, _ACCELERATE
    if _INITIALIZED:
        return
    from ._logging import configure_package_logging
    configure_package_logging(logging.getLevelName(log_level) if isinstance(log_level, int) else log_level)
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                 "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[name] = "1"
    import numpy  # Load before runtime control.
    import scipy.linalg  # SciPy may load a separate BLAS.
    import pyarrow as pa
    from threadpoolctl import threadpool_limits
    for package in (numpy, scipy):
        _BUILD_BACKENDS[package.__name__] = getattr(package.__config__, "CONFIG", {}).get(
            "Build Dependencies", {}).get("blas", {}).get("name", "unknown")
    for api in ("blas", "openmp"):
        try:
            _LIMITERS.append(threadpool_limits(limits=1, user_api=api))
        except Exception as error:
            _CONTROL_ERRORS.append(f"{api} control unavailable: {error}")
    if any("accelerate" in name.lower() for name in _BUILD_BACKENDS.values()):
        _ACCELERATE = _accelerate_state(configure=True)
    pa.set_cpu_count(1)
    pa.set_io_thread_count(1)
    _INITIALIZED = True


def native_thread_state():
    """Inspect loaded pools in the executing process without changing them."""
    import pyarrow as pa
    from threadpoolctl import threadpool_info
    errors = list(_CONTROL_ERRORS)
    try:
        pools = threadpool_info()
    except Exception as error:
        pools = []
        errors.append(f"native pool inspection unavailable: {error}")
    accelerate = _accelerate_state() if _ACCELERATE is not None else None
    if _ACCELERATE is not None and not _ACCELERATE["verified"]:
        accelerate = _ACCELERATE
    if not any(pool["user_api"] == "blas" for pool in pools) and accelerate is None:
        errors.append("No supported BLAS backend found; BLAS limit is unverified")
    if not any(pool["user_api"] == "blas" for pool in pools):
        for package, backend in _BUILD_BACKENDS.items():
            if "accelerate" not in backend.lower():
                errors.append(f"{package} BLAS backend {backend!r} has no supported runtime inspection")
    verified = (not errors and all(pool["num_threads"] == 1 for pool in pools)
                and (accelerate is None or accelerate["verified"]))
    return {"pid": os.getpid(), "blas_target": 1, "openmp_target": 1, "pools": pools,
            "build_backends": dict(_BUILD_BACKENDS), "accelerate": accelerate,
            "arrow_cpu": pa.cpu_count(), "arrow_io": pa.io_thread_count(),
            "verified": verified, "diagnostics": errors}


def log_native_thread_state(logger=LOGGER):
    """Record observed backends/settings, including unmet or unverified limits."""
    state = native_thread_state()
    logger.log(logging.INFO if state["verified"] else logging.WARNING,
               "Worker pid=%d: BLAS target=1; OpenMP target=1; numerical limits %s; "
               "NumPy/SciPy BLAS=%s; native pools=%s; Accelerate=%s; Arrow CPU=%d; Arrow I/O=%d; diagnostics=%s.",
               os.getpid(), "verified" if state["verified"] else "unverified or unmet (continuing)",
               state["build_backends"], state["pools"], state["accelerate"],
               state["arrow_cpu"], state["arrow_io"], state["diagnostics"])
