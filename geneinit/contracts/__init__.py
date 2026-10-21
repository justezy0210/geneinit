"""Input, output, and provenance contracts."""

from importlib import import_module

__all__ = [
    "PreflightError",
    "SampleContractError",
    "SampleRecord",
    "load_samples",
    "run_preflight",
]


def __getattr__(name):
    # Evidence readers use streaming without loading deployment/preflight code.
    if name in {"PreflightError", "run_preflight"}:
        module = ".preflight"
    elif name in {"SampleContractError", "SampleRecord", "load_samples"}:
        module = ".samples"
    else:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    value = getattr(import_module(module, __name__), name)
    globals()[name] = value
    return value


def __dir__():
    return sorted(set(globals()) | set(__all__))
