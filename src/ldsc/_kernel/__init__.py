"""Internal compute kernels for the refactored LDSC package."""

__all__ = ["ref_panel_builder"]


def __getattr__(name: str):
    if name == "ref_panel_builder":
        from importlib import import_module

        module = import_module(f"{__name__}.ref_panel_builder")
        globals()[name] = module
        return module
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
