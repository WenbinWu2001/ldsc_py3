"""CLI entrypoint for the new LD-score estimation workflow."""

try:
    from ldscore.ldscore_workflow import main
except ImportError:  # pragma: no cover - package import fallback
    from .ldscore.ldscore_workflow import main


if __name__ == "__main__":
    main()
