"""Infrastructure shared by every workflow.

Path resolution, input preflight, workflow logging, progress reporting,
worker and CPU policy, Python keyword checks, CLI help formatting and the
shared result-file writers. This is tier 2 of the package layout: these
modules import only ``ldsc._core`` and ``ldsc.errors``.
"""
