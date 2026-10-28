"""Plots of saved results, drawn without refitting any model.

``dispatch`` reads a result directory's metadata, chooses the plot family
and writes it under ``<result-dir>/plots/``; ``builders`` holds the
Matplotlib figure code, imported only when a figure is drawn. This is tier 8
of the package layout.
"""
