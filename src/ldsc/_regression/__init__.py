"""LD score regression: estimators, regression workflows and post-fit results.

``estimators`` holds the LDSC estimators with their jackknife and IRWLS
routines; ``runner`` assembles regression datasets from saved LD-score
directories and summary statistics and runs ``h2``, ``partitioned-h2`` and
``rg``; ``h2_scale`` and ``quantile_h2`` post-process saved fits;
``outputs`` writes the result directories. This is tier 7 of the package
layout: it reads LD-score directories and never recomputes LD scores.
"""
