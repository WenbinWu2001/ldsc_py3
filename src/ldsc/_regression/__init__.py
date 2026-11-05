"""LD score regression: estimators, regression workflows and post-fit results.

``estimators`` holds the LDSC estimators with their jackknife and IRWLS
routines; ``alignment`` assembles regression datasets from saved LD-score
directories and summary statistics; ``fitting`` fits them and builds result
rows; ``partitioned`` and ``rg`` run partitioned-h2 batches and rg pairs;
``commands`` holds the ``h2``, ``partitioned-h2`` and ``rg`` entry points and
``runner`` the public ``RegressionRunner``; ``h2_scale`` and ``quantile_h2``
post-process saved fits; ``results`` holds the result types and ``outputs``
writes the result directories. This is tier 7 of the package
layout: it reads LD-score directories and never recomputes LD scores.
"""
