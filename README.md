# Paternal leaves

Supplementary material for the paper *Adaptation to paternal leave policies in Finnish municipalities: changing gender norms and cross-border policy legacies*.

A preprint is available at https://doi.org/10.31235/osf.io/k27yw.

## Files

This repository contains the following supplementary files:
- `supplementary_material.pdf`: the online supplementary material including descriptions of the political and attitudinal backgrounds, data, results, and model comparisons.
- `data_leaves.rds`: the data, describing utilized paternal leave quotas and other demographic information on Finnish municipalities in 2009-2017. More detailed description of the data can be found from file `run_models.R`.
- `data_leaves.xlsx`, `data_leaves.csv`: the data in formats easily accessible for programs other than R. These datasets do not include the geometry information for drawing the maps.
- `model_bin_nophi.stan`, `model_bin_phi.stan`, `model_beta_bin_nophi.stan`, `model_beta_bin_phi.stan`: the stan codes to fit the models. Four alternatives. `model_bin_phi.stan` is the main model.
- `run_models.R`: the R codes to read the data and fit the models. Includes description of the data.
- `run_models_cv.R`: the R codes to read the data and fit the models using a 10-fold cross-validation for model comparison. The actual model comparison is done in file `comparison.R`.
- `comparison.R`: the R codes to compare the different model alternatives.
- `barometer.qmd`, `barometer.pdf`: files including the R codes and tables considering Swedish media consumption by age and gender.
