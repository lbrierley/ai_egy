# Spatial models of avian influenza in Egyptian wild birds

Supporting data and code associated with Nabil et al. (2026). Molecular surveillance and predictive risk modelling of avian influenza virus in wild birds in Egypt. PLOS Pathogens. 2026;107(6):002278. https://doi.org/10.1099/jgv.0.002278

## Scripts

`data_egy.R` loads and formats all geolocated avian influenza records for Egypt from the above citation, FAO's EMPRES-i, WOAH's WAHIS, and BV-BRC.
`data_sdm.R` loads and formats all necessary spatial covariate layers. This requires functions from `functions.R`.

`fit_bart.R` and `fit_bart_LOSO.R` will fit, choose, and extract variable importance from a BART machine learning to predict avian influenza presence from the spatial covariates, validating models by either a 75:25 train:test split in spatial site clusters, or leaving each site cluster out in turn (LOSO). These require functions from `bart_scripts/bart_functions.R`.

`bart_scripts/predict_bart` will generate a layer of predicted risk for Egypt from a given model and for a given grid-cell resolution and calendar date(s).

All scripts in `sdm scripts/` handle raw calculation and processing (and sometimes extraction/download) of each used spatial covariate layer. These require functions from `functions.R`.
