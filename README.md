# MERS-CoV - Global - 2019

This repo contains the functions necessary to implement the boosted regression tree (BRT) methodology utilized in _Mapping the potential geographic extent of Middle East respiratory syndrome coronavirus_ published in _The Lancet: Global Health_.

1. `data_prep_functions.R`: contains functions to process literature extraction datasheet (or similar) and produce the necessary inputs for model functions

2. `model_functions.R`: contains functions to process files output from `data_prep_functions.R` through the BRT methods resulting in model outputs

3. `optimizers.py`: a python script to optimize the choice of BRT hyperparameters

4. `spillover_functions.R`: contains functions to process model outputs and compute the spillover and detection potential of level 2 administrative units.

5. `get_fit_stats.R`: defines a function for calculating the fit statistics for the bootstrap aggregated final model
