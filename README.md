# split-predictive-checks

This ```spc``` repository  was used to reproduce the experiments for:

*Calibrated Model Criticism Using Split Predictive Checks*. 2022.


The package includes functionality to load data, compute single and divided SPCs for conjugate Poisson model, Gaussian location model, two-level hierarchical Gaussian models, run SPCs for real-data examples and compare performance of SPCs to other methods.


# Usage
1. Core method implementation ([spc/](https://github.com/TARPS-group/split-predictive-checks/tree/main/spc) folder)
- compute_spcs.R and compute_spcs_for_hierarchical.R implement the SPCs method.
- generate_pvals.R contains functions for generating p-values for the simulation studies.
- evaluation.R provides visualization functions.
- util.R includes discrepancy functions used in the simulation study.

2. Reproduction scripts ([script/](https://github.com/TARPS-group/split-predictive-checks/blob/main/script/) folder)

To reproduce results for simulation study in Section 5 and Appendix 2:
- The folder [/simulation-models](https://github.com/TARPS-group/split-predictive-checks/tree/main/script/simulation-models) contains all the models used in the simulation study.
- To reproduce all simulation figures and tables, first run [script/simulation_main.r](https://github.com/TARPS-group/split-predictive-checks/blob/main/script/simulation_main.R), then knit
plot_simulation_results.Rmd.

To reproduce results for real-data experiments in Section 6:
- Section 6.1 (Microarray data): run [~/experiments/microarray/ microarray_main.Rmd](https://github.com/TARPS-group/split-predictive-checks/blob/main/script/experiments/microarray/microarray_main.Rmd)
- Section 6.2 (Airline delays data): run [~/experiments/airlines/ airlines_main.Rmd](https://github.com/TARPS-group/split-predictive-checks/blob/main/script/experiments/airlines/airlines_main.Rmd)




