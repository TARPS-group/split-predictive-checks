# split-predictive-checks

The ```spc``` package was used to reproduce the experiments for:

*Calibrated Model Criticism Using Split Predictive Checks*. 2022.


The package includes functionality to load data, compute single and divided SPCs for conjugate Poisson model, Gaussian location model, two-level hierarchical Gaussian models, run SPCs for real-data examples and compare performance of SPCs to other methods.


# Usage
To reproduce results for section 5, see [script/simulation_main.r](https://github.com/TARPS-group/split-predictive-checks/blob/main/script/simulation_main.R)

To reproduce plots in section 5, see [script/plot_simulation_results.Rmd](https://github.com/TARPS-group/split-predictive-checks/blob/main/script/plot_simulation_results.Rmd)

For results and plots in section 6.1, see [~/experiments/microarray/ microarray_main.Rmd](https://github.com/TARPS-group/split-predictive-checks/blob/main/script/experiments/microarray/microarray_main.Rmd)

For results and plots in section 6.2 & 6.3, see [~/experiments/airlines/ airlines_main.Rmd](https://github.com/TARPS-group/split-predictive-checks/blob/main/script/experiments/airlines/airlines_main.Rmd)

For basic spc functions, see [spc/](https://github.com/TARPS-group/split-predictive-checks/tree/main/spc) directory.



