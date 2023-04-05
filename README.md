# Example code for non-separable space-time models

This repository contains code for the examples in https://arxiv.org/abs/2006.04917, v4

* `R/handle_packages.R` Helper function for checking whether needed packages are installed
* `R/S2C.R` Helper functions for using FFT to convert between spectral densities and covariances
* `R/covariance_computation.R` Compute covariances. Can be run as `source("R/covariance_computation.R")`
* `R/covariance_plots.R` Make covariance plots, from data saved by the computation code. Can be run as `source("R/covariance_plots.R")`
* `R/example_computation.R` Compute prediction example data. Can be run as `source("R/example_computation.R")`
* `R/example_plots.R` Make prediction plot, from data saved by the example computation code. Can be run as `source("R/example_plots.R")`
* `R/globaltavg` Folder with multiple files for the global temperature example. See
  [`R/globaltavg/README.txt`](R/globaltavg/README.txt)
* `Makefile` Makefile for the global temperature example. The full analysis
  can be run in a terminal by the following sequence:
```
# Download, check and prepare the data:
make data
# Create the analysis mesh:
make mesh
# Plot the data and mesh (runs plot_data and plot_mesh):
make plot_pre
# Fit the models (runs model_fit and model_predict):
make model_fitting
# Produce tables and plots:
make results
```

Main location
https://github.com/finnlindgren/spacetime-paper-code
