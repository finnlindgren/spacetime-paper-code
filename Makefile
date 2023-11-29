# For the global temperature case study, run as
#   make pkgcheck|data|mesh|plots_pre|model_fitting|results
# in the root folder.
#
# For the covariance plots, run as
#   make covariance_computation|covariance_plots
# in the root folder.
#
# For the simulation example, run as
#   make example_computation|example_plots
# in the root folder.
#
# To override the Rscript version, use
#   make RSCRIPT=/path/to/Rscript ...

RSCRIPT?=Rscript

RGT=R/globaltavg

default:
	@grep "^#" Makefile

pkgcheck:

data:
	$(RSCRIPT) --vanilla $(RGT)/get_data.R
	$(RSCRIPT) --vanilla $(RGT)/datacheck.R
	$(RSCRIPT) --vanilla $(RGT)/dataprepare.R

mesh:
	$(RSCRIPT) --vanilla $(RGT)/mesh_creation.R

plot_data:
	$(RSCRIPT) --vanilla $(RGT)/data_plots.R
plot_mesh:
	$(RSCRIPT) --vanilla $(RGT)/mesh_plots.R
plots_pre: plot_data plot_mesh

model_fit:
	$(RSCRIPT) --vanilla $(RGT)/modelsfitting.R
model_predict:
	$(RSCRIPT) --vanilla $(RGT)/mpredictions7d_only_u.R
model_fitting: model_fit model_predict

results:
	$(RSCRIPT) --vanilla $(RGT)/results_tables_plots.R


covariance_computation_1: covariance_data/image-all-dim-1.Rdata
covariance_data/image-all-dim-1.Rdata:
	$(RSCRIPT) -e 'DIM <- 1; source("R/covariance_computation.R")'
covariance_plots_1: covariance_data/image-all-dim-1.Rdata
	$(RSCRIPT) -e 'DIM <- 1; source("R/covariance_plots.R")'

covariance_computation_2: covariance_data/image-all-dim-2.Rdata
covariance_data/image-all-dim-2.Rdata:
	$(RSCRIPT) -e 'DIM <- 2; source("R/covariance_computation.R")'
covariance_plots_2: covariance_data/image-all-dim-2.Rdata
	$(RSCRIPT) -e 'DIM <- 2; source("R/covariance_plots.R")'

example_computation: example_data/example.RData
example_data/example.RData:
	$(RSCRIPT) -e 'source("R/example_computation.R")'
example_plots: example_data/example.RData
	$(RSCRIPT) -e 'source("R/example_plots.R")'

.PHONY: default data mesh plot_data plot_mesh plots_pre model_fit \
	model_predict model_fitting results \
	covariance_computation_1 covariance_computation_2 \
	covariance_plots_1 covariance_plots_2 \
	example_computation example_plots
