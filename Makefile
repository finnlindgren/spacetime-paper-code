# Run as
#   make data|meshes|plots_pre|model_fitting|results
# in the root folder.
# To override the Rscript version, use
#   make RSCRIPT=/path/to/Rscript ... 

RSCRIPT?=Rscript

RGT=R/globaltavg

default:

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


.PHONY: default data mesh plot_data plot_mesh plots_pre model_fit model_predict model_fitting results
