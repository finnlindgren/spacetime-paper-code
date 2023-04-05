
List of R script files.
We have it in the order it is to be executed,
which is also in the ../../Makefile
- get_data.R
- functions.R (helper file only)
- datacheck.R
- datacheck_details_2022.Rmd (show data checking details)
- dataprepare.R
- mesh_creation.R
- data_plots.R
- mesh_plots.R
- modelsfitting.R
- mpredictions7d_only_u.R
- results_tables_plots.R

The get_data.R does
 1. Download data
  stations information
  all the daily data for one year
  elevation data "ETOPO2.RData"
 2. select TMIN and TMAX data with qflag=""
  save wide format to "w2data0.rds" file in "data_files"
 3. load elevation data
 4. load stations and convert to SpatialPointsDataFrame
 5. fix the elevation data and save stations data

The datacheck.R does
1. load the wide shape data from the "w2data0.rds" file
2. select stations with at least 14 observations
3. for each station perform a check for
  - outliers as data outside 7 stdev
  - periods with stdev lower/higher
4. set NA for outliers and periods with too low
  or too high variability.
  See "datacheck_details_2022.Rmd" for details
5. save this wide shape data to "wtavg.rds" in data_files

The dataprepare.R does
1. load the wide shape data from the "wgavg.rds" file
2. define latitude and temporal basis meshes
  and evaluate it for the data
3. collect the stations locations and convert it to
  the spherical coordinates
4. merge it (and elevation)
  to merge it with the data
5. save the long shape dataset to "londata.rds"

mesh_creation.R:
1. Creates spatial and temporal meshes for a given resolution
 and save it to "stmeshes.RData"

The data_plots.R does
1. get world map and define a ocean map
2. load the wide format dataset and the stations information
3. define a grid to group the time series and
 do the grouped time series plot, Figure 4
  - figures/wdata_tsplot.png

The mesh_plots.R
1. get world map and define a ocean map
2. load the "stmeshes.RData"
3. identify the mesh edges that wrap around the earth
4. project the coordinates in Mollweide projection
5. draw the mesh for Figure 4
  - figures/wmesh.pdf

The modelsfitting.R does
1. load the long shape data and meshes
2. setup the model without the spacetime effect v and fit (model M0)
3. fit each one of the models with spatial and spacetime fields
 saving selected elements and into data_files/

The mpredictions7d_only_u.R does
1. set 12 prediction scenarios, one for each month
2. load each model and do a prediction considering
  the hyperparameters fixed as well the
  part of the linear predictor that does not includes u

The results_tables_plots.R does
1. load all the outputs (fitted results and predictions)
2. compute 'goodness-of-fit' statistics, for Table 3
3. transform posteriors and its summary, for Table 4
4. visualize the smoothed latitude temporal mean, for Figure 5
5. visualize the posterior mean for the v field, for Figure 6
6. visualize the posterior mean for the u field, for figure 7
7. visualize the forecast multihorizon errors, for figure 8
8. visualize detailed forecast multi-horizon and scenario errors,
 for the figures in the last appendix

