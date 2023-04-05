### mesh and grouped time series plots

R.dir <- here::here("R", "globaltavg")
data.dir <- here::here("data_files")
figures.dir <- here::here("figures")

year <- 2022

library(rgeos)
library(INLA)
library(INLAspacetime)
library(inlabru)
library(sp)
library(maptools)
library(rgeos)

sp::set_evolution_status(2L)

### load the wide format data
wtavg <- readRDS(file.path(data.dir, "wtavg.rds"))

### load stations
stations.longlat <- readRDS(file.path(data.dir, "stations.longlat.rds"))

### 1. get maps needed in the plots

### world map
world.ll <- worldMap("+proj=longlat +datum=WGS84")
world <- fm_transform(world.ll, "+proj=moll +units=km")

### world box
wbox0 <- cbind(
  long = c(
    seq(-1, 1, length = 201), rep(1, 101),
    seq(1, -1, length = 201), rep(-1, 101)
  ) * 179.99999,
  lat = c(
    rep(1, 201), seq(1, -1, length = 101),
    rep(-1, 201), seq(-1, 1, length = 101)
  ) * 89.99999
)
wbox <- SpatialPolygons(
  list(Polygons(list(Polygon(wbox0)), "0")),
  proj4string = CRS("+proj=longlat +datum=WGS84")
)

### oceans = box - world
oceans.ll <- gDifference(wbox, world.ll)
oceans <- fm_transform(oceans.ll, CRS("+proj=moll +units=km"))

### stations locations
stations <- fm_transform(stations.longlat, fm_CRS(world))

### bounding box
bb <- matrix(c(-2, -1, 2, 1) * 9020.048, 2)
bb

#### define group considering a SpatialGrid
dd <- c(20, 10)
pix.s <- apply(bb, 1, diff) / dd
GT <- GridTopology(bb[, 1] + pix.s / 2, pix.s, dd)
SG <- SpatialGrid(GT, world@proj4string)
sgroups <- over(stations, SG)

### a spatially grouped plot of some time series
png(file.path(figures.dir, "wdata_tsplot.png"),
  pointsize = 20,
  width = 16,
  heigh = 8.1,
  units = "in",
  res = 300,
  type = "cairo"
)

par(mfrow = c(1, 1), mar = c(0, 0, 1.5, 0), mgp = c(1, .5, 0), xaxs = "i", yaxs = "i", cex = 1.5)
plot(oceans, col = rgb(.3, .7, 1, .5), border = c(0.5, 0.5), xlim = bb[1, ], ylim = bb[2, ])
plot(SG, col = gray(0.5, 0.5), lty = 2, add = TRUE)
plot(world, col = gray(0.7, 0.5), border = gray(0.5, 0.5), add = TRUE)
system.time(stlines(t(wtavg), SG, sgroups, cex = 0.05, xscale = 3))
title(paste("Daily temperature at", nrow(wtavg), "stations for year 2022."))

dev.off()
