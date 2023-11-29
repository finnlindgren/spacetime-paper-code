### mesh and grouped time series plots

R.dir <- here::here("R", "globaltavg")
data.dir <- here::here("data_files")
figures.dir <- here::here("figures")

year <- 2022

source(here::here("R", "handle_packages.R"))
handle_packages(
  c(
    "INLA" = NA,
    "INLAspacetime" = NA,
    "inlabru" = NA,
    "fmesher" = NA,
    "sp" = NA,
    "sf" = NA,
    "rnaturalearth" = NA
  ),
  attach = TRUE
)

### load the wide format data
wtavg <- readRDS(file.path(data.dir, "wtavg.rds"))

### load stations
stations.longlat <- readRDS(file.path(data.dir, "stations.longlat.rds"))

### 1. get maps needed in the plots

### world map
world.ll <- ne_countries(scale = "medium", returnclass = "sf")
world <- fm_transform(world.ll, "+proj=moll +units=km")

### world box
wbox0 <- cbind(
  long = c(
    seq(-1, 1, length.out = 201), rep(1, 101),
    seq(1, -1, length.out = 201), rep(-1, 101)
  ) * 179.99999,
  lat = c(
    rep(1, 201), seq(1, -1, length.out = 101),
    rep(-1, 201), seq(-1, 1, length.out = 101)
  ) * 89.99999
)
wbox <- sf::st_as_sf(SpatialPolygons(
  list(Polygons(list(Polygon(wbox0)), "0")),
  proj4string = CRS("+proj=longlat +datum=WGS84")
))

### oceans = box - world
use_s2 <- sf::sf_use_s2(FALSE)
oceans.ll <- sf::st_difference(wbox, sf::st_union(world.ll))
sf::sf_use_s2(use_s2)
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
SG <- SpatialGrid(GT, proj4string = fm_CRS(world))
sgroups <- over(stations, SG)

### a spatially grouped plot of some time series
png(file.path(figures.dir, "wdata_tsplot.png"),
  pointsize = 30,
  width = 16,
  heigh = 8.1,
  units = "in",
  res = 300,
  type = "cairo"
)

par(mfrow = c(1, 1), mar = c(0, 0, 1.5, 0), mgp = c(1, .5, 0), xaxs = "i", yaxs = "i", cex = 1.5)
plot(oceans, col = rgb(.3, .7, 1, .5), border = c(0.5, 0.5), xlim = bb[1, ], ylim = bb[2, ])
plot(SG, col = gray(0.5, 0.5), lty = 2, add = TRUE)
plot(world$geometry, col = gray(0.7, 0.5), border = gray(0.5, 0.5), add = TRUE)
system.time(stlines(t(wtavg), SG, sgroups, cex = 0.05, xscale = 3))
title(paste0("Daily temperature at ", nrow(wtavg), " stations for year ", year, "."))

dev.off()
