### mesh and grouped time series plots

## R.dir <- here::here("R", "globaltavg")
data.dir <- here::here("data_files")
figures.dir <- here::here("figures")

source(here::here("R", "handle_packages.R"))
handle_packages(
  c(
    INLA = NA,
    INLAspacetime = NA,
    inlabru = NA),
  attach = TRUE)

dir.create(figures.dir, showWarnings = FALSE, recursive = TRUE)

stations.longlat <- readRDS(
  stations.file <- file.path(
    data.dir, "stations.longlat.rds"
  )
)
load(file.path(data.dir, "stmeshes.RData"))

### 1. get maps needed in the plots

### world map
world.ll <- worldMap("+proj=longlat +datum=WGS84", returnclass = "sf")
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


### projection
gmesh$crs <- fm_CRS("sphere")
wmesh <- fm_transform(gmesh, crs = fm_crs(world), passthrough = TRUE)

### transform back to longlat
lmesh <- fm_transform(gmesh, crs = CRS("+proj=longlat +datum=WGS84"), passthrough = TRUE)

### identify the turning edges
i180 <- apply(lmesh$graph$tv, 1, function(ii) {
  any(lmesh$loc[ii, 1] < (-100)) & any(lmesh$loc[ii, 1] > 100)
})

wmfig <- file.path(figures.dir, paste0("wmesh.pdf"))
wmfig

# png(wmfig,
#    pointsize = 20,
#    units = "in",
#    width=16,
#    heigh=8,
#    res=300,
#    type = "cairo")
pdf(wmfig,
  pointsize = 30,
  width = 16,
  height = 8.1
)

par(mfrow = c(1, 1), mar = c(0, 0, 1.5, 0), mgp = c(1, .5, 0), xaxs = "i", yaxs = "i", cex = 1.4)
plot(oceans, col = rgb(.3, .7, 1, .5), border = c(0.5, 0.5), xlim = bb[1, ], ylim = bb[2, ])
plot(world$geometry, col = gray(0.7, 0.5), border = gray(0.5, 0.5), add = TRUE)
points(stations, pch = 19, cex = 0.75, col = "green")
plot(wmesh, add = TRUE, t.sub = which(!i180))
title(paste0(
  "Spatial mesh with ", wmesh$n,
  " nodes, and the ", nrow(stations.longlat), " stations."
))

dev.off()
