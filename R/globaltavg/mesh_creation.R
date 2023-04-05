### The code to create the temporal and spatial meshes

### libraries
library(INLA)

### data directory
data.dir <- here::here("data_files")

### The meshes file name
meshes.filename <- file.path(data.dir, "stmeshes.RData")

### load stations and data
ldata <- readRDS(file.path(data.dir, "longdata.rds"))
stations.longlat <- readRDS(file.path(data.dir, "stations.longlat.rds"))

ls()
ndata <- nrow(ldata)

### The spacetime model resolution setup
### The run time and the model fit quality depends on
### the number of nodes over the spacetime domain

### Seting the temporal model resolution
h.t <- 1 ## higher defines less time knots

### Setting the spatial mesh parameters
## globe: smaller gives less mesh nodes
## cutoff: smaller gives more mesh nodes
sresol <- c(7, 0.05)

### resulting approximate edge length over the oceans:
0.2 / sresol[1]

### resulting approximate edge length (in km) over the oceans:
36080.2 * 0.2 / sresol[1]

### cutoff value in km
6371 * sresol[2]

### space time model domain setting
(ntimes <- max(ldata$time))
time0.h <- seq(1, ntimes + h.t * .99, h.t)
cat("time mesh domain:", range(time0.h), "\n")
tmesh <- inla.mesh.1d(time0.h)

### transform data locations to sphere
locs.sph <- inla.mesh.map(
  loc = coordinates(stations.longlat),
  projection = "longlat",
  inverse = TRUE
)

### mesh for the spacetime fields
gmesh <- inla.mesh.create(
  locs.sph,
  globe = sresol[1],
  cutoff = sresol[2],
  refine = list(
    max.edge = 20 * sresol[2]
  )
)
gmesh$n

### save the meshes for later use (e.g. in modelsfitting and outplot)
save(
  list = c("gmesh", "tmesh", "h.t", "sresol"),
  file = meshes.filename
)

(nst <- tmesh$n * gmesh$n)

cat("# gmesh nodes:  ", gmesh$n, "\n",
  "# time knots:   ", tmesh$n, "\n",
  "st model size:  ", nst, "\n",
  "ocean res (km): ", paste(round(36080.2 * 0.2 / sresol[1], 2), collapse = ", "), "\n",
  "obs res (km):   ", paste(round(6371 * sresol[2], 2), collapse = ", "), "\n",
  "data size:      ", ndata, "\n",
  sep = ""
)
