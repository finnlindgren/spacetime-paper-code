### load the checked dataset
### define basis functions, spacetime and variable id

library(INLA)
library(inlabru)

data.dir <- here::here("data_files")

### output files:
stations.file <- file.path(data.dir, "stations.longlat.rds")
ldata.file <- file.path(data.dir, "longdata.rds")
Bbasis.file <- file.path(data.dir, "B_meshes.RData")

### load the checked data (wide format)
wtavg <- readRDS(file.path(data.dir, "wtavg.rds"))
year <- attr(wtavg, "year")
cat("Data on", nrow(wtavg), "stations loaded\n")

### load the stations information data
all.stations.longlat <- readRDS(
  file.path(data.dir, "all.stations.longlat.rds")
)

### id to stations
id.to.stations <- pmatch(dimnames(wtavg)[[1]], all.stations.longlat$station)

### stations with data
stations.longlat <- all.stations.longlat[id.to.stations, ]
saveRDS(
  object = stations.longlat,
  file = stations.file
)
cat(
  "Saved data on", nrow(stations.longlat),
  "selected stations to", stations.file, "\n"
)

### data information resume
c(
  nstations <- dim(wtavg)[1],
  ntimes <- dim(wtavg)[2]
)

## latitude basis setting
slat.basis.knots <- seq(-1, 1, length = 3)
Blat.mesh <- inla.mesh.1d(
  loc = slat.basis.knots,
  interval = c(-1, 1),
  boundary = "free",
  degree = 2
)

### setting fixed effects temporal basis
time.basis <- bru_mapper_harmonics(
  order = 1, interval = c(1, 366)
)

### save
save(
  list = c("Blat.mesh", "time.basis"),
  file = Bbasis.file
)
cat("Saved temporal latitude mesh/basis to", Bbasis.file, "\n")

### stations longlat coordinates
locs.ll <- sp::coordinates(stations.longlat)

### evaluate latitude basis for the data
Blat.data0 <- as.matrix(inla.spde.make.A(
  Blat.mesh, sin(pi * locs.ll[, 2] / 180)
))
Blat.data <- Blat.data0[
  rep(1:nstations, ntimes),
  c(1, ncol(Blat.data0))
] ## select 1st and last
colnames(Blat.data) <- c("south", "north")

### data evaluation of the time basi
Btime.data0 <- as.matrix(
  ibm_jacobian(
    mapper = time.basis,
    input = 1:ntimes
  )
)

Btime.data <- Btime.data0[ ## (each station first!!!)
  rep(1:ntimes, each = nstations),
]
colnames(Btime.data) <- paste0("t", 1:ncol(Btime.data))

### stations coordinates in the unit sphere
locs.sph <- inla.mesh.map(locs.ll, "longlat", inverse = TRUE)

### start the data.frame
summary(elev0km <- stations.longlat$elevation / 1000)

### long data.frame, consider dim(wtavg) = c(nstations, ndays)
### thus the data will start with all the stations for 1st time
system.time(
  ldata <- data.frame(
    mu = rep(1, nstations * ntimes),
    elev = rep(elev0km, ntimes),
    south = as.matrix(Blat.data[, 1] * Btime.data),
    north = as.matrix(Blat.data[, 2] * Btime.data)
  )
)
print(colnames(ldata))

### add the variable
ldata$y <- as.vector(wtavg)

## add space and time (observe space and time order)
ldata$s1loc <- rep(locs.sph[, 1], ntimes)
ldata$s2loc <- rep(locs.sph[, 2], ntimes)
ldata$s3loc <- rep(locs.sph[, 3], ntimes)
ldata$time <- rep(1:ntimes, each = nstations)

cat("Data structure:\n")
print(str(ldata))
cat("Will use a total of", sum(!is.na(ldata$y)), "non-missing observations!\n")

### save the dataset
saveRDS(object = ldata, file = ldata.file)
cat("Saved longdata into", ldata.file, "\n")
