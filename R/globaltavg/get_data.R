### 1. Download data
###  stations information
###  all the daily data for one year
###  elevation data "ETOPO2.RData"
### 2. select TMIN and TMAX
###  save wide format
### 3. load elevation data
### 4. load stations and convert to SpatialPointsDataFrame
### 5. fix the elevation data and save stations data

R.dir <- here::here("R", "globaltavg")
data.dir <- here::here("data_files")
dir.create(data.dir, showWarnings = FALSE, recursive = TRUE)

### raw wide shape data file name
w2fl <- file.path(data.dir, "w2data0.rds")

### (TMIN + TMAX)/2 raw data in wide format file name
wfl0 <- file.path(data.dir, "wtavg0.rds")

### information on all stations rds file
all.stations.rds <- file.path(data.dir, "all.stations.longlat.rds")

### functions
source(file.path(R.dir, "functions.R"))

### Extend default download timeout
options(timeout = 3600)

### download the raw data files
year <- 2022
data.files <- dataDownload(
  data.dir = data.dir,
  year = year
)

cat("Local files downloaded:", data.files, "\n")

### select temperature variables with qflag=""
if (file.exists(w2fl)) {
  w2data <- readRDS(w2fl)
} else {
  system.time(w2data <- dataSelect(
    gzfile = data.files["daily.data"],
    variable = c("TMIN", "TMAX"),
    qflag = ""
  ))
  attr(w2data, "year") <- year

  ### saving the raw temperature data (integer 10 Celcius degree)
  saveRDS(
    object = w2data,
    file = w2fl
  )
}

print(str(w2data))
cat("Raw", dimnames(w2data)[[3]], "data dimension:", dim(w2data), "\n")

### Compute (TMIN + TMAX)/2 and save
### (consider that the original daa is in 10 Celcius degrees)
wtavg0 <- (w2data[, , 1] + w2data[, , 2]) / 20 ## average TMIN and TMAX
attr(wtavg0, "year") <- attr(w2data, "year")
saveRDS(
  object = wtavg0,
  file = wfl0
)
cat("Data saved to", wfl0, "\n")

### load the ETOPO  data
load(file.path(data.files["elevation"]))
cat("ETOPO data dimension:", dim(ETOPO2), "\n")

### extract the longitude (and fix) and latitude
elon <- attr(ETOPO2, "lon")
elon[elon >= 180] <- 180 - rev(elon[elon >= 180])
elat <- attr(ETOPO2, "lat")

### fix the order of the lines
ETOPO2 <- ETOPO2[, ncol(ETOPO2):1]

### read the informations on stations
all.stations.longlat <- read.fwf(
  file = file.path(data.dir, "ghcnd-stations.txt"),
  widths = diff(c(0, 11, 20, 30, 37, 40, 71, 75, 79, 85)),
  comment.char = ""
)

### set colnames
colnames(all.stations.longlat) <- c(
  "station", "latitude", "longitude",
  "elevation", "state", "name", "gsn", "hcn/crn", "wmo"
)

### convert into SpatialPoints
sp:::coordinates(all.stations.longlat) <- ~ longitude + latitude

### projection
sp::proj4string(all.stations.longlat) <- "+proj=longlat +datum=WGS84"

alocs.ll <- sp::coordinates(all.stations.longlat)

ij <- list(
  i = findInterval(alocs.ll[, 1], c(-180, elon + 1 / 60)),
  j = findInterval(alocs.ll[, 2], elat)
)

etopoll <- sapply(1:nrow(alocs.ll), function(i) ETOPO2[ij$i[i], ij$j[i]])

### index of the elevation data to fix
ii.to.fix <- which(all.stations.longlat$elevation < (-999))

### fix the elevation data
all.stations.longlat$elevation[ii.to.fix] <- etopoll[ii.to.fix]

### save stations.longlat
saveRDS(
  object = all.stations.longlat,
  file = all.stations.rds
)
cat("Saved data on", nrow(all.stations.longlat), "stations\n")
