### check against outlier and
### locally lower/hight variability

R.dir <- here::here("R", "globaltavg")
data.dir <- here::here("data_files")

### output file name
wtavg.rds <- file.path(data.dir, "wtavg.rds")

### functions
source(file.path(R.dir, "functions.R"))

### the raw data
w2tavg <- readRDS(file.path(data.dir, "w2data0.rds"))
year <- attr(w2tavg, "year")

### the (TMIN + TMAX)/2 wide format data
wtavg0 <- (w2tavg[, , 1] + w2tavg[, , 2]) / 20

### number of observations per station and variable
wna <- is.na(wtavg0)
nobs.s <- rowSums(!wna)

### select only stations with at least 14 observations
iissel <- which(nobs.s > 14)
nssel <- length(iissel)
wtavg <- wtavg0[iissel, ]
cat("Consider", sum(nobs.s[iissel]), "observations on", nssel, "stations\n")

### check outliers in the data
wEpan <- 1 - seq(0, 1, length.out = 30)^2
wEpan <- wEpan[2:(length(wEpan) - 1)]
wout <- parallel::mclapply(1:nssel, function(i) {
  return(outDetect(wtavg[i, ], weights = wEpan, ff = c(5, 5)))
}, mc.cores = 8L)
nout <- sapply(wout, sum, na.rm = TRUE)
ntot.out <- sum(nout)
iout <- which(nout > 0)
cat("Detected", ntot.out, "outliers on", length(iout), "stations\n")

c(ntot.out, length(iout))

### set NA to detected outliers
for (i in iout) {
  wtavg[i, which(wout[[i]])] <- NA
}

### Define subsegments
ss.setups <- list(
  s1 = list(jj = 1:357, nsub = 21, fs = 10),
  s2 = list(jj = 9:365, nsub = 21, fs = 10),
  s2 = list(jj = 1:360, nsub = 15, fs = 10),
  s2 = list(jj = 6:365, nsub = 15, fs = 10)
)

### compute standard deviation for segments of each time series
ss.results <- parallel::mclapply(ss.setups, function(s) {
  sapply(1:nrow(wtavg), function(i) {
    stdSubs(wtavg[i, s$jj], nsub = s$nsub, fs = s$fs)
  })
}, mc.cores = 4L)

### summarize results
stdsubs <- Reduce("+", ss.results) > 0
istdout <- which(stdsubs)
nstdout <- length(istdout)
cat("Detected", nstdout, "stations with locally low/high variance\n")

### manually define which data to remove
rm.periods <- vector("list", nstdout)

if (year == 2021) {
  removal <-
    list(
      "USR0000ACOT" = 262:365
    )
}

if (year == 2022) {
  removal <-
    list(
      "USC00091732" = 178:243,
      "USC00411048" = integer(0),
      "USC00412786" = integer(0),
      "USC00413507" = integer(0),
      "USC00518108" = 170:365,
      "USR0000ACOT" = 1:133,
      "USS0021B31S" = 235:340,
      "USS0022C12S" = 294:331,
      "USW00000230" = integer(0),
      "USW00053952" = integer(0)
    )
}

snames <- names(removal)
for (idx in seq_along(istdout)) {
  nm <- dimnames(wtavg)[[1]][istdout[idx]]
  if (nm %in% snames) {
    rm.periods[[idx]] <- removal[[nm]]
  }
}


cat("Manually set to NA data on days\n")
for (i in seq_len(nstdout)) {
  if (length(rm.periods[[i]]) > 0) {
    cat(rm.periods[[i]], "\nof station", dimnames(wtavg)[[1]][istdout[i]], "\n")
  }
}

if (FALSE) {
  par(mfrow = c(3, 4), mar = c(2, 2, 2, .1), mgp = c(2, 0.5, 0))
  for (i in 1:nstdout) {
    plot(wtavg[istdout[i], ],
      type = "l", xlab = "", ylab = "", axes = FALSE,
      main = dimnames(wtavg)[[1]][istdout[i]]
    )
    points(wtavg[istdout[i], ], pch = 8, col = 1 + (1:365 %in% rm.periods[[i]]))
    axis(1, 5 * (0:73), las = 2)
    axis(2)
  }
}

### remove the periods
for (i in 1:length(istdout)) {
  wtavg[istdout[i], rm.periods[[i]]] <- NA
}

natavg <- is.na(wtavg)
nobs.new <- rowSums(!natavg)
summary(nobs.new)
iissnew <- which(nobs.new > 14)
wtavg <- wtavg[iissnew, ]
print(dim(wtavg))

attr(wtavg, "year") <- year

cat(nrow(wtavg), "stations with a total of", sum(nobs.new), "observations\n")

cat("Save the data to be used in the analysis\n")

### save the checked dataset
saveRDS(
  object = wtavg,
  file = wtavg.rds
)
