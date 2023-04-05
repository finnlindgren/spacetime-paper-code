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
wtavg0 <- (w2tavg[,,1] + w2tavg[,,2])/20

### number of observations per station and variable 
wna <- is.na(wtavg0)
nobs.s <- rowSums(!wna)

### select only stations with at least 14 observations
iissel <- which(nobs.s>14)
nssel <- length(iissel)
wtavg <- wtavg0[iissel,]
cat("Consider", sum(nobs.s[iissel]), "observations on", nssel, "stations\n")

### check outliers in the data 
wEpan <- 1-seq(0, 1, length=30)^2
wEpan <- wEpan[2:(length(wEpan)-1)]
wout <- parallel::mclapply(1:nssel, function(i) {
  return(outDetect(wtavg[i,], weights = wEpan, ff = c(5,5)))
}, mc.cores=8L)
nout <- sapply(wout, sum, na.rm=TRUE)
ntot.out <- sum(nout) 
iout <- which(nout>0)
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
  s2 = list(jj = 6:365, nsub = 15, fs = 10))

### compute standard deviation for segments of each time series
ss.results <- parallel::mclapply(ss.setups, function(s) {
  sapply(1:nrow(wtavg), function(i) 
    stdSubs(wtavg[i, s$jj], nsub = s$nsub, fs = s$fs))
}, mc.cores = 4L)

### summarize results
stdsubs <- Reduce("+", ss.results)>0 
istdout <- which(stdsubs)
nstdout <- length(istdout)
cat("Detected", nstdout, "stations with locally low/high variance\n")

### manually define which data to remove
rm.periods <- vector('list', nstdout)

if(year==2021) {

	snames <- "USR0000ACOT"
	i.s.rm <- which(dimnames(wtavg)[[1]][istdout] %in% snames)
	stopifnot(length(i.s.rm)==1)
        rm.periods[[i.s.rm[1]]] <- 262:365
	
}

if(year==2022) {
	
    snames <- c("USC00091732", "USC00411048", "USC00412786", "USC00413507",
                "USC00518108", "USR0000ACOT", "USS0021B31S", "USS0022C12S",
                "USW00000230", "USW00053952")##, "USS0022G24S")
	i.s.rm <- which(dimnames(wtavg)[[1]][istdout] %in% snames)
	
	stopifnot(length(i.s.rm)==10)

    rm.periods[[i.s.rm[1]]] <- 178:243
    rm.periods[[i.s.rm[5]]] <- 170:365
    rm.periods[[i.s.rm[6]]] <- 1:133 
    rm.periods[[i.s.rm[7]]] <- 235:340
    rm.periods[[i.s.rm[8]]] <- 294:331
##	rm.periods[[i.s.rm[5]]] <- 210:365

} 


cat("Manually set to NA data on days\n")
for(i in 1:nstdout) {
    if(length(rm.periods[[i]])>0)
        cat(rm.periods[[i]], "\nof station", dimnames(wtavg)[[1]][istdout[i]], "\n")
}

if(FALSE) {

    par(mfrow=c(3,4), mar=c(2,2,2,.1), mgp=c(2,0.5,0))
    for(i in 1:nstdout) {
        plot(wtavg[istdout[i],], type='l', xlab='', ylab='', axes = FALSE, 
             main = dimnames(wtavg)[[1]][istdout[i]]) 
        points(wtavg[istdout[i],], pch=8, col = 1 + (1:365%in%rm.periods[[i]]))
        axis(1, 5*(0:73), las=2)
        axis(2)
    }

}

### remove the periods
for (i in 1:length(istdout))
  wtavg[istdout[i], rm.periods[[i]]] <- NA

natavg <- is.na(wtavg)
nobs.new <- rowSums(!natavg)
summary(nobs.new)
iissnew <- which(nobs.new>14)
wtavg <- wtavg[iissnew, ]
print(dim(wtavg))

attr(wtavg, "year") <- year 

cat(nrow(wtavg), "stations with a total of", sum(nobs.new), "observations\n")

cat("Save the data to be used in the analysis\n")

### save the checked dataset
saveRDS(
    object = wtavg,
    file = wtavg.rds)

