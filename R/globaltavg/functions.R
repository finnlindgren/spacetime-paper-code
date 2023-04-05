#' Donload the data files used 
#'
#' @param data.dir the folder to store the files.
#' @param year the year of the daily weather data.
#' @param force logical indicating if it is to force
#' the download. If FALSE each file will be downloaded
#' if it does not exists locally yet.
#' @return a named character vector with the local file names:
#' daily.data, stations.all, elevation.
dataDownload <- function(data.dir, year=2022, force=FALSE)
{
    ## function dataDownload() downloads 
    ## 1. daily weather data for one year
    ## 2. stations information
    ## 3. ETOPO2 elevation data

    
### base URL
    ghcnd <- "https://www.ncei.noaa.gov/pub/data/ghcn/daily/"
    
### daily weather data for a given year
    dfl <- paste0(year, ".csv.gz")
    loc.dfl <- file.path(data.dir, dfl)
    if(force | (!file.exists(loc.dfl)))
        utils::download.file(
                   url = paste0(ghcnd, "by_year/", dfl), 
                   destfile = loc.dfl)
    
### all the available stations information
    sfl <- "ghcnd-stations.txt"
    loc.sfl <- file.path(data.dir, sfl)
    if(force | (!file.exists(loc.sfl)))
        utils::download.file(
                   url = paste0(ghcnd, sfl),
                   destfile = loc.sfl)

### elevation data (the finner grid: ETOPO2)
    efl <- "ETOPO2.RData"
    loc.efl <- file.path(data.dir, efl)
    if(force | (!file.exists(loc.efl))) 
        utils::download.file(
                   url = paste0("http://leesj.sites.oasis.unc.edu/",
                                "FETCH/GRAB/RPACKAGES/", efl), 
                   destfile = loc.efl)

    return(c(daily.data = loc.dfl,
             stations.all = loc.sfl,
             elevation = loc.efl))
    
}
#' Select data from the daily dataset
#'
#' @param gzfile the local filename
#' @param variable string with the variable name(s) to be selected
#' @param qflag a string with quality control flag(s)
#' @param verbose logical indicating if progress is to be printed
#' @section The default selects TMIN, TAVG and TMAX and
#' return it as integer because the original data is also integer
#' with units in 10 Celcius degrees.
#' @references
#' Menne, M., Durre, I., Vose, R., Gleason, B. and Houston, T. (2012)
#' An overview of the global historical climatology network-daily database.
#' Journal of Atmospheric and Oceanic Technology, 897–910.
#' @return array [days, stations, variables] if more than one variable
#' or a matrix [days, stations] if one variable.

dataSelect <- function(gzfile,
                       variable = c("TMIN", "TAVG", "TMAX"),
                       qflag = "", 
                       verbose = TRUE,
                       astype=as.integer)
{

### this function selects `variable` from the daily dataset
### it select data with the given quality control `qfrag`
### it can return the selected data in long or wide format

    if (verbose) 
        t0 <- Sys.time()

### read the full dataset 
    if (requireNamespace("data.table", quietly = TRUE)) {
        d <- data.table::fread(gzfile, data.table = FALSE)
    } else {
        if (verbose) 
            warning("\"data.table\" is not available... it may take a while.")
        d <- utils::read.csv(gzfile)
    }
    
    if (verbose) {
        cat("readed ", nrow(d), "")
        t1 <- Sys.time()
        print(t1-t0)
    }

### select the variables and qflag
    ii <- which(d$V3 %in% variable)
    if (verbose) {
        cat("found ", length(ii), "observations on", variable, "")
        t2 <- Sys.time()
        print(t2-t1)
    }


    ii <- ii[which(d$V6[ii] %in% qflag)]
    d <- d[ii, ]
    
    if (verbose) {
        cat("selected ", length(ii), "observations. ")
        t3 <- Sys.time()
        print(t3-t2)
    }

    cnames <- c("day", "station")
    names(d)[2:1] <- cnames
    if(length(variable)==1) {
        d <- tapply(d[, 4], d[, cnames[2:1]], astype)
    } else {
        cnames <- c(cnames, "variable")
        names(d)[3] <- "variable"
        d <- tapply(d[,4], d[, cnames[c(2,1,3)]], astype)
        d <- d[, , pmatch(variable, dimnames(d)[[3]]), drop = FALSE]
    }
    if (verbose) {
        cat("wide data dim =", dim(d), "")
        t4 <- Sys.time()
        print(t4-t3)
    }
    
    return(d)    
}
#' Detect outliers in a time series considering the raw data
#' and a smoothed version of it.
#'
#' @param x numeric vector
#' @param weights non-increasing numeric vector used as weights for
#' computing a smoothed vector as a rooling window average.
#' Default is null and then \eqn{w_j} is proportional to j
#' in the equation in the Details below.
#' @param ff numeric length two vector with the factors
#' used to consider how many times the standard deviation
#' one data point is out to be considered as an outlier.
#' @return logical vector indicating if the data is an outlier
#' with attributes as detailed bellow.
#' @section attr(, 'm') is the mean of x.
#' @section attr(, 's') is the standard devation of x.
#' @section attr(, 'ss') is the standard deviation for
#' the smoothed data \eqn{y_t} that is defined as
#' \eqn{y_t = \sum_{k=j}^h w_j * (x_{t-j}+x_{t+j})/2}
#' Both `s` and `ss` are used to define outliers if
#'  \eqn{|x_t-m|/s>ff_1} or \eqn{|x_t-y_t|/ss>ff_2}
#' @section attr(, 'xs') the smoothed time series \eqn{y_t}
outDetect <- function(x, weights=NULL, ff = c(7,7))
{
### detect outliers in a time serifes and compute standard deviation for segments
    ## x: n length time series data 
    ## weights: half way weights for the smoothing
    ## ff: factors for detecting outliers (#stdev)
    ## jumps: jump size in the moving window to compute local sd
    n <- length(x)
    if(is.null(weights)) {
        h <- 7
        ## define triangular weights for the time smoothing
        weights <- h:1
    } else {
        stopifnot(all(diff(weights)<=0))
        h <- length(weights)
    }
    weights <- c(rev(weights), 0, weights)
    m <- mean(x, na.rm=TRUE)
    s <- sd(x, na.rm=TRUE)
    x <- x - m
    xx <- c(rep(NA, h), x, rep(NA, h))
    xs <- x*0
    ii <- which(complete.cases(xs))
    if(length(ii)>0) {
        for(i in ii) {### time smoothing
            xw <- weights * xx[-h:h + i + h]
            sw <- sum(weights[complete.cases(xw)])
            xs[i] <- sum(xw, na.rm = TRUE)/sw
        }
    }
    ss <- sd(x - xs, na.rm = TRUE)
### check for outliers in the centered and smoothed data
    r <- (abs(x / s) > ff[1]) | (abs((x - xs) / ss) > ff[2])
    attr(r, "m") <- m
    attr(r, "s") <- s
    attr(r, "ss") <- ss
    attr(r, "xs") <- xs
    return(r)
}
#' To check unusual low/high variance segments
#' @param x numeric vector
#' @param nsub number for the segments length
#' @param fs numeric to use for detecting too
#' hight or too low local standard deviations.
#' @return logical indicating if any of the
#' `st` are `fs` times lower/higher the average
#' of `st`, where is returned as an attribute
#' ad detailed below.
#' @section attr(, 'st') numeric vector with the
#' standard deviation at each segment of the data.
stdSubs <- function(x, nsub=12, fs=15)
{
    n <- length(x)
    ## compute stdev for each segment of the time series
    st <- sapply(split(x, 0:(n-1)%/%nsub), function(xw) {
        if(mean(is.na(xw))>0.5) return(NA)
        return(sd(xw, na.rm=TRUE))
    })
    st.m <- mean(st, na.rm=TRUE)
    r <- any((st/st.m)>fs, na.rm = TRUE) |
        any((st.m/st)>fs, na.rm = TRUE)
    attr(r, 'st') <- st
    return(r)
}
