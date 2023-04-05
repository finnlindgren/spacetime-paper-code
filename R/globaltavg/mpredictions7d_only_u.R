### code to condition on the available available data
### and on the fixed effect and v(s,t) terms
### for the first 14 days of each month to predict 7 days ahead

### libraries
library(INLA)
library(INLAspacetime)
library(inlabru)

data.dir <- here::here("data_files")

### load the stations, data and meshes
ldata <- readRDS(file.path(data.dir, "longdata.rds"))
(ndata <- nrow(ldata))
(ntimes <- max(ldata$time))

### set prior for the likelihood parameter
lprec <- list(theta=list(prior="pc.prec", param=c(1,0.5)))

### inla options
inla.setOption(
    inla.mode="compact",
    smtp="pardiso",
    pardiso.license="~/.pardiso.lic",
    num.threads="1:-4") 

### define control.compute options
ccomp <- list(dic=TRUE, cpo=TRUE, waic=TRUE,
              mlik=TRUE, return.marginals=FALSE)

stations <- readRDS(file.path(data.dir, "stations.longlat.rds"))

load(file.path(data.dir, "stmeshes.RData"))
load(file.path(data.dir, "B_meshes.RData"))

ls()
c(gmesh$n, tmesh$n)

### space time model domain setting
(nst <- tmesh$n*gmesh$n)

cat("# mesh nodes:", gmesh$n, "\n",
    "# time knots:", tmesh$n, "\n",
    "st model size:", nst, "\n",
    "data size:", ndata, "\n",
    sep = "")

### setting up the indexes for the observations
### to be used in each month/scenario
mm <- sprintf("%02d", 1:12); names(mm) <- mm
mm.times <- lapply(mm, function(m) {
    as.integer(difftime(as.Date(paste0("2022-", m, "-01")) + 0:20,
                        as.Date("2021-12-31"), units = "days"))
})
str(mm.times)
    
### for M0 is the same as the other case
system(paste("cp",
             file.path(data.dir, "tavg_m0_mpred.rds"),
             file.path(data.dir, "tavg_m0_mpred_u.rds")))

### the spacetime mapper for u on the full time period
### to be used to compute u(s,t) at the observations level
stmapper.full <-
    inlabru::bru_mapper_multi(list(
                 space = inlabru::bru_mapper(gmesh),
                 time = inlabru::bru_mapper(tmesh,
                                            indexed = TRUE)))

#################################################################################
#### models A, B, C and D
#################################################################################

inla.setOption(
    num.threads="1:1") 

### models id
models <- c("102", "121", "202", "220")

for(imodel in 1:4) {

    cat("Starting the forecasts for model", models[imodel], "\n")
    tt0 <- Sys.time()

    ## file name for the fitted results (each mode one file)
    rfl <- file.path(
        data.dir,
        paste0(
            "tavg_m",
            imodel,
            "_n",
            ndata,
            ".rds"
        ))

    ## file name for the output (each model one file)
    mpredfl <- file.path(
        data.dir,
        paste0(
            "tavg_m",
            imodel,
            "_mpred_u.rds"
        ))

    ### load the fitted model for the full dataset
    fit <- readRDS(rfl)

    ## select the fitted hypeparameters mode:
    ## likelihood precision and the u(s,t) parameters
    theta.ini <- fit$mode$theta[c(1, 4, 5, 6)]

    tt1 <- Sys.time()
    print(tt1 - tt0)
        
    cat("Doing for each scenario ... \n")

    mpred <- mclapply(mm, function(m) {
        
        ## select the data in the prediction scenario
        ii.m <- which(ldata$time %in% mm.times[[m]])
        sdata <- ldata[ii.m, ]
        
        ## evaluate the posterior mean for the u(s,t) field
        ## at the observation level scenario
        u.fit <- ibm_eval(
            mapper = stmapper.full,
            input = list(space = cbind(
                             sdata$s1loc, sdata$s1loc, sdata$s1loc),
                         time = sdata$time),
            state = fit$summary.random$spacetime$mean)
        
        ## The offset as the sum of the fixed effects and v(s,t)
        offset.mean <- fit$summary.fitted.values$mean[ii.m] - u.fit
        
        sdata$y0 <- sdata$y - offset.mean ### exclude the offset
        iipred <- which(sdata$time %in% tail(sort(unique(sdata$time)), 7))
        sdata$y0[iipred] <- NA

### define the likelihood part
        lhood <- like(
            formula = y0 ~ ., ## y0 for prediction
            family="gaussian",
            control.family = list(
                hyper = lprec),
            data = sdata)
        
        ## define the temporal mesh (in the prediction scenario)
        tmesh.m <- inla.mesh.1d(loc = mm.times[[m]])
        
        ## define the spacetime model (in the prediction scenario)
        rt0 <- c(1,2,1,2)*1 ## temporal range
        stmodel <- stModel.define(
            gmesh, tmesh.m, models[imodel],
            control.priors = list(
                prs = c(600/6371, 0.01),
                prt = c(rt0[imodel], 0.01),
                psigma = c(5, 0.01)),
            constr = TRUE) ### one spacetime constraint
        
        ## print number of non-zeros in Q_u
        cat("Number of non-zeros in Q_u:",
            stmodel$f$cgeneric$data$matrices$xx[2], "\n")
        
        ## define model components: only u(s,t)
        Mcomps.u <- ~ 0 + 
            spacetime(list(space = cbind(s1loc, s2loc, s3loc),
                           time = time),
                      model=stmodel)

        ## evaluate the prediction scenario given the hyperparameters keep fixed
        cat("Evaluating m", m, "...\n")
        pred <- bru(
            components = Mcomps.u,
            lhood,
            options = list(
                verbose = !TRUE,
                control.compute = ccomp,
                control.mode = list(theta = theta.ini, fixed = TRUE),
                control.inla = list(int.strategy = "eb")))
        cat("Selecting results for scenario m", m, "...\n")
        
        ## extract the predicted summary
        pred <- pred$summary.fitted.values[iipred, ]
        rownames(pred) <- NULL
        for(j in c(1, 3:6))
            pred[, j] <- pred[, j] + offset.mean[iipred]
        return(pred) 
        
    }, mc.cores = 4L)
    
    ## save the fitted values summary
    cat("Saving the forecasts for model", models[imodel], "\n")
    
    saveRDS(object = mpred, 
            file = mpredfl)

    cat("Finished the forecasts for model", models[imodel], "\n")
    print(Sys.time() - tt1)
    
}

