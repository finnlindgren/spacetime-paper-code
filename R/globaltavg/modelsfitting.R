### code to fit each one of the spacetime models

### libraries
library(INLA)
library(INLAspacetime)
library(inlabru)

data.dir <- here::here("data_files")

### load the data and meshes
ldata <- readRDS(file.path(data.dir, "longdata.rds"))
load(file.path(data.dir, "stmeshes.RData"))
load(file.path(data.dir, "B_meshes.RData"))

ls()

str(ldata)
(ndata <- nrow(ldata))

### define the output files
rfls <- file.path(
  data.dir,
  paste0(
    "tavg_m", 0:4,
    "_n", ndata,
    ".rds"
  )
)

### space time model domain setting
(ntimes <- max(ldata$time))
(nst <- tmesh$n * gmesh$n)

cat("# gmesh nodes:", gmesh$n, "\n",
  "# time knots:", tmesh$n, "\n",
  "st model size:", nst, "\n",
  "data size:", ndata, "\n",
  sep = ""
)

### define the spatial SPDE model for v
sspde <- inla.spde2.pcmatern(
  gmesh,
  alpha = 2,
  prior.range = c(1000 / 6371, 0.01),
  prior.sigma = c(5, 0.01),
  constr = TRUE
) ## not needed, here is the right place to set it

### define the mapper for the spatial model
smapper <- bru_mapper(gmesh)

### set prior for the likelihood parameter
lprec <- list(theta = list(
  prior = "pc.prec",
  param = c(5, 0.01)
))

### the likelihood part
lhood <- like(
  formula = y ~ .,
  family = "gaussian",
  control.family = list(
    hyper = lprec
  ),
  data = ldata
)

### Note:
### inla.mode="compact" is the "new avenue" in INLA (see https://arxiv.org/abs/2204.06797)
### smtp="pardiso" allows two levels of parallelism (see https://arxiv.org/abs/2204.04678)
### pardiso.license : please follow instructions from inla.pardiso()
### num.threads="A:B" : uses A parallel function evaluations each one with B threads
###   (the amount of memory scales approx. linearly in A)
###   (recommended to have A = 2p, or p+1, where p is the number of hyperparameters)
###   (increasing B only makes sense for large problems)
###   (see INLA documentation for more details)
inla.setOption(
  inla.mode = "compact",
  smtp = "pardiso",
  pardiso.license = "~/.pardiso.lic"
)

### define control.compute options
ccomp <- list(
  dic = TRUE, cpo = TRUE, waic = TRUE,
  mlik = TRUE, return.marginals = FALSE
)

### define model components
M0f <- ~ 0 + mu + elev +
  south.t1 + south.t2 + south.t3 +
  north.t1 + north.t2 + north.t3 +
  space(cbind(s1loc, s2loc, s3loc),
    model = sspde, mapper = smapper,
    replicate = time, replicate_mapper = time.basis
  )


### verify the mode components
print(summary(component_list(
  M0f,
  lhoods = list(lhood), verbose = FALSE
)))

### initial values for theta
theta.ini <- c(-2.4, -1, 2)

fit <- bru(
  components = M0f,
  lhood,
  options = list(
    verbose = TRUE,
    control.compute = ccomp,
    control.mode = list(theta = theta.ini, restart = TRUE),
    control.inla = list(int.strategy = "eb")
  )
)

### select elements of the output to save
rnams <- c(
  paste0(
    "summary.",
    c("fixed", "hyperpar", "random", "fitted.values")
  ),
  "internal.marginals.hyperpar",
  "mlik", "mode", "cpu.used"
)

sres <- fit[rnams]
sres$summary.fitted.values <-
  fit$summary.fitted.values[1:ndata, c(1, 2)]
rownames(sres$summary.fitted.values) <- NULL

### R^2 = 1 - ssres/sstot
### sstot = sum((y - mean(y))^2)
sqrtot <- sum((ldata$y - mean(ldata$y, na.rm = TRUE))^2, na.rm = TRUE)
cat("SQR Total:", sqrtot, "\n")

### ssres = sum(y - fit) ^2
sqrres <- sum((ldata$y - sres$summary.fitted.values$mean)^2, na.rm = TRUE)
cat("SQR Res:", sqrres, "\n")

### Some fit statistics
sres$stats <- c(
  R2 = 1 - sqrres / sqrtot,
  INLAspacetime::stats.inla(
    fit,
    y = ldata$y, ### same data used to fit
    fsummarize = function(x) {
      mean(x, na.rm = TRUE)
    }
  )
)

print(round(sres$stats, 4))

### save
saveRDS(
  object = sres,
  file = rfls[1]
)

rm(fit)
gc(reset = TRUE)

### define model components
Mcomps <- ~ 0 + mu + elev +
  south.t1 + south.t2 + south.t3 +
  north.t1 + north.t2 + north.t3 +
  space(cbind(s1loc, s2loc, s3loc),
    model = sspde, mapper = smapper,
    replicate = time, replicate_mapper = time.basis
  ) +
  spacetime(
    list(
      space = cbind(s1loc, s2loc, s3loc),
      time = time
    ),
    model = stmodel
  )

### models id
models <- c("102", "121", "202", "220")

for (imodel in c(1, 2, 3, 4)) {
  ### initial values for theta
  theta.ini <- c(
    -1.45, -1, 1.2,
    c(-1.5, -1, -1.5, -1.5)[imodel],
    c(2.0, 4, 1.5, 2.0)[imodel],
    1
  )

  cat("Fitting model", models[imodel], "\n")
  t0 <- Sys.time()
  rfl <- rfls[1 + imodel]

  ### define the spacetime model
  rt0 <- c(1, 2.0, 1, 2.0) * 1 ## temporal range
  stmodel <- stModel.define(
    gmesh, tmesh, models[imodel],
    control.priors = list(
      prs = c(600 / 6371, 0.01),
      prt = c(rt0[imodel], 0.01),
      psigma = c(5, 0.01)
    ),
    constr = TRUE
  )

  ### print number of non-zeros in Q_u
  cat(
    "Number of non-zeros in Q_u:",
    stmodel$f$cgeneric$data$matrices$xx[2], "\n"
  )

  ### define the mapper for the spacetime model
  stmapper <- bru_get_mapper(stmodel)

  ### verify the mode components
  print(summary(component_list(Mcomps, lhoods = list(lhood), verbose = FALSE)))

  ### set overall INLA options
  ## fit
  fit <- bru(
    components = Mcomps,
    lhood,
    options = list(
      verbose = TRUE,
      num.threads = "3:16", ### with 48 threads
      control.compute = ccomp,
      control.mode = list(
        theta = theta.ini,
        restart = TRUE
      ),
      control.inla = list(int.strategy = "eb")
    )
  )

  ### select elements of the output to save
  rnams <- c(
    paste0(
      "summary.",
      c("fixed", "hyperpar", "random", "fitted.values")
    ),
    "internal.marginals.hyperpar", "mlik", "mode", "cpu.used"
  )

  sres <- fit[rnams]
  sres$summary.fitted.values <-
    fit$summary.fitted.values[1:ndata, c("mean", "sd")]
  rownames(sres$summary.fitted.values) <- NULL ### to use less mem

  ### compute some "goodness-of-fit" statistics
  sqrres <- sum((ldata$y - sres$summary.fitted.values$mean)^2, na.rm = TRUE)
  cat("SQR Res:", sqrres, "\n")

  ### Some fit stats
  sres$stats <- c(
    R2 = 1 - sqrres / sqrtot,
    INLAspacetime::stats.inla(
      fit,
      y = ldata$y, ### same data used to fit
      fsummarize = function(x) {
        mean(x, na.rm = TRUE)
      }
    )
  )
  print(round(sres$stats, 4))

  ### save
  saveRDS(
    object = sres,
    file = rfl
  )

  print(Sys.time() - t0)
}
