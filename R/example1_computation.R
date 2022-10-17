## FL 2022-10-16:
## Code wasn't in repo, so downloaded the version from
##   https://haakonbakkagit.github.io/btopic132.html
## and updated to use the new rgeneric code, and inlabru

R.dir <- here::here("R")
data.dir <- here::here("example1.data")
dir.create(data.dir, showWarnings = FALSE, recursive = TRUE)

## ---- warning=FALSE, message=FALSE-----------------------------
library(INLA)
library(inlabru)
library(sp)
#library(ggplot2)
#library(fields)
#library(viridisLite)
## library(gtools)

## Load functions
source(file.path(R.dir, "rgeneric-models.R"))
## Seed
set.seed(20200119)


## --------------------------------------------------------------
## Range in time
## range.t = 4 for DEMF(1,0,2)
## range.t = 4*1.8 fo DEMF(1,2,1) for effective comparable range
range.t <- 4

## Range in space
range.s <- 2.5

## Max edge in spatial mesh
## Small numbers makes algorithm very slow
## Should be less than 1/4 of range.s, preferably much smaller
max.e <- c(0.25, 1)

## Number of timepoints used
## Must be 2 or greater for first order time models
## Must be 4 or greater for second order time models
t.max <- 3

## --------------------------------------------------------------
mesh.t <- inla.mesh.1d(seq_len(t.max))

## --------------------------------------------------------------
## Slightly nudge one corner to break the symmetry.
boundary <- SpatialPolygons(list(Polygons(list(Polygon(
  matrix(c(
    0, 0,
    10, 0,
    10 - 1e-6, 10 - 1e-6,
    0, 10
  ),
  nrow = 4, byrow = TRUE
  )
)), ID = "Boundary")))
mesh.s <- inla.mesh.2d(boundary = list(boundary), max.edge = max.e)

## --------------------------------------------------------------
ggplot() +
  gg(mesh.s) +
  coord_equal() +
  theme_bw()

## --------------------------------------------------------------
mco.space <- inla.spde2.pcmatern(
  mesh = mesh.s, prior.range = c(5, .5), prior.sigma = c(1, .5)
)


## --------------------------------------------------------------
Qsep.s <- inla.spde2.precision(spde = mco.space, theta = log(c(range.s, 1)))


## --------------------------------------------------------------
## Gaussian noise
sig.eps <- 0.01
## Seed used in C code
inla.seed <- sample(1E12, 1)
## Sample with INLA
u.sim <- inla.qsample(n = 1, Qsep.s, seed = inla.seed, num.threads = 1)[, 1]
u.sim <- u.sim - mean(u.sim) # Artificially centre the simulation
sim.noise <- rnorm(length(u.sim), 0, 1) * sig.eps

df1 <- data.frame(
  y = u.sim + sim.noise,
  u.sim = u.sim,
  sim.noise = sim.noise,
  year = 1,
  locx = mesh.s$loc[, 1],
  locy = mesh.s$loc[, 2]
)
summary(df1)


## --------------------------------------------------------------
## Rgeneric object containing needed variables
## Mesh in space and time
## Lambdas for exponential prior on transformed hyper-param (1/rt, 1/rs and sig)
rgen.obj <- list(
  mesh.space = mesh.s,
  mesh.time = mesh.t,
  lambdas = c(1, 1, 1)
)

# Construct mapping between location&time and the kronecker mesh:
mapper <- bru_mapper_multi(list(
  space = bru_mapper(mesh.s),
  time = bru_mapper(mesh.t, indexed = TRUE)
))

## Nonsep model definition
nm <- mesh.s$n * mesh.t$n

## The non-separable random effect / random field
## We use the function loaded in the beginning of the document
mco.nonsep <- inla.rgeneric.define(
  model = inla.stmodel121, debug = FALSE, n = nm, obj = rgen.obj
)


## --------------------------------------------------------------
ggplot() +
  gg(mesh.s, color = df1$u.sim, mask = boundary) +
  coord_equal()



## --------------------------------------------------------------
M <- list()
for (i in 1:2) M[[i]] <- list()


## --------------------------------------------------------------
## We want to fix the autocorrelation in time
## The theta2 as defined before is the INLA internal parametrisation of rho
## 0.13 = rho^range.t
## rho = 0.13^(1/range.t)
## theta2 = qlogis((1+rho)/(1-rho))
rho.t <- 0.13^(1 / range.t)
hyper.ar1.rho <- list(rho = list(initial = qlogis((1 + rho.t) / (1 - rho.t)), fixed = TRUE))
comp.sep <- ~ -1 + field(cbind(locx, locy),
  model = mco.space,
  mapper = bru_mapper(mesh.s),
  group = year,
  group_mapper = bru_mapper(mesh.t, indexed = TRUE),
  control.group = list(model = "ar1", hyper = hyper.ar1.rho)
)
M[[1]]$components <- comp.sep
M[[1]]$formula <- y ~ .


## --------------------------------------------------------------
## We need to fix the temporal range in the non-separable model
rgen.obj2 <- rgen.obj
rgen.obj2$fixed.theta <- c(log(range.t * 1.8), NA, NA)
mco.nonsep.fix <- inla.rgeneric.define(model = inla.stmodel121, debug = FALSE, n = nm, obj = rgen.obj2)

comp.nonsep <- ~ -1 + field(list(space = cbind(locx, locy), time = year),
  model = mco.nonsep.fix,
  mapper = mapper
)
M[[2]]$components <- comp.nonsep
M[[2]]$formula <- y ~ .


## --------------------------------------------------------------
print(M)


## --------------------------------------------------------------
## WARNING: Set these variables to NULL if you change the model in any way!
M[[1]]$init <- c(1.8, 0.1)
M[[2]]$init <- c(1.142, -0.035)

fits <- list()


## --------------------------------------------------------------
for (i in 1:length(M)) {
  print(paste("Running:  ", i))
  fits[[i]] <- bru(M[[i]]$components,
    like(
      formula = M[[i]]$formula,
      family = "gaussian",
      control.family = list(hyper = list(prec = list(
        initial = -2 * log(sig.eps), fixed = TRUE
      ))),
      data = df1
    ),
    options = list(
      control.predictor = list(compute = TRUE),
      verbose = TRUE,
      num.threads = 4,
      control.inla = list(int.strategy = "eb"),
      control.mode = list(restart = TRUE, theta = M[[i]]$init),
      control.compute = list(
        dic = TRUE, cpo = TRUE, waic = TRUE,
        mlik = TRUE, return.marginals = FALSE,
        openmp.strategy = "default", smtp = "taucs"
      )
    )
  )

  print(round(fits[[i]]$cpu.used[4], 2))
}


# --------------------------------------------------------------
## Check what we can set the initial values to
for (i in 1:length(M)) {
  print(paste(round(
    fits[[i]]$internal.summary.hyperpar$mean, 3
  ),
  collapse = " , "
  ))
}


## --------------------------------------------------------------
## Comparison
fits[[1]]$summary.hyperpar[, c(4, 3, 5)]

## This only works in this specific case, when we fixed the first hyper-parameter
data.frame(var = c("Range", "Stdev"), exp(fits[[2]]$summary.hyperpar[c(1, 2), c(4, 3, 5)]))


save.image(file.path(data.dir, "example1.RData"))
