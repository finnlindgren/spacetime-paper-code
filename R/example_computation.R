## ETK 2022-10-17, considering example1 from FL 2022-10-16
## FL updated 2022-12-03

R.dir <- here::here("R")
data.dir <- here::here("example_data")
dir.create(data.dir, showWarnings = FALSE, recursive = TRUE)

source(file.path(R.dir, "handle_packages.R"))
handle_packages(
  c(
    "INLA" = "22.11.28",
    "INLAspacetime" = "0.1.3",
    "inlabru" = "2.7.0",
    "sp" = NA,
    fmesher = NA
  ),
  attach = TRUE
)

## Seed
set.seed(20200119)

### --------------------------------------------------------------
### Range in time
### For nonseparable
### (we use a factor of 1.8 modification when comparing to separable, see paper)
range.t <- 4.0

### Range in space
range.s <- 2.5

### Max edge in spatial mesh
### Small numbers makes algorithm very slow
### Should be less than 1/4 of range.s, preferably much smaller
max.e <- c(0.25, 1)

### Number of timepoints used
## Must be 1 or greater for first order time models
## Must be 3 or greater for second order time models
## The first time point in this example is 0, so the nr of times = t.max + 1
t.max <- 3

### --------------------------------------------------------------
mesh.t <- fm_mesh_1d(seq(0, t.max, by = 1))

### --------------------------------------------------------------
pol <- matrix(
  c(
    0, 0,
    10, 0,
    10 - 1e-6, 10 - 1e-6,
    0, 10
  ),
  nrow = 4, byrow = TRUE
)
boundary <- SpatialPolygons(list(Polygons(list(Polygon(
  pol
)), ID = "Boundary")))

mesh.s <- fm_mesh_2d_inla(boundary = list(boundary), max.edge = max.e)

(ns <- mesh.s$n)

### --------------------------------------------------------------
mco.space <- inla.spde2.pcmatern(
  mesh = mesh.s, prior.range = c(5, .5), prior.sigma = c(1, .5)
)

### --------------------------------------------------------------
Qsep.s <- inla.spde2.precision(spde = mco.space, theta = log(c(range.s, 1)))

### --------------------------------------------------------------
### Gaussian noise
sig.eps <- 0.01
### Seed used in C code
inla.seed <- sample(1E12, 1)
### Sample with INLA
u.sim <- inla.qsample(n = 1, Qsep.s, seed = inla.seed, num.threads = 1)[, 1]
u.sim <- u.sim - mean(u.sim)
sim.noise <- rnorm(length(u.sim), 0, 1) * sig.eps

### st is spacetime index
df2 <- data.frame(
  y = u.sim + sim.noise,
  u.sim = u.sim,
  sim.noise = sim.noise,
  time = 0,
  locx = mesh.s$loc[, 1],
  locy = mesh.s$loc[, 2]
)
summary(df2)

### model size
(nm <- mesh.s$n * mesh.t$n)

####################################################################
##### Models
####################################################################
### The non-separable random effect / random field model definition
### using cgeneric, one for each model, containing needed variables
### using a mesh in space and a mesh in time.
### We fix the temporal range .
models <- c("102", "121", "202", "220")
names(models) <- paste0("M", 1:4)

### define the cgeneric models
cmodels <- vector("list", 4)
names(cmodels) <- names(models)

for (k in c(1, 3)) { ### separable models
  cmodels[[k]] <- stModel.define(
    mesh.s, mesh.t, models[k],
    control.priors = list(
      prs = c(range.s, 0.5),
      prt = c(range.t, 0),
      psigma = c(1, 0.5)
    )
  )
}
for (k in c(2, 4)) { ### non-separable models (with 1.8 factor for temporal range)
  cmodels[[k]] <- stModel.define(
    mesh.s, mesh.t, models[k],
    control.priors = list(
      prs = c(range.s, 0.5),
      prt = c(range.t * 1.8, 0),
      psigma = c(1, 0.5)
    )
  )
}

### Construct mapping between location&time and the kronecker mesh:
### Not needed since INLAspacetime 0.1.2 that has a bru_get_mapper() method
### defined for stModel_cgeneric objects.
# mapper <- bru_mapper_multi(
#     list(space = bru_mapper(mesh.s),
#          time = bru_mapper(mesh.t, indexed = TRUE)))
#
### From INLAspacetime 0.1.2, if we need the mapper ourselves, do
# mapper <- bru_get_mapper(cmodels[[k]])


### object to store the models
M <- list()
for (i in 1:length(models)) {
  M[[i]] <- list(formula = y ~ .)
  M[[i]]$components <-
    ~ -1 + field(list(space = cbind(locx, locy), time = time),
      model = cmodels[[i]],
      mapper = bru_get_mapper(cmodels[[i]])
    )
}


### --------------------------------------------------------------
print(M)


fits <- vector("list", length(models))

if (file.exists("~/.pardiso.lic")) {
  inla.setOption(
    pardiso.license = "~/.pardiso.lic",
    smtp = "pardiso",
    num.threads = "4:-1",
  )
}

### --------------------------------------------------------------
for (i in 1:length(M)) {
  print(paste("Running:  ", i))
  fits[[i]] <- bru(
    M[[i]]$components,
    like(
      formula = M[[i]]$formula,
      family = "gaussian",
      control.family = list(
        hyper = list(
          prec = list(
            initial = -2 * log(sig.eps),
            fixed = TRUE
          )
        )
      ),
      data = df2
    ),
    options = list(
      verbose = TRUE, safe = FALSE,
      control.inla = list(int.strategy = "eb")
    )
  )
  print(round(fits[[i]]$cpu.used[4], 2))
}

(cput <- sapply(fits, function(x) x$cpu.used[[4]]))
(nfn <- sapply(fits, function(x) x$misc$nfunc))
# nfn / cput

sapply(fits, function(x) x$mode$theta)

### Check what we can set the initial values to
for (i in 1:length(M)) {
  print(paste(
    round(
      fits[[i]]$internal.summary.hyperpar$mean, 3
    ),
    collapse = " , "
  ))
}

### --------------------------------------------------------------
### Comparison
### This only works in this specific case, when we fixed the first hyper-parameter
lapply(fits, function(m) {
  data.frame(var = c("Range", "Stdev"), exp(m$summary.hyperpar[c(1, 2), c(4, 3, 5)]))
})


### --------------------------------------------------------------
pred.mean <- sapply(fits, function(m) m$summary.random$field$mean)
summary(pred.mean)

save.image(file.path(data.dir, "example.RData"))
