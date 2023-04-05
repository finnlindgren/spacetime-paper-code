R.dir <- here::here("R")
# Must run in the directory above the R code.
stopifnot(file.exists(file.path(R.dir, "covariance_computation.R")))

DIM <- 2 # spatial dimension (only 1d and 2d implemented)

# number of points in each direction to evaluate the covariance function in
N <- 64
# use parallel computations?
use_parallel <- TRUE
# number of cores to use for parallel computations
n.cores <- 8
# Where to store the data, to be read by R/covariance_plots.R
data.dir <- here::here("covariance_data")

source(file.path(R.dir, "handle_packages.R"))
if (!handle_packages(c(
  "tidyverse" = NA,
  "fftwtools" = NA,
  if (use_parallel) {
    c(
      "parallel" = NA,
      "doSNOW" = NA
    )
  } else {
    NULL
  }
))) {
  stop("Packages not fully installed.")
}

dir.create(data.dir, showWarnings = FALSE, recursive = TRUE)
source(file.path(R.dir, "S2C.R"))

library(tidyverse)
if (use_parallel) {
  library(parallel)
  library(doSNOW)
}

#' @param w_s the norm of the spatial (angular) wave number vector
#' @param w_t temporal norm of the temporal (angular) frequency
#' @param gamma_v = [gamma_t gamma_s gamma_0]
#' @param alpha_v = [alpha_t alpha_s alpha_e]
#' @param DIM the spatial domain dimension
nonsep_spectrum1d <- function(w_s, w_t, gamma_v, alpha_v, DIM) {
  # alpha = [alpha_t alpha_s alpha_e]
  # gamma = [gamma_t gamma_s gamma_0]

  S_spat <- gamma_v[2]^2 + w_s^2
  S <- ((gamma_v[1]^2 * w_t^2 + S_spat^alpha_v[2])^alpha_v[1]) * S_spat^alpha_v[3]
  S <- 1 / ((2 * pi)^(1 + DIM) * gamma_v[3]^2 * S)

  return(S)
}

# w_s = norm of space-omega
# ht = time-lag
nonsep_space_cross_spectrum <- function(w_s, ht, gamma_v, alpha_v, DIM) {
  # alpha = [alpha_t alpha_s alpha_e]
  # gamma = [gamma_t gamma_s gamma_0]

  S_spat <- gamma_v[2]^2 + w_s^2
  # Time-integral parameter:
  kappa <- S_spat^(alpha_v[2] / 2) / gamma_v[1]
  # Correlation part of time-integral:
  S <- INLA::inla.matern.cov(nu = alpha_v[1] - 1 / 2, kappa = kappa, ht, corr = TRUE)
  # Variance part of time-integral:
  S <- S * gamma(alpha_v[1] - 1 / 2) / gamma(alpha_v[1]) /
    kappa^(2 * alpha_v[1] - 1) / sqrt(4 * pi)
  # Scaling constants and purely spatial part:
  S <- S / ((2 * pi)^(DIM) * gamma_v[3]^2 * S_spat^alpha_v[3] *
    gamma_v[1]^(2 * alpha_v[1]))

  return(S)
}


# Explicit integration for time, FFT for space
# hx must be an increasing regularly spaced vector, starting at 0
nonsep_covar <- function(hx, ht, gamma_v, alpha_v, tol, DIM,
                         expand_factor = 2) {
  # alpha = [alpha_t alpha_s alpha_e]
  # gamma = [gamma_t gamma_s gamma_0]
  stopifnot(hx[1] == 0)
  stopifnot(all(hx >= 0))
  ht <- abs(ht)

  S_fun <- function(omega, ht, gamma_v, alpha_v, DIM) {
    nonsep_space_cross_spectrum(rowSums(omega^2)^0.5,
      ht,
      gamma_v,
      alpha_v,
      DIM = DIM
    )
  }

  calc_dim <- 2^ceiling(log2(length(hx))) * 2 * expand_factor
  dim <- rep(1, DIM) * calc_dim
  h <- rep(hx[2] - hx[1], DIM)
  L <- dim * h
  #  x_ <- make_x(dim, L)
  omega_ <- make_omega(dim, L)
  omega <- as.matrix(expand.grid(omega_))
  #  S <- S_fun(omega = omega, ht, gamma_v, alpha_v, DIM)
  S <- fold_spectrum(
    omega = omega, S_fun = S_fun, h,
    ht = ht, gamma_v = gamma_v, alpha_v = alpha_v,
    DIM = DIM
  )
  C <- fftshift(S2C(S, dim, h), dim)
  if (DIM == 1) {
    C <- C[seq_along(hx)]
  } else if (DIM == 2) {
    C <- C[seq_along(hx), seq_along(hx)]
  } else {
    stop()
  }

  return(C)
}



## Sets smoothness parameters

## List of models in the paper example
model_defn_paper_example <-
  tibble::tribble(
    ~alpha_t, ~alpha_s, ~alpha_e, ~rt.factor, ~Type,
    ## Separable
    ## Model 1 in the paper figure
    1, 0, 2, 1, "A: Separable order 1",
    ## Model 2 in the paper figure
    1, 2, DIM / 2, 1, "B: Critical diffusion",
    ## Model 3 in the paper figure
    2, 0, 2, 1, "C: Separable order 2",
    ## Model 3 in the paper figure
    2, 2, 0, 1, "D: Iterated diffusion"
  )

# Fixed nu-values, varying beta_s
model_defn_fixed_nu <-
  tibble::tribble(
    ~alpha_t, ~alpha_s, ~alpha_e, ~rt.factor, ~Type,
    ## Separable
    ## Model 1 in the paper figure
    1, 0, 1 + DIM / 2, 1, "beta_s=0",
    ## Model 2 in the paper figure
    1, 2, DIM / 2, 1, "beta_s=beta_*",
    ## Fully nonsep
    ## Model 3 in the paper figure
    1 + DIM / 4, 2, 0, 1, "beta_s=1",
    ## Model 4 in the paper figure
    2, 2, 0, 1, "Iterated diffusion"
  )

# Fixed nu-values of higher order, varying beta_s
model_defn_higher_order <-
  tibble::tribble(
    ~nu_t, ~nu_s, ~beta_s, ~rt.factor, ~Type,
    ## Separable
    1, 3, 0, 1, "beta_s=0",
    ## Critical
    1, 3, 3 / (3 + DIM / 2), 1, "beta_s=beta_*",
    ## Fully nonsep
    1, 3, 1, 1, "beta_s=1",
    3, 1, 1, 1, "Special nu=1,1"
  ) %>%
  rowwise() %>%
  mutate(
    alpha_t = nu_t * max(1, beta_s * (nu_s + DIM / 2) / nu_s) + 0.5,
    alpha_s = 1 / nu_t * min(beta_s * (nu_s + DIM / 2), nu_s),
    alpha_e = (nu_s + DIM / 2) * (1 - beta_s)
  ) %>%
  ungroup()


model_defn <- as.data.frame(model_defn_paper_example)


M <- list()

for (model.id in seq_len(nrow(model_defn))) {
  alpha_t <- model_defn[model.id, "alpha_t"]
  alpha_s <- model_defn[model.id, "alpha_s"]
  alpha_e <- model_defn[model.id, "alpha_e"]
  rt.factor <- model_defn[model.id, "rt.factor"]

  alpha_v <- c(alpha_t, alpha_s, alpha_e)

  # decide gamma_s by fixing the spatial range
  alpha <- alpha_e + alpha_s * (alpha_t - 1 / 2)
  nu_spatial <- alpha - DIM / 2
  range_s <- 1
  gamma_s <- sqrt(8 * nu_spatial) / range_s

  nu_time <- min(alpha_t - 1 / 2, nu_spatial / alpha_s)
  print(paste("nu_time = ", nu_time, " is the min of ", alpha_t - 1 / 2, "and", nu_spatial / alpha_s))
  print(paste("smoothness in time and space", nu_time, nu_spatial))

  # decide gamma_t by fixing spatial range to 1 in the separable case
  range_t <- 1 * rt.factor

  ## alpha_t-1/2 matches the interpretation in the paper, relating to
  ## the temporal correlation range of a spatial constant
  gamma_t <- range_t * gamma_s^alpha_s / sqrt(8 * (alpha_t - 1 / 2))

  # decide gamma_0 by fixing variance to one
  sigma2 <- gamma(alpha_t - 1 / 2) * gamma(alpha - DIM / 2) / (gamma(alpha_t) * gamma(alpha) * (4 * pi)^((DIM + 1) / 2) * gamma_t * gamma_s^(2 * (alpha - DIM / 2)))
  gamma_0 <- sqrt(sigma2)
  gamma_v <- c(gamma_t, gamma_s, gamma_0)


  hs <- hx <- seq(from = 0, to = range_s * 1.25, length.out = N)
  ht <- seq(from = 0, to = range_t * 1.25, length.out = N)

  if (use_parallel) {
    cl <- parallel::makeCluster(n.cores)
    doSNOW::registerDoSNOW(cl)
  }

  if (use_parallel) {
    parallel::clusterExport(cl,
      varlist = c("hx", "ht", "gamma_v", "alpha_v", "tol", "DIM"),
      envir = environment()
    )
    par <- foreach(j = 1:N) %dopar% {
      res <- NA
      try({
        res <- nonsep_covar(hx, ht[j], gamma_v, alpha_v, tol, DIM)
        res <- res[, 1]
      })
      return(res)
    }
    C <- matrix(unlist(par), length(hx), length(ht))
  } else {
    C <- matrix(0, length(hx), length(ht))
    for (j in 1:N) {
      cat(".")
      res <- nonsep_covar(hx, ht[j], gamma_v, alpha_v, tol, DIM)
      C[, j] <- res[, 1]
    }
  }

  save.image(
    file = file.path(
      data.dir,
      paste0(
        "image-a",
        paste(round(alpha_v, 2), collapse = "-"), "dim-", DIM, ".Rdata"
      )
    )
  )

  M[[model.id]] <- list()
  M[[model.id]]$C <- C

  if (use_parallel) {
    try({
      parallel::stopCluster(cl)
    })
  }
}

save.image(file = file.path(
  data.dir,
  paste0("image-all-dim-", DIM, ".Rdata")
))
