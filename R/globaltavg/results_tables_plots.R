### tables and plots for the paper

source(here::here("R", "handle_packages.R"))
handle_packages(
  c(
    INLA = NA,
    INLAspacetime = NA,
    inlabru = NA,
    fmesher = NA,
    fields = NA,
    RColorBrewer = NA,
    ggplot2 = NA,
    terra = NA,
    tidyterra = NA,
    tidyverse = NA,
    "ggplot2" = NA,
    "sf" = NA,
    "rnaturalearth" = NA,
    "rnaturalearthdata" = NA,
    patchwork = NA),
  attach = TRUE)

data.dir <- here::here("data_files")
fig.dir <- here::here("figures")

model_names <- c("M0", LETTERS[1:4])
model_names_latex <- paste0("$\\M", c("o", letters[1:4]), "$")

### 2022 data
ldata <- readRDS(file.path(data.dir, "longdata.rds"))
str(ldata)
(ndata <- nrow(ldata))
(ntimes <- max(ldata$time))

### meshes used
load(file.path(data.dir, "stmeshes.RData"))
ls()

c(gmesh$n, tmesh$n)

gmesh$crs <- fm_CRS("sphere")

m_u <- fm_transform(gmesh, crs = fm_crs("globe"))
cat(paste0(
  "Average node distance (km): ",
  paste0(
    "v = ",
    round(
      c(
        u = median((fm_fem(m_u)$va / pi)^0.5 * 2)
      ),
      digits = 3
    ),
    collapse = ", "
  )
), "\n")

mtags <- model_names

### results for each model fitted
results <- vector("list", length = 5L)
for (i in 1:5) {
  cat("Reading results for model", mtags[i], "\n")
  results[[i]] <- readRDS(file.path(
    data.dir,
    paste0("tavg_m", i - 1, "_fit.rds")
  ))
  results[[i]]$mpred <- readRDS(file.path(
    data.dir,
    paste0("tavg_m", i - 1, "_mpred_u.rds")
  )) ## NEW: tavg_m[0-4]_mpred_u.rds
}
names(results) <- model_names

### some outputs
sapply(results, function(x) x$cpu.used)
sapply(results, function(x) x$cpu.used[["Total"]]) / 3600 ## in h
sapply(results, function(x) unname(x$mode$theta))
sapply(results, function(x) x$mlik[[2]])
sapply(results, function(x) x$summary.fixed$mean)

### gof stats
round(gof.tab <- as.data.frame(t(sapply(results, function(r) r$stats))), 5)

df <- t(as.matrix(gof.tab))
df[1:length(as.vector(df))] <- sprintf("%1.4f", df)
colnames(df) <- model_names_latex
rownames(df) <- toupper(rownames(df))
rownames(df)[1] <- "$R^2$"
for (k in seq_len(nrow(df))) {
  if (k == 1) {
    i <- which.max(gof.tab[, k])
  } else {
    i <- which.min(gof.tab[, k])
  }
  df[k, i] <- paste0("\\textbf{", df[k, i], "}")
}
df <- rbind(Model = colnames(df), df)
print(
  xtable::xtable(df, align = "lrrrrr", digits = 5),
  sanitize.text.function = function(x) x,
  include.colnames = FALSE,
  include.rownames = TRUE
)

# df <- as.matrix(gof.tab)
# df[1:length(as.vector(df))] <- sprintf("%1.4f", df)
# rownames(df) <- model_names_latex
# colnames(df) <- toupper(colnames(df))
# colnames(df)[1] <- "$R^2$"
# for (k in seq_len(ncol(df))) {
#   if (k == 1) {
#     i <- which.max(gof.tab[, k])
#   } else {
#     i <- which.min(gof.tab[, k])
#   }
#   df[i, k] <- paste0("\\textbf{", df[i, k], "}")
# }
# df <- cbind(Model = rownames(df), df)
# print(
#   xtable::xtable(df, align = "lrrrrrrrrrr", digits = 5),
#   sanitize.text.function = function(x) x,
#   include.rownames = FALSE
# )




crps.g <- function(y, m, s) {
  md <- y - m
  s / sqrt(pi) - 2 * s * dnorm(md / s) + md * (1 - 2 * pnorm(md / s))
}
scrps.g <- function(y, m, s) {
  md <- y - m
  -0.5 * log(2 * s / sqrt(pi)) - sqrt(pi) *
    (s * dnorm(md / s) - md / 2 + md * pnorm(md / s)) / s
}
stats.f <- function(y, m, v) {
  data.frame(
    ME = y - m,
    MAE = abs(y - m),
    MSE = (y - m)^2,
    DS = (y - m)^2 / v + log(v),
    CRPS = -crps.g(y, m, sqrt(v)),
    SCRPS = -scrps.g(y, m, sqrt(v))
  )
}

mm <- sprintf("%02d", 1:12)
names(mm) <- mm
mm.times <- lapply(mm, function(m) {
  as.integer(difftime(as.Date(paste0("2022-", m, "-01")) + 0:20,
    as.Date("2021-12-31"),
    units = "days"
  ))
})
str(mm.times)

mpred.ts <- do.call("rbind", lapply(1:length(mm), function(m) {
  iidpred <- which(ldata$time %in% tail(mm.times[[m]], 7))
  do.call("rbind", lapply(2:5, function(i) {
    r <- results[[i]]
    ss <- stats::aggregate(
      stats.f(
        y = ldata$y[iidpred],
        m = r$mpred[[m]]$mean,
        v = r$mpred[[m]]$sd^2 + 1 / r$summary.hyperpar$mean[1]
      ),
      ldata[iidpred, "time", drop = FALSE],
      mean,
      na.rm = TRUE
    )
    ns <- ncol(ss) - 1
    df <- data.frame(
      Month = month.abb[m],
      Model = model_names[i],
      Time = rep(ss$time, ns),
      Score = rep(colnames(ss)[-1], each = nrow(ss)),
      Statistic = unlist(ss[, -1])
    )
    df$ahead <- df$Time - min(df$Time) + 1L
    df$Score <- factor(df$Score, names(stats.f(1, 1, 1)))
    return(df)
  }))
}))
mpred.ts$Month <- factor(mpred.ts$Month, month.abb)

head(mpred.ts)
tail(mpred.ts)

gg.mpred <- ggplot(mpred.ts) +
  geom_line(aes(x = ahead, y = Statistic, color = Month, lty = Month)) +
  theme_minimal() +
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
  geom_abline(slope = 0, intercept = 0, linetype = 2) +
  xlab("Days ahead") +
  facet_grid(vars(Score), vars(Model), scales = "free")

gg.mpred.diff <-
  left_join(mpred.ts,
    mpred.ts %>% filter(Model == "D"),
    by = c("Month", "Time", "ahead", "Score"),
    suffix = c("", ".D")
  ) %>%
  filter(!(Score %in% "ME")) %>%
  filter(!(Model %in% "D")) %>%
  ggplot() +
  geom_line(aes(x = ahead, y = Statistic - Statistic.D, color = Month, lty = Month)) +
  theme_minimal() +
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
  geom_abline(slope = 0, intercept = 0, linetype = 2) +
  xlab("Days ahead") +
  ylab("Statistic difference to D") +
  facet_grid(vars(Score), vars(Model), scales = "free")

png(file.path(fig.dir, "forecast_error_months.png"),
  pointsize = 12,
  units = "in",
  width = 6,
  height = 10,
  res = 300,
  type = "cairo"
)

gg.mpred

dev.off()

png(file.path(fig.dir, "forecast_error_months_difference.png"),
  pointsize = 12,
  units = "in",
  width = 6,
  height = 10,
  res = 300,
  type = "cairo"
)

gg.mpred.diff

dev.off()

ggb1.mpred <-
  mpred.ts %>%
  filter(Score %in% c("MAE", "MSE", "DS", "SCRPS")) %>%
  ggplot() +
  geom_boxplot(aes(x = factor(ahead), y = Statistic, fill = Model)) +
  theme_minimal() +
  geom_abline(slope = 0, intercept = 0, linetype = 2) +
  xlab("Days ahead") +
  ylab("Score") +
  facet_wrap(~Score, scales = "free", nrow = 2)

ggb2.mpred <-
  left_join(mpred.ts,
    mpred.ts %>% filter(Model == "D"),
    by = c("Month", "Time", "ahead", "Score"),
    suffix = c("", ".D")
  ) %>%
  filter(Score %in% c("MAE", "MSE", "DS", "SCRPS")) %>%
  ggplot() +
  geom_boxplot(aes(x = factor(ahead), y = Statistic - Statistic.D, fill = Model)) +
  theme_minimal() +
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
  geom_abline(slope = 0, intercept = 0, linetype = 2) +
  xlab("Days ahead") +
  ylab("Score difference to D") +
  facet_wrap(~Score, scales = "free", nrow = 2)

png(file.path(fig.dir, "forecast_error.png"),
  pointsize = 12,
  units = "in",
  width = 12 * 2 / 3,
  height = 5 * 2 / 3,
  res = 300,
  type = "cairo"
)

ggb1.mpred

dev.off()

png(file.path(fig.dir, "forecast_error_difference.png"),
  pointsize = 12,
  units = "in",
  width = 12 * 2 / 3,
  height = 5 * 2 / 3,
  res = 300,
  type = "cairo"
)

ggb2.mpred

dev.off()

### Equator mean and elevation effect
lapply(results, function(r) r$summary.fixed[1:2, 1:2])

### User-interpretable posterior marginals
eR <- c(1, 6371, 1, 6371, 1, 1)
for (i in 1:5) {
  if (i == 1) {
    results[[i]]$marginals.user.hy <- vector("list", 3)
    for (j in 2:3) {
      results[[i]]$marginals.user.hy[[j]] <-
        inla.tmarginal(
          function(x) exp(x) * eR[j],
          results[[i]]$internal.marginals.hyperpar[[j]]
        )
    }
    names(results[[i]]$marginals.user.hy)[2:3] <-
      c("rs_v", "sigma_v")
  } else {
    results[[i]]$marginals.user.hy <- vector("list", 6)
    for (j in 2:(c(3, 6)[(i > 1) + 1])) {
      results[[i]]$marginals.user.hy[[j]] <-
        inla.tmarginal(
          function(x) exp(x) * eR[j],
          results[[i]]$internal.marginals.hyperpar[[j]]
        )
    }
    names(results[[i]]$marginals.user.hy)[2:6] <-
      c("rs_v", "sigma_v", "rs", "rt", "sigma")
  }
  names(results[[i]]$marginals.user.hy)[1] <- "sigma.e"
  results[[i]]$marginals.user.hy[[1]] <- inla.tmarginal(
    function(x) exp(-x / 2), results[[i]]$internal.marginals.hyperpar[[1]]
  )
}

hstabs <- lapply(results, function(r) {
  sapply(r$marginals.user.hy, function(m) {
    unlist(inla.zmarginal(m, TRUE))[c(1, 2, 3, 7)]
  })
})

lapply(hstabs, round, 4)

df <-
  sapply(hstabs[2:5], function(m) {
    paste0(
      sprintf("%2.2f", m[1, ]), " (",
      sprintf("%2.3f", m[2, ]), ")"
    )
  })
colnames(df) <- model_names_latex[-1]
rownames(df) <- paste0(
  "$",
  c(
    "\\sigma_e",
    "r_v",
    "\\sigma_v",
    "r_s",
    "r_t",
    "\\sigma"
  ),
  "$"
)
print(
  xtable::xtable(df, align = "lrrrr"),
  sanitize.text.function = function(x) x
)

sapply(hstabs, function(x) x[1, 1])
sapply(hstabs[2:5], function(x) x[1, ])

xtable::xtable(hstabs[[4]], digits = 4)

rbind(
  "mean(sigma)" = sapply(results[2:5], function(r) {
    round(exp(r$summary.hyperpar$mean[5]), 2)
  }),
  "sd(mean(u(s,t))" = sapply(results[2:5], function(r) {
    sd(r$summary.random$spacetime$mean)
  })
)

iparams <- 1:6
names(iparams) <- c("sigma_e", "r_v", "sigma_v", "rs", "rt", "sigma")

pmdata <- lapply(iparams, function(p) {
  Reduce(
    "rbind",
    lapply(2:5, function(i) {
      data.frame(
        model = LETTERS[i - 1],
        results[[i]]$marginals.user.hy[[p]]
      )
    })
  )
})

ggpars <- lapply(pmdata, function(d) {
  ggplot(d) +
    geom_line(aes(x = x, y = y, group = model, color = model)) +
    ylab("Density") +
    xlab("")
})

if (FALSE) {
  wrap_plots(ggpars)

  ggpars[[1]] + xlab(expression(sigma[e]))

  ggpars[[2]] + xlab("Spatial range of v")

  ggpars[[3]] + xlab("Marginal variance of v")

  ggpars[[4]] + xlab("Spatial range of u")

  ggpars[[5]] + xlab("Temporal range of u")

  ggpars[[6]] + xlab("Marginal variance of u")

  ### Histograms
  par(mfrow = c(2, 2), mar = c(3, 3, 2, 1), mgp = c(1.5, 0.5, 0), las = 1)
  for (imodel in 2:5) {
    hist(results[[imodel]]$summary.random$spacetime$mean, 100,
      xlab = "", ylab = "",
      main = paste0("Posterior mean of u in ", model_names[imodel])
    )
  }
}

### temporal basis for the fixed effect
load(file.path(data.dir, "B_meshes.RData"))

## for plotting time basis functions
Btime0 <- ibm_jacobian(time.basis, 1:ntimes)
timeDate0 <- as.Date("2021-12-31") + 1:ntimes

### for eval and plotting the latitude basis functions
lat0 <- -90:90
lat0sin <- sin(pi * lat0 / 180)
Blat0 <- inla.spde.make.A(Blat.mesh, lat0sin)

### Basis products:
iibb1 <- rep(1:length(lat0), each = length(timeDate0))
iibb2 <- rep(1:length(timeDate0), length(lat0))
BB0 <- cbind(
  Blat0[iibb1, 1] * Btime0[iibb2, ], ### south
  Blat0[iibb1, 4] * Btime0[iibb2, ]
) ### north

lapply(results, function(r) r$summary.fixed$mean)

seasonal.lat.fit <- lapply(results, function(r) {
  if (length(r$summary.fixed$mean) == 7) {
    return(matrix(BB0 %*% r$summary.fixed$mean[2:7], length(timeDate0)))
  }
  return(matrix(BB0 %*% r$summary.fixed$mean[3:8], length(timeDate0)))
})
str(seasonal.lat.fit)

mu0fl <- file.path(fig.dir, "seasonal_latitude.png")
mu0fl
b_rasters <- c(
  rast(t(seasonal.lat.fit[["A"]])[181:1, ]),
  rast(t(seasonal.lat.fit[["B"]])[181:1, ]),
  rast(t(seasonal.lat.fit[["C"]])[181:1, ]),
  rast(t(seasonal.lat.fit[["D"]])[181:1, ])
)
names(b_rasters) <- model_names[-1]
ext(b_rasters) <- c(
  xmin = 1,
  xmax = length(timeDate0),
  ymin = min(lat0),
  ymax = max(lat0)
)


png(mu0fl,
  pointsize = 12,
  units = "in",
  width = 10,
  height = 3,
  res = 600,
  type = "cairo"
)

ggplot() +
  geom_spatraster(data = b_rasters) +
  scale_fill_distiller(
    type = "seq",
    palette = "Blues",
    limits = range(values(b_rasters))
  ) +
  geom_spatraster_contour(data = b_rasters, binwidth = 10) +
  theme_minimal() +
  facet_wrap(~lyr, nrow = 1) +
  xlab("day of year") +
  ylab("latitude") +
  labs(fill = "b(s,t)") +
  guides()

dev.off()

## system(paste("eog", mu0fl, "&"))

bb <- matrix(c(-2, -1, 2, 1) * 9020.048, 2)
bb

apply(bb, 1, diff)
apply(bb, 1, diff) / 100

### define a grid over the full domain
world <- worldMap(returnclass = "sf")
wx0y0 <- list(
  x = seq(bb[1, 1] + 40, bb[1, 2], 80),
  y = seq(bb[2, 1] + 20, bb[2, 2], 80)
)
sapply(wx0y0, length)
wxy0 <- sf::st_as_sf(SpatialPoints(
  as.matrix(
    expand.grid(x = wx0y0$x, y = wx0y0$y)
  ),
  fm_CRS(world)
))

### world box
wbox0 <- cbind(
  long = c(
    seq(-1, 1, length.out = 201), rep(1, 101),
    seq(1, -1, length.out = 201), rep(-1, 101)
  ) * 179.99999,
  lat = c(
    rep(1, 201), seq(1, -1, length.out = 101),
    rep(-1, 201), seq(-1, 1, length.out = 101)
  ) * 89.99999
)
wbox <- SpatialPolygons(
  list(Polygons(list(Polygon(wbox0)), "0")),
  proj4string = CRS("+proj=longlat +datum=WGS84")
)
earth.moll <- fm_transform(wbox, "+proj=moll +units=km")

world_sf <- ne_countries(scale = "medium", returnclass = "sf")
world_sf <- fm_transform(world_sf, "+proj=moll +units=km")

### the range of the spatial coefficients v from each model (slow temporal varying)
sapply(results[2:5], function(r) range(r$summary.random$space$mean))

### wich times to visualize
jjl <- list(
  quarters = round(((seq_len(4) - 1) + 0.5) / 4 * tmesh$n)
)
jjl

for (el in 1:2) {
  ### which pixels fall in the whole Earth or only over land
  if (el == 1) { ### whole Earth
    system.time(ii.wxy0 <-
      unlist(sf::st_intersects(
        sf::st_as_sfc(earth.moll),
        sf::st_as_sfc(wxy0)
      )))
  } else { ### land
    system.time(ii.wxy0 <-
      unlist(sf::st_intersects(
        world_sf,
        sf::st_as_sfc(wxy0)
      )))
  }

  cat(
    "Pixels in", c("whole Earth", "land")[el], ":",
    length(ii.wxy0), "(",
    100 * length(ii.wxy0) / nrow(wxy0), "%)\n"
  )

  ### convert selected pixels to longlat for the projection
  wxy0ll <- fm_transform(
    wxy0[ii.wxy0, ],
    fm_crs(wbox)
  )

  ### projector based on gmesh
  pgridl <- inla.mesh.projector(gmesh, wxy0ll)
  stopifnot(all.equal(length(ii.wxy0), sum(pgridl$proj$A)))

  for (ijj in seq_along(jjl)) {
    jj1 <- jjl[[ijj]]
    print(timeDate0[tmesh$loc[jj1]])

    Btime0 <- ibm_jacobian(time.basis, jj1)
    dim(Btime0)

    ### project selected time to the grid
    v.grid <- lapply(results[-1], function(r) {
      v3m <- inla.mesh.project(
        pgridl, matrix(r$summary.random$space$mean, ncol = 3)
      )
      vv <- kronecker(Btime0[, 1], v3m[, 1])
      if (ncol(Btime0) > 1) {
        for (k in 2:ncol(Btime0)) {
          vv <- vv + kronecker(Btime0[, k], v3m[, k])
        }
      }
      return(matrix(vv, ncol = nrow(Btime0)))
    })
    (n.v <- length(v.grid))
    str(v.grid)

    ##        ab.st <- max(abs(sapply(v.grid, range)))
    ##      ab.st <- max(40, ab.st)
    ab.st <- 31
    cat("Scale range:", ab.st, "\n")

    figfl <- paste0(
      "slow_st", c("Earth", "Land")[el],
      "4m", gmesh$n, "_", tmesh$n,
      "t", paste0(jj1, collapse = "_"), ".png"
    )
    st4mfl <- file.path(fig.dir, figfl)
    cat("Figure file:", st4mfl, "\n")

    subfig <- list()
    for (j in 1:length(jj1)) {
      subfig[[j]] <- list()
      for (i in 1:n.v) {
        mtag <- mtags[1 + i]
        wzz <- matrix(NA, length(wx0y0$x), length(wx0y0$y))
        wzz[ii.wxy0] <- drop(v.grid[[i]][, j])
        cat("E(v|y,", model_names[i+1], ") range: ",
          paste(format(range(wzz, na.rm = TRUE), digits = 4),
            collapse = " to "
          ), "\n",
          sep = ""
        )
        wzz_raster <- rast(
          data.frame(
            x = rep(wx0y0$x, times = length(wx0y0$y)) * 1000,
            y = rep(wx0y0$y, each = length(wx0y0$x)) * 1000,
            u = as.vector(wzz)
          ),
          crs = "+proj=moll"
        )
        subfig[[j]][[i]] <-
          ggplot() +
          geom_spatraster(data = wzz_raster) +
          scale_fill_distiller(
            type = "div",
            palette = "RdBu",
            na.value = "transparent",
            limits = ab.st * c(-1, 1)
          ) +
          geom_sf(data = world_sf, alpha = 0, col = "grey") +
          scale_x_continuous(expand = c(0, 0)) +
          scale_y_continuous(expand = c(0, 0)) +
          theme_minimal() +
          labs(fill = "v(s,t)") +
          ggtitle(paste0(
            mtag,
            ": ",
            format(timeDate0[tmesh$loc[jj1[j]]], "%b %d")
          )) +
          theme(plot.title = element_text(margin = margin(b = -15)))
      }
    }

    png(st4mfl,
      pointsize = 10,
      units = "in",
      width = 2 * 2 * 4,
      height = 2 * 1 * length(jj1) * 1.07,
      res = 300,
      type = "cairo"
    )

    print(wrap_plots(do.call(c, subfig), ncol = 4) +
      plot_layout(guides = "collect") &
      theme(
        legend.position = "bottom",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.key.width = unit(3, "cm"),
        plot.title = element_text(
          size = 12
        )
      ))

    dev.off()
  }
}

### the range of the posterior mean for u from each model
sapply(results[2:5], function(r) range(r$summary.random$spacetime$mean))

for (el in 1:2) {
  ### which pixels fall in the whole Earth or only over land
  if (el == 1) { ### whole Earth
    system.time(ii.wxy0 <-
      unlist(sf::st_intersects(
        sf::st_as_sfc(earth.moll),
        sf::st_as_sfc(wxy0)
      )))
  } else { ### land
    system.time(ii.wxy0 <-
      unlist(sf::st_intersects(
        world_sf,
        sf::st_as_sfc(wxy0)
      )))
  }

  cat(
    "Pixels in", c("whole Earth", "land")[el], ":",
    length(ii.wxy0), "(",
    100 * length(ii.wxy0) / nrow(wxy0), "%)\n"
  )

  ### convert selected pixels to longlat for the projection
  wxy0ll <- fm_transform(
    wxy0[ii.wxy0, ],
    fm_crs(wbox)
  )

  ### projector based on gmesh
  pgridl <- inla.mesh.projector(gmesh, wxy0ll)
  stopifnot(all.equal(length(ii.wxy0), sum(pgridl$proj$A)))

  for (ijj in seq_along(jjl)) {
    jj1 <- jjl[[ijj]]
    print(timeDate0[tmesh$loc[jj1]])

    ### project selected time to the grid
    u.grid <- lapply(results[-1], function(r) {
      inla.mesh.project(
        pgridl, matrix(
          r$summary.random$spacetime$mean,
          gmesh$n
        )[, jj1]
      )
    })
    (n.u <- length(u.grid))

    #        ab.st <- max(abs(sapply(u.grid, range)))
    #       ab.st <- max(40, ab.st)
    ab.st <- 16
    cat("Scale range:", ab.st, "\n")

    figfl <- paste0(
      "st", c("Earth", "Land")[el],
      "4m", gmesh$n, "_", tmesh$n,
      "t", paste0(jj1, collapse = "_"), ".png"
    )
    st4mfl <- file.path(fig.dir, figfl)
    cat("Figure file:", st4mfl, "\n")

    subfig <- list()
    for (j in 1:length(jj1)) {
      subfig[[j]] <- list()
      for (i in 1:n.u) {
        mtag <- mtags[1 + i]
        wzz <- matrix(NA, length(wx0y0$x), length(wx0y0$y))
        wzz[ii.wxy0] <- drop(u.grid[[i]][, j])
        cat("E(u|y,", model_names[i+1], ") range: ",
          paste(format(range(wzz, na.rm = TRUE), digits = 4),
            collapse = " to "
          ), "\n",
          sep = ""
        )
        wzz_raster <- rast(
          data.frame(
            x = rep(wx0y0$x, times = length(wx0y0$y)) * 1000,
            y = rep(wx0y0$y, each = length(wx0y0$x)) * 1000,
            u = as.vector(wzz)
          ),
          crs = "+proj=moll"
        )
        subfig[[j]][[i]] <-
          ggplot() +
          geom_spatraster(data = wzz_raster) +
          scale_fill_distiller(
            type = "div",
            palette = "RdBu",
            na.value = "transparent",
            limits = ab.st * c(-1, 1)
          ) +
          geom_sf(data = world_sf, alpha = 0, col = "grey") +
          scale_x_continuous(expand = c(0, 0)) +
          scale_y_continuous(expand = c(0, 0)) +
          theme_minimal() +
          labs(fill = "u(s,t)") +
          ggtitle(paste0(
            mtag,
            ": ",
            format(timeDate0[tmesh$loc[jj1[j]]], "%b %d")
          )) +
          theme(plot.title = element_text(margin = margin(b = -15)))
      }
    }

    png(st4mfl,
      pointsize = 10,
      units = "in",
      width = 2 * 2 * 4,
      height = 2 * 1 * length(jj1) * 1.1,
      res = 300,
      type = "cairo"
    )

    print(wrap_plots(do.call(c, subfig), ncol = 4) +
      plot_layout(guides = "collect") &
      theme(
        legend.position = "bottom",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.key.width = unit(3, "cm"),
        plot.title = element_text(
          size = 12
        )
      ))

    dev.off()
  }
}

## system(paste("eog", st4mfl, "&"))
