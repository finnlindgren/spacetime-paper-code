## Partly based on old code from
##   https://haakonbakkagit.github.io/btopic132.html
## and updated to use newer rgeneric code, and inlabru
## 2022-10-17 Finn Lindgren

overwrite_plots <- FALSE

R.dir <- here::here("R")
data.dir <- here::here("example1.data")
figures.dir <- here::here("figures")
if (overwrite_plots) {
  dir.create(figures.dir, showWarnings = FALSE, recursive = TRUE)
}


## ---- warning=FALSE, message=FALSE-----------------------------
library(INLA)
library(inlabru)
library(sp)
library(ggplot2)
library(fields)
library(viridisLite)
library(patchwork)

data_filename <- file.path(data.dir, "example1.RData")
if (!file.exists(data_filename)) {
  stop(paste0("Data file ", data_filename, " is missing. Run R/example1_computation.R first."))
}
load(file = data_filename)

ibm_eval <- function(mapper, input, state) {
  A <- ibm_amatrix(mapper, input)
  A %*% state
}

## --------------------------------------------------------------
pred.sep <- fits[[1]]$summary.random$field$mean
pred.nonsep <- fits[[2]]$summary.random$field$mean


# The automated scales involve the values on the mesh extension, and is too wide
# zlim_abs <- range(c(pred.sep, pred.nonsep))
# zlim_diff <- c(-1, 1) * max(abs(c(pred.nonsep - pred.sep)))
zlim_abs <- c(-3.3, 3.3)
zlim_diff <- c(-0.7, 0.7)
scale_abs <- list(
  scale_fill_viridis(limits = zlim_abs, option = "A"),
  coord_equal()
)
scale_diff <- list(
  scale_fill_viridis(limits = zlim_diff, option = "A"),
  coord_equal()
)
pl <- list()
for (time in seq_len(3)) {
  ## --------------------------------------------------------------
  pl1 <- ggplot() +
    gg(mesh.s,
       mask = boundary,
       color = ibm_eval(mapper,
                        list(space = mesh.s$loc, time = time),
                        state = pred.sep
       )
    ) +
    scale_abs +
    ggtitle(paste0("DEMF(1,0,2), time = ", time))
  pl2 <- ggplot() +
    gg(mesh.s,
       mask = boundary,
       color = ibm_eval(mapper,
                        list(space = mesh.s$loc, time = time),
                        state = pred.nonsep
       )
    ) +
    scale_abs +
    ggtitle(paste0("DEMF(1,2,1), time = ", time))

  ## --------------------------------------------------------------
  ## Difference between predictions.
  pl3 <- ggplot() +
    gg(mesh.s,
       mask = boundary,
       color = ibm_eval(mapper,
                        list(space = mesh.s$loc, time = time),
                        state = pred.nonsep - pred.sep
       )
    ) +
    scale_diff +
    ggtitle(paste0("Difference, time = ", time))

  pl[[time]] <- list(pl1, pl2, pl3)
}

pl_1 <-
  ((pl[[1]][[1]] | pl[[1]][[2]]) /
     (pl[[2]][[1]] | pl[[2]][[2]]) /
     (pl[[3]][[1]] | pl[[3]][[2]])) +
  plot_layout(guides = "collect") &
  labs(fill = "Value") & xlab("") & ylab("") & theme_bw()

if (!overwrite_plots) {
  print(pl_1)
} else {
  png(file.path(figures.dir, "example-1.png"),
      pointsize = 20,
      units = "in",
      width = 6.5,
      height = 8,
      res = 600
  )
  print(pl_1)
  dev.off()
}
