## ETK 2022-10-17, considering example1 from FL 2022-10-16
## FL 2022-12-03/06 various updates

R.dir <- here::here("R")

### load packages
source(file.path(R.dir, "handle_packages.R"))
# Minimum versions:
if (!handle_packages(c(
  "INLA" = "22.11.28",
  "INLAspacetime" = "0.1.2",
  "inlabru" = "2.7.0",
  "fmesher" = NA,
  "ggplot2" = NA,
  "fields" = NA,
  "viridis" = NA,
  "viridisLite" = NA,
  "patchwork" = NA,
  "grid" = NA
))) {
  stop("Packages not fully installed.")
}
library(INLA)
library(inlabru)
library(fmesher)
library(ggplot2)
library(fields)
library(viridis)
library(viridisLite)
library(patchwork)
library(grid)

rerun <- FALSE

## source the computations file
if (rerun) {
  source(file.path(R.dir, "example_computation.R"))
}
data.dir <- here::here("example_data")
fig.dir <- here::here("figures")
dir.create(data.dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig.dir, showWarnings = FALSE, recursive = TRUE)

load(file.path(data.dir, "example.RData"))

### --------------------------------------------------------------
ggplot() +
  gg(mesh.s) +
  coord_equal() +
  theme_bw()

### Automated scales involve the values on the mesh extension, and would be too wide

zlim_abs <- c(-3.7, 3.7)

scale_abs <- list(
  viridis::scale_fill_viridis(limits = zlim_abs, option = "A"),
  coord_equal()
)

mapper <- bru_get_mapper(cmodels[[1]])

pl <- list()
for (time_idx in seq_along(mesh.t$loc)) {
  time <- mesh.t$loc[time_idx]
  ## --------------------------------------------------------------
  pl1 <- ggplot() +
    gg(mesh.s,
      mask = boundary,
      color = ibm_eval(mapper,
        list(space = mesh.s$loc, time = time),
        state = pred.mean[, 1]
      )
    ) +
    scale_abs
  pl2 <- ggplot() +
    gg(mesh.s,
      mask = boundary,
      color = ibm_eval(mapper,
        list(space = mesh.s$loc, time = time),
        state = pred.mean[, 2]
      )
    ) +
    scale_abs
  pl3 <- ggplot() +
    gg(mesh.s,
      mask = boundary,
      color = ibm_eval(mapper,
        list(space = mesh.s$loc, time = time),
        state = pred.mean[, 3]
      )
    ) +
    scale_abs
  pl4 <- ggplot() +
    gg(mesh.s,
      mask = boundary,
      color = ibm_eval(mapper,
        list(space = mesh.s$loc, time = time),
        state = pred.mean[, 4]
      )
    ) +
    scale_abs
  pl[[time_idx]] <- list(pl1, pl2, pl3, pl4)
}

pl_2 <-
  ((pl[[1]][[1]] | pl[[1]][[2]] | pl[[1]][[3]] | pl[[1]][[4]]) /
    (pl[[2]][[1]] | pl[[2]][[2]] | pl[[2]][[3]] | pl[[2]][[4]]) /
    (pl[[3]][[1]] | pl[[3]][[2]] | pl[[3]][[3]] | pl[[3]][[4]])) +
    plot_layout(guides = "collect") &
    labs(fill = "Value") &
    xlab("") & ylab("") &
    theme_bw() +
      theme(legend.position = "right",
            legend.title = element_text(size = 20),
            legend.text = element_text(size = 20),
            legend.key.height = unit(3, "cm"),
      ) +
      theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank()
      )

png(file.path(fig.dir, "example-predictions.png"),
  pointsize = 20,
  units = "in",
  width = 16,
  height = 12,
  res = 300
)
print(pl_2)
grid.text("A", x = unit(0.13, "npc"), y = unit(0.98+0.01, "npc"))
grid.text("B", x = unit(0.365, "npc"), y = unit(0.98+0.01, "npc"))
grid.text("C", x = unit(0.6, "npc"), y = unit(0.98+0.01, "npc"))
grid.text("D", x = unit(0.835, "npc"), y = unit(0.98+0.01, "npc"))
grid.text("t=0", x = unit(0.01, "npc"), y = unit(0.84, "npc"), rot = 90)
grid.text("t=1", x = unit(0.01, "npc"), y = unit(0.51, "npc"), rot = 90)
grid.text("t=2", x = unit(0.01, "npc"), y = unit(0.18, "npc"), rot = 90)
dev.off()

if (FALSE) {
  system(paste0("eog ", file.path(fig.dir, "example-predictions.png"), " &"))
}
