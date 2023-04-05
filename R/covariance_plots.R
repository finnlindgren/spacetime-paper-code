# Must run in the directory above the R code.
stopifnot(file.exists("R/covariance_plots.R"))

# Set TRUE to overwrite plots on disk
overwrite_plots <- FALSE
# Set to 1 or 2. The example in the paper is for DIM=2
DIM <- 2

# Where to read the computed covariances, from R/covariance_computation.R
data.dir <- here::here("covariance_data")
# Where to save the figures, if overwrite_plots is TRUE
figures.dir <- here::here("figures")

source("R/S2C.R")

library(tidyverse)
library(ggplot2)
library(patchwork)
library(metR) # For geom_contour_fill

data_filename <- file.path(data.dir, paste0("image-all-dim-", DIM, ".Rdata"))
if (!file.exists(data_filename)) {
  stop(paste0("Data file ", data_filename, " is missing. Run R/covariance_computation.R first."))
}
load(file = data_filename)

# Index to start with -h as the first lag
plot_idx <- c(2, seq_len(N))
plot_mult <- c(-1, rep(1, N))

data <- NULL
for (model.id in 1:length(M)) {
  if (is.null(M[[model.id]])) next
  data <- rbind(
    data,
    tibble(
      ID = as.character(model.id),
      Type = model_defn[model.id, "Type"],
      nu_s = with(model_defn[model.id, ], alpha_e + alpha_s * (alpha_t - 0.5) - DIM / 2),
      nu_t = min(
        model_defn[model.id, "alpha_t"] - 0.5,
        nu_s / model_defn[model.id, "alpha_s"]
      ),
      beta_s = 1 - model_defn[model.id, "alpha_e"] / (nu_s + DIM / 2),
      x = rep(hs[plot_idx] * plot_mult, times = N + 1),
      y = rep(ht[plot_idx] * plot_mult, each = N + 1),
      z = as.vector(M[[model.id]]$C[plot_idx, plot_idx])
    )
  )
}
# Compute relative temporal decay and relative spatial decay
data <- data %>%
  group_by(ID, x) %>%
  mutate(z_temporal_decay = z / z[y == 0]) %>%
  ungroup() %>%
  group_by(ID, y) %>%
  mutate(z_spatial_decay = z / z[x == 0]) %>%
  ungroup()

# All marginals
pl1 <- ggplot(data %>% filter(y == 0, x >= 0)) +
  geom_line(aes(x = x, y = z, color = Type, linetype = Type)) +
  xlab(expression("" ~ h[s])) +
  ylab("Covariance") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(color = "Model", linetype = "Model")
pl2 <- ggplot(data %>% filter(x == 0, y >= 0)) +
  geom_line(aes(x = y, y = z, color = Type, linetype = Type)) +
  xlab(expression("" ~ h[t])) +
  ylab("Covariance") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(color = "Model", linetype = "Model")
pl_1 <- (pl1 / pl2) +
  plot_layout(guides = "collect") &
  labs(color = "Model") &
  theme_minimal()

# Spacetime contours
getPalette <- colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))
n_levels <- 10
my_labels <- as.character(seq_len(n_levels - 1) / n_levels)
my_labels[seq(1, n_levels - 1, 2)] <- ""
pl_2 <- ggplot(data) +
  geom_contour_fill(aes(x = x, y = y, z = z, fill = after_stat(level)),
    breaks = (0:n_levels) / n_levels,
  ) +
  geom_contour(aes(x = x, y = y, z = z),
    color = "black",
    breaks = (0:n_levels) / n_levels
  ) +
  geom_text_contour(aes(x = x, y = y, z = z),
    breaks = (2:8) / 10,
    stroke = 0.2,
    label.placer = label_placer_flattest(ref_angle = -15)
  ) +
  scale_fill_distiller(
    super = ScaleDiscretised, type = "seq", palette = "Blues", direction = +1,
    labels = my_labels,
    guide = guide_colourbar(ticks = FALSE)
  ) +
  geom_contour(aes(x = x, y = y, z = z_temporal_decay),
    color = "grey",
    breaks = c(0.2, 0.4, 0.6, 0.8)
  ) +
  geom_contour(aes(x = x, y = y, z = z_spatial_decay),
    color = "grey",
    breaks = c(0.2, 0.4, 0.6, 0.8)
  ) +
  facet_wrap(vars(Type), nrow = 2) + # , labeller=label_parsed) +
  theme_minimal() +
  xlab(expression("" ~ h[s])) +
  ylab(expression("" ~ h[t]))

plotfile_all_margs <- file.path(figures.dir, paste0("covar-all-marginals", ifelse(DIM == 2, "", "-1"), ".pdf"))
plotfile_bivar <- file.path(figures.dir, paste0("covar-bivar", ifelse(DIM == 2, "", "-1"), ".pdf"))
if (!overwrite_plots &&
  file.exists(plotfile_all_margs)) {
  message("Displaying all marginal covariances")
  print(pl_1)
} else {
  message(paste0("Saving plot with all marginals to ", plotfile_all_margs))
  dir.create(figures.dir, showWarnings = FALSE, recursive = TRUE)
  pdf(plotfile_all_margs,
    pointsize = 20,
    width = 7,
    height = 4
  )
  print(pl_1)
  dev.off()
}
if (!overwrite_plots &&
  file.exists(plotfile_bivar)) {
  message("Displaying bivariate covariances")
  print(pl_2)
} else {
  message(paste0("Saving plot with bivariate covariances to ", plotfile_bivar))
  pdf(plotfile_bivar,
    pointsize = 20,
    width = 7,
    height = 4
  )
  print(pl_2)
  dev.off()
}
