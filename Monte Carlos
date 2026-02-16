############################################################
# 🌏 Monte Carlo NN distances for Mangroves, 2017–2024
############################################################

# =========================================================
# 📦 Libraries
# =========================================================
library(dplyr)
library(readr)
library(ggplot2)
library(sf)
library(cowplot)

# =========================================================
# 📂 Paths & years
# =========================================================
in_dir   <- "username/"
out_root <- "username/"
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

years <- 2017:2024

# =========================================================
# 🔧 Helper functions
# =========================================================
deg2rad <- function(d) d * pi / 180
R_earth <- 6371  # km

# ---------------------------------------------------------
# Function: compute Monte-Carlo NN panels for one year
# ---------------------------------------------------------
run_mc_for_year <- function(yr) {
  message("========== Year ", yr, " ==========")

  # ---------- 1) Load data ----------
  f <- file.path(in_dir, sprintf("sea_emb_samples_%d_all.csv", yr))
  if (!file.exists(f)) {
    stop("File not found for year ", yr, ": ", f)
  }

  d <- read_csv(f, show_col_types = FALSE)
  message("Loaded ", nrow(d), " points for ", yr)

  # sanity checks
  stopifnot("lon_x" %in% names(d))
  stopifnot("lat_y" %in% names(d))

  # identify class column
  class_col <- if ("land_class" %in% names(d)) "land_class" else "class"

  # keep only mangroves & valid coords
  df <- d %>%
    filter(
      !is.na(lon_x),
      !is.na(lat_y),
      !!sym(class_col) == 1
    ) %>%
    select(lon_x, lat_y)

  message("Mangrove-only points: ", nrow(df))

  if (nrow(df) < 10) {
    stop("Too few mangrove points for year ", yr, " to run Monte Carlo.")
  }

  # ---------- 2) 25 km grid and centroids ----------
  sf_pts  <- st_as_sf(df, coords = c("lon_x", "lat_y"), crs = 4326)
  sf_3857 <- st_transform(sf_pts, 3857)

  grid <- st_make_grid(sf_3857, cellsize = 25000, what = "polygons")
  grid_sf <- st_sf(id = seq_along(grid), geometry = grid)

  joined <- st_join(sf_3857, grid_sf, join = st_intersects, left = FALSE)
  agg <- joined %>%
    st_drop_geometry() %>%
    as.data.frame() %>%
    group_by(id) %>%
    summarise(n = n(), .groups = "drop")

  # keep centroids only for occupied cells
  centroids <- st_centroid(grid_sf) %>% st_transform(4326)
  agg_sf <- left_join(agg, centroids, by = "id") %>% st_as_sf()

  coords <- st_coordinates(agg_sf)
  agg_sf$lon_grid <- coords[, 1]
  agg_sf$lat_grid <- coords[, 2]

  N <- nrow(agg_sf)
  message("25 km grid cells with mangroves: ", N)

  if (N < 5) {
    stop("Not enough occupied grid cells for year ", yr, " to run Monte Carlo.")
  }

  # ---------- 3) Monte-Carlo NN distances ----------
  max_k <- min(1000, N - 1)   # up to 800 neighbors
  draws <- 1000

  lat_r <- deg2rad(agg_sf$lat_grid)
  lon_r <- deg2rad(agg_sf$lon_grid)

  tot_mat <- matrix(NA, draws, max_k)
  ns_mat  <- matrix(NA, draws, max_k)
  ew_mat  <- matrix(NA, draws, max_k)

  message("Running Monte Carlo (", draws,
          " draws, up to k = ", max_k, ") for ", yr, "...")

  for (i in 1:draws) {

    idx <- sample.int(N, 1)

    dlat <- lat_r - lat_r[idx]
    dlon <- lon_r - lon_r[idx]

    a <- sin(dlat / 2)^2 +
      cos(lat_r[idx]) * cos(lat_r) * sin(dlon / 2)^2
    gc <- R_earth * 2 * atan2(sqrt(a), sqrt(1 - a))

    ns <- abs(dlat) * R_earth
    ew <- abs(dlon) * R_earth * cos(lat_r[idx])

    # remove self
    gc[idx] <- NA
    ns[idx] <- NA
    ew[idx] <- NA

    ord <- order(gc, na.last = NA)[1:max_k]

    tot_mat[i, ] <- gc[ord]
    ns_mat[i, ]  <- ns[ord]
    ew_mat[i, ]  <- ew[ord]
  }

  df_nn <- data.frame(
    k        = 1:max_k,
    tot_mean = apply(tot_mat, 2, mean,     na.rm = TRUE),
    tot_low  = apply(tot_mat, 2, quantile, 0.10, na.rm = TRUE),
    tot_high = apply(tot_mat, 2, quantile, 0.90, na.rm = TRUE),
    ns_mean  = apply(ns_mat,  2, mean,     na.rm = TRUE),
    ns_low   = apply(ns_mat,  2, quantile, 0.10, na.rm = TRUE),
    ns_high  = apply(ns_mat,  2, quantile, 0.90, na.rm = TRUE),
    ew_mean  = apply(ew_mat,  2, mean,     na.rm = TRUE),
    ew_low   = apply(ew_mat,  2, quantile, 0.10, na.rm = TRUE),
    ew_high  = apply(ew_mat,  2, quantile, 0.90, na.rm = TRUE)
  )

  # ---------- 4) Plots for this year ----------
  p_tot_nn <- ggplot(df_nn, aes(k)) +
    geom_ribbon(aes(ymin = tot_low, ymax = tot_high), fill = "grey80") +
    geom_line(aes(y = tot_mean), size = 0.8) +
    labs(title = "Total distance", x = "Neighbor rank", y = "Distance (km)") +
    theme_bw(base_size = 10)

  p_ns_nn <- ggplot(df_nn, aes(k)) +
    geom_ribbon(aes(ymin = ns_low, ymax = ns_high),
                fill = "skyblue", alpha = .4) +
    geom_line(aes(y = ns_mean), color = "navy", size = 0.8) +
    labs(title = "North–South", x = "Neighbor rank", y = "Distance (km)") +
    theme_bw(base_size = 10)

  p_ew_nn <- ggplot(df_nn, aes(k)) +
    geom_ribbon(aes(ymin = ew_low, ymax = ew_high),
                fill = "salmon", alpha = .4) +
    geom_line(aes(y = ew_mean), color = "darkred", size = 0.8) +
    labs(title = "East–West", x = "Neighbor rank", y = "Distance (km)") +
    theme_bw(base_size = 10)

  # row of 3 panels for this year
  row_panels <- plot_grid(p_tot_nn, p_ns_nn, p_ew_nn, ncol = 3)

  # add a small year title above this row
  year_panel <- ggdraw() +
    draw_label(
      paste0("Year ", yr),
      fontface = "bold",
      size     = 11,
      x        = 0.02, y = 0.98,
      hjust    = 0
    ) +
    draw_plot(row_panels, x = 0, y = 0, width = 1, height = 0.94)

  # (optional) save per-year CSV of distances
  write.csv(
    df_nn,
    file.path(out_root, sprintf("Mangrove_NN_Distances_GridOnly_%d.csv", yr)),
    row.names = FALSE
  )

  return(year_panel)
}

# =========================================================
# 🚀 Run for all years & combine into 2×4 grid
# =========================================================
year_panels <- lapply(years, run_mc_for_year)

big_panel <- plot_grid(
  plotlist   = year_panels,
  ncol       = 2,      # 2 columns → 4 rows for 8 years
  rel_heights = rep(1, length(years))
)

# Add global title
final_plot <- ggdraw() +
  draw_label(
    "Monte Carlo Nearest-Neighbor Distances (Mangroves Only, 25 km Grid, 2017–2024)",
    fontface = "bold",
    size     = 16,
    x        = 0.5, y = 0.99, hjust = 0.5
  ) +
  draw_plot(big_panel, x = 0, y = 0, width = 1, height = 0.95)

ggsave(
  file.path(out_root, "MC_NN_Mangroves_2017_2024_2x4years_3panelsEach.png"),
  final_plot,
  width  = 18,
  height = 20,
  dpi    = 300
)

ggsave(
  file.path(out_root, "Monte Carlos Distance.png"),
  final_plot,
  width  = 24,   # wider
  height = 32,   # taller
  dpi    = 300
)

cat("🎉 Done:\n", out_root, "\n")
