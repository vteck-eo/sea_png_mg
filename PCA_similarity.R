############################################################
# 🌏 PCA RGB – Mangroves Only, 2017–2024
############################################################

# =========================================================
# 📦 Libraries
# =========================================================
library(dplyr)
library(readr)
library(ggplot2)
library(sf)
library(proxy)
library(irlba)
library(cowplot)

# =========================================================
# 📂 Paths & years
# =========================================================
in_dir  <- "username/"
out_dir <- "username/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

years <- 2017:2024

############################################################
# 🔧 Helper: compute PCA RGB map for ONE YEAR
############################################################
compute_pca_rgb_map <- function(year) {

  message("Processing year ", year)

  f <- file.path(in_dir, sprintf("sea_emb_samples_%d_all.csv", year))
  if (!file.exists(f)) return(NULL)

  # ---------- 1) Load & filter mangroves ----------
  d <- read_csv(f, show_col_types = FALSE)

  if (!all(c("lon_x","lat_y") %in% names(d))) return(NULL)

  class_col <- if ("land_class" %in% names(d)) "land_class" else "class"

  emb_cols <- grep("^a\\d{1,3}$", names(d), value = TRUE, ignore.case = TRUE)
  if (length(emb_cols) < 3) return(NULL)

  d[emb_cols] <- lapply(d[emb_cols], as.numeric)

  df <- d %>%
    filter(
      !is.na(lon_x),
      !is.na(lat_y),
      !!sym(class_col) == 1
    ) %>%
    select(lon_x, lat_y, all_of(emb_cols))

  if (nrow(df) < 10) return(NULL)

  # ---------- 2) 25 km grid aggregation ----------
  sf_pts  <- st_as_sf(df, coords = c("lon_x","lat_y"), crs = 4326)
  sf_3857 <- st_transform(sf_pts, 3857)

  grid <- st_make_grid(sf_3857, cellsize = 25000, what = "polygons")
  grid_sf <- st_sf(id = seq_along(grid), geometry = grid)

  joined <- st_join(sf_3857, grid_sf, left = FALSE)
  jdf <- joined %>% st_drop_geometry()

  agg <- jdf %>%
    group_by(id) %>%
    summarise(
      across(all_of(emb_cols), ~ median(.x, na.rm = TRUE)),
      .groups = "drop"
    )

  if (nrow(agg) < 2) return(NULL)

  centroids <- st_centroid(grid_sf) %>% st_transform(4326)
  agg_sf <- left_join(agg, centroids, by = "id") %>% st_as_sf()

  coords <- st_coordinates(agg_sf)
  agg_sf$lon <- coords[,1]
  agg_sf$lat <- coords[,2]

  # ---------- 3) Cosine similarity ----------
  emb_mat <- as.matrix(agg_sf %>% st_drop_geometry() %>% select(all_of(emb_cols)))
  sim_mat <- 1 - as.matrix(proxy::dist(emb_mat, method = "cosine"))
  sim_mat[!is.finite(sim_mat)] <- 0

  # ---------- 4) PCA ----------
  pca <- tryCatch(
    irlba::prcomp_irlba(sim_mat, n = 3, scale. = TRUE),
    error = function(e) NULL
  )
  if (is.null(pca)) return(NULL)

  # ---------- 5) Normalize to RGB ----------
  rgb_norm <- apply(pca$x[,1:3], 2, function(x) {
    r <- range(x, na.rm = TRUE)
    if (diff(r) == 0) rep(0.5, length(x)) else (x - r[1]) / diff(r)
  })

  agg_sf$rgb <- rgb(
    pmin(pmax(rgb_norm[,1],0),1),
    pmin(pmax(rgb_norm[,2],0),1),
    pmin(pmax(rgb_norm[,3],0),1)
  )

  # ---------- 6) Plot (NO TITLE) ----------
  p <- ggplot(agg_sf) +
    geom_point(
      aes(lon, lat, color = rgb),
      size = 1.1,
      alpha = 0.9,
      show.legend = FALSE
    ) +
    scale_color_identity() +
    coord_sf(
      xlim = c(88, 160),
      ylim = c(-14, 25),
      expand = FALSE
    ) +
    labs(
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank()
    )

  return(p)
}

############################################################
# 🌟 Generate plots for all years
############################################################
plots_list <- lapply(years, compute_pca_rgb_map)
plots_list <- plots_list[!sapply(plots_list, is.null)]

if (length(plots_list) == 0) stop("No plots created.")

############################################################
# ⭐ Combine into 2 × 4 panel with A–H labels
############################################################
final_plot <- plot_grid(
  plotlist = plots_list,
  ncol = 2,
  labels = LETTERS[1:length(plots_list)],
  label_size = 18,
  label_fontface = "bold",
  label_x = 0.97,
  label_y = 0.97,
  hjust = 1,
  vjust = 1
)

############################################################
# 💾 Save output
############################################################
ggsave(
  file.path(out_dir, "PCA_RGB_Mangroves_2017_2024.jpg"),
  final_plot,
  width = 20,
  height = 28,
  dpi = 300
)

cat("\n✅ DONE :\n", out_dir, "\n")
