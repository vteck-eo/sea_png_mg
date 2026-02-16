############################################################
# 🌏 Monte-Carlo Nearest-Neighbor Distances (2024)

############################################################

library(dplyr)
library(readr)
library(ggplot2)
library(sf)
library(cowplot)
library(jsonlite)

# ==========================================================
# 📂 Input / Output
# ==========================================================
in_file <- "username"
out_dir <- "username/"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

year <- 2024

# ==========================================================
# 🔧 Utility functions
# ==========================================================
deg2rad <- function(d) d * pi / 180
R_earth <- 6371  # km

extract_lonlat <- function(geo_json) {
  if (is.na(geo_json) || geo_json == "") return(c(NA, NA))
  g <- tryCatch(fromJSON(geo_json), error = function(e) NULL)
  if (is.null(g) || is.null(g$coordinates)) return(c(NA, NA))
  c(lon_x = g$coordinates[1], lat_y = g$coordinates[2])
}

# ==========================================================
# 🚀 MONTE-CARLO WORKFLOW (2024)
# ==========================================================
message("========== Processing year 2024 ==========")

# 1️⃣ Load data
df <- read_csv(in_file, show_col_types = FALSE)

# Extract lon/lat if needed
if (!("lon_x" %in% names(df)) || !("lat_y" %in% names(df))) {
  if (!(".geo" %in% names(df))) stop("No lon/lat or .geo column found.")
  coords <- t(sapply(df$.geo, extract_lonlat))
  coords <- as.data.frame(coords)
  df$lon_x <- coords$lon_x
  df$lat_y <- coords$lat_y
}

class_col <- if ("land_class" %in% names(df)) "land_class" else "class"

# Mangroves only
df_mg <- df %>%
  filter(!is.na(lon_x), !is.na(lat_y), !!sym(class_col) == 1) %>%
  select(lon_x, lat_y)

message("Mangrove points:", nrow(df_mg))
if (nrow(df_mg) < 10) stop("Too few mangrove samples")

# ==========================================================
# 2️⃣ 25 km grid in Web Mercator
# ==========================================================
sf_pts <- st_as_sf(df_mg, coords = c("lon_x", "lat_y"), crs = 4326)
sf_3857 <- st_transform(sf_pts, 3857)

grid <- st_make_grid(sf_3857, cellsize = 25000, what = "polygons")
grid_sf <- st_sf(id = seq_along(grid), geometry = grid)

joined <- st_join(sf_3857, grid_sf, left = FALSE)

agg <- joined %>%
  st_drop_geometry() %>%
  group_by(id) %>%
  summarise(n = n(), .groups = "drop")

# ==========================================================
# 3️⃣ Grid centroids
# ==========================================================
centroids <- st_centroid(grid_sf) %>% st_transform(4326)
agg_sf <- left_join(agg, centroids, by = "id") %>% st_as_sf()

coords <- st_coordinates(agg_sf)
agg_sf$lon <- coords[, 1]
agg_sf$lat <- coords[, 2]

N <- nrow(agg_sf)
message("Occupied grid cells:", N)
if (N < 5) stop("Too few grid cells")

# ==========================================================
# 4️⃣ Monte-Carlo settings
# ==========================================================
draws <- 1000
max_k <- min(1000, N - 1)

lat_r <- deg2rad(agg_sf$lat)
lon_r <- deg2rad(agg_sf$lon)

tot_mat <- matrix(NA, draws, max_k)
ns_mat  <- matrix(NA, draws, max_k)
ew_mat  <- matrix(NA, draws, max_k)

# ==========================================================
# 5️⃣ Monte-Carlo simulation
# ==========================================================
set.seed(42)

for (i in 1:draws) {

  idx <- sample.int(N, 1)

  dlat <- lat_r - lat_r[idx]
  dlon <- lon_r - lon_r[idx]

  a <- sin(dlat/2)^2 +
       cos(lat_r[idx]) * cos(lat_r) * sin(dlon/2)^2
  gc <- 2 * R_earth * atan2(sqrt(a), sqrt(1 - a))

  ns <- abs(dlat) * R_earth
  ew <- abs(dlon) * R_earth * cos(lat_r[idx])

  gc[idx] <- ns[idx] <- ew[idx] <- NA

  ord <- order(gc, na.last = NA)[1:max_k]

  tot_mat[i, ] <- gc[ord]
  ns_mat[i, ]  <- ns[ord]
  ew_mat[i, ]  <- ew[ord]
}

# ==========================================================
# 6️⃣ Summary statistics
# ==========================================================
df_nn <- data.frame(
  k = 1:max_k,

  tot_mean = colMeans(tot_mat, na.rm = TRUE),
  tot_low  = apply(tot_mat, 2, quantile, 0.10, na.rm = TRUE),
  tot_high = apply(tot_mat, 2, quantile, 0.90, na.rm = TRUE),

  ns_mean  = colMeans(ns_mat, na.rm = TRUE),
  ns_low   = apply(ns_mat, 2, quantile, 0.10, na.rm = TRUE),
  ns_high  = apply(ns_mat, 2, quantile, 0.90, na.rm = TRUE),

  ew_mean  = colMeans(ew_mat, na.rm = TRUE),
  ew_low   = apply(ew_mat, 2, quantile, 0.10, na.rm = TRUE),
  ew_high  = apply(ew_mat, 2, quantile, 0.90, na.rm = TRUE)
)

write_csv(
  df_nn,
  file.path(out_dir, "MC_NN_Distances_2024.csv")
)

# ==========================================================
# 7️⃣ Plots (3-panel)
# ==========================================================
p_tot <- ggplot(df_nn, aes(k)) +
  geom_ribbon(aes(ymin = tot_low, ymax = tot_high), fill = "grey80") +
  geom_line(aes(y = tot_mean), linewidth = 0.9) +
  labs(title = "Total distance", x = "Neighbor rank", y = "Distance (km)") +
  theme_bw(base_size = 11)

p_ns <- ggplot(df_nn, aes(k)) +
  geom_ribbon(aes(ymin = ns_low, ymax = ns_high), fill = "#8DA0CB", alpha = 0.4) +
  geom_line(aes(y = ns_mean), color = "#1F3A93", linewidth = 0.9) +
  labs(title = "North–South", x = "Neighbor rank", y = "Distance (km)") +
  theme_bw(base_size = 11)

p_ew <- ggplot(df_nn, aes(k)) +
  geom_ribbon(aes(ymin = ew_low, ymax = ew_high), fill = "#FC8D62", alpha = 0.4) +
  geom_line(aes(y = ew_mean), color = "#8B0000", linewidth = 0.9) +
  labs(title = "East–West", x = "Neighbor rank", y = "Distance (km)") +
  theme_bw(base_size = 11)

final_plot <- plot_grid(p_tot, p_ns, p_ew, ncol = 3)

# ==========================================================
# 💾 Save
# ==========================================================
ggsave(
  file.path(out_dir, "Monte Carlos_2024.png"),
  final_plot,
  width = 18,
  height = 6,
  dpi = 300
)

cat("\n✅ Done :\n", out_dir, "\n")
