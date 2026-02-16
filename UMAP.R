
############################################################
# 🌏 UMAP – Mangrove Embedding Structure (2017–2024)
############################################################

# =========================================================
# 📦 Libraries
# =========================================================
library(dplyr)
library(readr)
library(ggplot2)
library(viridis)
library(jsonlite)
library(umap)

# =========================================================
# 📂 Paths & Years
# =========================================================
in_dir  <- "username/"
out_dir <- "username/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

years <- 2017:2024

# =========================================================
# 🔧 Function: Extract lon/lat from .geo JSON
# =========================================================
extract_lonlat <- function(geo_json) {
  if (is.na(geo_json) || geo_json == "") return(c(NA, NA))
  g <- tryCatch(fromJSON(geo_json), error = function(e) NULL)
  if (is.null(g) || is.null(g$coordinates)) return(c(NA, NA))
  c(lon_x = g$coordinates[1], lat_y = g$coordinates[2])
}

# =========================================================
# 1️⃣ Load Data
# =========================================================
all_data <- lapply(years, function(yr) {

  f <- file.path(in_dir, sprintf("sea_emb_samples_%d_all.csv", yr))
  if (!file.exists(f)) return(NULL)

  d <- read_csv(f, show_col_types = FALSE)

  coords <- t(sapply(d$.geo, extract_lonlat))
  coords <- as.data.frame(coords)

  d$lon_x <- coords$lon_x
  d$lat_y <- coords$lat_y
  d$year  <- yr
  d

}) %>% bind_rows()

# =========================================================
# 2️⃣ Detect Embedding Columns
# =========================================================
numeric_cols <- names(all_data)[sapply(all_data, is.numeric)]
exclude <- c("lon_x","lat_y","year","class","land_class","class_num")
emb_cols <- setdiff(numeric_cols, exclude)

# =========================================================
# 3️⃣ Clean Dataset
# =========================================================
class_col <- if ("land_class" %in% names(all_data)) "land_class" else "class"

df <- all_data %>%
  filter(!is.na(lon_x),
         !is.na(lat_y),
         !!sym(class_col) %in% c(0,1)) %>%
  mutate(
    class_num = !!sym(class_col),
    year      = factor(year)
  ) %>%
  select(year, class_num, all_of(emb_cols))

# Optional subsampling for speed
set.seed(42)
n_max <- 60000
if (nrow(df) > n_max) df <- df %>% sample_n(n_max)

emb_mat <- as.matrix(df %>% select(all_of(emb_cols)))

# =========================================================
# 4️⃣ Run UMAP
# =========================================================
umap_config <- umap::umap.defaults
umap_config$n_neighbors <- 30
umap_config$min_dist    <- 0.1
umap_config$metric      <- "euclidean"

umap_res <- umap::umap(emb_mat, config = umap_config)

df$UMAP1 <- umap_res$layout[,1]
df$UMAP2 <- umap_res$layout[,2]

# =========================================================
# 5️⃣ Plot
# =========================================================
pal <- viridis(length(levels(df$year)), option = "D")

p_umap <- ggplot() +
  geom_point(
    data = df %>% filter(class_num == 0),
    aes(UMAP1, UMAP2),
    color = "grey90",
    alpha = 0.3,
    size = 0.5
  ) +
  geom_point(
    data = df %>% filter(class_num == 1),
    aes(UMAP1, UMAP2, color = year),
    size = 1.1,
    alpha = 0.9
  ) +
  scale_color_manual(values = pal, name = "Year") +
  theme_minimal(base_size = 14) +
  labs(x = "UMAP-1", y = "UMAP-2") +
  theme(
    axis.line = element_line(color="black", linewidth=1),
    axis.ticks = element_line(color="black"),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(face="bold")
  )

# =========================================================
# 💾 Save
# =========================================================
ggsave(
  file.path(out_dir, "UMAP_mg.png"),
  p_umap,
  width = 10,
  height = 7,
  dpi = 300,
  bg = "white"
)

cat("✅ UMAP export.\n")

