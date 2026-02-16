############################################################ 
# 🌏 UMAP – Mangrove Class Separation Over Years
############################################################

# =========================================================
# 📦 Libraries
# =========================================================
library(dplyr)
library(readr)
library(ggplot2)
library(viridis)
library(umap)
library(jsonlite)

# =========================================================
# 📂 Paths & years
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
# 1️⃣ Load all years
# =========================================================
all_data <- lapply(years, function(yr) {

  f <- file.path(in_dir, sprintf("sea_emb_samples_%d_all.csv", yr))
  if (!file.exists(f)) {
    warning("⚠ Missing file: ", f)
    return(NULL)
  }

  d <- read_csv(f, show_col_types = FALSE)

  if (!(".geo" %in% names(d)))
    stop("❌ '.geo' column not found in: ", f)

  coords <- t(sapply(d$.geo, extract_lonlat))
  coords <- as.data.frame(coords)

  d$lon_x <- coords$lon_x
  d$lat_y <- coords$lat_y
  d$year  <- yr

  d
}) %>% bind_rows()

cat("Total loaded:", nrow(all_data), "\n")

# =========================================================
# 2️⃣ Detect embedding columns
# =========================================================
numeric_cols <- names(all_data)[sapply(all_data, is.numeric)]

exclude <- c("lon_x","lat_y","year",
             "class","land_class","class_num")

emb_candidates <- setdiff(numeric_cols, exclude)

emb_cols <- emb_candidates[
  sapply(all_data[emb_candidates], function(x)
    length(unique(x)) > 50)
]

cat("Embedding dimensions detected:", length(emb_cols), "\n")

if (length(emb_cols) == 0)
  stop("❌ No embedding columns detected!")

# =========================================================
# 3️⃣ Clean dataset
# =========================================================
class_col <- if ("land_class" %in% names(all_data))
  "land_class" else "class"

df <- all_data %>%
  filter(
    !is.na(lon_x),
    !is.na(lat_y),
    !!sym(class_col) %in% c(0,1)
  ) %>%
  mutate(
    class_num   = !!sym(class_col),
    class_label = ifelse(class_num == 1,
                         "Mangrove",
                         "Non-mangrove"),
    year        = factor(year, levels = years)
  ) %>%
  select(year, class_num, class_label,
         lon_x, lat_y, all_of(emb_cols))

cat("Points for analysis:", nrow(df), "\n")

# =========================================================
# 4️⃣ Subsample (max 60k)
# =========================================================
set.seed(42)
n_max <- 60000

df_sub <- if (nrow(df) > n_max)
  df %>% sample_n(n_max)
else df

cat("Subsampled:", nrow(df_sub), "\n")

# =========================================================
# 5️⃣ Run UMAP
# =========================================================
cat("Running UMAP...\n")

emb_mat <- as.matrix(df_sub %>% select(all_of(emb_cols)))

# IMPORTANT: scale embeddings
emb_mat_scaled <- scale(emb_mat)

umap_cfg <- umap::umap.defaults
umap_cfg$n_neighbors <- 30
umap_cfg$min_dist    <- 0.1
umap_cfg$metric      <- "euclidean"

umap_res <- umap::umap(emb_mat_scaled,
                       config = umap_cfg)

df_sub$UMAP1 <- umap_res$layout[,1]
df_sub$UMAP2 <- umap_res$layout[,2]

cat("UMAP completed.\n")

# =========================================================
# 🎨 Colors
# =========================================================
pal <- viridis::viridis(length(years),
                        option = "D")

# =========================================================
# 📊 Plot UMAP
# =========================================================
p_umap <- ggplot() +

  # Non-mangrove background
  geom_point(
    data = df_sub %>% filter(class_num == 0),
    aes(UMAP1, UMAP2),
    color = "grey90",
    alpha = 0.3,
    size = 0.5
  ) +

  # Mangrove colored by year
  geom_point(
    data = df_sub %>% filter(class_num == 1),
    aes(UMAP1, UMAP2, color = year),
    size = 1.1,
    alpha = 0.9
  ) +

  scale_color_manual(
    name = "Year",
    values = pal,
    guide = guide_legend(
      override.aes = list(size = 4)
    )
  ) +

  theme_minimal(base_size = 14) +

  labs(
    x = "UMAP-1",
    y = "UMAP-2"
  ) +

  theme(
    axis.line = element_line(color="black",
                             size=0.5),
    axis.ticks = element_line(color="black",
                              size=0.4),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color="grey92",
                                    size=0.3),
    legend.position = "right"
  )

# =========================================================
# 💾 Save
# =========================================================
ggsave(
  file.path(out_dir,
            "UMAP_Mangrove_2017_2024.png"),
  p_umap,
  width = 10,
  height = 7,
  dpi = 300,
  bg = "white"
)

cat("✅ DONE.\n")
