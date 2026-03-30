# Paleobiology Database Cluster Analysis
# Felicia Selbst
#
# This project merges fossil occurrence datasets from multiple regions,
# engineers categorical geological features, and applies PAM clustering
# with Gower distance to identify patterns across fossil localities.

library(cluster)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tidyverse)

# -------------------- Read in Regional Datasets -------------------- #

xin      <- read.csv("data/dinoXinjiang.csv")
west     <- read.csv("data/dinoWesternCape.csv")
vic      <- read.csv("data/dinoVictoria.csv")
utah     <- read.csv("data/dinoUtah.csv")
ulaa     <- read.csv("data/dinoUlaanbaatar.csv")
tran     <- read.csv("data/dinoTransylvania.csv")
thai     <- read.csv("data/dinoThailand.csv")
tex      <- read.csv("data/dinoTexas.csv")
south    <- read.csv("data/dinoSouthIsland.csv")
somer    <- read.csv("data/dinoSomerset.csv")
sici     <- read.csv("data/dinoSicily.csv")
santa    <- read.csv("data/dinoSantaCruz.csv")
rio      <- read.csv("data/dinoRioGrande.csv")
queens   <- read.csv("data/dinoQueensland.csv")
norma    <- read.csv("data/dinoNormandy.csv")
moro     <- read.csv("data/dinoMorocco.csv")
montana  <- read.csv("data/dinoMontana.csv")
meta     <- read.csv("data/dinoMeta.csv")
liao     <- read.csv("data/dinoLiaoning.csv")
karo     <- read.csv("data/dinoKaronga.csv")
hokk     <- read.csv("data/dinoHokkaido.csv")
guj      <- read.csv("data/dinoGujarat.csv")
coa      <- read.csv("data/dinoCoahuila.csv")
chu      <- read.csv("data/dinoChubut.csv")
cata     <- read.csv("data/dinoCatalonia.csv")
briCol   <- read.csv("data/dinoBritColumbia.csv")
bavar    <- read.csv("data/dinoBavaria.csv")
aya      <- read.csv("data/dinoAyacucho.csv")
alb      <- read.csv("data/dinoAlberta.csv")
afar     <- read.csv("data/dinoAfar.csv")

locations <- list(
  xin, west, vic, utah, ulaa, tran, thai, tex, south, somer,
  sici, santa, rio, queens, norma, moro, montana, meta, liao, karo,
  hokk, guj, coa, chu, cata, briCol, bavar, aya, alb, afar
)

# -------------------- Merge and Reduce Data -------------------- #

merged_data <- do.call(rbind, locations)
merged_data <- unique(merged_data)

cols_reduced <- c(
  "occurrence_no", "accepted_name", "accepted_rank", "early_interval",
  "max_ma", "min_ma", "composition", "lng", "lat", "formation",
  "lithology1", "lithadj1", "lithology2", "environment"
)

reduced_data <- merged_data[, names(merged_data) %in% cols_reduced]

# -------------------- Clean Composition -------------------- #

reduced_data <- reduced_data %>%
  mutate(
    composition_clean = str_to_lower(str_trim(str_remove_all(composition, '"'))),
    composition_group = case_when(
      str_detect(composition_clean, "no hard parts") ~ "none",
      str_detect(composition_clean, "calcite|aragonite") ~ "carbonate",
      str_detect(composition_clean, "apatite|phosphatic") ~ "phosphate",
      str_detect(composition_clean, "silica|siliceous|opal") ~ "silica",
      str_detect(composition_clean, "chitin|lignin|sclero-protein|agglutinated") ~ "organic",
      str_detect(composition_clean, ",") ~ "mixed",
      TRUE ~ "unknown"
    )
  )

# -------------------- Clean Lithology -------------------- #

carbonate_terms <- c(
  "limestone", "wackestone", "packstone", "grainstone", "rudstone",
  "framestone", "floatstone", "bafflestone", "bindstone",
  "reef rocks", "marl", "dolomite"
)

clastic_terms <- c(
  "shale", "mudstone", "siltstone", "sandstone", "claystone",
  "conglomerate", "gravel"
)

silica_terms <- c(
  "chert", "diatomite", "phosphorite", "evaporite", "travertine"
)

volcanic_terms <- c(
  "tuff", "volcaniclastic", "ash"
)

organic_terms <- c(
  "coal", "lignite", "peat", "subbituminous coal"
)

categorize_lith <- function(lith) {
  case_when(
    lith %in% carbonate_terms ~ "carbonate",
    lith %in% clastic_terms   ~ "clastic",
    lith %in% silica_terms    ~ "silica_chemical",
    lith %in% volcanic_terms  ~ "volcanic",
    lith %in% organic_terms   ~ "organic",
    is.na(lith) | lith == ""  ~ "unknown",
    TRUE                      ~ "other"
  )
}

reduced_data <- reduced_data %>%
  mutate(
    lith1_clean = str_to_lower(str_trim(str_remove_all(lithology1, '"'))),
    lith2_clean = str_to_lower(str_trim(str_remove_all(lithology2, '"'))),
    lith1_group = categorize_lith(lith1_clean),
    lith2_group = categorize_lith(lith2_clean),
    lith_combined = case_when(
      lith1_group == lith2_group ~ lith1_group,
      lith2_group == "unknown" ~ lith1_group,
      lith1_group == "unknown" ~ lith2_group,
      TRUE ~ paste(lith1_group, lith2_group, sep = "_")
    )
  )

# -------------------- Clean Environment -------------------- #

reduced_data <- reduced_data %>%
  mutate(
    env_clean = str_to_lower(str_trim(environment)),
    environment_group = case_when(
      str_detect(env_clean, "marine|subtidal|shelf|offshore|basinal|reef|perireef|platform|slope|ramp") ~ "marine",
      str_detect(env_clean, "delta|shoreface|foreshore|estuary|bay|lagoon|paralic|sand shoal|prodelta|beach") ~ "marginal_marine",
      str_detect(env_clean, "fluvial|channel|crevasse|floodplain|alluvial|levee|braid") ~ "fluvial",
      str_detect(env_clean, "lacustrine|lake|pond") ~ "lacustrine",
      str_detect(env_clean, "terrestrial|eolian|dune|mire|swamp|soil|glacial|spring|sinkhole") ~ "terrestrial",
      TRUE ~ "unknown"
    )
  )

# -------------------- Rename and Select Analysis Columns -------------------- #

cols_mutated_reduced <- c(
  "accepted_name", "accepted_rank", "early_interval", "max_ma", "min_ma",
  "composition_clean", "lng", "lat", "composition_group",
  "lith1_clean", "lith2_clean", "lith1_group", "lith2_group",
  "lith_combined", "env_clean", "environment_group"
)

mutated_reduced_data <- reduced_data[, names(reduced_data) %in% cols_mutated_reduced]

mutated_reduced_data <- mutated_reduced_data %>%
  rename(
    taxon        = accepted_name,
    rank         = accepted_rank,
    interval     = early_interval,
    longitude    = lng,
    latitude     = lat,
    material     = composition_clean,
    material_grp = composition_group,
    lith1        = lith1_clean,
    lith2        = lith2_clean,
    lith1_grp    = lith1_group,
    lith2_grp    = lith2_group,
    lith_combo   = lith_combined,
    env          = env_clean,
    env_group    = environment_group
  )

cluster_ready <- mutated_reduced_data %>%
  select(latitude, longitude, lith_combo, material_grp, env_group)

cluster_ready <- cluster_ready %>%
  mutate(across(where(is.factor), as.character)) %>%
  mutate(across(where(is.character), ~ na_if(., "unknown"))) %>%
  mutate(across(where(is.character), as.factor))

complete_cluster_ready <- na.omit(cluster_ready)

# Optional: export cleaned dataset
# write_csv(complete_cluster_ready, "data/clusterReady_cleanCompleteNoDates.csv")

# -------------------- Read Cleaned Clustering File -------------------- #

df <- read.csv("data/clusterReady_clean.csv")

df <- df %>%
  mutate(
    lith_combo = as.factor(lith_combo),
    material_grp = as.factor(material_grp),
    env_group = as.factor(env_group)
  )

# -------------------- Sample Data for Clustering -------------------- #

set.seed(32)
df_sub <- df %>% sample_n(8000)

df_cluster <- df_sub %>%
  select(lith_combo, material_grp, env_group)

# -------------------- Compute Gower Distance -------------------- #

gower_dist <- daisy(df_cluster, metric = "gower")

# -------------------- Select Number of Clusters -------------------- #

sil_scores <- c()

for (k in 3:8) {
  cat("Fitting PAM with k =", k, "...\n")
  pam_fit <- pam(gower_dist, diss = TRUE, k = k)
  sil <- silhouette(pam_fit$clustering, gower_dist)
  mean_sil <- mean(sil[, 3])
  sil_scores[k] <- mean_sil
  cat("Silhouette score:", mean_sil, "\n\n")
}

best_k <- which.max(sil_scores)
final_pam <- pam(gower_dist, diss = TRUE, k = best_k)

df_sub$cluster <- as.factor(final_pam$clustering)

# -------------------- Add Geologic Period Labels -------------------- #

df_sub <- df_sub %>%
  mutate(
    periods = case_when(
      max_ma <= 541 & min_ma >= 485 ~ "Cambrian",
      max_ma <= 485 & min_ma >= 444 ~ "Ordovician",
      max_ma <= 444 & min_ma >= 419 ~ "Silurian",
      max_ma <= 419 & min_ma >= 359 ~ "Devonian",
      max_ma <= 359 & min_ma >= 299 ~ "Carboniferous",
      max_ma <= 299 & min_ma >= 252 ~ "Permian",
      max_ma <= 252 & min_ma >= 201 ~ "Triassic",
      max_ma <= 201 & min_ma >= 145 ~ "Jurassic",
      max_ma <= 145 & min_ma >= 66  ~ "Cretaceous",
      max_ma <= 66  & min_ma >= 23  ~ "Paleogene",
      max_ma <= 23  & min_ma >= 2.6 ~ "Neogene",
      max_ma <= 2.6 & min_ma >= 0   ~ "Quaternary"
    )
  )

# -------------------- Cluster Summary -------------------- #

cluster_summary <- df_sub %>%
  group_by(cluster) %>%
  summarize(
    count = n(),
    top_env = names(sort(table(env_group), decreasing = TRUE))[1],
    top_lith = names(sort(table(lith_combo), decreasing = TRUE))[1],
    top_material = names(sort(table(material_grp), decreasing = TRUE))[1],
    lat_min = min(latitude),
    lat_max = max(latitude),
    lon_min = min(longitude),
    lon_max = max(longitude)
  )

print(cluster_summary)

# -------------------- Plot Settings -------------------- #

my_colors <- c("lightgoldenrod", "thistle", "plum", "darkseagreen", "cornflowerblue")

# -------------------- Save Figures -------------------- #

png("figures/clusters_map.png", width = 900, height = 650)
ggplot(df_sub, aes(x = longitude, y = latitude, color = cluster)) +
  geom_point(alpha = 0.5, size = 1) +
  theme_minimal() +
  scale_color_manual(values = my_colors) +
  labs(
    title = "Clusters of Fossil Localities",
    x = "Longitude",
    y = "Latitude"
  )
dev.off()

png("figures/silhouette_plot.png", width = 900, height = 650)
par(mar = c(10, 5, 5.5, 6.5))
sil <- silhouette(final_pam$clustering, dist = gower_dist)
plot(sil, col = my_colors, border = NA, main = "Silhouette Plot")
dev.off()

png("figures/period_composition_by_cluster.png", width = 900, height = 650)
ggplot(df_sub, aes(x = cluster, fill = periods)) +
  geom_bar(position = "fill") +
  theme_minimal() +
  labs(
    title = "Period Composition by Cluster",
    y = "Proportion",
    x = "Cluster"
  )
dev.off()

png("figures/environment_composition_by_cluster.png", width = 900, height = 650)
ggplot(df_sub, aes(x = cluster, fill = env_group)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal() +
  labs(
    title = "Environment Composition by Cluster",
    y = "Percent",
    x = "Cluster"
  )
dev.off()

png("figures/lithology_composition_by_cluster.png", width = 900, height = 650)
ggplot(df_sub, aes(x = cluster, fill = lith_combo)) +
  geom_bar(position = "fill") +
  theme_minimal() +
  labs(
    title = "Lithology Composition by Cluster",
    y = "Proportion",
    x = "Cluster"
  )
dev.off()

png("figures/material_composition_by_cluster.png", width = 900, height = 650)
ggplot(df_sub, aes(x = cluster, fill = material_grp)) +
  geom_bar(position = "fill") +
  theme_minimal() +
  labs(
    title = "Fossil Material Composition by Cluster",
    y = "Proportion",
    x = "Cluster"
  )
dev.off()

png("figures/cluster5_map.png", width = 900, height = 650)
ggplot(
  df_sub %>% filter(cluster == 5),
  aes(x = longitude, y = latitude, color = material_grp)
) +
  geom_point(alpha = 0.7, size = 1.2) +
  theme_minimal() +
  labs(
    title = "Cluster 5 Localities by Fossil Material",
    x = "Longitude",
    y = "Latitude",
    color = "Material"
  )
dev.off()

# -------------------- Cluster 5 Follow-Up Analysis -------------------- #

cluster5 <- df_sub %>% filter(cluster == 5)

env_comp <- cluster5 %>%
  count(env_group) %>%
  mutate(percent = round((n / sum(n)) * 100, 2))

lith_comp <- cluster5 %>%
  count(lith_combo) %>%
  mutate(percent = round((n / sum(n)) * 100, 2))

material_comp <- cluster5 %>%
  count(material_grp) %>%
  mutate(percent = round((n / sum(n)) * 100, 2))

cluster5_composition <- list(
  environment = env_comp,
  lithology = lith_comp,
  material = material_comp
)

print(cluster5_composition)

cluster5_long <- cluster5 %>%
  pivot_longer(
    cols = c(env_group, lith_combo, material_grp),
    names_to = "factor",
    values_to = "category"
  ) %>%
  count(factor, category) %>%
  group_by(factor) %>%
  mutate(percent = round((n / sum(n)) * 100, 2)) %>%
  arrange(factor, desc(percent))

print(cluster5_long)