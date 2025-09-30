# =============================
# Step 0: Load libraries
# =============================
library(openxlsx)
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)

# =============================
# Step 1: Import dataset
# =============================
setwd("#File_directory")
df <- read_excel("allelo_dataset_all_parameters.xlsx", sheet = "pot")

# =============================
# Step 2: Fill missing values
# =============================
df_filled <- df %>%
  group_by(species, treatment) %>%
  mutate(across(c(dry_weight_gm,nodule_number, gp,GI,
                  h2o2_nmol_g_FW, total_chlorophyll_content_mg_g,
                  carotenoid_mg_g, MDA_nmol_g_FW, proline_umol_g_FW),
                ~ ifelse(is.na(.), mean(., na.rm = TRUE), .))) %>%
  ungroup() %>%
  mutate(
    species = as.factor(species),
    treatment = as.factor(treatment)
  )

# =============================
# Step 3: Function to calculate RI & SE + Heatmap
# =============================
calc_RI_SE_heatmap <- function(df, traits, SE_traits, group_vars, control_label = "1") {
  
  # Summarise mean values per group
  df_summary <- df %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(across(all_of(traits), mean, na.rm = TRUE), .groups = "drop")
  
  # Extract control means
  control_means <- df_summary %>%
    filter(treatment == control_label) %>%
    select(all_of(group_vars), all_of(traits))
  
  # Join control values to each row
  join_vars <- setdiff(group_vars, "treatment")
  df_joined <- df_summary %>%
    left_join(control_means, by = join_vars, suffix = c("", "_control"))
  
  # Calculate RI
  for (t in traits) {
    RI_col <- paste0("RI_", t)
    control_col <- paste0(t, "_control")
    df_joined[[RI_col]] <- ifelse(df_joined[[t]] >= df_joined[[control_col]],
                                  1 - df_joined[[control_col]] / df_joined[[t]],
                                  df_joined[[t]] / df_joined[[control_col]] - 1)
  }
  
  # Calculate SE (average of selected RI traits)
  RI_SE_cols <- paste0("RI_", SE_traits)
  df_joined$SE <- rowMeans(df_joined[, RI_SE_cols], na.rm = TRUE)
  
  # =============================
  # Step 3a: Improve labels
  # =============================
  df_joined <- df_joined %>%
    mutate(
      species_label = case_when(
        species == "A" ~ "Cajanus cajan",
        species == "M" ~ "Phaseolus mungo",
        species == "B" ~ "Vigna unguiculata",
        TRUE ~ as.character(species)
      ),
      treatment_label = case_when(
        treatment == "1"   ~ "Control",
        treatment == "2" ~ "Leaf Powder | 100:5",
        treatment == "3" ~ "Leaf Powder | 100:10",
        treatment == "4" ~ "Leaf Powder | 100:20",
        treatment == "5" ~ "Leaf Water Extract | 25 mg/ml",
        treatment == "6" ~ "Leaf Water Extract | 50 mg/ml",
        treatment == "7" ~ "Leaf Water Extract | 100 mg/ml",
        TRUE ~ as.character(treatment)
      ),
      group_label = paste(species_label, treatment_label, sep = " — ")
    )
  
  # =============================
  # Step 3b: Prepare long format
  # =============================
  df_long <- df_joined %>%
    select(group_label, all_of(paste0("RI_", traits)), SE) %>%
    pivot_longer(cols = -group_label, names_to = "Index", values_to = "Value")
  
  # Nice names for indices
  # Step 3g: Custom index order
  index_order <- c(
    # Germination
    "RI_gp", "RI_GI", "RI_GE", "RI_GS",
    # Morphology
    "RI_shoot_lentgh_cm", "RI_root_length_cm", "RI_dry_weight_gm", "RI_nodule_number",
    # Biochemical
    "RI_total_chlorophyll_content_mg_g", "RI_carotenoid_mg_g",
    "RI_MDA_nmol_g_FW", "RI_h2o2_nmol_g_FW", "RI_proline_umol_g_FW",
    # Synthetic
    "SE"
  )
  df_long$Index <- factor(df_long$Index, levels = index_order)
  
  index_labels <- c(
    "RI_root_length_cm" = "Root Length",
    "RI_shoot_lentgh_cm" = "Shoot Length",
    "RI_GI" = "Germination Index",
    "RI_gp" = "Germination %",
    "RI_h2o2_nmol_g_FW" = expression(H[2]*O[2]~"(nmol/g FW)"),
    "RI_MDA_nmol_g_FW" = "MDA (nmol/g FW)",
    "RI_proline_umol_g_FW" = "Proline (µmol/g FW)",
    "RI_total_chlorophyll_content_mg_g" = "Chlorophyll (mg/g)",
    "RI_carotenoid_mg_g" = "Carotenoid (mg/g)",
    "RI_dry_weight_gm" = " Biomass (g)",
    "SE" = "Synthetic Effect (SE)"
  )
  
  # Order y-axis: species grouped, treatments ordered
  treatment_order <- c("Control",
                       "Leaf Powder | 100:5", "Leaf Powder | 100:10", "Leaf Powder | 100:20",
                       "Leaf Water Extract | 25 mg/ml", "Leaf Water Extract | 50 mg/ml", "Leaf Water Extract | 100 mg/ml")
  
  df_long <- df_long %>%
    mutate(
      treatment_ordered = factor(sub(".*— ", "", group_label), levels = treatment_order),
      species_ordered = factor(sub(" —.*", "", group_label),
                               levels = c("Cajanus cajan", "Phaseolus mungo", "Vigna unguiculata")),
      group_label = fct_inorder(paste(species_ordered, treatment_ordered, sep = " — "))
    )
  
  # =============================
  # Step 3c: Plot heatmap
  # =============================
  heatmap_plot <- ggplot(df_long, aes(x = Index, y = group_label, fill = Value)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", Value),
                  fontface = ifelse(abs(Value) > 0.5, "bold", "plain"),
                  color = abs(Value) > 0.5), size = 3) +
    scale_color_manual(values = c("TRUE" = "white", "FALSE" = "black"), guide = "none") +
    scale_fill_gradient2(low = "#D73027", mid = "white", high = "#1A9850", midpoint = 0, name = "RI / SE") +
    scale_x_discrete(labels = index_labels) +
    scale_y_discrete(limits = rev) +
    labs(x = "Index", y = "Species × Treatment",
         title = "Allelopathic Response (RI & SE)") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1, size = 10),
      axis.text.y = element_text(size = 9, face = "italic"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      panel.grid = element_blank()
    )
  
  return(list(results = df_joined, heatmap = heatmap_plot))
}

# =============================
# Step 4: Run analysis for Pot dataset
# =============================
RI_traits_pot <- c("shoot_lentgh_cm","root_length_cm","dry_weight_gm", "gp", "GI",
                   "h2o2_nmol_g_FW", "total_chlorophyll_content_mg_g",
                   "carotenoid_mg_g", "MDA_nmol_g_FW", "proline_umol_g_FW")

SE_traits_pot <- c("shoot_lentgh_cm","root_length_cm","dry_weight_gm", "gp", "GI",
                   "h2o2_nmol_g_FW", "total_chlorophyll_content_mg_g",
                   "carotenoid_mg_g", "MDA_nmol_g_FW", "proline_umol_g_FW")

res_pot <- calc_RI_SE_heatmap(df_filled,
                              traits = RI_traits_pot,
                              SE_traits = SE_traits_pot,
                              group_vars = c("species", "treatment"),
                              control_label = "1")

# View results
#View(res_pot$results)

# Show heatmap
print(res_pot$heatmap)
ggsave("res_pot$heatmap.png", plot = res_pot$heatmap, dpi = 600, width = 8, height = 5, units = "in", bg = "white")

