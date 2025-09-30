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
df <- read_excel("allelo_dataset_all_parameters.xlsx", sheet = "petri")

# =============================
# Step 2: Fill missing values
# =============================
df_filled <- df %>%
  group_by(extract_type, species, conc) %>%
  mutate(across(c(seeds, gp, D4, D7, GI),
                ~ ifelse(is.na(.), mean(., na.rm = TRUE), .))) %>%
  ungroup() %>%
  mutate(
    species = as.factor(species),
    extract_type = as.factor(extract_type),
    conc = as.factor(conc)
  )

# =============================
# Step 3: Function to calculate RI & SE + heatmap
# =============================
calc_RI_SE_heatmap <- function(df, traits, SE_traits, group_vars, control_label = "1") {
  
  # Summarize mean per group
  df_summary <- df %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(across(all_of(traits), mean, na.rm = TRUE), .groups = "drop")
  
  # Extract control means
  control_means <- df_summary %>%
    filter(conc == control_label) %>%
    select(all_of(group_vars), all_of(traits))
  
  # Join control to each row
  join_vars <- setdiff(group_vars, "conc")
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
  
  # Calculate SE
  RI_SE_cols <- paste0("RI_", SE_traits)
  df_joined$SE <- rowMeans(df_joined[, RI_SE_cols], na.rm = TRUE)
  
  # =============================
  # Improve labels
  # =============================
  df_joined <- df_joined %>%
    mutate(
      species_label = case_when(
        species == "A" ~ "Cajanus cajan",
        species == "M" ~ "Phaseolus mungo",
        species == "B" ~ "Vigna unguiculata",
        TRUE ~ as.character(species)
      ),
      extract_label = case_when(
        extract_type == "akw"   ~ "Leaf Water Extract",
        extract_type == "ak50"  ~ "50% Ethanol Extract",
        extract_type == "ak100" ~ "100% Ethanol Extract",
        TRUE ~ as.character(extract_type)
      ),
      conc_label = case_when(
        conc == "1" ~ "Control",
        conc == "2" ~ "25 mg/ml",
        conc == "3" ~ "50 mg/ml",
        conc == "4" ~ "100 mg/ml",
        TRUE ~ as.character(conc)
      ),
      # For y-axis: species first, then extract, then conc
      group_label = paste(species_label, extract_label, conc_label, sep = " — ")
    )
  
  # =============================
  # Prepare long format for heatmap
  # =============================
  df_long <- df_joined %>%
    select(group_label, all_of(paste0("RI_", traits)), SE) %>%
    pivot_longer(cols = -group_label, names_to = "Index", values_to = "Value")
  
  # Index order and labels
  index_order <- c("RI_gp","RI_GI",  "RI_plumule_lentgh_cm", "RI_radicle_lentgh_cm","RI_fresh_biomass_gm","SE")
  df_long$Index <- factor(df_long$Index, levels = index_order)
  
  index_labels <- c(
    "RI_gp" = "Germination %",
     "RI_GI" = "Germination Index",
    "RI_plumule_lentgh_cm" = "Plumule Length (cm)",
    "RI_radicle_lentgh_cm" = "Radicle Length (cm)",
    "RI_fresh_biomass_gm" = "Fresh Biomass (g)",
    "SE" = "Synthetic Effect"
  )
  
  # Order y-axis: species -> extract -> conc
  species_levels <- c("Cajanus cajan", "Phaseolus mungo", "Vigna unguiculata")
  extract_levels <- c("Leaf Water Extract","50% Ethanol Extract", "100% Ethanol Extract" )
  conc_levels <- c("Control","25 mg/ml","50 mg/ml","100 mg/ml")
  
  df_long <- df_long %>%
    mutate(
      species_ordered = factor(sub(" —.*", "", group_label), levels = species_levels),
      extract_ordered = factor(sub(".*— (.*) — .*", "\\1", group_label), levels = extract_levels),
      conc_ordered = factor(sub(".*— .* — (.*)", "\\1", group_label), levels = conc_levels),
      group_label = fct_inorder(paste(species_ordered, extract_ordered, conc_ordered, sep = " — "))
    )
  
  # =============================
  # Plot heatmap
  # =============================
  heatmap_plot <- ggplot(df_long, aes(x = Index, y = group_label, fill = Value)) +
    geom_tile(color = "white") +
    geom_text(aes(
      label = ifelse(is.na(Value), "", sprintf("%.2f", Value)),
      fontface = ifelse(is.na(Value), "plain", ifelse(abs(Value) > 0.5, "bold", "plain")),
      color = ifelse(is.na(Value), "black", ifelse(abs(Value) > 0.5, "white", "black"))
    ), size = 3) +
    scale_color_manual(values = c("white" = "white", "black" = "black"), guide = "none") +
    scale_fill_gradient2(low = "#D73027", mid = "white", high = "#1A9850", midpoint = 0, name = "RI / SE") +
    scale_x_discrete(labels = index_labels) +
    scale_y_discrete(limits = rev) +
    labs(x = "Index", y = "Species — Extract — Concentration",
         title = "Allelopathic Response (RI & SE) - Petri Dish") +
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
# Step 4: Run analysis for Petri Dish dataset
# =============================
RI_traits_petri <- c("GI", "gp", "plumule_lentgh_cm", "radicle_lentgh_cm","fresh_biomass_gm")
SE_traits_petri <- RI_traits_petri

res_petri <- calc_RI_SE_heatmap(df_filled,
                                traits = RI_traits_petri,
                                SE_traits = SE_traits_petri,
                                group_vars = c("species", "extract_type", "conc"),
                                control_label = "1")

# =============================
# Step 5: View results and heatmap
# =============================
#View(res_petri$results)
print(res_petri$heatmap)
ggsave("res_petri_heatmap.png", plot = res_petri$heatmap, dpi = 600, width = 8, height = 6, units = "in", bg = "white")
