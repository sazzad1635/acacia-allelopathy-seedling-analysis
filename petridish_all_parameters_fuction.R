library(openxlsx)
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(multcompView)
library(FactoMineR)
library(factoextra)

# =============================
# Step 1: Import dataset
# =============================
setwd("#file_directory")

df <- read_excel("allelo_dataset_all_parameters.xlsx", sheet = "petri")

# Step 2: Fill missing values with group mean (per extract_type + species + conc)
df_filled <- df %>%
  group_by(extract_type, species, conc) %>%
  mutate(across(c(seeds, D4, D7,gp, GI),
                ~ ifelse(is.na(.), mean(., na.rm = TRUE), .))) %>%
  ungroup() %>%
  mutate(
    species = as.factor(species),
    extract_type = as.factor(extract_type),
    conc = as.factor(conc)   # treat concentration as categorical factor
  )

# =============================
# Helper: Extract significance letters
# =============================
get_significance_letters <- function(anova_model, tukey_res) {
  letters_obj <- multcompLetters4(anova_model, tukey_res)
  data.frame(
    conc = names(letters_obj$conc$Letters),
    Letters = letters_obj$conc$Letters,
    stringsAsFactors = FALSE
  )
}

# =============================
# General Function for Extract & Variable
# =============================
analyze_extract <- function(data, extract_code, response_var, title_text, subtitle_text, y_label) {
  
  # Subset for this extract type
  df_sub <- data %>% filter(extract_type == extract_code)
  
  # Formula dynamically (gp, GI, etc.)
  f <- as.formula(paste(response_var, "~ conc"))
  
  # Run ANOVA + Tukey for each species
  run_anova <- function(df) {
    if(nrow(df) > 0) {
      model <- aov(f, data = df)
      tukey <- TukeyHSD(model)
      letters <- get_significance_letters(model, tukey)
      return(list(model = model, letters = letters))
    }
    return(NULL)
  }
  
  res_A <- run_anova(df_sub %>% filter(species == "A"))
  res_B <- run_anova(df_sub %>% filter(species == "B"))
  res_M <- run_anova(df_sub %>% filter(species == "M"))
  
  # Collect letters
  letters_all <- bind_rows(
    if(!is.null(res_A)) mutate(res_A$letters, species = "A"),
    if(!is.null(res_B)) mutate(res_B$letters, species = "B"),
    if(!is.null(res_M)) mutate(res_M$letters, species = "M")
  )
  
  # =============================
  # Summary Data
  # =============================
  summary_df <- df_sub %>%
    group_by(species, conc) %>%
    summarize(mean_val = mean(.data[[response_var]], na.rm = TRUE),
              se_val = sd(.data[[response_var]], na.rm = TRUE) / sqrt(n()), .groups = "drop") %>%
    left_join(letters_all, by = c("species", "conc"))
  
  # Dynamic offset for labels
  y_offset <- (max(summary_df$mean_val + summary_df$se_val) - 
                 min(summary_df$mean_val - summary_df$se_val)) * 0.05
  
  summary_df <- summary_df %>%
    mutate(label_y = mean_val + se_val + y_offset)
  
  y_max <- max(summary_df$label_y) * 1.1
  
  # =============================
  # Plot
  # =============================
  p <- ggplot(summary_df, aes(x = species, y = mean_val, fill = conc)) + 
    geom_bar(stat = "identity", color = "black", 
             position = position_dodge(), width = 0.7) +
    scale_x_discrete(breaks = c("A", "B", "M"),
                     labels = c("Cajanus cajan", "Vigna unguiculata", "Phaseolus mungo")) +
    geom_errorbar(aes(ymin = mean_val - se_val, ymax = mean_val + se_val), 
                  width = 0.2, position = position_dodge(0.7)) +
    geom_text(aes(label = Letters, y = label_y), 
              position = position_dodge(0.7), vjust = -0.5, size = 5) +
    guides(fill = guide_legend(title = "Concentration")) +
    scale_fill_manual(values = c("1" = "#4F6F52", 
                                 "2" = "#739072",
                                 "3" = "#86A789",
                                 "4" = "#D2E3C8"),
                      labels = c("1" = "Control",
                                 "2" = "25 mg/ml",
                                 "3" = "50 mg/ml",
                                 "4" = "100 mg/ml")) +
    labs(title = title_text,
         subtitle = subtitle_text,
         y = y_label,
         x = "Species") +
    theme_minimal(base_size = 14) +
    theme(
      panel.background = element_blank(),
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(color = "black",size=12,face = "bold.italic"),
      axis.line = element_line(color = "black", size = 0.2),
      axis.ticks.y = element_line(color = "black", size = 0.5),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12, face = "bold"),
      axis.ticks.length = unit(6, "pt"),
      legend.position = "right",
      legend.key.size = unit(0.4, "cm"),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(hjust = 0, face = "bold", size = 14),
      aspect.ratio = 5/8
    ) +
    scale_y_continuous(
      limits = c(0, y_max),
      breaks = pretty(c(0, y_max)),
      expand = expansion(mult = c(0.02, 0.05))
    )
  
  return(p)
}

# =============================
# Example usage
# =============================

# Germination % (gp)

#View(df_filled)
head(df_filled)
gp_akw <- analyze_extract(df_filled, "akw", "gp", "Leaf Water Extract", "Petridish Experiment", "Germination (%)")
gp_ak50 <- analyze_extract(df_filled, "ak50", "gp", "50% Ethanol Extract", "Petridish Experiment", "Germination (%)")
gp_ak100 <- analyze_extract(df_filled, "ak100", "gp", "100% Ethanol Extract", "Petridish Experiment", "Germination (%)")
plu_akw <- analyze_extract(df_filled, "akw", "plumule_lentgh_cm", "Leaf Water Extract", "Petridish Experiment", "Plumule length (cm)")
plu_ak50 <- analyze_extract(df_filled, "ak50", "plumule_lentgh_cm", "50% Ethanol Extract", "Petridish Experiment", "Plumule length (cm)")
plu_ak100 <- analyze_extract(df_filled, "ak100", "plumule_lentgh_cm", "100% Ethanol Extract", "Petridish Experiment", "Plumule length (cm)")
rad_akw <- analyze_extract(df_filled, "akw", "radicle_lentgh_cm", "Leaf Water Extract", "Petridish Experiment", "Radicle length (cm)")
rad_ak50 <- analyze_extract(df_filled, "ak50", "radicle_lentgh_cm", "50% Ethanol Extract", "Petridish Experiment", "Radicle length (cm)")
rad_ak100 <- analyze_extract(df_filled, "ak100", "radicle_lentgh_cm", "100% Ethanol Extract", "Petridish Experiment", "Radicle length (cm)")
bio_akw <- analyze_extract(df_filled, "akw", "fresh_biomass_gm", "Leaf Water Extract", "Petridish Experiment", " Biomass (g/dish)")
bio_ak50 <- analyze_extract(df_filled, "ak50", "fresh_biomass_gm", "50% Ethanol Extract", "Petridish Experiment", "Biomass (g/dish)")
bio_ak100 <- analyze_extract(df_filled, "ak100", "fresh_biomass_gm", "100% Ethanol Extract", "Petridish Experiment", "Biomass (g/dish)")


# Print one plot to test
print(gp_akw)
print(gp_ak50)
print(gp_ak100)
print(plu_akw)
print(plu_ak50)
print(plu_ak100)
print(rad_akw)
print(rad_ak50)
print(rad_ak100)
print(bio_akw)
print(bio_ak50)
print(bio_ak100)

ggsave("gp_akw.png", plot = gp_akw, dpi = 600, width = 8, height = 5, units = "in", bg = "white")
ggsave("gp_ak50.png", plot = gp_ak50, dpi = 600, width = 8, height = 5, units = "in", bg = "white")
ggsave("gp_ak100.png", plot = gp_ak100, dpi = 600, width = 8, height = 5, units = "in", bg = "white")
ggsave("plu_akw.png", plot = plu_akw, dpi = 600, width = 8, height = 5, units = "in", bg = "white")
ggsave("plu_ak50.png", plot = plu_ak50, dpi = 600, width = 8, height = 5, units = "in", bg = "white")
ggsave("plu_ak100.png", plot = plu_ak100, dpi = 600, width = 8, height = 5, units = "in", bg = "white")
ggsave("rad_akw.png", plot = rad_akw, dpi = 600, width = 8, height = 5, units = "in", bg = "white")
ggsave("rad_ak50.png", plot = rad_ak50, dpi = 600, width = 8, height = 5, units = "in", bg = "white")
ggsave("rad_ak100.png", plot = rad_ak100, dpi = 600, width = 8, height = 5, units = "in", bg = "white")
ggsave("bio_akw.png", plot = bio_akw, dpi = 600, width = 8, height = 5, units = "in", bg = "white")
ggsave("bio_ak50.png", plot = bio_ak50, dpi = 600, width = 8, height = 5, units = "in", bg = "white")
ggsave("bio_ak100.png", plot = bio_ak100, dpi = 600, width = 8, height = 5, units = "in", bg = "white")




######################################PCA##############################
# =============================
# PCA Analysis
# =============================

# Select numeric variables for PCA
PCA_numeric <- df_filled %>%
  select(radicle_lentgh_cm,
         gp, GI)

# Drop constant/NA columns 
PCA_numeric_clean <- PCA_numeric %>% select(where(~ var(., na.rm = TRUE) > 0))

# Run PCA
pca_result <- prcomp(PCA_numeric_clean, scale. = TRUE)

# Explained variance
pca_var <- pca_result$sdev^2
pca_var_percent <- round(100 * pca_var / sum(pca_var), 1)

# =============================
# PCA Plot - by Species
# =============================
# Relabel species
df_filled$species <- factor(df_filled$species,
                            levels = c("A", "B", "M"),
                            labels = c("Cajanus cajan",
                                       "Vigna unguiculata",
                                       "Phaseolus mungo"))

# Relabel extract_type
df_filled$extract_type <- factor(df_filled$extract_type,
                                 levels = c("akw", "ak50", "ak100"),
                                 labels = c("Leaf Water Extract",
                                            "50% Ethanol Extract",
                                            "100% Ethanol Extract"))
pca_plot_species <- fviz_pca_ind(
  pca_result,
  geom.ind = "point",
  col.ind = df_filled$species,   # color by species
  addEllipses = TRUE,            # confidence ellipses
  legend.title = "Species",
  palette = "Set2"
) +
  labs(title = "Comparative Allelopathic Effects Across Species",
       subtitle = "Petridish Experiment",
       x = paste0("PC1 (", pca_var_percent[1], "%)"),
       y = paste0("PC2 (", pca_var_percent[2], "%)")) +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(face = "italic"),
    axis.line = element_line(color = "black", size = 0.2),
    axis.ticks.y = element_line(color = "black", size = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),
    axis.ticks.length = unit(6, "pt"),
    legend.position = "right",
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 12,face="italic"),
    legend.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(hjust = 0, face = "bold", size = 14),
    aspect.ratio = 5/6
    )

print(pca_plot_species)
ggsave("pca_plot_species.png", plot = pca_plot_species, dpi = 600, width = 8, height = 5, units = "in", bg = "white")


# =============================
# PCA Plot - by Extract Type
# =============================
pca_plot_extract <- fviz_pca_ind(
  pca_result,
  geom.ind = "point",
  col.ind = df_filled$extract_type,  # color by extract type
  addEllipses = TRUE,
  legend.title = "Extract Type",
  palette = "Dark2"
) +
  labs(title = "PCA of Allelopathy Effects by Extract Type",
       x = paste0("PC1 (", pca_var_percent[1], "%)"),
       y = paste0("PC2 (", pca_var_percent[2], "%)")) +
  theme_minimal(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12)
  )

print(pca_plot_extract)
ggsave("pca_plot_extract.png", plot = pca_plot_extract, dpi = 600, width = 8, height = 5, units = "in", bg = "white")

# =============================
# Variables contribution plot
# =============================
fviz_pca_var(
  pca_result,
  col.var = "contrib",
  gradient.cols = c("blue", "orange", "red"),
  repel = TRUE
)



