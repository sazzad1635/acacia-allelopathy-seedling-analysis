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

df <- read_excel("allelo_dataset_all_parameters.xlsx", sheet = "pot")

# =============================
# Step 2: Fill missing values (species + treatment)
# =============================
df_filled <- df %>%
  group_by(species, treatment) %>%
  mutate(across(c(nodule_number,seeds, D4, D7,D14,gp, GI,
                  h2o2_nmol_g_FW,total_chlorophyll_content_mg_g,carotenoid_mg_g,
                 MDA_nmol_g_FW,proline_umol_g_FW,proline_mg_g_FW),
                ~ ifelse(is.na(.), mean(., na.rm = TRUE), .))) %>%
  ungroup() %>%
  mutate(
    species = as.factor(species),
    treatment = as.factor(treatment)
  )

# =============================
# Helper: Extract significance letters
# =============================
get_significance_letters <- function(anova_model, tukey_res) {
  letters_obj <- multcompLetters4(anova_model, tukey_res)
  data.frame(
    treatment = names(letters_obj$treatment$Letters),
    Letters = letters_obj$treatment$Letters,
    stringsAsFactors = FALSE
  )
}

# =============================
# General function for pot experiment (like petridish code)
# =============================
analyze_pot <- function(data, response_var, title_text, subtitle_text, y_label) {
  
  # Formula dynamically
  f <- as.formula(paste(response_var, "~ treatment"))
  
  # Run ANOVA per species
  run_anova <- function(df) {
    if(nrow(df) > 0) {
      model <- aov(f, data = df)
      tukey <- TukeyHSD(model)
      letters <- get_significance_letters(model, tukey)
      return(letters)
    }
    return(NULL)
  }
  
  letters_A <- run_anova(data %>% filter(species == "A")) %>% mutate(species = "A")
  letters_B <- run_anova(data %>% filter(species == "B")) %>% mutate(species = "B")
  letters_M <- run_anova(data %>% filter(species == "M")) %>% mutate(species = "M")
  
  letters_all <- bind_rows(
    letters_A, letters_B, letters_M
  )
  
  # =============================
  # Summary Data
  # =============================
  summary_df <- data %>%
    group_by(species, treatment) %>%
    summarize(mean_val = mean(.data[[response_var]], na.rm = TRUE),
              se_val = sd(.data[[response_var]], na.rm = TRUE)/sqrt(n()), .groups = "drop") %>%
    left_join(letters_all, by = c("species", "treatment"))
  
  # Dynamic offset for labels
  y_offset <- (max(summary_df$mean_val + summary_df$se_val) - 
                 min(summary_df$mean_val - summary_df$se_val)) * 0.05
  summary_df <- summary_df %>% mutate(label_y = mean_val + se_val + y_offset)
  y_max <- max(summary_df$label_y) * 1.1
  
  # =============================
  # Plot
  # =============================
  p <- ggplot(summary_df, aes(x = species, y = mean_val, fill = treatment)) +
    geom_bar(stat = "identity", color = "black", position = position_dodge(), width = 0.7) +
    scale_x_discrete(breaks = c("A", "B", "M"),
                     labels = c("Cajanus cajan", "Vigna unguiculata", "Phaseolus mungo")) +
    geom_errorbar(aes(ymin = mean_val - se_val, ymax = mean_val + se_val),
                  width = 0.2, position = position_dodge(0.7)) +
    geom_text(aes(label = Letters, y = label_y),
              position = position_dodge(0.7), vjust = -0.5, size = 5) +
    guides(fill = guide_legend(title = "Treatment")) +
    scale_fill_manual(values = c("1" = "#4F6F52", "2" = "#6A8A69", "3" = "#86A789",
                                 "4" = "#A9C3A8", "5" = "#C4DDBC", "6" = "#DEEED0",
                                 "7" = "#F1FAE3"),
                      labels = c("1" = "Control", "2" = "100:5", "3" = "100:10",
                                 "4" = "100:20", "5" = "25 mg/ml",
                                 "6" = "50 mg/ml", "7" = "100 mg/ml")) +
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
      aspect.ratio = 5/10
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
#gp_pot <- analyze_pot(df_filled, "gp", "Germination Percentage", "Pot Experiment", "Germination (%)")
bio_pot <- analyze_pot(df_filled, "fresh_biomass_gm", "Biomass", "Pot Experiment", "Biomass (g)")
proline_plot <- analyze_pot(df_filled, "proline_mg_g_FW", "Proline Content", "Pot Experiment", "Proline (mg/g FW)")
nod_pot <- analyze_pot(df_filled, "nodule_number", "Nodule Number", "Pot Experiment", "Nodule Number")
rot_pot <- analyze_pot(df_filled, "root_length_cm", "Root Length", "Pot Experiment", "Root Length (cm)")
shoot_pot <- analyze_pot(df_filled, "shoot_lentgh_cm", "Shoot Length", "Pot Experiment", "Shoot Length (cm)")
# Print plots
print(gp_pot)
print(bio_pot)
print(proline_plot)
print(rot_pot)
print(shoot_pot)

#ggsave("gp_pot.png", plot = gp_pot, dpi = 600, width = 10, height = 5, units = "in", bg = "white")
ggsave("rot_pot.png", plot = rot_pot, dpi = 600, width = 10, height = 5, units = "in", bg = "white")
ggsave("shoot_pot.png", plot = shoot_pot, dpi = 600, width = 10, height = 5, units = "in", bg = "white")
ggsave("bio_pot.png", plot = bio_pot, dpi = 600, width = 10, height = 5, units = "in", bg = "white")


################PCA#############
######################################PCA##############################
# =============================
# PCA Analysis
# =============================

# Select numeric variables for PCA
PCA_numeric <- df_filled %>%
  select(root_length_cm,shoot_lentgh_cm,dry_weight_gm,gp,GE, GI, GS,
         h2o2_nmol_g_FW,total_chlorophyll_content_mg_g,carotenoid_mg_g,
         MDA_nmol_g_FW,proline_mg_g_FW)
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

pca_plot_species <- fviz_pca_ind(
  pca_result,
  geom.ind = "point",
  col.ind = df_filled$species,   # color by species
  addEllipses = TRUE,            # confidence ellipses
  legend.title = "Species",
  palette = "Set2"
) +
  labs(title = "Comparative Allelopathic Effects Across Species",
       subtitle = "Pot Experiment",
       x = paste0("PC1 (", pca_var_percent[1], "%)"),
       y = paste0("PC2 (", pca_var_percent[2], "%)")) +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(face = "bold"),
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
# Variables contribution plot
# =============================
fviz_pca_var(
  pca_result,
  col.var = "contrib",
  gradient.cols = c("blue", "orange", "red"),
  repel = TRUE
)



