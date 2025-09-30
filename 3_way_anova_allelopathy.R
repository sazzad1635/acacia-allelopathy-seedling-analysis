# Load required libraries
library(multcomp)
library(broom)
library(openxlsx)
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(multcompView)
library(FactoMineR)
library(factoextra)


# Set working directory and load data
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




# Drop unused combinations to prevent post-hoc errors
df_filled <- droplevels(df_filled)

# Rename biochar levels for legend clarity
levels(df_filled$species) <- c("A"="Cajanus cajan","B"= "Vigna unguiculata","M"="Phaseolus mungo")
levels(df_filled$extract_type) <- c("akw"="Leaf Water Extract","ak50"= "50% Ethanol Leaf Extract","ak100"= "100% Ethanol Leaf Extract")
levels(df_filled$conc) <- c("1"="Control","2"= "25 mg/ml","3"= "50 mg/ml","4"="100 mg/ml")


# ------------------------
# Three-way ANOVA for plumule length
# ------------------------
anova_plumule_lentgh_cm <- aov(plumule_lentgh_cm ~ species * extract_type * conc, data = df_filled)
anova_result_plumule <- summary(anova_plumule_lentgh_cm)
print(anova_result_plumule)

# ------------------------


# ------------------------
# Three-way ANOVA for radicle length
# ------------------------
anova_radicle_lentgh_cm <- aov(radicle_lentgh_cm ~ species * extract_type * conc, data = df_filled)
anova_result_radicle <- summary(anova_radicle_lentgh_cm)
print(anova_result_radicle)


# ------------------------
# ------------------------
# Three-way ANOVA for germination percent
# ------------------------
anova_gp <- aov(gp ~ species * extract_type * conc, data = df_filled)
anova_result_gp <- summary(anova_gp)
print(anova_result_gp)

# Three-way ANOVA for germination index
# ------------------------
anova_gi <- aov(GI ~ species * extract_type * conc, data = df_filled)
anova_result_gi <- summary(anova_gi)
print(anova_result_gi)

fresh_biomass_gm
# Three-way ANOVA for biomass
# ------------------------
anova_fresh_biomass_gm <- aov(fresh_biomass_gm ~ species * extract_type * conc, data = df_filled)
anova_result_fresh_biomass_gm <- summary(anova_fresh_biomass_gm)
print(anova_result_fresh_biomass_gm)
# ------------------------
#######################
##################pot##############################################################
# -----------------------------
# Load required libraries
# -----------------------------
library(readxl)
library(dplyr)

# -----------------------------
# Step 1: Read dataset
# -----------------------------
df <- read_excel("#file_directory", 
                 sheet = "pot")

# -----------------------------
# Step 2: Fill missing values within existing groups
# -----------------------------
df_filled <- df %>%
  group_by(species, treatment, extract_type) %>%
  mutate(across(c(dry_weight_gm, nodule_number, gp,  GI,
                  h2o2_nmol_g_FW, total_chlorophyll_content_mg_g,
                  carotenoid_mg_g, MDA_nmol_g_FW, proline_umol_g_FW),
                ~ ifelse(is.na(.), mean(., na.rm = TRUE), .))) %>%
  ungroup() %>%
  mutate(
    species = as.factor(species),
    treatment = as.factor(treatment),
    extract_type = as.factor(extract_type)
  )

# -----------------------------
# Step 3: Drop rows where key response variable is still NA
# -----------------------------
df_clean <- df_filled %>%
  filter(!is.na(fresh_biomass_gm)) %>%
  droplevels()

# -----------------------------
# Step 4: Rename factor levels for clarity
# -----------------------------
levels(df_clean$species) <- c("A"="Cajanus cajan","B"= "Vigna unguiculata","M"="Phaseolus mungo")
levels(df_clean$extract_type) <- c("LP"="Leaf Powder", "LW"="Leaf Water Extract")
levels(df_clean$treatment) <- c("1"="Control", "2"="100:5", "3"="100:10", "4"="100:20",
                                "5"="25 mg/ml", "6"="50 mg/ml","7"= "100 mg/ml")

# -----------------------------
# Step 5: Three-way ANOVA for shoot length etc.
# -----------------------------
anova_fresh_biomass_gm <- aov(fresh_biomass_gm~ species * extract_type + species * treatment %in% extract_type, 
                             data = df_clean)

summary(anova_fresh_biomass_gm)



