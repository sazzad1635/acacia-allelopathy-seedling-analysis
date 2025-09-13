
🌱 Allelopathy Seedling Analysis

This repository contains the R scripts, datasets, and analysis pipeline used to evaluate the allelopathic effects of Acacia auriculiformis leaf extracts on crop seedlings. The study examines germination, morphological, and biochemical responses under different treatments.

📂 Repository Contents

data/ – Raw and processed datasets (Excel/CSV).

scripts/ – R scripts for data preprocessing, statistical analysis, and visualization.

results/ – Generated figures (PCA plots, heatmaps, barplots, etc.).

README.md – Documentation of repository usage.

🔬 Analysis Overview

Software: MS Excel, RStudio (Version 2025.05.1)

R Packages:

ggplot2 (3.5.2) – Data visualization

dplyr (1.1.4) – Data wrangling

tidyr (1.3.1) – Data tidying

multcompView (0.1-10) – Post-hoc analysis visualization

Statistical Methods:

One-way ANOVA for treatment effects

TukeyHSD test for pairwise comparisons

Heatmap analysis for allelopathic response index (RI) and synthesis effect index (SE)

📊 Parameters Analyzed

Seedling Germination Index

Morphology: shoot length, root length, fresh biomass

Biochemical Traits: chlorophyll, carotenoid, proline, MDA, H₂O₂

🚀 How to Use

Clone the repository:

git clone https://github.com/sazzad1635/allelopathy-seedling-analysis.git
cd allelopathy-seedling-analysis


Open scripts/ in RStudio.

Run scripts in order (data preprocessing → statistical analysis → visualization).
