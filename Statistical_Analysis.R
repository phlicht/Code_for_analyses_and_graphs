                                              #############################################################
                                              ###        Statistical analyses and Plotting with R       ###
                                              ### Code shown in the order the respective plots/analyses ###
                                              ###               appear in the manuscript                ###
                                              #############################################################

###############################################
### Heatmap of Rarefied Metagenomic Samples ###
###                Figure 1                 ###
###############################################

#### Build the heatmap (language: R)
library(pheatmap)
library(tidyverse)
library(RColorBrewer)
library(ggplotify)

# load the metaphlan data
species_rarefied_heatmap <- read_tsv("rarefied_profiled_merged_species.tsv") %>%
    column_to_rownames("sample")
metadata <- read_tsv("metadata.tsv") %>%
    column_to_rownames("sample")
metadata$Patient <- factor(metadata$Patient, levels = c(
  "Pat1",
  "Pat2",
  "Pat3",
  "Pat4",
  "Pat5",
  "Pat6",
  "Pat14",
  "Pat15",
  "Pat17",
  "Pat18",
  "Pat20"
  )
  )

# Inspect how the data is structured 
x <- species_rarefied_heatmap %>% as.vector() %>%
    unlist
ggpubr::ggdensity(x)

# Set breaks for the color scale of the pheatmap legend (to increase the contrast in the heatmap)
q <- unique(
  quantile(x, 
  probs = c(
  0,  
  0.8,
  0.81,
  0.82,
  0.83,
  0.84,
  0.85,
  0.86,
  0.87,
  0.88,
  0.89,
  0.9,
  0.91,
  0.92,
  0.93,
  0.94,
  0.945,
  0.95,
  0.952,
  0.955,
  0.957,
  0.958,
  0.96,
  0.961,
  0.962,
  0.964,
  0.967,
  0.968,
  0.969,
  0.97,
  0.972,
  0.973,
  0.974,
  0.975,
  0.976,
  0.977,
  0.978,
  0.979,
  0.98,
  0.981,
  0.982,
  0.983,
  0.984,
  0.985,
  0.987,
  0.987,
  0.99,
  0.991,
  0.991,
  0.992,
  0.993,
  1
)))

breaks <- c(q, 101)

# Specify the number of color levels
num_colors <- length(q)

# Set color palette
colorPalette <- colorRampPalette(
  brewer.pal(n = 9, name = "YlOrRd"))(num_colors)
colorPalette[1] <- "white"

# Set Itacilized microbial species names
italic_species_names <- lapply(
  rownames(species_rarefied_heatmap),
  function(x) bquote(italic(.(x))))

# Plot the heatmap
heatmap_rarefied <-
pheatmap(
  mat = species_rarefied_heatmap,
  color = colorPalette,
  breaks = breaks,
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  cutree_cols = 3,
  treeheight_row = 0,
  treeheight_col = 25,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_col = metadata,
  labels_row = as.expression(italic_species_names),
  annotation_colors = ann_colors,
  display_numbers = FALSE
  )

#####################
### Shannon Index ###
###   Figure 2a   ###
#####################
library(vegan)
diversity(
  microbiome_rarefied_species,
  index = "shannon",
)

###################################
### Bray-Curtis Dissimilarities ###
###          Figure 2b          ###
###################################

# Calculate Bray-Curtis distances
library(vegan)
bray_species <- vegan::vegdist(
  t(feature_table_species_rarefied),
  method = "bray"
)

# Transform class dist into matrix
bray_species_matrix <- as.matrix(bray_species)
# Drop redudant and sel-comparing entries
bray_species_matrix[upper.tri(bray_species_matrix)] <- NA
diag(bray_species_matrix) <- NA
# Transform into long data.frame
bray_species_matrix <- melt(bray_species_matrix)
bray_species_matrix <- bray_species_matrix[complete.cases(bray_species_matrix$value), ]

# Select metadata needed
sd <- metadata %>%
  select(sample, stage) %>%
  mutate_if(is.factor, as.character)

# Combine distances with sample data
colnames(sd) <- c("Var1", "Stage1")
bray_species_matrix.sd <- left_join(bray_species_matrix, sd, by = "Var1")

colnames(sd) <- c("Var2", "Stage2")
bray_species_matrix.sd <- left_join(bray_species_matrix.sd, sd, by = "Var2")

# Add Labels
bray_species_matrix.sd$Label <-
  ifelse(bray_species_matrix.sd$Stage1 == "nonlesional" & bray_species_matrix.sd$Stage2 == "nonlesional", "Within nonlesional",
  ifelse(bray_species_matrix.sd$Stage1 == "Patch" & bray_species_matrix.sd$Stage2 == "Patch", "Within Patch",
  ifelse(bray_species_matrix.sd$Stage1 == "Plaque" & bray_species_matrix.sd$Stage2 == "Plaque", "Within Plaque",

  ifelse(bray_species_matrix.sd$Stage1 == "nonlesional" & bray_species_matrix.sd$Stage2 == "Patch", "Between nonlesional and Patch",
  ifelse(bray_species_matrix.sd$Stage1 == "Patch" & bray_species_matrix.sd$Stage2 == "nonlesional", "Between nonlesional and Patch",

  ifelse(bray_species_matrix.sd$Stage1 == "nonlesional" & bray_species_matrix.sd$Stage2 == "Plaque", "Between nonlesional and Plaque",
  ifelse(bray_species_matrix.sd$Stage1 == "Plaque" & bray_species_matrix.sd$Stage2 == "nonlesional", "Between nonlesional and Plaque",

  ifelse(bray_species_matrix.sd$Stage1 == "Patch" & bray_species_matrix.sd$Stage2 == "Plaque", "Between Patch and Plaque",
  ifelse(bray_species_matrix.sd$Stage1 == "Plaque" & bray_species_matrix.sd$Stage2 == "Patch", "Between Patch and Plaque", "")))))))))

# PERMANOVA
library(vegan)
set.seed(42)
adonis_species <- adonis(
  formula = bray_species ~ stage,
  data = metadata
  )
print(adonis_species$aov.tab)
# Permutation: free
# Number of permutations: 999
#
# Terms added sequentially (first to last)
#
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# stage      2    1.4018 0.70092  2.1364 0.12842  0.025 *
# Residuals 29    9.5144 0.32808         0.87158
# Total     31   10.9163                 1.00000
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Pairwise PERMANOVA
library(pairwiseAdonis)
set.seed(42)
pairwise_adonis_species <- pairwise.adonis2(
  bray_species ~ stage,
  data = metadata
)
print(pairwise_adonis_species)
# $parent_call
# [1] "bray_species ~ stage , strata = Null , permutations 999"
#
# $Patch_vs_nonlesional
#          Df SumOfSqs     R2      F Pr(>F)
# stage     1   0.3799 0.0541 1.2011  0.274
# Residual 21   6.6427 0.9459
# Total    22   7.0226 1.0000
#
# $Patch_vs_Plaque
#          Df SumOfSqs      R2      F Pr(>F)
# stage     1   1.0449 0.17151 3.1052  0.019 *
# Residual 15   5.0475 0.82849
# Total    16   6.0924 1.00000
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# $nonlesional_vs_Plaque
#          Df SumOfSqs      R2     F Pr(>F)
# stage     1   0.7649 0.09439 2.293  0.047 *
# Residual 22   7.3387 0.90561
# Total    23   8.1036 1.00000
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# PERMISP
library(vegan)
betadisper_species <- betadisper(
  d = bray_species,
  group = metadata$stage
) %>%
  anova()
print(betadisper_species)
# Analysis of Variance Table
#
# Response: Distances
#           Df  Sum Sq  Mean Sq F value Pr(>F)
# Groups     2 0.01927 0.009637  0.1805 0.8358
# Residuals 29 1.54822 0.053387


############################################################
### Differential Abundance Analysis of Microbial Species ###
###                      Figure 3a                       ###
############################################################

library(Maaslin2)
diff_abundance_analysis <- Maaslin2(
  input_data = feature_table,
  input_metadata = metadata,
  min_abundance = 0.0,
  min_prevalence = 0.25,
  normalization = "NONE",
  analysis_method = "LM",
  random_effects = "patient",
  fixed_effects = c(
    "stage",
    "Patch_seq_depth",
    "Plaque_seq_depth",
    "nonlesional_seq_depth"
    ),
  reference = c("stage,nonlesional"),
  output = "diff_abundance_analysis"
)

# Transform coef into log2fc
diff_abundance_analysis$results$fc <- exp(
  diff_abundance_analysis$results$coef
)
diff_abundance_analysis$results$log2fc <- log2(
  diff_abundance_analysis$results$fc
)

# Retain only results for the dependent variable "stage"
library(dplyr)
results <- filter(
  diff_abundance_analysis$results,
  grepl("stage", metadata)
  )

########################
### PCA of PanPhlAn  ###
###   Figure 6       ###
########################

# Shown is a representitive example for S. aureus
# The other bacteria were processed the same way


#load presence/absence matrix of gene families into R
# (1 means the gene family is present in the sample, 0 it is not)
library(tidyverse)
Saureus <- read.table(
  "Panphlan_Saureus.tsv",
  sep = "\t",
  strip.white = TRUE,
  stringsAsFactors = FALSE,
  row.names = 1,
  header = TRUE,
  check.names = FALSE
  )

Saureus <- t(Saureus)

# Remove  genefamilies that were entirely absent in all strains
# i.e., columns that are zero
Saureus <- Saureus[, colSums(Saureus != 0) > 0]

# Perform PCA
Saureus_pca <- prcomp(Saureus, center = T, scale. = T)
summary(Saureus_pca)
# Importance of components:
#                          PC1    PC2    PC3    PC4     PC5     PC6     PC7     PC8     PC9    PC10
# Standard deviation     11.5193 9.1731 8.9037 8.1906 6.82921 6.32964 6.11302 5.43338 4.39771 4.30506
# Proportion of Variance  0.2369 0.1502 0.1416 0.1198 0.08328 0.07154 0.06673 0.05271 0.03453 0.03309
# Cumulative Proportion   0.2369 0.3872 0.5287 0.6485 0.73181 0.80335 0.87007 0.92278 0.95732 0.99041
#                          PC11      PC12
# Standard deviation     2.31737 1.017e-13
# Proportion of Variance 0.00959 0.000e+00
# Cumulative Proportion  1.00000 1.000e+00

### Plot the PCA
# Get Variables to colour and shape the PCA for S. aureus
metadataSaureus <- filter(
  metadata,
  rowname %in% rownames(Saureus_pca[["x"]])
  )

# Make sure that the sample order in metadataSaureus
# is the same as in the pcobj object
sample_order <- Saureus_pca[["x"]]
metadataSaureus <- metadataSaureus[rownames(sample_order),]
all(rownames(sample_order) == rownames(metadataSaureus))
# [1] TRUE

# Define Shapes and Colours to plot Patient and Stage Variables
Stage_shape <- c(16, 15, 17)
names(Stage_shape) <- c("nonlesional", "Patch", "Plaque")

Patient_Colours <- c(
  "#F8766D",
  "#EA8331",
  "#D89000",
  "#C09B00",
  "#A3A500",
  "#7CAE00",
  "#39B600",
  "#00BB4E",
  "#00BF7D",
  "#00C1A3",
  "#00BFC4",
  "#00BAE0",
  "#00B0F6",
  "#35A2FF",
  "#9590FF",
  "#C77CFF",
  "#E76BF3",
  "#FA62DB",
  "#FF62BC",
  "#FF6A98"
  )
names(Patient_Colours) <- c(
  "Pat1",
  "Pat2",
  "Pat3",
  "Pat4",
  "Pat5",
  "Pat6",
  "Pat7",
  "Pat8",
  "Pat9",
  "Pat10",
  "Pat11",
  "Pat12",
  "Pat13",
  "Pat14",
  "Pat15",
  "Pat16",
  "Pat17",
  "Pat18",
  "Pat19",
  "Pat20"
  )

# Plot
library(ggplot2)
library(ggpubr)
Saureus.pca.plot <- autoplot(
  Saureus_least_stringent_no0.pca,
  data = metadataSaureus_less_stringent,
  colour = "Patient",
  shape = "Stage", size=4) +
  theme_bw() +
  ggtitle("S. aureus") +
  theme(plot.title = element_text(face = "italic")) +
  scale_shape_manual(values = Stage_shape) +
  scale_color_manual(values = Patient_Colours)

# Repeat the process for S. hominis, S. epidermidis, C. acnes
# Then arrange the plots as follows
ggarrange(
  Saureus.pca.plot,
  Shominis.pca.plot,
  Sepidermidis.pca.plot,
  Cacnes.pca.plot,
  common.legend = TRUE,
  legend = "bottom",
  labels = "auto"
  )

###############################################################################
### Differential Abundance Analysis of Virulence Factors found by ShortBRED ###
###                               Figure 8a                                 ###
###############################################################################

library(Maaslin2)
diff_abundance_vf <- Maaslin2(
  input_data = virulence_factors,
  input_metadata = metadata,
  min_abundance = 0.0,
  min_prevalence = 0.1,
  transform = "NONE",
  normalization = "NONE",
  analysis_method = "LM",
  random_effects = "patient",
  fixed_effects = "stage",
  reference = c("stage,nonlesional"),
  output = "diff_abundance_analysis"
)

# Transform coef into log2fc
diff_abundance_vf$results$fc <- exp(
  diff_abundance_vf$results$coef
)
diff_abundance_vf$results$log2fc <- log2(
  diff_abundance_vf$results$fc
)

##################################
### Heatmap of Control Samples ###
###      Suppl. Figure 7       ###
##################################

library(pheatmap)
library(tidyverse)
library(RColorBrewer)
library(ggplotify)

# Get the color annotation
controls <- read_csv("./zymoBiomics_controls.csv")
controls_col <- controls$color
names(controls_col) <- controls$Sample

ann_colors <- list(
  `Standard Type` = controls_col
)

# Get the data!
controls_species <- read_tsv("./controls_vs_zymoBIOMICS.tsv")
controls_species <- column_to_rownames(controls_species, "sample")
controls_species$sample103_1 <- NULL
controls_species$sample87_1 <- NULL
controls_species$sample70_1 <- NULL
metadata <- read_tsv("./metadata-controls.tsv")
metadata <- column_to_rownames(metadata, "Sample")
metadata$`Standard Type` <- factor(metadata$`Standard Type`, levels = c(
  "ZymoBIOMICS",
  "Cell-derived",
  "DNA-derived"
))

### Do the heatmap
# Inspect how the data is strucured
x <- controls_species %>% as.vector() %>% unlist
ggpubr::ggdensity(x)

q <- unique(
  quantile(x,
  probs = seq(0, 1, 0.001)
))

# Specify the number of color levels
num_colors <- length(q)

# Set the breaks for the color scale
breaks <- c(q, 35)

# Set color palette
option1 <- colorRampPalette(
  brewer.pal(n = 9, name = "YlOrRd"))(num_colors)
option1[1] <- "white"

# Set Itacilized microbial species names
italic_species_names <- lapply(
  rownames(controls_species),
  function(x) bquote(italic(.(x))))

# Plot the heatmap
heatmap_controls <-
pheatmap(
  mat = controls_species,
  color = option1,
  breaks = breaks,
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  cutree_cols = 3,
  treeheight_row = 0,
  treeheight_col = 25,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_col = metadata,
  labels_row = as.expression(italic_species_names),
  annotation_colors = ann_colors,
  display_numbers = FALSE
  )