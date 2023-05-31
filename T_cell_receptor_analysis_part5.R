                                            ################################
                                            ### T cell receptor analysis ###
                                            ################################

                                                      ##############
                                                      ### Part 5 ###
                                                      ##############

### Extract the meaningful results (language: R)
# Background: We now profiled a nigh number of TCR-epitope pairs (n = 3950)
# for every sample. Most of it is noise.
# To dynamically extract the TCRs that are reactive to one of
# the S. aureus epitopes per sample, we calculate the point at
# which the binding score drops the most from one TCR-epitope-pair
# to the next. Only the most reactive TCR-epitope pairs will be retained.

# load data
library(tidyverse)
ERGOII_results <- list.files(path = "path/to/ERGOII/results/Saureus_IEDB_predicted/", pattern = "\\.csv$")
score_list_Saureus <- list()

for (ERGOII in ERGOII_results) {

  sample_name <- str_split(ERGOII, pattern = "\\.")
  sample_name <- sample_name[[1]][[1]]
  sample_name <- str_split(sample_name, "_")
  sample_name <- sample_name[[1]][[1]]

  filepath_ERGOII <- paste0("path/to/ERGOII/results/Saureus_IEDB_predicted/", ERGOII)
  df_ERGOII <- read_csv(filepath_ERGOII, col_names = T, col_select = c(2,5,6,8,10)) %>% suppressMessages()
  df_ERGOII <- df_ERGOII[order(df_ERGOII$Score, decreasing = T),]
  df_ERGOII$row_num <- as.numeric(rownames(df_ERGOII))
  score_list_Saureus[[sample_name]] <- df_ERGOII
}

# Function to find the Peptide where
# the slope decreases the most and subset data.frame there
find_decline_diff <- function(df) {
  # Calculate first and second derivative
  diffs <- diff(df$Score)
  slopes <- diffs / diff(seq_along(df$row_num))
  diff_slopes <- diff(slopes) / diff(seq_along(df$row_num[-1]))

  # Find where second derivative is smallest
  max_decline_index <- which.min(diff_slopes)
  cat("The point where the slope decreases the most is at peptide", df$row_num[max_decline_index + 1], "\n")

  # Get row_index
  row_index = max_decline_index+1

  # Subset data.frame
  df_subset = df[1:row_index,]
  return(df_subset)
}

# Apply the function to my ERGOII list (score_list)
df_ERGOII_MaxSlope <- lapply(df_ERGOII, find_decline_diff)

# Add Epitope information
for (sample in names(df_ERGOII_MaxSlope)){
  df_ERGOII_MaxSlope[[sample]] <-
  left_join(df_ERGOII_MaxSlope[[sample]], Saureus_epitopes_IEDB_predicted, by = "Peptide")
}