                                            ################################
                                            ### T cell receptor analysis ###
                                            ################################

                                                      ##############
                                                      ### Part 3 ###
                                                      ##############

# Load the data
# language: R
TCRs <- list.files(path = "/path/to/TCR_clean/", pattern = "\\.csv$")
TCR_list = list()

for (TCR in TCRs){

  sample_name = str_split(TCR, pattern = "\\.")
  sample_name = sample_name[[1]][[1]]

  filepath_TCR = paste0("/path/to/TCR_clean/", TCR)
  df_TCR = read_tsv(
    filepath_TCR, col_names = T, col_select = c("cdr3aa", "v", "j")
    ) %>% suppressMessages()
  df_TCR = rename(df_TCR, TRB = cdr3aa, TRBV = v, TRBJ = j)

  TCR_list[[sample_name]] = df_TCR
}

# Load epitopes and Check for duplicates between predicted and emperic epitopes
# language: R
Saureus_epitopes_IEDB_predicted = read_csv(
    "path/to/Saureus_epitopes_IEDB_and_predicted.csv")

table(duplicated(Saureus_epitopes_IEDB_predicted$Peptide)) # 1 duplicate found!
Saureus_epitopes_IEDB_predicted <-
Saureus_epitopes_IEDB_predicted[
    !duplicated(Saureus_epitopes_IEDB_predicted$Peptide), ]
table(duplicated(Saureus_epitopes_IEDB_predicted$Peptide)) # 0 duplicates found!

# Build the TCR-epitope-pairs
# language: R
TCR_epitopeSaureus_IEDB_predicted_pairs = list()

for (TCR in names(TCR_list)){

  df <- expand_grid(
    TCR_list[[TCR]][1:50,], Saureus_epitopes_IEDB_predicted$Peptide)
  # Need to bring the df into the format that ERGO-II requires
  library(tidyverse)
  df <- rename(df, Peptide = `Saureus_epitopes_IEDB_predicted$Peptide`)
  df$TRA <- NA
  df$TRAV <- NA
  df$TRAJ <- NA
  df$`T-Cell-Type` <- NA
  df$MHC <- NA
  df <- relocate(
    df, c(
        "TRA",
        "TRB",
        "TRAV",
        "TRAJ",
        "TRBV",
        "TRBJ",
        "T-Cell-Type",
        "Peptide",
        "MHC"
        ))

  TCR_epitopeSaureus_IEDB_predicted_pairs[[TCR]] = df
}

#export to csv
for (TCR in names(TCR_epitopeSaureus_IEDB_predicted_pairs)) {
  filename = paste0(
    "path/to/ERGOII/epitope_pairs/Saureus_IEDB_predicted/", TCR, ".csv")
  write_csv(
    TCR_epitopeTcellHTLV_neoplasms_pairs[[TCR]], file = filename, na = "")
}

###########################
### PROCEED WITH PART 4 ###
###########################