###### This script is for merging pg-matrices from DIA-NN for the further processing in AlphPeptStats ###########

library(tidyverse)
library(stringr)
library(openxlsx)


# Files loading -----------------------------------------------------------

# read pg_matric files into list
# ensure all pg_matrix files in a separate folder and having a .tsv extension
matrix_files <- list.files(path = 'D:/O2TM_proteomics_processed/DIANN_otp_pg_mtrcs/', pattern="*.tsv", full.names=TRUE) # change a path
pg_mtrx_dfs <- lapply(matrix_files, read_tsv)


# Files processing --------------------------------------------------------

#### This script loop over dfs in a list, assigns 0 to NA, then parse path names to patient codes, sort columns with intensities,
#### assigns colnames as patient codes to the columns with intensities, then separate df to two fractions df and one fraction df, sum
#### two fraction df, cbind it with one fraction df, adding first 5 columns from original pg_matrix

for(idx in c(1:length(pg_mtrx_dfs))){
  
  idx_needed_clmns <- c(6:ncol(pg_mtrx_dfs[[idx]]))
  
  pg_mtrx_dfs[[idx]][, idx_needed_clmns][is.na(pg_mtrx_dfs[[idx]][, idx_needed_clmns])] <- 0
  
  patient_codes <- str_extract(colnames(pg_mtrx_dfs[[idx]][idx_needed_clmns]),
                               '(?<=[:digit:]_)[:alnum:]+(\\*?)(?=(_20%|_100%))')
  
  sorted_pat_codes <- sort(patient_codes)
  ordered_pat_codes <- order(patient_codes)
  
  pg_mtrx_dfs[[idx]][idx_needed_clmns] <- pg_mtrx_dfs[[idx]][idx_needed_clmns][, ordered_pat_codes]
  
  colnames(pg_mtrx_dfs[[idx]])[idx_needed_clmns] <- sorted_pat_codes
  
  duplicated_names <- colnames(pg_mtrx_dfs[[idx]][idx_needed_clmns])[duplicated(colnames(pg_mtrx_dfs[[idx]][idx_needed_clmns])) | duplicated(colnames(pg_mtrx_dfs[[idx]][idx_needed_clmns]),
                                                                                                                 fromLast=TRUE)]
  idx_of_duplicated <- which(colnames(pg_mtrx_dfs[[idx]][idx_needed_clmns]) %in% duplicated_names)
  lfq_intens_col_ordered_two_frac <- pg_mtrx_dfs[[idx]][, idx_needed_clmns][, idx_of_duplicated]
  lfq_intens_col_ordered_one_frac <- pg_mtrx_dfs[[idx]][, idx_needed_clmns][, -idx_of_duplicated]
  
  # sum intensities in files that have two fractions
  first_columns <- seq(1, ncol(lfq_intens_col_ordered_two_frac) - 1, 2)
  second_columns <- seq(2, ncol(lfq_intens_col_ordered_two_frac), 2)
  lfq_intens_sumed <- lfq_intens_col_ordered_two_frac[, first_columns] + lfq_intens_col_ordered_two_frac[, second_columns]
  
  # add samples with one fraction to summed
  pg_mtrx_dfs[[idx]] <- cbind(pg_mtrx_dfs[[idx]][, -idx_needed_clmns], lfq_intens_col_ordered_one_frac, lfq_intens_sumed)
  
}


# Full join of dfs in a list ----------------------------------------------

#### We need a full join to add all proteins accross files and all patient codes accros the files

new_df <- pg_mtrx_dfs %>%
  reduce(full_join, by = c('Protein.Group', 'Protein.Ids', 'Protein.Names', 'Genes', 'First.Protein.Description'))

# fill in NA after full join
new_df[is.na(new_df)] <- 0


# Write data in a file ----------------------------------------------------


write.xlsx(new_df, 'D:/O2TM_proteomics_processed/pg_matrices_o2td_diann.xlsx')
