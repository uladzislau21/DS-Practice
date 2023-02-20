library(devtools)
library(diann)
library(tidyverse)
library(stringr)
library(openxlsx)

###### I manually made some modifications to DIA-NN output files as well as to the names of RAW files (PC and server) #######
###### I added patient codes to the 2nd group (File.Name and Run columns of DIA-NN output) and to RAW file names      #######
###### I added * to the duplicated patient codes (AK123, AK112), 21st, 9th, 1st, 5th groups. Again both to file names #######
###### and to the columns in DIA-NN output)                                                                           #######


# Data loading ------------------------------------------------------------


# read samples description
sample_info <- read.xlsx('D:/O2TM_proteomics_processed/SERUM_GROUPS.xlsx') # change a path
sample_info$Group <- paste(sample_info$Group, sample_info$Number, sep = '')
sample_info <- sample_info[, -c(3)]

# read files final output from DIA-NN into list
# ensure all output files in a separate folder and having a .tsv extension
diann_files <- list.files(path = 'D:/O2TM_proteomics_processed/DIANN_output/', pattern="*.tsv", full.names=TRUE) # change a path
diann_dfs <- lapply(diann_files, read_tsv)

# Data processing ---------------------------------------------------------


# filter protein groups identified with minimum of 2 peptides
for(df_idx in c(1:length(diann_dfs))){
  
  # identify which prot groups were identified with min 2 peptides
  diann_2pept <- diann_dfs[[df_idx]] %>%
    group_by(Protein.Group) %>%
    summarise(count = n_distinct(Stripped.Sequence))
  
  # get protein groups names identified with min 2 pept
  diann_2pept <- diann_2pept[diann_2pept$count > 1,]
  protein_groups <- diann_2pept$Protein.Group
  
  # filter out identified with less than 2 peptides
  diann_dfs[[df_idx]] <- diann_dfs[[df_idx]][diann_dfs[[df_idx]]$Protein.Group %in% protein_groups, ]
}

# merge dfs in a list
diann_df <- bind_rows(diann_dfs)

## apply LFQ algorithm from diann R package (!!! no 0.01 FDR filtrattion !!!)
lfq_intens <- diann_maxlfq(diann_df, group.header="Protein.Group",
                                  id.header = "Precursor.Id",
                                  quantity.header = "Precursor.Normalised")


# process rownames, make them with one identifier
rownames(lfq_intens) <- str_extract(rownames(lfq_intens), '\\w+|[:punctuation:](?=;)')

# remove NAs, change for 0
lfq_intens[is.na(lfq_intens)] <- 0

## sum 20 and 100 fractions
patient_codes <- str_extract(colnames(lfq_intens),
                             '(?<=[:digit:]_)[:alnum:]+(\\*?)(?=(_20%|_100%))')

sorted_pat_codes <- sort(patient_codes)
ordered_pat_codes <- order(patient_codes)

# sort columns according to sorted patient codes in order to sum columns pair by pair
lfq_intens_col_ordered <- lfq_intens[,ordered_pat_codes]

# rename columns to patient codes
colnames(lfq_intens_col_ordered) <- sorted_pat_codes

# filter out not duplicated columns (sometimes 0 proteins are found in a run so DIA-NN does not put these files
# names in the final output)
duplicated_names <- colnames(lfq_intens_col_ordered)[duplicated(colnames(lfq_intens_col_ordered)) | duplicated(colnames(lfq_intens_col_ordered),
                                                                                           fromLast=TRUE)]
idx_of_duplicated <- which(colnames(lfq_intens_col_ordered) %in% duplicated_names)
lfq_intens_col_ordered_two_frac <- lfq_intens_col_ordered[, idx_of_duplicated]
lfq_intens_col_ordered_one_frac <- lfq_intens_col_ordered[, -idx_of_duplicated]

# sum intensities in files that have two fractions
first_columns <- seq(1, ncol(lfq_intens_col_ordered_two_frac) - 1, 2)
second_columns <- seq(2, ncol(lfq_intens_col_ordered_two_frac), 2)
lfq_intens_sumed <- lfq_intens_col_ordered_two_frac[, first_columns] + lfq_intens_col_ordered_two_frac[, second_columns]

# add samples with one fraction to summed
lfq_intens_sumed <- cbind(lfq_intens_col_ordered_one_frac, lfq_intens_sumed)

## rename columns to D1, C1, O1, D2,......
group_ids <- match(colnames(lfq_intens_sumed), sample_info$ID)
colnames(lfq_intens_sumed) <- sample_info[group_ids, 'Group']
lfq_intens_sumed <- lfq_intens_sumed[, sort(colnames(lfq_intens_sumed))]


# Data processing to export for Miklos ------------------------------------

##### This section is to assign TAJ card ids instead of patient codes.
##### This will be used to introduce MS measurments to database
##### Use lfq_intens_summed df before renaming columns to D1, C1, O1, D2,......

# make patient codes as row names
lfq_intens_sumed_T <- t(lfq_intens_sumed)

# Data export -------------------------------------------------------------


## write final table
write.xlsx(as.data.frame(lfq_intens_sumed), 'D:/O2TM_proteomics_processed/final_table_o2td_diann.xlsx', rowNames = TRUE) # change a path
