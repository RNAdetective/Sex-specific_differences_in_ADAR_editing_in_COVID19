# Overall ADAR editing levels

library(dplyr)

#=============== Pre-females ===================
input_folder <- "G:\\Chapter4\\ch4_ctlfemale.89_subset\\ch4_ctlfemale.89_subset\\ch4_ctlfemale_89.snpEff_counts_csv\\snpEff_counts_csv"
output_folder <- "G:\\Chapter4\\Overall_ADAR_editing\\Pre_infection_females"
#==============================================


files <- list.files(input_folder, full.names = TRUE, pattern = "\\.csv$|\\.tsv$")

for (file in files) {
  
  is_csv <- grepl("\\.csv$", file)
  
  df <- if (is_csv) {
    read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    read.table(file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  }
  
  # remove columns if present
  df <- dplyr::select(df, -dplyr::any_of(c("misc1_Rank", "misc2_HGVS.c", "misc3_HGVS.p")))
  
  # create new col
  df$ref_alt <- paste0(df$REF, df$ALT)
  
  # filter AG and TC
  df <- dplyr::filter(df, ref_alt %in% c("AG", "TC"))
  
  # write output
  out_file <- file.path(output_folder, basename(file))
  
  if (is_csv) {
    write.csv(df, out_file, row.names = FALSE)
  } else {
    write.table(df, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

#=============== Pre-males ===================
input_folder <- "G:\\Chapter4\\ch4_ctlmale.100_subset\\ch4_ctlmale.100_subset\\ch4_ctlmale_100_215.snpEff_counts_csv\\snpEff_counts_csv"
output_folder <- "G:\\Chapter4\\Overall_ADAR_editing\\Pre_infection_males"
#==============================================


files <- list.files(input_folder, full.names = TRUE, pattern = "\\.csv$|\\.tsv$")

for (file in files) {
  
  is_csv <- grepl("\\.csv$", file)
  
  df <- if (is_csv) {
    read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    read.table(file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  }
  
  # remove columns 
  df <- dplyr::select(df, -dplyr::any_of(c("misc1_Rank", "misc2_HGVS.c", "misc3_HGVS.p")))
  
  # create new col
  df$ref_alt <- paste0(df$REF, df$ALT)
  
  # filter AG and TC
  df <- dplyr::filter(df, ref_alt %in% c("AG", "TC"))
  
  # write
  out_file <- file.path(output_folder, basename(file))
  
  if (is_csv) {
    write.csv(df, out_file, row.names = FALSE)
  } else {
    write.table(df, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

#=============== Mid-females ===================
input_folder <- "G:\\Chapter4\\ch4_midfemale.89_subset\\ch4_midfemale.89_subset\\ch4_midfemale_89.snpEff_counts_csv\\snpEff_counts_csv"
output_folder <- "G:\\Chapter4\\Overall_ADAR_editing\\Mid_infection_females"
#==============================================


files <- list.files(input_folder, full.names = TRUE, pattern = "\\.csv$|\\.tsv$")

for (file in files) {
  
  is_csv <- grepl("\\.csv$", file)
  
  df <- if (is_csv) {
    read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    read.table(file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  }
  
  # remove columns
  df <- dplyr::select(df, -dplyr::any_of(c("misc1_Rank", "misc2_HGVS.c", "misc3_HGVS.p")))
  
  # create new col
  df$ref_alt <- paste0(df$REF, df$ALT)
  
  # filter AG and TC
  df <- dplyr::filter(df, ref_alt %in% c("AG", "TC"))
  
  # write 
  out_file <- file.path(output_folder, basename(file))
  
  if (is_csv) {
    write.csv(df, out_file, row.names = FALSE)
  } else {
    write.table(df, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
  }
}


#=============== Mid-males ===================
input_folder <- "G:\\Chapter4\\ch4_midmale.100_subset\\ch4_midmale.100_subset\\ch4_midmale_100_345.snpEff_counts_csv\\snpEff_counts_csv"
output_folder <- "G:\\Chapter4\\Overall_ADAR_editing\\Mid_infection_males"
#==============================================


files <- list.files(input_folder, full.names = TRUE, pattern = "\\.csv$|\\.tsv$")

for (file in files) {
  
  is_csv <- grepl("\\.csv$", file)
  
  df <- if (is_csv) {
    read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    read.table(file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  }
  
  # remove columns if present
  df <- dplyr::select(df, -dplyr::any_of(c("misc1_Rank", "misc2_HGVS.c", "misc3_HGVS.p")))
  
  # create new col
  df$ref_alt <- paste0(df$REF, df$ALT)
  
  # filter AG and TC
  df <- dplyr::filter(df, ref_alt %in% c("AG", "TC"))
  
  # write output
  out_file <- file.path(output_folder, basename(file))
  
  if (is_csv) {
    write.csv(df, out_file, row.names = FALSE)
  } else {
    write.table(df, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
  }
}



#=============== Post-females ===================
input_folder <- "G:\\Chapter4\\ch4_postfemale.56_subset\\ch4_postfemale.56_subset\\ch4_postfemale_56.snpEff_counts_csv\\snpEff_counts_csv"
output_folder <- "G:\\Chapter4\\Overall_ADAR_editing\\Post_infection_females"
#==============================================


files <- list.files(input_folder, full.names = TRUE, pattern = "\\.csv$|\\.tsv$")

for (file in files) {
  
  is_csv <- grepl("\\.csv$", file)
  
  df <- if (is_csv) {
    read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    read.table(file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  }
  
  # remove columns
  df <- dplyr::select(df, -dplyr::any_of(c("misc1_Rank", "misc2_HGVS.c", "misc3_HGVS.p")))
  
  # create new col
  df$ref_alt <- paste0(df$REF, df$ALT)
  
  # filter AG and TC
  df <- dplyr::filter(df, ref_alt %in% c("AG", "TC"))
  
  # write output
  out_file <- file.path(output_folder, basename(file))
  
  if (is_csv) {
    write.csv(df, out_file, row.names = FALSE)
  } else {
    write.table(df, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

#=============== Post-males ===================
input_folder <- "G:\\Chapter4\\ch4_postmale.100_subset\\ch4_postmale.100_subset\\ch4_postmale_100_295.snpEff_counts_csv\\snpEff_counts_csv"
output_folder <- "G:\\Chapter4\\Overall_ADAR_editing\\Post_infection_males"
#==============================================


files <- list.files(input_folder, full.names = TRUE, pattern = "\\.csv$|\\.tsv$")

for (file in files) {
  
  is_csv <- grepl("\\.csv$", file)
  
  df <- if (is_csv) {
    read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    read.table(file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  }
  
  # remove columns if present
  df <- dplyr::select(df, -dplyr::any_of(c("misc1_Rank", "misc2_HGVS.c", "misc3_HGVS.p")))
  
  # create new col
  df$ref_alt <- paste0(df$REF, df$ALT)
  
  # filter AG and TC
  df <- dplyr::filter(df, ref_alt %in% c("AG", "TC"))
  
  # write 
  out_file <- file.path(output_folder, basename(file))
  
  if (is_csv) {
    write.csv(df, out_file, row.names = FALSE)
  } else {
    write.table(df, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

#===============================================================================
#                           Overall ADAR editing levels - pre_females
#=================================================================================


## Overall ADAR editing levels

library(dplyr)
library(tidyr)
library(data.table)
library(gridExtra)
library(patchwork)
setwd("G:\\Chapter4\\Overall_ADAR_editing\\")
directoryPath <- "G:\\Chapter4\\Overall_ADAR_editing\\"

#Read Redi file
redi_data <- fread("G:\\Chapter4\\Overall_ADAR_editing\\REDI\\TABLE1_hg38.txt", header = TRUE, fill = TRUE, stringsAsFactors = FALSE)

# Converting position column in redi_data to character (to match with coordinate in sample_file)
redi_data$Position <- as.character(redi_data$Position)
redi_data$Region <- as.character(redi_data$Region)

# Remove "chr" from the value in Region column of redi_data (to match with chromosome in sample_file)
redi_data$Region <- gsub("chr", "", redi_data$Region)
head(redi_data)

# Filter for reads >=100
apply_filter <- function(df) {
  df[df$Total >= 100, ]
}

# Function
calculate_values <- function(input_csv_file) {
  
  sample_file <- fread(input_csv_file)
  setDT(sample_file)
  
  sample_file[, POS := as.character(POS)]
  sample_file[, CHROM := as.character(CHROM)]
  
  
  sample_POS <- unique(sample_file$POS)
  redi_POS <- unique(redi_data$Position)
  
    
  # Merge CHROM and POS
  filtered_sample_file <- sample_file[
    .(redi_data$Region, redi_data$Position),
    on = c("CHROM", "POS"),
    nomatch = 0
  ]
  
  
  # Coverage filter
  filtered_sample_file <- apply_filter(filtered_sample_file)
  
  # Convert count columns to numeric
  filtered_sample_file <- filtered_sample_file %>%
    mutate(across(c(G, C, Total), as.numeric))
  
  # Add editing level
  filtered_sample_file <- filtered_sample_file %>%
    mutate(editinglevel = ifelse(REF == "A", G / Total, C / Total))
  
  # Save 
  filteredFileName = paste0(directoryPath,"Result_Pre_females\\",basename(input_csv_file),"_filtered.csv")
  write.csv(filtered_sample_file, filteredFileName, row.names = FALSE) 

  
  # Summary values
  A_count <- sum(filtered_sample_file[REF == "A", G])
  T_count <- sum(filtered_sample_file[REF == "T", C])
  coverage_sum <- sum(filtered_sample_file$Total, na.rm = TRUE)
  
  overall_editing_level <- (A_count + T_count) / coverage_sum
  
  new_df <- data.frame(
    sample = basename(input_csv_file),
    As = A_count,
    Ts = T_count,
    coverage = coverage_sum,
    overall_editing_level = overall_editing_level
  )
  
    return(new_df)
}

# path to CSV files
csv_files <- list.files(
  "G:/Chapter4/Overall_ADAR_editing/Pre_infection_females",
  pattern = "\\.csv$",
  full.names = TRUE
)

# Run over all files
list_of_new_dfs <- lapply(csv_files, calculate_values)
new_df <- do.call(rbind, list_of_new_dfs)
infection_stage_sex <- rep("Pre_infection_females", nrow(new_df))
new_df$infection_stage_sex <- infection_stage_sex
# Save 
finalResult <- file.path(directoryPath, "Result_Pre_females", "overall_editing_pre_infection_females.csv")
write.csv(new_df, finalResult, row.names = FALSE)

#===============================================================================
#                           Overall ADAR editing levels - pre_males
#=================================================================================


## Overall ADAR editing levels
library(dplyr)
library(tidyr)
library(data.table)
library(gridExtra)
library(patchwork)
setwd("G:\\Chapter4\\Overall_ADAR_editing\\")
directoryPath <- "G:\\Chapter4\\Overall_ADAR_editing\\"

#Read Redi file
redi_data <- fread("G:\\Chapter4\\Overall_ADAR_editing\\REDI\\TABLE1_hg38.txt", header = TRUE, fill = TRUE, stringsAsFactors = FALSE)

# Converting position column in redi_data to character (to match with coordinate in sample_file)
redi_data$Position <- as.character(redi_data$Position)
redi_data$Region <- as.character(redi_data$Region)

# Remove "chr" from the value in Region column of redi_data (to match with chromosome in sample_file)
redi_data$Region <- gsub("chr", "", redi_data$Region)
head(redi_data)

# Filter for reads >=100
apply_filter <- function(df) {
  df[df$Total >= 100, ]
}

# Function
calculate_values <- function(input_csv_file) {
  
    sample_file <- fread(input_csv_file)
  setDT(sample_file)
  
  sample_file[, POS := as.character(POS)]
  sample_file[, CHROM := as.character(CHROM)]
  
 
  sample_POS <- unique(sample_file$POS)
  redi_POS <- unique(redi_data$Position)
  
    
  # Merge using CHROM and POS
  filtered_sample_file <- sample_file[
    .(redi_data$Region, redi_data$Position),
    on = c("CHROM", "POS"),
    nomatch = 0
  ]
  
  
  # Apply coverage filter
  filtered_sample_file <- apply_filter(filtered_sample_file)
  
  # Convert count columns to numeric
  filtered_sample_file <- filtered_sample_file %>%
    mutate(across(c(G, C, Total), as.numeric))
  
  # Add editing level
  filtered_sample_file <- filtered_sample_file %>%
    mutate(editinglevel = ifelse(REF == "A", G / Total, C / Total))
  
  # Save filtered file
  filteredFileName = paste0(directoryPath,"Result_pre_males\\",basename(input_csv_file),"_filtered.csv")
  write.csv(filtered_sample_file, filteredFileName, row.names = FALSE)
  
   
  # Calculate summary values
  A_count <- sum(filtered_sample_file[REF == "A", G])
  T_count <- sum(filtered_sample_file[REF == "T", C])
  coverage_sum <- sum(filtered_sample_file$Total, na.rm = TRUE)
  
  overall_editing_level <- (A_count + T_count) / coverage_sum
  
  new_df <- data.frame(
    sample = basename(input_csv_file),
    As = A_count,
    Ts = T_count,
    coverage = coverage_sum,
    overall_editing_level = overall_editing_level
  )
  
   return(new_df)
}

# Path containing filtered CSV files
csv_files <- list.files(
  "G:/Chapter4/Overall_ADAR_editing/Pre_infection_males",
  pattern = "\\.csv$",
  full.names = TRUE
)

# Run over all files
list_of_new_dfs <- lapply(csv_files, calculate_values)
new_df <- do.call(rbind, list_of_new_dfs)
infection_stage_sex <- rep("Pre_infection_males", nrow(new_df))
new_df$infection_stage_sex <- infection_stage_sex
# Save 
finalResult <- file.path(directoryPath, "Result_pre_males", "overall_editing_pre_infection_males.csv")
write.csv(new_df, finalResult, row.names = FALSE)


#===============================================================================
#                           Overall ADAR editing levels - Mid_females
#=================================================================================


## Overall ADAR editing levels
library(dplyr)
library(tidyr)
library(data.table)
library(gridExtra)
library(patchwork)
setwd("G:\\Chapter4\\Overall_ADAR_editing\\")
directoryPath <- "G:\\Chapter4\\Overall_ADAR_editing\\"

# Read Redi file
redi_data <- fread("G:\\Chapter4\\Overall_ADAR_editing\\REDI\\TABLE1_hg38.txt", header = TRUE, fill = TRUE, stringsAsFactors = FALSE)

# Converting position column in redi_data to character (to match with coordinate in sample_file)
redi_data$Position <- as.character(redi_data$Position)
redi_data$Region <- as.character(redi_data$Region)

# Remove "chr" from the value in Region column of redi_data (to match with chromosome in sample_file)
redi_data$Region <- gsub("chr", "", redi_data$Region)
head(redi_data)

# Filter for reads >=100
apply_filter <- function(df) {
  df[df$Total >= 100, ]
}

# Function
calculate_values <- function(input_csv_file) {
  
   sample_file <- fread(input_csv_file)
  setDT(sample_file)
  
  sample_file[, POS := as.character(POS)]
  sample_file[, CHROM := as.character(CHROM)]
  
  
  sample_POS <- unique(sample_file$POS)
  redi_POS <- unique(redi_data$Position)
  
  
  
  # Merge using CHROM and POS
  filtered_sample_file <- sample_file[
    .(redi_data$Region, redi_data$Position),
    on = c("CHROM", "POS"),
    nomatch = 0
  ]
  
  
  # Apply coverage filter
  filtered_sample_file <- apply_filter(filtered_sample_file)
  
  # Convert count columns to numeric
  filtered_sample_file <- filtered_sample_file %>%
    mutate(across(c(G, C, Total), as.numeric))
  
  # Add editing level
  filtered_sample_file <- filtered_sample_file %>%
    mutate(editinglevel = ifelse(REF == "A", G / Total, C / Total))
  
  # Save filtered file
  filteredFileName = paste0(directoryPath,"Result_Mid_females\\",basename(input_csv_file),"_filtered.csv")
  write.csv(filtered_sample_file, filteredFileName, row.names = FALSE)
  
   
  # Calculate summary values
  A_count <- sum(filtered_sample_file[REF == "A", G])
  T_count <- sum(filtered_sample_file[REF == "T", C])
  coverage_sum <- sum(filtered_sample_file$Total, na.rm = TRUE)
  
  overall_editing_level <- (A_count + T_count) / coverage_sum
  
  new_df <- data.frame(
    sample = basename(input_csv_file),
    As = A_count,
    Ts = T_count,
    coverage = coverage_sum,
    overall_editing_level = overall_editing_level
  )
  
   return(new_df)
}

# Path containing your filtered CSV files
csv_files <- list.files(
  "G:/Chapter4/Overall_ADAR_editing/Mid_infection_females",
  pattern = "\\.csv$",
  full.names = TRUE
)

# Run over all files
list_of_new_dfs <- lapply(csv_files, calculate_values)
new_df <- do.call(rbind, list_of_new_dfs)
infection_stage_sex <- rep("Mid_infection_females", nrow(new_df))
new_df$infection_stage_sex <- infection_stage_sex
# Save 
finalResult <- file.path(directoryPath, "Result_Mid_females", "overall_editing_Mid_infection_females.csv")
write.csv(new_df, finalResult, row.names = FALSE)


#===============================================================================
#                           Overall ADAR editing levels - Mid_males
#=================================================================================


## Overall ADAR editing levels
library(dplyr)
library(tidyr)
library(data.table)
library(gridExtra)
library(patchwork)
setwd("G:\\Chapter4\\Overall_ADAR_editing\\")
directoryPath <- "G:\\Chapter4\\Overall_ADAR_editing\\"

#Read Redi file
redi_data <- fread("G:\\Chapter4\\Overall_ADAR_editing\\REDI\\TABLE1_hg38.txt", header = TRUE, fill = TRUE, stringsAsFactors = FALSE)

# Converting position column in redi_data to character (to match with coordinate in sample_file)
redi_data$Position <- as.character(redi_data$Position)
redi_data$Region <- as.character(redi_data$Region)

# Remove "chr" from the value in Region column of redi_data (to match with chromosome in sample_file)
redi_data$Region <- gsub("chr", "", redi_data$Region)
head(redi_data)

# Filter for reads >=100
apply_filter <- function(df) {
  df[df$Total >= 100, ]
}

# Main function
calculate_values <- function(input_csv_file) {
  
  
  sample_file <- fread(input_csv_file)
  setDT(sample_file)
  
  sample_file[, POS := as.character(POS)]
  sample_file[, CHROM := as.character(CHROM)]
  
 
  sample_POS <- unique(sample_file$POS)
  redi_POS <- unique(redi_data$Position)
  
    
  # Merge using CHROM and POS
  filtered_sample_file <- sample_file[
    .(redi_data$Region, redi_data$Position),
    on = c("CHROM", "POS"),
    nomatch = 0
  ]
  
  
  # Apply coverage filter
  filtered_sample_file <- apply_filter(filtered_sample_file)
  
  # Convert count columns to numeric
  filtered_sample_file <- filtered_sample_file %>%
    mutate(across(c(G, C, Total), as.numeric))
  
  # Add editing level
  filtered_sample_file <- filtered_sample_file %>%
    mutate(editinglevel = ifelse(REF == "A", G / Total, C / Total))
  
  # Save 
  filteredFileName = paste0(directoryPath,"Result_Mid_males\\",basename(input_csv_file),"_filtered.csv")
  write.csv(filtered_sample_file, filteredFileName, row.names = FALSE)
  
  
  # Calculate summary values
  A_count <- sum(filtered_sample_file[REF == "A", G])
  T_count <- sum(filtered_sample_file[REF == "T", C])
  coverage_sum <- sum(filtered_sample_file$Total, na.rm = TRUE)
  
  overall_editing_level <- (A_count + T_count) / coverage_sum
  
  new_df <- data.frame(
    sample = basename(input_csv_file),
    As = A_count,
    Ts = T_count,
    coverage = coverage_sum,
    overall_editing_level = overall_editing_level
  )
  
   return(new_df)
}

# Path 
csv_files <- list.files(
  "G:/Chapter4/Overall_ADAR_editing/Mid_infection_males",
  pattern = "\\.csv$",
  full.names = TRUE
)

# Run over all files
list_of_new_dfs <- lapply(csv_files, calculate_values)
new_df <- do.call(rbind, list_of_new_dfs)
infection_stage_sex <- rep("Mid_infection_males", nrow(new_df))
new_df$infection_stage_sex <- infection_stage_sex
# Save 
finalResult <- file.path(directoryPath, "Result_Mid_males", "overall_editing_Mid_infection_males.csv")
write.csv(new_df, finalResult, row.names = FALSE)


#===============================================================================
#                           Overall ADAR editing levels - Post_females
#=================================================================================


## Overall ADAR editing levels
library(dplyr)
library(tidyr)
library(data.table)
library(gridExtra)
library(patchwork)
setwd("G:\\Chapter4\\Overall_ADAR_editing\\")
directoryPath <- "G:\\Chapter4\\Overall_ADAR_editing\\"

#Read Redi file
redi_data <- fread("G:\\Chapter4\\Overall_ADAR_editing\\REDI\\TABLE1_hg38.txt", header = TRUE, fill = TRUE, stringsAsFactors = FALSE)

# Converting position column in redi_data to character (to match with coordinate in sample_file)
redi_data$Position <- as.character(redi_data$Position)
redi_data$Region <- as.character(redi_data$Region)

# Remove "chr" from the value in Region column of redi_data (to match with chromosome in sample_file)
redi_data$Region <- gsub("chr", "", redi_data$Region)
head(redi_data)

# Filter for reads >=100
apply_filter <- function(df) {
  df[df$Total >= 100, ]
}

# Function
calculate_values <- function(input_csv_file) {
  
   sample_file <- fread(input_csv_file)
  setDT(sample_file)
  
  sample_file[, POS := as.character(POS)]
  sample_file[, CHROM := as.character(CHROM)]
  
   
  sample_POS <- unique(sample_file$POS)
  redi_POS <- unique(redi_data$Position)
  
    
  # Merge using CHROM and POS
  filtered_sample_file <- sample_file[
    .(redi_data$Region, redi_data$Position),
    on = c("CHROM", "POS"),
    nomatch = 0
  ]
  
   
  # Apply coverage filter
  filtered_sample_file <- apply_filter(filtered_sample_file)
  
  # Convert count columns to numeric
  filtered_sample_file <- filtered_sample_file %>%
    mutate(across(c(G, C, Total), as.numeric))
  
  # Add editing level
  filtered_sample_file <- filtered_sample_file %>%
    mutate(editinglevel = ifelse(REF == "A", G / Total, C / Total))
  
  # Save filtered file
  filteredFileName = paste0(directoryPath,"Result_Post_females\\",basename(input_csv_file),"_filtered.csv")
  write.csv(filtered_sample_file, filteredFileName, row.names = FALSE)
  
   
  # Calculate summary values
  A_count <- sum(filtered_sample_file[REF == "A", G])
  T_count <- sum(filtered_sample_file[REF == "T", C])
  coverage_sum <- sum(filtered_sample_file$Total, na.rm = TRUE)
  
  overall_editing_level <- (A_count + T_count) / coverage_sum
  
  new_df <- data.frame(
    sample = basename(input_csv_file),
    As = A_count,
    Ts = T_count,
    coverage = coverage_sum,
    overall_editing_level = overall_editing_level
  )
  
    return(new_df)
}

# Path to filterd CSV
csv_files <- list.files(
  "G:/Chapter4/Overall_ADAR_editing/Post_infection_females",
  pattern = "\\.csv$",
  full.names = TRUE
)

# Run over all files
list_of_new_dfs <- lapply(csv_files, calculate_values)
new_df <- do.call(rbind, list_of_new_dfs)
infection_stage_sex <- rep("Post_infection_females", nrow(new_df))
new_df$infection_stage_sex <- infection_stage_sex
# Save 
finalResult <- file.path(directoryPath, "Result_Post_females", "overall_editing_Post_infection_females.csv")
write.csv(new_df, finalResult, row.names = FALSE)


#===============================================================================
#                           Overall ADAR editing levels - Post_males
#=================================================================================


## Overall ADAR editing levels
library(dplyr)
library(tidyr)
library(data.table)
library(gridExtra)
library(patchwork)
setwd("G:\\Chapter4\\Overall_ADAR_editing\\")
directoryPath <- "G:\\Chapter4\\Overall_ADAR_editing\\"

#Read Redi file
redi_data <- fread("G:\\Chapter4\\Overall_ADAR_editing\\REDI\\TABLE1_hg38.txt", header = TRUE, fill = TRUE, stringsAsFactors = FALSE)

# Converting position column in redi_data to character (to match with coordinate in sample_file)
redi_data$Position <- as.character(redi_data$Position)
redi_data$Region <- as.character(redi_data$Region)

# Remove "chr" from the value in Region column of redi_data (to match with chromosome in sample_file)
redi_data$Region <- gsub("chr", "", redi_data$Region)
head(redi_data)

# Filter for reads >=100
apply_filter <- function(df) {
  df[df$Total >= 100, ]
}

# Function
calculate_values <- function(input_csv_file) {
  
 
  sample_file <- fread(input_csv_file)
  setDT(sample_file)
  
  sample_file[, POS := as.character(POS)]
  sample_file[, CHROM := as.character(CHROM)]
  
    sample_POS <- unique(sample_file$POS)
  redi_POS <- unique(redi_data$Position)
  
  
  
  
  # Merge using CHROM and POS
  filtered_sample_file <- sample_file[
    .(redi_data$Region, redi_data$Position),
    on = c("CHROM", "POS"),
    nomatch = 0
  ]
  
   
  # Apply coverage filter
  filtered_sample_file <- apply_filter(filtered_sample_file)
  
  # Convert count columns to numeric
  filtered_sample_file <- filtered_sample_file %>%
    mutate(across(c(G, C, Total), as.numeric))
  
  # Add editing level
  filtered_sample_file <- filtered_sample_file %>%
    mutate(editinglevel = ifelse(REF == "A", G / Total, C / Total))
  
  # Save filtered file
  filteredFileName = paste0(directoryPath,"Result_Post_males\\",basename(input_csv_file),"_filtered.csv")
  write.csv(filtered_sample_file, filteredFileName, row.names = FALSE)
  
    
  # Calculate summary values
  A_count <- sum(filtered_sample_file[REF == "A", G])
  T_count <- sum(filtered_sample_file[REF == "T", C])
  coverage_sum <- sum(filtered_sample_file$Total, na.rm = TRUE)
  
  overall_editing_level <- (A_count + T_count) / coverage_sum
  
  new_df <- data.frame(
    sample = basename(input_csv_file),
    As = A_count,
    Ts = T_count,
    coverage = coverage_sum,
    overall_editing_level = overall_editing_level
  )
  
    return(new_df)
}

# Path containing your filtered CSV files
csv_files <- list.files(
  "G:/Chapter4/Overall_ADAR_editing/Post_infection_males",
  pattern = "\\.csv$",
  full.names = TRUE
)

# Run over all files
list_of_new_dfs <- lapply(csv_files, calculate_values)
new_df <- do.call(rbind, list_of_new_dfs)
infection_stage_sex <- rep("Post_infection_males", nrow(new_df))
new_df$infection_stage_sex <- infection_stage_sex
# Save 
finalResult <- file.path(directoryPath, "Result_Post_males", "overall_editing_Post_infection_males.csv")
write.csv(new_df, finalResult, row.names = FALSE)

#===============================================================================
#===============================================================================
#===============================================================================

library(tidyverse)
library(ggpubr)

# Define the function to read a file and add the group labels
read_and_label_data <- function(file_path, stage, sex) {
  data <- read.csv(file_path)
  data <- data %>%
    mutate(
      Stage = stage, 
      Sex = sex
    ) 
  return(data)
}

# Path to .csv files with overall editing levels calulated for each sample
file_paths <- list(
  pre_male = "G:\\Chapter4\\Overall_ADAR_editing\\Result_pre_males\\overall_editing_pre_infection_males.csv",
  pre_female = "G:\\Chapter4\\Overall_ADAR_editing\\Result_Pre_females\\overall_editing_pre_infection_females.csv",
  mid_male = "G:\\Chapter4\\Overall_ADAR_editing\\Result_Mid_males\\overall_editing_Mid_infection_males.csv",
  mid_female = "G:\\Chapter4\\Overall_ADAR_editing\\Result_Mid_females\\overall_editing_Mid_infection_females.csv",
  post_male = "G:\\Chapter4\\Overall_ADAR_editing\\Result_Post_males\\overall_editing_Post_infection_males.csv",
  post_female = "G:\\Chapter4\\Overall_ADAR_editing\\Result_Post_females\\overall_editing_Post_infection_females.csv"
)

# Import and combine the abpve files to a single DF
full_data <- bind_rows(
  read_and_label_data(file_paths$pre_male, "pre", "male"),
  read_and_label_data(file_paths$pre_female, "pre", "female"),
  read_and_label_data(file_paths$mid_male, "mid", "male"),
  read_and_label_data(file_paths$mid_female, "mid", "female"),
  read_and_label_data(file_paths$post_male, "post", "male"),
  read_and_label_data(file_paths$post_female, "post", "female")
)

# Order stages 
full_data <- full_data %>%
  mutate(
    # Set the correct order for the stages
    Stage = factor(Stage, levels = c("pre", "mid", "post")),
    Sex = factor(Sex, levels = c("male", "female"))
  )
sex_colors <- c("male" = "royalblue", "female" = "hotpink")

# Stat
# ANOVA within sex across different stages of infection (stage effect grouped by sex)
anova_pvals <- full_data %>%
  group_by(Sex) %>%
  summarise(
   
    p = tryCatch(summary(aov(overall_editing_level ~ Stage, data = cur_data()))[[1]][["Pr(>F)"]][1],  # will andle the input file errors if any
                 error = function(e) NA_real_),
    y = max(overall_editing_level, na.rm = TRUE) * 1.05, 
    .groups = "drop"
  ) %>%
  mutate(label = paste0("ANOVA p = ", signif(p, 3)))


# plot
wilcox_y_pos <- max(full_data$overall_editing_level, na.rm = TRUE) * 1.0
p_combined_adar <- ggplot(full_data, aes(x = Stage, y = overall_editing_level, fill = Sex)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, position = position_dodge(width = 0.8)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), alpha = 0.5) +
  scale_fill_manual(values = sex_colors) +
  labs(
    y = "Overall ADAR editing levels", 
    x = "Stages of SARS-CoV-2 infection"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.position = c(0.9, 0.9) 
  ) +
  
  # 1.Comparison between sexes at each stage of infection
  stat_compare_means(aes(group = Sex),
                     method = "wilcox.test", 
                     label = "p.format",
                     label.y = wilcox_y_pos, 
                     size = 4,
                     vjust = -0.5) +
  
  # 2.Comparison within sex across different stages of infection
  geom_bracket(data = anova_pvals,
               aes(xmin = 1, xmax = 3, y.position = y, label = label, color = Sex),
               inherit.aes = FALSE,
               tip.length = 0.03, 
               size = 0.8) +
  scale_color_manual(values = sex_colors, guide = "none") 
print(p_combined_adar)

ggsave("overall_ADAR_editing_combined_plot.png", plot = p_combined_adar, width = 15, height = 10, dpi = 400)
