library(data.table)
#install.packages("readxl")
library(readxl)
library(dplyr)
library(ggpubr)



# Supplementary excel file (SUpp file 13) from Taru Tukiainen et al - contains list of XCI genes
xci_file <- "G:/Chapter4/XCI_analysis/NIHMS905235-supplement-supp_table13.xlsx"
xci <- read_excel(xci_file)

# Extract XCI gene 
xci_genes <- unique(na.omit(xci$`Gene name`))
length(xci_genes)
head(xci_genes)

# Folder of VCF files
edit_dir <- "G:/Chapter4/ch4_ctlfemale.89_subset/ch4_ctlfemale.89_subset/ch4_ctlfemale_89.snpEff_counts_csv/snpEff_counts_csv"
edit_files <- list.files(edit_dir, pattern = "\\.csv$", full.names = TRUE)
length(edit_files)  

# Output dir
out_dir <- file.path(edit_dir, "XCI_filtered")
dir.create(out_dir, showWarnings = FALSE)

for (f in edit_files) {
   df <- fread(f)
  df_xci <- df %>%
    filter(snpEff_geneName %in% xci_genes)
  out_file <- file.path(out_dir, basename(f))
  fwrite(df_xci, out_file)
}

#============= XCI ADAR edits - pre-females ==================================================
library(data.table)
library(dplyr)
library(readxl)

# # Supplementary excel file (SUpp file 13) from Taru Tukiainen et al - contains list of XCI genes
xci_file <- "G:/Chapter4/XCI_analysis/NIHMS905235-supplement-supp_table13.xlsx"

# Folder of VCF files
edit_dir <- "G:/Chapter4/ch4_ctlfemale.89_subset/ch4_ctlfemale.89_subset/ch4_ctlfemale_89.snpEff_counts_csv/snpEff_counts_csv"

# XCI list 
xci <- read_excel(xci_file)
xci_genes <- unique(na.omit(xci$`Gene name`))

# CSV files
edit_files <- list.files(edit_dir, pattern = "\\.csv$", full.names = TRUE)

# o/p dir
out_dir <- file.path(edit_dir, "XCI_ADAR_only")
dir.create(out_dir, showWarnings = FALSE)

# ---- Process files ----
for (f in edit_files) {
  
  df <- fread(f)
  df_filtered <- df %>%
    
    filter(snpEff_geneName %in% xci_genes) %>% # keep only XCI genes
    mutate(
      ADAR_edits = paste0(toupper(REF), toupper(ALT)) # create ADAR_edits column
    ) %>%
    filter(ADAR_edits %in% c("AG", "TC")) ## Filter for ADAR edits
    fwrite(df_filtered,
         file.path(out_dir, basename(f)))
}

# Summary file
summarise_ADAR_per_gene <- function(df, sample_id,
                                    sex = NA,
                                    stage = NA) {
  
  df %>%
    mutate(
      snpEff_geneName = as.character(snpEff_geneName),
      edited_reads = case_when(
        REF == "A" & ALT == "G" ~ G,
        REF == "T" & ALT == "C" ~ C,
        TRUE ~ NA_real_
      ),
      editing_level = edited_reads / Total
    ) %>%
    group_by(snpEff_geneName) %>%
    summarise(
      sample = sample_id,
      sex = sex,
      stage = stage,
      n_ADAR_sites = n(),
      total_depth = sum(Total, na.rm = TRUE),
      edited_reads = sum(edited_reads, na.rm = TRUE),
      mean_editing_level = mean(editing_level, na.rm = TRUE),
      weighted_editing_level = edited_reads / total_depth,
      .groups = "drop"
    )
}


# Directory containing XCI_ADAR_only filtered CSVs
in_dir <- "G:/Chapter4/ch4_ctlfemale.89_subset/ch4_ctlfemale.89_subset/ch4_ctlfemale_89.snpEff_counts_csv/snpEff_counts_csv/XCI_ADAR_only"

files <- list.files(in_dir, pattern = "\\.csv$", full.names = TRUE)
adar_gene_summary <- lapply(files, function(f) {
  
  df <- fread(f)
  sample_id <- tools::file_path_sans_ext(basename(f))
  
  summarise_ADAR_per_gene(
    df = df,
    sample_id = sample_id,
    sex = "female",
    stage = "pre-infection"
  )
}) %>% bind_rows()
table(adar_gene_summary$sex, adar_gene_summary$stage)
fwrite(adar_gene_summary,
       file.path(in_dir, "ADAR_per_gene_summary.csv"))

#============= XCI ADAR edits - pre-males ==================================================
library(data.table)
library(dplyr)
library(readxl)

# # Supplementary excel file (SUpp file 13) from Taru Tukiainen et al - contains list of XCI genes
xci_file <- "G:/Chapter4/XCI_analysis/NIHMS905235-supplement-supp_table13.xlsx"

# Folder of VCF files
edit_dir <- "G:\\Chapter4\\ch4_ctlmale.100_subset\\ch4_ctlmale.100_subset\\ch4_ctlmale_100_215.snpEff_counts_csv\\snpEff_counts_csv"

# XCI list 
xci <- read_excel(xci_file)
xci_genes <- unique(na.omit(xci$`Gene name`))

# CSV files
edit_files <- list.files(edit_dir, pattern = "\\.csv$", full.names = TRUE)

# o/p dir
out_dir <- file.path(edit_dir, "XCI_ADAR_only")
dir.create(out_dir, showWarnings = FALSE)

# ---- Process files ----
for (f in edit_files) {
  
  df <- fread(f)
  df_filtered <- df %>%
    
    filter(snpEff_geneName %in% xci_genes) %>% # keep only XCI genes
    mutate(
      ADAR_edits = paste0(toupper(REF), toupper(ALT)) # create ADAR_edits column
    ) %>%
    filter(ADAR_edits %in% c("AG", "TC")) ## Filter for ADAR edits
  fwrite(df_filtered,
         file.path(out_dir, basename(f)))
}

# Summary file
summarise_ADAR_per_gene <- function(df, sample_id,
                                    sex = NA,
                                    stage = NA) {
  
  df %>%
    mutate(
      snpEff_geneName = as.character(snpEff_geneName),
      edited_reads = case_when(
        REF == "A" & ALT == "G" ~ G,
        REF == "T" & ALT == "C" ~ C,
        TRUE ~ NA_real_
      ),
      editing_level = edited_reads / Total
    ) %>%
    group_by(snpEff_geneName) %>%
    summarise(
      sample = sample_id,
      sex = sex,
      stage = stage,
      n_ADAR_sites = n(),
      total_depth = sum(Total, na.rm = TRUE),
      edited_reads = sum(edited_reads, na.rm = TRUE),
      mean_editing_level = mean(editing_level, na.rm = TRUE),
      weighted_editing_level = edited_reads / total_depth,
      .groups = "drop"
    )
}


# Directory containing XCI_ADAR_only filtered CSVs
in_dir <- "G:\\Chapter4\\ch4_ctlmale.100_subset\\ch4_ctlmale.100_subset\\ch4_ctlmale_100_215.snpEff_counts_csv\\snpEff_counts_csv\\XCI_ADAR_only"

files <- list.files(in_dir, pattern = "\\.csv$", full.names = TRUE)
adar_gene_summary <- lapply(files, function(f) {
  
  df <- fread(f)
  sample_id <- tools::file_path_sans_ext(basename(f))
  
  summarise_ADAR_per_gene(
    df = df,
    sample_id = sample_id,
    sex = "female",
    stage = "pre-infection"
  )
}) %>% bind_rows()
table(adar_gene_summary$sex, adar_gene_summary$stage)
fwrite(adar_gene_summary,
       file.path(in_dir, "ADAR_per_gene_summary_pre_males.csv"))


#============= XCI ADAR edits - mid-females ==================================================
library(data.table)
library(dplyr)
library(readxl)

# # Supplementary excel file (SUpp file 13) from Taru Tukiainen et al - contains list of XCI genes
xci_file <- "G:/Chapter4/XCI_analysis/NIHMS905235-supplement-supp_table13.xlsx"

# Folder of VCF files
edit_dir <- "G:\\Chapter4\\ch4_midfemale.89_subset\\ch4_midfemale.89_subset\\ch4_midfemale_89.snpEff_counts_csv\\snpEff_counts_csv"

# XCI list 
xci <- read_excel(xci_file)
xci_genes <- unique(na.omit(xci$`Gene name`))

# CSV files
edit_files <- list.files(edit_dir, pattern = "\\.csv$", full.names = TRUE)

# o/p dir
out_dir <- file.path(edit_dir, "XCI_ADAR_only")
dir.create(out_dir, showWarnings = FALSE)

# ---- Process files ----
for (f in edit_files) {
  
  df <- fread(f)
  df_filtered <- df %>%
    
    filter(snpEff_geneName %in% xci_genes) %>% # keep only XCI genes
    mutate(
      ADAR_edits = paste0(toupper(REF), toupper(ALT)) # create ADAR_edits column
    ) %>%
    filter(ADAR_edits %in% c("AG", "TC")) ## Filter for ADAR edits
  fwrite(df_filtered,
         file.path(out_dir, basename(f)))
}

# Summary file
summarise_ADAR_per_gene <- function(df, sample_id,
                                    sex = NA,
                                    stage = NA) {
  
  df %>%
    mutate(
      snpEff_geneName = as.character(snpEff_geneName),
      edited_reads = case_when(
        REF == "A" & ALT == "G" ~ G,
        REF == "T" & ALT == "C" ~ C,
        TRUE ~ NA_real_
      ),
      editing_level = edited_reads / Total
    ) %>%
    group_by(snpEff_geneName) %>%
    summarise(
      sample = sample_id,
      sex = sex,
      stage = stage,
      n_ADAR_sites = n(),
      total_depth = sum(Total, na.rm = TRUE),
      edited_reads = sum(edited_reads, na.rm = TRUE),
      mean_editing_level = mean(editing_level, na.rm = TRUE),
      weighted_editing_level = edited_reads / total_depth,
      .groups = "drop"
    )
}


# Directory containing XCI_ADAR_only filtered CSVs
in_dir <- "G:\\Chapter4\\ch4_midfemale.89_subset\\ch4_midfemale.89_subset\\ch4_midfemale_89.snpEff_counts_csv\\snpEff_counts_csv\\XCI_ADAR_only"

files <- list.files(in_dir, pattern = "\\.csv$", full.names = TRUE)
adar_gene_summary <- lapply(files, function(f) {
  
  df <- fread(f)
  sample_id <- tools::file_path_sans_ext(basename(f))
  
  summarise_ADAR_per_gene(
    df = df,
    sample_id = sample_id,
    sex = "female",
    stage = "pre-infection"
  )
}) %>% bind_rows()
table(adar_gene_summary$sex, adar_gene_summary$stage)
fwrite(adar_gene_summary,
       file.path(in_dir, "ADAR_per_gene_summary_mid_females.csv"))


#============= XCI ADAR edits - mid-males ==================================================
library(data.table)
library(dplyr)
library(readxl)

# # Supplementary excel file (SUpp file 13) from Taru Tukiainen et al - contains list of XCI genes
xci_file <- "G:/Chapter4/XCI_analysis/NIHMS905235-supplement-supp_table13.xlsx"

# Folder of VCF files
edit_dir <- "G:\\Chapter4\\ch4_midmale.100_subset\\ch4_midmale.100_subset\\ch4_midmale_100_345.snpEff_counts_csv\\snpEff_counts_csv"

# XCI list 
xci <- read_excel(xci_file)
xci_genes <- unique(na.omit(xci$`Gene name`))

# CSV files
edit_files <- list.files(edit_dir, pattern = "\\.csv$", full.names = TRUE)

# o/p dir
out_dir <- file.path(edit_dir, "XCI_ADAR_only")
dir.create(out_dir, showWarnings = FALSE)

# ---- Process files ----
for (f in edit_files) {
  
  df <- fread(f)
  df_filtered <- df %>%
    
    filter(snpEff_geneName %in% xci_genes) %>% # keep only XCI genes
    mutate(
      ADAR_edits = paste0(toupper(REF), toupper(ALT)) # create ADAR_edits column
    ) %>%
    filter(ADAR_edits %in% c("AG", "TC")) ## Filter for ADAR edits
  fwrite(df_filtered,
         file.path(out_dir, basename(f)))
}

# Summary file
summarise_ADAR_per_gene <- function(df, sample_id,
                                    sex = NA,
                                    stage = NA) {
  
  df %>%
    mutate(
      snpEff_geneName = as.character(snpEff_geneName),
      edited_reads = case_when(
        REF == "A" & ALT == "G" ~ G,
        REF == "T" & ALT == "C" ~ C,
        TRUE ~ NA_real_
      ),
      editing_level = edited_reads / Total
    ) %>%
    group_by(snpEff_geneName) %>%
    summarise(
      sample = sample_id,
      sex = sex,
      stage = stage,
      n_ADAR_sites = n(),
      total_depth = sum(Total, na.rm = TRUE),
      edited_reads = sum(edited_reads, na.rm = TRUE),
      mean_editing_level = mean(editing_level, na.rm = TRUE),
      weighted_editing_level = edited_reads / total_depth,
      .groups = "drop"
    )
}


# Directory containing XCI_ADAR_only filtered CSVs
in_dir <- "G:\\Chapter4\\ch4_midmale.100_subset\\ch4_midmale.100_subset\\ch4_midmale_100_345.snpEff_counts_csv\\snpEff_counts_csv\\XCI_ADAR_only"

files <- list.files(in_dir, pattern = "\\.csv$", full.names = TRUE)
adar_gene_summary <- lapply(files, function(f) {
  
  df <- fread(f)
  sample_id <- tools::file_path_sans_ext(basename(f))
  
  summarise_ADAR_per_gene(
    df = df,
    sample_id = sample_id,
    sex = "Male",
    stage = "Mid-infection"
  )
}) %>% bind_rows()
table(adar_gene_summary$sex, adar_gene_summary$stage)
fwrite(adar_gene_summary,
       file.path(in_dir, "ADAR_per_gene_summary_mid_males.csv"))



#============= XCI ADAR edits - post-females ==================================================
library(data.table)
library(dplyr)
library(readxl)

# # Supplementary excel file (SUpp file 13) from Taru Tukiainen et al - contains list of XCI genes
xci_file <- "G:/Chapter4/XCI_analysis/NIHMS905235-supplement-supp_table13.xlsx"

# Folder of VCF files
edit_dir <- "G:\\Chapter4\\ch4_postfemale.56_subset\\ch4_postfemale.56_subset\\ch4_postfemale_56.snpEff_counts_csv\\snpEff_counts_csv"

# XCI list 
xci <- read_excel(xci_file)
xci_genes <- unique(na.omit(xci$`Gene name`))

# CSV files
edit_files <- list.files(edit_dir, pattern = "\\.csv$", full.names = TRUE)

# o/p dir
out_dir <- file.path(edit_dir, "XCI_ADAR_only")
dir.create(out_dir, showWarnings = FALSE)

# ---- Process files ----
for (f in edit_files) {
  
  df <- fread(f)
  df_filtered <- df %>%
    
    filter(snpEff_geneName %in% xci_genes) %>% # keep only XCI genes
    mutate(
      ADAR_edits = paste0(toupper(REF), toupper(ALT)) # create ADAR_edits column
    ) %>%
    filter(ADAR_edits %in% c("AG", "TC")) ## Filter for ADAR edits
  fwrite(df_filtered,
         file.path(out_dir, basename(f)))
}

# Summary file
summarise_ADAR_per_gene <- function(df, sample_id,
                                    sex = NA,
                                    stage = NA) {
  
  df %>%
    mutate(
      snpEff_geneName = as.character(snpEff_geneName),
      edited_reads = case_when(
        REF == "A" & ALT == "G" ~ G,
        REF == "T" & ALT == "C" ~ C,
        TRUE ~ NA_real_
      ),
      editing_level = edited_reads / Total
    ) %>%
    group_by(snpEff_geneName) %>%
    summarise(
      sample = sample_id,
      sex = sex,
      stage = stage,
      n_ADAR_sites = n(),
      total_depth = sum(Total, na.rm = TRUE),
      edited_reads = sum(edited_reads, na.rm = TRUE),
      mean_editing_level = mean(editing_level, na.rm = TRUE),
      weighted_editing_level = edited_reads / total_depth,
      .groups = "drop"
    )
}


# Directory containing XCI_ADAR_only filtered CSVs
in_dir <- "G:\\Chapter4\\ch4_postfemale.56_subset\\ch4_postfemale.56_subset\\ch4_postfemale_56.snpEff_counts_csv\\snpEff_counts_csv\\XCI_ADAR_only"

files <- list.files(in_dir, pattern = "\\.csv$", full.names = TRUE)
adar_gene_summary <- lapply(files, function(f) {
  
  df <- fread(f)
  sample_id <- tools::file_path_sans_ext(basename(f))
  
  summarise_ADAR_per_gene(
    df = df,
    sample_id = sample_id,
    sex = "Male",
    stage = "Mid-infection"
  )
}) %>% bind_rows()
table(adar_gene_summary$sex, adar_gene_summary$stage)
fwrite(adar_gene_summary,
       file.path(in_dir, "ADAR_per_gene_summary_mid_males.csv"))



#============= XCI ADAR edits - post-males ==================================================
library(data.table)
library(dplyr)
library(readxl)

# # Supplementary excel file (SUpp file 13) from Taru Tukiainen et al - contains list of XCI genes
xci_file <- "G:/Chapter4/XCI_analysis/NIHMS905235-supplement-supp_table13.xlsx"

# Folder of VCF files
edit_dir <- "G:\\Chapter4\\ch4_postmale.100_subset\\ch4_postmale.100_subset\\ch4_postmale_100_295.snpEff_counts_csv\\snpEff_counts_csv"

# XCI list 
xci <- read_excel(xci_file)
xci_genes <- unique(na.omit(xci$`Gene name`))

# CSV files
edit_files <- list.files(edit_dir, pattern = "\\.csv$", full.names = TRUE)

# o/p dir
out_dir <- file.path(edit_dir, "XCI_ADAR_only")
dir.create(out_dir, showWarnings = FALSE)

# ---- Process files ----
for (f in edit_files) {
  
  df <- fread(f)
  df_filtered <- df %>%
    
    filter(snpEff_geneName %in% xci_genes) %>% # keep only XCI genes
    mutate(
      ADAR_edits = paste0(toupper(REF), toupper(ALT)) # create ADAR_edits column
    ) %>%
    filter(ADAR_edits %in% c("AG", "TC")) ## Filter for ADAR edits
  fwrite(df_filtered,
         file.path(out_dir, basename(f)))
}

# Summary file
summarise_ADAR_per_gene <- function(df, sample_id,
                                    sex = NA,
                                    stage = NA) {
  
  df %>%
    mutate(
      snpEff_geneName = as.character(snpEff_geneName),
      edited_reads = case_when(
        REF == "A" & ALT == "G" ~ G,
        REF == "T" & ALT == "C" ~ C,
        TRUE ~ NA_real_
      ),
      editing_level = edited_reads / Total
    ) %>%
    group_by(snpEff_geneName) %>%
    summarise(
      sample = sample_id,
      sex = sex,
      stage = stage,
      n_ADAR_sites = n(),
      total_depth = sum(Total, na.rm = TRUE),
      edited_reads = sum(edited_reads, na.rm = TRUE),
      mean_editing_level = mean(editing_level, na.rm = TRUE),
      weighted_editing_level = edited_reads / total_depth,
      .groups = "drop"
    )
}


# Directory containing XCI_ADAR_only filtered CSVs
in_dir <- "G:\\Chapter4\\ch4_postmale.100_subset\\ch4_postmale.100_subset\\ch4_postmale_100_295.snpEff_counts_csv\\snpEff_counts_csv\\XCI_ADAR_only"

files <- list.files(in_dir, pattern = "\\.csv$", full.names = TRUE)
adar_gene_summary <- lapply(files, function(f) {
  
  df <- fread(f)
  sample_id <- tools::file_path_sans_ext(basename(f))
  
  summarise_ADAR_per_gene(
    df = df,
    sample_id = sample_id,
    sex = "Male",
    stage = "Mid-infection"
  )
}) %>% bind_rows()
table(adar_gene_summary$sex, adar_gene_summary$stage)
fwrite(adar_gene_summary,
       file.path(in_dir, "ADAR_per_gene_summary_post_males.csv"))


# ================== plot======================================================

Pre_females <- read.csv("G:\\Chapter4\\ch4_ctlfemale.89_subset\\ch4_ctlfemale.89_subset\\ch4_ctlfemale_89.snpEff_counts_csv\\snpEff_counts_csv\\XCI_ADAR_only\\ADAR_per_gene_summary.csv")
Pre_males <- read.csv("G:\\Chapter4\\ch4_ctlmale.100_subset\\ch4_ctlmale.100_subset\\ch4_ctlmale_100_215.snpEff_counts_csv\\snpEff_counts_csv\\XCI_ADAR_only\\ADAR_per_gene_summary_pre_males.csv")

library(ggplot2)
library(dplyr)

# Assuming you have combined table: adar_combined
adar_combined <- rbind(Pre_females, Pre_males)

# Global weighted editing per sample
global_summary <- adar_combined %>%
  group_by(sample, sex) %>%
  summarise(global_XCI_editing = sum(edited_reads)/sum(total_depth),
            total_XCI_ADAR_sites = sum(n_ADAR_sites),
            .groups = "drop")
library(ggpubr)


library(effsize)
cd <- cliff.delta(
  global_XCI_editing ~ sex,
  data = global_summary
)

cd
delta_value <- round(cd$estimate, 3)

ggplot(global_summary, aes(x = sex, y = global_XCI_editing, fill = sex)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.7) +
  scale_fill_manual(values = c("Male" = "royalblue", "female" = "hotpink")) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format"
  ) +
  annotate(
    "text",
    x = 1.5,
    y = max(global_summary$global_XCI_editing, na.rm = TRUE) * 1.05,
    label = paste0("Cliff's \u03B4 = ", delta_value),
    size = 5
  ) +
  labs(
    y = "Weighted ADAR editing pre-infection (XCI genes)",
    x = "Sex"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_blank(),
    legend.position = c(0.9, 0.9)
  )

ggsave("G:\\Chapter4\\XCI_analysis\\Pre_males_females_weighted_average.png", width = 8, height = 8)


#=========================

library(dplyr)

adar_pre <- adar_combined %>%
  filter(stage == "pre-infection")

gene_mean_edits <- adar_pre %>%
  group_by(snpEff_geneName, sex) %>%
  summarise(
    mean_ADAR_sites = mean(n_ADAR_sites, na.rm = TRUE),
    .groups = "drop"
  )
library(tidyr)

gene_scatter <- gene_mean_edits %>%
  pivot_wider(
    names_from = sex,
    values_from = mean_ADAR_sites,
    values_fill = 0
  )

library(ggplot2)
library(effsize)


ggplot(gene_scatter,
       aes(x = Male, y = female)) +
  geom_point(alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "red") +
  labs(x = "Mean number of ADAR edits per gene (pre-infection males)",
    y = "Mean number of ADAR edits per gene (Mid-infection females)"      ) +
  theme_bw() +
theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  
      axis.text.y = element_text(size = 14),  
      axis.title.x = element_text(size = 16),  
      axis.title.y = element_text(size = 16),  
      legend.text = element_text(size = 14), 
      legend.title = element_blank()) # Remove legend title
      #legend.position = c(0.9, 0.9))

ggsave("G:\\Chapter4\\XCI_analysis\\Pre_males_females.png", width = 8, height = 8)


tlr8_df <- adar_combined_post %>%
  filter(snpEff_geneName == "TLR8")


#-------------------------------------------------------------------------------
Mid_females <- read.csv("G:\\Chapter4\\ch4_midfemale.89_subset\\ch4_midfemale.89_subset\\ch4_midfemale_89.snpEff_counts_csv\\snpEff_counts_csv\\XCI_ADAR_only\\ADAR_per_gene_summary_mid_females.csv")
Mid_males <- read.csv("G:\\Chapter4\\ch4_midmale.100_subset\\ch4_midmale.100_subset\\ch4_midmale_100_345.snpEff_counts_csv\\snpEff_counts_csv\\XCI_ADAR_only\\ADAR_per_gene_summary_mid_males.csv")

# Assuming you have combined table: adar_combined
adar_combined_mid <- rbind(Mid_females, Mid_males)

global_summary <- adar_combined_mid %>%
  group_by(sample, sex) %>%
  summarise(global_XCI_editing = sum(edited_reads)/sum(total_depth),
            total_XCI_ADAR_sites = sum(n_ADAR_sites),
            .groups = "drop")
cd <- cliff.delta(
  global_XCI_editing ~ sex,
  data = global_summary
)

cd
delta_value <- round(cd$estimate, 3)


ggplot(global_summary, aes(x = sex, y = global_XCI_editing, fill = sex)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.7) +
  scale_fill_manual(values = c("Male" = "royalblue", "female" = "hotpink")) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format"
  ) +
  annotate(
    "text",
    x = 1.5,
    y = max(global_summary$global_XCI_editing, na.rm = TRUE) * 1.05,
    label = paste0("Cliff's \u03B4 = ", delta_value),
    size = 5
  ) +
  labs(
    y = "Weighted ADAR editing mid-infection (XCI genes)",
    x = "Sex"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_blank(),
    legend.position = c(0.9, 0.9)
  )

ggsave("G:\\Chapter4\\XCI_analysis\\mid_males_females_weighted_average.png", width = 8, height = 8)
library(dplyr)

adar_mid <- adar_combined_mid %>%
  filter(stage == "Mid-infection")

gene_mean_edits <- adar_mid %>%
  group_by(snpEff_geneName, sex) %>%
  summarise(
    mean_ADAR_sites = mean(n_ADAR_sites, na.rm = TRUE),
    .groups = "drop"
  )
library(tidyr)

gene_scatter <- gene_mean_edits %>%
  pivot_wider(
    names_from = sex,
    values_from = mean_ADAR_sites,
    values_fill = 0
  )

library(ggplot2)

ggplot(gene_scatter,
       aes(x = Male, y = female)) +
  geom_point(alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "red") +
  labs(x = "Mean number of ADAR edits per gene (mid-infection males)",
       y = "Mean number of ADAR edits per gene (mid-infection females)"      ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  
        axis.text.y = element_text(size = 14),  
        axis.title.x = element_text(size = 16),  
        axis.title.y = element_text(size = 16),  
        legend.text = element_text(size = 14), 
        legend.title = element_blank()) # Remove legend title
#legend.position = c(0.9, 0.9))

ggsave("G:\\Chapter4\\XCI_analysis\\mid_males_females.png", width = 8, height = 8)



#-------------------------------------------------------------------------------
Post_females <- read.csv("G:\\Chapter4\\ch4_postfemale.56_subset\\ch4_postfemale.56_subset\\ch4_postfemale_56.snpEff_counts_csv\\snpEff_counts_csv\\XCI_ADAR_only\\ADAR_per_gene_summary_post_females.csv")
Post_males <- read.csv("G:\\Chapter4\\ch4_postmale.100_subset\\ch4_postmale.100_subset\\ch4_postmale_100_295.snpEff_counts_csv\\snpEff_counts_csv\\XCI_ADAR_only\\ADAR_per_gene_summary_post_males.csv")

# Assuming you have combined table: adar_combined
adar_combined_post <- rbind(Post_females, Post_males)
View(adar_combined_post)
global_summary <- adar_combined_post %>%
  group_by(sample, sex) %>%
  summarise(global_XCI_editing = sum(edited_reads)/sum(total_depth),
            total_XCI_ADAR_sites = sum(n_ADAR_sites),
            .groups = "drop")

cd <- cliff.delta(
  global_XCI_editing ~ sex,
  data = global_summary
)

cd
delta_value <- round(cd$estimate, 3)


ggplot(global_summary, aes(x = sex, y = global_XCI_editing, fill = sex)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.7) +
  scale_fill_manual(values = c("Male" = "royalblue", "Female" = "hotpink")) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format"
  ) +
  annotate(
    "text",
    x = 1.5,
    y = max(global_summary$global_XCI_editing, na.rm = TRUE) * 1.05,
    label = paste0("Cliff's \u03B4 = ", delta_value),
    size = 5
  ) +
  labs(
    y = "Weighted ADAR editing post-infection (XCI genes)",
    x = "Sex"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_blank(),
    legend.position = c(0.9, 0.9)
  )

ggsave("G:\\Chapter4\\XCI_analysis\\post_males_females_weighted_average.png", width = 8, height = 8)


library(dplyr)

adar_post <- adar_combined_post %>%
  filter(stage == "post-infection")
View(adar_post)
gene_mean_edits <- adar_combined_post %>%
  group_by(snpEff_geneName, sex) %>%
  summarise(
    mean_ADAR_sites = mean(n_ADAR_sites, na.rm = TRUE),
    .groups = "drop"
  )
View(gene_mean_edits)
library(tidyr)

gene_scatter <- gene_mean_edits %>%
  pivot_wider(
    names_from = sex,
    values_from = mean_ADAR_sites,
    values_fill = 0
  )

library(ggplot2)

ggplot(gene_scatter,
       aes(x = Male, y = Female)) +
  geom_point(alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "red") +
  labs(x = "Mean number of ADAR edits per gene (post-infection males)",
       y = "Mean number of ADAR edits per gene (post-infection females)"      ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  
        axis.text.y = element_text(size = 14),  
        axis.title.x = element_text(size = 16),  
        axis.title.y = element_text(size = 16),  
        legend.text = element_text(size = 14), 
        legend.title = element_blank()) # Remove legend title
#legend.position = c(0.9, 0.9))

ggsave("G:\\Chapter4\\XCI_analysis\\Post_males_females.png", width = 8, height = 8)















