# ADAR expression
library(ggplot2)
library(dplyr)
library(ggpubr)
library(tidyverse)
library(data.table)

#WD
setwd("G:\\Chapter4")

# Path to all the counts folder
folders <- list(
  pre_male   = "G:\\Chapter4\\ch4_ctlmale.100_subset\\ch4_ctlmale.100_subset\\ch4_ctlmale_100_215.counts-etc\\counts",
  pre_female = "G:\\Chapter4\\ch4_ctlfemale.89_subset\\ch4_ctlfemale.89_subset\\ch4_ctlfemale_89.counts-etc\\counts",
  mid_male   = "G:\\Chapter4\\ch4_midmale.100_subset\\ch4_midmale.100_subset\\ch4_midmale_100_345.counts-etc\\counts",
  mid_female = "G:\\Chapter4\\ch4_midfemale.89_subset\\ch4_midfemale.89_subset\\ch4_midfemale_89.counts-etc\\counts",
  post_male  = "G:\\Chapter4\\ch4_postmale.100_subset\\ch4_postmale.100_subset\\ch4_postmale_100_295.counts-etc\\counts",
  post_female= "G:\\Chapter4\\ch4_postfemale.56_subset\\ch4_postfemale.56_subset\\ch4_postfemale_56.counts-etc\\counts"
)


##------------- ADAR1 ----------------------------------------------------------
# Function to read one file and extract ADAR1 expression
extract_adar1 <- function(file, stage, sex) {
  df <- fread(file)
  df %>%
    filter(`Gene Name` == "ADAR") %>%   
    mutate(
      sample_id = basename(file),
      stage = stage,
      sex = sex
    ) %>%
    select(sample_id, stage, sex, TPM)
}


# Loop all the folders
adar1_list <- list()

for (nm in names(folders)) {
  path <- folders[[nm]]
  files <- list.files(path, pattern = "\\.tab", full.names = TRUE)  
  # split folder name to get stage + sex
  parts <- strsplit(nm, "_")[[1]]
  stage <- parts[1]
  sex <- parts[2]
  tmp <- map_dfr(files, extract_adar1, stage = stage, sex = sex)
  adar1_list[[nm]] <- tmp
}

adar1_expr <- bind_rows(adar1_list)
head(adar1_expr)




# Order stages
adar1_expr <- adar1_expr %>%
  mutate(
    stage = factor(stage, levels = c("pre", "mid", "post")),
    sex = factor(sex, levels = c("male", "female"))
  )

# Colours
sex_colors <- c("male" = "royalblue", "female" = "hotpink")

# Run ANOVA for stage specific differences within each sex
anova_pvals <- adar1_expr %>%
  group_by(sex) %>%
  summarise(
    p = tryCatch(summary(aov(TPM ~ stage))[[1]][["Pr(>F)"]][1], error = function(e) NA_real_),
    y = max(TPM, na.rm = TRUE) * 1.05,
    .groups = "drop"
  ) %>%
  mutate(label = paste0("p = ", signif(p, 3)))



# Plot
ggplot(adar1_expr, aes(x = stage, y = TPM, fill = sex)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, position = position_dodge(width = 0.8)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), alpha = 0.5) +
  scale_fill_manual(values = sex_colors) +
  labs(y = "ADAR1 expression in TPM", x = "Infection stage") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.position = c(0.9, 0.9)) +
  # Sex specific differences at each stage between males and female
  stat_compare_means(aes(group = sex),
                     method = "t.test",
                     label = "p.format",
                     label.y = 600) +
  # ANOVA for differences across stages within each sex
  geom_bracket(data = anova_pvals,
               aes(xmin = 1, xmax = 3, y.position = y, label = label, color = sex),
               inherit.aes = FALSE,
               tip.length = 0.03, size = 0.8) +
  scale_color_manual(values = sex_colors, guide = "none") 
ggsave("adar1.png", width = 8, height = 6, dpi = 350)  

##------------- ADAR2 ----------------------------------------------------------
# Function to read one file and extract ADAR1 expression
extract_adar2 <- function(file, stage, sex) {
  df <- fread(file)
  df %>%
    filter(`Gene Name` == "ADARB1") %>%   
    mutate(
      sample_id = basename(file),
      stage = stage,
      sex = sex
    ) %>%
    select(sample_id, stage, sex, TPM)
}


# Loop all the folders
adar2_list <- list()

for (nm in names(folders)) {
  path <- folders[[nm]]
  files <- list.files(path, pattern = "\\.tab", full.names = TRUE)  
  # split folder name to get stage + sex
  parts <- strsplit(nm, "_")[[1]]
  stage <- parts[1]
  sex <- parts[2]
  tmp <- map_dfr(files, extract_adar2, stage = stage, sex = sex)
  adar2_list[[nm]] <- tmp
}

adar2_expr <- bind_rows(adar2_list)
head(adar2_expr)




# Order stages
adar2_expr <- adar2_expr %>%
  mutate(
    stage = factor(stage, levels = c("pre", "mid", "post")),
    sex = factor(sex, levels = c("male", "female"))
  )

# Colours
sex_colors <- c("male" = "royalblue", "female" = "hotpink")

# Run ANOVA for stage specific differences within each sex
anova_pvals <- adar2_expr %>%
  group_by(sex) %>%
  summarise(
    p = tryCatch(summary(aov(TPM ~ stage))[[1]][["Pr(>F)"]][1], error = function(e) NA_real_),
    y = max(TPM, na.rm = TRUE) * 1.5,
    .groups = "drop"
  ) %>%
  mutate(label = paste0("p = ", signif(p, 3)))



# Plot
ggplot(adar2_expr, aes(x = stage, y = TPM, fill = sex)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, position = position_dodge(width = 0.8)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), alpha = 0.5) +
  scale_fill_manual(values = sex_colors) +
  labs(y = "ADAR2 expression in TPM", x = "Infection stage") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.position = c(0.9, 0.9)) +
  # Sex specific differences at each stage between males and female
  stat_compare_means(aes(group = sex),
                     method = "t.test",
                     label = "p.format",
                     label.y = 14) +
  # ANOVA for differences across stages within each sex
  geom_bracket(data = anova_pvals,
               aes(xmin = 1, xmax = 3, y.position = y, label = label, color = sex),
               inherit.aes = FALSE,
               tip.length = 0.03, size = 0.8) +
  scale_color_manual(values = sex_colors, guide = "none") 
ggsave("adar2.png", width = 8, height = 7, dpi = 350)  

##------------- adar3 ----------------------------------------------------------
# Function to read one file and extract ADAR1 expression
extract_adar3 <- function(file, stage, sex) {
  df <- fread(file)
  df %>%
    filter(`Gene Name` == "ADARB2") %>%   
    mutate(
      sample_id = basename(file),
      stage = stage,
      sex = sex
    ) %>%
    select(sample_id, stage, sex, TPM)
}


# Loop all the folders
adar3_list <- list()

for (nm in names(folders)) {
  path <- folders[[nm]]
  files <- list.files(path, pattern = "\\.tab", full.names = TRUE)  
  # split folder name to get stage + sex
  parts <- strsplit(nm, "_")[[1]]
  stage <- parts[1]
  sex <- parts[2]
  tmp <- map_dfr(files, extract_adar3, stage = stage, sex = sex)
  adar3_list[[nm]] <- tmp
}

adar3_expr <- bind_rows(adar3_list)
head(adar3_expr)


# Order stages
adar3_expr <- adar3_expr %>%
  mutate(
    stage = factor(stage, levels = c("pre", "mid", "post")),
    sex = factor(sex, levels = c("male", "female"))
  )

# Colours
sex_colors <- c("male" = "royalblue", "female" = "hotpink")

# Run ANOVA for stage specific differences within each sex
anova_pvals <- adar3_expr %>%
  group_by(sex) %>%
  summarise(
    p = tryCatch(summary(aov(TPM ~ stage))[[1]][["Pr(>F)"]][1], error = function(e) NA_real_),
    y = max(TPM, na.rm = TRUE) * 1.5,
    .groups = "drop"
  ) %>%
  mutate(label = paste0("p = ", signif(p, 3)))



# Plot
ggplot(adar3_expr, aes(x = stage, y = TPM, fill = sex)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, position = position_dodge(width = 0.8)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), alpha = 0.5) +
  scale_fill_manual(values = sex_colors) +
  labs(y = "ADAR3 expression in TPM", x = "Infection stage") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.position = c(0.9, 0.9)) +
  # Sex specific differences at each stage between males and female
  stat_compare_means(aes(group = sex),
                     method = "t.test",
                     label = "p.format",
                     label.y = 12) +
  # ANOVA for differences across stages within each sex
  geom_bracket(data = anova_pvals,
               aes(xmin = 1, xmax = 3, y.position = y, label = label, color = sex),
               inherit.aes = FALSE,
               tip.length = 0.03, size = 0.8) +
  scale_color_manual(values = sex_colors, guide = "none") 
ggsave("adar3.png", width = 8, height = 6, dpi = 350)  



