library(tidyverse)
library(ggpubr)
library(ggpmisc)
library(reshape2)
library(patchwork)
library(RColorBrewer)
setwd("G:\\Chapter4\\correlations_ADAR_editing")

# ADAR1
# Folders with ADAR TPM counts
pre_male   <- "G:/Chapter4/ch4_ctlmale.100_subset/ch4_ctlmale.100_subset/ch4_ctlmale_100_215.counts-etc/counts/"
pre_female <- "G:/Chapter4/ch4_ctlfemale.89_subset/ch4_ctlfemale.89_subset/ch4_ctlfemale_89.counts-etc/counts/"
mid_male   <- "G:/Chapter4/ch4_midmale.100_subset/ch4_midmale.100_subset/ch4_midmale_100_345.counts-etc/counts/"
mid_female <- "G:/Chapter4/ch4_midfemale.89_subset/ch4_midfemale.89_subset/ch4_midfemale_89.counts-etc/counts/"
post_male  <- "G:/Chapter4/ch4_postmale.100_subset/ch4_postmale.100_subset/ch4_postmale_100_295.counts-etc/counts/"
post_female<- "G:/Chapter4/ch4_postfemale.56_subset/ch4_postfemale.56_subset/ch4_postfemale_56.counts-etc/counts/"

folder_list <- list(
  Pre_Male = pre_male,
  Pre_Female = pre_female,
  Mid_Male = mid_male,
  Mid_Female = mid_female,
  Post_Male = post_male,
  Post_Female = post_female
)

# Extract ADAR TPM

extract_ADAR1 <- function(folder, group_name) {
  files <- list.files(folder, pattern = "\\.tab$", full.names = TRUE)
  
  map_df(files, function(f) {
    tab <- read.delim(f, stringsAsFactors = FALSE)
    
    adar_row <- tab %>% filter(Gene.Name %in% c("ADAR", "ADAR1"))
    if (nrow(adar_row) == 0) return(NULL)
    
    tibble(
      sample = gsub(".tab$", "", basename(f)),
      ADAR1_TPM = adar_row$TPM,
      group = group_name
    )
  })
}

# ADAR exp
adar_expression <- map2_df(folder_list, names(folder_list), extract_ADAR1)

# .csv files with pre caluclated overall ADAR editing values
pre_female  <- "G:/Chapter4/Overall_ADAR_editing/Result_Pre_females/overall_editing_pre_infection_females.csv"
pre_male    <- "G:/Chapter4/Overall_ADAR_editing/Result_pre_males/overall_editing_pre_infection_males.csv"

mid_female  <- "G:/Chapter4/Overall_ADAR_editing/Result_Mid_females/overall_editing_Mid_infection_females.csv"
mid_male    <- "G:/Chapter4/Overall_ADAR_editing/Result_Mid_males/overall_editing_Mid_infection_males.csv"

post_female <- "G:/Chapter4/Overall_ADAR_editing/Result_Post_females/overall_editing_Post_infection_females.csv"
post_male   <- "G:/Chapter4/Overall_ADAR_editing/Result_Post_males/overall_editing_Post_infection_males.csv"

editing_list <- list(
  Pre_Female  = read.csv(pre_female),
  Pre_Male    = read.csv(pre_male),
  Mid_Female  = read.csv(mid_female),
  Mid_Male    = read.csv(mid_male),
  Post_Female = read.csv(post_female),
  Post_Male   = read.csv(post_male)
)

editing_df <- bind_rows(editing_list, .id = "group")

### Clean sample names
editing_df$sample_clean <- editing_df$sample %>% 
  gsub("\\.csv$|\\.txt$", "", .) %>% 
  gsub("_.*$", "", .)

# Merge files with expression and editing levels

merged <- left_join(
  adar_expression,
  editing_df,
  by = c("sample" = "sample_clean"),
  suffix = c("_expr", "_edit")
)

# Modify/add colums
merged <- merged %>%
  mutate(
    stage = case_when(
      grepl("Pre", group_expr) ~ "Pre-infection",
      grepl("Mid", group_expr) ~ "Mid-infection",
      grepl("Post", group_expr) ~ "Post-infection"
    ),
    sex = case_when(
      grepl("Male", group_expr) ~ "Male",
      grepl("Female", group_expr) ~ "Female"
    )
  )

# Order stages
merged$stage <- factor(merged$stage, 
                       levels = c("Pre-infection", "Mid-infection", "Post-infection"))

# spearman correlations
corr_table <- merged %>%
  group_by(stage, sex) %>%
  summarise(
    rho = cor(ADAR1_TPM, overall_editing_level, method = "spearman", use = "complete.obs"),
    p_value = cor.test(ADAR1_TPM, overall_editing_level, method = "spearman")$p.value,
    n = n()
  )

write.csv(corr_table, "Spearman_correlation_results.csv", row.names = FALSE)

# plot with spearman

# p <- ggplot(merged, aes(x = ADAR1_TPM, y = overall_editing_level)) +
#   geom_point(size = 3, alpha = 0.8) +
#   geom_smooth(method = "lm", se = FALSE, color = "black") +
#   # Remove the 'aes(label = ...)' line entirely.
#   # The default label is usually sufficient and avoids parsing issues.
#   stat_cor(
#     method = "spearman",
#     label.x.npc = "left",
#     label.y.npc = "top",
#     size = 4
#   ) +
#   facet_grid(stage ~ sex) +
#   theme_bw(base_size = 14) +
#   labs(
#     x = "ADAR1 TPM",
#     y = "Overall ADAR Editing Level",
#     title = "Correlation between ADAR1 expression and ADAR editing"
#   )
# 
# ggsave("ADAR1_correlation_facet_plot.png", p, width = 10, height = 8)
# 


# plot with spearman  
p <- ggplot(merged, aes(x = ADAR1_TPM, y = overall_editing_level)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), # This will give R2
    method = "spearman", 
    label.x.npc = "left",
    label.y.npc = "top",
    size = 4
  ) +
  facet_grid(stage ~ sex) +
  theme_bw(base_size = 14) +
  labs(
    x = "ADAR1 TPM",
    y = "Overall ADAR Editing Level",
    title = "Correlation between ADAR1 expression and ADAR editing"
  )
# Correlation heatmap
heatmap_df <- corr_table %>%
  dcast(stage ~ sex, value.var = "rho")

rownames(heatmap_df) <- heatmap_df$stage
heatmap_df$stage <- NULL

p_heat <- ggplot(melt(as.matrix(heatmap_df)), 
                 aes(Var2, Var1, fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 2)), size = 6) +
  scale_fill_distiller(palette = "RdYlBu", direction = -1) +
  labs(
    x = "Sex",
    y = "Stage"
  ) +
  theme_minimal(base_size = 14)

ggsave("Spearman_heatmap.png", p_heat, width = 6, height = 5)

write.csv(merged, "ADAR1_TPM_and_editing_merged.csv", row.names = FALSE)

# Combine the two plots 
combined_plot <- p | p_heat # this is using lib patchwork
combined_plot

# Save 
ggsave("Combined_Correlation_Plot.png", combined_plot, width = 13, height = 5, dpi = 350)
  #=========================================== ADAR2==============================
  # ADAR2
  # Folders with ADAR TPM counts
  pre_male   <- "G:/Chapter4/ch4_ctlmale.100_subset/ch4_ctlmale.100_subset/ch4_ctlmale_100_215.counts-etc/counts/"
  pre_female <- "G:/Chapter4/ch4_ctlfemale.89_subset/ch4_ctlfemale.89_subset/ch4_ctlfemale_89.counts-etc/counts/"
  mid_male   <- "G:/Chapter4/ch4_midmale.100_subset/ch4_midmale.100_subset/ch4_midmale_100_345.counts-etc/counts/"
  mid_female <- "G:/Chapter4/ch4_midfemale.89_subset/ch4_midfemale.89_subset/ch4_midfemale_89.counts-etc/counts/"
  post_male  <- "G:/Chapter4/ch4_postmale.100_subset/ch4_postmale.100_subset/ch4_postmale_100_295.counts-etc/counts/"
  post_female<- "G:/Chapter4/ch4_postfemale.56_subset/ch4_postfemale.56_subset/ch4_postfemale_56.counts-etc/counts/"
  
  folder_list <- list(
    Pre_Male = pre_male,
    Pre_Female = pre_female,
    Mid_Male = mid_male,
    Mid_Female = mid_female,
    Post_Male = post_male,
    Post_Female = post_female
  )
  
  # Extract ADAR TPM
  
  extract_ADAR2 <- function(folder, group_name) {
    files <- list.files(folder, pattern = "\\.tab$", full.names = TRUE)
    
    map_df(files, function(f) {
      tab <- read.delim(f, stringsAsFactors = FALSE)
      
      adar_row <- tab %>% filter(Gene.Name %in% c("ADAR2", "ADARB1"))
      if (nrow(adar_row) == 0) return(NULL)
      
      tibble(
        sample = gsub(".tab$", "", basename(f)),
        ADAR2_TPM = adar_row$TPM,
        group = group_name
      )
    })
  }
  
  # ADAR exp
  adar_expression <- map2_df(folder_list, names(folder_list), extract_ADAR2)
  
  # .csv files with pre-caluclated overall ADAR editing values
  pre_female  <- "G:/Chapter4/Overall_ADAR_editing/Result_Pre_females/overall_editing_pre_infection_females.csv"
  pre_male    <- "G:/Chapter4/Overall_ADAR_editing/Result_pre_males/overall_editing_pre_infection_males.csv"
  
  mid_female  <- "G:/Chapter4/Overall_ADAR_editing/Result_Mid_females/overall_editing_Mid_infection_females.csv"
  mid_male    <- "G:/Chapter4/Overall_ADAR_editing/Result_Mid_males/overall_editing_Mid_infection_males.csv"
  
  post_female <- "G:/Chapter4/Overall_ADAR_editing/Result_Post_females/overall_editing_Post_infection_females.csv"
  post_male   <- "G:/Chapter4/Overall_ADAR_editing/Result_Post_males/overall_editing_Post_infection_males.csv"
  
  editing_list <- list(
    Pre_Female  = read.csv(pre_female),
    Pre_Male    = read.csv(pre_male),
    Mid_Female  = read.csv(mid_female),
    Mid_Male    = read.csv(mid_male),
    Post_Female = read.csv(post_female),
    Post_Male   = read.csv(post_male)
  )
  
  editing_df <- bind_rows(editing_list, .id = "group")
  
  ### Clean sample names
  editing_df$sample_clean <- editing_df$sample %>% 
    gsub("\\.csv$|\\.txt$", "", .) %>% 
    gsub("_.*$", "", .)
  
  # Merge files with expression and editing levels
  merged <- left_join(
    adar_expression,
    editing_df,
    by = c("sample" = "sample_clean"),
    suffix = c("_expr", "_edit")
  )
  
  # Modify/add colums
  merged <- merged %>%
    mutate(
      stage = case_when(
        grepl("Pre", group_expr) ~ "Pre-infection",
        grepl("Mid", group_expr) ~ "Mid-infection",
        grepl("Post", group_expr) ~ "Post-infection"
      ),
      sex = case_when(
        grepl("Male", group_expr) ~ "Male",
        grepl("Female", group_expr) ~ "Female"
      )
    )
  
  # Order stages
  merged$stage <- factor(merged$stage, 
                         levels = c("Pre-infection", "Mid-infection", "Post-infection"))
  
  # spearman correlations
  corr_table <- merged %>%
    group_by(stage, sex) %>%
    summarise(
      rho = cor(ADAR2_TPM, overall_editing_level, method = "spearman", use = "complete.obs"),
      p_value = cor.test(ADAR2_TPM, overall_editing_level, method = "spearman")$p.value,
      n = n()
    )
  
  write.csv(corr_table, "Spearman_correlation_results_adar2.csv", row.names = FALSE)
  
  # plot with spearman
  
  # p <- ggplot(merged, aes(x = adar2_TPM, y = overall_editing_level)) +
  #   geom_point(size = 3, alpha = 0.8) +
  #   geom_smooth(method = "lm", se = FALSE, color = "black") +
  #   # Remove the 'aes(label = ...)' line entirely.
  #   # The default label is usually sufficient and avoids parsing issues.
  #   stat_cor(
  #     method = "spearman",
  #     label.x.npc = "left",
  #     label.y.npc = "top",
  #     size = 4
  #   ) +
  #   facet_grid(stage ~ sex) +
  #   theme_bw(base_size = 14) +
  #   labs(
  #     x = "adar2 TPM",
  #     y = "Overall ADAR Editing Level",
  #     title = "Correlation between adar2 expression and ADAR editing"
  #   )
  # 
  # ggsave("adar2_correlation_facet_plot.png", p, width = 10, height = 8)
  # 
  
  
  # plot with spearman  ##r2 
  p <- ggplot(merged, aes(x = ADAR2_TPM, y = overall_editing_level)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    stat_cor(
      aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), # This will give R2
      method = "spearman", 
      label.x.npc = "left",
      label.y.npc = "top",
      size = 4
    ) +
    facet_grid(stage ~ sex) +
    theme_bw(base_size = 14) +
    labs(
      x = "ADAR2 TPM",
      y = "Overall ADAR Editing Level",
      title = "Correlation between ADAR2 expression and ADAR editing"
    )
  # Correlation heatmap
  heatmap_df <- corr_table %>%
    dcast(stage ~ sex, value.var = "rho")
  
  rownames(heatmap_df) <- heatmap_df$stage
  heatmap_df$stage <- NULL
  
  p_heat <- ggplot(melt(as.matrix(heatmap_df)), 
                   aes(Var2, Var1, fill = value)) +
    geom_tile() +
    geom_text(aes(label = round(value, 2)), size = 6) +
    scale_fill_distiller(palette = "RdYlBu", direction = -1) +
    labs(
      x = "Sex",
      y = "Stage"
    ) +
    theme_minimal(base_size = 14)
  
  ggsave("Spearman_heatmap_adar2.png", p_heat, width = 6, height = 5)
  
  write.csv(merged, "adar2_TPM_and_editing_merged.csv", row.names = FALSE)
  
  # Combine the two plots 
  combined_plot <- p | p_heat # this is using lib patchwork
  combined_plot
  
  # Save 
  ggsave("Combined_Correlation_Plot_ADAR2.png", combined_plot, width = 13, height = 5, dpi = 350)
  
  #=========================================== ADAR3==============================
  # ADAR3
  # Folders with ADAR TPM counts
  pre_male   <- "G:/Chapter4/ch4_ctlmale.100_subset/ch4_ctlmale.100_subset/ch4_ctlmale_100_215.counts-etc/counts/"
  pre_female <- "G:/Chapter4/ch4_ctlfemale.89_subset/ch4_ctlfemale.89_subset/ch4_ctlfemale_89.counts-etc/counts/"
  mid_male   <- "G:/Chapter4/ch4_midmale.100_subset/ch4_midmale.100_subset/ch4_midmale_100_345.counts-etc/counts/"
  mid_female <- "G:/Chapter4/ch4_midfemale.89_subset/ch4_midfemale.89_subset/ch4_midfemale_89.counts-etc/counts/"
  post_male  <- "G:/Chapter4/ch4_postmale.100_subset/ch4_postmale.100_subset/ch4_postmale_100_295.counts-etc/counts/"
  post_female<- "G:/Chapter4/ch4_postfemale.56_subset/ch4_postfemale.56_subset/ch4_postfemale_56.counts-etc/counts/"
  
  folder_list <- list(
    Pre_Male = pre_male,
    Pre_Female = pre_female,
    Mid_Male = mid_male,
    Mid_Female = mid_female,
    Post_Male = post_male,
    Post_Female = post_female
  )
  
  # Extract ADAR TPM
  
  extract_ADAR3 <- function(folder, group_name) {
    files <- list.files(folder, pattern = "\\.tab$", full.names = TRUE)
    
    map_df(files, function(f) {
      tab <- read.delim(f, stringsAsFactors = FALSE)
      
      adar_row <- tab %>% filter(Gene.Name %in% c("ADAR3", "ADARB2"))
      if (nrow(adar_row) == 0) return(NULL)
      
      tibble(
        sample = gsub(".tab$", "", basename(f)),
        ADAR3_TPM = adar_row$TPM,
        group = group_name
      )
    })
  }
  
  # ADAR exp
  adar_expression <- map2_df(folder_list, names(folder_list), extract_ADAR3)
  
  # .csv files with pre-caluclated overall ADAR editing values
  pre_female  <- "G:/Chapter4/Overall_ADAR_editing/Result_Pre_females/overall_editing_pre_infection_females.csv"
  pre_male    <- "G:/Chapter4/Overall_ADAR_editing/Result_pre_males/overall_editing_pre_infection_males.csv"
  
  mid_female  <- "G:/Chapter4/Overall_ADAR_editing/Result_Mid_females/overall_editing_Mid_infection_females.csv"
  mid_male    <- "G:/Chapter4/Overall_ADAR_editing/Result_Mid_males/overall_editing_Mid_infection_males.csv"
  
  post_female <- "G:/Chapter4/Overall_ADAR_editing/Result_Post_females/overall_editing_Post_infection_females.csv"
  post_male   <- "G:/Chapter4/Overall_ADAR_editing/Result_Post_males/overall_editing_Post_infection_males.csv"
  
  editing_list <- list(
    Pre_Female  = read.csv(pre_female),
    Pre_Male    = read.csv(pre_male),
    Mid_Female  = read.csv(mid_female),
    Mid_Male    = read.csv(mid_male),
    Post_Female = read.csv(post_female),
    Post_Male   = read.csv(post_male)
  )
  
  editing_df <- bind_rows(editing_list, .id = "group")
  
  ### Clean sample names
  editing_df$sample_clean <- editing_df$sample %>% 
    gsub("\\.csv$|\\.txt$", "", .) %>% 
    gsub("_.*$", "", .)
  
  # Merge files with expression and editing levels
  merged <- left_join(
    adar_expression,
    editing_df,
    by = c("sample" = "sample_clean"),
    suffix = c("_expr", "_edit")
  )
  
  # Modify/add colums
  merged <- merged %>%
    mutate(
      stage = case_when(
        grepl("Pre", group_expr) ~ "Pre-infection",
        grepl("Mid", group_expr) ~ "Mid-infection",
        grepl("Post", group_expr) ~ "Post-infection"
      ),
      sex = case_when(
        grepl("Male", group_expr) ~ "Male",
        grepl("Female", group_expr) ~ "Female"
      )
    )
  
  # Order stages
  merged$stage <- factor(merged$stage, 
                         levels = c("Pre-infection", "Mid-infection", "Post-infection"))
  
  # spearman correlations
  corr_table <- merged %>%
    group_by(stage, sex) %>%
    summarise(
      rho = cor(ADAR3_TPM, overall_editing_level, method = "spearman", use = "complete.obs"),
      p_value = cor.test(ADAR3_TPM, overall_editing_level, method = "spearman")$p.value,
      n = n()
    )
  
  write.csv(corr_table, "Spearman_correlation_results_ADAR3.csv", row.names = FALSE)
  
  # plot with spearman
  
  # p <- ggplot(merged, aes(x = ADAR3_TPM, y = overall_editing_level)) +
  #   geom_point(size = 3, alpha = 0.8) +
  #   geom_smooth(method = "lm", se = FALSE, color = "black") +
  #   # Remove the 'aes(label = ...)' line entirely.
  #   # The default label is usually sufficient and avoids parsing issues.
  #   stat_cor(
  #     method = "spearman",
  #     label.x.npc = "left",
  #     label.y.npc = "top",
  #     size = 4
  #   ) +
  #   facet_grid(stage ~ sex) +
  #   theme_bw(base_size = 14) +
  #   labs(
  #     x = "ADAR3 TPM",
  #     y = "Overall ADAR Editing Level",
  #     title = "Correlation between ADAR3 expression and ADAR editing"
  #   )
  # 
  # ggsave("ADAR3_correlation_facet_plot.png", p, width = 10, height = 8)
  # 
  
  
  # plot with spearman  ##r2 
  p <- ggplot(merged, aes(x = ADAR3_TPM, y = overall_editing_level)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    stat_cor(
      aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), # This will give R2
      method = "spearman", 
      label.x.npc = "left",
      label.y.npc = "top",
      size = 4
    ) +
    facet_grid(stage ~ sex) +
    theme_bw(base_size = 14) +
    labs(
      x = "ADAR3 TPM",
      y = "Overall ADAR Editing Level",
      title = "Correlation between ADAR3 expression and ADAR editing"
    )
  # Correlation heatmap
  heatmap_df <- corr_table %>%
    dcast(stage ~ sex, value.var = "p_value")
  
  rownames(heatmap_df) <- heatmap_df$stage
  heatmap_df$stage <- NULL
  
  p_heat <- ggplot(melt(as.matrix(heatmap_df)), 
                   aes(Var2, Var1, fill = value)) +
    geom_tile() +
    geom_text(aes(label = round(value, 2)), size = 6) +
    scale_fill_distiller(palette = "RdYlBu", direction = -1) +
    labs(
      x = "Sex",
      y = "Stage"
    ) +
    theme_minimal(base_size = 14)
  
  ggsave("Spearman_heatmap_ADAR3.png", p_heat, width = 6, height = 5)
  
write.csv(merged, "ADAR3_TPM_and_editing_merged.csv", row.names = FALSE)
  
  # Combine the two plots 
  combined_plot <- p | p_heat # this is using lib patchwork
  combined_plot
  
  # Save 
  ggsave("Combined_Correlation_Plot_ADAR3.png", combined_plot, width = 13, height = 5, dpi = 350)
  
  
  