# Differential expression analyses

###############################################
### Helper Functions
###############################################

## 1. Mean NPQ
summarise_mean_NPQ <- function(data, matrix_type,adjusted = FALSE) {
  if(adjusted){
    data %>%
      filter(SampleMatrixType == matrix_type, !is.na(NPQ)) %>%
      group_by(type, Target, UniProtID) %>%
      summarise(mean_NPQ = mean(NPQ_adj), .groups = "drop") %>%
      pivot_wider(
        names_from = type,
        values_from = mean_NPQ,
        names_prefix = "mean_NPQ_"
      )
  } else{
  data %>%
    filter(SampleMatrixType == matrix_type, !is.na(NPQ)) %>%
    group_by(type, Target, UniProtID) %>%
    summarise(mean_NPQ = mean(NPQ), .groups = "drop") %>%
    pivot_wider(
      names_from = type,
      values_from = mean_NPQ,
      names_prefix = "mean_NPQ_"
    )
    }
}

## 2. add log2FC for groups and PGMC subgroups
add_log2fc_standard <- function(df) {
  df %>%
    mutate(
      log2FC_ALS_CTR = mean_NPQ_ALS - mean_NPQ_CTR,
      log2FC_ALS_PGMC = mean_NPQ_ALS - mean_NPQ_PGMC,
      log2FC_ALS_mimic = mean_NPQ_ALS - mean_NPQ_mimic,
      log2FC_mimic_CTR = mean_NPQ_mimic - mean_NPQ_CTR,
      log2FC_PGMC_mimic = mean_NPQ_PGMC - mean_NPQ_mimic,
      log2FC_PGMC_CTR = mean_NPQ_PGMC - mean_NPQ_CTR
    )
}

add_log2fc_pgmc <- function(df) {
  df %>%
    mutate(
      log2FC_C9orf72_CTR = mean_NPQ_C9orf72 - mean_NPQ_CTR,
      log2FC_SOD1_CTR = mean_NPQ_SOD1 - mean_NPQ_CTR,
      log2FC_TARDBP_CTR = mean_NPQ_TARDBP - mean_NPQ_CTR,
      log2FC_C9orf72_SOD1 = mean_NPQ_C9orf72 - mean_NPQ_SOD1,
      log2FC_C9orf72_TARDBP = mean_NPQ_C9orf72 - mean_NPQ_TARDBP
    )
}

## 3. Final table
finalise_DE_table <- function(df,PGMC_groups = FALSE) {
  if(PGMC_groups){
    df %>%
      select(
        Target, UniProtID, Fluid,
        contains("C9orf72_CTR"),contains("CTR_C9orf72"),
        contains("TARDBP_CTR"),contains("CTR_TARDBP"),
        contains("SOD1_CTR"),contains("CTR_SOD1"),
        contains("C9orf72_SOD1"),contains("C9orf72_TARDBP")) %>%
      rename_with(~ gsub("CTR_C9orf72", "C9orf72_CTR", .x), .cols = contains("CTR_C9orf72")) %>%
      rename_with(~ gsub("CTR_TARDBP", "TARDBP_CTR", .x), .cols = contains("CTR_TARDBP")) %>%
      rename_with(~ gsub("CTR_SOD1", "SOD1_CTR", .x), .cols = contains("CTR_SOD1")) %>%
      rename_with(~ gsub("CTR_TARDBP", "TARDBP_CTR", .x), .cols = contains("CTR_TARDBP"))
  }else{
    df %>%
      select(
        Target, UniProtID, Fluid,
        contains("ALS_CTR"),contains("CTR_ALS"),
        contains("ALS_PGMC"),contains("PGMC_ALS"),
        contains("ALS_mimic"),contains("mimic_ALS"),
        contains("mimic_CTR"),contains("CTR_mimic"),
        contains("PGMC_CTR"),contains("CTR_PGMC"),
        contains("PGMC_mimic"),contains("mimic_PGMC")) %>%
      rename_with(~ gsub("CTR_ALS", "ALS_CTR", .x), .cols = contains("CTR_ALS")) %>%
      rename_with(~ gsub("PGMC_ALS", "ALS_PGMC", .x), .cols = contains("PGMC_ALS")) %>%
      rename_with(~ gsub("CTR_mimic", "mimic_CTR", .x), .cols = contains("CTR_mimic")) %>%
      rename_with(~ gsub("CTR_PGMC", "PGMC_CTR", .x), .cols = contains("CTR_PGMC")) %>%
      rename_with(~ gsub("mimic_PGMC", "PGMC_mimic", .x), .cols = contains("mimic_PGMC"))
  }
}

## 4. add DEx info in data
add_DE_flag <- function(df, lfc, padj, up, down, alpha) {
  df %>%
    mutate(
      log10_padj = -log10(.data[[padj]]),
      DE = case_when(
        .data[[padj]] < alpha & .data[[lfc]] > 0 ~ up,
        .data[[padj]] < alpha & .data[[lfc]] < 0 ~ down,
        TRUE ~ "ns"
      )
    )
}

## 5. Volcano plot
volcano_plot <- function(
    df, log2fc, title,
    fdr = 0.05,
    colors,
    label_targets = NULL,
    add_x = 0,
    add_y = 0
) {
  
  
  ggplot(df, aes_string(log2fc, "log10_padj")) +
    geom_point(aes(color = DE, size = log10_padj)) +
    geom_text_repel(
      data = df %>% filter(DE != "ns" | Target %in% label_targets),
      aes(label = Target),
      max.overlaps = 40,
      size = 5
    ) +
    geom_hline(yintercept = -log10(fdr), linetype = "dashed", color = "darkgrey") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey") +
    scale_color_manual(values = colors) +
    annotate(
      "text",
      x = add_x,
      y = -log10(fdr) + add_y,
      label = paste0("FDR = ", fdr * 100, "%"),
      color = "darkgrey",
      size = 5
    ) +
    theme_minimal() +
    ggtitle(title) +
    scale_size_continuous(range = c(2, 7)) +
    ylab(expression("-log"[10]*"(adjusted p-value)")) +
    xlab(expression("log"[2]*"(fold-change)")) + 
    labs(size = expression("-log"[10]*"(adjusted p-value)")) +
    theme(
      axis.title.x = element_blank(),
      text = element_text(size = 15),               # Base font size for everything
      axis.title = element_text(size = 18),         # Axis titles
      axis.text = element_text(size = 16),          # Axis tick labels
      plot.title = element_text(size = 18, hjust = 0.5,face = "bold"),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 15))
}


###############################################################################
# Run pipeline
###############################################################################

### PLASMA
## -> standard subgroups (not adjusted)
DE_subgroups_plasma = summarise_mean_NPQ(protein_data_IDs,matrix_type = "PLASMA") %>%
  add_log2fc_standard() %>%
  left_join(results_ALL$PLASMA$pvals) %>%
  mutate(Fluid = rep("PLASMA"))

DE_subgroups_plasma_final = finalise_DE_table(DE_subgroups_plasma)
writexl::write_xlsx(DE_subgroups_plasma_final,"results/DE_subgroups_plasma.xlsx")

## -> standard subgroups (adjusted)
protein_data_IDs_covariates = protein_data_IDs %>%
  mutate(center = dplyr::case_when(
    grepl("TR", ParticipantCode) ~ "Turkey",
    grepl("CH", ParticipantCode) ~ "Switzerland",
    grepl("DE", ParticipantCode) ~ "Germany",
    grepl("SK", ParticipantCode) ~ "Slovakia",
    grepl("FR", ParticipantCode) ~ "France",
    grepl("IL", ParticipantCode) ~ "Israel",
    TRUE                 ~ NA_character_)) %>%
  left_join(Sex_age_all_participants %>% dplyr::rename(PatientID = Pseudonyme))
DE_subgroups_plasma_adj = summarise_mean_NPQ(adjust_dataset(protein_data_IDs_covariates),
                                             matrix_type = "PLASMA",adjusted = TRUE) %>%
  add_log2fc_standard() %>%
  left_join(results_ALL$PLASMA$pvals_adjusted) %>%
  mutate(Fluid = rep("PLASMA"))

DE_subgroups_plasma_adj_final = finalise_DE_table(DE_subgroups_plasma_adj)
writexl::write_xlsx(DE_subgroups_plasma_adj_final,"results/DE_subgroups_plasma_adjusted.xlsx")

## -> PGMC subgroups (not adjusted)
DE_PGMC_subgroups_plasma = summarise_mean_NPQ(protein_data_CTR_PGMC_IDs,matrix_type = "PLASMA") %>%
  add_log2fc_pgmc() %>%
  left_join(results_PGMC$PLASMA$pvals) %>%
  mutate(Fluid = rep("PLASMA"))

DE_PGMC_subgroups_plasma_final = finalise_DE_table(DE_PGMC_subgroups_plasma,PGMC_groups = TRUE)
writexl::write_xlsx(DE_PGMC_subgroups_plasma_final,"results/DE_PGMC_subgroups_plasma.xlsx")

## -> PGMC subgroups ( adjusted)
protein_data_CTR_PGMC_IDs_covariates = protein_data_CTR_PGMC_IDs %>% select(-type) %>%
  left_join(samples_PGMC_CTR_ID_type %>% dplyr::rename(SampleName = `Sample ID`)) %>%
  mutate(center = dplyr::case_when(
    grepl("TR", ParticipantCode) ~ "Turkey",
    grepl("CH", ParticipantCode) ~ "Switzerland",
    grepl("DE", ParticipantCode) ~ "Germany",
    grepl("SK", ParticipantCode) ~ "Slovakia",
    grepl("FR", ParticipantCode) ~ "France",
    grepl("IL", ParticipantCode) ~ "Israel",
    TRUE                 ~ NA_character_)) %>%
  left_join(Sex_age_all_participants %>% dplyr::rename(PatientID = Pseudonyme))
DE_PGMC_subgroups_plasma_adj = summarise_mean_NPQ(adjust_dataset(protein_data_CTR_PGMC_IDs_covariates),
                                              matrix_type = "PLASMA",adjusted = TRUE) %>%
  add_log2fc_pgmc() %>%
  left_join(results_PGMC$PLASMA$pvals_adjusted) %>%
  mutate(Fluid = rep("PLASMA"))

DE_PGMC_subgroups_plasma_adj_final = finalise_DE_table(DE_PGMC_subgroups_plasma_adj,PGMC_groups = TRUE)
writexl::write_xlsx(DE_PGMC_subgroups_plasma_adj_final,"results/DE_PGMC_subgroups_plasma_adjusted.xlsx")

### SERUM
## -> standard subgroups (not adjusted)
DE_subgroups_SERUM = summarise_mean_NPQ(protein_data_IDs,matrix_type = "SERUM") %>%
  add_log2fc_standard() %>%
  left_join(results_ALL$SERUM$pvals) %>%
  mutate(Fluid = rep("SERUM"))

DE_subgroups_SERUM_final = finalise_DE_table(DE_subgroups_SERUM)
writexl::write_xlsx(DE_subgroups_SERUM_final,"results/DE_subgroups_SERUM.xlsx")

## -> standard subgroups (adjusted)
DE_subgroups_SERUM_adj = summarise_mean_NPQ(adjust_dataset(protein_data_IDs_covariates),
                                             matrix_type = "SERUM",adjusted = TRUE) %>%
  add_log2fc_standard() %>%
  left_join(results_ALL$SERUM$pvals_adjusted) %>%
  mutate(Fluid = rep("SERUM"))

DE_subgroups_SERUM_adj_final = finalise_DE_table(DE_subgroups_SERUM_adj)
writexl::write_xlsx(DE_subgroups_SERUM_adj_final,"results/DE_subgroups_SERUM_adjusted.xlsx")

## -> PGMC subgroups (not adjusted)
DE_PGMC_subgroups_SERUM = summarise_mean_NPQ(protein_data_CTR_PGMC_IDs,matrix_type = "SERUM") %>%
  add_log2fc_pgmc() %>%
  left_join(results_PGMC$SERUM$pvals) %>%
  mutate(Fluid = rep("SERUM"))

DE_PGMC_subgroups_SERUM_final = finalise_DE_table(DE_PGMC_subgroups_SERUM,PGMC_groups = TRUE)
writexl::write_xlsx(DE_PGMC_subgroups_SERUM_final,"results/DE_PGMC_subgroups_SERUM.xlsx")

## -> PGMC subgroups ( adjusted)
DE_PGMC_subgroups_SERUM_adj = summarise_mean_NPQ(adjust_dataset(protein_data_CTR_PGMC_IDs_covariates),
                                                  matrix_type = "SERUM",adjusted = TRUE) %>%
  add_log2fc_pgmc() %>%
  left_join(results_PGMC$SERUM$pvals_adjusted) %>%
  mutate(Fluid = rep("SERUM"))

DE_PGMC_subgroups_SERUM_adj_final = finalise_DE_table(DE_PGMC_subgroups_SERUM_adj,PGMC_groups = TRUE)
writexl::write_xlsx(DE_PGMC_subgroups_SERUM_adj_final,"results/DE_PGMC_subgroups_SERUM_adjusted.xlsx")

### CSF
## -> standard subgroups (not adjusted)
DE_subgroups_CSF = summarise_mean_NPQ(protein_data_IDs,matrix_type = "CSF") %>%
  add_log2fc_standard() %>%
  left_join(results_ALL$CSF$pvals) %>%
  mutate(Fluid = rep("CSF"))

DE_subgroups_CSF_final = finalise_DE_table(DE_subgroups_CSF)
writexl::write_xlsx(DE_subgroups_CSF_final,"results/DE_subgroups_CSF.xlsx")

## -> standard subgroups (adjusted)
DE_subgroups_CSF_adj = summarise_mean_NPQ(adjust_dataset(protein_data_IDs_covariates),
                                            matrix_type = "CSF",adjusted = TRUE) %>%
  add_log2fc_standard() %>%
  left_join(results_ALL$CSF$pvals_adjusted) %>%
  mutate(Fluid = rep("CSF"))

DE_subgroups_CSF_adj_final = finalise_DE_table(DE_subgroups_CSF_adj)
writexl::write_xlsx(DE_subgroups_CSF_adj_final,"results/DE_subgroups_CSF_adjusted.xlsx")

## -> PGMC subgroups (not adjusted)
DE_PGMC_subgroups_CSF = summarise_mean_NPQ(protein_data_CTR_PGMC_IDs,matrix_type = "CSF") %>%
  add_log2fc_pgmc() %>%
  left_join(results_PGMC$CSF$pvals) %>%
  mutate(Fluid = rep("CSF"))

DE_PGMC_subgroups_CSF_final = finalise_DE_table(DE_PGMC_subgroups_CSF,PGMC_groups = TRUE)
writexl::write_xlsx(DE_PGMC_subgroups_CSF_final,"results/DE_PGMC_subgroups_CSF.xlsx")

## -> PGMC subgroups ( adjusted)
DE_PGMC_subgroups_CSF_adj = summarise_mean_NPQ(adjust_dataset(protein_data_CTR_PGMC_IDs_covariates),
                                                 matrix_type = "CSF",adjusted = TRUE) %>%
  add_log2fc_pgmc() %>%
  left_join(results_PGMC$CSF$pvals_adjusted) %>%
  mutate(Fluid = rep("CSF"))

DE_PGMC_subgroups_CSF_adj_final = finalise_DE_table(DE_PGMC_subgroups_CSF_adj,PGMC_groups = TRUE)
writexl::write_xlsx(DE_PGMC_subgroups_CSF_adj_final,"results/DE_PGMC_subgroups_CSF_adjusted.xlsx")

###### VOLCANO PLOTS 
## 1: ALS vs CTR 
# -> plasma (not adjusted)
DE_plasma_ALS_CTR = DE_subgroups_plasma_final %>%
  select(Target,UniProtID,log2FC_ALS_CTR,pvalue_ALS_CTR,padj_ALS_CTR) %>%
  add_DE_flag(lfc = "log2FC_ALS_CTR",padj = "padj_ALS_CTR",up = "ALS",down = "CTR",alpha = 0.05)

volcano_ALS_CTR_plasma = volcano_plot(DE_plasma_ALS_CTR,"log2FC_ALS_CTR", 
                                      title = "ALS vs CTR in plasma",
                                      colors = c('CTR' = "#C25F3D", 
                                                 'ns' = 'lightgrey', 
                                                 'ALS'= "#3da0c2"),
                                      add_x = 4.9,add_y = 0.6)

pdf("plots/volcano_plots/volcano_ALS_CTR_plasma.pdf",width = 8, height = 6)
volcano_ALS_CTR_plasma
dev.off()

# -> plasma (adjusted)
DE_plasma_ALS_CTR_adj = DE_subgroups_plasma_adj_final %>%
  select(Target,UniProtID,log2FC_ALS_CTR,pvalue_ALS_CTR,padj_ALS_CTR) %>%
  add_DE_flag(lfc = "log2FC_ALS_CTR",padj = "padj_ALS_CTR",up = "ALS",down = "CTR",alpha = 0.05)

volcano_ALS_CTR_plasma_adj = volcano_plot(DE_plasma_ALS_CTR_adj,"log2FC_ALS_CTR", 
                                      title = "ALS vs CTR in plasma (adjusted)",
                                      colors = c('CTR' = "#C25F3D", 
                                                 'ns' = 'lightgrey', 
                                                 'ALS'= "#3da0c2"),
                                      add_x = 2.5,add_y = 0.2)

pdf("plots/volcano_plots/volcano_ALS_CTR_plasma_adj.pdf",width = 8, height = 6)
volcano_ALS_CTR_plasma_adj
dev.off()

# -> CSF (not adjusted)
DE_CSF_ALS_CTR = DE_subgroups_CSF_final %>%
  select(Target,UniProtID,log2FC_ALS_CTR,pvalue_ALS_CTR,padj_ALS_CTR) %>%
  add_DE_flag(lfc = "log2FC_ALS_CTR",padj = "padj_ALS_CTR",up = "ALS",down = "CTR",alpha = 0.05)

volcano_ALS_CTR_CSF = volcano_plot(DE_CSF_ALS_CTR,"log2FC_ALS_CTR", 
                                      title = "ALS vs CTR in CSF",
                                      colors = c('CTR' = "#C25F3D", 
                                                 'ns' = 'lightgrey', 
                                                 'ALS'= "#3da0c2"),
                                      add_x = 3.4,add_y = 0.6)

pdf("plots/volcano_plots/volcano_ALS_CTR_CSF.pdf",width = 8, height = 6)
volcano_ALS_CTR_CSF
dev.off()

# -> CSF (adjusted)
DE_CSF_ALS_CTR_adj = DE_subgroups_CSF_adj_final %>%
  select(Target,UniProtID,log2FC_ALS_CTR,pvalue_ALS_CTR,padj_ALS_CTR) %>%
  add_DE_flag(lfc = "log2FC_ALS_CTR",padj = "padj_ALS_CTR",up = "ALS",down = "CTR",alpha = 0.05)

volcano_ALS_CTR_CSF_adj = volcano_plot(DE_CSF_ALS_CTR_adj,"log2FC_ALS_CTR", 
                                          title = "ALS vs CTR in CSF (adjusted)",
                                          colors = c('CTR' = "#C25F3D", 
                                                     'ns' = 'lightgrey', 
                                                     'ALS'= "#3da0c2"),
                                          add_x = 1.8,add_y = 0.2)

pdf("plots/volcano_plots/volcano_ALS_CTR_CSF_adj.pdf",width = 8, height = 6)
volcano_ALS_CTR_CSF_adj
dev.off()

# -> SERUM (not adjusted)
DE_SERUM_ALS_CTR = DE_subgroups_SERUM_final %>%
  select(Target,UniProtID,log2FC_ALS_CTR,pvalue_ALS_CTR,padj_ALS_CTR) %>%
  add_DE_flag(lfc = "log2FC_ALS_CTR",padj = "padj_ALS_CTR",up = "ALS",down = "CTR",alpha = 0.05)

volcano_ALS_CTR_SERUM = volcano_plot(DE_SERUM_ALS_CTR,"log2FC_ALS_CTR", 
                                   title = "ALS vs CTR in SERUM",
                                   colors = c('CTR' = "#C25F3D", 
                                              'ns' = 'lightgrey', 
                                              'ALS'= "#3da0c2"),
                                   add_x = 4.3,add_y = 0.6)

pdf("plots/volcano_plots/volcano_ALS_CTR_SERUM.pdf",width = 8, height = 6)
volcano_ALS_CTR_SERUM
dev.off()

# -> SERUM (adjusted)
DE_SERUM_ALS_CTR_adj = DE_subgroups_SERUM_adj_final %>%
  select(Target,UniProtID,log2FC_ALS_CTR,pvalue_ALS_CTR,padj_ALS_CTR) %>%
  add_DE_flag(lfc = "log2FC_ALS_CTR",padj = "padj_ALS_CTR",up = "ALS",down = "CTR",alpha = 0.05)

volcano_ALS_CTR_SERUM_adj = volcano_plot(DE_SERUM_ALS_CTR_adj,"log2FC_ALS_CTR", 
                                       title = "ALS vs CTR in SERUM (adjusted)",
                                       colors = c('CTR' = "#C25F3D", 
                                                  'ns' = 'lightgrey', 
                                                  'ALS'= "#3da0c2"),
                                       add_x = 2.2,add_y = 0.2)

pdf("plots/volcano_plots/volcano_ALS_CTR_SERUM_adj.pdf",width = 8, height = 6)
volcano_ALS_CTR_SERUM_adj
dev.off()

# -> all together
pdf("plots/volcano_plots/volcano_ALS_CTR_all.pdf", onefile = TRUE, width = 14, height = 15)
grid.arrange(volcano_ALS_CTR_plasma, 
             volcano_ALS_CTR_plasma_adj,
             volcano_ALS_CTR_SERUM,
             volcano_ALS_CTR_SERUM_adj,
             volcano_ALS_CTR_CSF,
             volcano_ALS_CTR_CSF_adj,nrow = 3, ncol = 2)
dev.off()

## 2: ALS vs PGMC 
# -> plasma (not adjusted)
DE_plasma_ALS_PGMC = DE_subgroups_plasma_final %>%
  select(Target,UniProtID,log2FC_ALS_PGMC,pvalue_ALS_PGMC,padj_ALS_PGMC) %>%
  add_DE_flag(lfc = "log2FC_ALS_PGMC",padj = "padj_ALS_PGMC",up = "ALS",down = "PGMC",alpha = 0.05)

volcano_ALS_PGMC_plasma = volcano_plot(DE_plasma_ALS_PGMC,"log2FC_ALS_PGMC", 
                                      title = "ALS vs PGMC in plasma",
                                      colors = c('PGMC' = "#C25F3D", 
                                                 'ns' = 'lightgrey', 
                                                 'ALS'= "#3da0c2"),
                                      add_x = 3.3,add_y = 0.6)

pdf("plots/volcano_plots/volcano_ALS_PGMC_plasma.pdf",width = 8, height = 6)
volcano_ALS_PGMC_plasma
dev.off()

# -> plasma (adjusted)
DE_plasma_ALS_PGMC_adj = DE_subgroups_plasma_adj_final %>%
  select(Target,UniProtID,log2FC_ALS_PGMC,pvalue_ALS_PGMC,padj_ALS_PGMC) %>%
  add_DE_flag(lfc = "log2FC_ALS_PGMC",padj = "padj_ALS_PGMC",up = "ALS",down = "PGMC",alpha = 0.05)

volcano_ALS_PGMC_plasma_adj = volcano_plot(DE_plasma_ALS_PGMC_adj,"log2FC_ALS_PGMC", 
                                          title = "ALS vs PGMC in plasma (adjusted)",
                                          colors = c('PGMC' = "#C25F3D", 
                                                     'ns' = 'lightgrey', 
                                                     'ALS'= "#3da0c2"),
                                          add_x = 1,add_y = 0.2)

pdf("plots/volcano_plots/volcano_ALS_PGMC_plasma_adj.pdf",width = 8, height = 6)
volcano_ALS_PGMC_plasma_adj
dev.off()

# -> CSF (not adjusted)
DE_CSF_ALS_PGMC = DE_subgroups_CSF_final %>%
  select(Target,UniProtID,log2FC_ALS_PGMC,pvalue_ALS_PGMC,padj_ALS_PGMC) %>%
  add_DE_flag(lfc = "log2FC_ALS_PGMC",padj = "padj_ALS_PGMC",up = "ALS",down = "PGMC",alpha = 0.05)

volcano_ALS_PGMC_CSF = volcano_plot(DE_CSF_ALS_PGMC,"log2FC_ALS_PGMC", 
                                   title = "ALS vs PGMC in CSF",
                                   colors = c('PGMC' = "#C25F3D", 
                                              'ns' = 'lightgrey', 
                                              'ALS'= "#3da0c2"),
                                   add_x = 3.4,add_y = 0.6)

pdf("plots/volcano_plots/volcano_ALS_PGMC_CSF.pdf",width = 8, height = 6)
volcano_ALS_PGMC_CSF
dev.off()

# -> CSF (adjusted)
DE_CSF_ALS_PGMC_adj = DE_subgroups_CSF_adj_final %>%
  select(Target,UniProtID,log2FC_ALS_PGMC,pvalue_ALS_PGMC,padj_ALS_PGMC) %>%
  add_DE_flag(lfc = "log2FC_ALS_PGMC",padj = "padj_ALS_PGMC",up = "ALS",down = "PGMC",alpha = 0.05)

volcano_ALS_PGMC_CSF_adj = volcano_plot(DE_CSF_ALS_PGMC_adj,"log2FC_ALS_PGMC", 
                                       title = "ALS vs PGMC in CSF (adjusted)",
                                       colors = c('PGMC' = "#C25F3D", 
                                                  'ns' = 'lightgrey', 
                                                  'ALS'= "#3da0c2"),
                                       add_x = 0.6,add_y = 0.1)

pdf("plots/volcano_plots/volcano_ALS_PGMC_CSF_adj.pdf",width = 8, height = 6)
volcano_ALS_PGMC_CSF_adj
dev.off()

# -> SERUM (not adjusted)
DE_SERUM_ALS_PGMC = DE_subgroups_SERUM_final %>%
  select(Target,UniProtID,log2FC_ALS_PGMC,pvalue_ALS_PGMC,padj_ALS_PGMC) %>%
  add_DE_flag(lfc = "log2FC_ALS_PGMC",padj = "padj_ALS_PGMC",up = "ALS",down = "PGMC",alpha = 0.05)

volcano_ALS_PGMC_SERUM = volcano_plot(DE_SERUM_ALS_PGMC,"log2FC_ALS_PGMC", 
                                     title = "ALS vs PGMC in SERUM",
                                     colors = c('PGMC' = "#C25F3D", 
                                                'ns' = 'lightgrey', 
                                                'ALS'= "#3da0c2"),
                                     add_x = 3.5,add_y = 0.6)

pdf("plots/volcano_plots/volcano_ALS_PGMC_SERUM.pdf",width = 8, height = 6)
volcano_ALS_PGMC_SERUM
dev.off()

# -> SERUM (adjusted)
DE_SERUM_ALS_PGMC_adj = DE_subgroups_SERUM_adj_final %>%
  select(Target,UniProtID,log2FC_ALS_PGMC,pvalue_ALS_PGMC,padj_ALS_PGMC) %>%
  add_DE_flag(lfc = "log2FC_ALS_PGMC",padj = "padj_ALS_PGMC",up = "ALS",down = "PGMC",alpha = 0.05)

volcano_ALS_PGMC_SERUM_adj = volcano_plot(DE_SERUM_ALS_PGMC_adj,"log2FC_ALS_PGMC", 
                                         title = "ALS vs PGMC in SERUM (adjusted)",
                                         colors = c('PGMC' = "#C25F3D", 
                                                    'ns' = 'lightgrey', 
                                                    'ALS'= "#3da0c2"),
                                         add_x = 1.1,add_y = 0.1)

pdf("plots/volcano_plots/volcano_ALS_PGMC_SERUM_adj.pdf",width = 8, height = 6)
volcano_ALS_PGMC_SERUM_adj
dev.off()

# -> all together
pdf("plots/volcano_plots/volcano_ALS_PGMC_all.pdf", onefile = TRUE, width = 14, height = 15)
grid.arrange(volcano_ALS_PGMC_plasma, 
             volcano_ALS_PGMC_plasma_adj,
             volcano_ALS_PGMC_SERUM,
             volcano_ALS_PGMC_SERUM_adj,
             volcano_ALS_PGMC_CSF,
             volcano_ALS_PGMC_CSF_adj,nrow = 3, ncol = 2)
dev.off()

## 3: ALS vs mimic 
# -> plasma (not adjusted)
DE_plasma_ALS_mimic = DE_subgroups_plasma_final %>%
  select(Target,UniProtID,log2FC_ALS_mimic,pvalue_ALS_mimic,padj_ALS_mimic) %>%
  add_DE_flag(lfc = "log2FC_ALS_mimic",padj = "padj_ALS_mimic",up = "ALS",down = "mimic",alpha = 0.05)

volcano_ALS_mimic_plasma = volcano_plot(DE_plasma_ALS_mimic,"log2FC_ALS_mimic", 
                                       title = "ALS vs mimic in plasma",
                                       colors = c('mimic' = "#C25F3D", 
                                                  'ns' = 'lightgrey', 
                                                  'ALS'= "#3da0c2"),
                                       add_x = 2.8,add_y = 0.4)

pdf("plots/volcano_plots/volcano_ALS_mimic_plasma.pdf",width = 8, height = 6)
volcano_ALS_mimic_plasma
dev.off()

# -> plasma (adjusted)
DE_plasma_ALS_mimic_adj = DE_subgroups_plasma_adj_final %>%
  select(Target,UniProtID,log2FC_ALS_mimic,pvalue_ALS_mimic,padj_ALS_mimic) %>%
  add_DE_flag(lfc = "log2FC_ALS_mimic",padj = "padj_ALS_mimic",up = "ALS",down = "mimic",alpha = 0.05)

volcano_ALS_mimic_plasma_adj = volcano_plot(DE_plasma_ALS_mimic_adj,"log2FC_ALS_mimic", 
                                           title = "ALS vs mimic in plasma (adjusted)",
                                           colors = c('mimic' = "#C25F3D", 
                                                      'ns' = 'lightgrey', 
                                                      'ALS'= "#3da0c2"),
                                           add_x = 2,add_y = 0.2)

pdf("plots/volcano_plots/volcano_ALS_mimic_plasma_adj.pdf",width = 8, height = 6)
volcano_ALS_mimic_plasma_adj
dev.off()

# -> CSF (not adjusted)
DE_CSF_ALS_mimic = DE_subgroups_CSF_final %>%
  select(Target,UniProtID,log2FC_ALS_mimic,pvalue_ALS_mimic,padj_ALS_mimic) %>%
  add_DE_flag(lfc = "log2FC_ALS_mimic",padj = "padj_ALS_mimic",up = "ALS",down = "mimic",alpha = 0.05)

volcano_ALS_mimic_CSF = volcano_plot(DE_CSF_ALS_mimic,"log2FC_ALS_mimic", 
                                    title = "ALS vs mimic in CSF",
                                    colors = c('mimic' = "#C25F3D", 
                                               'ns' = 'lightgrey', 
                                               'ALS'= "#3da0c2"),
                                    add_x = 1.5,add_y = 0.1)

pdf("plots/volcano_plots/volcano_ALS_mimic_CSF.pdf",width = 8, height = 6)
volcano_ALS_mimic_CSF
dev.off()

# -> CSF (adjusted)
DE_CSF_ALS_mimic_adj = DE_subgroups_CSF_adj_final %>%
  select(Target,UniProtID,log2FC_ALS_mimic,pvalue_ALS_mimic,padj_ALS_mimic) %>%
  add_DE_flag(lfc = "log2FC_ALS_mimic",padj = "padj_ALS_mimic",up = "ALS",down = "mimic",alpha = 0.05)

volcano_ALS_mimic_CSF_adj = volcano_plot(DE_CSF_ALS_mimic_adj,"log2FC_ALS_mimic", 
                                        title = "ALS vs mimic in CSF (adjusted)",
                                        colors = c('mimic' = "#C25F3D", 
                                                   'ns' = 'lightgrey', 
                                                   'ALS'= "#3da0c2"),
                                        add_x = 1.1,add_y = 0.1)

pdf("plots/volcano_plots/volcano_ALS_mimic_CSF_adj.pdf",width = 8, height = 6)
volcano_ALS_mimic_CSF_adj
dev.off()

# -> SERUM (not adjusted)
DE_SERUM_ALS_mimic = DE_subgroups_SERUM_final %>%
  select(Target,UniProtID,log2FC_ALS_mimic,pvalue_ALS_mimic,padj_ALS_mimic) %>%
  add_DE_flag(lfc = "log2FC_ALS_mimic",padj = "padj_ALS_mimic",up = "ALS",down = "mimic",alpha = 0.05)

volcano_ALS_mimic_SERUM = volcano_plot(DE_SERUM_ALS_mimic,"log2FC_ALS_mimic", 
                                      title = "ALS vs mimic in SERUM",
                                      colors = c('mimic' = "#C25F3D", 
                                                 'ns' = 'lightgrey', 
                                                 'ALS'= "#3da0c2"),
                                      add_x = 3,add_y = 0.3)

pdf("plots/volcano_plots/volcano_ALS_mimic_SERUM.pdf",width = 8, height = 6)
volcano_ALS_mimic_SERUM
dev.off()

# -> SERUM (adjusted)
DE_SERUM_ALS_mimic_adj = DE_subgroups_SERUM_adj_final %>%
  select(Target,UniProtID,log2FC_ALS_mimic,pvalue_ALS_mimic,padj_ALS_mimic) %>%
  add_DE_flag(lfc = "log2FC_ALS_mimic",padj = "padj_ALS_mimic",up = "ALS",down = "mimic",alpha = 0.05)

volcano_ALS_mimic_SERUM_adj = volcano_plot(DE_SERUM_ALS_mimic_adj,"log2FC_ALS_mimic", 
                                          title = "ALS vs mimic in SERUM (adjusted)",
                                          colors = c('mimic' = "#C25F3D", 
                                                     'ns' = 'lightgrey', 
                                                     'ALS'= "#3da0c2"),
                                          add_x = 2.1,add_y = 0.1)

pdf("plots/volcano_plots/volcano_ALS_mimic_SERUM_adj.pdf",width = 8, height = 6)
volcano_ALS_mimic_SERUM_adj
dev.off()

# -> all together
pdf("plots/volcano_plots/volcano_ALS_mimic_all.pdf", onefile = TRUE, width = 14, height = 15)
grid.arrange(volcano_ALS_mimic_plasma, 
             volcano_ALS_mimic_plasma_adj,
             volcano_ALS_mimic_SERUM,
             volcano_ALS_mimic_SERUM_adj,
             volcano_ALS_mimic_CSF,
             volcano_ALS_mimic_CSF_adj,nrow = 3, ncol = 2)
dev.off()

## 4: mimic vs CTR 
# -> plasma (not adjusted)
DE_plasma_mimic_CTR = DE_subgroups_plasma_final %>%
  select(Target,UniProtID,log2FC_mimic_CTR,pvalue_mimic_CTR,padj_mimic_CTR) %>%
  add_DE_flag(lfc = "log2FC_mimic_CTR",padj = "padj_mimic_CTR",up = "mimic",down = "CTR",alpha = 0.05)

volcano_mimic_CTR_plasma = volcano_plot(DE_plasma_mimic_CTR,"log2FC_mimic_CTR", 
                                      title = "mimic vs CTR in plasma",
                                      colors = c('CTR' = "#C25F3D", 
                                                 'ns' = 'lightgrey', 
                                                 'mimic'= "#3da0c2"),
                                      add_x = 2.3,add_y = 0.1)

pdf("plots/volcano_plots/volcano_mimic_CTR_plasma.pdf",width = 8, height = 6)
volcano_mimic_CTR_plasma
dev.off()

# -> plasma (adjusted)
DE_plasma_mimic_CTR_adj = DE_subgroups_plasma_adj_final %>%
  select(Target,UniProtID,log2FC_mimic_CTR,pvalue_mimic_CTR,padj_mimic_CTR) %>%
  add_DE_flag(lfc = "log2FC_mimic_CTR",padj = "padj_mimic_CTR",up = "mimic",down = "CTR",alpha = 0.05)

volcano_mimic_CTR_plasma_adj = volcano_plot(DE_plasma_mimic_CTR_adj,"log2FC_mimic_CTR", 
                                          title = "mimic vs CTR in plasma (adjusted)",
                                          colors = c('CTR' = "#C25F3D", 
                                                     'ns' = 'lightgrey', 
                                                     'mimic'= "#3da0c2"),
                                          add_x = 1.1,add_y = 0.1)

pdf("plots/volcano_plots/volcano_mimic_CTR_plasma_adj.pdf",width = 8, height = 6)
volcano_mimic_CTR_plasma_adj
dev.off()

# -> CSF (not adjusted)
DE_CSF_mimic_CTR = DE_subgroups_CSF_final %>%
  select(Target,UniProtID,log2FC_mimic_CTR,pvalue_mimic_CTR,padj_mimic_CTR) %>%
  add_DE_flag(lfc = "log2FC_mimic_CTR",padj = "padj_mimic_CTR",up = "mimic",down = "CTR",alpha = 0.05)

volcano_mimic_CTR_CSF = volcano_plot(DE_CSF_mimic_CTR,"log2FC_mimic_CTR", 
                                   title = "mimic vs CTR in CSF",
                                   colors = c('CTR' = "#C25F3D", 
                                              'ns' = 'lightgrey', 
                                              'mimic'= "#3da0c2"),
                                   add_x = 3.4,add_y = 0.1)

pdf("plots/volcano_plots/volcano_mimic_CTR_CSF.pdf",width = 8, height = 6)
volcano_mimic_CTR_CSF
dev.off()

# -> CSF (adjusted)
DE_CSF_mimic_CTR_adj = DE_subgroups_CSF_adj_final %>%
  select(Target,UniProtID,log2FC_mimic_CTR,pvalue_mimic_CTR,padj_mimic_CTR) %>%
  add_DE_flag(lfc = "log2FC_mimic_CTR",padj = "padj_mimic_CTR",up = "mimic",down = "CTR",alpha = 0.05)

volcano_mimic_CTR_CSF_adj = volcano_plot(DE_CSF_mimic_CTR_adj,"log2FC_mimic_CTR", 
                                       title = "mimic vs CTR in CSF (adjusted)",
                                       colors = c('CTR' = "#C25F3D", 
                                                  'ns' = 'lightgrey', 
                                                  'mimic'= "#3da0c2"),
                                       add_x = 2.8,add_y = 0.1)

pdf("plots/volcano_plots/volcano_mimic_CTR_CSF_adj.pdf",width = 8, height = 6)
volcano_mimic_CTR_CSF_adj
dev.off()

# -> SERUM (not adjusted)
DE_SERUM_mimic_CTR = DE_subgroups_SERUM_final %>%
  select(Target,UniProtID,log2FC_mimic_CTR,pvalue_mimic_CTR,padj_mimic_CTR) %>%
  add_DE_flag(lfc = "log2FC_mimic_CTR",padj = "padj_mimic_CTR",up = "mimic",down = "CTR",alpha = 0.05)

volcano_mimic_CTR_SERUM = volcano_plot(DE_SERUM_mimic_CTR,"log2FC_mimic_CTR", 
                                     title = "mimic vs CTR in SERUM",
                                     colors = c('CTR' = "#C25F3D", 
                                                'ns' = 'lightgrey', 
                                                'mimic'= "#3da0c2"),
                                     add_x = 2,add_y = 0.2)

pdf("plots/volcano_plots/volcano_mimic_CTR_SERUM.pdf",width = 8, height = 6)
volcano_mimic_CTR_SERUM
dev.off()

# -> SERUM (adjusted)
DE_SERUM_mimic_CTR_adj = DE_subgroups_SERUM_adj_final %>%
  select(Target,UniProtID,log2FC_mimic_CTR,pvalue_mimic_CTR,padj_mimic_CTR) %>%
  add_DE_flag(lfc = "log2FC_mimic_CTR",padj = "padj_mimic_CTR",up = "mimic",down = "CTR",alpha = 0.05)

volcano_mimic_CTR_SERUM_adj = volcano_plot(DE_SERUM_mimic_CTR_adj,"log2FC_mimic_CTR", 
                                         title = "mimic vs CTR in SERUM (adjusted)",
                                         colors = c('CTR' = "#C25F3D", 
                                                    'ns' = 'lightgrey', 
                                                    'mimic'= "#3da0c2"),
                                         add_x = 1.2,add_y = 0.1)

pdf("plots/volcano_plots/volcano_mimic_CTR_SERUM_adj.pdf",width = 8, height = 6)
volcano_mimic_CTR_SERUM_adj
dev.off()

# -> all together
pdf("plots/volcano_plots/volcano_mimic_CTR_all.pdf", onefile = TRUE, width = 14, height = 15)
grid.arrange(volcano_mimic_CTR_plasma, 
             volcano_mimic_CTR_plasma_adj,
             volcano_mimic_CTR_SERUM,
             volcano_mimic_CTR_SERUM_adj,
             volcano_mimic_CTR_CSF,
             volcano_mimic_CTR_CSF_adj,nrow = 3, ncol = 2)
dev.off()

## 5: PGMC vs CTR 
# -> plasma (not adjusted)
DE_plasma_PGMC_CTR = DE_subgroups_plasma_final %>%
  select(Target,UniProtID,log2FC_PGMC_CTR,pvalue_PGMC_CTR,padj_PGMC_CTR) %>%
  add_DE_flag(lfc = "log2FC_PGMC_CTR",padj = "padj_PGMC_CTR",up = "PGMC",down = "CTR",alpha = 0.05)

volcano_PGMC_CTR_plasma = volcano_plot(DE_plasma_PGMC_CTR,"log2FC_PGMC_CTR", 
                                        title = "PGMC vs CTR in plasma",
                                        colors = c('CTR' = "#C25F3D", 
                                                   'ns' = 'lightgrey', 
                                                   'PGMC'= "#3da0c2"),
                                        add_x = 1.2,add_y = 0.1,
                                       label_targets = c("MAPT","NEFH","pTau-217","TAFA5",
                                                         "IL7","SOD1","HTT"))

pdf("plots/volcano_plots/volcano_PGMC_CTR_plasma.pdf",width = 8, height = 6)
volcano_PGMC_CTR_plasma
dev.off()

# -> plasma (adjusted)
DE_plasma_PGMC_CTR_adj = DE_subgroups_plasma_adj_final %>%
  select(Target,UniProtID,log2FC_PGMC_CTR,pvalue_PGMC_CTR,padj_PGMC_CTR) %>%
  add_DE_flag(lfc = "log2FC_PGMC_CTR",padj = "padj_PGMC_CTR",up = "PGMC",down = "CTR",alpha = 0.05)

volcano_PGMC_CTR_plasma_adj = volcano_plot(DE_plasma_PGMC_CTR_adj,"log2FC_PGMC_CTR", 
                                            title = "PGMC vs CTR in plasma (adjusted)",
                                            colors = c('CTR' = "#C25F3D", 
                                                       'ns' = 'lightgrey', 
                                                       'PGMC'= "#3da0c2"),
                                            add_x = 1.5,add_y = 0.1,
                                           label_targets = c("MAPT","NEFH","pTau-217","TAFA5",
                                                             "IL7","SOD1","HTT"))

pdf("plots/volcano_plots/volcano_PGMC_CTR_plasma_adj.pdf",width = 8, height = 6)
volcano_PGMC_CTR_plasma_adj
dev.off()

# -> CSF (not adjusted)
DE_CSF_PGMC_CTR = DE_subgroups_CSF_final %>%
  select(Target,UniProtID,log2FC_PGMC_CTR,pvalue_PGMC_CTR,padj_PGMC_CTR) %>%
  add_DE_flag(lfc = "log2FC_PGMC_CTR",padj = "padj_PGMC_CTR",up = "PGMC",down = "CTR",alpha = 0.05)

volcano_PGMC_CTR_CSF = volcano_plot(DE_CSF_PGMC_CTR,"log2FC_PGMC_CTR", 
                                     title = "PGMC vs CTR in CSF",
                                     colors = c('CTR' = "#C25F3D", 
                                                'ns' = 'lightgrey', 
                                                'PGMC'= "#3da0c2"),
                                     add_x = 1,add_y = 0.1,
                                    label_targets = c("MAPT","NEFH","pTau-217","TAFA5",
                                                      "IL7","SOD1","HTT"))

pdf("plots/volcano_plots/volcano_PGMC_CTR_CSF.pdf",width = 8, height = 6)
volcano_PGMC_CTR_CSF
dev.off()

# -> CSF (adjusted)
DE_CSF_PGMC_CTR_adj = DE_subgroups_CSF_adj_final %>%
  select(Target,UniProtID,log2FC_PGMC_CTR,pvalue_PGMC_CTR,padj_PGMC_CTR) %>%
  add_DE_flag(lfc = "log2FC_PGMC_CTR",padj = "padj_PGMC_CTR",up = "PGMC",down = "CTR",alpha = 0.05)

volcano_PGMC_CTR_CSF_adj = volcano_plot(DE_CSF_PGMC_CTR_adj,"log2FC_PGMC_CTR", 
                                         title = "PGMC vs CTR in CSF (adjusted)",
                                         colors = c('CTR' = "#C25F3D", 
                                                    'ns' = 'lightgrey', 
                                                    'PGMC'= "#3da0c2"),
                                         add_x = 1.3,add_y = 0.1,
                                        label_targets = c("MAPT","NEFH","pTau-217","TAFA5",
                                                          "IL7","SOD1","HTT"))

pdf("plots/volcano_plots/volcano_PGMC_CTR_CSF_adj.pdf",width = 8, height = 6)
volcano_PGMC_CTR_CSF_adj
dev.off()

# -> SERUM (not adjusted)
DE_SERUM_PGMC_CTR = DE_subgroups_SERUM_final %>%
  select(Target,UniProtID,log2FC_PGMC_CTR,pvalue_PGMC_CTR,padj_PGMC_CTR) %>%
  add_DE_flag(lfc = "log2FC_PGMC_CTR",padj = "padj_PGMC_CTR",up = "PGMC",down = "CTR",alpha = 0.05)

volcano_PGMC_CTR_SERUM = volcano_plot(DE_SERUM_PGMC_CTR,"log2FC_PGMC_CTR", 
                                       title = "PGMC vs CTR in SERUM",
                                       colors = c('CTR' = "#C25F3D", 
                                                  'ns' = 'lightgrey', 
                                                  'PGMC'= "#3da0c2"),
                                       add_x = 0.7,add_y = 0.1,
                                      label_targets = c("MAPT","NEFH","pTau-217","TAFA5",
                                                        "IL7","SOD1","HTT"))

pdf("plots/volcano_plots/volcano_PGMC_CTR_SERUM.pdf",width = 8, height = 6)
volcano_PGMC_CTR_SERUM
dev.off()

# -> SERUM (adjusted)
DE_SERUM_PGMC_CTR_adj = DE_subgroups_SERUM_adj_final %>%
  select(Target,UniProtID,log2FC_PGMC_CTR,pvalue_PGMC_CTR,padj_PGMC_CTR) %>%
  add_DE_flag(lfc = "log2FC_PGMC_CTR",padj = "padj_PGMC_CTR",up = "PGMC",down = "CTR",alpha = 0.05)

volcano_PGMC_CTR_SERUM_adj = volcano_plot(DE_SERUM_PGMC_CTR_adj,"log2FC_PGMC_CTR", 
                                           title = "PGMC vs CTR in SERUM (adjusted)",
                                           colors = c('CTR' = "#C25F3D", 
                                                      'ns' = 'lightgrey', 
                                                      'PGMC'= "#3da0c2"),
                                           add_x = 1,add_y = 0.1,
                                          label_targets = c("MAPT","NEFH","pTau-217","TAFA5",
                                                            "IL7","SOD1","HTT"))

pdf("plots/volcano_plots/volcano_PGMC_CTR_SERUM_adj.pdf",width = 8, height = 6)
volcano_PGMC_CTR_SERUM_adj
dev.off()

# -> all together
pdf("plots/volcano_plots/volcano_PGMC_CTR_all.pdf", onefile = TRUE, width = 14, height = 15)
grid.arrange(volcano_PGMC_CTR_plasma, 
             volcano_PGMC_CTR_plasma_adj,
             volcano_PGMC_CTR_SERUM,
             volcano_PGMC_CTR_SERUM_adj,
             volcano_PGMC_CTR_CSF,
             volcano_PGMC_CTR_CSF_adj,nrow = 3, ncol = 2)
dev.off()

## 6: PGMC vs mimic 
# -> plasma (not adjusted)
DE_plasma_PGMC_mimic = DE_subgroups_plasma_final %>%
  select(Target,UniProtID,log2FC_PGMC_mimic,pvalue_PGMC_mimic,padj_PGMC_mimic) %>%
  add_DE_flag(lfc = "log2FC_PGMC_mimic",padj = "padj_PGMC_mimic",up = "PGMC",down = "mimic",alpha = 0.05)

volcano_PGMC_mimic_plasma = volcano_plot(DE_plasma_PGMC_mimic,"log2FC_PGMC_mimic", 
                                       title = "PGMC vs mimic in plasma",
                                       colors = c('mimic' = "#C25F3D", 
                                                  'ns' = 'lightgrey', 
                                                  'PGMC'= "#3da0c2"),
                                       add_x = 0.8,add_y = 0.1)

pdf("plots/volcano_plots/volcano_PGMC_mimic_plasma.pdf",width = 8, height = 6)
volcano_PGMC_mimic_plasma
dev.off()

# -> plasma (adjusted)
DE_plasma_PGMC_mimic_adj = DE_subgroups_plasma_adj_final %>%
  select(Target,UniProtID,log2FC_PGMC_mimic,pvalue_PGMC_mimic,padj_PGMC_mimic) %>%
  add_DE_flag(lfc = "log2FC_PGMC_mimic",padj = "padj_PGMC_mimic",up = "PGMC",down = "mimic",alpha = 0.05)

volcano_PGMC_mimic_plasma_adj = volcano_plot(DE_plasma_PGMC_mimic_adj,"log2FC_PGMC_mimic", 
                                           title = "PGMC vs mimic in plasma (adjusted)",
                                           colors = c('mimic' = "#C25F3D", 
                                                      'ns' = 'lightgrey', 
                                                      'PGMC'= "#3da0c2"),
                                           add_x = 1,add_y = 0.1)

pdf("plots/volcano_plots/volcano_PGMC_mimic_plasma_adj.pdf",width = 8, height = 6)
volcano_PGMC_mimic_plasma_adj
dev.off()

# -> CSF (not adjusted)
DE_CSF_PGMC_mimic = DE_subgroups_CSF_final %>%
  select(Target,UniProtID,log2FC_PGMC_mimic,pvalue_PGMC_mimic,padj_PGMC_mimic) %>%
  add_DE_flag(lfc = "log2FC_PGMC_mimic",padj = "padj_PGMC_mimic",up = "PGMC",down = "mimic",alpha = 0.05)

volcano_PGMC_mimic_CSF = volcano_plot(DE_CSF_PGMC_mimic,"log2FC_PGMC_mimic", 
                                    title = "PGMC vs mimic in CSF",
                                    colors = c('mimic' = "#C25F3D", 
                                               'ns' = 'lightgrey', 
                                               'PGMC'= "#3da0c2"),
                                    add_x = -2.5,add_y = 0.1)

pdf("plots/volcano_plots/volcano_PGMC_mimic_CSF.pdf",width = 8, height = 6)
volcano_PGMC_mimic_CSF
dev.off()

# -> CSF (adjusted)
DE_CSF_PGMC_mimic_adj = DE_subgroups_CSF_adj_final %>%
  select(Target,UniProtID,log2FC_PGMC_mimic,pvalue_PGMC_mimic,padj_PGMC_mimic) %>%
  add_DE_flag(lfc = "log2FC_PGMC_mimic",padj = "padj_PGMC_mimic",up = "PGMC",down = "mimic",alpha = 0.05)

volcano_PGMC_mimic_CSF_adj = volcano_plot(DE_CSF_PGMC_mimic_adj,"log2FC_PGMC_mimic", 
                                        title = "PGMC vs mimic in CSF (adjusted)",
                                        colors = c('mimic' = "#C25F3D", 
                                                   'ns' = 'lightgrey', 
                                                   'PGMC'= "#3da0c2"),
                                        add_x = -2.5,add_y = 0.1)

pdf("plots/volcano_plots/volcano_PGMC_mimic_CSF_adj.pdf",width = 8, height = 6)
volcano_PGMC_mimic_CSF_adj
dev.off()

# -> SERUM (not adjusted)
DE_SERUM_PGMC_mimic = DE_subgroups_SERUM_final %>%
  select(Target,UniProtID,log2FC_PGMC_mimic,pvalue_PGMC_mimic,padj_PGMC_mimic) %>%
  add_DE_flag(lfc = "log2FC_PGMC_mimic",padj = "padj_PGMC_mimic",up = "PGMC",down = "mimic",alpha = 0.05)

volcano_PGMC_mimic_SERUM = volcano_plot(DE_SERUM_PGMC_mimic,"log2FC_PGMC_mimic", 
                                      title = "PGMC vs mimic in SERUM",
                                      colors = c('mimic' = "#C25F3D", 
                                                 'ns' = 'lightgrey', 
                                                 'PGMC'= "#3da0c2"),
                                      add_x = 0.7,add_y = 0.1)

pdf("plots/volcano_plots/volcano_PGMC_mimic_SERUM.pdf",width = 8, height = 6)
volcano_PGMC_mimic_SERUM
dev.off()

# -> SERUM (adjusted)
DE_SERUM_PGMC_mimic_adj = DE_subgroups_SERUM_adj_final %>%
  select(Target,UniProtID,log2FC_PGMC_mimic,pvalue_PGMC_mimic,padj_PGMC_mimic) %>%
  add_DE_flag(lfc = "log2FC_PGMC_mimic",padj = "padj_PGMC_mimic",up = "PGMC",down = "mimic",alpha = 0.05)

volcano_PGMC_mimic_SERUM_adj = volcano_plot(DE_SERUM_PGMC_mimic_adj,"log2FC_PGMC_mimic", 
                                          title = "PGMC vs mimic in SERUM (adjusted)",
                                          colors = c('mimic' = "#C25F3D", 
                                                     'ns' = 'lightgrey', 
                                                     'PGMC'= "#3da0c2"),
                                          add_x = 1,add_y = 0.1)

pdf("plots/volcano_plots/volcano_PGMC_mimic_SERUM_adj.pdf",width = 8, height = 6)
volcano_PGMC_mimic_SERUM_adj
dev.off()

# -> all together
pdf("plots/volcano_plots/volcano_PGMC_mimic_all.pdf", onefile = TRUE, width = 14, height = 15)
grid.arrange(volcano_PGMC_mimic_plasma, 
             volcano_PGMC_mimic_plasma_adj,
             volcano_PGMC_mimic_SERUM,
             volcano_PGMC_mimic_SERUM_adj,
             volcano_PGMC_mimic_CSF,
             volcano_PGMC_mimic_CSF_adj,nrow = 3, ncol = 2)
dev.off()

## 7: C9orf72 vs CTR 
# -> plasma (not adjusted)
DE_plasma_C9orf72_CTR = DE_PGMC_subgroups_plasma_final %>%
  select(Target,UniProtID,log2FC_C9orf72_CTR,pvalue_C9orf72_CTR,padj_C9orf72_CTR) %>%
  add_DE_flag(lfc = "log2FC_C9orf72_CTR",padj = "padj_C9orf72_CTR",up = "C9orf72",down = "CTR",alpha = 0.05)

volcano_C9orf72_CTR_plasma = volcano_plot(DE_plasma_C9orf72_CTR,"log2FC_C9orf72_CTR", 
                                        title = "C9orf72 vs CTR in plasma",
                                        colors = c('CTR' = "#C25F3D", 
                                                   'ns' = 'lightgrey', 
                                                   'C9orf72'= "#3da0c2"),
                                        add_x = 2.3,add_y = -0.1)

pdf("plots/volcano_plots/volcano_C9orf72_CTR_plasma.pdf",width = 8, height = 6)
volcano_C9orf72_CTR_plasma
dev.off()

# -> plasma (adjusted)
DE_plasma_C9orf72_CTR_adj = DE_PGMC_subgroups_plasma_adj_final %>%
  select(Target,UniProtID,log2FC_C9orf72_CTR,pvalue_C9orf72_CTR,padj_C9orf72_CTR) %>%
  add_DE_flag(lfc = "log2FC_C9orf72_CTR",padj = "padj_C9orf72_CTR",up = "C9orf72",down = "CTR",alpha = 0.05)

volcano_C9orf72_CTR_plasma_adj = volcano_plot(DE_plasma_C9orf72_CTR_adj,"log2FC_C9orf72_CTR", 
                                            title = "C9orf72 vs CTR in plasma (adjusted)",
                                            colors = c('CTR' = "#C25F3D", 
                                                       'ns' = 'lightgrey', 
                                                       'C9orf72'= "#3da0c2"),
                                            add_x = 2.5,add_y = -0.1)

pdf("plots/volcano_plots/volcano_C9orf72_CTR_plasma_adj.pdf",width = 8, height = 6)
volcano_C9orf72_CTR_plasma_adj
dev.off()

# -> CSF (not adjusted)
DE_CSF_C9orf72_CTR = DE_PGMC_subgroups_CSF_final %>%
  select(Target,UniProtID,log2FC_C9orf72_CTR,pvalue_C9orf72_CTR,padj_C9orf72_CTR) %>%
  add_DE_flag(lfc = "log2FC_C9orf72_CTR",padj = "padj_C9orf72_CTR",up = "C9orf72",down = "CTR",alpha = 0.05)

volcano_C9orf72_CTR_CSF = volcano_plot(DE_CSF_C9orf72_CTR,"log2FC_C9orf72_CTR", 
                                     title = "C9orf72 vs CTR in CSF",
                                     colors = c('CTR' = "#C25F3D", 
                                                'ns' = 'lightgrey', 
                                                'C9orf72'= "#3da0c2"),
                                     add_x = 2,add_y = 0.1)

pdf("plots/volcano_plots/volcano_C9orf72_CTR_CSF.pdf",width = 8, height = 6)
volcano_C9orf72_CTR_CSF
dev.off()

# -> CSF (adjusted)
DE_CSF_C9orf72_CTR_adj = DE_PGMC_subgroups_CSF_adj_final %>%
  select(Target,UniProtID,log2FC_C9orf72_CTR,pvalue_C9orf72_CTR,padj_C9orf72_CTR) %>%
  add_DE_flag(lfc = "log2FC_C9orf72_CTR",padj = "padj_C9orf72_CTR",up = "C9orf72",down = "CTR",alpha = 0.05)

volcano_C9orf72_CTR_CSF_adj = volcano_plot(DE_CSF_C9orf72_CTR_adj,"log2FC_C9orf72_CTR", 
                                         title = "C9orf72 vs CTR in CSF (adjusted)",
                                         colors = c('CTR' = "#C25F3D", 
                                                    'ns' = 'lightgrey', 
                                                    'C9orf72'= "#3da0c2"),
                                         add_x = 2,add_y = 0.1)

pdf("plots/volcano_plots/volcano_C9orf72_CTR_CSF_adj.pdf",width = 8, height = 6)
volcano_C9orf72_CTR_CSF_adj
dev.off()

# -> SERUM (not adjusted)
DE_SERUM_C9orf72_CTR = DE_PGMC_subgroups_SERUM_final %>%
  select(Target,UniProtID,log2FC_C9orf72_CTR,pvalue_C9orf72_CTR,padj_C9orf72_CTR) %>%
  add_DE_flag(lfc = "log2FC_C9orf72_CTR",padj = "padj_C9orf72_CTR",up = "C9orf72",down = "CTR",alpha = 0.05)

volcano_C9orf72_CTR_SERUM = volcano_plot(DE_SERUM_C9orf72_CTR,"log2FC_C9orf72_CTR", 
                                       title = "C9orf72 vs CTR in SERUM",
                                       colors = c('CTR' = "#C25F3D", 
                                                  'ns' = 'lightgrey', 
                                                  'C9orf72'= "#3da0c2"),
                                       add_x = 2,add_y = 0.1)

pdf("plots/volcano_plots/volcano_C9orf72_CTR_SERUM.pdf",width = 8, height = 6)
volcano_C9orf72_CTR_SERUM
dev.off()

# -> SERUM (adjusted)
DE_SERUM_C9orf72_CTR_adj = DE_PGMC_subgroups_SERUM_adj_final %>%
  select(Target,UniProtID,log2FC_C9orf72_CTR,pvalue_C9orf72_CTR,padj_C9orf72_CTR) %>%
  add_DE_flag(lfc = "log2FC_C9orf72_CTR",padj = "padj_C9orf72_CTR",up = "C9orf72",down = "CTR",alpha = 0.05)

volcano_C9orf72_CTR_SERUM_adj = volcano_plot(DE_SERUM_C9orf72_CTR_adj,"log2FC_C9orf72_CTR", 
                                           title = "C9orf72 vs CTR in SERUM (adjusted)",
                                           colors = c('CTR' = "#C25F3D", 
                                                      'ns' = 'lightgrey', 
                                                      'C9orf72'= "#3da0c2"),
                                           add_x = 2.8,add_y = 0.1)

pdf("plots/volcano_plots/volcano_C9orf72_CTR_SERUM_adj.pdf",width = 8, height = 6)
volcano_C9orf72_CTR_SERUM_adj
dev.off()

# -> all together
pdf("plots/volcano_plots/volcano_C9orf72_CTR_all.pdf", onefile = TRUE, width = 14, height = 15)
grid.arrange(volcano_C9orf72_CTR_plasma, 
             volcano_C9orf72_CTR_plasma_adj,
             volcano_C9orf72_CTR_SERUM,
             volcano_C9orf72_CTR_SERUM_adj,
             volcano_C9orf72_CTR_CSF,
             volcano_C9orf72_CTR_CSF_adj,nrow = 3, ncol = 2)
dev.off()

## 8: SOD1 vs CTR 
# -> plasma (not adjusted)
DE_plasma_SOD1_CTR = DE_PGMC_subgroups_plasma_final %>%
  select(Target,UniProtID,log2FC_SOD1_CTR,pvalue_SOD1_CTR,padj_SOD1_CTR) %>%
  add_DE_flag(lfc = "log2FC_SOD1_CTR",padj = "padj_SOD1_CTR",up = "SOD1",down = "CTR",alpha = 0.05)

volcano_SOD1_CTR_plasma = volcano_plot(DE_plasma_SOD1_CTR,"log2FC_SOD1_CTR", 
                                          title = "SOD1 vs CTR in plasma",
                                          colors = c('CTR' = "#C25F3D", 
                                                     'ns' = 'lightgrey', 
                                                     'SOD1'= "#3da0c2"),
                                          add_x = 1,add_y = 0.1)

pdf("plots/volcano_plots/volcano_SOD1_CTR_plasma.pdf",width = 8, height = 6)
volcano_SOD1_CTR_plasma
dev.off()

# -> plasma (adjusted)
DE_plasma_SOD1_CTR_adj = DE_PGMC_subgroups_plasma_adj_final %>%
  select(Target,UniProtID,log2FC_SOD1_CTR,pvalue_SOD1_CTR,padj_SOD1_CTR) %>%
  add_DE_flag(lfc = "log2FC_SOD1_CTR",padj = "padj_SOD1_CTR",up = "SOD1",down = "CTR",alpha = 0.05)

volcano_SOD1_CTR_plasma_adj = volcano_plot(DE_plasma_SOD1_CTR_adj,"log2FC_SOD1_CTR", 
                                              title = "SOD1 vs CTR in plasma (adjusted)",
                                              colors = c('CTR' = "#C25F3D", 
                                                         'ns' = 'lightgrey', 
                                                         'SOD1'= "#3da0c2"),
                                              add_x = 1.2,add_y = 0.1)

pdf("plots/volcano_plots/volcano_SOD1_CTR_plasma_adj.pdf",width = 8, height = 6)
volcano_SOD1_CTR_plasma_adj
dev.off()

# -> CSF (not adjusted)
DE_CSF_SOD1_CTR = DE_PGMC_subgroups_CSF_final %>%
  select(Target,UniProtID,log2FC_SOD1_CTR,pvalue_SOD1_CTR,padj_SOD1_CTR) %>%
  add_DE_flag(lfc = "log2FC_SOD1_CTR",padj = "padj_SOD1_CTR",up = "SOD1",down = "CTR",alpha = 0.05)

volcano_SOD1_CTR_CSF = volcano_plot(DE_CSF_SOD1_CTR,"log2FC_SOD1_CTR", 
                                       title = "SOD1 vs CTR in CSF",
                                       colors = c('CTR' = "#C25F3D", 
                                                  'ns' = 'lightgrey', 
                                                  'SOD1'= "#3da0c2"),
                                       add_x = -2.5,add_y = 0.1)

pdf("plots/volcano_plots/volcano_SOD1_CTR_CSF.pdf",width = 8, height = 6)
volcano_SOD1_CTR_CSF
dev.off()

# -> CSF (adjusted)
DE_CSF_SOD1_CTR_adj = DE_PGMC_subgroups_CSF_adj_final %>%
  select(Target,UniProtID,log2FC_SOD1_CTR,pvalue_SOD1_CTR,padj_SOD1_CTR) %>%
  add_DE_flag(lfc = "log2FC_SOD1_CTR",padj = "padj_SOD1_CTR",up = "SOD1",down = "CTR",alpha = 0.05)

volcano_SOD1_CTR_CSF_adj = volcano_plot(DE_CSF_SOD1_CTR_adj,"log2FC_SOD1_CTR", 
                                           title = "SOD1 vs CTR in CSF (adjusted)",
                                           colors = c('CTR' = "#C25F3D", 
                                                      'ns' = 'lightgrey', 
                                                      'SOD1'= "#3da0c2"),
                                           add_x = 1,add_y = 0.1)

pdf("plots/volcano_plots/volcano_SOD1_CTR_CSF_adj.pdf",width = 8, height = 6)
volcano_SOD1_CTR_CSF_adj
dev.off()

# -> SERUM (not adjusted)
DE_SERUM_SOD1_CTR = DE_PGMC_subgroups_SERUM_final %>%
  select(Target,UniProtID,log2FC_SOD1_CTR,pvalue_SOD1_CTR,padj_SOD1_CTR) %>%
  add_DE_flag(lfc = "log2FC_SOD1_CTR",padj = "padj_SOD1_CTR",up = "SOD1",down = "CTR",alpha = 0.05)

volcano_SOD1_CTR_SERUM = volcano_plot(DE_SERUM_SOD1_CTR,"log2FC_SOD1_CTR", 
                                         title = "SOD1 vs CTR in SERUM",
                                         colors = c('CTR' = "#C25F3D", 
                                                    'ns' = 'lightgrey', 
                                                    'SOD1'= "#3da0c2"),
                                         add_x = 0.7,add_y = 0.1)

pdf("plots/volcano_plots/volcano_SOD1_CTR_SERUM.pdf",width = 8, height = 6)
volcano_SOD1_CTR_SERUM
dev.off()

# -> SERUM (adjusted)
DE_SERUM_SOD1_CTR_adj = DE_PGMC_subgroups_SERUM_adj_final %>%
  select(Target,UniProtID,log2FC_SOD1_CTR,pvalue_SOD1_CTR,padj_SOD1_CTR) %>%
  add_DE_flag(lfc = "log2FC_SOD1_CTR",padj = "padj_SOD1_CTR",up = "SOD1",down = "CTR",alpha = 0.05)

volcano_SOD1_CTR_SERUM_adj = volcano_plot(DE_SERUM_SOD1_CTR_adj,"log2FC_SOD1_CTR", 
                                             title = "SOD1 vs CTR in SERUM (adjusted)",
                                             colors = c('CTR' = "#C25F3D", 
                                                        'ns' = 'lightgrey', 
                                                        'SOD1'= "#3da0c2"),
                                             add_x = 0.8,add_y = 0.1)

pdf("plots/volcano_plots/volcano_SOD1_CTR_SERUM_adj.pdf",width = 8, height = 6)
volcano_SOD1_CTR_SERUM_adj
dev.off()

# -> all together
pdf("plots/volcano_plots/volcano_SOD1_CTR_all.pdf", onefile = TRUE, width = 14, height = 15)
grid.arrange(volcano_SOD1_CTR_plasma, 
             volcano_SOD1_CTR_plasma_adj,
             volcano_SOD1_CTR_SERUM,
             volcano_SOD1_CTR_SERUM_adj,
             volcano_SOD1_CTR_CSF,
             volcano_SOD1_CTR_CSF_adj,nrow = 3, ncol = 2)
dev.off()

## 9: TARDBP vs CTR 
# -> plasma (not adjusted)
DE_plasma_TARDBP_CTR = DE_PGMC_subgroups_plasma_final %>%
  select(Target,UniProtID,log2FC_TARDBP_CTR,pvalue_TARDBP_CTR,padj_TARDBP_CTR) %>%
  add_DE_flag(lfc = "log2FC_TARDBP_CTR",padj = "padj_TARDBP_CTR",up = "TARDBP",down = "CTR",alpha = 0.05)

volcano_TARDBP_CTR_plasma = volcano_plot(DE_plasma_TARDBP_CTR,"log2FC_TARDBP_CTR", 
                                       title = "TARDBP vs CTR in plasma",
                                       colors = c('CTR' = "#C25F3D", 
                                                  'ns' = 'lightgrey', 
                                                  'TARDBP'= "#3da0c2"),
                                       add_x = 1.7,add_y = 0.1)

pdf("plots/volcano_plots/volcano_TARDBP_CTR_plasma.pdf",width = 8, height = 6)
volcano_TARDBP_CTR_plasma
dev.off()

# -> plasma (adjusted)
DE_plasma_TARDBP_CTR_adj = DE_PGMC_subgroups_plasma_adj_final %>%
  select(Target,UniProtID,log2FC_TARDBP_CTR,pvalue_TARDBP_CTR,padj_TARDBP_CTR) %>%
  add_DE_flag(lfc = "log2FC_TARDBP_CTR",padj = "padj_TARDBP_CTR",up = "TARDBP",down = "CTR",alpha = 0.05)

volcano_TARDBP_CTR_plasma_adj = volcano_plot(DE_plasma_TARDBP_CTR_adj,"log2FC_TARDBP_CTR", 
                                           title = "TARDBP vs CTR in plasma (adjusted)",
                                           colors = c('CTR' = "#C25F3D", 
                                                      'ns' = 'lightgrey', 
                                                      'TARDBP'= "#3da0c2"),
                                           add_x = 1.7,add_y = 0.05)

pdf("plots/volcano_plots/volcano_TARDBP_CTR_plasma_adj.pdf",width = 8, height = 6)
volcano_TARDBP_CTR_plasma_adj
dev.off()

# -> CSF (not adjusted)
DE_CSF_TARDBP_CTR = DE_PGMC_subgroups_CSF_final %>%
  select(Target,UniProtID,log2FC_TARDBP_CTR,pvalue_TARDBP_CTR,padj_TARDBP_CTR) %>%
  add_DE_flag(lfc = "log2FC_TARDBP_CTR",padj = "padj_TARDBP_CTR",up = "TARDBP",down = "CTR",alpha = 0.05)

volcano_TARDBP_CTR_CSF = volcano_plot(DE_CSF_TARDBP_CTR,"log2FC_TARDBP_CTR", 
                                    title = "TARDBP vs CTR in CSF",
                                    colors = c('CTR' = "#C25F3D", 
                                               'ns' = 'lightgrey', 
                                               'TARDBP'= "#3da0c2"),
                                    add_x = 5.5,add_y = 0.1)

pdf("plots/volcano_plots/volcano_TARDBP_CTR_CSF.pdf",width = 8, height = 6)
volcano_TARDBP_CTR_CSF
dev.off()

# -> CSF (adjusted)
DE_CSF_TARDBP_CTR_adj = DE_PGMC_subgroups_CSF_adj_final %>%
  select(Target,UniProtID,log2FC_TARDBP_CTR,pvalue_TARDBP_CTR,padj_TARDBP_CTR) %>%
  add_DE_flag(lfc = "log2FC_TARDBP_CTR",padj = "padj_TARDBP_CTR",up = "TARDBP",down = "CTR",alpha = 0.05)

volcano_TARDBP_CTR_CSF_adj = volcano_plot(DE_CSF_TARDBP_CTR_adj,"log2FC_TARDBP_CTR", 
                                        title = "TARDBP vs CTR in CSF (adjusted)",
                                        colors = c('CTR' = "#C25F3D", 
                                                   'ns' = 'lightgrey', 
                                                   'TARDBP'= "#3da0c2"),
                                        add_x = 5.5,add_y = 0.1)

pdf("plots/volcano_plots/volcano_TARDBP_CTR_CSF_adj.pdf",width = 8, height = 6)
volcano_TARDBP_CTR_CSF_adj
dev.off()

# -> SERUM (not adjusted)
DE_SERUM_TARDBP_CTR = DE_PGMC_subgroups_SERUM_final %>%
  select(Target,UniProtID,log2FC_TARDBP_CTR,pvalue_TARDBP_CTR,padj_TARDBP_CTR) %>%
  add_DE_flag(lfc = "log2FC_TARDBP_CTR",padj = "padj_TARDBP_CTR",up = "TARDBP",down = "CTR",alpha = 0.05)

volcano_TARDBP_CTR_SERUM = volcano_plot(DE_SERUM_TARDBP_CTR,"log2FC_TARDBP_CTR", 
                                      title = "TARDBP vs CTR in SERUM",
                                      colors = c('CTR' = "#C25F3D", 
                                                 'ns' = 'lightgrey', 
                                                 'TARDBP'= "#3da0c2"),
                                      add_x = 2,add_y = 0.1)

pdf("plots/volcano_plots/volcano_TARDBP_CTR_SERUM.pdf",width = 8, height = 6)
volcano_TARDBP_CTR_SERUM
dev.off()

# -> SERUM (adjusted)
DE_SERUM_TARDBP_CTR_adj = DE_PGMC_subgroups_SERUM_adj_final %>%
  select(Target,UniProtID,log2FC_TARDBP_CTR,pvalue_TARDBP_CTR,padj_TARDBP_CTR) %>%
  add_DE_flag(lfc = "log2FC_TARDBP_CTR",padj = "padj_TARDBP_CTR",up = "TARDBP",down = "CTR",alpha = 0.05)

volcano_TARDBP_CTR_SERUM_adj = volcano_plot(DE_SERUM_TARDBP_CTR_adj,"log2FC_TARDBP_CTR", 
                                          title = "TARDBP vs CTR in SERUM (adjusted)",
                                          colors = c('CTR' = "#C25F3D", 
                                                     'ns' = 'lightgrey', 
                                                     'TARDBP'= "#3da0c2"),
                                          add_x = -3.7,add_y = 0.1)

pdf("plots/volcano_plots/volcano_TARDBP_CTR_SERUM_adj.pdf",width = 8, height = 6)
volcano_TARDBP_CTR_SERUM_adj
dev.off()

# -> all together
pdf("plots/volcano_plots/volcano_TARDBP_CTR_all.pdf", onefile = TRUE, width = 14, height = 15)
grid.arrange(volcano_TARDBP_CTR_plasma, 
             volcano_TARDBP_CTR_plasma_adj,
             volcano_TARDBP_CTR_SERUM,
             volcano_TARDBP_CTR_SERUM_adj,
             volcano_TARDBP_CTR_CSF,
             volcano_TARDBP_CTR_CSF_adj,nrow = 3, ncol = 2)
dev.off()


## 10: C9orf72 vs SOD1 
# -> plasma (not adjusted)
DE_plasma_C9orf72_SOD1 = DE_PGMC_subgroups_plasma_final %>%
  select(Target,UniProtID,log2FC_C9orf72_SOD1,pvalue_C9orf72_SOD1,padj_C9orf72_SOD1) %>%
  add_DE_flag(lfc = "log2FC_C9orf72_SOD1",padj = "padj_C9orf72_SOD1",up = "C9orf72",down = "SOD1",alpha = 0.05)

volcano_C9orf72_SOD1_plasma = volcano_plot(DE_plasma_C9orf72_SOD1,"log2FC_C9orf72_SOD1", 
                                          title = "C9orf72 vs SOD1 in plasma",
                                          colors = c('SOD1' = "#C25F3D", 
                                                     'ns' = 'lightgrey', 
                                                     'C9orf72'= "#3da0c2"),
                                          add_x = 1.5,add_y = 0.1)

pdf("plots/volcano_plots/volcano_C9orf72_SOD1_plasma.pdf",width = 8, height = 6)
volcano_C9orf72_SOD1_plasma
dev.off()

# -> plasma (adjusted)
DE_plasma_C9orf72_SOD1_adj = DE_PGMC_subgroups_plasma_adj_final %>%
  select(Target,UniProtID,log2FC_C9orf72_SOD1,pvalue_C9orf72_SOD1,padj_C9orf72_SOD1) %>%
  add_DE_flag(lfc = "log2FC_C9orf72_SOD1",padj = "padj_C9orf72_SOD1",up = "C9orf72",down = "SOD1",alpha = 0.05)

volcano_C9orf72_SOD1_plasma_adj = volcano_plot(DE_plasma_C9orf72_SOD1_adj,"log2FC_C9orf72_SOD1", 
                                              title = "C9orf72 vs SOD1 in plasma (adjusted)",
                                              colors = c('SOD1' = "#C25F3D", 
                                                         'ns' = 'lightgrey', 
                                                         'C9orf72'= "#3da0c2"),
                                              add_x = 1.5,add_y = 0.1)

pdf("plots/volcano_plots/volcano_C9orf72_SOD1_plasma_adj.pdf",width = 8, height = 6)
volcano_C9orf72_SOD1_plasma_adj
dev.off()

# -> CSF (not adjusted)
DE_CSF_C9orf72_SOD1 = DE_PGMC_subgroups_CSF_final %>%
  select(Target,UniProtID,log2FC_C9orf72_SOD1,pvalue_C9orf72_SOD1,padj_C9orf72_SOD1) %>%
  add_DE_flag(lfc = "log2FC_C9orf72_SOD1",padj = "padj_C9orf72_SOD1",up = "C9orf72",down = "SOD1",alpha = 0.05)

volcano_C9orf72_SOD1_CSF = volcano_plot(DE_CSF_C9orf72_SOD1,"log2FC_C9orf72_SOD1", 
                                       title = "C9orf72 vs SOD1 in CSF",
                                       colors = c('SOD1' = "#C25F3D", 
                                                  'ns' = 'lightgrey', 
                                                  'C9orf72'= "#3da0c2"),
                                       add_x = 1.5,add_y = 0.1)

pdf("plots/volcano_plots/volcano_C9orf72_SOD1_CSF.pdf",width = 8, height = 6)
volcano_C9orf72_SOD1_CSF
dev.off()

# -> CSF (adjusted)
DE_CSF_C9orf72_SOD1_adj = DE_PGMC_subgroups_CSF_adj_final %>%
  select(Target,UniProtID,log2FC_C9orf72_SOD1,pvalue_C9orf72_SOD1,padj_C9orf72_SOD1) %>%
  add_DE_flag(lfc = "log2FC_C9orf72_SOD1",padj = "padj_C9orf72_SOD1",up = "C9orf72",down = "SOD1",alpha = 0.05)

volcano_C9orf72_SOD1_CSF_adj = volcano_plot(DE_CSF_C9orf72_SOD1_adj,"log2FC_C9orf72_SOD1", 
                                           title = "C9orf72 vs SOD1 in CSF (adjusted)",
                                           colors = c('SOD1' = "#C25F3D", 
                                                      'ns' = 'lightgrey', 
                                                      'C9orf72'= "#3da0c2"),
                                           add_x = 1.5,add_y = 0.1)

pdf("plots/volcano_plots/volcano_C9orf72_SOD1_CSF_adj.pdf",width = 8, height = 6)
volcano_C9orf72_SOD1_CSF_adj
dev.off()

# -> SERUM (not adjusted)
DE_SERUM_C9orf72_SOD1 = DE_PGMC_subgroups_SERUM_final %>%
  select(Target,UniProtID,log2FC_C9orf72_SOD1,pvalue_C9orf72_SOD1,padj_C9orf72_SOD1) %>%
  add_DE_flag(lfc = "log2FC_C9orf72_SOD1",padj = "padj_C9orf72_SOD1",up = "C9orf72",down = "SOD1",alpha = 0.05)

volcano_C9orf72_SOD1_SERUM = volcano_plot(DE_SERUM_C9orf72_SOD1,"log2FC_C9orf72_SOD1", 
                                         title = "C9orf72 vs SOD1 in SERUM",
                                         colors = c('SOD1' = "#C25F3D", 
                                                    'ns' = 'lightgrey', 
                                                    'C9orf72'= "#3da0c2"),
                                         add_x = 1.7,add_y = 0.1)

pdf("plots/volcano_plots/volcano_C9orf72_SOD1_SERUM.pdf",width = 8, height = 6)
volcano_C9orf72_SOD1_SERUM
dev.off()

# -> SERUM (adjusted)
DE_SERUM_C9orf72_SOD1_adj = DE_PGMC_subgroups_SERUM_adj_final %>%
  select(Target,UniProtID,log2FC_C9orf72_SOD1,pvalue_C9orf72_SOD1,padj_C9orf72_SOD1) %>%
  add_DE_flag(lfc = "log2FC_C9orf72_SOD1",padj = "padj_C9orf72_SOD1",up = "C9orf72",down = "SOD1",alpha = 0.05)

volcano_C9orf72_SOD1_SERUM_adj = volcano_plot(DE_SERUM_C9orf72_SOD1_adj,"log2FC_C9orf72_SOD1", 
                                             title = "C9orf72 vs SOD1 in SERUM (adjusted)",
                                             colors = c('SOD1' = "#C25F3D", 
                                                        'ns' = 'lightgrey', 
                                                        'C9orf72'= "#3da0c2"),
                                             add_x = 1.8,add_y = 0.1)

pdf("plots/volcano_plots/volcano_C9orf72_SOD1_SERUM_adj.pdf",width = 8, height = 6)
volcano_C9orf72_SOD1_SERUM_adj
dev.off()

# -> all together
pdf("plots/volcano_plots/volcano_C9orf72_SOD1_all.pdf", onefile = TRUE, width = 14, height = 15)
grid.arrange(volcano_C9orf72_SOD1_plasma, 
             volcano_C9orf72_SOD1_plasma_adj,
             volcano_C9orf72_SOD1_SERUM,
             volcano_C9orf72_SOD1_SERUM_adj,
             volcano_C9orf72_SOD1_CSF,
             volcano_C9orf72_SOD1_CSF_adj,nrow = 3, ncol = 2)
dev.off()






## signed p-values 
# ALS vs CTR + ALS vs PGMC
colors = c('ALS' = '#32A55E',
  'PGMC' = '#96AA9A',
  'ALS sign. both' = '#5E718B',
  'ns' = '#B4BFC5',
  'CTR' = '#CF7041',
  "CTR and PGMC" = "#E5C5BD")

signed_pvalue_ALS_PGMC_CTR <- function(data_to_plot,title_plot){
  ggplot(data_to_plot,aes(signed_ALS_CTR,signed_ALS_PGMC)) +
    geom_point(size = 2.25, color = "grey", alpha = 0.5) +
    geom_point(aes(colour = DE_new), size = 2.25, alpha = 0.8) +
    scale_colour_manual(values = colors) +
    theme(panel.background = element_blank(),
          plot.background = element_blank(), 
          axis.line = element_line()) +
    ggtitle(title_plot)+
    theme(text = element_text(size = 15),               # Base font size for everything
          axis.title = element_text(size = 18),         # Axis titles
          axis.text = element_text(size = 16),          # Axis tick labels
          plot.title = element_text(size = 18, hjust = 0.5,face = "bold"),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 15))+
    scale_x_continuous(breaks = seq(-50, 20, by = 5))+
    scale_y_continuous(breaks = seq(-30, 140, by = 10))+
    xlab(paste0("ALS vs CTR (sign(LFC) x -log10(p-adj))"))+
    theme(axis.title = element_text(size = 15), axis.text=element_text(size=10))+
    ylab(paste0("ALS vs PGMC (sign(LFC) x -log10(p-adj))"))+
    theme(axis.title = element_text(size = 15), axis.text=element_text(size=10))+
    geom_vline(xintercept = log10(0.05), linetype="dashed",color = "grey",alpha = 0.6)+
    geom_vline(xintercept = -log10(0.05), linetype="dashed",color = "grey",alpha = 0.6)+
    geom_hline(yintercept = log10(0.05), linetype="dashed",color = "grey",alpha = 0.6)+
    geom_hline(yintercept = -log10(0.05), linetype="dashed",color = "grey",alpha = 0.6)+
    #geom_abline(linetype="dashed",color = "grey",alpha = 0.6) +
    annotate(geom = "text", label = "FDR = 5%",x= -log10(0.05)+13,y= -log10(0.05)-3.15,
             color = "grey",size = 5) +
    annotate(geom = "text", label = "FDR = 5%",x=-log10(0.05)+13,y= -log10(0.05)+0.5,
             color = "grey",size = 5) +
    geom_text_repel(data = subset(data_to_plot, DE_new != "ns"), 
                     aes(label=Target),  color="black",
                     nudge_y = 0.5,
                     nudge_x = 0.5,  
                     min.segment.length = unit(0.1, "lines"), 
                     max.overlaps = 40, size = 5)  }

# -> Plasma
signed_subgroups_ALS_CTR_and_ALS_PGMC = merge(DE_subgroups_plasma_ALS_CTR %>%
  mutate(signed_ALS_CTR = sign(log2FC_ALS_CTR)*log10_padj) %>%
  rename(DE_ALS_CTR = DE) %>%
  select(Target,UniProtID,signed_ALS_CTR,DE_ALS_CTR),
  DE_subgroups_plasma_ALS_PGMC %>%
    mutate(signed_ALS_PGMC = sign(log2FC_ALS_PGMC)*log10_padj) %>%
    rename(DE_ALS_PGMC = DE) %>%
    select(Target,UniProtID,signed_ALS_PGMC,DE_ALS_PGMC)
  ) %>%
  mutate(DE_new = ifelse(DE_ALS_CTR == "ALS" & DE_ALS_PGMC == "ALS","ALS sign. both",
                         ifelse(DE_ALS_CTR == "CTR" & DE_ALS_PGMC == "PGMC","CTR and PGMC",
                         ifelse(DE_ALS_CTR == "CTR","CTR",
                                ifelse(DE_ALS_CTR == "ALS","ALS",
                                       ifelse(DE_ALS_PGMC == "ALS","ALS",
                                              ifelse(DE_ALS_PGMC == "PGMC","PGMC","ns")))))))

colors = c('ALS' = '#32A55E',
           'ALS sign. both' = '#5E718B',
           'ns' = '#B4BFC5')
pdf("plots/signed_ALSvsCTR_ALSvsPGMC_plasma.pdf")
signed_pvalue_ALS_PGMC_CTR(signed_subgroups_ALS_CTR_and_ALS_PGMC,"ALS vs CTR + ALS vs PGMC (PLASMA)")
dev.off()

# -> CSF
colors = c('ALS' = '#32A55E',
           'ALS sign. both' = '#5E718B',
           'ns' = '#B4BFC5',
           'PGMC' = '#96AA9A')
signed_subgroups_ALS_CTR_and_ALS_PGMC = merge(DE_subgroups_CSF_ALS_CTR %>%
                                                mutate(signed_ALS_CTR = sign(log2FC_ALS_CTR)*log10_padj) %>%
                                                rename(DE_ALS_CTR = DE) %>%
                                                select(Target,UniProtID,signed_ALS_CTR,DE_ALS_CTR),
                                              DE_subgroups_CSF_ALS_PGMC %>%
                                                mutate(signed_ALS_PGMC = sign(log2FC_ALS_PGMC)*log10_padj) %>%
                                                rename(DE_ALS_PGMC = DE) %>%
                                                select(Target,UniProtID,signed_ALS_PGMC,DE_ALS_PGMC)
) %>%
  mutate(DE_new = ifelse(DE_ALS_CTR == "ALS" & DE_ALS_PGMC == "ALS","ALS sign. both",
                         ifelse(DE_ALS_CTR == "CTR" & DE_ALS_PGMC == "PGMC","CTR and PGMC",
                                ifelse(DE_ALS_CTR == "CTR","CTR",
                                       ifelse(DE_ALS_CTR == "ALS","ALS",
                                              ifelse(DE_ALS_PGMC == "ALS","ALS",
                                                     ifelse(DE_ALS_PGMC == "PGMC","PGMC","ns")))))))

pdf("plots/signed_ALSvsCTR_ALSvsPGMC_CSF.pdf")
signed_pvalue_ALS_PGMC_CTR(signed_subgroups_ALS_CTR_and_ALS_PGMC,"ALS vs CTR + ALS vs PGMC (CSF)")
dev.off()

# -> Serum
colors = c('ALS' = '#32A55E',
           'PGMC' = '#96AA9A',
           'ALS sign. both' = '#5E718B',
           'ns' = '#B4BFC5',
           'CTR' = '#CF7041',
           "CTR and PGMC" = "#E5C5BD")
signed_subgroups_ALS_CTR_and_ALS_PGMC = merge(DE_subgroups_SERUM_ALS_CTR %>%
                                                mutate(signed_ALS_CTR = sign(log2FC_ALS_CTR)*log10_padj) %>%
                                                rename(DE_ALS_CTR = DE) %>%
                                                select(Target,UniProtID,signed_ALS_CTR,DE_ALS_CTR),
                                              DE_subgroups_SERUM_ALS_PGMC %>%
                                                mutate(signed_ALS_PGMC = sign(log2FC_ALS_PGMC)*log10_padj) %>%
                                                rename(DE_ALS_PGMC = DE) %>%
                                                select(Target,UniProtID,signed_ALS_PGMC,DE_ALS_PGMC)
) %>%
  mutate(DE_new = ifelse(DE_ALS_CTR == "ALS" & DE_ALS_PGMC == "ALS","ALS sign. both",
                         ifelse(DE_ALS_CTR == "CTR" & DE_ALS_PGMC == "PGMC","CTR and PGMC",
                                ifelse(DE_ALS_CTR == "CTR","CTR",
                                       ifelse(DE_ALS_CTR == "ALS","ALS",
                                              ifelse(DE_ALS_PGMC == "ALS","ALS",
                                                     ifelse(DE_ALS_PGMC == "PGMC","PGMC","ns")))))))

pdf("plots/signed_ALSvsCTR_ALSvsPGMC_SERUM.pdf")
signed_pvalue_ALS_PGMC_CTR(signed_subgroups_ALS_CTR_and_ALS_PGMC,"ALS vs CTR + ALS vs PGMC (SERUM)")
dev.off()

# ALS vs CTR + PGMC1 vs CTR 
colors = c('ALS' = '#32A55E',
           'mimic' = '#96AA9A',
           'CTR sign. both' = '#5E718B',
           'ns' = '#B4BFC5',
           'CTR' = '#CF7041',
           "ALS and mimic" = "#E5C5BD")
signed_pvalue_ALS_CTR_mimic <- function(data_to_plot,title_plot){
  ggplot(data_to_plot,aes(signed_ALS_CTR,signed_mimic_CTR)) +
    geom_point(aes(colour = DE_new),size = 2.25, color = "grey", alpha = 0.5) +
    geom_point(aes(colour = DE_new),size = 2.25, alpha = 0.5) +
    theme(panel.background = element_blank(),
          plot.background = element_blank(), 
          axis.line = element_line()) +
    ggtitle(title_plot)+
    scale_x_continuous(breaks = seq(-50, 20, by = 5))+
    scale_y_continuous(breaks = seq(-30, 140, by = 10))+
    xlab(paste0("ALS vs CTR (sign(LFC) x -log10(p-adj))"))+
    theme(axis.title = element_text(size = 15), axis.text=element_text(size=10))+
    ylab(paste0("mimic vs CTR (sign(LFC) x -log10(p-adj))"))+
    theme(axis.title = element_text(size = 18),         # Axis titles
          axis.text = element_text(size = 16),          # Axis tick labels
          plot.title = element_text(size = 18, hjust = 0.5,face = "bold"),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 15))+
    geom_vline(xintercept = log10(0.05), linetype="dashed",color = "grey",alpha = 0.6)+
    geom_vline(xintercept = -log10(0.05), linetype="dashed",color = "grey",alpha = 0.6)+
    geom_hline(yintercept = log10(0.05), linetype="dashed",color = "grey",alpha = 0.6)+
    geom_hline(yintercept = -log10(0.05), linetype="dashed",color = "grey",alpha = 0.6)+
    #geom_abline(linetype="dashed",color = "grey",alpha = 0.6) +
    annotate(geom = "text", label = "FDR = 5%",x= -log10(0.05)+0.9,y= log10(0.05)-0.15,
             color = "grey", size = 5) +
    annotate(geom = "text", label = "FDR = 5%",x=log10(0.05)-0.9,y= - log10(0.05)+0.15,
             color = "grey", size = 5) +
    geom_text_repel(data = subset(data_to_plot, DE_new != "ns"), 
                    aes(label=Target),  color="black",
                    nudge_y = 0.5,
                    nudge_x = 0.5, 
                    min.segment.length = unit(0.1, "lines"), 
                    max.overlaps = 40, size = 5) +
    scale_colour_manual(values = colors) }
# -> Plasma
signed_subgroups_ALS_CTR_and_mimic_CTR = merge(DE_subgroups_plasma_ALS_CTR %>%
                                                mutate(signed_ALS_CTR = sign(log2FC_ALS_CTR)*log10_padj) %>%
                                                rename(DE_ALS_CTR = DE) %>%
                                                select(Target,UniProtID,signed_ALS_CTR,DE_ALS_CTR),
                                              DE_subgroups_plasma_mimic_CTR %>%
                                                mutate(signed_mimic_CTR = sign(log2FC_mimic_CTR)*log10_padj) %>%
                                                rename(DE_mimic_CTR = DE) %>%
                                                select(Target,UniProtID,signed_mimic_CTR,DE_mimic_CTR)
) %>%
  mutate(DE_new = ifelse(DE_ALS_CTR == "CTR" & DE_mimic_CTR == "CTR","CTR sign. both",
                         ifelse(DE_ALS_CTR == "ALS" & DE_mimic_CTR == "mimic","ALS and mimic",
                         ifelse(DE_ALS_CTR == "CTR","CTR",
                                ifelse(DE_ALS_CTR == "ALS","ALS",
                                       ifelse(DE_mimic_CTR == "mimic","mimic",
                                              ifelse(DE_mimic_CTR == "CTR","CTR","ns")))))))

pdf("plots/signed_ALSvsCTR_mimicvsCTR_plasma.pdf")
signed_pvalue_ALS_CTR_mimic(signed_subgroups_ALS_CTR_and_mimic_CTR,"ALS vs CTR + mimic vs CTR (PLASMA)")
dev.off()

# -> CSF
signed_subgroups_ALS_CTR_and_mimic_CTR = merge(DE_subgroups_CSF_ALS_CTR %>%
                                                 mutate(signed_ALS_CTR = sign(log2FC_ALS_CTR)*log10_padj) %>%
                                                 rename(DE_ALS_CTR = DE) %>%
                                                 select(Target,UniProtID,signed_ALS_CTR,DE_ALS_CTR),
                                               DE_subgroups_CSF_mimic_CTR %>%
                                                 mutate(signed_mimic_CTR = sign(log2FC_mimic_CTR)*log10_padj) %>%
                                                 rename(DE_mimic_CTR = DE) %>%
                                                 select(Target,UniProtID,signed_mimic_CTR,DE_mimic_CTR)
) %>%
  mutate(DE_new = ifelse(DE_ALS_CTR == "CTR" & DE_mimic_CTR == "CTR","CTR sign. both",
                         ifelse(DE_ALS_CTR == "ALS" & DE_mimic_CTR == "mimic","ALS and mimic",
                                ifelse(DE_ALS_CTR == "CTR","CTR",
                                       ifelse(DE_ALS_CTR == "ALS","ALS",
                                              ifelse(DE_mimic_CTR == "mimic","mimic",
                                                     ifelse(DE_mimic_CTR == "CTR","CTR","ns")))))))

pdf("plots/signed_ALSvsCTR_mimicvsCTR_CSF.pdf")
signed_pvalue_ALS_CTR_mimic(signed_subgroups_ALS_CTR_and_mimic_CTR,"ALS vs CTR + mimic vs CTR (CSF)")
dev.off()

# -> Serum
signed_subgroups_ALS_CTR_and_mimic_CTR = merge(DE_subgroups_SERUM_ALS_CTR %>%
                                                 mutate(signed_ALS_CTR = sign(log2FC_ALS_CTR)*log10_padj) %>%
                                                 rename(DE_ALS_CTR = DE) %>%
                                                 select(Target,UniProtID,signed_ALS_CTR,DE_ALS_CTR),
                                               DE_subgroups_SERUM_mimic_CTR %>%
                                                 mutate(signed_mimic_CTR = sign(log2FC_mimic_CTR)*log10_padj) %>%
                                                 rename(DE_mimic_CTR = DE) %>%
                                                 select(Target,UniProtID,signed_mimic_CTR,DE_mimic_CTR)
) %>%
  mutate(DE_new = ifelse(DE_ALS_CTR == "CTR" & DE_mimic_CTR == "CTR","CTR sign. both",
                         ifelse(DE_ALS_CTR == "ALS" & DE_mimic_CTR == "mimic","ALS and mimic",
                                ifelse(DE_ALS_CTR == "CTR","CTR",
                                       ifelse(DE_ALS_CTR == "ALS","ALS",
                                              ifelse(DE_mimic_CTR == "mimic","mimic",
                                                     ifelse(DE_mimic_CTR == "CTR","CTR","ns")))))))

pdf("plots/signed_ALSvsCTR_mimicvsCTR_SERUM.pdf")
signed_pvalue_ALS_CTR_mimic(signed_subgroups_ALS_CTR_and_mimic_CTR,"ALS vs CTR + mimic vs CTR (SERUM)")
dev.off()

# C9orf72 vs CTR + SOD1 vs CTR
colors = c('CTR' = '#32A55E',
           'C9orf72' = '#96AA9A',
           'CTR sign. both' = '#5E718B',
           'ns' = '#B4BFC5',
           'SOD1' = '#CF7041',
           "C9orf72 and SOD1" = "#E5C5BD")
signed_pvalue_C9orf72_SOD1_CTR <- function(data_to_plot,title_plot){
  ggplot(data_to_plot,aes(signed_C9orf72_CTR,signed_SOD1_CTR)) +
    geom_point(aes(colour = DE_new),size = 2.25, color = "grey", alpha = 0.5) +
    geom_point(aes(colour = DE_new),size = 2.25, alpha = 0.5) +
    theme(panel.background = element_blank(),
          plot.background = element_blank(), 
          axis.line = element_line()) +
    ggtitle(title_plot)+
    theme(plot.title = element_text(size = 15, face = "bold"))+
    scale_x_continuous(breaks = seq(-50, 20, by = 5))+
    scale_y_continuous(breaks = seq(-30, 140, by = 10))+
    xlab(paste0("C9orf72 vs CTR (sign(LFC) x -log10(p-adj))"))+
    theme(axis.title = element_text(size = 15), axis.text=element_text(size=10))+
    ylab(paste0("SOD1 vs CTR (sign(LFC) x -log10(p-adj))"))+
    theme(text = element_text(size = 15),               # Base font size for everything
          axis.title = element_text(size = 18),         # Axis titles
          axis.text = element_text(size = 16),          # Axis tick labels
          plot.title = element_text(size = 18, hjust = 0.5,face = "bold"),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 15))+
    geom_vline(xintercept = log10(0.05), linetype="dashed",color = "grey",alpha = 0.6)+
    geom_vline(xintercept = -log10(0.05), linetype="dashed",color = "grey",alpha = 0.6)+
    geom_hline(yintercept = log10(0.05), linetype="dashed",color = "grey",alpha = 0.6)+
    geom_hline(yintercept = -log10(0.05), linetype="dashed",color = "grey",alpha = 0.6)+
    #geom_abline(linetype="dashed",color = "grey",alpha = 0.6) +
    annotate(geom = "text", label = "FDR = 5%",x= -log10(0.05)+0.9,y= log10(0.05)-0.15,
             color = "grey", size = 5) +
    annotate(geom = "text", label = "FDR = 5%",x=log10(0.05)-0.9,y= - log10(0.05)+0.15,
             color = "grey", size = 5) +
    geom_text_repel(data = subset(data_to_plot, DE_new != "ns"), 
                    aes(label=Target), color="black",
                    nudge_y = 0.5,
                    nudge_x = 0.5,  
                    min.segment.length = unit(0.1, "lines"), 
                    max.overlaps = 40, size = 5) +
    scale_colour_manual(values = colors) }
# -> Plasma
signed_subgroups_C9orf72_CTR_and_SOD1_CTR = merge(DE_subgroups_PLASMA_C9orf72_CTR %>%
                                                mutate(signed_C9orf72_CTR = sign(log2FC_C9orf72_CTR)*log10_padj) %>%
                                                rename(DE_C9orf72_CTR = DE) %>%
                                                select(Target,UniProtID,signed_C9orf72_CTR,DE_C9orf72_CTR),
                                              DE_subgroups_PLASMA_SOD1_CTR %>%
                                                mutate(signed_SOD1_CTR = sign(log2FC_SOD1_CTR)*log10_padj) %>%
                                                rename(DE_SOD1_CTR = DE) %>%
                                                select(Target,UniProtID,signed_SOD1_CTR,DE_SOD1_CTR)
) %>%
  mutate(DE_new = ifelse(DE_C9orf72_CTR == "CTR" & DE_SOD1_CTR == "CTR","CTR sign. both",
                         ifelse(DE_C9orf72_CTR == "C9orf72" & DE_SOD1_CTR == "SOD1","C9orf72 and SOD1",
                                ifelse(DE_C9orf72_CTR == "CTR","CTR",
                                       ifelse(DE_C9orf72_CTR == "C9orf72","C9orf72",
                                              ifelse(DE_SOD1_CTR == "SOD1","SOD1",
                                                     ifelse(DE_SOD1_CTR == "CTR","CTR","ns")))))))

pdf("plots/signed_C9orf72vsCTR_SOD1vsCTR_plasma.pdf")
signed_pvalue_C9orf72_SOD1_CTR(signed_subgroups_C9orf72_CTR_and_SOD1_CTR,"C9orf72 vs CTR + SOD1 vs CTR (PLASMA)")
dev.off()

# -> CSF
signed_subgroups_C9orf72_CTR_and_SOD1_CTR = merge(DE_subgroups_CSF_C9orf72_CTR %>%
                                                    mutate(signed_C9orf72_CTR = sign(log2FC_C9orf72_CTR)*log10_padj) %>%
                                                    rename(DE_C9orf72_CTR = DE) %>%
                                                    select(Target,UniProtID,signed_C9orf72_CTR,DE_C9orf72_CTR),
                                                  DE_subgroups_CSF_SOD1_CTR %>%
                                                    mutate(signed_SOD1_CTR = sign(log2FC_SOD1_CTR)*log10_padj) %>%
                                                    rename(DE_SOD1_CTR = DE) %>%
                                                    select(Target,UniProtID,signed_SOD1_CTR,DE_SOD1_CTR)
) %>%
  mutate(DE_new = ifelse(DE_C9orf72_CTR == "CTR" & DE_SOD1_CTR == "CTR","CTR sign. both",
                         ifelse(DE_C9orf72_CTR == "C9orf72" & DE_SOD1_CTR == "SOD1","C9orf72 and SOD1",
                                ifelse(DE_C9orf72_CTR == "CTR","CTR",
                                       ifelse(DE_C9orf72_CTR == "C9orf72","C9orf72",
                                              ifelse(DE_SOD1_CTR == "SOD1","SOD1",
                                                     ifelse(DE_SOD1_CTR == "CTR","CTR","ns")))))))

pdf("plots/signed_C9orf72vsCTR_SOD1vsCTR_CSF.pdf")
signed_pvalue_C9orf72_SOD1_CTR(signed_subgroups_C9orf72_CTR_and_SOD1_CTR,"C9orf72 vs CTR + SOD1 vs CTR (CSF)")
dev.off()

# -> SERUM
signed_subgroups_C9orf72_CTR_and_SOD1_CTR = merge(DE_subgroups_SERUM_C9orf72_CTR %>%
                                                    mutate(signed_C9orf72_CTR = sign(log2FC_C9orf72_CTR)*log10_padj) %>%
                                                    rename(DE_C9orf72_CTR = DE) %>%
                                                    select(Target,UniProtID,signed_C9orf72_CTR,DE_C9orf72_CTR),
                                                  DE_subgroups_SERUM_SOD1_CTR %>%
                                                    mutate(signed_SOD1_CTR = sign(log2FC_SOD1_CTR)*log10_padj) %>%
                                                    rename(DE_SOD1_CTR = DE) %>%
                                                    select(Target,UniProtID,signed_SOD1_CTR,DE_SOD1_CTR)
) %>%
  mutate(DE_new = ifelse(DE_C9orf72_CTR == "CTR" & DE_SOD1_CTR == "CTR","CTR sign. both",
                         ifelse(DE_C9orf72_CTR == "C9orf72" & DE_SOD1_CTR == "SOD1","C9orf72 and SOD1",
                                ifelse(DE_C9orf72_CTR == "CTR","CTR",
                                       ifelse(DE_C9orf72_CTR == "C9orf72","C9orf72",
                                              ifelse(DE_SOD1_CTR == "SOD1","SOD1",
                                                     ifelse(DE_SOD1_CTR == "CTR","CTR","ns")))))))

pdf("plots/signed_C9orf72vsCTR_SOD1vsCTR_SERUM.pdf")
signed_pvalue_C9orf72_SOD1_CTR(signed_subgroups_C9orf72_CTR_and_SOD1_CTR,"C9orf72 vs CTR + SOD1 vs CTR (SERUM)")
dev.off()


