# Differential expression analyses

###############################################
### Helper Functions
###############################################

## 1. Mean NPQ
summarise_mean_NPQ <- function(data, matrix_type,adjusted = FALSE) {
  val_col <- if(adjusted) "NPQ_adj" else "NPQ"
  
  data %>%
    filter(SampleMatrixType == matrix_type, !is.na(.data[[val_col]])) %>%
    group_by(type, Target, UniProtID) %>%
    summarise(mean_NPQ = mean(.data[[val_col]]), .groups = "drop") %>%
    pivot_wider(
      names_from = type,
      values_from = mean_NPQ,
      names_prefix = "mean_NPQ_"
    )
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
      log2FC_others_CTR = mean_NPQ_others - mean_NPQ_CTR,
      log2FC_C9orf72_SOD1 = mean_NPQ_C9orf72 - mean_NPQ_SOD1,
      log2FC_C9orf72_TARDBP = mean_NPQ_C9orf72 - mean_NPQ_TARDBP,
      log2FC_C9orf72_others = mean_NPQ_C9orf72 - mean_NPQ_others,
      log2FC_SOD1_TARDBP = mean_NPQ_SOD1 - mean_NPQ_TARDBP,
      log2FC_SOD1_others = mean_NPQ_SOD1 - mean_NPQ_others,
      log2FC_TARDBP_others = mean_NPQ_TARDBP - mean_NPQ_others)
}

## 3. Final table
finalise_DE_table <- function(df,PGMC_groups = FALSE) {
  if(PGMC_groups){
    df %>%
      select(
        Target, UniProtID, Fluid,pvalue_anova,
        contains("C9orf72_CTR"),contains("CTR_C9orf72"),
        contains("TARDBP_CTR"),contains("CTR_TARDBP"),
        contains("SOD1_CTR"),contains("CTR_SOD1"),
        contains("others_CTR"),contains("CTR_others"),
        contains("C9orf72_SOD1"),contains("C9orf72_TARDBP"),
        contains("C9orf72_others"),contains("SOD1_TARDBP"),
        contains("others_SOD1"),contains("SOD1_others"),
        contains("others_TARDBP"),contains("TARDBP_others"))  %>%
      rename_with(~ gsub("CTR_C9orf72", "C9orf72_CTR", .x), .cols = contains("CTR_C9orf72")) %>%
      rename_with(~ gsub("CTR_TARDBP", "TARDBP_CTR", .x), .cols = contains("CTR_TARDBP")) %>%
      rename_with(~ gsub("CTR_SOD1", "SOD1_CTR", .x), .cols = contains("CTR_SOD1")) %>%
      rename_with(~ gsub("CTR_others","others_CTR", .x), .cols = contains("CTR_others")) %>%
      rename_with(~ gsub("others_SOD1", "SOD1_others", .x), .cols = contains("others_SOD1")) %>%
      rename_with(~ gsub("others_TARDBP", "TARDBP_others", .x), .cols = contains("others_TARDBP")) 
  }else{
    df %>%
      select(
        Target, UniProtID, Fluid,pvalue_anova,
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
    label_targets = NULL
) {
  
  x_max <- max(df[[log2fc]], na.rm = TRUE)
  x_min <- min(df[[log2fc]], na.rm = TRUE)
  y_max <- max(df$log10_padj, na.rm = TRUE)
  
  # define x and y positions for the significance label
  label_x_pos <- x_min + (x_max - x_min) * 0.05  # 5% from the left edge
  label_y_pos <- -log10(fdr) - (y_max * 0.02)    # Just slightly beloew the FDR line
  
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
      x = label_x_pos,
      y = label_y_pos,
      label = paste0("FDR = ", fdr * 100, "%"),
      color = "darkgrey",
      size = 5, hjust = 0) +
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

volcano_plot_pvalue <- function(
    df, log2fc, title,
    fdr = 0.05,
    colors,
    label_targets = NULL,
    target_color = "#B22222"
) {
  
  # Calculate auto-coordinates for the significance label
  x_max <- max(df[[log2fc]], na.rm = TRUE)
  x_min <- min(df[[log2fc]], na.rm = TRUE)
  y_max <- max(df$log10_padj, na.rm = TRUE)
  
  # define x and y positions for the significance label
  label_x_pos <- x_min + (x_max - x_min) * 0.05  # 5% from the left edge
  label_y_pos <- -log10(fdr) - (y_max * 0.02)    # Just slightly beloew the FDR line
  
  ggplot(df, aes_string(log2fc, "log10_padj")) +
    geom_point(aes(color = DE, size = log10_padj), alpha = 0.6) +
    scale_color_manual(
      name = "P-value < 0.05", 
      values = colors,
      guide = guide_legend(order = 1)) +
    
    # Switch to new scale
    new_scale_color() + 
    
    # Highlight specific points
    geom_point(
      data = df %>% filter(Target %in% label_targets),
      aes(size = log10_padj, color = "Sign. Proteins"), 
      shape = 19) +
    scale_color_manual(
      name = "ANOVA & P-value < 0.05", 
      values = c("Sign. Proteins" = target_color),
      guide = guide_legend(order = 2)) +
    
    # Label specific points
    geom_text_repel(
      data = df %>% filter(Target %in% label_targets),
      aes(label = Target, size = 4),
      color = target_color,
      max.overlaps = 40,
      box.padding = 0.5,
      show.legend = FALSE) +
    geom_hline(yintercept = -log10(fdr), linetype = "dashed", color = "darkgrey") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey") +
    annotate("text", x = label_x_pos, y = label_y_pos,
             label = paste0("P-value = ", fdr * 100, "%"),
             color = "darkgrey", size = 5,hjust = 0) +
    theme_minimal() +
    ggtitle(title) +
    scale_size_continuous(range = c(2, 7), guide = guide_legend(order = 3)) + 
    ylab(expression("-log"[10]*"(p-value)")) +
    xlab(expression("log"[2]*"(fold-change)")) + 
    labs(size = expression("-log"[10]*"(p-value)")) +
    theme(
      axis.title.x = element_blank(),
      text = element_text(size = 15),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 16),
      plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 14),
      legend.spacing.y = unit(0.3, 'cm'))
}


###############################################################################
# Run pipeline
###############################################################################

###### Differential expression datasets for each fluid
tissues = c("PLASMA", "CSF", "SERUM")

# Initialize a list to store all final dataframes
DE_results <- list()

for (tissue in tissues) {
  message(paste("Processing tissue:", tissue))
  DE_results[[tissue]] <- list()
  
  # --- 1. Standard Subgroups (Unadjusted) ---
  DE_results[[tissue]]$standard <- summarise_mean_NPQ(protein_data_IDs, matrix_type = tissue) %>%
    add_log2fc_standard() %>%
    left_join(results_ALL[[tissue]]$pvals) %>%
    mutate(Fluid = tissue) %>%
    finalise_DE_table()
  
  # --- 2. Standard Subgroups (Adjusted) ---
  DE_results[[tissue]]$standard_adj <- summarise_mean_NPQ(results_ALL[[tissue]]$data_adjusted, 
                                                          matrix_type = tissue, adjusted = TRUE) %>%
    add_log2fc_standard() %>%
    left_join(results_ALL[[tissue]]$pvals_adjusted) %>%
    mutate(Fluid = tissue) %>%
    finalise_DE_table()
  
  # --- 3. PGMC Subgroups (Unadjusted) ---
  DE_results[[tissue]]$pgmc <- summarise_mean_NPQ(results_PGMC[[tissue]]$data, matrix_type = tissue) %>%
    add_log2fc_pgmc() %>%
    left_join(results_PGMC[[tissue]]$pvals) %>%
    mutate(Fluid = tissue) %>%
    finalise_DE_table(PGMC_groups = TRUE)
  
  # --- 4. PGMC Subgroups (Adjusted) ---
  DE_results[[tissue]]$pgmc_adj <- summarise_mean_NPQ(results_PGMC[[tissue]]$data_adjusted, 
                                                      matrix_type = tissue, adjusted = TRUE) %>%
    add_log2fc_pgmc() %>%
    left_join(results_PGMC[[tissue]]$pvals_adjusted) %>%
    mutate(Fluid = tissue) %>%
    finalise_DE_table(PGMC_groups = TRUE)
  
  # --- Export Files ---
  writexl::write_xlsx(DE_results[[tissue]]$standard,     paste0("results/DE_subgroups_", tissue, ".xlsx"))
  writexl::write_xlsx(DE_results[[tissue]]$standard_adj, paste0("results/DE_subgroups_", tissue, "_adjusted.xlsx"))
  writexl::write_xlsx(DE_results[[tissue]]$pgmc,         paste0("results/DE_PGMC_subgroups_", tissue, ".xlsx"))
  writexl::write_xlsx(DE_results[[tissue]]$pgmc_adj,     paste0("results/DE_PGMC_subgroups_", tissue, "_adjusted.xlsx"))
}


###### VOLCANO PLOTS 
## 1: ALS vs CTR 
volcano_plots <- list()
colors_v <- c('CTR' = "#C25F3D", 'ns' = 'lightgrey', 'ALS'= "#3da0c2")

for (tissue in tissues) {
  message(paste("Generating volcano plots for:", tissue))
  
  # --- 1. Standard (Not Adjusted) ---
  # Pulling from our structured list: DE_results[[tissue]]$standard
  df_std <- DE_results[[tissue]]$standard %>%
    select(Target, UniProtID, log2FC_ALS_CTR, pvalue_ALS_CTR, padj_ALS_CTR) %>%
    add_DE_flag(lfc = "log2FC_ALS_CTR", padj = "padj_ALS_CTR", up = "ALS", down = "CTR", alpha = 0.1)
  
  volcano_plots[[tissue]]$std <- volcano_plot(
    df_std, "log2FC_ALS_CTR", fdr = 0.1,
    title = paste("ALS vs CTR in", tissue),
    colors = colors_v)
  
  # --- 2. Standard (Adjusted) ---
  df_adj <- DE_results[[tissue]]$standard_adj %>%
    select(Target, UniProtID, log2FC_ALS_CTR, pvalue_ALS_CTR, padj_ALS_CTR) %>%
    add_DE_flag(lfc = "log2FC_ALS_CTR", padj = "padj_ALS_CTR", up = "ALS", down = "CTR", alpha = 0.1)
  
  volcano_plots[[tissue]]$adj <- volcano_plot(
    df_adj, "log2FC_ALS_CTR", fdr = 0.1,
    title = paste("ALS vs CTR in", tissue, "(adjusted)"),
    colors = colors_v
  )
  
  # --- 3. Adjusted + P-Value (with special labels) ---
  df_pv <- DE_results[[tissue]]$standard_adj %>%
    select(Target, UniProtID, log2FC_ALS_CTR, pvalue_ALS_CTR, padj_ALS_CTR, pvalue_anova) %>%
    add_DE_flag(lfc = "log2FC_ALS_CTR", padj = "pvalue_ALS_CTR", up = "ALS", down = "CTR", alpha = 0.05)
  
  # Auto-calculate targets based on criteria
  targets <- df_pv %>% 
    filter(pvalue_anova < 0.05 & pvalue_ALS_CTR < 0.05) %>% 
    pull(Target)
  
  volcano_plots[[tissue]]$pv_adj <- volcano_plot_pvalue(
    df_pv, "log2FC_ALS_CTR", 
    title = paste("ALS vs CTR in", tissue, "(adj + p-value)"),
    colors = colors_v,
    label_targets = targets,
    target_color = "#702963")
}

# PDF 1: Standard vs Adjusted (FDR 10%)
pdf("plots/volcano_plots/volcano_ALS_CTR_all_FDR10.pdf", width = 14, height = 15)
grid.arrange(
  volcano_plots$PLASMA$std, volcano_plots$PLASMA$adj,
  volcano_plots$SERUM$std,  volcano_plots$SERUM$adj,
  volcano_plots$CSF$std,    volcano_plots$CSF$adj,
  nrow = 3, ncol = 2)
dev.off()

# PDF 2: P-value Highlight vs Adjusted
pdf("plots/volcano_plots/volcano_ALS_CTR_all_pvalue5_FDR10.pdf", width = 14, height = 15)
grid.arrange(
  volcano_plots$PLASMA$pv_adj, volcano_plots$PLASMA$adj,
  volcano_plots$SERUM$pv_adj,  volcano_plots$SERUM$adj,
  volcano_plots$CSF$pv_adj,    volcano_plots$CSF$adj,
  nrow = 3, ncol = 2)
dev.off()
volcano_plots_ALS_CTR = volcano_plots

## 2: ALS vs PGMC 
volcano_plots <- list()
colors_v <- c('PGMC' = "#C25F3D", 'ns' = 'lightgrey', 'ALS'= "#3da0c2")

for (tissue in tissues) {
  message(paste("Generating volcano plots for:", tissue))
  
  # --- 1. Standard (Not Adjusted) ---
  # Pulling from our structured list: DE_results[[tissue]]$standard
  df_std <- DE_results[[tissue]]$standard %>%
    select(Target, UniProtID, log2FC_ALS_PGMC, pvalue_ALS_PGMC, padj_ALS_PGMC) %>%
    add_DE_flag(lfc = "log2FC_ALS_PGMC", padj = "padj_ALS_PGMC", up = "ALS", down = "PGMC", alpha = 0.1)
  
  volcano_plots[[tissue]]$std <- volcano_plot(
    df_std, "log2FC_ALS_PGMC", fdr = 0.1,
    title = paste("ALS vs PGMC in", tissue),
    colors = colors_v)
  
  # --- 2. Standard (Adjusted) ---
  df_adj <- DE_results[[tissue]]$standard_adj %>%
    select(Target, UniProtID, log2FC_ALS_PGMC, pvalue_ALS_PGMC, padj_ALS_PGMC) %>%
    add_DE_flag(lfc = "log2FC_ALS_PGMC", padj = "padj_ALS_PGMC", up = "ALS", down = "PGMC", alpha = 0.1)
  
  volcano_plots[[tissue]]$adj <- volcano_plot(
    df_adj, "log2FC_ALS_PGMC", fdr = 0.1,
    title = paste("ALS vs PGMC in", tissue, "(adjusted)"),
    colors = colors_v
  )
  
  # --- 3. Adjusted + P-Value (with special labels) ---
  df_pv <- DE_results[[tissue]]$standard_adj %>%
    select(Target, UniProtID, log2FC_ALS_PGMC, pvalue_ALS_PGMC, padj_ALS_PGMC, pvalue_anova) %>%
    add_DE_flag(lfc = "log2FC_ALS_PGMC", padj = "pvalue_ALS_PGMC", up = "ALS", down = "PGMC", alpha = 0.05)
  
  # Auto-calculate targets based on criteria
  targets <- df_pv %>% 
    filter(pvalue_anova < 0.05 & pvalue_ALS_PGMC < 0.05) %>% 
    pull(Target)
  
  volcano_plots[[tissue]]$pv_adj <- volcano_plot_pvalue(
    df_pv, "log2FC_ALS_PGMC", 
    title = paste("ALS vs PGMC in", tissue, "(adj + p-value)"),
    colors = colors_v,
    label_targets = targets,
    target_color = "#702963")
}

# PDF 1: Standard vs Adjusted (FDR 10%)
pdf("plots/volcano_plots/volcano_ALS_PGMC_all_FDR10.pdf", width = 14, height = 15)
grid.arrange(
  volcano_plots$PLASMA$std, volcano_plots$PLASMA$adj,
  volcano_plots$SERUM$std,  volcano_plots$SERUM$adj,
  volcano_plots$CSF$std,    volcano_plots$CSF$adj,
  nrow = 3, ncol = 2)
dev.off()

# PDF 2: P-value Highlight vs Adjusted
pdf("plots/volcano_plots/volcano_ALS_PGMC_all_pvalue5_FDR10.pdf", width = 14, height = 15)
grid.arrange(
  volcano_plots$PLASMA$pv_adj, volcano_plots$PLASMA$adj,
  volcano_plots$SERUM$pv_adj,  volcano_plots$SERUM$adj,
  volcano_plots$CSF$pv_adj,    volcano_plots$CSF$adj,
  nrow = 3, ncol = 2)
dev.off()

volcano_plots_ALS_PGMC = volcano_plots

## 3: ALS vs mimic 
volcano_plots <- list()
colors_v <- c('mimic' = "#C25F3D", 'ns' = 'lightgrey', 'ALS'= "#3da0c2")

for (tissue in tissues) {
  message(paste("Generating volcano plots for:", tissue))
  
  # --- 1. Standard (Not Adjusted) ---
  # Pulling from our structured list: DE_results[[tissue]]$standard
  df_std <- DE_results[[tissue]]$standard %>%
    select(Target, UniProtID, log2FC_ALS_mimic, pvalue_ALS_mimic, padj_ALS_mimic) %>%
    add_DE_flag(lfc = "log2FC_ALS_mimic", padj = "padj_ALS_mimic", up = "ALS", down = "mimic", alpha = 0.1)
  
  volcano_plots[[tissue]]$std <- volcano_plot(
    df_std, "log2FC_ALS_mimic", fdr = 0.1,
    title = paste("ALS vs mimic in", tissue),
    colors = colors_v)
  
  # --- 2. Standard (Adjusted) ---
  df_adj <- DE_results[[tissue]]$standard_adj %>%
    select(Target, UniProtID, log2FC_ALS_mimic, pvalue_ALS_mimic, padj_ALS_mimic) %>%
    add_DE_flag(lfc = "log2FC_ALS_mimic", padj = "padj_ALS_mimic", up = "ALS", down = "mimic", alpha = 0.1)
  
  volcano_plots[[tissue]]$adj <- volcano_plot(
    df_adj, "log2FC_ALS_mimic", fdr = 0.1,
    title = paste("ALS vs mimic in", tissue, "(adjusted)"),
    colors = colors_v
  )
  
  # --- 3. Adjusted + P-Value (with special labels) ---
  df_pv <- DE_results[[tissue]]$standard_adj %>%
    select(Target, UniProtID, log2FC_ALS_mimic, pvalue_ALS_mimic, padj_ALS_mimic, pvalue_anova) %>%
    add_DE_flag(lfc = "log2FC_ALS_mimic", padj = "pvalue_ALS_mimic", up = "ALS", down = "mimic", alpha = 0.05)
  
  # Auto-calculate targets based on criteria
  targets <- df_pv %>% 
    filter(pvalue_anova < 0.05 & pvalue_ALS_mimic < 0.05) %>% 
    pull(Target)
  
  volcano_plots[[tissue]]$pv_adj <- volcano_plot_pvalue(
    df_pv, "log2FC_ALS_mimic", 
    title = paste("ALS vs mimic in", tissue, "(adj + p-value)"),
    colors = colors_v,
    label_targets = targets,
    target_color = "#702963")
}

# PDF 1: Standard vs Adjusted (FDR 10%)
pdf("plots/volcano_plots/volcano_ALS_mimic_all_FDR10.pdf", width = 14, height = 15)
grid.arrange(
  volcano_plots$PLASMA$std, volcano_plots$PLASMA$adj,
  volcano_plots$SERUM$std,  volcano_plots$SERUM$adj,
  volcano_plots$CSF$std,    volcano_plots$CSF$adj,
  nrow = 3, ncol = 2)
dev.off()

# PDF 2: P-value Highlight vs Adjusted
pdf("plots/volcano_plots/volcano_ALS_mimic_all_pvalue5_FDR10.pdf", width = 14, height = 15)
grid.arrange(
  volcano_plots$PLASMA$pv_adj, volcano_plots$PLASMA$adj,
  volcano_plots$SERUM$pv_adj,  volcano_plots$SERUM$adj,
  volcano_plots$CSF$pv_adj,    volcano_plots$CSF$adj,
  nrow = 3, ncol = 2)
dev.off()

volcano_plots_ALS_mimic = volcano_plots

## 4: mimic vs CTR 
volcano_plots <- list()
colors_v <- c('CTR' = "#C25F3D", 'ns' = 'lightgrey', 'mimic'= "#3da0c2")

for (tissue in tissues) {
  message(paste("Generating volcano plots for:", tissue))
  
  # --- 1. Standard (Not Adjusted) ---
  # Pulling from our structured list: DE_results[[tissue]]$standard
  df_std <- DE_results[[tissue]]$standard %>%
    select(Target, UniProtID, log2FC_mimic_CTR, pvalue_mimic_CTR, padj_mimic_CTR) %>%
    add_DE_flag(lfc = "log2FC_mimic_CTR", padj = "padj_mimic_CTR", up = "mimic", down = "CTR", alpha = 0.1)
  
  volcano_plots[[tissue]]$std <- volcano_plot(
    df_std, "log2FC_mimic_CTR", fdr = 0.1,
    title = paste("mimic vs CTR in", tissue),
    colors = colors_v)
  
  # --- 2. Standard (Adjusted) ---
  df_adj <- DE_results[[tissue]]$standard_adj %>%
    select(Target, UniProtID, log2FC_mimic_CTR, pvalue_mimic_CTR, padj_mimic_CTR) %>%
    add_DE_flag(lfc = "log2FC_mimic_CTR", padj = "padj_mimic_CTR", up = "mimic", down = "CTR", alpha = 0.1)
  
  volcano_plots[[tissue]]$adj <- volcano_plot(
    df_adj, "log2FC_mimic_CTR", fdr = 0.1,
    title = paste("mimic vs CTR in", tissue, "(adjusted)"),
    colors = colors_v
  )
  
  # --- 3. Adjusted + P-Value (with special labels) ---
  df_pv <- DE_results[[tissue]]$standard_adj %>%
    select(Target, UniProtID, log2FC_mimic_CTR, pvalue_mimic_CTR, padj_mimic_CTR, pvalue_anova) %>%
    add_DE_flag(lfc = "log2FC_mimic_CTR", padj = "pvalue_mimic_CTR", up = "mimic", down = "CTR", alpha = 0.05)
  
  # Auto-calculate targets based on criteria
  targets <- df_pv %>% 
    filter(pvalue_anova < 0.05 & pvalue_mimic_CTR < 0.05) %>% 
    pull(Target)
  
  volcano_plots[[tissue]]$pv_adj <- volcano_plot_pvalue(
    df_pv, "log2FC_mimic_CTR", 
    title = paste("mimic vs CTR in", tissue, "(adj + p-value)"),
    colors = colors_v,
    label_targets = targets,
    target_color = "#702963")
}

# PDF 1: Standard vs Adjusted (FDR 10%)
pdf("plots/volcano_plots/volcano_mimic_CTR_all_FDR10.pdf", width = 14, height = 15)
grid.arrange(
  volcano_plots$PLASMA$std, volcano_plots$PLASMA$adj,
  volcano_plots$SERUM$std,  volcano_plots$SERUM$adj,
  volcano_plots$CSF$std,    volcano_plots$CSF$adj,
  nrow = 3, ncol = 2)
dev.off()

# PDF 2: P-value Highlight vs Adjusted
pdf("plots/volcano_plots/volcano_mimic_CTR_all_pvalue5_FDR10.pdf", width = 14, height = 15)
grid.arrange(
  volcano_plots$PLASMA$pv_adj, volcano_plots$PLASMA$adj,
  volcano_plots$SERUM$pv_adj,  volcano_plots$SERUM$adj,
  volcano_plots$CSF$pv_adj,    volcano_plots$CSF$adj,
  nrow = 3, ncol = 2)
dev.off()

volcano_plots_mimic_CTR = volcano_plots

## 5: PGMC vs CTR 
volcano_plots <- list()
colors_v <- c('CTR' = "#C25F3D", 'ns' = 'lightgrey', 'PGMC'= "#3da0c2")

for (tissue in tissues) {
  message(paste("Generating volcano plots for:", tissue))
  
  # --- 1. Standard (Not Adjusted) ---
  # Pulling from our structured list: DE_results[[tissue]]$standard
  df_std <- DE_results[[tissue]]$standard %>%
    select(Target, UniProtID, log2FC_PGMC_CTR, pvalue_PGMC_CTR, padj_PGMC_CTR) %>%
    add_DE_flag(lfc = "log2FC_PGMC_CTR", padj = "padj_PGMC_CTR", up = "PGMC", down = "CTR", alpha = 0.1)
  
  volcano_plots[[tissue]]$std <- volcano_plot(
    df_std, "log2FC_PGMC_CTR", fdr = 0.1,
    title = paste("PGMC vs CTR in", tissue),
    colors = colors_v)
  
  # --- 2. Standard (Adjusted) ---
  df_adj <- DE_results[[tissue]]$standard_adj %>%
    select(Target, UniProtID, log2FC_PGMC_CTR, pvalue_PGMC_CTR, padj_PGMC_CTR) %>%
    add_DE_flag(lfc = "log2FC_PGMC_CTR", padj = "padj_PGMC_CTR", up = "PGMC", down = "CTR", alpha = 0.1)
  
  volcano_plots[[tissue]]$adj <- volcano_plot(
    df_adj, "log2FC_PGMC_CTR", fdr = 0.1,
    title = paste("PGMC vs CTR in", tissue, "(adjusted)"),
    colors = colors_v
  )
  
  # --- 3. Adjusted + P-Value (with special labels) ---
  df_pv <- DE_results[[tissue]]$standard_adj %>%
    select(Target, UniProtID, log2FC_PGMC_CTR, pvalue_PGMC_CTR, padj_PGMC_CTR, pvalue_anova) %>%
    add_DE_flag(lfc = "log2FC_PGMC_CTR", padj = "pvalue_PGMC_CTR", up = "PGMC", down = "CTR", alpha = 0.05)
  
  # Auto-calculate targets based on criteria
  targets <- df_pv %>% 
    filter(pvalue_anova < 0.05 & pvalue_PGMC_CTR < 0.05) %>% 
    pull(Target)
  
  volcano_plots[[tissue]]$pv_adj <- volcano_plot_pvalue(
    df_pv, "log2FC_PGMC_CTR", 
    title = paste("PGMC vs CTR in", tissue, "(adj + p-value)"),
    colors = colors_v,
    label_targets = targets,
    target_color = "#702963")
}

# PDF 1: Standard vs Adjusted (FDR 10%)
pdf("plots/volcano_plots/volcano_PGMC_CTR_all_FDR10.pdf", width = 14, height = 15)
grid.arrange(
  volcano_plots$PLASMA$std, volcano_plots$PLASMA$adj,
  volcano_plots$SERUM$std,  volcano_plots$SERUM$adj,
  volcano_plots$CSF$std,    volcano_plots$CSF$adj,
  nrow = 3, ncol = 2)
dev.off()

# PDF 2: P-value Highlight vs Adjusted
pdf("plots/volcano_plots/volcano_PGMC_CTR_all_pvalue5_FDR10.pdf", width = 14, height = 15)
grid.arrange(
  volcano_plots$PLASMA$pv_adj, volcano_plots$PLASMA$adj,
  volcano_plots$SERUM$pv_adj,  volcano_plots$SERUM$adj,
  volcano_plots$CSF$pv_adj,    volcano_plots$CSF$adj,
  nrow = 3, ncol = 2)
dev.off()

volcano_plots_PGMC_CTR = volcano_plots

## 6: PGMC vs mimic 
volcano_plots <- list()
colors_v <- c('mimic' = "#C25F3D", 'ns' = 'lightgrey', 'PGMC'= "#3da0c2")

for (tissue in tissues) {
  message(paste("Generating volcano plots for:", tissue))
  
  # --- 1. Standard (Not Adjusted) ---
  # Pulling from our structured list: DE_results[[tissue]]$standard
  df_std <- DE_results[[tissue]]$standard %>%
    select(Target, UniProtID, log2FC_PGMC_mimic, pvalue_PGMC_mimic, padj_PGMC_mimic) %>%
    add_DE_flag(lfc = "log2FC_PGMC_mimic", padj = "padj_PGMC_mimic", up = "PGMC", down = "mimic", alpha = 0.1)
  
  volcano_plots[[tissue]]$std <- volcano_plot(
    df_std, "log2FC_PGMC_mimic", fdr = 0.1,
    title = paste("PGMC vs mimic in", tissue),
    colors = colors_v)
  
  # --- 2. Standard (Adjusted) ---
  df_adj <- DE_results[[tissue]]$standard_adj %>%
    select(Target, UniProtID, log2FC_PGMC_mimic, pvalue_PGMC_mimic, padj_PGMC_mimic) %>%
    add_DE_flag(lfc = "log2FC_PGMC_mimic", padj = "padj_PGMC_mimic", up = "PGMC", down = "mimic", alpha = 0.1)
  
  volcano_plots[[tissue]]$adj <- volcano_plot(
    df_adj, "log2FC_PGMC_mimic", fdr = 0.1,
    title = paste("PGMC vs mimic in", tissue, "(adjusted)"),
    colors = colors_v
  )
  
  # --- 3. Adjusted + P-Value (with special labels) ---
  df_pv <- DE_results[[tissue]]$standard_adj %>%
    select(Target, UniProtID, log2FC_PGMC_mimic, pvalue_PGMC_mimic, padj_PGMC_mimic, pvalue_anova) %>%
    add_DE_flag(lfc = "log2FC_PGMC_mimic", padj = "pvalue_PGMC_mimic", up = "PGMC", down = "mimic", alpha = 0.05)
  
  # Auto-calculate targets based on criteria
  targets <- df_pv %>% 
    filter(pvalue_anova < 0.05 & pvalue_PGMC_mimic < 0.05) %>% 
    pull(Target)
  
  volcano_plots[[tissue]]$pv_adj <- volcano_plot_pvalue(
    df_pv, "log2FC_PGMC_mimic", 
    title = paste("PGMC vs mimic in", tissue, "(adj + p-value)"),
    colors = colors_v,
    label_targets = targets,
    target_color = "#702963")
}

# PDF 1: Standard vs Adjusted (FDR 10%)
pdf("plots/volcano_plots/volcano_PGMC_mimic_all_FDR10.pdf", width = 14, height = 15)
grid.arrange(
  volcano_plots$PLASMA$std, volcano_plots$PLASMA$adj,
  volcano_plots$SERUM$std,  volcano_plots$SERUM$adj,
  volcano_plots$CSF$std,    volcano_plots$CSF$adj,
  nrow = 3, ncol = 2)
dev.off()

# PDF 2: P-value Highlight vs Adjusted
pdf("plots/volcano_plots/volcano_PGMC_mimic_all_pvalue5_FDR10.pdf", width = 14, height = 15)
grid.arrange(
  volcano_plots$PLASMA$pv_adj, volcano_plots$PLASMA$adj,
  volcano_plots$SERUM$pv_adj,  volcano_plots$SERUM$adj,
  volcano_plots$CSF$pv_adj,    volcano_plots$CSF$adj,
  nrow = 3, ncol = 2)
dev.off()

volcano_plots_PGMC_mimic = volcano_plots

## 7: C9orf72 vs CTR 
volcano_plots <- list()
colors_v <- c('CTR' = "#C25F3D", 'ns' = 'lightgrey', 'C9orf72'= "#3da0c2")

for (tissue in tissues) {
  message(paste("Generating volcano plots for:", tissue))
  
  # --- 1. Standard (Not Adjusted) ---
  df_std <- DE_results[[tissue]]$pgmc %>%
    select(Target, UniProtID, log2FC_C9orf72_CTR, pvalue_C9orf72_CTR, padj_C9orf72_CTR) %>%
    add_DE_flag(lfc = "log2FC_C9orf72_CTR", padj = "padj_C9orf72_CTR", up = "C9orf72", down = "CTR", alpha = 0.1)
  
  volcano_plots[[tissue]]$std <- volcano_plot(
    df_std, "log2FC_C9orf72_CTR", fdr = 0.1,
    title = paste("C9orf72 vs CTR in", tissue),
    colors = colors_v)
  
  # --- 2. Standard (Adjusted) ---
  df_adj <- DE_results[[tissue]]$pgmc_adj %>%
    select(Target, UniProtID, log2FC_C9orf72_CTR, pvalue_C9orf72_CTR, padj_C9orf72_CTR) %>%
    add_DE_flag(lfc = "log2FC_C9orf72_CTR", padj = "padj_C9orf72_CTR", up = "C9orf72", down = "CTR", alpha = 0.1)
  
  volcano_plots[[tissue]]$adj <- volcano_plot(
    df_adj, "log2FC_C9orf72_CTR", fdr = 0.1,
    title = paste("C9orf72 vs CTR in", tissue, "(adjusted)"),
    colors = colors_v
  )
  
  # --- 3. Adjusted + P-Value (with special labels) ---
  df_pv <- DE_results[[tissue]]$pgmc_adj %>%
    select(Target, UniProtID, log2FC_C9orf72_CTR, pvalue_C9orf72_CTR, padj_C9orf72_CTR, pvalue_anova) %>%
    add_DE_flag(lfc = "log2FC_C9orf72_CTR", padj = "pvalue_C9orf72_CTR", up = "C9orf72", down = "CTR", alpha = 0.05)
  
  # Auto-calculate targets based on criteria
  targets <- df_pv %>% 
    filter(pvalue_anova < 0.05 & pvalue_C9orf72_CTR < 0.05) %>% 
    pull(Target)
  
  volcano_plots[[tissue]]$pv_adj <- volcano_plot_pvalue(
    df_pv, "log2FC_C9orf72_CTR", 
    title = paste("C9orf72 vs CTR in", tissue, "(adj + p-value)"),
    colors = colors_v,
    label_targets = targets,
    target_color = "#702963")
}

# PDF 1: Standard vs Adjusted (FDR 10%)
pdf("plots/volcano_plots/volcano_C9orf72_CTR_all_FDR10.pdf", width = 14, height = 15)
grid.arrange(
  volcano_plots$PLASMA$std, volcano_plots$PLASMA$adj,
  volcano_plots$SERUM$std,  volcano_plots$SERUM$adj,
  volcano_plots$CSF$std,    volcano_plots$CSF$adj,
  nrow = 3, ncol = 2)
dev.off()

# PDF 2: P-value Highlight vs Adjusted
pdf("plots/volcano_plots/volcano_C9orf72_CTR_all_pvalue5_FDR10.pdf", width = 14, height = 15)
grid.arrange(
  volcano_plots$PLASMA$pv_adj, volcano_plots$PLASMA$adj,
  volcano_plots$SERUM$pv_adj,  volcano_plots$SERUM$adj,
  volcano_plots$CSF$pv_adj,    volcano_plots$CSF$adj,
  nrow = 3, ncol = 2)
dev.off()

volcano_plots_C9orf72_CTR = volcano_plots

## 8: SOD1 vs CTR 
volcano_plots <- list()
colors_v <- c('CTR' = "#C25F3D", 'ns' = 'lightgrey', 'SOD1'= "#3da0c2")

for (tissue in tissues) {
  message(paste("Generating volcano plots for:", tissue))
  
  # --- 1. Standard (Not Adjusted) ---
  df_std <- DE_results[[tissue]]$pgmc %>%
    select(Target, UniProtID, log2FC_SOD1_CTR, pvalue_SOD1_CTR, padj_SOD1_CTR) %>%
    add_DE_flag(lfc = "log2FC_SOD1_CTR", padj = "padj_SOD1_CTR", up = "SOD1", down = "CTR", alpha = 0.1)
  
  volcano_plots[[tissue]]$std <- volcano_plot(
    df_std, "log2FC_SOD1_CTR", fdr = 0.1,
    title = paste("SOD1 vs CTR in", tissue),
    colors = colors_v)
  
  # --- 2. Standard (Adjusted) ---
  df_adj <- DE_results[[tissue]]$pgmc_adj %>%
    select(Target, UniProtID, log2FC_SOD1_CTR, pvalue_SOD1_CTR, padj_SOD1_CTR) %>%
    add_DE_flag(lfc = "log2FC_SOD1_CTR", padj = "padj_SOD1_CTR", up = "SOD1", down = "CTR", alpha = 0.1)
  
  volcano_plots[[tissue]]$adj <- volcano_plot(
    df_adj, "log2FC_SOD1_CTR", fdr = 0.1,
    title = paste("SOD1 vs CTR in", tissue, "(adjusted)"),
    colors = colors_v
  )
  
  # --- 3. Adjusted + P-Value (with special labels) ---
  df_pv <- DE_results[[tissue]]$pgmc_adj %>%
    select(Target, UniProtID, log2FC_SOD1_CTR, pvalue_SOD1_CTR, padj_SOD1_CTR, pvalue_anova) %>%
    add_DE_flag(lfc = "log2FC_SOD1_CTR", padj = "pvalue_SOD1_CTR", up = "SOD1", down = "CTR", alpha = 0.05)
  
  # Auto-calculate targets based on criteria
  targets <- df_pv %>% 
    filter(pvalue_anova < 0.05 & pvalue_SOD1_CTR < 0.05) %>% 
    pull(Target)
  
  volcano_plots[[tissue]]$pv_adj <- volcano_plot_pvalue(
    df_pv, "log2FC_SOD1_CTR", 
    title = paste("SOD1 vs CTR in", tissue, "(adj + p-value)"),
    colors = colors_v,
    label_targets = targets,
    target_color = "#702963")
}

# PDF 1: Standard vs Adjusted (FDR 10%)
pdf("plots/volcano_plots/volcano_SOD1_CTR_all_FDR10.pdf", width = 14, height = 15)
grid.arrange(
  volcano_plots$PLASMA$std, volcano_plots$PLASMA$adj,
  volcano_plots$SERUM$std,  volcano_plots$SERUM$adj,
  volcano_plots$CSF$std,    volcano_plots$CSF$adj,
  nrow = 3, ncol = 2)
dev.off()

# PDF 2: P-value Highlight vs Adjusted
pdf("plots/volcano_plots/volcano_SOD1_CTR_all_pvalue5_FDR10.pdf", width = 14, height = 15)
grid.arrange(
  volcano_plots$PLASMA$pv_adj, volcano_plots$PLASMA$adj,
  volcano_plots$SERUM$pv_adj,  volcano_plots$SERUM$adj,
  volcano_plots$CSF$pv_adj,    volcano_plots$CSF$adj,
  nrow = 3, ncol = 2)
dev.off()

volcano_plots_SOD1_CTR = volcano_plots

## 9: TARDBP vs CTR 
volcano_plots <- list()
colors_v <- c('CTR' = "#C25F3D", 'ns' = 'lightgrey', 'TARDBP'= "#3da0c2")

for (tissue in tissues) {
  message(paste("Generating volcano plots for:", tissue))
  
  # --- 1. Standard (Not Adjusted) ---
  df_std <- DE_results[[tissue]]$pgmc %>%
    select(Target, UniProtID, log2FC_TARDBP_CTR, pvalue_TARDBP_CTR, padj_TARDBP_CTR) %>%
    add_DE_flag(lfc = "log2FC_TARDBP_CTR", padj = "padj_TARDBP_CTR", up = "TARDBP", down = "CTR", alpha = 0.1)
  
  volcano_plots[[tissue]]$std <- volcano_plot(
    df_std, "log2FC_TARDBP_CTR", fdr = 0.1,
    title = paste("TARDBP vs CTR in", tissue),
    colors = colors_v)
  
  # --- 2. Standard (Adjusted) ---
  df_adj <- DE_results[[tissue]]$pgmc_adj %>%
    select(Target, UniProtID, log2FC_TARDBP_CTR, pvalue_TARDBP_CTR, padj_TARDBP_CTR) %>%
    add_DE_flag(lfc = "log2FC_TARDBP_CTR", padj = "padj_TARDBP_CTR", up = "TARDBP", down = "CTR", alpha = 0.1)
  
  volcano_plots[[tissue]]$adj <- volcano_plot(
    df_adj, "log2FC_TARDBP_CTR", fdr = 0.1,
    title = paste("TARDBP vs CTR in", tissue, "(adjusted)"),
    colors = colors_v
  )
  
  # --- 3. Adjusted + P-Value (with special labels) ---
  df_pv <- DE_results[[tissue]]$pgmc_adj %>%
    select(Target, UniProtID, log2FC_TARDBP_CTR, pvalue_TARDBP_CTR, padj_TARDBP_CTR, pvalue_anova) %>%
    add_DE_flag(lfc = "log2FC_TARDBP_CTR", padj = "pvalue_TARDBP_CTR", up = "TARDBP", down = "CTR", alpha = 0.05)
  
  # Auto-calculate targets based on criteria
  targets <- df_pv %>% 
    filter(pvalue_anova < 0.05 & pvalue_TARDBP_CTR < 0.05) %>% 
    pull(Target)
  
  volcano_plots[[tissue]]$pv_adj <- volcano_plot_pvalue(
    df_pv, "log2FC_TARDBP_CTR", 
    title = paste("TARDBP vs CTR in", tissue, "(adj + p-value)"),
    colors = colors_v,
    label_targets = targets,
    target_color = "#702963")
}

# PDF 1: Standard vs Adjusted (FDR 10%)
pdf("plots/volcano_plots/volcano_TARDBP_CTR_all_FDR10.pdf", width = 14, height = 15)
grid.arrange(
  volcano_plots$PLASMA$std, volcano_plots$PLASMA$adj,
  volcano_plots$SERUM$std,  volcano_plots$SERUM$adj,
  volcano_plots$CSF$std,    volcano_plots$CSF$adj,
  nrow = 3, ncol = 2)
dev.off()

# PDF 2: P-value Highlight vs Adjusted
pdf("plots/volcano_plots/volcano_TARDBP_CTR_all_pvalue5_FDR10.pdf", width = 14, height = 15)
grid.arrange(
  volcano_plots$PLASMA$pv_adj, volcano_plots$PLASMA$adj,
  volcano_plots$SERUM$pv_adj,  volcano_plots$SERUM$adj,
  volcano_plots$CSF$pv_adj,    volcano_plots$CSF$adj,
  nrow = 3, ncol = 2)
dev.off()

volcano_plots_TARDBP_CTR = volcano_plots

## 10: C9orf72 vs SOD1 
volcano_plots <- list()
colors_v <- c('SOD1' = "#C25F3D", 'ns' = 'lightgrey', 'C9orf72'= "#3da0c2")

for (tissue in tissues) {
  message(paste("Generating volcano plots for:", tissue))
  
  # --- 1. Standard (Not Adjusted) ---
  df_std <- DE_results[[tissue]]$pgmc %>%
    select(Target, UniProtID, log2FC_C9orf72_SOD1, pvalue_C9orf72_SOD1, padj_C9orf72_SOD1) %>%
    add_DE_flag(lfc = "log2FC_C9orf72_SOD1", padj = "padj_C9orf72_SOD1", up = "C9orf72", down = "SOD1", alpha = 0.1)
  
  volcano_plots[[tissue]]$std <- volcano_plot(
    df_std, "log2FC_C9orf72_SOD1", fdr = 0.1,
    title = paste("C9orf72 vs SOD1 in", tissue),
    colors = colors_v)
  
  # --- 2. Standard (Adjusted) ---
  df_adj <- DE_results[[tissue]]$pgmc_adj %>%
    select(Target, UniProtID, log2FC_C9orf72_SOD1, pvalue_C9orf72_SOD1, padj_C9orf72_SOD1) %>%
    add_DE_flag(lfc = "log2FC_C9orf72_SOD1", padj = "padj_C9orf72_SOD1", up = "C9orf72", down = "SOD1", alpha = 0.1)
  
  volcano_plots[[tissue]]$adj <- volcano_plot(
    df_adj, "log2FC_C9orf72_SOD1", fdr = 0.1,
    title = paste("C9orf72 vs SOD1 in", tissue, "(adjusted)"),
    colors = colors_v
  )
  
  # --- 3. Adjusted + P-Value (with special labels) ---
  df_pv <- DE_results[[tissue]]$pgmc_adj %>%
    select(Target, UniProtID, log2FC_C9orf72_SOD1, pvalue_C9orf72_SOD1, padj_C9orf72_SOD1, pvalue_anova) %>%
    add_DE_flag(lfc = "log2FC_C9orf72_SOD1", padj = "pvalue_C9orf72_SOD1", up = "C9orf72", down = "SOD1", alpha = 0.05)
  
  # Auto-calculate targets based on criteria
  targets <- df_pv %>% 
    filter(pvalue_anova < 0.05 & pvalue_C9orf72_SOD1 < 0.05) %>% 
    pull(Target)
  
  volcano_plots[[tissue]]$pv_adj <- volcano_plot_pvalue(
    df_pv, "log2FC_C9orf72_SOD1", 
    title = paste("C9orf72 vs SOD1 in", tissue, "(adj + p-value)"),
    colors = colors_v,
    label_targets = targets,
    target_color = "#702963")
}

# PDF 1: Standard vs Adjusted (FDR 10%)
pdf("plots/volcano_plots/volcano_C9orf72_SOD1_all_FDR10.pdf", width = 14, height = 15)
grid.arrange(
  volcano_plots$PLASMA$std, volcano_plots$PLASMA$adj,
  volcano_plots$SERUM$std,  volcano_plots$SERUM$adj,
  volcano_plots$CSF$std,    volcano_plots$CSF$adj,
  nrow = 3, ncol = 2)
dev.off()

# PDF 2: P-value Highlight vs Adjusted
pdf("plots/volcano_plots/volcano_C9orf72_SOD1_all_pvalue5_FDR10.pdf", width = 14, height = 15)
grid.arrange(
  volcano_plots$PLASMA$pv_adj, volcano_plots$PLASMA$adj,
  volcano_plots$SERUM$pv_adj,  volcano_plots$SERUM$adj,
  volcano_plots$CSF$pv_adj,    volcano_plots$CSF$adj,
  nrow = 3, ncol = 2)
dev.off()

volcano_plots_C9orf72_SOD1 = volcano_plots

## Save all in one pdf per fluid
master_volcano_list <- list(
  "ALS_CTR"     = volcano_plots_ALS_CTR,
  "ALS_PGMC"    = volcano_plots_ALS_PGMC,
  "ALS_mimic"   = volcano_plots_ALS_mimic,
  "mimic_CTR"   = volcano_plots_mimic_CTR,
  "PGMC_CTR"    = volcano_plots_PGMC_CTR,
  "PGMC_mimic"  = volcano_plots_PGMC_mimic,
  "C9orf72_CTR" = volcano_plots_C9orf72_CTR,
  "SOD1_CTR"    = volcano_plots_SOD1_CTR,
  "TARDBP_CTR"  = volcano_plots_TARDBP_CTR,
  "C9orf72_SOD1" = volcano_plots_C9orf72_SOD1)

tissues <- c("PLASMA", "CSF", "SERUM")

plot_versions <- list(
  "unadjusted"  = "std",
  "adjusted"    = "adj",
  "adj_PVALUE5" = "pv_adj"
)

for (tissue in tissues) {
  for (version_label in names(plot_versions)) {
    
    file_path <- paste0("plots/volcano_plots/volcano_", tissue, "_", version_label, ".pdf")
    version_key <- plot_versions[[version_label]]
    
    message("Creating: ", file_path)
    
    pdf(file_path, width = 10, height = 8)
    
    # Loop through each comparison in the master list
    for (comp_name in names(master_volcano_list)) {
      p <- master_volcano_list[[comp_name]][[tissue]][[version_key]]
      if (!is.null(p)) {
        print(p)
      } else {
        warning(paste("Plot missing for:", comp_name, tissue, version_label))
      }
    }
    dev.off()
  }
}

## signed p-values 
compute_signed <- function(df, logfc_col, comp_name, up, down, alpha_threshold = 0.1, use_fdr = TRUE){
  
  # Explicitly choose padj or pvalue
  sig_col <- if(use_fdr) paste0("padj_", comp_name) else paste0("pvalue_", comp_name)
  
  #if(!sig_col %in% colnames(df)) sig_col <- paste0("pvalue_", comp_name)
  
  df %>%
    add_DE_flag(lfc = logfc_col, padj = sig_col, up = up, down = down, alpha = alpha_threshold) %>%
    mutate(signed = sign(.data[[logfc_col]]) * log10_padj) %>%
    #mutate(signed = ifelse(is.na(signed), 0, signed)) %>%
    select(Target, UniProtID, signed, DE, any_of("pvalue_anova")) %>% 
    rename(!!paste0("signed_", comp_name) := signed,
           !!paste0("DE_", comp_name) := DE)}

classify_ALS_PGMC <- function(df){
  
  df %>%
    mutate(
      DE_new = case_when(
        DE_ALS_CTR == "ALS" & DE_ALS_PGMC == "ALS" ~ "ALS sign. both",
        DE_ALS_CTR == "CTR" & DE_ALS_PGMC == "PGMC" ~ "CTR and PGMC",
        DE_ALS_CTR == "CTR" ~ "CTR",
        DE_ALS_CTR == "ALS" ~ "ALS",
        DE_ALS_PGMC == "ALS" ~ "ALS",
        DE_ALS_PGMC == "PGMC" ~ "PGMC",
        TRUE ~ "ns"
      )
    )
  
}

classify_ALS_PGMC_CTR <- function(df){
  
  df %>%
    mutate(
      DE_new = case_when(
        DE_ALS_CTR == "ALS" & DE_PGMC_CTR == "PGMC" ~ "ALS and PGMC",
        DE_ALS_CTR == "CTR" & DE_PGMC_CTR == "CTR" ~ "CTR sign. both",
        DE_ALS_CTR == "CTR" ~ "CTR",
        DE_ALS_CTR == "ALS" ~ "ALS",
        DE_PGMC_CTR == "CTR" ~ "CTR",
        DE_PGMC_CTR == "PGMC" ~ "PGMC",
        TRUE ~ "ns"
      )
    )
  
}

classify_ALS_mimic <- function(df){
  
  df %>%
    mutate(
      DE_new = case_when(
        DE_ALS_CTR == "CTR" & DE_mimic_CTR == "CTR" ~ "CTR sign. both",
        DE_ALS_CTR == "ALS" & DE_mimic_CTR == "mimic" ~ "ALS and mimic",
        DE_ALS_CTR == "CTR" ~ "CTR",
        DE_ALS_CTR == "ALS" ~ "ALS",
        DE_mimic_CTR == "mimic" ~ "mimic",
        DE_mimic_CTR == "CTR" ~ "CTR",
        TRUE ~ "ns"
      )
    )
  
}

classify_C9_SOD1 <- function(df){
  
  df %>%
    mutate(
      DE_new = case_when(
        DE_C9orf72_CTR == "CTR" & DE_SOD1_CTR == "CTR" ~ "CTR sign. both",
        DE_C9orf72_CTR == "C9orf72" & DE_SOD1_CTR == "SOD1" ~ "C9orf72 and SOD1",
        DE_C9orf72_CTR == "CTR" ~ "CTR",
        DE_C9orf72_CTR == "C9orf72" ~ "C9orf72",
        DE_SOD1_CTR == "SOD1" ~ "SOD1",
        DE_SOD1_CTR == "CTR" ~ "CTR",
        TRUE ~ "ns"
      )
    )
  
}


signed_plot <- function(data, xvar, yvar, title, colors, label_type = "default", 
                        alpha = 0.1, highlight_color = "#702963"){
  
  line_pos <- -log10(alpha)
  anova_col <- intersect(colnames(data), c("pvalue_anova", "pvalue_anova.x", "pvalue_anova.y"))[1]
  
  if(label_type == "anova" & !is.na(anova_col)){
    label_subset <- data %>% filter(.data[[anova_col]] < 0.05 & DE_new != "ns")
  } else {
    label_subset <- data %>% filter(DE_new != "ns")
  }
  
  get_groups <- function(v) { strsplit(gsub("signed_", "", v), "_")[[1]] }
  groups_x <- get_groups(xvar)
  groups_y <- get_groups(yvar)
  
  x_title <- paste0("signed -log10 p-value (", groups_x[1], " vs ", groups_x[2], ")")
  y_title <- paste0("signed -log10 p-value (", groups_y[1], " vs ", groups_y[2], ")")
  
  limit_val <- max(abs(c(data[[xvar]], data[[yvar]])), na.rm = TRUE) * 1.1
  
  # Calibrated for a tighter "Interpretation Zone"
  x_arrow_pos <- -limit_val * 1.40 
  y_arrow_pos <- -limit_val * 1.42 
  
  p <- ggplot(data, aes(x = .data[[xvar]], y = .data[[yvar]])) +
    geom_point(aes(colour = DE_new), size = 2.25, alpha = 0.4)
  
  if(nrow(label_subset) > 0){
    if(label_type == "anova"){
      p <- p + geom_point(data = label_subset, color = highlight_color, size = 3, shape = 19)
    } else {
      p <- p + geom_point(data = label_subset, aes(colour = DE_new), size = 3, alpha = 1)
    }
    p <- p + geom_text_repel(data = label_subset, aes(label = Target), color = "black",
                             max.overlaps = 40, size = 3.5, fontface = "bold", box.padding = 0.5)
  }
  
  p <- p + 
    geom_vline(xintercept = c(line_pos, -line_pos), linetype = "dashed", colour = "grey", alpha = 0.6) + 
    geom_hline(yintercept = c(line_pos, -line_pos), linetype = "dashed", colour = "grey", alpha = 0.6) +
    scale_colour_manual(values = colors, name = "Significance") + 
    theme_classic(base_size = 12) + 
    labs(title = title, x = x_title, y = y_title) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13, margin = margin(b = 10)),
      axis.title = element_text(size = 11, face = "plain"),
      axis.title.x = element_text(margin = margin(t = 8, b = 32)), 
      axis.title.y = element_text(margin = margin(r = 8, l = 38)), 
      legend.title = element_text(face = "bold", size = 10),
      plot.margin = margin(t = 5, r = 5, b = 10, l = 10) 
    ) + 
    coord_cartesian(clip = "off", xlim = c(-limit_val, limit_val), ylim = c(-limit_val, limit_val)) 
  
  p <- p +
    annotate("segment", x = line_pos, xend = limit_val*0.8, y = x_arrow_pos, yend = x_arrow_pos, 
             arrow = arrow(length = unit(0.12, "cm")), colour = "grey30", size = 0.6) +
    annotate("text", x = (line_pos + limit_val*0.8)/2, y = x_arrow_pos, label = groups_x[1], 
             size = 3.8, fontface = "bold.italic", vjust = 1.4) + 
    
    annotate("segment", x = -line_pos, xend = -limit_val*0.8, y = x_arrow_pos, yend = x_arrow_pos, 
             arrow = arrow(length = unit(0.12, "cm")), colour = "grey30", size = 0.6) +
    annotate("text", x = (-line_pos - limit_val*0.8)/2, y = x_arrow_pos, label = groups_x[2], 
             size = 3.8, fontface = "bold.italic", vjust = 1.4) + 
    
    annotate("segment", x = y_arrow_pos, xend = y_arrow_pos, y = line_pos, yend = limit_val*0.8, 
             arrow = arrow(length = unit(0.12, "cm")), colour = "grey30", size = 0.6) +
    annotate("text", x = y_arrow_pos, y = (line_pos + limit_val*0.8)/2, label = groups_y[1], 
             angle = 90, size = 3.8, fontface = "bold.italic", vjust = -1.1) + 
    
    annotate("segment", x = y_arrow_pos, xend = y_arrow_pos, y = -line_pos, yend = -limit_val*0.8, 
             arrow = arrow(length = unit(0.12, "cm")), colour = "grey30", size = 0.6) +
    annotate("text", x = y_arrow_pos, y = (-line_pos - limit_val*0.8)/2, label = groups_y[2], 
             angle = 90, size = 3.8, fontface = "bold.italic", vjust = -1.1)
  
  return(p)
}

# Define the components 
tissues <- c("PLASMA", "CSF", "SERUM")
comparisons <- c("ALS_CTR", "ALS_PGMC", "PGMC_CTR", "mimic_CTR", "C9orf72_CTR", "SOD1_CTR")

datasets <- list()

for (t in tissues) {
  
  # 1. Unadjusted
  raw_combined <- inner_join(
    DE_results[[t]]$standard, 
    DE_results[[t]]$pgmc %>% select(-any_of(c("Fluid", "pvalue_anova"))), 
    by = c("Target", "UniProtID")
  )
  
  # 2. ADJUSTED
  adj_combined <- inner_join(
    DE_results[[t]]$standard_adj, 
    DE_results[[t]]$pgmc_adj %>% select(-any_of(c("Fluid", "pvalue_anova"))), 
    by = c("Target", "UniProtID")
  )
  
  # Assign to the datasets structure
  datasets[[tolower(t)]] <- list(
    raw = lapply(comparisons, function(comp) raw_combined),
    adjusted = lapply(comparisons, function(comp) adj_combined),
    adjusted_pvalue = lapply(comparisons, function(comp) adj_combined)
  )
  
  # Name the internal comparison lists
  names(datasets[[tolower(t)]]$raw) <- comparisons
  names(datasets[[tolower(t)]]$adjusted) <- comparisons
  names(datasets[[tolower(t)]]$adjusted_pvalue) <- comparisons
}

# colors
colors <- c(
  'ALS'='#32A55E',
  'PGMC'='#96AA9A',
  'mimic'='#96AA9A',
  'C9orf72'='#96AA9A',
  'SOD1'='#CF7041',
  'CTR'='#CF7041',
  'ALS sign. both'='#5E718B',
  'CTR sign. both' = '#E5C5BD',
  'ALS and mimic'='#E5C5BD',
  'CTR and PGMC' = '#E5C5BD',
  'ALS and PGMC' = '#CF7041',
  'C9orf72 and SOD1'='#E5C5BD',
  'ns'='#B4BFC5')

# run all signed plots
for(tissue in names(datasets)){
  
  for(version in names(datasets[[tissue]])){
    
    tissue_data <- datasets[[tissue]][[version]]
    
    # Logic for significance selection
    is_pvalue_version <- (version == "adjusted_pvalue")
    use_fdr_flag <- !is_pvalue_version
    alpha <- ifelse(is_pvalue_version, 0.05, 0.1)
    current_label_type <- ifelse(is_pvalue_version,"anova","default")
    
    ### ALS vs CTR + ALS vs PGMC
    
    df1 <- compute_signed(tissue_data$ALS_CTR,"log2FC_ALS_CTR","ALS_CTR",
                          up = "ALS",down = "CTR",
                          alpha_threshold = alpha,
                          use_fdr = use_fdr_flag)
    df2 <- compute_signed(tissue_data$ALS_PGMC,"log2FC_ALS_PGMC","ALS_PGMC",
                          up = "ALS",down = "PGMC",
                          alpha_threshold = alpha,
                          use_fdr = use_fdr_flag)
    
    merged <- inner_join(df1,df2) %>%
      classify_ALS_PGMC()
    
    p <- signed_plot(
      merged,
      "signed_ALS_CTR",
      "signed_ALS_PGMC",
      paste0("ALS vs CTR + ALS vs PGMC (",toupper(tissue),")"),
      colors,
      label_type = current_label_type,alpha = alpha)
    
    ggsave(
      paste0(
        "plots/signed_plots/signed_ALSvsCTR_ALSvsPGMC_",
        tissue,"_",
        version,
        ".pdf"
      ),
      p,
      width=7,
      height=6
    )
    
    ### ALS vs CTR + PGMC vs CTR
    
    df1 <- compute_signed(tissue_data$ALS_CTR,"log2FC_ALS_CTR","ALS_CTR",
                          up = "ALS",down = "CTR",
                          alpha_threshold = alpha,
                          use_fdr = use_fdr_flag)
    df2 <- compute_signed(tissue_data$PGMC_CTR,"log2FC_PGMC_CTR","PGMC_CTR",
                          up = "PGMC",down = "CTR",
                          alpha_threshold = alpha,
                          use_fdr = use_fdr_flag)
    
    merged <- inner_join(df1,df2) %>%
      classify_ALS_PGMC_CTR()
    
    p <- signed_plot(
      merged,
      "signed_ALS_CTR",
      "signed_PGMC_CTR",
      paste0("ALS vs CTR + PGMC vs CTR (",toupper(tissue),")"),
      colors,
      label_type = current_label_type,alpha = alpha)
    
    ggsave(
      paste0(
        "plots/signed_plots/signed_ALSvsCTR_PGMCvsCTR_",
        tissue,"_",
        version,
        ".pdf"
      ),
      p,
      width=7,
      height=6
    )
    
    
    ### ALS vs CTR + mimic vs CTR
    
    df1 <- compute_signed(tissue_data$ALS_CTR,"log2FC_ALS_CTR","ALS_CTR",
                          up = "ALS",down = "CTR",
                          alpha_threshold = alpha,
                          use_fdr = use_fdr_flag)
    df2 <- compute_signed(tissue_data$mimic_CTR,"log2FC_mimic_CTR","mimic_CTR",
                          up = "mimic",down = "CTR",
                          alpha_threshold = alpha,
                          use_fdr = use_fdr_flag)
    
    merged <- inner_join(df1,df2) %>%
      classify_ALS_mimic()
    
    p <- signed_plot(
      merged,
      "signed_ALS_CTR",
      "signed_mimic_CTR",
      paste0("ALS vs CTR + mimic vs CTR (",toupper(tissue),")"),
      colors,
      label_type = current_label_type,alpha = alpha)
    
    ggsave(
      paste0(
        "plots/signed_plots/signed_ALSvsCTR_mimicvsCTR_",
        tissue,"_",
        version,
        ".pdf"
      ),
      p,
      width=7,
      height=6
    )
    
    
    ### C9orf72 vs CTR + SOD1 vs CTR
    
    df1 <- compute_signed(tissue_data$C9orf72_CTR, "log2FC_C9orf72_CTR", "C9orf72_CTR", 
                          up = "C9orf72", down = "CTR", alpha_threshold = alpha,
                          use_fdr = use_fdr_flag)
    df2 <- compute_signed(tissue_data$SOD1_CTR, "log2FC_SOD1_CTR", "SOD1_CTR", 
                          up = "SOD1", down = "CTR", alpha_threshold = alpha,
                          use_fdr = use_fdr_flag)
    
    merged <- inner_join(df1,df2) %>%
      classify_C9_SOD1()
    
    p <- signed_plot(
      merged,
      "signed_C9orf72_CTR",
      "signed_SOD1_CTR",
      paste0("C9orf72 vs CTR + SOD1 vs CTR (",toupper(tissue),")"),
      colors,
      label_type = current_label_type,alpha = alpha)
    
    ggsave(
      paste0(
        "plots/signed_plots/signed_C9orf72vsCTR_SOD1vsCTR_",
        tissue,"_",
        version,
        ".pdf"
      ),
      p,
      width=7,
      height=6
    )
  }
}