source("00_initialization.R")
source("01_data_mining.R")

###############################################
### Helper Functions
###############################################

## 1. Clean datasets (all samples OR PGMC subset)
prepare_dataset <- function(protein_data, sample_map) {
  protein_data %>%
    filter(SampleType == "Sample") %>%
    select(SampleName, SampleMatrixType, Target, UniProtID, ProteinName, NPQ) %>%
    left_join(sample_map %>% rename(SampleName = `Sample ID`), by = "SampleName") %>%
    filter(!is.na(type)) %>%
    mutate(type = factor(type))
}

## 2. Run ANOVA + pairwise tests for one fluid (CSF/SERUM/PLASMA)
run_stats_for_fluid <- function(df, fluid) {
  
  df_f <- df %>%
    filter(SampleMatrixType == fluid) %>%
    filter(!is.na(NPQ))
  
  # ANOVA
  anova_res <- df_f %>%
    group_by(Target) %>%
    rstatix::anova_test(NPQ ~ type, effect.size = "partial_eta_squared") %>%
    as_tibble() %>%
    mutate(Fluid = fluid, Model = "ANOVA")
  
  
  if(nrow(anova_res) == 0) {
    warning("No ANOVA results for fluid: ", fluid)
    anova_res <- NULL
  }
  
  # Welch t-tests
  pairwise_res <- df_f %>%
    filter(type %in% c("ALS","CTR","PGMC","mimic","C9orf72","SOD1","TARDBP")) %>%
    group_by(Target) %>%
    rstatix::t_test(NPQ ~ type, var.equal = FALSE) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance() %>%
    mutate(Fluid = fluid, Model = "Pairwise")
  
  if(nrow(pairwise_res) == 0) {
    warning("No pairwise results for fluid: ", fluid)
    pairwise_res <- NULL
  }
  
  list(
    anova = anova_res,
    pairwise = pairwise_res,
    data = df_f
  )
}

## 3. Extract p-values for all targets
extract_pvalues <- function(stats_list) {
  
  pw <- stats_list$pairwise %>%
    mutate(
      comparison = paste(group1, group2, sep = "_")
    ) 
  
  pval_wide = pw %>%
    select(Fluid, Target, comparison, p) %>%
    pivot_wider(
      names_from = comparison,
      values_from = p,
      names_glue = "pvalue_{comparison}"
    )
  
  padj_wide = pw %>%
    select(Fluid, Target, comparison, p.adj) %>%
    pivot_wider(
      names_from = comparison,
      values_from = p.adj,
      names_glue = "padj_{comparison}"
    )
  
  signif_wide = pw %>%
    select(Fluid, Target, comparison, p.adj.signif) %>%
    pivot_wider(
      names_from = comparison,
      values_from = p.adj.signif,
      names_glue = "padj_signif_{comparison}"
    )
  
  combined <- pval_wide %>%
    left_join(padj_wide, by = c("Fluid","Target")) %>%
    left_join(signif_wide, by = c("Fluid","Target")) 
  
  comparisons <- unique(pw$comparison)
  
  ordered_cols <- c(
    "Fluid", "Target",
    unlist(lapply(comparisons, function(cmp) {
      c(
        paste0("pvalue_", cmp),
        paste0("padj_", cmp),
        paste0("padj_signif_", cmp)
      )
    }))
  )
  
  final_table <- combined %>% select(any_of(ordered_cols))
  
  pw_final = stats_list$anova %>% select(Target,pvalue_anova = p) %>%
    left_join(final_table)  %>%
    arrange(pvalue_anova)
  
  pw_final = pw_final[,c(1,3,2,4:ncol(pw_final))]
 
  return(pw_final)
}

## 4. Get LOD for a given protein and fluid 
get_lod <- function(protein, fluid, td) {
  lod_val <- td %>%
    filter(Target == protein, SampleMatrixType == fluid) %>%
    pull(TargetLOD_NPQ)
  if(length(lod_val) == 0) return(NA) else return(lod_val)
}

## 5. Get pairwise p-values filtered by cutoff for plotting
get_pairwise_sig <- function(stats_list, fluid, proteins, df_top, p_cutoff = 0.1) {
  
  pairwise_df <- stats_list$pairwise %>%
    filter(Target %in% proteins, p <= p_cutoff)
  
  if(nrow(pairwise_df) == 0) return(pairwise_df)
  
  # Max NPQ per protein 
  max_vals <- df_top %>%
    group_by(Target) %>%
    summarise(max_y = max(NPQ, na.rm = TRUE), .groups = "drop")
  
  # Assign incremental y positions
  pairwise_df <- pairwise_df %>%
    left_join(max_vals, by = "Target") %>%
    group_by(Target) %>%
    mutate(
      comparison_index = row_number(),
      y.position = max_y * (0.9 + 0.1 * (comparison_index - 1))
    ) %>%
    ungroup()
  
  return(pairwise_df)
}

## 6. Top N significantly changing proteins 
plot_top_proteins_violin <- function(df, stats_list, td, fluid, top_n = 15) {
  
  # Check ANOVA results
  anova_res <- stats_list$anova
  
  if(is.null(anova_res) || nrow(anova_res) == 0) {
    warning("No ANOVA results for ", fluid, "; skipping plot.")
    return(NULL)
  }
  
  # Select top proteins
  top_proteins_df <- anova_res %>%
    arrange(p) %>%
    slice(1:min(top_n, nrow(.)))
  
  top_proteins <- top_proteins_df$Target
  df_top <- df %>% filter(SampleMatrixType == fluid, 
                          Target %in% top_proteins,
                          type %in% c("ALS","CTR","PGMC","mimic","C9orf72","SOD1","TARDBP"))
  
  if(nrow(df_top) == 0) {
    warning("No data for top proteins in ", fluid, "; skipping plot.")
    return(NULL)
  }
  
  # Pairwise p-values
  pairwise_res <- get_pairwise_sig(stats_list, fluid, top_proteins, df_top)
  
  # LOD per protein
  lod_df <- df_top %>%
    group_by(Target) %>%
    summarise(
      LOD_val = get_lod(Target[1], fluid, td),
      max_val = max(NPQ, na.rm = TRUE)
    ) %>%
    mutate(y_position_text = LOD_val * 1.02) %>%
    ungroup()
  
  pgmc_mutation_proteins <- c("C9orf72", "SOD1", "TARDBP")
  df_top <- df_top %>%
    mutate(type = case_when(
      type %in% pgmc_mutation_proteins ~ factor(type, levels = c("CTR", "C9orf72", "SOD1", "TARDBP")),
      TRUE ~ factor(type, levels = c("CTR", "PGMC", "ALS", "mimic"))
    ))
  
  # Build plot
  p <- ggplot(df_top, aes(x = type, y = NPQ, fill = type)) +
    geom_violin(trim = FALSE, alpha = 0.4) +
    geom_boxplot(width = 0.2, outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.5) +
    facet_wrap(~Target, scales = "free_y") +
    scale_fill_manual(values  = c('CTR' = '#6F8EB2',  
                                  'ALS' = '#B2936F',
                                  'PGMC' = '#ad5291',
                                  'mimic' = '#62cda9',
                                  'other' = '#ad5291',
                                  'C9orf72' = '#55aa82',
                                  'SOD1' = '#4661b9',
                                  'TARDBP' = '#B99E46')) +
    labs(
      x = "Group",
      y = "NPQ",
      title = paste("Top", top_n, "Protein Expression -", fluid)
      ) +
    theme(
      panel.background  = element_rect(fill = "white", color = NA),
      plot.background   = element_rect(fill = "white", color = NA),
      panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.title.x = element_blank(),
      text = element_text(size = 15),               # Base font size for everything
      axis.title = element_text(size = 18),         # Axis titles
      axis.text = element_text(size = 14),          # Axis tick labels
      strip.text = element_text(size = 16, face = "bold"),  # Facet labels
      plot.title = element_text(size = 18, hjust = 0.5),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14),
      legend.position = "none"
    ) 
  
  
  # Add LOD lines and labels
  p <- p +
    geom_hline(
      data = lod_df %>% tidyr::unnest(LOD_val),
      aes(yintercept = LOD_val),
      linetype = "dashed",
      color = "gray55"
    ) 
  # +
  #   geom_text(
  #     data = lod_df %>% mutate(label = paste("LOD =", paste(signif(LOD_val, 3), collapse = ", "))),
  #     aes(x = 4.2, y = y_position_text, label = label),
  #     color = "gray55",
  #     size = 4,
  #     hjust = 0,
  #     inherit.aes = FALSE
  #   )
    
  # Add pairwise p-values only if available
  if (!is.null(pairwise_res) && nrow(pairwise_res) > 0) {
    p <- p + stat_pvalue_manual(
      pairwise_res,
      label = "p",
      size = 5,
      bracket.nudge.y = 0.04 * max(df_top$NPQ, na.rm = TRUE),
      # xmin = "group1",
      # xmax = "group2",
      y.position = "y.position"
      # tip.length = 0.03
    )
  }
  
  return(p)
}

## 7. Single proteins
plot_single_protein_violin <- function(df, stats_list, td, fluid, protein) {
  
  df_prot <- df %>% filter(SampleMatrixType == fluid, 
                           Target == protein,
                           type %in% c("ALS","CTR","PGMC","mimic","C9orf72","SOD1","TARDBP"))
  anova_p <- stats_list$anova %>% filter(Target == protein) %>% pull(p)
  pairwise_res <- get_pairwise_sig(stats_list, fluid, protein, df_top = df_prot)
  
  # Determine type order
  pgmc_mutation_proteins <- c("C9orf72", "SOD1", "TARDBP")
  
  df_prot = df_prot %>%
    mutate(type = case_when(
      type %in% pgmc_mutation_proteins ~ factor(type, levels = c("CTR", "C9orf72", "SOD1", "TARDBP")),
      TRUE ~ factor(type, levels = c("CTR", "PGMC", "ALS", "mimic"))
    ))
  
  # LOD for this protein
  LOD_val <- get_lod(protein, fluid, td) %>% unique()
  max_val <- max(df_prot$NPQ, na.rm = TRUE)
  y_text <- max_val * 1.02
  
  p <- ggplot(df_prot, aes(x = type, y = NPQ, fill = type)) +
    geom_violin(trim = FALSE, alpha = 0.4) +
    geom_jitter(width = 0.15, alpha = 0.5) +
    scale_fill_manual(values  = c('CTR' = '#6F8EB2',  
                                  'ALS' = '#B2936F',
                                  'PGMC' = '#ad5291',
                                  'mimic' = '#62cda9',
                                  'other' = '#ad5291',
                                  'C9orf72' = '#55aa82',
                                  'SOD1' = '#4661b9',
                                  'TARDBP' = '#B99E46')) +
    labs(
      x = "Group",
      y = "NPQ",
      title = paste(protein),
      subtitle = paste("ANOVA p =", signif(anova_p, 3),
                       "; LOD1 =", signif(LOD_val[1],3), 
                       "; LOD2 =", signif(LOD_val[2],3))
    ) +
    theme_test() + 
    theme(
      panel.background  = element_rect(fill = "white", color = NA),
      plot.background   = element_rect(fill = "white", color = NA),
      axis.title.x = element_blank(),
      text = element_text(size = 15),               # Base font size for everything
      axis.title = element_text(size = 18),         # Axis titles
      axis.text = element_text(size = 14),          # Axis tick labels
      strip.text = element_text(size = 16, face = "bold"),  # Facet labels
      plot.title = element_text(size = 18, hjust = 0.5,face = "bold"),
      plot.subtitle = element_text(size = 16, hjust = 0.5),
      legend.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 14),
      legend.position = "none"
    ) +
    geom_hline(yintercept = LOD_val, linetype = "dashed", color = "gray55") +
    annotate("text",
             x = 3.5,
             y = LOD_val * 1.02,
             label = paste0("LOD = ", signif(LOD_val, 3)),
             color = "gray55",
             hjust = 0,
             size = 5)
  
  # p <- p +
  #   annotation_custom(
  #     grob = grid::rectGrob(
  #       gp = grid::gpar(fill = "white", col = "black", lwd = 1)
  #     ),
  #     ymin = Inf, ymax = Inf, xmin = -Inf, xmax = Inf
  #   ) +
  #   annotation_custom(
  #     grob = grid::textGrob(
  #       label = protein,
  #       gp = grid::gpar(fontsize = 16, fontface = "bold")
  #     ),
  #     ymin = Inf, ymax = Inf, xmin = -Inf, xmax = Inf
  #   ) +
  #   theme(
  #     plot.margin = margin(t = 30),  # extra room for strip
  #     strip.background = element_blank(),
  #     strip.text = element_blank()
  #   )
    
  if (nrow(pairwise_res) > 0) {
    p <- p + stat_pvalue_manual(
      pairwise_res,
      label = "p",
      bracket.nudge.y = 0.1 * max_val,
      size = 6,
      # xmin = "group1",
      # xmax = "group2",
      y.position = "y.position"
      # tip.length = 0.03
    )
  }
  
  p = p + expand_limits(y = max_val * 1.25)
  
  return(p)
}

## 8. Full pipeline for one dataset (all samples OR PGMC subset)
run_full_pipeline <- function(protein_data, sample_map, td, prefix = "ALL") {
  
  message("Preparing dataset...")
  df <- prepare_dataset(protein_data, sample_map)
  
  fluids <- c("CSF", "SERUM", "PLASMA")
  results <- list()
  
  for (fluid in fluids) {
    message("Running stats for ", fluid, " ...")
    stats <- run_stats_for_fluid(df, fluid)
    pvals <- extract_pvalues(stats)
    
    results[[fluid]] <- list(
      stats = stats,
      pvals = pvals
    )
    
    # Save p-values
    write_xlsx(pvals,
               paste0("results/", prefix, "_pvalues_", fluid, ".xlsx"))
    
    # Top 15 plot (violin)
    p_top15 <- plot_top_proteins_violin(df, stats, td, fluid)
    ggsave(paste0("plots/boxplots_", fluid, "/", prefix, "_Top15_", fluid, ".pdf"),
           p_top15, width = 15, height = 18)
    
    # All proteins (multi-page PDF)
    targets <- unique(df$Target)
    pdf(paste0("plots/boxplots_", fluid,"/", prefix, "_ALLproteins_", fluid, ".pdf"),
        width = 6, height = 5.7)
    for (t in targets) {
      p <- plot_single_protein_violin(df, stats, td, fluid, t)
      print(p)
    }
    dev.off()
  }
  
  return(results)
}

## 9. PCA per fluid labelled based on status
run_pca <- function(df, matrix_type) {
  df_matrix <- df %>%
    filter(SampleMatrixType == matrix_type) %>%
    select(SampleName, Target, NPQ, type, subtype)
  
  # Reshape wide
  df_wide <- df_matrix %>%
    pivot_wider(
      names_from = Target,
      values_from = NPQ,
      values_fill = 0
    )
  
  meta <- df_wide %>% select(SampleName, type, subtype)
  X <- df_wide %>% select(-SampleName, -type, -subtype)
  
  X <- X %>% select(where(~ {
    s <- sd(.x, na.rm = TRUE)
    !is.na(s) && s > 0
  }))
  
  # PCA
  pca <- prcomp(X, scale. = TRUE)
  
  scores <- as_tibble(pca$x[, 1:2]) %>%
    bind_cols(meta)
  
  list(scores = scores, pca = pca)
}
# Step 3: Define color palette
my_colors <- c(
  'CTR'    = '#6F8EB2',
  'ALS'    = '#B2936F',
  'PGMC'   = '#ad5291',
  'mimic'  = '#62cda9',
  'other'  = '#ad5291',
  'C9orf72'= '#55aa82',
  'SOD1'   = '#4661b9',
  'TARDBP' = '#B99E46'
)

# Step 5: Plot function
plot_pca <- function(pca_res, title) {
  pca <- pca_res$pca
  scores <- pca_res$scores %>%
    mutate(
      subtype = factor(
        subtype,
        levels = c("No subtype", "C9orf72", "SOD1", "TARDBP","FUS","other")
      )
    )
  
  ggplot(scores, aes(x = PC1, y = PC2, color = type, shape = subtype)) +
    geom_point(size = 5, alpha = 0.8) +
    scale_color_manual(values = my_colors) +
    theme_minimal(base_size = 16) +
    labs(
      title = paste("PCA -", title),
      x = paste0("PC1 (", round(100 * summary(pca)$importance[2,1], 1), "%)"),
      y = paste0("PC2 (", round(100 * summary(pca)$importance[2,2], 1), "%)")
    ) +
    theme(
      text = element_text(size = 16),               # Base font size for everything
      axis.title = element_text(size = 18),         # Axis titles
      axis.text = element_text(size = 16),          # Axis tick labels
      plot.title = element_text(size = 18, hjust = 0.5,face = "bold"),
      legend.title = element_text(size = 17),
      legend.text = element_text(size = 16)
    )
}

# PCA all together
run_pca_all <- function(df) {
  
  # Ensure NPQ is numeric
  df <- df %>%
    mutate(NPQ = suppressWarnings(as.numeric(NPQ)))
  
  # Reshape wide
  df_wide <- df %>%
    tidyr::pivot_wider(
      names_from = Target,
      values_from = NPQ,
      values_fill = 0
    )
  
  # Metadata
  meta <- df_wide %>%
    dplyr::select(SampleName, SampleMatrixType)
  
  # Numeric matrix for PCA
  X <- df_wide %>%
    dplyr::select(-SampleName, -SampleMatrixType)
  
  # Force numeric (some columns may still be character/factor)
  X <- X %>% mutate(across(everything(), ~ suppressWarnings(as.numeric(.x))))
  
  # X <- X %>% select(where(~ {
  #   s <- sd(.x, na.rm = TRUE)
  #   !is.na(s) && s > 0
  # }))
  
  # Drop columns with all NA or zero variance
  X <- X %>%
    dplyr::select(where(~ {
      all_na <- all(is.na(.x))
      s <- sd(.x, na.rm = TRUE)
      !all_na && !is.na(s) && s > 0
    }))
  
  # Replace remaining NA with 0 
  X[is.na(X)] <- 0
  X <- X[,!names(X) %in% c("APOE4","CRP")]
  
  print(X)
  
  # PCA
  pca <- prcomp(X, scale. = TRUE)
  
  scores <- as_tibble(pca$x[, 1:2]) %>%
    bind_cols(meta)
  
  list(scores = scores, pca = pca)
}

my_colors <- c(
  'CSF'    = '#1B9E77',
  'PLASMA' = '#D95F02',
  'SERUM'  = '#7570B3'
)

###############################################################################
# Run pipeline
###############################################################################

## 1. All samples
results_ALL <- run_full_pipeline(
  protein_data    = protein_data,
  sample_map      = samples_ID_type,
  td              = td,
  prefix          = "ALLsamples"
)

## 2. PGMC mutation analysis
results_PGMC <- run_full_pipeline(
  protein_data = protein_data,
  sample_map   = samples_PGMC_CTR_ID_type,
  td           = td,
  prefix       = "PGMCvsCTR"
)

## 3. PCA by fluid and subtypes
protein_data_PCA <- protein_data_IDs %>%
  left_join(samples_PGMC_CTR_ID_type %>% rename(subtype = type,
                                                SampleName = `Sample ID`)) %>%
  mutate(
    subtype = ifelse(subtype == "CTR" | is.na(subtype), "No subtype",
                     ifelse(subtype %in% c("C9orf72","SOD1","TARDBP","FUS"), subtype, "other"))
  )

fluids <- c("SERUM", "PLASMA", "CSF")

subtypes_counts <- bind_rows(
  lapply(fluids, function(f) count_subtypes_per_fluid(protein_data_PCA, f))
)

sample_counts = rbind(sample_counts,subtypes_counts %>% rename(type = subtype)) %>%
  arrange(biofluid) %>%
  filter(type != "No subtype")

writexl::write_xlsx(sample_counts, "results/samples_biofluid_overview.xlsx")

protein_data_clean <- protein_data_PCA %>% filter(!is.na(type))

matrices <- unique(protein_data_clean$SampleMatrixType)
pca_results_subtype <- lapply(matrices, function(m) run_pca(protein_data_clean, m))
names(pca_results_subtype) <- matrices

plots_subtype <- lapply(names(pca_results_subtype), 
                        function(m) plot_pca(pca_results_subtype[[m]], m))

pdf("plots/PCA_SERUM_subtype.pdf", width = 8, height = 6.5) 
plots_subtype[[which(names(pca_results_subtype)=="SERUM")]] 
dev.off()
pdf("plots/PCA_PLASMA_subtype.pdf", width = 8, height = 6.5)
plots_subtype[[which(names(pca_results_subtype)=="PLASMA")]]
dev.off()
pdf("plots/PCA_CSF_subtype.pdf", width = 8, height = 6.5)
plots_subtype[[which(names(pca_results_subtype)=="CSF")]]
dev.off()

## 4. PCA all biofluids together
protein_data_PCA_all <- protein_data_IDs %>%
  filter(SampleMatrixType %in% c("CSF","SERUM","PLASMA")) %>%
  select(SampleName, Target, NPQ, SampleMatrixType)

pca_results_all <- run_pca_all(protein_data_PCA_all)
pdf("plots/PCA_all_fluids.pdf", width = 8, height = 6.5)
plot_pca_all(pca_results_all)
dev.off()


