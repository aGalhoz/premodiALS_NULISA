## compare CNS and immune panels

protein_data_IDs_CNS = protein_data_IDs
protein_data_IDs_immune = read_excel("~/Documents/HMGU/premodiALS/NULISA - immune panel/data input/P004_BSHRI_NULISAseq_InflammationPanel_NPQCounts_2025_10_10.xlsx") %>%
  filter(SampleType == "Sample") %>%
  select(SampleName, SampleMatrixType, Target, UniProtID, ProteinName, NPQ) %>%
  left_join(samples_ID_type %>% rename(SampleName = `Sample ID`), by = "SampleName")

protein_data_IDs_CNS = protein_data_IDs_CNS %>%
  mutate(panel = "CNS") %>%
  filter(!is.na(type))
proteins_data_IDs_immune = protein_data_IDs_immune %>%
  mutate(panel = "IMMUNE")  %>%
  filter(!is.na(type))

## Combine into one long dataframe
all_data <- bind_rows(protein_data_IDs_CNS, proteins_data_IDs_immune) %>%
  mutate(
    panel            = factor(panel, levels = c("CNS", "IMMUNE")),
    SampleMatrixType = factor(SampleMatrixType, levels = c("SERUM", "PLASMA", "CSF"))
  ) %>%
  filter(!is.na(NPQ)) %>%
  filter(!Target %in% c("APOE4","APOE","CRP","KNG1"))

## =========================================================================
## PART 1: HEATMAPS -- one per fluid, rows split by panel (CNS vs Immune)
## =========================================================================

## Build a protein x sample matrix of NPQ values for a given fluid,
## and a matching vector of panel labels (in the same row order)
build_matrix_for_fluid <- function(data, fluid) {
  df_fluid <- data %>% filter(SampleMatrixType == fluid)
  
  ## collapse duplicate protein/sample combos (if any) by mean NPQ
  df_wide <- df_fluid %>%
    group_by(panel, Target, SampleName) %>%
    summarise(NPQ = mean(NPQ, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = SampleName, values_from = NPQ)
  
  panel_vec <- df_wide$panel
  mat <- df_wide %>% select(-panel, -Target) %>% as.matrix()
  rownames(mat) <- df_wide$Target
  
  list(mat = mat, panel = panel_vec)
}

## Draw + return a ComplexHeatmap object for one fluid
make_fluid_heatmap <- function(data, fluid) {
  built <- build_matrix_for_fluid(data, fluid)
  mat   <- built$mat
  panel_vec <- built$panel
  
  q <- quantile(mat, probs = c(0.01, 0.5, 0.99), na.rm = TRUE)
  col_fun <- colorRamp2(
    q,
    c("#440154", "#21918C", "#FDE725")   # viridis: dark purple -> teal -> yellow
  )
  
  Heatmap(
    mat,
    name             = "NPQ",
    col              = col_fun,
    row_split        = panel_vec,          
    row_title        = c("CNS", "IMMUNE"),
    row_title_gp     = gpar(fontsize = 12, fontface = "bold"),
    cluster_rows     = TRUE,
    cluster_row_slices = FALSE,
    cluster_columns  = TRUE,
    show_row_names   = TRUE,
    show_column_names= FALSE,              # set TRUE if you want sample IDs
    row_names_gp     = gpar(fontsize = 6),
    column_title     = paste0(fluid, " - CNS & Immune panels"),
    heatmap_legend_param = list(title = "NPQ")
  )
}

## Save one heatmap per fluid into a single multi-page PDF
fluids <- levels(all_data$SampleMatrixType)

pdf("plots/CNS_IMMUNE_panels/heatmaps_by_fluid.pdf", width = 10, height = 12)
for (fl in fluids) {
  if (nrow(filter(all_data, SampleMatrixType == fl)) == 0) next
  ht <- make_fluid_heatmap(all_data, fl)
  draw(ht)
}
dev.off()

## =========================================================================
## PART 2: PROTEINS COMMON TO BOTH PANELS
## =========================================================================

cns_ids    <- protein_data_IDs_CNS    %>% 
  filter(!Target %in% c("APOE4","APOE","CRP","KNG1")) %>% 
  distinct(UniProtID, Target)
immune_ids <- proteins_data_IDs_immune  %>% 
  filter(!Target %in% c("APOE4","APOE","CRP","KNG1")) %>% 
  distinct(UniProtID, Target)

## Proteins measured in both panels (matched by UniProtID)
common_proteins <- inner_join(cns_ids, immune_ids, by = "UniProtID", suffix = c("_cns", "_immune"))

cat("Number of proteins common to both panels:", nrow(common_proteins), "\n")
print(common_proteins)

## =========================================================================
## PART 3: BOXPLOTS -- CNS vs Immune NPQ, per common protein, per fluid
## =========================================================================

## Subset to just the common proteins, for both panels
common_data <- all_data %>%
  filter(UniProtID %in% common_proteins$UniProtID)

## One boxplot per protein: x = panel (CNS/Immune), facet by fluid
plot_protein_boxplot <- function(data, uniprot_id) {
  df <- data %>% filter(UniProtID == uniprot_id)
  protein_label <- unique(df$Target)[1]
  
  ggplot(df, aes(x = panel, y = NPQ, fill = panel)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8) +
    geom_jitter(width = 0.15, size = 0.8, alpha = 0.4) +
    facet_wrap(~ SampleMatrixType, nrow = 1) +
    scale_fill_manual(values = c("CNS" = "#4477AA", "IMMUNE" = "#CC6677")) +
    labs(
      title = paste0(protein_label, " (", uniprot_id, ")"),
      x = NULL, y = "NPQ"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold"))
}

## 15 proteins per page
plot_list <- lapply(unique(common_proteins$UniProtID), function(uid) {
  plot_protein_boxplot(common_data, uid) +
    theme(
      plot.title  = element_text(size = 9, face = "bold"),
      axis.text   = element_text(size = 6),
      axis.title  = element_text(size = 7),
      strip.text  = element_text(size = 7)
    )
})

n_per_page <- 15   
n_row <- 5
n_col <- 3

ggsave_pages <- marrangeGrob(
  grobs = plot_list,
  nrow  = n_row,
  ncol  = n_col,
  top   = NULL
)

ggsave(filename = "plots/CNS_IMMUNE_panels/common_protein_boxplots.pdf",
  plot = ggsave_pages, width    = 14, height   = 16)

## ==========================================================
## PART 4: ML -- CNS and Immune panels together, per fluid
## ==========================================================
## Find highly correlated protein pairs across the CNS and immune panels

## combined data of CNS and immune panels
combined_long <- common_data %>%
  filter(panel %in% c("CNS", "IMMUNE")) %>%         
  mutate(protein_panel = paste(Target, panel, sep = "__")) %>%
  group_by(SampleName, protein_panel) %>%
  summarise(NPQ = mean(NPQ, na.rm = TRUE), .groups = "drop")

mat <- combined_long %>%
  pivot_wider(names_from = protein_panel, values_from = NPQ) %>%
  column_to_rownames("SampleName")

## correlation matrix
corr_mat <- cor(mat, use = "pairwise.complete.obs", method = "pearson")  # or "spearman"
corr_long <- corr_mat
corr_long[upper.tri(corr_long, diag = TRUE)] <- NA

corr_pairs <- as_tibble(corr_long, rownames = "Protein1") %>%
  pivot_longer(-Protein1, names_to = "Protein2", values_to = "r") %>%
  filter(!is.na(r)) %>%
  separate(Protein1, into = c("Target1", "Panel1"), sep = "__") %>%
  separate(Protein2, into = c("Target2", "Panel2"), sep = "__") %>%
  filter(Target1 == Target2) %>%
  arrange(desc(abs(r)))

write_csv(corr_pairs, "results/high_correlation_protein_pairs.csv")

high_corr_pairs = corr_pairs %>% 
  filter(r>0.8)
high_corr_proteins = high_corr_pairs %>% pull(Target1)

## ------------------------------------------------------------------
## Per-protein correlation: CNS panel vs Immune panel NPQ

threshold <- 0.8

fluid_counts <- common_data %>% distinct(SampleName, SampleMatrixType) %>% count(SampleMatrixType)
fluids <- fluid_counts %>% pull(SampleMatrixType) %>% as.character()

## same-protein CNS vs Immune correlation for each fluid
correlate_panels_by_fluid <- function(data, fluid_name) {
  
  fluid_data <- data %>% filter(SampleMatrixType == fluid_name, panel %in% c("CNS", "IMMUNE"))
  
  ## keep only proteins measured on both panels within this fluid
  shared_proteins <- fluid_data %>%
    distinct(Target, panel) %>%
    count(Target) %>%
    filter(n == 2) %>%
    pull(Target)
  
  if (length(shared_proteins) == 0) {
    message("No shared CNS/Immune proteins found for fluid: ", fluid_name)
    return(NULL)
  }
  
  paired_data <- fluid_data %>%
    filter(Target %in% shared_proteins) %>%
    group_by(SampleName, Target, panel) %>%
    summarise(NPQ = mean(NPQ, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = panel, values_from = NPQ) %>%
    drop_na(CNS, IMMUNE)
  
  protein_corr <- paired_data %>%
    group_by(Target) %>%
    summarise(
      r = cor(CNS, IMMUNE, use = "pairwise.complete.obs", method = "pearson"),
      n = n(),
      .groups = "drop"
    ) %>%
    filter(!is.na(r)) %>%
    arrange(r) %>%
    mutate(
      Target = factor(Target, levels = Target),
      agreement = if_else(r < threshold, "Low agreement", "High agreement"),
      fluid = fluid_name
    )
  
  protein_corr
}

all_fluid_results <- map(fluids, ~ correlate_panels_by_fluid(common_data, .x)) %>%
  set_names(fluids) %>%
  compact()   

all_fluid_corr <- bind_rows(all_fluid_results)
write_csv(all_fluid_corr, "results/protein_CNS_vs_Immune_correlation_by_fluid.csv")

# correlation plots per fluid
plots <- imap(all_fluid_results, function(df, fluid_name) {
  
  p <- ggplot(df, aes(x = r, y = Target, fill = agreement)) +
    geom_col() +
    geom_vline(xintercept = threshold, linetype = "dashed", color = "black") +
    scale_fill_manual(values = c("Low agreement" = "#D55E00", "High agreement" = "#0072B2")) +
    labs(
      title = paste0("CNS vs Immune panel correlation — ", fluid_name),
      subtitle = paste0("Dashed line = threshold (r = ", threshold, ")"),
      x = "Pearson r (CNS vs Immune)",
      y = NULL,
      fill = NULL
    ) +
    theme_minimal(base_size = 18) +
    theme(axis.text.y = element_text(size = 15))
  
  ggsave(
    filename = paste0("plots/CNS_IMMUNE_panels/protein_correlation_", fluid_name, ".png"),
    plot = p,
    width = 5, height = 6, dpi = 300, limitsize = FALSE
  )
  
  p
})

## ------------------------------------------------------------------
## ML of both CNS and immune panel together 

## z-score wide format
zscore_proteins <- function(wide_df) {
  wide_df %>%
    mutate(across(where(is.numeric), ~ as.numeric(scale(.x))))
}

## z-score matrice
build_panel_wide <- function(data, panel_name, fluid, types, npq_col, suffix) {
  data %>%
    filter(panel == panel_name, SampleMatrixType == fluid, type %in% types) %>%
    select(SampleName, Target, NPQ = all_of(npq_col), type) %>%
    pivot_wider(names_from = Target, values_from = NPQ) %>%
    zscore_proteins() %>%
    rename_with(~ paste0(.x, suffix), .cols = -c(SampleName, type))
}

## identify shared proteins (both panels) with correlation above threshold, per fluid ----
get_high_corr_targets <- function(data, fluid, npq_col, threshold = 0.8) {
  
  fluid_data <- data %>% filter(panel %in% c("CNS", "IMMUNE"), SampleMatrixType == fluid)
  
  ## Match by UniProtID
  shared_proteins <- fluid_data %>%
    distinct(UniProtID, Target, panel) %>%
    count(UniProtID) %>%
    filter(n == 2) %>%
    pull(UniProtID)
  
  if (length(shared_proteins) == 0) return(character(0))
  
  ## Keep a lookup so we know which Target name belongs to which panel, per UniProtID
  target_lookup <- fluid_data %>%
    filter(UniProtID %in% shared_proteins) %>%
    distinct(UniProtID, Target, panel)
  
  paired <- fluid_data %>%
    filter(UniProtID %in% shared_proteins) %>%
    select(SampleName, UniProtID, panel, NPQ = all_of(npq_col)) %>%
    group_by(SampleName, UniProtID, panel) %>%
    summarise(NPQ = mean(NPQ, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = panel, values_from = NPQ) %>%
    drop_na(CNS, IMMUNE)
  
  high_corr_ids <- paired %>%
    group_by(UniProtID) %>%
    summarise(r = cor(CNS, IMMUNE, use = "pairwise.complete.obs"), .groups = "drop") %>%
    filter(!is.na(r), r >= threshold) %>%
    pull(UniProtID)
  
  target_lookup %>%
    filter(UniProtID %in% high_corr_ids) %>%
    pivot_wider(names_from = panel, values_from = Target, names_prefix = "target_") %>%
    mutate(merged_name = target_CNS)  
}

## merge datasets and average high correlating proteins
merge_panels <- function(cns_wide, immune_wide, corr_lookup) {
  
  merged <- full_join(cns_wide, immune_wide, by = c("SampleName", "type"))
  
  if (nrow(corr_lookup) == 0) return(merged)
  
  for (i in seq_len(nrow(corr_lookup))) {
    col_cns     <- paste0(corr_lookup$target_CNS[i],    "_CNS")
    col_immune  <- paste0(corr_lookup$target_IMMUNE[i], "_IMMUNE")
    merged_name <- corr_lookup$merged_name[i]   # no panel suffix 
    
    if (all(c(col_cns, col_immune) %in% names(merged))) {
      merged[[merged_name]] <- rowMeans(merged[, c(col_cns, col_immune)], na.rm = TRUE)
      merged[[col_cns]]    <- NULL
      merged[[col_immune]] <- NULL
    } else {
      message("Skipping merge for protein '", merged_name,
              "' — one or both columns not found (", col_cns, " / ", col_immune, ")")
    }
  }
  
  merged
}

## build final ML datasets 
build_ml_dataset <- function(data, fluid, types, npq_col, status_positive) {
  
  corr_lookup <- get_high_corr_targets(data, fluid, npq_col)
  
  cns_wide    <- build_panel_wide(data, "CNS",    fluid, types, npq_col, suffix = "_CNS")
  immune_wide <- build_panel_wide(data, "IMMUNE", fluid, types, npq_col, suffix = "_IMMUNE")
  
  merged <- merge_panels(cns_wide, immune_wide, corr_lookup)
  
  merged %>%
    select(-SampleName) %>%
    rename(status = type) %>%
    mutate(status = ifelse(status == status_positive, 1, 0))
}

build_ml_dataset_type <- function(data, fluid, types, npq_col, status_positive) {
  
  corr_lookup <- get_high_corr_targets(data, fluid, npq_col)
  
  cns_wide    <- build_panel_wide(data, "CNS",    fluid, types, npq_col, suffix = "_CNS")
  immune_wide <- build_panel_wide(data, "IMMUNE", fluid, types, npq_col, suffix = "_IMMUNE")
  
  merged <- merge_panels(cns_wide, immune_wide, corr_lookup)
  
  merged
}

## ------------------------------------------------------------------
## Main loop to run all the functions
for (fluid in fluids) {
  
  ## Prepare merged, z-scored, correlation-averaged datasets
  protein_data_PGMCvsCTR_new <- build_ml_dataset(all_data, fluid, c("PGMC", "CTR"), "NPQ",     "PGMC")
  protein_data_ALSvsCTR_new  <- build_ml_dataset(all_data, fluid, c("ALS",  "CTR"), "NPQ", "ALS")
  protein_data_ALSvsPGMC_new <- build_ml_dataset(all_data, fluid, c("ALS",  "PGMC"), "NPQ", "ALS")
  
  ## Run Lasso and ROC curve with 5-fold cv and 500 bootstrap iterations
  lm_PGMC_CTR <- runML_with_lasso(protein_data_PGMCvsCTR_new, bs_count = 500)
  final_roc_plot_PGMC_CTR <- calculateROC(lm_PGMC_CTR,
                                          paste0("plots/CNS_IMMUNE_panels/ML/Lasso/ROC_lm_CNS_IMMUNE_", fluid, "_PGMC_CTR.pdf"))

  lm_ALS_CTR <- runML_with_lasso(protein_data_ALSvsCTR_new, bs_count = 500)
  final_roc_plot_ALS_CTR <- calculateROC(lm_ALS_CTR,
                                         paste0("plots/CNS_IMMUNE_panels/ML/Lasso/ROC_lm_CNS_IMMUNE_", fluid, "_ALS_CTR.pdf"))

  lm_ALS_PGMC <- runML_with_lasso(protein_data_ALSvsPGMC_new, bs_count = 500)
  final_roc_plot_ALS_PGMC <- calculateROC(lm_ALS_PGMC,
                                          paste0("plots/CNS_IMMUNE_panels/ML/Lasso/ROC_lm_CNS_IMMUNE_", fluid, "_ALS_PGMC.pdf"))

  ## Extract weights + plot averaged
  lm_weights_PGMC_CTR <- feature_importance(lm_PGMC_CTR,
                                            plot_path = paste0("plots/CNS_IMMUNE_panels/ML/Lasso/protein_selection_CNS_IMMUNE_", fluid, "_PGMC_CTR.pdf"))
  write_xlsx(lm_weights_PGMC_CTR, path = paste0("results/weights_lm_CNS_IMMUNE_", fluid, "_PGMC_CTR.xlsx"))

  lm_weights_ALS_CTR <- feature_importance(lm_ALS_CTR,
                                           plot_path = paste0("plots/CNS_IMMUNE_panels/ML/Lasso/protein_selection_CNS_IMMUNE_", fluid, "_ALS_CTR.pdf"))
  write_xlsx(lm_weights_ALS_CTR, path = paste0("results/weights_lm_CNS_IMMUNE_", fluid, "_ALS_CTR.xlsx"))

  lm_weights_ALS_PGMC <- feature_importance(lm_ALS_PGMC,
                                            plot_path = paste0("plots/CNS_IMMUNE_panels/ML/Lasso/protein_selection_CNS_IMMUNE_", fluid, "_ALS_PGMC.pdf"))
  write_xlsx(lm_weights_ALS_PGMC, path = paste0("results/weights_lm_CNS_IMMUNE_", fluid, "_ALS_PGMC.xlsx"))
  
  # Run Elastic Net and ROC curve with 5-fold cv and 500 bootstrap iterations
  lm_enet_PGMC_CTR <- runML_with_elasticnet(protein_data_PGMCvsCTR_new,bs_count = 500)
  final_roc_plot_PGMC_CTR_enet = calculateROC(lm_enet_PGMC_CTR,
                                              paste0("plots/CNS_IMMUNE_panels/ML/Elastic Net/ROC_lm_CNS_IMMUNE_",fluid,"_PGMC_CTR.pdf"))
  lm_enet_ALS_CTR <- runML_with_elasticnet(protein_data_ALSvsCTR_new,bs_count = 500)
  final_roc_plot_ALS_CTR_enet = calculateROC(lm_enet_ALS_CTR,
                                             paste0("plots/CNS_IMMUNE_panels/ML/Elastic Net/ROC_lm_CNS_IMMUNE_",fluid,"_ALS_CTR.pdf"))
  lm_enet_ALS_PGMC <- runML_with_elasticnet(protein_data_ALSvsPGMC_new,bs_count = 500)
  final_roc_plot_ALS_PGMC_enet = calculateROC(lm_enet_ALS_PGMC,
                                              paste0("plots/CNS_IMMUNE_panels/ML/Elastic Net/ROC_lm_CNS_IMMUNE_",fluid,"_ALS_PGMC.pdf"))
  
  # extract weights + plot averaged
  lm_weights_PGMC_CTR_enet = feature_importance(lm_enet_PGMC_CTR,
                                                plot_path = paste0("plots/CNS_IMMUNE_panels/ML/Elastic Net/protein_selection_CNS_IMMUNE_",fluid,"_PGMC_CTR.pdf"))
  write_xlsx(lm_weights_PGMC_CTR_enet, path = paste0("results/weights_lm_enet_CNS_IMMUNE_",fluid,"_PGMC_CTR.xlsx"))
  lm_weights_ALS_CTR_enet = feature_importance(lm_enet_ALS_CTR,
                                               plot_path = paste0("plots/CNS_IMMUNE_panels/ML/Elastic Net/protein_selection_CNS_IMMUNE_",fluid,"_ALS_CTR.pdf"))
  write_xlsx(lm_weights_ALS_CTR_enet, path = paste0("results/weights_lm_enet_CNS_IMMUNE_",fluid,"_ALS_CTR.xlsx"))
  lm_weights_ALS_PGMC_enet = feature_importance(lm_enet_ALS_PGMC,
                                                plot_path = paste0("plots/CNS_IMMUNE_panels/ML/Elastic Net/protein_selection_CNS_IMMUNE_",fluid,"_ALS_PGMC.pdf"))
  write_xlsx(lm_weights_ALS_PGMC_enet, path = paste0("results/weights_lm_enet_CNS_IMMUNE_",fluid,"_ALS_PGMC.xlsx"))
}

# =============================================================
# Find the most optimal set of proteins for each ML and fluid
find_optimal_signature_CNS_IMMUNE <- function(data_all, 
                                   ranked_proteins, 
                                   fluid = "SERUM",
                                   output_prefix = "optimization_",
                                   inner_folds = 5,
                                   seed = 123) {
  
  set.seed(seed)
  
  min_proteins = 3
  
  X <- as.matrix(data_all[, ranked_proteins])
  y <- data_all$status
  
  X <- scale(X)
  y_factor <- as.factor(y)
  
  # performance curve
  auc_list <- vector("list", length(ranked_proteins))
  folds_inner_full <- createFolds(y_factor, k = inner_folds)
  
  for (k in 1:length(ranked_proteins)) {
    
    vars <- ranked_proteins[1:k]
    auc_inner <- c()
    
    for (j in seq_along(folds_inner_full)) {
      
      val_idx <- folds_inner_full[[j]]
      tr_idx <- setdiff(seq_along(y), val_idx)
      
      X_tr <- X[tr_idx, vars, drop = FALSE]
      y_tr <- y[tr_idx]
      
      X_val <- X[val_idx, vars, drop = FALSE]
      y_val <- y[val_idx]
      
      if (length(unique(y_tr)) < 2 || length(unique(y_val)) < 2) next
      
      if (k == 1) {
        df_tr <- data.frame(y = y_tr, x = X_tr[, 1])
        df_val <- data.frame(x = X_val[, 1])
        fit <- glm(y ~ x, data = df_tr, family = "binomial")
        probs <- predict(fit, newdata = df_val, type = "response")
      } else {
        fit <- cv.glmnet(X_tr, y_tr,
                         family = "binomial",
                         alpha = 0,
                         nfolds = inner_folds)
        probs <- predict(fit, newx = X_val, s = "lambda.min", type = "response")
      }
      
      auc_inner <- c(auc_inner,
                     as.numeric(auc(roc(y_val, probs, quiet = TRUE,
                                        levels=c(0,1), direction="<"))))
    }
    
    auc_list[[k]] <- auc_inner
  }
  
  perf_results <- data.frame(
    n_proteins = 1:length(ranked_proteins),
    protein_added = ranked_proteins,
    auc_mean = sapply(auc_list, function(x) mean(x, na.rm = TRUE)),
    auc_sd   = sapply(auc_list, function(x) sd(x, na.rm = TRUE))
  )
  
  # selection of best number of proteins
  valid_idx <- min_proteins:length(ranked_proteins)
  
  best_idx <- which.max(perf_results$auc_mean[valid_idx]) + (min_proteins - 1)
  best_auc <- perf_results$auc_mean[best_idx]
  best_sd  <- perf_results$auc_sd[best_idx]
  threshold <- best_auc - 0.5 * (best_sd/sqrt(inner_folds))
  candidate_k <- valid_idx[perf_results$auc_mean[valid_idx] >= threshold]
  
  best_n <- min(candidate_k) 
  optimal_proteins <- ranked_proteins[1:best_n]
  selected_auc <- perf_results$auc_mean[best_n]
  
  message("Selected signature size: ", best_n)
  
  # plot
  selected_auc <- perf_results$auc_mean[perf_results$n_proteins == best_n]
  
  p <- ggplot(perf_results, aes(x = n_proteins, y = auc_mean)) +
    
    # SD error bars
    geom_errorbar(aes(ymin = auc_mean - auc_sd/sqrt(inner_folds),
                      ymax = auc_mean + auc_sd/sqrt(inner_folds)),
                  width = 0.2,
                  color = "grey50") +
    
    # mean line
    geom_line(color = "#ad5291", size = 1) +
    
    # selected (robust) signature
    geom_point(aes(color = (n_proteins == best_n)), size = 3) +
    geom_vline(xintercept = best_n, linetype = "dashed", color = "gray50") +
    
    # best AUC point
    #geom_point(data = perf_results[best_idx, ],
    #          aes(x = n_proteins, y = auc_mean),
    #               color = "blue", size = 3) +
    
    scale_x_continuous(breaks = 1:length(ranked_proteins),
                       labels = ranked_proteins) +
    
    scale_color_manual(values = c("black", "red"), guide = "none") +
    
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    
    labs(
      title = "Selection of best protein signature based on performance",
      subtitle = paste0(
        "Selected signature (red): ", best_n, " proteins with ",
        "AUC = ", round(selected_auc, 3), " (± ", 
        round(perf_results$auc_sd[best_n]/sqrt(inner_folds), 3), ")."),
      x = "Proteins added based on feature importance",
      y = "Cross-validated AUC (mean ± SE)"
    )
  
  ggsave(paste0(output_prefix, "protein_signature_", fluid, ".pdf"),
         p, width = 8, height = 6)
  
  return(list(
    optimal_n = best_n,
    optimal_proteins = optimal_proteins,
    performance_table = perf_results,
    plot = p
  ))
}

run_ALS_signature_workflow_CNS_IMMUNE <- function(
    data_all,
    proteins,
    fluid = "SERUM",
    model_type = c("lasso", "elastic_net"),
    alpha = NULL,
    output_prefix = "output",
    make_plots = TRUE
) {
  
  model_type <- match.arg(model_type)
  
  # define alpha
  if (is.null(alpha)) {
    alpha <- ifelse(model_type == "lasso", 1, 0.5)
  }
  alpha = 0
  
  message("Running model: ", model_type, " (alpha = ", alpha, ") on ", fluid)
  
  # prepare data
  data_train <- build_ml_dataset_type(all_data, 
                                 fluid, 
                                 c("ALS",  "CTR"), 
                                 "NPQ", "ALS")
  
  data_model <- data_train %>%
    select(-SampleName) %>%
    rename(status = type) %>%
    mutate(status = ifelse(status == "ALS", 1, 0))

  X <- as.matrix(data_model[, proteins])
  y <- data_model$status
  
  # scale
  X_scaled <- scale(X)
  scaling_center <- attr(X_scaled, "scaled:center")
  scaling_scale  <- attr(X_scaled, "scaled:scale")
  
  # model with protein signature
  set.seed(123)
  cv_model <- cv.glmnet(X_scaled, y, family = "binomial", alpha = alpha)
  
  coef_df <- as.data.frame(as.matrix(coef(cv_model, s = "lambda.min")))
  coef_df$feature <- rownames(coef_df)
  
  # prediction
  pred_train <- predict(cv_model,
                        newx = X_scaled,
                        s = "lambda.min",
                        type = "link")
  
  ALS_CTR_results <- data_train %>%
    mutate(ALS_risk_score = as.numeric(pred_train)) %>%
    select(-type)
  
  # check prediction in PGMCs
  data_pgmc <- build_ml_dataset_type(all_data, 
                                fluid, 
                                c("PGMC"), 
                                "NPQ", "PGMC") 
  
  X_pgmc <- as.matrix(data_pgmc[, proteins])
  X_pgmc_scaled <- scale(X_pgmc,
                         center = scaling_center,
                         scale = scaling_scale)
  
  pred_pgmc <- predict(cv_model,
                       newx = X_pgmc_scaled,
                       s = "lambda.min",
                       type = "link")
  
  PGMC_results <- data_pgmc %>%
    mutate(ALS_risk_score = as.numeric(pred_pgmc)) %>%
    select(-type)
  
  # combine both
  results_all <- bind_rows(PGMC_results, ALS_CTR_results)
  
  if (make_plots) {
    
    plot_data <- results_all %>%
      left_join(samples_ID_type %>% rename(SampleName = `Sample ID`)) %>%
      mutate(type = factor(type, levels = c("ALS","PGMC","CTR")),
             label = ifelse(type == "PGMC" & ALS_risk_score > 0,
                            ParticipantCode, NA))
    
    group_colors <- c("CTR"="#6F8EB2","ALS"="#B2936F","PGMC"="#ad5291")
    
    pdf(paste0(output_prefix, "ALS_risk_score_",fluid,".pdf"), height = 5, width = 7)
    
    print(
      ggplot(plot_data, aes(x = ALS_risk_score, y = type, fill = type)) +
        geom_violin(trim = FALSE, alpha = 0.4) +
        geom_jitter(aes(color = type), height = 0.08, size = 1.2) +
        ggrepel::geom_text_repel(
          data = plot_data %>% filter(!is.na(label)),
          aes(label = label),
          size = 3,
          segment.color = "black",
          max.overlaps = 50,
          direction = "y") +
        geom_vline(xintercept = 0, linetype = "dashed") +
        scale_fill_manual(values = group_colors) +
        scale_color_manual(values = group_colors) +
        theme_classic(base_size = 12) +
        theme(
          axis.line = element_line(size = 0.4),
          axis.ticks = element_line(size = 0.3),
          axis.text = element_text(color = "black"),
          axis.title = element_text(size = 12),
          plot.title = element_text(size = 13, face = "bold"),
          legend.position = "none"
        ) +
        labs(
          x = "ALS Risk Score",
          y = NULL,
          title = paste0("ALS risk score (", model_type, ", ", fluid, ")")))
    
    dev.off()
  }
  
  return(list(
    results = results_all,
    model = cv_model,
    coefficients = coef_df,
    scaling = list(center = scaling_center, scale = scaling_scale)
  ))
}

run_heatmap_signature_CNS_IMMUNE <- function(
    data_all,
    fluid,
    proteins,
    risk_scores_df,   
    group_colors,
    output_prefix = "heatmap",
    highlight_ids = NULL
) {
  
  message("Running heatmap for ", fluid)
  
  # heatmap data
  data_heatmap <- build_ml_dataset_type(data_all,
                                        fluid, c("PGMC","ALS","CTR"),
                                        "NPQ","ALS") %>%
    left_join(samples_ID_type %>% rename(SampleName = `Sample ID`))  %>%
    left_join(Sex_age_all_participants %>% rename(PatientID = Pseudonyme)) %>%
    left_join(ALSFRS_NULISA_V0 %>%
                filter(type == "ALS") %>%
                select(PatientID, ALSFRS_1),
              by = "PatientID") %>%
    left_join(site_onset %>%
                filter(type == "ALS") %>%
                select(ParticipantCode, SiteDiseaseOnset),
              by = "ParticipantCode") %>%
    left_join(disease_duration %>%
                filter(type == "ALS", Visit == "V0") %>%
                select(PatientID, `Disease duration`) %>%
                distinct(),
              by = "PatientID") %>%
    left_join(clinical_table_extended %>% select(ParticipantCode,MutationType)) %>%
    distinct() %>%
    select(proteins, ParticipantCode, type, sex, age,
           ALSFRS_1, SiteDiseaseOnset, `Disease duration`,SampleName,MutationType) %>%
    left_join(risk_scores_df %>%
                select(SampleName, ALS_risk_score),
              by = "SampleName")
  
  # heatmap matrix
  mat <- data_heatmap %>%
    select(-c(ParticipantCode, type, ALS_risk_score,
              sex, age, ALSFRS_1,
              SiteDiseaseOnset, `Disease duration`,SampleName,MutationType)) %>%
    as.matrix()
  
  rownames(mat) <- data_heatmap$ParticipantCode
  
  # scale rows (z-score per sample)
  mat_scaled <- t(scale(t(mat)))
  
  # definition of colors
  score_col <- colorRamp2(
    c(min(data_heatmap$ALS_risk_score, na.rm = TRUE),
      0,
      max(data_heatmap$ALS_risk_score, na.rm = TRUE)),
    c("lightgreen", "white", "darkred")
  )
  
  disease_duration_col <- colorRamp2(
    range(data_heatmap$`Disease duration`, na.rm = TRUE),
    c("#eae115","#e69619")
  )
  
  age_col <- colorRamp2(
    range(data_heatmap$age, na.rm = TRUE),
    c("lightgreen","darkgreen")
  )
  
  ALSFRS_col <- colorRamp2(
    c(0,22,48),
    c("#0c1cf3", "#d1d1f0", "#f4f5f9")
  )
  
  # row annotations
  ha <- rowAnnotation(
    Group = data_heatmap$type,
    ALS_score = data_heatmap$ALS_risk_score,
    IDs_interest = ifelse(data_heatmap$ParticipantCode %in% highlight_ids,
                          "yes","no"),
    Disease_onset = data_heatmap$SiteDiseaseOnset,
    Disease_duration = data_heatmap$`Disease duration`,
    Mutation = data_heatmap$MutationType,
    ALSFRS_R = data_heatmap$ALSFRS_1,
    Age = data_heatmap$age,
    Sex = data_heatmap$sex,
    na_col = "white",
    col = list(
      Group = group_colors,
      ALS_score = score_col,
      Disease_onset = c("spinal"="#3498db",
                        "bulbar"="#e74c3c",
                        "other"="grey"),
      Mutation = c("C9orf72" = "#66c2a5",
                   "SOD1" = "#fc8d62",
                   "TARDBP" = "#8da0cb",
                   "FUS" = "#e78ac3",
                   "FIG4" = "#a6d854",
                   "Other" = "#b3b3b3",
                   "UBQLN2" = "#ffd92f",
                   "ANG" = "#e5c494"),
      Disease_duration = disease_duration_col,
      ALSFRS_R = ALSFRS_col,
      Age = age_col,
      Sex = c("M"="#4c72b0","F"="#dd8452"),
      IDs_interest = c("yes"="black","no"="white")
    )
  )
  
  # clustered
  pdf(paste0(output_prefix, "_clustered.pdf"), height = 14)
  print(
    Heatmap(
      mat_scaled,
      name = "Z-score",
      left_annotation = ha,
      cluster_rows = TRUE,
      cluster_columns = TRUE,
      show_row_names = TRUE,
      show_column_names = TRUE,
      col = colorRamp2(c(-2,0,2), c("blue","white","red"))
    )
  )
  dev.off()
  
  # grouped
  pdf(paste0(output_prefix, "_grouped.pdf"), height = 14)
  print(
    Heatmap(
      mat_scaled,
      name = "Z-score",
      left_annotation = ha,
      row_split = data_heatmap$type,
      cluster_rows = TRUE,
      cluster_columns = TRUE,
      show_row_names = TRUE,
      show_column_names = TRUE,
      col = colorRamp2(c(-2,0,2), c("blue","white","red"))
    )
  )
  dev.off()
  
  return(data_heatmap)
}

## -> Lasso + Serum
proteins_serum_ALS_CTR = c("NEFL_CNS","pTau-181_CNS","CLEC4A_IMMUNE","HAVCR1_IMMUNE",
                           "BASP1_CNS","ACHE_CNS","TAFA5_CNS","pTau-231_CNS","SELE_IMMUNE")
lasso_optimal_protein_signature_serum = find_optimal_signature_CNS_IMMUNE(data_all = build_ml_dataset(all_data, 
                                                                                                      "SERUM", 
                                                                                                      c("ALS",  "CTR"), 
                                                                                                      "NPQ", "ALS"),
                                                               ranked_proteins = proteins_serum_ALS_CTR,
                                                               fluid = "SERUM",
                                                               output_prefix = "plots/CNS_IMMUNE_panels/ML/Lasso/performance_")

# Elastic Net + Serum
proteins_serum_ALS_CTR = c("NEFL_CNS","pTau-181_CNS","NEFH_CNS","pTau-231_CNS",
                           "CLEC4A_IMMUNE","HAVCR1_IMMUNE","FABP3_CNS","GDNF_CNS",
                           "SELE_IMMUNE","TNFRSF13C_IMMUNE","IL1RL1_IMMUNE","TAFA5_CNS")
enet_optimal_protein_signature_serum = find_optimal_signature_CNS_IMMUNE(data_all = build_ml_dataset(all_data, 
                                                                                                     "SERUM", 
                                                                                                     c("ALS",  "CTR"), 
                                                                                                     "NPQ", "ALS"),
                                                              ranked_proteins = proteins_serum_ALS_CTR,
                                                              fluid = "SERUM",
                                                              output_prefix = "plots/CNS_IMMUNE_panels/ML/Elastic Net/performance_")

# Elastic Net + Plasma
proteins_plasma_ALS_CTR = c("NEFL_CNS","NEFH_CNS","pTau-181_CNS","FABP3_CNS",
                            "pTau-231_CNS","TNFRSF9_IMMUNE","pTau-217_CNS",
                            "CD3E_IMMUNE","SCG2_IMMUNE","SFRP1_CNS","IL1B_CNS")
enet_optimal_protein_signature_plasma = find_optimal_signature_CNS_IMMUNE(data_all = build_ml_dataset(all_data, 
                                                                                                      "PLASMA", 
                                                                                                      c("ALS",  "CTR"), 
                                                                                                      "NPQ", "ALS"),
                                                               ranked_proteins = proteins_plasma_ALS_CTR,
                                                               fluid = "PLASMA",
                                                               output_prefix = "plots/CNS_IMMUNE_panels/ML/Elastic Net/performance_")

# Elastic Net + CSF
proteins_CSF_ALS_CTR =  c("NEFL_CNS","NEFH_CNS","CHI3L1","IL33_IMMUNE","CHIT1_CNS",
                          "CCL2","CCL3_IMMUNE","IFNG_IMMUNE","IL18_CNS","CCL3_CNS",
                          "IL23_IMMUNE")

enet_optimal_protein_signature_CSF = find_optimal_signature_CNS_IMMUNE(data_all = build_ml_dataset(all_data, 
                                                                                                   "CSF", 
                                                                                                   c("ALS",  "CTR"), 
                                                                                                   "NPQ", "ALS"),
                                                            ranked_proteins = proteins_CSF_ALS_CTR,
                                                            fluid = "CSF",
                                                            output_prefix = "plots/CNS_IMMUNE_panels/ML/Elastic Net/performance_")


# ====================================================
# Compute an ALS risk score based on ALS vs CTR model 
## -> Lasso + Serum
proteins_serum_ALS_CTR = lasso_optimal_protein_signature_serum$optimal_proteins

lasso_PGMC_serum_signature = run_ALS_signature_workflow_CNS_IMMUNE(all_data,
                                                        proteins_serum_ALS_CTR,
                                                        fluid = "SERUM",
                                                        model_type = "lasso",
                                                        output_prefix = "plots/CNS_IMMUNE_panels/ML/Lasso/")

participant_code_label = lasso_PGMC_serum_signature$results %>%
  left_join(protein_data_IDs %>% select(SampleName,ParticipantCode,type)) %>%
  distinct() %>%
  filter(type == "PGMC" & ALS_risk_score > 0) %>%
  pull(ParticipantCode)

# ================================================================================
# Unsupervised visualisation (heatmap) of PGMC, ALS, CTR based on 4-protein signature

group_colors <- c(
  "CTR"  = "#6F8EB2",
  "ALS"  = "#B2936F",
  "PGMC" = "#ad5291")

run_heatmap_signature_CNS_IMMUNE(all_data,
                      "SERUM",
                      lasso_optimal_protein_signature_serum$optimal_proteins,
                      lasso_PGMC_serum_signature$results,
                      group_colors = group_colors,
                      output_prefix = "plots/CNS_IMMUNE_panels/ML/Lasso/heatmap_SERUM",
                      highlight_ids = participant_code_label)

##### ----
# Elastic Net

## -> Elastic Net + Serum
EN_PGMC_serum_signature = run_ALS_signature_workflow_CNS_IMMUNE(all_data,
                                                     enet_optimal_protein_signature_serum$optimal_proteins,
                                                     fluid = "SERUM",
                                                     model_type = "elastic_net",
                                                     output_prefix = "plots/CNS_IMMUNE_panels/ML/Elastic Net/")

participant_code_label = EN_PGMC_serum_signature$results %>%
  left_join(protein_data_IDs %>% select(SampleName,ParticipantCode,type)) %>%
  distinct() %>%
  filter(type == "PGMC" & ALS_risk_score > 0) %>%
  pull(ParticipantCode)

# ================================================================================
# Unsupervised visualisation (heatmap) of PGMC, ALS, CTR based on 9-protein signature

run_heatmap_signature_CNS_IMMUNE(all_data,
                      "SERUM",
                      enet_optimal_protein_signature_serum$optimal_proteins,
                      EN_PGMC_serum_signature$results,
                      group_colors = group_colors,
                      output_prefix = "plots/CNS_IMMUNE_panels/ML/Elastic Net/heatmap_SERUM",
                      highlight_ids = participant_code_label)

## -> Elastic Net + Plasma
EN_PGMC_plasma_signature = run_ALS_signature_workflow_CNS_IMMUNE(all_data,
                                                      enet_optimal_protein_signature_plasma$optimal_proteins,
                                                      fluid = "PLASMA",
                                                      model_type = "elastic_net",
                                                      output_prefix = "plots/CNS_IMMUNE_panels/ML/Elastic Net/")


participant_code_label = EN_PGMC_plasma_signature$results %>%
  left_join(protein_data_IDs %>% select(SampleName,ParticipantCode,type)) %>%
  distinct() %>%
  filter(type == "PGMC" & ALS_risk_score > 0) %>%
  pull(ParticipantCode)

# ================================================================================
# Unsupervised visualisation (heatmap) of PGMC, ALS, CTR based on 3-protein signature

run_heatmap_signature_CNS_IMMUNE(all_data,
                      "PLASMA",
                      enet_optimal_protein_signature_plasma$optimal_proteins,
                      EN_PGMC_plasma_signature$results,
                      group_colors = group_colors,
                      output_prefix = "plots/CNS_IMMUNE_panels/ML/Elastic Net/heatmap_PLASMA",
                      highlight_ids = participant_code_label)

## -> Elastic Net + CSF
EN_PGMC_CSF_signature = run_ALS_signature_workflow_CNS_IMMUNE(all_data,
                                                   enet_optimal_protein_signature_CSF$optimal_proteins,
                                                   fluid = "CSF",
                                                   model_type = "elastic_net",
                                                   output_prefix = "plots/CNS_IMMUNE_panels/ML/Elastic Net/")

participant_code_label = EN_PGMC_CSF_signature$results %>%
  left_join(protein_data_IDs %>% select(SampleName,ParticipantCode,type)) %>%
  distinct() %>%
  filter(type == "PGMC" & ALS_risk_score > 0) %>%
  pull(ParticipantCode)

# ================================================================================
# Unsupervised visualisation (heatmap) of PGMC, ALS, CTR based on 9-protein signature

run_heatmap_signature_CNS_IMMUNE(all_data,
                      "CSF",
                      enet_optimal_protein_signature_CSF$optimal_proteins,
                      EN_PGMC_CSF_signature$results,
                      group_colors = group_colors,
                      output_prefix = "plots/CNS_IMMUNE_panels/ML/Elastic Net/heatmap_CSF",
                      highlight_ids = participant_code_label)

