###############################################
### Helper Functions
###############################################

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

# make covariate table for MAXOMOD data
build_external_covariate_table <- function(ids) {
  
  dropped <- ids %>% filter(!disease %in% c("als", "ctrl"))
  if (nrow(dropped) > 0) {
    message("Dropping ", nrow(dropped), " samples with disease label(s): ",
            paste(unique(dropped$disease), collapse = ", "),
            " (not als/ctrl)")
  }
  
  ids %>%
    filter(disease %in% c("als", "ctrl")) %>%
    distinct(SampleName = Tube_ID, age, sex, disease) %>%
    mutate(
      sex_male = ifelse(sex == "Male", 1, 0),   # matches internal coding; flag if external sex has other levels
      status   = ifelse(disease == "als", 1, 0)
    ) %>%
    select(SampleName, age, sex_male, status)
}

## z-score wide format
zscore_proteins <- function(wide_df, protein_cols, center = NULL, scale = NULL) {
  
  mat <- as.matrix(wide_df[, protein_cols, drop = FALSE])
  
  if (is.null(center) || is.null(scale)) {
    ## discovery mode: compute fresh from this data, and expose the params
    scaled_mat <- scale(mat)
    center <- attr(scaled_mat, "scaled:center")
    scale  <- attr(scaled_mat, "scaled:scale")
  } else {
    ## validation mode: reuse frozen training params -- never recompute
    center <- center[protein_cols]
    scale  <- scale[protein_cols]
    scaled_mat <- scale(mat, center = center, scale = scale)
  }
  
  wide_df[, protein_cols] <- as.data.frame(scaled_mat)
  attr(wide_df, "zscore_center") <- center
  attr(wide_df, "zscore_scale")  <- scale
  wide_df
}

## mak data for ML analysis
build_panel_wide <- function(data, panel_name, fluid, types, npq_col, suffix,
                             center = NULL, scale = NULL) {
  
  wide <- data %>%
    filter(panel == panel_name, SampleMatrixType == fluid, type %in% types) %>%
    select(SampleName, Target, NPQ = all_of(npq_col), type) %>%
    pivot_wider(names_from = Target, values_from = NPQ)
  
  protein_cols <- setdiff(names(wide), c("SampleName", "type"))
  
  wide <- zscore_proteins(wide, protein_cols, center = center, scale = scale)
  
  ## capture center and scale
  panel_center <- attr(wide, "zscore_center")
  panel_scale  <- attr(wide, "zscore_scale")
  
  wide <- wide %>% rename_with(~ paste0(.x, suffix), .cols = -c(SampleName, type))
  
  attr(wide, "zscore_center") <- panel_center
  attr(wide, "zscore_scale")  <- panel_scale
  wide
}

## identify shared proteins (both panels) with correlation above threshold, per fluid 
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

# make covariate table for ML
build_covariate_table <- function(metadata) {
  metadata %>%
    distinct(SampleName, age, sex, center) %>%   
    mutate(
      sex_male       = ifelse(sex == "M", 1, 0),  # female as reference
      center_Turkey = ifelse(center == "Turkey", 1, 0) # Germany = reference (0)
    ) %>%
    select(SampleName, age, sex_male, center_Turkey)
}

## build final ML datasets 
build_ml_dataset <- function(data, fluid, types, npq_col, status_positive,
                             covariate_table,
                             panel_mode = "both") {
  
  if (panel_mode == "both") {
    corr_lookup <- get_high_corr_targets(data, fluid, npq_col)
    cns_wide    <- build_panel_wide(data, "CNS",    fluid, types, npq_col, suffix = "_CNS")
    immune_wide <- build_panel_wide(data, "IMMUNE", fluid, types, npq_col, suffix = "_IMMUNE")
    merged <- merge_panels(cns_wide, immune_wide, corr_lookup)
  } else {
    merged <- build_panel_wide(data, panel_mode, fluid, types, npq_col, suffix = paste0("_", panel_mode))
  }
  
  merged <- merged %>% left_join(covariate_table, by = "SampleName")
  
  merged %>%
    select(-SampleName) %>%
    rename(status = type) %>%
    mutate(status = ifelse(status == status_positive, 1, 0))
}

build_ml_dataset_type <- function(data, fluid, types, npq_col, status_positive,
                                  panel_zscore_params = NULL,
                                  covariate_table,
                                  panel_mode = "both") {
  
  if (panel_mode == "both") {
    corr_lookup <- get_high_corr_targets(data, fluid, npq_col)
    
    cns_wide <- build_panel_wide(data, "CNS", fluid, types, npq_col, suffix = "_CNS",
                                 center = panel_zscore_params$CNS$center,
                                 scale  = panel_zscore_params$CNS$scale)
    
    immune_wide <- build_panel_wide(data, "IMMUNE", fluid, types, npq_col, suffix = "_IMMUNE",
                                    center = panel_zscore_params$IMMUNE$center,
                                    scale  = panel_zscore_params$IMMUNE$scale)
    
    merged <- merge_panels(cns_wide, immune_wide, corr_lookup)
    zparams <- list(
      CNS    = list(center = attr(cns_wide, "zscore_center"),    scale = attr(cns_wide, "zscore_scale")),
      IMMUNE = list(center = attr(immune_wide, "zscore_center"), scale = attr(immune_wide, "zscore_scale"))
    )
  } else {
    ## single-panel mode (only for CNS or immune)
    corr_lookup <- tibble(target_CNS = character(0), target_IMMUNE = character(0), merged_name = character(0))
    suffix <- paste0("_", panel_mode)
    
    panel_wide <- build_panel_wide(data, panel_mode, fluid, types, npq_col, suffix = suffix,
                                   center = panel_zscore_params[[panel_mode]]$center,
                                   scale  = panel_zscore_params[[panel_mode]]$scale)
    merged <- panel_wide
    zparams <- setNames(
      list(list(center = attr(panel_wide, "zscore_center"), scale = attr(panel_wide, "zscore_scale"))),
      panel_mode
    )
  }
  
  merged <- merged %>% left_join(covariate_table, by = "SampleName")
  
  attr(merged, "corr_lookup") <- corr_lookup
  attr(merged, "panel_zscore_params") <- zparams
  merged
}

# for the MAXOMOD data
build_external_ml_dataset <- function(external_long, fluid, npq_col,
                                      covariate_table, corr_lookup_discovery,
                                      panel_zscore_params,
                                      panel_mode = "both") {
  
  external_tagged <- external_long %>% mutate(type = "EXTERNAL")
  
  if (panel_mode == "both") {
    cns_wide <- build_panel_wide(external_tagged, "CNS", fluid, "EXTERNAL", npq_col, suffix = "_CNS",
                                 center = panel_zscore_params$CNS$center,
                                 scale  = panel_zscore_params$CNS$scale)
    immune_wide <- build_panel_wide(external_tagged, "IMMUNE", fluid, "EXTERNAL", npq_col, suffix = "_IMMUNE",
                                    center = panel_zscore_params$IMMUNE$center,
                                    scale  = panel_zscore_params$IMMUNE$scale)
    merged <- merge_panels(cns_wide, immune_wide, corr_lookup_discovery)
  } else {
    suffix <- paste0("_", panel_mode)
    merged <- build_panel_wide(external_tagged, panel_mode, fluid, "EXTERNAL", npq_col, suffix = suffix,
                               center = panel_zscore_params[[panel_mode]]$center,
                               scale  = panel_zscore_params[[panel_mode]]$scale)
  }
  
  merged %>% inner_join(covariate_table, by = "SampleName")
}

validate_external <- function(cv_model, external_wide, proteins, covariate_cols,
                              scaling_center, scaling_scale, fluid_label = "",
                              plot_path = NULL) {
  
  missing_cols <- setdiff(c(proteins, covariate_cols), names(external_wide))
  if (length(missing_cols) > 0) {
    stop("External dataset is missing required column(s): ",
         paste(missing_cols, collapse = ", "),
         "\n(check panel coverage / merge mapping above before proceeding)")
  }
  
  X_ext <- as.matrix(external_wide[, c(proteins, covariate_cols)])
  
  ## apply the scales from premodiALS
  X_ext_scaled <- X_ext
  X_ext_scaled[, proteins] <- scale(X_ext[, proteins],
                                    center = scaling_center,
                                    scale  = scaling_scale)
  beta <- coef(cv_model, s = "lambda.min")
  
  n_na <- sum(!complete.cases(X_ext_scaled))
  if (n_na > 0) message(n_na, " external samples have missing protein/covariate values and will be dropped from AUC")
  
  keep <- complete.cases(X_ext_scaled)
  X_ext_scaled <- X_ext_scaled[keep, , drop = FALSE]
  y_ext        <- external_wide$status[keep]
  
  probs <- predict(cv_model, newx = X_ext_scaled, s = "lambda.min", type = "response")
  
  roc_obj <- pROC::roc(y_ext, as.numeric(probs), levels = c(0, 1), direction = "<", quiet = TRUE)
  auc_val <- as.numeric(pROC::auc(roc_obj))
  ci_val  <- as.numeric(pROC::ci.auc(roc_obj))
  
  message("External validation AUC (", fluid_label, "): ",
          round(auc_val, 3), " (95% CI: ", round(ci_val[1], 3), "-", round(ci_val[3], 3), ")")
  
  ## ------------------------------------------------------------------
  ## ROC plot, styled consistently as calculateROC()
  p <- NULL
  if (!is.null(plot_path)) {
    
    roc_df <- data.frame(
      fpr = 1 - roc_obj$specificities,
      tpr = roc_obj$sensitivities
    )
    roc_df <- roc_df[order(roc_df$fpr, roc_df$tpr), ]
    
    auc_text <- sprintf("AUC = %.3f (95%% CI: %.3f\u2013%.3f)",
                        auc_val, ci_val[1], ci_val[3])
    
    p <- ggplot(roc_df, aes(x = fpr, y = tpr)) +
      geom_line(color = "darkblue", linewidth = 1.2) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                  color = "darkgrey", linewidth = 0.6) +
      scale_x_continuous(limits = c(0, 1), expand = c(0.01, 0.01)) +
      scale_y_continuous(limits = c(0, 1), expand = c(0.01, 0.01)) +
      coord_equal() +
      theme_minimal(base_size = 13) +
      theme(
        panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA),
        panel.grid.major = element_line(color = "grey90"),
        axis.line  = element_line(color = "grey30"),
        axis.ticks = element_line(color = "grey30")
      ) +
      labs(
        x = "1 - Specificity",
        y = "Sensitivity",
        title = paste0("External validation ROC curve (", fluid_label, ")"),
        subtitle = paste0(auc_text, "  |  n = ", nrow(X_ext_scaled))
      )
    
    ggsave(filename = plot_path, plot = p, width = 7, height = 7, dpi = 300)
    message("Saved external ROC plot to: ", plot_path)
  }
  
  list(
    auc = auc_val,
    auc_ci = ci_val,
    roc = roc_obj,
    plot = p,
    n = nrow(X_ext_scaled),
    predictions = data.frame(SampleName = external_wide$SampleName[keep],
                             status = y_ext,
                             risk_score = as.numeric(probs))
  )
}

# Z-score scaling
scale_manual <- function(df) {
  status_col <- df$status
  numeric_df <- df %>% dplyr::select(-status)
  scaled_df <- as.data.frame(apply(numeric_df, 2, function(x) (x - mean(x)) / sd(x)))
  scaled_df$status <- status_col
  return(scaled_df)
}

# =============================
# Main ML function with robust bootstrap + CV
# =============================
runML = function(data_frame,algorithm, cv = 10, BS_number = 100, seed = 123){
  
  set.seed(seed)
  
  #check input
  stopifnot('Last argument (BS_number) is not a number.' = is.numeric(BS_number),
            'Number given for cross validation is not a number ,' = is.numeric(cv),
            'Data is not a dataframe.' = is.data.frame(data_frame),
            'status column of df not correct' = all(data_frame$status == 1 | data_frame$status == 0),
            'unknown algorithm. please select \'lm\',\'enet\',\'rf\',\'svm rad\' or \'svm lin\' ' = algorithm %in% c('lm','enet','rf','svm rad','svm lin'))
  
  key_output = c('linear regression (lasso)',
                 'elastic net',
                 'random forest', 
                 'Support Vector Machine (radial kernel)', 
                 'Support Vector Machine (linear kernel)')
  names(key_output) = c('lm','enet','rf','svm rad','svm lin')
  # Output parameters for user check
  print(paste0('Running a ',BS_number,'x boot strap with a ',cv,
               ' fold cross validation. Algorithm is ', key_output[[algorithm]]))
  
  data_frame$status = as.factor(make.names(data_frame$status)) # caret needs prediction variable to have name; turns 0 -> X0 and 1 -> X1
  covariate_cols <- intersect(c("age", "sex_male", "center_Turkey"), names(data_frame))
  predictor_cols <- setdiff(names(data_frame), c("status"))
  penalty_vec <- ifelse(predictor_cols %in% covariate_cols, 0, 1)
  names(penalty_vec) <- predictor_cols
  
  
  # create lists to save results from bs
  models = list()
  importance = list()
  predictions = list()
  indices = list()
  predictions_raw = list()
  actuals_list = list()
  smp_size = floor(0.8 * nrow(data_frame)) # use 80% of data for training, 20% for testing
  
  #set the right parameters for each algorithm
  if(algorithm == 'lm'| algorithm == 'enet'| algorithm == 'svm rad' | algorithm == 'svm lin'){
    ctrl = trainControl(method = "cv",
                        number = cv,
                        classProbs = TRUE,
                        summaryFunction = twoClassSummary,
                        savePredictions = TRUE)
  }
  else if(algorithm == 'rf'){
    
    n_features = ncol(data_frame[ ,!names(data_frame) == 'status'])
    ctrl = trainControl(method = "cv", 
                        number = cv, 
                        search = 'grid',classProbs = TRUE, savePredictions = TRUE, 
                        summaryFunction=twoClassSummary )
  }
  
  # run machine learning
  for(i in 1:BS_number){
    
    train_ind = createDataPartition(data_frame$status, p = 0.8, list = FALSE)
    
    if(length(unique(data_frame[-train_ind, "status"])) < 2) {
      message(paste0("Bootstrap iteration ", i, " skipped: test set has only one class."))
      next
    }
    
    train = data_frame[train_ind, ]
    test = data_frame[ -train_ind,!names(data_frame) == "status"]
    
    fit_result <- tryCatch({
      
      if (algorithm == 'lm'){
        model = train(x=train[ , !names(train) == "status"],
                      y= train$status,
                      method = "glmnet", family = "binomial", tuneLength = 5, metric = "ROC",
                      trControl = ctrl,
                      tuneGrid=expand.grid(.alpha=1, .lambda=10^seq(-4, 0, length.out = 20)),
                      penalty.factor = penalty_vec[names(train)[!names(train) == "status"]])
      }
      else if (algorithm == 'enet'){
        model = train(x=train[ , !names(train) == "status"],
                      y= train$status,
                      method = "glmnet", family = "binomial", metric = "ROC",
                      trControl = ctrl,
                      tuneGrid = expand.grid(.alpha = seq(0, 1, length.out = 10),
                                             .lambda = 10^seq(-4, 0, length.out = 20)),
                      penalty.factor = penalty_vec[names(train)[!names(train) == "status"]])
      }
    
    else if (algorithm == 'rf'){
      tunegrid = expand.grid(
        .mtry = c(2, 3, 4, 7, 11, 17, 27, floor(sqrt(n_features)), 41, 64, 99, 154, 237, 367, 567, 876),
        .splitrule = c("extratrees","gini"),
        .min.node.size = c(1,2,3,4,5)
      )
      
      model = train(x=train[ , !names(train) == "status"],
                    y= train$status,
                    tuneGrid = tunegrid, 
                    method = "ranger",  tuneLength = 15, metric = "ROC",
                    num.trees = 500,
                    trControl = ctrl,
                    importance = 'impurity'
      )
    }
    
    else if (algorithm == 'svm lin'){
      model = train(x=train[ , !names(train) == "status"],
                    y= train$status,
                    method = "svmLinear", tuneLength = 5, metric = "ROC",
                    trControl = ctrl,
      )
    }
    
    else if (algorithm == 'svm rad'){
      model = train(x=train[ , !names(train) == "status"],
                    y= train$status,
                    method = "svmRadial", tuneLength = 5, metric = "ROC",
                    trControl = ctrl,
      )
    }
      list(model = model, imp = varImp(model))
    }, error = function(e) {
      message("Bootstrap iteration ", i, " FAILED and was skipped. Reason: ", conditionMessage(e))
      NULL
    })
    if (is.null(fit_result)) next   
    
    models[[i]] = fit_result$model
    importance[[i]] = fit_result$imp
    predictions[[i]] = predict(fit_result$model, newdata = test, type = "prob")
    predictions_raw[[i]] = predict(fit_result$model, newdata = test, type = "raw")
    indices[[i]] = train_ind
    actuals_list[[i]] = data_frame[-train_ind,]$status
    print(paste0('finished loop #',i))
  }
  keep <- !sapply(models, is.null)
  models           <- models[keep]
  importance       <- importance[keep]
  predictions      <- predictions[keep]
  predictions_raw  <- predictions_raw[keep]
  indices          <- indices[keep]
  actuals_list     <- actuals_list[keep]
  
  n_failed <- BS_number - sum(keep)
  if (n_failed > 0) {
    message(n_failed, " of ", BS_number, " bootstrap iterations failed and were excluded from downstream results.")
  }
  return_list = list(models,importance,predictions,predictions_raw,indices,actuals_list)
  names(return_list) = c('models','importance','predictions','predictions_raw','indices','actuals')
  return(return_list)
}

runML_with_lasso <- function(df_ml, cv = 5, bs_count = 500, seed = 123) {
  
  # Scale data
  df_scaled <- scale_manual(df_ml)
  
  # Run the existing ML function 
  ml_results <- runML(df_scaled, 'lm', cv, BS_number = bs_count, seed = seed)
  
  return(ml_results)
}

runML_with_elasticnet <- function(df_ml, cv = 5, bs_count = 500, seed = 123) {
  
  # Scale data
  df_scaled <- scale_manual(df_ml)
  
  # Run the existing ML function 
  ml_results <- runML(df_scaled, 'enet', cv, BS_number = bs_count, seed = seed)
  
  return(ml_results)
}

# =============================
# Make ROC curve
# =============================
calculateROC <- function(ml_results,
                         plot_path = NULL,
                         positive_class = "X1",
                         n_grid = 200) {
  
  if (!requireNamespace("pROC", quietly = TRUE)) {
    stop("Package 'pROC' is required.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }
  
  preds_list <- ml_results$predictions
  actuals_list <- ml_results$actuals
  n_boot <- length(actuals_list)
  fpr_grid <- seq(0, 1, length.out = n_grid)
  
  tpr_mat <- matrix(NA, nrow = length(preds_list), ncol = n_grid)
  aucs <- rep(NA_real_, length(preds_list))
  
  for (i in seq_along(preds_list)) {
    pred_df <- preds_list[[i]]
    actual <- actuals_list[[i]]
    
    if (is.null(pred_df) || is.null(actual)) next
    if (!positive_class %in% colnames(pred_df)) next
    
    probs <- pred_df[[positive_class]]
    actual <- as.factor(actual)
    
    if (length(unique(actual)) < 2) next
    
    roc_obj <- pROC::roc(actual, probs, levels = c("X0", "X1"), direction = "<", quiet = TRUE)
    aucs[i] <- as.numeric(pROC::auc(roc_obj))
    
    fpr <- 1 - roc_obj$specificities
    tpr <- roc_obj$sensitivities
    
    ord <- order(fpr, tpr)
    fpr <- fpr[ord]
    tpr <- tpr[ord]
    
    tmp <- aggregate(tpr, by = list(fpr = fpr), FUN = max)
    fpr_u <- tmp$fpr
    tpr_u <- tmp$x
    
    if (min(fpr_u) > 0) {
      fpr_u <- c(0, fpr_u)
      tpr_u <- c(0, tpr_u)
    }
    if (max(fpr_u) < 1) {
      fpr_u <- c(fpr_u, 1)
      tpr_u <- c(tpr_u, 1)
    }
    
    tpr_mat[i, ] <- approx(fpr_u, tpr_u, xout = fpr_grid, rule = 2)$y
  }
  
  valid_rows <- complete.cases(tpr_mat)
  tpr_mat <- tpr_mat[valid_rows, , drop = FALSE]
  aucs <- aucs[valid_rows]
  
  n_valid <- nrow(tpr_mat)
  if (n_valid == 0) {
    stop("No valid bootstrap iterations with both classes present.")
  }
  
  message("Computed ROC curves for ", n_valid, " valid bootstrap iterations out of ", n_boot)
  
  # Compute statistics
  mean_tpr <- colMeans(tpr_mat)
  lower_tpr <- apply(tpr_mat, 2, quantile, probs = 0.025, na.rm = TRUE)
  upper_tpr <- apply(tpr_mat, 2, quantile, probs = 0.975, na.rm = TRUE)
  mean_auc <- mean(aucs)
  auc_ci <- quantile(aucs, probs = c(0.025, 0.975), na.rm = TRUE)
  sd_auc <- sd(aucs)
  
  # Prepare data for plotting
  ci_ribbon_data <- data.frame(fpr = fpr_grid, ymin = lower_tpr, ymax = upper_tpr)
  mean_data <- data.frame(fpr = fpr_grid, tpr = mean_tpr)
  
  auc_text <- sprintf("AUC = %.3f (95%% CI: %.3f\u2013%.3f)", 
                      mean_auc, auc_ci[1], auc_ci[2])
  
  p <-ggplot() +
    geom_ribbon(data = ci_ribbon_data, aes(x = fpr, ymin = ymin, ymax = ymax, fill = "95% CI"), 
                alpha = 0.6, color = NA) +
    geom_line(data = mean_data, aes(x = fpr, y = tpr, color = "Mean ROC"), 
              linewidth = 1.2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
                color = "darkgrey", linewidth = 0.6) +
    scale_x_continuous(limits = c(0, 1), expand = c(0.01, 0.01)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0.01, 0.01)) +
    scale_fill_manual(values = c("95% CI" = "grey85")) +
    scale_color_manual(values = c("Mean ROC" = "darkblue")) +
    labs(
      x = "1 - Specificity",
      y = "Sensitivity",
      title = "Mean of LASSO Bootstrap ROC curve",
      subtitle = paste0(
        "Mean AUC: ", round(mean_auc, 3),
        " (\u00B1 ", round(sd_auc, 3), "), 95% CI: ",
        round(auc_ci[1], 3), "-",
        round(auc_ci[2], 3)
      )
    ) +
    coord_equal() +
    theme_minimal(base_size = 13) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      #panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90"),
      axis.line = element_line(color = "grey30"),
      axis.ticks = element_line(color = "grey30"),
      legend.position="none"
    ) 
  
  if (!is.null(plot_path)) {
    ggsave(filename = plot_path, plot = p, width = 8, height = 7, dpi = 300)
  }
  
  # Return
  list(
    plot = p,
    roc_data = list(
      fpr_grid = fpr_grid,
      mean_tpr = mean_tpr,
      lower_tpr = lower_tpr,
      upper_tpr = upper_tpr,
      tpr_mat = tpr_mat
    ),
    auc_values = aucs
  )
}

# =============================
# Plot feature importance
# =============================
feature_importance = function(ml_object, 
                              plot_path = "plots/protein_stability_selection.pdf", 
                              scale_importance = FALSE,
                              number = 20,
                              covariate_cols = c("age","sex_male","center_Turkey")) {
  
  if (!isFALSE(plot_path) && !dir.exists(dirname(plot_path))) {
    dir.create(dirname(plot_path), recursive = TRUE)
  }
  
  imp_list     <- ml_object$importance
  num_runs     <- length(imp_list)
  all_proteins <- unique(unlist(lapply(imp_list, function(x) rownames(x$importance))))
  all_proteins <- setdiff(all_proteins, covariate_cols) 
  
  imp_matrix           <- matrix(NA, nrow = length(all_proteins), ncol = num_runs)
  rownames(imp_matrix) <- all_proteins
  
  for (i in seq_len(num_runs)) {
    current_df        <- as.data.frame(imp_list[[i]]$importance)
    vals              <- current_df$Overall
    names(vals)       <- rownames(current_df)
    common_names      <- intersect(names(vals), all_proteins)
    imp_matrix[common_names, i] <- vals[common_names]
  }
  
  freq_prop <- rowSums(imp_matrix > 0, na.rm = TRUE) / num_runs
  
  weight <- data.frame(
    avg_coef = rowMeans(imp_matrix, na.rm = TRUE),
    sd_coef  = apply(imp_matrix, 1, sd, na.rm = TRUE),
    sem_coef = apply(imp_matrix, 1, function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))),
    freq     = freq_prop * 100,
    row.names = all_proteins
  )
  
  # sort by selection frequency, take top N
  weight   <- weight[order(-weight$freq), ]
  number   <- min(number, nrow(weight))
  plot_df  <- weight[seq_len(number), ]
  plot_df$protein <- rownames(plot_df)
  
  p <- ggplot2::ggplot(plot_df,
                       ggplot2::aes(x = reorder(protein, avg_coef),
                                    y = avg_coef,
                                    fill = freq)) +
    ggplot2::geom_bar(stat = "identity", width = 0.7, color = "white") +
    ggplot2::scale_fill_gradient(low = "#deebf7", high = "#084594") +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = avg_coef - sem_coef,
        ymax = avg_coef + sem_coef
      ),
      width = 0.2
    ) + 
    #ggplot2::scale_y_continuous(limits = c(0, 115), breaks = seq(0, 100, 20)) +
    ggplot2::labs(
      x       = "Protein",
      y       = "Importance (+/- Standard Mean Error)",
      title   = paste0("Protein stability across ", num_runs, " bootstraps"),
      fill    = "Selection Frequency"
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::coord_flip() +
    ggplot2::theme(
      axis.title       = ggplot2::element_text(face = "bold"),
      plot.title       = ggplot2::element_text(size = 14, face = "bold"),
      plot.caption     = ggplot2::element_text(size = 8, color = "grey50", hjust = 0),
      legend.position  = "right"
    )
  
  if (!isFALSE(plot_path)) {
    pdf(plot_path, paper = "a4r", width = 11, height = 8)
    print(p)
    dev.off()
    message("Plot saved to ", plot_path)
  }
  
  return(plot_df)
}

# ====================================================
# Check the best protein signature based on nested CV
# ====================================================
find_optimal_signature_CNS_IMMUNE <- function(data_all, 
                                              ranked_proteins, 
                                              fluid = "SERUM",
                                              output_prefix = "optimization_",
                                              inner_folds = 5,
                                              seed = 123,
                                              covariate_cols = c("age","sex_male")) {
  
  set.seed(seed)
  
  min_proteins = 3
  
  X <- as.matrix(data_all[, c(ranked_proteins, covariate_cols)])
  y <- data_all$status
  
  X[, ranked_proteins] <- scale(X[, ranked_proteins])
  y_factor <- as.factor(y)
  
  # performance curve
  auc_list <- vector("list", length(ranked_proteins))
  folds_inner_full <- createFolds(y_factor, k = inner_folds)
  
  for (k in 1:length(ranked_proteins)) {
    
    vars <- c(ranked_proteins[1:k], covariate_cols)
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
        df_tr <- data.frame(y = y_tr, x = X_tr[, ranked_proteins[1]], age = X_tr[,"age"], sex_male = X_tr[,"sex_male"])
        df_val <- data.frame(x = X_val[, ranked_proteins[1]], age = X_val[,"age"], sex_male = X_val[,"sex_male"])
        fit <- glm(y ~ x + age + sex_male, data = df_tr, family = "binomial")
        probs <- predict(fit, newdata = df_val, type = "response")
      } else {
        penalty_vec <- ifelse(colnames(X_tr) %in% covariate_cols, 0, 1)
        fit <- cv.glmnet(X_tr, y_tr, family = "binomial", alpha = 0,
                         nfolds = inner_folds, penalty.factor = penalty_vec)
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
  threshold <- best_auc - 0.4 * (best_sd/sqrt(inner_folds))
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

# =============================
# Get ALS risk score in PGMC
# =============================
run_ALS_signature_workflow_CNS_IMMUNE <- function(
    data_all,
    proteins,
    fluid = "SERUM",
    model_type = c("lasso", "elastic_net"),
    alpha = NULL,
    output_prefix = "output",
    make_plots = TRUE,
    panel_mode = "both"
) {
  
  model_type <- match.arg(model_type)
  
  # define alpha
  if (is.null(alpha)) {
    alpha <- ifelse(model_type == "lasso", 1, 0.5)
  }
  alpha = 0
  
  message("Running model: ", model_type, " (alpha = ", alpha, ") on ", fluid)
  
  covariate_cols <- c("age", "sex_male") 
  
  # prepare data
  data_train <- build_ml_dataset_type(all_data, 
                                      fluid, 
                                      c("ALS",  "CTR"), 
                                      "NPQ", "ALS",
                                      covariate_table = covariate_table,
                                      panel_mode = panel_mode)
  
  panel_zscore_params_discovery <- attr(data_train, "panel_zscore_params") 
  
  # get which proteins used for the premodiALS discovery dataset
  corr_lookup_discovery <- attr(data_train, "corr_lookup")
  
  data_model <- data_train %>%
    select(-SampleName) %>%
    rename(status = type) %>%
    mutate(status = ifelse(status == "ALS", 1, 0))
  
  X <- as.matrix(data_model[, c(proteins,"age","sex_male")])
  y <- data_model$status
  
  # scale
  X_scaled <- scale(X[,proteins])
  scaling_center <- attr(X_scaled, "scaled:center")
  scaling_scale  <- attr(X_scaled, "scaled:scale")
  
  X_scaled <- cbind(X_scaled, X[, covariate_cols, drop = FALSE])
  penalty_vec <- ifelse(colnames(X_scaled) %in% covariate_cols, 0, 1)
  
  # model with protein signature
  set.seed(123)
  cv_model <- cv.glmnet(X_scaled, y, family = "binomial", alpha = alpha,
                        penalty.factor = penalty_vec)
  
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
                                     "NPQ", "PGMC",
                                     panel_zscore_params = panel_zscore_params_discovery,
                                     covariate_table = covariate_table,
                                     panel_mode = panel_mode) 
  
  X_pgmc <- as.matrix(data_pgmc[, c(proteins, covariate_cols)])
  X_pgmc_scaled <- X_pgmc
  X_pgmc_scaled[, proteins] <- scale(X_pgmc[, proteins],
                                     center = scaling_center,
                                     scale  = scaling_scale)
  
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
    scaling = list(center = scaling_center, scale = scaling_scale),
    corr_lookup = corr_lookup_discovery,
    panel_zscore_params = panel_zscore_params_discovery 
  ))
}

# =============================
# Heatmap of protein signature
# =============================
run_heatmap_signature_CNS_IMMUNE <- function(
    data_all,
    fluid,
    proteins,
    risk_scores_df,   
    group_colors,
    output_prefix = "heatmap",
    highlight_ids = NULL,
    panel_mode = "both"
) {
  
  message("Running heatmap for ", fluid)
  
  # heatmap data
  data_heatmap <- build_ml_dataset_type(data_all,
                                        fluid, c("PGMC","ALS","CTR"),
                                        "NPQ","ALS",
                                        covariate_table = covariate_table,
                                        panel_mode = panel_mode) %>%
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


# ------------------------------------------------------------------------------
# ====================================================
## Compare CNS and immune panels
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

## External dataset: MAXOMOD NULISA
MAXOMOD_NULISA_CNS = read_excel("~/Documents/HMGU/premodiALS/MAXOMOD NULISA data/P005_BSHRI_NULISAseq_CNSDiseasePanel_NPQCounts_2025_03_10.xlsx")
MAXOMOD_NULISA_IMMUNE = read_excel("~/Documents/HMGU/premodiALS/MAXOMOD NULISA data/P005_BSHRI_NULISAseq_InflammationPanel_NPQ_03022026.xlsx")
MAXOMOD_IDs = read_excel("~/Documents/HMGU/premodiALS/MAXOMOD NULISA data/all_participants_IDs.xlsx")

external_cns <- MAXOMOD_NULISA_CNS %>%
  filter(SampleMatrixType %in% c("CSF","PLASMA","SERUM")) %>%
  select(SampleName, SampleMatrixType, Target, UniProtID, ProteinName, NPQ) %>%
  mutate(panel = "CNS")

external_immune <- MAXOMOD_NULISA_IMMUNE %>%
  filter(SampleMatrixType %in% c("CSF","PLASMA","SERUM")) %>%
  select(SampleName, SampleMatrixType, Target, UniProtID, ProteinName, NPQ) %>%
  mutate(panel = "IMMUNE")

external_all_data <- bind_rows(external_cns, external_immune) %>%
  mutate(
    panel            = factor(panel, levels = c("CNS", "IMMUNE")),
    SampleMatrixType = factor(SampleMatrixType, levels = c("SERUM", "PLASMA", "CSF"))
  )

external_covariate_table <- build_external_covariate_table(MAXOMOD_IDs)


## =========================================================================
## PART 1: HEATMAPS -- one per fluid, rows split by panel (CNS vs Immune)
## =========================================================================
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
corr_mat <- cor(mat, use = "pairwise.complete.obs", method = "pearson")  
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
## Main loop to run all the functions
for (fluid in fluids) {
  
  meta_data = all_data %>% 
    select(SampleName,PatientID,ParticipantCode) %>%
    left_join(Sex_age_all_participants %>% dplyr::rename(PatientID = Pseudonyme)) %>%
    mutate(
      center = dplyr::case_when(
        grepl("TR", ParticipantCode) ~ "Turkey",
        grepl("CH", ParticipantCode) ~ "Switzerland",
        grepl("DE", ParticipantCode) ~ "Germany",
        grepl("SK", ParticipantCode) ~ "Slovakia",
        grepl("FR", ParticipantCode) ~ "France",
        grepl("IL", ParticipantCode) ~ "Israel",
        TRUE                 ~ NA_character_
      ))
  
  covariate_table = build_covariate_table(meta_data)
  
  ## Prepare merged, z-scored, correlation-averaged datasets
  protein_data_PGMCvsCTR_new <- build_ml_dataset(all_data, fluid, c("PGMC", "CTR"), "NPQ",     "PGMC",
                                                 covariate_table)
  protein_data_ALSvsCTR_new  <- build_ml_dataset(all_data, fluid, c("ALS",  "CTR"), "NPQ", "ALS",
                                                 covariate_table)
  protein_data_ALSvsPGMC_new <- build_ml_dataset(all_data, fluid, c("ALS",  "PGMC"), "NPQ", "ALS",
                                                 covariate_table)
  
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
# Protein signature
## -> Lasso + Serum
proteins_serum_ALS_CTR = c("NEFL_CNS","pTau-181_CNS","BASP1_CNS","CLEC4A_IMMUNE",
                           "ACHE_CNS","CXCL14_IMMUNE","HAVCR1_IMMUNE","SELE_IMMUNE",
                           "TNFRSF13C_IMMUNE","GDNF_CNS","NGF_CNS")
lasso_optimal_protein_signature_serum_both = find_optimal_signature_CNS_IMMUNE(data_all = build_ml_dataset(all_data, 
                                                                                                      "SERUM", 
                                                                                                      c("ALS",  "CTR"), 
                                                                                                      "NPQ", "ALS",
                                                                                                      covariate_table),
                                                               ranked_proteins = proteins_serum_ALS_CTR,
                                                               fluid = "SERUM",
                                                               output_prefix = "plots/CNS_IMMUNE_panels/ML/Lasso/performance_")

## -> Lasso + Plasma
proteins_plasma_ALS_CTR = c("NEFL_CNS","pTau-181_CNS","TEK_CNS","TNFSF9_IMMUNE",
                            "IL1B_CNS","CD3E_IMMUNE","GDNF_CNS","CXCL14_IMMUNE",
                            "VSNL1_CNS","TNFRSF13C_IMMUNE","TSLP_IMMUNE")
lasso_optimal_protein_signature_plasma_both = find_optimal_signature_CNS_IMMUNE(data_all = build_ml_dataset(all_data, 
                                                                                                      "PLASMA", 
                                                                                                      c("ALS",  "CTR"), 
                                                                                                      "NPQ", "ALS",
                                                                                                      covariate_table),
                                                                          ranked_proteins = proteins_plasma_ALS_CTR,
                                                                          fluid = "PLASMA",
                                                                          output_prefix = "plots/CNS_IMMUNE_panels/ML/Lasso/performance_")

# Elastic Net + Serum
proteins_serum_ALS_CTR = c("NEFL_CNS","pTau-181_CNS","CLEC4A_IMMUNE","GDNF_CNS",
                           "HAVCR1_IMMUNE","BASP1_CNS","TNFRSF13C_IMMUNE","SELE_IMMUNE",
                           "pTau-231_CNS","CXCL14_IMMUNE","ACHE_CNS")
enet_optimal_protein_signature_serum_both = find_optimal_signature_CNS_IMMUNE(data_all = build_ml_dataset(all_data, 
                                                                                                     "SERUM", 
                                                                                                     c("ALS",  "CTR"), 
                                                                                                     "NPQ", "ALS",
                                                                                                     covariate_table),
                                                              ranked_proteins = proteins_serum_ALS_CTR,
                                                              fluid = "SERUM",
                                                              output_prefix = "plots/CNS_IMMUNE_panels/ML/Elastic Net/performance_")

# Elastic Net + Plasma
proteins_plasma_ALS_CTR = c("NEFL_CNS","pTau-181_CNS","TEK_CNS","IL1B_CNS","TNFRSF9_IMMUNE",
                            "CD3E_IMMUNE","GDNF_CNS","CXCL14_IMMUNE","VSNL1_CNS",
                            "CCL25_IMMUNE","IL1RL1_IMMUNE")
enet_optimal_protein_signature_plasma_both = find_optimal_signature_CNS_IMMUNE(data_all = build_ml_dataset(all_data, 
                                                                                                      "PLASMA", 
                                                                                                      c("ALS",  "CTR"), 
                                                                                                      "NPQ", "ALS",
                                                                                                      covariate_table),
                                                               ranked_proteins = proteins_plasma_ALS_CTR,
                                                               fluid = "PLASMA",
                                                               output_prefix = "plots/CNS_IMMUNE_panels/ML/Elastic Net/performance_")

# Elastic Net + CSF
proteins_CSF_ALS_CTR =  c("NEFL_CNS","NEFH_CNS","CCL13_CNS","TNFRSF13C_IMMUNE",
                          "IL18_IMMUNE","IFNG_IMMUNE","PDCD1LG2_IMMUNE",
                          "IL4_CNS","EPO_IMMUNE","IL27_IMMUNE","CSF3R_IMMUNE",
                          "IL33_IMMUNE")

enet_optimal_protein_signature_CSF_both = find_optimal_signature_CNS_IMMUNE(data_all = build_ml_dataset(all_data, 
                                                                                                   "CSF", 
                                                                                                   c("ALS",  "CTR"), 
                                                                                                   "NPQ", "ALS",
                                                                                                   covariate_table),
                                                            ranked_proteins = proteins_CSF_ALS_CTR,
                                                            fluid = "CSF",
                                                            output_prefix = "plots/CNS_IMMUNE_panels/ML/Elastic Net/performance_")

# ====================================================
# Compute an ALS risk score based on ALS vs CTR model 
## -> Lasso + Serum
proteins_serum_ALS_CTR_both = lasso_optimal_protein_signature_serum_both$optimal_proteins

lasso_PGMC_serum_signature_both = run_ALS_signature_workflow_CNS_IMMUNE(all_data,
                                                        proteins_serum_ALS_CTR_both,
                                                        fluid = "SERUM",
                                                        model_type = "lasso",
                                                        output_prefix = "plots/CNS_IMMUNE_panels/ML/Lasso/")

external_lasso_serum_wide_both <- build_external_ml_dataset(
  external_all_data, "SERUM", "NPQ",
  external_covariate_table,
  corr_lookup_discovery = lasso_PGMC_serum_signature_both$corr_lookup,
  panel_zscore_params = lasso_PGMC_serum_signature_both$panel_zscore_params
)

external_validation_lasso_serum_both <- validate_external(
  cv_model        = lasso_PGMC_serum_signature_both$model,
  external_wide   = external_lasso_serum_wide_both,
  proteins        = lasso_optimal_protein_signature_serum_both$optimal_proteins,
  covariate_cols  = c("age", "sex_male"),
  scaling_center  = lasso_PGMC_serum_signature_both$scaling$center,
  scaling_scale   = lasso_PGMC_serum_signature_both$scaling$scale,
  fluid_label     = "SERUM",
  plot_path       = "plots/CNS_IMMUNE_panels/ML/Lasso/ROC_external_validation_SERUM.pdf"
)

participant_code_label = lasso_PGMC_serum_signature_both$results %>%
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
                      lasso_optimal_protein_signature_serum_both$optimal_proteins,
                      lasso_PGMC_serum_signature_both$results,
                      group_colors = group_colors,
                      output_prefix = "plots/CNS_IMMUNE_panels/ML/Lasso/heatmap_SERUM",
                      highlight_ids = participant_code_label)

## -> Lasso + plasma
proteins_plasma_ALS_CTR_both = lasso_optimal_protein_signature_plasma_both$optimal_proteins

lasso_PGMC_plasma_signature_both = run_ALS_signature_workflow_CNS_IMMUNE(all_data,
                                                                   proteins_plasma_ALS_CTR_both,
                                                                   fluid = "PLASMA",
                                                                   model_type = "lasso",
                                                                   output_prefix = "plots/CNS_IMMUNE_panels/ML/Lasso/")

external_lasso_plasma_wide_both <- build_external_ml_dataset(
  external_all_data, "PLASMA", "NPQ",
  external_covariate_table,
  corr_lookup_discovery = lasso_PGMC_plasma_signature_both$corr_lookup,
  panel_zscore_params = lasso_PGMC_plasma_signature_both$panel_zscore_params
)

external_validation_lasso_plasma_both <- validate_external(
  cv_model        = lasso_PGMC_plasma_signature_both$model,
  external_wide   = external_lasso_plasma_wide_both,
  proteins        = lasso_optimal_protein_signature_plasma_both$optimal_proteins,
  covariate_cols  = c("age", "sex_male"),
  scaling_center  = lasso_PGMC_plasma_signature_both$scaling$center,
  scaling_scale   = lasso_PGMC_plasma_signature_both$scaling$scale,
  fluid_label     = "PLASMA",
  plot_path       = "plots/CNS_IMMUNE_panels/ML/Lasso/ROC_external_validation_plasma.pdf"
)

participant_code_label = lasso_PGMC_plasma_signature_both$results %>%
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
                                 "PLASMA",
                                 lasso_optimal_protein_signature_plasma_both$optimal_proteins,
                                 lasso_PGMC_plasma_signature_both$results,
                                 group_colors = group_colors,
                                 output_prefix = "plots/CNS_IMMUNE_panels/ML/Lasso/heatmap_plasma",
                                 highlight_ids = participant_code_label)


##### ----
# Elastic Net

## -> Elastic Net + Serum
EN_PGMC_serum_signature_both = run_ALS_signature_workflow_CNS_IMMUNE(all_data,
                                                     enet_optimal_protein_signature_serum_both$optimal_proteins,
                                                     fluid = "SERUM",
                                                     model_type = "elastic_net",
                                                     output_prefix = "plots/CNS_IMMUNE_panels/ML/Elastic Net/")

external_enet_serum_wide_both <- build_external_ml_dataset(
  external_all_data, "SERUM", "NPQ",
  external_covariate_table,
  corr_lookup_discovery = EN_PGMC_serum_signature_both$corr_lookup,
  panel_zscore_params = EN_PGMC_serum_signature_both$panel_zscore_params
)

external_validation_enet_serum_both <- validate_external(
  cv_model        = EN_PGMC_serum_signature_both$model,
  external_wide   = external_enet_serum_wide_both,
  proteins        = enet_optimal_protein_signature_serum_both$optimal_proteins,
  covariate_cols  = c("age", "sex_male"),
  scaling_center  = EN_PGMC_serum_signature_both$scaling$center,
  scaling_scale   = EN_PGMC_serum_signature_both$scaling$scale,
  fluid_label     = "SERUM",
  plot_path       = "plots/CNS_IMMUNE_panels/ML/Elastic Net/ROC_external_validation_SERUM.pdf"
)


participant_code_label = EN_PGMC_serum_signature_both_both$results %>%
  left_join(protein_data_IDs %>% select(SampleName,ParticipantCode,type)) %>%
  distinct() %>%
  filter(type == "PGMC" & ALS_risk_score > 0) %>%
  pull(ParticipantCode)

# ================================================================================
# Unsupervised visualisation (heatmap) of PGMC, ALS, CTR based on 9-protein signature

run_heatmap_signature_CNS_IMMUNE(all_data,
                      "SERUM",
                      enet_optimal_protein_signature_serum_both$optimal_proteins,
                      EN_PGMC_serum_signature_both$results,
                      group_colors = group_colors,
                      output_prefix = "plots/CNS_IMMUNE_panels/ML/Elastic Net/heatmap_SERUM",
                      highlight_ids = participant_code_label)

## -> Elastic Net + Plasma
EN_PGMC_plasma_signature_both = run_ALS_signature_workflow_CNS_IMMUNE(all_data,
                                                      enet_optimal_protein_signature_plasma_both$optimal_proteins,
                                                      fluid = "PLASMA",
                                                      model_type = "elastic_net",
                                                      output_prefix = "plots/CNS_IMMUNE_panels/ML/Elastic Net/")

external_enet_plasma_wide_both <- build_external_ml_dataset(
  external_all_data, "PLASMA", "NPQ",
  external_covariate_table,
  corr_lookup_discovery = EN_PGMC_plasma_signature_both$corr_lookup,
  panel_zscore_params = EN_PGMC_plasma_signature_both$panel_zscore_params
)

external_validation_enet_plasma_both <- validate_external(
  cv_model        = EN_PGMC_plasma_signature_both$model,
  external_wide   = external_enet_plasma_wide_both,
  proteins        = enet_optimal_protein_signature_plasma_both$optimal_proteins,
  covariate_cols  = c("age", "sex_male"),
  scaling_center  = EN_PGMC_plasma_signature_both$scaling$center,
  scaling_scale   = EN_PGMC_plasma_signature_both$scaling$scale,
  fluid_label     = "PLASMA",
  plot_path       = "plots/CNS_IMMUNE_panels/ML/Elastic Net/ROC_external_validation_plasma.pdf"
)


participant_code_label = EN_PGMC_plasma_signature_both$results %>%
  left_join(protein_data_IDs %>% select(SampleName,ParticipantCode,type)) %>%
  distinct() %>%
  filter(type == "PGMC" & ALS_risk_score > 0) %>%
  pull(ParticipantCode)

# ================================================================================
# Unsupervised visualisation (heatmap) of PGMC, ALS, CTR based on 3-protein signature

run_heatmap_signature_CNS_IMMUNE(all_data,
                      "PLASMA",
                      enet_optimal_protein_signature_plasma_both$optimal_proteins,
                      EN_PGMC_plasma_signature_both$results,
                      group_colors = group_colors,
                      output_prefix = "plots/CNS_IMMUNE_panels/ML/Elastic Net/heatmap_PLASMA",
                      highlight_ids = participant_code_label)

## -> Elastic Net + CSF
EN_PGMC_CSF_signature_both = run_ALS_signature_workflow_CNS_IMMUNE(all_data,
                                                   enet_optimal_protein_signature_CSF_both$optimal_proteins,
                                                   fluid = "CSF",
                                                   model_type = "elastic_net",
                                                   output_prefix = "plots/CNS_IMMUNE_panels/ML/Elastic Net/")

external_enet_csf_wide_both <- build_external_ml_dataset(
  external_all_data, "CSF", "NPQ",
  external_covariate_table,
  corr_lookup_discovery = EN_PGMC_CSF_signature_both$corr_lookup,
  panel_zscore_params = EN_PGMC_CSF_signature_both$panel_zscore_params
)

external_validation_enet_csf_both <- validate_external(
  cv_model        = EN_PGMC_CSF_signature_both$model,
  external_wide   = external_enet_csf_wide_both,
  proteins        = enet_optimal_protein_signature_CSF_both$optimal_proteins,
  covariate_cols  = c("age", "sex_male"),
  scaling_center  = EN_PGMC_CSF_signature_both$scaling$center,
  scaling_scale   = EN_PGMC_CSF_signature_both$scaling$scale,
  fluid_label     = "CSF",
  plot_path       = "plots/CNS_IMMUNE_panels/ML/Elastic Net/ROC_external_validation_CSF.pdf"
)

participant_code_label = EN_PGMC_CSF_signature_both$results %>%
  left_join(protein_data_IDs %>% select(SampleName,ParticipantCode,type)) %>%
  distinct() %>%
  filter(type == "PGMC" & ALS_risk_score > 0) %>%
  pull(ParticipantCode)

# ================================================================================
# Unsupervised visualisation (heatmap) of PGMC, ALS, CTR based on 9-protein signature

run_heatmap_signature_CNS_IMMUNE(all_data,
                      "CSF",
                      enet_optimal_protein_signature_CSF_both$optimal_proteins,
                      EN_PGMC_CSF_signature_both$results,
                      group_colors = group_colors,
                      output_prefix = "plots/CNS_IMMUNE_panels/ML/Elastic Net/heatmap_CSF",
                      highlight_ids = participant_code_label)


## ==========================================================
## PART 5: ML -- CNS panel only per fluid
## ==========================================================

for (fluid in fluids) {
  
  meta_data = all_data %>% 
    select(SampleName,PatientID,ParticipantCode) %>%
    left_join(Sex_age_all_participants %>% dplyr::rename(PatientID = Pseudonyme)) %>%
    mutate(
      center = dplyr::case_when(
        grepl("TR", ParticipantCode) ~ "Turkey",
        grepl("CH", ParticipantCode) ~ "Switzerland",
        grepl("DE", ParticipantCode) ~ "Germany",
        grepl("SK", ParticipantCode) ~ "Slovakia",
        grepl("FR", ParticipantCode) ~ "France",
        grepl("IL", ParticipantCode) ~ "Israel",
        TRUE                 ~ NA_character_
      ))
  
  covariate_table = build_covariate_table(meta_data)
  
  ## Prepare merged, z-scored, correlation-averaged datasets
  protein_data_PGMCvsCTR_CNS <- build_ml_dataset(all_data, fluid, c("PGMC", "CTR"), "NPQ",     "PGMC",
                                                 covariate_table,
                                                 panel_mode = "CNS")
  protein_data_ALSvsCTR_CNS  <- build_ml_dataset(all_data, fluid, c("ALS",  "CTR"), "NPQ", "ALS",
                                                 covariate_table,
                                                 panel_mode = "CNS")
  protein_data_ALSvsPGMC_CNS <- build_ml_dataset(all_data, fluid, c("ALS",  "PGMC"), "NPQ", "ALS",
                                                 covariate_table,
                                                 panel_mode = "CNS")
  
  ## Run Lasso and ROC curve with 5-fold cv and 500 bootstrap iterations
  lm_PGMC_CTR_CNS <- runML_with_lasso(protein_data_PGMCvsCTR_CNS, bs_count = 500)
  final_roc_plot_PGMC_CTR_CNS <- calculateROC(lm_PGMC_CTR_CNS,
                                          paste0("plots/ML/Lasso/ROC_lm_CNS_", fluid, "_PGMC_CTR.pdf"))
  
  lm_ALS_CTR_CNS <- runML_with_lasso(protein_data_ALSvsCTR_CNS, bs_count = 500)
  final_roc_plot_ALS_CTR_CNS <- calculateROC(lm_ALS_CTR_CNS,
                                         paste0("plots/ML/Lasso/ROC_lm_CNS_", fluid, "_ALS_CTR.pdf"))
  
  lm_ALS_PGMC_CNS <- runML_with_lasso(protein_data_ALSvsPGMC_CNS, bs_count = 500)
  final_roc_plot_ALS_PGMC_CNS <- calculateROC(lm_ALS_PGMC_CNS,
                                          paste0("plots/ML/Lasso/ROC_lm_CNS_", fluid, "_ALS_PGMC.pdf"))
  
  ## Extract weights + plot averaged
  lm_weights_PGMC_CTR_CNS <- feature_importance(lm_PGMC_CTR_CNS,
                                            plot_path = paste0("plots/ML/Lasso/protein_selection_CNS_", fluid, "_PGMC_CTR.pdf"))
  write_xlsx(lm_weights_PGMC_CTR_CNS, path = paste0("results/weights_lm_CNS_", fluid, "_PGMC_CTR.xlsx"))
  
  lm_weights_ALS_CTR_CNS <- feature_importance(lm_ALS_CTR_CNS,
                                           plot_path = paste0("plots/ML/Lasso/protein_selection_CNS_", fluid, "_ALS_CTR.pdf"))
  write_xlsx(lm_weights_ALS_CTR_CNS, path = paste0("results/weights_lm_CNS_", fluid, "_ALS_CTR.xlsx"))
  
  lm_weights_ALS_PGMC_CNS <- feature_importance(lm_ALS_PGMC_CNS,
                                            plot_path = paste0("plots/ML/Lasso/protein_selection_CNS_", fluid, "_ALS_PGMC.pdf"))
  write_xlsx(lm_weights_ALS_PGMC_CNS, path = paste0("results/weights_lm_CNS_", fluid, "_ALS_PGMC.xlsx"))
  
  # Run Elastic Net and ROC curve with 5-fold cv and 500 bootstrap iterations
  lm_enet_PGMC_CTR_CNS <- runML_with_elasticnet(protein_data_PGMCvsCTR_CNS,bs_count = 500)
  final_roc_plot_PGMC_CTR_enet_CNS = calculateROC(lm_enet_PGMC_CTR_CNS,
                                              paste0("plots/ML/Elastic Net/ROC_lm_CNS_",fluid,"_PGMC_CTR.pdf"))
  lm_enet_ALS_CTR_CNS <- runML_with_elasticnet(protein_data_ALSvsCTR_CNS,bs_count = 500)
  final_roc_plot_ALS_CTR_enet_CNS = calculateROC(lm_enet_ALS_CTR_CNS,
                                             paste0("plots/ML/Elastic Net/ROC_lm_CNS_",fluid,"_ALS_CTR.pdf"))
  lm_enet_ALS_PGMC_CNS <- runML_with_elasticnet(protein_data_ALSvsPGMC_CNS,bs_count = 500)
  final_roc_plot_ALS_PGMC_enet_CNS = calculateROC(lm_enet_ALS_PGMC_CNS,
                                              paste0("plots/ML/Elastic Net/ROC_lm_CNS_",fluid,"_ALS_PGMC.pdf"))
  
  # extract weights + plot averaged
  lm_weights_PGMC_CTR_enet_CNS = feature_importance(lm_enet_PGMC_CTR_CNS,
                                                plot_path = paste0("plots/ML/Elastic Net/protein_selection_CNS_",fluid,"_PGMC_CTR.pdf"))
  write_xlsx(lm_weights_PGMC_CTR_enet_CNS, path = paste0("results/weights_lm_enet_CNS_",fluid,"_PGMC_CTR.xlsx"))
  lm_weights_ALS_CTR_enet_CNS = feature_importance(lm_enet_ALS_CTR_CNS,
                                               plot_path = paste0("plots/ML/Elastic Net/protein_selection_CNS_",fluid,"_ALS_CTR.pdf"))
  write_xlsx(lm_weights_ALS_CTR_enet_CNS, path = paste0("results/weights_lm_enet_CNS_",fluid,"_ALS_CTR.xlsx"))
  lm_weights_ALS_PGMC_enet_CNS = feature_importance(lm_enet_ALS_PGMC_CNS,
                                                plot_path = paste0("plots/ML/Elastic Net/protein_selection_CNS_",fluid,"_ALS_PGMC.pdf"))
  write_xlsx(lm_weights_ALS_PGMC_enet_CNS, path = paste0("results/weights_lm_enet_CNS_",fluid,"_ALS_PGMC.xlsx"))
}

# ================================================================================
# Protein signature (CNS)
proteins_serum_ALS_CTR = c("NEFL_CNS","pTau-181_CNS","ACHE_CNS","BASP1_CNS",
                           "GDNF_CNS","NGF_CNS","pTau-231_CNS","IL33_CNS",
                           "CXCL8_CNS","NRGN_CNS","PDLIM5_CNS")
lasso_optimal_protein_signature_serum_CNS = find_optimal_signature_CNS_IMMUNE(data_all = build_ml_dataset(all_data, 
                                                                                                      "SERUM", 
                                                                                                      c("ALS",  "CTR"), 
                                                                                                      "NPQ", "ALS",
                                                                                                      covariate_table,
                                                                                                      panel_mode = "CNS"),
                                                                          ranked_proteins = proteins_serum_ALS_CTR,
                                                                          fluid = "SERUM",
                                                                          output_prefix = "plots/ML/Lasso/performance_CNS_")

## -> Lasso + Plasma
proteins_plasma_ALS_CTR = c("NEFL_CNS","pTau-181_CNS","TEK_CNS",
                            "IL1B_CNS","CXCL8_CNS","GDNF_CNS","PDLIM5_CNS",
                            "CX3CL1_CNS","VSNL1_CNS","FCN2_CNS","HBA1_CNS")
lasso_optimal_protein_signature_plasma_CNS = find_optimal_signature_CNS_IMMUNE(data_all = build_ml_dataset(all_data, 
                                                                                                       "PLASMA", 
                                                                                                       c("ALS",  "CTR"), 
                                                                                                       "NPQ", "ALS",
                                                                                                       covariate_table,
                                                                                                       panel_mode = "CNS"),
                                                                           ranked_proteins = proteins_plasma_ALS_CTR,
                                                                           fluid = "PLASMA",
                                                                           output_prefix = "plots/ML/Lasso/performance_CNS_")

# Elastic Net + Serum
proteins_serum_ALS_CTR = c("NEFL_CNS","pTau-181_CNS","GDNF_CNS","ACHE_CNS",
                           "BASP1_CNS","pTau-231_CNS","NGF_CNS","IL33_CNS",
                           "CXCL8_CNS","NRGN_CNS","PDLIM5_CNS")
enet_optimal_protein_signature_serum_CNS = find_optimal_signature_CNS_IMMUNE(data_all = build_ml_dataset(all_data, 
                                                                                                     "SERUM", 
                                                                                                     c("ALS",  "CTR"), 
                                                                                                     "NPQ", "ALS",
                                                                                                     covariate_table,
                                                                                                     panel_mode = "CNS"),
                                                                         ranked_proteins = proteins_serum_ALS_CTR,
                                                                         fluid = "SERUM",
                                                                         output_prefix = "plots/ML/Elastic Net/performance_CNS_")

# Elastic Net + Plasma
proteins_plasma_ALS_CTR = c("NEFL_CNS","pTau-181_CNS","TEK_CNS","IL1B_CNS",
                            "CXCL8_CNS","GDNF_CNS","VSNL1_CNS","PDLIM5_CNS",
                            "CX3CL1_CNS","FCN2_CNS","SFRP1_CNS")
enet_optimal_protein_signature_plasma_CNS = find_optimal_signature_CNS_IMMUNE(data_all = build_ml_dataset(all_data, 
                                                                                                      "PLASMA", 
                                                                                                      c("ALS",  "CTR"), 
                                                                                                      "NPQ", "ALS",
                                                                                                      covariate_table,
                                                                                                      panel_mode = "CNS"),
                                                                          ranked_proteins = proteins_plasma_ALS_CTR,
                                                                          fluid = "PLASMA",
                                                                          output_prefix = "plots/ML/Elastic Net/performance_CNS_")

# Elastic Net + CSF
proteins_CSF_ALS_CTR =  c("NEFH_CNS","NEFL_CNS","CCL13_CNS","IL4_CNS",
                          "UCHL1_CNS","IL6R_CNS","pTDP43-409_CNS","CCL2_CNS",
                          "CALB2_CNS","Aβ38_CNS","ICAM1_CNS")

enet_optimal_protein_signature_CSF_CNS = find_optimal_signature_CNS_IMMUNE(data_all = build_ml_dataset(all_data, 
                                                                                                   "CSF", 
                                                                                                   c("ALS",  "CTR"), 
                                                                                                   "NPQ", "ALS",
                                                                                                   covariate_table,
                                                                                                   panel_mode = "CNS"),
                                                                       ranked_proteins = proteins_CSF_ALS_CTR,
                                                                       fluid = "CSF",
                                                                       output_prefix = "plots/ML/Elastic Net/performance_CNS_")

# ====================================================
# Compute an ALS risk score based on ALS vs CTR model 
## -> Lasso + Serum
proteins_serum_ALS_CTR_CNS = lasso_optimal_protein_signature_serum_CNS$optimal_proteins

lasso_PGMC_serum_signature_CNS = run_ALS_signature_workflow_CNS_IMMUNE(all_data,
                                                                   proteins_serum_ALS_CTR_CNS,
                                                                   fluid = "SERUM",
                                                                   model_type = "lasso",
                                                                   output_prefix = "plots/ML/Lasso/CNS_",
                                                                   panel_mode = "CNS")

external_lasso_serum_wide_CNS <- build_external_ml_dataset(
  external_all_data, "SERUM", "NPQ",
  external_covariate_table,
  corr_lookup_discovery = lasso_PGMC_serum_signature_CNS$corr_lookup,
  panel_zscore_params = lasso_PGMC_serum_signature_CNS$panel_zscore_params,
  panel_mode = "CNS"
)

external_validation_lasso_serum_CNS <- validate_external(
  cv_model        = lasso_PGMC_serum_signature_CNS$model,
  external_wide   = external_lasso_serum_wide_CNS,
  proteins        = lasso_optimal_protein_signature_serum_CNS$optimal_proteins,
  covariate_cols  = c("age", "sex_male"),
  scaling_center  = lasso_PGMC_serum_signature_CNS$scaling$center,
  scaling_scale   = lasso_PGMC_serum_signature_CNS$scaling$scale,
  fluid_label     = "SERUM",
  plot_path       = "plots/ML/Lasso/ROC_external_validation_SERUM_CNS.pdf"
)

participant_code_label = lasso_PGMC_serum_signature_CNS$results %>%
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
                                 lasso_optimal_protein_signature_serum_CNS$optimal_proteins,
                                 lasso_PGMC_serum_signature_CNS$results,
                                 group_colors = group_colors,
                                 output_prefix = "plots/ML/Lasso/heatmap_SERUM_CNS",
                                 highlight_ids = participant_code_label,
                                 panel_mode = "CNS")

## -> Lasso + plasma
proteins_plasma_ALS_CTR_CNS = lasso_optimal_protein_signature_plasma_CNS$optimal_proteins

lasso_PGMC_plasma_signature_CNS = run_ALS_signature_workflow_CNS_IMMUNE(all_data,
                                                                    proteins_plasma_ALS_CTR_CNS,
                                                                    fluid = "PLASMA",
                                                                    model_type = "lasso",
                                                                    output_prefix = "plots/ML/Lasso/CNS_",
                                                                    panel_mode = "CNS")

external_lasso_plasma_wide_CNS <- build_external_ml_dataset(
  external_all_data, "PLASMA", "NPQ",
  external_covariate_table,
  corr_lookup_discovery = lasso_PGMC_plasma_signature_CNS$corr_lookup,
  panel_zscore_params = lasso_PGMC_plasma_signature_CNS$panel_zscore_params,
  panel_mode = "CNS"
)

external_validation_lasso_plasma_CNS <- validate_external(
  cv_model        = lasso_PGMC_plasma_signature_CNS$model,
  external_wide   = external_lasso_plasma_wide_CNS,
  proteins        = lasso_optimal_protein_signature_plasma_CNS$optimal_proteins,
  covariate_cols  = c("age", "sex_male"),
  scaling_center  = lasso_PGMC_plasma_signature_CNS$scaling$center,
  scaling_scale   = lasso_PGMC_plasma_signature_CNS$scaling$scale,
  fluid_label     = "PLASMA",
  plot_path       = "plots/ML/Lasso/ROC_external_validation_plasma_CNS.pdf"
)

participant_code_label = lasso_PGMC_plasma_signature_CNS$results %>%
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
                                 "PLASMA",
                                 lasso_optimal_protein_signature_plasma_CNS$optimal_proteins,
                                 lasso_PGMC_plasma_signature_CNS$results,
                                 group_colors = group_colors,
                                 output_prefix = "plots/ML/Lasso/heatmap_plasma_CNS",
                                 highlight_ids = participant_code_label,
                                 panel_mode = "CNS"
                                 )


##### ----
# Elastic Net

## -> Elastic Net + Serum
EN_PGMC_serum_signature_CNS = run_ALS_signature_workflow_CNS_IMMUNE(all_data,
                                                                enet_optimal_protein_signature_serum_CNS$optimal_proteins,
                                                                fluid = "SERUM",
                                                                model_type = "elastic_net",
                                                                output_prefix = "plots/ML/Elastic Net/CNS_",
                                                                panel_mode = "CNS")

external_enet_serum_wide_CNS <- build_external_ml_dataset(
  external_all_data, "SERUM", "NPQ",
  external_covariate_table,
  corr_lookup_discovery = EN_PGMC_serum_signature_CNS$corr_lookup,
  panel_zscore_params = EN_PGMC_serum_signature_CNS$panel_zscore_params,
  panel_mode = "CNS"
)

external_validation_enet_serum_CNS <- validate_external(
  cv_model        = EN_PGMC_serum_signature_CNS$model,
  external_wide   = external_enet_serum_wide_CNS,
  proteins        = enet_optimal_protein_signature_serum_CNS$optimal_proteins,
  covariate_cols  = c("age", "sex_male"),
  scaling_center  = EN_PGMC_serum_signature_CNS$scaling$center,
  scaling_scale   = EN_PGMC_serum_signature_CNS$scaling$scale,
  fluid_label     = "SERUM",
  plot_path       = "plots/ML/Elastic Net/ROC_external_validation_SERUM_CNS.pdf"
)


participant_code_label = EN_PGMC_serum_signature_CNS$results %>%
  left_join(protein_data_IDs %>% select(SampleName,ParticipantCode,type)) %>%
  distinct() %>%
  filter(type == "PGMC" & ALS_risk_score > 0) %>%
  pull(ParticipantCode)

# ================================================================================
# Unsupervised visualisation (heatmap) of PGMC, ALS, CTR based on 9-protein signature

run_heatmap_signature_CNS_IMMUNE(all_data,
                                 "SERUM",
                                 enet_optimal_protein_signature_serum_CNS$optimal_proteins,
                                 EN_PGMC_serum_signature_CNS$results,
                                 group_colors = group_colors,
                                 output_prefix = "plots/ML/Elastic Net/heatmap_SERUM_CNS",
                                 highlight_ids = participant_code_label,
                                 panel_mode ="CNS")

## -> Elastic Net + Plasma
EN_PGMC_plasma_signature_CNS = run_ALS_signature_workflow_CNS_IMMUNE(all_data,
                                                                 enet_optimal_protein_signature_plasma_CNS$optimal_proteins,
                                                                 fluid = "PLASMA",
                                                                 model_type = "elastic_net",
                                                                 output_prefix = "plots/ML/Elastic Net/CNS_",
                                                                 panel_mode = "CNS")

external_enet_plasma_wide_CNS <- build_external_ml_dataset(
  external_all_data, "PLASMA", "NPQ",
  external_covariate_table,
  corr_lookup_discovery = EN_PGMC_plasma_signature_CNS$corr_lookup,
  panel_zscore_params = EN_PGMC_plasma_signature_CNS$panel_zscore_params,
  panel_mode = "CNS"
)

external_validation_enet_plasma_CNS <- validate_external(
  cv_model        = EN_PGMC_plasma_signature_CNS$model,
  external_wide   = external_enet_plasma_wide_CNS,
  proteins        = enet_optimal_protein_signature_plasma_CNS$optimal_proteins,
  covariate_cols  = c("age", "sex_male"),
  scaling_center  = EN_PGMC_plasma_signature_CNS$scaling$center,
  scaling_scale   = EN_PGMC_plasma_signature_CNS$scaling$scale,
  fluid_label     = "PLASMA",
  plot_path       = "plots/ML/Elastic Net/ROC_external_validation_plasma_CNS.pdf"
)


participant_code_label = EN_PGMC_plasma_signature_CNS$results %>%
  left_join(protein_data_IDs %>% select(SampleName,ParticipantCode,type)) %>%
  distinct() %>%
  filter(type == "PGMC" & ALS_risk_score > 0) %>%
  pull(ParticipantCode)

# ================================================================================
# Unsupervised visualisation (heatmap) of PGMC, ALS, CTR based on 3-protein signature

run_heatmap_signature_CNS_IMMUNE(all_data,
                                 "PLASMA",
                                 enet_optimal_protein_signature_plasma_CNS$optimal_proteins,
                                 EN_PGMC_plasma_signature_CNS$results,
                                 group_colors = group_colors,
                                 output_prefix = "plots/ML/Elastic Net/heatmap_PLASMA_CNS",
                                 highlight_ids = participant_code_label,
                                 panel_mode = "CNS")

## -> Elastic Net + CSF
EN_PGMC_CSF_signature_CNS = run_ALS_signature_workflow_CNS_IMMUNE(all_data,
                                                              enet_optimal_protein_signature_CSF_CNS$optimal_proteins,
                                                              fluid = "CSF",
                                                              model_type = "elastic_net",
                                                              output_prefix = "plots/ML/Elastic Net/CNS_",
                                                              panel_mode = "CNS")

external_enet_csf_wide_CNS <- build_external_ml_dataset(
  external_all_data, "CSF", "NPQ",
  external_covariate_table,
  corr_lookup_discovery = EN_PGMC_CSF_signature_CNS$corr_lookup,
  panel_zscore_params = EN_PGMC_CSF_signature_CNS$panel_zscore_params,
  panel_mode = "CNS"
)

external_validation_enet_csf_CNS <- validate_external(
  cv_model        = EN_PGMC_CSF_signature_CNS$model,
  external_wide   = external_enet_csf_wide_CNS,
  proteins        = enet_optimal_protein_signature_CSF_CNS$optimal_proteins,
  covariate_cols  = c("age", "sex_male"),
  scaling_center  = EN_PGMC_CSF_signature_CNS$scaling$center,
  scaling_scale   = EN_PGMC_CSF_signature_CNS$scaling$scale,
  fluid_label     = "CSF",
  plot_path       = "plots/ML/Elastic Net/ROC_external_validation_CSF_CNS.pdf"
)

participant_code_label = EN_PGMC_CSF_signature_CNS$results %>%
  left_join(protein_data_IDs %>% select(SampleName,ParticipantCode,type)) %>%
  distinct() %>%
  filter(type == "PGMC" & ALS_risk_score > 0) %>%
  pull(ParticipantCode)

# ================================================================================
# Unsupervised visualisation (heatmap) of PGMC, ALS, CTR based on 9-protein signature

run_heatmap_signature_CNS_IMMUNE(all_data,
                                 "CSF",
                                 enet_optimal_protein_signature_CSF_CNS$optimal_proteins,
                                 EN_PGMC_CSF_signature_CNS$results,
                                 group_colors = group_colors,
                                 output_prefix = "plots/ML/Elastic Net/heatmap_CSF_CNS",
                                 highlight_ids = participant_code_label,
                                 panel_mode = "CNS")



## ==========================================================
## PART 6: ML -- Immune panel only per fluid
## ==========================================================

for (fluid in fluids) {
  
  meta_data = all_data %>% 
    select(SampleName,PatientID,ParticipantCode) %>%
    left_join(Sex_age_all_participants %>% dplyr::rename(PatientID = Pseudonyme)) %>%
    mutate(
      center = dplyr::case_when(
        grepl("TR", ParticipantCode) ~ "Turkey",
        grepl("CH", ParticipantCode) ~ "Switzerland",
        grepl("DE", ParticipantCode) ~ "Germany",
        grepl("SK", ParticipantCode) ~ "Slovakia",
        grepl("FR", ParticipantCode) ~ "France",
        grepl("IL", ParticipantCode) ~ "Israel",
        TRUE                 ~ NA_character_
      ))
  
  covariate_table = build_covariate_table(meta_data)
  
  ## Prepare merged, z-scored, correlation-averaged datasets
  protein_data_PGMCvsCTR_IMMUNE <- build_ml_dataset(all_data, fluid, c("PGMC", "CTR"), "NPQ",     "PGMC",
                                                 covariate_table,
                                                 panel_mode = "IMMUNE")
  protein_data_ALSvsCTR_IMMUNE  <- build_ml_dataset(all_data, fluid, c("ALS",  "CTR"), "NPQ", "ALS",
                                                 covariate_table,
                                                 panel_mode = "IMMUNE")
  protein_data_ALSvsPGMC_IMMUNE <- build_ml_dataset(all_data, fluid, c("ALS",  "PGMC"), "NPQ", "ALS",
                                                 covariate_table,
                                                 panel_mode = "IMMUNE")
  
  ## Run Lasso and ROC curve with 5-fold cv and 500 bootstrap iterations
  lm_PGMC_CTR_IMMUNE <- runML_with_lasso(protein_data_PGMCvsCTR_IMMUNE, bs_count = 500)
  final_roc_plot_PGMC_CTR_IMMUNE <- calculateROC(lm_PGMC_CTR_IMMUNE,
                                              paste0("plots/ML/Lasso/ROC_lm_IMMUNE_", fluid, "_PGMC_CTR.pdf"))
  
  lm_ALS_CTR_IMMUNE <- runML_with_lasso(protein_data_ALSvsCTR_IMMUNE, bs_count = 500)
  final_roc_plot_ALS_CTR_IMMUNE <- calculateROC(lm_ALS_CTR_IMMUNE,
                                             paste0("plots/ML/Lasso/ROC_lm_IMMUNE_", fluid, "_ALS_CTR.pdf"))
  
  lm_ALS_PGMC_IMMUNE <- runML_with_lasso(protein_data_ALSvsPGMC_IMMUNE, bs_count = 500)
  final_roc_plot_ALS_PGMC_IMMUNE <- calculateROC(lm_ALS_PGMC_IMMUNE,
                                              paste0("plots/ML/Lasso/ROC_lm_IMMUNE_", fluid, "_ALS_PGMC.pdf"))
  
  ## Extract weights + plot averaged
  lm_weights_PGMC_CTR_IMMUNE <- feature_importance(lm_PGMC_CTR_IMMUNE,
                                                plot_path = paste0("plots/ML/Lasso/protein_selection_IMMUNE_", fluid, "_PGMC_CTR.pdf"))
  write_xlsx(lm_weights_PGMC_CTR_IMMUNE, path = paste0("results/weights_lm_IMMUNE_", fluid, "_PGMC_CTR.xlsx"))
  
  lm_weights_ALS_CTR_IMMUNE <- feature_importance(lm_ALS_CTR_IMMUNE,
                                               plot_path = paste0("plots/ML/Lasso/protein_selection_IMMUNE_", fluid, "_ALS_CTR.pdf"))
  write_xlsx(lm_weights_ALS_CTR_IMMUNE, path = paste0("results/weights_lm_IMMUNE_", fluid, "_ALS_CTR.xlsx"))
  
  lm_weights_ALS_PGMC_IMMUNE <- feature_importance(lm_ALS_PGMC_IMMUNE,
                                                plot_path = paste0("plots/ML/Lasso/protein_selection_IMMUNE_", fluid, "_ALS_PGMC.pdf"))
  write_xlsx(lm_weights_ALS_PGMC_IMMUNE, path = paste0("results/weights_lm_IMMUNE_", fluid, "_ALS_PGMC.xlsx"))
  
  # Run Elastic Net and ROC curve with 5-fold cv and 500 bootstrap iterations
  lm_enet_PGMC_CTR_IMMUNE <- runML_with_elasticnet(protein_data_PGMCvsCTR_IMMUNE,bs_count = 500)
  final_roc_plot_PGMC_CTR_enet_IMMUNE = calculateROC(lm_enet_PGMC_CTR_IMMUNE,
                                                  paste0("plots/ML/Elastic Net/ROC_lm_IMMUNE_",fluid,"_PGMC_CTR.pdf"))
  lm_enet_ALS_CTR_IMMUNE <- runML_with_elasticnet(protein_data_ALSvsCTR_IMMUNE,bs_count = 500)
  final_roc_plot_ALS_CTR_enet_IMMUNE = calculateROC(lm_enet_ALS_CTR_IMMUNE,
                                                 paste0("plots/ML/Elastic Net/ROC_lm_IMMUNE_",fluid,"_ALS_CTR.pdf"))
  lm_enet_ALS_PGMC_IMMUNE <- runML_with_elasticnet(protein_data_ALSvsPGMC_IMMUNE,bs_count = 500)
  final_roc_plot_ALS_PGMC_enet_IMMUNE = calculateROC(lm_enet_ALS_PGMC_IMMUNE,
                                                  paste0("plots/ML/Elastic Net/ROC_lm_IMMUNE_",fluid,"_ALS_PGMC.pdf"))
  
  # extract weights + plot averaged
  lm_weights_PGMC_CTR_enet_IMMUNE = feature_importance(lm_enet_PGMC_CTR_IMMUNE,
                                                    plot_path = paste0("plots/ML/Elastic Net/protein_selection_IMMUNE_",fluid,"_PGMC_CTR.pdf"))
  write_xlsx(lm_weights_PGMC_CTR_enet_IMMUNE, path = paste0("results/weights_lm_enet_IMMUNE_",fluid,"_PGMC_CTR.xlsx"))
  lm_weights_ALS_CTR_enet_IMMUNE = feature_importance(lm_enet_ALS_CTR_IMMUNE,
                                                   plot_path = paste0("plots/ML/Elastic Net/protein_selection_IMMUNE_",fluid,"_ALS_CTR.pdf"))
  write_xlsx(lm_weights_ALS_CTR_enet_IMMUNE, path = paste0("results/weights_lm_enet_IMMUNE_",fluid,"_ALS_CTR.xlsx"))
  lm_weights_ALS_PGMC_enet_IMMUNE = feature_importance(lm_enet_ALS_PGMC_IMMUNE,
                                                    plot_path = paste0("plots/ML/Elastic Net/protein_selection_IMMUNE_",fluid,"_ALS_PGMC.pdf"))
  write_xlsx(lm_weights_ALS_PGMC_enet_IMMUNE, path = paste0("results/weights_lm_enet_IMMUNE_",fluid,"_ALS_PGMC.xlsx"))
}

# ================================================================================
# Protein signature (IMMUNE)
proteins_serum_ALS_CTR = c("CLEC4A_IMMUNE","HAVCR1_IMMUNE","IL18R1_IMMUNE","SELE_IMMUNE",
                           "TNFRSF13C_IMMUNE","CD200_IMMUNE", "CTLA4_IMMUNE",
                           "SIRPA_IMMUNE","NGF_IMMUNE","GZMB_IMMUNE","BST2_IMMUNE")
lasso_optimal_protein_signature_serum_IMMUNE = find_optimal_signature_CNS_IMMUNE(data_all = build_ml_dataset(all_data, 
                                                                                                      "SERUM", 
                                                                                                      c("ALS",  "CTR"), 
                                                                                                      "NPQ", "ALS",
                                                                                                      covariate_table,
                                                                                                      panel_mode = "IMMUNE"),
                                                                          ranked_proteins = proteins_serum_ALS_CTR,
                                                                          fluid = "SERUM",
                                                                          output_prefix = "plots/ML/Lasso/performance_IMMUNE_")

## -> Lasso + Plasma
proteins_plasma_ALS_CTR = c("IL1RL1_IMMUNE","CD200_IMMUNE","HAVCR1_IMMUNE","CXCL10_IMMUNE",
                            "CCL25_IMMUNE","TNFSF9_IMMUNE","TSLP_IMMUNE","TNFRSF13C_IMMUNE",
                            "AGRP_IMMUNE","CD3E_IMMUNE","SIRPA_IMMUNE")
lasso_optimal_protein_signature_plasma_IMMUNE = find_optimal_signature_CNS_IMMUNE(data_all = build_ml_dataset(all_data, 
                                                                                                       "PLASMA", 
                                                                                                       c("ALS",  "CTR"), 
                                                                                                       "NPQ", "ALS",
                                                                                                       covariate_table,
                                                                                                       panel_mode = "IMMUNE"),
                                                                           ranked_proteins = proteins_plasma_ALS_CTR,
                                                                           fluid = "PLASMA",
                                                                           output_prefix = "plots/ML/Lasso/performance_IMMUNE_")

# Elastic Net + Serum
proteins_serum_ALS_CTR = c("CLEC4A_IMMUNE","HAVCR1_IMMUNE","SELE_IMMUNE","TNFRSF13C_IMMUNE",
                           "IL18R1_IMMUNE", "CTLA4_IMMUNE","CD200_IMMUNE","CD276_IMMUNE",
                           "SIRPA_IMMUNE","NGF_IMMUNE","GZMB_IMMUNE")
enet_optimal_protein_signature_serum_IMMUNE = find_optimal_signature_CNS_IMMUNE(data_all = build_ml_dataset(all_data, 
                                                                                                     "SERUM", 
                                                                                                     c("ALS",  "CTR"), 
                                                                                                     "NPQ", "ALS",
                                                                                                     covariate_table,
                                                                                                     panel_mode = "IMMUNE"),
                                                                         ranked_proteins = proteins_serum_ALS_CTR,
                                                                         fluid = "SERUM",
                                                                         output_prefix = "plots/ML/Elastic Net/performance_IMMUNE_")

# Elastic Net + Plasma
proteins_plasma_ALS_CTR = c("IL1RL1_IMMUNE","CXCL10_IMMUNE","CCL25_IMMUNE","CD200_IMMUNE",
                            "TNFSF9_IMMUNE","SELE_IMMUNE","HAVCR1_IMMUNE","BST2_IMMUNE",
                            "PGF_IMMUNE","TSLP_IMMUNE","AGRP_IMMUNE")
enet_optimal_protein_signature_plasma_IMMUNE = find_optimal_signature_CNS_IMMUNE(data_all = build_ml_dataset(all_data, 
                                                                                                      "PLASMA", 
                                                                                                      c("ALS",  "CTR"), 
                                                                                                      "NPQ", "ALS",
                                                                                                      covariate_table,
                                                                                                      panel_mode = "IMMUNE"),
                                                                          ranked_proteins = proteins_plasma_ALS_CTR,
                                                                          fluid = "PLASMA",
                                                                          output_prefix = "plots/ML/Elastic Net/performance_IMMUNE_")

# Elastic Net + CSF
proteins_CSF_ALS_CTR =  c("CCL28_IMMUNE","EPO_IMMUNE","PDCD1LG2_IMMUNE","IL18_IMMUNE",
                          "IL27_IMMUNE","IL11_IMMUNE","CSF3R_IMMUNE","TNFRSF13C_IMMUNE",
                          "IL19_IMMUNE","CXCL11_IMMUNE","CHI3L1_IMMUNE")

enet_optimal_protein_signature_CSF_IMMUNE = find_optimal_signature_CNS_IMMUNE(data_all = build_ml_dataset(all_data, 
                                                                                                   "CSF", 
                                                                                                   c("ALS",  "CTR"), 
                                                                                                   "NPQ", "ALS",
                                                                                                   covariate_table,
                                                                                                   panel_mode = "IMMUNE"),
                                                                       ranked_proteins = proteins_CSF_ALS_CTR,
                                                                       fluid = "CSF",
                                                                       output_prefix = "plots/ML/Elastic Net/performance_IMMUNE_")

# ====================================================
# Compute an ALS risk score based on ALS vs CTR model 
## -> Lasso + Serum
proteins_serum_ALS_CTR_IMMUNE = lasso_optimal_protein_signature_serum_IMMUNE$optimal_proteins

lasso_PGMC_serum_signature_IMMUNE = run_ALS_signature_workflow_CNS_IMMUNE(all_data,
                                                                   proteins_serum_ALS_CTR_IMMUNE,
                                                                   fluid = "SERUM",
                                                                   model_type = "lasso",
                                                                   output_prefix = "plots/ML/Lasso/IMMUNE_",
                                                                   panel_mode = "IMMUNE")

external_lasso_serum_wide_IMMUNE <- build_external_ml_dataset(
  external_all_data, "SERUM", "NPQ",
  external_covariate_table,
  corr_lookup_discovery = lasso_PGMC_serum_signature_IMMUNE$corr_lookup,
  panel_zscore_params = lasso_PGMC_serum_signature_IMMUNE$panel_zscore_params,
  panel_mode = "IMMUNE"
)

external_validation_lasso_serum_IMMUNE <- validate_external(
  cv_model        = lasso_PGMC_serum_signature_IMMUNE$model,
  external_wide   = external_lasso_serum_wide_IMMUNE,
  proteins        = lasso_optimal_protein_signature_serum_IMMUNE$optimal_proteins,
  covariate_cols  = c("age", "sex_male"),
  scaling_center  = lasso_PGMC_serum_signature_IMMUNE$scaling$center,
  scaling_scale   = lasso_PGMC_serum_signature_IMMUNE$scaling$scale,
  fluid_label     = "SERUM",
  plot_path       = "plots/ML/Lasso/ROC_external_validation_SERUM_IMMUNE.pdf"
)

participant_code_label = lasso_PGMC_serum_signature_IMMUNE$results %>%
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
                                 lasso_optimal_protein_signature_serum_IMMUNE$optimal_proteins,
                                 lasso_PGMC_serum_signature_IMMUNE$results,
                                 group_colors = group_colors,
                                 output_prefix = "plots/ML/Lasso/heatmap_SERUM_IMMUNE",
                                 highlight_ids = participant_code_label,
                                 panel_mode = "IMMUNE")

## -> Lasso + plasma
proteins_plasma_ALS_CTR_IMMUNE = lasso_optimal_protein_signature_plasma_IMMUNE$optimal_proteins

lasso_PGMC_plasma_signature_IMMUNE = run_ALS_signature_workflow_CNS_IMMUNE(all_data,
                                                                    proteins_plasma_ALS_CTR_IMMUNE,
                                                                    fluid = "PLASMA",
                                                                    model_type = "lasso",
                                                                    output_prefix = "plots/ML/Lasso/IMMUNE_",
                                                                    panel_mode = "IMMUNE")

external_lasso_plasma_wide_IMMUNE <- build_external_ml_dataset(
  external_all_data, "PLASMA", "NPQ",
  external_covariate_table,
  corr_lookup_discovery = lasso_PGMC_plasma_signature_IMMUNE$corr_lookup,
  panel_zscore_params = lasso_PGMC_plasma_signature_IMMUNE$panel_zscore_params,
  panel_mode = "IMMUNE"
)

external_validation_lasso_plasma_IMMUNE <- validate_external(
  cv_model        = lasso_PGMC_plasma_signature_IMMUNE$model,
  external_wide   = external_lasso_plasma_wide_IMMUNE,
  proteins        = lasso_optimal_protein_signature_plasma_IMMUNE$optimal_proteins,
  covariate_cols  = c("age", "sex_male"),
  scaling_center  = lasso_PGMC_plasma_signature_IMMUNE$scaling$center,
  scaling_scale   = lasso_PGMC_plasma_signature_IMMUNE$scaling$scale,
  fluid_label     = "PLASMA",
  plot_path       = "plots/ML/Lasso/ROC_external_validation_plasma_IMMUNE.pdf"
)

participant_code_label = lasso_PGMC_plasma_signature_IMMUNE$results %>%
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
                                 "PLASMA",
                                 lasso_optimal_protein_signature_plasma_IMMUNE$optimal_proteins,
                                 lasso_PGMC_plasma_signature_IMMUNE$results,
                                 group_colors = group_colors,
                                 output_prefix = "plots/ML/Lasso/heatmap_plasma_IMMUNE",
                                 highlight_ids = participant_code_label,
                                 panel_mode = "IMMUNE"
)


##### ----
# Elastic Net

## -> Elastic Net + Serum
EN_PGMC_serum_signature_IMMUNE = run_ALS_signature_workflow_CNS_IMMUNE(all_data,
                                                                enet_optimal_protein_signature_serum_IMMUNE$optimal_proteins,
                                                                fluid = "SERUM",
                                                                model_type = "elastic_net",
                                                                output_prefix = "plots/ML/Elastic Net/IMMUNE_",
                                                                panel_mode = "IMMUNE")

external_enet_serum_wide_IMMUNE <- build_external_ml_dataset(
  external_all_data, "SERUM", "NPQ",
  external_covariate_table,
  corr_lookup_discovery = EN_PGMC_serum_signature_IMMUNE$corr_lookup,
  panel_zscore_params = EN_PGMC_serum_signature_IMMUNE$panel_zscore_params,
  panel_mode = "IMMUNE"
)

external_validation_enet_serum_IMMUNE <- validate_external(
  cv_model        = EN_PGMC_serum_signature_IMMUNE$model,
  external_wide   = external_enet_serum_wide_IMMUNE,
  proteins        = enet_optimal_protein_signature_serum_IMMUNE$optimal_proteins,
  covariate_cols  = c("age", "sex_male"),
  scaling_center  = EN_PGMC_serum_signature_IMMUNE$scaling$center,
  scaling_scale   = EN_PGMC_serum_signature_IMMUNE$scaling$scale,
  fluid_label     = "SERUM",
  plot_path       = "plots/ML/Elastic Net/ROC_external_validation_SERUM_IMMUNE.pdf"
)


participant_code_label = EN_PGMC_serum_signature_IMMUNE$results %>%
  left_join(protein_data_IDs %>% select(SampleName,ParticipantCode,type)) %>%
  distinct() %>%
  filter(type == "PGMC" & ALS_risk_score > 0) %>%
  pull(ParticipantCode)

# ================================================================================
# Unsupervised visualisation (heatmap) of PGMC, ALS, CTR based on 9-protein signature

run_heatmap_signature_CNS_IMMUNE(all_data,
                                 "SERUM",
                                 enet_optimal_protein_signature_serum_IMMUNE$optimal_proteins,
                                 EN_PGMC_serum_signature_IMMUNE$results,
                                 group_colors = group_colors,
                                 output_prefix = "plots/ML/Elastic Net/heatmap_SERUM_IMMUNE",
                                 highlight_ids = participant_code_label,
                                 panel_mode ="IMMUNE")

## -> Elastic Net + Plasma
EN_PGMC_plasma_signature_IMMUNE = run_ALS_signature_workflow_CNS_IMMUNE(all_data,
                                                                 enet_optimal_protein_signature_plasma_IMMUNE$optimal_proteins,
                                                                 fluid = "PLASMA",
                                                                 model_type = "elastic_net",
                                                                 output_prefix = "plots/ML/Elastic Net/IMMUNE_",
                                                                 panel_mode = "IMMUNE")

external_enet_plasma_wide_IMMUNE <- build_external_ml_dataset(
  external_all_data, "PLASMA", "NPQ",
  external_covariate_table,
  corr_lookup_discovery = EN_PGMC_plasma_signature_IMMUNE$corr_lookup,
  panel_zscore_params = EN_PGMC_plasma_signature_IMMUNE$panel_zscore_params,
  panel_mode = "IMMUNE"
)

external_validation_enet_plasma_IMMUNE <- validate_external(
  cv_model        = EN_PGMC_plasma_signature_IMMUNE$model,
  external_wide   = external_enet_plasma_wide_IMMUNE,
  proteins        = enet_optimal_protein_signature_plasma_IMMUNE$optimal_proteins,
  covariate_cols  = c("age", "sex_male"),
  scaling_center  = EN_PGMC_plasma_signature_IMMUNE$scaling$center,
  scaling_scale   = EN_PGMC_plasma_signature_IMMUNE$scaling$scale,
  fluid_label     = "PLASMA",
  plot_path       = "plots/ML/Elastic Net/ROC_external_validation_plasma_IMMUNE.pdf"
)


participant_code_label = EN_PGMC_plasma_signature_IMMUNE$results %>%
  left_join(protein_data_IDs %>% select(SampleName,ParticipantCode,type)) %>%
  distinct() %>%
  filter(type == "PGMC" & ALS_risk_score > 0) %>%
  pull(ParticipantCode)

# ================================================================================
# Unsupervised visualisation (heatmap) of PGMC, ALS, CTR based on 3-protein signature

run_heatmap_signature_CNS_IMMUNE(all_data,
                                 "PLASMA",
                                 enet_optimal_protein_signature_plasma_IMMUNE$optimal_proteins,
                                 EN_PGMC_plasma_signature_IMMUNE$results,
                                 group_colors = group_colors,
                                 output_prefix = "plots/ML/Elastic Net/heatmap_PLASMA_IMMUNE",
                                 highlight_ids = participant_code_label,
                                 panel_mode = "IMMUNE")

## -> Elastic Net + CSF
EN_PGMC_CSF_signature_IMMUNE = run_ALS_signature_workflow_CNS_IMMUNE(all_data,
                                                              enet_optimal_protein_signature_CSF_IMMUNE$optimal_proteins,
                                                              fluid = "CSF",
                                                              model_type = "elastic_net",
                                                              output_prefix = "plots/ML/Elastic Net/IMMUNE_",
                                                              panel_mode = "IMMUNE")

external_enet_csf_wide_IMMUNE <- build_external_ml_dataset(
  external_all_data, "CSF", "NPQ",
  external_covariate_table,
  corr_lookup_discovery = EN_PGMC_CSF_signature_IMMUNE$corr_lookup,
  panel_zscore_params = EN_PGMC_CSF_signature_IMMUNE$panel_zscore_params,
  panel_mode = "IMMUNE"
)

external_validation_enet_csf_IMMUNE <- validate_external(
  cv_model        = EN_PGMC_CSF_signature_IMMUNE$model,
  external_wide   = external_enet_csf_wide_IMMUNE,
  proteins        = enet_optimal_protein_signature_CSF_IMMUNE$optimal_proteins,
  covariate_cols  = c("age", "sex_male"),
  scaling_center  = EN_PGMC_CSF_signature_IMMUNE$scaling$center,
  scaling_scale   = EN_PGMC_CSF_signature_IMMUNE$scaling$scale,
  fluid_label     = "CSF",
  plot_path       = "plots/ML/Elastic Net/ROC_external_validation_CSF_IMMUNE.pdf"
)

participant_code_label = EN_PGMC_CSF_signature_IMMUNE$results %>%
  left_join(protein_data_IDs %>% select(SampleName,ParticipantCode,type)) %>%
  distinct() %>%
  filter(type == "PGMC" & ALS_risk_score > 0) %>%
  pull(ParticipantCode)

# ================================================================================
# Unsupervised visualisation (heatmap) of PGMC, ALS, CTR based on 9-protein signature

run_heatmap_signature_CNS_IMMUNE(all_data,
                                 "CSF",
                                 enet_optimal_protein_signature_CSF_IMMUNE$optimal_proteins,
                                 EN_PGMC_CSF_signature_IMMUNE$results,
                                 group_colors = group_colors,
                                 output_prefix = "plots/ML/Elastic Net/heatmap_CSF_IMMUNE",
                                 highlight_ids = participant_code_label,
                                 panel_mode = "IMMUNE")

########################################  
## EXTERNAL BOXPLOTS of MAXOMOD

get_external_raw_values <- function(external_long, fluid, proteins, panel_mode = "both",
                                    corr_lookup_discovery = NULL) {
  
  fluid_data <- external_long %>% filter(SampleMatrixType == fluid)
  
  raw_panel_wide <- function(panel_name, suffix) {
    fluid_data %>%
      filter(panel == panel_name) %>%
      select(SampleName, Target, NPQ) %>%
      group_by(SampleName, Target) %>%
      summarise(NPQ = mean(NPQ, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = Target, values_from = NPQ) %>%
      rename_with(~ paste0(.x, suffix), .cols = -SampleName)
  }
  
  if (panel_mode == "both") {
    cns_raw    <- raw_panel_wide("CNS", "_CNS")
    immune_raw <- raw_panel_wide("IMMUNE", "_IMMUNE")
    wide <- full_join(cns_raw, immune_raw, by = "SampleName")
    
    if (!is.null(corr_lookup_discovery) && nrow(corr_lookup_discovery) > 0) {
      for (i in seq_len(nrow(corr_lookup_discovery))) {
        col_cns     <- paste0(corr_lookup_discovery$target_CNS[i],    "_CNS")
        col_immune  <- paste0(corr_lookup_discovery$target_IMMUNE[i], "_IMMUNE")
        merged_name <- corr_lookup_discovery$merged_name[i]
        
        if (all(c(col_cns, col_immune) %in% names(wide))) {
          wide[[merged_name]] <- rowMeans(wide[, c(col_cns, col_immune)], na.rm = TRUE)
        }
      }
    }
  } else {
    suffix <- paste0("_", panel_mode)
    wide <- raw_panel_wide(panel_mode, suffix)
  }
  
  available <- intersect(proteins, names(wide))
  missing   <- setdiff(proteins, names(wide))
  if (length(missing) > 0) {
    message("No raw NPQ found for: ", paste(missing, collapse = ", "),
            " in external ", fluid, " data -- check panel coverage / naming.")
  }
  
  wide %>% select(SampleName, all_of(available))
}


## Long-format (protein, NPQ, group) data ready for ggplot
build_boxplot_data <- function(external_long, fluid, proteins, panel_mode,
                               corr_lookup_discovery, covariate_table) {
  
  raw_wide <- get_external_raw_values(external_long, fluid, proteins, panel_mode, corr_lookup_discovery)
  value_cols <- intersect(proteins, names(raw_wide))
  
  raw_wide %>%
    inner_join(covariate_table %>% select(SampleName, status), by = "SampleName") %>%
    mutate(group = ifelse(status == 1, "ALS", "CTR")) %>%
    pivot_longer(cols = all_of(value_cols), names_to = "protein", values_to = "NPQ")
}


## Draw the boxplots, one facet per protein
plot_external_signature_boxplots <- function(external_long, fluid, proteins, panel_mode,
                                             corr_lookup_discovery, covariate_table,
                                             plot_path = NULL,
                                             group_colors = c("CTR" = "#6F8EB2", "ALS" = "#B2936F"),
                                             max_ncol = 6) {
  
  plot_df <- build_boxplot_data(external_long, fluid, proteins, panel_mode,
                                corr_lookup_discovery, covariate_table)
  
  n_proteins <- length(unique(plot_df$protein))
  
  ncol_facet <- min(max_ncol, ceiling(sqrt(n_proteins)))
  nrow_facet <- ceiling(n_proteins / ncol_facet)
  
  ## wilcoxon test
  stats_table <- plot_df %>%
    group_by(protein) %>%
    summarise(
      p_value = tryCatch(wilcox.test(NPQ ~ group)$p.value, error = function(e) NA_real_),
      .groups = "drop"
    ) %>%
    mutate(
      p_adj = p.adjust(p_value, method = "BH"),
      signif_label = case_when(
        p_adj < 0.001 ~ "***",
        p_adj < 0.01  ~ "**",
        p_adj < 0.05  ~ "*",
        TRUE          ~ "ns"
      )
    )
  
  use_ggpubr <- requireNamespace("ggpubr", quietly = TRUE)
  
  p <- ggplot(plot_df, aes(x = group, y = NPQ, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.6) +
    geom_jitter(width = 0.15, size = 0.9, alpha = 0.4, color = "black") +
    facet_wrap(~ protein, scales = "free_y", ncol = ncol_facet)
  
  if (use_ggpubr) {
    p <- p + ggpubr::stat_compare_means(method = "wilcox.test", label = "p.format",
                                        label.x.npc = "center", size = 3)
  } else {
    label_df <- plot_df %>%
      group_by(protein) %>%
      summarise(y_pos = max(NPQ, na.rm = TRUE) * 1.05, .groups = "drop") %>%
      left_join(stats_table, by = "protein")
    p <- p + geom_text(data = label_df,
                       aes(x = 1.5, y = y_pos, label = paste0("p = ", signif(p_value, 2))),
                       inherit.aes = FALSE, size = 3)
  }
  
  p <- p +
    scale_fill_manual(values = group_colors) +
    labs(
      title = paste0("External validation cohort (MAXOMOD) in ", fluid),
      x = NULL, y = "NPQ (raw)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold")
    )
  
  if (!is.null(plot_path)) {
    ggsave(plot_path, p, width = 3.2 * ncol_facet, height = 3.4 * nrow_facet, limitsize = FALSE)
    write.csv(stats_table, gsub("\\.pdf$", "_stats.csv", plot_path), row.names = FALSE)
  }
  
  list(plot = p, stats = stats_table)
}


## gt protein signature
signature_registry <- list(
  list(label = "Lasso_SERUM_both",   proteins = lasso_optimal_protein_signature_serum_both$optimal_proteins,   model_result = lasso_PGMC_serum_signature_both,   fluid = "SERUM",  panel_mode = "both",   out_dir = "plots/ML/Lasso"),
  list(label = "Lasso_PLASMA_both",  proteins = lasso_optimal_protein_signature_plasma_both$optimal_proteins,  model_result = lasso_PGMC_plasma_signature_both,  fluid = "PLASMA", panel_mode = "both",   out_dir = "plots/ML/Lasso"),
  list(label = "EN_SERUM_both",      proteins = enet_optimal_protein_signature_serum_both$optimal_proteins,    model_result = EN_PGMC_serum_signature_both,      fluid = "SERUM",  panel_mode = "both",   out_dir = "plots/ML/Elastic Net"),
  list(label = "EN_PLASMA_both",     proteins = enet_optimal_protein_signature_plasma_both$optimal_proteins,   model_result = EN_PGMC_plasma_signature_both,     fluid = "PLASMA", panel_mode = "both",   out_dir = "plots/ML/Elastic Net"),
  list(label = "EN_CSF_both",        proteins = enet_optimal_protein_signature_CSF_both$optimal_proteins,      model_result = EN_PGMC_CSF_signature_both,        fluid = "CSF",    panel_mode = "both",   out_dir = "plots/ML/Elastic Net"),
  
  list(label = "Lasso_SERUM_CNS",    proteins = lasso_optimal_protein_signature_serum_CNS$optimal_proteins,    model_result = lasso_PGMC_serum_signature_CNS,    fluid = "SERUM",  panel_mode = "CNS",    out_dir = "plots/ML/Lasso"),
  list(label = "Lasso_PLASMA_CNS",   proteins = lasso_optimal_protein_signature_plasma_CNS$optimal_proteins,   model_result = lasso_PGMC_plasma_signature_CNS,   fluid = "PLASMA", panel_mode = "CNS",    out_dir = "plots/ML/Lasso"),
  list(label = "EN_SERUM_CNS",       proteins = enet_optimal_protein_signature_serum_CNS$optimal_proteins,     model_result = EN_PGMC_serum_signature_CNS,       fluid = "SERUM",  panel_mode = "CNS",    out_dir = "plots/ML/Elastic Net"),
  list(label = "EN_PLASMA_CNS",      proteins = enet_optimal_protein_signature_plasma_CNS$optimal_proteins,    model_result = EN_PGMC_plasma_signature_CNS,      fluid = "PLASMA", panel_mode = "CNS",    out_dir = "plots/ML/Elastic Net"),
  list(label = "EN_CSF_CNS",         proteins = enet_optimal_protein_signature_CSF_CNS$optimal_proteins,       model_result = EN_PGMC_CSF_signature_CNS,         fluid = "CSF",    panel_mode = "CNS",    out_dir = "plots/ML/Elastic Net"),
  
  list(label = "Lasso_SERUM_IMMUNE", proteins = lasso_optimal_protein_signature_serum_IMMUNE$optimal_proteins, model_result = lasso_PGMC_serum_signature_IMMUNE, fluid = "SERUM",  panel_mode = "IMMUNE", out_dir = "plots/ML/Lasso"),
  list(label = "Lasso_PLASMA_IMMUNE",proteins = lasso_optimal_protein_signature_plasma_IMMUNE$optimal_proteins,model_result = lasso_PGMC_plasma_signature_IMMUNE,fluid = "PLASMA", panel_mode = "IMMUNE", out_dir = "plots/ML/Lasso"),
  list(label = "EN_SERUM_IMMUNE",    proteins = enet_optimal_protein_signature_serum_IMMUNE$optimal_proteins,  model_result = EN_PGMC_serum_signature_IMMUNE,    fluid = "SERUM",  panel_mode = "IMMUNE", out_dir = "plots/ML/Elastic Net"),
  list(label = "EN_PLASMA_IMMUNE",   proteins = enet_optimal_protein_signature_plasma_IMMUNE$optimal_proteins, model_result = EN_PGMC_plasma_signature_IMMUNE,   fluid = "PLASMA", panel_mode = "IMMUNE", out_dir = "plots/ML/Elastic Net"),
  list(label = "EN_CSF_IMMUNE",      proteins = enet_optimal_protein_signature_CSF_IMMUNE$optimal_proteins,    model_result = EN_PGMC_CSF_signature_IMMUNE,      fluid = "CSF",    panel_mode = "IMMUNE", out_dir = "plots/ML/Elastic Net")
)

## one pdf pr fluid
union_proteins_by_fluid <- function(fluid_name) {
  entries <- Filter(function(e) e$fluid == fluid_name, signature_registry)
  unique(unlist(lapply(entries, function(e) e$proteins)))
}

## correlation of protins across panels
get_both_corr_lookup <- function(fluid_name) {
  entry <- Filter(function(e) e$fluid == fluid_name && e$panel_mode == "both", signature_registry)[[1]]
  entry$model_result$corr_lookup
}

fluid_boxplot_results <- list()
for (fl in c("SERUM", "PLASMA", "CSF")) {
  
  union_proteins <- union_proteins_by_fluid(fl)
  message(fl, ": ", length(union_proteins), " unique proteins across all signatures")
  
  fluid_boxplot_results[[fl]] <- tryCatch({
    plot_external_signature_boxplots(
      external_all_data, fl,
      proteins = union_proteins,
      panel_mode = "both",   # superset -- correctly resolves suffixed AND merged names
      corr_lookup_discovery = get_both_corr_lookup(fl),
      covariate_table = external_covariate_table,
      plot_path = paste0("plots/ML/external_boxplots_", fl, "_all_signature_proteins.pdf")
    )
  }, error = function(e) {
    message("  FAILED for ", fl, ": ", conditionMessage(e))
    NULL
  })
}

fluid_boxplot_results$SERUM$stats
fluid_boxplot_results$PLASMA$stats
fluid_boxplot_results$CSF$stats