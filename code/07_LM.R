###############################################
### Helper Functions
###############################################

# Min-Max Scaling: brings everything to a 0-1 range
scale_manual <- function(df) {
  # We exclude 'status' from scaling
  status_col <- df$status
  numeric_df <- df %>% dplyr::select(-status)
  scaled_df <- as.data.frame(apply(numeric_df, 2, function(x) (x - min(x)) / diff(range(x))))
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
    actuals_list[[i]] = data_frame[ -train_ind,]$status
    
    if (algorithm == 'lm'){
      model = train(x=train[ , !names(train) == "status"],
                    y= train$status,
                    method = "glmnet", family = "binomial", tuneLength = 5, metric = "ROC",
                    trControl = ctrl,
                    tuneGrid=expand.grid(
                      .alpha=1, # alpha 1 == lasso
                      .lambda=10^seq(-4, 0, length.out = 20))
      )
    }
    
    else if (algorithm == 'enet'){
      model = train(x=train[ , !names(train) == "status"],
                    y= train$status,
                    method = "glmnet",
                    family = "binomial",
                    metric = "ROC",
                    trControl = ctrl,
                    tuneGrid = expand.grid(
                      .alpha = seq(0, 1, length.out = 10),  # elastic net mix
                      .lambda = 10^seq(-4, 0, length.out = 20))
      )
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
    
    models[[i]] = model
    importance[[i]] = varImp(model)
    predictions[[i]] = predict(model, newdata = test,type = "prob")
    predictions_raw[[i]] = predict(model, newdata = test,type = "raw")
    indices[[i]] = train_ind
    print(paste0('finshed loop #',i)) # keep track of what is happening
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
calculateROC_old <- function(list_from_ML, plot_path = FALSE) {
 
  # Process ROC points from each bootstrap "test" set
  all_roc_points <- list()
  
  for(i in seq_along(list_from_ML$predictions)) {
    # Extract actuals for the samples not used in training 
    actuals <- list_from_ML$actuals[[i]]
    # Get probabilities (Class X1)
    probs <- list_from_ML$predictions[[i]]$X1
    
    if(length(unique(actuals)) > 1) {
      all_roc_points[[i]] <- get_roc_points(actuals, probs, i)
    }
  }
  
  roc_df <- bind_rows(all_roc_points)
  
  # Summarize for Plotting
  summary_roc <- roc_df %>%
    dplyr::mutate(FPR = round(FPR, 4)) %>%
    group_by(FPR) %>%
    dplyr::summarise(
      mean_TPR  = mean(TPR),
      lower_TPR = quantile(TPR, 0.025),
      upper_TPR = quantile(TPR, 0.975),
      .groups = "drop"
    )
  
  auc_vec <- roc_df %>%
    dplyr::group_by(run) %>%
    dplyr::summarise(a = dplyr::first(auc), .groups = "drop") %>%
    dplyr::pull(a)
   
  p <- ggplot(summary_roc, aes(x = FPR, y = mean_TPR)) + 
    geom_ribbon(aes(ymin = lower_TPR, ymax = upper_TPR), fill = "grey70", alpha = 0.4) +
    geom_line(color = "darkblue", size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color="darkgrey") +
    theme_bw() +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(
      x = "1 - Specificity",
      y = "Sensitivity",
      title = "Mean of LASSO Bootstrap ROC curve",
      subtitle = paste0(
        "Mean AUC: ", round(mean(auc_vec), 3),
        " (\u00B1 ", round(sd(auc_vec), 3), "), 95% CI: ",
        round(quantile(auc_vec, 0.025), 3), "-",
        round(quantile(auc_vec, 0.975), 3)
      )
    ) +
    coord_equal()
  
  if(plot_path != FALSE){
    ggsave(plot_path, p, width = 7, height = 7, device = "pdf")
  }
  
  return(list(plot = p, roc_data = summary_roc, auc_values = auc_vec))
}

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
                              number = 20) {
  
  if (!dir.exists(dirname(plot_path))) dir.create(dirname(plot_path), recursive = TRUE)
  
  imp_list     <- ml_object$importance
  num_runs     <- length(imp_list)
  all_proteins <- unique(unlist(lapply(imp_list, function(x) rownames(x$importance))))
  
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
  
  return(weight)
}

# ================================================
# Check the best protein signature based on LOOCV
# ================================================
find_optimal_signature <- function(results_ALL, 
                                   ranked_proteins, 
                                   fluid = "SERUM", 
                                   output_prefix = "optimization_",
                                   nfolds = NULL) {
  
  # data prep
  data_train <- results_ALL[[fluid]]$data_adjusted %>%
    filter(type %in% c("ALS", "CTR"), Target %in% ranked_proteins) %>%
    select(SampleName, Target, NPQ_adj, type) %>%
    pivot_wider(names_from = Target, values_from = NPQ_adj)
  
  X <- as.matrix(data_train[, ranked_proteins])
  y <- ifelse(data_train$type == "ALS", 1, 0)
  X_scaled <- scale(X)
  
  if (is.null(nfolds)) nfolds <- length(y)
  
  # performance for each protein in the list
  perf_results <- data.frame(n_proteins = 1:length(ranked_proteins), 
                             protein_added = ranked_proteins,
                             auc = NA)
  
  message("Optimizing signature size...")
  
  for (i in 1:length(ranked_proteins)) {
    current_vars <- ranked_proteins[1:i]
    X_sub <- X_scaled[, current_vars, drop = FALSE]
    
    if (i == 1) {
      probs_vec <- numeric(length(y)) 
      
      for (j in 1:length(y)) {
        train_X <- X_sub[-j, , drop = FALSE]
        train_y <- y[-j]
        test_X  <- X_sub[j, , drop = FALSE]
        
        # Convert to data frame for glm
        df_train <- data.frame(status = train_y, protein = as.numeric(train_X))
        df_test  <- data.frame(protein = as.numeric(test_X))
        
        fit <- suppressWarnings(glm(status ~ protein, data = df_train, family = "binomial"))
        probs_vec[j] <- predict(fit, newdata = df_test, type = "response")
      }
      
      # Calculate AUC
      perf_results$auc[i] <- as.numeric(auc(roc(y, probs_vec, quiet = TRUE, levels=c(0,1), direction="<")))
      
    } else {
      set.seed(123)
      cv_fit <- cv.glmnet(X_sub, y, family = "binomial", alpha = 0, 
                          type.measure = "auc", nfolds = nfolds, keep = TRUE)
      best_lambda_idx <- which(cv_fit$lambda == cv_fit$lambda.min)
      cv_probs <- cv_fit$fit.preval[, best_lambda_idx]
      perf_results$auc[i] <- as.numeric(auc(roc(y, cv_probs, quiet = TRUE, levels=c(0,1), direction="<")))
    }
  }
  
  # check the oprimal set of proteins
  if(nrow(perf_results) > 1) {
    best_n <- perf_results$n_proteins[which.max(perf_results$auc[-1]) + 1]
  } else {
    best_n <- 1 
  }
  
  optimal_proteins <- ranked_proteins[1:best_n]
  
  # plot number of proteins versus performance
  p <- ggplot(perf_results, aes(x = n_proteins, y = auc)) +
    geom_line(color = "#ad5291", size = 1) +
    geom_point(aes(color = (n_proteins == best_n)), size = 3) +
    geom_vline(xintercept = best_n, linetype = "dashed", color = "gray50") +
    scale_x_continuous(breaks = 1:length(ranked_proteins), labels = ranked_proteins) +
    scale_color_manual(values = c("black", "red"), guide = "none") +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Selection of best protein signature based on model performance",
         subtitle = paste0("Max AUC: ", round(max(perf_results$auc), 3), " at ", best_n, " proteins"),
         x = "Proteins added based on feature importance", y = "AUC of LOOCV")
  
  ggsave(paste0(output_prefix, "protein_signature_", fluid, ".pdf"), p, width = 8, height = 6)
  
  return(list(
    optimal_n = best_n,
    optimal_proteins = optimal_proteins,
    performance_table = perf_results
  ))
}

# =============================
# Get ALS risk score in PGMC
# =============================
run_ALS_signature_workflow <- function(
    results_ALL,
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
  data_train <- results_ALL[[fluid]]$data_adjusted %>%
    filter(type %in% c("ALS","CTR"), Target %in% proteins) %>%
    select(SampleName, Target, NPQ_adj, type) %>%
    pivot_wider(names_from = Target, values_from = NPQ_adj)
  
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
  data_pgmc <- results_ALL[[fluid]]$data_adjusted %>%
    filter(type == "PGMC", Target %in% proteins) %>%
    select(SampleName, Target, NPQ_adj) %>%
    pivot_wider(names_from = Target, values_from = NPQ_adj)
  
  X_pgmc <- as.matrix(data_pgmc[, proteins])
  X_pgmc_scaled <- scale(X_pgmc,
                         center = scaling_center,
                         scale = scaling_scale)
  
  pred_pgmc <- predict(cv_model,
                       newx = X_pgmc_scaled,
                       s = "lambda.min",
                       type = "link")
  
  PGMC_results <- data_pgmc %>%
    mutate(ALS_risk_score = as.numeric(pred_pgmc))
  
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

# =============================
# Heatmap of protein signature
# =============================
run_heatmap_signature <- function(
    results_ALL,
    fluid,
    proteins,
    risk_scores_df,   
    group_colors,
    output_prefix = "heatmap",
    highlight_ids = NULL
) {
  
  message("Running heatmap for ", fluid)
  
  # heatmap data
  data_heatmap <- results_ALL[[fluid]]$data_adjusted %>%
    filter(type %in% c("PGMC","ALS","CTR"),
           Target %in% proteins) %>%
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
    distinct() %>%
    select(Target, NPQ_adj, ParticipantCode, type, sex, age,
           ALSFRS_1, SiteDiseaseOnset, `Disease duration`,SampleName) %>%
    pivot_wider(names_from = Target, values_from = NPQ_adj) %>%
    left_join(risk_scores_df %>%
                select(SampleName, ALS_risk_score),
              by = "SampleName")
  
  # heatmap matrix
  mat <- data_heatmap %>%
    select(-c(ParticipantCode, type, ALS_risk_score,
              sex, age, ALSFRS_1,
              SiteDiseaseOnset, `Disease duration`,SampleName)) %>%
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
# Perform ML models in all fluids for all comparisons

tissues = c("PLASMA", "CSF", "SERUM")

#### -> Lasso modeling
for (tissue in tissues) {
  
  # Prepara data for ML
  # -> PGMC vs CTR
  protein_data_PGMCvsCTR = results_ALL[[tissue]]$data_adjusted %>%
    filter(type %in% c("PGMC","CTR")) %>%
    select(SampleName,Target,NPQ_adj,type) %>%
    pivot_wider(names_from = Target,
                values_from = NPQ_adj)
  
  protein_data_PGMCvsCTR_new = protein_data_PGMCvsCTR %>%
    select(-SampleName) %>%
    rename(status = type) %>%
    mutate(status = ifelse(status == "PGMC",1,0))
  
  # -> ALS vs CTR
  protein_data_ALSvsCTR = results_ALL[[tissue]]$data_adjusted %>%
    filter(type %in% c("ALS","CTR")) %>%
    select(SampleName,Target,NPQ_adj,type) %>%
    pivot_wider(names_from = Target,
                values_from = NPQ_adj)
  
  protein_data_ALSvsCTR_new = protein_data_ALSvsCTR %>%
    select(-SampleName) %>%
    rename(status = type) %>%
    mutate(status = ifelse(status == "ALS",1,0))
  
  # -> ALS vs PGMC
  protein_data_ALSvsPGMC = results_ALL[[tissue]]$data_adjusted %>%
    filter(type %in% c("ALS","PGMC")) %>%
    select(SampleName,Target,NPQ_adj,type) %>%
    pivot_wider(names_from = Target,
                values_from = NPQ_adj)
  
  protein_data_ALSvsPGMC_new = protein_data_ALSvsPGMC %>%
    select(-SampleName) %>%
    rename(status = type) %>%
    mutate(status = ifelse(status == "ALS",1,0))
  
  # Run Lasso and ROC curve with 5-fold cv and 500 bootstrap iterations
  lm_PGMC_CTR <- runML_with_lasso(protein_data_PGMCvsCTR_new,bs_count = 500)
  final_roc_plot_PGMC_CTR = calculateROC(lm_PGMC_CTR,
                                         paste0("plots/ML/Lasso/ROC_lm_",tissue,"_PGMC_CTR.pdf"))
  lm_ALS_CTR <- runML_with_lasso(protein_data_ALSvsCTR_new,bs_count = 500)
  final_roc_plot_ALS_CTR = calculateROC(lm_ALS_CTR,
                                         paste0("plots/ML/Lasso/ROC_lm_",tissue,"_ALS_CTR.pdf"))
  lm_ALS_PGMC <- runML_with_lasso(protein_data_ALSvsPGMC_new,bs_count = 500)
  final_roc_plot_ALS_PGMC = calculateROC(lm_ALS_PGMC,
                                        paste0("plots/ML/Lasso/ROC_lm_",tissue,"_ALS_PGMC.pdf"))
  
  # extract weights + plot averaged
  lm_weights_PGMC_CTR = feature_importance(lm_PGMC_CTR,
                                                plot_path = paste0("plots/ML/Lasso/protein_selection_",tissue,"_PGMC_CTR.pdf"))
  write_xlsx(lm_weights_PGMC_CTR, path = paste0("results/weights_lm_",tissue,"_PGMC_CTR.xlsx"))
  lm_weights_ALS_CTR = feature_importance(lm_ALS_CTR,
                                                plot_path = paste0("plots/ML/Lasso/protein_selection_",tissue,"_ALS_CTR.pdf"))
  write_xlsx(lm_weights_ALS_CTR, path = paste0("results/weights_lm_",tissue,"_ALS_CTR.xlsx"))
  lm_weights_ALS_PGMC = feature_importance(lm_ALS_PGMC,
                                               plot_path = paste0("plots/ML/Lasso/protein_selection_",tissue,"_ALS_PGMC.pdf"))
  write_xlsx(lm_weights_ALS_PGMC, path = paste0("results/weights_lm_",tissue,"_ALS_PGMC.xlsx"))
}

# for high detectable proteins
# for (tissue in tissues) {
#   
#   high_detected_proteins = detectability_summary %>% 
#     filter(SampleMatrixType == tissue) %>%
#     filter(detectability == "high") %>%
#     pull(Target)
#   
#   # Prepara data for ML
#   # -> PGMC vs CTR
#   protein_data_PGMCvsCTR = results_ALL[[tissue]]$data_adjusted %>%
#     filter(type %in% c("PGMC","CTR")) %>%
#     filter(Target %in% high_detected_proteins) %>%
#     select(SampleName,Target,NPQ_adj,type) %>%
#     pivot_wider(names_from = Target,
#                 values_from = NPQ_adj)
#   
#   protein_data_PGMCvsCTR_new = protein_data_PGMCvsCTR %>%
#     select(-SampleName) %>%
#     rename(status = type) %>%
#     mutate(status = ifelse(status == "PGMC",1,0))
#   
#   # -> ALS vs CTR
#   protein_data_ALSvsCTR = results_ALL[[tissue]]$data_adjusted %>%
#     filter(type %in% c("ALS","CTR")) %>%
#     filter(Target %in% high_detected_proteins) %>%
#     select(SampleName,Target,NPQ_adj,type) %>%
#     pivot_wider(names_from = Target,
#                 values_from = NPQ_adj)
#   
#   protein_data_ALSvsCTR_new = protein_data_ALSvsCTR %>%
#     select(-SampleName) %>%
#     rename(status = type) %>%
#     mutate(status = ifelse(status == "ALS",1,0))
#   
#   # -> ALS vs PGMC
#   protein_data_ALSvsPGMC = results_ALL[[tissue]]$data_adjusted %>%
#     filter(type %in% c("ALS","PGMC")) %>%
#     filter(Target %in% high_detected_proteins) %>%
#     select(SampleName,Target,NPQ_adj,type) %>%
#     pivot_wider(names_from = Target,
#                 values_from = NPQ_adj)
#   
#   protein_data_ALSvsPGMC_new = protein_data_ALSvsPGMC %>%
#     select(-SampleName) %>%
#     rename(status = type) %>%
#     mutate(status = ifelse(status == "ALS",1,0))
#   
#   # Run Lasso and ROC curve with 5-fold cv and 500 bootstrap iterations
#   lm_PGMC_CTR <- runML_with_lasso(protein_data_PGMCvsCTR_new,bs_count = 500)
#   final_roc_plot_PGMC_CTR = calculateROC(lm_PGMC_CTR,
#                                          paste0("plots/ML/Lasso/high_detectable/ROC_lm_",tissue,"_PGMC_CTR_high_detected.pdf"))
#   lm_ALS_CTR <- runML_with_lasso(protein_data_ALSvsCTR_new,bs_count = 500)
#   final_roc_plot_ALS_CTR = calculateROC(lm_ALS_CTR,
#                                         paste0("plots/ML/Lasso/high_detectable/ROC_lm_",tissue,"_ALS_CTR_high_detected.pdf"))
#   lm_ALS_PGMC <- runML_with_lasso(protein_data_ALSvsPGMC_new,bs_count = 500)
#   final_roc_plot_ALS_PGMC = calculateROC(lm_ALS_PGMC,
#                                          paste0("plots/ML/Lasso/high_detectable/ROC_lm_",tissue,"_ALS_PGMC_high_detected.pdf"))
#   
#   # extract weights + plot averaged
#   lm_weights_PGMC_CTR = feature_importance(lm_PGMC_CTR,
#                                                 plot_path = paste0("plots/ML/Lasso/high_detectable/protein_selection_",
#                                                                    tissue,"_PGMC_CTR_high_detected.pdf"))
#   write_xlsx(lm_weights_PGMC_CTR, path = paste0("results/weights_lm_",
#                                                 tissue,"_PGMC_CTR_high_detected.xlsx"))
#   lm_weights_ALS_CTR = feature_importance(lm_ALS_CTR,
#                                                plot_path = paste0("plots/ML/Lasso/high_detectable/protein_selection_",
#                                                                   tissue,"_ALS_CTR_high_detected.pdf"))
#   write_xlsx(lm_weights_ALS_CTR, path = paste0("results/weights_lm_",
#                                                tissue,
#                                                "_ALS_CTR_high_detected.xlsx"))
#   lm_weights_ALS_PGMC = feature_importance(lm_ALS_PGMC,
#                                                 plot_path = paste0("plots/ML/Lasso/high_detectable/protein_selection_",tissue,"_ALS_PGMC_high_detected.pdf"))
#   write_xlsx(lm_weights_ALS_PGMC, path = paste0("results/weights_lm_",
#                                                 tissue,"_ALS_PGMC_high_detected.xlsx"))
# }

# for all proteins except NEFL
# for (tissue in tissues) {
#   
#   # Prepara data for ML
#   # -> PGMC vs CTR
#   protein_data_PGMCvsCTR = results_ALL[[tissue]]$data_adjusted %>%
#     filter(type %in% c("PGMC","CTR")) %>%
#     filter(Target != "NEFL") %>%
#     select(SampleName,Target,NPQ_adj,type) %>%
#     pivot_wider(names_from = Target,
#                 values_from = NPQ_adj)
#   
#   protein_data_PGMCvsCTR_new = protein_data_PGMCvsCTR %>%
#     select(-SampleName) %>%
#     rename(status = type) %>%
#     mutate(status = ifelse(status == "PGMC",1,0))
#   
#   # -> ALS vs CTR
#   protein_data_ALSvsCTR = results_ALL[[tissue]]$data_adjusted %>%
#     filter(type %in% c("ALS","CTR")) %>%
#     filter(Target != "NEFL") %>%
#     select(SampleName,Target,NPQ_adj,type) %>%
#     pivot_wider(names_from = Target,
#                 values_from = NPQ_adj)
#   
#   protein_data_ALSvsCTR_new = protein_data_ALSvsCTR %>%
#     select(-SampleName) %>%
#     rename(status = type) %>%
#     mutate(status = ifelse(status == "ALS",1,0))
#   
#   # -> ALS vs PGMC
#   protein_data_ALSvsPGMC = results_ALL[[tissue]]$data_adjusted %>%
#     filter(type %in% c("ALS","PGMC")) %>%
#     filter(Target != "NEFL") %>%
#     select(SampleName,Target,NPQ_adj,type) %>%
#     pivot_wider(names_from = Target,
#                 values_from = NPQ_adj)
#   
#   protein_data_ALSvsPGMC_new = protein_data_ALSvsPGMC %>%
#     select(-SampleName) %>%
#     rename(status = type) %>%
#     mutate(status = ifelse(status == "ALS",1,0))
#   
#   # Run Lasso and ROC curve with 5-fold cv and 500 bootstrap iterations
#   # lm_PGMC_CTR <- runML_with_lasso(protein_data_PGMCvsCTR_new,bs_count = 500)
#   # final_roc_plot_PGMC_CTR = calculateROC(lm_PGMC_CTR,
#   #                                        paste0("plots/ML/Lasso/without_NEFL/ROC_lm_",tissue,"_PGMC_CTR_noNEFL.pdf"))
#   lm_ALS_CTR <- runML_with_lasso(protein_data_ALSvsCTR_new,bs_count = 500)
#   final_roc_plot_ALS_CTR = calculateROC(lm_ALS_CTR,
#                                         paste0("plots/ML/Lasso/without_NEFL/ROC_lm_",tissue,"_ALS_CTR_noNEFL.pdf"))
#   # lm_ALS_PGMC <- runML_with_lasso(protein_data_ALSvsPGMC_new,bs_count = 500)
#   # final_roc_plot_ALS_PGMC = calculateROC(lm_ALS_PGMC,
#   #                                        paste0("plots/ML/Lasso/without_NEFL/ROC_lm_",tissue,"_ALS_PGMC_noNEFL.pdf"))
#   # 
#   # extract weights + plot averaged
#   # lm_weights_PGMC_CTR = feature_importance(lm_PGMC_CTR,
#   #                                               plot_path = paste0("plots/ML/Lasso/without_NEFL/protein_selection_",
#   #                                                                  tissue,"_PGMC_CTR_noNEFL.pdf"))
#   # write_xlsx(lm_weights_PGMC_CTR, path = paste0("results/weights_lm_",
#   #                                               tissue,"_PGMC_CTR_noNEFL.xlsx"))
#   lm_weights_ALS_CTR = feature_importance(lm_ALS_CTR,
#                                                plot_path = paste0("plots/ML/Lasso/without_NEFL/protein_selection_",
#                                                                   tissue,"_ALS_CTR_noNEFL.pdf"))
#   write_xlsx(lm_weights_ALS_CTR, path = paste0("results/weights_lm_",
#                                                tissue,
#                                                "_ALS_CTR_noNEFL.xlsx"))
#   # lm_weights_ALS_PGMC = feature_importance(lm_ALS_PGMC,
#   #                                               plot_path = paste0("plots/ML/Lasso/without_NEFL/protein_selection_",tissue,"_ALS_PGMC_noNEFL.pdf"))
#   # write_xlsx(lm_weights_ALS_PGMC, path = paste0("results/weights_lm_",
#   #                                               tissue,"_ALS_PGMC_noNEFL.xlsx"))
# }

#### -> Elastic Net modeling
for (tissue in tissues) {
  
  # Prepara data for ML
  # -> PGMC vs CTR
  protein_data_PGMCvsCTR = results_ALL[[tissue]]$data_adjusted %>%
    filter(type %in% c("PGMC","CTR")) %>%
    select(SampleName,Target,NPQ_adj,type) %>%
    pivot_wider(names_from = Target,
                values_from = NPQ_adj)
  
  protein_data_PGMCvsCTR_new = protein_data_PGMCvsCTR %>%
    select(-SampleName) %>%
    rename(status = type) %>%
    mutate(status = ifelse(status == "PGMC",1,0))
  
  # -> ALS vs CTR
  protein_data_ALSvsCTR = results_ALL[[tissue]]$data_adjusted %>%
    filter(type %in% c("ALS","CTR")) %>%
    select(SampleName,Target,NPQ_adj,type) %>%
    pivot_wider(names_from = Target,
                values_from = NPQ_adj)
  
  protein_data_ALSvsCTR_new = protein_data_ALSvsCTR %>%
    select(-SampleName) %>%
    rename(status = type) %>%
    mutate(status = ifelse(status == "ALS",1,0))
  
  # -> ALS vs PGMC
  protein_data_ALSvsPGMC = results_ALL[[tissue]]$data_adjusted %>%
    filter(type %in% c("ALS","PGMC")) %>%
    select(SampleName,Target,NPQ_adj,type) %>%
    pivot_wider(names_from = Target,
                values_from = NPQ_adj)
  
  protein_data_ALSvsPGMC_new = protein_data_ALSvsPGMC %>%
    select(-SampleName) %>%
    rename(status = type) %>%
    mutate(status = ifelse(status == "ALS",1,0))
  
  # Run Lasso and ROC curve with 5-fold cv and 500 bootstrap iterations
  lm_enet_PGMC_CTR <- runML_with_elasticnet(protein_data_PGMCvsCTR_new,bs_count = 500)
  final_roc_plot_PGMC_CTR_enet = calculateROC(lm_enet_PGMC_CTR,
                                         paste0("plots/ML/Elastic Net/ROC_lm_",tissue,"_PGMC_CTR.pdf"))
  lm_enet_ALS_CTR <- runML_with_elasticnet(protein_data_ALSvsCTR_new,bs_count = 500)
  final_roc_plot_ALS_CTR_enet = calculateROC(lm_enet_ALS_CTR,
                                        paste0("plots/ML/Elastic Net/ROC_lm_",tissue,"_ALS_CTR.pdf"))
  lm_enet_ALS_PGMC <- runML_with_elasticnet(protein_data_ALSvsPGMC_new,bs_count = 500)
  final_roc_plot_ALS_PGMC_enet = calculateROC(lm_enet_ALS_PGMC,
                                         paste0("plots/ML/Elastic Net/ROC_lm_",tissue,"_ALS_PGMC.pdf"))
  
  # extract weights + plot averaged
  lm_weights_PGMC_CTR_enet = feature_importance(lm_enet_PGMC_CTR,
                                                plot_path = paste0("plots/ML/Elastic Net/protein_selection_",tissue,"_PGMC_CTR.pdf"))
  write_xlsx(lm_weights_PGMC_CTR_enet, path = paste0("results/weights_lm_enet_",tissue,"_PGMC_CTR.xlsx"))
  lm_weights_ALS_CTR_enet = feature_importance(lm_enet_ALS_CTR,
                                               plot_path = paste0("plots/ML/Elastic Net/protein_selection_",tissue,"_ALS_CTR.pdf"))
  write_xlsx(lm_weights_ALS_CTR_enet, path = paste0("results/weights_lm_enet_",tissue,"_ALS_CTR.xlsx"))
  lm_weights_ALS_PGMC_enet = feature_importance(lm_enet_ALS_PGMC,
                                                plot_path = paste0("plots/ML/Elastic Net/protein_selection_",tissue,"_ALS_PGMC.pdf"))
  write_xlsx(lm_weights_ALS_PGMC_enet, path = paste0("results/weights_lm_enet_",tissue,"_ALS_PGMC.xlsx"))
}


# =============================================================
# Find the most optimal set of proteins for each ML and fluid

## -> Lasso + Serum
proteins_serum_ALS_CTR = c("NEFL","GDNF","pTau-181","TAFA5","VEGFD","NGF","BASP1","IL6","SOD1")
lasso_optimal_protein_signature_serum = find_optimal_signature(results_ALL = results_ALL,
                                                         ranked_proteins = proteins_serum_ALS_CTR,
                                                         fluid = "SERUM",
                                                         output_prefix = "plots/ML/Lasso/performance_")

# Elastic Net + Serum
proteins_serum_ALS_CTR = c("NEFL","GDNF","pTau-181","NEFH","VEGFD","TAFA5",
                           "pTau-231","NGF","FABP3","CALB2")
enet_optimal_protein_signature_serum = find_optimal_signature(results_ALL = results_ALL,
                                                               ranked_proteins = proteins_serum_ALS_CTR,
                                                               fluid = "SERUM",
                                                               output_prefix = "plots/ML/Elastic Net/performance_")

# Elastic Net + Plasma
proteins_plasma_ALS_CTR = c("NEFL","NEFH","pTau-181","GDNF","FABP3","IL16","TEK")
enet_optimal_protein_signature_plasma = find_optimal_signature(results_ALL = results_ALL,
                                                              ranked_proteins = proteins_plasma_ALS_CTR,
                                                              fluid = "PLASMA",
                                                              output_prefix = "plots/ML/Elastic Net/performance_")


# Elastic Net + CSF
proteins_CSF_ALS_CTR =  c("NEFL","NEFH","CHIT1","IL6","MSLN","TNF","CHI3L1",
                          "IL12p70","UCHL1","CCL2","CCL3")

enet_optimal_protein_signature_CSF = find_optimal_signature(results_ALL = results_ALL,
                                                               ranked_proteins = proteins_CSF_ALS_CTR,
                                                               fluid = "CSF",
                                                               output_prefix = "plots/ML/Elastic Net/performance_")


# ====================================================
# Compute an ALS risk score based on ALS vs CTR model 
## -> Lasso + Serum
proteins_serum_ALS_CTR = lasso_optimal_protein_signature_serum$optimal_proteins

lasso_PGMC_serum_signature = run_ALS_signature_workflow(results_ALL,
                                                        proteins_serum_ALS_CTR,
                                                        fluid = "SERUM",
                                                        model_type = "lasso",
                                                        output_prefix = "plots/ML/Lasso/")

# ================================================================================
# Unsupervised visualisation (PCA) of PGMC, ALS, CTR based on 6-protein signature
protein_data_PCA_serum = protein_data_clean %>%
  filter(Target %in% lasso_optimal_protein_signature_serum$optimal_proteins) %>%
  filter(type %in% c("PGMC","ALS","CTR"))

pca_results_serum_adj <- run_pca_adjusted(remove_effect_covariates(protein_data_PCA_serum,keep = "type"), 
                                          "SERUM")

plots_subtype_adj <- plot_pca(pca_results_serum_adj,"SERUM based on 6-protein signature")

pdf("plots/ML/Lasso/PCA_SERUM_protein_signature.pdf", width = 8, height = 6.5) 
plots_subtype_adj
dev.off()

pca_results_serum_adj$scores <- pca_results_serum_adj$scores %>%
  mutate(ParticipantCode = ifelse(type == "PGMC",ParticipantCode,NA))

plots_subtype_adj_label <- plot_pca(pca_results_serum_adj,"SERUM based on 6-protein signature",
                                                       label = TRUE)

pdf("plots/ML/Lasso/PCA_SERUM_protein_signature_label.pdf", width = 8, height = 6.5) 
plots_subtype_adj_label
dev.off()


# ================================================================================
# Unsupervised visualisation (heatmap) of PGMC, ALS, CTR based on 5-protein signature

group_colors <- c(
  "CTR"  = "#6F8EB2",
  "ALS"  = "#B2936F",
  "PGMC" = "#ad5291")

run_heatmap_signature(results_ALL,
                      "SERUM",
                      lasso_optimal_protein_signature_serum$optimal_proteins,
                      lasso_PGMC_serum_signature$results,
                      group_colors = group_colors,
                      output_prefix = "plots/ML/Lasso/heatmap_SERUM",
                      highlight_ids = c("DE101","TR120","TR112"))

##### ----
# Elastic Net

## -> Elastic Net + Serum
EN_PGMC_serum_signature = run_ALS_signature_workflow(results_ALL,
                                                     enet_optimal_protein_signature_serum$optimal_proteins,
                                                        fluid = "SERUM",
                                                        model_type = "elastic_net",
                                                        output_prefix = "plots/ML/Elastic Net/")


# ================================================================================
# Unsupervised visualisation (PCA) of PGMC, ALS, CTR based on 9-protein signature
protein_data_PCA_serum = protein_data_clean %>%
  filter(Target %in% enet_optimal_protein_signature_serum$optimal_proteins) %>%
  filter(type %in% c("PGMC","ALS","CTR"))

pca_results_serum_adj <- run_pca_adjusted(remove_effect_covariates(protein_data_PCA_serum,keep = "type"), 
                                          "SERUM")

plots_subtype_adj <- plot_pca(pca_results_serum_adj,"SERUM based on 9-protein signature")

pdf("plots/ML/Elastic Net/PCA_SERUM_protein_signature.pdf", width = 8, height = 6.5) 
plots_subtype_adj
dev.off()

pca_results_serum_adj$scores <- pca_results_serum_adj$scores %>%
  mutate(ParticipantCode = ifelse(type == "PGMC",ParticipantCode,NA))

plots_subtype_adj_label <- plot_pca(pca_results_serum_adj,"SERUM based on 9-protein signature",
                                    label = TRUE)

pdf("plots/ML/Elastic Net/PCA_SERUM_protein_signature_label.pdf", width = 8, height = 6.5) 
plots_subtype_adj_label
dev.off()


# ================================================================================
# Unsupervised visualisation (heatmap) of PGMC, ALS, CTR based on 9-protein signature

run_heatmap_signature(results_ALL,
                      "SERUM",
                      enet_optimal_protein_signature_serum$optimal_proteins,
                      EN_PGMC_serum_signature$results,
                      group_colors = group_colors,
                      output_prefix = "plots/ML/Elastic Net/heatmap_SERUM",
                      highlight_ids = c("DE101","TR120","TR114","TR112"))

## -> Elastic Net + Plasma
EN_PGMC_plasma_signature = run_ALS_signature_workflow(results_ALL,
                                                     enet_optimal_protein_signature_plasma$optimal_proteins,
                                                     fluid = "PLASMA",
                                                     model_type = "elastic_net",
                                                     output_prefix = "plots/ML/Elastic Net/")


# ================================================================================
# Unsupervised visualisation (PCA) of PGMC, ALS, CTR based on 3-protein signature
protein_data_PCA_PLASMA = protein_data_clean %>%
  filter(Target %in% enet_optimal_protein_signature_plasma$optimal_proteins) %>%
  filter(type %in% c("PGMC","ALS","CTR"))

pca_results_PLASMA_adj <- run_pca_adjusted(remove_effect_covariates(protein_data_PCA_PLASMA,keep = "type"), 
                                          "PLASMA")

plots_subtype_adj <- plot_pca(pca_results_PLASMA_adj,"PLASMA based on 3-protein signature")

pdf("plots/ML/Elastic Net/PCA_PLASMA_protein_signature.pdf", width = 8, height = 6.5) 
plots_subtype_adj
dev.off()

pca_results_PLASMA_adj$scores <- pca_results_PLASMA_adj$scores %>%
  mutate(ParticipantCode = ifelse(type == "PGMC",ParticipantCode,NA))

plots_subtype_adj_label <- plot_pca(pca_results_PLASMA_adj,"PLASMA based on 3-protein signature",
                                    label = TRUE)

pdf("plots/ML/Elastic Net/PCA_PLASMA_protein_signature_label.pdf", width = 8, height = 6.5) 
plots_subtype_adj_label
dev.off()


# ================================================================================
# Unsupervised visualisation (heatmap) of PGMC, ALS, CTR based on 3-protein signature

run_heatmap_signature(results_ALL,
                      "PLASMA",
                      enet_optimal_protein_signature_plasma$optimal_proteins,
                      EN_PGMC_plasma_signature$results,
                      group_colors = group_colors,
                      output_prefix = "plots/ML/Elastic Net/heatmap_PLASMA",
                      highlight_ids = c("DE101","TR120","DE102","TR113","TR128"))

## -> Elastic Net + CSF
EN_PGMC_CSF_signature = run_ALS_signature_workflow(results_ALL,
                                                      enet_optimal_protein_signature_CSF$optimal_proteins,
                                                      fluid = "CSF",
                                                      model_type = "elastic_net",
                                                      output_prefix = "plots/ML/Elastic Net/")


# ================================================================================
# Unsupervised visualisation (PCA) of PGMC, ALS, CTR based on 2-protein signature
protein_data_PCA_CSF = protein_data_clean %>%
  filter(Target %in% enet_optimal_protein_signature_CSF$optimal_proteins) %>%
  filter(type %in% c("PGMC","ALS","CTR"))

pca_results_CSF_adj <- run_pca_adjusted(remove_effect_covariates(protein_data_PCA_CSF,keep = "type"), 
                                           "CSF")

plots_subtype_adj <- plot_pca(pca_results_CSF_adj,"CSF based on 2-protein signature")

pdf("plots/ML/Elastic Net/PCA_CSF_protein_signature.pdf", width = 8, height = 6.5) 
plots_subtype_adj
dev.off()

pca_results_CSF_adj$scores <- pca_results_CSF_adj$scores %>%
  mutate(ParticipantCode = ifelse(type == "PGMC",ParticipantCode,NA))

plots_subtype_adj_label <- plot_pca(pca_results_CSF_adj,"CSF based on 2-protein signature",
                                    label = TRUE)

pdf("plots/ML/Elastic Net/PCA_CSF_protein_signature_label.pdf", width = 8, height = 6.5) 
plots_subtype_adj_label
dev.off()


# ================================================================================
# Unsupervised visualisation (heatmap) of PGMC, ALS, CTR based on 9-protein signature

run_heatmap_signature(results_ALL,
                      "CSF",
                      enet_optimal_protein_signature_CSF$optimal_proteins,
                      EN_PGMC_CSF_signature$results,
                      group_colors = group_colors,
                      output_prefix = "plots/ML/Elastic Net/heatmap_CSF",
                      highlight_ids = c("DE101"))



