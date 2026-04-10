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

# Extracts standardized ROC coordinates for plotting
get_roc_points <- function(y_true, y_score, run_id) {
  roc_obj <- pROC::roc(y_true, y_score, direction = "<", quiet = TRUE)
  coords <- pROC::coords(roc_obj, x = seq(0, 1, length.out = 2000), 
                         input = "specificity", ret = c("specificity", "sensitivity"))
  
  data.frame(
    FPR = 1 - coords$specificity,
    TPR = coords$sensitivity,
    run = run_id,
    auc = as.numeric(pROC::auc(roc_obj))
  )
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
            'unknown algorithm. please select \'lm\',\'rf\',\'svm rad\' or \'svm lin\' ' = algorithm %in% c('lm','rf','svm rad','svm lin'))
  
  key_output = c('linear regression (lasso)','random forest', 'Support Vector Machine (radial kernel)', 'Support Vector Machine (linear kernel)')
  names(key_output) = c('lm','rf','svm rad','svm lin')
  # Output parameters for user check
  print(paste0('Running a ',BS_number,'x boot strap with a ',cv,' fold cross validation. Algorithm is ', key_output[[algorithm]]))
  
  data_frame$status = as.factor(make.names(data_frame$status)) # caret needs prediction variable to have name; turns 0 -> X0 and 1 -> X1
  
  # create lists to save results from bs
  models = list()
  importance = list()
  predictions = list()
  indices = list()
  predictions_raw = list()
  actuals_list = list()
  smp_size = floor(0.8 * nrow(data_frame)) # use 80% of data for training, 20% for testing
  
  #set the right parameters for each algo
  if(algorithm == 'lm'| algorithm == 'svm rad' | algorithm == 'svm lin'){
    ctrl = trainControl(method="cv",   
                        number = cv,        
                        summaryFunction=twoClassSummary,   # Use AUC to pick the best model
                        classProbs=TRUE,savePredictions = TRUE)
  }
  else if(algorithm == 'rf'){
    
    n_features = ncol(data_frame[ ,!names(data_frame) == 'status'])
    ctrl = trainControl(method = "cv", 
                        number = cv, 
                        search = 'grid',classProbs = TRUE, savePredictions = TRUE, summaryFunction=twoClassSummary )
  }
  
  # run machine learning
  for(i in 1:BS_number){
    
    train_ind = createDataPartition(data_frame$status, p = 0.8, list = FALSE)
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

# =============================
# Make ROC curve
# =============================
calculateROC <- function(list_from_ML, plot_path = FALSE) {
  
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

# =============================
# Plot feature importance
# =============================
analyze_lasso_stability = function(lasso_analysis_object, plot_path = "plots/protein_stability_selection.pdf", number = 15) {
  
  if (!dir.exists(dirname(plot_path))) dir.create(dirname(plot_path), recursive = TRUE)
  
  imp_list <- lasso_analysis_object$importance
  num_runs <- length(imp_list)
  all_proteins <- unique(unlist(lapply(imp_list, function(x) rownames(x$importance))))
  
  imp_matrix <- matrix(0, nrow = length(all_proteins), ncol = num_runs)
  rownames(imp_matrix) <- all_proteins
  
  for(i in 1:num_runs) {
    current_df <- as.data.frame(imp_list[[i]]$importance)
    vals <- current_df$Overall
    names(vals) <- rownames(current_df)
    common_names <- intersect(names(vals), all_proteins)
    imp_matrix[common_names, i] <- vals[common_names]
  }
  
  weight <- data.frame(
    avg_coef = rowMeans(imp_matrix,na.rm = T),
    sd_coef = apply(imp_matrix, 1, function(x) sd(x,na.rm = T)),
    sem_coef = apply(imp_matrix, 1, function(x) sd(x,na.rm = T)/sqrt(length(x))),
    freq = (rowSums(imp_matrix > 0,na.rm = T) / num_runs) * 100)
  
  # Sort by Selection Frequency
  weight <- weight[order(-weight$freq), ]
  number <- min(number, nrow(weight))
  plot_df <- weight[1:number, ]
  plot_df$protein <- rownames(plot_df)
  
  if(plot_path != FALSE){
    pdf(plot_path, paper="a4r", width = 11, height = 8)
    
    p <- ggplot(plot_df, 
                aes(x = reorder(protein, freq), 
                    y = freq, 
                    fill = avg_coef)) +
      geom_bar(stat = "identity", width = 0.7, color = "white") +
      geom_text(aes(label = paste0(round(freq, 1), "%")), 
                hjust = -0.2, size = 3.5, color = "grey30") +
      scale_fill_gradient(low = "#deebf7", high = "#084594") +
      labs(x = "Protein",
           y = "Selection Frequency (%)",
           title = "Protein frequency across 500 bootstraps",
           fill = "Mean importance score") +
      theme_classic() +
      coord_flip() + 
      scale_y_continuous(limits = c(0, 110), breaks = seq(0, 100, 20)) +
      theme(
        axis.title = element_text(face = "bold"),
        title = element_text(size = 14)
      )
    
    print(p)
    dev.off()
  }
  
  message("Plot saved to ", plot_path)
  return(weight)
}

tissues = c("PLASMA", "CSF", "SERUM")

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
                                         paste0("plots/ML/ROC_lm_",tissue,"_PGMC_CTR.pdf"))
  lm_ALS_CTR <- runML_with_lasso(protein_data_ALSvsCTR_new,bs_count = 500)
  final_roc_plot_ALS_CTR = calculateROC(lm_ALS_CTR,
                                        paste0("plots/ML/ROC_lm_",tissue,"_ALS_CTR.pdf"))
  lm_ALS_PGMC <- runML_with_lasso(protein_data_ALSvsPGMC_new,bs_count = 500)
  final_roc_plot_ALS_PGMC = calculateROC(lm_ALS_PGMC,
                                         paste0("plots/ML/ROC_lm_",tissue,"_ALS_PGMC.pdf"))
  
  # extract weights + plot averaged
  lm_weights_PGMC_CTR = analyze_lasso_stability(lm_PGMC_CTR,
                                                plot_path = paste0("plots/ML/protein_selection_",tissue,"_PGMC_CTR.pdf"))
  write_xlsx(lm_weights_PGMC_CTR, path = paste0("results/weights_lm_",tissue,"_PGMC_CTR.xlsx"))
  lm_weights_ALS_CTR = analyze_lasso_stability(lm_ALS_CTR,
                                               plot_path = paste0("plots/ML/protein_selection_",tissue,"_ALS_CTR.pdf"))
  write_xlsx(lm_weights_ALS_CTR, path = paste0("results/weights_lm_",tissue,"_ALS_CTR.xlsx"))
  lm_weights_ALS_PGMC = analyze_lasso_stability(lm_ALS_PGMC,
                                                plot_path = paste0("plots/ML/protein_selection_",tissue,"_ALS_PGMC.pdf"))
  write_xlsx(lm_weights_ALS_PGMC, path = paste0("results/weights_lm_",tissue,"_ALS_PGMC.xlsx"))
}