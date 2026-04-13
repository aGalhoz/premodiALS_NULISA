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
analyze_lasso_stability = function(lasso_analysis_object, 
                                   plot_path = "plots/protein_stability_selection.pdf", 
                                   number = 15) {
  
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

# ------------------------------------------------------------------------------
# ====================================================
# Perform ML models in all fluids for all comparisons

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

# for high detectable proteins
for (tissue in tissues) {
  
  high_detected_proteins = detectability_summary %>% 
    filter(SampleMatrixType == tissue) %>%
    filter(detectability == "high") %>%
    pull(Target)
  
  # Prepara data for ML
  # -> PGMC vs CTR
  protein_data_PGMCvsCTR = results_ALL[[tissue]]$data_adjusted %>%
    filter(type %in% c("PGMC","CTR")) %>%
    filter(Target %in% high_detected_proteins) %>%
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
    filter(Target %in% high_detected_proteins) %>%
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
    filter(Target %in% high_detected_proteins) %>%
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
                                         paste0("plots/ML/high_detectable/ROC_lm_",tissue,"_PGMC_CTR_high_detected.pdf"))
  lm_ALS_CTR <- runML_with_lasso(protein_data_ALSvsCTR_new,bs_count = 500)
  final_roc_plot_ALS_CTR = calculateROC(lm_ALS_CTR,
                                        paste0("plots/ML/high_detectable/ROC_lm_",tissue,"_ALS_CTR_high_detected.pdf"))
  lm_ALS_PGMC <- runML_with_lasso(protein_data_ALSvsPGMC_new,bs_count = 500)
  final_roc_plot_ALS_PGMC = calculateROC(lm_ALS_PGMC,
                                         paste0("plots/ML/high_detectable/ROC_lm_",tissue,"_ALS_PGMC_high_detected.pdf"))
  
  # extract weights + plot averaged
  lm_weights_PGMC_CTR = analyze_lasso_stability(lm_PGMC_CTR,
                                                plot_path = paste0("plots/ML/high_detectable/protein_selection_",
                                                                   tissue,"_PGMC_CTR_high_detected.pdf"))
  write_xlsx(lm_weights_PGMC_CTR, path = paste0("results/weights_lm_",
                                                tissue,"_PGMC_CTR_high_detected.xlsx"))
  lm_weights_ALS_CTR = analyze_lasso_stability(lm_ALS_CTR,
                                               plot_path = paste0("plots/ML/high_detectable/protein_selection_",
                                                                  tissue,"_ALS_CTR_high_detected.pdf"))
  write_xlsx(lm_weights_ALS_CTR, path = paste0("results/weights_lm_",
                                               tissue,
                                               "_ALS_CTR_high_detected.xlsx"))
  lm_weights_ALS_PGMC = analyze_lasso_stability(lm_ALS_PGMC,
                                                plot_path = paste0("plots/ML/high_detectable/protein_selection_",tissue,"_ALS_PGMC_high_detected.pdf"))
  write_xlsx(lm_weights_ALS_PGMC, path = paste0("results/weights_lm_",
                                                tissue,"_ALS_PGMC_high_detected.xlsx"))
}

# for all proteins except NEFL
for (tissue in tissues) {
  
  # Prepara data for ML
  # -> PGMC vs CTR
  protein_data_PGMCvsCTR = results_ALL[[tissue]]$data_adjusted %>%
    filter(type %in% c("PGMC","CTR")) %>%
    filter(Target != "NEFL") %>%
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
    filter(Target != "NEFL") %>%
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
    filter(Target != "NEFL") %>%
    select(SampleName,Target,NPQ_adj,type) %>%
    pivot_wider(names_from = Target,
                values_from = NPQ_adj)
  
  protein_data_ALSvsPGMC_new = protein_data_ALSvsPGMC %>%
    select(-SampleName) %>%
    rename(status = type) %>%
    mutate(status = ifelse(status == "ALS",1,0))
  
  # Run Lasso and ROC curve with 5-fold cv and 500 bootstrap iterations
  # lm_PGMC_CTR <- runML_with_lasso(protein_data_PGMCvsCTR_new,bs_count = 500)
  # final_roc_plot_PGMC_CTR = calculateROC(lm_PGMC_CTR,
  #                                        paste0("plots/ML/without_NEFL/ROC_lm_",tissue,"_PGMC_CTR_noNEFL.pdf"))
  lm_ALS_CTR <- runML_with_lasso(protein_data_ALSvsCTR_new,bs_count = 500)
  final_roc_plot_ALS_CTR = calculateROC(lm_ALS_CTR,
                                        paste0("plots/ML/without_NEFL/ROC_lm_",tissue,"_ALS_CTR_noNEFL.pdf"))
  # lm_ALS_PGMC <- runML_with_lasso(protein_data_ALSvsPGMC_new,bs_count = 500)
  # final_roc_plot_ALS_PGMC = calculateROC(lm_ALS_PGMC,
  #                                        paste0("plots/ML/without_NEFL/ROC_lm_",tissue,"_ALS_PGMC_noNEFL.pdf"))
  # 
  # extract weights + plot averaged
  # lm_weights_PGMC_CTR = analyze_lasso_stability(lm_PGMC_CTR,
  #                                               plot_path = paste0("plots/ML/without_NEFL/protein_selection_",
  #                                                                  tissue,"_PGMC_CTR_noNEFL.pdf"))
  # write_xlsx(lm_weights_PGMC_CTR, path = paste0("results/weights_lm_",
  #                                               tissue,"_PGMC_CTR_noNEFL.xlsx"))
  lm_weights_ALS_CTR = analyze_lasso_stability(lm_ALS_CTR,
                                               plot_path = paste0("plots/ML/without_NEFL/protein_selection_",
                                                                  tissue,"_ALS_CTR_noNEFL.pdf"))
  write_xlsx(lm_weights_ALS_CTR, path = paste0("results/weights_lm_",
                                               tissue,
                                               "_ALS_CTR_noNEFL.xlsx"))
  # lm_weights_ALS_PGMC = analyze_lasso_stability(lm_ALS_PGMC,
  #                                               plot_path = paste0("plots/ML/without_NEFL/protein_selection_",tissue,"_ALS_PGMC_noNEFL.pdf"))
  # write_xlsx(lm_weights_ALS_PGMC, path = paste0("results/weights_lm_",
  #                                               tissue,"_ALS_PGMC_noNEFL.xlsx"))
}

# ====================================================
# Compute an ALS risk score based on ALS vs CTR model

proteins_serum_ALS_CTR = c("NEFL","GDNF","pTau-181","TAFA5","VEGFD")

# get data of ALS vs CTR
protein_data_ALSvsCTR_serum = results_ALL[["SERUM"]]$data_adjusted %>%
  filter(type %in% c("ALS","CTR")) %>%
  filter(Target %in% proteins_serum_ALS_CTR) %>%
  select(SampleName,Target,NPQ_adj,type) %>%
  pivot_wider(names_from = Target,
              values_from = NPQ_adj)
protein_data_ALSvsCTR_serum_new = protein_data_ALSvsCTR_serum %>%
  select(-SampleName) %>%
  rename(status = type) %>%
  mutate(status = ifelse(status == "ALS",1,0))

# prediction of ALS vs CTR using protein signature
y_pred = protein_data_ALSvsCTR_serum_new$status 
X_serum = as.matrix(protein_data_ALSvsCTR_serum_new[,proteins_serum_ALS_CTR])

# scaling
X_scaled <- scale(X_serum)
scaling_center <- attr(X_scaled, "scaled:center")
scaling_scale  <- attr(X_scaled, "scaled:scale")

# model 
set.seed(123)
cv_model <- cv.glmnet(X_scaled, y_pred,family = "binomial",alpha = 1)
final_model <- cv_model$glmnet.fit
lambda_opt <- cv_model$lambda.min

# model coeficients
coef_vec <- coef(cv_model, s = "lambda.min")
coef_df <- as.data.frame(as.matrix(coef_vec))
coef_df$feature <- rownames(coef_df)
print(coef_df)

# get PGMC only data
protein_data_PGMC_serum = results_ALL[["SERUM"]]$data_adjusted %>%
  filter(type == "PGMC") %>%
  filter(Target %in% proteins_serum_ALS_CTR) %>%
  select(Target,NPQ_adj,SampleName) %>%
  pivot_wider(names_from = Target,
              values_from = NPQ_adj) 

PGMC_samples = protein_data_PGMC_serum$SampleName
X_PGMC = protein_data_PGMC_serum %>%
  select(all_of(proteins_serum_ALS_CTR))

# same scaling as before on the protein data
X_PGMC = as.matrix(X_PGMC)
X_PGMC_scaled = scale(X_PGMC,center = scaling_center,scale = scaling_scale)

# Risk scores of PGMC patients
ALS_risk_score_PGMC = predict(cv_model,
                              newx = X_PGMC_scaled,
                              s = "lambda.min",
                              type = "link") 

PGMC_results = protein_data_PGMC_serum %>%
  mutate(ALS_risk_score = as.numeric(ALS_risk_score_PGMC)) 

# Risk scores of ALS and CTR 
ALS_risk_scores_ALS_CTR = predict(cv_model,
                                  newx = X_scaled,
                                  s = "lambda.min",
                                  type = "link") 

ALS_CTR_results = protein_data_ALSvsCTR_serum %>%
  mutate(ALS_risk_score = as.numeric(ALS_risk_scores_ALS_CTR)) %>%
  select(-type)

# Risk scores of all and visualisation
ALS_risk_scores_all = rbind(PGMC_results,
                            ALS_CTR_results) %>%
  left_join(samples_ID_type %>% rename(SampleName = `Sample ID`))

plot_data_risk_scores = ALS_risk_scores_all %>%
  mutate(label = ifelse(type == "PGMC" & ALS_risk_score > 0,ParticipantCode,NA))
plot_data_risk_scores$type <- factor(plot_data_risk_scores$type, 
                                     levels = c("ALS", "PGMC", "CTR"))
group_colors <- c(
  "CTR"  = "#6F8EB2",
  "ALS"  = "#B2936F",
  "PGMC" = "#ad5291")

pdf("plots/ML/ALS_risk_score_serum.pdf",height = 5,width=7)
ggplot(plot_data_risk_scores, aes(x = ALS_risk_score, y = type, fill = type)) +
  geom_violin(trim = FALSE, alpha = 0.4) +
  #geom_boxplot(width = 0.15, outlier.shape = NA, size = 0.3) +
  geom_jitter(aes(color = type),width = 0,height = 0.08,size = 1.2,alpha = 0.8) +
  geom_text_repel(data = plot_data_risk_scores %>% filter(!is.na(label)),
    aes(label = label),
    size = 3,segment.color = "black",
    segment.size = 0.4,
    segment.alpha = 0.6,
    segment.curvature = 0.2,
    segment.angle = 20,
    box.padding = 0.7,
    point.padding = 0.4,
    max.overlaps = 50, direction = "y") +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black",size = 0.4) +
  theme_classic(base_size = 12) +
  theme(
    axis.line = element_line(size = 0.4),
    axis.ticks = element_line(size = 0.3),
    axis.text = element_text(color = "black"),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 13, face = "bold"),
    legend.position = "none"
  ) +
  labs(x = "ALS Risk Score",y = NULL,
    title = "ALS risk score based on 5-protein signature in Serum")
dev.off()


# ================================================================================
# Unsupervised visualisation (PCA) of PGMC, ALS, CTR based on 5-protein signature
protein_data_PCA_serum = protein_data_clean %>%
  filter(Target %in% proteins_serum_ALS_CTR) %>%
  filter(type %in% c("PGMC","ALS","CTR"))

pca_results_serum_adj <- run_pca_adjusted(remove_effect_covariates(protein_data_PCA_serum,keep = "type"), 
                                          "SERUM")

plots_subtype_adj <- plot_pca(pca_results_serum_adj,"SERUM based on 5-protein signature")

pdf("plots/ML/PCA_SERUM_protein_signature.pdf", width = 8, height = 6.5) 
plots_subtype_adj
dev.off()

pca_results_serum_adj$scores <- pca_results_serum_adj$scores %>%
  mutate(ParticipantCode = ifelse(type == "PGMC",ParticipantCode,NA))

plots_subtype_adj_label <- plot_pca(pca_results_serum_adj,"SERUM based on 5-protein signature",
                                                       label = TRUE)

pdf("plots/ML/PCA_SERUM_protein_signature_label.pdf", width = 8, height = 6.5) 
plots_subtype_adj_label
dev.off()


# ================================================================================
# Unsupervised visualisation (heatmap) of PGMC, ALS, CTR based on 5-protein signature

# data for heatmap
data_heatmap_serum_groups = results_ALL[["SERUM"]]$data_adjusted %>%
  filter(type %in% c("PGMC","ALS","CTR")) %>%
  filter(Target %in% proteins_serum_ALS_CTR) %>%
  select(Target,NPQ_adj,ParticipantCode,type) %>%
  pivot_wider(names_from = Target,
              values_from = NPQ_adj) %>%
  left_join(plot_data_risk_scores %>% select(ALS_risk_score,ParticipantCode))

heatmap_matrix_serum_groups = data_heatmap_serum_groups %>%
  select(-c(ParticipantCode,type,ALS_risk_score)) %>%
  as.matrix()

rownames(heatmap_matrix_serum_groups) = data_heatmap_serum_groups$ParticipantCode

# scale data
heatmap_matrix_serum_groups_scaled = t(scale(t(heatmap_matrix_serum_groups)))

# ALS score info
score_column_colors = colorRamp2(c(min(data_heatmap_serum_groups$ALS_risk_score),0,
                            max(data_heatmap_serum_groups$ALS_risk_score)),
                            c("lightgrey", "red","darkred"))

# annotation for participants
ha = rowAnnotation(Group = data_heatmap_serum_groups$type,
                   ALS_score = data_heatmap_serum_groups$ALS_risk_score,
                   IDs_interest = ifelse(data_heatmap_serum_groups$ParticipantCode %in% c("DE101", "TR120", "TR114"),"yes","no"),
                   col = list(Group = group_colors,
                              ALS_score = score_column_colors,
                              IDs_interest = c("yes" = "black","no" = "white")))

# plot heatmap
pdf("plots/ML/heatmap_protein_signature_clustered.pdf",height = 14)
Heatmap(
  heatmap_matrix_serum_groups_scaled,
  name = "Z-score",
  left_annotation  = ha,
  #row_split = data_heatmap_serum_groups$type,
  show_row_names = TRUE,
  show_column_names = TRUE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  col = colorRamp2(c(-2, 0, 2), c("#28d778", "white", "purple")))
dev.off()

pdf("plots/ML/heatmap_protein_signature_group_organized.pdf",height = 14)
Heatmap(
  heatmap_matrix_serum_groups_scaled,
  name = "Z-score",
  left_annotation  = ha,
  row_split = data_heatmap_serum_groups$type,
  show_row_names = TRUE,
  show_column_names = TRUE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  col = colorRamp2(c(-2, 0, 2), c("#28d778", "white", "purple")))
dev.off()


