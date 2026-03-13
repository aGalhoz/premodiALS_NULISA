source("00_initialization.R")

###############################################
### Helper Functions
###############################################

## 1. Assign group labels (ALS / CTR / PGMC / mimic)
assign_group_labels <- function(sample_df, ...) {
  bind_rows(...) %>%
    rename(ParticipantCode = `Optional Informarion (patient ID)`) %>%
    left_join(sample_df, by = "ParticipantCode")
}

## 2. Collect sample counts per fluid
count_samples_per_fluid <- function(df, fluid) {
  df %>%
    filter(SampleMatrixType == fluid, !is.na(type)) %>%
    distinct(SampleName, type) %>%
    group_by(type) %>%
    summarise(nr_samples = n(), .groups = "drop") %>%
    mutate(biofluid = fluid)
}

count_subtypes_per_fluid <- function(df, fluid) {
  df %>%
    filter(SampleMatrixType == fluid, !is.na(subtype)) %>%
    distinct(SampleName, subtype) %>%
    group_by(subtype) %>%
    summarise(nr_samples = n(), .groups = "drop") %>%
    mutate(biofluid = fluid)
}

## 3. Get mean per fluid
get_mean_per_fluid <- function(df, fluid) {
  df %>%
    filter(SampleMatrixType == fluid) %>%
    distinct(Target, UniProtID, NPQ) %>%
    group_by(Target, UniProtID) %>%
    summarise(NPQ = mean(NPQ), .groups = "drop") %>%
    rename_with(~ paste0("NPQ_", fluid), "NPQ")
}

## 4. Correlation plot generator
corr_plot <- function(data, x, y, title) {
  
  cor_test <- cor.test(data[[x]], data[[y]], use = "complete.obs")
  r_val <- round(cor_test$estimate, 2)
  p_val <- ifelse(cor_test$p.value < 0.001,
                  formatC(cor_test$p.value, format = "e", digits = 2),
                  round(cor_test$p.value, 3))
  
  cor_text <- paste0("italic(R) == ", r_val,
                     " * ',' ~ italic(p) == ", p_val)
  
  p <- ggscatter(data, x = x, y = y,
                 add = "reg.line",
                 add.params = list(color = "blue",fill = "lightblue"),
                 conf.int = TRUE) +
    annotate("text",
             x = min(data[[x]], na.rm = TRUE),
             y = max(data[[y]], na.rm = TRUE),
             hjust = 0,
             label = cor_text,
             parse = TRUE, size = 5) +
    labs(title = title, x = x, y = y) +
    theme( plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
           axis.title = element_text(size = 16), 
           axis.text = element_text(size = 14), 
           legend.title = element_text(size = 14), 
           legend.text = element_text(size = 14))
  
  ggMarginal(p, type = "histogram", color = "lightblue", fill = "lightblue")
}

make_corr <- function(df1, df2, label1, label2, filename) {
  merged <- inner_join(df1, df2, by = c("Target", "UniProtID"))
  p <- corr_plot(merged,
                 paste0("NPQ_", label1),
                 paste0("NPQ_", label2),
                 paste("Correlation", label1, "vs", label2))
  
  pdf(filename)
  print(p)
  dev.off()
  p
}

## 5. Detectability by fluid and plate
plot_detectability <- function(df, fluid, by_plate = FALSE) {
  
  df_f <- df %>% filter(SampleMatrixType == fluid)
  
  if (by_plate) {
    df_f <- df_f %>%
      group_by(PlateID, Target) %>%
      mutate(mean_value = mean(TargetDetectability_value, na.rm = TRUE),
             mean_pct   = mean(TargetDetectability, na.rm = TRUE)) %>%
      ungroup()
  } else {
    df_f <- df_f %>%
      group_by(Target) %>%
      mutate(mean_value = mean(TargetDetectability_value, na.rm = TRUE),
             mean_pct   = mean(TargetDetectability, na.rm = TRUE)) %>%
      ungroup()
  }
  
  df_f <- df_f %>%
    mutate(Target = factor(Target,
                           levels = unique(Target[order(mean_value, decreasing = TRUE)])
    ))
  
  p <- ggplot(df_f,
              aes(x = Target, y = TargetDetectability_value,
                  fill = mean_pct)) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_gradient(low = "blue", high = "orange") +
    theme_bw(base_size = 14) +
    labs(
      title = paste("Target Detectability -", fluid),
      x = "Target",
      y = "Detectability (NPQ - LOD)",
      fill = "Detectability (%)"
    ) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
  
  if (by_plate) p <- p + facet_wrap(~PlateID, ncol = 1)
  
  p
}

## 6. Get median per fluid and plate
get_median_per_fluid_plate <- function(df, fluid) {
  df %>%
    filter(SampleMatrixType == fluid) %>%
    group_by(Target, UniProtID,PlateID) %>%
    summarise(median_NPQ = median(NPQ), .groups = "drop")
}

## 7. Get median delta based on ref plate
delta_median_per_fluid <- function(df,plate_ref){
  df %>%
    group_by(Target) %>%
    mutate(ref_median = median_NPQ[PlateID == plate_ref],
           delta = ref_median - median_NPQ) %>%
    ungroup() %>%
    select(Target,PlateID,delta)
}


## 8. Get variance by fluid in respect to the covariates
variance_by_fluid = function(df, value_col = "NPQ", covariates = c("age","sex","center")) {
  
  results <- lapply(unique(df$SampleMatrixType), function(fluid) {
    df_f <- df %>% filter(SampleMatrixType == fluid)
    
    lapply(unique(df_f$Target), function(tgt) {
      sub_df <- df_f %>% filter(Target == tgt)
      
      # Ensure 'type' is factor
      sub_df$type <- as.factor(sub_df$type)
      
      # Skip if 'type' has < 2 levels
      if(nlevels(sub_df$type) < 2){
        warning("Skipping Target ", tgt, " in Fluid ", fluid, " because 'type' has < 2 levels")
        return(NULL)
      }
      
      # Keep covariates with ≥2 unique values
      valid_covs <- covariates[sapply(sub_df[covariates], function(x) n_distinct(x) > 1)]
      
      # Convert character covariates to factors if needed
      for(cv in valid_covs){
        if(is.character(sub_df[[cv]]) || is.logical(sub_df[[cv]])) sub_df[[cv]] <- as.factor(sub_df[[cv]])
        if(is.factor(sub_df[[cv]]) && nlevels(sub_df[[cv]]) < 2) valid_covs <- setdiff(valid_covs, cv)
      }
      
      # Build full formula
      fmla_full <- as.formula(
        paste(value_col, "~ type",
              if(length(valid_covs) > 0) paste("+", paste(valid_covs, collapse = " + ")))
      )
      
      fit_full <- try(lm(fmla_full, data = sub_df), silent = TRUE)
      if(inherits(fit_full, "try-error")) return(NULL)
      
      # Type-only variance
      fmla_type <- as.formula(
        paste(value_col, "~", if(length(valid_covs) > 0) paste(valid_covs, collapse = " + ") else "1")
      )
      fit_type <- try(lm(fmla_type, data = sub_df), silent = TRUE)
      if(inherits(fit_type, "try-error")) return(NULL)
      
      var_type <- sum(resid(fit_type)^2) - sum(resid(fit_full)^2)
      
      # Covariate-specific variance
      var_covs <- sapply(valid_covs, function(cov) {
        covs_reduced <- setdiff(valid_covs, cov)
        fmla_reduced <- as.formula(
          paste(value_col, "~ type",
                if(length(covs_reduced) > 0) paste("+", paste(covs_reduced, collapse = " + ")))
        )
        fit_reduced <- try(lm(fmla_reduced, data = sub_df), silent = TRUE)
        if(inherits(fit_reduced, "try-error")) return(NA_real_)
        sum(resid(fit_reduced)^2) - sum(resid(fit_full)^2)
      })
      
      tibble(
        Fluid = fluid,
        Target = tgt,
        variable = c("type", names(var_covs)),
        variance_explained = c(var_type, var_covs)
      )
      
    }) %>% bind_rows()
    
  }) %>% bind_rows()
  
  # Summarize per variable per fluid
  summary_df <- results %>%
    group_by(Fluid, variable) %>%
    summarize(
      median_variance_pct = median(variance_explained, na.rm = TRUE),
      mean_variance_pct   = mean(variance_explained, na.rm = TRUE),
      .groups = "drop"
    )
  
  list(
    raw = results,
    summary = summary_df
  )
}

## 9. Plot variance of each covariate per fluid
plot_variance_by_fluid = function(summary_df) {
  
  # Ensure variable is a factor with a consistent order
  summary_df <- summary_df %>%
    mutate(variable = factor(variable, levels = c("type", "age", "sex", "center")))
  
  ggplot(summary_df, aes(x = variable, y = median_variance_pct, fill = variable)) +
    geom_col(width = 0.6, show.legend = FALSE) +
    geom_point(aes(y = mean_variance_pct, color = variable), size = 3, show.legend = FALSE) +
    facet_wrap(~ Fluid, scales = "free_y") +
    labs(
      x = "Variable",
      y = "Variance Explained (%)",
      title = "Variance explained by each covariate",
      subtitle = "Bars = median, points = mean"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    scale_fill_brewer(palette = "Set2") +
    scale_color_brewer(palette = "Set2") # match dots to bars
}

###############################################
### 1. Assign Group Labels
###############################################

samples_ID_type <- sample_ID_info %>%
  rename(ParticipantCode = `Optional Informarion (patient ID)`) %>%
  left_join(bind_rows(ALS_ID, CTR_ID, PGMC_ID, mimic_ID), by = "ParticipantCode") %>%
  filter(ParticipantCode != "DE320")

samples_PGMC_CTR_ID_type = sample_ID_info %>%
  rename(ParticipantCode = `Optional Informarion (patient ID)`) %>%
  left_join(bind_rows(CTR_ID,PGMC_mutations_ID), by = "ParticipantCode") %>%
  filter(ParticipantCode != "DE320")

###############################################
### 2. Merge protein data with sample type
###############################################

protein_data_IDs <- protein_data %>%
  filter(SampleType == "Sample") %>%
  select(SampleName, SampleMatrixType, Target, UniProtID, ProteinName, NPQ) %>%
  left_join(samples_ID_type %>% rename(SampleName = `Sample ID`), by = "SampleName")

protein_data_CTR_PGMC_IDs <- protein_data %>%
  filter(SampleType == "Sample") %>%
  select(SampleName, SampleMatrixType, Target, UniProtID, ProteinName, NPQ) %>%
  left_join(samples_PGMC_CTR_ID_type %>% rename(SampleName = `Sample ID`), by = "SampleName")

# get count table of APOE
table_APOE = protein_data_IDs %>%
  left_join(samples_PGMC_CTR_ID_type %>% rename(subtype = type,
                                                SampleName = `Sample ID`)) %>%
  mutate(
    subtype = ifelse(subtype == "CTR" | is.na(subtype), "No subtype",
                     ifelse(subtype %in% c("FUS","UBQLN2","FIG4","other"),"others",subtype))) %>%
  filter(Target == "APOE4") %>%
  mutate(APOE_status =  ifelse(NPQ > 10, "Carrier",ifelse(NPQ <= 10, "Non-carrier",NA)))

table_APOE_notNA <- table_APOE %>% filter(!is.na(type))
table_APOE_NA <- table_APOE %>% filter(is.na(type)) %>%
  select(-PatientID,-type) %>%
  left_join(clinical_table_slim %>% select(PatientID,ParticipantCode,type))
table_APOE <- rbind(table_APOE_notNA,table_APOE_NA[,match(colnames(table_APOE_NA),colnames(table_APOE_notNA))])

carrier_table <- table_APOE %>%
  group_by(SampleMatrixType, type, subtype, APOE_status) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = APOE_status,
    values_from = n,
    values_fill = 0
  ) %>%
  mutate(
    Total = Carrier + `Non-carrier`
  )
writexl::write_xlsx(carrier_table,"results/APOE4_carriers.xlsx")

###############################################
### 3. Check effect of each covariate on the NPQ
###############################################

covariate_missing = protein_data_IDs %>%
  filter(!is.na(type)) %>%
  left_join(Sex_age_all_participants %>% dplyr::rename(PatientID = Pseudonyme)) %>%
  mutate(center = dplyr::case_when(
    grepl("TR", ParticipantCode) ~ "Turkey",
    grepl("CH", ParticipantCode) ~ "Switzerland",
    grepl("DE", ParticipantCode) ~ "Germany",
    grepl("SK", ParticipantCode) ~ "Slovakia",
    grepl("FR", ParticipantCode) ~ "France",
    grepl("IL", ParticipantCode) ~ "Israel",
    TRUE                 ~ NA_character_
  )) %>%
  group_by(Target) %>%
  summarise(
    n = n(),
    age_missing = sum(is.na(age)),
    sex_missing = sum(is.na(sex)),
    center_missing = sum(is.na(center))
  )

df_f = protein_data_IDs %>%
  filter(!is.na(type)) %>%
  left_join(Sex_age_all_participants %>% dplyr::rename(PatientID = Pseudonyme)) %>%
  mutate(center = dplyr::case_when(
    grepl("TR", ParticipantCode) ~ "Turkey",
    grepl("CH", ParticipantCode) ~ "Switzerland",
    grepl("DE", ParticipantCode) ~ "Germany",
    grepl("SK", ParticipantCode) ~ "Slovakia",
    grepl("FR", ParticipantCode) ~ "France",
    grepl("IL", ParticipantCode) ~ "Israel",
    TRUE                 ~ NA_character_
  ))

lm_full <- lm(NPQ ~ type + age + sex + center, data = df_f)
summary(lm_full) # the covariates are significantly effecting the NPQ

lm_reduced <- lm(NPQ ~ type + age + sex, data = df_f)
anova(lm_reduced, lm_full)

# check variance by fluid and each covariate
var_fluid_df <- variance_by_fluid(df_f, value_col = "NPQ", covariates = c("age","sex","center"))

writexl::write_xlsx(var_fluid_df$summary,"results/variance_by_covariate_fluid.xlsx")

pdf("plots/variance_covariate.pdf")
plot_variance_by_fluid(var_fluid_df$summary)
dev.off()

###############################################
### 5. Sample Counts Across Fluids
###############################################

fluids <- c("SERUM", "PLASMA", "CSF")

sample_counts <- bind_rows(
  lapply(fluids, function(f) count_samples_per_fluid(protein_data_IDs, f))
)

writexl::write_xlsx(sample_counts, "results/samples_biofluid_overview.xlsx")

###############################################
### 6. Correlation Between Fluids
###############################################

df_SERUM  <- get_mean_per_fluid(protein_data_IDs %>% filter(!is.na(type)), "SERUM")
df_PLASMA <- get_mean_per_fluid(protein_data_IDs %>% filter(!is.na(type)), "PLASMA")
df_CSF    <- get_mean_per_fluid(protein_data_IDs %>% filter(!is.na(type)), "CSF")

p1 = make_corr(df_PLASMA, df_SERUM,  "PLASMA", "SERUM", "plots/correlation_plasma_serum.pdf")
p2 = make_corr(df_SERUM,  df_CSF,    "SERUM",  "CSF",   "plots/correlation_serum_CSF.pdf")
p3 = make_corr(df_PLASMA, df_CSF,    "PLASMA", "CSF",   "plots/correlation_plasma_CSF.pdf")

# -> Correlation plots all together 
combined <- plot_grid(p2, p1, p3, nrow = 1, align = "h") 
pdf("plots/correlation_combined.pdf", width = 18, height = 6) 
print(combined) 
dev.off()

###############################################
### 7. Target Detectability per Fluid & Plate
###############################################

td <- protein_data  %>%
  select(PlateID, SampleMatrixType, Target, NPQ) %>%
  left_join(target_detectability %>% rename(Target = TargetName), by = c("Target","PlateID")) %>%
  mutate(TargetDetectability_value = NPQ - TargetLOD_NPQ,
         TargetDetectability =
           as.numeric(TargetDetectability) * 100) %>%
  mutate(TargetDetectability_value =
           ifelse(Target %in% c("APOE", "CRP"),
                  abs(TargetDetectability_value),
                  TargetDetectability_value))

fluids <- c("CSF", "PLASMA", "SERUM")

# Save plots
for (fluid in fluids) {
  p1 <- plot_detectability(td, fluid, by_plate = FALSE)
  p2 <- plot_detectability(td, fluid, by_plate = TRUE)
  
  Cairo::CairoPDF(paste0("plots/TargetDetectability_", fluid, ".pdf"),
                  width = 28, height = 6)
  print(p1)
  dev.off()
  
  Cairo::CairoPDF(paste0("plots/TargetDetectability_byPlate_", fluid, ".pdf"),
                  width = 28, height = 8)
  print(p2)
  dev.off()
}

###############################################
### 8. Explore Normalization
###############################################
protein_data_median_all_together = protein_data %>%
  group_by(Target, UniProtID,PlateID) %>%
  summarise(median_NPQ = median(NPQ), .groups = "drop")

df_delta_all_together = delta_median_per_fluid(protein_data_median_all_together,"P004_Plate1_090225_CB1.xml")

protein_data_norm_all_together = protein_data %>%
  left_join(df_delta_all_together) %>%
  mutate(NPQ_norm = NPQ + delta) %>%
  select(-delta)

#################################################################
### 9. NPQ Distributions by Plate original and normalized data
#################################################################
protein_long_alltogether <- protein_data_norm_all_together %>%
  mutate(
    SampleMatrixType = factor(
      SampleMatrixType,
      levels = c("PLASMA", "SERUM", "CSF", "CONTROL")
    )
  ) %>%
  select(
    Target, SampleMatrixType, PlateID, NPQ, NPQ_norm
  ) %>%
  pivot_longer(
    cols = c(NPQ, NPQ_norm),
    names_to = "Normalization",
    values_to = "value"
  ) %>%
  mutate(
    Normalization = factor(recode(
      Normalization,
      NPQ = "Before normalization",
      NPQ_norm = "After normalization"
    ),levels = c("Before normalization", "After normalization")
    )
  )

unique_targets <- unique(protein_data$Target)

for (target in unique_targets) {
  
  df_t_alltogether <- protein_long_alltogether %>%
    filter(Target == target)
  
  p_alltogether <- ggboxplot(
    df_t_alltogether,
    x = "SampleMatrixType",
    y = "value",
    color = "PlateID",
    add = "jitter"
  ) +
    facet_wrap(~ Normalization, ncol = 2) +
    labs(
      title = target,
      y = "NPQ"
    ) +
    theme_bw()
  
  pdf(paste0("plots/NPQ_fluid_plate/NPQ_", target, "_before_after_alltogether.pdf"),
      width = 10, height = 5)
  print(p_alltogether)
  dev.off()
}


#################################################################
### 10. Project-LOD of original and  normalized data
#################################################################

##############
# Project-LOD normalized data
df_reads_norm <- protein_data_norm_all_together %>%
  filter(SampleType == "NC") %>%
  mutate(
    reads = 2^NPQ_norm - 1
  )

lod_linear_norm <- df_reads_norm %>%
  group_by(Target) %>%
  summarise(
    mean_reads = mean(reads, na.rm = TRUE),
    sd_reads   = sd(reads, na.rm = TRUE),
    lod_reads  = mean_reads + 3 * sd_reads,
    .groups = "drop"
  )

lod_project_norm <- lod_linear_norm %>%
  mutate(
    LOD_NPQ_norm = log2(lod_reads + 1)
  ) %>%
  select(Target, LOD_NPQ_norm)

# Project-LOD non-normalized data
df_reads <- protein_data %>%
  filter(SampleType == "NC") %>%
  mutate(
    reads = 2^NPQ - 1
  )

lod_linear <- df_reads %>%
  group_by(Target) %>%
  summarise(
    mean_reads = mean(reads, na.rm = TRUE),
    sd_reads   = sd(reads, na.rm = TRUE),
    lod_reads  = mean_reads + 3 * sd_reads,
    .groups = "drop"
  )

lod_project <- lod_linear %>%
  mutate(
    LOD_NPQ = log2(lod_reads + 1)
  ) %>%
  select(Target, LOD_NPQ)

## Attached to original target detectability
target_detectability_extra = target_detectability %>%
  dplyr::rename(Target = TargetName) %>%
  left_join(lod_project_norm) %>%
  left_join(lod_project) %>%
  dplyr::rename(original_TargetLOD = TargetLOD_NPQ,
                ProjectLOD_norm = LOD_NPQ_norm,
                ProjectLOD = LOD_NPQ) %>%
  arrange(Target)

#################################################################
### 11. Check targets' detectability across samples
#################################################################
protein_data_with_lod <- protein_data %>%
  left_join(target_detectability_extra %>%
              select(Target,ProjectLOD) %>% 
              distinct()) %>%
  mutate(below_lod = ifelse(Target %in% c("APOE","CRP"),
                            NPQ > ProjectLOD,
                            NPQ < ProjectLOD))

detectability_summary <- protein_data_with_lod %>%
  group_by(SampleMatrixType, Target) %>%
  summarise(
    n_samples = n(),
    n_below_lod = sum(below_lod, na.rm = TRUE),
    frac_below_lod = n_below_lod / n_samples,
    .groups = "drop"
  ) %>%
  filter(SampleMatrixType != "CONTROL") %>%
  arrange(Target,frac_below_lod) %>%
  mutate(percent_below_lod = frac_below_lod * 100,
         detectability = ifelse(frac_below_lod>0.5,"low","high"))

writexl::write_xlsx(detectability_summary,"results/detectability_summary.xlsx")

targets_below_LOD = detectability_summary %>%
  select(SampleMatrixType, Target, percent_below_lod) %>%
  tidyr::pivot_wider(
    names_from = SampleMatrixType,
    values_from = percent_below_lod
  )

writexl::write_xlsx(targets_below_LOD,"results/targets_below_LOD.xlsx")